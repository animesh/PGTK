#!/usr/bin/env python3
"""Deep, read-only audit of a completed PGTK results tree.

Required environment:
  PGTK_RESULTS      Completed results directory.
  PGTK_SOURCE       PGTK source checkout used to generate the results.

Optional environment:
  PGTK_AUDIT_OUT    Output directory. Default: ./PGTK-deep-audit
  PGTK_JOB_ID       Job ID. Auto-detected from pipeline_trace-*.tsv when unique.
  PGTK_MAX_EVENTS   Maximum stratified events for BAM/display audit. Default: 1200.
  PGTK_TARGET_IDS   Comma-separated finding IDs always included.
  PGTK_STRICT       true to exit nonzero on ERROR findings. Default: false.

Run this script inside the exact pinned Pysam container used by PGTK.
The script never modifies the results tree.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import html
import json
import os
import re
import shutil
import statistics
import sys
import tarfile
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

try:
    import pysam
except Exception as exc:
    raise SystemExit(f"ERROR: pysam is required in the executing container: {exc}")

RESULTS = Path(os.environ.get("PGTK_RESULTS", "")).expanduser().resolve()
SOURCE = Path(os.environ.get("PGTK_SOURCE", "")).expanduser().resolve()
OUT = Path(os.environ.get("PGTK_AUDIT_OUT", "PGTK-deep-audit")).expanduser().resolve()
MAX_EVENTS = int(os.environ.get("PGTK_MAX_EVENTS", "1200"))
TARGET_IDS = {
    value.strip()
    for value in os.environ.get(
        "PGTK_TARGET_IDS", "ALDH16A1_TK14_19_49464415_A_G_A_G"
    ).split(",")
    if value.strip()
}
STRICT = os.environ.get("PGTK_STRICT", "false").lower() in {"1", "true", "yes"}

if not RESULTS.is_dir():
    raise SystemExit(f"ERROR: PGTK_RESULTS is not a directory: {RESULTS}")
if not SOURCE.is_dir():
    raise SystemExit(f"ERROR: PGTK_SOURCE is not a directory: {SOURCE}")
if OUT == RESULTS or RESULTS in OUT.parents:
    raise SystemExit("ERROR: PGTK_AUDIT_OUT must be outside the results tree")

if OUT.exists():
    shutil.rmtree(OUT)
OUT.mkdir(parents=True)

STARTED = time.time()
ISSUES: list[dict[str, Any]] = []
METRICS: dict[str, Any] = {}


def issue(level: str, code: str, message: str, path: str = "", details: Any = "") -> None:
    ISSUES.append(
        {
            "level": level,
            "code": code,
            "message": message,
            "path": path,
            "details": details if isinstance(details, str) else json.dumps(details, sort_keys=True),
        }
    )


def write_tsv(path: Path, rows: Iterable[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fields})


def first_value(row: dict[str, Any], aliases: Iterable[str], default: Any = "") -> Any:
    lower = {str(key).lower(): value for key, value in row.items()}
    for alias in aliases:
        value = lower.get(alias.lower())
        if value not in (None, ""):
            return value
    return default


def as_int(value: Any, default: int = 0) -> int:
    try:
        return int(float(str(value).replace(",", "")))
    except Exception:
        return default


def find_index(bam: Path) -> Path | None:
    for candidate in (Path(str(bam) + ".bai"), bam.with_suffix(".bai"), Path(str(bam) + ".csi")):
        if candidate.is_file() and candidate.stat().st_size > 0:
            return candidate
    return None


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(RESULTS))
    except ValueError:
        return str(path)


def cigar_event_state(read: pysam.AlignedSegment, pos0: int) -> str:
    """Classify one 0-based reference coordinate against one alignment record."""
    if read.is_unmapped or read.reference_start is None:
        return "unmapped"
    ref_pos = read.reference_start
    for operation, length in read.cigartuples or []:
        if operation in (0, 7, 8):  # M, =, X
            if ref_pos <= pos0 < ref_pos + length:
                return "aligned_base"
            ref_pos += length
        elif operation == 2:  # D
            if ref_pos <= pos0 < ref_pos + length:
                return "deletion"
            ref_pos += length
        elif operation == 3:  # N
            if ref_pos <= pos0 < ref_pos + length:
                return "skipped_region"
            ref_pos += length
        elif operation == 1:  # I, anchored after previous reference base
            if pos0 == ref_pos - 1:
                return "insertion_anchor"
        elif operation in (4, 5, 6):
            continue
    return "outside_aligned_blocks"


def alignment_key(read: pysam.AlignedSegment) -> tuple[Any, ...]:
    try:
        hi = read.get_tag("HI")
    except Exception:
        hi = None
    return (
        read.query_name,
        read.reference_id,
        read.reference_start,
        read.flag,
        read.cigarstring,
        hi,
    )


# 1. Source identity and required files
source_rows = []
manifest = SOURCE / "pipeline_required_files.txt"
required = []
if manifest.is_file():
    required = [line.strip() for line in manifest.read_text().splitlines() if line.strip() and not line.startswith("#")]
else:
    issue("ERROR", "SOURCE_MANIFEST_MISSING", "pipeline_required_files.txt is absent", str(manifest))

for name in required:
    path = SOURCE / name
    exists = path.is_file()
    size = path.stat().st_size if exists else 0
    source_rows.append({"path": name, "exists": exists, "bytes": size})
    if not exists or size == 0:
        issue("ERROR", "SOURCE_FILE_MISSING", f"Required source file missing or empty: {name}", str(path))
write_tsv(OUT / "source_manifest_audit.tsv", source_rows, ["path", "exists", "bytes"])

# 2. Result inventory
inventory_rows = []
family_counts: Counter[str] = Counter()
for path in sorted(RESULTS.rglob("*")):
    if not path.is_file():
        continue
    relative = rel(path)
    family_counts[Path(relative).parts[0]] += 1
    inventory_rows.append(
        {
            "path": relative,
            "bytes": path.stat().st_size,
            "mtime_ns": path.stat().st_mtime_ns,
            "suffix": "".join(path.suffixes),
        }
    )
write_tsv(OUT / "results_inventory.tsv", inventory_rows, ["path", "bytes", "mtime_ns", "suffix"])
write_tsv(
    OUT / "result_family_counts.tsv",
    [{"family": key, "files": value} for key, value in sorted(family_counts.items())],
    ["family", "files"],
)
METRICS["result_files"] = len(inventory_rows)

# 3. Job and trace audit
trace_candidates = sorted(RESULTS.glob("pipeline_trace-*.tsv"), key=lambda p: p.stat().st_mtime_ns)
job_id = os.environ.get("PGTK_JOB_ID", "").strip()
if job_id:
    trace = RESULTS / f"pipeline_trace-{job_id}.tsv"
else:
    trace = trace_candidates[-1] if trace_candidates else Path()
    if trace:
        match = re.search(r"pipeline_trace-(.+)\.tsv$", trace.name)
        job_id = match.group(1) if match else "unknown"

trace_summary: Counter[str] = Counter()
trace_rows: list[dict[str, Any]] = []
if trace and trace.is_file():
    with trace.open(encoding="utf-8") as handle:
        trace_rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in trace_rows:
        status = first_value(row, ["status"])
        trace_summary[status] += 1
        if status not in {"COMPLETED", "CACHED"}:
            issue("ERROR", "TRACE_BAD_STATUS", f"Trace contains status {status}", rel(trace), row)
        if as_int(first_value(row, ["attempt"], 1), 1) > 1:
            issue("WARN", "TRACE_RETRY", "Task attempt above 1", rel(trace), row)
    METRICS["trace_rows"] = len(trace_rows)
else:
    issue("ERROR", "TRACE_MISSING", "No pipeline trace found", str(RESULTS))
write_tsv(
    OUT / "trace_status_counts.tsv",
    [{"status": key, "count": value} for key, value in sorted(trace_summary.items())],
    ["status", "count"],
)

# 4. Failure ledger
ledger = RESULTS / "failure_logs" / job_id / "failure_ledger.tsv"
if ledger.is_file():
    with ledger.open(encoding="utf-8") as handle:
        ledger_rows = list(csv.DictReader(handle, delimiter="\t"))
    if ledger_rows:
        issue("ERROR", "FAILURE_LEDGER_NONEMPTY", f"Failure ledger has {len(ledger_rows)} rows", rel(ledger))
else:
    issue("WARN", "FAILURE_LEDGER_MISSING", "Job-specific failure ledger is absent", rel(ledger))

# 5. Samplesheet and expected samples
samplesheet_candidates = [SOURCE / "samples.csv", RESULTS.parent / "samples.csv"]
samples: list[str] = []
for candidate in samplesheet_candidates:
    if candidate.is_file():
        with candidate.open(encoding="utf-8-sig") as handle:
            sample_rows = list(csv.DictReader(handle))
        samples = [str(first_value(row, ["sample"])) for row in sample_rows if first_value(row, ["sample"])]
        break
if not samples:
    issue("WARN", "SAMPLESHEET_UNAVAILABLE", "Could not determine expected samples")
METRICS["samples"] = samples

# 6. VCF structure and counts
vcf_rows = []
for path in sorted(RESULTS.rglob("*.vcf.gz")):
    records = pass_records = 0
    header_samples: list[str] = []
    error = ""
    try:
        with pysam.VariantFile(str(path)) as vcf:
            header_samples = list(vcf.header.samples)
            for record in vcf:
                records += 1
                if not record.filter.keys() or "PASS" in record.filter.keys():
                    pass_records += 1
    except Exception as exc:
        error = str(exc).replace("\t", " ").replace("\n", " ")
        issue("ERROR", "VCF_READ_ERROR", error, rel(path))
    index = Path(str(path) + ".tbi")
    if not index.is_file():
        index = Path(str(path) + ".csi")
    indexed = index.is_file() and index.stat().st_size > 0
    if not indexed:
        issue("ERROR", "VCF_INDEX_MISSING", "Compressed VCF index missing", rel(path))
    vcf_rows.append(
        {
            "path": rel(path),
            "bytes": path.stat().st_size,
            "records": records,
            "pass_records": pass_records,
            "samples": ",".join(header_samples),
            "index": rel(index) if indexed else "",
            "error": error,
        }
    )
write_tsv(OUT / "vcf_audit.tsv", vcf_rows, ["path", "bytes", "records", "pass_records", "samples", "index", "error"])
METRICS["vcfs"] = len(vcf_rows)

# 7. Global BAM integrity and CIGAR audit
bam_rows = []
display_bams: dict[tuple[str, str], Path] = {}
for bam in sorted(RESULTS.rglob("*.bam")):
    index = find_index(bam)
    quickcheck = "FAIL"
    index_readable = False
    total = primary = secondary = supplementary = duplicate = 0
    cigar_n = cigar_d = cigar_i = malformed = missing_seq = 0
    unique_keys: set[tuple[Any, ...]] = set()
    duplicate_keys = 0
    error = ""
    try:
        if index is None:
            raise RuntimeError("missing index")
        check = pysam.quickcheck(str(bam))
        if check:
            raise RuntimeError("quickcheck: " + "; ".join(check))
        quickcheck = "PASS"
        with pysam.AlignmentFile(str(bam), "rb", require_index=True) as handle:
            handle.get_index_statistics()
            index_readable = True
        with pysam.AlignmentFile(str(bam), "rb", check_sq=False) as handle:
            for read in handle.fetch(until_eof=True):
                total += 1
                secondary += int(read.is_secondary)
                supplementary += int(read.is_supplementary)
                duplicate += int(read.is_duplicate)
                primary += int(not read.is_secondary and not read.is_supplementary)
                if read.query_sequence is None:
                    missing_seq += 1
                key = alignment_key(read)
                if key in unique_keys:
                    duplicate_keys += 1
                else:
                    unique_keys.add(key)
                q_consumed = 0
                for operation, length in read.cigartuples or []:
                    cigar_n += int(operation == 3)
                    cigar_d += int(operation == 2)
                    cigar_i += int(operation == 1)
                    if operation in (0, 1, 4, 7, 8):
                        q_consumed += length
                if read.query_length is not None and q_consumed != read.query_length:
                    malformed += 1
        if malformed:
            issue("ERROR", "BAM_CIGAR_QUERY_LENGTH", f"{malformed} malformed CIGAR/query-length rows", rel(bam))
        if duplicate_keys:
            issue("WARN", "BAM_DUPLICATE_ALIGNMENT_KEYS", f"{duplicate_keys} exact duplicate alignment keys", rel(bam))
    except Exception as exc:
        error = str(exc).replace("\t", " ").replace("\n", " ")
        issue("ERROR", "BAM_READ_ERROR", error, rel(bam))

    match = re.match(r"(.+)\.(event_display|exact_alt_display|reference_display)\.bam$", bam.name)
    if match and "finding_reviews" in bam.parts:
        display_bams[(match.group(1), match.group(2))] = bam

    bam_rows.append(
        {
            "path": rel(bam),
            "bytes": bam.stat().st_size,
            "index": rel(index) if index else "",
            "index_not_older": bool(index and index.stat().st_mtime_ns >= bam.stat().st_mtime_ns),
            "quickcheck": quickcheck,
            "index_readable": index_readable,
            "alignments": total,
            "primary": primary,
            "secondary": secondary,
            "supplementary": supplementary,
            "duplicates": duplicate,
            "cigar_N_ops": cigar_n,
            "cigar_D_ops": cigar_d,
            "cigar_I_ops": cigar_i,
            "missing_sequence": missing_seq,
            "malformed_query_length": malformed,
            "duplicate_alignment_keys": duplicate_keys,
            "error": error,
        }
    )
write_tsv(
    OUT / "bam_audit.tsv",
    bam_rows,
    [
        "path", "bytes", "index", "index_not_older", "quickcheck", "index_readable",
        "alignments", "primary", "secondary", "supplementary", "duplicates",
        "cigar_N_ops", "cigar_D_ops", "cigar_I_ops", "missing_sequence",
        "malformed_query_length", "duplicate_alignment_keys", "error",
    ],
)
METRICS["bams"] = len(bam_rows)

# 8. Load published findings
partition_candidates = sorted((RESULTS / "igv" / "findings" / "finding_explorer" / "partitions").glob("all.jsonl.gz"))
findings: list[dict[str, Any]] = []
if partition_candidates:
    partition = partition_candidates[0]
    with gzip.open(partition, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            try:
                row = json.loads(line)
                if isinstance(row, dict):
                    findings.append(row)
            except Exception as exc:
                issue("ERROR", "FINDING_JSON_ERROR", f"Line {line_number}: {exc}", rel(partition))
else:
    issue("ERROR", "FINDINGS_PARTITION_MISSING", "all.jsonl.gz is absent")
METRICS["findings"] = len(findings)

# 9. Finding arithmetic invariants and stratification
finding_summary: Counter[tuple[str, str, str]] = Counter()
parsed_findings: list[dict[str, Any]] = []
for row in findings:
    finding_id = str(first_value(row, ["EventID", "event_id", "FindingID", "finding_id", "unique_id", "id"]))
    sample = str(first_value(row, ["Sample", "sample"]))
    chromosome = str(first_value(row, ["Chrom", "chrom", "chromosome", "Chr"]))
    position = as_int(first_value(row, ["Start", "start", "Pos", "pos", "position"]))
    candidate_class = str(first_value(row, ["CandidateClass", "candidate_class", "class", "FindingClass"]))
    status = str(first_value(row, ["ReadValidationStatus", "read_validation_status", "status"]))
    impact = str(first_value(row, ["PredictedImpact", "predicted_impact", "impact"]))
    alt = as_int(first_value(row, ["ExactAltReads", "exact_alt_reads", "AltReads", "alt_reads"]))
    ref = as_int(first_value(row, ["CleanReferenceReads", "clean_reference_reads", "ReferenceReads", "ref_reads"]))
    excluded = as_int(first_value(row, ["ExcludedReads", "excluded_reads"]))
    callable_count = as_int(first_value(row, ["CallableAlignments", "callable_alignments", "CallableReads"]))
    unique_count = as_int(first_value(row, ["UniqueAlignments", "unique_alignments", "PrimaryAlignments"]))

    if callable_count != alt + ref:
        issue("ERROR", "FINDING_CALLABLE_INVARIANT", "CallableAlignments != ExactAltReads + CleanReferenceReads", finding_id, row)
    if unique_count and unique_count != callable_count + excluded:
        issue("ERROR", "FINDING_UNIQUE_INVARIANT", "UniqueAlignments != CallableAlignments + ExcludedReads", finding_id, row)
    if callable_count == 0:
        fraction = first_value(row, ["AltFractionAmongClean", "alt_fraction_among_clean", "AltFraction"])
        if str(fraction).strip().lower() not in {"", "na", "nan", "none", "null", "."}:
            issue("ERROR", "ZERO_DENOMINATOR_FRACTION", "Callable depth is zero but ALT fraction is numeric", finding_id, row)

    finding_summary[(candidate_class, status, impact)] += 1
    parsed_findings.append(
        {
            "id": finding_id,
            "sample": sample,
            "chrom": chromosome,
            "pos": position,
            "class": candidate_class,
            "status": status,
            "impact": impact,
            "alt": alt,
            "ref": ref,
            "excluded": excluded,
            "callable": callable_count,
            "unique": unique_count,
            "raw": row,
        }
    )
write_tsv(
    OUT / "finding_strata.tsv",
    [
        {"candidate_class": key[0], "status": key[1], "impact": key[2], "findings": value}
        for key, value in sorted(finding_summary.items())
    ],
    ["candidate_class", "status", "impact", "findings"],
)

# 10. Deterministic stratified event selection
selected: list[dict[str, Any]] = []
selected_ids: set[str] = set()
for row in parsed_findings:
    if row["id"] in TARGET_IDS:
        selected.append(row)
        selected_ids.add(row["id"])

strata: dict[tuple[str, str, str, str], list[dict[str, Any]]] = defaultdict(list)
for row in parsed_findings:
    strata[(row["sample"], row["class"], row["status"], row["impact"])].append(row)

per_stratum = max(1, MAX_EVENTS // max(1, len(strata)))
for key in sorted(strata):
    candidates = sorted(strata[key], key=lambda row: (row["chrom"], row["pos"], row["id"]))
    if len(candidates) <= per_stratum:
        picks = candidates
    else:
        picks = [candidates[round(i * (len(candidates) - 1) / (per_stratum - 1))] for i in range(per_stratum)] if per_stratum > 1 else [candidates[len(candidates) // 2]]
    for row in picks:
        if row["id"] not in selected_ids and len(selected) < MAX_EVENTS:
            selected.append(row)
            selected_ids.add(row["id"])

# 11. Event-level display BAM semantics
# This catches splice-span-only contamination, same-name expansion, wrong sample tracks,
# mismatches between finding counts and display-read counts, and malformed event BAMs.
event_rows = []
for number, finding in enumerate(selected, 1):
    sample = finding["sample"]
    chrom = finding["chrom"]
    pos1 = finding["pos"]
    if not sample or not chrom or pos1 <= 0:
        issue("WARN", "EVENT_COORDINATE_UNUSABLE", "Finding lacks sample/chromosome/position", finding["id"])
        continue
    pos0 = pos1 - 1
    for category in ("event_display", "exact_alt_display", "reference_display"):
        bam = display_bams.get((sample, category))
        if bam is None:
            issue("ERROR", "DISPLAY_BAM_MISSING", f"Missing {sample}.{category}.bam", finding["id"])
            continue
        state_counts: Counter[str] = Counter()
        query_names: Counter[str] = Counter()
        alignment_keys: set[tuple[Any, ...]] = set()
        exact_duplicate_keys = 0
        records = 0
        try:
            with pysam.AlignmentFile(str(bam), "rb", require_index=True) as handle:
                for read in handle.fetch(chrom, max(0, pos0), pos0 + 1):
                    records += 1
                    state = cigar_event_state(read, pos0)
                    state_counts[state] += 1
                    query_names[read.query_name] += 1
                    key = alignment_key(read)
                    if key in alignment_keys:
                        exact_duplicate_keys += 1
                    else:
                        alignment_keys.add(key)
        except Exception as exc:
            issue("ERROR", "DISPLAY_QUERY_ERROR", str(exc), f"{finding['id']}:{category}")
            continue

        observable = state_counts["aligned_base"] + state_counts["deletion"] + state_counts["insertion_anchor"]
        skip_only = state_counts["skipped_region"]
        duplicated_names = sum(1 for value in query_names.values() if value > 1)

        if skip_only:
            issue(
                "ERROR",
                "DISPLAY_SPLICE_SPAN_CONTAMINATION",
                f"{category} includes {skip_only} records whose CIGAR N skips the event coordinate",
                finding["id"],
                {"sample": sample, "chrom": chrom, "pos": pos1, "records": records, "states": dict(state_counts)},
            )
        if exact_duplicate_keys:
            issue("ERROR", "DISPLAY_DUPLICATE_ALIGNMENT", f"{exact_duplicate_keys} duplicate alignment identities", finding["id"])
        if category == "exact_alt_display" and observable < finding["alt"]:
            issue("ERROR", "EXACT_ALT_DISPLAY_UNDERCOUNT", f"Observable display records {observable} < reported exact ALT {finding['alt']}", finding["id"])
        if category == "exact_alt_display" and observable > finding["alt"]:
            issue("WARN", "EXACT_ALT_DISPLAY_OVERCOUNT", f"Observable display records {observable} > reported exact ALT {finding['alt']}", finding["id"])
        if category == "reference_display" and observable < finding["ref"]:
            issue("ERROR", "REFERENCE_DISPLAY_UNDERCOUNT", f"Observable display records {observable} < reported clean REF {finding['ref']}", finding["id"])

        event_rows.append(
            {
                "finding_id": finding["id"],
                "sample": sample,
                "class": finding["class"],
                "status": finding["status"],
                "impact": finding["impact"],
                "chrom": chrom,
                "pos": pos1,
                "category": category,
                "reported_alt": finding["alt"],
                "reported_ref": finding["ref"],
                "records_returned": records,
                "observable_at_event": observable,
                "aligned_base": state_counts["aligned_base"],
                "deletion": state_counts["deletion"],
                "insertion_anchor": state_counts["insertion_anchor"],
                "skipped_region": skip_only,
                "outside_blocks": state_counts["outside_aligned_blocks"],
                "duplicated_query_names": duplicated_names,
                "duplicate_alignment_keys": exact_duplicate_keys,
            }
        )

write_tsv(
    OUT / "display_event_audit.tsv",
    event_rows,
    [
        "finding_id", "sample", "class", "status", "impact", "chrom", "pos", "category",
        "reported_alt", "reported_ref", "records_returned", "observable_at_event", "aligned_base",
        "deletion", "insertion_anchor", "skipped_region", "outside_blocks",
        "duplicated_query_names", "duplicate_alignment_keys",
    ],
)
METRICS["display_events_selected"] = len(selected)
METRICS["display_event_track_queries"] = len(event_rows)

# 12. HTML/report/link checks
html_rows = []
for path in sorted(RESULTS.rglob("*.html")):
    text = path.read_text(encoding="utf-8", errors="replace")
    lower = text.lower()
    obvious_error = any(token in lower for token in ("traceback (most recent call last)", "internal server error", "cannot read properties of undefined"))
    if obvious_error:
        issue("ERROR", "HTML_CONTAINS_ERROR", "HTML contains an explicit error signature", rel(path))
    html_rows.append(
        {
            "path": rel(path),
            "bytes": path.stat().st_size,
            "has_igv": "igv" in lower,
            "has_traceback": "traceback (most recent call last)" in lower,
            "has_js_undefined_error": "cannot read properties of undefined" in lower,
        }
    )
write_tsv(OUT / "html_audit.tsv", html_rows, ["path", "bytes", "has_igv", "has_traceback", "has_js_undefined_error"])

# 13. FASTA structure
fasta_rows = []
for path in sorted(RESULTS.rglob("*.fa")) + sorted(RESULTS.rglob("*.fasta")) + sorted(RESULTS.rglob("*.fa.gz")) + sorted(RESULTS.rglob("*.fasta.gz")):
    opener = gzip.open if path.suffix == ".gz" else open
    sequences = duplicate_headers = invalid_sequences = 0
    headers: set[str] = set()
    current_sequence = []
    try:
        with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if current_sequence and not re.fullmatch(r"[A-Za-z*.-]+", "".join(current_sequence)):
                        invalid_sequences += 1
                    current_sequence = []
                    header = line[1:].split()[0]
                    duplicate_headers += int(header in headers)
                    headers.add(header)
                    sequences += 1
                else:
                    current_sequence.append(line)
        if current_sequence and not re.fullmatch(r"[A-Za-z*.-]+", "".join(current_sequence)):
            invalid_sequences += 1
    except Exception as exc:
        issue("ERROR", "FASTA_READ_ERROR", str(exc), rel(path))
    if duplicate_headers:
        issue("WARN", "FASTA_DUPLICATE_HEADERS", f"{duplicate_headers} duplicate FASTA headers", rel(path))
    fasta_rows.append({"path": rel(path), "sequences": sequences, "duplicate_headers": duplicate_headers, "invalid_sequences": invalid_sequences})
write_tsv(OUT / "fasta_audit.tsv", fasta_rows, ["path", "sequences", "duplicate_headers", "invalid_sequences"])

# 14. Final reports
level_counts = Counter(entry["level"] for entry in ISSUES)
code_counts = Counter(entry["code"] for entry in ISSUES)
METRICS["issues"] = dict(level_counts)
METRICS["issue_codes"] = dict(code_counts)
METRICS["job_id"] = job_id
METRICS["elapsed_seconds"] = round(time.time() - STARTED, 3)
METRICS["pysam_version"] = pysam.__version__
METRICS["samtools_version"] = getattr(pysam, "__samtools_version__", "unknown")

write_tsv(OUT / "issues.tsv", ISSUES, ["level", "code", "message", "path", "details"])
(OUT / "summary.json").write_text(json.dumps(METRICS, indent=2, sort_keys=True) + "\n", encoding="utf-8")

markdown = [
    "# PGTK deep audit",
    "",
    f"Results: `{RESULTS}`",
    f"Source: `{SOURCE}`",
    f"Job ID: `{job_id}`",
    f"Pysam: `{METRICS['pysam_version']}`",
    f"Samtools through Pysam: `{METRICS['samtools_version']}`",
    "",
    "## Overall status",
    "",
    f"ERROR: {level_counts.get('ERROR', 0)}  ",
    f"WARN: {level_counts.get('WARN', 0)}  ",
    f"INFO: {level_counts.get('INFO', 0)}",
    "",
    "## Scope",
    "",
    f"Result files: {METRICS.get('result_files', 0)}  ",
    f"Trace rows: {METRICS.get('trace_rows', 0)}  ",
    f"VCFs: {METRICS.get('vcfs', 0)}  ",
    f"BAMs: {METRICS.get('bams', 0)}  ",
    f"Published findings: {METRICS.get('findings', 0)}  ",
    f"Stratified findings selected for display-BAM audit: {METRICS.get('display_events_selected', 0)}  ",
    f"Display track queries: {METRICS.get('display_event_track_queries', 0)}",
    "",
    "## Issue counts",
    "",
]
if code_counts:
    for code, count in code_counts.most_common():
        markdown.append(f"- `{code}`: {count}")
else:
    markdown.append("No issues detected.")
markdown.extend(
    [
        "",
        "## Important interpretation",
        "",
        "This audit validates file structure, tabular invariants, BAM/VCF readability, CIGAR consistency, and a deterministic stratified sample of real finding-to-display-BAM relationships. It does not execute browser JavaScript. Any display-BAM contamination detected here is sufficient to reject IGV report readiness.",
        "",
        "See `issues.tsv`, `display_event_audit.tsv`, `bam_audit.tsv`, and `vcf_audit.tsv` for details.",
    ]
)
(OUT / "REPORT.md").write_text("\n".join(markdown) + "\n", encoding="utf-8")

# Checksums and compact archive
checksum_path = OUT / "checksums.sha256"
with checksum_path.open("w", encoding="utf-8") as handle:
    for path in sorted(OUT.rglob("*")):
        if path.is_file() and path != checksum_path:
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            handle.write(f"{digest}  {path.relative_to(OUT)}\n")

archive = OUT.parent / f"{OUT.name}.tar.gz"
if archive.exists():
    archive.unlink()
with tarfile.open(archive, "w:gz", compresslevel=9) as handle:
    handle.add(OUT, arcname=OUT.name)
archive_hash = hashlib.sha256(archive.read_bytes()).hexdigest()
(archive.parent / f"{archive.name}.sha256").write_text(f"{archive_hash}  {archive.name}\n", encoding="utf-8")

print(f"PGTK deep audit complete")
print(f"Results: {RESULTS}")
print(f"Findings: {METRICS.get('findings', 0)}")
print(f"Display events selected: {METRICS.get('display_events_selected', 0)}")
print(f"ERROR: {level_counts.get('ERROR', 0)}")
print(f"WARN: {level_counts.get('WARN', 0)}")
print(f"Report: {OUT / 'REPORT.md'}")
print(f"Archive: {archive}")
print(f"SHA256: {archive_hash}")

if STRICT and level_counts.get("ERROR", 0):
    raise SystemExit(2)
