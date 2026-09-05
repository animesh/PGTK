#!/usr/bin/env python3
"""PGTK local result auditor.

Runs read-only checks against a completed PGTK project/results directory and creates
an evidence-only shareable audit bundle. Uses only the Python standard library;
optional external tools improve BAM/CRAM and shell validation.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable


@dataclass
class Check:
    category: str
    name: str
    status: str
    detail: str
    path: str = ""


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            h.update(chunk)
    return h.hexdigest()


def human_size(value: int) -> str:
    units = ["B", "KiB", "MiB", "GiB", "TiB"]
    n = float(value)
    for unit in units:
        if n < 1024 or unit == units[-1]:
            return f"{n:.1f} {unit}"
        n /= 1024
    return f"{value} B"


def run_command(argv: list[str], cwd: Path | None = None, timeout: int = 300) -> tuple[int, str]:
    try:
        completed = subprocess.run(
            argv, cwd=str(cwd) if cwd else None, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            timeout=timeout, check=False,
        )
        return completed.returncode, completed.stdout[-20000:]
    except (OSError, subprocess.TimeoutExpired) as exc:
        return 125, str(exc)


def iter_files(root: Path) -> Iterable[Path]:
    for path in sorted(root.rglob("*")):
        if path.is_file() and not path.is_symlink():
            yield path


def rel(path: Path, root: Path) -> str:
    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path)


def add(checks: list[Check], category: str, name: str, status: str, detail: str, path: str = "") -> None:
    checks.append(Check(category, name, status, detail.replace("\t", " ").replace("\n", " ")[:4000], path))


def read_text_sample(path: Path, limit: int = 2_000_000) -> str:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        return handle.read(limit)


def check_gzip(path: Path) -> tuple[bool, str]:
    try:
        with gzip.open(path, "rb") as handle:
            while handle.read(8 * 1024 * 1024):
                pass
        return True, "gzip stream readable to EOF"
    except Exception as exc:
        return False, str(exc)


def check_tabular(path: Path) -> tuple[str, dict]:
    opener = gzip.open if path.suffix == ".gz" else open
    delimiter = "," if path.name.lower().endswith((".csv", ".csv.gz")) else "\t"
    rows = 0
    malformed = 0
    duplicate_header_names: list[str] = []
    header: list[str] = []
    expected = None
    try:
        with opener(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.reader(handle, delimiter=delimiter)
            for row in reader:
                if not row or (len(row) == 1 and row[0].startswith("#")):
                    continue
                if expected is None:
                    header = row
                    expected = len(row)
                    duplicate_header_names = sorted(k for k, v in Counter(header).items() if k and v > 1)
                    continue
                rows += 1
                if len(row) != expected:
                    malformed += 1
        status = "PASS" if malformed == 0 and not duplicate_header_names else "FAIL"
        return status, {
            "columns": expected or 0, "rows": rows, "malformed_rows": malformed,
            "duplicate_header_names": duplicate_header_names,
        }
    except Exception as exc:
        return "FAIL", {"error": str(exc), "rows": rows}


def parse_trace(path: Path) -> dict:
    with path.open("rt", encoding="utf-8", errors="replace", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    status_counts = Counter((row.get("status") or "UNKNOWN").upper() for row in rows)
    processes = Counter(row.get("process") or "UNKNOWN" for row in rows)
    failed = [row for row in rows if (row.get("status") or "").upper() not in {"COMPLETED", "CACHED"}]
    duplicate_keys = Counter((row.get("task_id"), row.get("hash"), row.get("process"), row.get("tag")) for row in rows)
    duplicates = sum(1 for count in duplicate_keys.values() if count > 1)
    return {
        "rows": len(rows), "status_counts": dict(status_counts),
        "processes": len(processes), "failed_or_other": len(failed),
        "duplicate_task_keys": duplicates,
    }


def parse_checksum_file(checksum_file: Path) -> list[tuple[str, str]]:
    entries = []
    for line in checksum_file.read_text(encoding="utf-8", errors="replace").splitlines():
        match = re.match(r"^([0-9a-fA-F]{64})\s+[* ]?(.+?)\s*$", line)
        if match:
            entries.append((match.group(1).lower(), match.group(2)))
    return entries


def safe_tar_members(path: Path) -> tuple[bool, str, int]:
    try:
        with tarfile.open(path, "r:*") as archive:
            members = archive.getmembers()
            unsafe = [m.name for m in members if Path(m.name).is_absolute() or ".." in Path(m.name).parts]
            if unsafe:
                return False, f"unsafe member paths: {unsafe[:5]}", len(members)
            for member in members:
                if member.isfile():
                    extracted = archive.extractfile(member)
                    if extracted:
                        while extracted.read(8 * 1024 * 1024):
                            pass
            return True, "archive readable; member paths safe", len(members)
    except Exception as exc:
        return False, str(exc), 0


def audit(project: Path, results: Path, job_id: str, output_root: Path, hash_large: bool,
          max_hash_gib: float, run_source_checks: bool) -> tuple[Path, dict]:
    started = time.time()
    checks: list[Check] = []
    project = project.resolve()
    results = results.resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    report_dir = output_root / f"PGTK-independent-audit-{job_id}"
    if report_dir.exists():
        shutil.rmtree(report_dir)
    report_dir.mkdir()

    add(checks, "layout", "project directory", "PASS" if project.is_dir() else "FAIL", str(project), str(project))
    add(checks, "layout", "results directory", "PASS" if results.is_dir() else "FAIL", str(results), str(results))
    if not project.is_dir() or not results.is_dir():
        raise SystemExit("Project or results directory does not exist")

    files = list(iter_files(results))
    add(checks, "inventory", "result files discovered", "PASS" if files else "FAIL", str(len(files)), str(results))

    inventory_path = report_dir / "results_inventory.tsv"
    manifest_path = report_dir / "results_checksums.sha256"
    extension_counts: Counter[str] = Counter()
    total_bytes = 0
    hash_limit = int(max_hash_gib * 1024**3)
    with inventory_path.open("w", encoding="utf-8", newline="") as inv, manifest_path.open("w", encoding="utf-8") as man:
        writer = csv.writer(inv, delimiter="\t", lineterminator="\n")
        writer.writerow(["relative_path", "size_bytes", "mtime_epoch", "sha256", "hash_status"])
        for path in files:
            size = path.stat().st_size
            total_bytes += size
            suffix = "".join(path.suffixes[-2:]).lower() or "[none]"
            extension_counts[suffix] += 1
            digest = ""
            hash_status = "SKIPPED_SIZE_LIMIT"
            if hash_large or size <= hash_limit:
                try:
                    digest = sha256_file(path)
                    hash_status = "HASHED"
                    man.write(f"{digest}  {rel(path, results)}\n")
                except OSError as exc:
                    hash_status = f"ERROR:{exc}"
            writer.writerow([rel(path, results), size, int(path.stat().st_mtime), digest, hash_status])
    add(checks, "inventory", "result bytes inventoried", "PASS", f"{total_bytes} ({human_size(total_bytes)})")

    empty_files = [path for path in files if path.stat().st_size == 0]
    allowed_empty_patterns = ("failure_ledger.tsv", ".gitkeep")
    allowed_empty_picard = {
        "picard_MarkIlluminaAdapters_histogram.txt",
        "picard_MeanQualityByCycle_histogram.txt",
        "picard_MeanQualityByCycle_histogram_1.txt",
        "picard_QualityScoreDistribution_histogram.txt",
    }
    suspicious_empty = [
        path for path in empty_files
        if not path.name.endswith(allowed_empty_patterns)
        and not ("qc/multiqc_report_data/" in rel(path, results) and path.name in allowed_empty_picard)
    ]
    add(checks, "inventory", "unexpected empty files", "PASS" if not suspicious_empty else "WARN",
        "none" if not suspicious_empty else ", ".join(rel(p, results) for p in suspicious_empty[:100]))

    broken_links = [p for p in results.rglob("*") if p.is_symlink() and not p.exists()]
    add(checks, "inventory", "broken symlinks", "PASS" if not broken_links else "FAIL",
        "none" if not broken_links else ", ".join(rel(p, results) for p in broken_links[:100]))

    for path in files:
        name = path.name.lower()
        rpath = rel(path, results)
        if name.endswith(".gz"):
            ok, detail = check_gzip(path)
            add(checks, "format", "gzip integrity", "PASS" if ok else "FAIL", detail, rpath)
        if name.endswith((".tar.gz", ".tgz", ".tar")):
            ok, detail, count = safe_tar_members(path)
            add(checks, "archive", "tar integrity", "PASS" if ok else "FAIL", f"{detail}; members={count}", rpath)
        if name.endswith((".tsv", ".csv", ".tsv.gz", ".csv.gz")) and path.stat().st_size:
            status, detail = check_tabular(path)
            add(checks, "tabular", "rectangular table", status, json.dumps(detail, sort_keys=True), rpath)
        if name.endswith(".json") and path.stat().st_size:
            try:
                json.loads(path.read_text(encoding="utf-8", errors="strict"))
                add(checks, "format", "JSON parse", "PASS", "valid JSON", rpath)
            except Exception as exc:
                add(checks, "format", "JSON parse", "FAIL", str(exc), rpath)
        if name.endswith(".xml") and path.stat().st_size:
            try:
                ET.parse(path)
                add(checks, "format", "XML parse", "PASS", "valid XML", rpath)
            except Exception as exc:
                add(checks, "format", "XML parse", "FAIL", str(exc), rpath)

    checksum_files = [p for p in files if p.name.endswith((".sha256", "checksums.sha256"))]
    for checksum_file in checksum_files:
        entries = parse_checksum_file(checksum_file)
        if not entries:
            add(checks, "checksum", "checksum manifest syntax", "WARN", "no standard SHA-256 entries", rel(checksum_file, results))
            continue
        passed = missing = mismatched = 0
        for expected, listed in entries:
            candidate = (checksum_file.parent / listed).resolve()
            if not candidate.is_file():
                missing += 1
            elif sha256_file(candidate) == expected:
                passed += 1
            else:
                mismatched += 1
        status = "PASS" if missing == 0 and mismatched == 0 else "FAIL"
        add(checks, "checksum", "checksum manifest verification", status,
            f"entries={len(entries)} passed={passed} missing={missing} mismatched={mismatched}", rel(checksum_file, results))

    trace_candidates = sorted(results.glob(f"*trace*{job_id}*.tsv")) or sorted(results.glob("*trace*.tsv"))
    if trace_candidates:
        for trace in trace_candidates:
            try:
                detail = parse_trace(trace)
                status = "PASS" if detail["failed_or_other"] == 0 and detail["duplicate_task_keys"] == 0 else "FAIL"
                add(checks, "execution", "Nextflow trace", status, json.dumps(detail, sort_keys=True), rel(trace, results))
            except Exception as exc:
                add(checks, "execution", "Nextflow trace", "FAIL", str(exc), rel(trace, results))
    else:
        add(checks, "execution", "Nextflow trace", "WARN", "trace file not found")

    log_candidates = [project / ".nextflow.log", project / f"pgtk-wrapper-{job_id}.log"]
    for log in log_candidates:
        if log.is_file():
            text = log.read_text(encoding="utf-8", errors="replace")
            if log.name == ".nextflow.log":
                complete = "Execution complete -- Goodbye" in text
                stats_ok = bool(re.search(r"failedCount=0.*abortedCount=0", text))
                add(checks, "execution", "Nextflow completion log", "PASS" if complete and stats_ok else "FAIL",
                    f"execution_complete={complete}; zero_failed_aborted={stats_ok}", rel(log, project))
            else:
                normal = "Job exited normally" in text or re.search(r"completed at", text, re.I)
                add(checks, "execution", "wrapper completion log", "PASS" if normal else "WARN", f"normal_completion={bool(normal)}", rel(log, project))
        else:
            add(checks, "execution", log.name, "WARN", "not found", rel(log, project))

    failure_ledgers = sorted(results.glob(f"failure_logs/{job_id}/failure_ledger.tsv"))
    if failure_ledgers:
        for ledger in failure_ledgers:
            with ledger.open("rt", encoding="utf-8", errors="replace") as handle:
                nonempty = [line for line in handle if line.strip()]
            rows = max(0, len(nonempty) - 1)
            add(checks, "execution", "failure ledger", "PASS" if rows == 0 else "FAIL", f"failure_rows={rows}", rel(ledger, results))
    else:
        add(checks, "execution", "failure ledger", "WARN", "not found")

    validation_archives = sorted(results.glob(f"validation/*{job_id}*.tar.gz"))
    add(checks, "validation", "validation archive present", "PASS" if validation_archives else "FAIL", str(len(validation_archives)))

    vcf_files = [p for p in files if p.name.endswith((".vcf", ".vcf.gz"))]
    for vcf in vcf_files:
        rpath = rel(vcf, results)
        try:
            sample = read_text_sample(vcf, 3_000_000)
            has_fileformat = "##fileformat=VCF" in sample
            has_columns = any(line.startswith("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO") for line in sample.splitlines())
            data_seen = any(line and not line.startswith("#") for line in sample.splitlines())
            add(checks, "VCF", "VCF header", "PASS" if has_fileformat and has_columns else "FAIL",
                f"fileformat={has_fileformat}; columns={has_columns}; data_seen_in_sample={data_seen}", rpath)
            if vcf.name.endswith(".vcf.gz"):
                indexes = [Path(str(vcf) + ".tbi"), Path(str(vcf) + ".csi")]
                add(checks, "VCF", "compressed VCF index", "PASS" if any(p.is_file() for p in indexes) else "FAIL",
                    "index present" if any(p.is_file() for p in indexes) else "missing .tbi/.csi", rpath)
        except Exception as exc:
            add(checks, "VCF", "VCF readability", "FAIL", str(exc), rpath)

    bam_files = [p for p in files if p.name.endswith((".bam", ".cram"))]
    samtools = shutil.which("samtools")
    for alignment in bam_files:
        rpath = rel(alignment, results)
        if alignment.suffix == ".bam":
            indexes = [Path(str(alignment) + ".bai"), alignment.with_suffix(".bai")]
        else:
            indexes = [Path(str(alignment) + ".crai"), alignment.with_suffix(".crai")]
        existing_indexes = [path for path in indexes if path.is_file()]
        add(checks, "alignment", "alignment index", "PASS" if existing_indexes else "FAIL",
            "index present" if existing_indexes else "index missing", rpath)
        if existing_indexes:
            newest_index = max(existing_indexes, key=lambda path: path.stat().st_mtime_ns)
            index_current = newest_index.stat().st_mtime_ns >= alignment.stat().st_mtime_ns
            freshness = "index timestamp current" if index_current else "index timestamp older than alignment; regenerate the index after rewriting the alignment"
            add(checks, "alignment", "alignment index freshness", "PASS" if index_current else "WARN", freshness, rpath)
        if samtools:
            rc, output = run_command([samtools, "quickcheck", "-v", str(alignment)], timeout=600)
            add(checks, "alignment", "samtools quickcheck", "PASS" if rc == 0 else "FAIL", output or "OK", rpath)
        else:
            add(checks, "alignment", "samtools quickcheck", "WARN", "samtools not available", rpath)

    bed_like = [p for p in files if p.name.endswith((".bed", ".bedpe"))]
    for bed in bed_like:
        invalid = 0
        rows = 0
        try:
            with bed.open("rt", encoding="utf-8", errors="replace") as handle:
                for line in handle:
                    if not line.strip() or line.startswith(("#", "track", "browser")):
                        continue
                    fields = line.rstrip("\n").split("\t")
                    rows += 1
                    required = 6 if bed.suffix == ".bedpe" else 3
                    try:
                        valid = len(fields) >= required and int(fields[1]) >= 0 and int(fields[2]) > int(fields[1])
                        if bed.suffix == ".bedpe":
                            valid = valid and int(fields[4]) >= 0 and int(fields[5]) > int(fields[4])
                    except (ValueError, IndexError):
                        valid = False
                    invalid += int(not valid)
            add(checks, "coordinates", "BED/BEDPE coordinates", "PASS" if invalid == 0 else "FAIL",
                f"rows={rows}; invalid={invalid}", rel(bed, results))
        except Exception as exc:
            add(checks, "coordinates", "BED/BEDPE readability", "FAIL", str(exc), rel(bed, results))

    if run_source_checks:
        py_files = sorted(project.glob("*.py"))
        for source in py_files:
            rc, output = run_command([sys.executable, "-m", "py_compile", str(source)], cwd=project)
            add(checks, "source", "Python syntax", "PASS" if rc == 0 else "FAIL", output or "OK", rel(source, project))
        for source in sorted(project.glob("*.sh")) + ([project / "scratch.slurm"] if (project / "scratch.slurm").is_file() else []):
            rc, output = run_command(["bash", "-n", str(source)], cwd=project)
            add(checks, "source", "shell syntax", "PASS" if rc == 0 else "FAIL", output or "OK", rel(source, project))

    severity = Counter(check.status for check in checks)
    overall = "FAIL" if severity["FAIL"] else ("WARN" if severity["WARN"] else "PASS")
    summary = {
        "audit_version": "1.0.0", "job_id": job_id, "project": str(project), "results": str(results),
        "overall_status": overall, "counts": dict(severity), "result_files": len(files),
        "result_bytes": total_bytes, "extension_counts": dict(extension_counts),
        "started_epoch": started, "finished_epoch": time.time(),
        "limitations": [
            "This is an automated structural, integrity, provenance, and cross-format audit.",
            "It does not establish clinical validity or biological truth.",
            "Read-level claims require BAM/CRAM files and indexes to be present.",
            "Hashes skipped by the configured size limit are explicitly marked in results_inventory.tsv.",
        ],
    }

    with (report_dir / "checks.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["category", "name", "status", "detail", "path"])
        for check in checks:
            writer.writerow([check.category, check.name, check.status, check.detail, check.path])
    (report_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    failures = [c for c in checks if c.status == "FAIL"]
    warnings = [c for c in checks if c.status == "WARN"]
    report = [
        f"# PGTK independent audit for job {job_id}", "",
        f"Overall status: **{overall}**", "",
        f"Project: `{project}`", f"Results: `{results}`", "",
        "## Summary", "",
        f"- Result files: {len(files)}", f"- Result size: {human_size(total_bytes)}",
        f"- PASS: {severity['PASS']}", f"- WARN: {severity['WARN']}", f"- FAIL: {severity['FAIL']}", "",
        "## Failures", "",
    ]
    report += [f"- `{c.path}`: {c.category} / {c.name}: {c.detail}" for c in failures] or ["- None"]
    report += ["", "## Warnings", ""]
    report += [f"- `{c.path}`: {c.category} / {c.name}: {c.detail}" for c in warnings] or ["- None"]
    report += ["", "## Scope and limitations", ""] + [f"- {item}" for item in summary["limitations"]]
    report += ["", "## Included evidence", "", "- checks.tsv", "- summary.json", "- results_inventory.tsv", "- results_checksums.sha256", ""]
    (report_dir / "REPORT.md").write_text("\n".join(report), encoding="utf-8")

    evidence_manifest = report_dir / "audit_bundle_checksums.sha256"
    with evidence_manifest.open("w", encoding="utf-8") as handle:
        for path in sorted(report_dir.iterdir()):
            if path.is_file() and path != evidence_manifest:
                handle.write(f"{sha256_file(path)}  {path.name}\n")

    bundle = output_root / f"PGTK-independent-audit-{job_id}.tar.gz"
    if bundle.exists():
        bundle.unlink()
    with tarfile.open(bundle, "w:gz", format=tarfile.PAX_FORMAT) as archive:
        archive.add(report_dir, arcname=report_dir.name)
    bundle_sha = output_root / f"PGTK-independent-audit-{job_id}.tar.gz.sha256"
    bundle_sha.write_text(f"{sha256_file(bundle)}  {bundle.name}\n", encoding="utf-8")
    return bundle, summary


def main() -> int:
    parser = argparse.ArgumentParser(description="Audit a completed PGTK run and create a shareable evidence bundle.")
    parser.add_argument("--project-dir", required=True, type=Path)
    parser.add_argument("--results-dir", required=True, type=Path)
    parser.add_argument("--job-id", required=True)
    parser.add_argument("--output-dir", type=Path, default=Path.cwd())
    parser.add_argument("--hash-large-files", action="store_true", help="Hash files larger than --max-hash-gib.")
    parser.add_argument("--max-hash-gib", type=float, default=5.0)
    parser.add_argument("--skip-source-checks", action="store_true")
    args = parser.parse_args()
    bundle, summary = audit(
        args.project_dir, args.results_dir, args.job_id, args.output_dir.resolve(),
        args.hash_large_files, args.max_hash_gib, not args.skip_source_checks,
    )
    print(f"OVERALL: {summary['overall_status']}")
    print(f"REPORT: {bundle.with_suffix('').with_suffix('') / 'REPORT.md'}")
    print(f"BUNDLE: {bundle}")
    print(f"CHECKSUM: {bundle}.sha256")
    return 1 if summary["overall_status"] == "FAIL" else 0


if __name__ == "__main__":
    raise SystemExit(main())
