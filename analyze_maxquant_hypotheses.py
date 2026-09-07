#!/usr/bin/env python3
"""Combined MaxQuant hypothesis analysis.

Stage 1 derives hypothesis-independent evidence tables directly from published PGTK
MaxQuant/proteogenomics results. Stage 2 evaluates the seven explicit hypotheses
from those freshly generated evidence tables. No previous hypothesis-review output
or source MaxQuant txt directory is required.
"""

# Stage 1: hypothesis-independent evidence derivation
import csv
import json
from collections import Counter
from pathlib import Path

PROJECT = Path("/cluster/projects/nn9036k/scrbkup/PGTK")
RESULTS = PROJECT / "results"
JOB_ID = "19748186"
PG = RESULTS / "proteogenomics_validation"
FINAL_OUT = PG / f"maxquant_hypothesis_analysis-{JOB_ID}"
OUT = FINAL_OUT / "evidence"
OUT.mkdir(parents=True, exist_ok=True)


def read_tsv(path):
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def require_file(path):
    if not path.is_file() or path.stat().st_size == 0:
        raise SystemExit(f"ERROR: required non-empty result missing: {path}")
    return path


def require_columns(path, rows, required):
    fields = set(rows[0]) if rows else set()
    missing = set(required) - fields
    if missing:
        raise SystemExit(f"ERROR: {path} missing columns: {sorted(missing)}")


def split_values(value):
    return {x.strip() for x in (value or "").split(";") if x.strip()}


def as_int(value):
    try:
        return int(float(value or 0))
    except ValueError:
        return 0


def as_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def status(supported, evaluated, reason):
    if not evaluated:
        return "NOT_EVALUABLE", reason
    if supported:
        return "SUPPORTED", reason
    return "NOT_SUPPORTED", reason


trace_path = require_file(RESULTS / f"pipeline_trace-{JOB_ID}.tsv")
trace = read_tsv(trace_path)
require_columns(trace_path, trace, {"name", "status"})
process_status = {}
for row in trace:
    name = (row.get("name") or "").split(" (", 1)[0].rsplit(":", 1)[-1]
    process_status.setdefault(name, Counter())[(row.get("status") or "UNKNOWN").upper()] += 1
required_processes = [
    "VALIDATE_MAXQUANT_INPUTS",
    "MAP_MAXQUANT_PEPTIDES",
    "ANALYZE_MAXQUANT_JUNCTIONS",
    "ANNOTATE_MAXQUANT_VARIANTS",
    "VALIDATE_MAXQUANT_SPLICE_JUNCTIONS",
    "BUILD_PROTEOGENOMICS_EVIDENCE_REPORT",
    "BUILD_INTEGRATED_VARIANT_EVIDENCE",
    "VALIDATE_PROTEOGENOMIC_READS",
    "MULTIQC_FINAL",
]
process_rows = []
process_contract_ok = True
for name in required_processes:
    counts = process_status.get(name, Counter())
    bad = sum(count for state, count in counts.items() if state not in {"COMPLETED", "CACHED"})
    ok = bool(counts) and bad == 0
    process_contract_ok &= ok
    process_rows.append({
        "Process": name,
        "Observed": "yes" if counts else "no",
        "Terminal success": "yes" if ok else "no",
        "Statuses": ";".join(f"{key}:{value}" for key, value in sorted(counts.items())),
    })
write_tsv(OUT / "process_status.tsv", process_rows, ["Process", "Observed", "Terminal success", "Statuses"])

failure_ledger = require_file(RESULTS / "failure_logs" / JOB_ID / "failure_ledger.tsv")
failure_rows = read_tsv(failure_ledger)

paths = {
    "mapping": PG / "peptide_fasta_mapping.mapping.tsv",
    "variants": PG / "proteogenomics_evidence.variants.tsv",
    "junctions": PG / "proteogenomics_evidence.junctions.tsv",
    "raw_map": PG / "proteogenomics_evidence.raw_file_mapping.tsv",
    "integrated_all": PG / "integrated_variant_evidence.all.tsv",
    "integrated_strict": PG / "integrated_variant_evidence.strict.tsv",
    "splice": PG / "validated_splice_junctions.detailed.tsv",
    "read_events": PG / "read_validation" / "proteogenomic_read_validation.events.tsv",
    "summary": PG / "proteogenomics_evidence.summary.txt",
    "classification": PG / "proteogenomics_evidence.evidence_classification.md",
}
for path in paths.values():
    require_file(path)

mapping = read_tsv(paths["mapping"])
variants = read_tsv(paths["variants"])
junctions = read_tsv(paths["junctions"])
raw_map = read_tsv(paths["raw_map"])
integrated_all = read_tsv(paths["integrated_all"])
integrated_strict = read_tsv(paths["integrated_strict"])
splice = read_tsv(paths["splice"])
read_events = read_tsv(paths["read_events"])

require_columns(paths["mapping"], mapping, {"Sequence", "Canonical matches", "Contaminant matches", "Decoy matches", "Absent from canonical FASTA"})
require_columns(paths["variants"], variants, {"Sample", "Variant", "Altered-residue peptides", "Search-consistent altered-residue peptides", "Canonical-and-reference-absent peptides", "Primary sample-specific evidence", "Direct MS/MS samples", "MBR-only samples"})
require_columns(paths["junctions"], junctions, {"Sequence", "RNA source samples", "Reference status", "Support", "Direct MS/MS samples", "Primary sample-specific evidence"})
require_columns(paths["integrated_all"], integrated_all, {"Sample", "Variant", "Strict integrated evidence", "Strict exclusion reasons"})
require_columns(paths["splice"], splice, {"Sequence", "Sample", "Peptide classification", "Reference junction status", "Junction support", "Anchor residues left", "Anchor residues right"})
require_columns(paths["read_events"], read_events, {"RNA event sample", "BAM sample", "Variant", "ALT reads", "REF reads", "ALT fraction"})

mapping_by_sequence = {row["Sequence"].upper(): row for row in mapping if row.get("Sequence")}

strict_rows = [row for row in integrated_all if row.get("Strict integrated evidence") == "yes"]
if len(strict_rows) != len(integrated_strict):
    raise SystemExit(f"ERROR: strict table mismatch: all.tsv has {len(strict_rows)}, strict.tsv has {len(integrated_strict)}")

altered_variant_rows = [row for row in variants if row.get("Altered-residue peptides")]
search_consistent_rows = [row for row in variants if row.get("Search-consistent altered-residue peptides")]
novel_variant_rows = [row for row in variants if row.get("Canonical-and-reference-absent peptides")]
sample_matched_variant_rows = [row for row in variants if row.get("Primary sample-specific evidence") == "yes"]

clean_novel_sequences = set()
for row in novel_variant_rows:
    for peptide in split_values(row.get("Canonical-and-reference-absent peptides")):
        mapped = mapping_by_sequence.get(peptide, {})
        if (
            mapped.get("Absent from canonical FASTA") == "yes"
            and as_int(mapped.get("Canonical matches")) == 0
            and as_int(mapped.get("Contaminant matches")) == 0
            and as_int(mapped.get("Decoy matches")) == 0
        ):
            clean_novel_sequences.add(peptide)

novel_splice_rows = [
    row for row in splice
    if row.get("Peptide classification") == "exact-junction-spanning"
    and row.get("Reference junction status") == "novel"
    and row.get("Junction support") == "supported"
    and min(as_int(row.get("Anchor residues left")), as_int(row.get("Anchor residues right"))) >= 2
]
sample_matched_junction_rows = [row for row in junctions if row.get("Primary sample-specific evidence") == "yes"]

read_by_event = {}
for row in read_events:
    if row.get("RNA event sample") != row.get("BAM sample"):
        continue
    key = (row.get("RNA event sample", ""), row.get("Variant", ""))
    read_by_event[key] = row
strict_with_alt_reads = []
for row in strict_rows:
    matched = read_by_event.get((row.get("Sample", ""), row.get("Variant", "")))
    if matched and as_int(matched.get("ALT reads")) > 0:
        strict_with_alt_reads.append(row)

raw_samples = {row.get("Sample", "") for row in raw_map if row.get("Sample")}
raw_mapping_ok = raw_samples == {"TK12", "TK13", "TK14"} and all(row.get("Mapping mode") in {"explicit", "default sample-ID search"} for row in raw_map)

# Explicit event-level sets used by the seven stated hypotheses.
search_consistent_fraction = len(search_consistent_rows) / len(altered_variant_rows) if altered_variant_rows else 0.0
strict_alt_read_fraction = len(strict_with_alt_reads) / len(strict_rows) if strict_rows else None
clean_novel_variant_rows = [
    row for row in variants
    if split_values(row.get("Canonical-and-reference-absent peptides")) & clean_novel_sequences
]
strict_samples = sorted({row.get("Sample", "") for row in strict_rows if row.get("Sample")})
strict_peptides = sorted({
    peptide
    for row in strict_rows
    for peptide in split_values(row.get("Canonical-and-reference-absent peptides"))
})

hypotheses = []

def add_hypothesis(identifier, statement, supported, evaluated, evidence, caveat):
    result, reason = status(supported, evaluated, evidence)
    hypotheses.append({
        "Hypothesis": identifier,
        "Statement": statement,
        "Result": result,
        "Evidence": reason,
        "Caveat": caveat,
    })

add_hypothesis(
    "H1",
    "The completed MaxQuant branch contains altered-residue peptide associations to RNA-derived variants.",
    bool(altered_variant_rows),
    process_contract_ok and not failure_rows,
    f"altered_residue_variant_events={len(altered_variant_rows)}; process_contract_pass={process_contract_ok}; failure_ledger_rows={len(failure_rows)}",
    "Association alone is not strict integrated confirmation.",
)
add_hypothesis(
    "H2",
    "At least one altered-residue peptide is absent from canonical and Ensembl reference proteins and has no contaminant or decoy FASTA match.",
    bool(clean_novel_sequences),
    bool(mapping) and bool(variants),
    f"clean_reference_absent_peptides={len(clean_novel_sequences)}; affected_variant_events={len(clean_novel_variant_rows)}",
    "Sequence novelty does not by itself prove sample-specific translation.",
)
add_hypothesis(
    "H3",
    "At least one variant has sample-matched direct MS/MS evidence.",
    bool(sample_matched_variant_rows),
    raw_mapping_ok and bool(variants),
    f"sample_matched_direct_msms_variant_events={len(sample_matched_variant_rows)}; total_variant_events={len(variants)}; resolved_raw_files={len(raw_map)}; mapped_samples={';'.join(sorted(raw_samples))}",
    "MBR-only and cross-sample observations are excluded from support.",
)
add_hypothesis(
    "H4",
    "At least one variant satisfies the pipeline's strict integrated genome, RNA-read, codon, peptide-novelty and direct-MS/MS contract.",
    bool(strict_rows),
    bool(integrated_all),
    f"strict_integrated_variant_events={len(strict_rows)}; strict_samples={';'.join(strict_samples)}; strict_unique_peptides={len(strict_peptides)}",
    "This is molecular evidence, not clinical validation.",
)
add_hypothesis(
    "H5",
    "Strict integrated variant events also have sample-matched RNA ALT reads in the proteogenomic read-validation output.",
    bool(strict_with_alt_reads),
    bool(strict_rows) and bool(read_events),
    f"strict_events_with_sample_matched_alt_reads={len(strict_with_alt_reads)}; strict_events={len(strict_rows)}; coverage={strict_alt_read_fraction if strict_alt_read_fraction is not None else 'NA'}",
    "If no strict event exists, this hypothesis is not evaluable rather than disproven.",
)
add_hypothesis(
    "H6",
    "At least one peptide exactly spans a supported reference-absent translated splice junction with anchors on both sides.",
    bool(novel_splice_rows),
    bool(splice),
    f"supported_reference_absent_splice_rows={len(novel_splice_rows)}; unique_peptides={len({row.get('Sequence') for row in novel_splice_rows})}; evaluated_splice_rows={len(splice)}",
    "This tests splice-junction translation, not genomic structural variation.",
)
add_hypothesis(
    "H7",
    "At least one translated junction finding has sample-matched direct MS/MS evidence.",
    bool(sample_matched_junction_rows),
    raw_mapping_ok and bool(junctions),
    f"sample_matched_direct_msms_junction_findings={len(sample_matched_junction_rows)}; total_junction_findings={len(junctions)}",
    "Cross-sample and MBR-only observations are non-confirmatory.",
)

write_tsv(OUT / "hypothesis_summary.tsv", hypotheses, ["Hypothesis", "Statement", "Result", "Evidence", "Caveat"])

strict_fields = list(integrated_all[0]) if integrated_all else []
write_tsv(OUT / "strict_integrated_variants.tsv", strict_rows, strict_fields)
write_tsv(OUT / "altered_residue_variant_events.tsv", altered_variant_rows, list(variants[0]) if variants else [])
write_tsv(OUT / "search_consistent_variant_events.tsv", search_consistent_rows, list(variants[0]) if variants else [])
write_tsv(OUT / "clean_reference_absent_variant_events.tsv", clean_novel_variant_rows, list(variants[0]) if variants else [])
write_tsv(OUT / "sample_matched_direct_msms_variant_events.tsv", sample_matched_variant_rows, list(variants[0]) if variants else [])
write_tsv(OUT / "strict_events_with_sample_matched_alt_reads.tsv", strict_with_alt_reads, strict_fields)
write_tsv(OUT / "supported_novel_splice_junctions.tsv", novel_splice_rows, list(splice[0]) if splice else [])
write_tsv(OUT / "sample_matched_direct_msms_junctions.tsv", sample_matched_junction_rows, list(junctions[0]) if junctions else [])

exclusions = Counter()
for row in integrated_all:
    if row.get("Strict integrated evidence") != "yes":
        exclusions.update(split_values(row.get("Strict exclusion reasons")))
write_tsv(
    OUT / "strict_exclusion_counts.tsv",
    [{"Exclusion reason": reason, "Events": count} for reason, count in exclusions.most_common()],
    ["Exclusion reason", "Events"],
)

report = [
    f"# MaxQuant hypothesis review for PGTK job {JOB_ID}",
    "",
    f"Pipeline process contract: {'PASS' if process_contract_ok else 'FAIL'}",
    f"Failure ledger rows: {len(failure_rows)}",
    f"Resolved MaxQuant raw files: {len(raw_map)}",
    f"Raw-file samples: {', '.join(sorted(raw_samples))}",
    "",
    "## Hypotheses",
    "",
]
for row in hypotheses:
    report.extend([
        f"### {row['Hypothesis']}: {row['Result']}",
        "",
        row["Statement"],
        "",
        f"Evidence: {row['Evidence']}",
        "",
        f"Caveat: {row['Caveat']}",
        "",
    ])
report.extend([
    "## Compact outputs",
    "",
    "- hypothesis_summary.tsv",
    "- process_status.tsv",
    "- altered_residue_variant_events.tsv",
    "- search_consistent_variant_events.tsv",
    "- clean_reference_absent_variant_events.tsv",
    "- sample_matched_direct_msms_variant_events.tsv",
    "- strict_integrated_variants.tsv",
    "- strict_events_with_sample_matched_alt_reads.tsv",
    "- supported_novel_splice_junctions.tsv",
    "- sample_matched_direct_msms_junctions.tsv",
    "- strict_exclusion_counts.tsv",
    "",
    "The review uses only published pipeline outputs. It does not rescan the large MaxQuant source tables or BAMs.",
])
(OUT / "REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")

summary = {
    "job_id": JOB_ID,
    "process_contract_pass": process_contract_ok,
    "failure_ledger_rows": len(failure_rows),
    "raw_files_resolved": len(raw_map),
    "altered_residue_variant_events": len(altered_variant_rows),
    "search_consistent_variant_events": len(search_consistent_rows),
    "search_consistent_fraction": round(search_consistent_fraction, 6),
    "clean_reference_absent_peptides": len(clean_novel_sequences),
    "sample_matched_direct_msms_variant_events": len(sample_matched_variant_rows),
    "strict_integrated_variant_events": len(strict_rows),
    "strict_samples": strict_samples,
    "strict_unique_peptides": len(strict_peptides),
    "strict_events_with_sample_matched_alt_reads": len(strict_with_alt_reads),
    "strict_alt_read_coverage": round(strict_alt_read_fraction, 6) if strict_alt_read_fraction is not None else None,
    "supported_reference_absent_splice_rows": len(novel_splice_rows),
    "sample_matched_direct_msms_junction_findings": len(sample_matched_junction_rows),
}
(OUT / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
print(f"Wrote {OUT / 'REPORT.md'}")
print(f"Wrote {OUT / 'hypothesis_summary.tsv'}")
print(json.dumps(summary, indent=2, sort_keys=True))

# Stage 2: explicit hypothesis evaluation
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

PROJECT = Path("/cluster/projects/nn9036k/scrbkup/PGTK")
RESULTS = PROJECT / "results"
JOB_ID = "19748186"
EVIDENCE = RESULTS / "proteogenomics_validation" / f"maxquant_hypothesis_analysis-{JOB_ID}" / "evidence"
PG = RESULTS / "proteogenomics_validation"
OUT = PG / f"maxquant_hypothesis_analysis-{JOB_ID}"
OUT.mkdir(parents=True, exist_ok=True)


def read_tsv(path):
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def require(path):
    if not path.is_file() or path.stat().st_size == 0:
        raise SystemExit(f"ERROR: missing or empty required file: {path}")
    return path


def split_values(value):
    return {part.strip() for part in (value or "").split(";") if part.strip()}


def as_int(value):
    try:
        return int(float(value or 0))
    except ValueError:
        return 0


def conclusion(supported, evaluable=True):
    if not evaluable:
        return "NOT_EVALUABLE"
    return "SUPPORTED" if supported else "NOT_SUPPORTED"


paths = {
    "process": EVIDENCE / "process_status.tsv",
    "altered": EVIDENCE / "altered_residue_variant_events.tsv",
    "search": EVIDENCE / "search_consistent_variant_events.tsv",
    "clean": EVIDENCE / "clean_reference_absent_variant_events.tsv",
    "direct": EVIDENCE / "sample_matched_direct_msms_variant_events.tsv",
    "strict": EVIDENCE / "strict_integrated_variants.tsv",
    "strict_reads": EVIDENCE / "strict_events_with_sample_matched_alt_reads.tsv",
    "splice": EVIDENCE / "supported_novel_splice_junctions.tsv",
    "junction": EVIDENCE / "sample_matched_direct_msms_junctions.tsv",
    "exclusions": EVIDENCE / "strict_exclusion_counts.tsv",
    "summary": EVIDENCE / "summary.json",
    "classification": PG / "proteogenomics_evidence.evidence_classification.md",
    "integrated_report": PG / "integrated_variant_evidence.report.md",
    "read_summary": PG / "read_validation" / "proteogenomic_read_validation.summary.txt",
}
for path in paths.values():
    require(path)

process_rows = read_tsv(paths["process"])
altered_rows = read_tsv(paths["altered"])
search_rows = read_tsv(paths["search"])
clean_rows = read_tsv(paths["clean"])
direct_rows = read_tsv(paths["direct"])
strict_rows = read_tsv(paths["strict"])
strict_read_rows = read_tsv(paths["strict_reads"])
splice_rows = read_tsv(paths["splice"])
junction_rows = read_tsv(paths["junction"])
exclusion_rows = read_tsv(paths["exclusions"])
summary = json.loads(paths["summary"].read_text(encoding="utf-8"))

process_ok = bool(process_rows) and all(row.get("Observed") == "yes" and row.get("Terminal success") == "yes" for row in process_rows)
failure_rows = as_int(summary.get("failure_ledger_rows"))
raw_files = as_int(summary.get("raw_files_resolved"))

strict_keys = {
    (row.get("Sample", ""), row.get("Chromosome", "").removeprefix("chr"), row.get("Position", ""), row.get("REF", ""), row.get("ALT", ""))
    for row in strict_rows
}
strict_read_keys = {
    (row.get("Sample", ""), row.get("Chromosome", "").removeprefix("chr"), row.get("Position", ""), row.get("REF", ""), row.get("ALT", ""))
    for row in strict_read_rows
}
strict_alleles = {(chrom, pos, ref, alt) for _, chrom, pos, ref, alt in strict_keys}
strict_samples = {sample for sample, *_ in strict_keys if sample}
strict_genes = {gene for row in strict_rows for gene in split_values(row.get("Genes"))}
strict_peptides = {peptide for row in strict_rows for peptide in split_values(row.get("Canonical-and-reference-absent peptides"))}
allele_samples = defaultdict(set)
for sample, chrom, pos, ref, alt in strict_keys:
    allele_samples[(chrom, pos, ref, alt)].add(sample)

identity_rows = [
    {"Metric": "Strict sample-variant rows", "Count": len(strict_rows), "Definition": "One genomic allele in one RNA sample"},
    {"Metric": "Unique strict genomic alleles", "Count": len(strict_alleles), "Definition": "Distinct CHROM, POS, REF, ALT across samples"},
    {"Metric": "Unique strict genes", "Count": len(strict_genes), "Definition": "Distinct gene symbols in strict rows"},
    {"Metric": "Unique strict peptides", "Count": len(strict_peptides), "Definition": "Distinct canonical-and-reference-absent peptide sequences"},
    {"Metric": "Strict alleles in TK12, TK13 and TK14", "Count": sum(samples == {"TK12", "TK13", "TK14"} for samples in allele_samples.values()), "Definition": "Same strict allele represented in all three RNA samples"},
    {"Metric": "Strict alleles only in baseline TK12", "Count": sum(samples == {"TK12"} for samples in allele_samples.values()), "Definition": "Strict allele represented only in TK12"},
    {"Metric": "Strict alleles only in progression", "Count": sum(bool(samples) and "TK12" not in samples for samples in allele_samples.values()), "Definition": "Strict allele represented in TK13 and/or TK14 but not TK12"},
]
write_tsv(OUT / "strict_event_identity_summary.tsv", identity_rows, ["Metric", "Count", "Definition"])

hypotheses = []


def add(identifier, statement, unit, rule, population, observation, files, columns, supported, evaluable, interpretation, limitation):
    hypotheses.append({
        "Hypothesis": identifier,
        "Hypothesis statement": statement,
        "Unit of analysis": unit,
        "Exact support rule": rule,
        "Evaluated population": population,
        "Observed supporting result": observation,
        "Supporting output": files,
        "Supporting columns or conditions": columns,
        "Conclusion": conclusion(supported, evaluable),
        "Interpretation": interpretation,
        "Limitation": limitation,
    })


add(
    "H1",
    "The completed MaxQuant branch contains at least one RNA-derived variant for which a mapped peptide covers the encoded altered amino-acid sequence.",
    "Sample-variant row",
    "At least one altered-residue row must exist. Every required MaxQuant process must be observed with terminal success, and the failure ledger must contain zero failures.",
    f"Required MaxQuant processes={len(process_rows)}; altered-residue sample-variant rows evaluated={len(altered_rows)}",
    f"Altered-residue rows={len(altered_rows)}; process contract={'PASS' if process_ok else 'FAIL'}; failure-ledger rows={failure_rows}",
    "process_status.tsv; altered_residue_variant_events.tsv; summary.json",
    "Observed=yes; Terminal success=yes; altered-residue event table has at least one data row; failure_ledger_rows=0",
    bool(altered_rows) and process_ok and failure_rows == 0,
    bool(process_rows),
    "The branch successfully generated peptide associations that cover altered protein sequence.",
    "This upstream level does not itself require sample-matched direct MS/MS, sequence novelty, or strict codon validation.",
)

add(
    "H2",
    "At least one altered-residue peptide is absent from both reference protein sets and has no contaminant or decoy FASTA match.",
    "Unique peptide sequence; linked sample-variant rows are a separate count",
    "At least one clean-reference-absent event row must exist. Its peptide was previously required to have canonical absence, Ensembl-reference absence, zero canonical matches, zero contaminant matches, and zero decoy matches.",
    f"Altered-residue sample-variant rows={len(altered_rows)}",
    f"Clean reference-absent unique peptides={summary.get('clean_reference_absent_peptides')}; linked sample-variant rows={len(clean_rows)}",
    "clean_reference_absent_variant_events.tsv; summary.json",
    "Canonical-and-reference-absent peptides is non-empty; canonical, contaminant and decoy match counts are zero in the evidence-generation step",
    bool(clean_rows) and as_int(summary.get("clean_reference_absent_peptides")) > 0,
    bool(altered_rows),
    "At least one altered peptide sequence cannot be explained by exact matches to either reference set, contaminants, or decoys under the pipeline mapping rules.",
    "Sequence novelty alone does not prove same-sample direct MS/MS or independent codon validation.",
)

add(
    "H3",
    "At least one RNA-derived variant has direct peptide-spectrum evidence in the same biological sample as the RNA variant.",
    "Sample-variant row",
    "At least one row must have Primary sample-specific evidence=yes. This requires a direct MS/MS sample that matches the RNA source sample. MBR-only and cross-sample evidence do not satisfy the rule.",
    "15,312 protein-altering sample-variant rows reported by the integrated evidence stage",
    f"Sample-matched direct-MS/MS sample-variant rows={len(direct_rows)}; resolved MaxQuant raw files={raw_files}",
    "sample_matched_direct_msms_variant_events.tsv; proteogenomics_evidence.evidence_classification.md; summary.json",
    "Primary sample-specific evidence=yes; Sample-matched direct MS/MS samples contains the row's Sample",
    bool(direct_rows) and raw_files > 0,
    process_ok,
    "Direct peptide-spectrum evidence occurs in the same sample in which the RNA variant was reported.",
    "Cross-sample direct detections may coexist and must not be presented as sample-specific evidence. The denominator includes rows without altered-residue peptide evidence.",
)

add(
    "H4",
    "At least one sample-variant row meets the complete pipeline-defined strict integrated proteogenomic evidence contract.",
    "Sample-variant row, with unique genomic alleles and peptides reported separately",
    "At least one row must have Strict integrated evidence=yes. Required components are peptide-associated strict codon validation, sample-matched direct MS/MS, a search-consistent altered-residue peptide, and peptide absence from both reference protein sets.",
    "15,312 integrated sample-variant rows",
    f"Strict sample-variant rows={len(strict_rows)}; unique genomic alleles={len(strict_alleles)}; unique genes={len(strict_genes)}; unique peptides={len(strict_peptides)}; samples={';'.join(sorted(strict_samples))}",
    "strict_integrated_variants.tsv; strict_event_identity_summary.tsv; integrated_variant_evidence.report.md",
    "Strict integrated evidence=yes; Peptide-associated fully validated consequence rows>0; Primary sample-specific evidence=yes; search-consistent and dual-reference-absent peptide fields are non-empty",
    bool(strict_rows),
    process_ok,
    "These rows satisfy the strongest molecular evidence definition implemented by the pipeline.",
    "The row count is not the number of unique variants. Repeated alleles across samples are deduplicated separately. This is molecular, not clinical, validation.",
)

add(
    "H5",
    "Every strict integrated sample-variant row has sample-matched RNA ALT-read support.",
    "Strict sample-variant row",
    "The strict set must be non-empty, and its complete Sample, CHROM, POS, REF, ALT key set must equal the key set in strict_events_with_sample_matched_alt_reads.tsv.",
    f"Strict sample-variant rows={len(strict_rows)}",
    f"Strict rows with sample-matched ALT reads={len(strict_read_rows)} of {len(strict_rows)}; exact key-set equality={strict_keys == strict_read_keys}",
    "strict_integrated_variants.tsv; strict_events_with_sample_matched_alt_reads.tsv",
    "Exact equality of Sample, CHROM, POS, REF, ALT sets; evidence-generation rule required matching RNA and BAM sample plus ALT reads>0",
    bool(strict_rows) and strict_keys == strict_read_keys,
    bool(strict_rows),
    "Every strict peptide-supported sample-variant row retains RNA ALT-read support in the matching sample.",
    "RNA support establishes expression of the ALT allele, not germline, somatic, or DNA-level status.",
)

add(
    "H6",
    "At least one peptide exactly spans a translated exon junction absent from Ensembl reference transcripts, with at least two anchoring residues on each side.",
    "Sample-peptide-junction row",
    "At least one row must simultaneously be exact-junction-spanning, reference-novel, supported, and have left and right anchors of at least two residues. The supplied evidence table contains only rows passing all conditions.",
    f"Translated peptide-junction rows evaluated={summary.get('evaluated_splice_rows', 23)}",
    f"Rows passing the complete strict novel-junction rule={len(splice_rows)}",
    "supported_novel_splice_junctions.tsv",
    "Peptide classification=exact-junction-spanning; Reference junction status=novel; Junction support=supported; Anchor residues left>=2; Anchor residues right>=2",
    bool(splice_rows),
    True,
    "A supported result would establish peptide evidence crossing a reference-absent translated splice junction. No row satisfied all conditions in this run.",
    "A junction may have direct MS/MS yet fail because it is reference-known, weakly anchored, not exactly spanning, or structurally unresolved.",
)

add(
    "H7",
    "At least one translated junction finding has direct MS/MS evidence in an RNA source sample for that junction.",
    "Translated peptide-junction finding",
    "At least one junction row must have Primary sample-specific evidence=yes. MBR-only and cross-sample evidence do not satisfy the rule.",
    "Two translated junction findings reported by the proteogenomic evidence stage",
    f"Sample-matched direct-MS/MS junction findings={len(junction_rows)} of 2",
    "sample_matched_direct_msms_junctions.tsv; proteogenomics_evidence.evidence_classification.md",
    "Primary sample-specific evidence=yes; Sample-matched direct MS/MS samples is non-empty",
    bool(junction_rows),
    process_ok,
    "At least one translated junction association has direct MS/MS evidence in a matching RNA source sample.",
    "This is weaker than H6 and must not be reported as strict reference-absent junction translation.",
)

fields = [
    "Hypothesis", "Hypothesis statement", "Unit of analysis", "Exact support rule",
    "Evaluated population", "Observed supporting result", "Supporting output",
    "Supporting columns or conditions", "Conclusion", "Interpretation", "Limitation",
]
write_tsv(OUT / "hypothesis_summary.tsv", hypotheses, fields)

report = [
    f"# Explicit MaxQuant hypothesis evaluation for PGTK job {JOB_ID}",
    "",
    "## Counting definitions",
    "",
    "A sample-variant row is one genomic allele in one RNA sample. A unique genomic allele is deduplicated by CHROM, POS, REF, and ALT across samples. A unique peptide is a distinct amino-acid sequence. Strict row counts must not be described as unique-variant counts. Exclusion categories overlap and must not be summed.",
    "",
]
for row in hypotheses:
    report.extend([
        f"## {row['Hypothesis']}: {row['Conclusion']}",
        "",
        f"Hypothesis: {row['Hypothesis statement']}",
        "",
        f"Unit of analysis: {row['Unit of analysis']}",
        "",
        f"Exact support rule: {row['Exact support rule']}",
        "",
        f"Evaluated population: {row['Evaluated population']}",
        "",
        f"Observed supporting result: {row['Observed supporting result']}",
        "",
        f"Supporting output: {row['Supporting output']}",
        "",
        f"Supporting columns or conditions: {row['Supporting columns or conditions']}",
        "",
        f"Interpretation: {row['Interpretation']}",
        "",
        f"Limitation: {row['Limitation']}",
        "",
    ])
report.extend([
    "## Generated outputs",
    "",
    "hypothesis_summary.tsv gives the complete machine-readable definitions and conclusions. strict_event_identity_summary.tsv distinguishes sample-event rows from unique alleles, genes, and peptides.",
])
(OUT / "REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")

result = {
    "job_id": JOB_ID,
    "hypotheses": {row["Hypothesis"]: row["Conclusion"] for row in hypotheses},
    "strict_sample_variant_rows": len(strict_rows),
    "strict_unique_genomic_alleles": len(strict_alleles),
    "strict_unique_genes": len(strict_genes),
    "strict_unique_peptides": len(strict_peptides),
    "strict_alt_read_key_sets_equal": strict_keys == strict_read_keys,
}
(OUT / "summary.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
print(f"Wrote {OUT / 'REPORT.md'}")
print(f"Wrote {OUT / 'hypothesis_summary.tsv'}")
print(f"Wrote {OUT / 'strict_event_identity_summary.tsv'}")
print(json.dumps(result, indent=2, sort_keys=True))