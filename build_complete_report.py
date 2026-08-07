import argparse
import csv
import gzip
import sys
from collections import Counter
from pathlib import Path


def open_text(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def iter_tsv(path):
    with open_text(path) as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def count_records(path):
    with open_text(path) as handle:
        return sum(1 for line in handle if line and not line.startswith("#"))


def count_fasta(path):
    with open_text(path) as handle:
        return sum(1 for line in handle if line.startswith(">"))


def count_gtf(path):
    with open_text(path) as handle:
        return sum(
            1
            for line in handle
            if not line.startswith("#") and "\ttranscript\t" in line
        )


def sample_from_path(path):
    return Path(path).name.split(".", 1)[0]


def escape_markdown(value):
    return str(value or "").replace("|", "\\|").replace("\n", " ")


def write_tsv(path, rows, fields):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def read_multiqc_general_stats(directory):
    candidates = list(Path(directory).rglob("multiqc_general_stats.txt"))
    if not candidates:
        return [], "multiqc_general_stats.txt"
    with candidates[0].open(
        encoding="utf-8",
        errors="replace",
        newline="",
    ) as handle:
        return list(csv.DictReader(handle, delimiter="\t")), candidates[0].name


def audit_fields():
    return [
        "Class",
        "Sample",
        "Event type",
        "Event",
        "Initial finding",
        "Validation rule",
        "Observed evidence",
        "Status",
        "Failure codes",
        "Failure explanation",
        "Required resolution",
        "Depth",
        "ALT reads",
        "ALT fraction",
        "Consequences",
        "Split reads",
        "Discordant mates",
        "Confidence",
        "Reading frame",
        "Junctions",
        "Junction read counts",
        "Minimum junction reads",
        "Source",
    ]


def stream_audits(args, audit_path, failure_path):
    fields = audit_fields()
    class_counts = Counter()
    status_counts = Counter()
    failure_code_counts = Counter()
    total = 0
    failed = 0

    with Path(audit_path).open(
        "w", encoding="utf-8", newline=""
    ) as audit_handle, Path(failure_path).open(
        "w", encoding="utf-8", newline=""
    ) as failure_handle:
        audit_writer = csv.DictWriter(
            audit_handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        failure_writer = csv.DictWriter(
            failure_handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        audit_writer.writeheader()
        failure_writer.writeheader()

        for event_class, paths in (
            ("variant", args.variant_audit),
            ("fusion", args.fusion_audit),
            ("splice", args.splice_audit),
        ):
            for path in paths:
                for source_row in iter_tsv(path):
                    row = {"Class": event_class, **source_row}
                    audit_writer.writerow(row)
                    total += 1
                    class_counts[event_class] += 1
                    status = row.get("Status", "UNSPECIFIED")
                    status_counts[(event_class, status)] += 1
                    if status != "RNA_VALIDATED":
                        failure_writer.writerow(row)
                        failed += 1
                        codes = row.get("Failure codes", "") or "UNSPECIFIED"
                        for code in filter(None, codes.split(";")):
                            failure_code_counts[(event_class, code)] += 1

        for path in args.arriba_discarded:
            for source_row in iter_tsv(path):
                gene1 = source_row.get("#gene1", source_row.get("gene1", ""))
                gene2 = source_row.get("gene2", "")
                filters = source_row.get("filters", "") or "ARRIBA_DISCARDED"
                row = {
                    "Class": "fusion",
                    "Sample": sample_from_path(path),
                    "Event type": "fusion",
                    "Event": f"{gene1}--{gene2}",
                    "Initial finding": (
                        f"Arriba discarded fusion candidate {gene1}--{gene2}"
                    ),
                    "Validation rule": (
                        "Arriba internal filtering and RNA fusion prevalidation"
                    ),
                    "Observed evidence": (
                        f"BREAKPOINTS={source_row.get('breakpoint1', '')}|"
                        f"{source_row.get('breakpoint2', '')};"
                        f"CONFIDENCE={source_row.get('confidence', '')};"
                        f"FILTERS={filters}"
                    ),
                    "Status": "REPORTED_DISCARDED",
                    "Failure codes": filters,
                    "Failure explanation": (
                        "Arriba discarded this initial RNA fusion finding."
                    ),
                    "Required resolution": (
                        "Retained for audit only; excluded from translation "
                        "and proteomic claims."
                    ),
                    "Source": Path(path).name,
                }
                audit_writer.writerow(row)
                failure_writer.writerow(row)
                total += 1
                failed += 1
                class_counts["fusion"] += 1
                status_counts[("fusion", "REPORTED_DISCARDED")] += 1
                for code in filter(None, filters.split(";")):
                    failure_code_counts[("fusion", code)] += 1

    return total, failed, class_counts, status_counts, failure_code_counts


def build_inventory(args, multiqc_rows, multiqc_name):
    inventory = [
        {
            "Software": "FastQC/Trim Galore/STAR/SAMtools/MultiQC",
            "Finding class": "combined QC metrics",
            "Sample": "all",
            "Count": len(multiqc_rows),
            "File": multiqc_name,
        }
    ]
    specifications = [
        ("VEP", "all RNA variant findings", args.vep_vcf, count_records),
        (
            "Arriba",
            "accepted RNA fusion findings",
            args.arriba,
            lambda path: sum(1 for _ in iter_tsv(path)),
        ),
        (
            "Arriba",
            "discarded RNA fusion findings",
            args.arriba_discarded,
            lambda path: sum(1 for _ in iter_tsv(path)),
        ),
        ("StringTie", "assembled transcripts", args.assembled_gtf, count_gtf),
        (
            "GFFCompare",
            "annotated transcripts",
            args.annotated_gtf,
            count_gtf,
        ),
        (
            "GFFCompare",
            "selected novel transcripts",
            args.novel_gtf,
            count_gtf,
        ),
        (
            "pypgatk",
            "validated variant proteins",
            args.variant_fasta,
            count_fasta,
        ),
        (
            "pVACfuse",
            "validated fusion proteins",
            args.fusion_fasta,
            count_fasta,
        ),
        (
            "TransDecoder",
            "validated splice proteins",
            args.splice_fasta,
            count_fasta,
        ),
        (
            "bcftools",
            "progression-specific validated variants",
            args.progression_vcf,
            count_records,
        ),
    ]
    for tool, finding_class, paths, counter in specifications:
        for path in paths:
            inventory.append(
                {
                    "Software": tool,
                    "Finding class": finding_class,
                    "Sample": sample_from_path(path),
                    "Count": counter(path),
                    "File": Path(path).name,
                }
            )
    return inventory


def write_failure_report(path, failed, failure_code_counts, failure_tsv):
    with Path(path).open("w", encoding="utf-8") as handle:
        handle.write("# RNA validation failures\n\n")
        handle.write(f"Rejected or discarded initial findings: {failed}\n\n")
        handle.write(
            "Every failure is retained row-by-row in "
            f"`{Path(failure_tsv).name}` with the initial finding, validation "
            "rule, observed evidence, failure code, explanation and required "
            "resolution.\n\n"
        )
        handle.write("| Class | Failure code | Count |\n")
        handle.write("|---|---|---:|\n")
        for (event_class, code), count in sorted(
            failure_code_counts.items(),
            key=lambda item: (-item[1], item[0]),
        ):
            handle.write(
                f"| {escape_markdown(event_class)} | "
                f"{escape_markdown(code)} | {count} |\n"
            )


def write_variant_explanation_report(path, status_counts, failure_code_counts):
    meanings={
        "VCF_FILTER_NOT_PASS":"Candidate failed upstream VCF filtering.",
        "REFERENCE_ALLELE_MISMATCH":"VCF REF disagrees with GRCh38.",
        "REFERENCE_LOOKUP_FAILED":"GRCh38 sequence could not be retrieved.",
        "RNA_DEPTH_BELOW_MINIMUM":"VCF RNA depth is below the configured minimum.",
        "RNA_ALT_READS_BELOW_MINIMUM":"VCF ALT-supporting reads are below the configured minimum.",
        "RNA_ALT_FRACTION_BELOW_MINIMUM":"VCF ALT fraction is below the configured minimum.",
        "NO_SUPPORTED_PROTEIN_ALTERING_CONSEQUENCE":"VEP reported no supported protein-altering consequence.",
        "MULTIALLELIC_REQUIRES_NORMALIZATION":"Multiallelic record was excluded to prevent incorrect allele assignment.",
    }
    reported=sum(count for (event_class,_status),count in status_counts.items() if event_class=="variant")
    validated=status_counts[("variant","RNA_VALIDATED")]
    with Path(path).open("w",encoding="utf-8") as handle:
        handle.write("# Why RNA variant candidates fail validation\n\n")
        handle.write("The upstream count contains GATK candidates retained by SelectVariants and annotated by VEP. The RNA validator rechecks VCF evidence and GRCh38 consistency; rejection does not mean a second caller proved the site absent. Independent BAM and codon checks are reported separately.\n\n")
        handle.write(f"Reported: {reported}\n\nValidated: {validated}\n\nRejected: {reported-validated}\n\n")
        handle.write("| Failure code | Count | Meaning |\n|---|---:|---|\n")
        for (event_class,code),count in sorted(failure_code_counts.items(),key=lambda item:(-item[1],item[0])):
            if event_class=="variant": handle.write(f"| {escape_markdown(code)} | {count} | {escape_markdown(meanings.get(code,'See row-level observed evidence.'))} |\n")

def write_main_report(path, inventory, class_counts, status_counts, failed):
    with Path(path).open("w", encoding="utf-8") as handle:
        handle.write("# Complete RNA-seq and proteogenomics findings\n\n")
        handle.write(
            "All reported RNA findings remain visible. Only RNA_VALIDATED "
            "events are translated into the proteomics search database. Only "
            "post-search sequence-validated peptide associations are called "
            "proteogenomically validated.\n\n"
        )
        handle.write("## Software findings inventory\n\n")
        handle.write(
            "| Software | Finding class | Sample | Count | File |\n"
            "|---|---|---|---:|---|\n"
        )
        for row in inventory:
            handle.write(
                f"| {escape_markdown(row['Software'])} | "
                f"{escape_markdown(row['Finding class'])} | "
                f"{escape_markdown(row['Sample'])} | {row['Count']} | "
                f"{escape_markdown(row['File'])} |\n"
            )
        handle.write("\n## RNA prevalidation\n\n")
        handle.write(
            "| Class | Reported | RNA validated | Rejected or discarded |\n"
            "|---|---:|---:|---:|\n"
        )
        for event_class in ("variant", "fusion", "splice"):
            reported = class_counts[event_class]
            validated = status_counts[(event_class, "RNA_VALIDATED")]
            handle.write(
                f"| {event_class} | {reported} | {validated} | "
                f"{reported - validated} |\n"
            )
        handle.write("\n## RNA validation failures\n\n")
        handle.write(
            f"Rejected or discarded findings: {failed}. See "
            "`complete_findings.rna_validation_failures.tsv` for the complete "
            "row-level audit and `.md` for failure-code counts.\n\n"
        )
        handle.write("## Proteogenomics\n\n")
        handle.write(
            "Detailed direct-MS/MS, MBR-only, validated and rejected "
            "peptide-event associations are in "
            "`results/proteogenomics_validation/`. RNA-only findings are never "
            "relabelled as proteomically supported.\n"
        )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", required=True)
    parser.add_argument("--multiqc-data", required=True)
    parser.add_argument("--variant-audit", nargs="*", default=[])
    parser.add_argument("--fusion-audit", nargs="*", default=[])
    parser.add_argument("--splice-audit", nargs="*", default=[])
    parser.add_argument("--vep-vcf", nargs="*", default=[])
    parser.add_argument("--arriba", nargs="*", default=[])
    parser.add_argument("--arriba-discarded", nargs="*", default=[])
    parser.add_argument("--assembled-gtf", nargs="*", default=[])
    parser.add_argument("--annotated-gtf", nargs="*", default=[])
    parser.add_argument("--novel-gtf", nargs="*", default=[])
    parser.add_argument("--variant-fasta", nargs="*", default=[])
    parser.add_argument("--fusion-fasta", nargs="*", default=[])
    parser.add_argument("--splice-fasta", nargs="*", default=[])
    parser.add_argument("--progression-vcf", nargs="*", default=[])
    parser.add_argument("--prefix", default="complete_findings")
    args = parser.parse_args()

    prefix = Path(args.prefix)
    audit_path = f"{prefix}.rna_event_audit.tsv"
    failure_path = f"{prefix}.rna_validation_failures.tsv"
    failure_md = f"{prefix}.rna_validation_failures.md"
    inventory_path = f"{prefix}.software_inventory.tsv"
    multiqc_path = f"{prefix}.multiqc_general_stats.tsv"
    report_path = f"{prefix}.report.md"
    variant_explanation_path = f"{prefix}.rna_variant_validation_explanations.md"

    multiqc_rows, multiqc_name = read_multiqc_general_stats(args.multiqc_data)
    multiqc_fields = list(multiqc_rows[0]) if multiqc_rows else ["Sample"]
    write_tsv(multiqc_path, multiqc_rows, multiqc_fields)

    total, failed, class_counts, status_counts, failure_code_counts = (
        stream_audits(args, audit_path, failure_path)
    )
    inventory = build_inventory(args, multiqc_rows, multiqc_name)
    write_tsv(
        inventory_path,
        inventory,
        ["Software", "Finding class", "Sample", "Count", "File"],
    )
    write_failure_report(
        failure_md,
        failed,
        failure_code_counts,
        failure_path,
    )
    write_variant_explanation_report(variant_explanation_path, status_counts, failure_code_counts)
    write_main_report(
        report_path,
        inventory,
        class_counts,
        status_counts,
        failed,
    )

    print(f"Wrote {report_path}")
    print(f"Streamed RNA audit rows: {total}")
    print(f"Failure rows: {failed}")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
