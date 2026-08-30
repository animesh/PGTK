#!/usr/bin/env python3
import argparse
import csv
import html
import json
from pathlib import Path
from report_legend import HTML_LEGEND


def rows(path):
    source = Path(path)
    if not source.exists() or source.stat().st_size == 0:
        return []
    with source.open(encoding="utf-8", errors="replace", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def number(value):
    try:
        return float(str(value).replace(",", ""))
    except (TypeError, ValueError):
        return 0.0


def write_json(output_dir, identifier, title, description, data, plot_type="bargraph", pconfig=None, headers=None):
    if not data:
        return
    payload = {
        "id": identifier,
        "section_name": title,
        "description": "Colors distinguish categories only; color intensity is not confidence. " + description,
        "plot_type": plot_type,
        "pconfig": {"id": identifier, "title": title, **(pconfig or {})},
        "data": data,
    }
    if headers:
        payload["headers"] = headers
    (output_dir / f"{identifier}_mqc.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")


def html_section(identifier, title, body):
    return f"---\nid: {identifier}\nsection_name: {title}\n---\n{body}\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--complete-report", required=True)
    parser.add_argument("--rna-failure-report", required=True)
    parser.add_argument("--rna-variant-explanations", required=True)
    parser.add_argument("--comparative-report", required=True)
    parser.add_argument("--progression-report", required=True)
    parser.add_argument("--variant-inventory", required=True)
    parser.add_argument("--fasta-inventory", required=True)
    parser.add_argument("--rna-inventory", required=True)
    parser.add_argument("--external-comparison", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--progression-summary", required=True)
    parser.add_argument("--progression-enrichment", required=True)
    parser.add_argument("--progression-pairwise", required=True)
    parser.add_argument("--variant-landscape-summary", required=True)
    parser.add_argument("--variant-landscape-genes", required=True)
    parser.add_argument("--pairwise-alleles", required=True)
    parser.add_argument("--pairwise-genes", required=True)
    parser.add_argument("--candidate-priority", required=True)
    parser.add_argument("--maxquant-enabled", choices=["true", "false"], default="false")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    variant_data = {}
    for row in rows(args.variant_inventory):
        sample = row.get("Sample") or row.get("sample") or row.get("Name") or "unknown"
        stage = row.get("Stage") or row.get("stage") or row.get("Category") or "value"
        count = next((number(row[key]) for key in row if key.lower() in {"records", "count", "variants", "alleles"}), 0)
        variant_data.setdefault(sample, {})[stage] = count
    write_json(
        output_dir,
        "pgtk_variant_stages",
        "Variant stages",
        "Counts retained across calling, filtering, annotation and RNA-validation stages.",
        variant_data,
        pconfig={"ylab": "Records"},
    )

    rna_data = {}
    for row in rows(args.rna_inventory):
        sample = row.get("Sample") or row.get("sample") or "all"
        category = row.get("Event class") or row.get("Category") or row.get("EventType") or row.get("Type") or "events"
        count = next((number(row[key]) for key in row if key.lower() in {"data rows", "records", "count", "events", "findings"}), 0)
        rna_data.setdefault(sample, {})[category] = count
    write_json(
        output_dir,
        "pgtk_rna_evidence",
        "RNA evidence",
        "RNA variants, fusions and validated splice evidence by sample. Event classes use different counting units.",
        rna_data,
        pconfig={"ylab": "Rows"},
    )

    progression_data = {}
    for row in rows(args.progression_summary):
        sample = row.get("Sample") or row.get("sample") or row.get("Set") or "progression"
        for key, value in row.items():
            if key.lower() in {"sample", "subject", "group", "set", "comparison"} or not str(value).strip():
                continue
            try:
                progression_data.setdefault(sample, {})[key] = float(value)
            except (TypeError, ValueError):
                pass
    write_json(
        output_dir,
        "pgtk_progression_summary",
        "Progression summary",
        "Baseline-subtracted RNA evidence. Absence from baseline RNA does not establish DNA-level absence.",
        progression_data,
        pconfig={"ylab": "Count"},
    )

    significant_go = []
    for row in rows(args.progression_enrichment):
        fdr = number(row.get("FDR", 1))
        overlap = number(row.get("ProgressionGenesInTerm") or row.get("Overlap") or 0)
        if fdr <= 0.1 and overlap > 0:
            significant_go.append((fdr, -overlap, row))
    significant_go.sort(key=lambda item: (item[0], item[1], item[2].get("GO_ID", "")))
    go_table = {}
    for index, (_, _, row) in enumerate(significant_go[:20], 1):
        key = f"{row.get('Sample', '')} | {row.get('GO_Name', '')} | {index}"
        go_table[key] = {
            "GO ID": row.get("GO_ID", ""),
            "FDR": number(row.get("FDR")),
            "Odds ratio": number(row.get("OddsRatio")),
            "Genes": number(row.get("ProgressionGenesInTerm") or row.get("Overlap")),
        }
    write_json(
        output_dir,
        "pgtk_progression_go",
        "Progression GO enrichment",
        "Top progression-associated terms passing FDR 0.1. Full tables remain linked below.",
        go_table,
        plot_type="table",
    )

    landscape_rows = rows(args.variant_landscape_summary)
    stage_labels = {
        "raw_genotyped": "Raw genotyped",
        "normalized": "Normalized",
        "hard_filtered_all": "Hard filtered",
        "hard_filter_pass": "PASS",
        "vep_pass": "VEP PASS",
        "rna_validated": "RNA validated",
        "progression_nonbaseline_only": "Progression only",
        "progression_baseline_only": "Baseline only",
        "progression_shared_with_baseline": "Shared with baseline",
    }

    def landscape_plot(identifier, title, description, category, metrics, ylab="Alleles"):
        data = {}
        for row in landscape_rows:
            if row.get("Category") != category or row.get("Metric") not in metrics:
                continue
            key = f"{row.get('Sample', '')} | {stage_labels.get(row.get('Stage', ''), row.get('Stage', ''))}"
            data.setdefault(key, {})[row.get("Metric", "")] = number(row.get("Count"))
        write_json(output_dir, identifier, title, description, data, pconfig={"ylab": ylab})

    landscape_plot(
        "pgtk_variant_types", "Variant types by sample and stage",
        "SNVs, MNVs, insertions, deletions and complex alleles across raw, filtered, RNA-validated and progression stages.",
        "variant_type", {"SNV", "MNV", "insertion", "deletion", "complex_allele"}
    )
    landscape_plot(
        "pgtk_variant_genotypes", "Variant genotypes by sample and stage",
        "Heterozygous, homozygous-alternate and missing genotypes across variant-processing stages.",
        "variant_type", {"heterozygous", "homozygous_alt", "missing_genotype"}
    )
    landscape_plot(
        "pgtk_substitution_classes", "Transition and transversion counts",
        "Substitution classes by sample and stage. These counts are meaningful for SNV-containing stages.",
        "variant_type", {"transition", "transversion"}
    )
    landscape_plot(
        "pgtk_variant_impact", "Predicted variant impact",
        "VEP HIGH, MODERATE, LOW and MODIFIER transcript consequences by sample and stage.",
        "VEP_impact", {"HIGH", "MODERATE", "LOW", "MODIFIER", "unannotated"}, "Annotations"
    )

    consequence_totals = {}
    for row in landscape_rows:
        if row.get("Category") != "VEP_consequence" or row.get("Stage") not in {"vep_pass", "rna_validated", "progression_nonbaseline_only"}:
            continue
        key = f"{row.get('Sample', '')} | {stage_labels.get(row.get('Stage', ''), row.get('Stage', ''))}"
        consequence_totals.setdefault(key, {})[row.get("Metric", "")] = number(row.get("Count"))
    consequence_top = {}
    for key, values in consequence_totals.items():
        consequence_top[key] = dict(sorted(values.items(), key=lambda item: (-item[1], item[0]))[:15])
    write_json(
        output_dir, "pgtk_variant_consequences", "Top functional consequences",
        "Top 15 VEP consequence classes per sample and selected stage. Full consequence counts remain in variant_landscape.summary.tsv.",
        consequence_top, pconfig={"ylab": "Annotations"}
    )

    pairwise = rows(args.pairwise_alleles)
    pair_counts = {}
    chromosome_counts = {}
    shared_gene_counts = {}
    shared_high = []
    chromosome_order = [str(value) for value in range(1, 23)] + ["X", "Y", "MT", "Other"]
    for row in pairwise:
        subject = row.get("Subject", "")
        comparison = f"{row.get('SampleA', '')} vs {row.get('SampleB', '')}"
        contrast = row.get("ContrastClass", "unclassified")
        key = f"{subject} | {comparison}"
        pair_counts.setdefault(key, {})[contrast] = pair_counts.setdefault(key, {}).get(contrast, 0) + 1
        if contrast == "shared":
            chromosome = str(row.get("Chrom", "")).removeprefix("chr")
            if chromosome not in chromosome_order:
                chromosome = "Other"
            chromosome_counts.setdefault(key, {})[chromosome] = chromosome_counts.setdefault(key, {}).get(chromosome, 0) + 1
            gene = row.get("Gene", "").strip()
            if gene:
                shared_gene_counts[gene] = shared_gene_counts.get(gene, 0) + 1
            if row.get("Impact") == "HIGH":
                shared_high.append(row)
    write_json(
        output_dir, "pgtk_progression_pair_overlap", "Progression overlap",
        "Alleles shared between progression samples, or exclusive to either progression sample, after baseline subtraction.",
        pair_counts, pconfig={"ylab": "Alleles"}
    )
    write_json(
        output_dir, "pgtk_shared_progression_chromosomes", "Shared progression variants by chromosome",
        "Chromosomal distribution of alleles present in both progression samples and absent from the baseline callset. RNA callable-space differences are not normalized.",
        chromosome_counts, pconfig={"ylab": "Shared alleles", "categories": chromosome_order}
    )
    top_shared_genes = dict(sorted(shared_gene_counts.items(), key=lambda item: (-item[1], item[0]))[:30])
    write_json(
        output_dir, "pgtk_shared_progression_genes", "Genes shared by progression samples",
        "Top genes carrying alleles found in both progression samples but absent from the baseline callset. Full alleles are linked from the report catalogue.",
        {"Shared progression": top_shared_genes}, pconfig={"ylab": "Alleles"}
    )
    shared_high.sort(key=lambda row: (row.get("Gene", ""), row.get("Chrom", ""), number(row.get("Pos"))))
    high_table = {}
    for index, row in enumerate(shared_high[:100], 1):
        high_table[str(index)] = {
            "Gene": row.get("Gene", ""), "Chromosome": row.get("Chrom", ""),
            "Position": row.get("Pos", ""), "Ref": row.get("Ref", ""), "Alt": row.get("Alt", ""),
            "Consequence": row.get("Consequence", ""), "Impact": row.get("Impact", ""),
        }
    write_json(
        output_dir, "pgtk_shared_progression_high_impact", "Shared progression HIGH-impact candidates",
        "HIGH-impact alleles observed in both progression samples and absent from the baseline RNA callset. These are candidates, not DNA-confirmed acquired mutations.",
        high_table, plot_type="table"
    )

    priority_rows = rows(args.candidate_priority)
    priority_table = {}
    for index, row in enumerate(priority_rows[:100], 1):
        priority_table[str(index)] = row
    write_json(
        output_dir, "pgtk_progression_priority", "Progression candidate priorities",
        "Top progression candidates from the complete candidate-priority table.",
        priority_table, plot_type="table"
    )

    external_rows = [
        row for row in rows(args.external_comparison)
        if any(str(value).strip() not in {"", "0", "NA", "N/A", "not_run", "disabled"} for value in row.values())
    ]
    if external_rows:
        external_table = {str(index): row for index, row in enumerate(external_rows[:100], 1)}
        write_json(
            output_dir,
            "pgtk_external_validation",
            "External caller validation",
            "Overlap and concordance with the configured external caller.",
            external_table,
            plot_type="table",
        )

    reports = [
        ("Complete findings", "../reports/complete_findings.report.md"),
        ("RNA validation failures", "../reports/complete_findings.rna_validation_failures.md"),
        ("RNA variant explanations", "../reports/complete_findings.rna_variant_validation_explanations.md"),
        ("Comparative biological evidence", "../comparative_advantage/comparative_advantage.report.md"),
        ("Progression biology", "../progression_biology/progression_biology.report.md"),
        ("Expression and GO", "../expression/go/"),
        ("Variant landscape tables and methods", "../variant_landscape/"),
        ("Shared and exclusive progression alleles", "../progression_biology/progression_biology.pairwise_allele_contrasts.tsv"),
        ("Shared and exclusive progression genes", "../progression_biology/progression_biology.pairwise_gene_contrasts.tsv"),
        ("Progression candidate priorities", "../progression_biology/progression_biology.candidate_priority.tsv"),
        ("Finding explorer", "../igv/findings/finding_explorer/index.html"),
    ]
    if external_rows:
        reports.append(("External caller comparison", "../comparison/external_vcf/"))
    if args.maxquant_enabled == "true":
        reports.append(("MaxQuant validation", "../proteogenomics_validation/"))

    body = (HTML_LEGEND + 
        "<div class=\"alert alert-info\"><strong>Interpretation:</strong> RNA-observed events are exploratory evidence and are not DNA-confirmed somatic mutations.</div>"
        "<p>This lightweight dashboard contains summaries only. Full reports remain separate files.</p><ul>"
        + "".join(f'<li><a href="{link}">{html.escape(title)}</a></li>' for title, link in reports)
        + "</ul><p><a href=\"../qc/multiqc_report.html\">Open the detailed sequencing-QC report</a>.</p>"
    )
    (output_dir / "pgtk_results_overview_mqc.html").write_text(
        html_section("pgtk_results_overview", "PGTK results overview", body), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
