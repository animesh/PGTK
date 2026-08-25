#!/usr/bin/env python3
import argparse
import csv
import html
import json
import math
import re
from collections import defaultdict
from pathlib import Path


def rows(path):
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def number(value, default=0.0):
    try:
        parsed = float(str(value).replace(",", ""))
        return parsed if math.isfinite(parsed) else default
    except (TypeError, ValueError):
        return default


def pick(row, *names, default=""):
    lowered = {key.lower(): value for key, value in row.items()}
    for name in names:
        value = lowered.get(name.lower())
        if value is not None and str(value).strip():
            return value
    return default


def emit(output_dir, identifier, title, data, description, ylab, ymax=None):
    if not data:
        return
    config = {"id": identifier, "title": title, "ylab": ylab}
    if ymax is not None:
        config["ymax"] = ymax
    payload = {
        "id": identifier,
        "section_name": title,
        "description": description,
        "plot_type": "bargraph",
        "pconfig": config,
        "data": data,
    }
    (output_dir / f"{identifier}_mqc.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )


def metrics(path):
    extracted = {}
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        match = re.match(
            r"\s*[-*]?\s*([^:\t|]{2,100})\s*[:\t|]\s*([0-9][0-9,]*(?:\.\d+)?)\s*$",
            line,
        )
        if match:
            extracted[match.group(1).strip()] = number(match.group(2))
    return extracted


def markdown_table_value(path, row_name, column_name):
    table_rows = []
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        if line.strip().startswith("|"):
            table_rows.append([cell.strip() for cell in line.strip().strip("|").split("|")])
    if not table_rows:
        raise SystemExit(f"no markdown table found in {path}")
    header = table_rows[0]
    if column_name not in header:
        raise SystemExit(f"column {column_name!r} missing from markdown table in {path}")
    index = header.index(column_name)
    for values in table_rows[2:]:
        if values and values[0] == row_name:
            return number(values[index])
    raise SystemExit(f"row {row_name!r} missing from markdown table in {path}")


def require_metric(metric_map, name):
    if name not in metric_map:
        raise SystemExit(f"required compact-dashboard metric missing: {name}")
    return metric_map[name]


def nonredundant_go(source_rows, mode, maximum=10):
    selected = []
    seen_names = set()
    for row in sorted(
        source_rows,
        key=lambda item: (number(item.get("FDR"), 1.0), item.get("GO_ID", "")),
    ):
        fdr = number(row.get("FDR"), 1.0)
        if fdr > 0.1:
            continue
        if mode == "ranked":
            if abs(number(row.get("MeanScore"))) < 0.25 or abs(number(row.get("ZScore"))) < 2.0:
                continue
        else:
            if number(row.get("OddsRatio")) < 1.5 or number(row.get("ForegroundGenesInTerm")) < 5:
                continue
        term = (row.get("GO_Name") or row.get("GO_ID") or "term").strip()
        normalized = re.sub(r"[^a-z0-9]+", " ", term.lower()).strip()
        if normalized in seen_names:
            continue
        seen_names.add(normalized)
        selected.append((fdr, row, term))
        if len(selected) == maximum:
            break
    return selected


def build_go_plot(output_dir, identifier, title, path, mode):
    grouped = defaultdict(list)
    for row in rows(path):
        grouped[pick(row, "Analysis", "Sample", "Comparison", default="analysis")].append(row)
    data = {}
    for analysis in sorted(grouped):
        for fdr, row, term in nonredundant_go(grouped[analysis], mode):
            label = f"{analysis} | {term[:70]}"
            data[label] = {"-log10 FDR": min(50.0, -math.log10(max(fdr, 1e-50)))}
    description = (
        "Top nonredundant terms after dashboard filters: FDR <= 0.1, absolute mean "
        "log2-TPM score >= 0.25 and absolute Z score >= 2. Full ranked table is linked."
        if mode == "ranked"
        else "Top nonredundant terms after dashboard filters: FDR <= 0.1, odds ratio >= 1.5 "
        "and at least 5 foreground genes. Full enrichment table is linked."
    )
    emit(output_dir, identifier, title, data, description, "-log10 FDR")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", required=True)
    for name in [
        "variant-inventory", "rna-inventory", "progression-summary", "external-comparison",
        "expression-ora", "ranked-go", "variant-set-go", "proteogenomics-summary",
        "integrated-report", "evidence-classification", "read-summary", "codon-summary", "provenance-summary",
    ]:
        parser.add_argument(f"--{name}", required=True)
    parser.add_argument("--maxquant-enabled", choices=["true", "false"], required=True)
    args = parser.parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    stages = {}
    for row in rows(args.variant_inventory):
        sample = pick(row, "Sample", "Name", default="unknown")
        stage = pick(row, "Stage", "Category", default="variants")
        stages.setdefault(sample, {})[stage] = number(
            pick(row, "Alleles", "Records", "Count", "Variants")
        )
    emit(output_dir, "pgtk_variant_attrition", "Variant attrition", stages,
         "Allele counts retained from raw calls through PASS, RNA validation and progression subtraction.",
         "Alleles")

    evidence = {}
    for row in rows(args.rna_inventory):
        sample = pick(row, "Sample", default="all")
        event_class = pick(row, "Event class", "Category", "EventType", "Type", default="events")
        count = number(pick(row, "Data rows", "Records", "Count", "Events", "Findings"))
        evidence.setdefault(sample, {})[event_class] = count
    if not evidence or not any(value > 0 for sample in evidence.values() for value in sample.values()):
        raise SystemExit("RNA evidence inventory produced no positive counts")
    emit(output_dir, "pgtk_rna_evidence", "RNA evidence", evidence,
         "Validated fusion rows and validated-splice audit rows by sample. Different event classes use different counting units.",
         "Rows")

    progression = {}
    for row in rows(args.progression_summary):
        sample = pick(row, "Sample", "Set", "Comparison", default="progression")
        for key, value in row.items():
            if key.lower() not in {"sample", "subject", "group", "set", "comparison"}:
                try:
                    progression.setdefault(sample, {})[key] = float(value)
                except (TypeError, ValueError):
                    pass
    emit(output_dir, "pgtk_progression", "Progression evidence", progression,
         "Nonbaseline-only progression evidence. RNA absence at baseline is not DNA confirmation.", "Count")

    build_go_plot(output_dir, "pgtk_expression_go", "Expression GO enrichment",
                  args.expression_ora, "ora")
    build_go_plot(output_dir, "pgtk_ranked_go", "Ranked expression GO",
                  args.ranked_go, "ranked")
    build_go_plot(output_dir, "pgtk_variant_set_go", "Progression variant-set GO",
                  args.variant_set_go, "ora")

    overlap, concordance = {}, {}
    for row in rows(args.external_comparison):
        label = f"{pick(row, 'Sample', default='unknown')} | {pick(row, 'Stage', default='unknown')}"
        metric = pick(row, "Metric").lower()
        value = number(pick(row, "Value"))
        if "overlap" in metric and ("percent" in metric or "%" in metric):
            overlap.setdefault(label, {})["PGTK overlap %"] = value
        if "concordance" in metric and ("percent" in metric or "%" in metric):
            concordance.setdefault(label, {})["Genotype concordance %"] = value
    if overlap or concordance:
        if len(overlap) != 9 or len(concordance) != 9:
            raise SystemExit(f"expected 9 external-caller stage comparisons, found overlap={len(overlap)}, concordance={len(concordance)}")
        emit(output_dir, "pgtk_external_overlap", "External caller overlap", overlap,
             "Percentage of each PGTK stage present in the configured external callset.", "Percent", 100)
        emit(output_dir, "pgtk_external_concordance", "External caller genotype concordance", concordance,
             "Genotype concordance among alleles shared with the configured external caller.", "Percent", 100)

    proteogenomics = metrics(args.proteogenomics_summary)
    integrated = metrics(args.integrated_report)
    funnel = {
        "Protein-altering events evaluated": require_metric(proteogenomics, "Protein-altering sample-specific genomic variants"),
        "With mapped variant-protein peptide": require_metric(proteogenomics, "Variants with any mapped variant-protein peptide"),
        "Altered-residue events": require_metric(proteogenomics, "Variants with altered-residue peptide association"),
        "Search-consistent events": require_metric(proteogenomics, "Search-consistent variant events"),
        "Sample-matched direct MS/MS events": markdown_table_value(args.evidence_classification, "SAMPLE_MATCHED_DIRECT_MSMS", "Variant events"),
        "Absent from canonical and Ensembl": require_metric(proteogenomics, "Canonical-and-reference-absent variant events"),
        "Strict integrated events": require_metric(integrated, "Strict integrated events"),
    }
    emit(output_dir, "pgtk_maxquant_evidence", "MaxQuant evidence funnel",
         {"Evidence": funnel},
         "Explicit event-level evidence funnel. Stages are selected by metric name and never by report-text order.",
         "Events")

    read_metrics = metrics(args.read_summary)
    read_summary = {
        key: require_metric(read_metrics, key)
        for key in [
            "Altered-residue genomic events", "Read-level sample-event comparisons",
            "Accepted Arriba events", "Translated junction findings", "SNV events",
            "complex_allele events", "insertion events",
        ]
    }
    emit(output_dir, "pgtk_read_validation", "Read-level validation",
         {"Observed": read_summary}, "Compact read-evidence event summary.", "Count")

    codon = metrics(args.codon_summary)
    provenance = metrics(args.provenance_summary)
    independent = dict(list(codon.items())[:8] + list(provenance.items())[:8])
    emit(output_dir, "pgtk_independent_validation", "Independent validation",
         {"Validated": independent}, "Codon and read-provenance validation summary.", "Count")

    links = [
        ("Variant landscape and nonsynonymous GO", "../variant_landscape/"),
        ("Complete findings", "../reports/complete_findings.report.md"),
        ("Expression and GO", "../expression/go/"),
        ("Progression biology", "../progression_biology/"),
        ("Sarek comparison", "../comparison/external_vcf/"),
        ("MaxQuant validation", "../proteogenomics_validation/"),
        ("Offline variant explorer", "../igv/findings/finding_explorer/index.html"),
    ]
    guide = """<h3>PGTK Results Guide and Navigation</h3><p><b>Start here.</b> Raw, normalized, hard-filtered, PASS, RNA-validated and progression calls are distinct stages.</p><ul><li><b>RNA validated:</b> RNA-supported and protein-altering by VEP, not DNA-confirmed.</li><li><b>Progression nonbaseline-only:</b> absent from baseline RNA callset, not proof of DNA acquisition.</li><li><b>Nonsynonymous GO:</b> over-representation among unique genes with protein-altering VEP consequences.</li><li><b>Offline explorer:</b> direct-open compact finding metadata and precomputed evidence counts; server mode is optional for full IGV Reports.</li></ul><h4>Variant landscape files</h4><ul><li><code>variant_landscape.summary.tsv</code>: sample-stage variant counts.</li><li><code>variant_landscape.nonsynonymous_genes.tsv</code>: unique protein-altering genes.</li><li><code>variant_landscape.go_significant.tsv</code>: FDR-significant terms only.</li><li><code>variant_landscape.go_top.tsv</code>: top 100 ranked terms per sample-stage.</li><li><code>variant_landscape.go_summary.tsv</code>: tested and significant term counts.</li></ul>"""
    body = '<div class="alert alert-info">' + guide + '</div><h4>Open detailed outputs</h4><ul>' + "".join(f'<li><a href="{url}">{html.escape(label)}</a></li>' for label, url in links) + "</ul>"
    (output_dir / "pgtk_results_guide_mqc.html").write_text("---\nid: pgtk_results_guide\nsection_name: PGTK results overview\n---\n" + body, encoding="utf-8")


if __name__ == "__main__":
    main()
