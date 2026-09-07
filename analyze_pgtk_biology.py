#!/usr/bin/env python3
"""Read-only biological and hypothesis review of completed PGTK results.

The program preserves the original biological-review outputs and adds a focused,
offline assessment of five progression and ATX-101 hypotheses using only local
PGTK results. It never modifies pipeline outputs and requires only Python's
standard library.
"""

import argparse
import csv
import gzip
import json
import math
import re
import statistics
from collections import Counter, defaultdict
from pathlib import Path

NULLS = {"", ".", "NA", "N/A", "NULL", "NONE", "NAN"}
IMPACT_RANK = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1, "": 0}
DIRECT_STATUSES = {"ALT_SUPPORTED", "MIXED_ALT_AND_REFERENCE"}
STRUCTURAL_TYPES = {"FUSION", "SPLICE_JUNCTION", "CONTEXT_EVENT"}

# Fixed, offline gene sets. These are intentionally transparent and editable.
GENE_SETS = {
    "translation_ribosome": {
        "RPL3", "RPL4", "RPL5", "RPL6", "RPL7", "RPL7A", "RPL8", "RPL9", "RPL10", "RPL11", "RPL12", "RPL13", "RPL13A", "RPL14", "RPL15", "RPL17", "RPL18", "RPL18A", "RPL19", "RPL21", "RPL22", "RPL23", "RPL23A", "RPL24", "RPL26", "RPL27", "RPL27A", "RPL28", "RPL29", "RPL30", "RPL31", "RPL32", "RPL34", "RPL35", "RPL35A", "RPL36", "RPL36A", "RPL37", "RPL37A", "RPL38", "RPL39", "RPLP0", "RPLP1", "RPLP2", "RPS2", "RPS3", "RPS3A", "RPS4X", "RPS5", "RPS6", "RPS7", "RPS8", "RPS9", "RPS10", "RPS11", "RPS12", "RPS13", "RPS14", "RPS15", "RPS15A", "RPS16", "RPS17", "RPS18", "RPS19", "RPS20", "RPS21", "RPS23", "RPS24", "RPS25", "RPS26", "RPS27", "RPS27A", "RPS28", "RPS29", "EEF1A1", "EEF1B2", "EEF1D", "EEF1G", "EEF2", "EIF2S1", "EIF2S2", "EIF2S3", "EIF3A", "EIF3B", "EIF3C", "EIF3D", "EIF3E", "EIF3F", "EIF3G", "EIF3H", "EIF3I", "EIF3J", "EIF4A1", "EIF4E", "EIF4G1"
    },
    "proteostasis_er_upr": {
        "XBP1", "ATF4", "ATF6", "ERN1", "EIF2AK3", "HSPA5", "HSP90B1", "CALR", "CANX", "PDIA3", "PDIA4", "PDIA6", "P4HB", "DNAJB9", "DNAJC3", "DDIT3", "HERPUD1", "SEL1L", "SYVN1", "EDEM1", "EDEM2", "EDEM3", "DERL1", "DERL2", "DERL3", "VCP", "NPLOC4", "UFD1", "SEC61A1", "SEC61B", "SEC61G", "SEC62", "SEC63", "RRBP1", "ERO1A", "ERO1B", "CRELD2"
    },
    "proteasome_ubiquitin": {
        "PSMA1", "PSMA2", "PSMA3", "PSMA4", "PSMA5", "PSMA6", "PSMA7", "PSMB1", "PSMB2", "PSMB3", "PSMB4", "PSMB5", "PSMB6", "PSMB7", "PSMB8", "PSMB9", "PSMB10", "PSMC1", "PSMC2", "PSMC3", "PSMC4", "PSMC5", "PSMC6", "PSMD1", "PSMD2", "PSMD3", "PSMD4", "PSMD5", "PSMD6", "PSMD7", "PSMD8", "PSMD9", "PSMD10", "PSMD11", "PSMD12", "PSMD13", "PSMD14", "PSME1", "PSME2", "PSME3", "UBB", "UBC", "UBE2D1", "UBE2D2", "UBE2D3", "UBE2S", "UBE2I", "UBE2O"
    },
    "rna_processing_surveillance": {
        "UPF1", "UPF2", "UPF3A", "UPF3B", "SMG1", "SMG5", "SMG6", "SMG7", "SMG8", "SMG9", "RBM8A", "MAGOH", "EIF4A3", "CASC3", "RNPS1", "SRSF1", "SRSF2", "SRSF3", "SRSF4", "SRSF5", "SRSF6", "SRSF7", "SRSF9", "SF3B1", "SF3B2", "SF3B3", "U2AF1", "U2AF2", "PRPF8", "PRPF19", "PRPF38B", "DDX5", "DDX17", "DDX20", "DDX39B", "HNRNPA1", "HNRNPC", "HNRNPK", "RBFOX2", "SON", "ACIN1", "THOC5", "NXF1"
    },
    "surface_adhesion_trafficking": {
        "TNFRSF17", "SDC1", "CD38", "CD44", "ITGA4", "ITGB1", "ITGB4", "ITGB7", "CXCR4", "NCAM1", "ICAM1", "MICA", "HLA-A", "HLA-B", "HLA-C", "HLA-G", "LILRB1", "LILRB4", "TSG101", "COG1", "COG3", "COG4", "SEC23A", "SEC24A", "SEC31A", "SEC31B", "GOLGA3", "GOLGA5", "STX8", "VPS33B", "VPS35", "HGS", "RAB5A", "RAB7A", "RAB11A", "B4GALT1", "B4GALNT1", "ST3GAL4", "ST6GAL1", "FUT2", "MAN2A2", "MAN2B1", "POMGNT1", "ALG1", "ALG6", "ALG12"
    },
    "glycolysis_ppp": {
        "SLC2A1", "SLC2A3", "HK1", "HK2", "GPI", "PFKM", "PFKP", "PFKL", "ALDOA", "ALDOB", "TPI1", "GAPDH", "PGK1", "PGAM1", "ENO1", "ENO2", "ENO3", "PKM", "LDHA", "LDHB", "G6PD", "PGLS", "PGD", "RPE", "RPIA", "TKT", "TALDO1"
    },
    "mitochondrial_oxphos_tca": {
        "CS", "ACO2", "IDH2", "IDH3A", "IDH3B", "IDH3G", "OGDH", "DLST", "SUCLA2", "SUCLG1", "SUCLG2", "SDHA", "SDHB", "SDHC", "SDHD", "FH", "MDH2", "NDUFS1", "NDUFS2", "NDUFS3", "NDUFS7", "NDUFV1", "NDUFV2", "UQCRC1", "UQCRC2", "COX4I1", "COX5A", "COX5B", "ATP5F1A", "ATP5F1B", "ATP5PD", "ATP5PO"
    },
    "redox_nad_glutathione": {
        "NADK", "NADK2", "NADSYN1", "NAMPT", "NMNAT1", "NMNAT2", "NMNAT3", "NMRK1", "QPRT", "NAPRT", "GCLC", "GCLM", "GSS", "GSR", "GPX1", "GPX2", "GPX4", "PRDX1", "PRDX2", "PRDX3", "PRDX5", "TXN", "TXNRD1", "TXNRD2", "GLRX", "GLRX3", "SLC7A11", "SLC25A39", "GSTM1", "GSTP1", "GSTO1", "GSTZ1"
    },
    "mitophagy_autophagy": {
        "PINK1", "PRKN", "BNIP3", "BNIP3L", "FUNDC1", "ULK1", "ULK2", "ATG2A", "ATG2B", "ATG5", "ATG7", "ATG12", "ATG16L1", "BECN1", "BCL2", "MAP1LC3A", "MAP1LC3B", "SQSTM1", "OPTN", "NDP52", "TFEB", "LAMP1", "LAMP2", "VPS11", "UVRAG"
    },
    "proliferation_replication_stress": {
        "PCNA", "MKI67", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "CHEK1", "ATR", "ATM", "BRCA1", "BRCA2", "FANCA", "FANCI", "FANCM", "RAD51", "RPA1", "RPA2", "RPA3", "DTL", "CDT1", "TICRR", "AURKA", "AURKB", "BIRC5", "CENPF", "CEP55"
    },
}

HYPOTHESES = {
    "H1_translation_proteostasis": {
        "title": "Progression increases translation and proteostasis load",
        "positive_sets": ["translation_ribosome", "proteostasis_er_upr", "proteasome_ubiquitin"],
        "negative_sets": [],
        "go_keywords": ["translation", "ribosome", "protein folding", "unfolded protein", "endoplasmic reticulum", "proteasome", "ubiquitin", "er-associated degradation"],
    },
    "H2_rna_processing_quality_control": {
        "title": "Progression remodels RNA processing and quality control",
        "positive_sets": ["rna_processing_surveillance"],
        "negative_sets": [],
        "go_keywords": ["rna splicing", "spliceosome", "rna transport", "rna surveillance", "nonsense-mediated decay", "mrna processing"],
    },
    "H3_surface_adhesion_remodeling": {
        "title": "Progression remodels surface, adhesion, glycosylation and trafficking",
        "positive_sets": [],
        "negative_sets": ["surface_adhesion_trafficking"],
        "go_keywords": ["cell surface", "cell adhesion", "integrin", "extracellular matrix", "glycosylation", "golgi", "vesicle", "membrane trafficking"],
    },
    "H4_metabolic_redox_stress": {
        "title": "Progression creates a stressed metabolic and redox state",
        "positive_sets": ["glycolysis_ppp", "redox_nad_glutathione", "proliferation_replication_stress"],
        "negative_sets": ["mitochondrial_oxphos_tca", "mitophagy_autophagy"],
        "go_keywords": ["glycolysis", "pentose phosphate", "nad", "glutathione", "oxidative phosphorylation", "tricarboxylic acid", "mitochond", "redox", "reactive oxygen", "mitophagy"],
    },
    "H5_composite_atx101_vulnerability": {
        "title": "ATX-101 sensitivity reflects a composite stress state rather than one mutation",
        "positive_sets": ["translation_ribosome", "proteostasis_er_upr", "proteasome_ubiquitin", "rna_processing_surveillance", "glycolysis_ppp", "redox_nad_glutathione", "proliferation_replication_stress"],
        "negative_sets": ["mitophagy_autophagy"],
        "go_keywords": ["replication stress", "dna repair", "protein folding", "proteasome", "ribosome", "rna surveillance", "glycolysis", "redox"],
    },
}

ATX_CANDIDATE_GENES = {
    "PCNA", "ENO1", "PGD", "GAPDH", "PFKP", "LDHA", "UPF3B", "TNFRSF17", "LILRB4", "TSG101", "COG4"
}


def open_text(path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if str(path).endswith(".gz") else path.open("r", encoding="utf-8", errors="replace")


def clean(value):
    return "" if value is None else str(value).strip()


def as_int(value, default=0):
    try:
        return int(float(clean(value)))
    except (TypeError, ValueError):
        return default


def as_float(value):
    text = clean(value)
    if text.upper() in NULLS:
        return None
    try:
        number = float(text)
        return number if math.isfinite(number) else None
    except ValueError:
        return None


def first(row, names, default=""):
    lower = {str(k).lower(): v for k, v in row.items()}
    for name in names:
        if name in row and clean(row[name]):
            return row[name]
        value = lower.get(name.lower())
        if clean(value):
            return value
    return default


def read_table(path):
    with open_text(path) as handle:
        sample = handle.read(8192)
        handle.seek(0)
        first_line = sample.splitlines()[0] if sample.splitlines() else ""
        delimiter = "\t" if "\t" in first_line else ","
        yield from csv.DictReader(handle, delimiter=delimiter)


def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sample_from_name(path):
    return path.name.split(".", 1)[0]


def variant_type(ref, alt):
    if len(ref) == len(alt) == 1:
        return "SNV"
    if len(ref) == len(alt):
        return "MNV"
    if len(alt) > len(ref) and alt.startswith(ref):
        return "INSERTION"
    if len(ref) > len(alt) and ref.startswith(alt):
        return "DELETION"
    return "COMPLEX_ALLELE"


def parse_info(text):
    data = {}
    for item in text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            data[key] = value
        elif item:
            data[item] = True
    return data


def parse_csq(info, csq_fields):
    if "CSQ" not in info or not csq_fields:
        return []
    return [dict(zip(csq_fields, entry.split("|"))) for entry in str(info["CSQ"]).split(",")]


def scan_vcf(path, stage, samples, warnings):
    sample = sample_from_name(path)
    rows = []
    csq_fields = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("##INFO=<ID=CSQ"):
                match = re.search(r'Format: ([^">]+)', line)
                if match:
                    csq_fields = match.group(1).strip().split("|")
                continue
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                warnings.append(f"Malformed VCF row in {path}")
                continue
            chrom, pos, _, ref, alts, qual, filt, info_text = fields[:8]
            info = parse_info(info_text)
            csq = parse_csq(info, csq_fields)
            for alt in alts.split(","):
                if alt.startswith("<") or alt == "*":
                    continue
                vep_alleles = {alt, alt[len(ref):] if alt.startswith(ref) else alt}
                allele_csq = [x for x in csq if x.get("Allele") in vep_alleles]
                impacts = [clean(x.get("IMPACT")) for x in allele_csq]
                consequences = sorted({c for x in allele_csq for c in clean(x.get("Consequence")).split("&") if c})
                genes = sorted({clean(x.get("SYMBOL")) for x in allele_csq if clean(x.get("SYMBOL"))})
                proteins = sorted({clean(x.get("HGVSp")) for x in allele_csq if clean(x.get("HGVSp"))})
                impact = max(impacts, key=lambda x: IMPACT_RANK.get(x, 0), default="")
                rows.append({
                    "sample": sample, "stage": stage, "chrom": chrom.removeprefix("chr"), "position": as_int(pos),
                    "ref": ref, "alt": alt, "variant_type": variant_type(ref, alt), "filter": filt,
                    "quality": as_float(qual), "depth": as_int(info.get("DP")), "genes": ";".join(genes),
                    "impact": impact, "consequences": ";".join(consequences), "protein_changes": ";".join(proteins),
                    "key": f"{chrom.removeprefix('chr')}:{pos}:{ref}:{alt}", "source": str(path),
                })
    if sample not in samples:
        warnings.append(f"VCF sample not present in samplesheet: {sample} ({path})")
    return rows


def load_samples(results, project, warnings):
    candidates = [project / "samples.csv", results / "qc/samplesheet/samplesheet_design.tsv"]
    table = next((p for p in candidates if p.is_file()), None)
    if table is None:
        raise FileNotFoundError("samples.csv or qc/samplesheet/samplesheet_design.tsv was not found")
    rows = []
    for row in read_table(table):
        sample = clean(first(row, ["sample", "Sample"]))
        if not sample:
            continue
        rows.append({
            "sample": sample,
            "srr": clean(first(row, ["srr", "SRA", "accession"])),
            "subject": clean(first(row, ["TK", "tk", "subject", "Subject"], sample)),
            "group": clean(first(row, ["Group", "group"], sample)),
            "baseline": clean(first(row, ["baseline", "Baseline"])).lower() in {"true", "1", "yes"},
        })
    if not rows:
        raise ValueError(f"No samples found in {table}")
    by_subject = defaultdict(list)
    for row in rows:
        by_subject[row["subject"]].append(row)
    for subject, group in by_subject.items():
        count = sum(x["baseline"] for x in group)
        if count != 1:
            warnings.append(f"Subject {subject} has {count} baseline samples")
    return rows, table


def scan_explorer(results, samples, warnings):
    path = results / "igv/findings/finding_explorer/partitions/all.jsonl.gz"
    if not path.is_file():
        warnings.append(f"Finding Explorer records missing: {path}")
        return [], {}

    geometry_path = path.parent.parent / "event_geometry.json"
    geometry = {}
    if geometry_path.is_file():
        try:
            geometry = json.loads(geometry_path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError) as exc:
            warnings.append(f"Cannot read Explorer geometry: {exc}")

    records = []
    seen = set()
    per_sample = Counter()
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, 1):
            try:
                row = json.loads(line)
            except json.JSONDecodeError as exc:
                warnings.append(f"Invalid Explorer JSON line {line_number}: {exc}")
                continue
            event_id = clean(row.get("EventID"))
            if event_id in seen:
                warnings.append(f"Duplicate Explorer EventID: {event_id}")
            seen.add(event_id)
            sample = clean(row.get("Sample"))
            per_sample[sample] += 1
            alt = as_int(row.get("ExactAltReads"))
            ref_count = as_int(row.get("CleanReferenceReads"))
            excluded = as_int(row.get("ExcludedReads"))
            callable_count = as_int(row.get("CallableReads", row.get("CallableAlignments")))
            examined = as_int(row.get("TotalReadsExamined", row.get("UniqueAlignments")))
            if callable_count != alt + ref_count or examined != callable_count + excluded:
                warnings.append(f"Explorer count invariant failed: {event_id}")

            geometry_item = geometry.get(event_id, {})
            regions = []
            for region in geometry_item.get("regions", []):
                chrom = clean(region.get("chrom")).removeprefix("chr")
                start0 = as_int(region.get("start0"), -1)
                end0 = as_int(region.get("end0"), -1)
                role = clean(region.get("role"))
                if chrom and start0 >= 0 and end0 > start0:
                    regions.append((chrom, start0, end0, role))
            regions.sort()
            geometry_key = ";".join(f"{chrom}:{start0}:{end0}:{role}" for chrom, start0, end0, role in regions)

            records.append({
                "event_id": event_id, "sample": sample, "gene": clean(row.get("Gene")).upper(),
                "event_type": clean(row.get("EventType")), "impact": clean(row.get("PredictedImpact")),
                "consequence": clean(row.get("PredictedConsequence")), "chrom": clean(row.get("Chrom")).removeprefix("chr"),
                "position": as_int(row.get("Position")), "ref": clean(row.get("REF")), "alt": clean(row.get("ALT")),
                "evidence_classes": clean(row.get("EvidenceClasses")), "status": clean(row.get("VisualEvidenceStatus", row.get("ReadValidationStatus"))),
                "read_validation_status": clean(row.get("ReadValidationStatus")), "exact_alt": alt, "clean_reference": ref_count,
                "excluded": excluded, "callable": callable_count, "examined": examined,
                "alt_fraction": as_float(row.get("AltFractionAmongClean")), "context_alignments": as_int(row.get("ContextAlignments")),
                "protein_change": clean(row.get("ProteinChange")), "transcript": clean(row.get("Transcript")),
                "source_events": clean(row.get("SourceEvents")), "geometry_key": geometry_key,
            })

    geometry_count = len(geometry) if geometry else None
    if geometry_count is not None and geometry_count != len(records):
        warnings.append(f"Explorer geometry count differs: {geometry_count} geometry vs {len(records)} findings")
    return records, {"path": str(path), "findings": len(records), "geometry_count": geometry_count, "per_sample": dict(per_sample)}

def select_top_findings(records, limit):
    def score(row):
        structural = row["event_type"] in STRUCTURAL_TYPES
        evidence = row["context_alignments"] if structural else row["exact_alt"]
        fraction = row["alt_fraction"] if row["alt_fraction"] is not None else -1
        return IMPACT_RANK.get(row["impact"], 0), evidence, fraction, row["callable"]
    ordered = sorted(records, key=score, reverse=True)
    return ordered if limit <= 0 else ordered[:limit]


def scan_generic_events(results, subdir, kind, warnings):
    root = results / subdir
    rows = []
    if not root.is_dir():
        return rows
    for path in sorted(root.rglob("*.tsv")):
        if any(token in path.name.lower() for token in ("audit", "rejected", "failed", "summary")):
            continue
        try:
            for row in read_table(path):
                sample = clean(first(row, ["Sample", "sample"], sample_from_name(path)))
                gene = clean(first(row, ["Gene", "gene", "gene1", "Gene1", "gene_name", "Name"])).upper()
                partner = clean(first(row, ["gene2", "Gene2", "partner_gene", "PartnerGene"])).upper()
                status = clean(first(row, ["Status", "status", "ValidationStatus", "ReadValidationStatus", "Confidence"]))
                support = as_int(first(row, ["support", "Support", "split_reads1", "SplitReads", "junction_reads", "JunctionReads", "read_support"]))
                rows.append({"kind": kind, "sample": sample, "gene": gene, "partner": partner, "status": status, "support": support, "source": str(path)})
        except Exception as exc:
            warnings.append(f"Cannot parse {path}: {exc}")
    return rows


def scan_progression(results, warnings):
    root = results / "progression_vcf"
    rows = []
    variants = []
    if not root.is_dir():
        return rows, variants
    categories = {"shared_with_baseline": "shared", "non_baseline_only": "followup_only", "nonbaseline_only": "followup_only", "baseline_only": "baseline_only"}
    for path in sorted(root.glob("*.vcf.gz")):
        category = next((value for token, value in categories.items() if token in path.name.lower()), "other")
        sample = sample_from_name(path)
        genes = Counter()
        try:
            parsed = scan_vcf(path, f"progression_{category}", {sample}, warnings)
            for row in parsed:
                row["category"] = category
                variants.append(row)
                for gene in row["genes"].split(";"):
                    if gene:
                        genes[gene.upper()] += 1
        except Exception as exc:
            warnings.append(f"Cannot parse progression VCF {path}: {exc}")
            parsed = []
        rows.append({"sample": sample, "category": category, "variants": len(parsed), "top_genes": ";".join(g for g, _ in genes.most_common(10)), "source": str(path)})
    return rows, variants


def discover_expression_matrix(results, sample_names, warnings):
    """Load the published gene-expression matrix, preferring TPM values."""
    root = results / "expression"
    if not root.is_dir():
        return {}, "", []

    preferred = root / "gene_expression.gene_expression.tsv"
    candidates = [preferred] if preferred.is_file() else []
    candidates.extend(
        path for path in sorted(root.glob("*.tsv"))
        if path != preferred and not any(token in path.name.lower() for token in ("go", "summary", "metadata", "multiqc"))
    )

    ranked = []
    for path in candidates:
        try:
            with open_text(path) as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                fields = reader.fieldnames or []
            sample_columns = {}
            for sample in sorted(sample_names):
                choices = [f"{sample}_TPM", f"{sample}_CPM", f"{sample}_raw_count", sample]
                column = next((choice for choice in choices if choice in fields), None)
                if column:
                    sample_columns[sample] = column
            if sample_columns:
                metric_priority = sum(column.endswith("_TPM") for column in sample_columns.values())
                ranked.append((len(sample_columns), metric_priority, path == preferred, path.stat().st_size, path, fields, sample_columns))
        except Exception:
            continue

    if not ranked:
        warnings.append("No gene-by-sample expression matrix with sample TPM, CPM, raw-count, or exact sample columns was discovered")
        return {}, "", []

    _, _, _, _, path, fields, sample_columns = max(ranked)
    matrix = defaultdict(lambda: defaultdict(float))
    gene_columns = ["Gene", "gene", "GeneSymbol", "symbol", "gene_name", "Gene_Name", "Gene_ID", "gene_id", "Geneid"]
    try:
        for row in read_table(path):
            gene = clean(first(row, gene_columns)).split(".")[0].upper()
            if not gene:
                continue
            for sample, column in sample_columns.items():
                value = as_float(row.get(column))
                if value is not None and value >= 0:
                    matrix[gene][sample] += value
    except Exception as exc:
        warnings.append(f"Cannot parse expression matrix {path}: {exc}")
        return {}, str(path), fields

    if not matrix:
        warnings.append(f"Expression matrix had no usable gene rows: {path}")
    missing = sorted(set(sample_names) - set(sample_columns))
    if missing:
        warnings.append(f"Expression matrix is missing sample columns for: {', '.join(missing)}")
    return dict(matrix), str(path), fields

def scan_go(results, warnings):
    roots = [results / "expression/go", results / "progression_biology", results / "variant_landscape"]
    go_rows = []
    seen_paths = set()
    for root in roots:
        if not root.is_dir():
            continue
        for path in sorted(root.rglob("*.tsv")):
            if path in seen_paths or any(token in path.name.lower() for token in ("summary", "metadata")):
                continue
            seen_paths.add(path)
            try:
                for row in read_table(path):
                    go_id = clean(first(row, ["GO_ID", "go_id", "GO"]))
                    go_name = clean(first(row, ["GO_Name", "go_name", "Description", "description", "Term", "term", "name"]))
                    analysis = clean(first(row, ["Analysis", "analysis", "Comparison", "comparison"], path.stem))
                    fdr = as_float(first(row, ["FDR", "fdr", "padj", "adjusted_p_value", "qvalue", "q_value"]))
                    pvalue = as_float(first(row, ["PValue", "pvalue", "p_value", "pval"]))
                    mean_score = as_float(first(row, ["MeanScore", "mean_score"]))
                    z_score = as_float(first(row, ["ZScore", "z_score"]))
                    leading_genes = clean(first(row, ["LeadingGenes", "leading_genes", "OverlapGenes", "overlap_genes"]))
                    direction = ""
                    if mean_score is not None:
                        direction = "UP" if mean_score > 0 else "DOWN" if mean_score < 0 else "NEUTRAL"
                    if go_id or go_name:
                        go_rows.append({
                            "analysis": analysis, "sample_or_comparison": analysis.split(":", 1)[0],
                            "term": go_id, "description": go_name, "direction": direction,
                            "fdr": fdr, "pvalue": pvalue, "mean_score": mean_score,
                            "z_score": z_score, "leading_genes": leading_genes, "source": str(path),
                        })
            except Exception as exc:
                warnings.append(f"Cannot parse GO table {path}: {exc}")
    go_rows.sort(key=lambda x: (x["fdr"] is None, x["fdr"] if x["fdr"] is not None else 1, x["pvalue"] if x["pvalue"] is not None else 1))
    return go_rows

def aggregate(records, fields):
    counts = Counter(tuple(row.get(field, "") for field in fields) for row in records)
    return [dict(zip(fields, key), count=value) for key, value in sorted(counts.items())]


def log2_ratio(followup, baseline, pseudocount=0.5):
    return math.log2((followup + pseudocount) / (baseline + pseudocount))


def pathway_expression_scores(matrix, samples):
    rows = []
    if not matrix:
        return rows
    by_subject = defaultdict(list)
    for row in samples:
        by_subject[row["subject"]].append(row)
    for subject, group in by_subject.items():
        baseline_rows = [x for x in group if x["baseline"]]
        if len(baseline_rows) != 1:
            continue
        baseline = baseline_rows[0]["sample"]
        for followup_row in [x for x in group if not x["baseline"]]:
            followup = followup_row["sample"]
            for set_name, genes in GENE_SETS.items():
                gene_rows = []
                for gene in sorted(genes):
                    if gene not in matrix or baseline not in matrix[gene] or followup not in matrix[gene]:
                        continue
                    base_value = matrix[gene][baseline]
                    follow_value = matrix[gene][followup]
                    lfc = log2_ratio(follow_value, base_value)
                    gene_rows.append({"gene": gene, "baseline_value": base_value, "followup_value": follow_value, "log2_fold_change": lfc})
                changes = [x["log2_fold_change"] for x in gene_rows]
                detected = len(gene_rows)
                median_lfc = statistics.median(changes) if changes else None
                mean_lfc = statistics.mean(changes) if changes else None
                up_fraction = sum(value > 0 for value in changes) / detected if detected else None
                down_fraction = sum(value < 0 for value in changes) / detected if detected else None
                rows.append({
                    "subject": subject, "baseline": baseline, "followup": followup, "gene_set": set_name,
                    "configured_genes": len(genes), "detected_genes": detected,
                    "median_log2_fold_change": median_lfc, "mean_log2_fold_change": mean_lfc,
                    "up_fraction": up_fraction, "down_fraction": down_fraction,
                    "top_increased_genes": ";".join(f"{x['gene']}:{x['log2_fold_change']:.3f}" for x in sorted(gene_rows, key=lambda x: x["log2_fold_change"], reverse=True)[:15]),
                    "top_decreased_genes": ";".join(f"{x['gene']}:{x['log2_fold_change']:.3f}" for x in sorted(gene_rows, key=lambda x: x["log2_fold_change"])[:15]),
                })
    return rows


def gene_expression_evidence(matrix, samples):
    rows = []
    if not matrix:
        return rows
    gene_to_sets = defaultdict(list)
    for set_name, genes in GENE_SETS.items():
        for gene in genes:
            gene_to_sets[gene].append(set_name)
    by_subject = defaultdict(list)
    for row in samples:
        by_subject[row["subject"]].append(row)
    for subject, group in by_subject.items():
        baselines = [x for x in group if x["baseline"]]
        if len(baselines) != 1:
            continue
        baseline = baselines[0]["sample"]
        for followup_row in [x for x in group if not x["baseline"]]:
            followup = followup_row["sample"]
            for gene, sets in sorted(gene_to_sets.items()):
                if gene not in matrix or baseline not in matrix[gene] or followup not in matrix[gene]:
                    continue
                rows.append({
                    "subject": subject, "baseline": baseline, "followup": followup, "gene": gene,
                    "gene_sets": ";".join(sorted(sets)), "baseline_value": matrix[gene][baseline],
                    "followup_value": matrix[gene][followup],
                    "log2_fold_change": log2_ratio(matrix[gene][followup], matrix[gene][baseline]),
                })
    return rows


def splice_progression_evidence(explorer, samples, matrix):
    rows = []
    by_subject = defaultdict(list)
    for row in samples:
        by_subject[row["subject"]].append(row)

    supported_by_sample_gene = Counter()
    event_keys_by_sample = defaultdict(set)
    for row in explorer:
        if row["event_type"] != "SPLICE_JUNCTION":
            continue
        geometry_key = clean(row.get("geometry_key"))
        if not geometry_key:
            continue
        event_keys_by_sample[row["sample"]].add(geometry_key)
        if row["status"] == "CONTEXT_ALIGNMENTS_AVAILABLE":
            supported_by_sample_gene[(row["sample"], row["gene"])] += 1

    for subject, group in by_subject.items():
        baselines = [x for x in group if x["baseline"]]
        if len(baselines) != 1:
            continue
        baseline = baselines[0]["sample"]
        baseline_keys = event_keys_by_sample[baseline]
        for followup_row in [x for x in group if not x["baseline"]]:
            followup = followup_row["sample"]
            followup_keys = event_keys_by_sample[followup]
            rows.append({
                "subject": subject, "baseline": baseline, "followup": followup, "gene": "ALL",
                "baseline_events": len(baseline_keys), "followup_events": len(followup_keys),
                "shared_events": len(baseline_keys & followup_keys), "followup_only_events": len(followup_keys - baseline_keys),
                "baseline_only_events": len(baseline_keys - followup_keys), "baseline_expression": "", "followup_expression": "",
                "expression_normalized_log2_ratio": "",
            })
            genes = sorted({row["gene"] for row in explorer if row["event_type"] == "SPLICE_JUNCTION" and row["sample"] in {baseline, followup} and row["gene"]})
            for gene in genes:
                base_count = supported_by_sample_gene[(baseline, gene)]
                follow_count = supported_by_sample_gene[(followup, gene)]
                base_expr = matrix.get(gene, {}).get(baseline) if matrix else None
                follow_expr = matrix.get(gene, {}).get(followup) if matrix else None
                normalized = None
                if base_expr is not None and follow_expr is not None:
                    normalized = log2_ratio(follow_count / (base_expr + 0.5), base_count / (base_expr + 0.5)) if follow_expr is None else log2_ratio(follow_count / (follow_expr + 0.5), base_count / (base_expr + 0.5))
                rows.append({
                    "subject": subject, "baseline": baseline, "followup": followup, "gene": gene,
                    "baseline_events": base_count, "followup_events": follow_count,
                    "shared_events": "", "followup_only_events": "", "baseline_only_events": "",
                    "baseline_expression": base_expr if base_expr is not None else "",
                    "followup_expression": follow_expr if follow_expr is not None else "",
                    "expression_normalized_log2_ratio": normalized if normalized is not None else "",
                })
    return rows

def progression_variant_evidence(progression_variants, explorer, samples):
    rows = []
    explorer_by_sample_key = defaultdict(list)
    for row in explorer:
        if row["ref"] and row["alt"]:
            key = f"{row['chrom']}:{row['position']}:{row['ref']}:{row['alt']}"
            explorer_by_sample_key[(row["sample"], key)].append(row)
    gene_to_hypotheses = defaultdict(set)
    for hypothesis_id, spec in HYPOTHESES.items():
        for set_name in spec["positive_sets"] + spec["negative_sets"]:
            for gene in GENE_SETS[set_name]:
                gene_to_hypotheses[gene].add(hypothesis_id)
    for row in progression_variants:
        genes = [gene.upper() for gene in row["genes"].split(";") if gene]
        relevant = sorted({hyp for gene in genes for hyp in gene_to_hypotheses.get(gene, set())})
        if not relevant and not (set(genes) & ATX_CANDIDATE_GENES):
            continue
        explorer_matches = explorer_by_sample_key.get((row["sample"], row["key"]), [])
        direct = [x for x in explorer_matches if x["status"] in DIRECT_STATUSES and x["exact_alt"] > 0]
        rows.append({
            "sample": row["sample"], "category": row.get("category", ""), "key": row["key"],
            "genes": row["genes"], "impact": row["impact"], "consequences": row["consequences"],
            "protein_changes": row["protein_changes"], "hypotheses": ";".join(relevant),
            "atx_candidate_gene": ";".join(sorted(set(genes) & ATX_CANDIDATE_GENES)),
            "direct_explorer_support": bool(direct),
            "exact_alt_reads": max((x["exact_alt"] for x in direct), default=0),
            "max_alt_fraction": max((x["alt_fraction"] for x in direct if x["alt_fraction"] is not None), default=None),
            "source": row["source"],
        })
    return rows


def go_hypothesis_evidence(go_rows):
    output = []
    for row in go_rows:
        text = f"{row['term']} {row['description']}".lower()
        for hypothesis_id, spec in HYPOTHESES.items():
            matched = sorted({keyword for keyword in spec["go_keywords"] if keyword in text})
            if matched:
                output.append(dict(row, hypothesis=hypothesis_id, matched_keywords=";".join(matched), significant=row["fdr"] is not None and row["fdr"] <= 0.1))
    return output


def classify_pathway_row(row, expected_direction):
    detected = row["detected_genes"]
    median_lfc = row["median_log2_fold_change"]
    direction_fraction = row["up_fraction"] if expected_direction == "up" else row["down_fraction"]
    if detected < 5 or median_lfc is None:
        return "INSUFFICIENT_DATA"
    directional_median = median_lfc if expected_direction == "up" else -median_lfc
    if directional_median >= 0.25 and direction_fraction >= 0.60:
        return "SUPPORT"
    if directional_median <= -0.25 and direction_fraction <= 0.40:
        return "CONTRADICT"
    return "MIXED"


def evaluate_hypotheses(pathway_rows, go_evidence, splice_evidence, variant_evidence, samples):
    results = []
    pathway_lookup = {(row["followup"], row["gene_set"]): row for row in pathway_rows}
    followups = [row["sample"] for row in samples if not row["baseline"]]
    for hypothesis_id, spec in HYPOTHESES.items():
        per_followup = []
        for followup in followups:
            components = []
            for set_name in spec["positive_sets"]:
                row = pathway_lookup.get((followup, set_name))
                if row:
                    components.append((set_name, "up", classify_pathway_row(row, "up")))
            for set_name in spec["negative_sets"]:
                row = pathway_lookup.get((followup, set_name))
                if row:
                    components.append((set_name, "down", classify_pathway_row(row, "down")))
            informative = [status for _, _, status in components if status != "INSUFFICIENT_DATA"]
            support_count = informative.count("SUPPORT")
            contradict_count = informative.count("CONTRADICT")
            if not informative:
                verdict = "INSUFFICIENT_DATA"
            elif support_count >= max(1, math.ceil(len(informative) / 2)) and contradict_count == 0:
                verdict = "SUPPORTED"
            elif contradict_count >= max(1, math.ceil(len(informative) / 2)) and support_count == 0:
                verdict = "CONTRADICTED"
            else:
                verdict = "MIXED"
            per_followup.append({"followup": followup, "verdict": verdict, "components": components})
        informative_followups = [x["verdict"] for x in per_followup if x["verdict"] != "INSUFFICIENT_DATA"]
        if len(informative_followups) < len(followups):
            overall = "INSUFFICIENT_DATA" if not informative_followups else "PARTIAL_DATA"
        elif all(x == "SUPPORTED" for x in informative_followups):
            overall = "SUPPORTED"
        elif all(x == "CONTRADICTED" for x in informative_followups):
            overall = "CONTRADICTED"
        else:
            overall = "MIXED"
        matched_go = [x for x in go_evidence if x["hypothesis"] == hypothesis_id and x["significant"]]
        relevant_variants = [x for x in variant_evidence if hypothesis_id in x["hypotheses"].split(";")]
        details = []
        for entry in per_followup:
            component_text = ", ".join(f"{name}:{status}" for name, _, status in entry["components"]) or "no expression components"
            details.append(f"{entry['followup']}={entry['verdict']} ({component_text})")
        if hypothesis_id == "H2_rna_processing_quality_control":
            all_rows = [x for x in splice_evidence if x["gene"] == "ALL"]
            details.append("splice comparisons: " + "; ".join(f"{x['followup']} followup-only={x['followup_only_events']} shared={x['shared_events']}" for x in all_rows))
        results.append({
            "hypothesis": hypothesis_id, "title": spec["title"], "verdict": overall,
            "followup_results": "; ".join(details), "significant_go_rows": len(matched_go),
            "relevant_progression_variants": len(relevant_variants),
            "directly_supported_progression_variants": sum(bool(x["direct_explorer_support"]) for x in relevant_variants),
            "interpretation": "Expression-based verdict; variants, splice events and GO rows are supporting context and do not establish causality.",
        })
    # Refine H5 with the single-mutation check without changing H1-H4.
    h5 = next(row for row in results if row["hypothesis"] == "H5_composite_atx101_vulnerability")
    candidate_followup_only = [x for x in variant_evidence if x["category"] == "followup_only" and x["atx_candidate_gene"] and x["direct_explorer_support"]]
    shared_candidate_keys = Counter(x["key"] for x in candidate_followup_only)
    common_unique_candidates = [key for key, count in shared_candidate_keys.items() if count >= 2]
    h5["single_candidate_mutations_shared_by_followups"] = len(common_unique_candidates)
    h5["single_candidate_mutation_keys"] = ";".join(common_unique_candidates)
    if common_unique_candidates:
        h5["interpretation"] += " A shared follow-up-only directly supported candidate-gene allele exists, so a single-allele model cannot be rejected."
    else:
        h5["interpretation"] += " No shared follow-up-only directly supported allele in the predefined ATX-101 candidate genes was found; this weakens a simple single-mutation explanation."
    return results


def biological_markdown(results, samples, explorer, top, recurring, progression, fusions, splicing, go_rows, warnings, inventory):
    lines = [
        "# PGTK biological findings review", "", f"Results directory: `{results}`", "",
        "## Interpretation boundary", "",
        "This report summarizes RNA-derived research evidence. It does not establish DNA confirmation, somatic status, pathogenicity, treatment relevance, or a clinical diagnosis.", "",
        "## Dataset", "", f"- Samples: {len(samples)}", f"- Explorer findings: {len(explorer):,}",
        f"- Fusion table records discovered: {len(fusions):,}", f"- Splice table records discovered: {len(splicing):,}",
        f"- Progression VCF summaries: {len(progression):,}", f"- GO records discovered: {len(go_rows):,}", "",
        "## Sample design", "", "| Sample | Subject | Group | Baseline |", "|---|---|---|---|",
    ]
    for row in samples:
        lines.append(f"| {row['sample']} | {row['subject']} | {row['group']} | {str(row['baseline']).lower()} |")
    lines += ["", "## Finding counts", ""]
    for row in aggregate(explorer, ["sample", "event_type", "status"]):
        lines.append(f"- {row['sample']} / {row['event_type']} / {row['status']}: {row['count']:,}")
    lines += ["", "## Highest-priority RNA findings", "", "Ranking is for review triage only.", "", "| Sample | Gene | Type | Impact | Status | ALT/context | REF | ALT fraction | EventID |", "|---|---|---|---|---|---:|---:|---:|---|"]
    for row in top[:30]:
        evidence = row["context_alignments"] if row["event_type"] in STRUCTURAL_TYPES else row["exact_alt"]
        fraction = "NA" if row["alt_fraction"] is None else f"{row['alt_fraction']:.4f}"
        lines.append(f"| {row['sample']} | {row['gene']} | {row['event_type']} | {row['impact']} | {row['status']} | {evidence} | {row['clean_reference']} | {fraction} | {row['event_id']} |")
    lines += ["", "## Recurrent allele findings across samples", ""]
    for row in recurring[:30]:
        lines.append(f"- {row['key']}: {row['sample_count']} samples ({row['samples']}), genes={row['genes'] or 'NA'}")
    lines += ["", "## Progression summaries", ""]
    for row in progression:
        lines.append(f"- {row['sample']} / {row['category']}: {row['variants']:,} variants; top genes: {row['top_genes'] or 'NA'}")
    lines += ["", "## GO enrichment", ""]
    significant = [x for x in go_rows if x["fdr"] is not None and x["fdr"] <= 0.1]
    lines.append(f"- GO rows with FDR <= 0.1: {len(significant):,}")
    for row in significant[:20]:
        lines.append(f"- {row['analysis']}: {row['term']} {row['description']} (FDR={row['fdr']:.4g})")
    lines += ["", "## Data integrity and review warnings", ""]
    lines.extend(f"- {item}" for item in warnings[:200]) if warnings else lines.append("- No review-program warnings.")
    lines += ["", "## Output inventory", ""]
    lines.extend(f"- {key}: {value}" for key, value in inventory.items())
    return "\n".join(lines) + "\n"


def hypothesis_markdown(hypothesis_summary, pathway_rows, gene_rows, splice_rows, variant_rows, go_rows, expression_source):
    lines = [
        "# PGTK progression and ATX-101 hypothesis review", "",
        "## Interpretation boundary", "",
        "The verdicts are research-level assessments from RNA-derived PGTK outputs. They do not establish DNA-level clonality, protein activity, metabolic flux, drug causality, clinical actionability, or treatment response.", "",
        f"Expression matrix used: `{expression_source or 'not discovered'}`", "",
        "## Hypothesis verdicts", "", "| Hypothesis | Verdict | Significant matched GO rows | Relevant progression variants | Directly supported variants |", "|---|---|---:|---:|---:|",
    ]
    for row in hypothesis_summary:
        lines.append(f"| {row['title']} | {row['verdict']} | {row['significant_go_rows']} | {row['relevant_progression_variants']} | {row['directly_supported_progression_variants']} |")
    for row in hypothesis_summary:
        lines += ["", f"## {row['hypothesis']}: {row['title']}", "", f"Verdict: **{row['verdict']}**", "", row["followup_results"], "", row["interpretation"]]
        if row["hypothesis"] == "H5_composite_atx101_vulnerability":
            lines.append(f"Shared follow-up-only directly supported predefined candidate alleles: {row.get('single_candidate_mutations_shared_by_followups', 0)}")
    lines += ["", "## Pathway scoring rules", "", "A pathway supports an expected direction when at least five configured genes are detected, the median absolute directional log2 fold change is at least 0.25, and at least 60% of detected genes move in the expected direction. Opposite movement is classified as contradiction; other informative results are mixed.", ""]
    for row in pathway_rows:
        median_text = "NA" if row["median_log2_fold_change"] is None else f"{row['median_log2_fold_change']:.4f}"
        up_text = "NA" if row["up_fraction"] is None else f"{row['up_fraction']:.3f}"
        lines.append(f"- {row['followup']} vs {row['baseline']} / {row['gene_set']}: detected={row['detected_genes']}/{row['configured_genes']}, median log2FC={median_text}, up fraction={up_text}")
    lines += ["", "## Evidence inventories", "", f"- Hypothesis-set gene comparison rows: {len(gene_rows):,}", f"- Splice comparison rows: {len(splice_rows):,}", f"- Relevant progression-variant rows: {len(variant_rows):,}", f"- Matched GO rows: {len(go_rows):,}"]
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description="Generate a read-only biological and five-hypothesis review of completed PGTK results.")
    parser.add_argument("--results-dir", required=True, type=Path)
    parser.add_argument("--project-dir", type=Path, default=Path.cwd())
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--top-findings", type=int, default=200, help="Number of ranked findings to write; use 0 for all findings")
    args = parser.parse_args()
    results = args.results_dir.resolve()
    project = args.project_dir.resolve()
    output = args.output_dir.resolve()
    if not results.is_dir():
        parser.error(f"results directory does not exist: {results}")
    output.mkdir(parents=True, exist_ok=True)

    warnings = []
    samples, sample_source = load_samples(results, project, warnings)
    sample_names = {x["sample"] for x in samples}
    explorer, explorer_meta = scan_explorer(results, sample_names, warnings)
    top = select_top_findings(explorer, args.top_findings)

    vcf_stages = {"raw": results / "vcf_raw", "pass": results / "vcf_pass", "rna_validated": results / "rna_validation/variants"}
    vcf_records = []
    for stage, root in vcf_stages.items():
        if not root.is_dir():
            warnings.append(f"VCF stage directory missing: {root}")
            continue
        pattern = "*.raw.vcf.gz" if stage == "raw" else "*.pass.vcf.gz" if stage == "pass" else "*.rna.validated.vcf.gz"
        for path in sorted(root.glob(pattern)):
            vcf_records.extend(scan_vcf(path, stage, sample_names, warnings))

    recurrence_map = defaultdict(lambda: {"samples": set(), "genes": set()})
    for row in explorer:
        if row["ref"] and row["alt"]:
            key = f"{row['chrom']}:{row['position']}:{row['ref']}:{row['alt']}"
            recurrence_map[key]["samples"].add(row["sample"])
            if row["gene"]:
                recurrence_map[key]["genes"].add(row["gene"])
    recurring = [{"key": key, "sample_count": len(value["samples"]), "samples": ";".join(sorted(value["samples"])), "genes": ";".join(sorted(value["genes"]))} for key, value in recurrence_map.items() if len(value["samples"]) > 1]
    recurring.sort(key=lambda x: (-x["sample_count"], x["key"]))

    fusions = scan_generic_events(results, "rna_validation/fusions", "fusion", warnings)
    splicing = scan_generic_events(results, "rna_validation/splicing", "splice", warnings)
    progression, progression_variants = scan_progression(results, warnings)
    expression_summary_path = results / "expression/gene_expression.summary.tsv"
    expression_summary = list(read_table(expression_summary_path)) if expression_summary_path.is_file() else []
    go_rows = scan_go(results, warnings)
    expression_matrix, expression_source, expression_fields = discover_expression_matrix(results, sample_names, warnings)

    pathway_rows = pathway_expression_scores(expression_matrix, samples)
    gene_rows = gene_expression_evidence(expression_matrix, samples)
    splice_rows = splice_progression_evidence(explorer, samples, expression_matrix)
    variant_rows = progression_variant_evidence(progression_variants, explorer, samples)
    go_hypothesis_rows = go_hypothesis_evidence(go_rows)
    hypothesis_summary = evaluate_hypotheses(pathway_rows, go_hypothesis_rows, splice_rows, variant_rows, samples)

    finding_counts = aggregate(explorer, ["sample", "event_type", "impact", "status"])
    vcf_counts = aggregate(vcf_records, ["sample", "stage", "variant_type", "impact", "filter"])

    write_tsv(output / "samples.tsv", samples, ["sample", "srr", "subject", "group", "baseline"])
    write_tsv(output / "finding_counts.tsv", finding_counts, ["sample", "event_type", "impact", "status", "count"])
    top_fields = ["event_id", "sample", "gene", "event_type", "impact", "consequence", "status", "chrom", "position", "ref", "alt", "exact_alt", "clean_reference", "excluded", "callable", "alt_fraction", "context_alignments", "protein_change", "transcript", "evidence_classes"]
    write_tsv(output / "top_findings.tsv", top, top_fields)
    write_tsv(output / "recurrent_findings.tsv", recurring, ["key", "sample_count", "samples", "genes"])
    write_tsv(output / "vcf_stage_counts.tsv", vcf_counts, ["sample", "stage", "variant_type", "impact", "filter", "count"])
    write_tsv(output / "progression_summary.tsv", progression, ["sample", "category", "variants", "top_genes", "source"])
    write_tsv(output / "fusion_records.tsv", fusions, ["kind", "sample", "gene", "partner", "status", "support", "source"])
    write_tsv(output / "splice_records.tsv", splicing, ["kind", "sample", "gene", "partner", "status", "support", "source"])
    write_tsv(output / "go_top.tsv", go_rows[:1000], ["analysis", "sample_or_comparison", "term", "description", "direction", "fdr", "pvalue", "mean_score", "z_score", "leading_genes", "source"])

    pathway_fields = ["subject", "baseline", "followup", "gene_set", "configured_genes", "detected_genes", "median_log2_fold_change", "mean_log2_fold_change", "up_fraction", "down_fraction", "top_increased_genes", "top_decreased_genes"]
    write_tsv(output / "hypothesis_pathway_scores.tsv", pathway_rows, pathway_fields)
    write_tsv(output / "hypothesis_gene_evidence.tsv", gene_rows, ["subject", "baseline", "followup", "gene", "gene_sets", "baseline_value", "followup_value", "log2_fold_change"])
    write_tsv(output / "hypothesis_splice_evidence.tsv", splice_rows, ["subject", "baseline", "followup", "gene", "baseline_events", "followup_events", "shared_events", "followup_only_events", "baseline_only_events", "baseline_expression", "followup_expression", "expression_normalized_log2_ratio"])
    write_tsv(output / "hypothesis_variant_evidence.tsv", variant_rows, ["sample", "category", "key", "genes", "impact", "consequences", "protein_changes", "hypotheses", "atx_candidate_gene", "direct_explorer_support", "exact_alt_reads", "max_alt_fraction", "source"])
    write_tsv(output / "hypothesis_go_evidence.tsv", go_hypothesis_rows, ["hypothesis", "analysis", "sample_or_comparison", "term", "description", "direction", "fdr", "pvalue", "mean_score", "z_score", "leading_genes", "matched_keywords", "significant", "source"])
    summary_fields = ["hypothesis", "title", "verdict", "followup_results", "significant_go_rows", "relevant_progression_variants", "directly_supported_progression_variants", "single_candidate_mutations_shared_by_followups", "single_candidate_mutation_keys", "interpretation"]
    write_tsv(output / "hypothesis_summary.tsv", hypothesis_summary, summary_fields)
    write_tsv(output / "review_warnings.tsv", [{"warning": x} for x in warnings], ["warning"])

    inventory = {
        "sample_source": str(sample_source), "explorer_path": explorer_meta.get("path", ""),
        "explorer_findings": len(explorer), "geometry_count": explorer_meta.get("geometry_count"),
        "vcf_records_scanned": len(vcf_records), "fusion_rows_scanned": len(fusions), "splice_rows_scanned": len(splicing),
        "progression_files_scanned": len(progression), "progression_variants_scanned": len(progression_variants),
        "expression_summary_rows": len(expression_summary), "expression_matrix_source": expression_source,
        "expression_matrix_genes": len(expression_matrix), "expression_matrix_columns": expression_fields,
        "go_rows_scanned": len(go_rows), "hypothesis_pathway_rows": len(pathway_rows),
        "hypothesis_gene_rows": len(gene_rows), "hypothesis_splice_rows": len(splice_rows),
        "hypothesis_variant_rows": len(variant_rows), "hypothesis_go_rows": len(go_hypothesis_rows), "warnings": len(warnings),
    }

    biological_summary = {
        "results_dir": str(results), "project_dir": str(project), "samples": samples, "inventory": inventory,
        "explorer": explorer_meta, "finding_counts": finding_counts, "vcf_stage_counts": vcf_counts,
        "progression": progression, "top_findings": top, "recurrent_findings": recurring, "warnings": warnings,
    }
    hypothesis_json = {
        "hypotheses": hypothesis_summary, "gene_sets": {key: sorted(value) for key, value in GENE_SETS.items()},
        "pathway_scores": pathway_rows, "splice_evidence": splice_rows, "variant_evidence": variant_rows,
        "go_evidence": go_hypothesis_rows, "expression_source": expression_source,
    }
    (output / "biological_review.json").write_text(json.dumps(biological_summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (output / "hypothesis_review.json").write_text(json.dumps(hypothesis_json, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (output / "biological_review.md").write_text(biological_markdown(results, samples, explorer, top, recurring, progression, fusions, splicing, go_rows, warnings, inventory), encoding="utf-8")
    (output / "hypothesis_review.md").write_text(hypothesis_markdown(hypothesis_summary, pathway_rows, gene_rows, splice_rows, variant_rows, go_hypothesis_rows, expression_source), encoding="utf-8")

    manifest = []
    for path in sorted(output.iterdir()):
        if path.is_file() and path.name != "output_manifest.tsv":
            manifest.append({"file": path.name, "bytes": path.stat().st_size})
    write_tsv(output / "output_manifest.tsv", manifest, ["file", "bytes"])

    print(f"Biological review written to: {output}")
    print(f"Samples: {len(samples)}")
    print(f"Explorer findings: {len(explorer)}")
    print(f"VCF records scanned: {len(vcf_records)}")
    print(f"Expression matrix genes: {len(expression_matrix)}")
    print(f"Hypotheses evaluated: {len(hypothesis_summary)}")
    print(f"Warnings: {len(warnings)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())