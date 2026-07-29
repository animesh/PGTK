#!/usr/bin/env python3

import argparse
import csv
import gzip
import re
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

PROTEIN_ALTERING = {
    "missense_variant", "frameshift_variant", "stop_gained", "stop_lost",
    "start_lost", "splice_donor_variant", "splice_acceptor_variant",
    "inframe_insertion", "inframe_deletion",
}

def open_text(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def read_tsv(path):
    with open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_delimited(path):
    with open_text(path) as handle:
        first = handle.readline()
        delimiter = "\t" if "\t" in first else ","
        handle.seek(0)
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_tsv(path, rows, fields):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def split_ids(value):
    return [item.strip() for item in re.split(r"[;,]", value or "") if item.strip()]


def safe_float(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def safe_int(value, default=0):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def normalize_transcript_id(value):
    return (value or "").strip().split(".", 1)[0]


def normalize_chromosome(value):
    value = (value or "").strip()
    return value[3:] if value.lower().startswith("chr") else value


def sample_from_name(value):
    match = re.search(r"(TK12|TK13|TK14)", value or "", re.IGNORECASE)
    return match.group(1).upper() if match else ""


def sorted_join(values):
    return ";".join(sorted({str(value).strip() for value in values if str(value).strip()}))


def read_samples(path):
    samples = {}
    for row in read_delimited(path):
        sample = (row.get("sample") or "").strip()
        if sample:
            samples[sample] = {
                "srr": (row.get("srr") or "").strip(),
                "group": (row.get("Group") or "").strip(),
                "baseline": (row.get("baseline") or "").strip(),
            }
    return samples


def parse_mqpar(path):
    root = ET.parse(path).getroot()
    text = lambda name, default="": (root.findtext(name) or default).strip()
    fastas = [
        (node.text or "").strip()
        for node in root.findall("./fastaFiles/FastaFileInfo/fastaFilePath")
        if (node.text or "").strip()
    ]
    return {
        "version": text("maxQuantVersion", "unknown"),
        "include_contaminants": text("includeContaminants", "unknown"),
        "match_between_runs": text("matchBetweenRuns", "unknown"),
        "peptide_fdr": text("peptideFdr", "unknown"),
        "protein_fdr": text("proteinFdr", "unknown"),
        "min_peptide_length": text("minPeptideLength", "unknown"),
        "fasta_paths": fastas,
    }


def load_mapping(path):
    mapping = {}
    for row in read_tsv(path):
        sequence = (row.get("Sequence") or "").strip().upper()
        if sequence:
            mapping[sequence] = row
    return mapping


def load_evidence(path):
    by_sequence = defaultdict(list)
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"id", "Sequence", "Raw file", "Experiment", "PEP", "Score", "MS/MS IDs"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"evidence.txt missing columns: {sorted(missing)}")
        for row in reader:
            sequence = (row.get("Sequence") or "").strip().upper()
            if not sequence or (row.get("Decoy") or "").strip() == "+":
                continue
            by_sequence[sequence].append(row)
    return by_sequence


def load_msms(path):
    by_id = {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"id", "Raw file", "Scan number", "Sequence", "Score", "PEP"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"msms.txt missing columns: {sorted(missing)}")
        for row in reader:
            by_id[(row.get("id") or "").strip()] = row
    return by_id


def validate_protein_groups(path):
    with open_text(path) as handle:
        fields = set(csv.DictReader(handle, delimiter="\t").fieldnames or [])
    missing = {"id", "Protein IDs"} - fields
    if missing:
        raise ValueError(f"proteinGroups.txt missing columns: {sorted(missing)}")


def parse_vep_vcfs(paths):
    events = defaultdict(list)
    for path in paths:
        sample = sample_from_name(Path(path).name)
        with open_text(path) as handle:
            csq_fields = None
            for line in handle:
                if line.startswith("##INFO=<ID=CSQ"):
                    match = re.search(r"Format: ([^\">]+)", line)
                    if not match:
                        raise ValueError(f"Cannot parse CSQ header in {path}")
                    csq_fields = match.group(1).strip().split("|")
                    continue
                if line.startswith("#"):
                    continue
                if csq_fields is None:
                    raise ValueError(f"CSQ definition missing in {path}")
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 8:
                    continue
                chrom, pos, variant_id, ref, alts, _qual, filt, info = fields[:8]
                info_map = {}
                for item in info.split(";"):
                    if "=" in item:
                        key, value = item.split("=", 1)
                        info_map[key] = value
                for csq_text in info_map.get("CSQ", "").split(","):
                    if not csq_text:
                        continue
                    values = csq_text.split("|")
                    csq = dict(zip(csq_fields, values + [""] * (len(csq_fields) - len(values))))
                    consequences = set((csq.get("Consequence") or "").split("&"))
                    if not consequences.intersection(PROTEIN_ALTERING):
                        continue
                    alt = csq.get("Allele") or alts
                    event_key = (sample, normalize_chromosome(chrom), pos, ref, alt)
                    events[event_key].append({
                        "Variant ID": variant_id,
                        "FILTER": filt,
                        "Gene": csq.get("SYMBOL", ""),
                        "Gene ID": csq.get("Gene", ""),
                        "Transcript": csq.get("Feature", ""),
                        "Protein ID": csq.get("ENSP", "") or csq.get("Protein", ""),
                        "Consequence": csq.get("Consequence", ""),
                        "IMPACT": csq.get("IMPACT", ""),
                        "HGVSc": csq.get("HGVSc", ""),
                        "HGVSp": csq.get("HGVSp", ""),
                        "Protein position": csq.get("Protein_position", ""),
                        "Amino acids": csq.get("Amino_acids", ""),
                        "VEP VCF": Path(path).name,
                    })
    return events


def annotation_join_key(row):
    return (
        (row.get("FASTA sample") or "").upper(),
        normalize_chromosome(row.get("Chromosome")),
        str(row.get("Position") or ""),
        row.get("REF") or "",
        normalize_transcript_id(row.get("Transcript")),
    )


def peptide_has_clean_msms(sequence, evidence_by_sequence):
    rows = evidence_by_sequence.get(sequence, [])
    return any(
        split_ids(row.get("MS/MS IDs"))
        and (row.get("Potential contaminant") or "").strip() != "+"
        and (row.get("Decoy") or "").strip() != "+"
        for row in rows
    )


def passes_user_filters(row, mapping_row, filters):
    if not filters["enabled"]:
        return False
    sequence = (row.get("Sequence") or "").upper()
    if filters["max_pep"] is not None and safe_float(row.get("PEP"), 1.0) > filters["max_pep"]:
        return False
    if filters["min_score"] is not None and safe_float(row.get("Score"), 0.0) < filters["min_score"]:
        return False
    if filters["min_msms_count"] is not None and safe_int(row.get("MS/MS Count"), 0) < filters["min_msms_count"]:
        return False
    if filters["min_peptide_length"] is not None and len(sequence) < filters["min_peptide_length"]:
        return False
    if filters["require_canonical_absence"] and mapping_row.get("Absent from canonical FASTA") != "yes":
        return False
    if filters["require_reference_absence"] and row.get("Peptide found in Ensembl reference protein") not in {"no", "unknown"}:
        return False
    if filters["exclude_contaminant_matches"] and safe_int(mapping_row.get("Contaminant matches"), 0) != 0:
        return False
    if filters["exclude_decoy_matches"] and safe_int(mapping_row.get("Decoy matches"), 0) != 0:
        return False
    return True


def evidence_summary(peptides, evidence_by_sequence, msms_by_id):
    raw_files, experiments, evidence_ids, msms_ids, scans = set(), set(), set(), set(), set()
    peps, scores = [], []
    contaminant = False
    for peptide in sorted(peptides):
        for row in evidence_by_sequence.get(peptide, []):
            raw = (row.get("Raw file") or "").strip()
            experiment = (row.get("Experiment") or "").strip()
            if raw:
                raw_files.add(raw)
            if experiment:
                experiments.add(experiment)
            if row.get("id"):
                evidence_ids.add(row["id"])
            contaminant = contaminant or (row.get("Potential contaminant") or "").strip() == "+"
            peps.append(safe_float(row.get("PEP"), 1.0))
            scores.append(safe_float(row.get("Score"), 0.0))
            for msms_id in split_ids(row.get("MS/MS IDs")):
                msms_ids.add(msms_id)
                msms = msms_by_id.get(msms_id)
                if msms and msms.get("Scan number"):
                    scans.add(f"{msms.get('Raw file') or raw}:{msms['Scan number']}")
    detection_samples = {sample_from_name(raw) for raw in raw_files}
    detection_samples.discard("")
    number_key = lambda value: (0, int(value)) if value.isdigit() else (1, value)
    return {
        "MS evidence": "yes" if raw_files else "no",
        "MS detection samples": sorted_join(detection_samples),
        "Raw files": sorted_join(raw_files),
        "Experiments": sorted_join(experiments),
        "Evidence IDs": ";".join(sorted(evidence_ids, key=number_key)),
        "MS/MS IDs": ";".join(sorted(msms_ids, key=number_key)),
        "Scan numbers": sorted_join(scans),
        "Best PEP": f"{min(peps):.6g}" if peps else "",
        "Best score": f"{max(scores):.6g}" if scores else "",
        "PSM count": str(len(msms_ids)),
        "Contaminant evidence": "yes" if contaminant else "no",
    }


def build_variant_rows(events, annotations, mapping, evidence_by_sequence, msms_by_id, samples, mqpar, filters):
    annotations_by_key = defaultdict(list)
    for row in annotations:
        annotations_by_key[annotation_join_key(row)].append(row)

    mq_min_length = safe_int(mqpar.get("min_peptide_length"), 0)
    output = []
    for event_key, consequences in events.items():
        sample, chrom, pos, ref, alt = event_key
        event_annotations = []
        for consequence in consequences:
            key = (sample, chrom, pos, ref, normalize_transcript_id(consequence.get("Transcript")))
            event_annotations.extend(annotations_by_key.get(key, []))

        unique_annotations = {}
        for annotation in event_annotations:
            identity = (
                annotation.get("Sequence", ""), annotation.get("Transcript", ""),
                annotation.get("HGVSp", ""), annotation.get("Variant FASTA header", ""),
            )
            unique_annotations[identity] = annotation
        event_annotations = list(unique_annotations.values())

        all_peptides = {(row.get("Sequence") or "").upper() for row in event_annotations if row.get("Sequence")}
        altered_rows = [
            row for row in event_annotations
            if row.get("Peptide spans VEP protein position") == "yes" and row.get("Sequence")
        ]
        altered_peptides = {(row.get("Sequence") or "").upper() for row in altered_rows}

        search_consistent_peptides = {
            (row.get("Sequence") or "").upper()
            for row in altered_rows
            if len((row.get("Sequence") or "")) >= mq_min_length
            and safe_int(row.get("MS/MS Count"), 0) >= 1
            and peptide_has_clean_msms((row.get("Sequence") or "").upper(), evidence_by_sequence)
        }
        canonical_absent_peptides = {
            sequence for sequence in altered_peptides
            if mapping.get(sequence, {}).get("Absent from canonical FASTA") == "yes"
            and safe_int(mapping.get(sequence, {}).get("Canonical matches"), 0) == 0
        }
        reference_absent_peptides = {
            (row.get("Sequence") or "").upper()
            for row in altered_rows
            if row.get("Peptide found in Ensembl reference protein") in {"no", "unknown"}
        }
        canonical_and_reference_absent = canonical_absent_peptides & reference_absent_peptides
        user_filtered_peptides = {
            (row.get("Sequence") or "").upper()
            for row in altered_rows
            if passes_user_filters(row, mapping.get((row.get("Sequence") or "").upper(), {}), filters)
        }

        evidence = evidence_summary(altered_peptides, evidence_by_sequence, msms_by_id)
        meta = samples.get(sample, {})
        output.append({
            "Sample": sample,
            "SRA": meta.get("srr", ""),
            "Variant": f"{chrom}:{pos}:{ref}>{alt}",
            "Chromosome": chrom,
            "Position": pos,
            "REF": ref,
            "ALT": alt,
            "Genes": sorted_join(row.get("Gene") for row in consequences),
            "Gene IDs": sorted_join(row.get("Gene ID") for row in consequences),
            "Transcripts": sorted_join(row.get("Transcript") for row in consequences),
            "Protein IDs": sorted_join(row.get("Protein ID") for row in consequences),
            "Consequences": sorted_join(row.get("Consequence") for row in consequences),
            "IMPACT": sorted_join(row.get("IMPACT") for row in consequences),
            "HGVSc": sorted_join(row.get("HGVSc") for row in consequences),
            "HGVSp": sorted_join(row.get("HGVSp") for row in consequences),
            "All mapped variant-protein peptides": sorted_join(all_peptides),
            "Altered-residue peptides": sorted_join(altered_peptides),
            "Search-consistent altered-residue peptides": sorted_join(search_consistent_peptides),
            "Canonical-absent altered-residue peptides": sorted_join(canonical_absent_peptides),
            "Ensembl-reference-absent altered-residue peptides": sorted_join(reference_absent_peptides),
            "Canonical-and-reference-absent peptides": sorted_join(canonical_and_reference_absent),
            "User-filtered altered-residue peptides": sorted_join(user_filtered_peptides),
            **evidence,
            "Evidence category": (
                "altered-residue association" if altered_peptides
                else "variant-protein peptide not spanning altered residue" if all_peptides
                else "RNA-only variant"
            ),
            "VEP VCF": sorted_join(row.get("VEP VCF") for row in consequences),
        })

    chrom_key = lambda value: (0, int(value)) if value.isdigit() else (1, value)
    output.sort(key=lambda row: (row["Sample"], chrom_key(row["Chromosome"]), int(row["Position"]), row["REF"], row["ALT"]))
    return output


def build_junction_rows(splice_rows, evidence_by_sequence, msms_by_id, samples):
    grouped = defaultdict(list)
    for row in splice_rows:
        if row.get("Peptide classification") != "exact-junction-spanning":
            continue
        key = (
            (row.get("Sequence") or "").upper(),
            row.get("Genomic junction", ""),
            row.get("Reference junction status", ""),
            row.get("Junction support", ""),
            row.get("Anchor residues left", ""),
            row.get("Anchor residues right", ""),
        )
        grouped[key].append(row)

    output = []
    for key, rows in grouped.items():
        sequence, junction, reference_status, support, left, right = key
        source_samples = {row.get("Sample", "") for row in rows if row.get("Sample")}
        evidence = evidence_summary({sequence}, evidence_by_sequence, msms_by_id)
        output.append({
            "Sequence": sequence,
            "RNA source samples": sorted_join(source_samples),
            "RNA source SRAs": sorted_join(samples.get(sample, {}).get("srr", "") for sample in source_samples),
            "Genomic junction": junction,
            "Reference status": reference_status,
            "Support": support,
            "Anchors": f"{left}+{right}",
            "MS detection samples": evidence["MS detection samples"],
            "Raw files": evidence["Raw files"],
            "Best PEP": evidence["Best PEP"],
            "Best score": evidence["Best score"],
            "PSM count": evidence["PSM count"],
        })
    output.sort(key=lambda row: (row["Reference status"], row["Sequence"], row["Genomic junction"]))
    return output


def md_escape(value):
    return str(value or "").replace("|", "\\|").replace("\n", " ")


def unique_peptides(rows, field):
    return {peptide for row in rows for peptide in split_ids(row.get(field, ""))}


def write_event_table(handle, rows, peptide_field, empty_message):
    columns = ["RNA sample", "SRA", "Variant", "Genes", "HGVSp", "Altered peptide", "MS samples", "Raw files", "Best PEP", "Best score", "PSMs"]
    handle.write("| " + " | ".join(columns) + " |\n")
    handle.write("|" + "---|" * len(columns) + "\n")
    selected = [row for row in rows if row.get(peptide_field)]
    for row in selected:
        values = [
            row["Sample"], row["SRA"], row["Variant"], row["Genes"], row["HGVSp"],
            row[peptide_field], row["MS detection samples"], row["Raw files"],
            row["Best PEP"], row["Best score"], row["PSM count"],
        ]
        handle.write("| " + " | ".join(md_escape(value) for value in values) + " |\n")
    if not selected:
        handle.write(f"| {md_escape(empty_message)} |  |  |  |  |  |  |  |  |  |  |\n")


def write_markdown(path, variant_rows, junction_rows, mqpar, samples, fasta_status, filters):
    altered_events = [row for row in variant_rows if row["Altered-residue peptides"]]
    search_events = [row for row in variant_rows if row["Search-consistent altered-residue peptides"]]
    canonical_events = [row for row in variant_rows if row["Canonical-absent altered-residue peptides"]]
    reference_events = [row for row in variant_rows if row["Ensembl-reference-absent altered-residue peptides"]]
    both_events = [row for row in variant_rows if row["Canonical-and-reference-absent peptides"]]
    user_events = [row for row in variant_rows if row["User-filtered altered-residue peptides"]]
    mapped_events = [row for row in variant_rows if row["All mapped variant-protein peptides"]]
    novel_junctions = [row for row in junction_rows if row["Reference status"] == "reference-absent"]

    with Path(path).open("w", encoding="utf-8") as handle:
        handle.write("# PGTK Proteogenomics Evidence Report\n\n")
        handle.write("## Evidence overview\n\n")
        handle.write(f"- Protein-altering sample-specific genomic variants: {len(variant_rows)}\n")
        handle.write(f"- Variants with any mapped variant-protein peptide: {len(mapped_events)}\n")
        handle.write(f"- Variants with an altered-residue peptide association: {len(altered_events)}\n")
        handle.write(f"- Unique altered-residue peptides: {len(unique_peptides(altered_events, 'Altered-residue peptides'))}\n")
        handle.write(f"- Search-consistent variant events: {len(search_events)}\n")
        handle.write(f"- Canonical-absent variant events: {len(canonical_events)}\n")
        handle.write(f"- Ensembl-reference-absent variant events: {len(reference_events)}\n")
        handle.write(f"- Canonical- and reference-absent variant events: {len(both_events)}\n")
        handle.write(f"- Unique translated junction-spanning peptides: {len({row['Sequence'] for row in junction_rows})}\n")
        handle.write(f"- Reference-absent translated junctions: {len(novel_junctions)}\n\n")
        handle.write("`PEP`, score, and PSM count are reported as evidence attributes. They are not used as default analyst-defined thresholds.\n\n")

        handle.write("## Block A: all altered-residue associations\n\n")
        handle.write("No additional analyst threshold is applied. The peptide must map to the variant protein and span the altered amino-acid position.\n\n")
        write_event_table(handle, altered_events, "Altered-residue peptides", "No altered-residue associations")

        handle.write("\n## Block B: search-consistent altered-residue evidence\n\n")
        handle.write(f"Uses only search-derived or MaxQuant-output conditions: peptide length >= {mqpar['min_peptide_length']} from `mqpar.xml`, at least one MS/MS identification, not reverse/decoy, and not marked as a potential contaminant by MaxQuant. No additional PEP, score, or PSM-count threshold is imposed.\n\n")
        write_event_table(handle, search_events, "Search-consistent altered-residue peptides", "No search-consistent altered-residue evidence")

        handle.write("\n## Block C: sequence-novelty subsets\n\n")
        handle.write("These are classification subsets, not confidence thresholds. Canonical absence refers to the searched canonical FASTA; reference absence refers to the corresponding Ensembl reference protein.\n\n")
        handle.write("### Canonical-absent associations\n\n")
        write_event_table(handle, canonical_events, "Canonical-absent altered-residue peptides", "No canonical-absent associations")
        handle.write("\n### Ensembl-reference-absent associations\n\n")
        write_event_table(handle, reference_events, "Ensembl-reference-absent altered-residue peptides", "No reference-absent associations")
        handle.write("\n### Canonical- and reference-absent associations\n\n")
        write_event_table(handle, both_events, "Canonical-and-reference-absent peptides", "No canonical- and reference-absent associations")

        handle.write("\n## Block D: optional user-filtered associations\n\n")
        if filters["enabled"]:
            handle.write("User-selected criteria were enabled. Values and sources:\n\n")
            for key, value in filters.items():
                if key != "enabled":
                    handle.write(f"- `{key}`: `{value}` (user parameter)\n")
            handle.write("\n")
            write_event_table(handle, user_events, "User-filtered altered-residue peptides", "No associations passed the selected user filters")
        else:
            handle.write("Disabled by default. No user-selected PEP, score, PSM-count, peptide-length, canonical-absence, reference-absence, contaminant-match, or decoy-match filter was applied.\n")

        handle.write("\n## Block E: translated splice-junction evidence\n\n")
        columns = ["Peptide", "RNA source samples", "RNA source SRAs", "Junction", "Reference status", "Support", "Anchors", "MS samples", "Raw files"]
        handle.write("| " + " | ".join(columns) + " |\n")
        handle.write("|" + "---|" * len(columns) + "\n")
        for row in junction_rows:
            values = [row["Sequence"], row["RNA source samples"], row["RNA source SRAs"], row["Genomic junction"], row["Reference status"], row["Support"], row["Anchors"], row["MS detection samples"], row["Raw files"]]
            handle.write("| " + " | ".join(md_escape(value) for value in values) + " |\n")
        if not junction_rows:
            handle.write("| No translated junction-spanning peptide evidence |  |  |  |  |  |  |  |  |\n")

        handle.write("\n## Search provenance and threshold sources\n\n")
        handle.write(f"- MaxQuant version: `{md_escape(mqpar['version'])}`\n")
        handle.write(f"- Minimum peptide length: `{md_escape(mqpar['min_peptide_length'])}` (mqpar.xml)\n")
        handle.write(f"- Peptide FDR: `{md_escape(mqpar['peptide_fdr'])}` (mqpar.xml; not converted into a PEP cutoff)\n")
        handle.write(f"- Protein FDR: `{md_escape(mqpar['protein_fdr'])}` (mqpar.xml)\n")
        handle.write(f"- Contaminants enabled: `{md_escape(mqpar['include_contaminants'])}` (mqpar.xml)\n")
        handle.write(f"- Match between runs: `{md_escape(mqpar['match_between_runs'])}` (mqpar.xml)\n\n")
        for fasta in mqpar["fasta_paths"]:
            handle.write(f"- `{md_escape(fasta)}`: {md_escape(fasta_status.get(fasta, 'historical path unavailable'))}\n")

        handle.write("\n## Sample and SRA mapping\n\n")
        handle.write("| Sample | SRA accession | Group | Baseline |\n|---|---|---|---|\n")
        for sample in sorted(samples):
            meta = samples[sample]
            handle.write(f"| {md_escape(sample)} | {md_escape(meta['srr'])} | {md_escape(meta['group'])} | {md_escape(meta['baseline'])} |\n")


def main():
    parser = argparse.ArgumentParser(description="Generate deduplicated PGTK variant and splice evidence reports.")
    parser.add_argument("--samples", required=True)
    parser.add_argument("--mqpar", required=True)
    parser.add_argument("--evidence", required=True)
    parser.add_argument("--msms", required=True)
    parser.add_argument("--protein-groups", required=True)
    parser.add_argument("--vep-vcf", nargs="+", required=True)
    parser.add_argument("--variant-annotation", required=True)
    parser.add_argument("--peptide-mapping", required=True)
    parser.add_argument("--splice-validation", required=True)
    parser.add_argument("--searched-fasta", nargs="*", default=[])
    parser.add_argument("--output-prefix", default="proteogenomics_evidence")
    parser.add_argument("--user-max-pep", type=float)
    parser.add_argument("--user-min-score", type=float)
    parser.add_argument("--user-min-msms-count", type=int)
    parser.add_argument("--user-min-peptide-length", type=int)
    parser.add_argument("--user-require-canonical-absence", action="store_true")
    parser.add_argument("--user-require-reference-absence", action="store_true")
    parser.add_argument("--user-exclude-contaminant-matches", action="store_true")
    parser.add_argument("--user-exclude-decoy-matches", action="store_true")
    args = parser.parse_args()

    samples = read_samples(args.samples)
    mqpar = parse_mqpar(args.mqpar)
    mapping = load_mapping(args.peptide_mapping)
    evidence_by_sequence = load_evidence(args.evidence)
    msms_by_id = load_msms(args.msms)
    validate_protein_groups(args.protein_groups)
    events = parse_vep_vcfs(args.vep_vcf)
    annotations = read_tsv(args.variant_annotation)
    splice_rows = read_tsv(args.splice_validation)

    filters = {
        "max_pep": args.user_max_pep,
        "min_score": args.user_min_score,
        "min_msms_count": args.user_min_msms_count,
        "min_peptide_length": args.user_min_peptide_length,
        "require_canonical_absence": args.user_require_canonical_absence,
        "require_reference_absence": args.user_require_reference_absence,
        "exclude_contaminant_matches": args.user_exclude_contaminant_matches,
        "exclude_decoy_matches": args.user_exclude_decoy_matches,
    }
    filters["enabled"] = any(value is not None and value is not False for value in filters.values())

    variant_rows = build_variant_rows(events, annotations, mapping, evidence_by_sequence, msms_by_id, samples, mqpar, filters)
    junction_rows = build_junction_rows(splice_rows, evidence_by_sequence, msms_by_id, samples)

    variant_fields = [
        "Sample", "SRA", "Variant", "Chromosome", "Position", "REF", "ALT", "Genes", "Gene IDs",
        "Transcripts", "Protein IDs", "Consequences", "IMPACT", "HGVSc", "HGVSp",
        "All mapped variant-protein peptides", "Altered-residue peptides", "Search-consistent altered-residue peptides",
        "Canonical-absent altered-residue peptides", "Ensembl-reference-absent altered-residue peptides",
        "Canonical-and-reference-absent peptides", "User-filtered altered-residue peptides",
        "MS evidence", "MS detection samples", "Raw files", "Experiments", "Evidence IDs", "MS/MS IDs",
        "Scan numbers", "Best PEP", "Best score", "PSM count", "Contaminant evidence", "Evidence category", "VEP VCF",
    ]
    junction_fields = [
        "Sequence", "RNA source samples", "RNA source SRAs", "Genomic junction", "Reference status",
        "Support", "Anchors", "MS detection samples", "Raw files", "Best PEP", "Best score", "PSM count",
    ]

    prefix = Path(args.output_prefix)
    variants_path = Path(f"{prefix}.variants.tsv")
    junctions_path = Path(f"{prefix}.junctions.tsv")
    report_path = Path(f"{prefix}.report.md")
    summary_path = Path(f"{prefix}.summary.txt")

    write_tsv(variants_path, variant_rows, variant_fields)
    write_tsv(junctions_path, junction_rows, junction_fields)

    supplied = {str(Path(path).resolve()): path for path in args.searched_fasta}
    fasta_status = {}
    for fasta in mqpar["fasta_paths"]:
        resolved = str(Path(fasta).resolve())
        if resolved in supplied:
            fasta_status[fasta] = "supplied and exact path matched"
        elif Path(fasta).exists():
            fasta_status[fasta] = "historical search path exists"
        elif any(Path(path).name == Path(fasta).name for path in args.searched_fasta):
            fasta_status[fasta] = "supplied by matching filename"
        else:
            fasta_status[fasta] = "historical path unavailable"

    write_markdown(report_path, variant_rows, junction_rows, mqpar, samples, fasta_status, filters)

    altered_events = [row for row in variant_rows if row["Altered-residue peptides"]]
    search_events = [row for row in variant_rows if row["Search-consistent altered-residue peptides"]]
    canonical_events = [row for row in variant_rows if row["Canonical-absent altered-residue peptides"]]
    reference_events = [row for row in variant_rows if row["Ensembl-reference-absent altered-residue peptides"]]
    both_events = [row for row in variant_rows if row["Canonical-and-reference-absent peptides"]]
    user_events = [row for row in variant_rows if row["User-filtered altered-residue peptides"]]
    mapped_events = [row for row in variant_rows if row["All mapped variant-protein peptides"]]
    novel_junctions = [row for row in junction_rows if row["Reference status"] == "reference-absent"]

    with summary_path.open("w", encoding="utf-8") as handle:
        handle.write(f"Protein-altering sample-specific genomic variants: {len(variant_rows)}\n")
        handle.write(f"Variants with any mapped variant-protein peptide: {len(mapped_events)}\n")
        handle.write(f"Variants with altered-residue peptide association: {len(altered_events)}\n")
        handle.write(f"Unique altered-residue peptides: {len(unique_peptides(altered_events, 'Altered-residue peptides'))}\n")
        handle.write(f"Search-consistent variant events: {len(search_events)}\n")
        handle.write(f"Canonical-absent variant events: {len(canonical_events)}\n")
        handle.write(f"Ensembl-reference-absent variant events: {len(reference_events)}\n")
        handle.write(f"Canonical-and-reference-absent variant events: {len(both_events)}\n")
        handle.write(f"User-filtered variant events: {len(user_events) if filters['enabled'] else 'disabled'}\n")
        handle.write(f"Unique translated junction-spanning peptides: {len({row['Sequence'] for row in junction_rows})}\n")
        handle.write(f"Reference-absent translated junctions: {len(novel_junctions)}\n")
        handle.write(f"MaxQuant minimum peptide length: {mqpar['min_peptide_length']}\n")
        handle.write(f"MaxQuant version: {mqpar['version']}\n")
        handle.write(f"Match between runs: {mqpar['match_between_runs']}\n")

    for output in (variants_path, junctions_path, report_path, summary_path):
        print(f"Wrote {output}")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
