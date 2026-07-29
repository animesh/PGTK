#wget -c https:/ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz
#python annotate_variant_peptides.py --candidates peptide_fasta_mapping.candidates.tsv --vep-vcf /cluster/home/ash022/scripts/pgtk/results/vep/TK12.vep.vcf.gz  /cluster/home/ash022/scripts/pgtk/results/vep/TK13.vep.vcf.gz /cluster/home/ash022/scripts/pgtk/results/vep/TK14.vep.vcf.gz --variant-fasta   /cluster/home/ash022/scripts/pgtk/results/combined_fasta/TK12.exploratory_proteogenomics.fasta        /cluster/home/ash022/scripts/pgtk/results/combined_fasta/TK13.exploratory_proteogenomics.fasta         /cluster/home/ash022/scripts/pgtk/results/combined_fasta/TK14.exploratory_proteogenomics.fasta --ensembl-pep /cluster/home/ash022/scripts/pgtk/Homo_sapiens.GRCh38.pep.all.fa.gz --il-equivalent --output-prefix variant_peptide_annotation
#tar cvzf variant_peptide_annotation_results.zip.txt variant_peptide_annotation.*
import argparse
import csv
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path


def open_text(path):
    path = Path(path)
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if path.suffix == ".gz" else path.open("rt", encoding="utf-8", errors="replace")


def norm_aa(seq, il=False):
    seq = re.sub(r"[^A-Za-z*]", "", seq).upper().rstrip("*")
    return seq.replace("I", "J").replace("L", "J") if il else seq


def parse_fasta(path):
    header = None
    parts = []
    with open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header = line[1:].strip()
                parts = []
            else:
                if header is None:
                    raise ValueError(f"Sequence before header in {path}")
                parts.append(line)
    if header is not None:
        yield header, "".join(parts)


def strip_version(identifier):
    return identifier.split(".", 1)[0] if identifier else ""


def parse_variant_header(header):
    token = header.split(None, 1)[0]
    match = re.match(r"^(TK\d+)\|var_([^|]+)$", token, re.I)
    if not match:
        return None
    sample = match.group(1).upper()
    body = match.group(2)
    # pypgatk format: ._CHROM.POS.REF.ALT_ENST.VERSION_FRAME
    match2 = re.match(r"^\.?_([^._]+)\.(\d+)\.([^._]+)\.([^_]+)_(ENST\d+(?:\.\d+)?)_(\d+)$", body)
    if not match2:
        return None
    chrom, pos, ref, alt, transcript, product = match2.groups()
    return {
        "sample": sample,
        "chrom": chrom.removeprefix("chr"),
        "pos": int(pos),
        "ref": ref,
        "alt": alt,
        "transcript": transcript,
        "transcript_base": strip_version(transcript),
        "product": product,
        "header": header,
        "token": token,
    }


def parse_csqs(vcf_paths):
    records = defaultdict(list)
    csq_fields_seen = None
    for vcf_path in vcf_paths:
        sample = Path(vcf_path).name.split(".", 1)[0]
        with open_text(vcf_path) as handle:
            csq_fields = None
            for line in handle:
                if line.startswith("##INFO=<ID=CSQ"):
                    match = re.search(r"Format: ([^\">]+)", line)
                    if not match:
                        raise ValueError(f"Cannot parse CSQ format in {vcf_path}")
                    csq_fields = match.group(1).strip().split("|")
                    csq_fields_seen = csq_fields
                    continue
                if line.startswith("#"):
                    continue
                if csq_fields is None:
                    raise ValueError(f"CSQ header not found before records in {vcf_path}")
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 8:
                    continue
                chrom, pos, _vid, ref, alts, _qual, _flt, info = fields[:8]
                info_map = {}
                for item in info.split(";"):
                    if "=" in item:
                        key, value = item.split("=", 1)
                        info_map[key] = value
                for csq in info_map.get("CSQ", "").split(","):
                    if not csq:
                        continue
                    values = csq.split("|")
                    annotation = dict(zip(csq_fields, values + [""] * (len(csq_fields) - len(values))))
                    transcript = annotation.get("Feature", "")
                    allele = annotation.get("Allele", "")
                    key = (sample.upper(), chrom.removeprefix("chr"), int(pos), ref, transcript)
                    annotation.update({"VCF_ALT": alts, "VCF_ALLELE": allele, "VCF_SAMPLE": sample.upper()})
                    records[key].append(annotation)
                    # Also index without sample for recovery if filenames are not TK-prefixed.
                    records[("", chrom.removeprefix("chr"), int(pos), ref, transcript)].append(annotation)
    if csq_fields_seen is None:
        raise ValueError("No VEP CSQ definition found in any VCF")
    return records


def choose_csq(records, variant):
    candidate_keys = []
    for sample in (variant["sample"], ""):
        for transcript in (variant["transcript"], variant["transcript_base"]):
            candidate_keys.append((sample, variant["chrom"], variant["pos"], variant["ref"], transcript))
    candidates = []
    seen = set()
    for key in candidate_keys:
        for csq in records.get(key, []):
            marker = tuple(sorted(csq.items()))
            if marker not in seen:
                candidates.append(csq)
                seen.add(marker)
    if not candidates:
        return None
    alt = variant["alt"]
    exact = [x for x in candidates if x.get("VCF_ALLELE") == alt or alt in x.get("VCF_ALT", "").split(",")]
    pool = exact or candidates
    # Prefer protein-changing records with explicit protein coordinates.
    pool.sort(key=lambda x: (not bool(x.get("Protein_position")), not bool(x.get("HGVSp")), not bool(x.get("Amino_acids"))))
    return pool[0]


def parse_protein_position(value):
    if not value:
        return None, None
    value = value.split("/", 1)[0]
    match = re.match(r"(\d+)(?:-(\d+))?", value)
    if not match:
        return None, None
    start = int(match.group(1))
    end = int(match.group(2) or start)
    return start, end


def load_reference_proteins(path):
    by_transcript = {}
    by_protein = {}
    for header, seq in parse_fasta(path):
        transcript = re.search(r"\btranscript:(ENST\d+(?:\.\d+)?)", header)
        protein = re.search(r"^(ENSP\d+(?:\.\d+)?)\b", header)
        if transcript:
            by_transcript[transcript.group(1)] = seq
            by_transcript[strip_version(transcript.group(1))] = seq
        if protein:
            by_protein[protein.group(1)] = seq
            by_protein[strip_version(protein.group(1))] = seq
    return by_transcript, by_protein


def load_variant_proteins(paths):
    proteins = defaultdict(list)
    for path in paths:
        for header, seq in parse_fasta(path):
            variant = parse_variant_header(header)
            if variant:
                proteins[variant["token"]].append((variant, seq, Path(path).name))
    return proteins


def load_candidates(path):
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Sequence", "All matching FASTA headers"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"Candidate TSV missing columns: {sorted(required - set(reader.fieldnames or []))}")
        return list(reader)


def main():
    parser = argparse.ArgumentParser(description="Annotate noncanonical MaxQuant peptides with VEP CSQ and test whether they span encoded protein changes.")
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--vep-vcf", nargs="+", required=True)
    parser.add_argument("--variant-fasta", nargs="+", required=True)
    parser.add_argument("--ensembl-pep", required=True)
    parser.add_argument("--output-prefix", default="variant_peptide_annotation")
    parser.add_argument("--il-equivalent", action="store_true", help="Also use I/L-equivalent matching if exact peptide matching fails")
    args = parser.parse_args()

    candidates = load_candidates(args.candidates)
    csq_records = parse_csqs(args.vep_vcf)
    ref_by_tx, ref_by_protein = load_reference_proteins(args.ensembl_pep)
    variant_proteins = load_variant_proteins(args.variant_fasta)

    output_rows = []
    unresolved = []
    for candidate in candidates:
        peptide = candidate["Sequence"].upper()
        headers = [x.strip() for x in candidate.get("All matching FASTA headers", "").split(";") if "|var_" in x]
        tokens = []
        for header in headers:
            token = header.split(None, 1)[0]
            if token not in tokens:
                tokens.append(token)
        if not tokens:
            continue
        for token in tokens:
            protein_entries = variant_proteins.get(token, [])
            if not protein_entries:
                parsed = parse_variant_header(token)
                unresolved.append((peptide, token, "variant header not found in supplied FASTAs"))
                continue
            for variant, alt_seq, fasta_name in protein_entries:
                csq = choose_csq(csq_records, variant)
                if csq is None:
                    unresolved.append((peptide, token, "matching VEP CSQ not found"))
                    continue
                peptide_start = alt_seq.find(peptide)
                match_mode = "exact"
                if peptide_start < 0 and args.il_equivalent:
                    peptide_start = norm_aa(alt_seq, True).find(norm_aa(peptide, True))
                    match_mode = "I/L-equivalent"
                peptide_start_1 = peptide_start + 1 if peptide_start >= 0 else None
                peptide_end_1 = peptide_start + len(peptide) if peptide_start >= 0 else None
                protein_start, protein_end = parse_protein_position(csq.get("Protein_position", ""))
                spans_change = "unknown"
                if peptide_start_1 is not None and protein_start is not None:
                    spans_change = "yes" if peptide_start_1 <= protein_end and peptide_end_1 >= protein_start else "no"

                protein_id = csq.get("ENSP", "") or csq.get("Protein", "")
                ref_seq = ref_by_tx.get(variant["transcript"]) or ref_by_tx.get(variant["transcript_base"])
                if not ref_seq and protein_id:
                    ref_seq = ref_by_protein.get(protein_id) or ref_by_protein.get(strip_version(protein_id))
                ref_peptide_exact = "unknown"
                if ref_seq:
                    ref_peptide_exact = "yes" if peptide in ref_seq else "no"
                    if ref_peptide_exact == "no" and args.il_equivalent and norm_aa(peptide, True) in norm_aa(ref_seq, True):
                        ref_peptide_exact = "I/L-equivalent"

                consequence = csq.get("Consequence", "")
                output_rows.append({
                    "Sequence": peptide,
                    "Observed pattern": candidate.get("Observed pattern", ""),
                    "PEP": candidate.get("PEP", ""),
                    "Score": candidate.get("Score", ""),
                    "MS/MS Count": candidate.get("MS/MS Count", ""),
                    "FASTA sample": variant["sample"],
                    "FASTA file": fasta_name,
                    "Chromosome": variant["chrom"],
                    "Position": variant["pos"],
                    "REF": variant["ref"],
                    "ALT": variant["alt"],
                    "Transcript": variant["transcript"],
                    "Protein ID": protein_id,
                    "Consequence": consequence,
                    "IMPACT": csq.get("IMPACT", ""),
                    "SYMBOL": csq.get("SYMBOL", ""),
                    "Gene": csq.get("Gene", ""),
                    "HGVSc": csq.get("HGVSc", ""),
                    "HGVSp": csq.get("HGVSp", ""),
                    "Protein position": csq.get("Protein_position", ""),
                    "Amino acids": csq.get("Amino_acids", ""),
                    "Codons": csq.get("Codons", ""),
                    "Peptide start in variant protein": peptide_start_1 or "",
                    "Peptide end in variant protein": peptide_end_1 or "",
                    "Peptide match mode": match_mode if peptide_start >= 0 else "not-found",
                    "Peptide spans VEP protein position": spans_change,
                    "Peptide found in Ensembl reference protein": ref_peptide_exact,
                    "Reference protein available": "yes" if ref_seq else "no",
                    "Variant FASTA header": variant["header"],
                })

    fields = [
        "Sequence", "Observed pattern", "PEP", "Score", "MS/MS Count", "FASTA sample", "FASTA file",
        "Chromosome", "Position", "REF", "ALT", "Transcript", "Protein ID", "Consequence", "IMPACT", "SYMBOL", "Gene",
        "HGVSc", "HGVSp", "Protein position", "Amino acids", "Codons",
        "Peptide start in variant protein", "Peptide end in variant protein", "Peptide match mode",
        "Peptide spans VEP protein position", "Peptide found in Ensembl reference protein", "Reference protein available", "Variant FASTA header"
    ]
    prefix = Path(args.output_prefix)
    detailed = prefix.with_suffix(".detailed.tsv")
    prioritized = prefix.with_suffix(".prioritized.tsv")
    summary = prefix.with_suffix(".summary.txt")
    unresolved_path = prefix.with_suffix(".unresolved.tsv")

    with detailed.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(output_rows)

    # This compatibility output contains all altered-residue associations that are
    # absent from the corresponding Ensembl reference protein. No analyst-defined
    # PEP, score, PSM-count, or peptide-length threshold is applied here.
    priority_rows = [
        row for row in output_rows
        if row.get("Peptide spans VEP protein position") == "yes"
        and row.get("Peptide found in Ensembl reference protein") in {"no", "unknown"}
    ]
    priority_rows.sort(key=lambda row: (row["Sequence"], row["FASTA sample"], row["Chromosome"], int(row["Position"]), row["Transcript"]))
    with prioritized.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(priority_rows)

    with unresolved_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Sequence", "Variant header", "Reason"]); writer.writerows(unresolved)

    unique_peptides = {r["Sequence"] for r in output_rows}
    spanning = {r["Sequence"] for r in output_rows if r["Peptide spans VEP protein position"] == "yes"}
    priority = {r["Sequence"] for r in priority_rows}
    with summary.open("w", encoding="utf-8") as handle:
        handle.write(f"Candidate rows supplied: {len(candidates)}\n")
        handle.write(f"Variant annotation rows: {len(output_rows)}\n")
        handle.write(f"Unique annotated peptides: {len(unique_peptides)}\n")
        handle.write(f"Unique peptides spanning VEP protein position: {len(spanning)}\n")
        handle.write(f"Unique reference-absent altered-residue peptides: {len(priority)}\n")
        handle.write(f"Unresolved mappings: {len(unresolved)}\n")
        handle.write("\nReference-absent altered-residue peptides (no analyst thresholds):\n")
        for peptide in sorted(priority):
            annotations = [r for r in priority_rows if r["Sequence"] == peptide]
            top = annotations[0]
            handle.write(f"  {peptide}\t{top['Observed pattern']}\t{top['SYMBOL']}\t{top['HGVSp']}\t{top['Amino acids']}\tPEP={top['PEP']}\tMSMS={top['MS/MS Count']}\n")

    print(f"Wrote {detailed}")
    print(f"Wrote {prioritized}")
    print(f"Wrote {summary}")
    print(f"Wrote {unresolved_path}")
    print(f"Unique reference-absent altered-residue peptides: {len(priority)}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
