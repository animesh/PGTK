import argparse
import csv
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path


def open_text(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def normalize_sequence(sequence, il_equivalent=True):
    sequence = re.sub(r"[^A-Za-z]", "", sequence).upper()
    if il_equivalent:
        sequence = sequence.replace("I", "J").replace("L", "J")
    return sequence


def sample_from_path(path):
    sample = Path(path).name.split(".", 1)[0].strip()
    if not sample:
        raise ValueError(f"Cannot determine sample from GTF filename: {path}")
    return sample


def parse_gtf_attributes(text):
    attributes = {}
    for item in text.strip().strip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item and " " not in item.split("=", 1)[0]:
            key, value = item.split("=", 1)
        else:
            parts = item.split(None, 1)
            if len(parts) != 2:
                continue
            key, value = parts
        attributes[key] = value.strip().strip('"')
    return attributes


def add_gtf_record(records, key, sample, path, fields, attributes):
    chrom, _, _, start, end, _, strand, _, _ = fields
    record = records.get(key)
    if record is None:
        record = {
            "sample": sample,
            "chrom": chrom,
            "strand": strand,
            "exons": [],
            "attrs": dict(attributes),
            "sources": set(),
        }
        records[key] = record
    elif record["chrom"] != chrom or record["strand"] != strand:
        raise ValueError(
            f"Conflicting transcript model for {key}: "
            f"{record['chrom']}:{record['strand']} versus {chrom}:{strand}"
        )

    record["exons"].append((int(start), int(end)))
    record["sources"].add(Path(path).name)
    record["attrs"].update({k: v for k, v in attributes.items() if v})


def finalize_gtf_records(records):
    for record in records.values():
        record["exons"] = sorted(
            set(record["exons"]),
            reverse=(record["strand"] == "-"),
        )
    return records


def parse_sample_gtfs(paths):
    records = {}
    for path in paths:
        sample = sample_from_path(path)
        with open_text(path) as handle:
            for line in handle:
                if not line or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 9 or fields[2] != "exon":
                    continue
                attributes = parse_gtf_attributes(fields[8])
                transcript = (
                    attributes.get("transcript_id")
                    or attributes.get("transcriptId")
                    or attributes.get("ID")
                )
                if not transcript:
                    continue
                add_gtf_record(
                    records,
                    (sample, transcript),
                    sample,
                    path,
                    fields,
                    attributes,
                )
    return finalize_gtf_records(records)


def parse_reference_gtf(path):
    records = {}
    with open_text(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "exon":
                continue
            attributes = parse_gtf_attributes(fields[8])
            transcript = (
                attributes.get("transcript_id")
                or attributes.get("transcriptId")
                or attributes.get("ID")
            )
            if not transcript:
                continue
            add_gtf_record(
                records,
                transcript,
                "REFERENCE",
                path,
                fields,
                attributes,
            )
    return finalize_gtf_records(records)


def build_reference_junctions(reference_transcripts):
    junctions = set()
    for record in reference_transcripts.values():
        genomic_exons = sorted(record["exons"])
        for left, right in zip(genomic_exons, genomic_exons[1:]):
            junctions.add(
                (
                    record["chrom"].removeprefix("chr"),
                    left[1],
                    right[0],
                    record["strand"],
                )
            )
    return junctions


def parse_fasta(path):
    header = None
    sequence_parts = []
    with open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence_parts)
                header = line[1:].strip()
                sequence_parts = []
            else:
                if header is None:
                    raise ValueError(f"Sequence before FASTA header in {path}")
                sequence_parts.append(line)
    if header is not None:
        yield header, "".join(sequence_parts)


def parse_splice_header(header):
    token = header.split(None, 1)[0]
    if "|SPLICE|" not in token:
        return None

    sample = token.split("|", 1)[0].upper()
    protein_id = token.split("|SPLICE|", 1)[1]
    transcript_id = protein_id.rsplit(".p", 1)[0]
    coordinate = re.search(
        r"(?:^|\s)(STRG\.\d+\.\d+):(\d+)-(\d+)\(([+-])\)",
        header,
    )
    if not coordinate:
        return None

    coord_transcript, start, end, orientation = coordinate.groups()
    return {
        "token": token,
        "sample": sample,
        "protein_id": protein_id,
        "transcript_id": transcript_id,
        "coord_transcript": coord_transcript,
        "orf_start": int(start),
        "orf_end": int(end),
        "orf_orientation": orientation,
        "header": header,
    }


def load_splice_fastas(paths):
    by_token = defaultdict(list)
    by_sample_protein = defaultdict(list)

    for path in paths:
        filename_sample = sample_from_path(path)
        for header, sequence in parse_fasta(path):
            info = parse_splice_header(header)
            if info is None:
                continue
            if info["sample"] != filename_sample:
                raise ValueError(
                    f"FASTA sample mismatch in {path}: header says "
                    f"{info['sample']}, filename says {filename_sample}"
                )
            info["sequence"] = sequence.upper().rstrip("*")
            info["file"] = Path(path).name
            by_token[info["token"]].append(info)
            by_sample_protein[(info["sample"], info["protein_id"])].append(info)

    return by_token, by_sample_protein


def load_candidates(path):
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Sequence", "Novel headers"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Candidate table missing columns: {sorted(missing)}")
        return list(reader)


def transcript_layout(record):
    offset = 1
    layout = []
    for genomic_start, genomic_end in record["exons"]:
        length = genomic_end - genomic_start + 1
        layout.append(
            {
                "g_start": genomic_start,
                "g_end": genomic_end,
                "t_start": offset,
                "t_end": offset + length - 1,
            }
        )
        offset += length
    return layout


def genomic_junction_coordinates(record, left_exon, right_exon):
    if record["strand"] == "+":
        return left_exon["g_end"], right_exon["g_start"]
    return right_exon["g_end"], left_exon["g_start"]


def _all_positions(sequence, peptide):
    positions = []
    offset = 0
    while peptide:
        position = sequence.find(peptide, offset)
        if position < 0:
            break
        positions.append(position)
        offset = position + 1
    return positions


def find_peptide(sequence, peptide, il_equivalent):
    positions = _all_positions(sequence, peptide)
    mode = "exact"
    if not positions and il_equivalent:
        positions = _all_positions(
            normalize_sequence(sequence, True),
            normalize_sequence(peptide, True),
        )
        mode = "I/L-equivalent"
    if len(positions) == 1:
        return positions[0], mode
    if len(positions) > 1:
        return -2, "ambiguous-multiple-occurrences"
    return -1, "not-found"


def analyze_protein(
    candidate,
    protein,
    transcript_record,
    reference_junctions,
    il_equivalent,
):
    peptide = candidate["Sequence"].upper()
    position, match_mode = find_peptide(
        protein["sequence"], peptide, il_equivalent
    )
    if position == -2:
        return [], "peptide has multiple occurrences in splice protein"
    if position < 0:
        return [], "peptide not found in splice protein"
    if transcript_record is None:
        return [], "sample-specific transcript not found in supplied GTF"

    peptide_start_aa = position + 1
    peptide_end_aa = position + len(peptide)
    layout = transcript_layout(transcript_record)

    orf_start, orf_end = sorted((protein["orf_start"], protein["orf_end"]))
    expected_nt = len(protein["sequence"]) * 3
    observed_nt = orf_end - orf_start + 1
    if abs(observed_nt - expected_nt) <= 3:
        orf_length_status = "consistent"
    else:
        orf_length_status = f"mismatch:{observed_nt}nt_vs_{expected_nt}nt"

    rows = []
    for junction_index in range(len(layout) - 1):
        left_exon = layout[junction_index]
        right_exon = layout[junction_index + 1]
        boundary_after_nt = left_exon["t_end"]

        if not (orf_start <= boundary_after_nt < orf_end):
            continue

        coding_nt_left = boundary_after_nt - orf_start + 1
        phase = coding_nt_left % 3

        if phase == 0:
            junction_aa_left = coding_nt_left // 3
            junction_aa_right = junction_aa_left + 1
            crosses = (
                peptide_start_aa <= junction_aa_left
                and peptide_end_aa >= junction_aa_right
            )
            left_anchor = (
                max(0, junction_aa_left - peptide_start_aa + 1)
                if crosses
                else 0
            )
            right_anchor = (
                max(0, peptide_end_aa - junction_aa_right + 1)
                if crosses
                else 0
            )
            boundary_type = "between-codons"
        else:
            junction_aa = coding_nt_left // 3 + 1
            junction_aa_left = junction_aa
            junction_aa_right = junction_aa
            crosses = peptide_start_aa <= junction_aa <= peptide_end_aa
            left_anchor = (
                max(0, junction_aa - peptide_start_aa) if crosses else 0
            )
            right_anchor = (
                max(0, peptide_end_aa - junction_aa) if crosses else 0
            )
            boundary_type = f"within-codon-phase-{phase}"

        donor, acceptor = genomic_junction_coordinates(
            transcript_record, left_exon, right_exon
        )
        junction_key = (
            transcript_record["chrom"].removeprefix("chr"),
            donor,
            acceptor,
            transcript_record["strand"],
        )
        reference_status = (
            "known" if junction_key in reference_junctions else "novel"
        )
        classification = (
            "exact-junction-spanning"
            if crosses
            else "same-ORF-not-this-junction"
        )
        if not crosses:
            support = "not-spanning"
        elif min(left_anchor, right_anchor) < 2:
            support = "weak-anchor"
        else:
            support = "supported"

        rows.append(
            {
                "Sequence": peptide,
                "Observed pattern": candidate.get("Observed pattern", ""),
                "Observed samples": candidate.get("Observed samples", ""),
                "PEP": candidate.get("PEP", ""),
                "Score": candidate.get("Score", ""),
                "MS/MS Count": candidate.get("MS/MS Count", ""),
                "Canonical match count": candidate.get(
                    "Canonical match count", ""
                ),
                "Sample": protein["sample"],
                "Splice protein": protein["protein_id"],
                "Transcript": protein["transcript_id"],
                "FASTA file": protein["file"],
                "Chromosome": transcript_record["chrom"],
                "Genomic strand": transcript_record["strand"],
                "TransDecoder ORF orientation": protein["orf_orientation"],
                "ORF cDNA start": orf_start,
                "ORF cDNA end": orf_end,
                "ORF length status": orf_length_status,
                "Protein length": len(protein["sequence"]),
                "Peptide start AA": peptide_start_aa,
                "Peptide end AA": peptide_end_aa,
                "Peptide match mode": match_mode,
                "Junction number": junction_index + 1,
                "Junction cDNA boundary after nt": boundary_after_nt,
                "Junction AA left": junction_aa_left,
                "Junction AA right": junction_aa_right,
                "Junction boundary type": boundary_type,
                "Left genomic exon": (
                    f"{left_exon['g_start']}-{left_exon['g_end']}"
                ),
                "Right genomic exon": (
                    f"{right_exon['g_start']}-{right_exon['g_end']}"
                ),
                "Genomic junction": (
                    f"{transcript_record['chrom']}:{donor}-{acceptor}"
                    f"({transcript_record['strand']})"
                ),
                "Reference junction status": reference_status,
                "Peptide classification": classification,
                "Anchor residues left": left_anchor,
                "Anchor residues right": right_anchor,
                "Junction support": support,
                "GTF class_code": transcript_record["attrs"].get(
                    "class_code", ""
                ),
                "GTF reference_id": transcript_record["attrs"].get(
                    "ref_id",
                    transcript_record["attrs"].get("reference_id", ""),
                ),
                "GTF gene_id": transcript_record["attrs"].get(
                    "gene_id", ""
                ),
                "GTF sources": ";".join(
                    sorted(transcript_record["sources"])
                ),
                "Splice FASTA header": protein["header"],
            }
        )

    if not rows:
        return [], "ORF contains no internal exon boundary"
    return rows, ""


def write_tsv(path, rows, fieldnames):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def output_path(prefix, suffix):
    return Path(f"{prefix}{suffix}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Validate whether candidate splice peptides cross translated "
            "sample-specific exon junctions and whether those junctions are "
            "absent from the matching Ensembl reference annotation."
        )
    )
    parser.add_argument(
        "--candidates",
        required=True,
        help="junction_peptide_analysis.splice_candidates.tsv",
    )
    parser.add_argument("--splice-fasta", nargs="+", required=True)
    parser.add_argument(
        "--transcript-gtf",
        nargs="+",
        required=True,
        help="Sample-specific StringTie assembled GTF files",
    )
    parser.add_argument(
        "--reference-gtf",
        required=True,
        help="Matching Ensembl reference GTF",
    )
    parser.add_argument(
        "--output-prefix",
        default="validated_splice_junctions_v3",
    )
    parser.add_argument("--exact-il", action="store_true")
    args = parser.parse_args()

    il_equivalent = not args.exact_il
    candidates = load_candidates(args.candidates)
    proteins_by_token, proteins_by_sample_id = load_splice_fastas(
        args.splice_fasta
    )
    sample_transcripts = parse_sample_gtfs(args.transcript_gtf)
    reference_transcripts = parse_reference_gtf(args.reference_gtf)
    reference_junctions = build_reference_junctions(reference_transcripts)

    detailed_rows = []
    unresolved_rows = []
    seen_entries = set()

    for candidate in candidates:
        raw_headers = candidate.get("Novel headers", "")
        tokens = []
        for raw_header in raw_headers.split(";"):
            raw_header = raw_header.strip()
            if "|SPLICE|" not in raw_header:
                continue
            token = raw_header.split(None, 1)[0]
            if token not in tokens:
                tokens.append(token)

        if not tokens:
            unresolved_rows.append(
                {
                    "Sequence": candidate["Sequence"],
                    "Splice header": "",
                    "Reason": "candidate row contains no splice FASTA header",
                }
            )
            continue

        for token in tokens:
            entries = proteins_by_token.get(token, [])
            if not entries:
                sample = token.split("|", 1)[0].upper()
                protein_id = token.split("|SPLICE|", 1)[-1]
                entries = proteins_by_sample_id.get((sample, protein_id), [])

            if not entries:
                unresolved_rows.append(
                    {
                        "Sequence": candidate["Sequence"],
                        "Splice header": token,
                        "Reason": "header not found in supplied splice FASTAs",
                    }
                )
                continue

            for protein in entries:
                marker = (
                    candidate["Sequence"],
                    protein["sample"],
                    protein["protein_id"],
                    protein["file"],
                )
                if marker in seen_entries:
                    continue
                seen_entries.add(marker)

                transcript_record = sample_transcripts.get(
                    (protein["sample"], protein["transcript_id"])
                ) or sample_transcripts.get(
                    (protein["sample"], protein["coord_transcript"])
                )

                rows, reason = analyze_protein(
                    candidate,
                    protein,
                    transcript_record,
                    reference_junctions,
                    il_equivalent,
                )
                detailed_rows.extend(rows)
                if reason:
                    unresolved_rows.append(
                        {
                            "Sequence": candidate["Sequence"],
                            "Splice header": protein["token"],
                            "Reason": reason,
                        }
                    )

    fieldnames = [
        "Sequence",
        "Observed pattern",
        "Observed samples",
        "PEP",
        "Score",
        "MS/MS Count",
        "Canonical match count",
        "Sample",
        "Splice protein",
        "Transcript",
        "FASTA file",
        "Chromosome",
        "Genomic strand",
        "TransDecoder ORF orientation",
        "ORF cDNA start",
        "ORF cDNA end",
        "ORF length status",
        "Protein length",
        "Peptide start AA",
        "Peptide end AA",
        "Peptide match mode",
        "Junction number",
        "Junction cDNA boundary after nt",
        "Junction AA left",
        "Junction AA right",
        "Junction boundary type",
        "Left genomic exon",
        "Right genomic exon",
        "Genomic junction",
        "Reference junction status",
        "Peptide classification",
        "Anchor residues left",
        "Anchor residues right",
        "Junction support",
        "GTF class_code",
        "GTF reference_id",
        "GTF gene_id",
        "GTF sources",
        "Splice FASTA header",
    ]

    junction_spanning_rows = [
        row
        for row in detailed_rows
        if row["Peptide classification"] == "exact-junction-spanning"
    ]
    novel_junction_rows = [
        row
        for row in junction_spanning_rows
        if row["Reference junction status"] == "novel"
    ]

    unique_prioritized = {}
    for row in novel_junction_rows:
        key = (
            row["Sequence"],
            row["Sample"],
            row["Genomic junction"],
        )
        old = unique_prioritized.get(key)
        balance = min(
            int(row["Anchor residues left"]),
            int(row["Anchor residues right"]),
        )
        old_balance = (
            -1
            if old is None
            else min(
                int(old["Anchor residues left"]),
                int(old["Anchor residues right"]),
            )
        )
        if old is None or balance > old_balance:
            unique_prioritized[key] = row

    prioritized_rows = sorted(
        unique_prioritized.values(),
        key=lambda row: (
            float(row["PEP"] or 1),
            -float(row["MS/MS Count"] or 0),
            row["Sequence"],
            row["Sample"],
        ),
    )

    prefix = args.output_prefix
    detailed_path = output_path(prefix, ".detailed.tsv")
    spanning_path = output_path(prefix, ".junction_spanning.tsv")
    prioritized_path = output_path(
        prefix, ".prioritized_novel_junctions.tsv"
    )
    unresolved_path = output_path(prefix, ".unresolved.tsv")
    summary_path = output_path(prefix, ".summary.txt")

    write_tsv(detailed_path, detailed_rows, fieldnames)
    write_tsv(spanning_path, junction_spanning_rows, fieldnames)
    write_tsv(prioritized_path, prioritized_rows, fieldnames)
    write_tsv(
        unresolved_path,
        unresolved_rows,
        ["Sequence", "Splice header", "Reason"],
    )

    unique_analyzed = {row["Sequence"] for row in detailed_rows}
    unique_spanning = {row["Sequence"] for row in junction_spanning_rows}
    unique_novel = {row["Sequence"] for row in novel_junction_rows}

    with summary_path.open("w", encoding="utf-8") as handle:
        handle.write(f"Candidate peptide rows: {len(candidates)}\n")
        handle.write(
            "Splice FASTA proteins loaded: "
            f"{sum(len(v) for v in proteins_by_token.values())}\n"
        )
        handle.write(
            f"Sample-specific transcript models loaded: "
            f"{len(sample_transcripts)}\n"
        )
        handle.write(
            f"Reference transcript models loaded: "
            f"{len(reference_transcripts)}\n"
        )
        handle.write(
            f"Detailed peptide-junction rows: {len(detailed_rows)}\n"
        )
        handle.write(
            f"Unique candidate peptides analyzed: {len(unique_analyzed)}\n"
        )
        handle.write(
            "Unique peptides spanning any translated exon junction: "
            f"{len(unique_spanning)}\n"
        )
        handle.write(
            "Unique peptides spanning a reference-absent junction: "
            f"{len(unique_novel)}\n"
        )
        handle.write(
            "Prioritized sample-specific peptide/junction events: "
            f"{len(prioritized_rows)}\n"
        )
        handle.write(f"Unresolved entries: {len(unresolved_rows)}\n")
        handle.write("\nPrioritized events:\n")
        for row in prioritized_rows:
            handle.write(
                f"  {row['Sequence']}\t{row['Observed pattern']}\t"
                f"{row['Sample']}\t{row['Genomic junction']}\t"
                f"anchors={row['Anchor residues left']}+"
                f"{row['Anchor residues right']}\tPEP={row['PEP']}\t"
                f"MSMS={row['MS/MS Count']}\n"
            )

    for path in (
        detailed_path,
        spanning_path,
        prioritized_path,
        unresolved_path,
        summary_path,
    ):
        print(f"Wrote {path}")
    print(f"Prioritized novel junction events: {len(prioritized_rows)}")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
