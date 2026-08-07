import argparse
import csv
import gzip
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path


def open_text(path):
    path = Path(path)
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def normalize_sequence(sequence, il_equivalent=True):
    sequence = re.sub(r"[^A-Za-z]", "", sequence).upper()
    if il_equivalent:
        sequence = sequence.replace("I", "J").replace("L", "J")
    return sequence


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
                    raise ValueError(f"Sequence found before first FASTA header in {path}")
                sequence_parts.append(line)
    if header is not None:
        yield header, "".join(sequence_parts)


def sample_from_header(header):
    match = re.match(r"^([^|]+)\|", header)
    return match.group(1).upper() if match else ""


def classify_header(header):
    upper = header.upper()
    if "|FUSION|" in upper:
        return "FUSION"
    if "|SPLICE|" in upper:
        return "SPLICE"
    return "OTHER"


def extract_gene_pair(header):
    if "|FUSION|" not in header:
        return "", ""
    token = header.split(None, 1)[0]
    body = token.split("|FUSION|", 1)[1]
    fields = body.split(".")
    candidate = fields[1] if len(fields) > 1 else fields[0]
    if "-" not in candidate:
        return candidate, ""
    left, right = candidate.split("-", 1)
    left = re.sub(r"\(.*?\)", "", left)
    right = re.sub(r"\(.*?\)", "", right)
    return left, right


def parse_group_map(items):
    mapping = []
    for item in items:
        if "=" not in item:
            raise ValueError(f"Invalid --group-map value {item!r}; expected PREFIX=SAMPLE")
        prefix, sample = item.split("=", 1)
        mapping.append((prefix.strip(), sample.strip()))
    return mapping


def numeric(value):
    try:
        return float(value) if value not in (None, "") else 0.0
    except ValueError:
        return 0.0


def load_peptides(path, group_map):
    rows = []
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        if "Sequence" not in fields:
            raise ValueError("peptides.txt must contain a Sequence column")
        intensity_columns = [field for field in fields if field.startswith("Intensity ") and field != "Intensity"]
        experiment_columns = [field for field in fields if field.startswith("Experiment ")]
        for row in reader:
            sequence = row.get("Sequence", "").strip().upper()
            if not sequence:
                continue
            observed_samples = []
            sample_intensities = {}
            for prefix, sample in group_map:
                intensity = sum(
                    numeric(row.get(column))
                    for column in intensity_columns
                    if column.startswith(f"Intensity {prefix}_")
                )
                evidence = sum(
                    1
                    for column in experiment_columns
                    if column.startswith(f"Experiment {prefix}_") and row.get(column, "").strip()
                )
                sample_intensities[sample] = intensity
                if intensity > 0 or evidence > 0:
                    observed_samples.append(sample)
            if len(group_map) == 3:
                expected = [sample for _, sample in group_map]
                present = [sample in observed_samples for sample in expected]
                if present == [True, False, False]:
                    pattern = f"{expected[0]}-only"
                elif present == [False, True, False]:
                    pattern = f"{expected[1]}-only"
                elif present == [False, False, True]:
                    pattern = f"{expected[2]}-only"
                elif present == [False, True, True]:
                    pattern = f"{expected[1]}+{expected[2]}"
                elif present == [True, True, True]:
                    pattern = "all-three"
                else:
                    pattern = "+".join(observed_samples) if observed_samples else "not-observed"
            else:
                pattern = "+".join(observed_samples) if observed_samples else "not-evaluated"
            rows.append(
                {
                    "Sequence": sequence,
                    "PEP": row.get("PEP", ""),
                    "Score": row.get("Score", ""),
                    "MS/MS Count": row.get("MS/MS Count", ""),
                    "Leading razor protein": row.get("Leading razor protein", ""),
                    "Observed samples": ";".join(observed_samples),
                    "Observed pattern": pattern,
                    "Sample intensities": ";".join(
                        f"{sample}:{sample_intensities.get(sample, 0):.0f}"
                        for _, sample in group_map
                    ),
                }
            )
    return rows


def load_proteins(paths, kind, il_equivalent):
    proteins = []
    for path in paths:
        for header, sequence in parse_fasta(path):
            if kind != "CANONICAL" and classify_header(header) != kind:
                continue
            proteins.append(
                {
                    "header": header,
                    "sequence": sequence.upper(),
                    "normalized": normalize_sequence(sequence, il_equivalent),
                    "sample": sample_from_header(header),
                    "file": Path(path).name,
                }
            )
    return proteins


def parse_arriba(paths):
    events = defaultdict(list)
    for path in paths:
        sample = Path(path).name.split(".", 1)[0].upper()
        with open_text(path) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames:
                continue
            normalized = {field.lstrip("#"): field for field in reader.fieldnames}
            gene1_field = normalized.get("gene1")
            gene2_field = normalized.get("gene2")
            if not gene1_field or not gene2_field:
                continue
            for row in reader:
                gene1 = row.get(gene1_field, "")
                gene2 = row.get(gene2_field, "")
                key = (sample, gene1.split("(", 1)[0], gene2.split("(", 1)[0])
                event = {
                    "breakpoint1": row.get(normalized.get("breakpoint1", ""), ""),
                    "breakpoint2": row.get(normalized.get("breakpoint2", ""), ""),
                    "confidence": row.get(normalized.get("confidence", ""), ""),
                    "type": row.get(normalized.get("type", ""), ""),
                    "reading_frame": row.get(normalized.get("reading_frame", ""), ""),
                    "split_reads1": row.get(normalized.get("split_reads1", ""), ""),
                    "split_reads2": row.get(normalized.get("split_reads2", ""), ""),
                    "discordant_mates": row.get(normalized.get("discordant_mates", ""), ""),
                    "file": Path(path).name,
                }
                events[key].append(event)
    return events


def peptide_matches(peptide_norm, proteins):
    matches = []
    for protein in proteins:
        offset = 0
        while peptide_norm:
            start = protein["normalized"].find(peptide_norm, offset)
            if start < 0:
                break
            matches.append((protein, start + 1))
            offset = start + 1
    return matches


def make_fragment_lookup(canonical_proteins):
    cache = {}

    def lookup(fragment):
        if fragment not in cache:
            cache[fragment] = [
                protein for protein in canonical_proteins if fragment in protein["normalized"]
            ]
        return cache[fragment]

    return lookup


def infer_split(peptide_norm, lookup, min_anchor):
    candidates = []
    for split in range(min_anchor, len(peptide_norm) - min_anchor + 1):
        left = peptide_norm[:split]
        right = peptide_norm[split:]
        left_hits = lookup(left)
        right_hits = lookup(right)
        if not left_hits or not right_hits:
            continue
        distinct_pair = None
        for left_hit in left_hits[:50]:
            for right_hit in right_hits[:50]:
                if left_hit["header"] != right_hit["header"]:
                    distinct_pair = (left_hit, right_hit)
                    break
            if distinct_pair:
                break
        if distinct_pair:
            candidates.append(
                {
                    "split": split,
                    "left": left,
                    "right": right,
                    "left_header": distinct_pair[0]["header"],
                    "right_header": distinct_pair[1]["header"],
                    "balance": min(len(left), len(right)),
                }
            )
    if not candidates:
        return None
    candidates.sort(key=lambda item: (-item["balance"], abs(len(peptide_norm) / 2 - item["split"])))
    return candidates[0]


def write_tsv(path, rows, fields):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Find canonical-absent peptides in fusion and splice FASTAs and infer candidate junction-spanning peptides."
    )
    parser.add_argument("--peptides", required=True)
    parser.add_argument("--canonical-fasta", nargs="+", required=True)
    parser.add_argument("--fusion-fasta", nargs="+", required=True)
    parser.add_argument("--splice-fasta", nargs="+", required=True)
    parser.add_argument("--arriba", nargs="*", default=[])
    parser.add_argument("--group-map", action="append", default=[])
    parser.add_argument("--output-prefix", default="junction_peptide_analysis")
    parser.add_argument("--min-anchor", type=int, default=4)
    parser.add_argument("--exact-il", action="store_true")
    parser.add_argument("--max-header-matches", type=int, default=100)
    args = parser.parse_args()

    il_equivalent = not args.exact_il
    group_map = parse_group_map(args.group_map)
    peptides = load_peptides(args.peptides, group_map)
    canonical = load_proteins(args.canonical_fasta, "CANONICAL", il_equivalent)
    fusion = load_proteins(args.fusion_fasta, "FUSION", il_equivalent)
    splice = load_proteins(args.splice_fasta, "SPLICE", il_equivalent)
    arriba_events = parse_arriba(args.arriba)
    lookup = make_fragment_lookup(canonical)

    common_fields = [
        "Sequence", "Class", "Observed pattern", "Observed samples", "Sample intensities",
        "PEP", "Score", "MS/MS Count", "Leading razor protein", "Canonical match count", "Canonical absent",
        "Novel match count", "FASTA samples", "FASTA files", "Novel headers",
        "Peptide positions", "Inferred junction split after peptide residue",
        "Left peptide segment", "Right peptide segment", "Left canonical support",
        "Right canonical support", "Junction inference", "Direct noncanonical evidence",
        "Arriba gene1", "Arriba gene2", "Arriba breakpoint1", "Arriba breakpoint2",
        "Arriba confidence", "Arriba type", "Arriba reading frame", "Arriba support",
    ]

    all_rows = []
    fusion_rows = []
    splice_rows = []
    for peptide in peptides:
        sequence = peptide["Sequence"]
        normalized = normalize_sequence(sequence, il_equivalent)
        canonical_hits = peptide_matches(normalized, canonical)
        for kind, proteins in (("FUSION", fusion), ("SPLICE", splice)):
            novel_hits = peptide_matches(normalized, proteins)
            if not novel_hits:
                continue
            split = None if canonical_hits else infer_split(normalized, lookup, args.min_anchor)
            headers = [match[0]["header"] for match in novel_hits[: args.max_header_matches]]
            positions = [
                f"{match[0]['header'].split(None, 1)[0]}:{match[1]}-{match[1] + len(sequence) - 1}"
                for match in novel_hits[: args.max_header_matches]
            ]
            samples = sorted({match[0]["sample"] for match in novel_hits if match[0]["sample"]})
            files = sorted({match[0]["file"] for match in novel_hits})
            direct = "yes" if not canonical_hits else "no"
            inference = "split-anchor-supported" if split else ("canonical-shared" if canonical_hits else "noncanonical-no-split-anchor")

            arriba_gene1 = arriba_gene2 = ""
            event_values = defaultdict(set)
            if kind == "FUSION":
                for match in novel_hits:
                    gene1, gene2 = extract_gene_pair(match[0]["header"])
                    if gene1 and not arriba_gene1:
                        arriba_gene1, arriba_gene2 = gene1, gene2
                    for event in arriba_events.get((match[0]["sample"], gene1, gene2), []):
                        for key, value in event.items():
                            if value:
                                event_values[key].add(value)

            output = {
                **peptide,
                "Class": kind,
                "Canonical match count": len(canonical_hits),
                "Canonical absent": "yes" if not canonical_hits else "no",
                "Novel match count": len(novel_hits),
                "FASTA samples": ";".join(samples),
                "FASTA files": ";".join(files),
                "Novel headers": ";".join(headers),
                "Peptide positions": ";".join(positions),
                "Inferred junction split after peptide residue": split["split"] if split else "",
                "Left peptide segment": sequence[: split["split"]] if split else "",
                "Right peptide segment": sequence[split["split"] :] if split else "",
                "Left canonical support": split["left_header"] if split else "",
                "Right canonical support": split["right_header"] if split else "",
                "Junction inference": inference,
                "Direct noncanonical evidence": direct,
                "Arriba gene1": arriba_gene1,
                "Arriba gene2": arriba_gene2,
                "Arriba breakpoint1": ";".join(sorted(event_values["breakpoint1"])),
                "Arriba breakpoint2": ";".join(sorted(event_values["breakpoint2"])),
                "Arriba confidence": ";".join(sorted(event_values["confidence"])),
                "Arriba type": ";".join(sorted(event_values["type"])),
                "Arriba reading frame": ";".join(sorted(event_values["reading_frame"])),
                "Arriba support": ";".join(
                    filter(
                        None,
                        [
                            "split1=" + ",".join(sorted(event_values["split_reads1"])) if event_values["split_reads1"] else "",
                            "split2=" + ",".join(sorted(event_values["split_reads2"])) if event_values["split_reads2"] else "",
                            "discordant=" + ",".join(sorted(event_values["discordant_mates"])) if event_values["discordant_mates"] else "",
                        ],
                    )
                ),
            }
            all_rows.append(output)
            if kind == "FUSION":
                fusion_rows.append(output)
            else:
                splice_rows.append(output)

    def is_high_confidence(row):
        return (
            row["Canonical absent"] == "yes"
            and numeric(row["PEP"]) <= 0.01
            and numeric(row["Score"]) >= 40
            and numeric(row["MS/MS Count"]) >= 2
            and len(row["Sequence"]) >= 8
        )

    fusion_candidates = [row for row in fusion_rows if is_high_confidence(row) and row["Junction inference"] == "split-anchor-supported" and row["Arriba breakpoint1"] and row["Arriba breakpoint2"]]
    splice_candidates = [row for row in splice_rows if is_high_confidence(row)]
    inferred_fusion = [row for row in fusion_candidates if row["Junction inference"] == "split-anchor-supported"]
    inferred_splice = [row for row in splice_candidates if row["Junction inference"] == "split-anchor-supported"]

    prefix = Path(args.output_prefix)
    write_tsv(prefix.with_suffix(".all_mappings.tsv"), all_rows, common_fields)
    write_tsv(prefix.with_suffix(".fusion_candidates.tsv"), fusion_candidates, common_fields)
    write_tsv(prefix.with_suffix(".splice_candidates.tsv"), splice_candidates, common_fields)
    write_tsv(prefix.with_suffix(".inferred_junctions.tsv"), inferred_fusion + inferred_splice, common_fields)

    summary_path = prefix.with_suffix(".summary.txt")
    with summary_path.open("w", encoding="utf-8") as handle:
        handle.write(f"Peptides analyzed: {len(peptides)}\n")
        handle.write(f"Canonical proteins scanned: {len(canonical)}\n")
        handle.write(f"Fusion proteins scanned: {len(fusion)}\n")
        handle.write(f"Splice proteins scanned: {len(splice)}\n")
        handle.write(f"I/L treated as equivalent: {'yes' if il_equivalent else 'no'}\n")
        handle.write(f"Peptides mapping to fusion proteins: {len(fusion_rows)}\n")
        handle.write(f"Peptides mapping to splice proteins: {len(splice_rows)}\n")
        handle.write(f"High-confidence canonical-absent fusion peptides: {len(fusion_candidates)}\n")
        handle.write(f"High-confidence canonical-absent splice peptides: {len(splice_candidates)}\n")
        handle.write(f"Fusion peptides with split-anchor junction inference: {len(inferred_fusion)}\n")
        handle.write(f"Splice peptides with split-anchor junction inference: {len(inferred_splice)}\n")
        handle.write("\nInterpretation:\n")
        handle.write("  canonical-absent = full peptide not found in supplied canonical FASTA\n")
        handle.write("  split-anchor-supported = peptide can be split into left and right segments supported by distinct canonical proteins\n")
        handle.write("  this is a computational junction inference, not definitive breakpoint proof\n")

    print(f"Wrote {prefix.with_suffix('.all_mappings.tsv')}")
    print(f"Wrote {prefix.with_suffix('.fusion_candidates.tsv')}")
    print(f"Wrote {prefix.with_suffix('.splice_candidates.tsv')}")
    print(f"Wrote {prefix.with_suffix('.inferred_junctions.tsv')}")
    print(f"Wrote {summary_path}")
    print(f"High-confidence fusion candidates: {len(fusion_candidates)}")
    print(f"High-confidence splice candidates: {len(splice_candidates)}")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
