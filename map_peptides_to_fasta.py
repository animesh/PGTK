
import argparse
import csv
import gzip
import re
import sys
from collections import Counter, defaultdict, deque
from pathlib import Path


def open_text(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def natural_key(value):
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", str(value))]


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
                    raise ValueError(f"Sequence before first FASTA header in {path}")
                sequence_parts.append(line)
    if header is not None:
        yield header, "".join(sequence_parts)


def classify_header(header):
    upper = header.upper()
    if re.search(r"(?:^|[| ])PROGRESSION(?:[| ])", upper):
        return "PROGRESSION"
    if "|FUSION|" in upper:
        return "FUSION"
    if "|SPLICE|" in upper:
        return "SPLICE"
    if "|VAR_" in upper:
        return "VARIANT"
    if upper.startswith("CON__") or " CONTAMINANT" in upper:
        return "CONTAMINANT"
    if upper.startswith("REV__") or upper.startswith("DECOY_"):
        return "DECOY"
    return "CANONICAL"


def sample_from_header(header):
    match = re.match(r"^([^|]+)\|", header)
    return match.group(1).upper() if match else ""


def extract_accession(header):
    token = header.split(None, 1)[0]
    fields = token.split("|")
    if len(fields) >= 3 and fields[0] in {"sp", "tr"}:
        return fields[1]
    return token


def build_automaton(patterns):
    transitions = [{}]
    failure = [0]
    outputs = [[]]
    for pattern_id, pattern in enumerate(patterns):
        state = 0
        for character in pattern:
            next_state = transitions[state].get(character)
            if next_state is None:
                next_state = len(transitions)
                transitions[state][character] = next_state
                transitions.append({})
                failure.append(0)
                outputs.append([])
            state = next_state
        outputs[state].append(pattern_id)

    queue = deque()
    for state in transitions[0].values():
        queue.append(state)
    while queue:
        state = queue.popleft()
        for character, child in transitions[state].items():
            queue.append(child)
            fallback = failure[state]
            while fallback and character not in transitions[fallback]:
                fallback = failure[fallback]
            failure[child] = transitions[fallback].get(character, 0)
            outputs[child].extend(outputs[failure[child]])
    return transitions, failure, outputs


def find_matches(sequence, automaton, pattern_lengths):
    transitions, failure, outputs = automaton
    state = 0
    found = {}
    for index, character in enumerate(sequence):
        while state and character not in transitions[state]:
            state = failure[state]
        state = transitions[state].get(character, 0)
        for pattern_id in outputs[state]:
            found.setdefault(pattern_id, index - pattern_lengths[pattern_id] + 2)
    return found


def parse_group_map(items):
    mapping = {}
    for item in items:
        if "=" not in item:
            raise ValueError(f"Invalid --group-map value: {item}; expected PREFIX=SAMPLE")
        prefix, sample = item.split("=", 1)
        mapping[prefix.strip()] = sample.strip()
    return mapping


def discover_fastas(inputs):
    files = []
    accepted = {".fa", ".faa", ".fasta", ".fas", ".fa.gz", ".faa.gz", ".fasta.gz", ".fas.gz"}
    for item in inputs:
        path = Path(item)
        if path.is_dir():
            for candidate in path.rglob("*"):
                name = candidate.name.lower()
                if candidate.is_file() and any(name.endswith(extension) for extension in accepted):
                    files.append(candidate)
        elif path.is_file():
            files.append(path)
        else:
            raise FileNotFoundError(f"FASTA input does not exist: {path}")
    unique = sorted(set(path.resolve() for path in files), key=natural_key)
    if not unique:
        raise FileNotFoundError("No FASTA files were found")
    return unique


def read_peptides(path, il_equivalent, min_length):
    rows = []
    normalized_to_indices = defaultdict(list)
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "Sequence" not in reader.fieldnames:
            raise ValueError("peptides.txt must contain a Sequence column")
        for row_index, row in enumerate(reader):
            sequence = row.get("Sequence", "").strip().upper()
            normalized = normalize_sequence(sequence, il_equivalent)
            if len(normalized) < min_length:
                continue
            row["_row_index"] = row_index
            row["_sequence"] = sequence
            row["_normalized"] = normalized
            rows.append(row)
            normalized_to_indices[normalized].append(len(rows) - 1)
    if not rows:
        raise ValueError("No peptides remained after filtering")
    return rows, normalized_to_indices


def numeric(value):
    try:
        return float(value) if value not in (None, "") else 0.0
    except ValueError:
        return 0.0


def main():
    parser = argparse.ArgumentParser(
        description="Map MaxQuant peptides back to every searched FASTA and report canonical/novel uniqueness and sample patterns."
    )
    parser.add_argument("--peptides", required=True, help="MaxQuant peptides.txt")
    parser.add_argument("--fasta", nargs="+", required=True, help="FASTA files and/or directories")
    parser.add_argument("--output-prefix", default="peptide_fasta_mapping", help="Output prefix")
    parser.add_argument("--group-map", action="append", default=[], metavar="PREFIX=SAMPLE", help="Map MaxQuant experiment prefix to biological sample, e.g. runA=sampleA")
    parser.add_argument("--min-length", type=int, default=7)
    parser.add_argument("--exact-il", action="store_true", help="Treat I and L as distinct. Default treats them as mass-spectrometry equivalent.")
    parser.add_argument("--max-headers", type=int, default=0, help="Maximum headers stored per peptide; 0 means all")
    args = parser.parse_args()

    group_map = parse_group_map(args.group_map)
    il_equivalent = not args.exact_il
    fasta_files = discover_fastas(args.fasta)
    peptide_rows, normalized_to_indices = read_peptides(args.peptides, il_equivalent, args.min_length)

    patterns = sorted(normalized_to_indices, key=lambda value: (-len(value), value))
    pattern_to_id = {pattern: index for index, pattern in enumerate(patterns)}
    pattern_lengths = [len(pattern) for pattern in patterns]
    automaton = build_automaton(patterns)

    matches = [dict(headers=[], class_counts=Counter(), sample_counts=Counter(), fasta_counts=Counter(), positions=[]) for _ in peptide_rows]
    proteins_scanned = 0

    for fasta_path in fasta_files:
        fasta_label = fasta_path.name
        for header, raw_sequence in parse_fasta(fasta_path):
            proteins_scanned += 1
            sequence = normalize_sequence(raw_sequence, il_equivalent)
            found = find_matches(sequence, automaton, pattern_lengths)
            if not found:
                continue
            sequence_class = classify_header(header)
            sequence_sample = sample_from_header(header)
            accession = extract_accession(header)
            for pattern_id, start in found.items():
                pattern = patterns[pattern_id]
                for row_index in normalized_to_indices[pattern]:
                    result = matches[row_index]
                    result["class_counts"][sequence_class] += 1
                    if sequence_sample:
                        result["sample_counts"][sequence_sample] += 1
                    result["fasta_counts"][fasta_label] += 1
                    if args.max_headers == 0 or len(result["headers"]) < args.max_headers:
                        result["headers"].append(header)
                        result["positions"].append(f"{accession}:{start}-{start + len(pattern) - 1}")

    fieldnames = []
    with open_text(args.peptides) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        original_fields = reader.fieldnames or []
        intensity_columns = [name for name in original_fields if name.startswith("Intensity ") and name != "Intensity"]
        experiment_columns = [name for name in original_fields if name.startswith("Experiment ")]

    output_prefix = Path(args.output_prefix)
    mapping_path = output_prefix.with_suffix(".mapping.tsv")
    candidate_path = output_prefix.with_suffix(".candidates.tsv")
    summary_path = output_prefix.with_suffix(".summary.txt")

    output_fields = [
        "Sequence", "Length", "PEP", "Score", "MS/MS Count", "Leading razor protein",
        "Canonical matches", "Variant matches", "Splice matches", "Fusion matches", "Progression matches",
        "Contaminant matches", "Decoy matches", "Other matches", "Total FASTA matches",
        "Absent from canonical FASTA", "Novel classes", "Novel samples in FASTA",
        "Observed samples", "Observed pattern", "FASTA files", "Match positions", "All matching FASTA headers"
    ]

    candidate_rows = []
    pattern_counter = Counter()
    class_counter = Counter()

    with mapping_path.open("w", encoding="utf-8", newline="") as output_handle:
        writer = csv.DictWriter(output_handle, fieldnames=output_fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row, result in zip(peptide_rows, matches):
            class_counts = result["class_counts"]
            total_matches = sum(class_counts.values())
            novel_classes = [name for name in ("VARIANT", "SPLICE", "FUSION", "PROGRESSION") if class_counts[name] > 0]
            observed_samples = []
            observed_values = {}
            for prefix, sample in group_map.items():
                intensity = sum(numeric(row.get(column)) for column in intensity_columns if column.startswith(f"Intensity {prefix}_"))
                evidence = sum(1 for column in experiment_columns if column.startswith(f"Experiment {prefix}_") and row.get(column, "").strip())
                observed_values[sample] = intensity
                if intensity > 0 or evidence > 0:
                    observed_samples.append(sample)
            expected_order = [sample for _, sample in group_map.items()]
            if len(expected_order) == 3:
                presence = [sample in observed_samples for sample in expected_order]
                if presence == [True, False, False]: observed_pattern = f"{expected_order[0]}-only"
                elif presence == [False, True, False]: observed_pattern = f"{expected_order[1]}-only"
                elif presence == [False, False, True]: observed_pattern = f"{expected_order[2]}-only"
                elif presence == [False, True, True]: observed_pattern = f"{expected_order[1]}+{expected_order[2]}"
                elif presence == [True, True, True]: observed_pattern = "all-three"
                else: observed_pattern = "+".join(observed_samples) if observed_samples else "not-observed"
            else:
                observed_pattern = "+".join(observed_samples) if observed_samples else "not-evaluated"

            output_row = {
                "Sequence": row["_sequence"],
                "Length": len(row["_sequence"]),
                "PEP": row.get("PEP", ""),
                "Score": row.get("Score", ""),
                "MS/MS Count": row.get("MS/MS Count", ""),
                "Leading razor protein": row.get("Leading razor protein", ""),
                "Canonical matches": class_counts["CANONICAL"],
                "Variant matches": class_counts["VARIANT"],
                "Splice matches": class_counts["SPLICE"],
                "Fusion matches": class_counts["FUSION"],
                "Progression matches": class_counts["PROGRESSION"],
                "Contaminant matches": class_counts["CONTAMINANT"],
                "Decoy matches": class_counts["DECOY"],
                "Other matches": class_counts["OTHER"],
                "Total FASTA matches": total_matches,
                "Absent from canonical FASTA": "yes" if class_counts["CANONICAL"] == 0 else "no",
                "Novel classes": ";".join(novel_classes),
                "Novel samples in FASTA": ";".join(sorted(result["sample_counts"], key=natural_key)),
                "Observed samples": ";".join(observed_samples),
                "Observed pattern": observed_pattern,
                "FASTA files": ";".join(f"{name}:{count}" for name, count in sorted(result["fasta_counts"].items(), key=lambda item: natural_key(item[0]))),
                "Match positions": ";".join(result["positions"]),
                "All matching FASTA headers": ";".join(result["headers"]),
            }
            writer.writerow(output_row)

            pattern_counter[observed_pattern] += 1
            for name in novel_classes:
                class_counter[name] += 1

            high_confidence = (
                class_counts["CANONICAL"] == 0
                and bool(novel_classes)
                and numeric(row.get("PEP")) <= 0.01
                and numeric(row.get("Score")) >= 40
                and numeric(row.get("MS/MS Count")) >= 2
                and len(row["_sequence"]) >= 8
            )
            if high_confidence:
                candidate_rows.append(output_row)

    with candidate_path.open("w", encoding="utf-8", newline="") as output_handle:
        writer = csv.DictWriter(output_handle, fieldnames=output_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(sorted(candidate_rows, key=lambda row: (row["Observed pattern"], numeric(row["PEP"]), -numeric(row["MS/MS Count"]))))

    with summary_path.open("w", encoding="utf-8") as output_handle:
        output_handle.write(f"Peptides analyzed: {len(peptide_rows)}\n")
        output_handle.write(f"Distinct peptide patterns: {len(patterns)}\n")
        output_handle.write(f"FASTA files scanned: {len(fasta_files)}\n")
        output_handle.write(f"Protein sequences scanned: {proteins_scanned}\n")
        output_handle.write(f"I/L treated as equivalent: {'yes' if il_equivalent else 'no'}\n")
        output_handle.write(f"High-confidence noncanonical candidates: {len(candidate_rows)}\n")
        output_handle.write("\nFASTA files:\n")
        for fasta_path in fasta_files:
            output_handle.write(f"  {fasta_path}\n")
        output_handle.write("\nObserved patterns:\n")
        for name, count in sorted(pattern_counter.items(), key=lambda item: (-item[1], item[0])):
            output_handle.write(f"  {name}: {count}\n")
        output_handle.write("\nNovel classes with peptide matches:\n")
        for name, count in sorted(class_counter.items()):
            output_handle.write(f"  {name}: {count}\n")

    print(f"Wrote {mapping_path}")
    print(f"Wrote {candidate_path}")
    print(f"Wrote {summary_path}")
    print(f"High-confidence noncanonical candidates: {len(candidate_rows)}")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
