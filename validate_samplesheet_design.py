#!/usr/bin/env python3
"""Validate PGTK samplesheet structure and longitudinal baseline design."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--validated-samplesheet", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    with args.input.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        rows = list(reader)

    missing = {"sample", "srr"} - set(fields)
    if missing:
        raise SystemExit(f"samplesheet missing columns: {sorted(missing)}")
    if not rows:
        raise SystemExit("samplesheet contains no data rows")

    seen_samples: set[str] = set()
    seen_srrs: set[str] = set()
    subjects: dict[str, list[tuple[str, str]]] = defaultdict(list)

    for line_number, row in enumerate(rows, 2):
        sample = (row.get("sample") or "").strip()
        srr = (row.get("srr") or "").strip()
        subject = (row.get("TK") or sample).strip()
        group = (row.get("Group") or sample).strip()
        baseline = (row.get("baseline") or "false").strip().lower()

        if not sample or not srr or not subject or not group:
            raise SystemExit(f"empty required value at samplesheet line {line_number}")
        if sample in seen_samples:
            raise SystemExit(f"duplicate sample: {sample}")
        if srr in seen_srrs:
            raise SystemExit(f"duplicate srr: {srr}")
        if baseline not in {"true", "false"}:
            raise SystemExit(f"baseline must be true or false for {sample}")

        seen_samples.add(sample)
        seen_srrs.add(srr)
        subjects[subject].append((sample, baseline))

    with args.report.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["TK", "Samples", "BaselineSamples", "NonBaselineSamples", "SubtractionStatus"]
        )
        for subject in sorted(subjects):
            members = subjects[subject]
            baselines = [sample for sample, value in members if value == "true"]
            non_baselines = [sample for sample, value in members if value == "false"]
            if len(baselines) > 1:
                raise SystemExit(f"multiple baselines for {subject}: {baselines}")
            if len(baselines) == 1 and non_baselines:
                status = "ENABLED"
            elif len(baselines) == 1:
                status = "NO_NONBASELINE"
            else:
                status = "SKIPPED_NO_BASELINE"
            writer.writerow(
                [subject, ",".join(sample for sample, _ in members), ",".join(baselines), ",".join(non_baselines), status]
            )

    args.validated_samplesheet.write_bytes(args.input.read_bytes())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
