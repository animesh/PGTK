#!/usr/bin/env python3
import json
import subprocess
import tempfile
from pathlib import Path

root = Path(__file__).resolve().parent
main = (root / "main.nf").read_text(encoding="utf-8")
builder = (root / "build_compact_multiqc_content.py").read_text(encoding="utf-8")

for token in [
    "${meta.sample}.raw.R1.fastq.gz",
    "${meta.sample}.raw.R2.fastq.gz",
    "${meta.sample}.trimmed.R1.fastq.gz",
    "${meta.sample}.trimmed.R2.fastq.gz",
    "--evidence-classification ${evidence_classification}",
    "evidence_report.classification_report,",
]:
    assert token in main, token

for token in [
    '"Event class"', '"Data rows"',
    '"Canonical-and-reference-absent variant events"',
    '"Strict integrated events"',
    '"SAMPLE_MATCHED_DIRECT_MSMS"',
    'abs(number(row.get("MeanScore"))) < 0.25',
    'number(row.get("OddsRatio")) < 1.5',
]:
    assert token in builder, token

print("PASS: reporting redesign static contract")
