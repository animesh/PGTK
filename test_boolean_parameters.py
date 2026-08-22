#!/usr/bin/env python3
from pathlib import Path
source = Path("main.nf").read_text(encoding="utf-8")
booleans = [
    "run_proteogenomic_validation", "hc_dont_use_soft_clipped_bases",
    "generate_priority_igv_reports", "run_external_vcf_comparison",
    "run_expression_go", "gene_count_count_read_pairs",
    "gene_count_require_both_ends", "gene_count_exclude_chimeric",
    "gene_count_primary_only", "gene_count_allow_multi_overlap",
    "gene_count_count_multimapping",
]
assert "def strictBooleanParam(value, String name)" in source
assert "['false','0','no','n','off']" in source
for name in booleans:
    uses = [line for line in source.splitlines() if f"params.{name}" in line and not line.startswith(f"params.{name} =")]
    assert uses, name
    for line in uses:
        assert "strictBooleanParam" in line or "booleanText" in line, (name, line)
assert "if (params.run_proteogenomic_validation)" not in source
assert "if (params.run_external_vcf_comparison)" not in source
print("PASS: all pipeline boolean parameters use strict value parsing")
