#!/usr/bin/env python3
from pathlib import Path
import ast
root=Path(__file__).parent
for name in ['build_pgtk_multiqc_content.py','build_expression_multiqc_content.py','build_final_multiqc_content.py']:
    ast.parse((root/name).read_text())
main=(root/'main.nf').read_text(); config=(root/'multiqc_config.yaml').read_text()
assert '<pre style="white-space:pre-wrap">' not in main
assert 'build_final_multiqc_content.py' in main
assert 'llms-full.txt' in main and '25000000' in main
assert '"_R1"' not in config and '"_R2"' not in config and '".trimmed"' not in config
assert main.count('process PREPARE_FINAL_MULTIQC_CONTENT {')==1
assert main.count('{')==main.count('}')
print('PASS: compact MultiQC design, report-size guard, and QC sample-name preservation')
