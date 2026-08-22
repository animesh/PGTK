#!/usr/bin/env python3
import ast
import re
from pathlib import Path

root = Path(__file__).resolve().parent
for path in root.glob('*.py'):
    ast.parse(path.read_text(encoding='utf-8'))

main = (root / 'main.nf').read_text(encoding='utf-8')
config = (root / 'nextflow.config').read_text(encoding='utf-8')
analysis = (root / 'expression_go_analysis.py').read_text(encoding='utf-8')
multiqc = (root / 'multiqc_config.yaml').read_text(encoding='utf-8')

assert main.count('{') == main.count('}')
assert 'show_hide_buttons:' not in multiqc and 'show_hide_mode:' not in multiqc
assert 'test "\\$(stat -c %s multiqc_report.html)"' in main
for line_number, line in enumerate(main.splitlines(), 1):
    if '$(' in line and '\\$(' not in line:
        raise AssertionError(f'unescaped shell substitution at {line_number}: {line}')

contracts = [
    ("args.output_prefix + '.expression_ora.tsv'", '${meta.sample}.expression_go.expression_ora.tsv'),
    ("args.output_prefix + '.ranked_go.tsv'", '${meta.sample}_vs_${baseline_sample}.expression_go.ranked_go.tsv'),
    ("args.output_prefix + '.expression_ora.tsv'", 'expression_go.expression_ora.tsv'),
    ("args.output_prefix + '.ranked_go.tsv'", 'expression_go.ranked_go.tsv'),
]
for script_suffix, nextflow_output in contracts:
    assert script_suffix in analysis, script_suffix
    assert nextflow_output in main, nextflow_output

for token in [
    "path 'expression_go.expression_ora.tsv', emit: ora",
    "path 'expression_go.ranked_go.tsv', emit: ranked",
    "path 'expression_go.summary.tsv', emit: summary",
    'expression_go.ora,',
    'expression_go.ranked,',
    'progression_variant_sets.enrichment,',
]:
    assert token in main, token

for forbidden in [
    '${meta.sample}.expression_go.ora.tsv',
    '${meta.sample}_vs_${baseline_sample}.expression_go.ranked.tsv',
    "path 'expression_go.ora.tsv', emit: ora",
    "path 'expression_go.ranked.tsv', emit: ranked",
    'expression_go.expression_ora,',
    'expression_go.ranked_go,',
]:
    assert forbidden not in main, forbidden

for observer in ['trace', 'timeline', 'report', 'dag']:
    block = re.search(rf'{observer}\s*\{{(?P<body>.*?)\n\}}', config, flags=re.DOTALL)
    assert block and re.search(r'overwrite\s*=\s*true', block.group('body')), observer

print('PASS: expression script filenames, process outputs, emitted channels, and observers are consistent')
