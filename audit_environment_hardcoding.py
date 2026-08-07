#!/usr/bin/env python3
import argparse
import re
import sys
from pathlib import Path

TEXT_SUFFIXES = {'.nf', '.config', '.groovy', '.sh', '.slurm', '.py', '.json', '.yaml', '.yml'}
SKIP_NAMES = {'audit_environment_hardcoding.py'}
RULES = {
    'absolute cluster path': re.compile(r'(?<![A-Za-z0-9_])/(?:cluster|home|Users|mnt|scratch|work)(?:/|\b)'),
    'Saga project account': re.compile(r'\bnn\d{4,}\w*\b', re.I),
    'local username': re.compile(r'\bash\d{3,}\b', re.I),
    'embedded SBATCH directive': re.compile(r'^\s*#SBATCH\b', re.M),
    'institutional email': re.compile(r'\b[A-Z0-9._%+-]+@(?:ntnu\.no|uio\.no|uib\.no|nmbu\.no)\b', re.I),
    'fixed Java module': re.compile(r'\bmodule\s+load\s+Java/[A-Za-z0-9._-]+'),
}
ALLOW_FILE_PATTERNS = {
    'validate_pipeline_commands.sh': ('/cluster/', 'nn9036k', '#SBATCH', 'Java/21'),
    'test_resource_configuration.py': ('/cluster/', 'nn9036k', '#SBATCH', 'Java/21', '@ntnu.no'),
}

def main():
    parser = argparse.ArgumentParser(description='Reject environment-specific hard-coding in pipeline code.')
    parser.add_argument('project_dir', nargs='?', default='.')
    args = parser.parse_args()
    root = Path(args.project_dir).resolve()
    failures = []
    checked = 0
    active_files = [
        root / "main.nf",
        root / "nextflow.config",
        root / "collect_pipeline_failures.py",
        root / "analyze_pipeline_trace.py",
        root / "validate_rna_events.py",
        root / "build_complete_report.py",
        root / "map_peptides_to_fasta.py",
        root / "annotate_variant_peptides.py",
        root / "analyze_chimeric_splice_peptides.py",
        root / "validate_splice_junction_peptides.py",
        root / "proteogenomics_evidence_report.py",
        root / "validate_proteogenomic_reads.py",
        root / "validate_variant_read_provenance.py",
        root / "validate_variant_codons.py",
        root / "merge_variant_validation.py",
        root / "analyze_codon_mismatches.py",
        root / "build_integrated_variant_evidence.py",
        root / "validate_runtime_inputs.py",
        root / "validate_haplotype_shards.py",
        root / "summarize_variant_stages.py",
        root / "compare_external_vcf.py",
        root / "build_comparative_advantage_report.py",
        root / "build_igv_evidence_bundle.py",
    ]
    for path in active_files:
        if not path.is_file() or path.name in SKIP_NAMES or path.suffix not in TEXT_SUFFIXES:
            continue
        checked += 1
        text = path.read_text(encoding='utf-8', errors='replace')
        rel = path.relative_to(root).as_posix()
        scrubbed = text
        for literal in ALLOW_FILE_PATTERNS.get(path.name, ()):
            scrubbed = scrubbed.replace(literal, '')
        for label, pattern in RULES.items():
            for match in pattern.finditer(scrubbed):
                line = scrubbed.count('\n', 0, match.start()) + 1
                failures.append(f'{rel}:{line}: {label}: {match.group(0)!r}')
    if failures:
        print('\n'.join(failures), file=sys.stderr)
        print(f'FAIL: {len(failures)} environment-specific hard-coded value(s)', file=sys.stderr)
        return 1
    print(f'PASS: checked {checked} code/configuration files; no environment-specific hard-coding found')
    return 0

raise SystemExit(main())
