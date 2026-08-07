#!/usr/bin/env python3
import argparse
import csv
import os
import sys
from pathlib import Path

REFERENCE_FILES = (
    'Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz',
    'Homo_sapiens.GRCh38.111.gtf.gz',
    'Homo_sapiens.GRCh38.cdna.all.fa.gz',
    'human_reviewed_isoforms.fasta.gz',
    'homo_sapiens_vep_111_GRCh38.tar.gz',
    'arriba_v2.4.0.tar.gz',
)
LOCAL_CONTAINER_FILES = (
    'pvactools-7.1.1.img',
    'stringtie-3.0.3.img',
    'gffcompare-0.12.10.img',
    'transdecoder-6.0.0.img',
)
REQUIRED_SAMPLE_COLUMNS = ('sample', 'srr')
OPTIONAL_SAMPLE_COLUMNS = ('TK', 'Group', 'baseline')

def require_file(path, label):
    if not path.is_file() or path.stat().st_size <= 0:
        raise ValueError(f'{label} is missing or empty: {path}')

def require_directory(path, label, writable=False):
    if not path.is_dir():
        raise ValueError(f'{label} is missing: {path}')
    if not os.access(path, os.R_OK):
        raise ValueError(f'{label} is not readable: {path}')
    if writable and not os.access(path, os.W_OK):
        raise ValueError(f'{label} is not writable: {path}')

def validate_samples(samplesheet, sra_dir):
    require_file(samplesheet, 'samplesheet')
    with samplesheet.open(encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle)
        missing = [name for name in REQUIRED_SAMPLE_COLUMNS if name not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f'samplesheet lacks columns: {missing}')
        rows = list(reader)
    if not rows:
        raise ValueError('samplesheet contains no samples')
    seen = set()
    for line_number, row in enumerate(rows, 2):
        for name in REQUIRED_SAMPLE_COLUMNS:
            if not str(row.get(name, '')).strip():
                raise ValueError(f'samplesheet line {line_number} has empty {name}')
        sample = row['sample'].strip()
        srr = row['srr'].strip()
        if sample in seen:
            raise ValueError(f'duplicate sample in samplesheet: {sample}')
        seen.add(sample)
        baseline = str(row.get('baseline', '') or '').strip().lower()
        if baseline and baseline not in {'true', 'false'}:
            raise ValueError(f'samplesheet line {line_number} baseline must be true or false')
        require_file(sra_dir / srr / f'{srr}.sra', f'SRA input for {sample}')

def main():
    parser = argparse.ArgumentParser(description='Validate all runtime files before Nextflow submits tasks.')
    parser.add_argument('--project-dir', required=True, type=Path)
    parser.add_argument('--reference-downloads', required=True, type=Path)
    parser.add_argument('--container-cache', required=True, type=Path)
    parser.add_argument('--sra-dir', required=True, type=Path)
    parser.add_argument('--samplesheet', required=True, type=Path)
    parser.add_argument('--ensembl-pep', required=True, type=Path)
    parser.add_argument('--work-dir', required=True, type=Path)
    parser.add_argument('--tmp-root', required=True, type=Path)
    parser.add_argument('--results-dir', required=True, type=Path)
    args = parser.parse_args()

    require_directory(args.project_dir, 'project directory')
    require_directory(args.reference_downloads, 'reference directory')
    require_directory(args.container_cache, 'container cache')
    require_directory(args.sra_dir, 'SRA directory')
    require_directory(args.work_dir, 'work directory', writable=True)
    require_directory(args.tmp_root, 'temporary root', writable=True)
    require_directory(args.results_dir, 'results directory', writable=True)

    for name in REFERENCE_FILES:
        require_file(args.reference_downloads / name, 'reference asset')
    require_file(args.ensembl_pep, 'Ensembl peptide FASTA')
    for name in LOCAL_CONTAINER_FILES:
        require_file(args.container_cache / name, 'local container image')
    validate_samples(args.samplesheet, args.sra_dir)
    print('PASS: runtime references, containers, SRA inputs, samplesheet and writable directories')

if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        print(f'ERROR: {error}', file=sys.stderr)
        raise SystemExit(1)
