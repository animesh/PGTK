#!/usr/bin/env python3
import argparse
import csv
import os
import subprocess
import sys
from pathlib import Path

REFERENCE_FILES = (
    'Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz',
    'Homo_sapiens.GRCh38.111.gtf.gz',
    'Homo_sapiens.GRCh38.cdna.all.fa.gz',
    'human_reviewed_isoforms.fasta.gz',
    'homo_sapiens_vep_111_GRCh38.tar.gz',
    'arriba_v2.4.0.tar.gz',
    'go-basic.obo',
    'goa_human.gaf.gz',
)
LOCAL_CONTAINER_FILES = (
    'pvactools-7.1.1.img',
    'stringtie-3.0.3.img',
    'gffcompare-0.12.10.img',
    'transdecoder-6.0.0.img',
    'subread-2.0.8.img',
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
    parser.add_argument('--host-python', required=True, type=Path)
    parser.add_argument('--apptainer', required=True, type=Path)
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
    require_file(args.host_python, 'host Python executable')
    require_file(args.apptainer, 'Apptainer executable')
    if not os.access(args.host_python, os.X_OK):
        raise ValueError(f'host Python is not executable: {args.host_python}')
    if not os.access(args.apptainer, os.X_OK):
        raise ValueError(f'Apptainer is not executable: {args.apptainer}')
    subread_image = args.container_cache / 'subread-2.0.8.img'
    subprocess.run([str(args.apptainer), 'exec', str(subread_image), 'featureCounts', '-v'], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    samtools_images = sorted(args.container_cache.rglob('*samtools-1.21*img'))
    if len(samtools_images) != 1:
        raise ValueError(f'expected exactly one samtools 1.21 image below {args.container_cache}; found {len(samtools_images)}')
    multiqc_images = sorted(args.container_cache.rglob('*multiqc-1.35*img'))
    if len(multiqc_images) != 1:
        raise ValueError(f'expected exactly one MultiQC 1.35 image below {args.container_cache}; found {len(multiqc_images)}')
    host_python_real = args.host_python.resolve()
    host_python_root = host_python_real.parent.parent
    bind_paths = ','.join(dict.fromkeys(map(str, (
        args.project_dir,
        args.reference_downloads,
        args.container_cache,
        args.sra_dir,
        args.work_dir,
        args.tmp_root,
        args.results_dir,
        host_python_root,
    ))))
    common = [str(args.apptainer), 'exec', '--no-home', '--pid', '-B', bind_paths]
    samtools_base = common + [str(samtools_images[0])]
    multiqc_base = common + [str(multiqc_images[0])]
    standard_library_probe = 'import csv,gzip,math,subprocess,xml.etree.ElementTree; print("HOST_PYTHON_STDLIB_OK")'
    subprocess.run(samtools_base + [str(host_python_real), '-c', standard_library_probe], check=True)
    subprocess.run(samtools_base + ['samtools', '--version'], check=True, stdout=subprocess.DEVNULL)
    subprocess.run(samtools_base + [str(host_python_real), str(args.project_dir / 'build_igv_evidence_bundle.py'), '--help'], check=True, stdout=subprocess.DEVNULL)
    subprocess.run(multiqc_base + [str(host_python_real), '-c', standard_library_probe], check=True)
    subprocess.run(multiqc_base + [str(host_python_real), str(args.project_dir / 'expression_go_analysis.py'), '--help'], check=True, stdout=subprocess.DEVNULL)
    subprocess.run(multiqc_base + [str(host_python_real), str(args.project_dir / 'analyze_progression_biology.py'), '--help'], check=True, stdout=subprocess.DEVNULL)
    validate_samples(args.samplesheet, args.sra_dir)
    print('PASS: runtime references, containers, SRA inputs, samplesheet, dependency-free host Python, featureCounts, samtools/Python and MultiQC/GO Python contracts')

if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        print(f'ERROR: {error}', file=sys.stderr)
        raise SystemExit(1)
