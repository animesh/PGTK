#!/usr/bin/env python3
import argparse
import csv
import os
import subprocess
import tempfile
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
    parser.add_argument('--pysam-image', required=True, type=Path)
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
    require_file(args.pysam_image, 'Pysam container image')
    if not os.access(args.host_python, os.X_OK):
        raise ValueError(f'host Python is not executable: {args.host_python}')
    if not os.access(args.apptainer, os.X_OK):
        raise ValueError(f'Apptainer is not executable: {args.apptainer}')
    subread_image = args.container_cache / 'subread-2.0.8.img'
    subprocess.run([str(args.apptainer), 'exec', str(subread_image), 'featureCounts', '-v'], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    igv_reports_image = args.container_cache / 'quay.io-biocontainers-igv-reports-1.16.0--pyh7e72e81_0.img'
    require_file(igv_reports_image, 'IGV Reports container image')
    with tempfile.TemporaryDirectory(prefix='.pgtk-igv-reports-smoke-', dir=args.project_dir) as smoke_dir:
        smoke = Path(smoke_dir)
        fasta = smoke / 'smoke.fa'
        fasta.write_text('>chr1\n' + 'A' * 200 + '\n', encoding='utf-8')
        (smoke / 'smoke.fa.fai').write_text('chr1\t200\t6\t200\t201\n', encoding='utf-8')
        bed = smoke / 'smoke.bed'
        bed.write_text('chr1\t9\t20\tSMOKE\n', encoding='utf-8')
        report = smoke / 'smoke.html'
        command_probe = 'command -v create_report || command -v create_reports'
        probe = subprocess.run(
            [str(args.apptainer), 'exec', '--cleanenv', str(igv_reports_image), 'sh', '-c', command_probe],
            check=True, text=True, capture_output=True,
        )
        report_command = probe.stdout.strip().splitlines()[-1]
        subprocess.run(
            [str(args.apptainer), 'exec', '--cleanenv', str(igv_reports_image), report_command, str(bed), '--fasta', str(fasta),
             '--tracks', str(bed), '--flanking', '25', '--title', 'PGTK offline smoke test',
             '--output', str(report)],
            check=True, timeout=180,
        )
        require_file(report, 'self-contained IGV Reports smoke-test HTML')
        html = report.read_text(encoding='utf-8', errors='replace').lower()
        if '<html' not in html or 'smoke' not in html:
            raise ValueError('IGV Reports smoke-test output lacks expected HTML content')
    samtools_images = sorted(args.container_cache.rglob('*samtools-1.21*img'))
    if len(samtools_images) != 1:
        raise ValueError(f'expected exactly one samtools 1.21 image below {args.container_cache}; found {len(samtools_images)}')
    multiqc_images = sorted(args.container_cache.rglob('*multiqc-1.35*img'))
    if len(multiqc_images) != 1:
        raise ValueError(f'expected exactly one MultiQC 1.35 image below {args.container_cache}; found {len(multiqc_images)}')
    host_python_real = args.host_python.resolve()
    native_probe = 'import csv,gzip,math,subprocess,xml.etree.ElementTree; print("CONTAINER_PYTHON_STDLIB_OK")'
    subprocess.run([str(host_python_real), '-c', 'import csv,gzip,math; print("WRAPPER_PYTHON_OK")'], check=True)
    subprocess.run([str(args.apptainer), 'exec', '--cleanenv', str(multiqc_images[0]), 'python3', '-c', native_probe], check=True)
    subprocess.run([str(args.apptainer), 'exec', '--cleanenv', str(multiqc_images[0]), 'python3', str(args.project_dir / 'expression_go_analysis.py'), '--help'], check=True, stdout=subprocess.DEVNULL)
    subprocess.run([str(args.apptainer), 'exec', '--cleanenv', str(multiqc_images[0]), 'python3', str(args.project_dir / 'analyze_progression_biology.py'), '--help'], check=True, stdout=subprocess.DEVNULL)
    pysam_probe = 'import pysam; assert pysam.__version__ == "0.24.0"; assert pysam.__samtools_version__ == "1.23.1"; print("PYSAM_RUNTIME_OK")'
    subprocess.run([str(args.apptainer), 'exec', '--cleanenv', str(args.pysam_image), 'python3', '-c', pysam_probe], check=True)
    pysam_scripts = (
        'build_igv_evidence_bundle.py', 'build_finding_igv_reviews.py',
        'validate_rna_events.py', 'validate_variant_codons.py',
        'validate_variant_read_provenance.py', 'validate_proteogenomic_reads.py',
    )
    for script_name in pysam_scripts:
        subprocess.run([str(args.apptainer), 'exec', '--cleanenv', str(args.pysam_image), 'python3', str(args.project_dir / script_name), '--help'], check=True, stdout=subprocess.DEVNULL)
    validate_samples(args.samplesheet, args.sra_dir)
    print('PASS: runtime references, containers, SRA inputs, samplesheet, wrapper Python, container-native Python, featureCounts, offline IGV Reports and Pysam/HTSlib contracts')

if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        print(f'ERROR: {error}', file=sys.stderr)
        raise SystemExit(1)
