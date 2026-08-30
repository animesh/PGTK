#!/usr/bin/env python3
import argparse, csv, gzip, os, subprocess, tarfile, xml.etree.ElementTree as ET
from pathlib import Path

REFERENCES = {
    'Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz': 'gzip',
    'Homo_sapiens.GRCh38.111.gtf.gz': 'gzip',
    'Homo_sapiens.GRCh38.cdna.all.fa.gz': 'gzip',
    'Homo_sapiens.GRCh38.pep.all.fa.gz': 'gzip',
    'homo_sapiens_vep_111_GRCh38.tar.gz': 'tar',
    'arriba_v2.4.0.tar.gz': 'tar',
    'go-basic.obo': 'obo',
    'goa_human.gaf.gz': 'gaf',
}
IMAGES = {
    'quay.io-biocontainers-sra-tools-3.2.1--h4304569_0.img': [('fasterq-dump','--version'),('vdb-validate','--version')],
    'quay.io-biocontainers-trim-galore-0.6.10--hdfd78af_0.img': [('trim_galore','--version')],
    'quay.io-biocontainers-fastqc-0.12.1--hdfd78af_0.img': [('fastqc','--version')],
    'quay.io-biocontainers-star-2.7.11b--h43eeafb_1.img': [('STAR','--version')],
    'quay.io-biocontainers-gatk4-4.6.1.0--py310hdfd78af_0.img': [('gatk','--version')],
    'quay.io-biocontainers-ensembl-vep-111.0--pl5321h2a3209d_0.img': [('vep','--help')],
    'quay.io-biocontainers-pypgatk-0.0.24--pyhdfd78af_0.img': [('pypgatk','--help')],
    'quay.io-biocontainers-arriba-2.4.0--h0033a41_2.img': [('arriba','-h')],
    'quay.io-biocontainers-bcftools-1.21--h8b25389_0.img': [('bcftools','--version')],
    'quay.io-biocontainers-samtools-1.21--h96c455f_1.img': [('samtools','--version')],
    'quay.io-biocontainers-multiqc-1.35--pyhdfd78af_1.img': [('multiqc','--version'),('python3','--version')],
    'subread-2.0.8.img': [('featureCounts','-v')],
    'stringtie-3.0.3.img': [('stringtie','--version')],
    'transdecoder-6.0.0.img': [('sh','-c','test -x /usr/local/opt/transdecoder/util/gtf_genome_to_cdna_fasta.pl && test -x /usr/local/opt/transdecoder/util/TransDecoder.LongOrfs && test -x /usr/local/opt/transdecoder/util/TransDecoder.Predict')],
    'gffcompare-0.12.10.img': [('gffcompare','--version')],
    'pvactools-7.1.1.img': [('pvacfuse','--help')],
    'quay.io-biocontainers-pysam-0.24.0--py312hf5ad864_1.img': [('python3','--version')],
    'quay.io-biocontainers-igv-reports-1.16.0--pyh7e72e81_0.img': [('sh','-c','command -v create_report || command -v create_reports')],
}
PYSAM_SCRIPTS = ('variant_read_evidence.py','build_igv_evidence_bundle.py','build_finding_igv_reviews.py','validate_rna_events.py','validate_variant_codons.py','validate_variant_read_provenance.py','validate_proteogenomic_reads.py')
MULTIQC_SCRIPTS = ('report_legend.py','validate_haplotype_shards.py','summarize_variant_stages.py','compare_external_vcf.py','expression_go_analysis.py','prepare_go_annotations.py','analyze_progression_biology.py','compare_progression_pair.py','merge_progression_biology.py','analyze_variant_landscape.py','build_finding_explorer.py','build_compact_multiqc_content.py','build_expression_multiqc_content.py','build_pgtk_multiqc_content.py','build_complete_report.py','build_comparative_advantage_report.py','map_peptides_to_fasta.py','annotate_variant_peptides.py','analyze_chimeric_splice_peptides.py','validate_splice_junction_peptides.py','proteogenomics_evidence_report.py','merge_variant_validation.py','analyze_codon_mismatches.py','build_integrated_variant_evidence.py')

def require_file(path, label):
    if not path.is_file() or path.stat().st_size == 0: raise RuntimeError(f'{label} missing or empty: {path}')
def require_dir(path, label, writable=False):
    if not path.is_dir(): raise RuntimeError(f'{label} missing: {path}')
    if not os.access(path, os.R_OK): raise RuntimeError(f'{label} unreadable: {path}')
    if writable and not os.access(path, os.W_OK): raise RuntimeError(f'{label} not writable: {path}')
def run(command, **kwargs):
    return subprocess.run([str(value) for value in command], check=True, **kwargs)
def describe_size(path):
    size = path.stat().st_size
    for unit in ('B','KiB','MiB','GiB','TiB'):
        if size < 1024 or unit == 'TiB': return f'{size:.1f} {unit}'
        size /= 1024
def exec_in(apptainer, image, *command, bind_paths=(), timeout=180, label=None):
    description = label or f'{Path(image).name}: {" ".join(map(str, command))}'
    print(f'  RUN   {description}', flush=True)
    try:
        bind_args=[item for path in bind_paths for item in ('--bind',f'{path}:{path}')]; run([apptainer,'exec','--cleanenv','--no-home',*bind_args,image,*command], stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True, timeout=timeout)
    except subprocess.TimeoutExpired as error:
        raise RuntimeError(f'{description} timed out after {timeout} seconds') from error
    except subprocess.CalledProcessError as error:
        diagnostic = (error.stderr or '').strip() or 'no stderr output'
        raise RuntimeError(f'{description} failed with exit {error.returncode}: {diagnostic}') from error
    print(f'  PASS  {description}', flush=True)
def strict_bool(value, name):
    normalized = str(value).strip().lower()
    if normalized in {'true','1','yes','y','on'}: return True
    if normalized in {'false','0','no','n','off'}: return False
    raise RuntimeError(f'{name} must be true or false: {value}')
def parse_pipeline_options(items):
    values = {}; index = 0
    while index < len(items):
        token = items[index]
        if token.startswith('--'):
            key = token[2:].replace('-','_')
            if index + 1 < len(items) and not items[index + 1].startswith('--'):
                values[key] = items[index + 1]; index += 2
            else:
                values[key] = 'true'; index += 1
        else: index += 1
    return values

def main():
    parser = argparse.ArgumentParser()
    for name in ('project-dir','reference-downloads','container-cache','sra-dir','samplesheet','work-dir','tmp-root','results-dir','host-python','apptainer','pysam-image','nextflow'):
        parser.add_argument('--' + name, required=True, type=Path)
    parser.add_argument('--bind-path', action='append', default=[], type=Path)
    parser.add_argument('pipeline_args', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    print('PGTK COMPLETE RUNTIME PREFLIGHT', flush=True)
    project = args.project_dir.resolve(); references = args.reference_downloads.resolve(); cache = args.container_cache.resolve(); apptainer = args.apptainer.resolve()
    print(f'Project: {project}', flush=True)
    print(f'References: {references}', flush=True)
    print(f'Containers: {cache}', flush=True)
    print(f'SRA directory: {args.sra_dir.resolve()}', flush=True)
    bind_paths=[]
    for value in args.bind_path:
        resolved=value.expanduser().resolve()
        if not resolved.exists(): raise RuntimeError(f'bind path does not exist: {resolved}')
        if any(char.isspace() for char in str(resolved)): raise RuntimeError(f'bind path contains whitespace: {resolved}')
        if resolved not in bind_paths: bind_paths.append(resolved)
    required_targets=(project,references,cache,args.sra_dir.resolve(),args.work_dir.resolve(),args.tmp_root.resolve(),args.results_dir.resolve(),args.samplesheet.resolve().parent)
    for target in required_targets:
        if not any(target == parent or parent in target.parents for parent in bind_paths):
            raise RuntimeError(f'bind contract does not cover required path: {target}')
    print('Container bind contract:', flush=True)
    for path in bind_paths: print(f'  BIND  {path}', flush=True)
    def require_bound(path, label):
        target=Path(path).expanduser().resolve()
        if not any(target == parent or parent in target.parents for parent in bind_paths):
            raise RuntimeError(f'bind contract does not cover {label}: {target}')
    print('\n[1/8] Directories and required source files', flush=True)
    for path,label,writable in ((project,'project',False),(references,'references',False),(cache,'container cache',False),(args.sra_dir,'SRA directory',False),(args.work_dir,'work directory',True),(args.tmp_root,'temporary root',True),(args.results_dir,'results directory',True)):
        require_dir(Path(path), label, writable)
        print(f'  PASS  {label}: {Path(path).resolve()}', flush=True)
    required_names=[line.strip() for line in (project/'pipeline_required_files.txt').read_text().splitlines() if line.strip()]
    for name in required_names: require_file(project/name, 'required source')
    print(f'  PASS  {len(required_names)} required source files', flush=True)
    print('\n[2/8] Host executables', flush=True)
    for executable,label in ((args.host_python,'host Python'),(apptainer,'Apptainer'),(args.nextflow,'Nextflow')):
        require_file(executable,label)
        if not os.access(executable,os.X_OK): raise RuntimeError(f'{label} is not executable: {executable}')
        print(f'  PASS  {label}: {executable}', flush=True)
    print(f'\n[3/8] Reference assets ({len(REFERENCES)})', flush=True)
    for name,kind in REFERENCES.items():
        path=references/name; require_file(path,'reference asset')
        print(f'  CHECK {name} ({describe_size(path)})', flush=True)
        if kind=='gzip':
            with gzip.open(path,'rb') as handle: handle.read(65536)
        elif kind=='tar':
            with tarfile.open(path,'r:gz') as handle:
                if not handle.getmembers(): raise RuntimeError(f'empty archive: {path}')
        elif kind=='obo':
            if '[Term]' not in path.read_text(errors='replace')[:200000]: raise RuntimeError(f'invalid OBO: {path}')
        elif kind=='gaf':
            with gzip.open(path,'rt',errors='replace') as handle:
                if not any(line.startswith('!gaf-version:') for _,line in zip(range(100),handle)): raise RuntimeError(f'invalid GAF: {path}')
        print(f'  PASS  {name}', flush=True)
    print('\n[4/8] Exact container path visibility', flush=True)
    import shlex
    visible=(project/'variant_read_evidence.py',references,cache,args.sra_dir.resolve(),args.work_dir.resolve(),args.tmp_root.resolve(),args.results_dir.resolve())
    command='set -eu; ' + '; '.join(f'test -e {shlex.quote(str(path))}' for path in visible)
    exec_in(apptainer,args.pysam_image.resolve(),'sh','-c',command,bind_paths=bind_paths,label='Pysam container: configured paths visible')
    print(f'\n[5/8] Exact container images and tool smoke tests ({len(IMAGES)})', flush=True)
    for name,commands in IMAGES.items():
        image=cache/name; require_file(image,'container image')
        print(f'  CHECK container {name} ({describe_size(image)})', flush=True)
        run([apptainer,'inspect',image],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,timeout=120)
        for command in commands: exec_in(apptainer,image,*command,bind_paths=bind_paths)
    print('\n[6/8] Pipeline Python scripts in target containers', flush=True)
    pysam_image=args.pysam_image.resolve(); require_file(pysam_image,'Pysam image')
    exec_in(apptainer,pysam_image,'python3','-c','import pysam; assert pysam.__version__=="0.24.0"; assert pysam.__samtools_version__=="1.23.1"',bind_paths=bind_paths)
    for script in PYSAM_SCRIPTS: exec_in(apptainer,pysam_image,'python3',project/script,'--help',bind_paths=bind_paths)
    multiqc=cache/'quay.io-biocontainers-multiqc-1.35--pyhdfd78af_1.img'
    for script in MULTIQC_SCRIPTS: exec_in(apptainer,multiqc,'python3',project/script,'--help',bind_paths=bind_paths)
    print('\n[7/8] Samplesheet and SRA integrity', flush=True)
    require_file(args.samplesheet,'samplesheet')
    with args.samplesheet.open(newline='',encoding='utf-8') as handle:
        reader=csv.DictReader(handle); fields=reader.fieldnames or []; rows=list(reader)
    if not {'sample','srr'} <= set(fields) or not rows: raise RuntimeError('samplesheet requires sample,srr and at least one row')
    samples=set(); srrs=set(); sra_image=cache/'quay.io-biocontainers-sra-tools-3.2.1--h4304569_0.img'
    print(f'  CHECK {len(rows)} samples from {args.samplesheet}', flush=True)
    for line,row in enumerate(rows,2):
        sample=(row.get('sample') or '').strip(); srr=(row.get('srr') or '').strip()
        if not sample or not srr: raise RuntimeError(f'empty sample or srr at line {line}')
        if sample in samples or srr in srrs: raise RuntimeError(f'duplicate sample or srr at line {line}')
        samples.add(sample); srrs.add(srr)
        if (row.get('baseline') or 'false').strip().lower() not in {'true','false'}: raise RuntimeError(f'invalid baseline at line {line}')
        archive=args.sra_dir/srr/f'{srr}.sra'; require_file(archive,f'SRA for {sample}'); print(f'  CHECK {sample}: {srr} ({describe_size(archive)})', flush=True); exec_in(apptainer,sra_image,'vdb-validate',archive,bind_paths=bind_paths,timeout=1800,label=f'vdb-validate {sample} {srr}')
    print('\n[8/8] Optional branches and Nextflow inspection', flush=True)
    options=parse_pipeline_options(args.pipeline_args[1:] if args.pipeline_args[:1]==['--'] else args.pipeline_args)
    external_enabled=strict_bool(options.get('run_external_vcf_comparison','false'),'run_external_vcf_comparison')
    proteogenomic_enabled=strict_bool(options.get('run_proteogenomic_validation','false'),'run_proteogenomic_validation')
    print(f'  External VCF comparison enabled: {external_enabled}', flush=True)
    print(f'  Proteogenomic validation enabled: {proteogenomic_enabled}', flush=True)
    if external_enabled:
        root=Path(options.get('external_vcf_dir','')); suffix=options.get('external_vcf_suffix',''); require_bound(root,'external VCF directory')
        require_dir(root,'external VCF directory')
        if not suffix: raise RuntimeError('--external_vcf_suffix is required when external comparison is enabled')
        for row in rows:
            matches=list(root.rglob(row['srr'].strip()+suffix))
            if len(matches)!=1: raise RuntimeError(f"expected one external VCF for {row['srr']}; found {len(matches)}")
            require_file(matches[0],'external VCF'); require_file(Path(str(matches[0])+'.tbi'),'external VCF index')
    if proteogenomic_enabled:
        ensembl=Path(options.get('ensembl_pep', references/'Homo_sapiens.GRCh38.pep.all.fa.gz')); require_file(ensembl,'Ensembl peptide FASTA')
        maxquant=Path(options.get('maxquant_txt','')); require_dir(maxquant,'MaxQuant txt directory'); require_bound(maxquant,'MaxQuant txt directory')
        for name in ('peptides.txt','evidence.txt','msms.txt','proteinGroups.txt'): require_file(maxquant/name,'MaxQuant input')
        override=options.get('maxquant_mqpar',''); candidates=(Path(override),) if override else (maxquant/'mqpar.xml',maxquant.parent/'mqpar.xml')
        mqpar=next((path for path in candidates if path.is_file() and path.stat().st_size),None)
        if mqpar is None: raise RuntimeError('mqpar.xml missing')
        require_bound(mqpar.parent,'MaxQuant mqpar directory')
        root=ET.parse(mqpar).getroot(); fasta_values=[(node.text or '').strip() for node in root.findall('./fastaFiles/FastaFileInfo/fastaFilePath') if (node.text or '').strip()]
        if not fasta_values: raise RuntimeError('mqpar.xml contains no FASTA paths')
        contaminant=options.get('maxquant_contaminants','')
        if contaminant: require_file(Path(contaminant),'MaxQuant contaminants FASTA'); require_bound(Path(contaminant).parent,'MaxQuant contaminants directory')
    environment=dict(os.environ,PGTK_ACCOUNT='validation',PGTK_NORMAL_PARTITION='normal',PGTK_BIGMEM_PARTITION='bigmem',PGTK_NORMAL_CPU_THRESHOLD='20',PGTK_NORMAL_MEMORY_THRESHOLD_GB='160',PGTK_MAX_CPUS='32',PGTK_MAX_MEMORY_GB='512',PGTK_QUEUE_SIZE='200',PGTK_SUBMIT_RATE_LIMIT='60/1min')
    print('  RUN   Nextflow inspect', flush=True)
    run([args.nextflow,'inspect',project/'main.nf'],env=environment,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,timeout=300)
    print('  PASS  Nextflow inspect', flush=True)
    print('\nPASS: COMPLETE PRE-SUBMISSION RUNTIME VALIDATION')
    print(f'PASS: {len(IMAGES)} exact containers and required executables')
    print(f'PASS: {len(REFERENCES)} required reference assets')
    print(f'PASS: {len(rows)} samples and vdb-validated SRA archives')
    print('PASS: enabled optional branches and exact inputs')

if __name__=='__main__':
    try: main()
    except Exception as error:
        print(f'ERROR: PRE-SUBMISSION VALIDATION FAILED: {error}',file=__import__('sys').stderr)
        raise SystemExit(1)
