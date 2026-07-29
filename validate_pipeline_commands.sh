#!/bin/bash
set -uo pipefail

PROJECT_DIR="${1:-$(pwd -P)}"
MAIN_NF="${PROJECT_DIR}/main.nf"
SLURM_FILE="${PROJECT_DIR}/scratch.slurm"
SAMPLES="${PROJECT_DIR}/samples.csv"
CONTAINER_DIR="${PROJECT_DIR}/singularity_cache"
REFERENCE_DIR="${PROJECT_DIR}/reference_downloads"
SRA_DIR="${PROJECT_DIR}/sra_cache"
WORK_DIR="/cluster/work/users/ash022/work"
TMP_ROOT="/cluster/work/users/ash022/tmp"
NEXTFLOW="/cluster/home/ash022/scripts/nextflow"
REPORT="${PROJECT_DIR}/pipeline_command_validation_$(date +%Y%m%d_%H%M%S).txt"
TEST_ROOT="${TMPDIR:-/tmp}/pgtk_validation_$$"

PASS=0
WARN=0
FAIL=0
mkdir -p "$TEST_ROOT"
trap 'rm -rf "$TEST_ROOT"' EXIT
exec > >(tee "$REPORT") 2>&1

pass() { PASS=$((PASS + 1)); printf 'PASS  %s\n' "$*"; }
warn() { WARN=$((WARN + 1)); printf 'WARN  %s\n' "$*"; }
fail() { FAIL=$((FAIL + 1)); printf 'FAIL  %s\n' "$*"; }
section() { printf '\n============================================================\n%s\n============================================================\n' "$*"; }

require_file() {
    [[ -s "$1" ]] && pass "File exists: $1" || fail "Missing or empty file: $1"
}

require_dir_writable() {
    local dir="$1" probe
    if [[ -d "$dir" && -w "$dir" ]]; then
        probe="$dir/.pgtk_write_test_$$"
        if : > "$probe" 2>/dev/null; then
            rm -f "$probe"
            pass "Directory writable: $dir"
        else
            fail "Write test failed: $dir"
        fi
    else
        fail "Directory missing or not writable: $dir"
    fi
}

require_terms() {
    local label="$1" file="$2" missing=0 term
    shift 2
    for term in "$@"; do
        if ! grep -Fq -- "$term" "$file"; then
            fail "$label missing: $term"
            missing=1
        fi
    done
    (( missing == 0 )) && pass "$label"
}

reject_terms() {
    local label="$1" file="$2" found=0 term
    shift 2
    for term in "$@"; do
        if grep -Fq -- "$term" "$file"; then
            fail "$label contains forbidden term: $term"
            found=1
        fi
    done
    (( found == 0 )) && pass "$label"
}

container_exec() {
    local image="$1"
    shift
    apptainer exec --cleanenv --bind "$PROJECT_DIR:$PROJECT_DIR" --pwd "$PROJECT_DIR" "$image" "$@"
}

check_container() {
    local image="$1"
    if [[ ! -s "$image" ]]; then
        fail "Container missing: $image"
    elif apptainer inspect "$image" >/dev/null 2>&1; then
        pass "Container valid: $(basename "$image")"
    else
        fail "Container inspection failed: $image"
    fi
}

check_executable() {
    local label="$1" image="$2" executable="$3"
    if container_exec "$image" sh -c "command -v '$executable' >/dev/null 2>&1 || test -x '$executable'"; then
        pass "$label executable"
    else
        fail "$label executable missing: $executable"
    fi
}

process_block() {
    local process_name="$1"
    awk -v name="$process_name" '
        $0 ~ "^process[[:space:]]+" name "[[:space:]]*\\{" { inside=1 }
        inside { print }
        inside && /^}/ { exit }
    ' "$MAIN_NF"
}

check_process_resources() {
    local name="$1" cpus="$2" memory="$3" queue="$4"
    local block="$TEST_ROOT/${name}.block"
    process_block "$name" > "$block"
    if [[ ! -s "$block" ]]; then
        fail "Process missing: $name"
        return
    fi
    require_terms "$name resources" "$block" "cpus $cpus" "memory '$memory'" "queue '$queue'"
}

check_process_terms() {
    local name="$1"
    shift
    local block="$TEST_ROOT/${name}.contract.block"
    process_block "$name" > "$block"
    if [[ ! -s "$block" ]]; then
        fail "Process missing: $name"
        return
    fi
    require_terms "$name contract" "$block" "$@"
}

section "Environment and filesystem"
printf 'Started: %s\nProject: %s\nReport: %s\n' "$(date --iso-8601=seconds)" "$PROJECT_DIR" "$REPORT"
require_file "$MAIN_NF"
require_file "$SLURM_FILE"
require_file "$SAMPLES"
require_file "$PROJECT_DIR/proteogenomics_evidence_report.py"
command -v apptainer >/dev/null 2>&1 && pass "Apptainer found: $(command -v apptainer)" || fail "Apptainer missing"
command -v java >/dev/null 2>&1 && pass "Java found: $(command -v java)" || fail "Java missing"
command -v python >/dev/null 2>&1 && pass "Python found: $(command -v python)" || fail "Python missing"
[[ -x "$NEXTFLOW" ]] && pass "Nextflow executable exists" || fail "Nextflow executable missing: $NEXTFLOW"
require_dir_writable "$WORK_DIR"
require_dir_writable "$TMP_ROOT"
[[ -d "$CONTAINER_DIR" ]] && pass "Container directory exists" || fail "Container directory missing"
[[ -d "$REFERENCE_DIR" ]] && pass "Reference directory exists" || fail "Reference directory missing"
[[ -d "$SRA_DIR" ]] && pass "SRA directory exists" || fail "SRA directory missing"
df -h "$PROJECT_DIR" "$WORK_DIR" "$TMP_ROOT" 2>/dev/null || true

section "Runtime versions"
apptainer --version || fail "Apptainer version command failed"
java -version 2>&1 | sed -n '1,3p'
"$NEXTFLOW" -version || fail "Nextflow version command failed"

section "Samplesheet contract"
if python - "$SAMPLES" "$SRA_DIR" <<'PY'
import csv
import pathlib
import re
import sys

samples = pathlib.Path(sys.argv[1])
sra_root = pathlib.Path(sys.argv[2])
with samples.open(newline='') as handle:
    reader = csv.DictReader(handle)
    fieldnames = reader.fieldnames
    rows = list(reader)
required = ['sample', 'srr', 'TK', 'Group', 'baseline']
if fieldnames != required:
    raise SystemExit(f'header must be exactly {required}; found {fieldnames}')
if not rows:
    raise SystemExit('samples.csv has no data rows')
seen = set()
keys = {}
for number, row in enumerate(rows, 2):
    for field in required[:4]:
        if not (row[field] or '').strip():
            raise SystemExit(f'row {number}: empty {field}')
    sample = row['sample'].strip()
    srr = row['srr'].strip()
    key = row['TK'].strip()
    baseline = (row['baseline'] or '').strip().lower()
    if sample in seen:
        raise SystemExit(f'duplicate sample: {sample}')
    seen.add(sample)
    if not re.fullmatch(r'SRR\d+', srr):
        raise SystemExit(f'row {number}: invalid SRR accession {srr}')
    if baseline not in {'true', 'false', ''}:
        raise SystemExit(f'row {number}: invalid baseline value {baseline!r}')
    data = keys.setdefault(key, {'baseline': 0, 'progression': 0})
    if baseline == 'true':
        data['baseline'] += 1
    elif baseline == 'false':
        data['progression'] += 1
    archive = sra_root / srr / f'{srr}.sra'
    if not archive.is_file() or archive.stat().st_size == 0:
        raise SystemExit(f'missing SRA archive: {archive}')
for key, data in keys.items():
    if data['progression'] and data['baseline'] != 1:
        raise SystemExit(f'{key}: expected one baseline for progression subtraction; found {data["baseline"]}')
print(f'validated {len(rows)} samples and local SRA paths')
PY
then
    pass "Samplesheet and SRA paths"
else
    fail "Samplesheet or SRA-path validation"
fi

section "Container integrity"
containers=(
quay.io-biocontainers-sra-tools-3.2.1--h4304569_0.img
quay.io-biocontainers-trim-galore-0.6.10--hdfd78af_0.img
quay.io-biocontainers-fastqc-0.12.1--hdfd78af_0.img
quay.io-biocontainers-star-2.7.11b--h43eeafb_1.img
quay.io-biocontainers-gatk4-4.6.1.0--py310hdfd78af_0.img
quay.io-biocontainers-samtools-1.21--h96c455f_1.img
quay.io-biocontainers-ensembl-vep-111.0--pl5321h2a3209d_0.img
quay.io-biocontainers-pypgatk-0.0.24--pyhdfd78af_0.img
quay.io-biocontainers-arriba-2.4.0--h0033a41_2.img
quay.io-biocontainers-bcftools-1.21--h8b25389_0.img
quay.io-biocontainers-multiqc-1.35--pyhdfd78af_1.img
stringtie-3.0.3.img
transdecoder-6.0.0.img
pvactools-7.1.1.img
gffcompare-0.12.10.img
)
for name in "${containers[@]}"; do check_container "$CONTAINER_DIR/$name"; done

section "Container executables"
check_executable "fasterq-dump" "$CONTAINER_DIR/${containers[0]}" fasterq-dump
check_executable "Trim Galore" "$CONTAINER_DIR/${containers[1]}" trim_galore
check_executable "FastQC" "$CONTAINER_DIR/${containers[2]}" fastqc
check_executable "STAR" "$CONTAINER_DIR/${containers[3]}" STAR
check_executable "GATK" "$CONTAINER_DIR/${containers[4]}" gatk
check_executable "SAMtools" "$CONTAINER_DIR/${containers[5]}" samtools
check_executable "VEP" "$CONTAINER_DIR/${containers[6]}" vep
check_executable "bgzip in VEP" "$CONTAINER_DIR/${containers[6]}" bgzip
check_executable "tabix in VEP" "$CONTAINER_DIR/${containers[6]}" tabix
check_executable "pypgatk" "$CONTAINER_DIR/${containers[7]}" pypgatk
check_executable "gzip in pypgatk" "$CONTAINER_DIR/${containers[7]}" gzip
check_executable "Arriba" "$CONTAINER_DIR/${containers[8]}" arriba
check_executable "bcftools" "$CONTAINER_DIR/${containers[9]}" bcftools
check_executable "MultiQC" "$CONTAINER_DIR/${containers[10]}" multiqc
check_executable "StringTie" "$CONTAINER_DIR/${containers[11]}" stringtie
check_executable "TransDecoder.LongOrfs" "$CONTAINER_DIR/${containers[12]}" /usr/local/opt/transdecoder/util/TransDecoder.LongOrfs
check_executable "TransDecoder.Predict" "$CONTAINER_DIR/${containers[12]}" /usr/local/opt/transdecoder/util/TransDecoder.Predict
check_executable "TransDecoder GTF converter" "$CONTAINER_DIR/${containers[12]}" /usr/local/opt/transdecoder/util/gtf_genome_to_cdna_fasta.pl
check_executable "pVACfuse" "$CONTAINER_DIR/${containers[13]}" pvacfuse
check_executable "gffcompare" "$CONTAINER_DIR/${containers[14]}" gffcompare

section "Reference archive integrity"
references=(
Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
Homo_sapiens.GRCh38.111.gtf.gz
Homo_sapiens.GRCh38.cdna.all.fa.gz
Homo_sapiens.GRCh38.pep.all.fa.gz
human_reviewed_isoforms.fasta.gz
homo_sapiens_vep_111_GRCh38.tar.gz
arriba_v2.4.0.tar.gz
)
for name in "${references[@]}"; do require_file "$REFERENCE_DIR/$name"; done
for name in "${references[@]:0:5}"; do gzip -t "$REFERENCE_DIR/$name" && pass "gzip integrity: $name" || fail "gzip integrity: $name"; done
tar -tzf "$REFERENCE_DIR/${references[5]}" >/dev/null 2>&1 && pass "tar integrity: ${references[5]}" || fail "tar integrity: ${references[5]}"
tar -tzf "$REFERENCE_DIR/${references[6]}" >/dev/null 2>&1 && pass "tar integrity: ${references[6]}" || fail "tar integrity: ${references[6]}"
VEP_ARCHIVE_LIST="$TEST_ROOT/vep_cache_archive.list"
ARRIBA_ARCHIVE_LIST="$TEST_ROOT/arriba_archive.list"
if tar -tzf "$REFERENCE_DIR/${references[5]}" > "$VEP_ARCHIVE_LIST"; then pass "VEP archive listing"; else fail "VEP archive listing"; fi
if tar -tzf "$REFERENCE_DIR/${references[6]}" > "$ARRIBA_ARCHIVE_LIST"; then pass "Arriba archive listing"; else fail "Arriba archive listing"; fi
grep -Eq '(^|/)homo_sapiens/111_GRCh38(/|$)' "$VEP_ARCHIVE_LIST" && pass "VEP cache layout" || fail "VEP cache layout"
grep -Fq 'blacklist_hg38_GRCh38' "$ARRIBA_ARCHIVE_LIST" && pass "Arriba blacklist asset" || fail "Arriba blacklist asset"
grep -Fq 'known_fusions_hg38_GRCh38' "$ARRIBA_ARCHIVE_LIST" && pass "Arriba known-fusions asset" || fail "Arriba known-fusions asset"
grep -Fq 'protein_domains_hg38_GRCh38' "$ARRIBA_ARCHIVE_LIST" && pass "Arriba protein-domains asset" || fail "Arriba protein-domains asset"

section "Workflow syntax and corruption checks"
HTML_ENTITY_REPORT="$TEST_ROOT/html_entities.txt"
if python - "$MAIN_NF" > "$HTML_ENTITY_REPORT" <<'PY_HTML'
from pathlib import Path
import sys

text = Path(sys.argv[1]).read_text()
amp = chr(38)
entities = (
    amp + 'gt;',
    amp + 'lt;',
    amp + 'amp;',
    '-' + amp + 'gt;',
)

found = []
for number, line in enumerate(text.splitlines(), 1):
    hits = [entity for entity in entities if entity in line]
    if hits:
        found.append((number, hits, line))

if found:
    for number, hits, line in found:
        print(f"{number}: {','.join(hits)}: {line}")
    raise SystemExit(1)
PY_HTML
then
    pass "No literal HTML entities"
else
    cat "$HTML_ENTITY_REPORT"
    fail "main.nf contains literal HTML entities"
fi
if grep -n $'\r' "$MAIN_NF" "$SLURM_FILE" >/dev/null; then fail "CRLF characters detected"; else pass "Unix line endings"; fi
processes=(DOWNLOAD_REFERENCES SRA_TO_FASTQ CAT_FASTQ FASTQC_RAW TRIM_GALORE FASTQC_TRIMMED STAR_INDEX STAR_ALIGN SORT_INDEX_BAM SAMTOOLS_FLAGSTAT REF_INDEX MARK_DUPLICATES SPLIT_N_CIGAR HAPLOTYPE_CALLER GENOTYPE_FILTER BCFTOOLS_STATS VEP_ANNOTATE PYPGATK_FASTA ARRIBA FUSION_FASTA STRINGTIE_ASSEMBLY GFFCOMPARE_NOVEL SPLICE_PROTEIN_FASTA COMBINE_PROTEIN_FASTA PROGRESSION_SUBTRACT PROGRESSION_FASTA MULTIQC VALIDATE_MAXQUANT_INPUTS MAP_MAXQUANT_PEPTIDES ANNOTATE_MAXQUANT_VARIANTS ANALYZE_MAXQUANT_JUNCTIONS VALIDATE_MAXQUANT_SPLICE_JUNCTIONS BUILD_PROTEOGENOMICS_EVIDENCE_REPORT)
[[ $(grep -c '^process ' "$MAIN_NF") -eq ${#processes[@]} ]] && pass "Expected ${#processes[@]} process declarations" || fail "Unexpected process count: $(grep -c '^process ' "$MAIN_NF")"
for name in "${processes[@]}"; do [[ $(grep -c "^process $name " "$MAIN_NF") -eq 1 ]] && pass "Process declared once: $name" || fail "Process declaration count invalid: $name"; done

section "Workflow command contracts"
require_terms "Task-local fasterq scratch" "$MAIN_NF" "--split-files" '--threads ${task.cpus}' "--temp fasterq_tmp" "--outdir ." "rm -rf fasterq_tmp"
require_terms "Reference indexing" "$MAIN_NF" "samtools faidx genome.fa" "samtools dict -o genome.dict genome.fa"
require_terms "STAR two-pass and Arriba BAM" "$MAIN_NF" "--twopassMode Basic" "--chimOutType WithinBAM HardClip" "--outSAMtype BAM Unsorted"
require_terms "GATK bounded heaps" "$MAIN_NF" "-Xmx56g" "-Xmx28g" "--MAX_RECORDS_IN_RAM 1000000"
require_terms "VEP all-transcript mode" "$MAIN_NF" "--cache_version 111" "--dir_cache" "--canonical" "--protein" "--hgvs" '--fork ${task.cpus}'
reject_terms "VEP selected-transcript options absent" "$MAIN_NF" "--pick" "--flag_pick"
require_terms "pypgatk 0.0.24 interface" "$MAIN_NF" "--vcf" "--input_fasta" "--gene_annotations_gtf" "--annotation_field_name CSQ" "--af_field AF" "--include_consequences" "--output_proteindb"
reject_terms "Legacy pypgatk options absent" "$MAIN_NF" "--input-vcf" "--protein-db-fasta" "--consequence-filter" "--output-proteindb"
require_terms "Zero-result protein FASTA guards" "$MAIN_NF" "no variant proteins generated" "no progression proteins generated"
require_terms "GFFCompare robust output handling" "$MAIN_NF" 'prefix=${meta.sample}.gffcompare' 'if [[ ! -s \${prefix}.annotated.gtf ]]' 'if [[ -s \${prefix}.stats ]]' 'No non-empty gffcompare statistics report' 'test -e ${meta.sample}.novel.gtf'
require_terms "TransDecoder absolute paths" "$MAIN_NF" "/usr/local/opt/transdecoder/util/gtf_genome_to_cdna_fasta.pl" "/usr/local/opt/transdecoder/util/TransDecoder.LongOrfs" "/usr/local/opt/transdecoder/util/TransDecoder.Predict"
require_terms "Reference routing" "$MAIN_NF" "PYPGATK_FASTA(annotated,refs.gtf,refs.cdna)" "PROGRESSION_FASTA(prog,refs.gtf,refs.cdna)" "COMBINE_PROTEIN_FASTA(combined_inputs,refs.proteome)"
require_terms "Novel class-code configuration" "$MAIN_NF" "params.splice_class_codes = 'j,u'" "-v allowed='\${params.splice_class_codes}'"
require_terms "MultiQC aggregation" "$MAIN_NF" "raw_qc.qc" "trimmed_qc.qc" "trim_result.reports" "star_result.logs" "md_result.metrics" "variant_stats" "MULTIQC(qc_files)"
check_process_terms MULTIQC "path 'multiqc_report.html'" "path 'multiqc_report_data'" "--filename multiqc_report.html" "--data-dir"
reject_terms "Obsolete MultiQC output absent" "$MAIN_NF" "path 'multiqc_data'"
reject_terms "No network downloads in main.nf" "$MAIN_NF" "curl " "wget " "prefetch "
reject_terms "No obsolete workflow scratch configuration" "$MAIN_NF" "params.shared_tmp_root" "singularity.runOptions"

require_terms "Optional MaxQuant branch" "$MAIN_NF" \
    "params.run_proteogenomic_validation = false" \
    "params.maxquant_txt = null" \
    "params.maxquant_mqpar = null" \
    "VALIDATE_MAXQUANT_INPUTS" \
    "MAP_MAXQUANT_PEPTIDES" \
    "ANNOTATE_MAXQUANT_VARIANTS" \
    "ANALYZE_MAXQUANT_JUNCTIONS" \
    "VALIDATE_MAXQUANT_SPLICE_JUNCTIONS" \
    "BUILD_PROTEOGENOMICS_EVIDENCE_REPORT" \
    "proteogenomics_evidence.report.md"
require_terms "MaxQuant raw-file evidence inputs" "$MAIN_NF" \
    'mq_evidence = file("${params.maxquant_txt}/evidence.txt"' \
    'mq_msms = file("${params.maxquant_txt}/msms.txt"' \
    'mq_protein_groups = file("${params.maxquant_txt}/proteinGroups.txt"'
require_terms "Ensembl peptide reference" "$PROJECT_DIR/download_assets.sh" \
    "Homo_sapiens.GRCh38.pep.all.fa.gz" \
    "Ensembl release 111 protein sequences"
section "Python post-processing scripts"
for script in map_peptides_to_fasta.py annotate_variant_peptides.py analyze_chimeric_splice_peptides.py validate_splice_junction_peptides.py proteogenomics_evidence_report.py; do
    require_file "$PROJECT_DIR/$script"
    if python -m py_compile "$PROJECT_DIR/$script"; then
        pass "Python syntax: $script"
    else
        fail "Python syntax: $script"
    fi
done
section "pypgatk VCF input contracts"
check_process_terms PYPGATK_FASTA \
    'gzip -t ${vcf}' \
    'gzip -dc ${vcf} > ${meta.sample}.pypgatk.vcf' \
    'test -s ${meta.sample}.pypgatk.vcf' \
    "grep -q '^##fileformat=VCF' \${meta.sample}.pypgatk.vcf" \
    "if ! grep -qv '^#' \${meta.sample}.pypgatk.vcf; then" \
    ': > ${meta.sample}.variant_proteins.fasta' \
    '--vcf ${meta.sample}.pypgatk.vcf'

check_process_terms PROGRESSION_FASTA \
    'gzip -t ${vcf}' \
    'gzip -dc ${vcf} > ${meta.sample}.progression.pypgatk.vcf' \
    'test -s ${meta.sample}.progression.pypgatk.vcf' \
    "grep -q '^##fileformat=VCF' \${meta.sample}.progression.pypgatk.vcf" \
    "if ! grep -qv '^#' \${meta.sample}.progression.pypgatk.vcf; then" \
    ': > ${meta.sample}.progression_proteins.fasta' \
    '--vcf ${meta.sample}.progression.pypgatk.vcf'

PYPGATK_BLOCK="$TEST_ROOT/PYPGATK_FASTA.disk.block"
PROGRESSION_BLOCK="$TEST_ROOT/PROGRESSION_FASTA.disk.block"
process_block PYPGATK_FASTA > "$PYPGATK_BLOCK"
process_block PROGRESSION_FASTA > "$PROGRESSION_BLOCK"
require_terms "PYPGATK_FASTA disk allocation" "$PYPGATK_BLOCK" "disk '80 GB'"
require_terms "PROGRESSION_FASTA disk allocation" "$PROGRESSION_BLOCK" "disk '80 GB'"

section "GATK task-local temporary storage"
for process_name in MARK_DUPLICATES SPLIT_N_CIGAR HAPLOTYPE_CALLER GENOTYPE_FILTER; do
    check_process_terms "$process_name" \
        "mkdir -p gatk_tmp" \
        "trap 'rm -rf gatk_tmp' EXIT" \
        '-Djava.io.tmpdir=\${PWD}/gatk_tmp'
done

section "SPLIT_N_CIGAR output contract"
check_process_terms SPLIT_N_CIGAR \
    'if [[ -s ${meta.sample}.split.bai && ! -e ${meta.sample}.split.bam.bai ]]; then' \
    'mv ${meta.sample}.split.bai ${meta.sample}.split.bam.bai' \
    'test -s ${meta.sample}.split.bam' \
    'test -s ${meta.sample}.split.bam.bai'

section "Resource contracts"
check_process_resources DOWNLOAD_REFERENCES 8 "32 GB" normal
check_process_resources SRA_TO_FASTQ 16 "32 GB" normal
check_process_resources STAR_INDEX 20 "64 GB" normal
check_process_resources STAR_ALIGN 32 "256 GB" bigmem
check_process_resources SORT_INDEX_BAM 20 "64 GB" normal
check_process_resources MARK_DUPLICATES 20 "64 GB" normal
check_process_resources SPLIT_N_CIGAR 20 "64 GB" normal
check_process_resources HAPLOTYPE_CALLER 20 "64 GB" normal
check_process_resources VEP_ANNOTATE 20 "64 GB" normal
check_process_resources STRINGTIE_ASSEMBLY 20 "64 GB" normal
check_process_resources SPLICE_PROTEIN_FASTA 20 "64 GB" normal

section "GFFCompare fixture smoke test"
FIXTURE="$TEST_ROOT/gffcompare_fixture"
mkdir -p "$FIXTURE"
printf 'chr1\tref\texon\t100\t199\t.\t+\t.\tgene_id "REF1"; transcript_id "REF1.1";\n' > "$FIXTURE/ref.gtf"
printf '%s\n' \
$'chr1\tquery\ttranscript\t100\t199\t.\t+\t.\tgene_id "Q1"; transcript_id "Q1.1";' \
$'chr1\tquery\texon\t100\t199\t.\t+\t.\tgene_id "Q1"; transcript_id "Q1.1";' \
$'chr1\tquery\ttranscript\t300\t399\t.\t+\t.\tgene_id "Q2"; transcript_id "Q2.1";' \
$'chr1\tquery\texon\t300\t399\t.\t+\t.\tgene_id "Q2"; transcript_id "Q2.1";' > "$FIXTURE/query.gtf"
if apptainer exec --cleanenv --bind "$FIXTURE:$FIXTURE" --pwd "$FIXTURE" "$CONTAINER_DIR/gffcompare-0.12.10.img" gffcompare -r ref.gtf -o fixture query.gtf >/dev/null 2>"$FIXTURE/stderr"; then
    [[ -s "$FIXTURE/fixture.annotated.gtf" ]] && pass "GFFCompare fixture annotated GTF" || fail "GFFCompare fixture annotated GTF missing"
    [[ -e "$FIXTURE/fixture.stats" ]] && pass "GFFCompare fixture stats path" || fail "GFFCompare fixture stats path missing"
    grep -q 'class_code' "$FIXTURE/fixture.annotated.gtf" && pass "GFFCompare fixture class codes" || fail "GFFCompare fixture class codes missing"
else
    fail "GFFCompare fixture execution"
    sed -n '1,80p' "$FIXTURE/stderr"
fi

section "Launcher contract"
bash -n "$SLURM_FILE" && pass "scratch.slurm Bash syntax" || fail "scratch.slurm Bash syntax"
require_terms "Launcher resources" "$SLURM_FILE" "#SBATCH --account=nn9036k" "#SBATCH --cpus-per-task=4" "#SBATCH --mem=16G" "#SBATCH --time=168:00:00"
require_terms "Storage placement" "$SLURM_FILE" 'WORK_DIR="/cluster/work/users/ash022/work"' 'TMP_ROOT="/cluster/work/users/ash022/tmp"' '-work-dir "${WORK_DIR}"'
require_terms "Native Apptainer configuration" "$SLURM_FILE" 'NXF_APPTAINER_CACHEDIR="${PROJECT_DIR}/singularity_cache"' 'APPTAINER_TMPDIR="${RUN_TMP}/apptainer"' "-with-apptainer"
require_terms "Java and Nextflow" "$SLURM_FILE" "module purge" "module load Java/21" 'NXF_OPTS="-Xms2g -Xmx12g"'
require_terms "SLURM child-job account" "$SLURM_FILE" '"-process.clusterOptions=--account=nn9036k"'
require_terms "Resume, trace, and optional argument forwarding" "$SLURM_FILE" "-resume" '-with-trace "${PROJECT_DIR}/results/pipeline_trace-${SLURM_JOB_ID}.tsv"' '"$@"'
reject_terms "No obsolete or suppressive launcher settings" "$SLURM_FILE" "module --force purge" "-with-singularity" "NXF_SINGULARITY_CACHEDIR" "SINGULARITY_TMPDIR" "MESSAGELEVEL" "-with-report" "-with-timeline" "singularity.runOptions"

section "SLURM scheduler syntax"
if command -v sbatch >/dev/null 2>&1; then
    if sbatch --test-only "$SLURM_FILE" >"$TEST_ROOT/sbatch_test.txt" 2>&1; then
        pass "sbatch --test-only"
        sed -n '1,3p' "$TEST_ROOT/sbatch_test.txt"
    else
        fail "sbatch --test-only"
        cat "$TEST_ROOT/sbatch_test.txt"
    fi
else
    warn "sbatch unavailable; scheduler validation skipped"
fi

section "Summary"
printf 'PASS: %d\nWARN: %d\nFAIL: %d\nREPORT: %s\n' "$PASS" "$WARN" "$FAIL" "$REPORT"
if (( FAIL > 0 )); then
    printf 'RESULT: FAILED\n'
    exit 1
fi
printf 'RESULT: PASSED\n'
