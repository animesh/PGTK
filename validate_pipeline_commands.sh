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
    local dir="$1"
    if [[ -d "$dir" && -w "$dir" ]]; then
        local probe="$dir/.pgtk_write_test_$$"
        if : > "$probe" 2>/dev/null; then rm -f "$probe"; pass "Directory writable: $dir"; else fail "Write test failed: $dir"; fi
    else
        fail "Directory missing or not writable: $dir"
    fi
}

require_terms() {
    local label="$1" file="$2"
    shift 2
    local missing=0 term
    for term in "$@"; do
        if ! grep -Fq -- "$term" "$file"; then fail "$label missing: $term"; missing=1; fi
    done
    (( missing == 0 )) && pass "$label"
}

reject_terms() {
    local label="$1" file="$2"
    shift 2
    local found=0 term
    for term in "$@"; do
        if grep -Fq -- "$term" "$file"; then fail "$label contains forbidden term: $term"; found=1; fi
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
    if [[ ! -s "$image" ]]; then fail "Container missing: $image"; return; fi
    if apptainer inspect "$image" >/dev/null 2>&1; then pass "Container valid: $(basename "$image")"; else fail "Container inspection failed: $image"; fi
}

check_executable() {
    local label="$1" image="$2" executable="$3"
    if container_exec "$image" sh -c "command -v '$executable' >/dev/null 2>&1 || test -x '$executable'"; then pass "$label executable"; else fail "$label executable missing: $executable"; fi
}

process_block() {
    local process_name="$1"
    awk -v name="$process_name" '
        $0 ~ "^process[[:space:]]+" name "[[:space:]]*\\{" {inside=1}
        inside {print}
        inside && /^}/ {exit}
    ' "$MAIN_NF"
}

check_process_resources() {
    local name="$1" cpus="$2" memory="$3" queue="$4"
    local block="$TEST_ROOT/${name}.block"
    process_block "$name" > "$block"
    if [[ ! -s "$block" ]]; then fail "Process missing: $name"; return; fi
    require_terms "$name resources" "$block" "cpus $cpus" "memory '$memory'" "queue '$queue'"
}

section "Environment and filesystem"
printf 'Started: %s\nProject: %s\nReport: %s\n' "$(date --iso-8601=seconds)" "$PROJECT_DIR" "$REPORT"
require_file "$MAIN_NF"
require_file "$SLURM_FILE"
require_file "$SAMPLES"
command -v apptainer >/dev/null 2>&1 && pass "Apptainer found: $(command -v apptainer)" || fail "Apptainer missing"
command -v java >/dev/null 2>&1 && pass "Java found: $(command -v java)" || fail "Java missing"
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
python - "$SAMPLES" "$SRA_DIR" <<'PY'
import csv, pathlib, re, sys
samples = pathlib.Path(sys.argv[1])
sra_root = pathlib.Path(sys.argv[2])
with samples.open(newline='') as handle:
    rows = list(csv.DictReader(handle))
required = ['sample', 'srr', 'TK', 'Group', 'baseline']
if not rows:
    raise SystemExit('samples.csv has no data rows')
if list(rows[0].keys()) != required:
    raise SystemExit(f'header must be exactly {required}; found {list(rows[0].keys())}')
seen = set()
baselines = {}
for number, row in enumerate(rows, 2):
    for key in required[:4]:
        if not (row[key] or '').strip():
            raise SystemExit(f'row {number}: empty {key}')
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
        raise SystemExit(f'row {number}: baseline must be true, false, or blank')
    if baseline == 'true':
        baselines[key] = baselines.get(key, 0) + 1
    archive = sra_root / srr / f'{srr}.sra'
    if not archive.is_file() or archive.stat().st_size == 0:
        raise SystemExit(f'missing SRA archive: {archive}')
for key, count in baselines.items():
    if count != 1:
        raise SystemExit(f'{key}: expected exactly one baseline, found {count}')
print(f'validated {len(rows)} samples and local SRA paths')
PY
[[ $? -eq 0 ]] && pass "Samplesheet and SRA paths" || fail "Samplesheet or SRA-path validation"

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
check_executable "pypgatk" "$CONTAINER_DIR/${containers[7]}" pypgatk
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
human_reviewed_isoforms.fasta.gz
homo_sapiens_vep_111_GRCh38.tar.gz
arriba_v2.4.0.tar.gz
)
for name in "${references[@]}"; do require_file "$REFERENCE_DIR/$name"; done
for name in "${references[@]:0:4}"; do gzip -t "$REFERENCE_DIR/$name" && pass "gzip integrity: $name" || fail "gzip integrity: $name"; done
tar -tzf "$REFERENCE_DIR/${references[4]}" >/dev/null 2>&1 && pass "tar integrity: ${references[4]}" || fail "tar integrity: ${references[4]}"
tar -tzf "$REFERENCE_DIR/${references[5]}" >/dev/null 2>&1 && pass "tar integrity: ${references[5]}" || fail "tar integrity: ${references[5]}"
VEP_ARCHIVE_LIST="$TEST_ROOT/vep_cache_archive.list"
ARRIBA_ARCHIVE_LIST="$TEST_ROOT/arriba_archive.list"

if tar -tzf "$REFERENCE_DIR/${references[4]}" > "$VEP_ARCHIVE_LIST"; then
    pass "VEP archive listing"
else
    fail "VEP archive listing"
fi

if tar -tzf "$REFERENCE_DIR/${references[5]}" > "$ARRIBA_ARCHIVE_LIST"; then
    pass "Arriba archive listing"
else
    fail "Arriba archive listing"
fi

if grep -Eq '(^|/)homo_sapiens/111_GRCh38(/|$)' "$VEP_ARCHIVE_LIST"; then
    pass "VEP cache layout"
else
    fail "VEP cache layout"
    grep -E 'homo_sapiens|111_GRCh38' "$VEP_ARCHIVE_LIST" | head -20 || true
fi

if grep -Fq 'blacklist_hg38_GRCh38' "$ARRIBA_ARCHIVE_LIST"; then
    pass "Arriba blacklist asset"
else
    fail "Arriba blacklist asset"
fi

if grep -Fq 'known_fusions_hg38_GRCh38' "$ARRIBA_ARCHIVE_LIST"; then
    pass "Arriba known-fusions asset"
else
    fail "Arriba known-fusions asset"
fi

if grep -Fq 'protein_domains_hg38_GRCh38' "$ARRIBA_ARCHIVE_LIST"; then
    pass "Arriba protein-domains asset"
else
    fail "Arriba protein-domains asset"
fi

section "Workflow syntax and corruption checks"
if grep -nE '&gt;|&lt;|&amp;|-&gt;' "$MAIN_NF"; then fail "main.nf contains HTML entities"; else pass "No HTML entities"; fi
if grep -n $'\r' "$MAIN_NF" "$SLURM_FILE" >/dev/null; then fail "CRLF characters detected"; else pass "Unix line endings"; fi
[[ $(grep -c '^process ' "$MAIN_NF") -eq 27 ]] && pass "Expected 27 process declarations" || fail "Unexpected process count: $(grep -c '^process ' "$MAIN_NF")"
processes=(DOWNLOAD_REFERENCES SRA_TO_FASTQ CAT_FASTQ FASTQC_RAW TRIM_GALORE FASTQC_TRIMMED STAR_INDEX STAR_ALIGN SORT_INDEX_BAM SAMTOOLS_FLAGSTAT REF_INDEX MARK_DUPLICATES SPLIT_N_CIGAR HAPLOTYPE_CALLER GENOTYPE_FILTER BCFTOOLS_STATS VEP_ANNOTATE PYPGATK_FASTA ARRIBA FUSION_FASTA STRINGTIE_ASSEMBLY GFFCOMPARE_NOVEL SPLICE_PROTEIN_FASTA COMBINE_PROTEIN_FASTA PROGRESSION_SUBTRACT PROGRESSION_FASTA MULTIQC)
for name in "${processes[@]}"; do [[ $(grep -c "^process $name " "$MAIN_NF") -eq 1 ]] && pass "Process declared once: $name" || fail "Process declaration count invalid: $name"; done

section "Workflow command contracts"
require_terms "Task-local fasterq scratch" "$MAIN_NF" "--split-files" "--threads \${task.cpus}" "--temp fasterq_tmp" "--outdir ." "rm -rf fasterq_tmp"
require_terms "Reference indexing" "$MAIN_NF" "samtools faidx genome.fa" "samtools dict -o genome.dict genome.fa"
require_terms "STAR two-pass and Arriba BAM" "$MAIN_NF" "--twopassMode Basic" "--chimOutType WithinBAM HardClip" "--outSAMtype BAM Unsorted"
require_terms "GATK bounded heaps" "$MAIN_NF" "-Xmx56g" "-Xmx28g" "--MAX_RECORDS_IN_RAM 1000000"
require_terms "VEP all-transcript mode" "$MAIN_NF" "--cache_version 111" "--dir_cache" "--canonical" "--protein" "--hgvs" "--fork \${task.cpus}"
reject_terms "VEP selected-transcript options absent" "$MAIN_NF" "--pick" "--flag_pick"
require_terms "pypgatk 0.0.24 interface" "$MAIN_NF" "--vcf" "--input_fasta" "--gene_annotations_gtf" "--annotation_field_name CSQ" "--af_field AF" "--include_consequences" "--output_proteindb"
reject_terms "Legacy pypgatk options absent" "$MAIN_NF" "--input-vcf" "--protein-db-fasta" "--consequence-filter" "--output-proteindb"
require_terms "Claude zero-result guards" "$MAIN_NF" "no variant proteins generated" "no progression proteins generated"
require_terms "GFFCompare robust output handling" "$MAIN_NF" 'prefix=${meta.sample}.gffcompare' 'if [[ ! -s \${prefix}.annotated.gtf ]]' 'if [[ -s \${prefix}.stats ]]' 'No non-empty gffcompare statistics report' 'test -e ${meta.sample}.novel.gtf'
require_terms "TransDecoder absolute paths" "$MAIN_NF" "/usr/local/opt/transdecoder/util/gtf_genome_to_cdna_fasta.pl" "/usr/local/opt/transdecoder/util/TransDecoder.LongOrfs" "/usr/local/opt/transdecoder/util/TransDecoder.Predict"
require_terms "Reference routing" "$MAIN_NF" "PYPGATK_FASTA(annotated,refs.gtf,refs.cdna)" "PROGRESSION_FASTA(prog,refs.gtf,refs.cdna)" "COMBINE_PROTEIN_FASTA(combined_inputs,refs.proteome)"
require_terms "Novel class-code configuration" "$MAIN_NF" "params.splice_class_codes = 'j,u'" "-v allowed='\${params.splice_class_codes}'"
require_terms "MultiQC aggregation" "$MAIN_NF" "raw_qc.qc" "trimmed_qc.qc" "trim_result.reports" "star_result.logs" "md_result.metrics" "variant_stats" "MULTIQC(qc_files)"
reject_terms "No network downloads in main.nf" "$MAIN_NF" "curl " "wget " "prefetch "
reject_terms "No obsolete workflow scratch configuration" "$MAIN_NF" "params.shared_tmp_root" "singularity.runOptions"

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
cat > "$FIXTURE/ref.gtf" <<'EOF'
chr1	ref	exon	100	199	.	+	.	gene_id "REF1"; transcript_id "REF1.1";
EOF
cat > "$FIXTURE/query.gtf" <<'EOF'
chr1	query	transcript	100	199	.	+	.	gene_id "Q1"; transcript_id "Q1.1";
chr1	query	exon	100	199	.	+	.	gene_id "Q1"; transcript_id "Q1.1";
chr1	query	transcript	300	399	.	+	.	gene_id "Q2"; transcript_id "Q2.1";
chr1	query	exon	300	399	.	+	.	gene_id "Q2"; transcript_id "Q2.1";
EOF
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
require_terms "Native Apptainer configuration" "$SLURM_FILE" 'NXF_APPTAINER_CACHEDIR="${PROJECT_DIR}/singularity_cache"' 'APPTAINER_TMPDIR="${RUN_TMP}/apptainer"' '-with-apptainer'
require_terms "Java and Nextflow" "$SLURM_FILE" "module purge" "module load Java/21" 'NXF_OPTS="-Xms2g -Xmx12g"'
require_terms "SLURM child-job account" "$SLURM_FILE" '"-process.clusterOptions=--account=nn9036k"'
require_terms "Resume and trace" "$SLURM_FILE" "-resume" '-with-trace "${PROJECT_DIR}/results/pipeline_trace-${SLURM_JOB_ID}.tsv"'
reject_terms "No obsolete or suppressive launcher settings" "$SLURM_FILE" "module --force purge" "-with-singularity" "NXF_SINGULARITY_CACHEDIR" "SINGULARITY_TMPDIR" "MESSAGELEVEL" "-with-report" "-with-timeline" "singularity.runOptions"

section "SLURM scheduler syntax"
if command -v sbatch >/dev/null 2>&1; then
    if sbatch --test-only "$SLURM_FILE" >"$TEST_ROOT/sbatch_test.txt" 2>&1; then pass "sbatch --test-only"; sed -n '1,3p' "$TEST_ROOT/sbatch_test.txt"; else fail "sbatch --test-only"; cat "$TEST_ROOT/sbatch_test.txt"; fi
else
    warn "sbatch unavailable; scheduler validation skipped"
fi

section "Summary"
printf 'PASS: %d\nWARN: %d\nFAIL: %d\nREPORT: %s\n' "$PASS" "$WARN" "$FAIL" "$REPORT"
if (( FAIL > 0 )); then printf 'RESULT: FAILED\n'; exit 1; fi
printf 'RESULT: PASSED\n'
