#!/bin/bash
set -uo pipefail

WORKDIR="${1:-$(pwd -P)}"
MAIN_NF="${2:-$WORKDIR/main.nf}"
CONTAINER_DIR="$WORKDIR/singularity_cache"
REFERENCE_DIR="$WORKDIR/reference_downloads"
SLURM_FILE="${3:-$WORKDIR/scratch.slurm}"
SCRATCH_ROOT="/cluster/work/users/ash022/tmp/pgtk"
STAMP=$(date +%Y%m%d_%H%M%S)
REPORT="$WORKDIR/pipeline_command_validation_${STAMP}.txt"
TMP_ROOT="${TMPDIR:-$WORKDIR/tmp}/pgtk_validation_$$"

mkdir -p "$TMP_ROOT"
trap 'rm -rf "$TMP_ROOT"' EXIT
exec > >(tee "$REPORT") 2>&1

PASS=0
WARN=0
FAIL=0

pass(){ PASS=$((PASS + 1)); printf 'PASS  %s\n' "$*"; }
warn(){ WARN=$((WARN + 1)); printf 'WARN  %s\n' "$*"; }
fail(){ FAIL=$((FAIL + 1)); printf 'FAIL  %s\n' "$*"; }
section(){ printf '\n============================================================\n%s\n============================================================\n' "$*"; }

run_image() {
    local image="$1"
    shift
    singularity exec --bind "$WORKDIR:$WORKDIR" --pwd "$WORKDIR" "$image" "$@"
}

capture() {
    local output="$1"
    local image="$2"
    shift 2
    set +e
    run_image "$image" "$@" >"$output" 2>&1
    local status=$?
    set -e
    [[ -s "$output" ]] || return "$status"
}

contains_all() {
    local label="$1"
    local file="$2"
    shift 2
    local missing=0
    local term
    for term in "$@"; do
        if ! grep -Fq -- "$term" "$file"; then
            fail "$label missing '$term'"
            missing=1
        fi
    done
    (( missing == 0 )) && pass "$label confirmed"
}

main_contains_all() {
    local label="$1"
    shift
    contains_all "$label" "$MAIN_NF" "$@"
}

check_image() {
    local image="$1"
    if [[ -s "$image" ]] && singularity inspect "$image" >/dev/null 2>&1; then
        pass "Container valid: $(basename "$image")"
    else
        fail "Container missing or invalid: $image"
    fi
}

SRA="$CONTAINER_DIR/quay.io-biocontainers-sra-tools-3.2.1--h4304569_0.img"
TRIM="$CONTAINER_DIR/quay.io-biocontainers-trim-galore-0.6.10--hdfd78af_0.img"
FASTQC="$CONTAINER_DIR/quay.io-biocontainers-fastqc-0.12.1--hdfd78af_0.img"
STAR="$CONTAINER_DIR/quay.io-biocontainers-star-2.7.11b--h43eeafb_1.img"
SAMTOOLS="$CONTAINER_DIR/quay.io-biocontainers-samtools-1.21--h96c455f_1.img"
GATK="$CONTAINER_DIR/quay.io-biocontainers-gatk4-4.6.1.0--py310hdfd78af_0.img"
BCFTOOLS="$CONTAINER_DIR/quay.io-biocontainers-bcftools-1.21--h8b25389_0.img"
VEP="$CONTAINER_DIR/quay.io-biocontainers-ensembl-vep-111.0--pl5321h2a3209d_0.img"
PYPGATK="$CONTAINER_DIR/quay.io-biocontainers-pypgatk-0.0.24--pyhdfd78af_0.img"
ARRIBA="$CONTAINER_DIR/quay.io-biocontainers-arriba-2.4.0--h0033a41_2.img"
MULTIQC="$CONTAINER_DIR/quay.io-biocontainers-multiqc-1.35--pyhdfd78af_1.img"
STRINGTIE="$CONTAINER_DIR/stringtie-3.0.3.img"
TRANSDECODER="$CONTAINER_DIR/transdecoder-6.0.0.img"
PVAC="$CONTAINER_DIR/pvactools-7.1.1.img"
GFFCOMPARE="$CONTAINER_DIR/gffcompare-0.12.10.img"

section "Environment"
printf 'Started: %s\nWorkdir: %s\nWorkflow: %s\nReport: %s\n' "$(date --iso-8601=seconds)" "$WORKDIR" "$MAIN_NF" "$REPORT"
command -v singularity >/dev/null 2>&1 && pass "Singularity found" || fail "singularity not found"
[[ -s "$MAIN_NF" ]] && pass "main.nf exists" || fail "main.nf missing"
[[ -s "$WORKDIR/samples.csv" ]] && pass "samples.csv exists" || fail "samples.csv missing"
[[ -s "$SLURM_FILE" ]] && pass "scratch.slurm exists" || fail "scratch.slurm missing"
[[ -d "$SCRATCH_ROOT" && -w "$SCRATCH_ROOT" ]] && pass "Shared Nextflow work root is writable" || fail "Shared Nextflow work root unavailable: $SCRATCH_ROOT"

section "Container integrity"
for image in "$SRA" "$TRIM" "$FASTQC" "$STAR" "$SAMTOOLS" "$GATK" "$BCFTOOLS" "$VEP" "$PYPGATK" "$ARRIBA" "$MULTIQC" "$STRINGTIE" "$TRANSDECODER" "$PVAC" "$GFFCOMPARE"; do
    check_image "$image"
done

section "Installed CLI interfaces"
capture "$TMP_ROOT/fasterq.txt" "$SRA" fasterq-dump --help || true
contains_all "fasterq-dump" "$TMP_ROOT/fasterq.txt" "--split-files" "--threads" "--temp" "--outdir"

capture "$TMP_ROOT/star.txt" "$STAR" STAR --help || true
contains_all "STAR" "$TMP_ROOT/star.txt" "runMode" "genomeDir" "genomeFastaFiles" "sjdbGTFfile" "sjdbOverhang" "readFilesIn" "twopassMode" "outSAMtype" "chimOutType"

capture "$TMP_ROOT/samtools.txt" "$SAMTOOLS" samtools --help || true
contains_all "SAMtools" "$TMP_ROOT/samtools.txt" "faidx" "dict" "sort" "index" "flagstat"

capture "$TMP_ROOT/gatk.txt" "$GATK" gatk --list || true
contains_all "GATK" "$TMP_ROOT/gatk.txt" "MarkDuplicates" "SplitNCigarReads" "HaplotypeCaller" "GenotypeGVCFs" "VariantFiltration" "SelectVariants"

capture "$TMP_ROOT/bcftools.txt" "$BCFTOOLS" bcftools --help || true
contains_all "bcftools" "$TMP_ROOT/bcftools.txt" "stats" "isec" "index"

capture "$TMP_ROOT/vep.txt" "$VEP" vep --help || true
contains_all "VEP basic interface" "$TMP_ROOT/vep.txt" "--input_file" "--output_file" "--offline" "--fork"

capture "$TMP_ROOT/pypgatk.txt" "$PYPGATK" pypgatk vcf-to-proteindb --help || true
contains_all "pypgatk 0.0.24" "$TMP_ROOT/pypgatk.txt" "--input_fasta" "--vcf" "--gene_annotations_gtf" "--output_proteindb" "--annotation_field_name" "--af_field" "--include_consequences"

capture "$TMP_ROOT/arriba.txt" "$ARRIBA" arriba -h || true
contains_all "Arriba" "$TMP_ROOT/arriba.txt" "-x" "-a" "-g" "-b" "-k" "-p" "-o" "-O"

capture "$TMP_ROOT/stringtie.txt" "$STRINGTIE" stringtie --help || true
contains_all "StringTie" "$TMP_ROOT/stringtie.txt" "-G" "-p" "-c" "-f" "-j" "-o"

capture "$TMP_ROOT/td_long.txt" "$TRANSDECODER" /usr/local/opt/transdecoder/util/TransDecoder.LongOrfs || true
contains_all "TransDecoder.LongOrfs" "$TMP_ROOT/td_long.txt" "-t <string>" "-m <int>" "--output_dir"

capture "$TMP_ROOT/td_predict.txt" "$TRANSDECODER" /usr/local/opt/transdecoder/util/TransDecoder.Predict || true
contains_all "TransDecoder.Predict" "$TMP_ROOT/td_predict.txt" "-t <string>" "--output_dir"

if run_image "$TRANSDECODER" test -x /usr/local/opt/transdecoder/util/gtf_genome_to_cdna_fasta.pl; then
    pass "TransDecoder GTF-to-cDNA utility executable"
else
    fail "TransDecoder GTF-to-cDNA utility unavailable"
fi

capture "$TMP_ROOT/pvac.txt" "$PVAC" pvacfuse generate_protein_fasta --help || true
contains_all "pVACfuse" "$TMP_ROOT/pvac.txt" "Arriba fusion.tsv" "flanking_sequence_length" "--downstream-sequence-length"
if grep -Fq -- '--input-type' "$TMP_ROOT/pvac.txt"; then warn "pVACfuse exposes --input-type unexpectedly"; else pass "pVACfuse has no --input-type option"; fi

capture "$TMP_ROOT/gffcompare.txt" "$GFFCOMPARE" gffcompare --help || true
contains_all "gffcompare" "$TMP_ROOT/gffcompare.txt" "-r" "-o"

capture "$TMP_ROOT/multiqc.txt" "$MULTIQC" multiqc --help || true
contains_all "MultiQC" "$TMP_ROOT/multiqc.txt" "--force" "--title" "--filename" "--data-dir" "--data-format"

section "Reference assets"
for asset in \
    "$REFERENCE_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" \
    "$REFERENCE_DIR/Homo_sapiens.GRCh38.111.gtf.gz" \
    "$REFERENCE_DIR/Homo_sapiens.GRCh38.cdna.all.fa.gz" \
    "$REFERENCE_DIR/homo_sapiens_vep_111_GRCh38.tar.gz" \
    "$REFERENCE_DIR/human_reviewed_isoforms.fasta.gz" \
    "$REFERENCE_DIR/arriba_v2.4.0.tar.gz"; do
    [[ -s "$asset" ]] && pass "Reference exists: $(basename "$asset")" || fail "Reference missing: $asset"
done

section "main.nf contract"
main_contains_all "cDNA reference" "params.cdna_url" "emit: cdna" "refs.cdna"
main_contains_all "Task-local fasterq scratch" \
    "--temp fasterq_tmp" \
    "--outdir ." \
    "test -s \${srr}_1.fastq" \
    "test -s \${srr}_2.fastq" \
    "rm -rf fasterq_tmp"

if grep -Fq -- 'params.shared_tmp_root' "$MAIN_NF"; then
    fail "main.nf still contains custom external scratch configuration"
else
    pass "main.nf uses the Nextflow task work directory for fasterq scratch"
fi

main_contains_all "Correct pypgatk interface" "--vcf \${vcf}" "--input_fasta \${cdna}" "--gene_annotations_gtf \${gtf}" "--annotation_field_name CSQ" "--af_field AF" "--include_consequences" "--output_proteindb"

main_contains_all "REF_INDEX SAMtools implementation" "process REF_INDEX" "container 'quay.io/biocontainers/samtools:1.21--h96c455f_1'" "samtools faidx genome.fa" "samtools dict -o genome.dict genome.fa"
if grep -Fq -- 'gatk CreateSequenceDictionary' "$MAIN_NF"; then fail "REF_INDEX still uses GATK CreateSequenceDictionary"; else pass "REF_INDEX does not use GATK CreateSequenceDictionary"; fi
if grep -Eq -- 'cp[[:space:]]+(\$\{genome\}|genome[.]fa)[[:space:]]+genome[.]fa' "$MAIN_NF"; then fail "REF_INDEX contains redundant genome.fa copy"; else pass "REF_INDEX has no redundant genome.fa copy"; fi

if grep -Eq -- '--input-vcf|--protein-db-fasta|--gene-annotations-gtf|--annotation-field-name|--af-field|--consequence-filter|--output-proteindb' "$MAIN_NF"; then
    fail "Unsupported legacy pypgatk options remain"
else
    pass "No unsupported legacy pypgatk options remain"
fi

main_contains_all "Absolute TransDecoder paths" "/usr/local/opt/transdecoder/util/gtf_genome_to_cdna_fasta.pl" "/usr/local/opt/transdecoder/util/TransDecoder.LongOrfs" "/usr/local/opt/transdecoder/util/TransDecoder.Predict"
main_contains_all "gffcompare novelty branch" "process GFFCOMPARE_NOVEL" "class_code" "params.splice_class_codes" "GFFCOMPARE_NOVEL(assembled,refs.gtf)"
main_contains_all "Splice threshold" "params.splice_min_protein_aa = 60"

if grep -Eq -- '(^|[[:space:]])--pick([[:space:]\\]|$)|--flag_pick' "$MAIN_NF"; then fail "VEP incorrectly uses --pick or --flag_pick"; else pass "VEP retains all transcript consequences"; fi
if grep -Fq -- '--input-type' "$MAIN_NF"; then fail "Unsupported pVACfuse --input-type remains"; else pass "pVACfuse has no unsupported --input-type"; fi

main_contains_all "MultiQC branch" "process MULTIQC" "MULTIQC(qc_files)" "raw_qc.qc" "trimmed_qc.qc" "trim_result.reports" "star_result.logs" "flagstat" "md_result.metrics" "variant_stats"
main_contains_all "Combined FASTA branch" "process COMBINE_PROTEIN_FASTA" "\${proteome}" "\${variant_fasta}" "\${fusion_fasta}" "\${splice_fasta}" "exploratory_proteogenomics.fasta"

section "scratch.slurm contract"
contains_all "Shared Nextflow work directory" "$SLURM_FILE" \
    'SCRATCH_ROOT="/cluster/work/users/ash022/tmp/pgtk"' \
    'NXF_WORK="${SCRATCH_ROOT}/work"' \
    '-work-dir "${NXF_WORK}"'
contains_all "Launcher temporary environment" "$SLURM_FILE" \
    'export TMPDIR="${NXF_TMP}"' \
    'export TMP="${NXF_TMP}"' \
    'export TEMP="${NXF_TMP}"'
if grep -Fq -- 'singularity.runOptions' "$SLURM_FILE"; then
    fail "scratch.slurm contains unnecessary custom Singularity bind options"
else
    pass "scratch.slurm relies on the Nextflow work directory without custom binds"
fi
if grep -Fq -- '-resume' "$SLURM_FILE"; then
    warn "scratch.slurm enables resume"
else
    pass "scratch.slurm is configured for a fresh run"
fi
bash -n "$SLURM_FILE" && pass "scratch.slurm shell syntax valid" || fail "scratch.slurm shell syntax invalid"

section "Summary"
printf 'PASS: %d\nWARN: %d\nFAIL: %d\nREPORT: %s\n' "$PASS" "$WARN" "$FAIL" "$REPORT"
if (( FAIL > 0 )); then printf 'RESULT: FAILED\n'; exit 1; fi
printf 'RESULT: PASSED\n'
