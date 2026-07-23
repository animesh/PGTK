#!/bin/bash
set -uo pipefail

WORKDIR="${1:-$(pwd -P)}"
MAIN_NF="${2:-$WORKDIR/main.nf}"
CONTAINER_DIR="$WORKDIR/singularity_cache"
REPORT="$WORKDIR/pipeline_command_validation_$(date +%Y%m%d_%H%M%S).txt"
TMP_ROOT="${TMPDIR:-$WORKDIR/tmp}/pgtk_command_validation_$$"

mkdir -p "$TMP_ROOT"
trap 'rm -rf "$TMP_ROOT"' EXIT

PASS=0
FAIL=0
WARN=0

exec > >(tee "$REPORT") 2>&1

pass() { PASS=$((PASS + 1)); printf 'PASS  %s\n' "$*"; }
fail() { FAIL=$((FAIL + 1)); printf 'FAIL  %s\n' "$*"; }
warn() { WARN=$((WARN + 1)); printf 'WARN  %s\n' "$*"; }

section() {
    printf '\n============================================================\n'
    printf '%s\n' "$*"
    printf '============================================================\n'
}

require_file() {
    if [[ -s "$1" ]]; then
        pass "File exists: $1"
    else
        fail "Missing or empty file: $1"
    fi
}

image_path() {
    printf '%s/%s' "$CONTAINER_DIR" "$1"
}

run_in_image() {
    local image="$1"
    shift
    singularity exec \
        --bind "$WORKDIR:$WORKDIR" \
        --pwd "$WORKDIR" \
        "$image" "$@"
}

check_image() {
    local image="$1"
    if [[ ! -s "$image" ]]; then
        fail "Missing container: $image"
        return 1
    fi
    if singularity inspect "$image" >/dev/null 2>&1; then
        pass "Container valid: $(basename "$image")"
        return 0
    fi
    fail "Container inspection failed: $image"
    return 1
}

capture_help() {
    local output_file="$1"
    local image="$2"
    shift 2
    set +e
    run_in_image "$image" "$@" >"$output_file" 2>&1
    local status=$?
    set -e
    if [[ -s "$output_file" ]]; then
        return 0
    fi
    return "$status"
}

check_help_terms() {
    local label="$1"
    local image="$2"
    local command_spec="$3"
    shift 3
    local output_file="$TMP_ROOT/${label//[^A-Za-z0-9_.-]/_}.txt"
    local -a command_array=()
    read -r -a command_array <<< "$command_spec"

    if ! capture_help "$output_file" "$image" "${command_array[@]}"; then
        fail "$label: command produced no usable help output"
        sed -n '1,40p' "$output_file" 2>/dev/null || true
        return
    fi

    local missing=0
    local term
    for term in "$@"; do
        if ! grep -Fq -- "$term" "$output_file"; then
            fail "$label: help does not contain '$term'"
            missing=1
        fi
    done
    if (( missing == 0 )); then
        pass "$label: executable and required arguments detected"
    else
        printf '      Help excerpt: %s\n' "$output_file"
    fi
}

check_main_terms() {
    local label="$1"
    shift
    local missing=0
    local term
    for term in "$@"; do
        if ! grep -Fq -- "$term" "$MAIN_NF"; then
            fail "$label: main.nf missing '$term'"
            missing=1
        fi
    done
    if (( missing == 0 )); then
        pass "$label: main.nf contains required command and arguments"
    fi
}

set -e

section "Environment"
printf 'Validation started: %s\n' "$(date --iso-8601=seconds)"
printf 'Work directory: %s\n' "$WORKDIR"
printf 'Workflow file: %s\n' "$MAIN_NF"
printf 'Report file: %s\n' "$REPORT"

if command -v singularity >/dev/null 2>&1; then
    pass "Singularity executable: $(command -v singularity)"
    singularity --version || true
else
    fail "singularity command not found"
fi

require_file "$MAIN_NF"
require_file "$WORKDIR/samples.csv"

SRA_IMG=$(image_path 'quay.io-biocontainers-sra-tools-3.2.1--h4304569_0.img')
TRIM_IMG=$(image_path 'quay.io-biocontainers-trim-galore-0.6.10--hdfd78af_0.img')
FASTQC_IMG=$(image_path 'quay.io-biocontainers-fastqc-0.12.1--hdfd78af_0.img')
STAR_IMG=$(image_path 'quay.io-biocontainers-star-2.7.11b--h43eeafb_1.img')
GATK_IMG=$(image_path 'quay.io-biocontainers-gatk4-4.6.1.0--py310hdfd78af_0.img')
SAMTOOLS_IMG=$(image_path 'quay.io-biocontainers-samtools-1.21--h96c455f_1.img')
VEP_IMG=$(image_path 'quay.io-biocontainers-ensembl-vep-111.0--pl5321h2a3209d_0.img')
PYPGATK_IMG=$(image_path 'quay.io-biocontainers-pypgatk-0.0.24--pyhdfd78af_0.img')
ARRIBA_IMG=$(image_path 'quay.io-biocontainers-arriba-2.4.0--h0033a41_2.img')
BCFTOOLS_IMG=$(image_path 'quay.io-biocontainers-bcftools-1.21--h8b25389_0.img')
MULTIQC_IMG=$(image_path 'quay.io-biocontainers-multiqc-1.35--pyhdfd78af_1.img')
STRINGTIE_IMG=$(image_path 'stringtie-3.0.3.img')
TRANSDECODER_IMG=$(image_path 'transdecoder-6.0.0.img')
PVACTOOLS_IMG=$(image_path 'pvactools-7.1.1.img')

section "Container integrity"
for image in \
    "$SRA_IMG" "$TRIM_IMG" "$FASTQC_IMG" "$STAR_IMG" "$GATK_IMG" \
    "$SAMTOOLS_IMG" "$VEP_IMG" "$PYPGATK_IMG" "$ARRIBA_IMG" \
    "$BCFTOOLS_IMG" "$MULTIQC_IMG" "$STRINGTIE_IMG" \
    "$TRANSDECODER_IMG" "$PVACTOOLS_IMG"
do
    check_image "$image" || true
done

section "Executable and argument checks"

check_help_terms "fasterq-dump" "$SRA_IMG" "fasterq-dump --help" \
    "--split-files" "--threads" "--temp" "--outdir"

check_help_terms "Trim Galore" "$TRIM_IMG" "trim_galore --help" \
    "--paired" "--quality" "--length" "--cores" "--gzip" "--basename"

check_help_terms "FastQC" "$FASTQC_IMG" "fastqc --help" \
    "--threads" "--outdir"

check_help_terms "STAR" "$STAR_IMG" "STAR --help" \
    "--runMode" "--genomeDir" "--genomeFastaFiles" "--sjdbGTFfile" \
    "--readFilesIn" "--readFilesCommand" "--twopassMode" \
    "--outSAMtype" "--chimOutType"

check_help_terms "samtools" "$SAMTOOLS_IMG" "samtools --help" \
    "faidx" "dict" "sort" "index" "flagstat"

check_help_terms "GATK tools" "$GATK_IMG" "gatk --list" \
    "MarkDuplicates" "SplitNCigarReads" "HaplotypeCaller" \
    "GenotypeGVCFs" "VariantFiltration" "SelectVariants"

check_help_terms "bcftools" "$BCFTOOLS_IMG" "bcftools --help" \
    "stats" "isec" "index"

check_help_terms "VEP" "$VEP_IMG" "vep --help" \
    "--input_file" "--output_file" "--cache" "--offline" \
    "--dir_cache" "--fasta" "--canonical" "--protein" "--hgvs" "--fork"

check_help_terms "pypgatk vcf-to-proteindb" "$PYPGATK_IMG" "pypgatk vcf-to-proteindb --help" \
    "--input-vcf" "--gene-annotations-gtf" "--protein-db-fasta" \
    "--af-field" "--annotation-field-name" "--consequence-filter" \
    "--output-proteindb"

check_help_terms "Arriba" "$ARRIBA_IMG" "arriba -h" \
    "-x" "-a" "-g" "-b" "-k" "-p" "-o" "-O"

check_help_terms "pVACfuse protein FASTA" "$PVACTOOLS_IMG" \
    "pvacfuse generate_protein_fasta --help" \
    "Arriba fusion.tsv" "flanking_sequence_length" \
    "--downstream-sequence-length"

check_help_terms "StringTie" "$STRINGTIE_IMG" "stringtie --help" \
    "-G" "-p" "-c" "-f" "-j" "-o"

check_help_terms "TransDecoder.LongOrfs" "$TRANSDECODER_IMG" \
    "TransDecoder.LongOrfs --help" "-t" "--output_dir"

check_help_terms "TransDecoder.Predict" "$TRANSDECODER_IMG" \
    "TransDecoder.Predict --help" "-t" "--output_dir"

check_help_terms "TransDecoder transcript extraction" "$TRANSDECODER_IMG" \
    "gtf_genome_to_cdna_fasta.pl" "usage"

check_help_terms "MultiQC" "$MULTIQC_IMG" "multiqc --help" \
    "--force" "--title" "--filename" "--data-dir" "--data-format"

section "Workflow command checks"

check_main_terms "SRA_TO_FASTQ" \
    "process SRA_TO_FASTQ" "fasterq-dump" "--split-files" "--threads"

if grep -Fq -- '--temp' "$MAIN_NF" && grep -Fq -- '--outdir' "$MAIN_NF"; then
    pass "SRA_TO_FASTQ explicitly sets temporary and output directories"
else
    warn "SRA_TO_FASTQ does not explicitly set both --temp and --outdir; fasterq-dump defaults to the task working directory"
fi

check_main_terms "TRIM_GALORE" \
    "process TRIM_GALORE" "--paired" "--quality 20" "--length 36" \
    "--cores \${task.cpus}" "--gzip" "--basename \${meta.sample}"

check_main_terms "STAR_INDEX" \
    "process STAR_INDEX" "--runMode genomeGenerate" "--genomeDir star_index" \
    "--genomeFastaFiles \${genome}" "--sjdbGTFfile \${gtf}" \
    "--sjdbOverhang \${params.read_length-1}" "--runThreadN \${task.cpus}"

check_main_terms "STAR_ALIGN" \
    "process STAR_ALIGN" "--twopassMode Basic" "--outSAMtype BAM Unsorted" \
    "--outSAMunmapped Within" "--outBAMcompression 0" \
    "--chimOutType WithinBAM HardClip" "--chimSegmentMin 10"

check_main_terms "SAMtools branch" \
    "samtools faidx genome.fa" "samtools dict -o genome.dict genome.fa" \
    "samtools sort" "samtools index" "samtools flagstat"

check_main_terms "GATK branch" \
    "gatk MarkDuplicates" "gatk SplitNCigarReads" "gatk HaplotypeCaller" \
    "gatk GenotypeGVCFs" "gatk VariantFiltration" "gatk SelectVariants"

check_main_terms "VEP all-transcript annotation" \
    "process VEP_ANNOTATE" "--cache" "--offline" "--canonical" \
    "--protein" "--symbol" "--numbers" "--biotype" "--total_length" "--hgvs"

if grep -Eq -- '(^|[[:space:]])--pick([[:space:]\\]|$)|--flag_pick' "$MAIN_NF"; then
    fail "VEP contains --pick or --flag_pick; exploratory all-transcript behavior is disabled"
else
    pass "VEP has no --pick or --flag_pick; all overlapping transcript consequences are retained"
fi

check_main_terms "pypgatk" \
    "pypgatk vcf-to-proteindb" "--annotation-field-name CSQ" \
    "--consequence-filter" "--output-proteindb"

check_main_terms "Arriba" \
    "process ARRIBA" "arriba -x" "-a \${genome}" "-g \${gtf}" \
    "-b \${blacklist}" "-k \${known}" "-p \${domains}"

check_main_terms "pVACfuse" \
    "process FUSION_FASTA" "pvacfuse generate_protein_fasta" \
    "--downstream-sequence-length full" "\${params.fusion_flank_aa}"

if grep -Fq -- '--input-type' "$MAIN_NF"; then
    fail "pVACfuse contains unsupported --input-type option"
else
    pass "pVACfuse does not contain unsupported --input-type option"
fi

check_main_terms "StringTie" \
    "process STRINGTIE_ASSEMBLY" "stringtie \${bam}" "-G \${gtf}" \
    "-c \${params.splice_min_coverage}" "-f \${params.splice_min_isoform_fraction}" \
    "-j \${params.splice_min_junction_reads}"

check_main_terms "TransDecoder" \
    "process SPLICE_PROTEIN_FASTA" "gtf_genome_to_cdna_fasta.pl" \
    "TransDecoder.LongOrfs" "TransDecoder.Predict" \
    "\${params.splice_min_protein_aa}"

check_main_terms "Combined FASTA" \
    "process COMBINE_PROTEIN_FASTA" "\${proteome}" "\${variant_fasta}" \
    "\${fusion_fasta}" "\${splice_fasta}" "exploratory_proteogenomics.fasta"

check_main_terms "MultiQC" \
    "process MULTIQC" "multiqc ." "--force" \
    "--filename multiqc_report.html" "--data-dir" "--data-format tsv" \
    "MULTIQC(qc_files)"

check_main_terms "MultiQC input channels" \
    "raw_qc.qc" "trimmed_qc.qc" "trim_result.reports" \
    "star_result.logs" "flagstat" "md_result.metrics" "variant_stats"

section "Container command versions"
for spec in \
    "$SRA_IMG|fasterq-dump --version" \
    "$TRIM_IMG|trim_galore --version" \
    "$FASTQC_IMG|fastqc --version" \
    "$STAR_IMG|STAR --version" \
    "$SAMTOOLS_IMG|samtools --version" \
    "$GATK_IMG|gatk --version" \
    "$BCFTOOLS_IMG|bcftools --version" \
    "$VEP_IMG|vep --version" \
    "$PYPGATK_IMG|pypgatk --version" \
    "$STRINGTIE_IMG|stringtie --version" \
    "$MULTIQC_IMG|multiqc --version"
do
    image=${spec%%|*}
    command_spec=${spec#*|}
    read -r -a command_array <<< "$command_spec"
    printf '\n[%s]\n' "$(basename "$image")"
    set +e
    run_in_image "$image" "${command_array[@]}" 2>&1 | sed -n '1,5p'
    set -e
done

section "Summary"
printf 'PASS: %d\n' "$PASS"
printf 'WARN: %d\n' "$WARN"
printf 'FAIL: %d\n' "$FAIL"
printf 'REPORT: %s\n' "$REPORT"

if (( FAIL > 0 )); then
    printf 'RESULT: FAILED\n'
    exit 1
fi

printf 'RESULT: PASSED\n'
exit 0
