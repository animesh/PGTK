#!/bin/bash
set -uo pipefail

PROJECT_DIR="${1:-$(pwd -P)}"
MAIN_NF="${PROJECT_DIR}/main.nf"
SLURM_FILE="${PROJECT_DIR}/scratch.slurm"
WORK_DIR="/cluster/work/users/ash022/work"
TMP_ROOT="/cluster/work/users/ash022/tmp"
CONTAINER_DIR="${PROJECT_DIR}/singularity_cache"
REFERENCE_DIR="${PROJECT_DIR}/reference_downloads"
REPORT="${PROJECT_DIR}/pipeline_command_validation_$(date +%Y%m%d_%H%M%S).txt"
PASS=0; WARN=0; FAIL=0
exec > >(tee "$REPORT") 2>&1
pass(){ PASS=$((PASS+1)); printf 'PASS  %s\n' "$*"; }
fail(){ FAIL=$((FAIL+1)); printf 'FAIL  %s\n' "$*"; }
section(){ printf '\n============================================================\n%s\n============================================================\n' "$*"; }
file_ok(){ [[ -s "$1" ]] && pass "File exists: $1" || fail "Missing: $1"; }
terms(){ local label="$1" file="$2"; shift 2; local bad=0 x; for x in "$@"; do grep -Fq -- "$x" "$file" || { fail "$label missing: $x"; bad=1; }; done; ((bad==0)) && pass "$label confirmed"; }

section "Environment"
file_ok "$MAIN_NF"; file_ok "$SLURM_FILE"; file_ok "$PROJECT_DIR/samples.csv"
command -v singularity >/dev/null && pass "Singularity found" || fail "Singularity missing"
[[ -d "$WORK_DIR" && -w "$WORK_DIR" ]] && pass "Work directory writable" || fail "Work directory unavailable: $WORK_DIR"
[[ -d "$TMP_ROOT" && -w "$TMP_ROOT" ]] && pass "Temporary root writable" || fail "Temporary root unavailable: $TMP_ROOT"

section "Containers"
for image in \
quay.io-biocontainers-sra-tools-3.2.1--h4304569_0.img \
quay.io-biocontainers-trim-galore-0.6.10--hdfd78af_0.img \
quay.io-biocontainers-fastqc-0.12.1--hdfd78af_0.img \
quay.io-biocontainers-star-2.7.11b--h43eeafb_1.img \
quay.io-biocontainers-gatk4-4.6.1.0--py310hdfd78af_0.img \
quay.io-biocontainers-samtools-1.21--h96c455f_1.img \
quay.io-biocontainers-ensembl-vep-111.0--pl5321h2a3209d_0.img \
quay.io-biocontainers-pypgatk-0.0.24--pyhdfd78af_0.img \
quay.io-biocontainers-arriba-2.4.0--h0033a41_2.img \
quay.io-biocontainers-bcftools-1.21--h8b25389_0.img \
quay.io-biocontainers-multiqc-1.35--pyhdfd78af_1.img \
stringtie-3.0.3.img transdecoder-6.0.0.img pvactools-7.1.1.img gffcompare-0.12.10.img; do
  path="$CONTAINER_DIR/$image"
  [[ -s "$path" ]] && singularity inspect "$path" >/dev/null 2>&1 && pass "Container valid: $image" || fail "Container invalid: $image"
done

section "References"
for asset in Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz Homo_sapiens.GRCh38.111.gtf.gz Homo_sapiens.GRCh38.cdna.all.fa.gz homo_sapiens_vep_111_GRCh38.tar.gz human_reviewed_isoforms.fasta.gz arriba_v2.4.0.tar.gz; do
  [[ -s "$REFERENCE_DIR/$asset" ]] && pass "Reference exists: $asset" || fail "Reference missing: $asset"
done

section "Workflow"
terms "fasterq" "$MAIN_NF" "--temp fasterq_tmp" "--outdir ." "rm -rf fasterq_tmp"
terms "MarkDuplicates" "$MAIN_NF" "process MARK_DUPLICATES" "cpus 20; memory '64 GB'" "-Xmx56g" "--MAX_RECORDS_IN_RAM 1000000"
terms "GATK" "$MAIN_NF" "SplitNCigarReads" "HaplotypeCaller" "GenotypeGVCFs" "VariantFiltration" "SelectVariants" "--java-options"
terms "STAR" "$MAIN_NF" "process STAR_ALIGN" "cpus 32; memory '256 GB'" "queue 'bigmem'"
terms "pypgatk" "$MAIN_NF" "--vcf" "--input_fasta" "--gene_annotations_gtf" "--annotation_field_name CSQ" "--af_field AF" "--include_consequences" "--output_proteindb"
terms "branches" "$MAIN_NF" "process ARRIBA" "process GFFCOMPARE_NOVEL" "process SPLICE_PROTEIN_FASTA" "process COMBINE_PROTEIN_FASTA" "process PROGRESSION_SUBTRACT" "process MULTIQC"

section "Launcher"
terms "storage" "$SLURM_FILE" 'WORK_DIR="/cluster/work/users/ash022/work"' 'TMP_ROOT="/cluster/work/users/ash022/tmp"' '-work-dir "${WORK_DIR}"'
terms "environment" "$SLURM_FILE" 'module purge' 'module load Java/21' 'NXF_SINGULARITY_CACHEDIR="${PROJECT_DIR}/singularity_cache"' '--sra_dir "${PROJECT_DIR}/sra_cache"'
terms "cluster submission" "$SLURM_FILE" '"-process.clusterOptions=--account=nn9036k"' '-resume'
terms "launcher resources" "$SLURM_FILE" '#SBATCH --cpus-per-task=4' '#SBATCH --mem=16G'
if grep -Eq -- 'module --force purge|-process[.]clusterOptions "--account|-with-report|-with-timeline|singularity[.]runOptions' "$SLURM_FILE"; then fail "Obsolete or malformed launcher option found"; else pass "Launcher options are clean"; fi
bash -n "$SLURM_FILE" && pass "Launcher shell syntax valid" || fail "Launcher syntax invalid"

section "Summary"
printf 'PASS: %d\nWARN: %d\nFAIL: %d\nREPORT: %s\n' "$PASS" "$WARN" "$FAIL" "$REPORT"
((FAIL==0)) && { printf 'RESULT: PASSED\n'; exit 0; }
printf 'RESULT: FAILED\n'; exit 1
