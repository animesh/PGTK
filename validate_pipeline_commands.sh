#!/usr/bin/env bash
set -uo pipefail
PROJECT_DIR= NEXTFLOW=
while [[ $# -gt 0 ]]; do
    case "$1" in
        --project-dir) PROJECT_DIR=${2:?}; shift 2 ;;
        --nextflow) NEXTFLOW=${2:?}; shift 2 ;;
        --help|-h) printf 'Usage: %s --project-dir PATH --nextflow PATH\n' "$0"; exit 0 ;;
        *) printf 'ERROR: unknown option: %s\n' "$1" >&2; exit 2 ;;
    esac
done
[[ -n $PROJECT_DIR && -n $NEXTFLOW ]] || { printf 'ERROR: --project-dir and --nextflow are required\n' >&2; exit 2; }
MAIN_NF="$PROJECT_DIR/main.nf"
REPORT="$PROJECT_DIR/pipeline_command_validation_$(date +%Y%m%d_%H%M%S).txt"
PASS=0; FAIL=0
exec > >(tee "$REPORT") 2>&1
pass(){ PASS=$((PASS+1)); printf 'PASS  %s\n' "$*"; }
fail(){ FAIL=$((FAIL+1)); printf 'FAIL  %s\n' "$*"; }
need(){ [[ -s $1 ]] && pass "File: $1" || fail "Missing or empty: $1"; }
has(){ grep -Fq -- "$2" "$1" && pass "$3" || fail "$3: missing $2"; }
lacks(){ grep -Fq -- "$2" "$1" && fail "$3: forbidden $2" || pass "$3"; }
FILES=(main.nf nextflow.config scratch.slurm download_assets.sh samples.csv collect_pipeline_failures.py analyze_pipeline_trace.py test_resource_configuration.py validate_rna_events.py build_complete_report.py map_peptides_to_fasta.py annotate_variant_peptides.py analyze_chimeric_splice_peptides.py validate_splice_junction_peptides.py proteogenomics_evidence_report.py validate_proteogenomic_reads.py validate_variant_read_provenance.py validate_variant_codons.py merge_variant_validation.py analyze_codon_mismatches.py build_integrated_variant_evidence.py PIPELINE_VALIDATION_SEMANTICS.md RESOURCE_CALIBRATION.md RESOURCE_RETRY_MATRIX.tsv HISTORY_REGRESSION_CHECKLIST.md maxquant_raw_file_map.none.tsv audit_environment_hardcoding.py validate_runtime_inputs.py validate_haplotype_shards.py summarize_variant_stages.py compare_external_vcf.py build_comparative_advantage_report.py build_igv_evidence_bundle.py test_igv_evidence_bundle.py prepare_go_annotations.py analyze_progression_biology.py compare_progression_pair.py merge_progression_biology.py test_progression_biology.py)
for file in "${FILES[@]}"; do need "$PROJECT_DIR/$file"; done
for file in "$PROJECT_DIR"/*.py; do python -m py_compile "$file" && pass "Python syntax: $(basename "$file")" || fail "Python syntax: $(basename "$file")"; done
bash -n "$PROJECT_DIR/scratch.slurm" && pass 'scratch.slurm syntax' || fail 'scratch.slurm syntax'
bash -n "$PROJECT_DIR/download_assets.sh" && pass 'download_assets.sh syntax' || fail 'download_assets.sh syntax'
bash -n "$0" && pass 'validator syntax' || fail 'validator syntax'
mapfile -t declared < <(awk '/^process / {print $2}' "$MAIN_NF")
[[ ${#declared[@]} -eq 57 ]] && pass '57 process declarations' || fail "Expected 57 processes; found ${#declared[@]}"
[[ $(printf '%s\n' "${declared[@]}" | sort | uniq -d | wc -l) -eq 0 ]] && pass 'No duplicate process names' || fail 'Duplicate process names'
for name in "${declared[@]}"; do [[ $(grep -c "^process $name " "$MAIN_NF") -eq 1 ]] && pass "Process declared once: $name" || fail "Process declaration error: $name"; done
has "$MAIN_NF" 'process DOWNLOAD_REFERENCES {' 'Reference preparation process present'
has "$MAIN_NF" 'path genome_archive' 'Reference genome staged into container'
has "$MAIN_NF" 'path vep_archive' 'VEP cache staged into container'
has "$MAIN_NF" 'reference_arriba_archive = file(' 'Arriba archive validated before scheduling'
has "$PROJECT_DIR/scratch.slurm" 'validate_runtime_inputs.py' 'Runtime preflight runs before Nextflow'
has "$PROJECT_DIR/scratch.slurm" 'APPTAINER_BINDPATH=' 'Runtime paths supplied to Apptainer dynamically'
has "$PROJECT_DIR/analyze_pipeline_trace.py" "text in {'', '-', 'NA', 'N/A'}" 'Incomplete trace values handled'
has "$MAIN_NF" 'process SPLIT_N_CIGAR {' 'Monolithic SplitNCigarReads retained'
has "$MAIN_NF" 'process VALIDATE_HAPLOTYPE_SHARDS {' 'Haplotype shard validation wired'
has "$MAIN_NF" 'process VARIANT_STAGE_QC {' 'Variant stage QC wired'
has "$MAIN_NF" 'process COMPARE_EXTERNAL_VCF {' 'External VCF comparison wired'
has "$MAIN_NF" 'process BUILD_COMPARATIVE_ADVANTAGE_REPORT {' 'Comparative biological report wired'
has "$MAIN_NF" 'process BUILD_IGV_EVIDENCE_BUNDLE {' 'RNA and progression IGV bundle wired'
has "$MAIN_NF" 'process PREPARE_GO_ANNOTATIONS {' 'Versioned GO annotation preparation wired'
has "$PROJECT_DIR/download_assets.sh" 'go-basic.obo' 'GO ontology download integrated'
has "$PROJECT_DIR/download_assets.sh" 'goa_human.gaf.gz' 'Human GO annotation download integrated'
has "$PROJECT_DIR/download_assets.sh" 'validate_gaf_gzip' 'GO annotation validation integrated'
has "$PROJECT_DIR/download_assets.sh" 'downloaded_assets.sha256' 'Reference asset checksums integrated'
has "$MAIN_NF" 'process ANALYZE_PROGRESSION_SAMPLE {' 'Parallel per-sample GO analysis wired'
has "$MAIN_NF" 'process COMPARE_PROGRESSION_PAIR {' 'Parallel pairwise GO contrast wired'
has "$MAIN_NF" 'process MERGE_PROGRESSION_BIOLOGY {' 'GO result merge wired'
has "$MAIN_NF" 'pairwise_go_contrasts.tsv' 'Pairwise GO contrasts wired'
has "$MAIN_NF" 'progression_biology.go_enrichment.tsv' 'RNA-callable-background GO enrichment wired'
lacks "$PROJECT_DIR/analyze_progression_biology.py" 'DEFAULT_CATEGORIES' 'No hard-coded biological categories'
has "$PROJECT_DIR/prepare_go_annotations.py" "v.startswith('part_of ')" 'GO part_of propagation implemented'
has "$MAIN_NF" '"${params.host_python}" ${bundle_script}' 'IGV bundle uses configured host Python inside samtools container'
has "$MAIN_NF" 'process PREPARE_COMPARATIVE_MULTIQC_CONTENT {' 'Comparative MultiQC content wired'
has "$MAIN_NF" 'shared_with_baseline.vep.vcf.gz' 'Shared baseline subtraction VCF wired'
has "$MAIN_NF" 'txtMQMBR' 'Default MaxQuant results folder wired'
has "$MAIN_NF" 'params.run_proteogenomic_validation = false' 'MaxQuant validation disabled by default'
has "$MAIN_NF" 'params.run_external_vcf_comparison = false' 'External comparison disabled by default'
has "$MAIN_NF" 'cat ${variant_fasta} ${fusion_fasta} ${splice_fasta}' 'Custom FASTA excludes canonical proteome'
lacks "$MAIN_NF" 'cat ${proteome} ${variant_fasta}' 'Canonical proteome not embedded in custom FASTA'
has "$MAIN_NF" '${projectDir}/sarek' 'Default Sarek results folder wired'
has "$MAIN_NF" 'publishDir "${params.outdir}/vcf_raw"' 'Raw VCF publication wired'
has "$MAIN_NF" 'emit: raw' 'Raw VCF output channel'
has "$MAIN_NF" 'row.TK ?: sample' 'Optional TK defaults to sample'
lacks "$MAIN_NF" 'SPLIT_N_CIGAR_SHARD' 'No scattered SplitNCigarReads'
lacks "$MAIN_NF" "queue 'normal'" 'No hard-coded normal queue in main.nf'
lacks "$MAIN_NF" "queue 'bigmem'" 'No hard-coded bigmem queue in main.nf'
has "$MAIN_NF" 'params.host_python = null' 'Host Python supplied as parameter'
has "$MAIN_NF" 'params.reference_downloads = null' 'Reference path supplied as parameter'
has "$MAIN_NF" 'params.container_cache = null' 'Container cache supplied as parameter'
has "$PROJECT_DIR/nextflow.config" 'task.cpus > pgtkEffectiveNormalCpuThreshold || task.memory > pgtkEffectiveNormalMemoryThresholdGb.GB' 'Dynamic partition routing'
has "$PROJECT_DIR/nextflow.config" "requiredEnv('PGTK_NORMAL_PARTITION')" 'Normal partition supplied at runtime'
has "$PROJECT_DIR/nextflow.config" "requiredEnv('PGTK_BIGMEM_PARTITION')" 'Bigmem partition supplied at runtime'
has "$PROJECT_DIR/nextflow.config" 'Math.min(pgtkMaxCpus' 'Dynamic CPU cap'
has "$PROJECT_DIR/nextflow.config" 'pgtkEffectiveMaxMemoryGb.GB' 'Dynamic memory cap'
has "$PROJECT_DIR/nextflow.config" 'pgtkAbsoluteMaxCpus = 32' 'Absolute CPU ceiling'
has "$PROJECT_DIR/nextflow.config" 'pgtkAbsoluteMaxMemoryGb = 512' 'Absolute memory ceiling'
has "$PROJECT_DIR/nextflow.config" 'pgtkAbsoluteNormalCpuThreshold = 20' 'Normal-partition CPU ceiling'
has "$PROJECT_DIR/nextflow.config" 'pgtkAbsoluteNormalMemoryThresholdGb = 160' 'Normal-partition memory ceiling'
has "$PROJECT_DIR/nextflow.config" 'maxRetries = 2' 'Three total attempts configured'
has "$PROJECT_DIR/nextflow.config" "task.exitStatus in [137, 140, 143] ? 'retry' : 'terminate'" 'Only resource-related exit codes are retried'
lacks "$PROJECT_DIR/nextflow.config" "queue = 'normal'" 'No fixed global partition'
lacks "$PROJECT_DIR/nextflow.config" "queue = 'bigmem'" 'No fixed bigmem partition'
has "$PROJECT_DIR/scratch.slurm" '#SBATCH --account=nn9036k' 'Default Saga account directive'
has "$PROJECT_DIR/scratch.slurm" 'PROJECT_DIR=/cluster/projects/nn9036k/scrbkup/pgtk' 'Default Saga project directory'
has "$PROJECT_DIR/scratch.slurm" 'ACCOUNT=nn9036k' 'Default Saga runtime account'
has "$PROJECT_DIR/scratch.slurm" 'JAVA_MODULE=Java/21' 'Default Saga Java module'
has "$PROJECT_DIR/scratch.slurm" 'PIPELINE_ARGS=("$@")' 'Pipeline arguments forwarded after separator'
has "$PROJECT_DIR/scratch.slurm" 'trap finalize EXIT' 'Failure collector always runs'
has "$PROJECT_DIR/collect_pipeline_failures.py" "'Partition':pick(row,'queue')" 'Failure ledger records partition'
has "$MAIN_NF" 'variant_codon_validation.partial.tsv' 'Partial codon output preserved'
has "$MAIN_NF" '17_validation_semantics_mqc.html' '17-section MultiQC preserved'
has "$MAIN_NF" 'task.memory.toGiga() * 0.80' 'GATK Java heap follows retry memory'
has "$MAIN_NF" '-Xmx${javaHeapGb}g' 'Dynamic GATK maximum heap wired'
lacks "$MAIN_NF" '-Xmx40g' 'No fixed MarkDuplicates heap'
lacks "$MAIN_NF" '-Xmx18g' 'No fixed SplitNCigarReads heap'
has "$MAIN_NF" 'ParallelGCThreads=${javaGcThreads}' 'MarkDuplicates GC threads follow allocated CPUs'
has "$PROJECT_DIR/nextflow.config" 'withName: SPLIT_N_CIGAR {' 'SplitNCigarReads resource selector present'
python "$PROJECT_DIR/audit_environment_hardcoding.py" "$PROJECT_DIR" && pass 'No environment-specific hard-coding' || fail 'Environment-specific hard-coding audit'
python "$PROJECT_DIR/test_resource_configuration.py" && pass 'Resource and wrapper fixtures' || fail 'Resource and wrapper fixtures'
python "$PROJECT_DIR/test_igv_evidence_bundle.py" && pass 'IGV evidence fixture' || fail 'IGV evidence fixture'
python "$PROJECT_DIR/test_progression_biology.py" && pass 'Scalable progression biology fixture' || fail 'Scalable progression biology fixture'
TMP_COMPARE=$(mktemp -d)
cat > "$TMP_COMPARE/pgtk.vcf" <<'VCF'
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
1	10	.	A	G	.	PASS	.	GT	0/1
1	20	.	C	T	.	PASS	.	GT	1/1
VCF
cat > "$TMP_COMPARE/external.vcf" <<'VCF'
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	10	.	A	G	.	PASS	.	GT	0/1
chr1	30	.	G	A	.	PASS	.	GT	0/1
VCF
python "$PROJECT_DIR/compare_external_vcf.py" --sample fixture --stage raw --pgtk "$TMP_COMPARE/pgtk.vcf" --external "$TMP_COMPARE/external.vcf" --output-prefix "$TMP_COMPARE/comparison" >/dev/null 2>&1 \
    && [[ -s "$TMP_COMPARE/comparison.summary.tsv" && -s "$TMP_COMPARE/comparison.shared.tsv" && -s "$TMP_COMPARE/comparison.pgtk_only.tsv" && -s "$TMP_COMPARE/comparison.external_only.tsv" && -s "$TMP_COMPARE/comparison.report.md" ]] \
    && pass 'External comparison output contract fixture' || fail 'External comparison output contract fixture'
rm -rf "$TMP_COMPARE"
if [[ -x $NEXTFLOW ]]; then
    PGTK_ACCOUNT=validation PGTK_NORMAL_PARTITION=normal PGTK_BIGMEM_PARTITION=bigmem PGTK_NORMAL_CPU_THRESHOLD=20 PGTK_NORMAL_MEMORY_THRESHOLD_GB=160 PGTK_MAX_CPUS=32 PGTK_MAX_MEMORY_GB=512 PGTK_QUEUE_SIZE=200 PGTK_SUBMIT_RATE_LIMIT=60/1min \
        "$NEXTFLOW" inspect "$MAIN_NF" >/tmp/pgtk_dynamic_inspect.txt 2>&1 \
        && pass 'Nextflow inspect' || { cat /tmp/pgtk_dynamic_inspect.txt; fail 'Nextflow inspect'; }
else
    fail "Nextflow executable missing: $NEXTFLOW"
fi
printf '\nPASS: %d\nFAIL: %d\nREPORT: %s\n' "$PASS" "$FAIL" "$REPORT"
(( FAIL == 0 )) || exit 1
printf 'RESULT: PASSED\n'
