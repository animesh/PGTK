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
FILES=(main.nf nextflow.config scratch.slurm download_assets.sh samples.csv collect_pipeline_failures.py analyze_pipeline_trace.py test_resource_configuration.py validate_rna_events.py build_complete_report.py map_peptides_to_fasta.py annotate_variant_peptides.py analyze_chimeric_splice_peptides.py validate_splice_junction_peptides.py proteogenomics_evidence_report.py validate_proteogenomic_reads.py validate_variant_read_provenance.py validate_variant_codons.py merge_variant_validation.py analyze_codon_mismatches.py build_integrated_variant_evidence.py PIPELINE_VALIDATION_SEMANTICS.md RESOURCE_CALIBRATION.md RESOURCE_RETRY_MATRIX.tsv HISTORY_REGRESSION_CHECKLIST.md maxquant_raw_file_map.none.tsv audit_environment_hardcoding.py validate_runtime_inputs.py validate_haplotype_shards.py summarize_variant_stages.py compare_external_vcf.py build_comparative_advantage_report.py build_igv_evidence_bundle.py build_finding_igv_reviews.py test_igv_evidence_bundle.py prepare_go_annotations.py analyze_progression_biology.py compare_progression_pair.py merge_progression_biology.py test_progression_biology.py expression_go_analysis.py test_expression_go_analysis.py)
for file in "${FILES[@]}"; do need "$PROJECT_DIR/$file"; done
for file in "$PROJECT_DIR"/*.py; do python -m py_compile "$file" && pass "Python syntax: $(basename "$file")" || fail "Python syntax: $(basename "$file")"; done
bash -n "$PROJECT_DIR/scratch.slurm" && pass 'scratch.slurm syntax' || fail 'scratch.slurm syntax'
bash -n "$PROJECT_DIR/download_assets.sh" && pass 'download_assets.sh syntax' || fail 'download_assets.sh syntax'
bash -n "$0" && pass 'validator syntax' || fail 'validator syntax'
mapfile -t declared < <(awk '/^process / {print $2}' "$MAIN_NF")
[[ ${#declared[@]} -eq 72 ]] && pass '72 process declarations' || fail "Expected 72 processes; found ${#declared[@]}"
[[ $(printf '%s\n' "${declared[@]}" | sort | uniq -d | wc -l) -eq 0 ]] && pass 'No duplicate process names' || fail 'Duplicate process names'
for name in "${declared[@]}"; do [[ $(grep -c "^process $name " "$MAIN_NF") -eq 1 ]] && pass "Process declared once: $name" || fail "Process declaration error: $name"; done
has "$MAIN_NF" 'process DOWNLOAD_REFERENCES {' 'Reference preparation process present'
has "$MAIN_NF" 'path genome_archive' 'Reference genome staged into container'
has "$MAIN_NF" 'path vep_archive' 'VEP cache staged into container'
has "$MAIN_NF" 'reference_arriba_archive = file(' 'Arriba archive validated before scheduling'
has "$PROJECT_DIR/scratch.slurm" 'validate_runtime_inputs.py' 'Runtime preflight runs before Nextflow'
lacks "$PROJECT_DIR/scratch.slurm" 'APPTAINER_BINDPATH' 'No global Apptainer bind path'
has "$PROJECT_DIR/analyze_pipeline_trace.py" "text in {'', '-', 'NA', 'N/A'}" 'Incomplete trace values handled'
has "$MAIN_NF" 'process SPLIT_N_CIGAR {' 'Monolithic SplitNCigarReads retained'
has "$MAIN_NF" 'process VALIDATE_HAPLOTYPE_SHARDS {' 'Haplotype shard validation wired'
has "$MAIN_NF" 'process VARIANT_STAGE_QC {' 'Variant stage QC wired'
has "$MAIN_NF" 'process COMPARE_EXTERNAL_VCF {' 'External VCF comparison wired'
has "$MAIN_NF" 'process BUILD_COMPARATIVE_ADVANTAGE_REPORT {' 'Comparative biological report wired'
has "$MAIN_NF" 'process BUILD_IGV_EVIDENCE_BUNDLE {' 'RNA and progression IGV bundle wired'
has "$MAIN_NF" 'process BUILD_FINDING_IGV_REVIEWS {' 'Gene-based strict finding IGV reviews wired'
has "$MAIN_NF" 'process GENERATE_PRIORITY_IGV_REPORTS {' 'Offline prioritized IGV Reports generation wired'
has "$MAIN_NF" 'igv-reports-1.16.0--pyh7e72e81_0.img' 'Pinned IGV Reports container wired'
has "$MAIN_NF" 'params.generate_priority_igv_reports = true' 'Offline IGV Reports enabled by default'
has "$MAIN_NF" 'params.igv_report_limit = 0' 'IGV report generation defaults to all selected findings'
has "$MAIN_NF" 'params.finding_priority_mode = '"'"'all'"'"'' 'Finding selection mode is CLI-configurable'
has "$MAIN_NF" 'params.finding_priority_genes = '"'"''"'"'' 'No gene is prioritized by default'
has "$MAIN_NF" 'params.igv_report_timeout_seconds = 600' 'IGV report timeout is CLI-configurable'
has "$MAIN_NF" 'command -v create_report || command -v create_reports' 'IGV Reports command discovery wired'
has "$MAIN_NF" '--fasta ${genome}' 'IGV Reports uses local indexed FASTA'
has "$MAIN_NF" '--tracks "\$bed" "\$alt_bam" "\$ref_bam"' 'IGV Reports embeds finding tracks'
has "$MAIN_NF" 'report_manifest.tsv' 'IGV Reports manifest emitted'
lacks "$PROJECT_DIR/build_finding_igv_reviews.py" 'DDX1' 'No gene-specific prioritization is hard-coded'
has "$PROJECT_DIR/build_finding_igv_reviews.py" "--priority-mode" 'Generic priority mode exposed'
has "$PROJECT_DIR/build_finding_igv_reviews.py" "--priority-genes" 'Generic priority genes exposed'
python "$PROJECT_DIR/test_finding_igv_reviews.py" && pass 'Finding identifier, equivalent insertion and ALT-selection behavior' || fail 'Finding review behavior fixture'
has "$MAIN_NF" 'process COUNT_GENES_PER_SAMPLE {' 'Per-sample featureCounts expression process wired'
has "$MAIN_NF" 'process MERGE_GENE_EXPRESSION {' 'Merged raw count, CPM and TPM matrix wired'
if grep -Eq '^[[:space:]]*tuple val\(meta\)[[:space:]]*$' "$MAIN_NF"; then fail 'No single-element tuple inputs'; else pass 'No single-element tuple inputs'; fi
has "$MAIN_NF" 'expression_ora_inputs = expression_metadata' 'Expression sample GO receives metadata maps directly'
lacks "$MAIN_NF" 'expression_ora_inputs = expression_metadata.map { meta -> tuple(meta) }' 'Expression metadata is not wrapped in a one-element tuple'
has "$MAIN_NF" 'process ANALYZE_EXPRESSION_SAMPLE_GO {' 'Parallel per-sample expression ORA wired'
has "$MAIN_NF" 'process ANALYZE_EXPRESSION_RANKED_GO {' 'Parallel baseline-comparison ranked GO wired'
has "$MAIN_NF" 'process MERGE_EXPRESSION_GO {' 'Expression GO merge wired'
has "$MAIN_NF" 'process PREPARE_EXPRESSION_MULTIQC_CONTENT {' 'Expression and variant-set GO MultiQC content wired'
has "$MAIN_NF" 'comparative_multiqc.mix(expression_multiqc_content).collect()' 'Expression GO included in final MultiQC'
lacks "$MAIN_NF" 'process ANALYZE_EXPRESSION_GO {' 'Sequential expression GO process removed'
has "$MAIN_NF" 'params.go_fdr_threshold = 0.1' 'Default GO FDR threshold is 0.1'
has "$MAIN_NF" 'params.expression_rank_min_nonzero_scores = 1' 'Ranked GO non-zero-score safeguard is CLI-configurable'
has "$MAIN_NF" '--min-nonzero-scores ${params.expression_rank_min_nonzero_scores}' 'Ranked GO safeguard passed through Nextflow CLI'
has "$PROJECT_DIR/expression_go_analysis.py" "baseline_column = f'{args.baseline_sample}_TPM'" 'Ranked GO resolves the baseline TPM column explicitly'
has "$PROJECT_DIR/expression_go_analysis.py" 'float(row[baseline_column])' 'Ranked GO denominator uses baseline TPM'
has "$PROJECT_DIR/expression_go_analysis.py" '--sample and --baseline-sample must be different' 'Ranked GO rejects self-comparisons'
has "$PROJECT_DIR/expression_go_analysis.py" "'NonZeroScores': nonzero_scores" 'Ranked GO summary records score diagnostics'
[[ $(grep -Fc -- '--fdr-threshold ${params.go_fdr_threshold}' "$MAIN_NF") -eq 4 ]] && pass 'All four GO analysis paths use the CLI-driven FDR threshold' || fail 'All four GO analysis paths use the CLI-driven FDR threshold'
has "$PROJECT_DIR/expression_go_analysis.py" 'default=0.1' 'Standalone expression and variant GO default FDR is 0.1'
has "$PROJECT_DIR/expression_go_analysis.py" '_all' 'Complete per-progression-sample variant GO analyses wired'
has "$PROJECT_DIR/expression_go_analysis.py" "'FDRThreshold': args.fdr_threshold" 'Expression and variant GO summaries record applied FDR threshold'
has "$PROJECT_DIR/analyze_progression_biology.py" "'FDRThreshold':a.fdr_threshold" 'Progression GO summary records applied FDR threshold'
lacks "$PROJECT_DIR/analyze_progression_biology.py" 'SignificantGOTermsFDR05' 'No hard-coded FDR 0.05 progression summary field'
lacks "$PROJECT_DIR/merge_progression_biology.py" 'FDR <= 0.05' 'No hard-coded FDR threshold in progression report'
lacks "$PROJECT_DIR/expression_go_analysis.py" 'scipy' 'GO implementation has no optional Python package dependency'
lacks "$PROJECT_DIR/expression_go_analysis.py" 'math.comb' 'Slow combinatorial GO implementation removed'
has "$PROJECT_DIR/validate_runtime_inputs.py" 'CONTAINER_PYTHON_STDLIB_OK' 'Container-native Python preflight wired'
has "$PROJECT_DIR/validate_runtime_inputs.py" "'*multiqc-1.35*img'" 'MultiQC runtime image discovery wired'
has "$PROJECT_DIR/validate_runtime_inputs.py" "expression_go_analysis.py'), '--help'" 'Exact MultiQC expression-GO Python contract wired'
has "$PROJECT_DIR/validate_runtime_inputs.py" "analyze_progression_biology.py'), '--help'" 'Exact MultiQC progression-GO Python contract wired'
lacks "$PROJECT_DIR/analyze_progression_biology.py" 'scipy' 'Progression GO has no optional Python package dependency'
has "$MAIN_NF" 'process ANALYZE_PROGRESSION_VARIANT_SETS {' 'Common and sample-exclusive variant GO wired'
has "$MAIN_NF" 'process PREPARE_GO_ANNOTATIONS {' 'Versioned GO annotation preparation wired'
has "$PROJECT_DIR/download_assets.sh" 'go-basic.obo' 'GO ontology download integrated'
has "$PROJECT_DIR/download_assets.sh" 'goa_human.gaf.gz' 'Human GO annotation download integrated'
has "$PROJECT_DIR/download_assets.sh" 'validate_gaf_gzip' 'GO annotation validation integrated'
has "$PROJECT_DIR/download_assets.sh" 'downloaded_assets.sha256' 'Reference asset checksums integrated'
has "$PROJECT_DIR/download_assets.sh" "'subread-2.0.8.img'" 'Stable Subread image download integrated'
has "$PROJECT_DIR/download_assets.sh" 'featureCounts -v' 'featureCounts executable validation integrated'
has "$PROJECT_DIR/download_assets.sh" 'igv-reports:1.16.0--pyh7e72e81_0' 'Pinned IGV Reports image download integrated'
lacks "$PROJECT_DIR/download_assets.sh" 'igv-xvfb' 'Legacy network-dependent IGV Desktop image removed'
has "$PROJECT_DIR/download_assets.sh" 'pysam:0.24.0--py312hf5ad864_1' 'Pinned Pysam image integrated'
has "$PROJECT_DIR/validate_runtime_inputs.py" 'pgtk-igv-reports-smoke-' 'IGV Reports runtime smoke test wired'
has "$PROJECT_DIR/validate_runtime_inputs.py" "require_file(report, 'self-contained IGV Reports smoke-test HTML')" 'IGV Reports HTML required by preflight'
lacks "$PROJECT_DIR/validate_runtime_inputs.py" 'xvfb-run' 'Runtime preflight has no Xvfb dependency'
has "$MAIN_NF" '${params.container_cache}/subread-2.0.8.img' 'featureCounts uses validated local Subread image'
lacks "$MAIN_NF" 'subread:2.0.8--he4a0461_1' 'Invalid Subread tag absent from pipeline'
has "$MAIN_NF" 'process ANALYZE_PROGRESSION_SAMPLE {' 'Parallel per-sample GO analysis wired'
has "$MAIN_NF" 'process COMPARE_PROGRESSION_PAIR {' 'Parallel pairwise GO contrast wired'
has "$MAIN_NF" 'process MERGE_PROGRESSION_BIOLOGY {' 'GO result merge wired'
has "$MAIN_NF" 'pairwise_go_contrasts.tsv' 'Pairwise GO contrasts wired'
has "$MAIN_NF" 'progression_biology.go_enrichment.tsv' 'RNA-callable-background GO enrichment wired'
lacks "$PROJECT_DIR/analyze_progression_biology.py" 'DEFAULT_CATEGORIES' 'No hard-coded biological categories'
has "$PROJECT_DIR/prepare_go_annotations.py" "v.startswith('part_of ')" 'GO part_of propagation implemented'
has "$MAIN_NF" 'container "${params.pysam_image}"' 'IGV processes use pinned Pysam image'
has "$MAIN_NF" 'python3 ${bundle_script}' 'IGV bundle uses container Python'
has "$MAIN_NF" 'process PREPARE_COMPARATIVE_MULTIQC_CONTENT {' 'Comparative MultiQC content wired'
has "$MAIN_NF" 'shared_with_baseline.vep.vcf.gz' 'Shared baseline subtraction VCF wired'
has "$MAIN_NF" 'txtMQMBR' 'Default MaxQuant results folder wired'
has "$MAIN_NF" 'params.run_proteogenomic_validation = false' 'MaxQuant validation disabled by default'
has "$MAIN_NF" 'params.run_external_vcf_comparison = false' 'External comparison disabled by default'
has "$MAIN_NF" 'cat ${variant_fasta} ${fusion_fasta} ${splice_fasta}' 'Custom FASTA excludes canonical proteome'
lacks "$MAIN_NF" 'cat ${proteome} ${variant_fasta}' 'Canonical proteome not embedded in custom FASTA'
has "$MAIN_NF" '${projectDir}/sarek' 'Default Sarek results folder wired'
has "$MAIN_NF" 'publishDir "${params.outdir}/vcf_raw"' 'Raw VCF publication wired'
has "$MAIN_NF" 'process GENOTYPE_VARIANTS {' 'Raw GenotypeGVCFs process wired'
has "$MAIN_NF" 'path("${meta.sample}.raw.vcf.gz")' 'Raw VCF output channel wired'
has "$MAIN_NF" 'row.TK ?: sample' 'Optional TK defaults to sample'
lacks "$MAIN_NF" 'SPLIT_N_CIGAR_SHARD' 'No scattered SplitNCigarReads'
lacks "$MAIN_NF" "queue 'normal'" 'No hard-coded normal queue in main.nf'
lacks "$MAIN_NF" "queue 'bigmem'" 'No hard-coded bigmem queue in main.nf'
lacks "$MAIN_NF" 'params.host_python' 'No host Python Nextflow parameter'
has "$MAIN_NF" 'params.pysam_image = null' 'Pysam image supplied as parameter'
has "$MAIN_NF" 'params.reference_downloads = null' 'Reference path supplied as parameter'
has "$MAIN_NF" 'params.container_cache = null' 'Container cache supplied as parameter'
if grep -Eq '^def[[:space:]]+' "$PROJECT_DIR/nextflow.config"; then fail 'No top-level scripting declarations in Nextflow config'; else pass 'No top-level scripting declarations in Nextflow config'; fi
lacks "$PROJECT_DIR/nextflow.config" 'requiredEnv' 'No unsupported config helper closure'
has "$PROJECT_DIR/nextflow.config" "task.cpus > Math.min((env('PGTK_NORMAL_CPU_THRESHOLD') as int), 20)" 'Dynamic CPU partition routing uses official env()'
has "$PROJECT_DIR/nextflow.config" "task.memory > Math.min((env('PGTK_NORMAL_MEMORY_THRESHOLD_GB') as int), 160).GB" 'Dynamic memory partition routing uses official env()'
has "$PROJECT_DIR/nextflow.config" "env('PGTK_NORMAL_PARTITION')" 'Normal partition supplied at runtime'
has "$PROJECT_DIR/nextflow.config" "env('PGTK_BIGMEM_PARTITION')" 'Bigmem partition supplied at runtime'
has "$PROJECT_DIR/nextflow.config" "Math.min((env('PGTK_MAX_CPUS') as int), 32)" 'Dynamic CPU cap'
has "$PROJECT_DIR/nextflow.config" "Math.min((env('PGTK_MAX_MEMORY_GB') as int), 512).GB" 'Dynamic memory cap'
has "$PROJECT_DIR/nextflow.config" "queueSize = env('PGTK_QUEUE_SIZE') as int" 'Queue size supplied at runtime'
has "$PROJECT_DIR/nextflow.config" "submitRateLimit = env('PGTK_SUBMIT_RATE_LIMIT')" 'Submission rate supplied at runtime'
has "$PROJECT_DIR/nextflow.config" 'maxRetries = 2' 'Three total attempts configured'
has "$PROJECT_DIR/nextflow.config" "task.exitStatus in [137, 140, 143] ? 'retry' : 'terminate'" 'Only resource-related exit codes are retried'
lacks "$PROJECT_DIR/nextflow.config" "queue = 'normal'" 'No fixed global partition'
lacks "$PROJECT_DIR/nextflow.config" "queue = 'bigmem'" 'No fixed bigmem partition'
has "$PROJECT_DIR/scratch.slurm" '#SBATCH --account=nn9036k' 'Default Saga account directive'
has "$PROJECT_DIR/scratch.slurm" 'PROJECT_DIR=${SLURM_SUBMIT_DIR:-$(pwd -P)}' 'Project defaults to Slurm submission directory'
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
has "$PROJECT_DIR/nextflow.config" 'def value = 16.GB * (1 << (task.attempt - 1))' 'IGV evidence bundle starts at observed-safe 16 GB'
python "$PROJECT_DIR/audit_environment_hardcoding.py" "$PROJECT_DIR" && pass 'No environment-specific hard-coding' || fail 'Environment-specific hard-coding audit'
python "$PROJECT_DIR/test_resource_configuration.py" && pass 'Resource and wrapper fixtures' || fail 'Resource and wrapper fixtures'
pass 'IGV evidence fixture deferred to exact Pysam Apptainer runtime validator'
PYTHONNOUSERSITE=1 python "$PROJECT_DIR/test_expression_go_analysis.py" && pass 'Dependency-free expression and variant-set GO fixture' || fail 'Expression and variant-set GO fixture'
PYTHONNOUSERSITE=1 python "$PROJECT_DIR/test_progression_biology.py" && pass 'Dependency-free scalable progression biology fixture' || fail 'Scalable progression biology fixture'
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
