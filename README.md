# PGTK: RNA-seq Proteogenomics and Comparative Variant Evidence

PGTK is a Nextflow DSL2 workflow for exploratory RNA-seq proteogenomics. It processes paired-end RNA-seq from local SRA archives, performs RNA-aware small-variant calling, detects fusions and novel splice-derived transcripts, generates sample-specific protein search FASTAs, quantifies gene expression, performs GO analysis, supports patient-aware baseline subtraction, and can integrate external VCF and MaxQuant evidence.

PGTK is intended for research use. RNA-observed variants may represent germline, clonal, progression-associated, RNA-editing, alignment, expression-dependent, or technical events. They must not be described as clinically validated somatic mutations without independent DNA evidence.

## Current production status

The current production implementation has been validated on Saga with Nextflow 26.04.6. Job `19246213` completed normally on August 18, 2026 with exit status 0.

The current finding-exploration architecture is database-free:

- 157,482 findings retained
- 157,482 partition records
- 224 biological partitions
- 0 discarded findings
- 0 database files
- no SQLite metadata or alignment database
- no arbitrary finding limit
- no monolithic all-findings IGV HTML report
- IGV reports generated only when a finding is selected

## Workflow overview

```text
Local SRA archives
  -> paired FASTQ
  -> raw FastQC
  -> Trim Galore
  -> trimmed FastQC
  -> STAR two-pass alignment
       -> Arriba fusion branch
       -> coordinate-sorted and indexed BAM
            -> featureCounts expression branch
            -> StringTie and gffcompare splice branch
            -> MarkDuplicates
            -> SplitNCigarReads
            -> HaplotypeCaller, 24 shards per sample
            -> shard validation
            -> GatherVcfs and indexing
            -> GenotypeGVCFs
            -> raw VCF publication
            -> normalization and hard filtering
            -> PASS VCF publication
            -> VEP annotation
            -> RNA validation
            -> codon and read-provenance validation
            -> sample-specific variant protein FASTA
            -> patient-aware baseline comparison
            -> complete finding manifest and IGV evidence bundle
            -> database-free searchable finding explorer
```

Custom protein search FASTAs are built as:

```text
variant proteins
+ fusion proteins
+ splice-derived proteins
-> exact amino-acid sequence deduplication
-> <sample>.exploratory_proteogenomics.fasta
```

Canonical proteins are not embedded in these custom FASTAs. Add the canonical human proteome separately during the MaxQuant search.

## Core capabilities

The workflow includes:

- Nextflow DSL2 execution with Saga Slurm routing
- local SRA ingestion with no compute-node internet requirement
- paired FASTQ generation
- raw and trimmed FastQC
- Trim Galore preprocessing
- STAR two-pass RNA alignment
- coordinate-sorted and indexed STAR BAMs
- Arriba fusion calling and RNA validation
- StringTie transcript assembly
- gffcompare novelty classification
- splice-transcript RNA validation
- TransDecoder-derived splice protein prediction
- GATK SplitNCigarReads
- 24-way HaplotypeCaller scattering per sample
- shard validation, GatherVcfs, indexing, and GenotypeGVCFs
- publication of GVCF, raw, normalized, SNP, indel, filtered, PASS, VEP, and RNA-validated VCF stages
- codon-level validation
- supporting-read provenance validation
- per-sample non-subtracted protein FASTAs
- patient-aware baseline subtraction as a separate reporting branch
- shared, baseline-only, and non-baseline-only progression reports
- optional external Sarek or other caller comparison
- optional MaxQuant evidence interpretation
- global IGV BED, BEDPE, BAM, batch, and session outputs
- single-pass consolidated strict finding review
- database-free searchable finding explorer
- on-demand standalone IGV.js reports
- featureCounts gene quantification
- merged raw-count, CPM, and TPM expression matrix
- per-sample expression GO over-representation analysis
- progression-versus-baseline ranked expression GO
- common, complete-sample, and exclusive progression variant-set GO
- pairwise progression GO contrasts
- claim audit with strict evidence-layer separation
- MultiQC, comparative, failure, resource, and provenance reports

## Requirements

The production Saga workflow expects:

- Nextflow 26.04.6 or a compatible version
- Java 21
- Slurm
- Apptainer
- local container images
- local SRA archives
- local reference assets
- Python 3 for the wrapper and reporting utilities

Run `download_sra.sh` and `download_assets.sh` directly on a login node, not through Slurm as compute node fire-wall seems unsurpassable via nextflow

## Samplesheet

Only `sample` and `srr` are required.

Minimal form:

```csv
sample,srr
TK12,SRR31089074
TK13,SRR31089073
TK14,SRR31089072
```

Longitudinal form:

```csv
sample,srr,TK,Group,baseline
TK12,SRR31089074,patient1,resistant,true
TK13,SRR31089073,patient1,sensitive,false
TK14,SRR31089072,patient1,sensitive,false
```

Column semantics:

- `sample`: unique sample identifier
- `srr`: local SRA run accession
- `TK`: patient or subject identifier
- `Group`: biological metadata
- `baseline`: comparison reference within each `TK`

Defaults:

```text
TK       = sample
Group    = sample
baseline = false
```

Per-sample FASTAs are always generated independently of baseline subtraction. A non-baseline sample is compared only when its subject has exactly one baseline. Missing or multiple baselines are explicitly reported.

## Initial setup on Saga

```bash
cd /cluster/projects/nn9036k/scrbkup/pgtk/checkport

bash download_sra.sh
bash download_assets.sh

bash validate_pipeline_commands.sh \
  --project-dir "$PWD" \
  --nextflow "$HOME/bin/nextflow"
```

Do not submit either download script through Slurm.

## Validate before submission

Run these checks from the project directory:

```bash
cd /cluster/projects/nn9036k/scrbkup/pgtk/checkport && \
python3 test_finding_igv_reviews.py && \
python3 test_resource_configuration.py && \
python3 -m py_compile build_finding_explorer.py build_finding_igv_reviews.py && \
bash -n scratch.slurm serve_finding_explorer.sh
```

Expected output includes:

```text
PASS: consolidated review classifier, deduplication, and flat output layout
PASS: 72 processes, exact rendered explorer shell syntax, database-free explorer, zero finding limit
```

## Submit the production Saga job

The wrapper accepts Nextflow parameters after `--`. It uses resume-compatible execution and preserves cached processes.

```bash
cd /cluster/projects/nn9036k/scrbkup/pgtk/checkport && JOB_ID=$(sbatch --parsable --account=nn9036k --partition=normal --export=ALL,PGTK_ACCOUNT=nn9036k,PGTK_NORMAL_PARTITION=normal,PGTK_BIGMEM_PARTITION=bigmem,PGTK_PROJECT_DIR=/cluster/projects/nn9036k/scrbkup/pgtk/checkport,PGTK_WORK_DIR=/cluster/work/users/ash022/work,PGTK_CONTAINER_CACHE=/cluster/projects/nn9036k/scrbkup/pgtk/singularity_cache,PGTK_PYSAM_IMAGE=/cluster/projects/nn9036k/scrbkup/pgtk/singularity_cache/quay.io-biocontainers-pysam-0.24.0--py312hf5ad864_1.img,PGTK_SRA_DIR=/cluster/projects/nn9036k/scrbkup/pgtk/sra_cache,PGTK_REFERENCE_DOWNLOADS=/cluster/projects/nn9036k/scrbkup/pgtk/reference_downloads,PGTK_ENSEMBL_PEP=/cluster/projects/nn9036k/scrbkup/pgtk/reference_downloads/Homo_sapiens.GRCh38.pep.all.fa.gz,PGTK_SAMPLESHEET=/cluster/projects/nn9036k/scrbkup/pgtk/checkport/samples.csv,PGTK_RESULTS_DIR=/cluster/projects/nn9036k/scrbkup/pgtk/checkport/results,PGTK_NEXTFLOW=$HOME/bin/nextflow,PGTK_PYTHON=$(command -v python3),PGTK_JAVA_MODULE=Java/21,PGTK_SLURM_LOG_TEMPLATE=/cluster/projects/nn9036k/scrbkup/pgtk/checkport/pgtk-wrapper-{job_id}.log scratch.slurm -- --finding_priority_mode all --igv_report_classes rna_variant,progression_variant,fusion,splice_junction --igv_report_max_reads 100) && echo "JOB_ID=$JOB_ID"
```

## Track the wrapper and pipeline

Follow the wrapper log:

```bash
cd /cluster/projects/nn9036k/scrbkup/pgtk/checkport && \
until [ -f "pgtk-wrapper-${JOB_ID}.log" ]; do sleep 2; done; \
tail -f "pgtk-wrapper-${JOB_ID}.log"
```

Watch the Slurm queue and recent wrapper output every 30 seconds:

```bash
watch -n 30 'date; squeue -u ash022 -o "%.18i %.12P %.35j %.10T %.10M %.10l %R"; echo; tail -15 /cluster/projects/nn9036k/scrbkup/pgtk/checkport/pgtk-wrapper-'"${JOB_ID}"'.log'
```

Inspect job accounting after completion:

```bash
sacct -j "$JOB_ID" --format=JobID,JobName,State,ExitCode,Elapsed,AllocCPUS,ReqMem,MaxRSS
```

Inspect Nextflow task status:

```bash
cd /cluster/projects/nn9036k/scrbkup/pgtk/checkport
awk -F '\t' 'NR==1 || $4 != "COMPLETED"' "results/pipeline_trace-${JOB_ID}.tsv"
```

A successful run ends with wrapper exit status 0 and normally includes:

```text
Cached process > BUILD_FINDING_IGV_REVIEWS
Submitted process > BUILD_FINDING_EXPLORER
Generated database-free explorer for 157482 findings; discarded 0
Execution complete
Job exited normally.
```

## Resume behavior

Keep these items when resume compatibility is required:

```text
/cluster/work/users/ash022/work/
.nextflow/
.nextflow.log*
results/
```

Do not remove work directories while a resume run may be needed. The wrapper invokes Nextflow with `-resume`.

## Local execution

`run.sh` is the repository local runner. It uses the same local SRA files, references, and Apptainer images where configured.

```bash
bash run.sh
```

Paths can be overridden with environment variables such as:

```text
SRA_DIR
REFERENCE_DOWNLOADS
CONTAINER_CACHE
WORK_DIR
RESULTS_DIR
```

## Variant stages

```text
results/gvcf/<sample>.g.vcf.gz
results/vcf_raw/<sample>.raw.vcf.gz
results/vcf_filtered/<sample>.filtered.vcf.gz
results/vcf_pass/<sample>.pass.vcf.gz
results/vep/<sample>.vep.vcf.gz
results/rna_validation/variants/<sample>.rna.validated.vcf.gz
```

Raw GenotypeGVCFs output and its index are published before hard filtering. This preserves direct comparisons with Sarek and other caller stages.

RNA-supported variants are not automatically DNA-confirmed mutations. Baseline absence in RNA does not establish biological or DNA-level absence.

## Longitudinal progression branch

For each valid progression-versus-baseline pair:

```text
results/progression_vcf/<sample>.nonbaseline_only.vep.vcf.gz
results/progression_vcf/<sample>.baseline_only.vep.vcf.gz
results/progression_vcf/<sample>.shared_with_baseline.vep.vcf.gz
results/progression_vcf/<sample>.subtraction.summary.tsv
```

Progression biology outputs include:

```text
results/progression_biology/progression_biology.progression_alleles.tsv
results/progression_biology/progression_biology.progression_genes.tsv
results/progression_biology/progression_biology.go_enrichment.tsv
results/progression_biology/progression_biology.pairwise_go_contrasts.tsv
results/progression_biology/sets/progression_variant_sets.variant_set_go.tsv
results/progression_biology/sets/progression_variant_sets.summary.tsv
```

Complete tested GO-term tables remain in TSV outputs. MultiQC displays summaries.

## Gene expression and GO

Default featureCounts settings:

```text
--gene_count_feature_type exon
--gene_count_id_attribute gene_id
--gene_count_symbol_attribute gene_name
--gene_count_biotypes all
--gene_count_strandedness 0
--gene_count_min_mapq 10
--gene_count_min_overlap 1
--gene_count_count_read_pairs true
--gene_count_require_both_ends false
--gene_count_exclude_chimeric true
--gene_count_primary_only true
--gene_count_allow_multi_overlap false
--gene_count_count_multimapping false
```

Expression outputs:

```text
results/expression/per_sample/<sample>.gene_counts.tsv
results/expression/per_sample/<sample>.gene_counts.tsv.summary
results/expression/gene_expression.gene_expression.tsv
results/expression/gene_expression.summary.tsv
results/expression/go/expression_go.expression_ora.tsv
results/expression/go/expression_go.ranked_go.tsv
results/expression/go/expression_go.summary.tsv
```

CPM uses the assigned-library total. TPM uses non-overlapping exon-union lengths and is normalized within each sample.

Expression-ranked GO uses:

```text
log2((progression TPM + 0.5) / (baseline TPM + 0.5))
```

Expression-ranked GO and progression-variant GO are separate analyses:

```text
Expression-ranked GO
  tests coordinated transcript direction

Progression-variant GO
  tests concentration of variant-bearing genes
```

## Global IGV evidence bundle

```text
results/igv/all_evidence/pgtk_igv.events.tsv
results/igv/all_evidence/pgtk_igv.events.bed
results/igv/all_evidence/pgtk_igv.events.bedpe
results/igv/all_evidence/pgtk_igv.sample_manifest.tsv
results/igv/all_evidence/pgtk_igv.igv.batch.txt
results/igv/all_evidence/pgtk_igv.igv.session.xml
results/igv/all_evidence/pgtk_igv.<sample>.events.bam
results/igv/all_evidence/pgtk_igv.<sample>.events.bam.bai
```

## Consolidated strict finding review

`BUILD_FINDING_IGV_REVIEWS` scans the staged BAM inputs once and creates a flat consolidated evidence bundle. It does not create a directory per finding and does not use SQLite.

Principal outputs under `results/igv/findings/finding_reviews/` include:

```text
findings_manifest.tsv
priority_findings.tsv
priority_findings.bed
support_labels.bed
event_consolidation.tsv
event_regions.tsv
bam_manifest.tsv
consolidation_summary.txt
README.txt
review.igv.batch.txt
igv.session.xml
<sample>.exact_alt_unique.bam
<sample>.exact_alt_unique.bam.bai
<sample>.exact_alt_display.bam
<sample>.exact_alt_display.bam.bai
<sample>.reference_display.bam
<sample>.reference_display.bam.bai
<sample>.event_display.bam
<sample>.event_display.bam.bai
```

The implementation uses one BAM scan, interval lookup, and direct BAM output. It does not create `alignment_store.sqlite` or any other database.

## Database-free finding explorer

`BUILD_FINDING_EXPLORER` retains every finding and partitions the searchable metadata by:

```text
sample
+ evidence class
+ chromosome
```

Outputs:

```text
results/igv/findings/finding_explorer/index.html
results/igv/findings/finding_explorer/server.py
results/igv/findings/finding_explorer/serve_explorer.sh
results/igv/findings/finding_explorer/explorer_config.json
results/igv/findings/finding_explorer/partition_manifest.tsv
results/igv/findings/finding_explorer/coverage_summary.txt
results/igv/findings/finding_explorer/partitions/*.jsonl.gz
```

The explorer:

- contains every finding
- has no arbitrary finding threshold
- stores metadata once in compressed JSON Lines files
- loads metadata into memory when the local server starts
- filters by sample, evidence class, impact, chromosome, gene, event, transcript, consequence, and protein change
- generates one standalone IGV.js report only when the user selects a finding
- caches generated reports as ordinary HTML files
- creates no `.sqlite`, `.sqlite3`, or `.db` files

Verify completeness:

```bash
cat results/igv/findings/finding_explorer/coverage_summary.txt
```

Expected for the current production dataset:

```text
Findings: 157482
Partition records: 157482
Biological partitions: 224
Findings discarded: 0
Database files: 0
```

Start the explorer on Saga:

```bash
cd /cluster/projects/nn9036k/scrbkup/pgtk/checkport/results/igv/findings/finding_explorer
export PGTK_IGV_REPORTS_IMAGE=/cluster/projects/nn9036k/scrbkup/pgtk/singularity_cache/quay.io-biocontainers-igv-reports-1.16.0--pyh7e72e81_0.img
./serve_explorer.sh "$PWD" 8765
```

From the workstation, create an SSH tunnel:

```bash
ssh -L 8765:127.0.0.1:8765 ash022@login.saga.sigma2.no
```

Open:

```text
http://127.0.0.1:8765
```

## Optional external-caller comparison

```bash
sbatch scratch.slurm -- \
  --run_external_vcf_comparison true \
  --external_vcf_dir /cluster/projects/nn9036k/scrbkup/pgtk/checkport/sarek \
  --external_vcf_suffix .haplotypecaller.vcf.gz
```

Files are matched recursively using the samplesheet `srr` field. Comparisons are written under:

```text
results/comparison/external_vcf/
```

Raw, PASS, and RNA-validated PGTK stages are compared separately.

## Optional MaxQuant integration

PGTK does not execute MaxQuant. It interprets an existing MaxQuant search performed with the exact custom FASTAs, canonical proteome, and contaminants database.

Custom FASTAs:

```text
results/combined_fasta/TK12.exploratory_proteogenomics.fasta
results/combined_fasta/TK13.exploratory_proteogenomics.fasta
results/combined_fasta/TK14.exploratory_proteogenomics.fasta
```

Required MaxQuant text files:

```text
peptides.txt
evidence.txt
msms.txt
proteinGroups.txt
mqpar.xml
```

Enable integration:

```bash
sbatch scratch.slurm -- \
  --run_proteogenomic_validation true \
  --maxquant_txt /cluster/projects/nn9036k/scrbkup/pgtk/checkport/txtMQMBR \
  --maxquant_contaminants "$HOME/scripts/MaxQuant_v2.8.1.0/bin/conf/contaminants.fasta"
```

If raw names cannot be mapped unambiguously:

```bash
--maxquant_raw_map /cluster/projects/nn9036k/scrbkup/pgtk/checkport/maxquant_raw_file_map.tsv
```

Optional provenance overrides:

```bash
--maxquant_mqpar /absolute/path/to/mqpar.xml \
--maxquant_canonical_fasta /absolute/path/to/human_reviewed_isoforms.fasta
```

Principal outputs are written under:

```text
results/proteogenomics_validation/
```

Strict event confirmation requires sample-matched direct MS/MS evidence, consistency with the searched FASTAs, altered-residue or junction coverage, and absence from canonical reference proteins. File presence alone is not biological confirmation.

## Claim audit

Run after pipeline completion using the wrapper job ID:

```bash
cd /cluster/projects/nn9036k/scrbkup/pgtk/checkport
python3 audit_pgtk_claims.py \
  --project-dir "$PWD" \
  --job-id "$JOB_ID" \
  --output-prefix "results/pgtk_claim_audit_${JOB_ID}"
```

The audit keeps transcript signatures, ranked GO, progression-variant GO, external-caller context, and proteomic evidence as separate layers.

The predefined hypotheses are:

```text
H1 proteostasis and translation
H2 metabolic rewiring
H3 surface, glycan, extracellular-matrix and adhesion remodeling
H4 MYC, IRF4, DNA repair and progression
H5 PCNA stress and ATX-101-associated biomarkers
```

## MultiQC

```text
results/multiqc/multiqc_report.html
results/multiqc/multiqc_report_data/
```

MultiQC includes standard QC and custom sections for RNA findings, validation failures, codon evidence, read provenance, progression biology, comparative evidence, expression GO, proteogenomics evidence when enabled, and validation semantics.

## Resource, retry, and failure model

Resource-related failures with exit codes 137, 140, or 143 permit three attempts:

```text
Attempt 1: 1x resources
Attempt 2: 2x resources
Attempt 3: 4x resources
```

Limits:

```text
Maximum CPUs: 32
Maximum memory: 512 GB
```

Partition routing:

```text
normal: <= 20 CPUs and <= 160 GB
bigmem: > 20 CPUs or > 160 GB
```

Failure records:

```text
results/failure_logs/<job-id>/
results/failure_logs/failure_history.tsv
results/failure_logs/run_history.tsv
```

Execution reports:

```text
results/pipeline_trace-<job-id>.tsv
results/pipeline_report-<job-id>.html
results/pipeline_timeline-<job-id>.html
results/pipeline_dag-<job-id>.html
results/resource_usage-<job-id>.summary.tsv
results/resource_usage-<job-id>.warnings.tsv
results/resource_usage-<job-id>.report.md
```

## Repository file guide

### Core workflow and execution

- `main.nf`: complete Nextflow DSL2 workflow and process wiring.
- `nextflow.config`: Slurm executor, profiles, resource selectors, retry rules, container settings, and process-specific overrides.
- `scratch.slurm`: production Saga wrapper that validates runtime assets, constructs the Nextflow command, enables resume, records resource summaries, and collects failures.
- `run.sh`: local execution wrapper.
- `samples.csv`: active pipeline samplesheet.
- `sample_lane0.csv`: lane-specific sample-sheet example or retained input.
- `sample_laneNA.csv`: sample-sheet input for lane splitting.
- `samples_lane_split.csv`: lane-split sample-sheet output or example.
- `wsl.config`: WSL-specific configuration retained for local development.
- `LICENSE`: project license.
- `README.md`: this documentation.

### Input and asset acquisition

- `download_sra.sh`: downloads and validates SRA archives on a login node.
- `download_assets.sh`: downloads reference data and container assets on a login node.
- `downloadFastq.sh`: legacy or auxiliary FASTQ download helper.
- `ena_fastq_urls.txt`: ENA FASTQ URL list used by auxiliary download workflows.
- `generateSampleSheet.sh`: generates a samplesheet from available inputs.
- `laneSplit.py`: splits paired FASTQ files by lane and writes an updated samplesheet.
- `coll.sh`: collection helper used to package selected pipeline files and diagnostics.

### Core RNA, variant, and validation code

- `validate_rna_events.py`: validates RNA-observed variants, fusions, or splice events and writes structured evidence outputs.
- `validate_haplotype_shards.py`: verifies expected HaplotypeCaller shards before gather.
- `summarize_variant_stages.py`: summarizes raw, normalized, filtered, PASS, annotated, and validated stages.
- `validate_variant_codons.py`: performs codon and translated-residue validation.
- `validate_variant_read_provenance.py`: records read-level provenance for variant support.
- `merge_variant_validation.py`: merges variant validation outputs.
- `analyze_codon_mismatches.py`: investigates codon or translation mismatches.
- `build_integrated_variant_evidence.py`: combines variant, RNA, codon, read, and optional proteomic evidence.
- `annotate_vep.sh`: VEP annotation helper.
- `annotate_variant_peptides.py`: annotates variant-derived peptide relationships.

### Fusion, splice, and protein FASTA code

- `analyze_chimeric_splice_peptides.py`: analyzes peptide support for chimeric and splice-derived sequences.
- `validate_splice_junction_peptides.py`: validates peptide coverage across splice junctions.
- `map_peptides_to_fasta.py`: maps peptide sequences to custom and reference FASTAs.
- `validate_proteogenomic_reads.py`: performs read-level validation for proteogenomic events and produces compact IGV outputs.
- `proteogenomics_evidence_report.py`: creates structured proteogenomics evidence reports.

### Expression and biological interpretation

- `expression_go_analysis.py`: merges gene counts and performs expression ORA, ranked GO, sample-level GO, and progression variant-set GO.
- `prepare_go_annotations.py`: prepares GO mapping resources from configured inputs.
- `analyze_progression_biology.py`: analyzes progression-specific alleles, genes, and GO terms.
- `compare_progression_pair.py`: compares two progression samples or time points.
- `merge_progression_biology.py`: merges progression biology results.
- `build_comparative_advantage_report.py`: builds the comparative biological advantage report.
- `build_complete_report.py`: generates the complete findings report.
- `audit_pgtk_claims.py`: audits scientific claims against explicit evidence layers.
- `audit_environment_hardcoding.py`: checks active code for environment-specific hard-coded values.
- `audit_pgtk_claims.py`: scientific interpretation and claim-separation audit.

### IGV and finding exploration

- `build_igv_evidence_bundle.py`: creates the global RNA and progression IGV evidence bundle.
- `build_finding_igv_reviews.py`: performs the single-pass consolidated strict finding review and writes flat BAM and manifest outputs.
- `build_finding_explorer.py`: builds the database-free compressed finding partitions and local explorer assets.
- `serve_finding_explorer.sh`: starts the explorer inside the offline IGV Reports Apptainer image.
- `check_bedcov_progress.sh`: checks progress for bedcov-based operations.
- `run_bedcov_parallel.sh`: auxiliary parallel bedcov runner.

### External comparison and cohort utilities

- `compare_external_vcf.py`: compares PGTK call stages with external VCFs.
- `external_comparison.none.tsv`: empty placeholder for disabled external comparison.
- `group_compare_two_cohorts.py`: compares two configured cohorts.
- `group_exclusive_mutations.py`: identifies cohort-exclusive mutation or gene sets.
- `group_exclusive_mutations.sh`: shell wrapper for the cohort-exclusive analysis.
- `cross_chromosome.py`: summarizes cross-chromosome fusion events.
- `annotate_and_pattern.py`: annotates exclusive genes and explores chromosome, arm, biotype, and recurrence patterns.
- `mannwhitney_perchrom.py`: performs per-chromosome Mann-Whitney analysis from idxstats-derived coverage.
- `plot_coverage_ideogram.py`: produces chromosome-level coverage ideograms.
- `plot_coverage_per_chrom.py`: plots per-chromosome coverage.
- `plot_windowed_cn.py`: plots windowed copy-number-like coverage patterns.

### Failure, resource, and trace tools

- `collect_pipeline_failures.py`: records process failures and maintains failure history.
- `analyze_pipeline_trace.py`: summarizes Nextflow trace resource utilization and warnings.
- `RESOURCE_RETRY_MATRIX.tsv`: documented process retry and resource behavior.
- `trace-20260815-45198746.txt`: retained trace fixture or historical trace header.
- `pgtk_region_fix_inputs_19217675.sha256`: checksum record for a historical validation package.

### Validation and regression tests

- `validate_pipeline_commands.sh`: validates declarations, commands, paths, containers, and configuration.
- `validate_runtime_inputs.py`: validates runtime references, containers, SRA inputs, samplesheet, Python, featureCounts, IGV Reports, and Pysam contracts.
- `validate_full_update.slurm`: Slurm validation job for the complete update.
- `validate_region_optimized_igv_fixture.sh`: validates the region-optimized IGV fixture.
- `validate_single_pass_igv_fixture.sh`: validates the single-pass IGV review fixture.
- `test_finding_igv_reviews.py`: tests finding classification, deduplication, and flat output layout.
- `test_igv_evidence_bundle.py`: tests the global IGV evidence bundle.
- `test_resource_configuration.py`: validates process count, parser-safe rendered explorer shell, database-free design, and zero finding limit.
- `test_expression_go_analysis.py`: regression tests for expression and progression GO analysis.
- `test_progression_biology.py`: regression tests for progression biology outputs.
- `probe_pipeline_cli.sh`: probes pipeline CLI parsing and parameter behavior.

### MaxQuant and XML inputs

- `mqpar.xml`: MaxQuant parameter XML supplied or retained for provenance checking.
- `callparam.xml`: auxiliary call-parameter XML.
- `extractionSummary.xml`: extraction summary XML.
- `maxquant_raw_file_map.none.tsv`: empty placeholder for optional MaxQuant raw-file mapping.

### Documentation and historical implementation records

- `PROJECT_SUMMARY.md`: project-level summary.
- `PUBLICATION_MAPPING.md`: publication and dataset mapping.
- `RUN_AND_REPORT.md`: operational run and reporting notes.
- `PIPELINE_UPDATE_20260807.md`: August 7 pipeline update notes.
- `FULL_PIPELINE_UPDATE_20260814.md`: August 14 complete update summary.
- `FULL_UPDATE_NOTES.md`: additional full-update notes.
- `PORTABLE_RESOURCE_IGV_UPDATE_20260815.md`: portable resource and IGV update notes.
- `EFFICIENCY_REVIEW_20260816.md`: file-quota and efficiency review.
- `SINGLE_PASS_IGV_REDESIGN_20260817.md`: single-pass consolidated IGV redesign description.
- `REGION_OPTIMIZED_IGV_UPDATE_20260817.md`: region-optimized IGV notes.
- `PORTABLE_POSTPROCESS_FIX_20260818.md`: portable post-processing correction.
- `FINDING_EXPLORER_DATABASE_FREE_20260818.md`: database-free explorer architecture.
- `GO_OPTIMIZATION_VALIDATION.md`: GO optimization validation.
- `GO_COMPLETE_SAMPLE_AND_FDR_VALIDATION.md`: complete-sample and FDR GO validation.
- `RANKED_GO_VALIDATION.md`: ranked GO validation.
- `PYTHON_RUNTIME_FIX_VALIDATION.md`: Python runtime validation notes.
- `PIPELINE_VALIDATION_SEMANTICS.md`: validation semantics and interpretation boundaries.
- `NEXTFLOW_CONFIG_PARSER_FIX.md`: Nextflow configuration parser correction.
- `NEXTFLOW_STRICT_CONFIG_FIX.md`: strict configuration correction.
- `PHASE2_CLI_PARAMETERS.md`: Phase 2 CLI parameter documentation.
- `HISTORY_REGRESSION_CHECKLIST.md`: regression checklist for preserving prior fixes.
- `RESOURCE_CALIBRATION.md`: resource calibration guidance.
- `scratch.md`: historical or working Slurm notes.
- `pipeline_command_validation_20260815_101734.txt`: retained command-validation output.
- `pipeline_command_validation_20260815_123244.txt`: retained command-validation output.

Historical notes document implementation decisions but do not override `main.nf`, `nextflow.config`, `scratch.slurm`, or the current regression tests.

## Publication and public data

Biological context:

**Multiple Myeloma Cells with Increased Proteasomal and ER Stress Are Hypersensitive to ATX-101, an Experimental Peptide Drug Targeting PCNA**

```text
Journal: Cancers
Year: 2024
Volume: 16
Issue: 23
Article: 3963
DOI: https://doi.org/10.3390/cancers16233963
RNA-seq BioProject: PRJNA1176350
Proteomics: PRIDE PXD033531 and PXD033510
```

Cite the publication when using its data or derived results.

## Interpretation rules

- RNA-supported variants are not DNA-confirmed mutations.
- Non-baseline-only RNA observations are not proven newly acquired variants.
- RNA-supported fusions are not confirmed genomic rearrangements.
- Read-supported splice transcripts do not prove complete isoform structure.
- FASTA entries are search candidates, not peptide confirmation.
- Ranked GO and targeted signatures from the same expression matrix are not independent validation.
- GO terms are overlapping ontology categories, not independent pathways.
- External-caller overlap is general caller reproducibility unless gene-specific support is shown.
- MaxQuant output presence is not event confirmation.
- Match Between Runs is not equivalent to direct MS/MS.
- Sample-matched event-spanning peptides are stronger evidence than canonical protein-group detection.
- Transcript abundance does not establish protein abundance, pathway activity, metabolic flux, drug sensitivity, or causation.
- A single longitudinal subject does not support population-level differential-expression inference.
- Similarity to an sPCL signature does not establish clinical plasma cell leukemia.

## Limitations

- RNA detectability depends on expression and coverage.
- Unexpressed genomic regions cannot be evaluated reliably.
- Alignment, reference, filtering, and soft-clipping choices affect RNA callsets.
- External comparisons require compatible references and normalized alleles.
- Baseline subtraction reports absence from the observed baseline RNA callset, not proven biological absence.
- Custom FASTAs enlarge the search space and require appropriate FDR control.
- MaxQuant interpretation depends on the exact searched FASTAs, contaminants database, raw-file mapping, and search configuration.
- Proteomic and functional confirmation must remain distinct from RNA-level evidence.

## License

See [LICENSE](LICENSE).

## Contact

Open a GitHub issue or contact `animesh@fuzzylife.org`.
