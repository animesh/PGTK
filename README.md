# PGTK: RNA-seq Proteogenomics and Comparative Variant Evidence

PGTK is a Nextflow DSL2 workflow for exploratory RNA-seq proteogenomics. It processes paired-end RNA-seq from local SRA archives, performs RNA-aware small-variant calling, detects fusions and novel splice-derived transcripts, creates sample-specific protein search FASTAs, quantifies gene expression, performs several distinct GO analyses, supports patient-aware baseline comparison, and can integrate external VCF and MaxQuant evidence.

PGTK is intended for research use. An RNA-observed event may represent germline variation, a clonal or progression-associated event, RNA editing, transcript-specific expression, alignment ambiguity, local assembly behavior, or technical noise. RNA evidence must not be described as a clinically validated somatic mutation without independent DNA evidence.

## Production status

The current implementation was validated on Saga with Nextflow 26.04.6.

Latest validated production run:

```text
Job: 19378463
Date: 2026-08-22
Exit status: 0
Validated processes: 73
Static and runtime contract checks: 337 PASS, 0 FAIL
Finding explorer records: 157,632
Database files: 0
```

The production finding explorer is database-free and has no arbitrary finding limit. It provides two modes:

- Direct-file mode opens `index.html` through `file:///` for local searching, filtering, read counts, and validation status.
- Server mode generates a standalone IGV.js alignment report only when a finding is selected.

## Evidence model

PGTK keeps candidate generation and strict validation separate.

```text
Upstream RNA candidate
  -> VCF filtering and VEP annotation
  -> RNA-event selection
  -> strict read-by-read review
  -> ALT-supported, mixed, reference-only, or non-callable status
  -> optional codon, provenance, external-VCF, and proteomic evidence
```

A candidate remains visible for audit even when strict read validation does not support the exact allele.

Strict read-validation statuses are:

- `ALT_SUPPORTED`: one or more exact ALT-supporting reads and no clean reference reads.
- `MIXED_ALT_AND_REFERENCE`: exact ALT and clean reference reads are both present.
- `REFERENCE_ONLY`: callable reads support the reference allele and no exact ALT read is present.
- `NO_CALLABLE_READS`: alignments overlap the candidate region, but none can determine the exact ALT-versus-reference allele.
- `NO_OVERLAPPING_READS`: no reviewed alignment overlaps the candidate locus.

For an insertion such as `C>CA`, a colored `A` mismatch at another reference position is not insertion evidence. Exact support requires the expected inserted sequence and an insertion CIGAR operation at the correct anchor.

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
            -> codon and supporting-read provenance validation
            -> sample-specific variant protein FASTA
            -> patient-aware baseline comparison
            -> progression and GO analyses
            -> complete finding manifest
            -> consolidated read-evidence BAMs
            -> database-free finding explorer
            -> on-demand per-finding IGV.js reports
            -> compact final MultiQC report
```

Custom protein search FASTAs contain:

```text
variant proteins
+ fusion proteins
+ splice-derived proteins
-> exact amino-acid sequence deduplication
-> <sample>.exploratory_proteogenomics.fasta
```

Canonical proteins are not embedded in the custom FASTAs. Add the canonical human proteome separately during the MaxQuant search.

## Main capabilities

- Nextflow DSL2 execution with Slurm and Apptainer.
- Local SRA ingestion with no compute-node internet dependency.
- Paired FASTQ generation, FastQC, trimming, and STAR two-pass alignment.
- Arriba fusion calling and RNA validation.
- StringTie assembly, gffcompare novelty classification, and splice-derived protein prediction.
- GATK RNA-aware small-variant calling with 24 HaplotypeCaller shards per sample.
- Publication of GVCF, raw, normalized, filtered, PASS, VEP, and RNA-validated VCF stages.
- Codon-level validation and supporting-read provenance.
- Per-sample non-subtracted custom protein FASTAs.
- Patient-aware baseline comparison as a separate reporting branch.
- Shared, baseline-only, and non-baseline-only progression reports.
- Optional external Sarek or other caller comparison.
- Optional interpretation of existing MaxQuant results.
- Global IGV BED, BEDPE, BAM, batch, and session outputs.
- Single-pass consolidated strict finding review.
- Database-free searchable explorer with explicit candidate-validation status.
- On-demand, sample-specific IGV.js reports containing all-overlapping, exact-ALT, and clean-reference tracks.
- featureCounts gene quantification with raw-count, CPM, and TPM matrices.
- Per-sample expression GO, ranked progression-versus-baseline expression GO, progression variant-set GO, pairwise GO contrasts, and protein-altering variant GO.
- Compact MultiQC, comparative, failure, resource, provenance, and interpretation reports.

## Requirements

The production Saga workflow expects:

- Nextflow 26.04.6 or a compatible version.
- Java 21.
- Slurm.
- Apptainer.
- Local container images.
- Local SRA archives.
- Local reference assets.
- Python 3 for wrappers and reporting utilities.

Run internet-dependent scripts directly on a login node:

```bash
bash download_sra.sh
bash download_assets.sh
```

Do not submit download scripts through Slurm. Saga compute nodes do not provide the required external network access.

## Samplesheet

Only `sample` and `srr` are required.

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

Column meanings:

- `sample`: unique sample identifier.
- `srr`: local SRA run accession.
- `TK`: patient or subject identifier.
- `Group`: biological metadata.
- `baseline`: comparison reference within each `TK`.

Defaults:

```text
TK       = sample
Group    = sample
baseline = false
```

Per-sample FASTAs are always generated independently of baseline subtraction. A non-baseline sample is compared only when the corresponding subject has exactly one baseline. Missing and multiple baselines are reported explicitly.

## Validate before submission

Run from the project directory:

```bash
python3 -m py_compile *.py
python3 test_full_pipeline_contract.py
python3 test_reporting_redesign.py
python3 test_resource_configuration.py
python3 audit_environment_hardcoding.py .

bash validate_pipeline_commands.sh \
  --project-dir "$PWD" \
  --nextflow "$HOME/bin/nextflow"
```

Expected final result:

```text
PASS: 337
FAIL: 0
RESULT: PASSED
```

## Submit on Saga

`scratch.slurm` supplies `-resume`. Do not pass a second `-resume` argument.

```bash
cd /cluster/projects/nn9036k/scrbkup/pgtk/checkport && MAXQUANT_DIR=/cluster/home/ash022/scripts/pgtk/ftp.pride.ebi.ac.uk/pride/data/archive/2024/11/PXD033510/combined/txt && SAREK_DIR=/cluster/home/ash022/scripts/pgtk/sarek/variant_calling/haplotypecaller && JOB_ID=$(sbatch --parsable --account=nn9036k --partition=normal --export=ALL,PGTK_ACCOUNT=nn9036k,PGTK_NORMAL_PARTITION=normal,PGTK_BIGMEM_PARTITION=bigmem,PGTK_PROJECT_DIR="$PWD",PGTK_WORK_DIR=/cluster/work/users/ash022/work,PGTK_CONTAINER_CACHE=/cluster/projects/nn9036k/scrbkup/pgtk/singularity_cache,PGTK_PYSAM_IMAGE=/cluster/projects/nn9036k/scrbkup/pgtk/singularity_cache/quay.io-biocontainers-pysam-0.24.0--py312hf5ad864_1.img,PGTK_SRA_DIR=/cluster/projects/nn9036k/scrbkup/pgtk/sra_cache,PGTK_REFERENCE_DOWNLOADS=/cluster/projects/nn9036k/scrbkup/pgtk/reference_downloads,PGTK_ENSEMBL_PEP=/cluster/projects/nn9036k/scrbkup/pgtk/reference_downloads/Homo_sapiens.GRCh38.pep.all.fa.gz,PGTK_SAMPLESHEET="$PWD/samples.csv",PGTK_RESULTS_DIR="$PWD/results",PGTK_NEXTFLOW="$HOME/bin/nextflow",PGTK_PYTHON="$(command -v python3)",PGTK_JAVA_MODULE=Java/21,PGTK_SLURM_LOG_TEMPLATE="$PWD/pgtk-wrapper-{job_id}.log" scratch.slurm -- --run_external_vcf_comparison true --external_vcf_dir "$SAREK_DIR" --external_vcf_suffix '.haplotypecaller.filtered.vcf.gz' --run_proteogenomic_validation true --maxquant_txt "$MAXQUANT_DIR") && printf '%s\n' "$JOB_ID" | tee .pgtk_current_job_id && echo "JOB_ID=$JOB_ID" && until [[ -f "pgtk-wrapper-${JOB_ID}.log" ]]; do sleep 2; done && tail -F "pgtk-wrapper-${JOB_ID}.log"
```

## Monitor and resume

Follow the wrapper log:

```bash
JOB_ID=$(cat .pgtk_current_job_id)
tail -F "pgtk-wrapper-${JOB_ID}.log"
```

Inspect Slurm accounting:

```bash
sacct -j "$JOB_ID" \
  --format=JobID,JobName,State,ExitCode,Elapsed,AllocCPUS,ReqMem,MaxRSS
```

Inspect incomplete Nextflow tasks:

```bash
awk -F '\t' 'NR==1 || $4 != "COMPLETED"' \
  "results/pipeline_trace-${JOB_ID}.tsv"
```

Keep the following for resume compatibility:

```text
/cluster/work/users/ash022/work/
.nextflow/
.nextflow.log*
results/
```

Do not delete the work directory while resume may be needed.

## Primary outputs

### Variant stages

```text
results/gvcf/<sample>.g.vcf.gz
results/vcf_raw/<sample>.raw.vcf.gz
results/vcf_normalized/<sample>.*
results/vcf_filtered/<sample>.filtered.vcf.gz
results/vcf_pass/<sample>.pass.vcf.gz
results/vep/<sample>.vep.vcf.gz
results/rna_validation/variants/<sample>.rna.validated.vcf.gz
```

Raw GenotypeGVCFs outputs and indexes are published before hard filtering to preserve direct caller-stage comparisons.

### Longitudinal progression

```text
results/progression_vcf/<sample>.nonbaseline_only.vep.vcf.gz
results/progression_vcf/<sample>.baseline_only.vep.vcf.gz
results/progression_vcf/<sample>.shared_with_baseline.vep.vcf.gz
results/progression_vcf/<sample>.subtraction.summary.tsv
```

```text
results/progression_biology/progression_biology.progression_alleles.tsv
results/progression_biology/progression_biology.progression_genes.tsv
results/progression_biology/progression_biology.go_enrichment.tsv
results/progression_biology/progression_biology.pairwise_go_contrasts.tsv
results/progression_biology/sets/progression_variant_sets.variant_set_go.tsv
results/progression_biology/sets/progression_variant_sets.summary.tsv
```

### Expression and GO

Default featureCounts settings include exon-level counting by `gene_id`, unstranded counting, minimum mapping quality 10, paired-read counting, primary alignments only, and exclusion of chimeric and multimapping reads.

```text
results/expression/per_sample/<sample>.gene_counts.tsv
results/expression/per_sample/<sample>.gene_counts.tsv.summary
results/expression/gene_expression.gene_expression.tsv
results/expression/gene_expression.summary.tsv
results/expression/go/expression_go.expression_ora.tsv
results/expression/go/expression_go.ranked_go.tsv
results/expression/go/expression_go.summary.tsv
```

Ranked expression GO uses:

```text
log2((progression TPM + 0.5) / (baseline TPM + 0.5))
```

The GO layers answer different questions:

- Expression GO tests gene abundance or coordinated transcript direction.
- Progression variant-set GO tests concentration of genes in baseline-subtracted RNA candidate sets.
- RNA-seq protein-altering variant GO tests genes with protein-altering VEP consequences at explicit RNA-derived stages.

Compact protein-altering variant GO outputs:

```text
results/variant_landscape/variant_landscape.summary.tsv
results/variant_landscape/variant_landscape.nonsynonymous_genes.tsv
results/variant_landscape/variant_landscape.go_significant.tsv
results/variant_landscape/variant_landscape.go_top.tsv
results/variant_landscape/variant_landscape.go_summary.tsv
results/variant_landscape/variant_landscape.report.md
```

The implementation writes significant terms and the top 100 terms per sample-stage. Row-count and file-size safeguards prevent accidental quadratic or unbounded output.

### Consolidated strict finding review

`BUILD_FINDING_IGV_REVIEWS` scans each staged BAM once and writes a flat evidence bundle under:

```text
results/igv/findings/finding_reviews/
```

Principal outputs:

```text
findings_manifest.tsv
priority_findings.tsv
priority_findings.bed
support_labels.bed
event_consolidation.tsv
event_regions.tsv
bam_manifest.tsv
read_classification.tsv
excluded_reads.tsv
consolidation_summary.txt
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

Track meanings:

- `exact_alt_unique`: all strict exact-ALT reads retained for analysis.
- `exact_alt_display`: capped exact-ALT subset for browser display.
- `reference_display`: capped clean-reference subset for browser display.
- `event_display`: browser-safe alignments overlapping candidate regions.

Analytical counts remain in `findings_manifest.tsv`. Browser-display BAMs may omit unrelated deletion-CIGAR records that trigger the bundled IGV.js deletion renderer. Genuine deletion-candidate evidence remains eligible for display.

Paired-end flags, mate coordinates, template length, and orientation remain in BAM records. A mate may not appear in a local report when the mate lies outside the displayed interval or was not retained in the selected evidence subset.

### Finding explorer and on-demand IGV.js

```text
results/igv/findings/finding_explorer/index.html
results/igv/findings/finding_explorer/server.py
results/igv/findings/finding_explorer/serve_explorer.sh
results/igv/findings/finding_explorer/explorer_config.json
results/igv/findings/finding_explorer/coverage_summary.txt
results/igv/findings/finding_explorer/partitions/all.jsonl.gz
results/igv/findings/finding_explorer/report_cache/
```

Direct-file mode:

```text
Open index.html through file:///
```

The direct page supports searching, filtering, validation status, ALT-supporting reads, reference-supporting reads, excluded or non-callable reads, total examined alignments, callable reads, and ALT fraction.

Server mode:

```bash
cd results/igv/findings/finding_explorer
export PGTK_IGV_REPORTS_IMAGE=/cluster/projects/nn9036k/scrbkup/pgtk/singularity_cache/quay.io-biocontainers-igv-reports-1.16.0--pyh7e72e81_0.img
./serve_explorer.sh "$PWD" 8765
```

Create a workstation tunnel:

```bash
ssh -L 8765:127.0.0.1:8765 ash022@login.saga.sigma2.no
```

Open:

```text
http://127.0.0.1:8765
```

Select a candidate and click **Open alignments**. The report loads only that sample’s `event_display`, `exact_alt_display`, and `reference_display` BAMs. Generated HTML reports are cached under `report_cache/`.

The default explorer filter is `ALT_SUPPORTED`. Select `All candidates` to inspect mixed, reference-only, and non-callable events.

### Global IGV bundle

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

### MultiQC

```text
results/multiqc/multiqc_report.html
results/multiqc/multiqc_report_data/
```

MultiQC includes standard QC plus compact custom sections for:

- Results navigation and evidence semantics.
- RNA findings and variant types by sample-stage.
- RNA-seq protein-altering variant GO summary and top-term plots.
- Expression GO and ranked baseline comparisons.
- Progression biology and progression variant-set GO.
- Codon evidence and read provenance.
- External-caller comparison when enabled.
- Proteogenomic evidence when enabled.
- Failure and resource summaries.

## Optional external VCF comparison

```bash
sbatch scratch.slurm -- \
  --run_external_vcf_comparison true \
  --external_vcf_dir /path/to/external/vcfs \
  --external_vcf_suffix .haplotypecaller.filtered.vcf.gz
```

Files are matched recursively using sample identifiers. Raw, PASS, and RNA-validated PGTK stages are compared separately under:

```text
results/comparison/external_vcf/
```

Ambiguous matches and empty external VCFs are rejected.

## Optional MaxQuant integration

PGTK does not execute MaxQuant. The workflow interprets an existing MaxQuant search performed with the exact custom FASTAs, canonical proteome, and contaminants database.

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
  --maxquant_txt /path/to/combined/txt
```

Principal outputs are written under:

```text
results/proteogenomics_validation/
```

File presence alone is not biological confirmation. Strong event confirmation requires sample-matched direct MS/MS evidence, consistency with the searched FASTAs, altered-residue or junction coverage, and absence from canonical reference proteins.

## Resource, retry, and failure model

Resource-related failures with exit codes 137, 140, or 143 permit three attempts:

```text
Attempt 1: 1x resources
Attempt 2: 2x resources
Attempt 3: 4x resources
```

Configured limits:

```text
Maximum CPUs: 32
Maximum memory: 512 GB
normal partition: <= 20 CPUs and <= 160 GB
bigmem partition: > 20 CPUs or > 160 GB
```

Failure and execution records:

```text
results/failure_logs/<job-id>/
results/failure_logs/failure_history.tsv
results/failure_logs/run_history.tsv
results/pipeline_trace-<job-id>.tsv
results/pipeline_report-<job-id>.html
results/pipeline_timeline-<job-id>.html
results/pipeline_dag-<job-id>.html
results/resource_usage-<job-id>.summary.tsv
results/resource_usage-<job-id>.warnings.tsv
results/resource_usage-<job-id>.report.md
```

## Repository file guide

The active implementation is defined by `main.nf`, `nextflow.config`, `scratch.slurm`, active scripts, and tests. Historical notes document decisions but do not override current code.

### Workflow, configuration, and execution

- `main.nf`: complete Nextflow DSL2 graph, process definitions, channels, optional branches, and output wiring.
- `nextflow.config`: executor profiles, Slurm resource routing, retry behavior, container settings, and process-specific resources.
- `scratch.slurm`: Saga submission wrapper, runtime preflight, resume execution, failure collection, and resource summaries.
- `run.sh`: local execution wrapper.
- `wsl.config`: WSL-specific local-development configuration.
- `samples.csv`: active samplesheet.
- `sample_lane0.csv`: lane-specific samplesheet fixture or example.
- `sample_laneNA.csv`: input fixture for lane handling.
- `samples_lane_split.csv`: lane-split output fixture or example.
- `multiqc_config.yaml`: MultiQC configuration and custom-content ordering.
- `LICENSE`: project license.
- `README.md`: current operational and repository documentation.

### Acquisition and input preparation

- `download_sra.sh`: downloads and validates SRA archives on a login node.
- `download_assets.sh`: downloads references, GO resources, and container assets on a login node.
- `downloadFastq.sh`: auxiliary or legacy FASTQ download helper.
- `ena_fastq_urls.txt`: ENA URL list for auxiliary downloads, when present.
- `generateSampleSheet.sh`: generates a samplesheet from available inputs.
- `laneSplit.py`: splits paired FASTQ data by lane and writes an updated samplesheet.

### Variant calling, annotation, and validation

- `validate_haplotype_shards.py`: verifies expected HaplotypeCaller shards before gather.
- `summarize_variant_stages.py`: summarizes raw, normalized, filtered, PASS, VEP, and validated call stages.
- `compare_external_vcf.py`: compares PGTK VCF stages with an external caller.
- `external_comparison.none.tsv`: placeholder emitted when external comparison is disabled.
- `annotate_vep.sh`: VEP annotation helper.
- `validate_rna_events.py`: validates RNA-observed variants, fusions, and splice events into structured evidence outputs.
- `validate_variant_codons.py`: validates reference and alternate codons and translated residues.
- `validate_variant_read_provenance.py`: records supporting-read provenance for variant calls.
- `merge_variant_validation.py`: merges variant validation tables.
- `analyze_codon_mismatches.py`: investigates codon and translation mismatches.
- `build_integrated_variant_evidence.py`: integrates RNA, codon, provenance, and optional proteomic evidence.
- `analyze_variant_landscape.py`: produces compact per-stage variant summaries, protein-altering genes, significant GO terms, top GO terms, and MultiQC content.
- `generate_pgtk_vcf_report_inputs.sh`: generates auxiliary VCF report-input summaries.

### Fusion, splice, peptide, and protein FASTA tools

- `analyze_chimeric_splice_peptides.py`: analyzes peptide support for fusion and splice-derived sequences.
- `validate_splice_junction_peptides.py`: validates peptides spanning splice junctions.
- `map_peptides_to_fasta.py`: maps peptides to custom and reference FASTAs.
- `annotate_variant_peptides.py`: annotates peptide-to-variant relationships.
- `validate_proteogenomic_reads.py`: validates read support for proteogenomic events and creates compact evidence BAMs.
- `proteogenomics_evidence_report.py`: writes structured proteogenomic evidence reports.

### Expression and biological interpretation

- `expression_go_analysis.py`: expression ORA, ranked GO, sample GO, and expression-related merging.
- `prepare_go_annotations.py`: builds validated gene-to-GO resources.
- `analyze_progression_biology.py`: progression allele, gene, and GO analysis.
- `compare_progression_pair.py`: pairwise progression comparison.
- `merge_progression_biology.py`: merges progression analysis outputs.
- `build_comparative_advantage_report.py`: creates the comparative biological report.
- `build_complete_report.py`: creates the complete findings report.
- `audit_pgtk_claims.py`: audits interpretation against explicit evidence layers.
- `audit_environment_hardcoding.py`: detects environment-specific hard-coding in active code and configuration.

### MultiQC content builders

- `build_compact_multiqc_content.py`: results guide, navigation, terminology, and compact integrated sections.
- `build_expression_multiqc_content.py`: expression and GO MultiQC content.
- `build_final_multiqc_content.py`: final report content assembly.
- `build_pgtk_multiqc_content.py`: core PGTK custom MultiQC modules.

### IGV evidence and explorer

- `build_igv_evidence_bundle.py`: global RNA and progression BED, BEDPE, BAM, batch, session, and manifest bundle.
- `build_finding_igv_reviews.py`: single-pass strict read classification and consolidated ALT, reference, and event BAM generation.
- `build_finding_explorer.py`: direct-open searchable HTML, compressed finding records, status logic, local report server, and on-demand IGV.js command generation.
- `serve_finding_explorer.sh`: starts the local explorer in the pinned offline IGV Reports container.
- `check_bedcov_progress.sh`: reports progress for auxiliary bedcov work.
- `run_bedcov_parallel.sh`: auxiliary parallel bedcov runner.

### Cohort and coverage utilities

- `group_compare_two_cohorts.py`: compares two configured cohorts.
- `group_exclusive_mutations.py`: identifies cohort-exclusive events or genes.
- `group_exclusive_mutations.sh`: shell wrapper for exclusive-event analysis.
- `cross_chromosome.py`: summarizes cross-chromosome fusion events.
- `annotate_and_pattern.py`: annotates exclusive genes and recurrence patterns.
- `mannwhitney_perchrom.py`: per-chromosome Mann-Whitney analysis from coverage summaries.
- `plot_coverage_ideogram.py`: chromosome-level coverage ideogram.
- `plot_coverage_per_chrom.py`: per-chromosome coverage plots.
- `plot_windowed_cn.py`: windowed copy-number-like coverage plots.

### Runtime, failure, and trace utilities

- `validate_runtime_inputs.py`: validates references, samples, containers, Python, featureCounts, IGV Reports, and Pysam runtime contracts.
- `collect_pipeline_failures.py`: writes failure ledgers and cumulative failure history.
- `analyze_pipeline_trace.py`: summarizes Nextflow trace resources and warnings.
- `RESOURCE_RETRY_MATRIX.tsv`: documented retry and resource behavior.
- `RESOURCE_CALIBRATION.md`: resource calibration guidance.
- `coll.sh`: local packaging or diagnostic collection helper. Not required by the workflow graph.

### Tests and validation

- `validate_pipeline_commands.sh`: complete static command, process, resource, output, and configuration contract validator.
- `test_full_pipeline_contract.py`: process, channel, observer, and script-name consistency tests.
- `test_reporting_redesign.py`: reporting and MultiQC redesign contract tests.
- `test_resource_configuration.py`: process count, resource configuration, explorer shell, database-free, and zero-limit tests.
- `test_boolean_parameters.py`: strict boolean CLI parsing tests.
- `test_compact_multiqc_design.py`: compact MultiQC design tests.
- `test_expression_go_analysis.py`: expression and variant-set GO regression tests.
- `test_progression_biology.py`: progression GO and propagation regression tests.
- `test_external_vcf_resolution.py`: external VCF matching and ambiguity tests.
- `test_finding_igv_reviews.py`: strict finding classification, identifier, deduplication, and flat-layout tests.
- `test_igv_evidence_bundle.py`: global IGV evidence bundle tests.
- `test_proteogenomic_read_runtime.py`: proteogenomic Pysam runtime contract test.
- `test_streaming_variant_validation.py`: streaming variant validation regression test.
- `validate_region_optimized_igv_fixture.sh`: region-optimized IGV fixture validation.
- `validate_single_pass_igv_fixture.sh`: single-pass review fixture validation.
- `validate_finding_explorer.sh`: explorer output validation.
- `probe_pipeline_cli.sh`: CLI parameter and parsing probe.

### Optional MaxQuant and placeholder inputs

- `maxquant_raw_file_map.none.tsv`: placeholder when no MaxQuant raw-file map is supplied.
- `mqpar.xml`: MaxQuant parameter XML retained for provenance, when present.
- `callparam.xml`: auxiliary call-parameter XML, when present.
- `extractionSummary.xml`: auxiliary extraction summary, when present.

### Documentation

Current high-level documentation:

- `PROJECT_SUMMARY.md`: project summary.
- `PUBLICATION_MAPPING.md`: publication and public-data mapping.
- `RUN_AND_REPORT.md`: operational run and reporting notes.
- `PIPELINE_VALIDATION_SEMANTICS.md`: evidence and validation semantics.
- `HISTORY_REGRESSION_CHECKLIST.md`: checklist for preserving prior fixes.

Implementation records:

- `PIPELINE_UPDATE_20260807.md`
- `FULL_PIPELINE_UPDATE_20260814.md`
- `FULL_UPDATE_NOTES.md`
- `PORTABLE_RESOURCE_IGV_UPDATE_20260815.md`
- `EFFICIENCY_REVIEW_20260816.md`
- `SINGLE_PASS_IGV_REDESIGN_20260817.md`
- `REGION_OPTIMIZED_IGV_UPDATE_20260817.md`
- `PORTABLE_POSTPROCESS_FIX_20260818.md`
- `FINDING_EXPLORER_DATABASE_FREE_20260818.md`
- `FINDING_EXPLORER_REDESIGN_20260818.md`
- `CURRENT_FIX_20260818.md`
- `CLEAN_RUNTIME_UPDATE.md`
- `MULTIQC_REDESIGN_20260821.md`
- `PIPELINE_UPDATE_20260821_PORTABLE_EXPLORER_VARIANT_GO.md`
- `PIPELINE_UPDATE_20260822_COMPACT_GO_EXPLORER.md`
- `PIPELINE_UPDATE_20260822_EXPLORER_STATUS_AND_GO.md`
- `PIPELINE_UPDATE_20260822_IGV_ALIGNMENT_BUTTON.md`

Focused validation records:

- `GO_OPTIMIZATION_VALIDATION.md`
- `GO_COMPLETE_SAMPLE_AND_FDR_VALIDATION.md`
- `RANKED_GO_VALIDATION.md`
- `PYTHON_RUNTIME_FIX_VALIDATION.md`
- `NEXTFLOW_CONFIG_PARSER_FIX.md`
- `NEXTFLOW_STRICT_CONFIG_FIX.md`
- `PHASE2_CLI_PARAMETERS.md`

Historical notes explain why changes were made. The current code and tests are authoritative.

## Interpretation rules

- RNA-supported variants are not DNA-confirmed mutations.
- Non-baseline-only RNA observations are not proven newly acquired variants.
- RNA-supported fusions are not confirmed genomic rearrangements.
- Read-supported splice transcripts do not prove a complete isoform structure.
- FASTA entries are search candidates, not peptide confirmation.
- Ranked GO and targeted signatures from the same expression matrix are not independent validation.
- GO terms overlap and are not independent pathways.
- External-caller overlap is general reproducibility unless gene-specific support is demonstrated.
- MaxQuant output presence is not event confirmation.
- Match Between Runs is not equivalent to direct MS/MS evidence.
- Transcript abundance does not establish protein abundance, pathway activity, metabolic flux, drug sensitivity, or causation.
- A single longitudinal subject does not support population-level differential-expression inference.

## Limitations

- RNA detectability depends on expression and coverage.
- Unexpressed genomic regions cannot be evaluated reliably.
- Alignment, reference, filtering, normalization, and soft-clipping choices affect RNA callsets.
- GATK local haplotype inference may produce a candidate even when the final BAM contains no individual alignment representing the exact allele.
- External comparisons require compatible references and normalized alleles.
- Baseline subtraction reports absence from the observed baseline RNA callset, not proven biological absence.
- Custom FASTAs enlarge the search space and require appropriate FDR control.
- Proteomic interpretation depends on the exact searched FASTAs, contaminants database, raw-file mapping, and search configuration.
- Proteomic and functional confirmation must remain distinct from RNA-level evidence.

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

## License

See [LICENSE](LICENSE).

## Contact

Open a GitHub issue or contact `animesh@fuzzylife.org`.
