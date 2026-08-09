# PGTK: RNA-seq Proteogenomics and Comparative Variant Evidence

PGTK is a Nextflow DSL2 workflow for exploratory RNA-seq proteogenomics. It processes paired-end RNA-seq data from local SRA archives, performs RNA-aware small-variant calling, detects fusions and novel splice-derived transcripts, generates sample-specific custom protein FASTAs, quantifies gene expression, performs GO analysis, supports longitudinal baseline subtraction, and optionally integrates external caller and MaxQuant evidence.

PGTK is intended for research use. RNA-observed variants can include germline, clonal, progression-associated, RNA-editing, alignment, expression-dependent, and technical events. They must not be described as clinically validated somatic mutations without independent DNA evidence.

## Quick start

```bash
git clone https://github.com/animesh/pgtk.git
cd pgtk

cat samples.csv
bash download_sra.sh
bash download_assets.sh

bash validate_pipeline_commands.sh \
  --project-dir . \
  --nextflow "$HOME/bin/nextflow"

sbatch scratch.slurm
```

The Slurm wrapper forwards arguments written after `--` to Nextflow and uses resume-compatible execution.

## Staged execution

### Stage 1: RNA pipeline, expression, GO, custom FASTAs, and reports

```bash
cd /cluster/projects/nn9036k/scrbkup/pgtk
sbatch scratch.slurm
```

### Stage 2: external Sarek comparison

```bash
cd /cluster/projects/nn9036k/scrbkup/pgtk
sbatch scratch.slurm -- \
  --run_external_vcf_comparison true \
  --external_vcf_dir /cluster/projects/nn9036k/scrbkup/pgtk/sarek \
  --external_vcf_suffix .haplotypecaller.vcf.gz
```

### Stage 3: search the PGTK FASTAs with MaxQuant

MaxQuant is not executed by Nextflow. Run MaxQuant separately using each sample-specific PGTK custom FASTA together with the canonical human proteome and the MaxQuant contaminants FASTA.

Custom FASTAs:

```text
results/combined_fasta/TK12.exploratory_proteogenomics.fasta
results/combined_fasta/TK13.exploratory_proteogenomics.fasta
results/combined_fasta/TK14.exploratory_proteogenomics.fasta
```

Copy or link the resulting MaxQuant text directory to:

```text
txtMQMBR/
```

It must contain:

```text
peptides.txt
evidence.txt
msms.txt
proteinGroups.txt
mqpar.xml
```

### Stage 4: integrate existing MaxQuant results

```bash
cd /cluster/projects/nn9036k/scrbkup/pgtk
sbatch scratch.slurm -- \
  --run_proteogenomic_validation true \
  --maxquant_txt /cluster/projects/nn9036k/scrbkup/pgtk/txtMQMBR \
  --maxquant_contaminants "$HOME/scripts/MaxQuant_v2.8.1.0/bin/conf/contaminants.fasta"
```

If raw-file names do not unambiguously map to samples, add:

```bash
--maxquant_raw_map /cluster/projects/nn9036k/scrbkup/pgtk/maxquant_raw_file_map.tsv
```

If needed, explicitly override other MaxQuant provenance inputs:

```bash
--maxquant_mqpar /absolute/path/to/mqpar.xml \
--maxquant_canonical_fasta /absolute/path/to/human_reviewed_isoforms.fasta
```

### Stage 5: run the scientific claim audit

Use the job ID from the completed pipeline invocation:

```bash
cd /cluster/projects/nn9036k/scrbkup/pgtk
python audit_pgtk_claims.py \
  --project-dir /cluster/projects/nn9036k/scrbkup/pgtk \
  --job-id JOB_ID \
  --output-prefix results/pgtk_claim_audit_JOB_ID
```

The audit intentionally keeps transcript signatures, ranked GO, progression-variant GO, Sarek context, and proteomic evidence as separate evidence layers. MaxQuant output presence is recognized, but detailed peptide-supported H1 to H5 aggregation remains separate from the transcript-level hypothesis audit.

## Current implementation

The workflow includes:

- Nextflow DSL2 execution with Saga Slurm routing
- local SRA ingestion and paired FASTQ generation
- raw and trimmed FastQC
- Trim Galore preprocessing
- STAR two-pass RNA alignment
- coordinate-sorted and indexed STAR BAMs
- Arriba fusion calling and RNA validation
- StringTie transcript assembly and gffcompare novelty classification
- splice-transcript RNA validation and TransDecoder protein prediction
- monolithic GATK SplitNCigarReads
- 24-way HaplotypeCaller scattering per sample
- shard validation, GatherVcfs, indexing, and GenotypeGVCFs
- publication of raw, filtered, PASS, VEP-annotated, and RNA-validated VCFs
- codon-level and supporting-read provenance validation
- per-sample non-subtracted custom protein FASTAs
- patient-aware baseline subtraction as a separate VCF/reporting branch
- shared, baseline-only, and non-baseline-only progression reports
- optional external Sarek or other caller comparison at raw, PASS, and RNA-validated stages
- optional MaxQuant evidence interpretation
- event-specific IGV BED, BEDPE, BAM, batch, and session outputs
- featureCounts gene quantification
- merged raw-count, CPM, and TPM expression matrix
- per-sample expression GO over-representation analysis
- progression-versus-baseline ranked expression GO
- common, complete-sample, and exclusive progression variant-set GO
- pairwise progression GO contrasts
- claim audit with strict separation of evidence layers
- MultiQC, comparative, failure, resource, and provenance reports

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

## Analysis model

```text
Local SRA archive
  -> paired FASTQ
  -> raw FastQC
  -> Trim Galore
  -> trimmed FastQC
  -> STAR two-pass alignment
       -> Arriba fusion branch
       -> coordinate-sorted BAM
            -> featureCounts expression branch
            -> StringTie/gffcompare splice branch
            -> MarkDuplicates
            -> SplitNCigarReads
            -> HaplotypeCaller, 24 shards
            -> shard validation
            -> GatherVcfs and indexing
            -> GenotypeGVCFs
            -> raw VCF publication
            -> hard filtering
            -> PASS VCF publication
            -> VEP annotation
            -> RNA validation
            -> codon/read-provenance validation
            -> variant protein FASTA
```

Custom FASTAs:

```text
variant proteins
+ fusion proteins
+ splice-derived proteins
-> exact amino-acid sequence deduplication
-> <sample>.exploratory_proteogenomics.fasta
```

Canonical proteins are not embedded in the custom FASTAs. Add the canonical proteome separately during the MaxQuant search.

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

Semantics:

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

Per-sample FASTAs are always generated independently of baseline subtraction. A non-baseline sample is compared only when its subject has exactly one baseline. Missing or multiple baselines are reported rather than silently resolved.

## Variant stages

```text
results/gvcf/<sample>.g.vcf.gz
results/vcf_raw/<sample>.raw.vcf.gz
results/vcf_filtered/<sample>.filtered.vcf.gz
results/vcf_pass/<sample>.pass.vcf.gz
results/vep/<sample>.vep.vcf.gz
results/rna_validation/variants/<sample>.rna.validated.vcf.gz
```

Raw GenotypeGVCFs output and its index are published before hard filtering. This preserves clean comparisons with Sarek and other caller stages.

RNA-supported variants are not automatically DNA-confirmed mutations. Baseline absence in RNA does not prove biological or DNA-level absence.

## Longitudinal progression branch

For each valid progression-versus-baseline pair:

```text
results/progression_vcf/<sample>.nonbaseline_only.vep.vcf.gz
results/progression_vcf/<sample>.baseline_only.vep.vcf.gz
results/progression_vcf/<sample>.shared_with_baseline.vep.vcf.gz
results/progression_vcf/<sample>.subtraction.summary.tsv
```

The branch reports non-baseline-only, baseline-only, and shared variants separately.

Progression biology outputs:

```text
results/progression_biology/progression_biology.progression_alleles.tsv
results/progression_biology/progression_biology.progression_genes.tsv
results/progression_biology/progression_biology.go_enrichment.tsv
results/progression_biology/progression_biology.pairwise_go_contrasts.tsv
results/progression_biology/sets/progression_variant_sets.variant_set_go.tsv
results/progression_biology/sets/progression_variant_sets.summary.tsv
```

These tables retain complete tested GO terms. MultiQC displays summaries only.

## Gene expression and GO

featureCounts parameters:

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

Expression ORA asks which GO terms are overrepresented among expressed genes. Ranked expression GO asks whether genes in a GO term shift upward or downward in progression versus baseline using:

```text
log2((progression TPM + 0.5) / (baseline TPM + 0.5))
```

Progression-variant GO is separate. It asks whether genes carrying non-baseline-only RNA-observed variants are overrepresented in GO categories. It does not test expression upregulation.

Important distinction:

```text
Expression-ranked GO
  tests coordinated transcript direction

Progression-variant GO
  tests concentration of variant-bearing genes
```

GO terms overlap extensively and must not be counted as independent pathways.

## Claim audit

Run:

```bash
python audit_pgtk_claims.py \
  --project-dir . \
  --job-id JOB_ID \
  --output-prefix results/pgtk_claim_audit_JOB_ID
```

Principal outputs:

```text
pgtk_claim_audit_JOB_ID.report.md
pgtk_claim_audit_JOB_ID.claims.tsv
pgtk_claim_audit_JOB_ID.supported_claims.tsv
pgtk_claim_audit_JOB_ID.unsupported_claims.tsv
pgtk_claim_audit_JOB_ID.insights.tsv
pgtk_claim_audit_JOB_ID.article_signatures.tsv
pgtk_claim_audit_JOB_ID.novel_findings.tsv
pgtk_claim_audit_JOB_ID.hypothesis_scores.tsv
pgtk_claim_audit_JOB_ID.hypothesis_evidence_matrix.tsv
pgtk_claim_audit_JOB_ID.hypothesis_variant_hits.tsv
pgtk_claim_audit_JOB_ID.hypothesis_novel_findings.tsv
pgtk_claim_audit_JOB_ID.checks.tsv
pgtk_claim_audit_JOB_ID.json
```

The audit evaluates five predefined hypotheses:

```text
H1 proteostasis and translation
H2 metabolic rewiring
H3 surface, glycan, extracellular-matrix and adhesion remodeling
H4 MYC, IRF4, DNA repair and progression
H5 PCNA stress and ATX-101-associated biomarkers
```

Evidence labels distinguish:

```text
SUPPORTED_WITHIN_EXPRESSION_LAYER
PARTIALLY_SUPPORTED_WITHIN_EXPRESSION_LAYER
NOT_SUPPORTED_OR_MIXED
NOT_TESTED
```

Targeted signatures and ranked GO use the same expression matrix. Their agreement is internal expression concordance, not independent validation. Sarek gives general callset-reproducibility context unless hypothesis-gene support is explicitly demonstrated. Proteomics and functional assays remain separate evidence layers.

The audit uses canonical protein-coding cytosolic and mitochondrial ribosome lists, excludes pseudogenes from those signatures, requires TPM of at least 1 in one compared sample, and labels olfactory-receptor GO findings for artifact review.

MHC class II reduction is reported only if at least two significant negative MHC class II GO terms are present. The audit records exact GO IDs, MeanScores, and FDR values.

The H3 terminology is deliberately `surface_glycan_adhesion_remodeling`; the audit does not infer glycan loss, marrow escape, or increased invasion from mixed transcript changes.

## Optional external-caller comparison

```bash
sbatch scratch.slurm -- \
  --run_external_vcf_comparison true \
  --external_vcf_dir ./sarek \
  --external_vcf_suffix .haplotypecaller.vcf.gz
```

Files are matched recursively using the samplesheet `srr` field. Results are written under:

```text
results/comparison/external_vcf/
```

Comparisons are performed separately for raw, PASS, and RNA-validated PGTK stages.

## Optional MaxQuant validation

MaxQuant validation is disabled by default:

```text
--run_proteogenomic_validation false
```

PGTK does not run MaxQuant. It interprets an existing MaxQuant search performed with the exact PGTK custom FASTAs, canonical proteome, and contaminants database.

Enable integration:

```bash
sbatch scratch.slurm -- \
  --run_proteogenomic_validation true \
  --maxquant_txt /cluster/projects/nn9036k/scrbkup/pgtk/txtMQMBR \
  --maxquant_contaminants "$HOME/scripts/MaxQuant_v2.8.1.0/bin/conf/contaminants.fasta"
```

The contaminants FASTA is passed through `params.maxquant_contaminants`, resolved by `resolveContaminants`, and staged automatically by Nextflow into processes that require it. It is not copied into the MaxQuant text directory. The supplied file should be the same contaminants database used in the MaxQuant search.

If `--maxquant_contaminants` is omitted, the workflow requires exactly one `contaminants.fasta` or `contaminants.fa` under one of these locations:

```text
--maxquant_txt directory
parent of --maxquant_txt
project directory
$HOME/scripts/MaxQuant_v2.8.1.0/bin/conf/contaminants.fasta
```

`mqpar.xml` must report `includeContaminants=True` and must contain searched FASTA paths. Canonical FASTAs are resolved from `mqpar.xml` unless overridden with `--maxquant_canonical_fasta`.

The branch executes:

```text
VALIDATE_MAXQUANT_INPUTS
MAP_MAXQUANT_PEPTIDES
ANNOTATE_MAXQUANT_VARIANTS
ANALYZE_MAXQUANT_JUNCTIONS
VALIDATE_MAXQUANT_SPLICE_JUNCTIONS
BUILD_PROTEOGENOMICS_EVIDENCE_REPORT
BUILD_INTEGRATED_VARIANT_EVIDENCE
VALIDATE_PROTEOGENOMIC_READS
PREPARE_FINAL_MULTIQC_CONTENT
MULTIQC_FINAL
```

It distinguishes:

- direct MS/MS from Match Between Runs
- sample-matched from cross-sample evidence
- canonical peptide support from event-specific support
- altered-residue variant peptides
- fusion-junction peptides
- splice-junction peptides
- unique from ambiguous mappings
- accepted, unresolved, and rejected associations

Principal outputs:

```text
results/proteogenomics_validation/maxquant_inputs.validated.txt
results/proteogenomics_validation/peptide_fasta_mapping.mapping.tsv
results/proteogenomics_validation/variant_peptide_annotation.detailed.tsv
results/proteogenomics_validation/validated_splice_junctions.detailed.tsv
results/proteogenomics_validation/proteogenomics_evidence.variants.tsv
results/proteogenomics_validation/proteogenomics_evidence.junctions.tsv
results/proteogenomics_validation/proteogenomics_evidence.direct_msms_variants.tsv
results/proteogenomics_validation/proteogenomics_evidence.direct_msms_junctions.tsv
results/proteogenomics_validation/proteogenomics_evidence.mbr_only_variants.tsv
results/proteogenomics_validation/proteogenomics_evidence.mbr_only_junctions.tsv
results/proteogenomics_validation/proteogenomics_evidence.sample_matched_direct_msms_variants.tsv
results/proteogenomics_validation/proteogenomics_evidence.sample_matched_direct_msms_junctions.tsv
results/proteogenomics_validation/integrated_variant_evidence.strict.tsv
results/proteogenomics_validation/proteogenomics_evidence.report.md
results/proteogenomics_validation/proteogenomics_evidence.summary.txt
```

Strict event confirmation requires sample-matched direct MS/MS evidence, consistency with the searched FASTAs, altered-residue or junction coverage, and absence from configured canonical reference proteins. File presence alone is not biological confirmation.

## IGV outputs

Global RNA/progression bundle:

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

MaxQuant-enabled read-validation outputs are written under:

```text
results/proteogenomics_validation/read_validation/
```

## MultiQC

```text
results/multiqc/multiqc_report.html
results/multiqc/multiqc_report_data/
```

MultiQC includes standard QC and custom sections for RNA findings, validation failures, codon evidence, read provenance, progression biology, comparative evidence, expression GO, proteogenomics evidence when enabled, and validation semantics.

Complete term-level GO tables remain in TSV outputs. MultiQC contains summaries and selected content.

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

## Project layout

```text
pgtk/
├── main.nf
├── nextflow.config
├── scratch.slurm
├── samples.csv
├── audit_pgtk_claims.py
├── reference_downloads/
├── sra_cache/
├── singularity_cache/
├── sarek/
├── txtMQMBR/
└── results/
```

Keep the Nextflow work directory and `.nextflow` metadata when resume compatibility is required.

## Validation

```bash
bash validate_pipeline_commands.sh \
  --project-dir . \
  --nextflow /cluster/home/ash022/bin/nextflow
```

Validation covers process declarations, resource selectors, Python and shell syntax, Nextflow inspection, runtime paths, raw VCF publication, shard completeness, retry behavior, progression outputs, custom FASTA composition, GO wiring, optional branch defaults, MaxQuant provenance checks, IGV reports, and MultiQC sections.


## Interpretation rules

- RNA-supported variants are not DNA-confirmed mutations.
- Non-baseline-only RNA observations are not proven newly acquired variants.
- RNA-supported fusions are not confirmed genomic rearrangements.
- Read-supported splice transcripts do not prove complete isoform structure.
- FASTA entries are search candidates, not peptide confirmation.
- Ranked GO and targeted signatures from the same expression matrix are not independent validation.
- GO terms are overlapping ontology categories, not independent pathways.
- Sarek overlap is general caller reproducibility unless gene-specific support is shown.
- MaxQuant output presence is not event confirmation.
- Match Between Runs is not equivalent to direct MS/MS.
- Sample-matched event-spanning peptides are stronger than canonical protein-group detection.
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
