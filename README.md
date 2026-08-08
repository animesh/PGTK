# PGTK: RNA-seq Proteogenomics and Comparative Variant Evidence

PGTK is a Nextflow DSL2 workflow for exploratory RNA-seq proteogenomics. It processes paired-end RNA-seq data from local SRA archives, performs RNA-aware small-variant calling, detects fusions and novel splice-derived transcripts, generates sample-specific custom protein FASTAs, and optionally integrates external caller results and MaxQuant evidence.

The workflow is designed for research use. RNA-derived variant calls may include germline, clonal, progression-associated, RNA-editing, alignment, and technical events. They must not be interpreted as clinically validated somatic variants without independent evidence.

## Current validated implementation

The current workflow has:

- 53 unique Nextflow processes
- Monolithic GATK `SplitNCigarReads`
- 24-way scattered HaplotypeCaller execution per sample
- Explicit GVCF gathering and indexing
- Published raw, filtered, PASS, VEP-annotated, and RNA-validated VCF stages
- RNA fusion and splice-transcript validation
- Codon-level and supporting-read provenance validation
- Patient-aware baseline subtraction when metadata is available
- Custom FASTAs containing only variant, fusion, and splice-derived proteins
- Optional Sarek or external-caller comparison
- Optional MaxQuant evidence integration
- Comparative evidence, IGV, resource, failure, and MultiQC reports
- Dynamic retry resources and normal/bigmem partition routing on Saga

The latest complete validation check  [passed](#validation).

## Publication and public data

The biological context is described in:

**Multiple Myeloma Cells with Increased Proteasomal and ER Stress Are Hypersensitive to ATX-101, an Experimental Peptide Drug Targeting PCNA**

- Journal: Cancers
- Year: 2024
- Volume: 16
- Issue: 23
- Article: 3963
- DOI: https://doi.org/10.3390/cancers16233963
- RNA-seq BioProject: `PRJNA1176350`
- Proteomics: PRIDE `PXD033531` and `PXD033510`

Please cite the study when using the associated data or derived results.

## Analysis model

The main RNA-seq path is:

```text
SRA archive
  -> paired FASTQ
  -> raw FastQC
  -> Trim Galore
  -> trimmed FastQC
  -> STAR two-pass alignment
       -> Arriba fusion calling
       -> coordinate-sorted BAM
            -> StringTie assembly
            -> gffcompare novelty classification
            -> validated splice transcripts
            -> TransDecoder protein prediction
            -> MarkDuplicates
            -> SplitNCigarReads
            -> HaplotypeCaller, 24 shards
            -> shard validation
            -> GatherVcfs and IndexFeatureFile
            -> GenotypeGVCFs
            -> raw VCF publication
            -> hard filtering
            -> PASS VCF publication
            -> VEP annotation
            -> RNA validation
            -> codon and read-provenance validation
            -> variant protein FASTA
```

Custom per-sample FASTAs are constructed from:

```text
variant proteins
+ fusion proteins
+ splice-derived proteins
-> exact amino-acid sequence deduplication
-> sample.exploratory_proteogenomics.fasta
```

Canonical proteins are not embedded in these custom FASTAs. Add the canonical proteome separately in MaxQuant.

## Samplesheet

Only `sample` and `srr` are required.

Minimal format:

```csv
sample,srr
TK12,SRR31089074
TK13,SRR31089073
TK14,SRR31089072
```

Full longitudinal format:

```csv
sample,srr,TK,Group,baseline
TK12,SRR31089074,patient1,resistant,true
TK13,SRR31089073,patient1,sensitive,false
TK14,SRR31089072,patient1,sensitive,false
```

Metadata semantics:

- `sample`: unique sample identifier used in filenames and outputs
- `srr`: SRA run accession
- `TK`: patient or subject identifier used for within-subject baseline matching
- `Group`: biological metadata used in reports
- `baseline`: `true` for the comparison reference and `false` otherwise

Defaults when optional columns are absent or empty:

```text
TK       = sample
Group    = sample
baseline = false
```

With a minimal samplesheet, all main analyses and per-sample FASTA generation run normally. Baseline subtraction is skipped because no baseline is defined.

When longitudinal metadata is supplied, each non-baseline sample is compared with the `baseline=true` sample having the same `TK` value. The pipeline reports missing or ambiguous baselines rather than silently performing an invalid comparison.

## Variant stages

The workflow preserves distinct calling and validation stages:

```text
results/gvcf/<sample>.g.vcf.gz
results/vcf_raw/<sample>.raw.vcf.gz
results/vcf_filtered/<sample>.filtered.vcf.gz
results/vcf_pass/<sample>.pass.vcf.gz
results/vep/<sample>.vep.vcf.gz
results/rna_validation/variants/<sample>.rna.validated.vcf.gz
```

This separation allows direct comparison with other callers without mixing calling, filtering, annotation, and RNA-validation effects.

## HaplotypeCaller strategy

PGTK uses an RNA-aware path:

```text
STAR BAM
-> MarkDuplicates
-> SplitNCigarReads
-> HaplotypeCaller in GVCF mode
-> GatherVcfs
-> IndexFeatureFile
-> GenotypeGVCFs
```

`SplitNCigarReads` is intentionally monolithic. Earlier interval-scattering experiments were removed because the tool is a context-sensitive, two-pass BAM transformation and scattered reconstruction introduced unnecessary correctness and maintenance risks.

HaplotypeCaller remains interval-scattered because that stage is naturally parallelizable. Every shard is checked before gathering.

Important caller settings are configurable, including calling confidence, soft-clipped-base handling, PCR indel model, thread count, and optional dbSNP resources.

## Baseline and progression analysis

Per-sample custom FASTAs are always non-subtracted. Baseline subtraction is a separate VCF and reporting branch.

For every valid non-baseline versus baseline pair, PGTK writes:

```text
results/progression_vcf/<sample>.nonbaseline_only.vep.vcf.gz
results/progression_vcf/<sample>.baseline_only.vep.vcf.gz
results/progression_vcf/<sample>.shared_with_baseline.vep.vcf.gz
results/progression_vcf/<sample>.subtraction.summary.tsv
```

These outputs distinguish:

- variants present only in the non-baseline sample
- variants present only in the baseline sample
- variants shared by both samples

This branch adds longitudinal biological interpretation that is not provided by independent per-sample HaplotypeCaller VCFs.

## RNA-event branches

### Variants

PGTK validates VEP-annotated calls using RNA evidence, including depth, alternate-read count, and alternate-allele fraction.

### Fusions

Arriba calls chimeric events from the STAR chimeric BAM. Accepted events are validated and translated into fusion-derived protein sequences with pVACfuse.

### Splicing

StringTie assembles expressed transcripts. gffcompare identifies selected transcript classes, RNA evidence is validated, and TransDecoder predicts splice-derived proteins.

These fusion and transcript-structure outputs are not representable as ordinary small-variant VCF records and therefore can complement DNA-oriented callers.

## Validation semantics

Variant-codon results use three top-level outcomes:

```text
VARIANT_CODON_VALIDATED
VARIANT_EVIDENCE_PARTIAL
VALIDATION_FAILED
```

`VARIANT_CODON_VALIDATED` requires reference agreement, sample-matched RNA alternate support, and a directly comparable coding consequence.

`VARIANT_EVIDENCE_PARTIAL` retains valid genome and RNA evidence when a consequence cannot safely be represented as a single codon comparison.

`VALIDATION_FAILED` requires one or more explicit failure codes.

Synonymous consequences are accepted as `SYNONYMOUS_CODON_TRANSLATION_CONFIRMED` when reference and alternate codons translate to the same amino acid. An empty VEP alternate-amino-acid field is accepted for synonymous VEP representations.

Primary failure categories are mutually exclusive. Overlapping failure occurrences are reported separately and are explicitly non-additive.

Strict integrated proteogenomic evidence requires:

- codon strict-integration eligibility
- sample-matched direct MS/MS evidence
- a search-consistent altered-residue peptide
- absence from the configured canonical reference protein sets

## Custom FASTA outputs

Main per-sample custom FASTAs:

```text
results/combined_fasta/<sample>.exploratory_proteogenomics.fasta
```

These contain only deduplicated non-canonical event-derived sequences:

```text
variant proteins
fusion proteins
splice-derived proteins
```

Additional component FASTAs are retained under:

```text
results/variant_fasta/
results/fusion_fasta/
results/splice_fasta/
```

For MaxQuant, search each custom FASTA together with:

- the canonical proteome FASTA
- the contaminants FASTA
- any other explicitly intended reference database

## Optional external-caller comparison

External comparison is disabled by default.

Enable it explicitly:

```bash
sbatch scratch.slurm -- \
  --run_external_vcf_comparison true \
  --external_vcf_dir /cluster/projects/nn9036k/scrbkup/pgtk/sarek \
  --external_vcf_suffix .haplotypecaller.vcf.gz
```

The external directory is searched recursively. Files are matched to samples using the samplesheet `srr` field, for example:

```text
SRR31089074.haplotypecaller.vcf.gz
```

Comparison outputs are written under:

```text
results/comparison/external_vcf/
```

Reports include:

- shared alleles
- PGTK-only alleles
- external-only alleles
- reciprocal overlap percentages
- Jaccard similarity
- SNP and indel counts
- genotype concordance
- discordant-call coordinates for review

The comparison runs separately for PGTK raw, PASS, and RNA-validated stages.

## Optional MaxQuant validation

MaxQuant validation is disabled by default.

Enable it only after MaxQuant has searched the newly generated custom FASTAs:

```bash
sbatch scratch.slurm -- \
  --run_proteogenomic_validation true \
  --maxquant_txt /cluster/projects/nn9036k/scrbkup/pgtk/txtMQMBR
```

The MaxQuant text directory must contain:

```text
peptides.txt
evidence.txt
msms.txt
proteinGroups.txt
```

`mqpar.xml` is resolved from that directory or its parent. The optional branch validates the actual searched FASTA provenance before interpreting peptide evidence.

The MaxQuant branch performs:

- input and search-provenance validation
- peptide mapping to every searched FASTA
- altered-residue variant-peptide annotation
- fusion and splice peptide analysis
- translated exon-junction validation
- sample-matched direct MS/MS classification
- Match Between Runs classification
- strict integrated variant evidence
- read-level validation and IGV output generation

## Recommended staged execution

### Stage 1: Generate RNA results and custom FASTAs

Both optional branches are disabled by default:

```bash
cd /cluster/projects/nn9036k/scrbkup/pgtk
sbatch scratch.slurm
```

### Stage 2: Add external-caller comparison

After Stage 1 finishes:

```bash
cd /cluster/projects/nn9036k/scrbkup/pgtk
sbatch scratch.slurm -- \
  --run_external_vcf_comparison true \
  --external_vcf_dir /cluster/projects/nn9036k/scrbkup/pgtk/sarek \
  --external_vcf_suffix .haplotypecaller.vcf.gz
```

MaxQuant remains disabled. Unchanged tasks are reused through `-resume`.

### Stage 3: Search custom FASTAs with MaxQuant

Use the custom FASTAs in `results/combined_fasta/` together with the canonical and contaminant FASTAs.

### Stage 4: Add MaxQuant evidence

After MaxQuant writes `txtMQMBR`:

```bash
cd /cluster/projects/nn9036k/scrbkup/pgtk
sbatch scratch.slurm -- \
  --run_proteogenomic_validation true \
  --maxquant_txt /cluster/projects/nn9036k/scrbkup/pgtk/txtMQMBR \
  --run_external_vcf_comparison true \
  --external_vcf_dir /cluster/projects/nn9036k/scrbkup/pgtk/sarek \
  --external_vcf_suffix .haplotypecaller.vcf.gz
```

Keep the work directory, `.nextflow` metadata, and result directory between stages so that compatible tasks can be resumed.

## Comparative evidence reports

The comparative report is written under:

```text
results/comparative_advantage/
```

Principal files include:

```text
comparative_advantage.report.md
comparative_advantage.variant_stage_inventory.tsv
comparative_advantage.fasta_inventory.tsv
comparative_advantage.rna_event_inventory.tsv
comparative_advantage.external_caller_comparison.tsv
comparative_advantage.multiqc_summary.tsv
```

These reports document advantages beyond simple variant discovery, including:

- RNA-supported call prioritization
- longitudinal baseline subtraction
- fusion discovery
- novel splice-transcript discovery
- custom protein-sequence generation
- codon and read-provenance validation
- MaxQuant peptide evidence
- IGV-ready event evidence

## IGV outputs

PGTK produces Nextflow-wired visualization bundles under:

```text
results/igv/all_evidence/
```

The bundle may include:

```text
pgtk_igv.events.tsv
pgtk_igv.events.bed
pgtk_igv.sample_manifest.tsv
pgtk_igv.igv.batch.txt
pgtk_igv.igv.session.xml
pgtk_igv.summary.txt
pgtk_igv.<sample>.events.bam
pgtk_igv.<sample>.events.bam.bai
```

It covers RNA-validated variants, progression-specific variants, fusion events, splice events, and proteogenomic evidence when available.

## MultiQC

The final MultiQC report is written to:

```text
results/multiqc/multiqc_report.html
results/multiqc/multiqc_report_data/
```

It includes standard QC plus custom sections for:

- complete RNA and proteogenomics findings
- RNA validation failures and explanations
- codon validation and mismatch investigation
- strict integrated variant evidence
- supporting-read provenance
- proteogenomics evidence and classification
- read-level validation
- validation and failure semantics
- comparative biological advantage
- comparative summary

The final MultiQC report is produced even when MaxQuant and external comparison are disabled. Optional sections report their disabled or unavailable status explicitly.

## Resource and retry model

The production resource profile is based on observed Saga execution. Important initial allocations include:

```text
FASTQC:                         2 CPUs, 2 GB
HaplotypeCaller shard:          2 CPUs, 4 GB
Gather Haplotype GVCF:          2 CPUs, 4 GB
Genotype and filtering:         2 CPUs, 4 GB
StringTie:                      2 CPUs, 3 GB
RNA splice validation:          2 CPUs, 1 GB
Variant codon validation:       2 CPUs, 1 GB
STAR alignment:                32 CPUs, 64 GB
MarkDuplicates:                 2 CPUs, 48 GB
SplitNCigarReads:               2 CPUs, 24 GB
```

Resource-related failures with exit codes 137, 140, or 143 allow three total attempts. Deterministic command and configuration failures terminate immediately:

```text
Attempt 1: 1x initial resources
Attempt 2: 2x initial resources
Attempt 3: 4x initial resources
```

Absolute retry limits:

```text
Maximum CPUs: 32
Maximum memory: 512 GB
```

Partition routing is recalculated per attempt:

```text
normal: request is at most 20 CPUs and at most 160 GB
bigmem: request exceeds 20 CPUs or 160 GB
```

GATK Java maximum heap is derived from 80 percent of effective task memory. CPU retry scaling is used only for tools with explicit thread controls. Effectively serial tools retain their calibrated CPU count while memory can escalate.

## Failure and provenance records

The wrapper records failures under:

```text
results/failure_logs/<job-id>/
```

Records include:

- failure ledger
- run summary
- pipeline trace
- Nextflow log
- controller log
- failed task commands
- task standard output and error
- exit codes
- effective resources and partition
- timestamps

Cumulative histories are maintained under:

```text
results/failure_logs/failure_history.tsv
results/failure_logs/run_history.tsv
```

A pipeline can complete successfully after a retry while retaining the earlier failed-attempt evidence.

Post-run execution artifacts include:

```text
results/pipeline_trace-<job-id>.tsv
results/pipeline_report-<job-id>.html
results/pipeline_timeline-<job-id>.html
results/pipeline_dag-<job-id>.html
results/resource_usage-<job-id>.summary.tsv
results/resource_usage-<job-id>.warnings.tsv
results/resource_usage-<job-id>.report.md
```

## Important implementation details

- SRA archives are local runtime inputs. Compute tasks do not require internet access.
- References and container paths are supplied at runtime and validated before task submission.
- Apptainer is used for process execution.
- STAR produces the chimeric BAM needed by Arriba.
- Sorted BAMs are used by StringTie, validation, and IGV branches.
- `SplitNCigarReads` remains monolithic for correctness.
- HaplotypeCaller is scattered and every shard is validated before gathering.
- `GatherVcfs` is followed by `IndexFeatureFile`.
- Raw GenotypeGVCFs output is published before filtering.
- VEP keeps all overlapping transcript consequences.
- Empty biological result sets generate valid empty outputs where appropriate.
- Nonzero tool failures remain fatal unless successfully resolved by the retry policy.
- Custom FASTAs exclude canonical sequences and deduplicate exact amino-acid sequences.
- Baseline subtraction never modifies the per-sample custom FASTAs.
- External-caller and MaxQuant branches require explicit `true` flags.

## Project layout

```text
pgtk/
├── main.nf
├── nextflow.config
├── scratch.slurm
├── samples.csv
├── reference_downloads/
├── sra_cache/
├── singularity_cache/
├── sarek/                 # optional external VCFs
├── txtMQMBR/              # optional MaxQuant text output
└── results/
```

Large Nextflow work files are kept outside the repository, for example:

```text
/cluster/work/users/ash022/work
```

Do not delete the work directory when resume compatibility is required.

## Validation

Run:

```bash
cd /cluster/projects/nn9036k/scrbkup/pgtk

bash validate_pipeline_commands.sh \
  --project-dir /cluster/projects/nn9036k/scrbkup/pgtk \
  --nextflow /cluster/home/ash022/bin/nextflow
```

The current validated version reports:

```text
PASS: 180
FAIL: 0
RESULT: PASSED
```

Validation covers:

- all 53 process declarations
- one resource selector per process
- Python and shell syntax
- Nextflow inspection
- reference and runtime preflight wiring
- raw VCF publication
- HaplotypeCaller shard completeness
- retry and partition semantics
- dynamic GATK heap wiring
- progression subtraction outputs
- custom FASTA composition
- optional external-caller and MaxQuant defaults
- comparative and IGV report wiring
- MultiQC custom sections
- known historical regression patterns

## Historical regressions removed

The production workflow no longer contains:

- scattered `SplitNCigarReads`
- custom split-BAM gathering
- boundary duplicate normalization
- split-BAM equivalence branches
- unsupported `samtools cat -@`
- Python execution assumptions inside samtools-only containers
- artificial global `maxForks` limits
- canonical proteins embedded in custom FASTAs
- default-enabled Sarek or MaxQuant validation

Historical failures are explicitly checked, including:

- multi-output process objects used as channels
- missing tools inside containers
- incomplete GVCF indexing
- incorrect BAM index names
- HTML entity corruption
- malformed or missing runtime inputs
- missing MaxQuant sentinel files
- incorrect synonymous-codon failure classification
- failure logging skipped because of shell exit behavior
- trace values represented as `-`
- resource selectors missing for newly added processes

## Interpreting comparisons with Sarek

Sarek and PGTK differ biologically and computationally. Sarek provides an independent HaplotypeCaller workflow, while PGTK applies RNA-aware preprocessing, filtering, RNA evidence, transcript translation, longitudinal subtraction, fusion and splice discovery, and proteogenomic integration.

Previous comparison of the same three samples showed:

```text
Approximately 80% to 83% of PGTK raw alleles were found by Sarek
Approximately 82% to 85% of PGTK PASS alleles were found by Sarek
Approximately 80% to 83% of PGTK RNA-validated alleles were found by Sarek
Approximately 97% to 98% genotype agreement among shared alleles
```

Sarek therefore provides strong independent support for most PGTK small variants. PGTK adds biological interpretation beyond the caller overlap:

- stringent RNA support
- progression-specific subtraction
- fusion and splice events
- protein-sequence generation
- codon consistency
- read provenance
- MaxQuant peptide evidence
- IGV-ready event bundles

## Limitations

- RNA-derived variants are not equivalent to validated DNA somatic variants.
- Transcript abundance strongly affects variant detectability.
- Regions without RNA expression cannot be evaluated reliably.
- Alternative alignment, reference, soft-clipping, and filtering choices can change callsets.
- Baseline subtraction means absence from the observed baseline RNA callset, not proven biological absence.
- External VCF comparison should use compatible references and normalized alleles.
- MaxQuant evidence depends on the exact searched FASTA set and search configuration.
- Custom FASTAs can increase the search space and require appropriate false-discovery-rate control.


## License

See [LICENSE](LICENSE) for details.

## Contact

For questions about the workflow or data, open an issue or contact `animesh@fuzzylife.org`.