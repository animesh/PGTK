# PGTK caller comparison and resource update

This update was prepared from the active pipeline package supplied after run 19029276.

Changed behavior:

- Publishes indexed raw GenotypeGVCFs calls under `results/vcf_raw`.
- Publishes indexed filtered calls under `results/vcf_filtered` and retains `results/vcf_pass`.
- Validates all HaplotypeCaller shards before gathering.
- Produces raw, PASS and RNA-validated QC and provenance under `results/qc/variant_stages`.
- Optionally compares PGTK raw calls with external VCFs under `results/comparison/external_vcf`.
- Accepts either the existing five-column samplesheet or a minimal `sample,srr` samplesheet.
- Retains independent per-sample exploratory proteogenomics FASTAs.
- Optimizes initial resources from observed run 19029276 while preserving retry scaling.

External comparison is optional. Files are matched by samplesheet `srr`:

```text
--external_vcf_dir /path/to/vcfs
--external_vcf_suffix .haplotypecaller.vcf.gz
```

Expected filenames include `SRR31089074.haplotypecaller.vcf.gz`.

Resource changes:

```text
FASTQC_RAW: 2 CPUs, 2 GB initial
FASTQC_TRIMMED: 2 CPUs, 2 GB initial
HAPLOTYPE_CALLER: 2 CPUs, 4 GB initial
GATHER_HAPLOTYPE_GVCF: 2 CPUs, 4 GB initial
GENOTYPE_FILTER: 2 CPUs, 4 GB initial
STRINGTIE_ASSEMBLY: 2 CPUs, 3 GB initial
VALIDATE_RNA_SPLICE_TRANSCRIPTS: 2 CPUs, 1 GB initial
VALIDATE_VARIANT_CODONS: 2 CPUs, 1 GB initial
```

All existing retry multipliers and hard caps remain active.

## Comparative evidence defaults

The run now expects these project-local directories by default:

```text
sarek/
txtMQMBR/
```

Sarek VCFs are found recursively by exact SRR filename, for example `SRR31089074.haplotypecaller.vcf.gz`.
MaxQuant validation is enabled by default and reads `peptides.txt`, `evidence.txt`, `msms.txt`, and `proteinGroups.txt` from `txtMQMBR`. `mqpar.xml` is resolved inside that directory or its parent. The contaminants FASTA is resolved from an explicit parameter, the MaxQuant results area, the project, or `$HOME/scripts/MaxQuant_v2.8.1.0/bin/conf/contaminants.fasta`.

Disable MaxQuant validation when needed with:

```text
--run_proteogenomic_validation false
```

## New evidence outputs

```text
results/progression_vcf/*.nonbaseline_only.vep.vcf.gz
results/progression_vcf/*.baseline_only.vep.vcf.gz
results/progression_vcf/*.shared_with_baseline.vep.vcf.gz
results/progression_vcf/*.subtraction.summary.tsv
results/comparative_advantage/
results/comparison/external_vcf/
results/igv/all_evidence/
results/qc/variant_stages/
results/qc/haplotype_shards/
results/multiqc/multiqc_report.html
```

The IGV bundle contains a combined event manifest, BED coordinates, event-specific indexed BAMs per sample, an IGV batch file, and an IGV session XML covering RNA-validated variants, progression variants, fusions, and splice events.

## Final default and FASTA corrections

- MaxQuant validation is disabled by default. Enable it only with `--run_proteogenomic_validation true`.
- Sarek/external caller comparison is disabled by default. Enable it only with `--run_external_vcf_comparison true`.
- `results/combined_fasta/<sample>.exploratory_proteogenomics.fasta` no longer embeds canonical proteins.
- The custom FASTA contains deduplicated variant, fusion, and splice-derived sequences only. Supply the canonical proteome separately in MaxQuant.
