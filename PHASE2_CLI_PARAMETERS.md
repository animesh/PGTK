# Phase 2 generic finding-review parameters

All finding selection and IGV Reports controls are Nextflow CLI parameters. No gene is prioritized by default.

## Defaults

```text
--finding_classes rna_variant,progression_variant
--finding_primary_class_order rna_variant,progression_variant
--finding_priority_mode all
--finding_priority_genes ''
--finding_priority_impacts ''
--finding_priority_consequences ''
--finding_review_mapq 20
--finding_review_baseq 20
--finding_review_reference_reads 20
--read_validation_padding 150
--generate_priority_igv_reports true
--igv_report_limit 0
--igv_report_timeout_seconds 600
--igv_report_title_prefix 'PGTK finding'
```

`--igv_report_limit 0` means all findings selected by the priority mode. A positive value limits output to the first N deterministically ranked findings. The previous value 500 was an operational cap, not a biological threshold.

`--finding_priority_mode all` includes every consolidated finding. `filter` includes a finding when its gene, impact, or consequence matches any corresponding nonempty CLI list.

## Examples

Generate reports for all consolidated findings:

```bash
sbatch scratch.slurm -- --finding_priority_mode all --igv_report_limit 0
```

Generate at most 100 reports for selected impacts and consequences:

```bash
sbatch scratch.slurm -- \
  --finding_priority_mode filter \
  --finding_priority_impacts HIGH,MODERATE \
  --finding_priority_consequences stop_gained,frameshift_variant \
  --igv_report_limit 100
```

Generate reports for a user-supplied gene list:

```bash
sbatch scratch.slurm -- \
  --finding_priority_mode filter \
  --finding_priority_genes GENE1,GENE2 \
  --igv_report_limit 0
```
