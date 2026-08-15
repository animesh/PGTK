# Run and report

## 1. Extract without replacing the active run directory

```bash
mkdir -p /cluster/projects/nn9036k/scrbkup/pgtk_update_20260814
tar -xzf pgtk_full_pipeline_update_20260814.tgz -C /cluster/projects/nn9036k/scrbkup/pgtk_update_20260814 --strip-components=1
cd /cluster/projects/nn9036k/scrbkup/pgtk_update_20260814
```

## 2. Preserve resume state only when intentionally updating the current run

Do not copy this package over an active Nextflow run. Keep the existing `work/`, `.nextflow/`, and `.nextflow.log*` with the active project when resume compatibility is required.

## 3. Validate on the Saga login node

```bash
bash validate_pipeline_commands.sh \
  --project-dir "$PWD" \
  --nextflow "$HOME/bin/nextflow"
```

Reference and container downloads must be performed directly on the login node, not through SLURM:

```bash
bash download_assets.sh
```

## 4. Submit the pipeline

```bash
sbatch scratch.slurm
```

The default Phase 2 review includes all consolidated RNA and progression variants and generates priority IGV Reports with no biological gene filter.

## 5. Optional filtered review

```bash
sbatch scratch.slurm -- \
  --finding_priority_mode filter \
  --finding_priority_impacts HIGH,MODERATE \
  --finding_priority_consequences stop_gained,frameshift_variant \
  --igv_report_limit 100
```

## 6. Report back

Return these files:

```text
results/pipeline_trace-<job-id>.tsv
results/resource_usage-<job-id>.summary.tsv
results/resource_usage-<job-id>.warnings.tsv
results/failure_logs/<job-id>/failure_ledger.tsv
results/igv/finding_reviews/consolidation_summary.txt
results/igv/finding_reviews/findings_manifest.tsv
results/multiqc/multiqc_report.html
.nextflow.log
slurm-<job-id>.out
```
