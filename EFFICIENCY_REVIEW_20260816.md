# Pipeline efficiency review, 2026-08-16

## Confirmed bottleneck

The failed run exhausted the shared project file quota, not the byte quota. `results/igv/findings/finding_reviews` contained 118,338 entries, while the complete project reached 2,000,000 entries. The previous strict-review implementation created a separate directory and roughly 17 files for every finding.

## Implemented correction

The strict finding review now produces one flat consolidated bundle. Output growth is proportional to sample count rather than finding count:

- shared findings, priority, classification, exclusion, region and BAM manifests
- one shared BED track and one priority BED track
- one IGV batch and one IGV session
- four BAM categories plus indexes per sample
- one optional consolidated IGV Reports HTML file

For three samples, the expected strict-review bundle is approximately 36 entries instead of 118,338. The Nextflow process enforces a sample-scaled entry bound and fails before publishing if per-finding expansion returns.

## Runtime evidence limitations

The supplied `trace-20260815-45198746.txt` contains only its header. The current failed run could not finish or write its resource summaries because the project was already at the hard file quota. Therefore CPU efficiency and peak-memory recalibration were not changed in this update.

The wrapper log still identifies long wall-clock stages including SRA conversion, trimming, STAR alignment, MarkDuplicates and SplitNCigarReads. These are expected candidates for later resource tuning, but allocations should not be reduced without a complete Nextflow trace containing `%cpu`, `peak_rss`, `realtime`, requested CPUs and requested memory.

## Recommended next measurement

After freeing project file quota, resume the run with this update. Preserve `results/pipeline_trace-<job-id>.tsv`, then run `analyze_pipeline_trace.py` on the completed trace. Compare at least two successful runs before changing production resources.
