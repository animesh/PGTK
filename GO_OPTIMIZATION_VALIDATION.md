# Expression and GO optimization

The sequential `ANALYZE_EXPRESSION_GO` task was replaced by independent Nextflow tasks:

- `ANALYZE_EXPRESSION_SAMPLE_GO`, one task per sample
- `ANALYZE_EXPRESSION_RANKED_GO`, one task per valid progression and baseline comparison
- `MERGE_EXPRESSION_GO`, one deterministic merge task
- `PREPARE_EXPRESSION_MULTIQC_CONTENT`, one reporting task

The GO over-representation test now uses `scipy.stats.hypergeom.sf` instead of repeated `math.comb` calculations. Ranked GO uses tied average ranks from `scipy.stats.rankdata` and a two-sided rank-sum normal approximation. The runtime preflight imports SciPy using the configured `--host_python` inside the Apptainer execution contract before Nextflow submits tasks.

Using the supplied real matrix containing 63,241 genes and the propagated mapping containing 1,995,057 gene-term rows:

- One sample ORA completed in 17 seconds.
- One ranked progression-versus-baseline GO test completed in 14 seconds.
- The cancelled sequential implementation was still running after 4 hours 48 minutes.

The benchmark was executed against the supplied real input files in the optimization archive. The output contained 7,529 tested GO terms per analysis.

Subjects with no baseline or multiple baselines are explicitly reported as skipped in the merged expression GO summary. No experiment-specific sample names, subjects, groups, genes, pathways, or category lists are embedded in the implementation.
