# Portable resource and IGV update, 2026-08-15

This release is based on the pipeline that completed Saga job 19199689.

Changes:

- `VALIDATE_VARIANT_READ_PROVENANCE` now uses 2 CPUs on all attempts and 12, 24, and 48 GB memory.
- Removed fixed Saga account, partitions, user paths, Mamba Python, Nextflow path, Java module, work path, and temporary path.
- Wrapper settings are supplied by command-line options or matching `PGTK_*` environment variables.
- The wrapper discovers `python3`, `nextflow`, and `apptainer` from `PATH` unless explicitly supplied.
- SLURM account and partitions are mandatory runtime settings.
- IGV HTML reports support RNA variants, progression variants, fusions, and splice junctions.
- Added report class, gene, sample, limit, displayed-read, timeout, and maximum-file-size controls.
- Full exact-ALT BAMs are retained. HTML reports use capped display BAMs.
- Report manifests record class, byte size, and generated or skipped status.

For SLURM, the wrapper allocation account and partition must be supplied to `sbatch`, while the same account and process partitions are supplied to the wrapper for Nextflow child jobs.
