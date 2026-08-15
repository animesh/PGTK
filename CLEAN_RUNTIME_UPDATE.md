# Clean container-native runtime update

This revision removes the redundant global Apptainer binding architecture.

Changes:

- Removed `APPTAINER_BINDPATH` and the global runtime bind list from `scratch.slurm`.
- Removed the `host_python` Nextflow parameter and stopped passing the Saga Python executable into containers.
- Python-only processes use each container's native `python3`.
- BAM, FASTA and HTSlib-aware Python processes use the pinned Pysam 0.24.0 container.
- `validate_rna_events.py`, `validate_variant_codons.py`, `validate_variant_read_provenance.py`, and `validate_proteogenomic_reads.py` no longer launch external samtools subprocesses. They use Pysam and bundled HTSlib APIs.
- Standalone samtools operations remain in the pinned samtools 1.21 container.
- IGV evidence and strict finding review remain in the pinned Pysam container.
- IGV Reports remain in the pinned IGV Reports 1.16.0 container.
- Runtime preflight validates native Python in the exact MultiQC and Pysam containers without manual bind arguments.
- Slurm stdout now defaults to `pgtk-wrapper-<job-id>.log` in the submission directory.
- GO ontology and GAF paths are derived from `--reference-downloads` automatically by the wrapper.

Expected result on Saga: no duplicate bind-mount warnings issued by this pipeline.
