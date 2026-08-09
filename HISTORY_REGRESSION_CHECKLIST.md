# Pipeline history regression checklist

Reviewed source archive SHA-256: `b5ce5c7aca9d729817ef6b5a72ecc042697c53e6082fdd1b3c4ffea474dc65ff`.

Checked historical source files, validation reports, Nextflow logs, Slurm logs and traces packaged on 2026-08-06.

Preserved:

- Validated RNA event, codon, read-provenance and proteogenomics reporting.
- Partial versus failed codon semantics.
- Seventeen-section final MultiQC content.
- GVCF gathering followed by explicit indexing.
- HaplotypeCaller interval scattering and GVCF gathering.
- Explicit Saga Python for scripts running in samtools containers.
- Nextflow resume, Apptainer cache and calibrated process resources.

Removed entirely:

- Scattered SplitNCigarReads.
- SplitNCigarReads shard gathering.
- Boundary duplicate normalization.
- SplitNCigarReads equivalence branch.
- SplitNCigarReads helper scripts from production wiring.
- Unsupported `samtools cat -@` usage.
- Python execution inside the samtools gather container.
- Artificial global `maxForks` limit.

Historical failures explicitly guarded against:

- Multi-output process object used as a data channel.
- Missing samtools inside the GATK container.
- Missing Python inside the samtools container.
- Unsupported samtools subcommand options.
- Old 47-process packages silently replacing newer scientific reporting files.
- HTML entity corruption in Nextflow source.
- Missing BAM indexes and incorrect BAM staging names.
- Missing MaxQuant sentinel mapping file.
- Incorrect partial-codon classification.
- Missing GVCF index after GatherVcfs.
- Failure logging skipped because `set -e` terminated the wrapper.


Dynamic-resource update checks:

- No `#SBATCH` settings remain in `scratch.slurm`.
- No cluster path, account, partition, Java module or executable is embedded in `scratch.slurm`.
- Wrapper options and Nextflow parameters are separated by a mandatory `--`.
- CPU and RAM are scaled 1x, 2x and 4x across three total attempts.
- Partition routing uses effective resources after retry scaling.
- Recommended routing is normal at up to 20 CPUs and 160 GB, bigmem above either threshold.
- Recommended absolute caps are 32 CPUs and 512 GB.
- Failure records include the effective partition.

Tool-resource wiring checks:

- Every GATK maximum heap is derived from 80 percent of effective `task.memory`.
- No fixed 40 GB, 18 GB, 6 GB or 4 GB GATK maximum heap remains.
- MarkDuplicates garbage-collector threads are bounded by `task.cpus`.
- HaplotypeCaller PairHMM threads remain wired to `task.cpus`.
- CPU retry scaling is limited to tools with an explicit thread-control argument.
- Single-threaded and effectively serial tools retain calibrated CPU counts while memory escalation remains active.
