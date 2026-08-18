# Saga resource calibration

The calibrated profile was derived from 171 successful task records in job 19002083.

Key changes:

- HaplotypeCaller: 8 CPUs and 20 GB to 4 CPUs and 6 GB per shard.
- STAR alignment: 20 CPUs and 128 GB to 32 CPUs and 64 GB, correcting observed CPU oversubscription.
- Trim Galore: retains 8 allocated CPUs but uses half as its internal core setting to avoid child-process oversubscription.
- SRA conversion: 12 CPUs and 32 GB to 4 CPUs and 4 GB.
- StringTie: 8 CPUs and 32 GB to 2 CPUs and 6 GB.
- Arriba: 8 CPUs and 32 GB to 1 CPU and 24 GB.
- Sort/index BAM: 12 CPUs and 48 GB to 8 CPUs and 24 GB.
- Memory-heavy GATK and read-validation processes retain conservative headroom.

Do not reduce `MARK_DUPLICATES`, `SPLIT_N_CIGAR`, or `VALIDATE_PROTEOGENOMIC_READS` further until at least two complete calibrated runs have been compared.

## Robust production retry policy

The `robust` profile permits three total attempts per failed task. Attempt 1 uses the calibrated allocation, attempt 2 uses twice the CPU and RAM, and attempt 3 uses four times the initial CPU and RAM. Allocations are capped at 64 CPUs and 256 GB. Runtime limits are not silently changed. Every failed attempt is preserved under `results/failure_logs/<job-id>/`, including trace metadata, work directory, command, stderr, stdout and exit code when available.

There is no `maxForks` limit. Nextflow may submit up to 200 dependency-ready tasks to Slurm. Biological dependencies still determine which tasks are eligible to run concurrently.
### Dynamic production retry and partition policy

The robust profile permits three total attempts per task. Attempt 1 uses calibrated resources, attempt 2 uses twice the initial CPU and RAM, and attempt 3 uses four times the initial CPU and RAM. Effective requests are capped at the command-line supplied limits. The recommended Saga values are 32 CPUs and 512 GB.

Partition selection is recalculated after scaling for every attempt. With the recommended Saga thresholds, requests above 20 CPUs or above 160 GB use the big-memory partition; all others use the normal partition. Partition names, thresholds, caps, account, queue size and submission rate are supplied through `scratch.slurm` command-line options and exported to `nextflow.config`.

Every failed attempt is preserved under `results/failure_logs/<job-id>/`, including partition, effective resources, work directory, command, stderr, stdout and exit code when available.

### Tool-level retry resource wiring

Memory escalation is consumed by every GATK command through a Java maximum heap equal to 80 percent of the effective Nextflow task memory. The remaining 20 percent is reserved for native libraries, compression buffers, temporary structures and operating-system overhead. MarkDuplicates garbage-collector threads are bounded by the effective task CPU allocation.

CPU escalation is enabled only for commands that explicitly consume `task.cpus`, including fasterq-dump, FastQC, Trim Galore, STAR, samtools sorting and indexing, HaplotypeCaller PairHMM, VEP, StringTie and threaded validation scripts. Processes without a supported tool-level thread control keep their calibrated CPU count while RAM still increases on retries. This prevents requesting idle CPUs merely because a task was retried.

## Completed run 19199689 calibration

`VALIDATE_VARIANT_READ_PROVENANCE` used 9.1 GB peak RSS for TK13 after the 8 GB attempt exited 137. The portable profile now keeps CPU fixed at 2 and requests 12 GB, 24 GB, and 48 GB across attempts. CPU is not escalated because the completed attempt used less than one effective core.

IGV report preparation now caps displayed ALT reads separately from the complete exact-ALT BAM. Full evidence BAMs remain published, while HTML reports use capped display BAMs to control embedded-report size.
