# Validation and resource semantics

## Validation outcomes

`VARIANT_CODON_VALIDATED` requires genome-reference agreement, sample-matched RNA ALT support, and a directly comparable codon consequence. `VARIANT_EVIDENCE_PARTIAL` retains genome and RNA evidence when the coding consequence is not safely comparable as one codon pair. `VALIDATION_FAILED` requires at least one explicit failure code.

Synonymous VEP consequences pass as `SYNONYMOUS_CODON_TRANSLATION_CONFIRMED` when reference and alternate codons translate to the same amino acid. An empty VEP alternate-amino-acid field is accepted for this synonymous representation.

Strict integrated evidence requires `Codon strict-integration eligible = yes`, sample-matched direct MS/MS, a search-consistent altered-residue peptide, and absence from both canonical reference sets.

## Resource profile

`robust` is the only configured Nextflow profile and is the production default. `scratch.slurm` selects it unless `PGTK_PROFILE` or `--profile` explicitly supplies another configured profile.

The `robust` profile applies process-specific CPU, memory, and time directives for Saga. Retry attempts increase resources for configured processes within the limits set by `PGTK_MAX_CPUS` and `PGTK_MAX_MEMORY_GB`. Retry eligibility and exit-status handling are defined by `nextflow.config`; this document does not define additional `calibrated`, `conservative`, or `recovery` profiles.

The Slurm wrapper writes Nextflow trace, report, timeline, DAG, failure-ledger, and post-run resource summaries under the configured results directory.
