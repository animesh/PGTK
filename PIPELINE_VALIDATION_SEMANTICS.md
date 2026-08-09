# Validation and resource semantics

## Validation outcomes

`VARIANT_CODON_VALIDATED` requires genome-reference agreement, sample-matched RNA ALT support, and a directly comparable codon consequence. `VARIANT_EVIDENCE_PARTIAL` retains genome and RNA evidence when the coding consequence is not safely comparable as one codon pair. `VALIDATION_FAILED` requires at least one explicit failure code.

Synonymous VEP consequences pass as `SYNONYMOUS_CODON_TRANSLATION_CONFIRMED` when reference and alternate codons translate to the same amino acid. An empty VEP alternate-amino-acid field is accepted for this synonymous representation.

Strict integrated evidence requires `Codon strict-integration eligible = yes`, sample-matched direct MS/MS, a search-consistent altered-residue peptide, and absence from both canonical reference sets.

## Resource profiles

`calibrated` is the production default and is based on successful Saga job 19002083. `conservative` uses the larger directives embedded in `main.nf`. `recovery` retries only exit codes 137, 140, and 143 and escalates resources for selected heavy processes.

The Slurm wrapper writes Nextflow trace, report, timeline, DAG, and post-run resource summaries under `results/`.
