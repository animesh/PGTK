# Ranked GO validation

The supplied current source already used the progression sample TPM in the numerator and matched baseline TPM in the denominator:

```text
log2((sample_TPM + pseudocount) / (baseline_TPM + pseudocount))
```

The uploaded completed results were therefore not zero-score artifacts. The real matrix validation reproduced every MeanScore, ZScore, PValue and FDR value in the supplied ranked GO output exactly.

This update hardens the implementation by:

- resolving sample and baseline TPM columns explicitly
- failing on missing or duplicate matrix columns
- rejecting self-comparisons
- rejecting zero-information comparisons through `--min-nonzero-scores`
- wiring `--expression_rank_min_nonzero_scores` through Nextflow, default 1
- recording rank metric, pseudocount, non-zero, positive and negative score counts, and score range in summary output
- adding tests that require negative pathway scores for baseline-enriched genes and positive scores for progression-enriched genes

## Real matrix validation

```text
TK13 versus TK12
Genes ranked: 61617
Non-zero scores: 36351
Positive scores: 21337
Negative scores: 15014
Minimum score: -6.787411816340121
Maximum score: 6.4290331853101765
Significant GO terms at FDR 0.1: 3859

TK14 versus TK12
Genes ranked: 61617
Non-zero scores: 35055
Positive scores: 20906
Negative scores: 14149
Minimum score: -5.955449844767238
Maximum score: 6.685715610848317
Significant GO terms at FDR 0.1: 5154
```

The GO FDR default remains CLI-driven at 0.1. Complete, common and exclusive progression variant GO analyses are preserved.
