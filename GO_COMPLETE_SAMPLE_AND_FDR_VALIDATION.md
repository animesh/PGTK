# Complete progression-sample GO and CLI-driven FDR validation

## Analyses generated per subject

For every subject, progression variant-set GO now emits:

- `common_all_progression`, intersection across all progression samples
- `<sample>_all`, the complete progression gene set for each progression sample
- `<sample>_exclusive`, genes found in that progression sample and absent from the subject's other progression samples

For TK13 and TK14, the output therefore contains exactly:

- `common_all_progression`
- `TK13_all`
- `TK13_exclusive`
- `TK14_all`
- `TK14_exclusive`

## SignificantGOTerms derivation

Every eligible GO term receives a raw one-sided over-representation P value. Benjamini-Hochberg correction is applied across all eligible GO terms within that analysis. `SignificantGOTerms` is then calculated as the number of rows satisfying:

```text
FDR <= configured fdr_threshold
```

The applied threshold is written to the `FDRThreshold` column in every summary row.

## Configuration

The Nextflow parameter is:

```text
--go_fdr_threshold
```

Its default is `0.1`. It is passed to all four GO analysis paths:

- per-sample expression ORA
- ranked expression GO
- complete/common/exclusive progression variant-set GO
- per-sample progression GO

The standalone Python commands expose:

```text
--fdr-threshold
```

Their default is also `0.1`. Values outside `[0, 1]` are rejected.

## Real-input validation

Using the supplied progression gene table, 63,241-gene expression matrix and 1,995,057-row propagated GO mapping, the complete progression variant-set analysis finished in 13 seconds and produced all five expected analyses.

At FDR threshold 0.1:

```text
common_all_progression: 0
TK13_all: 15
TK13_exclusive: 27
TK14_all: 5
TK14_exclusive: 5
```

With explicit CLI override `--fdr-threshold 0.05`:

```text
common_all_progression: 0
TK13_all: 9
TK13_exclusive: 14
TK14_all: 1
TK14_exclusive: 3
```

This confirms the significance counts are threshold-driven rather than hard-coded.
