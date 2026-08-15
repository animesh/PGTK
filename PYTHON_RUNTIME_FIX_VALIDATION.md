# Python runtime correction

This update removes SciPy and all other optional Python-package dependencies from the GO runtime scripts.

The optimized ORA uses a standard-library logarithmic hypergeometric recurrence based on `math.lgamma`. It does not use the previous repeated `math.comb` implementation. Ranked GO uses a standard-library tied-average-rank implementation and `math.erfc` for the two-sided normal approximation.

Runtime preflight now executes the configured host Python inside both actual target container classes with `--no-home` and the runtime bind set:

- samtools 1.21 container for host-Python consumers such as the IGV bundle
- MultiQC 1.35 container for expression and progression GO scripts

The preflight executes `expression_go_analysis.py --help` and `analyze_progression_biology.py --help` inside the MultiQC image before Nextflow submits any task.

Validation performed on the supplied source and real optimization inputs:

- no SciPy, NumPy, pandas, or user-site dependency in GO runtime files
- all Python files compile with `PYTHONNOUSERSITE=1`
- GO fixtures execute using Python `-S`, disabling site-package initialization
- 1,000 random hypergeometric tables match exact combinatorial probabilities within 2e-11
- all 7,529 real GO terms match the prior SciPy outputs within 1e-8
- real 63,241-gene matrix and 1,995,057-row GO mapping: sample ORA 11 seconds; ranked GO 12 seconds
- 64 unique Nextflow processes and matching selectors remain unchanged

The exact Saga Apptainer runtime contract must pass `validate_runtime_inputs.py` on Saga before the pipeline starts. The launcher performs this automatically.
