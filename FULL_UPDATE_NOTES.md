# PGTK grounded full update

Baseline: exact source bundle and job 19181644 artifacts supplied on 12 August 2026.

Implemented:

- Pysam 0.24.0 and HTSlib 1.23.1 runtime for BUILD_IGV_EVIDENCE_BUNDLE and BUILD_FINDING_IGV_REVIEWS.
- No external samtools subprocesses in either IGV Python script.
- ALT-specific VEP CSQ selection with PICK then CANONICAL preference.
- Deterministic consolidation of RNA-validated and progression duplicate findings.
- Preservation of EvidenceClasses, SourceEvents and Sources in strict-review outputs.
- Sorted, indexed and quickchecked evidence BAM generation.
- IGV batches terminate explicitly after snapshot generation.
- Exact job 19181644 regression fixture and validator.

Intentionally unchanged:

- SPLIT_N_CIGAR and the GATK calling branch.
- RNA splice, codon, provenance and proteogenomic validation implementations.
- Existing output layout outside the strict IGV review provenance additions.
