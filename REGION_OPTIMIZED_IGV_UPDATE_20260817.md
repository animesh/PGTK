# Region-optimized consolidated IGV update

The strict IGV review no longer performs read classification for every fusion and splice-junction event. Those non-allelic findings are converted to padded intervals, overlapping intervals are merged per sample, and each merged interval is scanned once. Variant and progression-variant findings retain allele-specific read classification.

Output bounds:

- `read_classification.tsv` contains at most 100,000 non-excluded diagnostic rows by default.
- `excluded_reads.tsv` contains at most 10,000 excluded diagnostic rows by default.
- fusion and splice-junction reads are written directly to one per-sample event-display BAM rather than duplicated into event-read TSV matrices or SQLite.
- exact-alt and clean-reference variant display BAMs remain globally deduplicated.
- progress is emitted every 1,000 findings and after every sample interval scan.

The complete findings manifest, event consolidation mapping, event regions, BED tracks, IGV session, IGV batch file, and priority manifest remain available.
