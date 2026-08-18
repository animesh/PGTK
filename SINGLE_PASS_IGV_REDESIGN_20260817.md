# Single-pass IGV redesign

The strict review stage now scans each per-sample event BAM exactly once.

- Variant positions are indexed in memory and matched to each alignment by binary search.
- Fusion and splice intervals are merged and matched during the same scan.
- Output BAM records are written directly to four unsorted BAM streams and coordinate-sorted once.
- SQLite has been removed.
- Repeated per-event BAM fetches have been removed.
- Read-classification diagnostics remain bounded.
- TSV remains the portable final manifest format. Parquet is not used because tabular serialization was not the dominant cost and downstream IGV/Nextflow consumers currently expect TSV.

This changes complexity from repeated random BAM queries plus duplicated SQLite inserts to one sequential pass over each of the three event BAMs.
