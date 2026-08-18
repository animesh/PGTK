# Complete finding explorer redesign

The monolithic IGV Reports process was replaced because embedding 157,482 loci into one HTML file cannot scale. No finding threshold is used.

The updated pipeline:

- indexes every finding in SQLite for fast search and filtering
- maps every finding to a biological sample, evidence-class, and chromosome partition
- verifies that indexed and mapped counts equal the full manifest count and that discarded findings are zero
- provides a local searchable web interface
- generates a standalone offline IGV.js report only when a finding is selected
- caches generated per-finding reports for repeat access
- retains all full TSV, BED, BAM, BAI, IGV session, and batch outputs
- removes the monolithic report timeout and report-size threshold from the execution path

The `igv_report_limit` parameter is no longer used to remove findings. `BUILD_FINDING_IGV_REVIEWS` always receives `--priority-limit 0`.
