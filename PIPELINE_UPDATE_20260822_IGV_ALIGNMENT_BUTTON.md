# IGV alignment-view button fix

The status-aware offline explorer now distinguishes direct-file mode from HTTP server mode.

In server mode, every candidate row contains an Open alignments button. The button calls the existing `/report/<EventID>` endpoint, generates a locus-specific standalone IGV.js report with the configured event, exact-ALT and reference BAM tracks, caches the report under `finding_explorer/report_cache/`, and opens it in a new browser tab.

In direct-file mode, the alignment control states that server mode must be started.

Validation performed:
- Python syntax validation
- JavaScript syntax validation with Node
- DDX1 status fixture validation
- exact HTTP `/report/<EventID>` route and cache validation using a deterministic create_report fixture
- full static pipeline contract, reporting redesign, resource configuration and hard-coding audits
