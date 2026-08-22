# Compact nonsynonymous GO and direct-open explorer update

Fixes the quadratic GO row duplication that produced a 146 GB table.

Changes:
- GO pairs and foreground genes are deduplicated.
- Hypergeometric survival probabilities use a stable recurrence rather than repeated combinatorial sums.
- Only FDR-significant terms, top 100 terms per sample-stage, and compact summary counts are written.
- GO output row-count and file-size safeguards are enforced.
- The direct-open explorer embeds compact finding arrays only; full records remain compressed for optional server mode.
- The explorer HTML has a 40 MB hard size limit.
- MultiQC custom filenames no longer expose numeric ordering prefixes.
- Variant-landscape runtime allocation is reduced from 8 hours to 2 hours.

Follow-up fix: BUILD_FINDING_EXPLORER now validates `Embedded compact records`, matching the compact explorer coverage summary.
