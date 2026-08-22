# Explorer validation status and visible variant GO update

Explorer changes:
- Separates upstream candidate presence from strict read-validation status.
- Defaults to ALT-supported candidates while preserving an All candidates view.
- Shows ALT-supporting, reference-supporting, excluded/uncallable, total examined, callable and ALT-fraction columns.
- Adds ALT_SUPPORTED, MIXED_ALT_AND_REFERENCE, REFERENCE_ONLY, NO_CALLABLE_READS and NO_OVERLAPPING_READS statuses.
- Aggregates retained exclusion reasons from finding_reviews/excluded_reads.tsv and explicitly reports when reason rows were unavailable or capped.

MultiQC changes:
- Replaces the link-only nonsynonymous GO section with a visible RNA-seq Protein-Altering Variant GO summary.
- Adds top-term plots for each RNA-validated sample and each progression nonbaseline-only sample.
- Reports protein-altering gene, overlapping GO-term and significant GO-term counts.
- Explains the distinction between expression GO, progression-set GO and protein-altering variant GO.
- Orders section assets immediately after the variant-type asset using hidden filename ordering while keeping human-readable section titles.
