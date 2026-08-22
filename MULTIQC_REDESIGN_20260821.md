# Compact MultiQC redesign

Changed only the final reporting layer. Upstream analytical outputs are unchanged.

- Replaces embedded full Markdown and audit reports with 14 compact dashboard sections.
- Adds plot-driven variant attrition, RNA evidence, progression, GO, Sarek, MaxQuant, read-validation and independent-validation summaries.
- Keeps full reports and TSV audits as relative links.
- Preserves raw/trimmed and R1/R2 FastQC and Cutadapt sample identities.
- Removes llms-full.txt from the final report data directory.
- Fails final MultiQC generation if HTML exceeds 25 MB.
