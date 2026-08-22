# Portable explorer, variant landscape and nonsynonymous GO update

Implemented from the complete source package dated 2026-08-21 21:02:35.

Changes:
- Finding explorer `index.html` is self-contained and works through `file:///` without fetch calls or an HTTP server.
- Direct-open mode embeds all finding metadata and retained diagnostic read observations, with sample, class, impact, chromosome and text filters.
- Optional server mode and IGV Reports generation remain available through `serve_explorer.sh`.
- Added per-sample variant landscape reporting for raw genotyped, normalized, all hard-filtered, PASS, VEP, RNA-validated, nonbaseline-only, baseline-only and shared-with-baseline stages.
- Added allele-class, Ti/Tv component, genotype, FILTER, VEP impact and VEP consequence counts.
- Added GO over-representation analysis restricted to genes with protein-altering VEP consequences.
- Added complete TSV, Markdown and MultiQC outputs under `results/variant_landscape`.
- Replaced the implementation-facing dashboard label with `PGTK Results Guide and Navigation` and added evidence/file explanations.

Interpretation constraints:
- Raw and normalized calls are pre-filter candidates.
- RNA-validated and progression calls are RNA evidence, not DNA-confirmed somatic mutations.
- GO analysis tests gene-set over-representation and does not establish pathway activation.
