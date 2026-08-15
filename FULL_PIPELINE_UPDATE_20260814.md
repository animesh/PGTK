# Full PGTK pipeline update 2026-08-14

Merge order:

1. GitHub repository snapshot (`PGTK-main.zip`)
2. Saga production snapshot (`pgtk.tgz`) as the operational source of truth
3. Phase 2 generic CLI update (`pgtk.generiCLIpipeline.zip`)

The Saga production implementation is preserved. Phase 2 changes add generic finding selection, consolidated per-finding read review, and optional self-contained IGV Reports without default gene prioritization.

Primary Phase 2 files:

- `main.nf`
- `build_finding_igv_reviews.py`
- `test_resource_configuration.py`
- `PHASE2_CLI_PARAMETERS.md`
- `validate_full_update.slurm`
- `validate_pipeline_commands.sh`

Default finding review behavior:

- classes: `rna_variant,progression_variant`
- priority mode: `all`
- gene, impact, and consequence filters: empty
- IGV report limit: `0`, meaning all selected consolidated findings

Before production submission, run the included login-node validation command documented in `RUN_AND_REPORT.md`.
