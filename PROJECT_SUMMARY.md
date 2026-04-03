# TK Project Summary

> Generated: 2026-04-03

---

## 🧬 Project Overview: Somatic Variant Calling in RNA-seq Data

This is a **somatic variant calling + downstream analysis project** on a set of cancer samples (the "TK" cohort), applied to **STAR-aligned RNA-seq BAMs** rather than DNA.

---

## 📁 The Samples

**28 tumor samples** organized into **9 groups**, each with 2–4 replicates:

| Group | Samples |
|-------|---------|
| TK10  | TK1049, TK1050, TK1051 |
| TK12  | TK12R2, TK12R3 |
| TK13  | TK131, TK132, TK133 |
| TK14  | TK141, TK142, TK143 |
| TK16  | TK16R1, TK16R2, TK16R3 |
| TK18  | TK18R1, TK18R2, TK18R3 |
| TK91  | TK91L003, TK91L005, TK91L006 |
| TK92  | TK92L002, TK92L003, TK92L005, TK92L006 |
| TK93  | TK93L002, TK93L003, TK93L005, TK93L006 |

All samples are pre-aligned, mark-duplicated, sorted BAMs (`*.markdup.sorted.bam`).

There is also a `fastq/` folder with 7 paired-end FASTQ files downloaded from SRA (SRR31089070–SRR31089076).

---

## ⚙️ The Nextflow Pipeline (`main.nf`)

A custom DSL2 pipeline for **tumor-only somatic variant calling** with two parallel arms:

| Arm | Tool | Purpose |
|-----|------|---------|
| SNV/Indel | **GATK Mutect2** | Tumor-only calling with F1R2 orientation bias correction |
| SV | **DELLY2** | Structural variant calling + genotyping |

### Key Pipeline Steps

1. **SAMTOOLS_FAIDX** — Index the reference genome (once)
2. **GATK_DICT** — Create sequence dictionary
3. **REMAP_MAPQ** — Remap MAPQ 255 → 60 *(RNA-seq specific fix: STAR sets MAPQ=255 for uniquely-mapped reads; GATK treats 255 as "unavailable" and drops them — remapping to 60 makes them visible to Mutect2)*
4. **MUTECT2** — Per-replicate tumor-only SNV/indel calling
5. **LEARN_ORIENTATION** — Learn F1R2 read orientation model
6. **FILTER_MUTECT2** — Apply Mutect2 filters with orientation priors
7. **MERGE_SNV_GROUP** — Merge per-replicate VCFs into per-group VCFs (bcftools)
8. **DELLY_CALL** — Per-replicate SV calling
9. **DELLY_MERGE / DELLY_GENOTYPE / DELLY_FILTER** — Merge, genotype, and filter SVs per group

### Infrastructure

- **Executor**: Local (12 CPUs, 78 GB RAM, `queueSize = 1` — one job at a time)
- **Containers**: Docker via Wave (`wave.seqera.io`)
- **Reference**: `genome.fa` (GRCh38), `Homo_sapiens.GRCh38.110.gtf`
- **Error strategy**: `retry` with `maxRetries = 2`
- **Strict syntax**: enabled (`nextflow.enable.strict = true`)

---

## 🏃 Pipeline Run History

**13 runs total** between March 25–31, 2026. All using `nextflow run main.nf`.

| Date | Run Name | Duration | Status | Notes |
|------|----------|----------|--------|-------|
| 2026-03-25 12:57 | `boring_austin` | 3.5s | ✅ OK | `-preview` dry run |
| 2026-03-25 13:00 | `distracted_agnesi` | 12m 20s | ❌ ERR | First full run attempt, failed |
| 2026-03-25 13:13 | `cranky_wing` | — | — | Crashed at startup (`-resume`) |
| 2026-03-25 13:30 | `hopeful_swirles` | — | — | Crashed at startup (`-resume`) |
| 2026-03-25 13:37 | `serene_nightingale` | — | — | Crashed at startup (`-resume`) |
| 2026-03-25 13:37 | `ridiculous_woese` | — | — | Crashed at startup (`-resume`) |
| 2026-03-25 14:12 | `grave_dalembert` | **19h 2m 8s** | ✅ OK | **Main successful run — full Mutect2 across all samples** |
| 2026-03-26 09:14 | `suspicious_meitner` | 5.9s | ✅ OK | Quick cache check resume |
| 2026-03-26 09:14 | `tiny_leakey` | — | — | Crashed at startup |
| 2026-03-26 15:51 | `angry_mcnulty` | 8h 27m 52s | ❌ ERR | Long DELLY run, errored |
| 2026-03-27 00:19 | `sad_knuth` | **3d 17h 15m** | ❌ ERR | Very long run (likely DELLY on large BAMs), errored |
| 2026-03-31 14:57 | `grave_kowalevski` | 1m 48s | ❌ ERR | Quick attempt, failed |
| 2026-03-31 15:04 | `silly_curran` | 5m 34s | ✅ OK | **Final successful run** ✅ |

### Arc Summary

- **March 25 afternoon**: Several startup crashes while debugging, then `grave_dalembert` ran for 19 hours and completed the bulk of Mutect2 work.
- **March 26–27**: Long DELLY runs (`angry_mcnulty`, `sad_knuth`) — both errored after many hours. Likely hitting issues with large RNA-seq BAMs in SV calling.
- **March 31**: Quick fix + final successful run (`silly_curran`) completed remaining steps.

---

## 📊 Results

### Variant Calls

| Directory | Contents |
|-----------|---------|
| `results/mutect2/` | Per-replicate raw + filtered Mutect2 VCFs |
| `results/mutect2_merged/` | Per-group merged SNV VCFs for all 9 groups (TK10–TK93) |
| `results/delly/` | Per-replicate DELLY SV calls |
| `results/delly_final/` | Per-group final SV VCFs for all 9 groups |
| `results/vep/` | VEP annotation (TK10 completed so far) |

### Downstream Analysis (`results/analysis/`)

| File | Description |
|------|-------------|
| `cohort_snv___indels_groupA_exclusive.tsv` | SNVs/indels exclusive to group A |
| `cohort_snv___indels_groupB_exclusive.tsv` | SNVs/indels exclusive to group B |
| `cohort_snv___indels_shared.tsv` | SNVs/indels shared across groups |
| `cohort_structural_variants_(sv)_groupA_exclusive.tsv` | SVs exclusive to group A |
| `cohort_structural_variants_(sv)_groupB_exclusive.tsv` | SVs exclusive to group B |
| `cohort_structural_variants_(sv)_shared.tsv` | SVs shared across groups |
| `gene_snv_matrix.tsv` | Gene × sample SNV mutation matrix |
| `gene_sv_matrix.tsv` | Gene × sample SV matrix |
| `exclusive_genes_annotated.tsv` | Annotated genes with exclusive mutations |
| `exclusive_genes_summary.tsv` | Summary of group-exclusive genes |
| `snv_exclusive_genes.tsv` | Genes with exclusive SNVs |
| `sv_exclusive_genes.tsv` | Genes with exclusive SVs |
| `cancer_genes_cn.tsv` | Copy number at known cancer genes |
| `genes_with_cn.tsv` | All genes with CN estimates |
| `gene_annotation_table.tsv` | Full gene annotation table |
| `perchrom_combined_stats.tsv` | Per-chromosome combined statistics |
| `perchrom_mannwhitney.tsv` | Mann-Whitney test results per chromosome |
| `perchrom_mannwhitney_idxstats.tsv` | idxstats-based per-chromosome stats |
| `perchrom_median_log2r.tsv` | Median log2 ratio per chromosome |
| `coverage_ideogram.png` | Genome-wide coverage ideogram |
| `coverage_per_chrom.png` | Per-chromosome coverage plot |
| `windowed_cn_profile.png` | Windowed copy number profile |
| `perchrom_mannwhitney_plot.png` | Per-chromosome Mann-Whitney plot |

### Coverage (`results/coverage/`)

| File | Description |
|------|-------------|
| `all_idxstats.tsv` | Merged samtools idxstats across all samples |
| `windows_1mb.bed` | 1 Mb genomic windows BED file |
| `bedcov/` | Per-window coverage (bedtools coverage) |

---

## 🐍 Python Analysis Scripts

Post-processing scripts run outside the Nextflow pipeline:

| Script | Purpose |
|--------|---------|
| `annotate_and_pattern.py` | VEP annotation loading + variant pattern analysis |
| `group_compare_two_cohorts.py` | Compare two cohorts (group A vs group B) |
| `group_exclusive_mutations.py` | Identify mutations unique to each group |
| `group_exclusive_mutations.sh` | Shell wrapper for exclusive mutation analysis |
| `mannwhitney_perchrom.py` | Per-chromosome Mann-Whitney U tests between groups |
| `plot_coverage_ideogram.py` | Genome-wide coverage ideogram plot |
| `plot_coverage_per_chrom.py` | Per-chromosome coverage visualization |
| `plot_windowed_cn.py` | Windowed copy number profile plotting |
| `run_bedcov_parallel.sh` | Parallel bedtools coverage across samples |
| `check_bedcov_progress.sh` | Monitor bedcov run progress |
| `scripts/annotate_vep.sh` | Shell script to run VEP annotation |

---

## 🗂️ Other Files

| File | Description |
|------|-------------|
| `genome.fa` / `genome.fa.fai` | GRCh38 reference genome + index |
| `Homo_sapiens.GRCh38.110.gtf` | Ensembl GRCh38 v110 GTF annotation |
| `samplesheet.csv` | Nextflow pipeline input (sample, group, bam, bai) |
| `nextflow.config` | Pipeline configuration (resources, Docker, profiles) |
| `combined.sorted.bam` | Combined/merged BAM file |
| `callparam.xml` / `extractionSummary.xml` | Possibly from upstream data extraction |
| `STDOUT_warnings.txt` | Captured stdout/warnings (~1.2 MB) |

---

## 📌 Current Status (as of 2026-04-03)

- ✅ **Variant calling complete** — All 9 groups have SNV and SV VCFs
- ✅ **Downstream analysis complete** — Full cohort comparison, gene matrices, CN profiles, statistical tests, and visualizations generated
- 🔄 **VEP annotation in progress** — TK10 annotated; remaining groups pending
- 🔄 **Coverage analysis ongoing** — Most recent file activity in `results/analysis/` and `results/coverage/`
