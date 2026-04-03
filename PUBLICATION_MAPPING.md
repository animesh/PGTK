# Mapping of Computational Results to Published Study

> Generated: 2026-04-03  
> **⚠️ Note:** The PMID you provided (39628151) appears to be a typo. That ID resolves to an unrelated Chinese environmental science paper. The correct PMID for your study is **39682151**.

---

## 📄 Publication Reference

**Title:** Multiple Myeloma Cells with Increased Proteasomal and ER Stress Are Hypersensitive to ATX-101, an Experimental Peptide Drug Targeting PCNA

**Journal:** Cancers (Basel). 2024 Nov 26;16(23):3963  
**DOI:** https://doi.org/10.3390/cancers16233963  
**PMID:** 39682151  
**PMCID:** PMC11640687

**Authors:** Olaisen C, Røst LM, **Sharma A**, Søgaard CK, Khong T, Berg S, Jang M, Nedal A, Spencer A, Bruheim P, Otterlei M  
*(Animesh Sharma — PROMEC, NTNU — is a co-author)*

---

## 🧪 Dataset Identity

The SRA accessions in `fastq/` (SRR31089070–SRR31089076) are linked to:
- **BioProject:** PRJNA1176350
- **SRA Study:** SRP540312
- **Institution:** Norwegian University of Science and Technology (NTNU)
- **Sample names in SRA:** TK18, matching `TK18R*` BAMs in this folder

The **TK cell lines** (TK9, 10, 12, 13, 14, 16, 18) are **patient-derived multiple myeloma cell lines** with a primary cell-like phenotype, developed and studied at NTNU/Otterlei lab. This folder contains 9 groups (TK10–TK93), where the TK9x/TK9xLxxx naming convention likely represents **longitudinal/replicate samples** from the same patients.

---

## 🔬 Publication Context vs. This Analysis

The paper is a **multi-omics study** of ATX-101 (a PCNA-targeting peptide drug) in multiple myeloma. The published work used:
- **Transcriptomics** (RNA-seq) ← *the BAMs in this folder are the source data*
- **Signallomics** (subproteomics)
- **Metabolomics**
- Traditional assays (viability, Western, ELISA)

The computational work in this folder represents **a new/extended genomic analysis layer** — somatic variant calling and copy number analysis — applied to the same RNA-seq data that was used for transcriptomics in the paper.

---

## 📊 Result-by-Result Mapping

### 1. Sample Groups → Cell Lines in Publication

| This Folder | Published Cell Line | Published Role |
|-------------|--------------------|----|
| TK10 (TK1049, TK1050, TK1051) | TK10 | Patient-derived MM cell line |
| TK12 (TK12R2, TK12R3) | TK12 | Patient-derived MM cell line |
| TK13 (TK131, TK132, TK133) | TK13 | **ATX-101 hypersensitive** (sensitive group) |
| TK14 (TK141, TK142, TK143) | TK14 | **ATX-101 hypersensitive** (sensitive group) |
| TK16 (TK16R1, TK16R2, TK16R3) | TK16 | **ATX-101 hypersensitive** (sensitive) — but placed in Group A in current analysis ⚠️ |
| TK18 (TK18R1, TK18R2, TK18R3) | TK18 | Patient-derived MM cell line |
| TK91 (TK91L003, TK91L005, TK91L006) | TK9 (replicate/timepoint) | Same patient as TK9 |
| TK92 (TK92L002, TK92L003, TK92L005, TK92L006) | TK9 (replicate/timepoint) | Same patient as TK9 |
| TK93 (TK93L002, TK93L003, TK93L005, TK93L006) | TK9 (replicate/timepoint) | Same patient as TK9 |

> **Key insight from the paper:** TK13, TK14, and TK16 are the **ATX-101-sensitive** cell lines (Group B in your analysis: `TK12, TK13, TK14`). TK10, TK16, TK18, TK91, TK92, TK93 map roughly to the **less sensitive** group (Group A). The paper notes that sensitivity correlates with elevated proteasomal/ER stress, elevated ribosomal gene expression, and low NAD+/NADH.

---

### 2. Cohort Grouping in Analysis Files

Your analysis splits samples into **Group A** and **Group B**:

| Group | Cell Lines | Publication Phenotype |
|-------|-----------|----------------------|
| **Group A** | TK10, TK16, TK18, TK91, TK92, TK93 | Mixed — contains TK16 (**sensitive**) alongside less-sensitive lines |
| **Group B** | TK12, TK13, TK14 | TK13 and TK14 are **ATX-101-hypersensitive**; TK12 sensitivity unclear |

> ⚠️ **Important mismatch:** The paper defines the **sensitive** lines as **TK13, TK14, and TK16**. However in your analysis, **TK16 is placed in Group A** (together with the less-sensitive lines), while TK13 and TK14 are in Group B. This means your Group A vs Group B split does **not** cleanly separate sensitive from resistant — TK16 is on the wrong side.
>
> The grouping may reflect a different biological question (e.g. patient of origin, disease stage, or treatment timepoint) rather than ATX-101 sensitivity per se. Worth revisiting the rationale for the current A/B split in `group_compare_two_cohorts.py`.

---

### 3. Variant Calling Results (`results/mutect2_merged/`, `results/delly_final/`)

| Result File | Connection to Paper |
|-------------|---------------------|
| `TK*.merged.snv.vcf.gz` | SNV/indel landscape of the MM cell lines — extends the paper's transcriptomic/proteomic view with genomic mutation data |
| `TK*.final.sv.vcf.gz` | Structural variant landscape — not covered in the publication; this is novel additional analysis |
| MAPQ remapping (255→60) | Necessary because BAMs are STAR-aligned RNA-seq (same data used for transcriptomics in the paper) |

> The paper characterized these cell lines at the **transcriptome, signalome, and metabolome** level. The somatic variant calls here add a **genomic/DNA-level** view of the same cells that is **not in the publication**.

---

### 4. Downstream Analysis Results (`results/analysis/`)

| Result File | Maps to Paper Section | Notes |
|-------------|----------------------|-------|
| `cohort_snv___indels_groupA_exclusive.tsv` | Not in paper | SNVs exclusive to Group A (TK10/16/18/91/92/93) — novel |
| `cohort_snv___indels_groupB_exclusive.tsv` | Not in paper | SNVs exclusive to Group B (TK12/13/14) — novel |
| `cohort_snv___indels_shared.tsv` | Not in paper | Mutations shared across all groups |
| `cohort_structural_variants_(sv)_*.tsv` | Not in paper | SV landscape — entirely novel analysis |
| `gene_snv_matrix.tsv` | Relates to transcriptomic data in paper | Per-gene mutation matrix across all TK lines |
| `gene_sv_matrix.tsv` | Not in paper | Per-gene SV matrix |
| `exclusive_genes_annotated.tsv` / `exclusive_genes_summary.tsv` | Partially relates to paper | Genes with group-exclusive mutations — connects to paper's biomarker discussion |
| `snv_exclusive_genes.tsv` | Not in paper | Genes mutated only in one group |
| `sv_exclusive_genes.tsv` | Not in paper | Genes with SVs only in one group |
| `cancer_genes_cn.tsv` | Relates to paper's discussion of MM biology | CN at known cancer genes (e.g., MCL1, BCMA/TNFRS17) |
| `perchrom_mannwhitney.tsv` + `.png` | Not in paper | Statistical comparison of per-chromosome CN between groups |
| `windowed_cn_profile.png` | Not in paper | Genome-wide CN landscape |
| `coverage_ideogram.png` + `coverage_per_chrom.png` | Not in paper | Coverage/CN visualizations |

---

### 5. Biomarker Genes from Paper — Check in Your Results

The paper identifies **11 proteins** specifically activated in sensitive MM lines (TK13, 14, 16, JJN3):

| Protein (Paper) | Gene | Check in Your Exclusive Mutations |
|-----------------|------|----------------------------------|
| TPD52 | TPD52 | Check `snv_exclusive_genes.tsv` / `sv_exclusive_genes.tsv` |
| TNFRS17/BCMA | TNFRSF17 | Present in `cancer_genes_cn.tsv` as TNFRS17 (Group B CN) |
| LILRB4/ILT3 | LILRB4 | Check exclusive genes |
| TSG101 | TSG101 | Check exclusive genes |
| ZNRF2 | ZNRF2 | Check exclusive genes |
| UPF3B | UPF3B | Check exclusive genes |
| FADS2 | FADS2 | Check exclusive genes |
| C11orf38/SMAP | SMAP | Check exclusive genes |
| CGREF1 | CGREF1 | Check exclusive genes |
| GAA | GAA | Check exclusive genes |
| COG4 | COG4 | Check exclusive genes |

> **Cross-reference opportunity:** Variants or CN alterations in these 11 genes across Group A vs B could directly extend the paper's biomarker findings at the genomic level.

---

### 6. VEP Annotation (`results/vep/`)

| Result | Maps to Paper |
|--------|--------------|
| `TK10.vep.tsv` | Functional annotation of TK10 SNVs — connects to paper's pathway analysis (glycolysis, PPP, ER stress) |
| Remaining groups pending | VEP would annotate all groups' variants for functional impact |

---

### 7. Coverage / Copy Number (`results/coverage/`)

| Result | Maps to Paper |
|--------|--------------|
| `all_idxstats.tsv` | Per-chromosome read counts across all 28 samples |
| `bedcov/` | High-resolution coverage for CN estimation |
| `windows_1mb.bed` | 1 Mb windowed CN — basis for `windowed_cn_profile.png` |

> The paper does not include copy number analysis. This is **entirely new work** extending the published dataset.

---

## 🧭 Summary: What's Published vs. What's New Here

| Analysis Type | In Publication (PMID 39682151) | In This Folder |
|---------------|-------------------------------|----------------|
| RNA-seq expression (transcriptomics) | ✅ Yes — central to the paper | Source BAMs present |
| Signallomics (subproteomics) | ✅ Yes | Not reproduced here |
| Metabolomics | ✅ Yes | Not reproduced here |
| Viability / Western / ELISA | ✅ Yes | Not reproduced here |
| **Somatic SNV/Indel calling** | ❌ Not in paper | ✅ **New — this folder** |
| **Structural variant calling** | ❌ Not in paper | ✅ **New — this folder** |
| **Copy number profiling** | ❌ Not in paper | ✅ **New — this folder** |
| **Group-exclusive mutation analysis** | ❌ Not in paper | ✅ **New — this folder** |
| **Per-chromosome statistical comparison** | ❌ Not in paper | ✅ **New — this folder** |
| **VEP functional annotation** | ❌ Not in paper | 🔄 In progress (TK10 done) |

---

## 💡 Potential for Follow-up / Extended Publication

This folder represents a substantial **genomic extension** of the published multi-omics study. Key angles for a follow-up:

1. **Genomic basis of ATX-101 sensitivity** — Do SNVs or SVs in the sensitive lines (TK13, TK14, TK16) explain the hypersensitivity phenotype? Check group-exclusive mutations in ER stress, proteasome, or DNA repair genes.
2. **TNFRS17/BCMA copy number** — Already showing up in `cancer_genes_cn.tsv`. BCMA is a key MM therapeutic target.
3. **MCL1 amplification** — Present in `cancer_genes_cn.tsv` with HIGH AMP in both groups. MCL1 is a canonical MM survival gene and ATX-101 resistance-related.
4. **Longitudinal genomic evolution** — TK91/92/93 are likely serial samples from the same TK9 patient. Comparing somatic variants across these could reveal clonal evolution under treatment.
5. **Cross-validation with transcriptomics** — Genes with somatic mutations AND expression changes (from the paper) would be particularly strong candidates for functional relevance.
