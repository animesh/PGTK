# PGTK - ProteoGenomics of Myeloma Cell Lines

## Overview

PGTK (ProteoGenomics TK) is a comprehensive multi-omics study of seven in-house Myeloma patient-derived cell lines (TK) with a more primary cell-like phenotype. This repository contains information and resources related to the integrated proteogenomics analysis of these cell lines.

## Publication

The detailed methodology and findings of this study are published in:

**Article:** [Cancers - MDPI](https://www.mdpi.com/2072-6694/16/23/3963)

Please cite this publication when using data or results from this project.

## Data Availability

### Transcriptome Data

The RNA-seq transcriptome data for the seven Myeloma cell lines is publicly available from NCBI BioProject:

- **NCBI BioProject ID:** [PRJNA1176350](https://www.ncbi.nlm.nih.gov/bioproject/1176350)
- **Access:** [https://www.ncbi.nlm.nih.gov/bioproject/1176350](https://www.ncbi.nlm.nih.gov/bioproject/1176350)

#### trancriptome https://www.ncbi.nlm.nih.gov/bioproject/1176350 download using https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit 

```
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.2.1/sratoolkit.3.2.1-ubuntu64.tar.gz
tar xvzf sratoolkit.3.2.1-ubuntu64.tar.gz
SRATOOL_DIR="$PWD$/sratoolkit.3.2.1-ubuntu64"; mkdir -p ./fasterq_tmp; curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=PRJNA1176350&usehistory=y" | grep -oP '(?<=<QueryKey>)\d+|(?<=<WebEnv>)[^<]+' | { read qk; read we; curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&query_key=$qk&WebEnv=$we&rettype=runinfo&retmode=text"; } | tail -n +2 | cut -d',' -f1 | xargs -P4 -I{} "$SRATOOL_DIR/bin/fasterq-dump" {} --split-files -e 4 -p -O . --temp ./fasterq_tmp
```


### Proteomics Data

The mass spectrometry proteomics data is available through the PRIDE Archive:

- **PRIDE Project 1:** [PXD033531](https://www.ebi.ac.uk/pride/archive/projects/PXD033531)
  - Access: [https://www.ebi.ac.uk/pride/archive/projects/PXD033531](https://www.ebi.ac.uk/pride/archive/projects/PXD033531)

- **PRIDE Project 2:** [PXD033510](https://www.ebi.ac.uk/pride/archive/projects/PXD033510/)
  - Access: [https://www.ebi.ac.uk/pride/archive/projects/PXD033510/](https://www.ebi.ac.uk/pride/archive/projects/PXD033510/)

  #### proteome https://www.ebi.ac.uk/pride/archive/projects/PXD033531 and https://www.ebi.ac.uk/pride/archive/projects/PXD033510/ dowload recursively from the ftp

```
wget -r "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2024/11/PXD033531/"
wget -r "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2024/11/PXD033510/"
```


## Project Description

This study focuses on seven patient-derived Myeloma cell lines (TK) that exhibit a more primary cell-like phenotype, making them valuable models for studying Multiple Myeloma biology. The integrated proteogenomics approach combines:

- **Transcriptomics:** Gene expression profiling through RNA-seq
- **Proteomics:** Protein abundance measurements through mass spectrometry

This multi-omics integration enables comprehensive characterization of the molecular landscape of these Myeloma cell lines.

## Citation

If you use this data or findings in your research, please cite:

```
[Citation details from the MDPI Cancers publication]
DOI: Available at https://www.mdpi.com/2072-6694/16/23/3963
```

---

## Nextflow Somatic Variant Calling Pipeline (`main.nf`)

### Overview

A DSL2 Nextflow pipeline for somatic variant calling from RNA-seq BAM files across multiple Myeloma cell line replicates. It calls SNVs/Indels with GATK Mutect2 and structural variants with DELLY2, then merges per-replicate calls into per-group VCFs.

### Pipeline steps

| Step | Tool | Notes |
|------|------|-------|
| FASTQ trimming | Trim Galore | Paired-end, quality ≥ 20, adapter auto-detection via `--paired` |
| Genome indexing | STAR 2.7.11b | Generated once, reused across samples |
| Alignment | STAR 2-pass | RNA-seq, outputs coordinate-sorted BAM |
| Mark duplicates | SAMtools | Removes PCR/optical duplicates |
| Reference indexing | SAMtools faidx + GATK CreateSequenceDictionary | One-time, shared across processes |
| MAPQ remap | SAMtools | STAR sets MAPQ=255 for unique reads; remapped to 60 so GATK doesn't drop them |
| SNV/Indel calling | GATK Mutect2 | Tumor-only mode; collects F1R2 read orientation counts |
| Orientation bias | GATK LearnReadOrientationModel | Corrects FFPE/oxidative damage artefacts |
| Variant filtering | GATK FilterMutectCalls | Applies orientation priors + default Mutect2 filters |
| Per-group SNV merge | BCFtools concat | Merges per-replicate filtered VCFs into one per group |
| SV calling | DELLY2 | Tumor-only; calls DEL, DUP, INV, BND |
| SV genotyping | DELLY2 genotype | Per-sample re-genotyping against merged site list |
| SV merge | DELLY2 merge + BCFtools | Merges per-replicate SV VCFs per group |

### Containers

All processes use stable `quay.io/biocontainers` images:

| Image | Used by |
|-------|---------|
| `quay.io/biocontainers/star:2.7.11b--h43eeafb_1` | STAR_GENOMEGENERATE, STAR_ALIGN |
| `quay.io/biocontainers/samtools:1.21--h50ea8bc_0` | SAMTOOLS_FAIDX, SAMTOOLS_MARKDUP, REMAP_MAPQ |
| `quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0` | GATK_DICT, MUTECT2, LEARN_ORIENTATION, FILTER_MUTECT2, MERGE_SNV_GROUP (header fix) |
| `quay.io/biocontainers/bcftools:1.21--h8b25389_0` | MERGE_SNV_GROUP, MERGE_SV_GROUP |
| `quay.io/biocontainers/delly:1.2.9--hd63ebec_1` | DELLY_CALL, DELLY_GENOTYPE, DELLY_MERGE |
| `wave.seqera.io/wt/…/trim-galore` | TRIM_GALORE |

> **Note:** `community.wave.seqera.io` images were replaced with `quay.io/biocontainers` equivalents after manifest resolution failures.

### Usage

```bash
# First run
nextflow run main.nf

# Resume after interruption / fix
nextflow run main.nf -resume
```

### Key design decisions

- **MAPQ remapping:** STAR uniquely-mapped reads carry MAPQ=255, which GATK treats as "unavailable" and silently drops. A SAMtools/awk one-liner remaps 255→60 before Mutect2 so all uniquely-mapped reads are used.
- **Tumor-only mode:** No matched normal available; Mutect2 runs in tumor-only mode with orientation bias correction to reduce artefacts.
- **Per-replicate → per-group:** Each replicate is called independently, then merged with BCFtools (SNVs) or DELLY merge (SVs) per biological group.
- **`--paired` is sufficient for adapter trimming:** Trim Galore auto-detects adapters when `--paired` is set; the now-removed `--detect_adapter_for_pe` flag was redundant and unsupported in current Trim Galore versions.

### Input samplesheet (`samplesheet.csv`)

```
sample,group,read1,read2
SRR31089070,TK12,fastq/SRR31089070_1.fastq.gz,fastq/SRR31089070_2.fastq.gz
...
```

### Output structure

```
results/
├── mutect2/<group>/          # Per-replicate unfiltered + filtered VCFs
├── mutect2_merged/           # Per-group merged SNV VCFs
├── delly/                    # Per-replicate SV calls
├── delly_final/              # Per-group merged SV VCFs
├── coverage/                 # Per-sample coverage BED files
├── vep/                      # VEP-annotated VCFs
├── analysis/                 # Downstream plots and tables
└── pipeline_info/            # Nextflow execution report, timeline, trace
```

---

## License

See [LICENSE](LICENSE) file for details.

## Contact

me 🤓