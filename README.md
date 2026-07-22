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

## Nextflow Somatic Variant Calling Pipeline (`main.nf`) updates

**Processes removed:**
- `MUTECT2`, `LEARN_ORIENTATION`, `FILTER_MUTECT2` -- replaced by HaplotypeCaller stack
- `DELLY_CALL`, `MERGE_SV_SITES`, `DELLY_GENOTYPE`, `MERGE_GENO_GROUP` -- replaced by Arriba

**Processes added:**

`SPLIT_N_CIGAR` (step 7) -- the critical missing step. Without it, GATK fires false positives at every splice junction because N-CIGAR operations look like indels. Goes between REMAP_MAPQ and BQSR.

`HAPLOTYPE_CALLER` + `GENOTYPE_GVCFS` -- diploid caller in GVCF mode. Gets you germline + somatic in one pass. `--dont-use-soft-clipped-bases` is mandatory for RNA-seq; `-stand-call-conf 20` is relaxed from 30 to compensate for 20M read depth.

`VARIANT_FILTRATION` + `SELECT_PASS` -- GATK RNA-seq hard filters (QD, FS, MQ, ReadPosRankSum thresholds). VQSR doesn't work on RNA-seq because annotation distributions are incompatible with WGS training data.

`VEP_ANNOTATE` -- produces CSQ field (consequence + HGVS protein notation) that pypgatk reads.

`PYPGATK_FASTA` -- VCF to mutant protein FASTA. Headers get a per-sample tag so MaxQuant/MSFragger results are traceable back to which sample contributed the variant.

`ARRIBA` -- STAR chimeric reads go here. STAR_ALIGN now emits both the main BAM and `Chimeric.out.sam` (via `WithinBAM SeparateSAMold`). The `outFilterMultimapNmax` was raised from 1 to 50 for Arriba; GATK ignores multimappers naturally via MAPQ.

`TK_PROGRESSION_SUB` -- `bcftools isec -C` keeps only variants private to the progression sample, absent from TK12. These are the variants that appeared between TK12 and TK13/TK14.

`PROGRESSION_FASTA` -- separate FASTA for progression-only variants, with `|PROGRESSION|` in the header so you can filter MS results by progression-acquired variants specifically.

**Samplesheet change needed:** add a `baseline` column. TK12 gets `true`, TK13/TK14 get `false`, all others get empty string. The workflow branches on this field.

**New params required:** `params.vep_cache`, `params.vep_assembly`, `params.reference_proteome`, `params.arriba_blacklist`, `params.arriba_known_fusions`, `params.arriba_protein_domains`. Arriba's database files are bundled with its GitHub release.

Caveat: `PROGRESSION_FASTA` currently calls pypgatk directly on the progression VCF without a fresh VEP run, which means it relies on CSQ annotations already present from the VEP step upstream. 

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

### Planned: Mutant FASTA generation

After variant calling, we plan to generate sample-specific FASTA files incorporating the expected somatic variants (SNPs and SNVs) identified by the pipeline. This will involve:

- Taking the per-group merged VCFs (SNVs from Mutect2, SVs from DELLY2)
- Applying high-confidence PASS variants to the reference genome (GRCh38) using a tool such as `bcftools consensus` or GATK `FastaAlternateReferenceMaker`
- Producing per-sample/per-group mutant reference FASTAs for downstream use in:
  - Neoantigen prediction (e.g. pVACseq / MHCflurry)
  - Proteogenomics database construction (custom protein FASTA for MS/MS searches)
  - Validation of variant-containing peptides against the PGTK proteomics data

---

## License

See [LICENSE](LICENSE) file for details.

## Contact

me 🤓