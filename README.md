# PGTK: Proteogenomics of Myeloma Cell Lines

## Overview

PGTK, ProteoGenomics TK, supports integrated transcriptomic and proteomic analysis of seven patient-derived multiple myeloma cell lines with a comparatively primary cell-like phenotype.

This repository contains a Nextflow DSL2 workflow for RNA-seq proteogenomics. The workflow builds sample-specific and progression-specific variant protein FASTA databases for downstream mass-spectrometry searches.

The current analysis focuses on the longitudinal TK12, TK13, and TK14 samples. TK12 is treated as the earliest available baseline, while TK13 and TK14 represent later progression samples from the same patient.

## Publication

The study methodology and biological findings are described in:

**Cancers 2024, Volume 16, Issue 23, Article 3963**

- Article: https://www.mdpi.com/2072-6694/16/23/3963
- PDF: https://www.mdpi.com/2072-6694/16/23/3963/pdf

Please cite the article when using the associated data, workflow, or results.

## Public Data

### RNA-seq

- NCBI BioProject: PRJNA1176350
- Access: https://www.ncbi.nlm.nih.gov/bioproject/1176350

The workflow uses the SRR accessions listed in `samples.csv`. SRA archives are downloaded on the Saga login node and converted to paired FASTQ files by the Nextflow workflow on compute nodes.

### Proteomics

- PRIDE PXD033531: https://www.ebi.ac.uk/pride/archive/projects/PXD033531
- PRIDE PXD033510: https://www.ebi.ac.uk/pride/archive/projects/PXD033510

Example recursive downloads:

```bash
wget -r "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2024/11/PXD033531/"
wget -r "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2024/11/PXD033510/"
```

## Analysis Strategy

The workflow is designed for proteogenomic database construction rather than clinical somatic-mutation classification.

Because no matched normal sample is available, Mutect2 tumor-normal calling is not used. Instead, the workflow creates a broad catalogue of expressed per-sample variants with GATK HaplotypeCaller after RNA-seq-specific preprocessing. The resulting calls may include germline, clonal, progression-associated, and RNA-specific events. They should not be interpreted as validated somatic mutations without additional evidence.

The principal outputs are:

- Per-sample PASS VCF files
- VEP-annotated VCF files
- Per-sample variant protein FASTA files
- Arriba fusion tables
- TK13-minus-TK12 and TK14-minus-TK12 progression VCF files
- Progression-specific protein FASTA files

## Workflow

The main processing path is:

```text
samples.csv
  |
  +--> local .sra archives
  |      |
  |      +--> fasterq-dump
  |      +--> FASTQ concatenation
  |      +--> Trim Galore
  |      +--> STAR two-pass alignment
  |             |
  |             +--> unsorted chimeric BAM --> Arriba
  |             |
  |             +--> sort and index BAM
  |                    +--> MarkDuplicates
  |                    +--> SplitNCigarReads
  |                    +--> HaplotypeCaller GVCF
  |                    +--> GenotypeGVCFs
  |                    +--> hard filtering and PASS selection
  |                    +--> VEP annotation
  |                    +--> pypgatk variant protein FASTA
  |
  +--> TK baseline/progression pairing
         +--> bcftools baseline subtraction
         +--> progression-specific protein FASTA
```

### Important implementation details

- `SplitNCigarReads` is run before variant calling to handle RNA-seq splice-junction alignments.
- STAR produces one Arriba-compatible unsorted BAM using `WithinBAM HardClip`.
- Arriba consumes the original STAR BAM directly.
- A separate `SORT_INDEX_BAM` branch creates the coordinate-sorted and indexed BAM required by the GATK branch.
- VEP uses `--pick` to retain one selected transcript consequence per variant and limit FASTA inflation.
- pypgatk 0.0.24 translates protein-altering VEP consequences into variant protein sequences.
- `bcftools isec -C -w 1` retains progression-sample variants absent from the matched baseline.
- Only `STAR_INDEX` and `STAR_ALIGN` use the Saga `bigmem` partition.

## Removed or Replaced Components

The following earlier processes are no longer part of the workflow:

- Mutect2, orientation modelling, and FilterMutectCalls
- DELLY structural-variant calling and group-level DELLY merging
- MAPQ 255-to-60 rewriting
- Redundant sorting of an already sorted STAR BAM
- Sensitive/resistant group-level VCF merging

They were replaced by:

- `SPLIT_N_CIGAR`
- `HAPLOTYPE_CALLER`
- `GENOTYPE_FILTER`
- `VEP_ANNOTATE`
- `PYPGATK_FASTA`
- `ARRIBA`
- `PROGRESSION_SUBTRACT`
- `PROGRESSION_FASTA`

## Repository Files

```text
main.nf                 Nextflow DSL2 workflow
download_assets.sh      Login-node downloader and validator for references
download_sra.sh         Login-node downloader and validator for SRA archives
scratch.slurm            Saga launcher for Nextflow
samples.csv              Sample metadata
singularity_cache/       Pre-downloaded Singularity images
reference_downloads/     Pre-downloaded reference archives
sra_cache/               Pre-downloaded SRA archives
work/                    Nextflow work directory
results/                 Published pipeline outputs
```

## Requirements

The tested Saga environment uses:

- Nextflow 26.04.6
- Java 21
- Singularity command backed by Apptainer 1.4.4
- SLURM account `nn9036k`
- Internet access on the login node for initial asset downloads
- No internet dependency for compute-node analysis

The workflow uses pre-downloaded BioContainers for SRA Tools, Trim Galore, STAR, GATK, VEP, pypgatk, Arriba, and bcftools.

## Input Samplesheet

The required file is `samples.csv`:

```csv
sample,srr,TK,Group,baseline
TK12,SRR31089074,patient1,resistant,true
TK13,SRR31089073,patient1,sensitive,false
TK14,SRR31089072,patient1,sensitive,false
```

Column definitions:

- `sample`: biological sample identifier
- `srr`: SRA run accession
- `TK`: longitudinal patient key used to pair baseline and progression samples
- `Group`: biological group label
- `baseline`: `true` for the earliest sample, `false` for a later sample, or blank for a sample without a matched baseline

## Setup on Saga

Run all commands from the repository directory.

### 1. Download and validate containers and references

```bash
chmod +x download_assets.sh
bash download_assets.sh
```

The downloader creates:

```text
singularity_cache/
reference_downloads/
```

Reference assets include:

- Ensembl GRCh38 primary assembly, release 111
- Ensembl release 111 GTF
- Ensembl VEP 111 GRCh38 cache
- Reviewed UniProt human proteome with isoforms
- Arriba 2.4.0 reference resources

Downloads are resumable. The VEP cache is large and may require repeated resume attempts through the Saga proxy.

### 2. Download and validate SRA archives

```bash
chmod +x download_sra.sh
bash download_sra.sh
```

Expected files:

```text
sra_cache/SRR31089072/SRR31089072.sra
sra_cache/SRR31089073/SRR31089073.sra
sra_cache/SRR31089074/SRR31089074.sra
```

Check them with:

```bash
find sra_cache -type f -name '*.sra' -printf '%p %s bytes\n'
```

`download_sra.sh` runs `prefetch` on the login node, explicitly binds the repository directory into Singularity, retries transient failures, and validates completed archives with `vdb-validate`.

### 3. Validate the SLURM launcher

```bash
bash -n scratch.slurm
sbatch --test-only scratch.slurm
```

### 4. Submit the pipeline

```bash
JOBID=$(sbatch --parsable scratch.slurm)
echo "$JOBID"
tail -f "resultsTKvep-${JOBID}.log"
```

Monitor launcher and child jobs:

```bash
watch -n 5 "squeue -u \$USER -o '%.18i %.38j %.10P %.10T %.6C %.12m %.20R'"
```

The launcher invokes Nextflow with `-resume`, so completed compatible tasks are reused after interruption or workflow correction.

## Saga Resource Routing

The current principal resource requests are:

```text
Process               Partition   CPUs   Memory
STAR_INDEX            bigmem       32     32 GB
STAR_ALIGN            bigmem       32     32 GB
HAPLOTYPE_CALLER      normal       10     16 GB
VEP_ANNOTATE          normal       10     16 GB
SRA_TO_FASTQ          normal        8     16 GB
SORT_INDEX_BAM        normal        8     16 GB
```

All other processes use the `normal` partition. The Nextflow launcher itself uses 2 CPUs and 8 GB on `normal`.

## Output Structure

```text
results/
├── references/             Prepared genome, annotation, proteome, VEP, and Arriba resources
│   └── star_index/         STAR genome index
├── qc/
│   └── trim_galore/        Trimming reports
├── bam/
│   └── star/               Coordinate-sorted and indexed STAR BAM files
├── gvcf/                   Per-sample GVCFs
├── vcf_pass/               Per-sample filtered PASS VCF files
├── vep/                    VEP-annotated VCF files
├── variant_fasta/          Per-sample variant protein FASTA files
├── fusions/                Arriba fusion and discarded-fusion tables
├── progression_vcf/        Variants absent from the matched baseline
├── progression_fasta/      Progression-specific variant protein FASTA files
├── pipeline_report-*.html  Nextflow execution reports
├── pipeline_timeline-*.html
└── pipeline_trace-*.tsv
```

## Interpretation and Limitations

- RNA-seq only detects variants in expressed and sufficiently covered regions.
- RNA editing, mapping artefacts, allele-specific expression, and transcript structure can affect calls.
- HaplotypeCaller outputs are not equivalent to validated tumor-only somatic calls.
- TK12 subtraction identifies variants absent from the TK12 RNA-seq callset. Absence may reflect biology or insufficient TK12 expression or coverage.
- `--pick` reduces database size but may omit protein consequences from non-selected isoforms.
- Arriba outputs fusion candidates. Fusion protein FASTA generation is not yet included.
- Alternative-splicing-derived protein FASTA generation is not yet included.
- The variant FASTA files should be combined with a canonical proteome and an appropriate contaminant database before MS/MS searching.
- Proteogenomic peptide validation should use strict target-decoy FDR control and inspection of variant-supporting spectra.

## Downstream Proteomics

Recommended downstream steps:

1. Combine the canonical reviewed human proteome, sample-specific variant proteins, progression proteins, and contaminants.
2. Search the corresponding PRIDE raw files with MSFragger/FragPipe or another proteogenomics-capable search engine.
3. Control peptide-spectrum match, peptide, and protein FDR.
4. Require peptides to span the altered residue or fusion junction.
5. Compare variant peptide detections with RNA read support, depth, and allele fraction.

## Citation

If you use this repository, its data, or derived results, cite the associated Cancers article:

```text
Cancers. 2024;16(23):3963.
https://www.mdpi.com/2072-6694/16/23/3963
```

## License

See [LICENSE](LICENSE) for details.

## Contact

For questions about the workflow or data, open an issue in this repository or contact me 🤓