# PGTK: Proteogenomics of Myeloma Cell Lines

## Overview

PGTK, ProteoGenomics TK, supports integrated transcriptomic and proteomic analysis of seven patient-derived multiple myeloma cell lines with a comparatively primary cell-like phenotype.

This repository contains a Nextflow DSL2 workflow for exploratory RNA-seq proteogenomics. It builds sample-specific and progression-specific variant protein FASTA databases for downstream mass-spectrometry searches.

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

The workflow uses the SRR accessions listed in `samples.csv`. SRA archives are downloaded and validated on the Saga login node, then converted into paired FASTQ files by the Nextflow workflow on compute nodes.

### Proteomics

- PRIDE PXD033531: https://www.ebi.ac.uk/pride/archive/projects/PXD033531
- PRIDE PXD033510: https://www.ebi.ac.uk/pride/archive/projects/PXD033510

Example recursive downloads:

```bash
wget -r "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2024/11/PXD033531/"
wget -r "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2024/11/PXD033510/"
```

## Analysis Strategy

The workflow is intended for exploratory proteogenomic database construction rather than clinical somatic-mutation classification.

Because no matched normal sample is available, Mutect2 tumor-normal calling is not used. Instead, GATK HaplotypeCaller creates a broad catalogue of expressed per-sample variants after RNA-seq-specific preprocessing. The resulting calls may include germline, clonal, progression-associated, RNA-edited, and technical events. They should not be interpreted as validated somatic mutations without additional evidence.

For exploratory sensitivity, VEP retains consequences for all overlapping transcripts. Neither `--pick` nor `--flag_pick` is used. pypgatk therefore receives the complete CSQ consequence set and can generate a wider transcript-aware variant protein database.

The principal outputs are:

- Per-sample PASS VCF files
- VEP-annotated VCF files
- Expanded per-sample variant protein FASTA files
- Arriba fusion candidate tables
- TK13-minus-TK12 and TK14-minus-TK12 progression VCF files
- Progression-specific protein FASTA files
- Raw and trimmed FastQC reports
- Alignment and variant summary statistics
- A consolidated MultiQC report

## Workflow

```text
samples.csv
  |
  +--> local SRA archives
         |
         +--> fasterq-dump --split-files
         +--> FASTQ concatenation
         +--> raw FastQC
         +--> Trim Galore
         +--> trimmed FastQC
         +--> STAR two-pass alignment
                |
                +--> Arriba-compatible unsorted chimeric BAM
                |      +--> Arriba fusion calling
                |
                +--> sort and index BAM
                       +--> samtools flagstat
                       +--> MarkDuplicates
                       +--> SplitNCigarReads
                       +--> HaplotypeCaller GVCF
                       +--> GenotypeGVCFs
                       +--> hard filtering and PASS selection
                       +--> bcftools stats
                       +--> VEP with all overlapping transcript consequences
                       +--> pypgatk expanded variant protein FASTA
                       +--> baseline/progression pairing
                              +--> bcftools baseline subtraction
                              +--> progression-specific FASTA

All QC logs and reports
  +--> MultiQC
```

## QC and Reporting

The workflow runs FastQC before and after trimming. MultiQC aggregates:

- Raw FastQC reports
- Trimmed FastQC reports
- Trim Galore/Cutadapt reports
- STAR `Log.final.out`
- `samtools flagstat`
- Picard/GATK MarkDuplicates metrics
- `bcftools stats`

Final reports:

```text
results/multiqc/multiqc_report.html
results/multiqc/multiqc_report_data/
```

FastQC and MultiQC use these pinned containers:

```text
quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0
quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1
```

## Important Implementation Details

- `fasterq-dump --split-files` converts each paired-end SRA archive into `_1.fastq` and `_2.fastq` files.
- Trim Galore automatically detects common adapter types, performs adapter trimming, trims low-quality ends at Phred 20, and removes pairs shorter than 36 bases after trimming.
- FastQC runs before and after trimming.
- STAR produces an Arriba-compatible unsorted BAM using `WithinBAM HardClip`.
- Arriba consumes the original STAR BAM directly.
- A separate `SORT_INDEX_BAM` branch creates the coordinate-sorted and indexed BAM required by GATK.
- `SAMTOOLS_FLAGSTAT` records alignment summary statistics for MultiQC.
- `SplitNCigarReads` handles RNA-seq splice-junction alignments before variant calling.
- VEP does not restrict output to a selected transcript.
- pypgatk 0.0.24 translates protein-altering consequences into sample-specific FASTA sequences.
- Both pypgatk branches validate and decompress BGZF VCF input to plain-text VCF before translation. Header-only VCFs produce an empty FASTA without invoking pypgatk.
- `bcftools isec -C -w 1` retains progression-sample variants absent from the matched baseline callset.
- `BCFTOOLS_STATS` records PASS-VCF summary statistics for MultiQC.
- Only `STAR_INDEX` and `STAR_ALIGN` use the Saga `bigmem` partition.

## Repository Files

```text
main.nf                 Nextflow DSL2 workflow
download_assets.sh      Login-node downloader and validator for references and containers
download_sra.sh         Login-node downloader and validator for SRA archives
scratch.slurm            Saga launcher for Nextflow
samples.csv              Sample metadata
singularity_cache/       Pre-downloaded Singularity images
reference_downloads/     Pre-downloaded reference archives
sra_cache/               Pre-downloaded SRA archives
work/                    Nextflow work directory
results/                 Published pipeline outputs
```

`download_assets.sh` and `download_sra.sh` perform download validation. `scratch.slurm` intentionally does not repeat validation of downloaded files. It launches Nextflow using the local paths prepared by the downloader scripts.

## Requirements

The tested Saga environment uses:

- Nextflow 26.04.6
- Java 21
- Singularity command backed by Apptainer 1.4.4
- SLURM account `nn9036k`
- Internet access on the login node for initial downloads
- No internet dependency for compute-node analysis

The workflow uses pre-downloaded BioContainers for:

- SRA Tools 3.2.1
- Trim Galore 0.6.10
- FastQC 0.12.1
- STAR 2.7.11b
- GATK 4.6.1.0
- Ensembl VEP 111
- pypgatk 0.0.24
- Arriba 2.4.0
- bcftools 1.21
- MultiQC 1.35

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
- `TK`: longitudinal key used to pair baseline and progression samples
- `Group`: biological group label
- `baseline`: `true` for the earliest sample, `false` for a later sample, or blank for a sample without a matched baseline

Exactly one baseline should be marked `true` for each longitudinal key used for progression subtraction.

## Setup on Saga

Run all commands from the repository directory.

### 1. Download and Validate Containers and References

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

Downloads are resumable where the remote server supports byte ranges. The UniProt streaming endpoint is downloaded afresh to a temporary file when needed because it does not support range-based resume.

### 2. Download and Validate SRA Archives

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
find sra_cache -type f -name '*.sra' -printf '%p %s bytes
'
```

`download_sra.sh` runs `prefetch` on the login node, binds the repository directory into Singularity, retries transient failures, and validates completed archives with `vdb-validate`.

### 3. Validate the SLURM Launcher

The SLURM launcher validates only Bash and SLURM syntax. It does not revalidate downloaded containers, references, or SRA files.

```bash
bash -n scratch.slurm
sbatch --test-only scratch.slurm
```

### 4. Submit the Pipeline

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

```text
Process               Partition   CPUs   Memory
STAR_INDEX            bigmem       32     32 GB
STAR_ALIGN            bigmem       32     32 GB
HAPLOTYPE_CALLER      normal       10     16 GB
VEP_ANNOTATE          normal       10     16 GB
SRA_TO_FASTQ          normal        8     16 GB
SORT_INDEX_BAM        normal        8     16 GB
FASTQC_RAW            normal        4      8 GB
FASTQC_TRIMMED        normal        4      8 GB
MULTIQC               normal        2      8 GB
```

All other processes use the `normal` partition. The Nextflow launcher uses 2 CPUs and 8 GB on `normal`.

## Output Structure

```text
results/
├── references/
│   └── star_index/
├── qc/
│   ├── fastqc_raw/
│   ├── fastqc_trimmed/
│   ├── trim_galore/
│   ├── flagstat/
│   └── bcftools/
├── multiqc/
│   ├── multiqc_report.html
│   └── multiqc_report_data/
├── bam/
│   └── star/
├── gvcf/
├── vcf_pass/
├── vep/
├── variant_fasta/
├── fusions/
├── progression_vcf/
├── progression_fasta/
├── pipeline_report-*.html
├── pipeline_timeline-*.html
└── pipeline_trace-*.tsv
```

## Interpretation and Limitations

- RNA-seq detects variants only in expressed and sufficiently covered regions.
- RNA editing, mapping artefacts, allele-specific expression, and transcript structure can affect calls.
- HaplotypeCaller outputs are not equivalent to validated tumor-only somatic calls.
- TK12 subtraction can reflect biological acquisition or insufficient TK12 expression or coverage.
- Retaining all transcript consequences substantially expands the database and its multiple-testing burden. Strict target-decoy FDR control and transcript-aware validation are required.
- Arriba outputs fusion candidates. Fusion protein FASTA generation is not yet included.
- Alternative-splicing-derived protein FASTA generation is not yet included.
- Variant FASTA files should be combined with a canonical proteome and an appropriate contaminant database before MS/MS searching.
- Proteogenomic peptide validation should include inspection of variant-supporting spectra and confirmation that identified peptides span the altered residue or junction.

## Downstream Proteomics

Recommended downstream steps:

1. Combine the canonical reviewed human proteome, sample-specific variant proteins, progression proteins, and contaminants.
2. Search the corresponding PRIDE raw files with MSFragger/FragPipe or another proteogenomics-capable search engine.
3. Control peptide-spectrum match, peptide, and protein FDR.
4. Require peptides to span the altered residue or fusion junction.
5. Compare variant peptide detections with RNA read support, depth, and allele fraction.

## Current Runtime Safeguards

GATK tasks use task-local temporary storage through `-Djava.io.tmpdir=${PWD}/gatk_tmp`. `SPLIT_N_CIGAR` normalizes `<sample>.split.bai` to `<sample>.split.bam.bai` and validates both outputs.

Both pypgatk branches receive an uncompressed plain-text VCF, validate the header, handle header-only VCFs, and request 80 GB task disk.

MultiQC 1.35 writes `multiqc_report.html` and `multiqc_report_data/`. The Nextflow output contract uses these exact names.

## Comprehensive Validation

The synchronized configuration passed on Saga on July 25, 2026:

```text
PASS: 146
WARN: 0
FAIL: 0
RESULT: PASSED
```

Validation report:

```text
pipeline_command_validation_20260725_205721.txt
```

Validation includes all 27 processes, 15 images, reference archives, pypgatk VCF handling, GATK temporary storage, split-index normalization, the gffcompare fixture, the MultiQC output contract, Apptainer settings, and `sbatch --test-only`.

The HTML detector checks literal encoded entities without misclassifying valid `>`, `<`, `&&`, or `->` operators.

## Citation

If you use this repository, its data, or derived results, cite:

```text
Cancers. 2024;16(23):3963.
https://www.mdpi.com/2072-6694/16/23/3963
```

## License

See [LICENSE](LICENSE) for details.

## Contact

For questions about the workflow or data, open an issue in this repository or contact the project maintainer.
