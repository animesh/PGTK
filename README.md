# PGTK: Proteogenomics of Myeloma Cell Lines

## Overview

PGTK, ProteoGenomics TK, supports integrated transcriptomic and proteomic analysis of seven patient-derived multiple myeloma cell lines with a comparatively primary cell-like phenotype.

This repository contains a Nextflow DSL2 workflow for exploratory RNA-seq proteogenomics. It generates sample-specific, progression-specific, fusion-derived, alternative-splicing-derived, and combined protein FASTA databases for downstream mass-spectrometry searches.

The current analysis focuses on the longitudinal TK12, TK13, and TK14 samples. TK12 is treated as the earliest available baseline, while TK13 and TK14 represent later progression samples from the same patient.

## Publication

The study methodology and biological findings are described in:

**Cancers 2024, Volume 16, Issue 23, Article 3963**

- Article: https://www.mdpi.com/2072-6694/16/23/3963
- PDF: https://www.mdpi.com/2072-6694/16/23/3963/pdf

Please cite the article when using the associated data, workflow, or results.

## Public Data

### RNA-seq

- NCBI BioProject: `PRJNA1176350`
- Access: https://www.ncbi.nlm.nih.gov/bioproject/1176350

The workflow uses the SRR accessions listed in `samples.csv`. SRA archives are downloaded and validated on the Saga login node, then converted into paired FASTQ files by the Nextflow workflow on compute nodes.

### Proteomics

- PRIDE `PXD033531`: https://www.ebi.ac.uk/pride/archive/projects/PXD033531
- PRIDE `PXD033510`: https://www.ebi.ac.uk/pride/archive/projects/PXD033510

Example recursive downloads:

```bash
wget -r "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2024/11/PXD033531/"
wget -r "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2024/11/PXD033510/"
```

## Analysis Strategy

The workflow is intended for exploratory proteogenomic database construction rather than clinical somatic-mutation classification.

Because no matched normal sample is available, Mutect2 tumor-normal calling is not used. Instead, GATK HaplotypeCaller creates a broad catalogue of expressed per-sample variants after RNA-seq-specific preprocessing. The resulting calls may include germline, clonal, progression-associated, RNA-edited, and technical events. They should not be interpreted as validated somatic mutations without additional evidence.

For exploratory sensitivity, VEP retains consequences for all overlapping transcripts. Neither `--pick` nor `--flag_pick` is used. pypgatk receives the complete CSQ consequence set and generates an expanded transcript-aware variant protein database.

The workflow also generates:

- Fusion-junction protein sequences from accepted Arriba candidates using pVACfuse
- Expressed transcript assemblies using StringTie
- Alternative-splicing and novel-transcript ORF predictions using TransDecoder
- Deduplicated combined databases containing the reviewed proteome, variant proteins, fusion proteins, and splice-derived proteins

## Principal Outputs

- Per-sample PASS VCF files
- VEP-annotated VCF files containing all overlapping transcript consequences
- Expanded per-sample variant protein FASTA files
- Arriba fusion candidate tables
- Fusion-junction protein FASTA files
- StringTie transcript assemblies
- Alternative-splicing-derived protein FASTA files
- TK13-minus-TK12 and TK14-minus-TK12 progression VCF files
- Progression-specific protein FASTA files
- Combined exploratory proteogenomics FASTA files
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
                |      +--> pVACfuse fusion protein FASTA
                |
                +--> coordinate-sorted and indexed BAM
                       +--> samtools flagstat
                       +--> StringTie transcript assembly
                       |      +--> transcript sequence extraction
                       |      +--> TransDecoder ORF prediction
                       |      +--> splice-derived protein FASTA
                       |
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
                              +--> progression-specific protein FASTA

reviewed proteome + variant FASTA + fusion FASTA + splice FASTA
  +--> exact protein-sequence deduplication
  +--> combined exploratory proteogenomics FASTA

all QC logs and reports
  +--> MultiQC
```

## Important Implementation Details

- `fasterq-dump --split-files` converts each paired-end SRA archive into `_1.fastq` and `_2.fastq` files.
- Trim Galore detects common adapter types, performs adapter trimming, trims low-quality ends at Phred 20, and removes pairs shorter than 36 bases after trimming.
- FastQC runs before and after trimming.
- STAR produces an Arriba-compatible unsorted BAM using `WithinBAM HardClip`.
- Arriba consumes the original STAR BAM directly.
- `SORT_INDEX_BAM` creates the coordinate-sorted and indexed BAM required by GATK and StringTie.
- `REF_INDEX`, `SORT_INDEX_BAM`, and `SAMTOOLS_FLAGSTAT` use a dedicated SAMtools 1.21 container.
- `SplitNCigarReads` handles RNA-seq splice-junction alignments before variant calling.
- VEP does not restrict output to a selected transcript.
- pypgatk 0.0.24 translates protein-altering consequences into sample-specific FASTA sequences.
- `bcftools isec -C -w 1` retains progression-sample variants absent from the matched baseline callset.
- `STAR_INDEX` uses the Saga `normal` partition with 8 CPUs and 64 GB memory.
- `STAR_ALIGN` uses the Saga `bigmem` partition with 32 CPUs and 256 GB memory.
- The larger `STAR_ALIGN` allocation follows an earlier exit-137 memory failure with 32 GB.
- Nextflow uses fail-fast behavior by default. An unhandled process failure stops the workflow and terminates other running tasks, including tasks on independent branches. Completed compatible tasks are recovered by the next `-resume` run.

## Exploratory Database Filters

The splice-derived database uses these defaults:

```text
Minimum StringTie coverage:         2.5
Minimum junction-supporting reads:  3
Minimum isoform fraction:           0.05
Minimum predicted protein length:   30 amino acids
```

Exact duplicate protein sequences are removed from splice-derived FASTA files and from the final combined database.

Fusion FASTA generation uses 50 amino acids around each fusion breakpoint where possible and retains the full downstream sequence for frameshift fusions. If Arriba reports no accepted candidates, an empty fusion FASTA is created so the remaining workflow can continue.

## QC and Reporting

MultiQC aggregates:

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
results/multiqc/multiqc_data/
```

## Removed or Replaced Components

The following earlier components are no longer part of the workflow:

- Mutect2, orientation modelling, and FilterMutectCalls
- DELLY structural-variant calling and group-level DELLY merging
- MAPQ 255-to-60 rewriting
- Sensitive/resistant group-level VCF merging

They were replaced by the current HaplotypeCaller, Arriba, progression, expanded FASTA, and QC branches.

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

`download_assets.sh` and `download_sra.sh` perform asset validation. `scratch.slurm` intentionally does not repeat downloaded-file validation and launches Nextflow using the prepared local paths.

## Requirements

The tested Saga environment uses:

- Nextflow 26.04.6
- Java 21
- Singularity command backed by Apptainer 1.4.4
- SLURM account `nn9036k`
- Internet access on the login node for initial downloads
- No internet dependency for compute-node analysis

The workflow uses 14 pre-downloaded container images:

- SRA Tools 3.2.1
- Trim Galore 0.6.10
- FastQC 0.12.1
- STAR 2.7.11b
- GATK 4.6.1.0
- SAMtools 1.21
- Ensembl VEP 111
- pypgatk 0.0.24
- Arriba 2.4.0
- bcftools 1.21
- MultiQC 1.35
- StringTie 3.0.3
- TransDecoder 6.0.0
- pVACtools 7.1.1

`download_assets.sh` resolves complete active BioContainers tags for StringTie and TransDecoder and saves stable local image names used by `main.nf`.

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

The downloader creates or updates:

```text
singularity_cache/
reference_downloads/
```

Reference assets include:

- Ensembl GRCh38 primary assembly, release 111
- Ensembl release 111 GTF
- Ensembl VEP 111 GRCh38 cache
- Reviewed UniProt human proteome with isoforms
- Arriba 2.4.0 resources

Downloads are resumable where supported. The UniProt streaming endpoint is downloaded afresh to a temporary file when needed because it does not support byte-range resume.

Expected final message:

```text
All 14 containers and all reference archives are valid.
```

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
find sra_cache -type f -name '*.sra' -printf '%p %s bytes\n'
```

### 3. Validate the SLURM Launcher

```bash
bash -n scratch.slurm
sbatch --test-only scratch.slurm
```

### 4. Submit or Resume the Pipeline

```bash
JOBID=$(sbatch --parsable scratch.slurm)
echo "$JOBID"
tail -f "resultsTKvep-${JOBID}.log"
```

Monitor launcher and child jobs:

```bash
watch -n 5 "squeue -u \$USER -o '%.18i %.38j %.10P %.10T %.6C %.12m %.20R'"
```

The launcher invokes Nextflow with `-resume`, so compatible completed tasks are reused after interruption or workflow correction.

## Saga Resource Routing

Process                 Partition   CPUs   Memory   Time    Disk
STAR_INDEX              normal        8     64 GB   12h    320 GB
STAR_ALIGN              bigmem       32    256 GB   24h    320 GB
HAPLOTYPE_CALLER        normal       10     16 GB   48h    120 GB
VEP_ANNOTATE            normal       10     16 GB   24h     60 GB
SRA_TO_FASTQ            normal        8     16 GB   24h    150 GB
SORT_INDEX_BAM          normal        8     16 GB   24h    200 GB
STRINGTIE_ASSEMBLY      normal        8     16 GB   24h     80 GB
SPLICE_PROTEIN_FASTA    normal        8     16 GB   24h     60 GB
FUSION_FASTA            normal        2      8 GB    8h     20 GB
FASTQC_RAW              normal        4      8 GB    8h    100 GB
FASTQC_TRIMMED          normal        4      8 GB    8h    100 GB
MULTIQC                 normal        2      8 GB    4h     40 GB

All processes except `STAR_ALIGN` use the `normal` partition. The Nextflow launcher uses 2 CPUs and 8 GB on `normal`.

### Observed STAR Resource Usage

`STAR_INDEX` was initially allocated 256 GB after an earlier 32 GB failure. A successful GRCh38 Ensembl release 111 index build used:

```text
CPUs allocated:       8
Peak memory:          38.06 GB
Wall-clock time:      1:09:06
CPU efficiency:       67.84%

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
│   └── multiqc_data/
├── bam/
│   └── star/
├── gvcf/
├── vcf_pass/
├── vep/
├── variant_fasta/
├── fusions/
├── fusion_fasta/
├── splicing/
│   └── stringtie/
├── splice_fasta/
├── combined_fasta/
├── progression_vcf/
├── progression_fasta/
├── pipeline_report-*.html
├── pipeline_timeline-*.html
└── pipeline_trace-*.tsv
```

The principal combined database is:

```text
results/combined_fasta/<sample>.exploratory_proteogenomics.fasta
```

## Interpretation and Limitations

- RNA-seq detects variants and isoforms only in expressed and sufficiently covered regions.
- RNA editing, mapping artefacts, allele-specific expression, transcript assembly errors, and incomplete ORFs can affect results.
- HaplotypeCaller outputs are not equivalent to validated tumor-only somatic calls.
- TK12 subtraction can reflect biological acquisition or insufficient TK12 expression or coverage.
- Retaining all transcript consequences, fusion proteins, and splice-derived ORFs substantially increases database size and the multiple-testing burden.
- Strict target-decoy FDR control is required.
- Variant peptides should span the altered residue, fusion peptides should span the fusion junction, and splice-derived peptides should support a novel exon junction or ORF.
- Fusion and splice-derived sequences remain exploratory and require orthogonal review of RNA read support.

## Downstream Proteomics

Recommended downstream steps:

1. Search each sample against its corresponding combined exploratory FASTA.
2. Include a contaminant database and target-decoy sequences.
3. Control peptide-spectrum match, peptide, and protein FDR.
4. Require variant, fusion, or splice-supporting peptides to cover the relevant altered residue or junction.
5. Compare peptide detections with RNA read support, depth, allele fraction, and transcript abundance.

## Citation

If you use this repository, its data, or derived results, cite:

```text
Cancers. 2024;16(23):3963.
https://www.mdpi.com/2072-6694/16/23/3963
```

## License

See [LICENSE](LICENSE) for details.

## Contact

For questions about the workflow or data, open an issue or contact `animesh@fuzzylife.org`.
