# PGTK: Proteogenomics of Myeloma Cell Lines

## Overview

PGTK, ProteoGenomics TK, supports integrated transcriptomic and proteomic analysis of seven patient-derived multiple myeloma cell lines with a comparatively primary cell-like phenotype.

This repository contains a Nextflow DSL2 workflow for exploratory RNA-seq proteogenomics. It generates sample-specific, progression-specific, fusion-derived, novel-splicing-derived, progression-specific, and combined protein FASTA databases for downstream mass-spectrometry searches.

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

The workflow uses the SRR accessions listed in `samples.csv`. SRA archives are downloaded and validated on the Saga login node, retained under `sra_cache/`, and converted into paired FASTQ files by Nextflow tasks on compute nodes.

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

Because no matched normal sample is available, Mutect2 tumor-normal calling is not used. GATK HaplotypeCaller instead creates a broad catalogue of expressed per-sample variants after RNA-seq-specific preprocessing. Calls may include germline, clonal, progression-associated, RNA-edited, and technical events. They must not be interpreted as validated somatic mutations without additional evidence.

For exploratory sensitivity, VEP retains consequences for all overlapping transcripts. Neither `--pick` nor `--flag_pick` is used. pypgatk receives the complete CSQ consequence set and uses the Ensembl release 111 cDNA transcript FASTA and matching GTF to generate an expanded transcript-aware variant protein database.

Additional branches generate:

- Fusion-junction proteins from accepted Arriba candidates using pVACfuse
- Expressed transcript assemblies using StringTie
- gffcompare classification against Ensembl release 111
- ORF predictions for selected novel transcript classes using TransDecoder
- Progression VCF and protein FASTA outputs after baseline subtraction
- Deduplicated combined FASTA databases containing reviewed, variant, fusion, and splice-derived proteins

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
                       |      +--> gffcompare classification
                       |      +--> retain class codes j and u
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
                       +--> VEP all-transcript annotation
                       +--> pypgatk variant protein FASTA
                       +--> baseline subtraction
                              +--> progression-specific VCF
                              +--> progression-specific protein FASTA

reviewed proteome + variant FASTA + fusion FASTA + splice FASTA
  +--> exact sequence deduplication
  +--> combined exploratory proteogenomics FASTA

all QC logs and reports
  +--> MultiQC
```

## Storage Model

Project files, downloaded assets, and published outputs remain under the repository directory:

```text
$PWD/
├── main.nf
├── scratch.slurm
├── samples.csv
├── singularity_cache/
├── reference_downloads/
├── sra_cache/
└── results/
```

Large execution data is placed on the Saga work filesystem:

```text
/cluster/work/users/ash022/
├── work/    Nextflow task directories, intermediates, and resume cache
└── tmp/     launcher, Java, Singularity, and Apptainer temporary files
```

The launcher passes:

```text
-work-dir /cluster/work/users/ash022/work
```

Nextflow creates a unique work directory for each task. Inputs are staged there, commands run there, outputs are cached there, and compatible completed tasks are reused with `-resume`.

`fasterq-dump` uses a `fasterq_tmp` subdirectory inside its task directory and writes FASTQ output to that same task directory. No external fasterq-specific path or custom Singularity bind is required.

Containers, reference archives, and SRA archives remain under `$PWD`. Compute-node analysis does not redownload them.

## Important Implementation Details

- `fasterq-dump --split-files` creates paired `_1.fastq` and `_2.fastq` files.
- Trim Galore trims adapters and Phred-20 low-quality ends, then removes pairs shorter than 36 bases.
- FastQC runs before and after trimming.
- STAR produces an Arriba-compatible unsorted BAM using `WithinBAM HardClip`.
- Arriba consumes the original STAR BAM directly.
- `SORT_INDEX_BAM` creates the coordinate-sorted and indexed BAM required by GATK and StringTie.
- `REF_INDEX`, `SORT_INDEX_BAM`, and `SAMTOOLS_FLAGSTAT` use SAMtools 1.21.
- `REF_INDEX` creates `genome.fa.fai` and `genome.dict` with `samtools faidx` and `samtools dict`.
- `SplitNCigarReads` handles RNA-seq splice-junction alignments before variant calling.
- VEP retains all overlapping transcript consequences.
- pypgatk 0.0.24 uses underscore-style options including `--input_fasta`, `--gene_annotations_gtf`, and `--output_proteindb`.
- The pypgatk nucleotide input is Ensembl release 111 cDNA, not the UniProt protein FASTA.
- The reviewed UniProt proteome is the canonical base of the final combined database.
- `bcftools isec -C -w 1` retains progression variants absent from the matched baseline callset.
- gffcompare retains class codes `j` and `u` by default.
- TransDecoder utilities are invoked using verified absolute paths under `/usr/local/opt/transdecoder/util/`.
- Large unpacked references and the STAR index remain in the Nextflow cache and are not duplicated into `results/`.
- Nextflow uses fail-fast behavior. An unhandled task failure stops the workflow and terminates running tasks. Completed compatible tasks are recovered by `-resume`.

## Full-Node Resource Strategy

Normal Saga nodes provide 20 CPUs and 80 GB memory. Compute-intensive normal-partition processes request up to 20 CPUs and 64 GB, leaving memory for the operating system and container overhead.

`MARK_DUPLICATES` previously failed with `OUT_OF_MEMORY` at 16 GB. It now requests 20 CPUs and 64 GB, uses a bounded 56 GB Java heap, and sets `--MAX_RECORDS_IN_RAM 1000000`.

`STAR_ALIGN` remains on `bigmem` with 32 CPUs and 256 GB after an earlier exit-137 failure at 32 GB.

```text
Process                 Partition   CPUs   Memory   Time    Disk
DOWNLOAD_REFERENCES     normal        8     32 GB   12h    100 GB
SRA_TO_FASTQ            normal       16     32 GB   24h    150 GB
CAT_FASTQ               normal        2      8 GB    4h    200 GB
FASTQC_RAW              normal        8     16 GB    8h    100 GB
TRIM_GALORE             normal        8     16 GB   12h    150 GB
FASTQC_TRIMMED          normal        8     16 GB    8h    100 GB
STAR_INDEX              normal       20     64 GB   12h    320 GB
STAR_ALIGN              bigmem       32    256 GB   24h    320 GB
SORT_INDEX_BAM          normal       20     64 GB   24h    200 GB
SAMTOOLS_FLAGSTAT       normal        4     16 GB    4h     20 GB
REF_INDEX               normal        4     16 GB    4h     30 GB
MARK_DUPLICATES         normal       20     64 GB   24h    200 GB
SPLIT_N_CIGAR           normal       20     64 GB   24h    200 GB
HAPLOTYPE_CALLER        normal       20     64 GB   48h    120 GB
GENOTYPE_FILTER         normal        8     32 GB   12h     40 GB
BCFTOOLS_STATS          normal        2      8 GB    4h     20 GB
VEP_ANNOTATE            normal       20     64 GB   24h     60 GB
PYPGATK_FASTA           normal        8     32 GB   12h     40 GB
ARRIBA                   normal        8     32 GB   12h     60 GB
FUSION_FASTA            normal        4     16 GB    8h     20 GB
STRINGTIE_ASSEMBLY      normal       20     64 GB   24h     80 GB
GFFCOMPARE_NOVEL        normal        4     16 GB    8h     30 GB
SPLICE_PROTEIN_FASTA    normal       20     64 GB   24h     60 GB
COMBINE_PROTEIN_FASTA   normal        2      8 GB    4h     30 GB
PROGRESSION_SUBTRACT    normal        2      8 GB    4h     20 GB
PROGRESSION_FASTA       normal        8     32 GB   12h     40 GB
MULTIQC                 normal        8     32 GB    4h     40 GB
```

The Nextflow launcher requests 4 CPUs and 16 GB on `normal`.

## Launcher Configuration

The launcher preserves Saga `StdEnv` and loads Java with:

```bash
module purge
module load Java/21
```

Do not use `module --force purge`, because it removes `StdEnv` and can make `Java/21` unavailable on compute nodes.

The SLURM account is passed to Nextflow as one argument:

```bash
"-process.clusterOptions=--account=nn9036k"
```

Do not use the two-argument form below, because Nextflow can parse `clusterOptions` as Boolean `true`:

```bash
-process.clusterOptions "--account=nn9036k"
```

The streamlined launcher enables `-resume` and creates only a task trace. HTML execution reports and timelines are not enabled by default.

## Exploratory Database Filters

```text
Minimum StringTie coverage:         2.5
Minimum junction-supporting reads:  3
Minimum isoform fraction:           0.05
Retained gffcompare class codes:    j,u
Minimum predicted protein length:   60 amino acids
```

Class code `j` represents a potentially novel splice-junction combination. Class code `u` represents an intergenic transcript. Both remain exploratory and require supporting evidence.

Exact duplicate protein sequences are removed from splice-derived FASTA files and from the final combined database.

Fusion FASTA generation uses 50 amino acids around each fusion breakpoint where possible and retains the full downstream sequence for frameshift fusions. If Arriba reports no accepted candidates, an empty fusion FASTA is created so the workflow can continue.

## QC and Reporting

MultiQC aggregates raw and trimmed FastQC reports, Trim Galore and Cutadapt reports, STAR `Log.final.out`, `samtools flagstat`, MarkDuplicates metrics, and `bcftools stats`.

```text
results/multiqc/multiqc_report.html
results/multiqc/multiqc_data/
results/pipeline_trace-<launcher-job-id>.tsv
```

## Validation

Run before submission:

```bash
dos2unix probe_pipeline_cli.sh validate_pipeline_commands.sh scratch.slurm
chmod +x probe_pipeline_cli.sh validate_pipeline_commands.sh scratch.slurm
bash -n probe_pipeline_cli.sh
bash -n validate_pipeline_commands.sh
bash -n scratch.slurm
bash probe_pipeline_cli.sh
bash validate_pipeline_commands.sh
sbatch --test-only scratch.slurm
```

The synchronized configuration currently passes 39 checks:

```text
PASS: 39
WARN: 0
FAIL: 0
RESULT: PASSED
```

The validator rejects force-purging modules, malformed `clusterOptions`, obsolete storage paths, custom container binds, and launcher report or timeline options.

## Requirements

- Nextflow 26.04.6
- Java 21
- Singularity command backed by Apptainer 1.4.4
- SLURM account `nn9036k`
- Internet access on the login node for initial downloads
- No internet dependency for compute-node analysis
- Fifteen pre-downloaded container images under `singularity_cache/`

## Input Samplesheet

```csv
sample,srr,TK,Group,baseline
TK12,SRR31089074,patient1,resistant,true
TK13,SRR31089073,patient1,sensitive,false
TK14,SRR31089072,patient1,sensitive,false
```

Exactly one baseline should be marked `true` for each longitudinal key used for progression subtraction.

## Setup on HPC/saga

Run from the repository directory.

### 1. Download assets

```bash
bash download_assets.sh
bash download_sra.sh
```

Do not remove `singularity_cache/`, `reference_downloads/`, or `sra_cache/` when resetting an execution.

### 2. Prepare storage

```bash
mkdir -p     /cluster/work/users/ash022/work     /cluster/work/users/ash022/tmp

test -w /cluster/work/users/ash022/work
test -w /cluster/work/users/ash022/tmp
```

For a completely fresh execution:

```bash
rm -rf /cluster/work/users/ash022/work/*
rm -rf /cluster/work/users/ash022/tmp/*
rm -rf .nextflow
rm -f .nextflow.log
```

### 3. Validate and submit

```bash
bash validate_pipeline_commands.sh
sbatch --test-only scratch.slurm

JOBID=$(sbatch --parsable scratch.slurm)
echo "$JOBID"

while [[ ! -f "resultsTKvep-${JOBID}.log" ]]; do
    sleep 2
done

tail -f "resultsTKvep-${JOBID}.log"
```

Monitor jobs and storage:

```bash
watch -n 60 '
squeue -u "$USER"   -o "%.18i %.42j %.10P %.10T %.6C %.12m %.20R"

echo

du -sh /cluster/work/users/ash022/work 2>/dev/null
du -sh /cluster/work/users/ash022/tmp 2>/dev/null
'
```

Monitor task resources:

```bash
watch -n 60 '
for job in $(squeue -h -u "$USER" -o "%A"); do
    sstat -j "${job}.batch"         --format=JobID,MaxRSS,AveRSS,AveCPU 2>/dev/null
done
'
```

## Principal Output Structure

```text
results/
├── qc/
├── multiqc/
├── bam/star/
├── gvcf/
├── vcf_pass/
├── vep/
├── variant_fasta/
├── fusions/
├── fusion_fasta/
├── splicing/stringtie/
├── splicing/gffcompare/
├── splice_fasta/
├── combined_fasta/
├── progression_vcf/
├── progression_fasta/
└── pipeline_trace-*.tsv
```

The principal combined database is:

```text
results/combined_fasta/<sample>.exploratory_proteogenomics.fasta
```

Large unpacked references and the STAR index remain cached under `/cluster/work/users/ash022/work` and are not duplicated into `results/`.

## Interpretation and Limitations

- RNA-seq detects variants and isoforms only in expressed and sufficiently covered regions.
- RNA editing, mapping artefacts, allele-specific expression, assembly errors, and incomplete ORFs can affect results.
- HaplotypeCaller outputs are not equivalent to validated tumor-only somatic calls.
- TK12 subtraction can reflect biological acquisition or insufficient TK12 expression or coverage.
- Retaining all transcript consequences, fusion proteins, and splice-derived ORFs increases database size and the multiple-testing burden.
- Strict target-decoy FDR control is required.
- Variant peptides should span the altered residue, fusion peptides the fusion junction, and splice-derived peptides a novel exon junction or ORF.
- Fusion and splice-derived sequences remain exploratory and require orthogonal review of RNA support.

## Downstream Proteomics

1. Search each sample against its corresponding combined exploratory FASTA.
2. Include contaminants and target-decoy sequences.
3. Control peptide-spectrum match, peptide, and protein FDR.
4. Require variant, fusion, or splice-supporting peptides to cover the relevant altered residue or junction.
5. Compare peptide detections with RNA depth, allele fraction, and transcript abundance.
6. Treat canonical, variant, fusion, and splice-derived identifications as separate evidence classes.

## Citation

```text
Cancers. 2024;16(23):3963.
https://www.mdpi.com/2072-6694/16/23/3963
```

## License

See [LICENSE](LICENSE) for details.

## Contact

For questions, open an issue or contact `animesh@fuzzylife.org`.
