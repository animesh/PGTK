# PGTK: Proteogenomics of Myeloma Cell Lines

## Overview

PGTK, ProteoGenomics TK, is a Nextflow DSL2 workflow for exploratory RNA-seq proteogenomics of patient-derived multiple myeloma cell lines. It generates sample-specific variant, progression-specific, fusion-derived, novel-splicing-derived, and combined protein FASTA databases for downstream mass-spectrometry searches.

The current analysis focuses on longitudinal TK12, TK13, and TK14 samples. TK12 is the earliest available baseline; TK13 and TK14 are later progression samples from the same patient.

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

The workflow uses SRR accessions from `samples.csv`. SRA archives are downloaded and validated on the Saga login node, retained under `sra_cache/`, and converted into paired FASTQ files by Nextflow tasks on compute nodes.

### Proteomics

- PRIDE `PXD033531`: https://www.ebi.ac.uk/pride/archive/projects/PXD033531
- PRIDE `PXD033510`: https://www.ebi.ac.uk/pride/archive/projects/PXD033510

```bash
wget -r "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2024/11/PXD033531/"
wget -r "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2024/11/PXD033510/"
```

## Analysis Strategy

The workflow constructs exploratory proteogenomic databases. It is not intended for clinical somatic-mutation classification.

No matched normal sample is available, so Mutect2 tumor-normal calling is not used. GATK HaplotypeCaller creates a broad catalogue of expressed per-sample variants after RNA-seq-specific preprocessing. Calls may include germline, clonal, progression-associated, RNA-edited, and technical events and must not be interpreted as validated somatic mutations without additional evidence.

VEP retains consequences for all overlapping transcripts. Neither `--pick` nor `--flag_pick` is used. pypgatk receives the complete CSQ consequence set and uses Ensembl release 111 cDNA and the matching GTF to generate an expanded transcript-aware variant protein database.

Additional branches generate:

- Fusion-junction proteins from accepted Arriba candidates using pVACfuse
- Expressed transcript assemblies using StringTie
- gffcompare classification against Ensembl release 111
- ORF predictions for selected novel transcript classes using TransDecoder
- Progression VCF and protein FASTA outputs after baseline subtraction
- Deduplicated combined databases containing reviewed, variant, fusion, and splice-derived proteins

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
                       |      +--> transcript extraction
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

## Storage and Runtime Model

Project-local assets and published outputs remain under the repository directory:

```text
$PWD/
├── main.nf
├── scratch.slurm
├── validate_pipeline_commands.sh
├── samples.csv
├── singularity_cache/
├── reference_downloads/
├── sra_cache/
└── results/
```

The existing image directory keeps its historical name, `singularity_cache/`, but execution uses native Apptainer support.

Large execution data uses the Saga work filesystem:

```text
/cluster/work/users/ash022/work    Nextflow tasks, intermediates, and resume cache
/cluster/work/users/ash022/tmp     Launcher and Apptainer temporary files
```

The launcher configures:

```bash
export NXF_APPTAINER_CACHEDIR="${PROJECT_DIR}/singularity_cache"
export APPTAINER_TMPDIR="${RUN_TMP}/apptainer"
```

Nextflow is invoked with:

```bash
-with-apptainer
"-process.clusterOptions=--account=nn9036k"
-work-dir "${WORK_DIR}"
-resume
```

The launcher does not use legacy Singularity environment variables, message-level suppression, custom bind options, HTML execution reports, or timelines.

`fasterq-dump` uses `fasterq_tmp` inside its normal Nextflow task directory. No external fasterq-specific path is required.

## Important Implementation Details

- `fasterq-dump --split-files` creates paired `_1.fastq` and `_2.fastq` files.
- Trim Galore trims adapters and Phred-20 low-quality ends and removes pairs shorter than 36 bases.
- FastQC runs before and after trimming.
- STAR produces an Arriba-compatible unsorted BAM using `WithinBAM HardClip`.
- `SORT_INDEX_BAM` creates the coordinate-sorted and indexed BAM required by GATK and StringTie.
- `REF_INDEX` creates `genome.fa.fai` and `genome.dict` with SAMtools 1.21.
- `SplitNCigarReads` handles RNA-seq splice-junction alignments before variant calling.
- pypgatk 0.0.24 uses underscore-style options including `--input_fasta`, `--gene_annotations_gtf`, and `--output_proteindb`.
- pypgatk uses Ensembl release 111 cDNA; the reviewed UniProt proteome is the canonical base of the final combined database.
- `bcftools isec -C -w 1` retains progression variants absent from the matched baseline.
- gffcompare retains class codes `j` and `u` by default.
- TransDecoder utilities use verified absolute paths under `/usr/local/opt/transdecoder/util/`.
- Unpacked references and the STAR index remain in the Nextflow cache and are not duplicated into `results/`.
- Nextflow uses fail-fast behavior. Compatible completed tasks are recovered with `-resume`.

## GATK Temporary Storage and Index Handling

GATK processes use task-local temporary storage inside their Nextflow work directories:

```bash
mkdir -p gatk_tmp
trap 'rm -rf gatk_tmp' EXIT
```

Java is directed to the task-local directory with:

```text
-Djava.io.tmpdir=${PWD}/gatk_tmp
```

This configuration is used by `MARK_DUPLICATES`, `SPLIT_N_CIGAR`, `HAPLOTYPE_CALLER`, and all commands in `GENOTYPE_FILTER`. It prevents large GATK sorting spills from using restricted compute-node temporary storage.

GATK may create the `SplitNCigarReads` index as `<sample>.split.bai`, while the declared Nextflow output is `<sample>.split.bam.bai`. The workflow normalizes the filename after a successful run:

```bash
if [[ -s ${meta.sample}.split.bai && ! -e ${meta.sample}.split.bam.bai ]]; then
    mv ${meta.sample}.split.bai ${meta.sample}.split.bam.bai
fi

test -s ${meta.sample}.split.bam
test -s ${meta.sample}.split.bam.bai
```

## Robust Empty-Result Handling

The workflow distinguishes legitimate empty biological results from command failures:

- `PYPGATK_FASTA` creates an empty FASTA when pypgatk succeeds but generates no variant proteins.
- `PROGRESSION_FASTA` creates an empty FASTA when subtraction yields no progression proteins.
- `FUSION_FASTA` creates an empty FASTA when Arriba reports no accepted fusions.
- `SPLICE_PROTEIN_FASTA` creates an empty FASTA when no selected novel transcripts or proteins are available.
- `GFFCOMPARE_NOVEL` accepts an empty `novel.gtf`, requires a non-empty annotated GTF, and creates an explanatory statistics file when no non-empty `.stats` file is produced.

A genuine nonzero tool exit remains fatal because process scripts use `set -euo pipefail`.

## Exploratory Database Filters

```text
Minimum StringTie coverage:         2.5
Minimum junction-supporting reads:  3
Minimum isoform fraction:           0.05
Retained gffcompare class codes:    j,u
Minimum predicted protein length:   60 amino acids
```

## Full-Node Resource Strategy

Normal Saga nodes provide 20 CPUs and 80 GB. Compute-intensive normal-partition processes request up to 20 CPUs and 64 GB, leaving headroom for the operating system and container runtime.

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

`MARK_DUPLICATES` uses a bounded 56 GB Java heap and `--MAX_RECORDS_IN_RAM 1000000` after a confirmed out-of-memory failure at 16 GB. `STAR_ALIGN` remains on `bigmem` after an earlier exit-137 failure at 32 GB.

The Nextflow launcher requests 4 CPUs and 16 GB on `normal`.

## Requirements

- Nextflow 26.04.6
- Java 21
- Apptainer 1.4.4
- SLURM account `nn9036k`
- Internet access on the login node for initial downloads
- No compute-node internet dependency
- Fifteen pre-downloaded images under `singularity_cache/`

## Input Samplesheet

```csv
sample,srr,TK,Group,baseline
TK12,SRR31089074,patient1,resistant,true
TK13,SRR31089073,patient1,sensitive,false
TK14,SRR31089072,patient1,sensitive,false
```

Exactly one baseline must be marked `true` for each longitudinal key used for progression subtraction.

## Comprehensive Validation

The synchronized configuration passed comprehensive validation on Saga on July 24, 2026:

```text
PASS: 136
WARN: 0
FAIL: 0
RESULT: PASSED
```

Validation covers:

- Filesystem paths and real write tests
- Java, Nextflow, and Apptainer versions
- Samplesheet schema, unique samples, SRR syntax, baseline values, and local SRA paths
- Inspection of all 15 images and executable checks inside each container
- Gzip and tar integrity plus VEP and Arriba archive layouts
- HTML corruption, CRLF characters, and all 27 expected process declarations
- Command, option, reference-routing, channel, and resource contracts
- Absence of legacy, unsupported, network-download, and message-suppression options
- Empty-result guards
- GATK task-local temporary-storage contracts
- `SplitNCigarReads` BAM-index normalization and output checks
- A real gffcompare fixture executed inside its container
- Apptainer-native launcher configuration
- Launcher Bash syntax and `sbatch --test-only`

Run before submission:

```bash
dos2unix main.nf scratch.slurm validate_pipeline_commands.sh
chmod +x scratch.slurm validate_pipeline_commands.sh
bash -n scratch.slurm
bash -n validate_pipeline_commands.sh
bash validate_pipeline_commands.sh
```

The successful validation report is:

```text
pipeline_command_validation_20260724_202805.txt
```

Check that `main.nf` contains no HTML entities:

```bash
grep -nE '&gt;|&lt;|&amp;|-&gt;' main.nf
```

Expected result: no output.

## Setup on HPC/Saga

### 1. Download assets

```bash
bash download_assets.sh
bash download_sra.sh
```

Do not remove `singularity_cache/`, `reference_downloads/`, or `sra_cache/` when resetting an execution.

### 2. Prepare storage

```bash
mkdir -p \
    /cluster/work/users/ash022/work \
    /cluster/work/users/ash022/tmp

test -w /cluster/work/users/ash022/work
test -w /cluster/work/users/ash022/tmp
```

For a completely fresh execution only:

```bash
rm -rf /cluster/work/users/ash022/work/*
rm -rf /cluster/work/users/ash022/tmp/*
rm -rf .nextflow
rm -f .nextflow.log
```

### 3. Submit or resume

```bash
JOBID=$(sbatch --parsable scratch.slurm)
echo "$JOBID"

while [[ ! -f "resultsTKvep-${JOBID}.log" ]]; do
    sleep 2
done

tail -f "resultsTKvep-${JOBID}.log"
```

Do not delete `/cluster/work/users/ash022/work` before a resume run.

Monitor jobs and storage:

```bash
watch -n 60 '
squeue -u "$USER" -o "%.18i %.42j %.10P %.10T %.6C %.12m %.20R"
echo
du -sh /cluster/work/users/ash022/work 2>/dev/null
du -sh /cluster/work/users/ash022/tmp 2>/dev/null
'
```

## QC and Reporting

MultiQC aggregates raw and trimmed FastQC reports, Trim Galore and Cutadapt reports, STAR `Log.final.out`, `samtools flagstat`, MarkDuplicates metrics, and `bcftools stats`.

```text
results/multiqc/multiqc_report.html
results/multiqc/multiqc_data/
results/pipeline_trace-<launcher-job-id>.tsv
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

## Interpretation and Limitations

- RNA-seq detects variants and isoforms only in expressed and sufficiently covered regions.
- RNA editing, mapping artefacts, allele-specific expression, assembly errors, and incomplete ORFs can affect results.
- HaplotypeCaller outputs are not validated tumor-only somatic calls.
- TK12 subtraction can reflect biological acquisition or insufficient TK12 expression or coverage.
- Retaining all transcript consequences, fusion proteins, and splice-derived ORFs increases database size and multiple-testing burden.
- Strict target-decoy FDR control is required.
- Variant peptides should span the altered residue, fusion peptides the fusion junction, and splice-derived peptides a novel junction or ORF.
- Fusion and splice-derived sequences remain exploratory and require orthogonal RNA evidence.

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
