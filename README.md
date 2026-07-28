# PGTK: Proteogenomics of Myeloma Cell Lines

## Overview

PGTK, ProteoGenomics TK, is a Nextflow DSL2 workflow for exploratory RNA-seq proteogenomics of patient-derived multiple myeloma cell lines. It generates sample-specific variant, fusion-derived, splice-derived, progression-specific, and combined protein FASTA databases for downstream mass-spectrometry searches.

The current longitudinal analysis focuses on TK12, TK13, and TK14. TK12 is treated as the earliest available baseline, while TK13 and TK14 represent later progression samples from the same patient.

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

The workflow constructs exploratory proteogenomic search databases. It is not intended for clinical somatic-mutation classification.

No matched normal sample is available, so Mutect2 tumor-normal calling is not used. GATK HaplotypeCaller instead creates a broad catalogue of expressed per-sample variants after RNA-seq-specific preprocessing. Calls may include germline, clonal, progression-associated, RNA-edited, and technical events. They must not be interpreted as validated somatic mutations without additional evidence.

VEP retains consequences for all overlapping transcripts. Neither `--pick` nor `--flag_pick` is used. pypgatk receives the complete CSQ consequence set and uses the Ensembl release 111 cDNA FASTA and matching GTF to construct transcript-aware variant protein sequences.

Additional branches generate:

- Fusion-junction protein sequences from accepted Arriba candidates using pVACfuse
- Expressed transcript assemblies using StringTie
- gffcompare classification against Ensembl release 111
- ORF predictions for selected transcript classes using TransDecoder
- Progression VCF and protein FASTA outputs after subtraction of the TK12 baseline
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
  +--> exact protein-sequence deduplication
  +--> combined exploratory proteogenomics FASTA

all QC logs and reports
  +--> MultiQC
```

## Principal Outputs

- Per-sample PASS and VEP-annotated VCF files
- Expanded per-sample variant protein FASTA files
- Accepted and discarded Arriba fusion tables
- Fusion-protein FASTA files
- StringTie transcript assemblies
- gffcompare annotated and novelty-filtered GTF files
- Splice-derived protein FASTA files
- TK13-minus-TK12 and TK14-minus-TK12 progression VCF files
- Progression-specific protein FASTA files
- Combined exploratory proteogenomics FASTA files
- Raw and trimmed FastQC reports
- Alignment, duplicate, and variant summary statistics
- A consolidated MultiQC report
- A Nextflow task trace for each launcher job

## Storage and Runtime Model

Project-local source files, downloaded assets, and published outputs remain under the repository directory:

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

The image directory retains its historical name, `singularity_cache/`, but execution uses native Apptainer support.

Large execution data is stored outside the repository:

```text
/cluster/work/users/ash022/work    Nextflow tasks, intermediates, and resume cache
/cluster/work/users/ash022/tmp     Launcher and Apptainer temporary files
```

The launcher configures:

```bash
export NXF_APPTAINER_CACHEDIR="${PROJECT_DIR}/singularity_cache"
export APPTAINER_TMPDIR="${RUN_TMP}/apptainer"
```

Nextflow is invoked with native Apptainer, SLURM execution, the external work directory, a task trace, and `-resume`. The launcher does not enable HTML execution reports or timelines by default.

## Important Implementation Details

- `fasterq-dump --split-files` creates paired `_1.fastq` and `_2.fastq` files.
- `fasterq-dump` uses `fasterq_tmp` inside its normal Nextflow task directory.
- Trim Galore removes adapters and Phred-20 low-quality ends, then removes pairs shorter than 36 bases.
- FastQC runs before and after trimming.
- STAR creates an Arriba-compatible unsorted BAM with `WithinBAM HardClip`.
- Arriba consumes the original STAR BAM directly.
- `SORT_INDEX_BAM` creates the coordinate-sorted and indexed BAM required by GATK and StringTie.
- `REF_INDEX` creates `genome.fa.fai` and `genome.dict` with SAMtools 1.21.
- `SplitNCigarReads` handles RNA-seq splice-junction alignments before variant calling.
- VEP retains all overlapping transcript consequences.
- pypgatk 0.0.24 uses underscore-style options, including `--input_fasta`, `--gene_annotations_gtf`, and `--output_proteindb`.
- pypgatk uses Ensembl release 111 cDNA. The reviewed UniProt proteome is the canonical base of each combined database.
- Both pypgatk branches validate and decompress BGZF VCF input to plain-text VCF before translation.
- Header-only VCFs produce an empty FASTA without invoking pypgatk.
- `bcftools isec -C -w 1` retains progression variants absent from the matched baseline.
- gffcompare retains class codes `j` and `u` by default.
- TransDecoder utilities use verified absolute paths under `/usr/local/opt/transdecoder/util/`.
- Compatible completed tasks are recovered with `-resume`.

## GATK Temporary Storage and Index Handling

GATK tasks use task-local temporary storage inside their Nextflow work directories:

```bash
mkdir -p gatk_tmp
trap 'rm -rf gatk_tmp' EXIT
```

Java is directed to the task-local directory with:

```text
-Djava.io.tmpdir=${PWD}/gatk_tmp
```

This is used by `MARK_DUPLICATES`, `SPLIT_N_CIGAR`, `HAPLOTYPE_CALLER`, and all commands in `GENOTYPE_FILTER`. It prevents large GATK sorting spills from using restricted compute-node temporary storage.

GATK can create the `SplitNCigarReads` index as `<sample>.split.bai`, while the declared Nextflow output is `<sample>.split.bam.bai`. The workflow normalizes and validates the filename:

```bash
if [[ -s ${meta.sample}.split.bai && ! -e ${meta.sample}.split.bam.bai ]]; then
    mv ${meta.sample}.split.bai ${meta.sample}.split.bam.bai
fi

test -s ${meta.sample}.split.bam
test -s ${meta.sample}.split.bam.bai
```

## Robust Empty-Result Handling

The workflow distinguishes legitimate empty biological results from command failures:

- `PYPGATK_FASTA` creates an empty FASTA when no qualifying variant proteins are generated.
- `PROGRESSION_FASTA` creates an empty FASTA when no progression proteins are generated.
- `FUSION_FASTA` creates an empty FASTA when Arriba reports no accepted fusions.
- `SPLICE_PROTEIN_FASTA` creates an empty FASTA when no selected transcript or predicted protein is available.
- `GFFCOMPARE_NOVEL` accepts an empty `novel.gtf`, requires a non-empty annotated GTF, and creates an explanatory statistics file if needed.

A genuine nonzero tool exit remains fatal because process scripts use `set -euo pipefail`.

## Exploratory Database Filters

```text
Minimum StringTie coverage:         2.5
Minimum junction-supporting reads:  3
Minimum isoform fraction:           0.05
Retained gffcompare class codes:    j,u
Minimum predicted protein length:   60 amino acids
```

## Resource Strategy

Normal Saga nodes provide 20 CPUs and 80 GB. Compute-intensive normal-partition tasks request up to 20 CPUs and 64 GB, leaving headroom for the operating system and container runtime. `STAR_ALIGN` uses `bigmem` with 32 CPUs and 256 GB.

`MARK_DUPLICATES` uses a bounded 56 GB Java heap and `--MAX_RECORDS_IN_RAM 1000000`. Both pypgatk processes request 80 GB task disk because each builds a large transcript annotation database and reads an uncompressed VCF.

The Nextflow launcher requests 4 CPUs and 16 GB on `normal`.

## Requirements

The validated Saga configuration uses:

- Nextflow 26.04.6
- Java 21
- Apptainer 1.4.4
- SLURM account `nn9036k`
- Internet access on the login node for initial downloads
- No compute-node internet dependency
- Fifteen pre-downloaded container images under `singularity_cache/`

The workflow includes containers for SRA Tools, Trim Galore, FastQC, STAR, GATK, SAMtools, VEP, pypgatk, Arriba, bcftools, MultiQC, StringTie, TransDecoder, pVACtools, and gffcompare.

## Input Samplesheet

```csv
sample,srr,TK,Group,baseline
TK12,SRR31089074,patient1,resistant,true
TK13,SRR31089073,patient1,sensitive,false
TK14,SRR31089072,patient1,sensitive,false
```

Exactly one baseline must be marked `true` for each longitudinal key used for progression subtraction.

## Setup on Saga

Run commands from the repository directory.

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

### 3. Validate

```bash
dos2unix main.nf scratch.slurm validate_pipeline_commands.sh
chmod +x scratch.slurm validate_pipeline_commands.sh
bash -n scratch.slurm
bash -n validate_pipeline_commands.sh
bash validate_pipeline_commands.sh
```

The synchronized configuration passed comprehensive validation on Saga on July 25, 2026:

```text
PASS: 146
WARN: 0
FAIL: 0
RESULT: PASSED
```

Validation covers all 27 processes, all 15 images, local SRA files, references, command contracts, pypgatk VCF conversion, GATK temporary storage, split-index normalization, a gffcompare fixture, the MultiQC output contract, native Apptainer settings, and `sbatch --test-only`.

### 4. Submit or resume

```bash
JOBID=$(sbatch --parsable scratch.slurm)
echo "$JOBID"

while [[ ! -f "resultsTKvep-${JOBID}.log" ]]; do
    sleep 2
done

tail -f "resultsTKvep-${JOBID}.log"
```

## Output Structure

```text
results/
├── qc/
├── multiqc/
│   ├── multiqc_report.html
│   └── multiqc_report_data/
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

The principal sample-specific search databases are:

```text
results/combined_fasta/TK12.exploratory_proteogenomics.fasta
results/combined_fasta/TK13.exploratory_proteogenomics.fasta
results/combined_fasta/TK14.exploratory_proteogenomics.fasta
```

## MaxQuant Search and Peptide Remapping

The generated databases were searched with MaxQuant together with the reviewed human UniProt proteome and the relevant contaminant database. Because MaxQuant groups proteins and can report only selected leading identifiers, peptide-to-protein assignments must be remapped against the exact FASTA files used in the search before interpreting noncanonical identifications.

The MaxQuant `txt/peptides.txt` table is the main input for the local post-processing tools described below. The observed TK pattern comes from the MaxQuant experiment-specific evidence, not from the sample prefix of a FASTA header. Identical variant sequences may occur under TK12, TK13, and TK14 headers.

Use the exact canonical and contaminant FASTA files supplied to MaxQuant. A peptide absent from the reviewed human proteome can still be a contaminant peptide if the contaminant database is omitted from remapping.

## Local Proteogenomic Validation Scripts

### `map_peptides_to_fasta.py`

Maps every MaxQuant peptide back to all searched FASTA files, preserving complete headers and all matching coordinates.

Inputs:

- MaxQuant `peptides.txt`
- Reviewed human UniProt FASTA
- MaxQuant contaminant FASTA, if used
- TK12, TK13, and TK14 combined FASTA files
- Experiment-to-sample mapping

Outputs:

```text
peptide_fasta_mapping.mapping.tsv
peptide_fasta_mapping.candidates.tsv
peptide_fasta_mapping.summary.txt
```

The candidate table retains peptides with no canonical match and at least one variant, splice, fusion, or progression match, subject to configurable quality thresholds. Isoleucine and leucine are treated as equivalent by default.

Example:

```bash
python map_peptides_to_fasta.py \
    --peptides /path/to/maxquant/txt/peptides.txt \
    --fasta \
        /path/to/reviewed_human.fasta \
        /path/to/maxquant_contaminants.fasta \
        results/combined_fasta/TK12.exploratory_proteogenomics.fasta \
        results/combined_fasta/TK13.exploratory_proteogenomics.fasta \
        results/combined_fasta/TK14.exploratory_proteogenomics.fasta \
    --group-map 2=TK12 \
    --group-map 3=TK13 \
    --group-map 4=TK14 \
    --output-prefix peptide_fasta_mapping
```

### `annotate_variant_peptides.py`

Combines remapped candidate peptides with VEP CSQ annotations, generated variant proteins, and the matching Ensembl release 111 protein FASTA.

It reports:

- Gene, transcript, and protein identifiers
- Genomic REF/ALT event
- Consequence, HGVSc, and HGVSp
- Reference and alternate amino acids
- Altered protein position
- Peptide position in the generated variant protein
- Whether the peptide spans the VEP-defined altered residue
- Whether the peptide occurs in the corresponding Ensembl reference protein

Outputs:

```text
variant_peptide_annotation.summary.txt
variant_peptide_annotation.prioritized.tsv
variant_peptide_annotation.detailed.tsv
variant_peptide_annotation.unresolved.tsv
```

The prioritized table contains peptides that overlap the VEP-defined protein change and are absent from the corresponding Ensembl reference protein.

### `analyze_chimeric_splice_peptides.py`

Screens MaxQuant peptides against fusion and splice FASTA files and the canonical proteome. It also associates fusion mappings with accepted Arriba events and performs an initial split-anchor screen.

Outputs:

```text
junction_peptide_analysis.summary.txt
junction_peptide_analysis.fusion_candidates.tsv
junction_peptide_analysis.splice_candidates.tsv
junction_peptide_analysis.inferred_junctions.tsv
junction_peptide_analysis.all_mappings.tsv
```

The split-anchor result is a screening classification. It does not replace exact translated breakpoint or exon-boundary reconstruction.

### `validate_splice_junction_peptides.py`

Performs definitive sample-aware validation of splice candidates. StringTie identifiers are indexed by `(sample, transcript_id)` so similarly named `STRG` transcripts from different samples are never merged.

The script:

- Uses the matching sample-specific assembled GTF for each splice protein
- Reconstructs transcript-oriented exon coordinates
- Converts exon boundaries into translated ORF positions
- Distinguishes genomic strand from TransDecoder ORF orientation
- Locates each peptide in the translated splice protein
- Determines whether the peptide crosses an exon boundary
- Reports amino-acid anchor lengths on both sides
- Tests whether the genomic junction exists in Ensembl release 111

Outputs:

```text
validated_splice_junctions_v3.summary.txt
validated_splice_junctions_v3.prioritized_novel_junctions.tsv
validated_splice_junctions_v3.junction_spanning.tsv
validated_splice_junctions_v3.detailed.tsv
validated_splice_junctions_v3.unresolved.tsv
```

Example:

```bash
python validate_splice_junction_peptides.py \
    --candidates junction_peptide_analysis.splice_candidates.tsv \
    --splice-fasta \
        results/splice_fasta/TK12.splice_proteins.fasta \
        results/splice_fasta/TK13.splice_proteins.fasta \
        results/splice_fasta/TK14.splice_proteins.fasta \
    --transcript-gtf \
        results/splicing/stringtie/TK12.assembled.gtf \
        results/splicing/stringtie/TK13.assembled.gtf \
        results/splicing/stringtie/TK14.assembled.gtf \
    --reference-gtf \
        reference_downloads/Homo_sapiens.GRCh38.111.gtf.gz \
    --output-prefix validated_splice_junctions_v3
```

## Interpretation and Limitations

- RNA-seq detects variants and isoforms only in expressed and sufficiently covered regions.
- RNA editing, mapping artefacts, allele-specific expression, assembly errors, and incomplete ORFs can affect results.
- HaplotypeCaller outputs are not validated tumor-only somatic calls.
- TK12 subtraction can reflect biological acquisition or insufficient TK12 expression or coverage.
- Retaining all transcript consequences, fusion products, and splice-derived ORFs increases the search database and multiple-testing burden.
- Strict target-decoy FDR control is required.
- A FASTA header or protein-group assignment is not direct evidence that a peptide contains the altered residue or junction.
- Variant peptides must overlap the VEP-defined altered residue.
- Fusion peptides must cross the translated breakpoint.
- Splice-junction peptides must cross the translated exon boundary and be checked against the reference annotation.
- Peptides mapped only to a novel ORF should be reported separately from junction-spanning peptides.
- MaxQuant leading protein identifiers may omit equivalent or shared FASTA mappings.
- Contaminants must be included during peptide remapping.
- Spectrum-level review and RNA-level confirmation remain necessary.

## Downstream Proteomics Recommendations

1. Search each sample against the corresponding combined exploratory FASTA plus canonical and contaminant sequences.
2. Use target-decoy searching and control peptide-spectrum match, peptide, and protein FDR.
3. Remap identified peptides to every FASTA used in the search.
4. Separate canonical, variant, fusion, splice-junction, and noncanonical-ORF evidence classes.
5. Require variant peptides to span the altered residue.
6. Require fusion peptides to cross the fusion breakpoint.
7. Require splice-junction peptides to cross the translated exon boundary.
8. Compare peptide detections with RNA depth, allele fraction, junction support, and transcript abundance.
9. Inspect supporting spectra manually before reporting high-value events.
10. Collapse overlapping peptides supporting the same genomic event.

## Citation

```text
Cancers. 2024;16(23):3963.
https://www.mdpi.com/2072-6694/16/23/3963
```

## License

See [LICENSE](LICENSE) for details.

## Contact

For questions about the workflow or data, open an issue or contact `animesh@fuzzylife.org`.
