# PGTK: Proteogenomics of Myeloma Cell Lines

## Overview

PGTK, ProteoGenomics TK, is a Nextflow DSL2 workflow for exploratory RNA-seq proteogenomics of patient-derived multiple myeloma cell lines. It generates sample-specific variant, fusion-derived, splice-derived, progression-specific, and combined protein FASTA databases for downstream mass-spectrometry searches. An optional MaxQuant evidence branch remaps identified peptides to the exact searched FASTA files and produces deduplicated variant and translated-junction evidence reports.

The current longitudinal analysis focuses on TK12, TK13, and TK14. TK12 is the earliest available baseline. TK13 and TK14 are later samples from the same patient.

This workflow supports exploratory research. It does not perform clinical somatic-mutation classification or establish clinical validity.

## Publication

The study methodology and biological context are described in:

**Multiple Myeloma Cells with Increased Proteasomal and ER Stress Are Hypersensitive to ATX-101, an Experimental Peptide Drug Targeting PCNA**

- Journal: Cancers
- Year: 2024
- Volume: 16
- Issue: 23
- Article: 3963
- DOI: https://doi.org/10.3390/cancers16233963
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

For the current validation, the MaxQuant text directory was:

```text
ftp.pride.ebi.ac.uk/pride/data/archive/2024/11/PXD033510/combined/txt
```

## Analysis Strategy

The workflow constructs exploratory proteogenomic search databases. No matched normal sample is available, so Mutect2 tumor-normal calling is not used. GATK HaplotypeCaller creates a broad catalogue of expressed per-sample variants after RNA-seq-specific preprocessing. Calls may include germline, clonal, progression-associated, RNA-edited, and technical events. They must not be interpreted as validated somatic mutations without independent evidence.

VEP retains consequences for all overlapping transcripts. Neither `--pick` nor `--flag_pick` is used. pypgatk receives the complete CSQ consequence set and uses the Ensembl release 111 cDNA FASTA and matching GTF to construct transcript-aware variant protein sequences.

Additional branches generate:

- Fusion-junction protein sequences from accepted Arriba candidates using pVACfuse
- Expressed transcript assemblies using StringTie
- gffcompare classification against Ensembl release 111
- ORF predictions for selected transcript classes using TransDecoder
- Progression VCF and protein FASTA outputs after subtraction of the TK12 baseline
- Deduplicated combined databases containing reviewed, variant, fusion, and splice-derived proteins
- Optional MaxQuant peptide remapping and deduplicated evidence reporting

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

optional MaxQuant branch
  +--> validate MaxQuant inputs and search provenance
  +--> remap peptides to every searched FASTA
  +--> annotate altered-residue peptide associations
  +--> analyze fusion and splice mappings
  +--> validate translated exon-junction crossings
  +--> build deduplicated evidence report
```

## Nextflow Processes

The validated workflow declares 33 processes:

```text
DOWNLOAD_REFERENCES
SRA_TO_FASTQ
CAT_FASTQ
FASTQC_RAW
TRIM_GALORE
FASTQC_TRIMMED
STAR_INDEX
STAR_ALIGN
SORT_INDEX_BAM
SAMTOOLS_FLAGSTAT
REF_INDEX
MARK_DUPLICATES
SPLIT_N_CIGAR
HAPLOTYPE_CALLER
GENOTYPE_FILTER
BCFTOOLS_STATS
VEP_ANNOTATE
PYPGATK_FASTA
ARRIBA
FUSION_FASTA
STRINGTIE_ASSEMBLY
GFFCOMPARE_NOVEL
SPLICE_PROTEIN_FASTA
COMBINE_PROTEIN_FASTA
PROGRESSION_SUBTRACT
PROGRESSION_FASTA
MULTIQC
VALIDATE_MAXQUANT_INPUTS
MAP_MAXQUANT_PEPTIDES
ANNOTATE_MAXQUANT_VARIANTS
ANALYZE_MAXQUANT_JUNCTIONS
VALIDATE_MAXQUANT_SPLICE_JUNCTIONS
BUILD_PROTEOGENOMICS_EVIDENCE_REPORT
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
- Consolidated MultiQC report
- Nextflow task trace for each launcher job
- Optional MaxQuant peptide-remapping and evidence-report outputs

## Storage and Runtime Model

Project-local source files, downloaded assets, and published outputs remain under the repository directory:

```text
$PWD/
├── main.nf
├── scratch.slurm
├── download_assets.sh
├── download_sra.sh
├── validate_pipeline_commands.sh
├── samples.csv
├── map_peptides_to_fasta.py
├── annotate_variant_peptides.py
├── analyze_chimeric_splice_peptides.py
├── validate_splice_junction_peptides.py
├── proteogenomics_evidence_report.py
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

Do not delete the work directory when a compatible `-resume` run is needed.

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
- `REF_INDEX` creates `genome.fa.fai` and `genome.dict` with SAMtools.
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
- The optional MaxQuant branch validates raw evidence tables and the recorded search FASTA paths before analysis.
- Transcript version suffixes are normalized only for joining. Full transcript annotations are retained in aggregated output columns.
- Variant outputs are deduplicated to one row per RNA sample, chromosome, position, REF, and ALT.
- Junction outputs are deduplicated by peptide sequence and genomic junction. RNA source samples and SRAs are aggregated separately from MS detection samples.

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

These filters apply to database construction. They are separate from MaxQuant evidence reporting.

## Resource Strategy

Normal Saga nodes provide 20 CPUs and 80 GB. Compute-intensive normal-partition tasks request up to 20 CPUs and 64 GB, leaving headroom for the operating system and container runtime. `STAR_ALIGN` uses `bigmem` with 32 CPUs and 256 GB.

`MARK_DUPLICATES` uses a bounded 56 GB Java heap and `--MAX_RECORDS_IN_RAM 1000000`. Both pypgatk processes request 80 GB task disk because each builds a large transcript annotation database and reads an uncompressed VCF.

The Nextflow launcher requests 4 CPUs and 16 GB on `normal`.

## Requirements

The final validated Saga configuration used:

- Nextflow 26.04.6
- Java 21.0.2
- Native Apptainer execution
- SLURM account `nn9036k`
- Internet access on the login node for initial downloads
- No compute-node internet dependency
- Fifteen pre-downloaded container images under `singularity_cache/`
- Seven validated reference archives under `reference_downloads/`

The validation environment reported Apptainer 1.4.4. A later compute-node launcher log reported 1.3.5, so exact patch-level availability can depend on the Saga node environment. The container and workflow contracts passed in both execution contexts.

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

### 1. Download and validate assets

```bash
bash download_assets.sh
bash download_sra.sh
```

The validated asset set contains 15 container images and seven reference archives. Do not remove `singularity_cache/`, `reference_downloads/`, or `sra_cache/` when resetting an execution.

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

The final synchronized configuration passed comprehensive validation on July 29, 2026:

```text
PASS: 168
WARN: 0
FAIL: 0
RESULT: PASSED
```

The final validation report was:

```text
pipeline_command_validation_20260729_185520.txt
```

Validation covers all 33 processes, all 15 images, local SRA files, seven references, command contracts, Python syntax, pypgatk VCF conversion, GATK temporary storage, split-index normalization, gffcompare fixture execution, MultiQC output contracts, optional MaxQuant contracts, native Apptainer settings, launcher forwarding, resume behavior, and `sbatch --test-only`.

### 4. Submit or resume the core workflow

```bash
JOBID=$(sbatch --parsable scratch.slurm)
echo "$JOBID"

while [[ ! -f "resultsTKvep-${JOBID}.log" ]]; do
    sleep 2
done

tail -f "resultsTKvep-${JOBID}.log"
```

### 5. Submit or resume with MaxQuant validation

```bash
MQTXT="$PWD/ftp.pride.ebi.ac.uk/pride/data/archive/2024/11/PXD033510/combined/txt"

JOBID=$(sbatch --parsable scratch.slurm \
    --run_proteogenomic_validation true \
    --maxquant_txt "$MQTXT" \
    --maxquant_mqpar "$PWD/mqpar.xml" \
    --maxquant_canonical_fasta \
        /cluster/home/ash022/FastaDB/uniprotkb_proteome_UP000005640_2026_06_25.fasta \
    --maxquant_contaminants \
        /cluster/home/ash022/scripts/MaxQuant_v2.8.1.0/bin/conf/contaminants.fasta)

echo "$JOBID"

while [[ ! -f "resultsTKvep-${JOBID}.log" ]]; do
    sleep 2
done

tail -f "resultsTKvep-${JOBID}.log"
```

The final evidence run was job `18955696`. It reused all compatible upstream tasks and reran `ANNOTATE_MAXQUANT_VARIANTS` and `BUILD_PROTEOGENOMICS_EVIDENCE_REPORT` after the reporting-script revision.

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
├── proteogenomics_validation/
└── pipeline_trace-*.tsv
```

The principal sample-specific search databases are:

```text
results/combined_fasta/TK12.exploratory_proteogenomics.fasta
results/combined_fasta/TK13.exploratory_proteogenomics.fasta
results/combined_fasta/TK14.exploratory_proteogenomics.fasta
```

## MaxQuant Search and Peptide Remapping

The generated databases were searched with MaxQuant together with the reviewed human UniProt proteome and the relevant contaminant database. Because MaxQuant groups proteins and can report selected leading identifiers, peptide-to-protein assignments must be remapped against every exact FASTA used in the search before noncanonical identifications are interpreted.

The optional evidence branch requires:

```text
peptides.txt
evidence.txt
msms.txt
proteinGroups.txt
mqpar.xml
canonical FASTA
contaminant FASTA
```

The observed TK detection pattern comes from MaxQuant evidence and raw-file metadata, not from the sample prefix of a FASTA header. Identical variant sequences may occur under TK12, TK13, and TK14 headers.

Use the exact canonical and contaminant FASTA files supplied to MaxQuant. A peptide absent from the reviewed human proteome can still be a contaminant peptide if the contaminant database is omitted during remapping.

The final search provenance was:

```text
MaxQuant version:          2.8.1.0
Match between runs:        False
Minimum peptide length:    7
Peptide FDR:               0.01
Protein FDR:               0.01
Contaminants enabled:      True
```

Peptide FDR is reported as search provenance. It is not converted into an arbitrary PEP cutoff. Score and PSM count are reported as evidence attributes and are not used as additional default thresholds.

## Local Proteogenomic Validation Scripts

### `map_peptides_to_fasta.py`

Maps every MaxQuant peptide back to all searched FASTA files, preserving complete headers and all matching coordinates.

Inputs:

- MaxQuant `peptides.txt`
- Reviewed human UniProt FASTA
- MaxQuant contaminant FASTA
- TK12, TK13, and TK14 combined FASTA files
- Experiment-to-sample mapping

Outputs:

```text
peptide_fasta_mapping.mapping.tsv
peptide_fasta_mapping.candidates.tsv
peptide_fasta_mapping.summary.txt
```

The mapping table contains every peptide. The candidate table is a convenience subset controlled by this script's explicit options. The final evidence report does not treat candidate-table thresholds as universal biological validation criteria.

Isoleucine and leucine are treated as equivalent by default.

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

Combines remapped peptides with VEP CSQ annotations, generated variant proteins, and the matching Ensembl release 111 protein FASTA.

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

For compatibility, `variant_peptide_annotation.prioritized.tsv` now means reference-absent altered-residue associations without analyst-defined PEP, score, PSM-count, or peptide-length thresholds. The filename is historical; the summary explicitly states the selection meaning.

### `analyze_chimeric_splice_peptides.py`

Screens MaxQuant peptides against fusion and splice FASTA files and the canonical proteome. It associates fusion mappings with accepted Arriba events and performs an initial split-anchor screen.

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
validated_splice_junctions.summary.txt
validated_splice_junctions.prioritized_novel_junctions.tsv
validated_splice_junctions.junction_spanning.tsv
validated_splice_junctions.detailed.tsv
validated_splice_junctions.unresolved.tsv
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
    --output-prefix validated_splice_junctions
```

### `proteogenomics_evidence_report.py`

Builds the final deduplicated evidence tables and Markdown report.

Variant rows are collapsed to one row per:

```text
RNA sample + chromosome + position + REF + ALT
```

Transcript-level VEP consequences are retained in aggregated columns. The report separates RNA source sample and SRA from the sample and raw file in which the peptide was detected by mass spectrometry.

The generated report contains five blocks:

```text
Block A: all altered-residue associations
Block B: search-consistent altered-residue evidence
Block C: sequence-novelty subsets
Block D: optional user-filtered associations
Block E: translated splice-junction evidence
```

Block A applies no analyst-defined PEP, score, PSM-count, or peptide-length threshold. The peptide must map to a generated variant protein and span the VEP-defined altered position.

Block B uses only search-derived and MaxQuant-output conditions:

- Minimum peptide length from `mqpar.xml`
- At least one MS/MS identification
- Not reverse/decoy
- Not marked as a potential contaminant by MaxQuant

Block C reports classification subsets:

- Canonical-absent altered-residue associations
- Ensembl-reference-absent altered-residue associations
- Canonical- and reference-absent altered-residue associations

These are sequence-novelty classifications, not confidence thresholds.

Block D is disabled by default. Optional user parameters are available when running the report script directly:

```text
--user-max-pep FLOAT
--user-min-score FLOAT
--user-min-msms-count INT
--user-min-peptide-length INT
--user-require-canonical-absence
--user-require-reference-absence
--user-exclude-contaminant-matches
--user-exclude-decoy-matches
```

Every enabled value is printed in the report as a user parameter.

Block E reports one row per peptide sequence and genomic junction. Equivalent sample-specific transcript mappings are aggregated instead of reported as duplicate biological junctions.

Outputs:

```text
proteogenomics_evidence.summary.txt
proteogenomics_evidence.report.md
proteogenomics_evidence.variants.tsv
proteogenomics_evidence.junctions.tsv
```

## Final Validation Output Directory

```text
results/proteogenomics_validation/
├── maxquant_inputs.validated.txt
├── peptide_fasta_mapping.mapping.tsv
├── peptide_fasta_mapping.candidates.tsv
├── peptide_fasta_mapping.summary.txt
├── variant_peptide_annotation.detailed.tsv
├── variant_peptide_annotation.prioritized.tsv
├── variant_peptide_annotation.summary.txt
├── variant_peptide_annotation.unresolved.tsv
├── junction_peptide_analysis.all_mappings.tsv
├── junction_peptide_analysis.fusion_candidates.tsv
├── junction_peptide_analysis.splice_candidates.tsv
├── junction_peptide_analysis.inferred_junctions.tsv
├── junction_peptide_analysis.summary.txt
├── validated_splice_junctions.detailed.tsv
├── validated_splice_junctions.junction_spanning.tsv
├── validated_splice_junctions.prioritized_novel_junctions.tsv
├── validated_splice_junctions.summary.txt
├── validated_splice_junctions.unresolved.tsv
├── proteogenomics_evidence.variants.tsv
├── proteogenomics_evidence.junctions.tsv
├── proteogenomics_evidence.summary.txt
└── proteogenomics_evidence.report.md
```


## Interpretation and Limitations

- RNA-seq detects variants and isoforms only in expressed and sufficiently covered regions.
- RNA editing, mapping artefacts, allele-specific expression, assembly errors, and incomplete ORFs can affect results.
- HaplotypeCaller outputs are not validated tumor-only somatic calls.
- TK12 subtraction can reflect biological acquisition or insufficient TK12 expression or coverage.
- Retaining all transcript consequences, fusion products, and splice-derived ORFs increases the search database and multiple-testing burden.
- Strict target-decoy FDR control is required during the MaxQuant search.
- A FASTA header or protein-group assignment is not direct evidence that a peptide contains an altered residue or junction.
- Variant peptides must overlap the VEP-defined altered position to enter the altered-residue blocks.
- Fusion peptides must cross the translated breakpoint to qualify as direct fusion-junction evidence.
- Splice-junction peptides must cross the translated exon boundary and be checked against the reference annotation.
- Peptides mapped only to a noncanonical ORF must be reported separately from junction-spanning peptides.
- MaxQuant leading protein identifiers may omit equivalent or shared FASTA mappings.
- Canonical absence and Ensembl-reference absence describe sequence novelty, not identification confidence.
- A peptide can be absent from one reference but present in another protein or isoform.
- RNA source sample and MS detection sample are different concepts and must not be conflated.
- Contaminants must be included during peptide remapping.
- The final report is exploratory. Spectrum-level review, RNA read-level inspection, and independent validation remain necessary before biological or clinical interpretation.

## Downstream Proteomics Recommendations

1. Search each sample against the corresponding combined exploratory FASTA plus canonical and contaminant sequences.
2. Use target-decoy searching and control peptide-spectrum match, peptide, and protein FDR.
3. Preserve `mqpar.xml`, all MaxQuant text tables, and the exact searched FASTA files.
4. Remap every identified peptide to every FASTA used in the search.
5. Separate canonical, altered-residue, sequence-novel, fusion, translated-junction, and noncanonical-ORF evidence classes.
6. Require variant peptides to span the altered residue before describing altered-residue peptide evidence.
7. Require fusion peptides to cross the translated fusion breakpoint.
8. Require splice-junction peptides to cross the translated exon boundary and compare the genomic junction with the reference annotation.
9. Compare peptide detections with RNA depth, allele fraction, junction support, and transcript abundance.
10. Inspect supporting spectra manually before reporting high-value events.
11. Collapse transcript duplicates and overlapping peptides supporting the same genomic event.
12. State every optional user-selected threshold and its source in the report.

## Reproducibility Checklist

Before interpreting or sharing results, retain:

```text
main.nf
scratch.slurm
samples.csv
mqpar.xml
all post-processing scripts
all exact searched FASTA files
MaxQuant txt directory
pipeline trace
pipeline validation report
final evidence report and TSV files
Nextflow work directory while resume compatibility is required
```

Recommended checks:

```bash
bash validate_pipeline_commands.sh

cat results/proteogenomics_validation/proteogenomics_evidence.summary.txt

grep -n '^## Block' \
    results/proteogenomics_validation/proteogenomics_evidence.report.md

wc -l \
    results/proteogenomics_validation/proteogenomics_evidence.variants.tsv \
    results/proteogenomics_validation/proteogenomics_evidence.junctions.tsv
```

## Citation

```text
Olaisen C, Røst LM, Sharma A, et al.
Multiple Myeloma Cells with Increased Proteasomal and ER Stress Are
Hypersensitive to ATX-101, an Experimental Peptide Drug Targeting PCNA.
Cancers. 2024;16(23):3963.
https://doi.org/10.3390/cancers16233963
```

## License

See [LICENSE](LICENSE) for details.

## Contact

For questions about the workflow or data, open an issue or contact `animesh@fuzzylife.org`.
