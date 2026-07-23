#!/bin/bash
#SBATCH --account=nn9036k
#SBATCH --job-name=resultsTKvep
#SBATCH --time=168:00:00
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --mail-user=animesh.sharma@ntnu.no
#SBATCH --mail-type=ALL
#SBATCH --output=resultsTKvep-%j.log

set -euo pipefail

PROJECT_ID="nn9036k"
NEXTFLOW="/cluster/home/ash022/scripts/nextflow"
WORKDIR="${SLURM_SUBMIT_DIR:-$(pwd -P)}"
RUN_ID="$SLURM_JOB_ID"
CONTAINER_DIR="$WORKDIR/singularity_cache"
REFERENCE_DIR="$WORKDIR/reference_downloads"
SRA_DIR="$WORKDIR/sra_cache"

cd "$WORKDIR"

export NXF_HOME="$WORKDIR/.nextflow_home"
export NXF_SINGULARITY_CACHEDIR="$CONTAINER_DIR"
export NXF_WORK="$WORKDIR/work"
export NXF_OPTS="-Xms1g -Xmx6g"
export SINGULARITY_TMPDIR="$WORKDIR/singularity_tmp"
export APPTAINER_TMPDIR="$SINGULARITY_TMPDIR"
export TMPDIR="$WORKDIR/tmp"

mkdir -p "$NXF_HOME" "$NXF_WORK" "$SINGULARITY_TMPDIR" "$TMPDIR" results

module purge
module load Java/21

command -v singularity
singularity --version
"$NEXTFLOW" -version

GENOME="$REFERENCE_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
GTF="$REFERENCE_DIR/Homo_sapiens.GRCh38.111.gtf.gz"
VEP="$REFERENCE_DIR/homo_sapiens_vep_111_GRCh38.tar.gz"
PROTEOME="$REFERENCE_DIR/human_reviewed_isoforms.fasta.gz"
ARRIBA="$REFERENCE_DIR/arriba_v2.4.0.tar.gz"

"$NEXTFLOW" run main.nf \
    -ansi-log false \
    -with-singularity \
    -process.executor slurm \
    "-process.clusterOptions=--account=$PROJECT_ID" \
    -work-dir "$NXF_WORK" \
    --sra_dir "$SRA_DIR" \
    --genome_url "file://$GENOME" \
    --gtf_url "file://$GTF" \
    --vep_cache_url "file://$VEP" \
    --proteome_url "file://$PROTEOME" \
    --arriba_url "file://$ARRIBA" \
    -with-report "results/pipeline_report-${RUN_ID}.html" \
    -with-timeline "results/pipeline_timeline-${RUN_ID}.html" \
    -with-trace "results/pipeline_trace-${RUN_ID}.tsv" \
    -resume
