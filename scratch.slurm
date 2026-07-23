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
RUN_ID="${SLURM_JOB_ID}"

CONTAINER_DIR="${WORKDIR}/singularity_cache"
REFERENCE_DIR="${WORKDIR}/reference_downloads"
SRA_DIR="${WORKDIR}/sra_cache"
SCRATCH_ROOT="/cluster/work/users/ash022/tmp/pgtk"
NXF_WORK="${SCRATCH_ROOT}/work"
NXF_TMP="${SCRATCH_ROOT}/tmp/launcher-${RUN_ID}"

cd "${WORKDIR}"

module purge
module load Java/21

export NXF_HOME="${WORKDIR}/.nextflow_home"
export NXF_SINGULARITY_CACHEDIR="${CONTAINER_DIR}"
export NXF_OPTS="-Xms1g -Xmx6g"
export TMPDIR="${NXF_TMP}"
export TMP="${NXF_TMP}"
export TEMP="${NXF_TMP}"
export SINGULARITY_TMPDIR="${NXF_TMP}/singularity"
export APPTAINER_TMPDIR="${NXF_TMP}/apptainer"

mkdir -p \
    "${NXF_HOME}" \
    "${NXF_WORK}" \
    "${TMPDIR}" \
    "${SINGULARITY_TMPDIR}" \
    "${APPTAINER_TMPDIR}" \
    "${WORKDIR}/results"

command -v singularity >/dev/null
command -v java >/dev/null
test -x "${NEXTFLOW}"
test -s "${WORKDIR}/main.nf"
test -s "${WORKDIR}/samples.csv"
test -d "${CONTAINER_DIR}"
test -d "${REFERENCE_DIR}"
test -d "${SRA_DIR}"
test -d "${SCRATCH_ROOT}"
test -w "${SCRATCH_ROOT}"
test -w "${NXF_WORK}"
test -w "${TMPDIR}"

printf 'Node: %s\n' "$(hostname)"
printf 'Started: %s\n' "$(date --iso-8601=seconds)"
printf 'Project directory: %s\n' "${WORKDIR}"
printf 'Nextflow work directory: %s\n' "${NXF_WORK}"
printf 'Temporary directory: %s\n' "${NXF_TMP}"

df -h "${WORKDIR}" "${SCRATCH_ROOT}"
df -ih "${WORKDIR}" "${SCRATCH_ROOT}"

command -v singularity
singularity --version
"${NEXTFLOW}" -version

"${NEXTFLOW}" run "${WORKDIR}/main.nf" \
    -ansi-log false \
    -with-singularity \
    -process.executor slurm \
    "-process.clusterOptions=--account=${PROJECT_ID}" \
    -work-dir "${NXF_WORK}" \
    --sra_dir "${SRA_DIR}" \
    -with-report "${WORKDIR}/results/pipeline_report-${RUN_ID}.html" \
    -with-timeline "${WORKDIR}/results/pipeline_timeline-${RUN_ID}.html" \
    -with-trace "${WORKDIR}/results/pipeline_trace-${RUN_ID}.tsv"

printf 'Finished: %s\n' "$(date --iso-8601=seconds)"
#dos2unix scratch.slurm validate_pipeline_commands.sh main.nf
#bash validate_pipeline_commands.sh
#sbatch --test-only scratch.slurm
#JOBID=$(sbatch --parsable scratch.slurm)
#echo "$JOBID"
#while [[ ! -f "resultsTKvep-${JOBID}.log" ]]; do sleep 2 done
#tail -f "resultsTKvep-${JOBID}.log"