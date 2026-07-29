#!/bin/bash
#SBATCH --account=nn9036k
#SBATCH --job-name=resultsTKvep
#SBATCH --time=168:00:00
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-user=animesh.sharma@ntnu.no
#SBATCH --mail-type=ALL
#SBATCH --output=resultsTKvep-%j.log

set -euo pipefail

NEXTFLOW="/cluster/home/ash022/scripts/nextflow"
PROJECT_DIR="${SLURM_SUBMIT_DIR:-$(pwd -P)}"
WORK_DIR="/cluster/work/users/ash022/work"
TMP_ROOT="/cluster/work/users/ash022/tmp"
RUN_TMP="${TMP_ROOT}/pgtk-${SLURM_JOB_ID}"

cd "${PROJECT_DIR}"
module purge
module load Java/21

export NXF_HOME="${PROJECT_DIR}/.nextflow_home"
export NXF_APPTAINER_CACHEDIR="${PROJECT_DIR}/singularity_cache"
export NXF_OPTS="-Xms2g -Xmx12g"
export TMPDIR="${RUN_TMP}"
export TMP="${RUN_TMP}"
export TEMP="${RUN_TMP}"
export APPTAINER_TMPDIR="${RUN_TMP}/apptainer"

mkdir -p \
    "${WORK_DIR}" \
    "${TMPDIR}" \
    "${APPTAINER_TMPDIR}" \
    "${NXF_HOME}" \
    "${PROJECT_DIR}/results"

test -x "${NEXTFLOW}"
test -s "${PROJECT_DIR}/main.nf"
test -s "${PROJECT_DIR}/samples.csv"
test -d "${PROJECT_DIR}/singularity_cache"
test -d "${PROJECT_DIR}/reference_downloads"
test -d "${PROJECT_DIR}/sra_cache"
test -w "${WORK_DIR}"
test -w "${TMPDIR}"
command -v apptainer >/dev/null

printf 'Node: %s\n' "$(hostname)"
printf 'Project: %s\n' "${PROJECT_DIR}"
printf 'Work: %s\n' "${WORK_DIR}"
printf 'Temporary: %s\n' "${TMPDIR}"
df -h "${PROJECT_DIR}" "${WORK_DIR}" "${TMPDIR}"

apptainer --version
"${NEXTFLOW}" -version

"${NEXTFLOW}" run "${PROJECT_DIR}/main.nf" \
    -ansi-log false \
    -with-apptainer \
    -process.executor slurm \
    "-process.clusterOptions=--account=nn9036k" \
    -work-dir "${WORK_DIR}" \
    --sra_dir "${PROJECT_DIR}/sra_cache" \
    -with-trace "${PROJECT_DIR}/results/pipeline_trace-${SLURM_JOB_ID}.tsv" \
    -resume \
    "$@"
