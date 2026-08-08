#!/usr/bin/env bash
#SBATCH --account=nn9036k
#SBATCH --partition=normal
#SBATCH --job-name=pgtk
#SBATCH --time=168:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/cluster/projects/nn9036k/scrbkup/pgtk/pgtk-wrapper-%j.log
set -uo pipefail

usage() {
    cat <<'EOF'
Usage: sbatch [SBATCH OPTIONS] scratch.slurm [WRAPPER OPTIONS] -- [PIPELINE OPTIONS]

Optional wrapper overrides (Saga defaults are built in):
  --project-dir PATH
  --nextflow PATH
  --python PATH
  --apptainer PATH
  --java-module NAME
  --work-dir PATH
  --tmp-root PATH
  --nxf-home PATH
  --container-cache PATH
  --sra-dir PATH
  --reference-downloads PATH
  --samplesheet PATH
  --results-dir PATH
  --ensembl-pep PATH
  --profile NAME
  --account NAME
  --normal-partition NAME
  --bigmem-partition NAME
  --normal-cpu-threshold INTEGER
  --normal-memory-threshold-gb INTEGER
  --max-cpus INTEGER
  --max-memory-gb INTEGER
  --queue-size INTEGER
  --submit-rate-limit VALUE
  --nxf-opts VALUE
  --slurm-log-template PATH_WITH_{job_id}

Everything following -- is forwarded unchanged to main.nf.
EOF
}

die() { printf 'ERROR: %s\n' "$*" >&2; exit 2; }
need_value() { [[ $# -ge 2 && -n "$2" ]] || die "$1 requires a value"; }

PROJECT_DIR=/cluster/projects/nn9036k/scrbkup/pgtk
NEXTFLOW=/cluster/home/ash022/bin/nextflow
HOST_PYTHON=/cluster/software/Mamba/4.14.0-0/bin/python
APPTAINER=$(command -v apptainer || true)
JAVA_MODULE=Java/21
WORK_DIR=/cluster/work/users/ash022/work
TMP_ROOT=/cluster/work/users/ash022/tmp
NXF_HOME_DIR=/cluster/projects/nn9036k/scrbkup/pgtk/.nextflow_home
CONTAINER_CACHE=/cluster/projects/nn9036k/scrbkup/pgtk/singularity_cache
SRA_DIR=/cluster/projects/nn9036k/scrbkup/pgtk/sra_cache
REFERENCE_DOWNLOADS=/cluster/projects/nn9036k/scrbkup/pgtk/reference_downloads
SAMPLESHEET=/cluster/projects/nn9036k/scrbkup/pgtk/samples.csv
RESULTS_DIR=/cluster/projects/nn9036k/scrbkup/pgtk/results
ENSEMBL_PEP=/cluster/projects/nn9036k/scrbkup/pgtk/reference_downloads/Homo_sapiens.GRCh38.pep.all.fa.gz
PROFILE=robust
ACCOUNT=nn9036k
NORMAL_PARTITION=normal
BIGMEM_PARTITION=bigmem
NORMAL_CPU_THRESHOLD=20
NORMAL_MEMORY_THRESHOLD_GB=160
MAX_CPUS=32
MAX_MEMORY_GB=512
QUEUE_SIZE=200
SUBMIT_RATE_LIMIT=60/1min
NXF_OPTS_VALUE='-Xms1g -Xmx6g'
SLURM_LOG_TEMPLATE=/cluster/projects/nn9036k/scrbkup/pgtk/pgtk-wrapper-{job_id}.log
PIPELINE_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --project-dir) need_value "$@"; PROJECT_DIR=$2; shift 2 ;;
        --nextflow) need_value "$@"; NEXTFLOW=$2; shift 2 ;;
        --python) need_value "$@"; HOST_PYTHON=$2; shift 2 ;;
        --apptainer) need_value "$@"; APPTAINER=$2; shift 2 ;;
        --java-module) need_value "$@"; JAVA_MODULE=$2; shift 2 ;;
        --work-dir) need_value "$@"; WORK_DIR=$2; shift 2 ;;
        --tmp-root) need_value "$@"; TMP_ROOT=$2; shift 2 ;;
        --nxf-home) need_value "$@"; NXF_HOME_DIR=$2; shift 2 ;;
        --container-cache) need_value "$@"; CONTAINER_CACHE=$2; shift 2 ;;
        --sra-dir) need_value "$@"; SRA_DIR=$2; shift 2 ;;
        --reference-downloads) need_value "$@"; REFERENCE_DOWNLOADS=$2; shift 2 ;;
        --samplesheet) need_value "$@"; SAMPLESHEET=$2; shift 2 ;;
        --results-dir) need_value "$@"; RESULTS_DIR=$2; shift 2 ;;
        --ensembl-pep) need_value "$@"; ENSEMBL_PEP=$2; shift 2 ;;
        --profile) need_value "$@"; PROFILE=$2; shift 2 ;;
        --account) need_value "$@"; ACCOUNT=$2; shift 2 ;;
        --normal-partition) need_value "$@"; NORMAL_PARTITION=$2; shift 2 ;;
        --bigmem-partition) need_value "$@"; BIGMEM_PARTITION=$2; shift 2 ;;
        --normal-cpu-threshold) need_value "$@"; NORMAL_CPU_THRESHOLD=$2; shift 2 ;;
        --normal-memory-threshold-gb) need_value "$@"; NORMAL_MEMORY_THRESHOLD_GB=$2; shift 2 ;;
        --max-cpus) need_value "$@"; MAX_CPUS=$2; shift 2 ;;
        --max-memory-gb) need_value "$@"; MAX_MEMORY_GB=$2; shift 2 ;;
        --queue-size) need_value "$@"; QUEUE_SIZE=$2; shift 2 ;;
        --submit-rate-limit) need_value "$@"; SUBMIT_RATE_LIMIT=$2; shift 2 ;;
        --nxf-opts) need_value "$@"; NXF_OPTS_VALUE=$2; shift 2 ;;
        --slurm-log-template) need_value "$@"; SLURM_LOG_TEMPLATE=$2; shift 2 ;;
        --) shift; PIPELINE_ARGS=("$@"); break ;;
        --help|-h) usage; exit 0 ;;
        *) die "unknown wrapper option before --: $1" ;;
    esac
done

for variable in PROJECT_DIR NEXTFLOW HOST_PYTHON APPTAINER JAVA_MODULE WORK_DIR TMP_ROOT NXF_HOME_DIR CONTAINER_CACHE SRA_DIR REFERENCE_DOWNLOADS SAMPLESHEET RESULTS_DIR ENSEMBL_PEP PROFILE ACCOUNT NORMAL_PARTITION BIGMEM_PARTITION NORMAL_CPU_THRESHOLD NORMAL_MEMORY_THRESHOLD_GB MAX_CPUS MAX_MEMORY_GB QUEUE_SIZE SUBMIT_RATE_LIMIT NXF_OPTS_VALUE SLURM_LOG_TEMPLATE; do
    [[ -n ${!variable} ]] || die "missing required wrapper option for $variable"
done
for variable in NORMAL_CPU_THRESHOLD NORMAL_MEMORY_THRESHOLD_GB MAX_CPUS MAX_MEMORY_GB QUEUE_SIZE; do
    [[ ${!variable} =~ ^[1-9][0-9]*$ ]] || die "$variable must be a positive integer"
done
(( NORMAL_CPU_THRESHOLD < MAX_CPUS )) || die 'normal CPU threshold must be below maximum CPUs'
(( NORMAL_MEMORY_THRESHOLD_GB < MAX_MEMORY_GB )) || die 'normal memory threshold must be below maximum memory'
[[ $SLURM_LOG_TEMPLATE == *'{job_id}'* ]] || die '--slurm-log-template must contain {job_id}'

JOB_ID=${SLURM_JOB_ID:-manual-$(date +%Y%m%d_%H%M%S)}
RUN_TMP="$TMP_ROOT/pgtk-$JOB_ID"
TRACE_FILE="$RESULTS_DIR/pipeline_trace-$JOB_ID.tsv"
RESOURCE_PREFIX="$RESULTS_DIR/resource_usage-$JOB_ID"
NEXTFLOW_LOG="$PROJECT_DIR/.nextflow.log"
SLURM_LOG=${SLURM_LOG_TEMPLATE//\{job_id\}/$JOB_ID}

finalize() {
    code=$?
    trap - EXIT
    set +e
    if [[ -s $TRACE_FILE ]]; then
        "$HOST_PYTHON" "$PROJECT_DIR/analyze_pipeline_trace.py" --trace "$TRACE_FILE" --output-prefix "$RESOURCE_PREFIX"
    fi
    "$HOST_PYTHON" "$PROJECT_DIR/collect_pipeline_failures.py" \
        --job-id "$JOB_ID" --exit-code "$code" --trace "$TRACE_FILE" \
        --nextflow-log "$NEXTFLOW_LOG" --slurm-log "$SLURM_LOG" --results-dir "$RESULTS_DIR"
    rm -rf "$RUN_TMP"
    exit "$code"
}
trap finalize EXIT

for executable in "$NEXTFLOW" "$HOST_PYTHON" "$APPTAINER"; do [[ -x $executable ]] || die "not executable: $executable"; done
for file in "$PROJECT_DIR/main.nf" "$PROJECT_DIR/nextflow.config" "$SAMPLESHEET" "$PROJECT_DIR/collect_pipeline_failures.py" "$PROJECT_DIR/analyze_pipeline_trace.py" "$ENSEMBL_PEP"; do [[ -s $file ]] || die "missing or empty file: $file"; done
for directory in "$SRA_DIR" "$REFERENCE_DOWNLOADS"; do [[ -d $directory && -r $directory ]] || die "input directory is missing or unreadable: $directory"; done
for directory in "$WORK_DIR" "$TMP_ROOT" "$NXF_HOME_DIR" "$CONTAINER_CACHE" "$RESULTS_DIR"; do mkdir -p "$directory" || die "cannot create directory: $directory"; [[ -d $directory && -w $directory ]] || die "directory is not writable: $directory"; done

cd "$PROJECT_DIR" || die "cannot enter project directory"
module purge || die 'module purge failed'
module load "$JAVA_MODULE" || die "cannot load Java module: $JAVA_MODULE"
export PATH="$(dirname "$APPTAINER"):$PATH"
export NXF_HOME="$NXF_HOME_DIR"
export NXF_APPTAINER_CACHEDIR="$CONTAINER_CACHE"
export NXF_OPTS="$NXF_OPTS_VALUE"
export TMPDIR="$RUN_TMP" TMP="$RUN_TMP" TEMP="$RUN_TMP" APPTAINER_TMPDIR="$RUN_TMP/apptainer"
export PGTK_ACCOUNT="$ACCOUNT"
export PGTK_NORMAL_PARTITION="$NORMAL_PARTITION"
export PGTK_BIGMEM_PARTITION="$BIGMEM_PARTITION"
export PGTK_NORMAL_CPU_THRESHOLD="$NORMAL_CPU_THRESHOLD"
export PGTK_NORMAL_MEMORY_THRESHOLD_GB="$NORMAL_MEMORY_THRESHOLD_GB"
export PGTK_MAX_CPUS="$MAX_CPUS"
export PGTK_MAX_MEMORY_GB="$MAX_MEMORY_GB"
export PGTK_QUEUE_SIZE="$QUEUE_SIZE"
export PGTK_SUBMIT_RATE_LIMIT="$SUBMIT_RATE_LIMIT"
mkdir -p "$TMPDIR" "$APPTAINER_TMPDIR" "$RESULTS_DIR/failure_logs" || die 'cannot create runtime directories'
HOST_PYTHON_REAL=$(readlink -f "$HOST_PYTHON") || die "cannot resolve host Python: $HOST_PYTHON"
HOST_PYTHON_ROOT=$(dirname "$(dirname "$HOST_PYTHON_REAL")")
RUNTIME_BINDS="$PROJECT_DIR,$REFERENCE_DOWNLOADS,$SRA_DIR,$CONTAINER_CACHE,$WORK_DIR,$TMP_ROOT,$HOST_PYTHON_ROOT"
export APPTAINER_BINDPATH="${APPTAINER_BINDPATH:+$APPTAINER_BINDPATH,}$RUNTIME_BINDS"
"$HOST_PYTHON" "$PROJECT_DIR/validate_runtime_inputs.py" \
    --project-dir "$PROJECT_DIR" \
    --reference-downloads "$REFERENCE_DOWNLOADS" \
    --container-cache "$CONTAINER_CACHE" \
    --sra-dir "$SRA_DIR" \
    --samplesheet "$SAMPLESHEET" \
    --ensembl-pep "$ENSEMBL_PEP" \
    --work-dir "$WORK_DIR" \
    --tmp-root "$TMP_ROOT" \
    --results-dir "$RESULTS_DIR" \
    --host-python "$HOST_PYTHON" \
    --apptainer "$APPTAINER" || die 'runtime input validation failed'

printf 'Job: %s\nProject: %s\nWork: %s\nTemporary: %s\nProfile: %s\nNormal partition: %s\nBigmem partition: %s\nThresholds: CPU>%s or memory>%sGB\nCaps: CPU=%s memory=%sGB\n' \
    "$JOB_ID" "$PROJECT_DIR" "$WORK_DIR" "$RUN_TMP" "$PROFILE" "$NORMAL_PARTITION" "$BIGMEM_PARTITION" "$NORMAL_CPU_THRESHOLD" "$NORMAL_MEMORY_THRESHOLD_GB" "$MAX_CPUS" "$MAX_MEMORY_GB"
"$APPTAINER" --version
"$NEXTFLOW" -version

"$NEXTFLOW" run "$PROJECT_DIR/main.nf" \
    -profile "$PROFILE" -ansi-log false -with-apptainer \
    -work-dir "$WORK_DIR" \
    -with-trace "$TRACE_FILE" \
    -with-report "$RESULTS_DIR/pipeline_report-$JOB_ID.html" \
    -with-timeline "$RESULTS_DIR/pipeline_timeline-$JOB_ID.html" \
    -with-dag "$RESULTS_DIR/pipeline_dag-$JOB_ID.html" \
    -resume \
    --samplesheet "$SAMPLESHEET" \
    --outdir "$RESULTS_DIR" \
    --sra_dir "$SRA_DIR" \
    --reference_downloads "$REFERENCE_DOWNLOADS" \
    --container_cache "$CONTAINER_CACHE" \
    --ensembl_pep "$ENSEMBL_PEP" \
    --host_python "$HOST_PYTHON" \
    "${PIPELINE_ARGS[@]}"
