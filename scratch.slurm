#!/usr/bin/env bash
#SBATCH --job-name=pgtk
#SBATCH --time=168:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=pgtk-wrapper-%j.log
set -uo pipefail

usage() {
cat <<'USAGE'
Usage:
  sbatch --account ACCOUNT --partition PARTITION scratch.slurm \
    --account ACCOUNT --normal-partition PARTITION --bigmem-partition PARTITION [OPTIONS] -- [NEXTFLOW PARAMS]

Wrapper options:
  --account NAME                    Required unless PGTK_ACCOUNT is set
  --normal-partition NAME           Required unless PGTK_NORMAL_PARTITION is set
  --bigmem-partition NAME           Required unless PGTK_BIGMEM_PARTITION is set
  --project-dir PATH                Default: SLURM_SUBMIT_DIR or current directory
  --nextflow PATH                   Default: nextflow found on PATH
  --python PATH                     Default: python3 found on PATH
  --apptainer PATH                  Default: apptainer found on PATH
  --java-module NAME                Optional; no module is loaded when omitted
  --work-dir PATH                   Default: PROJECT_DIR/work
  --tmp-root PATH                   Default: TMPDIR or PROJECT_DIR/tmp
  --nxf-home PATH                   Default: PROJECT_DIR/.nextflow_home
  --container-cache PATH            Default: PROJECT_DIR/singularity_cache
  --pysam-image PATH                Default: pinned image below container cache
  --sra-dir PATH                    Default: PROJECT_DIR/sra_cache
  --reference-downloads PATH        Default: PROJECT_DIR/reference_downloads
  --samplesheet PATH                Default: PROJECT_DIR/samples.csv
  --results-dir PATH                Default: PROJECT_DIR/results
  --ensembl-pep PATH                Default: under reference downloads
  --profile NAME                    Default: robust
  --normal-cpu-threshold INTEGER    Default: 20
  --normal-memory-threshold-gb INT  Default: 160
  --max-cpus INTEGER                Default: 32
  --max-memory-gb INTEGER           Default: 512
  --queue-size INTEGER              Default: 200
  --submit-rate-limit VALUE         Default: 60/1min
  --nxf-opts VALUE                  Default: -Xms1g -Xmx6g
  --slurm-log-template TEMPLATE     Default: PROJECT_DIR/pgtk-wrapper-{job_id}.log

Values may also be supplied through matching PGTK_* environment variables.
Everything after -- is forwarded unchanged to main.nf and overrides wrapper defaults.
USAGE
}

die(){ printf 'ERROR: %s\n' "$*" >&2; exit 2; }
need_value(){ [[ $# -ge 2 && -n $2 ]] || die "$1 requires a value"; }
PROJECT_DIR=${PGTK_PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$(pwd -P)}}
NEXTFLOW=${PGTK_NEXTFLOW:-$(command -v nextflow || true)}
HOST_PYTHON=${PGTK_PYTHON:-$(command -v python3 || true)}
APPTAINER=${PGTK_APPTAINER:-$(command -v apptainer || true)}
JAVA_MODULE=${PGTK_JAVA_MODULE:-}
WORK_DIR=${PGTK_WORK_DIR:-$PROJECT_DIR/work}
TMP_ROOT=${PGTK_TMP_ROOT:-${TMPDIR:-$PROJECT_DIR/tmp}}
NXF_HOME_DIR=${PGTK_NXF_HOME:-$PROJECT_DIR/.nextflow_home}
CONTAINER_CACHE=${PGTK_CONTAINER_CACHE:-$PROJECT_DIR/singularity_cache}
PYSAM_IMAGE=${PGTK_PYSAM_IMAGE:-}
SRA_DIR=${PGTK_SRA_DIR:-$PROJECT_DIR/sra_cache}
REFERENCE_DOWNLOADS=${PGTK_REFERENCE_DOWNLOADS:-$PROJECT_DIR/reference_downloads}
SAMPLESHEET=${PGTK_SAMPLESHEET:-$PROJECT_DIR/samples.csv}
RESULTS_DIR=${PGTK_RESULTS_DIR:-$PROJECT_DIR/results}
ENSEMBL_PEP=${PGTK_ENSEMBL_PEP:-}
PROFILE=${PGTK_PROFILE:-robust}
ACCOUNT=${PGTK_ACCOUNT:-}
NORMAL_PARTITION=${PGTK_NORMAL_PARTITION:-}
BIGMEM_PARTITION=${PGTK_BIGMEM_PARTITION:-}
NORMAL_CPU_THRESHOLD=${PGTK_NORMAL_CPU_THRESHOLD:-20}
NORMAL_MEMORY_THRESHOLD_GB=${PGTK_NORMAL_MEMORY_THRESHOLD_GB:-160}
MAX_CPUS=${PGTK_MAX_CPUS:-32}
MAX_MEMORY_GB=${PGTK_MAX_MEMORY_GB:-512}
QUEUE_SIZE=${PGTK_QUEUE_SIZE:-200}
SUBMIT_RATE_LIMIT=${PGTK_SUBMIT_RATE_LIMIT:-60/1min}
NXF_OPTS_VALUE=${PGTK_NXF_OPTS:-'-Xms1g -Xmx6g'}
SLURM_LOG_TEMPLATE=${PGTK_SLURM_LOG_TEMPLATE:-$PROJECT_DIR/pgtk-wrapper-{job_id}.log}
PIPELINE_ARGS=()
while [[ $# -gt 0 ]]; do
 case "$1" in
  --project-dir) need_value "$@"; PROJECT_DIR=$2; shift 2;; --nextflow) need_value "$@"; NEXTFLOW=$2; shift 2;;
  --python) need_value "$@"; HOST_PYTHON=$2; shift 2;; --apptainer) need_value "$@"; APPTAINER=$2; shift 2;;
  --java-module) need_value "$@"; JAVA_MODULE=$2; shift 2;; --work-dir) need_value "$@"; WORK_DIR=$2; shift 2;;
  --tmp-root) need_value "$@"; TMP_ROOT=$2; shift 2;; --nxf-home) need_value "$@"; NXF_HOME_DIR=$2; shift 2;;
  --container-cache) need_value "$@"; CONTAINER_CACHE=$2; shift 2;; --pysam-image) need_value "$@"; PYSAM_IMAGE=$2; shift 2;;
  --sra-dir) need_value "$@"; SRA_DIR=$2; shift 2;; --reference-downloads) need_value "$@"; REFERENCE_DOWNLOADS=$2; shift 2;;
  --samplesheet) need_value "$@"; SAMPLESHEET=$2; shift 2;; --results-dir) need_value "$@"; RESULTS_DIR=$2; shift 2;;
  --ensembl-pep) need_value "$@"; ENSEMBL_PEP=$2; shift 2;; --profile) need_value "$@"; PROFILE=$2; shift 2;;
  --account) need_value "$@"; ACCOUNT=$2; shift 2;; --normal-partition) need_value "$@"; NORMAL_PARTITION=$2; shift 2;;
  --bigmem-partition) need_value "$@"; BIGMEM_PARTITION=$2; shift 2;; --normal-cpu-threshold) need_value "$@"; NORMAL_CPU_THRESHOLD=$2; shift 2;;
  --normal-memory-threshold-gb) need_value "$@"; NORMAL_MEMORY_THRESHOLD_GB=$2; shift 2;; --max-cpus) need_value "$@"; MAX_CPUS=$2; shift 2;;
  --max-memory-gb) need_value "$@"; MAX_MEMORY_GB=$2; shift 2;; --queue-size) need_value "$@"; QUEUE_SIZE=$2; shift 2;;
  --submit-rate-limit) need_value "$@"; SUBMIT_RATE_LIMIT=$2; shift 2;; --nxf-opts) need_value "$@"; NXF_OPTS_VALUE=$2; shift 2;;
  --slurm-log-template) need_value "$@"; SLURM_LOG_TEMPLATE=$2; shift 2;; --) shift; PIPELINE_ARGS=("$@"); break;;
  --help|-h) usage; exit 0;; *) die "unknown wrapper option before --: $1";;
 esac
done
[[ -n $PYSAM_IMAGE ]] || PYSAM_IMAGE=$CONTAINER_CACHE/quay.io-biocontainers-pysam-0.24.0--py312hf5ad864_1.img
[[ -n $ENSEMBL_PEP ]] || ENSEMBL_PEP=$REFERENCE_DOWNLOADS/Homo_sapiens.GRCh38.pep.all.fa.gz
for v in PROJECT_DIR NEXTFLOW HOST_PYTHON APPTAINER WORK_DIR TMP_ROOT NXF_HOME_DIR CONTAINER_CACHE SRA_DIR REFERENCE_DOWNLOADS SAMPLESHEET RESULTS_DIR ENSEMBL_PEP PROFILE ACCOUNT NORMAL_PARTITION BIGMEM_PARTITION; do [[ -n ${!v} ]] || die "missing required value: $v"; done
for v in NORMAL_CPU_THRESHOLD NORMAL_MEMORY_THRESHOLD_GB MAX_CPUS MAX_MEMORY_GB QUEUE_SIZE; do [[ ${!v} =~ ^[1-9][0-9]*$ ]] || die "$v must be positive"; done
[[ $SLURM_LOG_TEMPLATE == *'{job_id}'* ]] || die '--slurm-log-template must contain {job_id}'
for exe in "$NEXTFLOW" "$HOST_PYTHON" "$APPTAINER"; do [[ -x $exe ]] || die "not executable: $exe"; done
JOB_ID=${SLURM_JOB_ID:-manual-$(date +%Y%m%d_%H%M%S)}
RUN_TMP=$TMP_ROOT/pgtk-$JOB_ID; TRACE_FILE=$RESULTS_DIR/pipeline_trace-$JOB_ID.tsv; RESOURCE_PREFIX=$RESULTS_DIR/resource_usage-$JOB_ID
NEXTFLOW_LOG=$PROJECT_DIR/.nextflow.log; SLURM_LOG=${SLURM_LOG_TEMPLATE//\{job_id\}/$JOB_ID}
finalize(){ code=$?; trap - EXIT; set +e; [[ -s $TRACE_FILE ]] && "$HOST_PYTHON" "$PROJECT_DIR/analyze_pipeline_trace.py" --trace "$TRACE_FILE" --output-prefix "$RESOURCE_PREFIX"; "$HOST_PYTHON" "$PROJECT_DIR/collect_pipeline_failures.py" --job-id "$JOB_ID" --exit-code "$code" --trace "$TRACE_FILE" --nextflow-log "$NEXTFLOW_LOG" --slurm-log "$SLURM_LOG" --results-dir "$RESULTS_DIR"; rm -rf "$RUN_TMP"; exit "$code"; }
trap finalize EXIT
for f in "$PROJECT_DIR/main.nf" "$PROJECT_DIR/nextflow.config" "$SAMPLESHEET" "$PYSAM_IMAGE" "$ENSEMBL_PEP"; do [[ -s $f ]] || die "missing file: $f"; done
for d in "$WORK_DIR" "$TMP_ROOT" "$NXF_HOME_DIR" "$CONTAINER_CACHE" "$RESULTS_DIR"; do mkdir -p "$d" || die "cannot create $d"; done
cd "$PROJECT_DIR" || die "cannot enter project directory"
if [[ -n $JAVA_MODULE ]]; then command -v module >/dev/null || die 'module command unavailable'; module purge; module load "$JAVA_MODULE" || die "cannot load $JAVA_MODULE"; fi
export NXF_HOME=$NXF_HOME_DIR NXF_APPTAINER_CACHEDIR=$CONTAINER_CACHE NXF_OPTS=$NXF_OPTS_VALUE
export TMPDIR=$RUN_TMP TMP=$RUN_TMP TEMP=$RUN_TMP APPTAINER_TMPDIR=$RUN_TMP/apptainer
export PGTK_ACCOUNT=$ACCOUNT PGTK_NORMAL_PARTITION=$NORMAL_PARTITION PGTK_BIGMEM_PARTITION=$BIGMEM_PARTITION
export PGTK_NORMAL_CPU_THRESHOLD=$NORMAL_CPU_THRESHOLD PGTK_NORMAL_MEMORY_THRESHOLD_GB=$NORMAL_MEMORY_THRESHOLD_GB PGTK_MAX_CPUS=$MAX_CPUS PGTK_MAX_MEMORY_GB=$MAX_MEMORY_GB PGTK_QUEUE_SIZE=$QUEUE_SIZE PGTK_SUBMIT_RATE_LIMIT=$SUBMIT_RATE_LIMIT
mkdir -p "$TMPDIR" "$APPTAINER_TMPDIR" "$RESULTS_DIR/failure_logs"
"$HOST_PYTHON" "$PROJECT_DIR/validate_runtime_inputs.py" --project-dir "$PROJECT_DIR" --reference-downloads "$REFERENCE_DOWNLOADS" --container-cache "$CONTAINER_CACHE" --sra-dir "$SRA_DIR" --samplesheet "$SAMPLESHEET" --ensembl-pep "$ENSEMBL_PEP" --work-dir "$WORK_DIR" --tmp-root "$TMP_ROOT" --results-dir "$RESULTS_DIR" --host-python "$HOST_PYTHON" --apptainer "$APPTAINER" --pysam-image "$PYSAM_IMAGE" || die 'runtime validation failed'
"$NEXTFLOW" run "$PROJECT_DIR/main.nf" -profile "$PROFILE" -ansi-log false -with-apptainer -work-dir "$WORK_DIR" -with-trace "$TRACE_FILE" -with-report "$RESULTS_DIR/pipeline_report-$JOB_ID.html" -with-timeline "$RESULTS_DIR/pipeline_timeline-$JOB_ID.html" -with-dag "$RESULTS_DIR/pipeline_dag-$JOB_ID.html" -resume --samplesheet "$SAMPLESHEET" --outdir "$RESULTS_DIR" --sra_dir "$SRA_DIR" --reference_downloads "$REFERENCE_DOWNLOADS" --container_cache "$CONTAINER_CACHE" --pysam_image "$PYSAM_IMAGE" --ensembl_pep "$ENSEMBL_PEP" --go_obo "$REFERENCE_DOWNLOADS/go-basic.obo" --go_gaf "$REFERENCE_DOWNLOADS/goa_human.gaf.gz" "${PIPELINE_ARGS[@]}"
