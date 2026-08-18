#!/usr/bin/env bash
#SBATCH --job-name=pgtk-full-update-check
#SBATCH --time=00:45:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=pgtk-full-update-check-%j.log
#SBATCH --error=pgtk-full-update-check-%j.log
set -euo pipefail
PROJECT_DIR=${PGTK_PROJECT_DIR:-$(pwd -P)}
NEXTFLOW=${PGTK_NEXTFLOW:-$(command -v nextflow || true)}
APPTAINER=${PGTK_APPTAINER:-$(command -v apptainer || true)}
IMAGE=${PGTK_PYSAM_IMAGE:-$PROJECT_DIR/singularity_cache/quay.io-biocontainers-pysam-0.24.0--py312hf5ad864_1.img}
EVENTS=$PROJECT_DIR/validation_fixtures/job_19181644/pgtk_igv.events.tsv
[[ -x $NEXTFLOW && -x $APPTAINER && -s $IMAGE && -s $EVENTS ]]
cd "$PROJECT_DIR"
bash validate_pipeline_commands.sh --project-dir "$PROJECT_DIR" --nextflow "$NEXTFLOW"
! grep -Fq 'DDX1' build_finding_igv_reviews.py
grep -Fq "params.finding_priority_mode = 'all'" main.nf
grep -Fq 'params.igv_report_limit = 0' main.nf
grep -Fq 'params.igv_report_timeout_seconds = 600' main.nf
rm -rf full_update_validation;mkdir full_update_validation
"$APPTAINER" exec --cleanenv -B "$PROJECT_DIR:/project" "$IMAGE" sh -c '
set -euo pipefail
cd /project
python3 - <<"PY"
import pysam
assert pysam.__version__=="0.24.0" and pysam.__samtools_version__=="1.23.1"
print("PYSAM_RUNTIME_OK")
PY
python3 test_igv_evidence_bundle.py
python3 test_finding_igv_reviews.py
python3 build_finding_igv_reviews.py --events validation_fixtures/job_19181644/pgtk_igv.events.tsv --genome hg38 --output-dir full_update_validation --plan-only
'
grep -qx 'Input variant rows: 17343' full_update_validation/consolidation_summary.txt
grep -qx 'Consolidated findings: 15321' full_update_validation/consolidation_summary.txt
grep -qx 'Merged duplicate rows: 2022' full_update_validation/consolidation_summary.txt
[[ $(tail -n +2 full_update_validation/event_consolidation.tsv | cut -f1 | sort | uniq -d | wc -l) -eq 0 ]]
echo GROUND_TRUTH_19181644_OK
echo PHASE2_GENERIC_CLI_PARAMETERS_OK
echo 'RESULT: PASSED'
