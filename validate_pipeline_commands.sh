#!/usr/bin/env bash
set -uo pipefail
PROJECT_DIR= NEXTFLOW= PYTHON=${PGTK_PYTHON:-$(command -v python3 || true)}
while [[ $# -gt 0 ]]; do
    case "$1" in
        --project-dir) PROJECT_DIR=${2:?}; shift 2 ;;
        --nextflow) NEXTFLOW=${2:?}; shift 2 ;;
        --python) PYTHON=${2:?}; shift 2 ;;
        --help|-h) printf 'Usage: %s --project-dir PATH --nextflow PATH [--python PATH]\n' "$0"; exit 0 ;;
        *) printf 'ERROR: unknown option: %s\n' "$1" >&2; exit 2 ;;
    esac
done
[[ -n $PROJECT_DIR && -n $NEXTFLOW && -n $PYTHON ]] || { echo 'ERROR: required arguments missing' >&2; exit 2; }
PROJECT_DIR=$(cd "$PROJECT_DIR" && pwd -P) || exit 2
REPORT="$PROJECT_DIR/pipeline_command_validation_$(date +%Y%m%d_%H%M%S).txt"
PASS=0 FAIL=0
exec > >(tee "$REPORT") 2>&1
pass(){ PASS=$((PASS+1)); printf 'PASS  %s\n' "$*"; }
fail(){ FAIL=$((FAIL+1)); printf 'FAIL  %s\n' "$*"; }

printf 'PGTK SOURCE PREFLIGHT\nProject: %s\n\n' "$PROJECT_DIR"
MANIFEST="$PROJECT_DIR/pipeline_required_files.txt"
if [[ -s $MANIFEST ]]; then
    pass 'Required-file manifest exists'
    while IFS= read -r file; do
        [[ -z $file ]] && continue
        [[ -s "$PROJECT_DIR/$file" ]] && pass "Required file: $file" || fail "Missing or empty: $file"
    done < "$MANIFEST"
else
    fail 'Required-file manifest missing'
fi
[[ -x $NEXTFLOW ]] && pass 'Nextflow executable' || fail "Nextflow executable: $NEXTFLOW"
[[ -x $PYTHON ]] && pass 'Python executable' || fail "Python executable: $PYTHON"
for file in "$PROJECT_DIR"/*.py; do "$PYTHON" -m py_compile "$file" && pass "Python syntax: $(basename "$file")" || fail "Python syntax: $(basename "$file")"; done
for file in "$PROJECT_DIR"/*.sh "$PROJECT_DIR"/*.slurm; do [[ -e $file ]] || continue; bash -n "$file" && pass "Shell syntax: $(basename "$file")" || fail "Shell syntax: $(basename "$file")"; done
mapfile -t processes < <(awk '/^process[[:space:]]+[A-Za-z0-9_]+[[:space:]]*\{/ {print $2}' "$PROJECT_DIR/main.nf")
[[ ${#processes[@]} -eq 73 ]] && pass '73 Nextflow processes' || fail "Expected 73 processes; found ${#processes[@]}"
[[ $(printf '%s\n' "${processes[@]}" | sort | uniq -d | wc -l) -eq 0 ]] && pass 'No duplicate process names' || fail 'Duplicate process names'
for token in human_reviewed_isoforms reference_proteome_archive proteome_archive; do
    if grep -Rq --exclude='pipeline_command_validation_*.txt' "$token" "$PROJECT_DIR/main.nf" "$PROJECT_DIR/download_assets.sh" "$PROJECT_DIR/validate_runtime_inputs.py"; then fail "Unused UniProt dependency remains: $token"; else pass "Unused dependency absent: $token"; fi
done
for token in 'import pysam' 'pysam.idxstats(' 'pysam.view(' 'pysam.index(' '__samtools_version__'; do grep -Fq "$token" "$PROJECT_DIR/validate_proteogenomic_reads.py" && pass "Pysam contract: $token" || fail "Missing Pysam contract: $token"; done
grep -Fq 'run(["samtools"' "$PROJECT_DIR/validate_proteogenomic_reads.py" && fail 'External samtools dependency remains' || pass 'No external samtools dependency'
grep -Fq 'validate_pipeline_commands.sh' "$PROJECT_DIR/scratch.slurm" && pass 'Source preflight wired before launch' || fail 'Source preflight not wired'
grep -Fq 'validate_runtime_inputs.py' "$PROJECT_DIR/scratch.slurm" && pass 'Runtime preflight wired before launch' || fail 'Runtime preflight not wired'
grep -Fq 'command -v apptainer' "$PROJECT_DIR/download_assets.sh" && pass 'Asset downloader uses Apptainer' || fail 'Asset downloader does not use Apptainer'
grep -Fq 'command -v apptainer' "$PROJECT_DIR/download_sra.sh" && pass 'SRA downloader uses Apptainer' || fail 'SRA downloader does not use Apptainer'
grep -Fq 'vdb-validate' "$PROJECT_DIR/download_sra.sh" && pass 'SRA validation wired' || fail 'SRA validation missing'
if [[ -x $NEXTFLOW ]]; then
    tmp=$(mktemp)
    if PGTK_ACCOUNT=validation PGTK_NORMAL_PARTITION=normal PGTK_BIGMEM_PARTITION=bigmem PGTK_NORMAL_CPU_THRESHOLD=20 PGTK_NORMAL_MEMORY_THRESHOLD_GB=160 PGTK_MAX_CPUS=32 PGTK_MAX_MEMORY_GB=512 PGTK_QUEUE_SIZE=200 PGTK_SUBMIT_RATE_LIMIT=60/1min "$NEXTFLOW" inspect "$PROJECT_DIR/main.nf" >"$tmp" 2>&1; then pass 'Nextflow inspect'; else cat "$tmp"; fail 'Nextflow inspect'; fi
    rm -f "$tmp"
fi
rm -rf "$PROJECT_DIR/__pycache__"
printf '\nPASS: %d\nFAIL: %d\nREPORT: %s\n' "$PASS" "$FAIL" "$REPORT"
(( FAIL == 0 )) || exit 1
printf 'RESULT: PASSED\n'
