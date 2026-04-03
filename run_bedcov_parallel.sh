#!/usr/bin/env bash
# run_bedcov_parallel.sh
# Runs samtools bedcov on all 28 BAMs in parallel (max 6 jobs at a time)
# Output: results/coverage/bedcov/<sample>.bedcov.tsv
# Logs:   results/coverage/bedcov/logs/<sample>.log
# On completion writes: results/coverage/bedcov/DONE

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BAM_DIR="$SCRIPT_DIR"
BED="$SCRIPT_DIR/results/coverage/windows_1mb.bed"
OUT_DIR="$SCRIPT_DIR/results/coverage/bedcov"
LOG_DIR="$OUT_DIR/logs"
MAX_PARALLEL=6

mkdir -p "$OUT_DIR" "$LOG_DIR"

BAMS=(
  TK1049 TK1050 TK1051
  TK12R2 TK12R3
  TK131 TK132 TK133
  TK141 TK142 TK143
  TK16R1 TK16R2 TK16R3
  TK18R1 TK18R2 TK18R3
  TK91L003 TK91L005 TK91L006
  TK92L002 TK92L003 TK92L005 TK92L006
  TK93L002 TK93L003 TK93L005 TK93L006
)

echo "[$(date '+%H:%M:%S')] Starting bedcov on ${#BAMS[@]} BAMs (max $MAX_PARALLEL parallel jobs)"
echo "[$(date '+%H:%M:%S')] Output dir: $OUT_DIR"
echo ""

run_one() {
  local PREFIX="$1"
  local BAM="$BAM_DIR/${PREFIX}.markdup.sorted.bam"
  local OUT="$OUT_DIR/${PREFIX}.bedcov.tsv"
  local LOG="$LOG_DIR/${PREFIX}.log"

  if [[ -f "$OUT" && -s "$OUT" ]]; then
    echo "  [SKIP] $PREFIX — already done ($(wc -l < "$OUT") lines)"
    return 0
  fi

  echo "  [START] $PREFIX"
  local T0=$SECONDS
  samtools bedcov "$BED" "$BAM" > "$OUT" 2>"$LOG"
  local RC=$?
  local ELAPSED=$(( SECONDS - T0 ))

  if [[ $RC -eq 0 ]]; then
    echo "  [OK]   $PREFIX — done in ${ELAPSED}s ($(wc -l < "$OUT") windows)"
  else
    echo "  [FAIL] $PREFIX — exit $RC after ${ELAPSED}s (see $LOG)"
    rm -f "$OUT"
  fi
}

export -f run_one
export BAM_DIR OUT_DIR LOG_DIR BED

# Use GNU parallel if available, otherwise xargs -P
if command -v parallel &>/dev/null; then
  echo "Using GNU parallel"
  printf '%s\n' "${BAMS[@]}" | parallel -j "$MAX_PARALLEL" run_one {}
else
  echo "Using xargs -P$MAX_PARALLEL"
  printf '%s\n' "${BAMS[@]}" | xargs -P "$MAX_PARALLEL" -I{} bash -c 'run_one "$@"' _ {}
fi

# Summary
N_OK=$(ls "$OUT_DIR"/*.bedcov.tsv 2>/dev/null | wc -l)
N_TOTAL=${#BAMS[@]}
echo ""
echo "[$(date '+%H:%M:%S')] Finished: $N_OK / $N_TOTAL BAMs completed"

if [[ "$N_OK" -eq "$N_TOTAL" ]]; then
  touch "$OUT_DIR/DONE"
  echo "[$(date '+%H:%M:%S')] All done. Wrote: $OUT_DIR/DONE"
  echo "[$(date '+%H:%M:%S')] Now run: python3 $SCRIPT_DIR/plot_windowed_cn.py"
else
  echo "[$(date '+%H:%M:%S')] WARNING: $(( N_TOTAL - N_OK )) BAMs failed. Check logs in $LOG_DIR"
  exit 1
fi
