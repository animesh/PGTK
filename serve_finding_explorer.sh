#!/usr/bin/env bash
set -euo pipefail

EXPLORER_DIR=${1:?Usage: serve_finding_explorer.sh EXPLORER_DIR [PORT]}
PORT=${2:-8765}
APPTAINER=${PGTK_APPTAINER:-$(command -v apptainer || true)}
PYSAM_IMAGE=${PGTK_PYSAM_IMAGE:-}
IGV_IMAGE=${PGTK_IGV_REPORTS_IMAGE:-}

EXPLORER_DIR=$(readlink -f "$EXPLORER_DIR")
test -d "$EXPLORER_DIR" || { echo "ERROR: missing explorer directory: $EXPLORER_DIR" >&2; exit 1; }
test -x "$APPTAINER" || { echo "ERROR: Apptainer executable not found: $APPTAINER" >&2; exit 1; }
test -s "$PYSAM_IMAGE" || { echo "ERROR: set PGTK_PYSAM_IMAGE to the pinned Pysam image" >&2; exit 1; }
test -s "$IGV_IMAGE" || { echo "ERROR: set PGTK_IGV_REPORTS_IMAGE to the pinned igv-reports image" >&2; exit 1; }
test -s "$EXPLORER_DIR/server.py" || { echo "ERROR: missing server.py" >&2; exit 1; }
test -s "$EXPLORER_DIR/prepare_event_igv_tracks.py" || { echo "ERROR: missing prepare_event_igv_tracks.py" >&2; exit 1; }
test -s "$EXPLORER_DIR/explorer_config.json" || { echo "ERROR: missing explorer_config.json" >&2; exit 1; }

"$APPTAINER" exec --cleanenv --no-home "$PYSAM_IMAGE" python3 -c 'import pysam; assert pysam.__version__ == "0.24.0"'
"$APPTAINER" exec --cleanenv --no-home "$IGV_IMAGE" sh -c 'command -v create_report >/dev/null'

export PGTK_APPTAINER="$APPTAINER"
export PGTK_PYSAM_IMAGE="$PYSAM_IMAGE"
export PGTK_IGV_REPORTS_IMAGE="$IGV_IMAGE"
exec python3 "$EXPLORER_DIR/server.py" "$PORT"
