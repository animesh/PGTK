#!/usr/bin/env bash
set -euo pipefail
PROJECT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
NEXTFLOW=${PGTK_NEXTFLOW:-$(command -v nextflow || true)}
[[ -n $NEXTFLOW && -x $NEXTFLOW ]] || { echo 'ERROR: supply PGTK_NEXTFLOW or put nextflow on PATH' >&2; exit 2; }
exec "$NEXTFLOW" run "$PROJECT_DIR/main.nf" "$@"
