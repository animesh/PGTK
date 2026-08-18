#!/usr/bin/env bash
set -euo pipefail
root=${1:-$(pwd -P)}; port=${2:-8765}; image=${PGTK_IGV_REPORTS_IMAGE:?Set PGTK_IGV_REPORTS_IMAGE}
exec apptainer exec --no-home --pid -B /cluster -B "$root:$root" --pwd "$root" "$image" python3 server.py "$port"
