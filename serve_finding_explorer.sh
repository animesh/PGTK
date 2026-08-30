#!/usr/bin/env bash
set -euo pipefail
root=${1:-$(pwd -P)}
port=${2:-8765}
image=${PGTK_IGV_REPORTS_IMAGE:?Set PGTK_IGV_REPORTS_IMAGE to the igv-reports image}
apptainer_bin=${PGTK_APPTAINER:-$(command -v apptainer || true)}
[[ -n $apptainer_bin && -x $apptainer_bin ]] || { echo "ERROR: apptainer is not executable" >&2; exit 2; }
root=$(cd "$root" && pwd -P)
config=$root/explorer_config.json
[[ -s $root/server.py && -s $config ]] || { echo "ERROR: incomplete explorer directory: $root" >&2; exit 2; }
[[ -s $image ]] || { echo "ERROR: IGV reports image not found: $image" >&2; exit 2; }
readarray -t resources < <(python3 - "$config" <<'PY2'
import json,os,sys
from pathlib import Path
root=Path(sys.argv[1]).resolve().parent
config=json.loads(Path(sys.argv[1]).read_text())
genome=Path(os.environ.get('PGTK_IGV_GENOME',config['genome'])).expanduser().resolve()
print(genome)
for track in config['tracks']:
 p=Path(track['path']).expanduser();print(p.resolve() if p.is_absolute() else (root/p).resolve())
PY2
)
((${#resources[@]})) || { echo "ERROR: explorer has no resources" >&2; exit 2; }
for resource in "${resources[@]}"; do [[ -s $resource ]] || { echo "ERROR: missing explorer resource: $resource" >&2; exit 2; }; done
findings_root=$(dirname "$root")
genome_dir=$(dirname "${resources[0]}")
binds=(--bind "$findings_root:$findings_root")
[[ $genome_dir == "$findings_root" || $genome_dir == "$findings_root"/* ]] || binds+=(--bind "$genome_dir:$genome_dir")
exec "$apptainer_bin" exec --cleanenv --no-home --pid "${binds[@]}" --pwd "$root" "$image" python3 "$root/server.py" "$port"
