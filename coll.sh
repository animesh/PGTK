#!/usr/bin/env bash
set -euo pipefail

project=/cluster/projects/nn9036k/scrbkup/pgtk/checkport
review=/cluster/work/users/ash022/work/95/3fd43550b0bdbb46dda914b2885309/finding_reviews
igv_image=/cluster/projects/nn9036k/scrbkup/pgtk/singularity_cache/quay.io-biocontainers-igv-reports-1.16.0--pyh7e72e81_0.img
output="$HOME/pgtk_explorer_redesign_inputs_small_19244152.tgz"
staging=$(mktemp -d)
trap 'rm -rf "$staging"' EXIT
mkdir -p "$staging/project" "$staging/production_review" "$staging/runtime"

cd "$project"
find . -maxdepth 1 -type f \
  ! -name '*.tgz' \
  ! -name '*.tar.gz' \
  ! -name '*.html' \
  ! -name '*.log' \
  ! -name '*.png' \
  ! -name '*.pdf' \
  ! -name '*.pyc' \
  ! -name 'main.nf.before_*' \
  ! -name 'SHA256SUMS*' \
  -print0 | tar --null -T - -cf - | tar -C "$staging/project" -xf -

for name in findings_manifest.tsv bam_manifest.tsv consolidation_summary.txt README.txt; do
  test -s "$review/$name"
  gzip -c -9 "$review/$name" > "$staging/production_review/$name.gz"
done

python3 - "$review/findings_manifest.tsv" "$staging/production_review/findings_stratified.tsv" <<'PY'
import csv
import sys
from collections import defaultdict

source, target = sys.argv[1:]
limit_per_group = 25
counts = defaultdict(int)
selected = []
with open(source, encoding='utf-8', newline='') as handle:
    reader = csv.DictReader(handle, delimiter='\t')
    fields = reader.fieldnames
    for row in reader:
        key = (row.get('Sample',''), row.get('EvidenceClasses',''), row.get('Chrom',''))
        if counts[key] < limit_per_group:
            selected.append(row)
            counts[key] += 1
with open(target, 'w', encoding='utf-8', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    writer.writerows(selected)
print(f'Stratified rows: {len(selected)}', file=sys.stderr)
PY

gzip -9 "$staging/production_review/findings_stratified.tsv"

{
  printf 'Path\tBytes\n'
  find -L "$review" -maxdepth 1 -type f -exec sh -c 'for f do printf "%s\t%s\n" "${f##*/}" "$(wc -c < "$f")"; done' sh {} +
} | sort > "$staging/production_review/inventory.tsv"

apptainer exec --no-home --pid -B /cluster "$igv_image" sh -c '
  set -eu
  echo "CREATE_REPORT=$(command -v create_report || command -v create_reports)"
  (create_report --help || create_reports --help) 2>&1 || true
  python3 - <<"PY"
import importlib.metadata as md
for dist in md.distributions():
    name = dist.metadata.get("Name", "")
    if "igv" in name.lower():
        print("DISTRIBUTION", name, dist.version, sep="\t")
PY
' > "$staging/runtime/igv_reports_runtime.txt"

find "$staging" -type f -print0 | sort -z | xargs -0 sha256sum > "$staging/SHA256SUMS"
tar -C "$staging" -czf "$output" .
gzip -t "$output"
sha256sum "$output" > "$output.sha256"
ls -lh "$output" "$output.sha256"
