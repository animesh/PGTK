#!/usr/bin/env bash
set -euo pipefail
project=${1:-$(pwd -P)}
review=${2:-$project/single_pass_igv_validation}
output=${3:-$project/finding_explorer_validation}
command -v create_report >/dev/null
test -s "$review/findings_manifest.tsv"
test -s "$review/bam_manifest.tsv"
rm -rf "$output"
python3 "$project/build_finding_explorer.py" --manifest "$review/findings_manifest.tsv" --bam-manifest "$review/bam_manifest.tsv" --genome genome.fa --output-dir "$output"
python3 - "$review/findings_manifest.tsv" "$output" <<'PY'
import csv,sqlite3,sys
from pathlib import Path
manifest=Path(sys.argv[1]);out=Path(sys.argv[2])
with manifest.open() as h: expected=sum(1 for _ in h)-1
con=sqlite3.connect(out/'findings.sqlite');indexed=con.execute('select count(*) from findings').fetchone()[0];mapped=con.execute('select count(*) from findings where partition_id is not null').fetchone()[0]
assert expected==indexed==mapped
assert 'Findings discarded: 0' in (out/'coverage_summary.txt').read_text()
for name in ('index.html','server.py','partition_manifest.tsv','finding_to_partition.tsv','explorer_config.json'):assert (out/name).stat().st_size>0
print(f'PASS: complete explorer findings={indexed} discarded=0')
PY
