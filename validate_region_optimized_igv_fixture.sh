#!/usr/bin/env bash
set -euo pipefail
project=${1:-$(pwd -P)}
fixture=${2:-$project/validation_fixtures/region_optimized_igv}
output=${3:-$project/region_optimized_igv_validation}
rm -rf "$output"
python3 "$project/build_finding_igv_reviews.py" \
  --events "$fixture/fixture.events.tsv" \
  --bam "TK12=$fixture/pgtk_igv.TK12.fixture.bam" \
  --bam "TK13=$fixture/pgtk_igv.TK13.fixture.bam" \
  --bam "TK14=$fixture/pgtk_igv.TK14.fixture.bam" \
  --genome genome.fa \
  --output-dir "$output" \
  --padding 150 \
  --mapq 20 \
  --baseq 20 \
  --reference-display-reads 20 \
  --alt-display-reads 100 \
  --finding-classes rna_variant,progression_variant,fusion,splice_junction \
  --priority-mode all \
  --priority-limit 0
python3 - "$output" <<'PY'
from pathlib import Path
import csv,sys
out=Path(sys.argv[1])
required=['findings_manifest.tsv','bam_manifest.tsv','support_labels.bed','priority_findings.bed','review.igv.batch.txt','igv.session.xml','consolidation_summary.txt']
for name in required:
    assert (out/name).is_file() and (out/name).stat().st_size>0,name
with (out/'findings_manifest.tsv').open() as handle:
    findings=list(csv.DictReader(handle,delimiter='\t'))
assert len(findings)>1000,len(findings)
assert any(row['REF'] and row['ALT'] for row in findings)
assert any(not row['REF'] and not row['ALT'] for row in findings)
assert not (out/'alignment_store.sqlite').exists()
assert not any(path.is_dir() for path in out.iterdir())
print(f'PASS: region-optimized fixture findings={len(findings)} files={sum(1 for _ in out.iterdir())}')
PY
