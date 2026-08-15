#!/usr/bin/env python3
import re
from pathlib import Path
main=Path('main.nf').read_text();config=Path('nextflow.config').read_text();slurm=Path('scratch.slurm').read_text();assets=Path('download_assets.sh').read_text()
processes=re.findall(r'^process\s+(\w+)\s*\{',main,re.M)
assert len(processes)==72 and len(set(processes))==72
assert 'GENERATE_PRIORITY_IGV_REPORTS' in processes and 'RENDER_PRIORITY_IGV_SNAPSHOTS' not in processes
assert 'params.generate_priority_igv_reports = true' in main
assert 'params.igv_report_limit = 0' in main
assert "params.finding_priority_mode = 'all'" in main
assert "params.finding_priority_genes = ''" in main
block=main.split('process GENERATE_PRIORITY_IGV_REPORTS {',1)[1].split('process BUILD_COMPARATIVE_ADVANTAGE_REPORT {',1)[0]
for required in ('igv-reports-1.16.0--pyh7e72e81_0.img','command -v create_report || command -v create_reports','--fasta ${genome}','report_manifest.tsv','igv_report_timeout_seconds'):
    assert required in block,required
assert 'xvfb' not in block.lower() and 'java -jar' not in block.lower()
assert 'withName: GENERATE_PRIORITY_IGV_REPORTS {' in config
assert 'withName: RENDER_PRIORITY_IGV_SNAPSHOTS {' not in config
assert '--pysam_image "$PYSAM_IMAGE"' in slurm
assert 'igv-reports:1.16.0--pyh7e72e81_0' in assets and 'igv-xvfb' not in assets
print('PASS: 72-process configuration uses pinned self-contained IGV Reports with no Desktop/Xvfb dependency')
