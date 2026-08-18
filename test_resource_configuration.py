#!/usr/bin/env python3
import re,subprocess,tempfile
from pathlib import Path
main=Path('main.nf').read_text();config=Path('nextflow.config').read_text();review=Path('build_finding_igv_reviews.py').read_text();explorer=Path('build_finding_explorer.py').read_text()
processes=re.findall(r'^process\s+(\w+)\s*\{',main,re.M);assert len(processes)==72 and len(set(processes))==72
assert 'process BUILD_FINDING_EXPLORER {' in main and 'GENERATE_PRIORITY_IGV_REPORTS' not in main
assert '--priority-limit 0' in main and 'params.igv_report_limit' not in main and 'timeout 600' not in main
assert 'withName: BUILD_FINDING_EXPLORER {' in config
assert 'fetch(until_eof=True)' in review and 'alignment_store.sqlite' not in review and 'sqlite3' not in review
assert 'import sqlite' not in explorer and 'jsonl.gz' in explorer
block=main.split('process BUILD_FINDING_EXPLORER {',1)[1].split('process BUILD_COMPARATIVE_ADVANTAGE_REPORT {',1)[0]
body=block.split('"""',2)[1]
subs={'${explorer_script}':'build_finding_explorer.py','${finding_reviews}':'finding_reviews','${genome}':'genome.fa','${params.read_validation_padding}':'150','${server_launcher}':'serve_finding_explorer.sh'}
for a,b in subs.items():body=body.replace(a,b)
body=body.replace('\\$','$')
with tempfile.NamedTemporaryFile('w',suffix='.sh') as h:
 h.write(body);h.flush();subprocess.run(['bash','-n',h.name],check=True)
assert '""' not in body.splitlines()[-2]
print('PASS: 72 processes, exact rendered explorer shell syntax, database-free explorer, zero finding limit')
