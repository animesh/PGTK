#!/usr/bin/env python3
from pathlib import Path
import ast
source=Path(__file__).with_name('validate_proteogenomic_reads.py').read_text()
ast.parse(source)
for text in ['import pysam','pysam.idxstats','pysam.view(str(bam), region)','catch_stdout=False','pysam.index(str(ob))','__samtools_version__']:
    assert text in source, text
assert '["samtools"' not in source
print('PASS: read validation uses the Pysam Python API without a samtools executable')
