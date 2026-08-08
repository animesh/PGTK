#!/usr/bin/env python3
import os
import subprocess
import sys
import tempfile
from pathlib import Path

root = Path(__file__).resolve().parent
with tempfile.TemporaryDirectory(prefix='pgtk_igv_fixture_') as temp:
    temp = Path(temp)
    bindir = temp / 'bin'
    bindir.mkdir()
    samtools = bindir / 'samtools'
    samtools.write_text("""#!/bin/sh
if [ \"$1\" = view ]; then
  out=''; prev=''
  for value in \"$@\"; do
    if [ \"$prev\" = -o ]; then out=$value; fi
    prev=$value
  done
  : > \"$out\"
elif [ \"$1\" = index ]; then
  : > \"$2.bai\"
fi
""", encoding='utf-8')
    samtools.chmod(0o755)
    (temp / 'genome.fa').write_text('>1\nA\n>2\nA\n>7\nA\n', encoding='utf-8')
    (temp / 'TK1.rna.validated.vcf').write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t10\t.\tA\tG\t.\tPASS\t.\n', encoding='utf-8')
    (temp / 'TK1.fusion.validated.tsv').write_text('#gene1\tgene2\tbreakpoint1\tbreakpoint2\nKIF1A\tDENND11\t2:20\t7:30\n', encoding='utf-8')
    (temp / 'TK1.splice.audit.tsv').write_text('Event\tStatus\tJunctions\nTX1\tRNA_VALIDATED\t1:40-50\n', encoding='utf-8')
    (temp / 'TK1.Aligned.sortedByCoord.out.bam').touch()
    prefix = temp / 'pgtk_igv'
    env = dict(os.environ)
    env['PATH'] = str(bindir) + os.pathsep + env.get('PATH', '')
    subprocess.run([
        sys.executable, str(root / 'build_igv_evidence_bundle.py'),
        '--genome', str(temp / 'genome.fa'),
        '--rna-vcf', str(temp / 'TK1.rna.validated.vcf'),
        '--fusion-table', str(temp / 'TK1.fusion.validated.tsv'),
        '--splice-table', str(temp / 'TK1.splice.audit.tsv'),
        '--bam', f'TK1={temp / "TK1.Aligned.sortedByCoord.out.bam"}',
        '--output-prefix', str(prefix),
    ], check=True, env=env)
    required = [
        prefix.with_suffix('.events.tsv'), prefix.with_suffix('.events.bed'), prefix.with_suffix('.events.bedpe'),
        prefix.with_suffix('.sample_manifest.tsv'), prefix.with_suffix('.igv.batch.txt'), prefix.with_suffix('.igv.session.xml'),
        prefix.with_suffix('.summary.txt'), Path(str(prefix) + '.TK1.events.bam'), Path(str(prefix) + '.TK1.events.bam.bai'),
    ]
    for path in required:
        if not path.exists():
            raise AssertionError(f'missing IGV fixture output: {path}')
    summary = prefix.with_suffix('.summary.txt').read_text(encoding='utf-8')
    for expected in ('RNA variants: 1', 'Fusions: 1', 'Splice junctions: 1'):
        if expected not in summary:
            raise AssertionError(f'missing summary value: {expected}')
    bedpe = prefix.with_suffix('.events.bedpe').read_text(encoding='utf-8')
    if not bedpe.startswith('2\t19\t20\t7\t29\t30\t'):
        raise AssertionError('interchromosomal BEDPE coordinates are incorrect')
print('PASS: IGV fixture covers RNA variants, splice junctions and both fusion breakpoints')
