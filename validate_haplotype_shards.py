#!/usr/bin/env python3
import argparse, csv, gzip, re, sys
from pathlib import Path

def open_text(path):
    return gzip.open(path, 'rt', encoding='utf-8', errors='replace')

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--sample',required=True); p.add_argument('--expected',required=True,type=int)
    p.add_argument('--gvcf',nargs='+',required=True); p.add_argument('--tbi',nargs='+',required=True)
    p.add_argument('--output-prefix',required=True); a=p.parse_args()
    gvcfs=sorted(map(Path,a.gvcf)); tbis=sorted(map(Path,a.tbi)); failures=[]; rows=[]
    if len(gvcfs)!=a.expected: failures.append(f'expected {a.expected} GVCFs, observed {len(gvcfs)}')
    if len(tbis)!=a.expected: failures.append(f'expected {a.expected} indexes, observed {len(tbis)}')
    index_names={x.name.removesuffix('.tbi') for x in tbis}
    seen_samples=set()
    for path in gvcfs:
        status=[]; records=0; sample_names=[]
        if not path.is_file() or path.stat().st_size==0: status.append('MISSING_OR_EMPTY')
        else:
            try:
                with open_text(path) as h:
                    for line in h:
                        if line.startswith('#CHROM'):
                            sample_names=line.rstrip().split('\t')[9:]
                        elif not line.startswith('#'): records+=1
            except Exception as e: status.append(f'UNREADABLE:{e}')
        if path.name not in index_names: status.append('INDEX_MISSING')
        seen_samples.update(sample_names)
        rows.append({'Sample':a.sample,'Shard':path.name,'Records':records,'VCF samples':';'.join(sample_names),'Status':';'.join(status) or 'PASS'})
        failures.extend(f'{path.name}: {x}' for x in status)
    if seen_samples and seen_samples!={a.sample}: failures.append(f'VCF sample mismatch: {sorted(seen_samples)}')
    out=Path(a.output_prefix)
    with Path(f'{out}.tsv').open('w',newline='',encoding='utf-8') as h:
        w=csv.DictWriter(h,fieldnames=list(rows[0]) if rows else ['Sample','Shard','Records','VCF samples','Status'],delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(rows)
    Path(f'{out}.summary.txt').write_text(f'Sample: {a.sample}\nExpected shards: {a.expected}\nObserved shards: {len(gvcfs)}\nStatus: {"FAIL" if failures else "PASS"}\n'+''.join(f'Failure: {x}\n' for x in failures),encoding='utf-8')
    if failures: raise SystemExit('; '.join(failures))
if __name__=='__main__': main()
