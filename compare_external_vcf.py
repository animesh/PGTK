#!/usr/bin/env python3
import argparse,csv,gzip,re
from pathlib import Path

def open_text(p): return gzip.open(p,'rt',encoding='utf-8',errors='replace') if str(p).endswith('.gz') else open(p,encoding='utf-8',errors='replace')
def norm_chrom(c): return c[3:] if c.lower().startswith('chr') else c
def trim(pos,ref,alt):
    while len(ref)>1 and len(alt)>1 and ref[-1]==alt[-1]: ref,alt=ref[:-1],alt[:-1]
    while len(ref)>1 and len(alt)>1 and ref[0]==alt[0]: ref,alt,pos=ref[1:],alt[1:],pos+1
    return pos,ref,alt
def read_vcf(path):
    out=set()
    with open_text(path) as h:
        for line in h:
            if line.startswith('#'): continue
            f=line.rstrip().split('\t');
            if len(f)<5: continue
            for alt in f[4].split(','):
                pos,ref,alt=trim(int(f[1]),f[3].upper(),alt.upper()); out.add((norm_chrom(f[0]),pos,ref,alt))
    return out
def kind(v): return 'SNP' if len(v[2])==len(v[3])==1 else 'INDEL'
def main():
    p=argparse.ArgumentParser(); p.add_argument('--sample',required=True); p.add_argument('--pgtk',required=True); p.add_argument('--external',required=True); p.add_argument('--output-prefix',required=True); a=p.parse_args()
    x=read_vcf(a.pgtk); y=read_vcf(a.external); shared=x&y; rows=[]
    for label,subset in [('PGTK',x),('External',y),('Shared',shared),('PGTK-only',x-y),('External-only',y-x)]: rows.append({'Sample':a.sample,'Set':label,'All':len(subset),'SNPs':sum(kind(v)=='SNP' for v in subset),'Indels':sum(kind(v)=='INDEL' for v in subset)})
    j=len(shared)/len(x|y) if x|y else 1.0; rows.append({'Sample':a.sample,'Set':'Metrics','All':f'Jaccard={j:.6f};PGTK_recall={len(shared)/len(x) if x else 1:.6f};External_recall={len(shared)/len(y) if y else 1:.6f}','SNPs':'','Indels':''})
    out=Path(a.output_prefix)
    with Path(f'{out}.summary.tsv').open('w',newline='',encoding='utf-8') as h: w=csv.DictWriter(h,fieldnames=['Sample','Set','All','SNPs','Indels'],delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(rows)
    with Path(f'{out}.shared.tsv').open('w',newline='',encoding='utf-8') as h:
        w=csv.writer(h,delimiter='\t',lineterminator='\n'); w.writerow(['CHROM','POS','REF','ALT']); w.writerows(sorted(shared))
    Path(f'{out}.report.md').write_text(f'# External VCF comparison\n\n- Sample: {a.sample}\n- PGTK alleles: {len(x)}\n- External alleles: {len(y)}\n- Shared: {len(shared)}\n- Jaccard: {j:.6f}\n\nComparison splits multiallelic records, harmonizes chr prefixes and trims common allele sequence. Full left alignment still requires the same reference FASTA.\n',encoding='utf-8')
if __name__=='__main__': main()
