#!/usr/bin/env python3
import argparse, csv, gzip, hashlib
from collections import Counter
from pathlib import Path

def open_text(path): return gzip.open(path,'rt',encoding='utf-8',errors='replace') if str(path).endswith('.gz') else open(path,encoding='utf-8',errors='replace')
def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
    return h.hexdigest()
def summarize(sample,stage,path):
    n=snps=indels=multi=0; filters=Counter(); titv=Counter(); contigs=set()
    pur={('A','G'),('G','A'),('C','T'),('T','C')}
    with open_text(path) as h:
        for line in h:
            if line.startswith('#'): continue
            f=line.rstrip().split('\t');
            if len(f)<7: continue
            chrom,ref,alts,flt=f[0],f[3],f[4].split(','),f[6]; contigs.add(chrom); filters[flt]+=1
            for alt in alts:
                n+=1
                if len(alts)>1: multi+=1
                if len(ref)==len(alt)==1:
                    snps+=1; titv['Ti' if (ref.upper(),alt.upper()) in pur else 'Tv']+=1
                else: indels+=1
    return {'Sample':sample,'Stage':stage,'File':Path(path).name,'Alleles':n,'SNPs':snps,'Indels':indels,'Multiallelic alleles':multi,'Ti':titv['Ti'],'Tv':titv['Tv'],'Ti/Tv':f"{titv['Ti']/titv['Tv']:.6f}" if titv['Tv'] else 'NA','Contigs':len(contigs),'Filter counts':';'.join(f'{k}:{v}' for k,v in sorted(filters.items())),'SHA256':sha256(path)}
def main():
    p=argparse.ArgumentParser(); p.add_argument('--sample',required=True); p.add_argument('--raw',required=True); p.add_argument('--pass-vcf',required=True); p.add_argument('--rna',required=True); p.add_argument('--genome',required=True); p.add_argument('--calling-confidence',required=True); p.add_argument('--soft-clipped-setting',required=True); p.add_argument('--pcr-indel-model',default='GATK_DEFAULT'); p.add_argument('--output-prefix',required=True); a=p.parse_args()
    rows=[summarize(a.sample,'raw_genotyped',a.raw),summarize(a.sample,'hard_filter_pass',a.pass_vcf),summarize(a.sample,'rna_validated',a.rna)]
    out=Path(a.output_prefix)
    with Path(f'{out}.tsv').open('w',newline='',encoding='utf-8') as h:
        w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(rows)
    Path(f'{out}.report.md').write_text('# Variant-stage QC\n\n'+'\n'.join(f"- {r['Stage']}: {r['Alleles']} alleles, {r['SNPs']} SNPs, {r['Indels']} indels, Ti/Tv {r['Ti/Tv']}" for r in rows)+'\n',encoding='utf-8')
    command_headers=[]
    with open_text(a.raw) as h:
        for line in h:
            if line.startswith('##GATKCommandLine='): command_headers.append(line.strip())
            elif line.startswith('#CHROM'): break
    provenance=[
        ('sample',a.sample),('reference_file',Path(a.genome).name),('reference_sha256',sha256(a.genome)),
        ('raw_vcf_sha256',sha256(a.raw)),('pass_vcf_sha256',sha256(a.pass_vcf)),('rna_vcf_sha256',sha256(a.rna)),
        ('hc_calling_confidence',a.calling_confidence),('hc_dont_use_soft_clipped_bases',a.soft_clipped_setting),
        ('hc_pcr_indel_model',a.pcr_indel_model),('gatk_command_headers',' || '.join(command_headers))]
    with Path(f'{out}.provenance.tsv').open('w',newline='',encoding='utf-8') as h:
        w=csv.writer(h,delimiter='\t',lineterminator='\n'); w.writerow(['Key','Value']); w.writerows(provenance)
if __name__=='__main__': main()
