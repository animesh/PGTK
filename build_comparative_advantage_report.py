#!/usr/bin/env python3
import argparse, csv, gzip, hashlib, re
from collections import Counter, defaultdict
from pathlib import Path

def ot(path):
    return gzip.open(path, 'rt', encoding='utf-8', errors='replace') if str(path).endswith('.gz') else open(path, encoding='utf-8', errors='replace')

def sample_from(path):
    return Path(path).name.split('.')[0]

def count_vcf(path):
    records = alleles = snps = indels = 0
    with ot(path) as h:
        for line in h:
            if line.startswith('#'): continue
            f = line.rstrip().split('\t')
            if len(f) < 5: continue
            records += 1
            for alt in f[4].split(','):
                if alt in {'*', '<NON_REF>', '.'}: continue
                alleles += 1
                if len(f[3]) == len(alt) == 1: snps += 1
                else: indels += 1
    return records, alleles, snps, indels

def fasta_stats(path):
    entries = residues = 0
    with open(path, encoding='utf-8', errors='replace') as h:
        for line in h:
            if line.startswith('>'): entries += 1
            else: residues += len(line.strip())
    return entries, residues

def table_rows(path):
    with open(path, encoding='utf-8', errors='replace', newline='') as h:
        return list(csv.DictReader(h, delimiter='\t'))

def count_data_rows(path):
    try:
        with open(path, encoding='utf-8', errors='replace') as h:
            return max(0, sum(1 for line in h if line.strip() and not line.startswith('#')) - 1)
    except Exception:
        return 0

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--samples', required=True)
    p.add_argument('--raw-vcf', nargs='*', default=[])
    p.add_argument('--pass-vcf', nargs='*', default=[])
    p.add_argument('--rna-vcf', nargs='*', default=[])
    p.add_argument('--progression-vcf', nargs='*', default=[])
    p.add_argument('--fusion-table', nargs='*', default=[])
    p.add_argument('--splice-detail', nargs='*', default=[])
    p.add_argument('--variant-fasta', nargs='*', default=[])
    p.add_argument('--fusion-fasta', nargs='*', default=[])
    p.add_argument('--splice-fasta', nargs='*', default=[])
    p.add_argument('--combined-fasta', nargs='*', default=[])
    p.add_argument('--stage-qc', nargs='*', default=[])
    p.add_argument('--external-comparison', nargs='*', default=[])
    p.add_argument('--codon-summary', nargs='*', default=[])
    p.add_argument('--provenance-summary', nargs='*', default=[])
    p.add_argument('--output-prefix', required=True)
    a=p.parse_args()

    with open(a.samples, encoding='utf-8', newline='') as h:
        samples=list(csv.DictReader(h))
    for r in samples:
        r['TK']=(r.get('TK') or r['sample']).strip(); r['Group']=(r.get('Group') or r['sample']).strip(); r['baseline']=(r.get('baseline') or 'false').strip().lower()

    vcf_rows=[]
    for stage, paths in [('raw',a.raw_vcf),('pass',a.pass_vcf),('rna_validated',a.rna_vcf),('progression_nonbaseline_only',a.progression_vcf)]:
        for path in paths:
            rec,alln,snp,indel=count_vcf(path); vcf_rows.append([sample_from(path),stage,Path(path).name,rec,alln,snp,indel])
    with open(a.output_prefix+'.variant_stage_inventory.tsv','w',newline='') as h:
        w=csv.writer(h,delimiter='\t',lineterminator='\n'); w.writerow(['Sample','Stage','File','Records','Alleles','SNPs','Indels']);w.writerows(vcf_rows)

    fasta_rows=[]
    for kind,paths in [('variant',a.variant_fasta),('fusion',a.fusion_fasta),('splice',a.splice_fasta),('combined',a.combined_fasta)]:
        for path in paths:
            entries,res=fasta_stats(path);fasta_rows.append([sample_from(path),kind,Path(path).name,entries,res])
    with open(a.output_prefix+'.fasta_inventory.tsv','w',newline='') as h:
        w=csv.writer(h,delimiter='\t',lineterminator='\n');w.writerow(['Sample','FASTA class','File','Sequences','Amino acid residues']);w.writerows(fasta_rows)

    event_rows=[]
    for kind,paths in [('fusion',a.fusion_table),('validated_splice',a.splice_detail)]:
        for path in paths:event_rows.append([sample_from(path),kind,Path(path).name,count_data_rows(path)])
    with open(a.output_prefix+'.rna_event_inventory.tsv','w',newline='') as h:
        w=csv.writer(h,delimiter='\t',lineterminator='\n');w.writerow(['Sample','Event class','File','Data rows']);w.writerows(event_rows)

    comp=[]
    for path in a.external_comparison:
        for row in table_rows(path): comp.append([sample_from(path),row.get('Metric',''),row.get('Value','')])
    with open(a.output_prefix+'.external_caller_comparison.tsv','w',newline='') as h:
        w=csv.writer(h,delimiter='\t',lineterminator='\n');w.writerow(['Sample','Metric','Value']);w.writerows(comp)

    progression=[r for r in vcf_rows if r[1]=='progression_nonbaseline_only']
    baseline_by_tk=defaultdict(list)
    for row in samples:
        if row['baseline']=='true': baseline_by_tk[row['TK']].append(row['sample'])
    warnings=[]
    for tk, members in defaultdict(list, {x['TK']:[] for x in samples}).items():
        bases=baseline_by_tk.get(tk,[])
        if not bases:warnings.append(f'{tk}: no baseline; progression subtraction skipped')
        elif len(bases)>1:warnings.append(f'{tk}: multiple baselines ({", ".join(bases)}); interpretation is ambiguous')
    if not warnings:warnings=['All subjects have exactly one baseline.']

    overview=['# PGTK comparative biological evidence report','',
              '## Why these outputs extend a basic small-variant caller','',
              '- Raw, hard-filtered and RNA-validated variant stages are retained separately.','- Baseline-aware progression subtraction is reported independently from per-sample FASTA generation.','- RNA fusions and novel splice-derived proteins are captured outside ordinary small-variant VCFs.','- Variant, fusion and splice proteins are assembled into independent per-sample MaxQuant FASTAs.','- Codon, read-provenance, external-caller and MaxQuant evidence can be connected without replacing the underlying caller stages.','',
              '## Samples','']
    overview += [f"- {r['sample']}: SRR={r['srr']}, subject={r['TK']}, group={r['Group']}, baseline={r['baseline']}" for r in samples]
    overview += ['', '## Baseline diagnostics',''] + [f'- {x}' for x in warnings]
    overview += ['', '## Progression outputs','']
    overview += [f'- {r[0]}: {r[4]} non-baseline-only alleles ({r[5]} SNPs, {r[6]} indels)' for r in progression] or ['- No progression VCF was generated.']
    overview += ['', '## External caller comparison','']
    for sample in sorted({x[0] for x in comp}):
        vals={m:v for s,m,v in comp if s==sample}; overview.append(f"- {sample}: shared={vals.get('Shared','NA')}, PGTK overlap={vals.get('PGTK overlap %','NA')}%, genotype concordance={vals.get('Genotype concordance %','NA')}%")
    if not comp:overview.append('- No external caller folder was supplied.')
    Path(a.output_prefix+'.report.md').write_text('\n'.join(overview)+'\n',encoding='utf-8')

    with open(a.output_prefix+'.multiqc_summary.tsv','w',newline='') as h:
        w=csv.writer(h,delimiter='\t',lineterminator='\n');w.writerow(['Sample','Raw alleles','PASS alleles','RNA-validated alleles','Progression alleles','Variant proteins','Fusion proteins','Splice proteins'])
        for s in [x['sample'] for x in samples]:
            va={(r[0],r[1]):r[4] for r in vcf_rows};fa={(r[0],r[1]):r[3] for r in fasta_rows}
            w.writerow([s,va.get((s,'raw'),0),va.get((s,'pass'),0),va.get((s,'rna_validated'),0),va.get((s,'progression_nonbaseline_only'),0),fa.get((s,'variant'),0),fa.get((s,'fusion'),0),fa.get((s,'splice'),0)])
if __name__=='__main__':main()
