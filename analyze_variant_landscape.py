#!/usr/bin/env python3
import argparse,csv,gzip,json,math,os
from collections import Counter,defaultdict
from pathlib import Path
ALTERING={'missense_variant','frameshift_variant','stop_gained','stop_lost','start_lost','splice_donor_variant','splice_acceptor_variant','inframe_insertion','inframe_deletion','protein_altering_variant'}
VARIANT_METRICS=('alleles','SNV','MNV','insertion','deletion','complex_allele','multiallelic_site','transition','transversion','heterozygous','homozygous_alt','missing_genotype')
def open_text(path):return gzip.open(path,'rt',encoding='utf-8',errors='replace') if str(path).endswith('.gz') else open(path,encoding='utf-8',errors='replace')
def bh(rows):
 order=sorted(range(len(rows)),key=lambda i:rows[i]['PValue']);q=[1.0]*len(rows);prev=1.0;m=len(rows)
 for rank in range(m-1,-1,-1):i=order[rank];prev=min(prev,rows[i]['PValue']*m/(rank+1));q[i]=prev
 for i,v in enumerate(q):rows[i]['FDR']=v
def hypergeom_sf(a,N,K,n):
 lo=max(a,n-(N-K));hi=min(n,K)
 if lo>hi or N<=0:return 1.0
 def logchoose(x,y):return math.lgamma(x+1)-math.lgamma(y+1)-math.lgamma(x-y+1)
 term=math.exp(logchoose(K,lo)+logchoose(N-K,n-lo)-logchoose(N,n));total=term
 for x in range(lo,hi):
  den=(x+1)*(N-K-n+x+1)
  if den<=0:break
  term*=((K-x)*(n-x))/den;total+=term
 return min(1.0,total)
def parse_vcf(stage,path):
 sample=Path(path).name.split('.',1)[0];csq=[];counts=Counter();filters=Counter();consequences=Counter();impacts=Counter();genes=set();pur={('A','G'),('G','A'),('C','T'),('T','C')}
 with open_text(path) as h:
  for line in h:
   if line.startswith('##INFO=<ID=CSQ') and 'Format:' in line:csq=line.split('Format:',1)[1].rsplit('"',1)[0].rstrip('>').split('|')
   if line.startswith('#'):continue
   f=line.rstrip('\n').split('\t')
   if len(f)<8:continue
   ref=f[3].upper();alts=f[4].upper().split(',');filters[f[6]]+=1
   fmt=dict(zip(f[8].split(':'),f[9].split(':'))) if len(f)>9 else {};gt=fmt.get('GT','')
   counts['heterozygous' if gt in {'0/1','1/0','0|1','1|0'} else 'homozygous_alt' if gt in {'1/1','1|1'} else 'missing_genotype' if '.' in gt or not gt else 'other_genotype']+=1
   if len(alts)>1:counts['multiallelic_site']+=1
   for alt in alts:
    counts['alleles']+=1
    if len(ref)==len(alt)==1:counts['SNV']+=1;counts['transition' if (ref,alt) in pur else 'transversion']+=1
    elif len(ref)==len(alt):counts['MNV']+=1
    elif len(ref)<len(alt):counts['insertion' if alt.startswith(ref) else 'complex_allele']+=1
    else:counts['deletion' if ref.startswith(alt) else 'complex_allele']+=1
   info=dict(x.split('=',1) for x in f[7].split(';') if '=' in x)
   for item in filter(None,info.get('CSQ','').split(',')):
    ann=dict(zip(csq,item.split('|')));terms=set(filter(None,ann.get('Consequence','').split('&')));consequences.update(terms);impacts[ann.get('IMPACT','') or 'unannotated']+=1
    if terms&ALTERING and ann.get('SYMBOL'):genes.add(ann['SYMBOL'])
 return sample,counts,filters,consequences,impacts,genes
def read_go(path):
 terms=defaultdict(set);names={};spaces={}
 with open(path,encoding='utf-8',newline='') as h:
  for r in csv.DictReader(h,delimiter='\t'):
   gene=r.get('Gene','');go=r.get('GO_ID','')
   if gene and go:terms[go].add(gene);names[go]=r.get('GO_Name','');spaces[go]=r.get('Namespace','')
 return terms,names,spaces
def write_tsv(path,rows,fields):
 with open(path,'w',encoding='utf-8',newline='') as h:w=csv.DictWriter(h,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def check_output(path,max_rows,max_bytes):
 size=os.path.getsize(path);rows=sum(1 for _ in open(path,encoding='utf-8',errors='replace'))-1
 if rows>max_rows or size>max_bytes:raise RuntimeError(f'Output safeguard exceeded for {path}: {rows} rows, {size} bytes')
 return rows,size
def main():
 p=argparse.ArgumentParser();p.add_argument('--vcf',action='append',required=True,help='STAGE=PATH');p.add_argument('--go-mapping',required=True);p.add_argument('--output-prefix',required=True);p.add_argument('--go-min-size',type=int,default=10);p.add_argument('--go-max-size',type=int,default=500);p.add_argument('--fdr-threshold',type=float,default=.1);p.add_argument('--go-top',type=int,default=100);p.add_argument('--max-significant-rows',type=int,default=100000);p.add_argument('--max-output-mb',type=int,default=100);a=p.parse_args()
 summary=[];gene_rows=[];significant=[];top=[];go_summary=[];terms,names,spaces=read_go(a.go_mapping);universe=set().union(*terms.values()) if terms else set();eligible={go:(genes&universe) for go,genes in terms.items() if a.go_min_size<=len(genes&universe)<=a.go_max_size}
 for spec in a.vcf:
  stage,path=spec.split('=',1);sample,c,filters,cons,imp,genes=parse_vcf(stage,path)
  for metric in VARIANT_METRICS:summary.append({'Sample':sample,'Stage':stage,'Category':'variant_type','Metric':metric,'Count':c[metric]})
  for k,v in sorted(filters.items()):summary.append({'Sample':sample,'Stage':stage,'Category':'filter','Metric':k,'Count':v})
  for k,v in sorted(imp.items()):summary.append({'Sample':sample,'Stage':stage,'Category':'VEP_impact','Metric':k,'Count':v})
  for k,v in sorted(cons.items()):summary.append({'Sample':sample,'Stage':stage,'Category':'VEP_consequence','Metric':k,'Count':v})
  for g in sorted(genes):gene_rows.append({'Sample':sample,'Stage':stage,'Gene':g})
  fg=genes&universe;rows=[]
  if fg:
   for go,term_genes in eligible.items():
    overlap=fg&term_genes
    if not overlap:continue
    pv=hypergeom_sf(len(overlap),len(universe),len(term_genes),len(fg));od=((len(overlap)+.5)*(len(universe)-len(fg)-len(term_genes)+len(overlap)+.5))/((len(fg)-len(overlap)+.5)*(len(term_genes)-len(overlap)+.5))
    rows.append({'Sample':sample,'Stage':stage,'GO_ID':go,'GO_Name':names[go],'Namespace':spaces[go],'ProteinAlteringGenes':len(fg),'BackgroundGenes':len(universe),'OverlapGenesInTerm':len(overlap),'TermGenesInBackground':len(term_genes),'OddsRatio':od,'PValue':pv,'FDR':1.0,'OverlapGenes':';'.join(sorted(overlap))})
   bh(rows);ranked=sorted(rows,key=lambda r:(r['FDR'],r['PValue'],-r['OddsRatio'],r['GO_ID']))
   significant.extend(r for r in ranked if r['FDR']<=a.fdr_threshold);top.extend(ranked[:a.go_top])
  go_summary.append({'Sample':sample,'Stage':stage,'ProteinAlteringGenes':len(fg),'EligibleGOTerms':len(eligible),'TermsWithOverlap':len(rows),'SignificantTerms':sum(r['FDR']<=a.fdr_threshold for r in rows),'TopTermsWritten':min(len(rows),a.go_top),'FDRThreshold':a.fdr_threshold})
 prefix=a.output_prefix;go_fields=['Sample','Stage','GO_ID','GO_Name','Namespace','ProteinAlteringGenes','BackgroundGenes','OverlapGenesInTerm','TermGenesInBackground','OddsRatio','PValue','FDR','OverlapGenes']
 write_tsv(prefix+'.summary.tsv',summary,['Sample','Stage','Category','Metric','Count']);write_tsv(prefix+'.nonsynonymous_genes.tsv',gene_rows,['Sample','Stage','Gene']);write_tsv(prefix+'.go_significant.tsv',significant,go_fields);write_tsv(prefix+'.go_top.tsv',top,go_fields);write_tsv(prefix+'.go_summary.tsv',go_summary,['Sample','Stage','ProteinAlteringGenes','EligibleGOTerms','TermsWithOverlap','SignificantTerms','TopTermsWritten','FDRThreshold'])
 checks=[]
 for path,limit in [(prefix+'.go_significant.tsv',a.max_significant_rows),(prefix+'.go_top.tsv',len(a.vcf)*a.go_top),(prefix+'.go_summary.tsv',len(a.vcf))]:checks.append((path,*check_output(path,limit,a.max_output_mb*1024*1024)))
 Path(prefix+'.report.md').write_text('# Variant landscape and nonsynonymous GO analysis\n\nCounts are separated by sample and stage. GO over-representation uses unique genes with protein-altering VEP consequences, a deduplicated GO-mapped human gene universe, one-sided hypergeometric tests and Benjamini-Hochberg FDR. Only significant terms and the top ranked terms are written. RNA-derived calls are not DNA-confirmed somatic mutations.\n\nOutput safeguards: '+', '.join(f'{Path(p).name}: {n} rows, {b} bytes' for p,n,b in checks)+'\n',encoding='utf-8')
 out=Path(prefix+'.multiqc');out.mkdir(exist_ok=True);data=defaultdict(dict)
 for r in summary:
  if r['Category']=='variant_type' and r['Metric'] in {'SNV','MNV','insertion','deletion','complex_allele'}:data[f"{r['Sample']} | {r['Stage']}"][r['Metric']]=r['Count']
 (out/'pgtk_variant_00_types_mqc.json').write_text(json.dumps({'id':'pgtk_variant_types','section_name':'Variant types by sample and stage','description':'Allele classes before and after filtering. Stages are not biologically equivalent.','plot_type':'bargraph','pconfig':{'id':'pgtk_variant_types','title':'Variant types by sample and stage','ylab':'Alleles'},'data':data},indent=2))
 summary_data={f"{r['Sample']} | {r['Stage']}":{'Protein-altering genes':r['ProteinAlteringGenes'],'GO terms with overlap':r['TermsWithOverlap'],'Significant GO terms':r['SignificantTerms']} for r in go_summary if r['Stage'] in {'rna_validated','progression_nonbaseline_only','progression_baseline_only','progression_shared_with_baseline'}}
 (out/'pgtk_variant_02_go_summary_mqc.json').write_text(json.dumps({'id':'pgtk_rna_variant_go_summary','section_name':'RNA-seq Protein-Altering Variant GO summary','description':'Protein-altering VEP genes tested by sample and RNA-derived stage. This is distinct from expression GO and progression-set GO. RNA-derived variants are not DNA-confirmed.','plot_type':'bargraph','pconfig':{'id':'pgtk_rna_variant_go_summary','title':'RNA-seq protein-altering variant GO summary','ylab':'Count'},'data':summary_data},indent=2))
 selected_stages={'rna_validated':'RNA-validated','progression_nonbaseline_only':'Progression nonbaseline-only'}
 for stage,label in selected_stages.items():
  for sample in sorted({r['Sample'] for r in top if r['Stage']==stage}):
   rows=[r for r in top if r['Stage']==stage and r['Sample']==sample][:15]
   plot_data={f"{r['GO_Name']} ({r['GO_ID']})":{'-log10 FDR':min(50.0,-math.log10(max(float(r['FDR']),1e-50)))} for r in reversed(rows)}
   ident=f"pgtk_variant_go_{stage}_{sample}".lower()
   (out/f'pgtk_variant_03_{ident}_mqc.json').write_text(json.dumps({'id':ident,'section_name':f'{label} GO: {sample}','description':'Top protein-altering variant GO terms ranked by FDR. Terms overlap and are not independent pathways.','plot_type':'bargraph','pconfig':{'id':ident,'title':f'{label} protein-altering variant GO: {sample}','ylab':'-log10(FDR)'},'data':plot_data},indent=2))
 guide='<h4>RNA-seq Protein-Altering Variant GO</h4><p><b>Different evidence layers:</b> expression GO uses gene abundance; progression-set GO uses baseline-subtracted candidate gene sets; this section uses genes carrying protein-altering VEP consequences in RNA-derived variant stages.</p><p><a href="../variant_landscape/variant_landscape.go_significant.tsv">Significant terms</a> | <a href="../variant_landscape/variant_landscape.go_top.tsv">Top terms</a> | <a href="../variant_landscape/variant_landscape.go_summary.tsv">Summary</a> | <a href="../variant_landscape/variant_landscape.report.md">Methods</a></p>'
 (out/'pgtk_variant_01_go_guide_mqc.html').write_text(guide,encoding='utf-8')
if __name__=='__main__':main()
