#!/usr/bin/env python3
import argparse,csv,gzip,math
from collections import defaultdict
IMPACT={'MODIFIER':0,'LOW':1,'MODERATE':2,'HIGH':3}
CR={'transcript_ablation':100,'splice_acceptor_variant':95,'splice_donor_variant':95,'stop_gained':90,'frameshift_variant':90,'stop_lost':85,'start_lost':85,'inframe_insertion':75,'inframe_deletion':75,'missense_variant':70,'protein_altering_variant':65,'splice_region_variant':60,'synonymous_variant':30}
def op(p):return gzip.open(p,'rt',encoding='utf-8',errors='replace') if str(p).endswith('.gz') else open(p,encoding='utf-8',errors='replace')
def fields(p):
 with op(p) as h:
  for l in h:
   if l.startswith('##INFO=<ID=CSQ') and 'Format:' in l:return l.split('Format:',1)[1].rsplit('"',1)[0].rstrip('>').split('|')
 raise SystemExit('CSQ header missing: '+p)
def best(xs):return max(xs,key=lambda x:(IMPACT.get(x.get('IMPACT',''),0),max([CR.get(y,0) for y in x.get('Consequence','').split('&')] or [0]),x.get('CANONICAL')=='YES',x.get('BIOTYPE')=='protein_coding'))
def vcf(p):
 fs=fields(p); out=[]
 with op(p) as h:
  for l in h:
   if l.startswith('#'):continue
   c=l.rstrip().split('\t'); info=dict(x.split('=',1) for x in c[7].split(';') if '=' in x); cs=[]
   for raw in info.get('CSQ','').split(','):
    v=raw.split('|')+['']*len(fs);cs.append(dict(zip(fs,v)))
   fmt=dict(zip(c[8].split(':'),c[9].split(':'))) if len(c)>9 else {}
   for alt in c[4].split(','):
    if alt in {'.','*','<NON_REF>'}:continue
    ms=[x for x in cs if x.get('Allele') in {alt,alt.lstrip(c[3][:1]) or alt}] or cs; bg=defaultdict(list)
    for x in ms:
     g=(x.get('SYMBOL') or x.get('Gene') or '').strip()
     if g:bg[g].append(x)
    anns=[best(x) for x in bg.values()]; pri=best(anns) if anns else {}
    out.append({'Chrom':c[0].removeprefix('chr'),'Pos':int(c[1]),'Ref':c[3],'Alt':alt,'Primary':pri,'AllGenes':';'.join(sorted(bg)),'GT':fmt.get('GT',''),'DP':fmt.get('DP',''),'AD':fmt.get('AD','')})
 return out
def fisher(a,b,c,d):
 n=a+b+c+d;r=a+b;k=a+c;den=math.comb(n,r);p=sum(math.comb(k,x)*math.comb(n-k,r-x)/den for x in range(a,min(r,k)+1) if 0<=r-x<=n-k);o=((a+.5)*(d+.5))/((b+.5)*(c+.5));return min(1,p),o
def bh(rows):
 order=sorted(range(len(rows)),key=lambda i:rows[i]['PValue']);prev=1.;adj=[1.]*len(rows)
 for z in range(len(order)-1,-1,-1):i=order[z];prev=min(prev,rows[i]['PValue']*len(rows)/(z+1));adj[i]=prev
 for r,x in zip(rows,adj):r['FDR']=x
def wr(p,rows,fs):
 with open(p,'w',encoding='utf-8',newline='') as h:w=csv.DictWriter(h,fieldnames=fs,delimiter='\t',lineterminator='\n',extrasaction='ignore');w.writeheader();w.writerows(rows)
def main():
 p=argparse.ArgumentParser();p.add_argument('--sample',required=True);p.add_argument('--subject',required=True);p.add_argument('--group',required=True);p.add_argument('--progression-vcf',required=True);p.add_argument('--background-vcf',required=True);p.add_argument('--go-mapping',required=True);p.add_argument('--go-min-size',type=int,default=10);p.add_argument('--go-max-size',type=int,default=500);p.add_argument('--output-prefix',required=True);a=p.parse_args()
 prog=vcf(a.progression_vcf);back=vcf(a.background_vcf);ars=[]
 for r in prog:
  x=r['Primary'];ars.append({'Sample':a.sample,'Subject':a.subject,'Group':a.group,'Chrom':r['Chrom'],'Pos':r['Pos'],'Ref':r['Ref'],'Alt':r['Alt'],'Gene':x.get('SYMBOL') or x.get('Gene',''),'AllGenes':r['AllGenes'],'Impact':x.get('IMPACT',''),'Consequence':x.get('Consequence',''),'Canonical':x.get('CANONICAL',''),'Feature':x.get('Feature',''),'HGVSc':x.get('HGVSc',''),'HGVSp':x.get('HGVSp',''),'GT':r['GT'],'DP':r['DP'],'AD':r['AD']})
 fg={r['Gene'] for r in ars if r['Gene']};univ={r['Primary'].get('SYMBOL') or r['Primary'].get('Gene','') for r in back};univ.discard('');fg &= univ
 go=defaultdict(set);names={};spaces={}
 with open(a.go_mapping,encoding='utf-8',newline='') as h:
  for r in csv.DictReader(h,delimiter='\t'):go[r['GO_ID']].add(r['Gene']);names[r['GO_ID']]=r['GO_Name'];spaces[r['GO_ID']]=r['Namespace']
 ers=[]
 for gid,genes in go.items():
  tested=genes&univ
  if not a.go_min_size<=len(tested)<=a.go_max_size:continue
  overlap=fg&tested;a1=len(overlap);b=len(fg-tested);c=len(tested-fg);d=len(univ-(fg|tested));pv,od=fisher(a1,b,c,d)
  ers.append({'Sample':a.sample,'Subject':a.subject,'Group':a.group,'GO_ID':gid,'GO_Name':names[gid],'Namespace':spaces[gid],'ProgressionGenesInTerm':a1,'ProgressionGenesTested':len(fg),'BackgroundGenesInTerm':len(tested),'BackgroundGenesTested':len(univ),'OddsRatio':od,'PValue':pv,'OverlapGenes':';'.join(sorted(overlap))})
 bh(ers);ers.sort(key=lambda r:(r['FDR'],r['PValue'],r['GO_ID']))
 gr=[]
 for g in sorted(fg):
  z=[r for r in ars if r['Gene']==g];gr.append({'Sample':a.sample,'Subject':a.subject,'Group':a.group,'Gene':g,'BestImpact':max((r['Impact'] for r in z),key=lambda x:IMPACT.get(x,0)),'Alleles':len(z)})
 cand=[{'PriorityScore':IMPACT.get(r['Impact'],0)*10,**r} for r in ars];cand.sort(key=lambda r:(-r['PriorityScore'],r['Gene'],int(r['Pos'])))
 pre=a.output_prefix;wr(pre+'.alleles.tsv',ars,['Sample','Subject','Group','Chrom','Pos','Ref','Alt','Gene','AllGenes','Impact','Consequence','Canonical','Feature','HGVSc','HGVSp','GT','DP','AD']);wr(pre+'.genes.tsv',gr,['Sample','Subject','Group','Gene','BestImpact','Alleles']);wr(pre+'.go_enrichment.tsv',ers,['Sample','Subject','Group','GO_ID','GO_Name','Namespace','ProgressionGenesInTerm','ProgressionGenesTested','BackgroundGenesInTerm','BackgroundGenesTested','OddsRatio','PValue','FDR','OverlapGenes']);wr(pre+'.candidates.tsv',cand,['PriorityScore','Sample','Subject','Group','Chrom','Pos','Ref','Alt','Gene','AllGenes','Impact','Consequence','Canonical','Feature','HGVSc','HGVSp','GT','DP','AD']);wr(pre+'.summary.tsv',[{'Sample':a.sample,'Subject':a.subject,'Group':a.group,'ProgressionAlleles':len(prog),'ProgressionGenes':len(fg),'GO_TermsTested':len(ers),'SignificantGOTermsFDR05':sum(r['FDR']<=.05 for r in ers)}],['Sample','Subject','Group','ProgressionAlleles','ProgressionGenes','GO_TermsTested','SignificantGOTermsFDR05'])
if __name__=='__main__':main()
