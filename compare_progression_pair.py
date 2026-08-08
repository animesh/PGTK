#!/usr/bin/env python3
import argparse,csv,math
def rd(p):
 with open(p,encoding='utf-8',newline='') as h:return list(csv.DictReader(h,delimiter='\t'))
def wr(p,rs,fs):
 with open(p,'w',encoding='utf-8',newline='') as h:w=csv.DictWriter(h,fieldnames=fs,delimiter='\t',lineterminator='\n',extrasaction='ignore');w.writeheader();w.writerows(rs)
def main():
 p=argparse.ArgumentParser()
 for x in ['subject','sample-a','sample-b','alleles-a','alleles-b','genes-a','genes-b','go-a','go-b','output-prefix']:p.add_argument('--'+x,required=True)
 a=p.parse_args();kf=('Chrom','Pos','Ref','Alt');ra=rd(a.alleles_a);rb=rd(a.alleles_b);aa={tuple(r[x] for x in kf):r for r in ra};ab={tuple(r[x] for x in kf):r for r in rb};ar=[]
 for lab,ks in [('shared',set(aa)&set(ab)),(a.sample_a+'_only',set(aa)-set(ab)),(a.sample_b+'_only',set(ab)-set(aa))]:
  for k in sorted(ks,key=lambda x:(x[0],int(x[1]),x[2],x[3])):r=aa.get(k) or ab[k];ar.append({'Subject':a.subject,'SampleA':a.sample_a,'SampleB':a.sample_b,'ContrastClass':lab,**dict(zip(kf,k)),'Gene':r.get('Gene',''),'Impact':r.get('Impact',''),'Consequence':r.get('Consequence','')})
 ga={r['Gene']:r for r in rd(a.genes_a)};gb={r['Gene']:r for r in rd(a.genes_b)};gr=[]
 for lab,gs in [('shared',set(ga)&set(gb)),(a.sample_a+'_only',set(ga)-set(gb)),(a.sample_b+'_only',set(gb)-set(ga))]:
  for g in sorted(gs):gr.append({'Subject':a.subject,'SampleA':a.sample_a,'SampleB':a.sample_b,'ContrastClass':lab,'Gene':g})
 ea={r['GO_ID']:r for r in rd(a.go_a)};eb={r['GO_ID']:r for r in rd(a.go_b)};er=[]
 for gid in sorted(set(ea)&set(eb)):
  x,y=ea[gid],eb[gid];oa=float(x['OddsRatio']);ob=float(y['OddsRatio']);z=math.log2(oa/ob)
  er.append({'Subject':a.subject,'SampleA':a.sample_a,'SampleB':a.sample_b,'GO_ID':gid,'GO_Name':x['GO_Name'],'Namespace':x['Namespace'],'SampleAOddsRatio':oa,'SampleAFDR':x['FDR'],'SampleBOddsRatio':ob,'SampleBFDR':y['FDR'],'Log2OddsRatioContrast':z,'Interpretation':a.sample_a if z>0 else a.sample_b if z<0 else 'equal'})
 su=[{'Subject':a.subject,'SampleA':a.sample_a,'SampleB':a.sample_b,'SharedAlleles':len(set(aa)&set(ab)),'SampleAOnlyAlleles':len(set(aa)-set(ab)),'SampleBOnlyAlleles':len(set(ab)-set(aa)),'SharedGenes':len(set(ga)&set(gb)),'SampleAOnlyGenes':len(set(ga)-set(gb)),'SampleBOnlyGenes':len(set(gb)-set(ga))}]
 pre=a.output_prefix;wr(pre+'.alleles.tsv',ar,['Subject','SampleA','SampleB','ContrastClass','Chrom','Pos','Ref','Alt','Gene','Impact','Consequence']);wr(pre+'.genes.tsv',gr,['Subject','SampleA','SampleB','ContrastClass','Gene']);wr(pre+'.go_contrasts.tsv',er,['Subject','SampleA','SampleB','GO_ID','GO_Name','Namespace','SampleAOddsRatio','SampleAFDR','SampleBOddsRatio','SampleBFDR','Log2OddsRatioContrast','Interpretation']);wr(pre+'.summary.tsv',su,['Subject','SampleA','SampleB','SharedAlleles','SampleAOnlyAlleles','SampleBOnlyAlleles','SharedGenes','SampleAOnlyGenes','SampleBOnlyGenes'])
if __name__=='__main__':main()
