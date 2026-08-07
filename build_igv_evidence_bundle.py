#!/usr/bin/env python3
import argparse,csv,gzip,html,subprocess,xml.etree.ElementTree as ET
from pathlib import Path

def ot(p):return gzip.open(p,'rt',errors='replace') if str(p).endswith('.gz') else open(p,errors='replace')
def sample(p):return Path(p).name.split('.')[0]
def main():
 p=argparse.ArgumentParser();p.add_argument('--genome',required=True);p.add_argument('--rna-vcf',nargs='*',default=[]);p.add_argument('--progression-vcf',nargs='*',default=[]);p.add_argument('--fusion-table',nargs='*',default=[]);p.add_argument('--splice-table',nargs='*',default=[]);p.add_argument('--bam',action='append',default=[]);p.add_argument('--padding',type=int,default=100);p.add_argument('--output-prefix',required=True);a=p.parse_args()
 bams=dict(x.split('=',1) for x in a.bam);events=[];eid=0
 for cls,paths in [('rna_variant',a.rna_vcf),('progression_variant',a.progression_vcf)]:
  for path in paths:
   s=sample(path)
   with ot(path) as h:
    for l in h:
     if l.startswith('#'):continue
     f=l.rstrip().split('\t');eid+=1;start=max(0,int(f[1])-1);end=start+max(1,len(f[3]));events.append([f'E{eid:08d}',s,cls,f[0],start,end,f'{f[3]}>{f[4]}',Path(path).name])
 for path in a.fusion_table:
  s=sample(path)
  with open(path,errors='replace',newline='') as h:
   for r in csv.DictReader(h,delimiter='\t'):
    c1=r.get('breakpoint1','');c2=r.get('breakpoint2','')
    if ':' not in c1 or ':' not in c2:continue
    ch1,p1=c1.rsplit(':',1);ch2,p2=c2.rsplit(':',1);eid+=1;events.append([f'E{eid:08d}',s,'fusion',ch1,max(0,int(p1)-1),int(p1),c2,Path(path).name])
 for path in a.splice_table:
  s=sample(path)
  with open(path,errors='replace',newline='') as h:
   for r in csv.DictReader(h,delimiter='\t'):
    chrom=r.get('chrom') or r.get('Chromosome') or r.get('seqname');start=r.get('start') or r.get('Start');end=r.get('end') or r.get('End')
    if not chrom or not start or not end:continue
    eid+=1;events.append([f'E{eid:08d}',s,'splice',chrom,max(0,int(float(start))-1),int(float(end)),r.get('transcript_id',''),Path(path).name])
 with open(a.output_prefix+'.events.tsv','w',newline='') as h:w=csv.writer(h,delimiter='\t',lineterminator='\n');w.writerow(['Event','Sample','Class','Chrom','Start0','End','Label','Source']);w.writerows(events)
 with open(a.output_prefix+'.events.bed','w',newline='') as h:
  w=csv.writer(h,delimiter='\t',lineterminator='\n');[w.writerow([x[3],x[4],x[5],x[0],0,'.']) for x in events if x[3]]
 regions={}
 for s in bams:
  rs=[]
  for x in events:
   if x[1]==s:rs.append(f'{x[3]}:{max(1,x[4]+1-a.padding)}-{x[5]+a.padding}')
  regions[s]=sorted(set(rs))
 manifest=[]
 for s,bam in bams.items():
  out=f'{a.output_prefix}.{s}.events.bam';reg=regions.get(s,[])
  if reg:
   subprocess.run(['samtools','view','-b','-o',out,bam,*reg],check=True);subprocess.run(['samtools','index',out],check=True)
  else:
   subprocess.run(['samtools','view','-H','-b','-o',out,bam],check=True);subprocess.run(['samtools','index',out],check=True)
  manifest.append([s,bam,out,out+'.bai',len(reg)])
 with open(a.output_prefix+'.sample_manifest.tsv','w',newline='') as h:w=csv.writer(h,delimiter='\t',lineterminator='\n');w.writerow(['Sample','Source BAM','Event BAM','Index','Regions']);w.writerows(manifest)
 with open(a.output_prefix+'.igv.batch.txt','w') as h:
  h.write('new\n');h.write(f'genome {Path(a.genome).resolve()}\n')
  for _,_,out,_,_ in manifest:h.write(f'load {out}\n')
  h.write(f'load {a.output_prefix}.events.bed\n')
  for x in events[:1000]:h.write(f'goto {x[3]}:{x[4]+1}-{x[5]}\n')
  h.write('exit\n')
 root=ET.Element('Session',genome=str(Path(a.genome).resolve()),version='8');res=ET.SubElement(root,'Resources')
 for _,_,out,_,_ in manifest:ET.SubElement(res,'Resource',path=out)
 ET.SubElement(res,'Resource',path=f'{a.output_prefix}.events.bed');ET.ElementTree(root).write(a.output_prefix+'.igv.session.xml',encoding='utf-8',xml_declaration=True)
 Path(a.output_prefix+'.summary.txt').write_text(f'Events: {len(events)}\nSamples: {len(bams)}\nRNA variants: {sum(x[2]=="rna_variant" for x in events)}\nProgression variants: {sum(x[2]=="progression_variant" for x in events)}\nFusions: {sum(x[2]=="fusion" for x in events)}\nSplice events: {sum(x[2]=="splice" for x in events)}\n')
if __name__=='__main__':main()
