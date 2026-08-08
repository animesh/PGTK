#!/usr/bin/env python3
import csv,gzip,subprocess,sys,tempfile
from pathlib import Path
R=Path(__file__).resolve().parent
CSQ='Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|HGVSc|HGVSp|CANONICAL'
def ann(alt,cons,imp,g):return '|'.join([alt,cons,imp,g,'ENSG'+g,'Transcript','ENST'+g,'protein_coding',g+':c.1A>G',g+':p.X1Y','YES'])
def vcf(p,s,rows):
 with gzip.open(p,'wt') as h:
  h.write('##fileformat=VCFv4.2\n##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: '+CSQ+'">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+s+'\n')
  for pos,g in rows:h.write(f'1\t{pos}\t.\tA\tG\t.\tPASS\tCSQ={ann("G","missense_variant","MODERATE",g)}\tGT:DP:AD\t0/1:20:12,8\n')
def rd(p):
 with open(p,newline='') as h:return list(csv.DictReader(h,delimiter='\t'))
with tempfile.TemporaryDirectory() as d:
 w=Path(d);obo=w/'go.obo';gaf=w/'go.gaf.gz'
 obo.write_text('format-version: 1.2\ndata-version: releases/fixture\n\n[Term]\nid: GO:0000001\nname: parent process\nnamespace: biological_process\n\n[Term]\nid: GO:0000002\nname: child process\nnamespace: biological_process\nis_a: GO:0000001 ! parent\n\n[Term]\nid: GO:0000003\nname: molecular activity\nnamespace: molecular_function\n',encoding='utf-8')
 with gzip.open(gaf,'wt') as h:
  h.write('!gaf-version: 2.2\n!generated-on: 2026-08-08\n')
  for g,go,asp in [('GENE1','GO:0000002','P'),('GENE2','GO:0000002','P'),('GENE3','GO:0000003','F'),('BG1','GO:0000001','P'),('BG2','GO:0000001','P')]:h.write(f'UniProtKB\tP\t{g}\t\t{go}\tPMID:1\tEXP\t\t{asp}\tX\t\tprotein\ttaxon:9606\t20260808\tTEST\t\t\n')
 subprocess.run([sys.executable,str(R/'prepare_go_annotations.py'),'--obo',str(obo),'--gaf',str(gaf),'--output-prefix',str(w/'go')],check=True)
 mapping=rd(w/'go.mapping.tsv');assert any(r['Gene']=='GENE1' and r['GO_ID']=='GO:0000001' for r in mapping)
 outs={}
 for s,rows in {'P1':[(10,'GENE1'),(20,'GENE3')],'P2':[(10,'GENE1'),(30,'GENE2')]}.items():
  pv=w/(s+'.p.vcf.gz');bv=w/(s+'.b.vcf.gz');vcf(pv,s,rows);vcf(bv,s,rows+[(40,'BG1'),(50,'BG2')]);pre=w/(s+'.x')
  subprocess.run([sys.executable,str(R/'analyze_progression_biology.py'),'--sample',s,'--subject','T1','--group',s,'--progression-vcf',str(pv),'--background-vcf',str(bv),'--go-mapping',str(w/'go.mapping.tsv'),'--go-min-size','1','--go-max-size','10','--output-prefix',str(pre)],check=True);outs[s]=pre
 e=rd(str(outs['P1'])+'.go_enrichment.tsv');assert {'GO:0000001','GO:0000002','GO:0000003'}<=set(r['GO_ID'] for r in e)
 pp=w/'pair';subprocess.run([sys.executable,str(R/'compare_progression_pair.py'),'--subject','T1','--sample-a','P1','--sample-b','P2','--alleles-a',str(outs['P1'])+'.alleles.tsv','--alleles-b',str(outs['P2'])+'.alleles.tsv','--genes-a',str(outs['P1'])+'.genes.tsv','--genes-b',str(outs['P2'])+'.genes.tsv','--go-a',str(outs['P1'])+'.go_enrichment.tsv','--go-b',str(outs['P2'])+'.go_enrichment.tsv','--output-prefix',str(pp)],check=True)
 assert rd(str(pp)+'.summary.tsv')[0]['SharedAlleles']=='1';assert len(rd(str(pp)+'.go_contrasts.tsv'))==2
 print('PASS: versioned GO OBO/GAF propagation, sample ORA and independent pair contrast')
