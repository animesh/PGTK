#!/usr/bin/env python3
import argparse,csv
from pathlib import Path
def many(ps):
 rs=[];fs=[]
 for p in ps:
  with open(p,encoding='utf-8',newline='') as h:r=csv.DictReader(h,delimiter='\t');fs=fs or list(r.fieldnames or []);rs.extend(r)
 return rs,fs
def wr(p,rs,fs):
 with open(p,'w',encoding='utf-8',newline='') as h:w=csv.DictWriter(h,fieldnames=fs,delimiter='\t',lineterminator='\n',extrasaction='ignore');w.writeheader();w.writerows(rs)
def main():
 p=argparse.ArgumentParser()
 for x in ['sample-alleles','sample-genes','sample-go','sample-candidates','sample-summary','pair-alleles','pair-genes','pair-go','pair-summary']:p.add_argument('--'+x,nargs='*',default=[])
 p.add_argument('--go-metadata',required=True);p.add_argument('--output-prefix',required=True);a=p.parse_args();defs=[('sample_alleles','progression_alleles.tsv'),('sample_genes','progression_genes.tsv'),('sample_go','go_enrichment.tsv'),('sample_candidates','candidate_priority.tsv'),('pair_alleles','pairwise_allele_contrasts.tsv'),('pair_genes','pairwise_gene_contrasts.tsv'),('pair_go','pairwise_go_contrasts.tsv'),('pair_summary','pairwise_summary.tsv')];m={}
 for at,su in defs:rs,fs=many(getattr(a,at));wr(a.output_prefix+'.'+su,rs,fs);m[at]=rs
 ss,sf=many(a.sample_summary);pc={}
 for r in m['pair_summary']:
  for x in ('SampleA','SampleB'):pc[r[x]]=pc.get(r[x],0)+1
 mq=[{**r,'PairwiseComparisons':pc.get(r['Sample'],0)} for r in ss];wr(a.output_prefix+'.multiqc_summary.tsv',mq,sf+['PairwiseComparisons'])
 Path(a.output_prefix+'.go_metadata.tsv').write_text(Path(a.go_metadata).read_text(),encoding='utf-8')
 lines=['# Progression biology with Gene Ontology over-representation analysis','','GO terms are derived from a versioned OBO ontology and GAF annotation file. Direct annotations were propagated through is_a and part_of relations. Fisher tests use each sample RNA-callable primary-gene set as background and Benjamini-Hochberg correction across all eligible GO terms.','','## Samples','']
 for r in sorted(ss,key=lambda x:(x['Subject'],x['Sample'])):lines.append(f"- {r['Sample']} ({r['Subject']}): {r['ProgressionAlleles']} alleles, {r['ProgressionGenes']} genes, {r['GO_TermsTested']} GO terms tested, {r['SignificantGOTermsFDR05']} significant at FDR <= 0.05")
 lines+=['','## Pairwise contrasts','']
 for r in sorted(m['pair_summary'],key=lambda x:(x['Subject'],x['SampleA'],x['SampleB'])):lines.append(f"- {r['Subject']}: {r['SampleA']} versus {r['SampleB']}: shared={r['SharedAlleles']}, {r['SampleA']}-only={r['SampleAOnlyAlleles']}, {r['SampleB']}-only={r['SampleBOnlyAlleles']}")
 lines+=['','## Limits','','- This is over-representation analysis, not a replicated classical difference-in-differences model.','- RNA coverage and allele-specific expression can create apparent baseline absence.','- GO result interpretation must use the recorded ontology and GAF versions.',''];Path(a.output_prefix+'.report.md').write_text('\n'.join(lines),encoding='utf-8')
if __name__=='__main__':main()
