#!/usr/bin/env python3
import argparse, html
from pathlib import Path
SECTIONS=[
("sequencing_qc","Sequencing QC outputs",{"qc"},"FASTQ quality was assessed before and after Trim Galore. Reads were aligned to GRCh38 with STAR. Alignment, duplicate and variant-call summaries were collected with Samtools, Picard and Bcftools.","SRA FASTQ files; GRCh38 primary assembly; Ensembl release 111 GTF.",["Trim quality: 20","Minimum trimmed length: 36","STAR unique-read MAPQ: 60"]),
("alignment_expression","Alignment and expression outputs",{"bam","expression"},"STAR coordinate-sorted BAMs were quantified with featureCounts. Per-sample counts were merged and used for per-sample over-representation and baseline-ranked GO analyses.","STAR BAMs; Ensembl 111 annotation; GO ontology; human GOA.",["featureCounts MAPQ: 10","Minimum overlap: 1","Expression CPM: 1.0","Pseudocount: 0.5","GO FDR: 0.1"]),
("raw_variants","Raw and filtered variant outputs",{"gvcf","vcf_raw","vcf_normalized","vcf_snp","vcf_indel","vcf_filtered","vcf_pass"},"RNA BAMs were processed with SplitNCigarReads, scattered GATK HaplotypeCaller, GenotypeGVCFs, normalization, SNP/INDEL separation, hard filtering and PASS selection.","Marked-duplicate split RNA BAMs; GRCh38.",["24 HaplotypeCaller shards","Calling confidence: 20","PCR indel model: CONSERVATIVE","SNP: QD<2, FS>60, SOR>3, MQ<40, MQRankSum<-12.5, ReadPosRankSum<-8","INDEL: QD<2, FS>200, SOR>10, ReadPosRankSum<-20"]),
("annotated_variants","Annotated and RNA-validated variants",{"vep","rna_validation"},"PASS variants were annotated with offline VEP 111. RNA support, codon consistency and supporting-read provenance were assessed with sample-matched BAMs.","PASS VCFs; VEP 111 cache; BAMs; GRCh38.",["RNA depth: 10","ALT reads: 3","ALT fraction: 0.05","MAPQ: 20","Base quality: 20"]),
("variant_landscape","Variant landscape and GO outputs",{"variant_landscape"},"Variant types, genotypes, substitution classes, impacts and consequences were summarized by stage. Protein-altering genes underwent hypergeometric GO testing with Benjamini-Hochberg correction.","Raw through progression-stage VCFs; GO mapping.",["GO term size: 10 to 500","GO FDR: 0.1","Top 100 terms per sample-stage"]),
("progression","Progression comparison outputs",{"progression_vcf","progression_biology"},"Each non-baseline sample was partitioned into nonbaseline-only, baseline-only and shared alleles. Progression samples from the same subject were compared for shared and exclusive alleles and genes, followed by GO analysis.","RNA-validated VCFs; TK subject and baseline fields; GO mapping.",["Exact CHROM:POS:REF:ALT comparison","Exactly one baseline per subject","GO FDR: 0.1"]),
("fusion_splice","Fusion and splice outputs",{"fusions","fusion_fasta","splicing","splice_fasta"},"Arriba fusion calls were support-filtered. StringTie transcripts were compared with Ensembl using gffcompare, splice support was validated, and translated FASTAs were generated.","STAR BAMs; Arriba resources; Ensembl GTF; GRCh38.",["Fusion split reads: 1","Fusion total support: 2","Splice junction reads: 3","StringTie coverage: 2.5","Isoform fraction: 0.05","Protein length: 60 aa"]),
("proteogenomic_fastas","Exploratory proteogenomic FASTAs",{"variant_fasta","combined_fasta"},"Protein-altering variants were translated with pyPGATK and combined with validated fusion and splice proteins. Identical sequences were deduplicated.","RNA-validated VCFs; Ensembl cDNA/GTF; fusion and splice FASTAs.",[]),
("igv","IGV evidence and finding explorer",{"igv"},"RNA variants, progression variants, fusions and splice events were converted to coordinates, event BAMs, IGV sessions and a database-free explorer.","Event tables; sample BAMs; GRCh38.",["MAPQ: 20","Base quality: 20","Padding: 150 bp","Maximum ALT display reads: 100"]),
("reports","Integrated reports",{"reports","comparative_advantage"},"Findings, validation audits, stage inventories and interpretation limits were consolidated into complete and comparative reports.","All core analytical branches.",[]),
("external","External caller comparison",{"comparison"},"When enabled, raw, PASS and RNA-validated alleles were compared with one indexed external VCF per accession.","PGTK VCFs; external VCFs.",[]),
("maxquant","MaxQuant proteogenomic validation",{"proteogenomics_validation"},"When enabled, MaxQuant outputs were mapped to custom and canonical FASTAs and integrated with RNA, codon and read support.","MaxQuant txt, mqpar.xml, FASTAs and PGTK evidence.",[]),
("execution","Execution records and failure logs",{"failure_logs","execution_root"},"Nextflow and wrapper outputs document status, resources, retries and failures.","Nextflow execution metadata.",[])]
def size(n):
 v=float(n)
 for u in ('B','KiB','MiB','GiB','TiB'):
  if v<1024 or u=='TiB': return f'{v:.1f} {u}'
  v/=1024
def main():
 p=argparse.ArgumentParser();p.add_argument('--results-dir',required=True);p.add_argument('--output-dir',required=True);p.add_argument('--external-enabled',choices=['true','false'],required=True);p.add_argument('--maxquant-enabled',choices=['true','false'],required=True);a=p.parse_args()
 root=Path(a.results_dir).resolve();out=Path(a.output_dir);out.mkdir(parents=True,exist_ok=True); groups={}
 for f in sorted(root.rglob("*")):
  if not f.is_file(): continue
  r=f.relative_to(root)
  if r.name=='multiqc_report.html' or r.name.startswith('.'): continue
  g='execution_root' if len(r.parts)==1 else r.parts[0]
  if g=='comparison' and a.external_enabled!='true': continue
  if g=='proteogenomics_validation' and a.maxquant_enabled!='true': continue
  groups.setdefault(g,[]).append((r,f.stat().st_size))
 made=[];covered=set();total=0
 for ident,title,kinds,method,inputs,params in SECTIONS:
  if ident=='external' and a.external_enabled!='true' or ident=='maxquant' and a.maxquant_enabled!='true': continue
  fs=sorted([x for k in kinds for x in groups.get(k,[])],key=lambda x:x[0].as_posix())
  if not fs: continue
  links=''.join(f'<li><a href="../{html.escape(r.as_posix(),quote=True)}">{html.escape(r.as_posix())}</a> <span class="text-muted">({size(n)})</span></li>' for r,n in fs)
  ps=''.join(f'<li>{html.escape(x)}</li>' for x in params) or '<li>See linked method and summary files.</li>'
  text=f'---\nid: pgtk_catalogue_{ident}\nsection_name: {title}\n---\n<h3>{html.escape(title)}</h3><p><strong>How produced:</strong> {html.escape(method)}</p><p><strong>Inputs:</strong> {html.escape(inputs)}</p><p><strong>Key thresholds and settings:</strong></p><ul>{ps}</ul><p><strong>Published files ({len(fs)}):</strong></p><ul>{links}</ul>\n'
  (out/f'pgtk_catalogue_{ident}_mqc.html').write_text(text);made.append((ident,title,len(fs)));covered|=kinds;total+=len(fs)
 other=sorted([x for k,v in groups.items() if k not in covered for x in v],key=lambda x:x[0].as_posix())
 if other:
  links=''.join(f'<li><a href="../{html.escape(r.as_posix(),quote=True)}">{html.escape(r.as_posix())}</a> ({size(n)})</li>' for r,n in other)
  (out/'pgtk_catalogue_other_mqc.html').write_text(f'---\nid: pgtk_catalogue_other\nsection_name: Other published outputs\n---\n<p><strong>Published files ({len(other)}):</strong></p><ul>{links}</ul>');made.append(('other','Other published outputs',len(other)));total+=len(other)
 links=''.join(f'<li><a href="#pgtk_catalogue_{i}">{html.escape(t)}</a> ({n} files)</li>' for i,t,n in made)
 (out/'pgtk_results_catalogue_mqc.html').write_text(f'---\nid: pgtk_results_catalogue\nsection_name: Results catalogue\n---\n<p><strong>{total} published result files</strong> are documented below. Each section explains steps, inputs, thresholds and direct output links.</p><ul>{links}</ul>')
 print(f'Wrote {len(made)} sections covering {total} files')
if __name__=='__main__':main()
