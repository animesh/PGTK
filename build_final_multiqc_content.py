#!/usr/bin/env python3
import argparse,html,json,re
from pathlib import Path

def text(path): return Path(path).read_text(encoding='utf-8',errors='replace') if Path(path).exists() else ''
def metrics(path):
    out={}
    for line in text(path).splitlines():
        m=re.match(r'\s*[-*]?\s*([^:\t|]{2,80})\s*[:\t|]\s*([0-9][0-9,]*(?:\.\d+)?)\s*$',line)
        if m: out[m.group(1).strip()]=float(m.group(2).replace(',',''))
    return out
def plot(out,ident,title,data,desc):
    if not data:return
    payload={'id':ident,'section_name':title,'description':desc,'plot_type':'bargraph','pconfig':{'id':ident,'title':title,'ylab':'Count'},'data':data}
    (out/f'{ident}_mqc.json').write_text(json.dumps(payload,indent=2),encoding='utf-8')
def main():
    p=argparse.ArgumentParser();p.add_argument('--output-dir',required=True)
    for n in ['variant-codon-summary','codon-mismatch-summary','integrated-report','provenance-summary','proteogenomics-summary','classification-report','read-summary','validation-semantics']:p.add_argument('--'+n,required=True)
    a=p.parse_args();out=Path(a.output_dir);out.mkdir(parents=True,exist_ok=True)
    pg=metrics(a.proteogenomics_summary); integ=metrics(a.integrated_report); read=metrics(a.read_summary); codon=metrics(a.variant_codon_summary); prov=metrics(a.provenance_summary)
    preferred=['Input peptides','Mapped peptides','Candidate peptides','Sample-matched direct MS/MS variant events','Sample-matched MBR-only variant events','Strict integrated events','Strict events']
    funnel={k:pg.get(k,integ.get(k)) for k in preferred if pg.get(k,integ.get(k)) is not None}
    if not funnel:
        funnel={k:v for k,v in list(pg.items())[:12]}
    plot(out,'10_pgtk_proteogenomics_funnel','MaxQuant evidence funnel',{'Evidence':funnel},'Summary of peptide mapping and strict proteogenomic evidence. Full evidence tables remain linked, not embedded.')
    plot(out,'11_pgtk_read_validation','Read-level validation',{'Observed':dict(list(read.items())[:12])},'Compact read-evidence summary. Full read-level TSV and IGV bundles remain outside MultiQC.')
    validation={**dict(list(codon.items())[:8]),**dict(list(prov.items())[:8])}
    plot(out,'12_pgtk_validation_summary','Independent validation summary',{'Validated':validation},'Codon and read-provenance validation summary.')
    links='<h3>Methods and detailed evidence</h3><p>This dashboard contains summaries only.</p><ul><li><a href="../proteogenomics_validation/proteogenomics_evidence.report.md">Proteogenomics report</a></li><li><a href="../proteogenomics_validation/integrated_variant_evidence.report.md">Strict integrated evidence</a></li><li><a href="../proteogenomics_validation/read_validation/proteogenomic_read_validation.report.md">Read validation</a></li><li><a href="../reports/complete_findings.rna_validation_failures.tsv">Full RNA failure audit</a></li></ul>'
    (out/'13_pgtk_detailed_evidence_mqc.html').write_text(links,encoding='utf-8')
if __name__=='__main__':main()
