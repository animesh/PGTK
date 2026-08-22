#!/usr/bin/env python3
import argparse, csv, html, json
from pathlib import Path

def rows(path):
    p=Path(path)
    if not p.exists() or p.stat().st_size==0: return []
    with p.open(encoding='utf-8', errors='replace', newline='') as h:
        return list(csv.DictReader(h, delimiter='\t'))

def num(v):
    try: return float(str(v).replace(',',''))
    except: return 0.0

def write_json(out, name, payload):
    (out/name).write_text(json.dumps(payload, indent=2), encoding='utf-8')

def plot(out, ident, title, description, data, pconfig, headers=None):
    payload={'id':ident,'section_name':title,'description':description,'plot_type':pconfig.pop('type'),'pconfig':pconfig,'data':data}
    if headers: payload['headers']=headers
    write_json(out, f'{ident}_mqc.json', payload)

def placeholder(out, ident, title, text, status):
    body=f'<div class="alert alert-info"><strong>Status: {html.escape(status)}</strong><br>{html.escape(text)}</div>'
    (out/f'{ident}_mqc.html').write_text(body, encoding='utf-8')

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--output-dir', required=True)
    ap.add_argument('--samples', required=True)
    ap.add_argument('--complete-report', required=True)
    ap.add_argument('--rna-failure-report', required=True)
    ap.add_argument('--rna-variant-explanations', required=True)
    ap.add_argument('--comparative-report', required=True)
    ap.add_argument('--progression-report', required=True)
    ap.add_argument('--variant-inventory', required=True)
    ap.add_argument('--fasta-inventory', required=True)
    ap.add_argument('--rna-inventory', required=True)
    ap.add_argument('--external-comparison', required=True)
    ap.add_argument('--summary', required=True)
    ap.add_argument('--progression-summary', required=True)
    ap.add_argument('--progression-enrichment', required=True)
    ap.add_argument('--progression-pairwise', required=True)
    ap.add_argument('--maxquant-enabled', choices=['true','false'], default='false')
    args=ap.parse_args(); out=Path(args.output_dir); out.mkdir(parents=True,exist_ok=True)

    vi=rows(args.variant_inventory)
    # Flexible detection of sample/stage/count columns.
    data={}
    for r in vi:
        sample=r.get('Sample') or r.get('sample') or r.get('Name') or 'unknown'
        stage=r.get('Stage') or r.get('stage') or r.get('Category') or 'value'
        count=next((num(r[k]) for k in r if k.lower() in {'records','count','variants','alleles'}),0)
        data.setdefault(sample,{})[stage]=count
    if data:
        plot(out,'pgtk_variant_stages','Variant stage attrition','Counts retained across calling, filtering, annotation and RNA-validation stages.',data,{'type':'bargraph','id':'pgtk_variant_stages','title':'Variants by sample and stage','ylab':'Records'})

    ri=rows(args.rna_inventory); rdata={}
    for r in ri:
        sample=r.get('Sample') or r.get('sample') or 'all'
        category=r.get('Category') or r.get('EventType') or r.get('Type') or 'events'
        count=next((num(r[k]) for k in r if k.lower() in {'records','count','events','findings'}),0)
        rdata.setdefault(sample,{})[category]=count
    if rdata:
        plot(out,'pgtk_rna_events','RNA evidence inventory','RNA variants, progression variants, fusions and splice-junction evidence by sample.',rdata,{'type':'bargraph','id':'pgtk_rna_events','title':'RNA evidence classes','ylab':'Events'})

    ps=rows(args.progression_summary); pdata={}
    for r in ps:
        sample=r.get('Sample') or r.get('sample') or r.get('Set') or 'progression'
        for k,v in r.items():
            if k.lower() not in {'sample','subject','group','set','comparison'} and str(v).strip():
                try: pdata.setdefault(sample,{})[k]=float(v)
                except: pass
    if pdata:
        plot(out,'pgtk_progression_sets','Progression variant sets','Shared and sample-exclusive progression evidence. RNA absence at baseline does not establish DNA-level absence.',pdata,{'type':'bargraph','id':'pgtk_progression_sets','title':'Progression set sizes','ylab':'Alleles'})

    go=rows(args.progression_enrichment)
    sig=[]
    for r in go:
        fdr=num(r.get('FDR',1)); overlap=num(r.get('ProgressionGenesInTerm') or r.get('Overlap') or 0)
        if fdr<=0.1 and overlap>0: sig.append((fdr,r,overlap))
    sig.sort(key=lambda x:(x[0],-x[2]))
    table={}
    for _,r,overlap in sig[:30]:
        key=f"{r.get('Sample','')} | {r.get('GO_ID','')} | {r.get('GO_Name','')}"
        table[key]={'FDR':num(r.get('FDR')),'OddsRatio':num(r.get('OddsRatio')),'Genes':overlap,'Namespace':r.get('Namespace','')}
    if table:
        plot(out,'pgtk_progression_go','Progression GO enrichment','Top significant progression-variant GO terms, sorted by FDR. Full tables remain in results/progression_biology/.',table,{'type':'table','id':'pgtk_progression_go','title':'Significant progression GO terms'},headers={'FDR':{'title':'FDR','format':'{:.3g}'},'OddsRatio':{'title':'Odds ratio','format':'{:.3g}'},'Genes':{'title':'Overlap genes'}})
    else:
        placeholder(out,'pgtk_progression_go','Progression GO enrichment','No progression GO terms passed FDR <= 0.1 with a non-zero overlap.','NO_SIGNIFICANT_TERMS')

    ext=rows(args.external_comparison)
    meaningful=[r for r in ext if any(str(v).strip() not in {'','0','NA','N/A','not_run','disabled'} for v in r.values())]
    if meaningful:
        table={str(i+1):r for i,r in enumerate(meaningful[:100])}
        plot(out,'pgtk_sarek_validation','Sarek / external caller validation','External-caller comparison populated from the resumed validation branch.',table,{'type':'table','id':'pgtk_sarek_validation','title':'External caller comparison'})
    else:
        placeholder(out,'pgtk_sarek_validation','Sarek / external caller validation','Placeholder. Resume with --run_external_vcf_comparison true and configured external VCF inputs to populate overlap, concordance and stage-specific plots.','NOT_RUN')

    if args.maxquant_enabled=='true':
        placeholder(out,'pgtk_maxquant_validation','MaxQuant proteogenomic validation','MaxQuant integration was enabled. Detailed peptide, event, direct-MS/MS and MBR sections are populated by the final evidence processes.','ENABLED')
    else:
        placeholder(out,'pgtk_maxquant_validation','MaxQuant proteogenomic validation','Placeholder. Resume with --run_proteogenomic_validation true and --maxquant_txt to populate peptide mapping, direct MS/MS, MBR, sample-matched variant and junction evidence.','NOT_RUN')

    # Stable navigation and interpretation section.
    links='''<h3>Detailed results</h3><ul>
<li><a href="../qc/">QC and per-process reports</a></li><li><a href="../vcf_raw/">Raw VCFs</a></li>
<li><a href="../vcf_pass/">PASS VCFs</a></li><li><a href="../rna_validation/">RNA validation</a></li>
<li><a href="../progression_vcf/">Progression subtraction</a></li><li><a href="../progression_biology/">Progression biology and GO</a></li>
<li><a href="../expression/">Expression and GO</a></li><li><a href="../igv/findings/finding_explorer/">Finding explorer</a></li>
<li><a href="../proteogenomics_validation/">MaxQuant validation outputs</a></li><li><a href="../comparison/external_vcf/">Sarek / external-caller comparison</a></li></ul>
<p>Links are relative to the published results bundle. A missing target means the optional branch was not run.</p>'''
    (out/'pgtk_result_navigation_mqc.html').write_text(links,encoding='utf-8')
    reports=[('Complete findings',args.complete_report,'../reports/complete_findings.report.md'),('RNA validation failures',args.rna_failure_report,'../reports/complete_findings.rna_validation_failures.md'),('RNA variant explanations',args.rna_variant_explanations,'../reports/complete_findings.rna_variant_validation_explanations.md'),('Comparative biological advantage',args.comparative_report,'../comparative_advantage/comparative_advantage.report.md'),('Progression biology',args.progression_report,'../progression_biology/progression_biology.report.md')]
    catalogue='<h2>Original report catalogue</h2><ul>'+''.join(f'<li><a href="{link}">{html.escape(title)}</a></li>' for title,source,link in reports)+'</ul><p><a href="../expression/go/">Expression and GO</a>; <a href="../igv/findings/finding_explorer/index.html">IGV finding explorer</a>; <a href="../comparison/external_vcf/">External caller comparison</a>; <a href="../proteogenomics_validation/">MaxQuant validation</a>.</p>'
    (out/'00_pgtk_report_catalogue_mqc.html').write_text(catalogue,encoding='utf-8')
    for index,(title,source,link) in enumerate(reports,1):
        text=Path(source).read_text(encoding='utf-8',errors='replace')
        (out/f'{index:02d}_core_report_mqc.html').write_text(f'<h2>{html.escape(title)}</h2><p><a href="{link}">Open original report</a></p><pre style="white-space:pre-wrap">{html.escape(text)}</pre>',encoding='utf-8')

if __name__=='__main__': main()
