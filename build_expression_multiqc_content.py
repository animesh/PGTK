#!/usr/bin/env python3
import argparse, csv, html, json, math
from collections import defaultdict
from pathlib import Path
from report_legend import HTML_LEGEND

def read_rows(path):
    with Path(path).open(encoding='utf-8', errors='replace', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))
def numeric(value, default=0.0):
    try:
        parsed=float(value); return parsed if math.isfinite(parsed) else default
    except (TypeError, ValueError): return default
def analysis(row):
    return row.get('Analysis') or row.get('Sample') or row.get('Comparison') or 'analysis'
def top_terms(rows, limit=10):
    grouped=defaultdict(list)
    for row in rows:
        fdr=numeric(row.get('FDR'),1.0)
        size=numeric(row.get('ProgressionGenesInTerm') or row.get('ForegroundGenesInTerm') or row.get('Overlap') or row.get('GenesInTerm'),0)
        if fdr <= 0.1 and (size > 0 or 'ZScore' in row or 'NES' in row): grouped[analysis(row)].append((fdr,-size,row))
    selected=[]
    for name in sorted(grouped): selected.extend(item[2] for item in sorted(grouped[name],key=lambda x:(x[0],x[1],x[2].get('GO_ID','')))[:limit])
    return selected
def html_table(rows):
    if not rows: return '<p>No rows.</p>'
    fields=list(rows[0]); head=''.join(f'<th>{html.escape(x)}</th>' for x in fields)
    body=''.join('<tr>'+''.join(f'<td>{html.escape(str(row.get(x,"")))}</td>' for x in fields)+'</tr>' for row in rows)
    return f'<table class="table table-sm"><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>'
def main():
    p=argparse.ArgumentParser(); p.add_argument('--output-dir',required=True); p.add_argument('--expression-ora',required=True); p.add_argument('--ranked-go',required=True); p.add_argument('--expression-summary',required=True); p.add_argument('--variant-set-go',required=True); p.add_argument('--variant-set-summary',required=True); a=p.parse_args(); out=Path(a.output_dir); out.mkdir(parents=True,exist_ok=True)
    overview='---\nid: pgtk_expression_go_overview\nsection_name: Expression and GO overview\n---\n<h2>Expression and Gene Ontology</h2><p>Complete source tables remain under <a href="../expression/go/">expression/go</a> and <a href="../progression_biology/sets/">progression_biology/sets</a>.</p>'+html_table(read_rows(a.expression_summary))+html_table(read_rows(a.variant_set_summary))
    (out/'pgtk_expression_go_overview_mqc.html').write_text(overview,encoding='utf-8')
    for ident,title,path in [('pgtk_expression_ora_top','Expression GO over-representation',a.expression_ora),('pgtk_expression_ranked_go_top','Ranked expression GO',a.ranked_go),('pgtk_progression_variant_set_go_top','Progression variant-set GO',a.variant_set_go)]:
        rows=top_terms(read_rows(path)); data={}
        for index,row in enumerate(rows,1):
            key=f"{analysis(row)} | {row.get('GO_ID','')} | {row.get('GO_Name','')} | {index}"
            data[key]={'Analysis':analysis(row),'GO ID':row.get('GO_ID',''),'GO term':row.get('GO_Name',''),'Namespace':row.get('Namespace',''),'FDR':numeric(row.get('FDR'),1),'P value':numeric(row.get('PValue'),1),'Effect':numeric(row.get('OddsRatio') or row.get('ZScore') or row.get('NES'),0),'Genes':numeric(row.get('ProgressionGenesInTerm') or row.get('ForegroundGenesInTerm') or row.get('Overlap') or row.get('GenesInTerm'),0)}
        if data:
            payload={'id':ident,'section_name':title,'description':'Top significant terms per analysis. Full tables are linked from the report catalogue.','plot_type':'table','pconfig':{'id':ident,'title':title},'data':data}; (out/f'{ident}_mqc.json').write_text(json.dumps(payload,indent=2),encoding='utf-8')
        else: (out/f'{ident}_mqc.html').write_text(f'<p>{title}: no terms passed FDR &le; 0.1.</p>',encoding='utf-8')
if __name__=='__main__': main()
