#!/usr/bin/env python3
import argparse
import csv
import gzip
import json
import os
import re
from collections import Counter, defaultdict
from pathlib import Path
from report_legend import HTML_LEGEND

FIELDS = (
    'EventID','EvidenceClasses','SourceEvents','Sample','Gene','PredictedConsequence',
    'PredictedImpact','Transcript','ProteinChange','Chrom','Position','REF','ALT',
    'ReadValidationStatus','ValidationExplanation','CountUnit','UniqueAlignments',
    'CallableAlignments','ExactAltReads','CleanReferenceReads','ExcludedReads',
    'AltFractionAmongClean','CallableFractionAmongExamined'
)
INTS = {'Position','UniqueAlignments','CallableAlignments','ExactAltReads','CleanReferenceReads','ExcludedReads'}

def safe(value):
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', value or 'none').strip('_') or 'none'

def allele_type(ref, alt):
    if not ref or not alt:
        return ''
    if len(ref) == len(alt) == 1:
        return 'SNV'
    if len(ref) == len(alt):
        return 'MNV'
    if len(alt) > len(ref) and alt.startswith(ref):
        return 'INSERTION'
    if len(ref) > len(alt) and ref.startswith(alt):
        return 'DELETION'
    return 'COMPLEX_ALLELE'

def load_source_events(path):
    events = {}
    with Path(path).open(newline='', encoding='utf-8', errors='replace') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            events[row['Event']] = row
    return events

def geometry_for(row, source_events):
    classes = set(filter(None, row['EvidenceClasses'].split(';')))
    sources = [source_events[event] for event in row.get('SourceEvents', '').split(';') if event in source_events]
    ref, alt = row['REF'], row['ALT']
    variant = allele_type(ref, alt)
    if 'fusion' in classes:
        event_type = 'FUSION'
    elif 'splice_junction' in classes:
        event_type = 'SPLICE_JUNCTION'
    elif variant:
        event_type = variant
    else:
        event_type = 'CONTEXT_EVENT'
    regions = []
    seen = set()
    if event_type in {'SNV','MNV','INSERTION','DELETION','COMPLEX_ALLELE'}:
        start0 = max(0, int(row['Position']) - 1)
        end0 = start0 + max(1, len(ref))
        regions.append({'chrom':row['Chrom'],'start0':start0,'end0':end0,'role':'TARGET_ALLELE'})
    else:
        for source in sources:
            candidates = [(source.get('Chrom',''), source.get('Start0',''), source.get('End',''), 'JUNCTION' if event_type == 'SPLICE_JUNCTION' else 'BREAKPOINT_1')]
            if source.get('Chrom2'):
                candidates.append((source['Chrom2'], source.get('Start2_0',''), source.get('End2',''), 'BREAKPOINT_2'))
            for chrom, start, end, role in candidates:
                if not chrom or start == '' or end == '':
                    continue
                key = (chrom, int(start), int(end), role)
                if key not in seen:
                    seen.add(key)
                    regions.append({'chrom':chrom,'start0':int(start),'end0':int(end),'role':role})
    return {'event_type':event_type,'regions':regions,'source_events':[s.get('Event','') for s in sources]}

SERVER = r'''#!/usr/bin/env python3
import csv,gzip,hashlib,html,json,os,re,shutil,subprocess,sys,urllib.parse
from http.server import ThreadingHTTPServer,SimpleHTTPRequestHandler
from pathlib import Path
from report_legend import HTML_LEGEND
R=Path(__file__).resolve().parent
C=json.loads((R/'explorer_config.json').read_text())
G=json.loads((R/C['geometry_file']).read_text())
K=R/'report_cache';K.mkdir(exist_ok=True)
REPORT_CACHE_VERSION='6-persistent-event-panel'
A=[]
def resolve(value):
 p=Path(value).expanduser();return p.resolve() if p.is_absolute() else (R/p).resolve()
def required(value,label):
 p=Path(value)
 if not p.is_file() or p.stat().st_size==0:raise RuntimeError(f'Missing {label}: {p}')
 return p
def binds(paths):
 roots=[]
 for value in paths:
  path=Path(value).resolve();path=path.parent if path.is_file() else path;text=str(path)
  if text not in roots:roots.append(text)
 result=[]
 for root in roots:result.extend(['--bind',root+':'+root])
 return result
def run(command,label):
 result=subprocess.run(command,cwd=R,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
 if result.returncode:raise RuntimeError(f"{label} failed: {(result.stderr or result.stdout).strip()}")
 return result
C['genome']=str(Path(os.environ.get('PGTK_IGV_GENOME',C['genome'])).expanduser().resolve())
C['display_manifest']=str(resolve(C['display_manifest']))
for track in C['tracks']:track['resolved_path']=str(resolve(track['path']))
required(C['genome'],'genome');required(C['display_manifest'],'display manifest')
for track in C['tracks']:
 path=required(track['resolved_path'],'track')
 index=Path(str(path)+'.bai') if Path(str(path)+'.bai').is_file() else path.with_suffix('.bai')
 required(index,'BAM index')
for part in C['partitions']:
 with gzip.open(R/part['File'],'rt') as handle:A.extend(json.loads(line) for line in handle if line.strip())
I={row['EventID']:row for row in A}
def search(query):
 word=query.get('q',[''])[0].lower();sample=query.get('sample',[''])[0];klass=query.get('class',[''])[0];impact=query.get('impact',[''])[0];chrom=query.get('chrom',[''])[0];event_type=query.get('event_type',[''])[0];rows=[]
 for row in A:
  if sample and row['Sample']!=sample:continue
  if klass and klass not in row['EvidenceClasses']:continue
  if impact and row['PredictedImpact']!=impact:continue
  if chrom and row['Chrom']!=chrom:continue
  if event_type and row.get('EventType')!=event_type:continue
  if word and word not in ' '.join((row['EventID'],row['Gene'],row['Transcript'],row['PredictedConsequence'],row['ProteinChange'])).lower():continue
  rows.append(row)
 offset=max(0,int(query.get('offset',['0'])[0]));limit=min(500,int(query.get('limit',['100'])[0]));return {'total':len(rows),'rows':rows[offset:offset+limit]}
def visual_status(row,geometry):
 if geometry['event_type'] in {'FUSION','SPLICE_JUNCTION','CONTEXT_EVENT'}:
  return 'CONTEXT_ALIGNMENTS_AVAILABLE' if row.get('ContextAlignments',0) else 'NO_CONTEXT_ALIGNMENTS'
 return row['ReadValidationStatus']
def coordinate(chrom,position):
 return f"chr{str(chrom).removeprefix('chr')}:{int(position):,}"
def event_description(row,geometry):
 event_type=geometry['event_type'];ref=row.get('REF','');alt=row.get('ALT','');chrom=row.get('Chrom','');position=int(row.get('Position') or 0)
 regions=geometry.get('regions',[])
 if event_type=='SNV':
  action=f"SUBSTITUTE {ref} WITH {alt} AT {coordinate(chrom,position)}";reference=ref;alternate=alt
 elif event_type=='MNV':
  end=position+len(ref)-1;action=f"REPLACE {ref} WITH {alt} AT {coordinate(chrom,position)}-{end:,}";reference=ref;alternate=alt
 elif event_type=='INSERTION' and alt.startswith(ref):
  inserted=alt[len(ref):];anchor_end=position+len(ref)-1;action=f"INSERT {inserted} AFTER {coordinate(chrom,anchor_end)} {ref}";reference=ref;alternate=ref+'['+inserted+']'
 elif event_type=='DELETION' and ref.startswith(alt):
  deleted=ref[len(alt):];anchor_end=position+len(alt)-1;action=f"DELETE {deleted} AFTER {coordinate(chrom,anchor_end)} {alt}";reference=alt+'['+deleted+']';alternate=alt
 elif event_type=='COMPLEX_ALLELE':
  end=position+len(ref)-1;action=f"REPLACE {ref} WITH {alt} AT {coordinate(chrom,position)}-{end:,}";reference=ref;alternate=alt
 elif event_type=='SPLICE_JUNCTION' and regions:
  region=regions[0];action=f"SPLICE {coordinate(region['chrom'],region['start0']+1)} TO {coordinate(region['chrom'],region['end0'])}";reference='genomic interval retained in reference';alternate='RNA junction joins the two boundaries'
 elif event_type=='FUSION' and len(regions)>=2:
  first,second=regions[0],regions[1];action=f"FUSE {coordinate(first['chrom'],first['start0']+1)} TO {coordinate(second['chrom'],second['start0']+1)}";reference='separate genomic loci';alternate='RNA fusion connects both breakpoints'
 elif regions:
  region=regions[0];action=f"CONTEXT EVENT AT {coordinate(region['chrom'],region['start0']+1)}-{region['end0']:,}";reference=ref or 'not applicable';alternate=alt or 'not applicable'
 else:
  action=f"{event_type} {row['EventID']}";reference=ref or 'not applicable';alternate=alt or 'not applicable'
 return {'action':action,'reference':reference,'alternate':alternate}
def event_summary_html(row,geometry,description):
 exact=int(row.get('ExactAltReads',0));clean=int(row.get('CleanReferenceReads',0));excluded=int(row.get('ExcludedReads',0));callable_count=int(row.get('CallableReads',exact+clean));fraction=row.get('AltFractionAmongClean');fraction_text='NA' if fraction is None else f"{float(fraction):.6f}"
 structural=geometry['event_type'] in {'FUSION','SPLICE_JUNCTION','CONTEXT_EVENT'}
 evidence_label='Context alignments' if structural else 'Exact ALT alignments'
 evidence_value=int(row.get('ContextAlignments',0)) if structural else exact
 note='Structural-event alignments provide context and are not allele-classified support.' if structural else 'Only the exact change at the red PGTK_TARGET marker defines ALT support. Other colored bases are unrelated alignment differences.'
 if geometry['event_type']=='INSERTION':
  following=int(row.get('Position') or 0)+len(row.get('REF',''))
  note+=f" The bracketed inserted sequence has no reference coordinate. The following reference base begins at {coordinate(row.get('Chrom',''),following)}. Purple I symbols at the target boundary in exact-alt reads represent the insertion."
 note+=' Counts are from the complete classification; displayed read tracks may be capped.'
 return '<section id="pgtk-event-summary" style="position:sticky;top:0;z-index:2147483647;font:14px Arial;margin:8px;padding:14px;border:3px solid #b00020;background:#fff7f7;box-shadow:0 2px 8px rgba(0,0,0,.22)">'+f'<h2 style="margin:0 0 8px;color:#8b0018">Selected event: {html.escape(description["action"])}</h2>'+f'<p><b>Event type:</b> {html.escape(geometry["event_type"])} &nbsp; <b>Gene:</b> {html.escape(row.get("Gene",""))} &nbsp; <b>Sample:</b> {html.escape(row.get("Sample",""))}</p>'+f'<p style="font-family:monospace;font-size:16px"><b>REFERENCE:</b> {html.escape(description["reference"])}<br><b>ALTERNATE:</b> {html.escape(description["alternate"])}</p>'+f'<p><b>{evidence_label}:</b> {evidence_value} &nbsp; <b>Clean reference alignments:</b> {clean} &nbsp; <b>Excluded:</b> {excluded} &nbsp; <b>Callable:</b> {callable_count} &nbsp; <b>ALT fraction:</b> {fraction_text} &nbsp; <b>Status:</b> {html.escape(row.get("VisualEvidenceStatus",row.get("ReadValidationStatus","")))}</p>'+f'<p><b>How to read this report:</b> {html.escape(note)}</p></section>'
def persistent_summary_script(summary):
 payload=json.dumps(summary).replace('</','<\\/')
 return '<script id="pgtk-event-summary-loader">(function(){const markup='+payload+';function ensure(){let panel=document.getElementById("pgtk-event-summary");if(panel)return;const template=document.createElement("template");template.innerHTML=markup.trim();panel=template.content.firstElementChild;const container=document.getElementById("container");if(container)container.insertBefore(panel,container.firstChild);else document.body.insertBefore(panel,document.body.firstChild);}if(document.readyState==="loading")document.addEventListener("DOMContentLoaded",ensure,{once:true});else ensure();new MutationObserver(ensure).observe(document.documentElement,{childList:true,subtree:true});})();</script>'
def report(event_id):
 if event_id not in I:raise KeyError(f'Unknown EventID: {event_id}')
 row=I[event_id];geometry=G[event_id];description=event_description(row,geometry);sample=row['Sample'];safe=re.sub(r'[^A-Za-z0-9_.-]+','_',event_id).strip('._') or 'event';safe+='-'+hashlib.sha256(event_id.encode()).hexdigest()[:12]
 output=K/(safe+'.html');event_dir=K/safe;signature_path=event_dir/'cache.signature'
 tracks={track['category']:track['resolved_path'] for track in C['tracks'] if track['sample']==sample}
 signature_data={'version':REPORT_CACHE_VERSION,'event':row,'geometry':geometry,'description':description,'genome':C['genome'],'manifest':C['display_manifest'],'tracks':{k:[v,Path(v).stat().st_size,Path(v).stat().st_mtime_ns] for k,v in tracks.items()}}
 signature=hashlib.sha256(json.dumps(signature_data,sort_keys=True,separators=(',',':')).encode()).hexdigest()
 if output.exists() and signature_path.is_file() and signature_path.read_text().strip()==signature:return output
 if event_dir.exists():shutil.rmtree(event_dir)
 event_dir.mkdir(parents=True);site=event_dir/'finding.tsv';marker=event_dir/'event_coordinate.bed'
 columns=['Chrom','Start','End','EventID','RegionRole','Gene','CandidateClass','EventType','REF','ALT','PredictedImpact','PredictedConsequence','ReadValidationStatus','VisualEvidenceStatus','ExactAltAlignments','CleanReferenceAlignments','ContextAlignments','ExcludedAlignments','TotalAlignmentsExamined','CallableAlignments','ALT_Fraction_Among_Callable','Callable_Fraction_Among_Examined','Interpretation']
 with site.open('w',newline='') as handle:
  writer=csv.writer(handle,delimiter='\t',lineterminator='\n');writer.writerow(columns)
  for region in geometry['regions']:
   writer.writerow([region['chrom'],region['start0']+1,region['end0'],row['EventID'],region['role'],row['Gene'],row['EvidenceClasses'],geometry['event_type'],row['REF'],row['ALT'],row['PredictedImpact'],row['PredictedConsequence'],row['ReadValidationStatus'],visual_status(row,geometry),row['ExactAltReads'],row['CleanReferenceReads'],row.get('ContextAlignments',0),row['ExcludedReads'],row['TotalReadsExamined'],row['CallableReads'],row['AltFractionAmongClean'] if row['AltFractionAmongClean'] is not None else 'NA',row['CallableFractionAmongExamined'] if row['CallableFractionAmongExamined'] is not None else 'NA',row['ValidationExplanation']])
 with marker.open('w') as handle:
  handle.write('track name="PGTK_TARGET" description="Exact event geometry" itemRgb="On"\n')
  for region in geometry['regions']:
   label=f"{description['action']} | {region['role']} | {event_id}"
   handle.write(f"{region['chrom']}\t{region['start0']}\t{region['end0']}\t{label}\t1000\t.\t{region['start0']}\t{region['end0']}\t220,0,0\n")
 apptainer=os.environ.get('PGTK_APPTAINER','apptainer');pysam_image=required(os.environ.get('PGTK_PYSAM_IMAGE',''),'Pysam image');igv_image=required(os.environ.get('PGTK_IGV_REPORTS_IMAGE',''),'igv-reports image');builder=required(R/'prepare_event_igv_tracks.py','event track builder')
 command=[apptainer,'exec','--cleanenv','--no-home',*binds([R,C['display_manifest'],*tracks.values()]),str(pysam_image),'python3',str(builder),'--event-id',event_id,'--sample',sample,'--display-manifest',C['display_manifest'],'--output-dir',str(event_dir)]
 for category,path in tracks.items():command.extend(['--bam',f'{category}={path}'])
 run(command,'event BAM extraction')
 with required(event_dir/'tracks.tsv','event track manifest').open(newline='') as handle:event_tracks=[entry['BAM'] for entry in csv.DictReader(handle,delimiter='\t')]
 command=[apptainer,'exec','--cleanenv','--no-home',*binds([R,event_dir,C['genome']]),str(igv_image),'create_report',str(site),'--fasta',C['genome'],'--sequence','1','--begin','2','--end','3','--tracks',str(marker),*event_tracks,'--flanking',str(C['flanking']),'--title','PGTK | '+description['action']+' | '+sample,'--standalone','--output',str(output)]
 run(command,'IGV report generation');required(output,'IGV report');signature_path.write_text(signature+'\n')
 text=output.read_text(encoding='utf-8',errors='replace');summary=event_summary_html(row,geometry,description);loader=persistent_summary_script(summary);insertion=summary+HTML_LEGEND+loader;rendered,count=re.subn(r'(<div\s+id=["\']container["\'][^>]*>)',lambda match:match.group(1)+insertion,text,count=1,flags=re.I);
 if not count:rendered,count=re.subn(r'(</body\s*>)',lambda match:insertion+match.group(1),text,count=1,flags=re.I);
 output.write_text(rendered if count else insertion+text,encoding='utf-8');return output
class H(SimpleHTTPRequestHandler):
 def do_GET(self):
  parsed=urllib.parse.urlparse(self.path)
  if parsed.path=='/api/findings':return self.json_response(search(urllib.parse.parse_qs(parsed.query)))
  if parsed.path=='/api/facets':return self.json_response(C['facets'])
  if parsed.path.startswith('/report/'):
   try:self.path='/report_cache/'+report(urllib.parse.unquote(parsed.path[8:])).name
   except Exception as error:return self.send_error(500,str(error))
  elif parsed.path=='/':self.path='/index.html'
  return super().do_GET()
 def json_response(self,value):
  body=json.dumps(value).encode();self.send_response(200);self.send_header('Content-Type','application/json');self.send_header('Content-Length',str(len(body)));self.end_headers();self.wfile.write(body)
def main():
 os.chdir(R);port=int(sys.argv[1]) if len(sys.argv)>1 else 8765;print(f'Loaded {len(A)} findings at http://127.0.0.1:{port}',flush=True);ThreadingHTTPServer(('127.0.0.1',port),H).serve_forever()
if __name__=='__main__':main()
'''

HTML_HEAD = '<!doctype html><html><meta charset="utf-8"><title>PGTK Offline Variant Explorer</title><style>body{font:14px Arial;margin:18px;color:#18202a}input,select,button{padding:7px;margin:3px}table{border-collapse:collapse;width:100%;font-size:12px}th,td{padding:5px;border-bottom:1px solid #ddd;text-align:left}.note{background:#eef6ff;padding:12px;border-left:4px solid #2878c8}</style><h1>PGTK Offline Variant Explorer</h1>'+HTML_LEGEND+'<div class="note"><b>Coordinate display:</b> reports include a red TARGET track. Allele events use normalized REF/ALT geometry; splice junctions use donor-to-acceptor intervals; fusions use both breakpoints. Structural-event event_display tracks are context alignments, not allele-classified support. Colored bases elsewhere are incidental read differences. <b>Direct-open mode:</b> counts work through file:///; server mode enables reports. <span id="mode"></span></div><p><input id="q" placeholder="Gene, event, consequence"><select id="status"><option value="ALT_SUPPORTED">ALT-supported only</option><option value="">All candidates</option></select><select id="sample"><option value="">All samples</option></select><select id="class"><option value="">All classes</option></select><select id="impact"><option value="">All impacts</option></select><select id="event_type"><option value="">All event types</option></select><input id="chrom" placeholder="Chromosome"><button id="prev">Previous</button><button id="next">Next</button> <span id="page"></span></p><table><thead><tr><th>Sample</th><th>Class</th><th>Event type</th><th>Gene</th><th>Impact</th><th>Consequence</th><th>Locus</th><th>ALT</th><th>REF</th><th>Context</th><th>Excluded</th><th>Callable</th><th>ALT fraction</th><th>Evidence status</th><th>Interpretation</th><th>Alignment view</th></tr></thead><tbody id="rows"></tbody></table><script>const D='
HTML_TAIL = r''';const $=x=>document.getElementById(x),E=x=>String(x??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));let Z=[],o=0,l=100;function init(){$('mode').innerHTML=location.protocol.startsWith('http')?'<b>Server mode:</b> use Open alignments.':'<b>Direct-file mode:</b> start serve_explorer.sh for reports.';for(const [id,n] of [['sample',0],['class',1],['impact',4],['event_type',2]])for(const v of [...new Set(D.map(r=>r[n]))].sort()){const x=document.createElement('option');x.value=x.textContent=v;$(id).appendChild(x)}for(const id of ['status','sample','class','impact','event_type','chrom'])$(id).onchange=()=>{o=0;load()};$('q').oninput=()=>{o=0;load()};$('prev').onclick=()=>{o=Math.max(0,o-l);draw()};$('next').onclick=()=>{o=Math.min(Math.max(0,Z.length-1),o+l);draw()};load()}function load(){const q=$('q').value.toLowerCase(),st=$('status').value,s=$('sample').value,c=$('class').value,i=$('impact').value,et=$('event_type').value,ch=$('chrom').value;Z=D.filter(r=>(!st||r[18]===st)&&(!s||r[0]===s)&&(!c||r[1].includes(c))&&(!i||r[4]===i)&&(!et||r[2]===et)&&(!ch||r[7]===ch)&&(!q||[r[3],r[5],r[11]].join(' ').toLowerCase().includes(q)));draw()}function draw(){$('page').textContent=`${Z.length?o+1:0}-${Math.min(o+l,Z.length)} of ${Z.length}`;$('rows').innerHTML=Z.slice(o,o+l).map(r=>`<tr><td>${E(r[0])}</td><td>${E(r[1])}</td><td>${E(r[2])}</td><td>${E(r[3])}</td><td>${E(r[4])}</td><td>${E(r[5])}</td><td>${E(r[6])}</td><td>${r[12]}</td><td>${r[13]}</td><td>${r[14]}</td><td>${r[15]}</td><td>${r[16]}</td><td>${r[17]===null?'NA':Number(r[17]).toFixed(4)}</td><td>${E(r[18])}</td><td>${E(r[19])}</td><td>${location.protocol.startsWith('http')?`<button onclick="window.open('/report/${encodeURIComponent(r[11])}','_blank')">Open alignments</button>`:'Start server mode'}</td></tr>`).join('')}init()</script></html>'''

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest', required=True)
    parser.add_argument('--events', required=True)
    parser.add_argument('--excluded-reads')
    parser.add_argument('--bam-manifest', required=True)
    parser.add_argument('--display-manifest', required=True)
    parser.add_argument('--genome', required=True)
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--flanking', type=int, default=150)
    parser.add_argument('--report-command', default='create_report')
    parser.add_argument('--max-html-mb', type=int, default=80)
    args = parser.parse_args()
    output = Path(args.output_dir); partitions = output/'partitions'; partitions.mkdir(parents=True,exist_ok=True)
    reasons=defaultdict(Counter)
    if args.excluded_reads and Path(args.excluded_reads).is_file():
        with open(args.excluded_reads,newline='',encoding='utf-8',errors='replace') as handle:
            for row in csv.DictReader(handle,delimiter='\t'):
                if row.get('EventID'):reasons[row['EventID']][row.get('Reason') or 'UNSPECIFIED']+=1
    source_events=load_source_events(args.events)
    display_counts=Counter()
    with gzip.open(args.display_manifest,'rt',encoding='utf-8',newline='') as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            display_counts[(row['EventID'],row['Category'])]+=1
    findings=[]; geometries={}; facets={'samples':set(),'classes':set(),'impacts':set(),'statuses':set(),'event_types':set()}
    with open(args.manifest,newline='',encoding='utf-8',errors='replace') as handle:
        for raw in csv.DictReader(handle,delimiter='\t'):
            row={key:raw.get(key,'') for key in FIELDS}
            for key in INTS:
                try:row[key]=int(float(row[key] or 0))
                except ValueError:row[key]=0
            text=(row.get('AltFractionAmongClean') or '').strip();row['AltFractionAmongClean']=None if text.upper() in {'','NA','N/A','NULL','NONE','.'} else float(text)
            callable_count=row['ExactAltReads']+row['CleanReferenceReads'];total=row['UniqueAlignments']
            if row['CallableAlignments'] != callable_count or total != callable_count + row['ExcludedReads']:
                raise ValueError(f"Manifest count invariant failed for {row['EventID']}")
            expected=(row['ExactAltReads']/callable_count) if callable_count else None
            if expected is None and row['AltFractionAmongClean'] not in (None,0.0):raise ValueError(f"Manifest ALT fraction invariant failed for {row['EventID']}")
            if expected is not None and (row['AltFractionAmongClean'] is None or f"{row['AltFractionAmongClean']:.6f}" != f"{expected:.6f}"):raise ValueError(f"Manifest ALT fraction invariant failed for {row['EventID']}")
            geometry=geometry_for(row,source_events);geometries[row['EventID']]=geometry
            context=display_counts[(row['EventID'],'event_display')]
            visual_status=('CONTEXT_ALIGNMENTS_AVAILABLE' if context else 'NO_CONTEXT_ALIGNMENTS') if geometry['event_type'] in {'FUSION','SPLICE_JUNCTION','CONTEXT_EVENT'} else row['ReadValidationStatus']
            detail='; '.join(f'{key}={value}' for key,value in reasons[row['EventID']].most_common())
            explanation=row.get('ValidationExplanation') or 'No authoritative validation explanation was supplied.'
            if detail:explanation+=' Retained exclusion reasons: '+detail+'.'
            elif row['ExcludedReads']>0:explanation+=' Per-reason diagnostic rows were unavailable or capped.'
            row.update({'AltFractionAmongClean':expected,'CallableFractionAmongExamined':(callable_count/total) if total else None,'TotalReadsExamined':total,'CallableReads':callable_count,'ValidationStatus':row['ReadValidationStatus'],'VisualEvidenceStatus':visual_status,'EventType':geometry['event_type'],'ContextAlignments':context,'ValidationExplanation':explanation})
            findings.append(row);facets['samples'].add(row['Sample']);facets['classes'].update(filter(None,row['EvidenceClasses'].split(';')));facets['impacts'].add(row['PredictedImpact']);facets['statuses'].add(visual_status);facets['event_types'].add(geometry['event_type'])
    compact=[[r['Sample'],r['EvidenceClasses'],r['EventType'],r['Gene'],r['PredictedImpact'],r['PredictedConsequence'],'; '.join(f"{g['chrom']}:{g['start0']+1}-{g['end0']}" for g in geometries[r['EventID']]['regions']),r['Chrom'],r['Position'],r['REF'],r['ALT'],r['EventID'],r['ExactAltReads'],r['CleanReferenceReads'],r['ContextAlignments'],r['ExcludedReads'],r['CallableReads'],r['AltFractionAmongClean'],r['VisualEvidenceStatus'],r['ValidationExplanation']] for r in findings]
    html=HTML_HEAD+json.dumps(compact,separators=(',',':')).replace('</','<\\/')+HTML_TAIL;html_bytes=len(html.encode());limit=args.max_html_mb*1024*1024
    if html_bytes>limit:raise RuntimeError(f'Direct-open HTML exceeds safeguard: {html_bytes} bytes > {limit}')
    (output/'index.html').write_text(html,encoding='utf-8')
    with gzip.open(partitions/'all.jsonl.gz','wt',encoding='utf-8',compresslevel=6) as handle:
        for row in findings:handle.write(json.dumps(row,separators=(',',':'))+'\n')
    (output/'event_geometry.json').write_text(json.dumps(geometries,indent=2,sort_keys=True)+'\n')
    tracks=[]
    with open(args.bam_manifest,newline='',encoding='utf-8',errors='replace') as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            if row.get('Category') in ('event_display','exact_alt_display','reference_display'):tracks.append({'sample':row['Sample'],'category':row['Category'],'path':'../finding_reviews/'+row['BAM'],'alignments':int(row.get('UniqueAlignments') or 0)})
    config={'genome':str(Path(args.genome).resolve()),'display_manifest':'../finding_reviews/'+Path(args.display_manifest).name,'geometry_file':'event_geometry.json','tracks':tracks,'flanking':args.flanking,'report_command':args.report_command,'findings':len(findings),'partitions':[{'PartitionID':'all','File':'partitions/all.jsonl.gz','Findings':len(findings)}],'facets':{key:sorted(value) for key,value in facets.items()}}
    (output/'explorer_config.json').write_text(json.dumps(config,indent=2)+'\n');(output/'server.py').write_text(SERVER);os.chmod(output/'server.py',0o755)
    (output/'README.txt').write_text('The red TARGET track is authoritative event geometry. Structural event_display tracks are context alignments, not allele-classified support.\n')
    (output/'coverage_summary.txt').write_text(f'Findings: {len(findings)}\nEmbedded compact records: {len(compact)}\nHTML bytes: {html_bytes}\nFindings discarded: 0\nDatabase files: 0\nDirect-open HTML: yes\n')
    print(f'Built compact direct-open explorer with {len(findings)} findings and {html_bytes} HTML bytes')
if __name__=='__main__':main()