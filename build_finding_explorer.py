#!/usr/bin/env python3
import argparse,csv,gzip,json,os,re
from collections import Counter,defaultdict
from pathlib import Path
from report_legend import HTML_LEGEND
FIELDS=('EventID','EvidenceClasses','Sample','Gene','PredictedConsequence','PredictedImpact','Transcript','ProteinChange','Chrom','Position','REF','ALT','ReadValidationStatus','ValidationExplanation','CountUnit','UniqueAlignments','CallableAlignments','ExactAltReads','CleanReferenceReads','ExcludedReads','AltFractionAmongClean','CallableFractionAmongExamined')
INTS={'Position','UniqueAlignments','CallableAlignments','ExactAltReads','CleanReferenceReads','ExcludedReads'}
def safe(x):return re.sub(r'[^A-Za-z0-9_.-]+','_',x or 'none').strip('_') or 'none'
SERVER='''#!/usr/bin/env python3
import gzip,json,os,re,subprocess,sys,urllib.parse
from http.server import ThreadingHTTPServer,SimpleHTTPRequestHandler
from pathlib import Path
from report_legend import HTML_LEGEND
R=Path(__file__).resolve().parent;C=json.loads((R/'explorer_config.json').read_text());K=R/'report_cache';K.mkdir(exist_ok=True);A=[]
def resolve_resource(value):
 p=Path(value).expanduser();return p.resolve() if p.is_absolute() else (R/p).resolve()
C['genome']=str(Path(os.environ.get('PGTK_IGV_GENOME',C['genome'])).expanduser().resolve())
for t in C['tracks']:t['resolved_path']=str(resolve_resource(t['path']))
missing=[]
if not Path(C['genome']).is_file():missing.append('genome='+C['genome'])
for t in C['tracks']:
 p=Path(t['resolved_path'])
 if not p.is_file():missing.append('track='+str(p))
 elif p.suffix=='.bam' and not (Path(str(p)+'.bai').is_file() or p.with_suffix('.bai').is_file() or Path(str(p)+'.csi').is_file()):missing.append('index='+str(p))
if missing:raise SystemExit('Missing explorer resources:\\n'+'\\n'.join(missing))
for p in C['partitions']:
 with gzip.open(R/p['File'],'rt') as h:A.extend(json.loads(x) for x in h if x.strip())
I={x['EventID']:x for x in A}
def search(q):
 w=q.get('q',[''])[0].lower();s=q.get('sample',[''])[0];c=q.get('class',[''])[0];i=q.get('impact',[''])[0];ch=q.get('chrom',[''])[0];z=[]
 for r in A:
  if s and r['Sample']!=s:continue
  if c and c not in r['EvidenceClasses']:continue
  if i and r['PredictedImpact']!=i:continue
  if ch and r['Chrom']!=ch:continue
  if w and w not in ' '.join((r['EventID'],r['Gene'],r['Transcript'],r['PredictedConsequence'],r['ProteinChange'])).lower():continue
  z.append(r)
 o=max(0,int(q.get('offset',['0'])[0]));n=min(500,int(q.get('limit',['100'])[0]));return {'total':len(z),'rows':z[o:o+n]}
def report(e):
 out=K/(e+'.html')
 if out.exists():return out
 r=I[e];site=K/(e+'.tsv');site.write_text('Chrom\\tStart\\tEnd\\tEventID\\tGene\\tCandidateClass\\tPredictedImpact\\tPredictedConsequence\\tReadValidationStatus\\tExactAltAlignments\\tCleanReferenceAlignments\\tExcludedAlignments\\tTotalAlignmentsExamined\\tCallableAlignments\\tALT_Fraction_Among_Callable\\tCallable_Fraction_Among_Examined\\tInterpretation\\n'+f"{r['Chrom']}\\t{r['Position']}\\t{r['Position']}\\t{r['EventID']}\\t{r['Gene']}\\t{r['EvidenceClasses']}\\t{r['PredictedImpact']}\\t{r['PredictedConsequence']}\\t{r['ReadValidationStatus']}\\t{r['ExactAltReads']}\\t{r['CleanReferenceReads']}\\t{r['ExcludedReads']}\\t{r['TotalReadsExamined']}\\t{r['CallableReads']}\\t{r['AltFractionAmongClean'] if r['AltFractionAmongClean'] is not None else 'NA'}\\t{r['CallableFractionAmongExamined'] if r['CallableFractionAmongExamined'] is not None else 'NA'}\\t{r['ValidationExplanation']}\\n")
 sample=r['Sample']
 tracks=[t['resolved_path'] for t in C['tracks'] if t['sample']==sample and t['category'] in ('event_display','exact_alt_display','reference_display')]
 if not tracks:raise RuntimeError(f'No alignment tracks configured for sample {sample}')
 cmd=[C['report_command'],str(site),'--fasta',C['genome'],'--sequence','1','--begin','2','--end','3','--tracks',*tracks,'--flanking',str(C['flanking']),'--title','PGTK '+e+' | '+sample+' | all overlapping, exact ALT, clean REF','--standalone','--output',str(out)]
 result=subprocess.run(cmd,cwd=R,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
 if result.returncode!=0:
  out.unlink(missing_ok=True)
  detail=(result.stderr or result.stdout or 'create_report failed without diagnostic output').strip()
  raise RuntimeError(f"IGV report generation failed (exit {result.returncode}): {detail}")
 if not out.is_file() or out.stat().st_size==0:raise RuntimeError('IGV report generation returned success but produced no HTML output')
 text=out.read_text(encoding='utf-8',errors='replace'); rendered,n=re.subn(r'(<body\b[^>]*>)',lambda m:m.group(1)+HTML_LEGEND,text,count=1,flags=re.I); out.write_text(rendered if n else HTML_LEGEND+text,encoding='utf-8')
 return out
class H(SimpleHTTPRequestHandler):
 def do_GET(self):
  u=urllib.parse.urlparse(self.path)
  if u.path=='/api/findings':return self.j(search(urllib.parse.parse_qs(u.query)))
  if u.path=='/api/facets':return self.j(C['facets'])
  if u.path.startswith('/report/'):
   try:self.path='/report_cache/'+report(urllib.parse.unquote(u.path[8:])).name
   except Exception as x:return self.send_error(500,str(x))
  elif u.path=='/':self.path='/index.html'
  return super().do_GET()
 def j(self,x):
  b=json.dumps(x).encode();self.send_response(200);self.send_header('Content-Type','application/json');self.send_header('Content-Length',str(len(b)));self.end_headers();self.wfile.write(b)
os.chdir(R);port=int(sys.argv[1]) if len(sys.argv)>1 else 8765;print(f'Loaded {len(A)} findings at http://127.0.0.1:{port}',flush=True);ThreadingHTTPServer(('127.0.0.1',port),H).serve_forever()
'''
HTML_HEAD='<!doctype html><html><meta charset="utf-8"><title>PGTK Offline Variant Explorer</title><style>body{font:14px Arial;margin:18px;color:#18202a}input,select,button{padding:7px;margin:3px}table{border-collapse:collapse;width:100%;font-size:12px}th,td{padding:5px;border-bottom:1px solid #ddd;text-align:left}.note{background:#eef6ff;padding:12px;border-left:4px solid #2878c8}</style><h1>PGTK Offline Variant Explorer</h1>'+HTML_LEGEND+'<div class="note"><b>Direct-open mode:</b> this file works through file:/// without a server. Every row is an upstream RNA/progression candidate; strict read-validation status is shown separately. <span id="mode"></span></div><p><input id="q" placeholder="Gene, event, consequence"><select id="status"><option value="ALT_SUPPORTED">ALT-supported only</option><option value="">All candidates</option></select><select id="sample"><option value="">All samples</option></select><select id="class"><option value="">All classes</option></select><select id="impact"><option value="">All impacts</option></select><input id="chrom" placeholder="Chromosome"><button id="prev">Previous</button><button id="next">Next</button> <span id="page"></span></p><table><thead><tr><th>Sample</th><th>Class</th><th>Gene</th><th>Impact</th><th>Consequence</th><th>Locus</th><th>ALT-supporting</th><th>REF-supporting</th><th>Excluded/uncallable</th><th>Total examined</th><th>Callable</th><th>ALT fraction among callable</th><th>Callable fraction among examined</th><th>Validation status</th><th>Reason / interpretation</th><th>Alignment view</th></tr></thead><tbody id="rows"></tbody></table><script>const D='
HTML_TAIL=';const $=x=>document.getElementById(x),E=x=>String(x??\'\').replace(/[&<>"\']/g,c=>({\'&\':\'&amp;\',\'<\':\'&lt;\',\'>\':\'&gt;\',\'"\':\'&quot;\',"\'":\'&#39;\'}[c]));let Z=[],o=0,l=100;function init(){$(\'mode\').innerHTML=location.protocol.startsWith(\'http\')?\'<b>Server mode:</b> use Open alignments to generate IGV.js reports.\':\'<b>Direct-file mode:</b> start serve_explorer.sh to enable alignment reports.\';for(const [id,n] of [[\'sample\',0],[\'class\',1],[\'impact\',3]])for(const v of [...new Set(D.map(r=>r[n]))].sort()){const x=document.createElement(\'option\');x.value=x.textContent=v;$(id).appendChild(x)}for(const id of [\'status\',\'sample\',\'class\',\'impact\',\'chrom\'])$(id).onchange=()=>{o=0;load()};$(\'q\').oninput=()=>{o=0;load()};$(\'prev\').onclick=()=>{o=Math.max(0,o-l);draw()};$(\'next\').onclick=()=>{o=Math.min(Math.max(0,Z.length-1),o+l);draw()};load()}function load(){const q=$(\'q\').value.toLowerCase(),st=$(\'status\').value,s=$(\'sample\').value,c=$(\'class\').value,i=$(\'impact\').value,ch=$(\'chrom\').value;Z=D.filter(r=>(!st||r[17]===st)&&(!s||r[0]===s)&&(!c||r[1]===c)&&(!i||r[3]===i)&&(!ch||r[5]===ch)&&(!q||[r[2],r[4],r[9]].join(\' \').toLowerCase().includes(q)));draw()}function draw(){$(\'page\').textContent=`${Z.length?o+1:0}-${Math.min(o+l,Z.length)} of ${Z.length}`;$(\'rows\').innerHTML=Z.slice(o,o+l).map(r=>`<tr><td>${E(r[0])}</td><td>${E(r[1])}</td><td>${E(r[2])}</td><td>${E(r[3])}</td><td>${E(r[4])}</td><td>${E(r[5])}:${r[6]} ${E(r[7])}&gt;${E(r[8])}</td><td>${r[10]}</td><td>${r[11]}</td><td>${r[13]}</td><td>${r[14]}</td><td>${r[15]}</td><td>${r[12]===null?\'NA\':Number(r[12]).toFixed(4)}</td><td>${r[16]===null?\'NA\':Number(r[16]).toFixed(4)}</td><td>${E(r[17])}</td><td>${E(r[18])}</td><td>${location.protocol.startsWith(\'http\')?`<button onclick=\"window.open(\'/report/${encodeURIComponent(r[9])}\',\'_blank\')\">Open alignments</button>`:\'Start server mode\'}</td></tr>`).join(\'\')}init()</script></html>'
def main():
 a=argparse.ArgumentParser();a.add_argument('--manifest',required=True);a.add_argument('--excluded-reads');a.add_argument('--bam-manifest',required=True);a.add_argument('--genome',required=True);a.add_argument('--output-dir',required=True);a.add_argument('--flanking',type=int,default=150);a.add_argument('--report-command',default='create_report');a.add_argument('--max-html-mb',type=int,default=60);x=a.parse_args();out=Path(x.output_dir);pd=out/'partitions';pd.mkdir(parents=True,exist_ok=True);findings=[];fac={'samples':set(),'classes':set(),'impacts':set(),'statuses':set()};reasons=defaultdict(Counter)
 if x.excluded_reads and Path(x.excluded_reads).is_file():
  with open(x.excluded_reads,newline='',encoding='utf-8',errors='replace') as eh:
   for er in csv.DictReader(eh,delimiter='\t'):
    if er.get('EventID'):reasons[er['EventID']][er.get('Reason') or 'UNSPECIFIED']+=1
 with open(x.manifest,newline='',encoding='utf-8',errors='replace') as h:
  for r in csv.DictReader(h,delimiter='\t'):
   z={k:r.get(k,'') for k in FIELDS}
   for k in INTS:
    try:z[k]=int(float(z[k] or 0))
    except ValueError:z[k]=0
   raw_fraction=(z.get('AltFractionAmongClean') or '').strip()
   z['AltFractionAmongClean']=None if raw_fraction.upper() in {'','NA','N/A','NULL','NONE','.'} else float(raw_fraction)
   callable_reads=z['ExactAltReads']+z['CleanReferenceReads'];total=z['UniqueAlignments']
   if z['CallableAlignments'] != callable_reads or total != callable_reads + z['ExcludedReads']:
    raise ValueError(f"Manifest count invariant failed for {z['EventID']}")
   expected_fraction=(z['ExactAltReads']/callable_reads) if callable_reads else None
   if expected_fraction is None:
    if z['AltFractionAmongClean'] not in (None,0.0):
     raise ValueError(f"Manifest ALT fraction invariant failed for {z['EventID']}: observed={z['AltFractionAmongClean']}, expected=NA")
   elif z['AltFractionAmongClean'] is None or f"{z['AltFractionAmongClean']:.6f}" != f"{expected_fraction:.6f}":
    observed='NA' if z['AltFractionAmongClean'] is None else f"{z['AltFractionAmongClean']:.6f}"
    raise ValueError(f"Manifest ALT fraction invariant failed for {z['EventID']}: observed={observed}, expected={expected_fraction:.6f}")
   z['AltFractionAmongClean']=expected_fraction
   z['CallableFractionAmongExamined']=(callable_reads/total) if total else None
   detail='; '.join(f'{k}={v}' for k,v in reasons[z['EventID']].most_common())
   status=z.get('ReadValidationStatus') or 'NO_CALLABLE_READS';explanation=z.get('ValidationExplanation') or 'No authoritative validation explanation was supplied.'
   if detail:explanation+=' Retained exclusion reasons: '+detail+'.'
   elif z['ExcludedReads']>0:explanation+=' Per-reason diagnostic rows were unavailable or capped.'
   z['TotalReadsExamined']=total;z['CallableReads']=callable_reads;z['ValidationStatus']=status;z['ValidationExplanation']=explanation
   findings.append(z);fac['samples'].add(z['Sample']);fac['classes'].update(filter(None,z['EvidenceClasses'].split(';')));fac['impacts'].add(z['PredictedImpact']);fac['statuses'].add(status)
 compact=[[r['Sample'],r['EvidenceClasses'],r['Gene'],r['PredictedImpact'],r['PredictedConsequence'],r['Chrom'],r['Position'],r['REF'],r['ALT'],r['EventID'],r['ExactAltReads'],r['CleanReferenceReads'],r['AltFractionAmongClean'],r['ExcludedReads'],r['TotalReadsExamined'],r['CallableReads'],r['CallableFractionAmongExamined'],r['ValidationStatus'],r['ValidationExplanation']] for r in findings]
 html=HTML_HEAD+json.dumps(compact,separators=(',',':')).replace('</','<\\/')+HTML_TAIL;html_bytes=len(html.encode('utf-8'));max_bytes=x.max_html_mb*1024*1024
 if html_bytes>max_bytes:raise RuntimeError(f'Direct-open HTML exceeds safeguard: {html_bytes} bytes > {max_bytes}')
 (out/'index.html').write_text(html,encoding='utf-8')
 with gzip.open(pd/'all.jsonl.gz','wt',encoding='utf-8',compresslevel=6) as h:
  for r in findings:h.write(json.dumps(r,separators=(',',':'))+'\n')
 tracks=[]
 with open(x.bam_manifest,newline='',encoding='utf-8',errors='replace') as h:
  for r in csv.DictReader(h,delimiter='\t'):
   if r.get('Category') in ('event_display','exact_alt_display','reference_display'):tracks.append({'sample':r['Sample'],'category':r['Category'],'path':'../finding_reviews/'+r['BAM'],'alignments':int(r.get('UniqueAlignments') or 0)})
 parts=[{'PartitionID':'all','File':'partitions/all.jsonl.gz','Findings':len(findings)}];cfg={'genome':str(Path(x.genome).resolve()),'tracks':tracks,'flanking':x.flanking,'report_command':x.report_command,'findings':len(findings),'partitions':parts,'facets':{k:sorted(v) for k,v in fac.items()}}
 (out/'explorer_config.json').write_text(json.dumps(cfg,indent=2),encoding='utf-8');(out/'server.py').write_text(SERVER,encoding='utf-8');os.chmod(out/'server.py',0o755);(out/'README.txt').write_text('Double-click index.html for offline counts. Run serve_explorer.sh, open the HTTP address, choose All candidates if needed, then click Open alignments for an IGV.js report.\n',encoding='utf-8');(out/'coverage_summary.txt').write_text(f'Findings: {len(findings)}\nEmbedded compact records: {len(compact)}\nHTML bytes: {html_bytes}\nFindings discarded: 0\nDatabase files: 0\nDirect-open HTML: yes\n',encoding='utf-8');print(f'Built compact direct-open explorer with {len(findings)} findings and {html_bytes} HTML bytes')
if __name__=='__main__':main()
