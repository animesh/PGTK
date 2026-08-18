#!/usr/bin/env python3
import argparse,csv,gzip,json,os,re
from collections import defaultdict
from pathlib import Path
FIELDS=('EventID','EvidenceClasses','Sample','Gene','Consequence','Impact','Transcript','ProteinChange','Chrom','Position','REF','ALT','UniqueAlignments','ExactAltReads','CleanReferenceReads','ExcludedReads','AltFractionAmongClean')
INTS={'Position','UniqueAlignments','ExactAltReads','CleanReferenceReads','ExcludedReads'}
def safe(x):return re.sub(r'[^A-Za-z0-9_.-]+','_',x or 'none').strip('_') or 'none'
SERVER='''#!/usr/bin/env python3
import gzip,json,os,subprocess,sys,urllib.parse
from http.server import ThreadingHTTPServer,SimpleHTTPRequestHandler
from pathlib import Path
R=Path(__file__).resolve().parent;C=json.loads((R/'explorer_config.json').read_text());K=R/'report_cache';K.mkdir(exist_ok=True);A=[]
for p in C['partitions']:
 with gzip.open(R/p['File'],'rt') as h:A.extend(json.loads(x) for x in h if x.strip())
I={x['EventID']:x for x in A}
def search(q):
 w=q.get('q',[''])[0].lower();s=q.get('sample',[''])[0];c=q.get('class',[''])[0];i=q.get('impact',[''])[0];ch=q.get('chrom',[''])[0];z=[]
 for r in A:
  if s and r['Sample']!=s:continue
  if c and c not in r['EvidenceClasses']:continue
  if i and r['Impact']!=i:continue
  if ch and r['Chrom']!=ch:continue
  if w and w not in ' '.join((r['EventID'],r['Gene'],r['Transcript'],r['Consequence'],r['ProteinChange'])).lower():continue
  z.append(r)
 o=max(0,int(q.get('offset',['0'])[0]));n=min(500,int(q.get('limit',['100'])[0]));return {'total':len(z),'rows':z[o:o+n]}
def report(e):
 out=K/(e+'.html')
 if out.exists():return out
 r=I[e];site=K/(e+'.tsv');site.write_text('Chrom\\tStart\\tEnd\\tEventID\\tGene\\tClass\\tImpact\\tConsequence\\n'+f"{r['Chrom']}\\t{r['Position']}\\t{r['Position']}\\t{r['EventID']}\\t{r['Gene']}\\t{r['EvidenceClasses']}\\t{r['Impact']}\\t{r['Consequence']}\\n")
 subprocess.run([C['report_command'],str(site),'--fasta',C['genome'],'--sequence','Chrom','--begin','Start','--end','End','--tracks',*C['tracks'],'--flanking',str(C['flanking']),'--title','PGTK '+e,'--standalone','--output',str(out)],check=True,cwd=R);return out
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
HTML='''<!doctype html><meta charset="utf-8"><title>PGTK Explorer</title><style>body{font:14px Arial;margin:15px}input,select,button{padding:7px;margin:3px}table{border-collapse:collapse;width:100%;font-size:12px}th,td{padding:5px;border-bottom:1px solid #ddd;text-align:left}iframe{width:100%;height:700px;border:0}</style><h2>PGTK Finding Explorer <span id="n"></span></h2><input id="q" placeholder="Search"><select id="sample"><option value="">All samples</option></select><select id="class"><option value="">All classes</option></select><select id="impact"><option value="">All impacts</option></select><input id="chrom" placeholder="Chromosome"><button id="prev">Previous</button><button id="next">Next</button><span id="page"></span><table><thead><tr><th>View</th><th>Sample</th><th>Class</th><th>Gene</th><th>Impact</th><th>Consequence</th><th>Locus</th><th>ALT</th><th>REF</th></tr></thead><tbody id="rows"></tbody></table><iframe id="viewer"></iframe><script>let o=0,l=100,$=x=>document.getElementById(x);async function init(){let f=await(await fetch('/api/facets')).json();for(let [i,k] of [['sample','samples'],['class','classes'],['impact','impacts']])for(let v of f[k]){let x=document.createElement('option');x.value=x.textContent=v;$(i).appendChild(x)}load()}async function load(){let p=new URLSearchParams({offset:o,limit:l});for(let i of ['q','sample','class','impact','chrom'])if($(i).value)p.set(i,$(i).value);let x=await(await fetch('/api/findings?'+p)).json();$('n').textContent=x.total;$('page').textContent=(x.total?o+1:0)+'-'+Math.min(o+l,x.total);$('rows').innerHTML=x.rows.map(r=>`<tr><td><button onclick="$('viewer').src='/report/${encodeURIComponent(r.EventID)}'">IGV</button></td><td>${r.Sample}</td><td>${r.EvidenceClasses}</td><td>${r.Gene}</td><td>${r.Impact}</td><td>${r.Consequence}</td><td>${r.Chrom}:${r.Position}</td><td>${r.ExactAltReads}</td><td>${r.CleanReferenceReads}</td></tr>`).join('')}for(let i of ['sample','class','impact','chrom'])$(i).onchange=()=>{o=0;load()};$('q').onkeyup=e=>{if(e.key==='Enter'){o=0;load()}};$('prev').onclick=()=>{o=Math.max(0,o-l);load()};$('next').onclick=()=>{o+=l;load()};init()</script>'''
def main():
 a=argparse.ArgumentParser();a.add_argument('--manifest',required=True);a.add_argument('--bam-manifest',required=True);a.add_argument('--genome',required=True);a.add_argument('--output-dir',required=True);a.add_argument('--flanking',type=int,default=150);a.add_argument('--report-command',default='create_report');x=a.parse_args();out=Path(x.output_dir);pd=out/'partitions';pd.mkdir(parents=True,exist_ok=True);hs={};counts=defaultdict(int);fac={'samples':set(),'classes':set(),'impacts':set()};total=0
 try:
  with open(x.manifest,newline='') as h:
   for r in csv.DictReader(h,delimiter='\t'):
    z={k:r.get(k,'') for k in FIELDS}
    for k in INTS:z[k]=int(z[k] or 0)
    z['AltFractionAmongClean']=float(z['AltFractionAmongClean'] or 0);pid='__'.join(safe(z[k]) for k in ('Sample','EvidenceClasses','Chrom'))
    if pid not in hs:hs[pid]=gzip.open(pd/(pid+'.jsonl.gz'),'wt',compresslevel=6)
    hs[pid].write(json.dumps(z,separators=(',',':'))+'\n');counts[pid]+=1;total+=1;fac['samples'].add(z['Sample']);fac['classes'].add(z['EvidenceClasses']);fac['impacts'].add(z['Impact'])
 finally:
  for h in hs.values():h.close()
 tracks=[]
 with open(x.bam_manifest,newline='') as h:
  for r in csv.DictReader(h,delimiter='\t'):
   if r['Category'] in ('exact_alt_display','reference_display','event_display'):tracks.append('../finding_reviews/'+r['BAM'])
 parts=[{'PartitionID':p,'File':'partitions/'+p+'.jsonl.gz','Findings':counts[p]} for p in sorted(counts)]
 with (out/'partition_manifest.tsv').open('w',newline='') as h:w=csv.DictWriter(h,fieldnames=('PartitionID','File','Findings'),delimiter='\t');w.writeheader();w.writerows(parts)
 cfg={'genome':str(Path(x.genome).resolve()),'tracks':tracks,'flanking':x.flanking,'report_command':x.report_command,'findings':total,'partitions':parts,'facets':{k:sorted(v) for k,v in fac.items()}};(out/'explorer_config.json').write_text(json.dumps(cfg,separators=(',',':')));(out/'index.html').write_text(HTML);(out/'server.py').write_text(SERVER);os.chmod(out/'server.py',0o755);(out/'coverage_summary.txt').write_text(f'Findings: {total}\nPartition records: {sum(counts.values())}\nBiological partitions: {len(parts)}\nFindings discarded: 0\nDatabase files: 0\n');print(f'Indexed {total} findings in {len(parts)} database-free partitions; discarded 0')
if __name__=='__main__':main()
