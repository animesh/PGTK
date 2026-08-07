#!/usr/bin/env python3
import argparse, csv, datetime as dt, fcntl, json, shutil
from pathlib import Path

p=argparse.ArgumentParser()
p.add_argument('--job-id',required=True); p.add_argument('--exit-code',required=True,type=int)
p.add_argument('--trace',required=True); p.add_argument('--nextflow-log',required=True)
p.add_argument('--slurm-log',required=True); p.add_argument('--results-dir',required=True)
a=p.parse_args()
now=dt.datetime.now().astimezone().isoformat(timespec='seconds')
results=Path(a.results_dir); run_dir=results/'failure_logs'/a.job_id; run_dir.mkdir(parents=True,exist_ok=True)
trace=Path(a.trace); failed=[]
if trace.is_file() and trace.stat().st_size:
    with trace.open(encoding='utf-8',errors='replace',newline='') as h:
        for row in csv.DictReader(h,delimiter='\t'):
            status=(row.get('status') or '').upper()
            exit_value=(row.get('exit') or '').strip()
            if status in {'FAILED','ABORTED'} or (exit_value and exit_value!='0'):
                failed.append(row)
fields=['Timestamp','Job ID','Task ID','Process','Tag','Attempt','Status','Exit','CPUs','Memory','Time','Partition','Realtime','Peak RSS','Work directory','Native ID','Error action']
def pick(row,*names):
    for name in names:
        if row.get(name) not in (None,''): return row[name]
    return ''
records=[]
failed_artifacts=[]
for row in failed:
    records.append({
      'Timestamp':now,'Job ID':a.job_id,'Task ID':pick(row,'task_id','task_id'),'Process':pick(row,'process','name').split(' (',1)[0],
      'Tag':pick(row,'tag','name'),'Attempt':pick(row,'attempt'),'Status':pick(row,'status'),'Exit':pick(row,'exit'),
      'CPUs':pick(row,'cpus'),'Memory':pick(row,'memory'),'Time':pick(row,'time'),'Partition':pick(row,'queue'),'Realtime':pick(row,'realtime'),
      'Peak RSS':pick(row,'peak_rss'),'Work directory':pick(row,'workdir'),'Native ID':pick(row,'native_id'),'Error action':pick(row,'error_action')})
    failed_artifacts.append((row, pick(row,'workdir')))
run_ledger=run_dir/'failure_ledger.tsv'
with run_ledger.open('w',encoding='utf-8',newline='') as h:
    w=csv.DictWriter(h,fieldnames=fields,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(records)
artifact_root=run_dir/'task_attempts'; artifact_root.mkdir(exist_ok=True)
for index,(row,workdir) in enumerate(failed_artifacts,1):
    safe=(pick(row,'name','process') or f'task_{index}').replace('/','_').replace(' ','_').replace('(','').replace(')','')
    target=artifact_root/f'{index:03d}_{safe}'; target.mkdir(exist_ok=True)
    work=Path(workdir)
    for name in ['.command.sh','.command.run','.command.err','.command.out','.exitcode']:
        try:
            source=work/name
            if source.is_file(): shutil.copy2(source,target/name.lstrip('.'))
        except OSError: pass
summary={'timestamp':now,'job_id':a.job_id,'pipeline_exit_code':a.exit_code,'failed_task_attempts':len(records),'trace':str(trace),'status':'SUCCESS' if a.exit_code==0 else 'FAILED'}
(run_dir/'run_summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
for source,name in [(Path(a.nextflow_log),'nextflow.log'),(trace,'pipeline_trace.tsv'),(Path(a.slurm_log),'slurm.log')]:
    try:
        if source.is_file(): shutil.copy2(source,run_dir/name)
    except OSError: pass
history=results/'failure_logs'/'failure_history.tsv'; history.parent.mkdir(parents=True,exist_ok=True)
with history.open('a+',encoding='utf-8',newline='') as h:
    fcntl.flock(h,fcntl.LOCK_EX); h.seek(0,2); empty=h.tell()==0
    w=csv.DictWriter(h,fieldnames=fields,delimiter='\t',lineterminator='\n')
    if empty: w.writeheader()
    w.writerows(records); h.flush(); fcntl.flock(h,fcntl.LOCK_UN)
run_history=results/'failure_logs'/'run_history.tsv'
run_fields=['Timestamp','Job ID','Pipeline exit code','Status','Failed task attempts']
with run_history.open('a+',encoding='utf-8',newline='') as h:
    fcntl.flock(h,fcntl.LOCK_EX); h.seek(0,2); empty=h.tell()==0
    w=csv.DictWriter(h,fieldnames=run_fields,delimiter='\t',lineterminator='\n')
    if empty:w.writeheader()
    w.writerow({'Timestamp':now,'Job ID':a.job_id,'Pipeline exit code':a.exit_code,'Status':summary['status'],'Failed task attempts':len(records)})
    h.flush(); fcntl.flock(h,fcntl.LOCK_UN)
print(f'Failure ledger: {run_ledger}')
print(f'Cumulative failure history: {history}')
