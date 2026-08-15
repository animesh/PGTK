#!/usr/bin/env python3
import argparse, csv, hashlib, re, shlex, sys
try:
    import pysam
except ImportError:
    pysam=None
from collections import defaultdict
from pathlib import Path
import xml.etree.ElementTree as ET
CIGAR_RE=re.compile(r'(\d+)([MIDNSHP=X])')

def sample_from_bam(path):
    m=re.match(r'pgtk_igv\.([^.]+)\.events\.bam$',Path(path).name)
    if not m: raise ValueError(f'Cannot derive sample from {path}')
    return m.group(1)
def parse_label(label):
    m=re.search(r'([ACGTN]+)>([ACGTN]+)',label.upper())
    return m.groups() if m else ('','')
def safe_id(event):
    ref,alt=parse_label(event['Label'])
    label=re.sub(r'[^A-Za-z0-9_.-]+','_',event.get('Label','event')).strip('_')
    gene=re.sub(r'[^A-Za-z0-9_.-]+','_',event.get('Gene','') or 'INTERGENIC').strip('_')
    return f"{gene}_{event['Sample']}_{event['Chrom']}_{int(event['Start0'])+1}_{ref}_{alt}_{label}"[:180]
def contig_for(bam, chrom):
    if pysam is None:raise RuntimeError('pysam is required for BAM operations')
    with pysam.AlignmentFile(bam,'rb') as alignment:names=set(alignment.references)
    plain=chrom[3:] if chrom.startswith('chr') else chrom
    for candidate in (chrom,plain,'chr'+plain):
        if candidate in names:return candidate
    raise RuntimeError(f'Contig {chrom} not present in {bam}')
def classify(fields,target,ref,alt,window_start,window_end,min_mapq,min_baseq):
    flag=int(fields[1]); start=int(fields[3]); mapq=int(fields[4]); cigar=fields[5]; seq=fields[9].upper(); qual=fields[10]
    if flag&4:return 'EXCLUDED','UNMAPPED','',''
    if flag&256:return 'EXCLUDED','SECONDARY','',''
    if flag&2048:return 'EXCLUDED','SUPPLEMENTARY','',''
    if flag&1024:return 'EXCLUDED','DUPLICATE_FLAG','',''
    if mapq<min_mapq:return 'EXCLUDED','LOW_MAPQ','',''
    rp=0; gp=start; bases={}; insertions={}; deletions=[]; disruptive=[]
    for ntext,op in CIGAR_RE.findall(cigar):
        n=int(ntext)
        if op in 'M=X':
            for i in range(n):
                pos=gp+i
                if window_start<=pos<=window_end and rp+i<len(seq):
                    q=ord(qual[rp+i])-33 if qual!='*' and rp+i<len(qual) else -1
                    bases[pos]=(seq[rp+i],q)
            rp+=n; gp+=n
        elif op=='I':
            anchor=gp-1; ins=seq[rp:rp+n]; iq=qual[rp:rp+n] if qual!='*' else ''
            insertions.setdefault(anchor,[]).append((ins,iq))
            if window_start-1<=anchor<=window_end: disruptive.append(f'{n}I@{anchor}:{ins}')
            rp+=n
        elif op=='D':
            deletions.append((gp,gp+n-1));
            if not(gp+n-1<window_start or gp>window_end):disruptive.append(f'{n}D@{gp}')
            gp+=n
        elif op=='N':
            if not(gp+n-1<window_start or gp>window_end):disruptive.append(f'{n}N@{gp}')
            gp+=n
        elif op=='S': rp+=n
    def qualities_ok(values): return values and min(values)>=min_baseq
    if len(ref)==1 and len(alt)==1:
        observed=bases.get(target)
        if not observed:return 'EXCLUDED','DOES_NOT_SPAN_TARGET','',''
        base,q=observed
        if q<min_baseq:return 'EXCLUDED','LOW_TARGET_BASE_QUALITY',base,str(q)
        if base==alt:return 'EXACT_ALT','EXACT_SNP',base,str(q)
        if base==ref:return 'CLEAN_REFERENCE','EXACT_REFERENCE',base,str(q)
        return 'EXCLUDED','OTHER_ALLELE',base,str(q)
    if len(alt)>len(ref) and alt.startswith(ref):
        ins=alt[len(ref):]
        candidate_anchors=(target,target-1,target+1)
        matches=[(anchor,x[1]) for anchor in candidate_anchors for x in insertions.get(anchor,[]) if x[0]==ins]
        if matches:
            matched_anchor,iq=matches[0]; qs=[ord(x)-33 for x in iq]
            if not qualities_ok(qs):return 'EXCLUDED','LOW_INSERTION_BASE_QUALITY',ins,iq
            others=[x for x in disruptive if not x.endswith(':'+ins)]
            if others:return 'EXCLUDED','ALT_WITH_OTHER_LOCAL_EVENT:'+','.join(others),ins,iq
            return 'EXACT_ALT',f'EXACT_INSERTION_ANCHOR_{matched_anchor}',ins,iq
        ref_obs=bases.get(target)
        if ref_obs and ref_obs[0]==ref and ref_obs[1]>=min_baseq and not disruptive:return 'CLEAN_REFERENCE','REFERENCE_NO_INSERTION',ref,str(ref_obs[1])
        return 'EXCLUDED','NO_CLEAN_INSERTION_COMPARISON','',''
    if len(ref)>len(alt) and ref.startswith(alt):
        deleted_start=target+len(alt); deleted_end=target+len(ref)-1
        if any(a==deleted_start and b==deleted_end for a,b in deletions):return 'EXACT_ALT','EXACT_DELETION',f'{deleted_start}-{deleted_end}',''
        needed=range(target,target+len(ref)); obs=[bases.get(x) for x in needed]
        if all(obs) and ''.join(x[0] for x in obs)==ref and qualities_ok([x[1] for x in obs]) and not disruptive:return 'CLEAN_REFERENCE','EXACT_REFERENCE_HAPLOTYPE',ref,''
        return 'EXCLUDED','NO_CLEAN_DELETION_COMPARISON','',''
    # normalized complex replacement
    needed=range(target,target+len(ref)); obs=[bases.get(x) for x in needed]
    if all(obs) and ''.join(x[0] for x in obs)==ref and qualities_ok([x[1] for x in obs]) and not disruptive:return 'CLEAN_REFERENCE','EXACT_REFERENCE_HAPLOTYPE',ref,''
    return 'EXCLUDED','COMPLEX_ALLELE_REQUIRES_MANUAL_REVIEW','',''
def write_bam(headers,records,path):
    if pysam is None:raise RuntimeError('pysam is required for BAM operations')
    path=Path(path);temporary=path.with_suffix('.tmp.bam');header=pysam.AlignmentHeader.from_text('\n'.join(headers)+'\n')
    try:
        with pysam.AlignmentFile(temporary,'wb',header=header) as output:
            for line in records:output.write(pysam.AlignedSegment.fromstring(line,header))
        pysam.sort('-o',str(path),str(temporary));pysam.index(str(path))
        if pysam.quickcheck(str(path))!='':raise RuntimeError(f'HTSlib quickcheck failed: {path}')
    finally:temporary.unlink(missing_ok=True)

def review(event,bam,outroot,genome,padding,min_mapq,min_baseq,display_n,script_path):
    ref,alt=parse_label(event['Label'])
    if not ref or not alt:return None
    event_id=safe_id(event); out=outroot/event_id; out.mkdir(parents=True,exist_ok=True)
    target=int(event['Start0'])+1; chrom=contig_for(bam,event['Chrom']); start=max(1,target-padding); end=target+max(len(ref),len(alt))+padding
    with pysam.AlignmentFile(bam,'rb') as alignment:
        headers=alignment.text.rstrip('\n').splitlines();records=[read.to_string() for read in alignment.fetch(chrom,start-1,end)]
    unique={tuple(x.split('\t')[:10]):x for x in records}; records=list(unique.values()); classes=defaultdict(list); rows=[]
    for line in records:
        f=line.split('\t'); cls,reason,observed,quality=classify(f,target,ref,alt,start,end,min_mapq,min_baseq); classes[cls].append(line); rows.append([event_id,event['Sample'],f[0],cls,reason,f[2],f[3],f[4],f[1],f[5],observed,quality])
    alt_bam=out/f'{event_id}.exact_alt.unique.bam'; ref_bam=out/f'{event_id}.reference.clean.bam'; display_bam=out/f'{event_id}.reference.display.bam'
    write_bam(headers,classes['EXACT_ALT'],alt_bam); write_bam(headers,classes['CLEAN_REFERENCE'],ref_bam)
    display=sorted(classes['CLEAN_REFERENCE'],key=lambda x:hashlib.sha256(x.split('\t',1)[0].encode()).hexdigest())[:display_n]; write_bam(headers,display,display_bam)
    with (out/f'{event_id}.read_classification.tsv').open('w',newline='') as h:w=csv.writer(h,delimiter='\t');w.writerow(['EventID','Sample','ReadName','Class','Reason','Contig','Start','MapQ','Flag','CIGAR','Observed','Quality']);w.writerows(rows)
    with (out/f'{event_id}.excluded.tsv').open('w',newline='') as h:w=csv.writer(h,delimiter='\t');w.writerow(['ReadName','Reason','Contig','Start','MapQ','Flag','CIGAR','Observed','Quality']);w.writerows([[r[2],r[4],*r[5:]] for r in rows if r[3]=='EXCLUDED'])
    alt_n=len(classes['EXACT_ALT']); ref_n=len(classes['CLEAN_REFERENCE']); den=alt_n+ref_n
    summary={'EventID':event_id,'EvidenceClasses':event.get('EvidenceClasses',event.get('Class','')),'SourceEvents':event.get('SourceEvents',event.get('Event','')),'Sources':event.get('Sources',event.get('Source','')),'Sample':event['Sample'],'Gene':event.get('Gene',''),'Consequence':event.get('Consequence',''),'Impact':event.get('Impact',''),'Transcript':event.get('Transcript',''),'ProteinChange':event.get('ProteinChange',''),'Chrom':chrom,'Position':target,'REF':ref,'ALT':alt,'UniqueAlignments':len(records),'ExactAltReads':alt_n,'CleanReferenceReads':ref_n,'ExcludedReads':len(classes['EXCLUDED']),'AltFractionAmongClean':f'{alt_n/den:.6f}' if den else '0.000000'}
    with (out/f'{event_id}.clean_summary.tsv').open('w',newline='') as h:w=csv.DictWriter(h,fieldnames=summary,delimiter='\t');w.writeheader();w.writerow(summary)
    bed=out/f'{event_id}.support_labels.bed'; bed.write_text(f'track name="{event_id}_support" itemRgb="On"\n{chrom}\t{target-1}\t{target}\t{event_id} ALT={alt_n} REF={ref_n}\t1000\t.\t{target-1}\t{target}\t255,80,80\n')
    batch=out/f'{event_id}.review.igv.batch.txt'; batch.write_text(f'new\ngenome {genome}\nload {bed.resolve()}\nload {alt_bam.resolve()}\nload {display_bam.resolve()}\ngoto {chrom}:{start}-{end}\nsort base\nexpand\nsnapshotDirectory {out.resolve()}\nsnapshot {event_id}.png\nexit\n')
    session=ET.Element('Session',genome=str(genome),version='8'); resources=ET.SubElement(session,'Resources'); [ET.SubElement(resources,'Resource',path=str(x.resolve())) for x in (bed,alt_bam,display_bam)]; ET.ElementTree(session).write(out/f'{event_id}.igv.session.xml',encoding='utf-8',xml_declaration=True)
    rerun=out/'review_event.sh'
    rerun.write_text(f'''#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd -P)
PROJECT_DIR=$(cd "$HERE/../../../../.." && pwd -P)
ALL_EVIDENCE=$(cd "$HERE/../../../all_evidence" && pwd -P)
PYTHON=${{HOST_PYTHON:-$(command -v python3)}}
"$PYTHON" "$PROJECT_DIR/build_finding_igv_reviews.py" \
  --events "$ALL_EVIDENCE/pgtk_igv.events.tsv" \
  --bam "{event['Sample']}=$ALL_EVIDENCE/pgtk_igv.{event['Sample']}.events.bam" \
  --genome "${{IGV_GENOME:-hg38}}" \
  --output-dir "$HERE/.." --event-id {shlex.quote(event['Event'])} \
  --padding {padding} --mapq {min_mapq} --baseq {min_baseq} --reference-display-reads {display_n}
''')
    rerun.chmod(0o755)
    (out/'README.txt').write_text(f"Event: {event_id}\nAllele: {ref}>{alt}\nRegion: {chrom}:{start}-{end}\nExact ALT reads: {alt_n}\nClean reference reads: {ref_n}\nExcluded reads: {len(classes['EXCLUDED'])}\nALT fraction among clean: {summary['AltFractionAmongClean']}\nRNA evidence requires cautious interpretation and matched DNA validation.\n")
    return summary

def main():
    p=argparse.ArgumentParser();p.add_argument('--events',required=True);p.add_argument('--bam',action='append',default=[]);p.add_argument('--genome',required=True);p.add_argument('--output-dir',required=True);p.add_argument('--event-id',default='');p.add_argument('--padding',type=int,default=100);p.add_argument('--mapq',type=int,default=20);p.add_argument('--baseq',type=int,default=20);p.add_argument('--reference-display-reads',type=int,default=20);p.add_argument('--finding-classes',default='rna_variant,progression_variant');p.add_argument('--primary-class-order',default='rna_variant,progression_variant');p.add_argument('--priority-genes',default='');p.add_argument('--priority-impacts',default='');p.add_argument('--priority-consequences',default='');p.add_argument('--priority-mode',choices=('all','filter'),default='all');p.add_argument('--plan-only',action='store_true');a=p.parse_args()
    bams=dict(x.split('=',1) for x in a.bam); out=Path(a.output_dir);out.mkdir(parents=True,exist_ok=True)
    with open(a.events,newline='') as h: events=list(csv.DictReader(h,delimiter='\t'))
    finding_classes={value.strip() for value in a.finding_classes.split(',') if value.strip()}
    primary_order=[value.strip() for value in a.primary_class_order.split(',') if value.strip()]
    priority_genes={value.strip() for value in a.priority_genes.split(',') if value.strip()}
    priority_impacts={value.strip() for value in a.priority_impacts.split(',') if value.strip()}
    priority_consequences={value.strip() for value in a.priority_consequences.split(',') if value.strip()}
    if not finding_classes: raise ValueError('--finding-classes must contain at least one class')
    selected=[]
    for event in events:
        if a.event_id and event['Event']!=a.event_id:continue
        if event['Class'] in finding_classes:selected.append((safe_id(event),event))
    grouped=defaultdict(list)
    for event_id,event in selected:grouped[event_id].append(event)
    consolidated=[]
    for event_id in sorted(grouped):
        members=grouped[event_id]
        primary=next((row for preferred in primary_order for row in members if row['Class']==preferred),members[0]).copy()
        primary['EvidenceClasses']=';'.join(sorted({row['Class'] for row in members}));primary['SourceEvents']=';'.join(row['Event'] for row in members)
        primary['Sources']=';'.join(sorted({row.get('Source','') for row in members if row.get('Source','')}));primary['_events']=a.events;consolidated.append(primary)
    duplicate_rows=len(selected)-len(consolidated)
    with (out/'event_consolidation.tsv').open('w',newline='') as handle:
        writer=csv.writer(handle,delimiter='\t');writer.writerow(['EventID','EvidenceClasses','SourceEvents','Sources'])
        for event in consolidated:writer.writerow([safe_id(event),event['EvidenceClasses'],event['SourceEvents'],event['Sources']])
    if a.plan_only:
        (out/'consolidation_summary.txt').write_text(f'Input variant rows: {len(selected)}\nConsolidated findings: {len(consolidated)}\nMerged duplicate rows: {duplicate_rows}\n');print(f'Planned {len(consolidated)} unique strict findings; merged {duplicate_rows} duplicate source rows');return
    summaries=[]
    for event in consolidated:
        if event['Sample'] not in bams:continue
        result=review(event,Path(bams[event['Sample']]),out,a.genome,a.padding,a.mapq,a.baseq,a.reference_display_reads,__file__)
        if result:summaries.append(result)
    fields=['EventID','EvidenceClasses','SourceEvents','Sources','Sample','Gene','Consequence','Impact','Transcript','ProteinChange','Chrom','Position','REF','ALT','UniqueAlignments','ExactAltReads','CleanReferenceReads','ExcludedReads','AltFractionAmongClean']
    with (out/'findings_manifest.tsv').open('w',newline='') as h:w=csv.DictWriter(h,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(summaries)
    priority=out/'priority_batches.txt'
    def priority_match(row):
        if a.priority_mode == 'all': return True
        consequences={value for value in re.split(r'[,&;]',row.get('Consequence','')) if value}
        return ((priority_genes and row.get('Gene') in priority_genes) or
                (priority_impacts and row.get('Impact') in priority_impacts) or
                (priority_consequences and bool(consequences & priority_consequences)))
    def priority_key(row):
        return (0 if priority_genes and row.get('Gene') in priority_genes else 1,
                0 if priority_impacts and row.get('Impact') in priority_impacts else 1,
                0 if priority_consequences and bool({value for value in re.split(r'[,&;]',row.get('Consequence','')) if value} & priority_consequences) else 1,
                row['EventID'])
    priority_rows=sorted((r for r in summaries if priority_match(r)),key=priority_key)
    priority.write_text(''.join(str((out/r['EventID']/f"{r['EventID']}.review.igv.batch.txt").resolve())+'\n' for r in priority_rows))
    launcher=out/'run_all_reviews.sh'; launcher.write_text('#!/usr/bin/env bash\nset -euo pipefail\nfind "$(cd "$(dirname "$0")" && pwd -P)" -mindepth 2 -maxdepth 2 -name review_event.sh -print0 | sort -z | xargs -0 -n1 bash\n');launcher.chmod(0o755)
    print(f'Generated {len(summaries)} strict finding reviews under {out}')
if __name__=='__main__':main()
