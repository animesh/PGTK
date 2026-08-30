#!/usr/bin/env python3
import argparse
import bisect
import csv
import hashlib
import re
from collections import defaultdict
from pathlib import Path
import xml.etree.ElementTree as ET
from variant_read_evidence import classify_sam_fields, evidence_status

try:
    import pysam
except ImportError:
    pysam = None

CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')
SUMMARY_FIELDS = ['EventID','EvidenceClasses','SourceEvents','Sources','Sample','Gene','Consequence','Impact','PredictedConsequence','PredictedImpact','Transcript','ProteinChange','Chrom','Position','REF','ALT','ReadValidationStatus','ValidationExplanation','CountUnit','UniqueAlignments','CallableAlignments','ExactAltReads','CleanReferenceReads','ExcludedReads','AltFractionAmongClean','CallableFractionAmongExamined']
READ_FIELDS = ['EventID','Sample','ReadName','Class','Reason','Contig','Start','MapQ','Flag','CIGAR','Observed','Quality']
BAM_CATEGORIES = ('exact_alt_unique','exact_alt_display','reference_display','event_display')


def parse_label(label):
    match = re.search(r'([ACGTN]+)>([ACGTN]+)', (label or '').upper())
    return match.groups() if match else ('','')


def safe_id(event):
    ref, alt = parse_label(event.get('Label',''))
    label = re.sub(r'[^A-Za-z0-9_.-]+','_',event.get('Label','event')).strip('_')
    gene = re.sub(r'[^A-Za-z0-9_.-]+','_',event.get('Gene','') or 'INTERGENIC').strip('_')
    return f"{gene}_{event['Sample']}_{event['Chrom']}_{int(event['Start0'])+1}_{ref}_{alt}_{label}"[:180]


def resolve_contig(alignment, chrom):
    names = set(alignment.references)
    plain = chrom[3:] if chrom.startswith('chr') else chrom
    for candidate in (chrom, plain, 'chr'+plain):
        if candidate in names:
            return candidate
    raise RuntimeError(f'Contig {chrom} is absent from BAM')


def classify(fields,target,ref,alt,window_start,window_end,min_mapq,min_baseq):
    del window_start,window_end
    return classify_sam_fields(fields,target,ref,alt,min_mapq,min_baseq)

def record_key(line):return tuple(line.split('\t')[:10])

def deterministic_subset(records,limit):
    values=sorted(records,key=lambda line:hashlib.sha256(line.split('\t',1)[0].encode()).hexdigest())
    return values if limit<=0 else values[:limit]


def write_bam(header_text,records,output_path):
    if pysam is None:raise RuntimeError('pysam is required')
    output_path=Path(output_path); temporary=output_path.with_suffix('.temporary.bam')
    header=pysam.AlignmentHeader.from_text(header_text.rstrip('\n')+'\n')
    try:
        with pysam.AlignmentFile(temporary,'wb',header=header) as output:
            for line in records:output.write(pysam.AlignedSegment.fromstring(line,header))
        pysam.sort('-o',str(output_path),str(temporary)); pysam.index(str(output_path))
        if pysam.quickcheck(str(output_path))!='':raise RuntimeError(f'HTSlib quickcheck failed: {output_path}')
    finally:temporary.unlink(missing_ok=True)


def write_table(path,fields,rows):
    with Path(path).open('w',encoding='utf-8',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t',lineterminator='\n'); writer.writeheader(); writer.writerows(rows)


def merge_intervals(intervals):
    merged=[]
    for chrom,start,end in sorted(intervals):
        if merged and merged[-1][0]==chrom and start<=merged[-1][2]: merged[-1]=(chrom,merged[-1][1],max(merged[-1][2],end))
        else: merged.append((chrom,start,end))
    return merged

def make_index(items, position_index):
    result={}
    grouped=defaultdict(list)
    for item in items: grouped[item[0]].append(item)
    for chrom,values in grouped.items():
        values.sort(key=lambda x:x[position_index]); result[chrom]=([x[position_index] for x in values],values)
    return result

def overlapping(index,chrom,start,end,pad=0):
    data=index.get(chrom)
    if not data:return ()
    positions,values=data
    return values[bisect.bisect_left(positions,start-pad):bisect.bisect_right(positions,end+pad)]

def browser_safe_display(read, ref='', alt=''):
    if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate or read.is_qcfail:
        return False
    if not read.cigarstring or not read.query_sequence:
        return False
    has_deletion = any(operation == 2 for operation, count in (read.cigartuples or ()))
    is_deletion_variant = bool(ref and alt and len(ref) > len(alt) and ref.startswith(alt))
    return not has_deletion or is_deletion_variant

def browser_display_copy(read, header):
    """Create a standalone primary alignment for a sparse browser display BAM.

    Evidence BAMs retain original records. Display BAMs deliberately remove pair,
    mate, secondary, supplementary, duplicate, and QC-fail state because sparse
    event extraction does not guarantee that a mate is present in the same BAM.
    """
    copied = pysam.AlignedSegment.fromstring(read.to_string(), header)
    copied.flag &= ~(1 | 2 | 4 | 8 | 32 | 64 | 128 | 256 | 512 | 1024 | 2048)
    copied.next_reference_id = -1
    copied.next_reference_start = -1
    copied.template_length = 0
    return copied

def finalize_bam(source,target):
    pysam.sort('-o',str(target),str(source));source.unlink(missing_ok=True)
    Path(str(target)+'.bai').unlink(missing_ok=True);target.with_suffix('.bai').unlink(missing_ok=True)
    pysam.index(str(target))
    if pysam.quickcheck(str(target))!='':raise RuntimeError(f'HTSlib quickcheck failed: {target}')

def main():
    p=argparse.ArgumentParser(description='Single-pass consolidated IGV review')
    p.add_argument('--events',required=True);p.add_argument('--bam',action='append',default=[]);p.add_argument('--genome',required=True);p.add_argument('--output-dir',required=True);p.add_argument('--event-id',default='')
    p.add_argument('--padding',type=int,default=100);p.add_argument('--mapq',type=int,default=20);p.add_argument('--baseq',type=int,default=20);p.add_argument('--reference-display-reads',type=int,default=20);p.add_argument('--alt-display-reads',type=int,default=100)
    p.add_argument('--finding-classes',default='rna_variant,progression_variant,fusion,splice_junction');p.add_argument('--primary-class-order',default='rna_variant,progression_variant');p.add_argument('--priority-mode',choices=('all','filter'),default='all');p.add_argument('--priority-genes',default='');p.add_argument('--priority-impacts',default='');p.add_argument('--priority-consequences',default='');p.add_argument('--priority-limit',type=int,default=0);p.add_argument('--gene-filter',default='');p.add_argument('--sample-filter',default='');p.add_argument('--diagnostic-read-limit',type=int,default=100000);p.add_argument('--excluded-read-limit',type=int,default=10000);p.add_argument('--progress-every-reads',type=int,default=1000000);p.add_argument('--plan-only',action='store_true');a=p.parse_args()
    bams=dict(x.split('=',1) for x in a.bam);out=Path(a.output_dir);out.mkdir(parents=True,exist_ok=True)
    with open(a.events,encoding='utf-8',newline='') as h:events=list(csv.DictReader(h,delimiter='\t'))
    classes={x.strip() for x in a.finding_classes.split(',') if x.strip()};genes={x.strip() for x in a.gene_filter.split(',') if x.strip()};samples={x.strip() for x in a.sample_filter.split(',') if x.strip()};order=[x.strip() for x in a.primary_class_order.split(',') if x.strip()]
    selected=[(safe_id(e),e) for e in events if e['Class'] in classes and (not a.event_id or e['Event']==a.event_id) and (not genes or e.get('Gene','') in genes) and (not samples or e.get('Sample','') in samples)]
    grouped=defaultdict(list)
    for eid,e in selected:grouped[eid].append(e)
    consolidated=[]
    for eid in sorted(grouped):
        members=grouped[eid];e=next((r for preferred in order for r in members if r['Class']==preferred),members[0]).copy();e['_id']=eid;e['EvidenceClasses']=';'.join(sorted({r['Class'] for r in members}));e['SourceEvents']=';'.join(r['Event'] for r in members);e['Sources']=';'.join(sorted({r.get('Source','') for r in members if r.get('Source','')}));consolidated.append(e)
    with (out/'event_consolidation.tsv').open('w',newline='') as h:w=csv.writer(h,delimiter='\t');w.writerow(['EventID','EvidenceClasses','SourceEvents','Sources']);w.writerows([[e['_id'],e['EvidenceClasses'],e['SourceEvents'],e['Sources']] for e in consolidated])
    if a.plan_only:(out/'consolidation_summary.txt').write_text(f'Input selected rows: {len(selected)}\nConsolidated findings: {len(consolidated)}\nMerged duplicate rows: {len(selected)-len(consolidated)}\n');print(f'Planned {len(consolidated)} unique strict findings; merged {len(selected)-len(consolidated)} duplicate source rows');return
    if pysam is None or not bams:raise RuntimeError('pysam and BAM inputs required')
    by_sample=defaultdict(list)
    for e in consolidated:
        if e['Sample'] in bams:by_sample[e['Sample']].append(e)
    counts={e['_id']:{'unique':0,'alt':0,'ref':0,'excluded':0} for e in consolidated};metadata={};regions=[];bed=['track name="PGTK_consolidated_findings" itemRgb="On"'];bam_manifest=[];stats=[]
    rh=(out/'read_classification.tsv').open('w',newline='');eh=(out/'excluded_reads.tsv').open('w',newline='');rw=csv.DictWriter(rh,fieldnames=READ_FIELDS,delimiter='\t');ew=csv.DictWriter(eh,fieldnames=READ_FIELDS,delimiter='\t');rw.writeheader();ew.writeheader();nr=ne=0
    for sample,bampath in sorted(bams.items()):
      with pysam.AlignmentFile(bampath,'rb') as bam:
        variants=[];nonvars=[]
        for e in by_sample[sample]:
          eid=e['_id'];ref,alt=parse_label(e.get('Label',''));target=int(e['Start0'])+1;chrom=resolve_contig(bam,e['Chrom']);start=max(1,target-a.padding);end=target+max(len(ref),len(alt),1)+a.padding;rs=[(chrom,start,end)]
          if e.get('Chrom2'):rs.append((resolve_contig(bam,e['Chrom2']),max(1,int(e['Start2_0'])+1-a.padding),int(e['End2'])+a.padding))
          for c,s,z in rs:regions.append({'EventID':eid,'Sample':sample,'Chrom':c,'Start':s,'End':z})
          metadata[eid]=(e,chrom,target,ref,alt)
          if ref and alt:
            variants.append((chrom,target,eid,e,start,end))
            nonvars.append((chrom,start-1,end))
            bed.append(f'{chrom}\t{target-1}\t{target}\t{eid}\t1000\t.\t{target-1}\t{target}\t255,80,80')
          else:
            for c,s,z in rs:nonvars.append((c,s-1,z));bed.append(f'{c}\t{s-1}\t{z}\t{eid} {e["Class"]}\t1000\t.\t{s-1}\t{z}\t80,120,255')
        vindex=make_index(variants,1);merged=merge_intervals(nonvars);mindex=make_index(merged,1)
        temp={c:out/f'{sample}.{c}.unsorted.bam' for c in BAM_CATEGORIES};writers={c:pysam.AlignmentFile(temp[c],'wb',template=bam) for c in BAM_CATEGORIES};scanned=eventn=exactn=altd=refd=0
        try:
          for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.reference_id<0:continue
            scanned+=1;c=bam.get_reference_name(read.reference_id);start=read.reference_start+1;end=read.reference_end or start
            overlaps=overlapping(mindex,c,read.reference_start,read.reference_end or read.reference_start)
            if any(x[1]<=read.reference_end and x[2]>=read.reference_start for x in overlaps) and browser_safe_display(read):writers['event_display'].write(browser_display_copy(read, bam.header));eventn+=1
            candidates=overlapping(vindex,c,start,end,2)
            if candidates:
              fields=read.to_string().split('\t');wa=wad=wr=False
              for target,eid,e,ws,we in [(x[1],x[2],x[3],x[4],x[5]) for x in candidates]:
                ref,alt=parse_label(e.get('Label',''));cls,reason,obs,qual=classify(fields,target,ref,alt,ws,we,a.mapq,a.baseq);co=counts[eid];co['unique']+=1
                if cls=='EXACT_ALT':co['alt']+=1;wa=True;wad=wad or co['alt']<=a.alt_display_reads
                elif cls=='CLEAN_REFERENCE':co['ref']+=1;wr=wr or co['ref']<=a.reference_display_reads
                else:co['excluded']+=1
                row=dict(zip(READ_FIELDS,[eid,sample,fields[0],cls,reason,fields[2],fields[3],fields[4],fields[1],fields[5],obs,qual]))
                if cls!='EXCLUDED' and nr<a.diagnostic_read_limit:rw.writerow(row);nr+=1
                elif cls=='EXCLUDED' and ne<a.excluded_read_limit:ew.writerow(row);ne+=1
              if wa:writers['exact_alt_unique'].write(read);exactn+=1
              display_ref,display_alt=parse_label(candidates[0][3].get('Label','')) if candidates else ('','')
              if wad and browser_safe_display(read,display_ref,display_alt):writers['exact_alt_display'].write(browser_display_copy(read, bam.header));altd+=1
              if wr and browser_safe_display(read,display_ref,display_alt):writers['reference_display'].write(browser_display_copy(read, bam.header));refd+=1
            if a.progress_every_reads and scanned%a.progress_every_reads==0:print(f'PROGRESS sample={sample} reads={scanned} event={eventn} exact_alt={exactn}',flush=True)
        finally:
          for w in writers.values():w.close()
        nums={'exact_alt_unique':exactn,'exact_alt_display':altd,'reference_display':refd,'event_display':eventn}
        for cat in BAM_CATEGORIES:
          final=out/f'{sample}.{cat}.bam';finalize_bam(temp[cat],final);bam_manifest.append({'Sample':sample,'Category':cat,'BAM':final.name,'Index':final.name+'.bai','UniqueAlignments':nums[cat]})
        stats.append((sample,scanned,len(variants),len(nonvars),len(merged),eventn,exactn));print(f'PROGRESS sample={sample} complete reads={scanned} variants={len(variants)} merged_intervals={len(merged)}',flush=True)
    rh.close();eh.close();write_table(out/'event_regions.tsv',['EventID','Sample','Chrom','Start','End'],regions)
    summaries=[]
    for e in consolidated:
      eid=e['_id']
      if eid not in metadata:continue
      e,chrom,target,ref,alt=metadata[eid];co=counts[eid];status,explanation,callable_count,fraction=evidence_status(co['alt'],co['ref'],co['excluded']);summaries.append({'EventID':eid,'EvidenceClasses':e['EvidenceClasses'],'SourceEvents':e['SourceEvents'],'Sources':e['Sources'],'Sample':e['Sample'],'Gene':e.get('Gene',''),'Consequence':e.get('Consequence',''),'Impact':e.get('Impact',''),'PredictedConsequence':e.get('Consequence',''),'PredictedImpact':e.get('Impact',''),'Transcript':e.get('Transcript',''),'ProteinChange':e.get('ProteinChange',''),'Chrom':chrom,'Position':target,'REF':ref,'ALT':alt,'ReadValidationStatus':status,'ValidationExplanation':explanation,'CountUnit':'primary_alignments','UniqueAlignments':co['unique'],'CallableAlignments':callable_count,'ExactAltReads':co['alt'],'CleanReferenceReads':co['ref'],'ExcludedReads':co['excluded'],'AltFractionAmongClean':f'{fraction:.6f}' if fraction is not None else 'NA','CallableFractionAmongExamined':f'{callable_count/co["unique"]:.6f}' if co['unique'] else 'NA'})
    write_table(out/'findings_manifest.tsv',SUMMARY_FIELDS,summaries);write_table(out/'bam_manifest.tsv',['Sample','Category','BAM','Index','UniqueAlignments'],bam_manifest);(out/'support_labels.bed').write_text('\n'.join(bed)+'\n')
    pg={x.strip() for x in a.priority_genes.split(',') if x.strip()};pi={x.strip() for x in a.priority_impacts.split(',') if x.strip()};pc={x.strip() for x in a.priority_consequences.split(',') if x.strip()}
    priority=[r for r in summaries if a.priority_mode=='all' or (pg and r['Gene'] in pg) or (pi and r['PredictedImpact'] in pi) or (pc and set(re.split(r'[,;&]',r['PredictedConsequence']))&pc)];priority=sorted(priority,key=lambda r:r['EventID']);priority=priority if a.priority_limit==0 else priority[:a.priority_limit];write_table(out/'priority_findings.tsv',SUMMARY_FIELDS,priority);ids={r['EventID'] for r in priority};(out/'priority_findings.bed').write_text('\n'.join([bed[0]]+[x for x in bed[1:] if x.split('\t')[3].split()[0] in ids])+'\n')
    visible=[out/r['BAM'] for r in bam_manifest if r['Category'] in {'exact_alt_display','reference_display','event_display'}];(out/'review.igv.batch.txt').write_text('\n'.join(['new',f'genome {Path(a.genome).resolve()}',f'load {(out/"support_labels.bed").resolve()}']+[f'load {x.resolve()}' for x in visible]+[f"goto {r['Chrom']}:{r['Position']}" for r in priority[:1000]]+['expand','exit'])+'\n');session=ET.Element('Session',genome=str(Path(a.genome).resolve()),version='8');resources=ET.SubElement(session,'Resources');ET.SubElement(resources,'Resource',path=str((out/'support_labels.bed').resolve()));[ET.SubElement(resources,'Resource',path=str(x.resolve())) for x in visible];ET.ElementTree(session).write(out/'igv.session.xml',encoding='utf-8',xml_declaration=True)
    text=[f'Input selected rows: {len(selected)}',f'Consolidated findings: {len(consolidated)}',f'Merged duplicate rows: {len(selected)-len(consolidated)}',f'Generated findings: {len(summaries)}',f'Diagnostic read rows: {nr}',f'Excluded diagnostic rows: {ne}']+[f'{s} scanned reads: {n}; variants: {v}; raw intervals: {ri}; merged intervals: {mi}; event reads: {er}; exact ALT reads: {ea}' for s,n,v,ri,mi,er,ea in stats];(out/'consolidation_summary.txt').write_text('\n'.join(text)+'\n');(out/'README.txt').write_text('Single-pass consolidated IGV review. Each event BAM is scanned once per sample. No SQLite and no repeated per-event BAM fetches.\n');print(f'Generated one single-pass consolidated IGV bundle containing {len(summaries)} findings under {out}')


if __name__=='__main__':main()
