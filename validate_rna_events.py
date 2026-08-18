import argparse, csv, gzip, re, sys
from collections import Counter, defaultdict
from pathlib import Path


def open_text(path):
    p=Path(path)
    return gzip.open(p,'rt',encoding='utf-8',errors='replace') if p.suffix=='.gz' else p.open(encoding='utf-8',errors='replace')

def write_tsv(path, rows, fields):
    with open(path,'w',encoding='utf-8',newline='') as h:
        w=csv.DictWriter(h,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore');w.writeheader();w.writerows(rows)

def parse_info(value):
    return dict(x.split('=',1) for x in value.split(';') if '=' in x)

def parse_format(fmt, sample):
    return dict(zip(fmt.split(':'), sample.split(':')))


def safe_int(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return 0


FAILURE_TEXT = {
    'VCF_FILTER_NOT_PASS': 'The VCF record did not pass upstream variant filtering.',
    'REFERENCE_ALLELE_MISMATCH': 'The VCF REF allele does not match the supplied GRCh38 reference.',
    'REFERENCE_LOOKUP_FAILED': 'The reference allele could not be retrieved from the supplied genome.',
    'RNA_DEPTH_BELOW_MINIMUM': 'RNA read depth is below the configured minimum.',
    'RNA_ALT_READS_BELOW_MINIMUM': 'ALT-supporting RNA reads are below the configured minimum.',
    'RNA_ALT_FRACTION_BELOW_MINIMUM': 'RNA ALT allele fraction is below the configured minimum.',
    'NO_SUPPORTED_PROTEIN_ALTERING_CONSEQUENCE': 'VEP reported no supported protein-altering consequence.',
    'MULTIALLELIC_REQUIRES_NORMALIZATION': "The record contains multiple ALT alleles and is excluded to prevent assigning one allele's support to another.",
    'MISSING_BREAKPOINT': 'One or both fusion breakpoints are missing.',
    'SPLIT_READS_BELOW_MINIMUM': 'Fusion split-read support is below the configured minimum.',
    'TOTAL_SUPPORT_BELOW_MINIMUM': 'Total fusion read support is below the configured minimum.',
    'LOW_CONFIDENCE': 'Arriba classified the fusion as low confidence.',
    'NOT_IN_FRAME': 'Arriba classified the fusion as out of frame, so it is not translated into the primary fusion FASTA.',
    'NO_EXON_JUNCTION': 'The transcript has no exon-exon junction to validate.',
    'JUNCTION_READS_BELOW_MINIMUM': 'At least one exact transcript junction has insufficient CIGAR-N read support.',
}


def explain(reasons, passed):
    if not reasons:
        return passed
    return '; '.join(FAILURE_TEXT.get(code, code) for code in reasons)

def ref_base(genome, chrom, pos, length):
    import pysam
    with pysam.FastaFile(genome) as fasta:
        names=set(fasta.references); plain=chrom[3:] if chrom.startswith('chr') else chrom
        resolved=next((x for x in (chrom,plain,'chr'+plain) if x in names),None)
        if resolved is None: raise ValueError(f'contig not found: {chrom}')
        return fasta.fetch(resolved,pos-1,pos-1+length).upper()

def variant_mode(a):
    headers=[]; rows=[]; accepted=[]; rejected=[]; csq_fields=[]
    with open_text(a.input) as h:
        for line in h:
            if line.startswith('#'):
                headers.append(line)
                if line.startswith('##INFO=<ID=CSQ'):
                    m=re.search(r'Format: ([^">]+)',line);csq_fields=m.group(1).split('|') if m else []
                continue
            f=line.rstrip('\n').split('\t'); chrom,pos,vid,ref,alt,qual,flt,info=f[:8];alts=alt.split(',')
            sample_data=parse_format(f[8],f[9]) if len(f)>9 else {}
            dp=safe_int(sample_data.get('DP'))
            ad=[safe_int(x) for x in sample_data.get('AD','').split(',') if x!='.']
            alt_depth=max(ad[1:] or [0]); af=alt_depth/dp if dp else 0.0
            reasons=[]
            if flt not in {'PASS','.'}:reasons.append('VCF_FILTER_NOT_PASS')
            if len(alts) != 1:reasons.append('MULTIALLELIC_REQUIRES_NORMALIZATION')
            try:
                observed=ref_base(a.genome,chrom,int(pos),len(ref))
                if observed!=ref.upper():reasons.append('REFERENCE_ALLELE_MISMATCH')
            except Exception:reasons.append('REFERENCE_LOOKUP_FAILED')
            if dp<a.min_depth:reasons.append('RNA_DEPTH_BELOW_MINIMUM')
            if alt_depth<a.min_alt_reads:reasons.append('RNA_ALT_READS_BELOW_MINIMUM')
            if af<a.min_alt_fraction:reasons.append('RNA_ALT_FRACTION_BELOW_MINIMUM')
            consequences=[]
            inf=parse_info(info)
            for item in filter(None,inf.get('CSQ','').split(',')):
                ann=dict(zip(csq_fields,item.split('|')))
                consequences += list(filter(None,ann.get('Consequence','').split('&')))
            allowed={'missense_variant','frameshift_variant','stop_gained','stop_lost','start_lost','splice_donor_variant','splice_acceptor_variant','inframe_insertion','inframe_deletion'}
            if not (set(consequences)&allowed):reasons.append('NO_SUPPORTED_PROTEIN_ALTERING_CONSEQUENCE')
            status='RNA_VALIDATED' if not reasons else 'REJECTED'
            row={'Sample':a.sample,'Event type':'variant','Event':f'{chrom}:{pos}:{ref}>{alt}','Initial finding':f'VEP protein-altering candidate {chrom}:{pos}:{ref}>{alt}','Validation rule':'PASS; GRCh38 REF match; supported protein consequence; depth/ALT-read/ALT-fraction thresholds','Observed evidence':f'FILTER={flt};DP={dp};ALT_READS={alt_depth};ALT_FRACTION={af:.6g};CONSEQUENCES={";".join(sorted(set(consequences)))}','Status':status,'Failure codes':';'.join(reasons),'Failure explanation':explain(reasons,'Passed all RNA variant validation rules'),'Required resolution':'none' if not reasons else 'Review the listed failed checks; rejected events are excluded from translation','Depth':dp,'ALT reads':alt_depth,'ALT fraction':f'{af:.6g}','Consequences':';'.join(sorted(set(consequences))),'Source':Path(a.input).name}
            rows.append(row)
            (accepted if status=='RNA_VALIDATED' else rejected).append(line)
    import pysam
    out=Path(a.output_prefix)
    validated_plain=Path(f'{out}.validated.vcf'); rejected_plain=Path(f'{out}.rejected.vcf')
    validated_plain.write_text(''.join(headers+accepted),encoding='utf-8'); rejected_plain.write_text(''.join(headers+rejected),encoding='utf-8')
    pysam.tabix_compress(str(validated_plain),f'{out}.validated.vcf.gz',force=True)
    pysam.tabix_index(f'{out}.validated.vcf.gz',preset='vcf',force=True)
    pysam.tabix_compress(str(rejected_plain),f'{out}.rejected.vcf.gz',force=True)
    validated_plain.unlink(); rejected_plain.unlink()
    fields=['Sample','Event type','Event','Initial finding','Validation rule','Observed evidence','Status','Failure codes','Failure explanation','Required resolution','Depth','ALT reads','ALT fraction','Consequences','Source']
    write_tsv(f'{out}.audit.tsv',rows,fields);write_tsv(f'{out}.validated.tsv',[x for x in rows if x['Status']=='RNA_VALIDATED'],fields);write_tsv(f'{out}.rejected.tsv',[x for x in rows if x['Status']=='REJECTED'],fields)

def fusion_mode(a):
    with open_text(a.input) as h:r=csv.DictReader(h,delimiter='\t');fields=r.fieldnames or [];source=list(r)
    audit=[];valid=[];reject=[]
    for row in source:
        split=safe_int(row.get('split_reads1'))+safe_int(row.get('split_reads2'));discord=safe_int(row.get('discordant_mates'));reasons=[]
        if not row.get('breakpoint1') or not row.get('breakpoint2'):reasons.append('MISSING_BREAKPOINT')
        if split<a.min_split_reads:reasons.append('SPLIT_READS_BELOW_MINIMUM')
        if split+discord<a.min_total_support:reasons.append('TOTAL_SUPPORT_BELOW_MINIMUM')
        if row.get('confidence','').lower()=='low':reasons.append('LOW_CONFIDENCE')
        if row.get('reading_frame','').lower()!='in-frame':reasons.append('NOT_IN_FRAME')
        status='RNA_VALIDATED' if not reasons else 'REJECTED'
        audit.append({'Sample':a.sample,'Event type':'fusion','Event':f"{row.get('#gene1',row.get('gene1',''))}--{row.get('gene2','')}:{row.get('breakpoint1','')}|{row.get('breakpoint2','')}",'Initial finding':f"Arriba fusion candidate {row.get('#gene1',row.get('gene1',''))}--{row.get('gene2','')}",'Validation rule':'Both breakpoints present; split and total support thresholds; confidence above low; in-frame','Observed evidence':f"BREAKPOINTS={row.get('breakpoint1','')}|{row.get('breakpoint2','')};SPLIT={split};DISCORDANT={discord};CONFIDENCE={row.get('confidence','')};FRAME={row.get('reading_frame','')}",'Status':status,'Failure codes':';'.join(reasons),'Failure explanation':explain(reasons,'Passed all RNA fusion validation rules'),'Required resolution':'none' if not reasons else 'Rejected fusion remains reported but is excluded from pVACfuse translation','Split reads':split,'Discordant mates':discord,'Confidence':row.get('confidence',''),'Reading frame':row.get('reading_frame',''),'Source':Path(a.input).name})
        (valid if status=='RNA_VALIDATED' else reject).append(row)
    write_tsv(f'{a.output_prefix}.validated.tsv',valid,fields);write_tsv(f'{a.output_prefix}.rejected.tsv',reject,fields)
    write_tsv(f'{a.output_prefix}.audit.tsv',audit,['Sample','Event type','Event','Initial finding','Validation rule','Observed evidence','Status','Failure codes','Failure explanation','Required resolution','Split reads','Discordant mates','Confidence','Reading frame','Source'])

def cigar_junctions(pos,cigar):
    ref=pos;out=[]
    for length,op in re.findall(r'(\d+)([MIDNSHP=X])',cigar):
        n=int(length)
        if op=='N':out.append((ref,ref+n-1));ref+=n
        elif op in 'MD=X':ref+=n
    return out

def transcript_models(path):
    models=defaultdict(list); lines=defaultdict(list); comments=[]
    with open_text(path) as h:
        for line in h:
            if line.startswith('#'):comments.append(line);continue
            f=line.rstrip('\n').split('\t')
            if len(f)<9:continue
            m=re.search(r'transcript_id "([^"]+)"',f[8]);tx=m.group(1) if m else ''
            if tx:lines[tx].append(line)
            if tx and f[2]=='exon':models[tx].append((f[0],int(f[3]),int(f[4]),f[6]))
    return models,lines,comments

def splice_mode(a):
    import pysam
    counts=Counter()
    with pysam.AlignmentFile(a.bam,'rb') as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate or not read.cigarstring or 'N' not in read.cigarstring: continue
            chrom=bam.get_reference_name(read.reference_id)
            for start,end in cigar_junctions(read.reference_start+1,read.cigarstring): counts[(chrom,start,end)]+=1
    models,lines,comments=transcript_models(a.input);audit=[];valid_lines=[];reject_lines=[]
    for tx,exons in models.items():
        exons=sorted(exons,key=lambda x:x[1]);junctions=[]
        for left,right in zip(exons,exons[1:]):junctions.append((left[0],left[2]+1,right[1]-1))
        support=[counts[j] for j in junctions];reasons=[]
        if not junctions:reasons.append('NO_EXON_JUNCTION')
        if support and min(support)<a.min_junction_reads:reasons.append('JUNCTION_READS_BELOW_MINIMUM')
        status='RNA_VALIDATED' if not reasons else 'REJECTED'
        audit.append({'Sample':a.sample,'Event type':'splice_transcript','Event':tx,'Initial finding':f'StringTie/GFFCompare candidate transcript {tx}','Validation rule':'At least two exons and every transcript junction supported by the minimum exact CIGAR-N read count','Observed evidence':f'JUNCTIONS={len(junctions)};COUNTS={";".join(map(str,support))}','Status':status,'Failure codes':';'.join(reasons),'Failure explanation':explain(reasons,'Passed all RNA splice validation rules'),'Required resolution':'none' if not reasons else 'Rejected transcript remains reported but is excluded from TransDecoder translation','Junctions':';'.join(f'{c}:{s}-{e}' for c,s,e in junctions),'Junction read counts':';'.join(map(str,support)),'Minimum junction reads':min(support) if support else 0,'Source':Path(a.input).name})
        (valid_lines if status=='RNA_VALIDATED' else reject_lines).extend(lines[tx])
    Path(f'{a.output_prefix}.validated.gtf').write_text(''.join(comments+valid_lines),encoding='utf-8')
    Path(f'{a.output_prefix}.rejected.gtf').write_text(''.join(comments+reject_lines),encoding='utf-8')
    fields=['Sample','Event type','Event','Initial finding','Validation rule','Observed evidence','Status','Failure codes','Failure explanation','Required resolution','Junctions','Junction read counts','Minimum junction reads','Source']
    write_tsv(f'{a.output_prefix}.audit.tsv',audit,fields);write_tsv(f'{a.output_prefix}.validated.tsv',[x for x in audit if x['Status']=='RNA_VALIDATED'],fields);write_tsv(f'{a.output_prefix}.rejected.tsv',[x for x in audit if x['Status']=='REJECTED'],fields)

def main():
    p=argparse.ArgumentParser();sub=p.add_subparsers(dest='mode',required=True)
    v=sub.add_parser('variant');v.add_argument('--input',required=True);v.add_argument('--genome',required=True);v.add_argument('--sample',required=True);v.add_argument('--output-prefix',required=True);v.add_argument('--min-depth',type=int,default=10);v.add_argument('--min-alt-reads',type=int,default=3);v.add_argument('--min-alt-fraction',type=float,default=0.05)
    f=sub.add_parser('fusion');f.add_argument('--input',required=True);f.add_argument('--sample',required=True);f.add_argument('--output-prefix',required=True);f.add_argument('--min-split-reads',type=int,default=1);f.add_argument('--min-total-support',type=int,default=2)
    s=sub.add_parser('splice');s.add_argument('--input',required=True);s.add_argument('--bam',required=True);s.add_argument('--sample',required=True);s.add_argument('--output-prefix',required=True);s.add_argument('--min-junction-reads',type=int,default=3)
    a=p.parse_args();{'variant':variant_mode,'fusion':fusion_mode,'splice':splice_mode}[a.mode](a)
if __name__=='__main__':
    try:main()
    except Exception as e:print(f'ERROR: {e}',file=sys.stderr);raise SystemExit(1)
