#!/usr/bin/env python3
"""Shared, allele-exact read classification for PGTK.

Coordinates are 1-based VCF coordinates. Counts are primary alignment counts, not
unique molecules. RNA splice skips (CIGAR N) advance the reference but are not
variant observations. Reference calls must be continuously aligned across the
entire comparison interval.
"""
import re
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
EXCLUDED_FLAGS = ((4,"UNMAPPED"),(256,"SECONDARY"),(2048,"SUPPLEMENTARY"),(1024,"DUPLICATE_FLAG"),(512,"QC_FAIL"))

def evidence_status(alt_reads, ref_reads, excluded_reads=0):
    alt_reads,ref_reads,excluded_reads=map(int,(alt_reads,ref_reads,excluded_reads))
    callable_reads=alt_reads+ref_reads
    if alt_reads and ref_reads:
        status='MIXED_ALT_AND_REFERENCE'; explanation=f'{alt_reads} exact ALT and {ref_reads} clean REF primary alignments.'
    elif alt_reads:
        status='ALT_SUPPORTED'; explanation=f'{alt_reads} exact ALT primary alignments and no clean REF alignments.'
    elif ref_reads:
        status='NO_EXACT_ALT_SUPPORT'; explanation=f'No exact ALT alignment; {ref_reads} clean REF primary alignments.'
    else:
        status='NO_CALLABLE_READS'; explanation=f'No alignment passed allele-specific callable criteria; {excluded_reads} excluded observations.'
    fraction=(alt_reads/callable_reads) if callable_reads else None
    return status,explanation,callable_reads,fraction

def classify_sam_fields(f,pos,ref,alt,min_mapq=20,min_baseq=20):
    if len(f)<11:return ('EXCLUDED','MALFORMED_SAM_RECORD','','')
    flag=int(f[1]); start=int(f[3]); mapq=int(f[4]); cigar=f[5]; seq=f[9].upper(); qual=f[10]
    for mask,reason in EXCLUDED_FLAGS:
        if flag&mask:return ('EXCLUDED',reason,'','')
    if mapq<min_mapq:return ('EXCLUDED','LOW_MAPQ','','')
    if cigar in ('','*') or seq in ('','*'):return ('EXCLUDED','MISSING_CIGAR_OR_SEQUENCE','','')
    bases={}; insertions={}; deletions=[]; rp=0; gp=start; block=0
    parsed=CIGAR_RE.findall(cigar)
    if not parsed or ''.join(n+op for n,op in parsed)!=cigar:return ('EXCLUDED','UNSUPPORTED_CIGAR','','')
    for ns,op in parsed:
        n=int(ns)
        if op in 'M=X':
            block+=1
            for i in range(n):
                j=rp+i
                if j<len(seq):
                    q=ord(qual[j])-33 if qual!='*' and j<len(qual) else -1
                    bases[gp+i]=(seq[j],q,block)
            rp+=n;gp+=n
        elif op=='I':
            qtext=qual[rp:rp+n] if qual!='*' else ''
            insertions.setdefault(gp-1,[]).append((seq[rp:rp+n],qtext))
            rp+=n;block+=1
        elif op=='D':deletions.append((gp,gp+n-1));gp+=n;block+=1
        elif op=='N':gp+=n;block+=1
        elif op=='S':rp+=n;block+=1
        elif op in 'HP':block+=1
    ref=ref.upper();alt=alt.upper()
    if not ref or not alt or ref==alt:return ('EXCLUDED','INVALID_ALLELES','','')
    def span(sequence_start,sequence):
        items=[bases.get(sequence_start+i) for i in range(len(sequence))]
        if not all(items):return None,'DOES_NOT_SPAN_COMPLETE_ALLELE'
        if len({x[2] for x in items})!=1:return None,'ALLELE_NOT_CONTIGUOUSLY_ALIGNED'
        observed=''.join(x[0] for x in items);qs=[x[1] for x in items]
        if min(qs)<min_baseq:return None,'LOW_ALLELE_BASE_QUALITY'
        return (observed,qs),'OK'
    # Substitutions and MNVs.
    if len(ref)==len(alt):
        call,reason=span(pos,ref)
        if not call:return ('EXCLUDED',reason,'','')
        observed,qs=call;quality=';'.join(map(str,qs))
        if observed==alt:return ('EXACT_ALT','EXACT_SUBSTITUTION',observed,quality)
        if observed==ref:return ('CLEAN_REFERENCE','EXACT_REFERENCE_HAPLOTYPE',observed,quality)
        return ('EXCLUDED','OTHER_ALLELE',observed,quality)
    # Left-anchored insertion. Require the complete REF anchor haplotype and exact insertion.
    if len(alt)>len(ref) and alt.startswith(ref):
        anchor_end=pos+len(ref)-1; anchor,reason=span(pos,ref)
        if not anchor:return ('EXCLUDED',reason,'','')
        observed_anchor,anchor_qs=anchor
        if observed_anchor!=ref:return ('EXCLUDED','OTHER_ANCHOR_ALLELE',observed_anchor,';'.join(map(str,anchor_qs)))
        expected=alt[len(ref):]
        matches=[(s,q) for s,q in insertions.get(anchor_end,[]) if s==expected]
        if matches:
            qtext=matches[0][1];qs=[ord(c)-33 for c in qtext]
            if not qs or min(qs)<min_baseq:return ('EXCLUDED','LOW_INSERTION_BASE_QUALITY',expected,';'.join(map(str,qs)))
            return ('EXACT_ALT','EXACT_INSERTION',expected,';'.join(map(str,qs)))
        if insertions.get(anchor_end):return ('EXCLUDED','OTHER_INSERTION_AT_BREAKPOINT',';'.join(x[0] for x in insertions[anchor_end]),'')
        following=bases.get(anchor_end+1)
        if not following:return ('EXCLUDED','DOES_NOT_SPAN_BOTH_SIDES_OF_INSERTION','','')
        if following[2]!=bases[anchor_end][2]:return ('EXCLUDED','INSERTION_BREAKPOINT_NOT_CONTIGUOUSLY_ALIGNED','','')
        if following[1]<min_baseq:return ('EXCLUDED','LOW_BREAKPOINT_BASE_QUALITY',following[0],str(following[1]))
        return ('CLEAN_REFERENCE','REFERENCE_CONTIGUOUS_ACROSS_INSERTION_BREAKPOINT',ref,';'.join(map(str,anchor_qs+[following[1]])))
    # Left-anchored deletion. Require the whole retained prefix and exact deleted interval.
    if len(ref)>len(alt) and ref.startswith(alt):
        retained,reason=span(pos,alt)
        if not retained:return ('EXCLUDED',reason,'','')
        observed_retained,retained_qs=retained
        if observed_retained!=alt:return ('EXCLUDED','OTHER_ANCHOR_ALLELE',observed_retained,';'.join(map(str,retained_qs)))
        deleted=(pos+len(alt),pos+len(ref)-1)
        if deleted in deletions:return ('EXACT_ALT','EXACT_DELETION',f'{deleted[0]}-{deleted[1]}',';'.join(map(str,retained_qs)))
        if any(not (d[1]<deleted[0] or d[0]>deleted[1]) for d in deletions):return ('EXCLUDED','OTHER_DELETION_AT_BREAKPOINT','','')
        reference,reason=span(pos,ref)
        if not reference:return ('EXCLUDED',reason,'','')
        observed,qs=reference;quality=';'.join(map(str,qs))
        if observed==ref:return ('CLEAN_REFERENCE','EXACT_REFERENCE_HAPLOTYPE',observed,quality)
        return ('EXCLUDED','OTHER_ALLELE',observed,quality)
    return ('EXCLUDED','UNSUPPORTED_COMPLEX_REPLACEMENT','','')
