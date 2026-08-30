#!/usr/bin/env python3
import argparse
import bisect
import csv
import gzip
import re
import pysam
import sys
import tempfile
from collections import Counter
from pathlib import Path
from variant_read_evidence import classify_sam_fields

CODE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S','TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T','AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A','GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}

def open_text(path):
    path = Path(path)
    return gzip.open(path, 'rt', encoding='utf-8', errors='replace') if path.suffix == '.gz' else path.open('rt', encoding='utf-8', errors='replace')

def parse_bams(items):
    result = {}
    for item in items:
        if '=' not in item:
            raise ValueError(f'Invalid --bam value: {item}; expected SAMPLE=FILE')
        sample, path = item.split('=', 1)
        result[sample] = Path(path)
    return result

def sample_for(path, samples):
    matches = [sample for sample in samples if re.search(rf'(?<![A-Za-z0-9]){re.escape(sample)}(?![A-Za-z0-9])', Path(path).name, re.I)]
    if len(matches) != 1:
        raise ValueError(f'Cannot map {path} uniquely to sample: {matches}')
    return matches[0]

def parse_info(value):
    return dict(item.split('=', 1) for item in value.split(';') if '=' in item)

def pileup_counts(text, reference):
    counts = Counter(); index = 0; reference = reference.upper()
    while index < len(text):
        char = text[index]
        if char == '^': index += 2; continue
        if char == '$': index += 1; continue
        if char in '.,': counts[reference] += 1; index += 1; continue
        if char in 'ACGTNacgtn':
            base = char.upper(); counts[base] += 1; index += 1; continue
        if char in '+-':
            sign = char; index += 1
            match = re.match(r'(\d+)', text[index:])
            if not match: continue
            length = int(match.group(1)); index += len(match.group(1))
            sequence = text[index:index + length].upper(); index += length
            counts[reference + sign + sequence] += 1
            continue
        index += 1
    return counts

def translate(codon):
    codon = re.sub('[^ACGT]', '', (codon or '').upper())
    return CODE.get(codon, '') if len(codon) == 3 else ''

def codon_check(csq):
    codons = (csq.get('Codons') or '').replace('-', '').upper()
    amino = (csq.get('Amino_acids') or '').upper()
    ref_codon, alt_codon = codons.split('/', 1) if '/' in codons else (codons, '')
    ref_aa, alt_aa = amino.split('/', 1) if '/' in amino else (amino, '')
    ref_translation, alt_translation = translate(ref_codon), translate(alt_codon)
    terms = set(filter(None, (csq.get('Consequence') or '').split('&')))
    if any(term in terms for term in {'frameshift_variant','inframe_insertion','inframe_deletion','protein_altering_variant'}):
        status, eligible = 'COMPLEX_CONSEQUENCE_NOT_SINGLE_CODON_COMPARABLE', 'no'
    elif len(ref_codon) != 3 or len(alt_codon) != 3 or not ref_translation or not alt_translation:
        status, eligible = 'CODON_NOT_SINGLE_TRIPLET', 'no'
    elif 'synonymous_variant' in terms and ref_translation == alt_translation and ref_aa in {'', ref_translation} and alt_aa in {'', alt_translation}:
        status, eligible = 'SYNONYMOUS_CODON_TRANSLATION_CONFIRMED', 'yes'
    elif ref_aa and alt_aa and (ref_translation, alt_translation) == (ref_aa, alt_aa):
        status, eligible = 'CODON_TRANSLATION_MATCH', 'yes'
    else:
        status, eligible = 'CODON_TRANSLATION_MISMATCH', 'no'
    return ref_codon, alt_codon, ref_aa, alt_aa, ref_translation, alt_translation, status, eligible

def cigar_ops(cigar):
    return [(int(length), operation) for length, operation in re.findall(r'(\d+)([MIDNSHP=X])', cigar)]

def alignment_observation(fields,position,ref,alt,min_base_quality):
    c,r,o,q=classify_sam_fields(fields,position,ref,alt,0,min_base_quality)
    return alt.upper() if c=='EXACT_ALT' else ref.upper() if c=='CLEAN_REFERENCE' else ''

def resolve_contig(alignment, chrom):
    names = set(alignment.references)
    plain = chrom[3:] if chrom.startswith('chr') else chrom
    resolved = next((name for name in (chrom, plain, 'chr' + plain) if name in names), None)
    if resolved is None:
        raise ValueError(f'contig not found in BAM: {chrom}')
    return resolved


def read_fields(read, alignment):
    sequence = read.query_sequence or '*'
    qualities = read.qual or '*'
    return [
        read.query_name or '', str(read.flag), alignment.get_reference_name(read.reference_id),
        str(read.reference_start + 1), str(read.mapping_quality), read.cigarstring or '*',
        '*', '0', '0', sequence, qualities,
    ]


def stream_alignment_counts(bam, records, minimum_base_quality, minimum_mapping_quality, threads, bed=None):
    del bed
    variants_by_chrom = {}
    for chrom, pos, _vid, ref, alts, _info in records:
        for alt in alts.split(','):
            variants_by_chrom.setdefault(chrom, {}).setdefault(pos, []).append((ref, alt))
    counts = {}
    with pysam.AlignmentFile(str(bam), 'rb', threads=max(1, threads)) as alignment:
        for chrom, position_map in variants_by_chrom.items():
            bam_chrom = resolve_contig(alignment, chrom)
            for position in sorted(position_map):
                seen = set()
                max_ref_length = max(len(ref) for ref, _alt in position_map[position])
                fetch_end = position - 1 + max(1, max_ref_length)
                for read in alignment.fetch(bam_chrom, position - 1, fetch_end):
                    if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate:
                        continue
                    if read.mapping_quality < minimum_mapping_quality or not read.cigarstring:
                        continue
                    read_key = (read.query_name, read.flag, read.reference_start, read.cigarstring)
                    if read_key in seen:
                        continue
                    seen.add(read_key)
                    fields = read_fields(read, alignment)
                    for ref, alt in position_map[position]:
                        observed = alignment_observation(fields, position, ref, alt, minimum_base_quality)
                        if observed:
                            counts.setdefault((chrom, position, ref, alt), Counter())[observed] += 1
    return counts

def write_tsv(path, rows, fields):
    with Path(path).open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fields, delimiter='\t', lineterminator='\n', extrasaction='ignore')
        writer.writeheader(); writer.writerows(rows)

def main():
    parser = argparse.ArgumentParser(description='Independently validate VCF REF against GRCh38, ALT in sample RNA reads, and VEP codon translation.')
    parser.add_argument('--vcf', nargs='+', required=True)
    parser.add_argument('--bam', nargs='+', required=True, metavar='SAMPLE=FILE')
    parser.add_argument('--genome', required=True)
    parser.add_argument('--min-base-quality', type=int, default=20)
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--min-mapping-quality', type=int, default=20)
    parser.add_argument('--min-alt-reads', type=int, default=3)
    parser.add_argument('--min-alt-fraction', type=float, default=0.05)
    parser.add_argument('--output-prefix', default='variant_codon_validation')
    args = parser.parse_args()
    bams = parse_bams(args.bam)
    rows = []
    for vcf in args.vcf:
        sample = sample_for(vcf, bams)
        records = []; csq_fields = None
        with open_text(vcf) as handle:
            for line in handle:
                if line.startswith('##INFO=<ID=CSQ'):
                    match = re.search(r'Format: ([^">]+)', line)
                    if match: csq_fields = match.group(1).strip().split('|')
                    continue
                if line.startswith('#'): continue
                fields = line.rstrip('\n').split('\t')
                if len(fields) >= 8:
                    records.append((fields[0], int(fields[1]), fields[2], fields[3], fields[4], parse_info(fields[7])))
        with tempfile.TemporaryDirectory(prefix='variant_codon_') as temporary:
            temporary = Path(temporary); bed = temporary / 'sites.bed'; regions = temporary / 'regions.txt'
            with bed.open('w') as bed_handle, regions.open('w') as region_handle:
                for chrom, pos, _vid, ref, _alts, _info in records:
                    bed_handle.write(f'{chrom}\t{pos-1}\t{pos}\n'); region_handle.write(f'{chrom}:{pos}-{pos+len(ref)-1}\n')
            references = []
            with pysam.FastaFile(args.genome) as fasta:
                names=set(fasta.references)
                for chrom,pos,_vid,ref,_alts,_info in records:
                    plain=chrom[3:] if chrom.startswith('chr') else chrom
                    resolved=next((x for x in (chrom,plain,'chr'+plain) if x in names),None)
                    if resolved is None: raise ValueError(f'contig not found: {chrom}')
                    references.append(fasta.fetch(resolved,pos-1,pos-1+len(ref)).upper())
            alignment_counts = stream_alignment_counts(
                bams[sample], records, args.min_base_quality, args.min_mapping_quality,
                args.threads, bed,
            )
        if len(references) != len(records):
            raise ValueError(f'faidx returned {len(references)} regions for {len(records)} VCF records in {vcf}')
        for record_index, (chrom, pos, variant_id, ref, alts, info) in enumerate(records):
            genome_ref = references[record_index]
            genome_status = 'GENOME_REF_VALIDATED' if genome_ref == ref.upper() else 'GENOME_REF_MISMATCH'
            for alt in alts.split(','):
                counts = alignment_counts.get((chrom, pos, ref, alt), Counter())
                alt_key = alt.upper() if (
                    len(ref) == len(alt) == 1
                    or (len(alt) > len(ref) and alt.startswith(ref))
                    or (len(ref) > len(alt) and ref.startswith(alt))
                ) else ''
                ref_reads = counts[ref.upper()]; alt_reads = counts[alt_key] if alt_key else 0
                depth = sum(counts.values())
                informative = ref_reads + alt_reads; alt_fraction = alt_reads / informative if informative else None
                if not alt_key: read_status = 'COMPLEX_ALLELE_NOT_COUNTED'
                elif alt_reads < args.min_alt_reads: read_status = 'ALT_READS_BELOW_MINIMUM'
                elif alt_fraction is not None and alt_fraction < args.min_alt_fraction: read_status = 'ALT_FRACTION_BELOW_MINIMUM'
                else: read_status = 'ALT_READS_VALIDATED'
                annotations = []
                for text in info.get('CSQ','').split(','):
                    values = text.split('|'); csq = dict(zip(csq_fields or [], values + [''] * max(0, len(csq_fields or []) - len(values))))
                    if csq.get('Allele') == alt or len(alts.split(',')) == 1: annotations.append(csq)
                if not annotations: annotations = [{}]
                for csq in annotations:
                    values = codon_check(csq); ref_codon, alt_codon, ref_aa, alt_aa, translated_ref, translated_alt, codon_status, codon_eligible = values
                    failures = []
                    if genome_status != 'GENOME_REF_VALIDATED': failures.append(genome_status)
                    if read_status != 'ALT_READS_VALIDATED': failures.append(read_status)
                    if codon_status == 'CODON_TRANSLATION_MISMATCH': failures.append(codon_status)
                    primary_failure = failures[0] if failures else ''
                    if failures:
                        overall_status = 'VALIDATION_FAILED'
                    elif codon_eligible == 'yes':
                        overall_status = 'VARIANT_CODON_VALIDATED'
                    else:
                        overall_status = 'VARIANT_EVIDENCE_PARTIAL'
                    rows.append({
                        'Sample':sample,'VCF':Path(vcf).name,'Chromosome':chrom,'Position':pos,'ID':variant_id,'REF':ref,'ALT':alt,
                        'Genome sequence at REF locus':genome_ref,'Genome REF status':genome_status,'BAM':bams[sample].name,'Pileup depth':depth,
                        'REF-supporting reads':ref_reads,'ALT-supporting reads':alt_reads,'Informative allele reads':informative,'ALT fraction':f'{alt_fraction:.6g}' if alt_fraction is not None else 'NA',
                        'Pileup ALT key':alt_key,'Pileup allele counts':';'.join(f'{key}:{value}' for key,value in sorted(counts.items())),'Read ALT status':read_status,
                        'Gene':csq.get('SYMBOL',''),'Transcript':csq.get('Feature',''),'Protein':csq.get('ENSP',''),'Consequence':csq.get('Consequence',''),
                        'HGVSc':csq.get('HGVSc',''),'HGVSp':csq.get('HGVSp',''),'VEP codons':csq.get('Codons',''),'Reference codon':ref_codon,
                        'Alternate codon':alt_codon,'VEP reference amino acid':ref_aa,'VEP alternate amino acid':alt_aa,
                        'Translated reference codon':translated_ref,'Translated alternate codon':translated_alt,'Codon translation status':codon_status,
                        'Codon strict-integration eligible':codon_eligible,'Overall status':overall_status,
                        'Primary failure code':primary_failure,'Failure codes':';'.join(failures),
                    })
    fields = list(rows[0]) if rows else ['Sample','Chromosome','Position','REF','ALT','Overall status','Failure codes']
    prefix = Path(args.output_prefix)
    write_tsv(f'{prefix}.all.tsv', rows, fields)
    write_tsv(f'{prefix}.validated.tsv', [row for row in rows if row['Overall status']=='VARIANT_CODON_VALIDATED'], fields)
    write_tsv(f'{prefix}.partial.tsv', [row for row in rows if row['Overall status']=='VARIANT_EVIDENCE_PARTIAL'], fields)
    write_tsv(f'{prefix}.failed.tsv', [row for row in rows if row['Overall status']=='VALIDATION_FAILED'], fields)
    category_rows = []
    for dimension, values in [('overall_status', Counter(row['Overall status'] for row in rows)), ('primary_failure', Counter(row['Primary failure code'] for row in rows if row['Primary failure code'])), ('failure_occurrence', Counter(code for row in rows for code in row['Failure codes'].split(';') if code)), ('codon_status', Counter(row['Codon translation status'] for row in rows))]:
        category_rows.extend({'Dimension':dimension,'Value':key,'Count':value} for key,value in values.most_common())
    write_tsv(f'{prefix}.category_summary.tsv', category_rows, ['Dimension','Value','Count'])
    failures = Counter(code for row in rows for code in row['Failure codes'].split(';') if code)
    with Path(f'{prefix}.summary.txt').open('w') as handle:
        handle.write(f'Validation rows: {len(rows)}\nValidated rows: {sum(row["Overall status"]=="VARIANT_CODON_VALIDATED" for row in rows)}\nFailed rows: {sum(row["Overall status"]!="VARIANT_CODON_VALIDATED" for row in rows)}\n')
        for code,count in sorted(failures.items(), key=lambda item:(-item[1],item[0])): handle.write(f'Failure {code}: {count}\n')
    with Path(f'{prefix}.report.md').open('w') as handle:
        handle.write('# Independent genome, RNA-read and codon validation\n\n')
        handle.write('This stage checks VCF REF against GRCh38, counts ALT directly in the matching RNA BAM, and independently translates VEP reference and alternate codons. It is independent of peptide evidence.\n\n')
        handle.write(f'- Validation rows: {len(rows)}\n- Validated rows: {sum(row["Overall status"]=="VARIANT_CODON_VALIDATED" for row in rows)}\n- Failed rows: {sum(row["Overall status"]!="VARIANT_CODON_VALIDATED" for row in rows)}\n\n')
        handle.write('| Failure code | Count | Meaning |\n|---|---:|---|\n')
        meanings={'GENOME_REF_MISMATCH':'VCF REF disagrees with GRCh38.','ALT_READS_BELOW_MINIMUM':'Fewer than the configured ALT-supporting RNA reads.','ALT_FRACTION_BELOW_MINIMUM':'ALT fraction is below the configured threshold.','COMPLEX_ALLELE_NOT_COUNTED':'Allele is not a simple SNV or normalized anchored indel.','CODON_TRANSLATION_MISMATCH':'VEP codons do not translate to the stated amino acids.'}
        for code,count in sorted(failures.items(), key=lambda item:(-item[1],item[0])): handle.write(f'| {code} | {count} | {meanings.get(code,"See row-level evidence.")} |\n')
        handle.write('\nCODON_NOT_SINGLE_TRIPLET is informational for events not representable as one reference and one alternate triplet.\n')
    for suffix in ('all.tsv','validated.tsv','partial.tsv','failed.tsv','category_summary.tsv','summary.txt','report.md'): print(f'Wrote {prefix}.{suffix}')

if __name__ == '__main__':
    try: main()
    except Exception as error:
        print(f'ERROR: {error}', file=sys.stderr); raise SystemExit(1)
