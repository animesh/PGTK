#!/usr/bin/env python3
import argparse
import csv
import gzip
import re
import pysam
import sys
from pathlib import Path


def open_text(path):
    path = Path(path)
    return gzip.open(path, 'rt', encoding='utf-8', errors='replace') if path.suffix == '.gz' else path.open('rt', encoding='utf-8', errors='replace')


def load_samples(path):
    result = {}
    with open_text(path) as handle:
        for row in csv.DictReader(handle):
            sample = (row.get('sample') or '').strip()
            if sample:
                result[sample] = {'SRA': (row.get('srr') or '').strip()}
    return result


def parse_bams(items):
    result = {}
    for item in items:
        if '=' not in item:
            raise ValueError(f'Invalid --bam value {item!r}; expected SAMPLE=FILE')
        sample, path = item.split('=', 1)
        result[sample] = Path(path)
    return result


def sample_for_file(path, samples):
    matches = [sample for sample in samples if re.search(rf'(?<![A-Za-z0-9]){re.escape(sample)}(?![A-Za-z0-9])', Path(path).name, re.I)]
    if len(matches) != 1:
        raise ValueError(f'Cannot map {path} uniquely to sample: {matches}')
    return matches[0]


def cigar_ops(cigar):
    return [(int(length), operation) for length, operation in re.findall(r'(\d+)([MIDNSHP=X])', cigar)]


def observe(fields, position, ref, alt):
    flag = int(fields[1]); start = int(fields[3]); mapq = int(fields[4])
    cigar = fields[5]; sequence = fields[9]; qualities = fields[10]
    operations = cigar_ops(cigar)
    ref_pos = start; read_pos = 0
    observed = ''; read_site = ''; base_quality = ''
    for index, (length, operation) in enumerate(operations):
        if operation in 'M=X':
            if ref_pos <= position < ref_pos + length:
                offset = position - ref_pos
                anchor = read_pos + offset
                read_site = anchor + 1
                if len(ref) == 1 and len(alt) == 1:
                    observed = sequence[anchor].upper()
                    if qualities != '*': base_quality = ord(qualities[anchor]) - 33
                elif len(alt) > len(ref) and alt.startswith(ref) and index + 1 < len(operations) and operations[index + 1][1] == 'I':
                    insertion_length = operations[index + 1][0]
                    observed = ref[0].upper() + sequence[anchor + 1:anchor + 1 + insertion_length].upper()
                break
            ref_pos += length; read_pos += length
        elif operation == 'I':
            read_pos += length
        elif operation in 'DN':
            if len(ref) > len(alt) and ref.startswith(alt) and position == ref_pos - 1 and length == len(ref) - len(alt):
                observed = alt.upper(); read_site = max(1, read_pos)
            ref_pos += length
        elif operation == 'S':
            read_pos += length
    alignment_end = start + sum(length for length, operation in operations if operation in 'MDN=X') - 1
    return {
        'Read name': fields[0],
        'Mate': 'R1' if flag & 64 else 'R2' if flag & 128 else 'unpaired',
        'Read strand': '-' if flag & 16 else '+',
        'Position in read (1-based)': read_site,
        'Observed allele': observed,
        'Base quality': base_quality,
        'Mapping quality': mapq,
        'CIGAR': cigar,
        'Genome alignment start': start,
        'Genome alignment end': alignment_end,
        'Genome mapping': f'{fields[2]}:{start}-{alignment_end}({"-" if flag & 16 else "+"})',
        'SAM flag': flag,
    }



def resolve_contig(alignment, chrom):
    names = set(alignment.references)
    plain = chrom[3:] if chrom.startswith('chr') else chrom
    resolved = next((name for name in (chrom, plain, 'chr' + plain) if name in names), None)
    if resolved is None:
        raise ValueError(f'contig not found in BAM: {chrom}')
    return resolved


def read_fields(read, alignment):
    return [
        read.query_name or '', str(read.flag), alignment.get_reference_name(read.reference_id),
        str(read.reference_start + 1), str(read.mapping_quality), read.cigarstring or '*',
        '*', '0', '0', read.query_sequence or '*', read.qual or '*',
    ]

def main():
    parser = argparse.ArgumentParser(description='List every RNA alignment supporting each VCF ALT allele, including read and genomic coordinates.')
    parser.add_argument('--vcf', nargs='+', required=True)
    parser.add_argument('--bam', nargs='+', required=True, metavar='SAMPLE=FILE')
    parser.add_argument('--samples', required=True)
    parser.add_argument('--min-base-quality', type=int, default=20)
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--min-mapping-quality', type=int, default=20)
    parser.add_argument('--output-prefix', default='variant_read_provenance')
    args = parser.parse_args()
    samples = load_samples(args.samples); bams = parse_bams(args.bam)
    output = []
    for vcf in args.vcf:
        sample = sample_for_file(vcf, bams)
        variants = []
        with open_text(vcf) as handle:
            for line in handle:
                if line.startswith('#'): continue
                fields = line.rstrip('\n').split('\t')
                if len(fields) >= 5:
                    chrom, position, variant_id, ref, alts = fields[:5]
                    for alt in alts.split(','):
                        variants.append((chrom, int(position), variant_id, ref, alt))
        variants_by_chrom = {}
        for variant in variants:
            variants_by_chrom.setdefault(variant[0], {}).setdefault(variant[1], []).append(variant)
        with pysam.AlignmentFile(str(bams[sample]), 'rb', threads=max(1, args.threads)) as alignment:
            for chrom, position_map in variants_by_chrom.items():
                bam_chrom = resolve_contig(alignment, chrom)
                for position in sorted(position_map):
                    seen = set()
                    max_ref_length = max(len(variant[3]) for variant in position_map[position])
                    fetch_end = position - 1 + max(1, max_ref_length)
                    for read in alignment.fetch(bam_chrom, position - 1, fetch_end):
                        if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate:
                            continue
                        if read.mapping_quality < args.min_mapping_quality or not read.cigarstring:
                            continue
                        read_key = (read.query_name, read.flag, read.reference_start, read.cigarstring)
                        if read_key in seen:
                            continue
                        seen.add(read_key)
                        fields = read_fields(read, alignment)
                        for _chrom, _position, variant_id, ref, alt in position_map[position]:
                            observation = observe(fields, position, ref, alt)
                            if observation['Observed allele'] != alt.upper():
                                continue
                            if observation['Base quality'] != '' and int(observation['Base quality']) < args.min_base_quality:
                                continue
                            fastq = f'{sample}_{observation["Mate"]}.fastq.gz' if observation['Mate'] in {'R1','R2'} else f'{sample}.fastq.gz'
                            output.append({
                                'Sample': sample, 'SRA': samples.get(sample, {}).get('SRA', ''), 'Source FASTQ': fastq,
                                'Variant': f'{chrom}:{position}:{ref}>{alt}', 'Variant ID': variant_id,
                                'Chromosome': chrom, 'Variant genome position': position, 'REF': ref, 'ALT': alt,
                                **observation,
                            })
    unique = {}
    for row in output:
        key = (row['Sample'], row['Variant'], row['Read name'], row['Mate'], row['Genome alignment start'], row['CIGAR'])
        unique[key] = row
    output = sorted(unique.values(), key=lambda row: (row['Sample'], row['Chromosome'], int(row['Variant genome position']), row['Read name'], row['Mate']))
    fields = ['Sample','SRA','Source FASTQ','Variant','Variant ID','Chromosome','Variant genome position','REF','ALT','Read name','Mate','Read strand','Position in read (1-based)','Observed allele','Base quality','Mapping quality','CIGAR','Genome alignment start','Genome alignment end','Genome mapping','SAM flag']
    prefix = Path(args.output_prefix)
    table = Path(f'{prefix}.supporting_reads.tsv'); summary = Path(f'{prefix}.summary.txt'); report = Path(f'{prefix}.report.md')
    with table.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n'); writer.writeheader(); writer.writerows(output)
    variants_supported = len({(row['Sample'], row['Variant']) for row in output})
    with summary.open('w', encoding='utf-8') as handle:
        handle.write(f'ALT-supporting read alignments: {len(output)}\nVariants with retained ALT-supporting reads: {variants_supported}\n')
    with report.open('w', encoding='utf-8') as handle:
        handle.write('# Variant-supporting RNA reads and mapping provenance\n\n')
        handle.write(f'- ALT-supporting read alignments: {len(output)}\n- Variants with retained supporting reads: {variants_supported}\n\n')
        handle.write('The row-level TSV identifies the sample, SRA accession, concatenated source FASTQ mate, original read name, one-based position of the allele in the sequenced read, strand, base quality, mapping quality, CIGAR string, and complete genomic alignment interval. Secondary, supplementary, duplicate, unmapped, low-MAPQ and low-base-quality alignments are excluded.\n\n')
        handle.write('The pipeline concatenates reads into sample-level FASTQs. Therefore, `Source FASTQ` identifies the sample R1/R2 file used by STAR; `SRA` and the preserved original read name provide provenance back to the imported run.\n')
    for path in (table, summary, report): print(f'Wrote {path}')


if __name__ == '__main__':
    try: main()
    except Exception as error:
        print(f'ERROR: {error}', file=sys.stderr); raise SystemExit(1)
