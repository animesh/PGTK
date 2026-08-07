#!/usr/bin/env python3
import argparse
import csv
from collections import Counter
from pathlib import Path


def merge_tsv(paths, output):
    fields = None
    count = 0
    with Path(output).open('w', encoding='utf-8', newline='') as out:
        writer = None
        for path in sorted(map(Path, paths), key=lambda p: p.name):
            with path.open(encoding='utf-8', errors='replace', newline='') as handle:
                reader = csv.DictReader(handle, delimiter='\t')
                if fields is None:
                    fields = reader.fieldnames or []
                    writer = csv.DictWriter(out, fieldnames=fields, delimiter='\t', lineterminator='\n', extrasaction='ignore')
                    writer.writeheader()
                elif (reader.fieldnames or []) != fields:
                    raise ValueError(f'Header mismatch in {path}')
                for row in reader:
                    writer.writerow(row)
                    count += 1
        if writer is None:
            out.write('')
    return count


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=['codon', 'provenance'], required=True)
    parser.add_argument('--inputs', nargs='+', required=True)
    parser.add_argument('--output-prefix', required=True)
    args = parser.parse_args()
    prefix = Path(args.output_prefix)
    by_suffix = {}
    for path in map(Path, args.inputs):
        for suffix in ('all.tsv','validated.tsv','partial.tsv','failed.tsv','category_summary.tsv','supporting_reads.tsv','summary.txt','report.md'):
            if path.name.endswith(suffix):
                by_suffix.setdefault(suffix, []).append(path)
                break
    if args.mode == 'codon':
        counts = {}
        for suffix in ('all.tsv','validated.tsv','partial.tsv','failed.tsv'):
            counts[suffix] = merge_tsv(by_suffix.get(suffix, []), f'{prefix}.{suffix}')
        failure_counts = Counter()
        for path in by_suffix.get('failed.tsv', []):
            with path.open(encoding='utf-8', errors='replace') as handle:
                for row in csv.DictReader(handle, delimiter='\t'):
                    for code in (row.get('Failure codes') or '').split(';'):
                        if code: failure_counts[code] += 1
        primary_counts = Counter()
        status_counts = Counter()
        codon_counts = Counter()
        with Path(f'{prefix}.all.tsv').open(encoding='utf-8', errors='replace') as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                status_counts[row.get('Overall status','')] += 1
                if row.get('Primary failure code'): primary_counts[row['Primary failure code']] += 1
                codon_counts[row.get('Codon translation status','')] += 1
        if sum(primary_counts.values()) != counts['failed.tsv']:
            raise ValueError('Primary failure categories must equal failed rows')
        with Path(f'{prefix}.category_summary.tsv').open('w', newline='') as handle:
            writer=csv.DictWriter(handle, fieldnames=['Dimension','Value','Count'], delimiter='\t', lineterminator='\n'); writer.writeheader()
            for dimension,counter in [('overall_status',status_counts),('primary_failure',primary_counts),('failure_occurrence',failure_counts),('codon_status',codon_counts)]:
                for key,value in counter.most_common(): writer.writerow({'Dimension':dimension,'Value':key,'Count':value})
        with Path(f'{prefix}.summary.txt').open('w') as handle:
            handle.write(f'Validation rows: {counts["all.tsv"]}\nValidated rows: {counts["validated.tsv"]}\nPartial rows: {counts["partial.tsv"]}\nFailed rows: {counts["failed.tsv"]}\n')
            for code, count in sorted(failure_counts.items(), key=lambda x: (-x[1], x[0])):
                handle.write(f'Failure {code}: {count}\n')
        with Path(f'{prefix}.report.md').open('w') as handle:
            handle.write('# Independent genome, RNA-read and codon validation\n\n')
            handle.write(f'- Validation rows: {counts["all.tsv"]}\n- Fully validated rows: {counts["validated.tsv"]}\n- Partial or non-comparable rows: {counts["partial.tsv"]}\n- Failed rows: {counts["failed.tsv"]}\n\n')
            handle.write('| Failure code | Count |\n|---|---:|\n')
            for code, count in sorted(failure_counts.items(), key=lambda x: (-x[1], x[0])):
                handle.write(f'| {code} | {count} |\n')
    else:
        count = merge_tsv(by_suffix.get('supporting_reads.tsv', []), f'{prefix}.supporting_reads.tsv')
        variants = set()
        with Path(f'{prefix}.supporting_reads.tsv').open(encoding='utf-8', errors='replace') as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                variants.add((row.get('Sample',''), row.get('Variant','')))
        with Path(f'{prefix}.summary.txt').open('w') as handle:
            handle.write(f'ALT-supporting read alignments: {count}\nVariants with retained ALT-supporting reads: {len(variants)}\n')
        with Path(f'{prefix}.report.md').open('w') as handle:
            handle.write('# Variant-supporting RNA reads and mapping provenance\n\n')
            handle.write(f'- ALT-supporting read alignments: {count}\n- Variants with retained supporting reads: {len(variants)}\n\n')
            handle.write('The row-level TSV records sample, SRA, source FASTQ mate, read name, one-based read position, allele, base and mapping quality, CIGAR, strand and genomic alignment interval.\n')

if __name__ == '__main__':
    main()
