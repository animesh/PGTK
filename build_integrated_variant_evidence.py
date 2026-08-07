#!/usr/bin/env python3
import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


def key(row):
    sample = row.get('Sample') or row.get('FASTA sample') or ''
    chrom = (row.get('Chromosome') or '').removeprefix('chr')
    return sample.upper(), chrom, str(row.get('Position','')), row.get('REF','').upper(), row.get('ALT','').upper()


def read_tsv(path):
    with open(path, encoding='utf-8', errors='replace', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write_tsv(path, rows, fields):
    with Path(path).open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n', extrasaction='ignore')
        writer.writeheader(); writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description='Build the strict intersection of genome/read/codon validation and sample-matched direct MS/MS evidence.')
    parser.add_argument('--variants', required=True)
    parser.add_argument('--codon-validation', required=True)
    parser.add_argument('--codon-mismatch-analysis', required=True)
    parser.add_argument('--output-prefix', default='integrated_variant_evidence')
    args = parser.parse_args()

    codon = defaultdict(list)
    for row in read_tsv(args.codon_validation): codon[key(row)].append(row)
    mismatch = defaultdict(list)
    for row in read_tsv(args.codon_mismatch_analysis): mismatch[key(row)].append(row)

    output = []
    for row in read_tsv(args.variants):
        event_key = key(row)
        codon_rows = codon.get(event_key, [])
        mismatch_rows = mismatch.get(event_key, [])
        validated_rows = [x for x in codon_rows if x.get('Overall status') == 'VARIANT_CODON_VALIDATED' and x.get('Codon strict-integration eligible') == 'yes']
        validated_consequences = len(validated_rows)
        evidence_hgvsp = {value for value in (row.get('HGVSp') or '').split(';') if value}
        matched_validated_consequences = sum(
            not evidence_hgvsp or x.get('HGVSp','') in evidence_hgvsp
            for x in validated_rows
        )
        partial_consequences = sum(x.get('Overall status') == 'VARIANT_EVIDENCE_PARTIAL' for x in codon_rows)
        failed_consequences = sum(x.get('Overall status') == 'VALIDATION_FAILED' for x in codon_rows)
        mismatch_categories = sorted({x.get('Mismatch diagnostic category','') for x in mismatch_rows if x.get('Mismatch diagnostic category')})
        direct = row.get('Primary sample-specific evidence') == 'yes'
        novel = bool(row.get('Canonical-and-reference-absent peptides'))
        search_consistent = bool(row.get('Search-consistent altered-residue peptides'))
        strict = bool(matched_validated_consequences and direct and novel and search_consistent)
        reasons = []
        if not matched_validated_consequences: reasons.append('NO_PEPTIDE_ASSOCIATED_VALIDATED_CODON_CONSEQUENCE')
        if not direct: reasons.append('NO_SAMPLE_MATCHED_DIRECT_MSMS')
        if not search_consistent: reasons.append('NO_SEARCH_CONSISTENT_ALTERED_RESIDUE_PEPTIDE')
        if not novel: reasons.append('PEPTIDE_NOT_ABSENT_FROM_BOTH_REFERENCE_SETS')
        output.append({
            **row,
            'Independent codon consequence rows': len(codon_rows),
            'Fully validated consequence rows': validated_consequences,
            'Peptide-associated fully validated consequence rows': matched_validated_consequences,
            'Partially assessed consequence rows': partial_consequences,
            'Independently failed consequence rows': failed_consequences,
            'Codon mismatch diagnostic categories': ';'.join(mismatch_categories),
            'Strict integrated evidence': 'yes' if strict else 'no',
            'Strict exclusion reasons': ';'.join(reasons),
        })

    prefix = Path(args.output_prefix)
    fields = list(output[0]) if output else ['Sample','Chromosome','Position','REF','ALT','Strict integrated evidence','Strict exclusion reasons']
    strict_rows = [row for row in output if row['Strict integrated evidence'] == 'yes']
    review_rows = [row for row in output if row['Strict integrated evidence'] != 'yes']
    write_tsv(f'{prefix}.all.tsv', output, fields)
    write_tsv(f'{prefix}.strict.tsv', strict_rows, fields)
    write_tsv(f'{prefix}.excluded.tsv', review_rows, fields)

    exclusions = Counter(code for row in review_rows for code in row['Strict exclusion reasons'].split(';') if code)
    with Path(f'{prefix}.report.md').open('w', encoding='utf-8') as handle:
        handle.write('# Strict integrated variant evidence\n\n')
        handle.write('The strict set requires an independently validated genome/read/codon consequence, sample-matched direct MS/MS evidence, a search-consistent altered-residue peptide, and peptide absence from both the canonical search proteome and Ensembl reference proteins.\n\n')
        handle.write(f'- Proteogenomic variant events evaluated: {len(output)}\n')
        handle.write(f'- Strict integrated events: {len(strict_rows)}\n')
        handle.write(f'- Events excluded from the strict set: {len(review_rows)}\n\n')
        handle.write('| Exclusion reason | Events |\n|---|---:|\n')
        for code, count in exclusions.most_common():
            handle.write(f'| {code} | {count} |\n')
        handle.write('\nThe detailed TSV retains every event and its exclusion reasons. Categories can overlap.\n')
    print(f'Wrote {prefix}.all.tsv')
    print(f'Wrote {prefix}.strict.tsv')
    print(f'Wrote {prefix}.excluded.tsv')
    print(f'Wrote {prefix}.report.md')

if __name__ == '__main__':
    main()
