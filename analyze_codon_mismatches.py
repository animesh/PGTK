#!/usr/bin/env python3
import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

STANDARD = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S','TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T','AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A','GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}
MITO = dict(STANDARD, ATA='M', TGA='W', AGA='*', AGG='*')
COMPLEMENT = str.maketrans('ACGT', 'TGCA')


def clean_codon(value):
    return re.sub('[^ACGT]', '', (value or '').upper())


def norm_aa(value):
    value = (value or '').strip().upper()
    return {'X':'*', 'TER':'*', 'STOP':'*'}.get(value, value)


def reverse_complement(value):
    return clean_codon(value).translate(COMPLEMENT)[::-1]


def mutation_orientation(ref, alt, ref_codon, alt_codon):
    if len(ref) != 1 or len(alt) != 1 or len(ref_codon) != 3 or len(alt_codon) != 3:
        return 'NON_SNV_OR_NON_TRIPLET'
    changes = [(a, b) for a, b in zip(ref_codon, alt_codon) if a != b]
    if len(changes) != 1:
        return f'CODON_DIFFERENCES_{len(changes)}'
    observed_ref, observed_alt = changes[0]
    if (observed_ref, observed_alt) == (ref.upper(), alt.upper()):
        return 'PLUS_STRAND_ALLELE_CONSISTENT'
    if (observed_ref, observed_alt) == (reverse_complement(ref), reverse_complement(alt)):
        return 'MINUS_STRAND_ALLELE_CONSISTENT'
    return 'CODON_ALLELE_INCONSISTENT'


def diagnosis(row):
    chrom = (row.get('Chromosome') or '').replace('chr', '').upper()
    ref_codon = clean_codon(row.get('Reference codon'))
    alt_codon = clean_codon(row.get('Alternate codon'))
    vep_ref = norm_aa(row.get('VEP reference amino acid'))
    vep_alt = norm_aa(row.get('VEP alternate amino acid'))
    standard_ref = STANDARD.get(ref_codon, '')
    standard_alt = STANDARD.get(alt_codon, '')
    mito_ref = MITO.get(ref_codon, '')
    mito_alt = MITO.get(alt_codon, '')
    orientation = mutation_orientation(row.get('REF',''), row.get('ALT',''), ref_codon, alt_codon)
    consequence = row.get('Consequence','')

    if len(ref_codon) != 3 or len(alt_codon) != 3:
        category = 'NON_SINGLE_TRIPLET_NOT_COMPARABLE'
    elif len(vep_ref) != 1 or len(vep_alt) != 1:
        category = 'VEP_AMINO_ACID_NOT_SINGLE_SYMBOL'
    elif chrom in {'M', 'MT'} and (mito_ref, mito_alt) == (vep_ref, vep_alt):
        category = 'MITOCHONDRIAL_CODE_RESOLVES_MISMATCH'
    elif vep_ref == 'U' and ref_codon == 'TGA':
        category = 'SELENOCYSTEINE_RECODING'
    elif (standard_ref, standard_alt) == (vep_ref, vep_alt):
        category = 'NO_REPRODUCIBLE_MISMATCH'
    elif standard_ref != vep_ref and standard_alt == vep_alt:
        category = 'REFERENCE_AMINO_ACID_MISMATCH_ONLY'
    elif standard_ref == vep_ref and standard_alt != vep_alt:
        category = 'ALTERNATE_AMINO_ACID_MISMATCH_ONLY'
    elif 'frameshift_variant' in consequence:
        category = 'FRAMESHIFT_NOT_SINGLE_CODON_COMPARABLE'
    elif orientation == 'CODON_ALLELE_INCONSISTENT':
        category = 'CODON_CHANGE_INCONSISTENT_WITH_VCF_ALLELES'
    else:
        category = 'UNEXPLAINED_TRANSLATION_MISMATCH'
    return category, orientation, standard_ref, standard_alt, mito_ref, mito_alt


def main():
    parser = argparse.ArgumentParser(description='Diagnose each codon-translation mismatch without changing the source validation evidence.')
    parser.add_argument('--input', required=True)
    parser.add_argument('--output-prefix', default='codon_mismatch_analysis')
    args = parser.parse_args()

    detailed = []
    with open(args.input, encoding='utf-8', errors='replace', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        for row in reader:
            if row.get('Codon translation status') != 'CODON_TRANSLATION_MISMATCH':
                continue
            category, orientation, standard_ref, standard_alt, mito_ref, mito_alt = diagnosis(row)
            detailed.append({
                **row,
                'Mismatch diagnostic category': category,
                'Codon versus VCF allele orientation': orientation,
                'Recomputed standard reference AA': standard_ref,
                'Recomputed standard alternate AA': standard_alt,
                'Recomputed mitochondrial reference AA': mito_ref,
                'Recomputed mitochondrial alternate AA': mito_alt,
                'Needs manual review': 'yes' if category in {
                    'CODON_CHANGE_INCONSISTENT_WITH_VCF_ALLELES',
                    'REFERENCE_AMINO_ACID_MISMATCH_ONLY',
                    'ALTERNATE_AMINO_ACID_MISMATCH_ONLY',
                    'UNEXPLAINED_TRANSLATION_MISMATCH',
                } else 'no',
            })

    prefix = Path(args.output_prefix)
    fields = list(detailed[0]) if detailed else [
        'Sample','Chromosome','Position','REF','ALT','Transcript','Consequence',
        'Mismatch diagnostic category','Codon versus VCF allele orientation','Needs manual review'
    ]
    with Path(f'{prefix}.detailed.tsv').open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n', extrasaction='ignore')
        writer.writeheader(); writer.writerows(detailed)

    categories = Counter(row['Mismatch diagnostic category'] for row in detailed)
    orientations = Counter(row['Codon versus VCF allele orientation'] for row in detailed)
    by_sample = defaultdict(Counter)
    by_consequence = Counter()
    by_transcript = Counter()
    for row in detailed:
        by_sample[row.get('Sample','')][row['Mismatch diagnostic category']] += 1
        for consequence in filter(None, row.get('Consequence','').split('&')):
            by_consequence[consequence] += 1
        by_transcript[row.get('Transcript','')] += 1

    summary_fields = ['Dimension','Value','Count']
    summary_rows = []
    summary_rows += [{'Dimension':'diagnostic_category','Value':key,'Count':value} for key,value in categories.most_common()]
    summary_rows += [{'Dimension':'orientation','Value':key,'Count':value} for key,value in orientations.most_common()]
    summary_rows += [{'Dimension':'consequence','Value':key,'Count':value} for key,value in by_consequence.most_common()]
    summary_rows += [{'Dimension':'transcript','Value':key,'Count':value} for key,value in by_transcript.most_common(50)]
    for sample in sorted(by_sample):
        summary_rows += [{'Dimension':f'sample:{sample}','Value':key,'Count':value} for key,value in by_sample[sample].most_common()]
    with Path(f'{prefix}.summary.tsv').open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields, delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(summary_rows)

    manual = [row for row in detailed if row['Needs manual review'] == 'yes']
    with Path(f'{prefix}.manual_review.tsv').open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n', extrasaction='ignore')
        writer.writeheader(); writer.writerows(manual)

    with Path(f'{prefix}.report.md').open('w', encoding='utf-8') as handle:
        handle.write('# Codon-translation mismatch investigation\n\n')
        handle.write('This analysis re-evaluates every `CODON_TRANSLATION_MISMATCH` row. It preserves the original result and adds diagnostic categories; it does not silently convert failures into passes.\n\n')
        handle.write(f'- Mismatch rows investigated: {len(detailed)}\n')
        handle.write(f'- Rows requiring manual review: {len(manual)}\n')
        handle.write(f'- Rows explained or marked non-comparable: {len(detailed)-len(manual)}\n\n')
        handle.write('## Diagnostic categories\n\n| Category | Rows |\n|---|---:|\n')
        for key, value in categories.most_common():
            handle.write(f'| {key} | {value} |\n')
        handle.write('\n## Allele-orientation checks\n\n| Result | Rows |\n|---|---:|\n')
        for key, value in orientations.most_common():
            handle.write(f'| {key} | {value} |\n')
        handle.write('\n## Interpretation\n\n')
        handle.write('- `MITOCHONDRIAL_CODE_RESOLVES_MISMATCH` indicates that the vertebrate mitochondrial code agrees with VEP.\n')
        handle.write('- `SELENOCYSTEINE_RECODING` flags TGA/U biology that cannot be represented by the standard code alone.\n')
        handle.write('- `CODON_CHANGE_INCONSISTENT_WITH_VCF_ALLELES` indicates that the codon nucleotide change matches neither the genomic allele nor its reverse complement.\n')
        handle.write('- Reference-only, alternate-only, and unexplained mismatches remain failures pending transcript, annotation-version, and sequence reconstruction review.\n')

    print(f'Wrote {prefix}.detailed.tsv')
    print(f'Wrote {prefix}.summary.tsv')
    print(f'Wrote {prefix}.manual_review.tsv')
    print(f'Wrote {prefix}.report.md')

if __name__ == '__main__':
    main()
