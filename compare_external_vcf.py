#!/usr/bin/env python3
import argparse
import csv
import gzip
from pathlib import Path


def open_text(path):
    return gzip.open(path, 'rt', encoding='utf-8', errors='replace') if str(path).endswith('.gz') else open(path, encoding='utf-8', errors='replace')


def norm_chrom(chrom):
    value = chrom[3:] if chrom.lower().startswith('chr') else chrom
    return 'MT' if value in {'M', 'MT'} else value


def trim_allele(pos, ref, alt):
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


def read_vcf(path):
    alleles = {}
    with open_text(path) as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            fields = line.rstrip().split('\t')
            if len(fields) < 5:
                continue
            chrom, pos, ref = norm_chrom(fields[0]), int(fields[1]), fields[3].upper()
            alts = fields[4].upper().split(',')
            gt = ''
            gt_indices = []
            if len(fields) > 9:
                fmt = fields[8].split(':')
                values = fields[9].split(':')
                sample_data = dict(zip(fmt, values))
                gt = sample_data.get('GT', '')
                gt_indices = [int(x) for x in gt.replace('|', '/').split('/') if x.isdigit()]
            for index, alt in enumerate(alts, 1):
                if alt in {'.', '*', '<NON_REF>'}:
                    continue
                npos, nref, nalt = trim_allele(pos, ref, alt)
                key = (chrom, npos, nref, nalt)
                alleles[key] = {'GT': gt, 'ALT dosage': gt_indices.count(index)}
    return alleles


def kind(key):
    return 'SNP' if len(key[2]) == len(key[3]) == 1 else 'INDEL'


def write_set(path, keys, left, right):
    fields = ['CHROM', 'POS', 'REF', 'ALT', 'Type', 'PGTK GT', 'External GT', 'PGTK ALT dosage', 'External ALT dosage', 'Genotype concordant']
    with Path(path).open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for key in sorted(keys):
            a, b = left.get(key, {}), right.get(key, {})
            comparable = a.get('GT', '') not in {'', '.', './.', '.|.'} and b.get('GT', '') not in {'', '.', './.', '.|.'}
            concordant = comparable and a.get('ALT dosage') == b.get('ALT dosage')
            writer.writerow({
                'CHROM': key[0], 'POS': key[1], 'REF': key[2], 'ALT': key[3], 'Type': kind(key),
                'PGTK GT': a.get('GT', ''), 'External GT': b.get('GT', ''),
                'PGTK ALT dosage': a.get('ALT dosage', ''), 'External ALT dosage': b.get('ALT dosage', ''),
                'Genotype concordant': 'yes' if concordant else 'no' if comparable else 'not_comparable',
            })


def main():
    parser = argparse.ArgumentParser(description='Compare normalized PGTK alleles with an external VCF.')
    parser.add_argument('--sample', required=True)
    parser.add_argument('--stage', required=True)
    parser.add_argument('--pgtk', required=True)
    parser.add_argument('--external', required=True)
    parser.add_argument('--output-prefix', required=True)
    args = parser.parse_args()

    pgtk = read_vcf(args.pgtk)
    external = read_vcf(args.external)
    pgtk_keys, external_keys = set(pgtk), set(external)
    shared = pgtk_keys & external_keys
    union = pgtk_keys | external_keys
    comparable = [key for key in shared if pgtk[key]['GT'] not in {'', '.', './.', '.|.'} and external[key]['GT'] not in {'', '.', './.', '.|.'}]
    concordant = [key for key in comparable if pgtk[key]['ALT dosage'] == external[key]['ALT dosage']]

    metrics = [
        ('PGTK alleles', len(pgtk_keys)),
        ('External alleles', len(external_keys)),
        ('Shared', len(shared)),
        ('PGTK-only', len(pgtk_keys - external_keys)),
        ('External-only', len(external_keys - pgtk_keys)),
        ('Shared SNPs', sum(kind(key) == 'SNP' for key in shared)),
        ('Shared indels', sum(kind(key) == 'INDEL' for key in shared)),
        ('Jaccard %', f'{100 * len(shared) / len(union):.6f}' if union else '100.000000'),
        ('PGTK overlap %', f'{100 * len(shared) / len(pgtk_keys):.6f}' if pgtk_keys else '100.000000'),
        ('External overlap %', f'{100 * len(shared) / len(external_keys):.6f}' if external_keys else '100.000000'),
        ('Genotypes comparable', len(comparable)),
        ('Genotypes concordant', len(concordant)),
        ('Genotype concordance %', f'{100 * len(concordant) / len(comparable):.6f}' if comparable else 'NA'),
    ]
    summary = Path(f'{args.output_prefix}.summary.tsv')
    with summary.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['Sample', 'Stage', 'Metric', 'Value'])
        for metric, value in metrics:
            writer.writerow([args.sample, args.stage, metric, value])

    write_set(f'{args.output_prefix}.shared.tsv', shared, pgtk, external)
    write_set(f'{args.output_prefix}.pgtk_only.tsv', pgtk_keys - external_keys, pgtk, external)
    write_set(f'{args.output_prefix}.external_only.tsv', external_keys - pgtk_keys, pgtk, external)

    metric_map = dict(metrics)
    Path(f'{args.output_prefix}.report.md').write_text(
        '# External VCF comparison\n\n'
        f'- Sample: {args.sample}\n- PGTK stage: {args.stage}\n'
        f'- PGTK alleles: {metric_map["PGTK alleles"]}\n- External alleles: {metric_map["External alleles"]}\n'
        f'- Shared alleles: {metric_map["Shared"]}\n- PGTK overlap: {metric_map["PGTK overlap %"]}%\n'
        f'- External overlap: {metric_map["External overlap %"]}%\n- Jaccard: {metric_map["Jaccard %"]}%\n'
        f'- Genotype concordance: {metric_map["Genotype concordance %"]}%\n\n'
        'Multiallelic records are split, chr prefixes are harmonized, and common allele prefixes and suffixes are trimmed. '
        'Reference-aware left alignment still requires normalization against the same FASTA.\n', encoding='utf-8')


if __name__ == '__main__':
    main()
