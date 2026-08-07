#!/usr/bin/env python3
import argparse
import csv
import math
import re
from collections import defaultdict
from pathlib import Path

UNITS = {'B':1, 'KB':1e3, 'MB':1e6, 'GB':1e9, 'TB':1e12}
TIME_UNITS = {'ms':0.001, 's':1, 'm':60, 'h':3600, 'd':86400}


def parse_bytes(value):
    match = re.fullmatch(r'\s*([0-9.]+)\s*([KMGT]?B)\s*', value or '', re.I)
    return float(match.group(1)) * UNITS[match.group(2).upper()] if match else 0.0


def parse_seconds(value):
    return sum(float(number) * TIME_UNITS[unit] for number, unit in re.findall(r'([0-9.]+)\s*(ms|d|h|m|s)', value or ''))


def parse_percent(value):
    text = str(value or '').strip().rstrip('%')
    if text in {'', '-', 'NA', 'N/A'}:
        return 0.0
    try:
        return float(text)
    except ValueError:
        return 0.0
def process_name(value):
    return (value or '').split(' (', 1)[0]


def main():
    parser = argparse.ArgumentParser(description='Summarize Nextflow trace resource efficiency and flag calibration problems.')
    parser.add_argument('--trace', required=True)
    parser.add_argument('--output-prefix', required=True)
    args = parser.parse_args()
    grouped = defaultdict(list)
    with Path(args.trace).open(encoding='utf-8', errors='replace', newline='') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            grouped[process_name(row.get('name'))].append(row)
    summary = []
    warnings = []
    for process, rows in sorted(grouped.items()):
        cpu_values = [parse_percent(row.get('%cpu')) for row in rows]
        rss_values = [parse_bytes(row.get('peak_rss', '')) for row in rows]
        runtimes = [parse_seconds(row.get('realtime', '')) for row in rows]
        summary.append({
            'Process': process,
            'Tasks': len(rows),
            'Median observed CPU cores': f'{sorted(cpu_values)[len(cpu_values)//2] / 100:.3f}',
            'Maximum peak RSS GB': f'{max(rss_values, default=0) / 1e9:.3f}',
            'Maximum runtime minutes': f'{max(runtimes, default=0) / 60:.3f}',
            'Failed tasks': sum(row.get('status') not in {'COMPLETED','CACHED'} for row in rows),
        })
        if max(cpu_values, default=0) > 3200:
            warnings.append({'Process':process, 'Warning':'Observed CPU exceeded 32 cores', 'Observed':f'{max(cpu_values)/100:.2f} cores'})
        if max(rss_values, default=0) > 0 and max(rss_values) > 48e9:
            warnings.append({'Process':process, 'Warning':'Peak RSS exceeded 48 GB', 'Observed':f'{max(rss_values)/1e9:.2f} GB'})
        if max(runtimes, default=0) > 12 * 3600:
            warnings.append({'Process':process, 'Warning':'Runtime exceeded 12 hours', 'Observed':f'{max(runtimes)/3600:.2f} h'})
    prefix = Path(args.output_prefix)
    for suffix, rows, fields in [
        ('summary.tsv', summary, ['Process','Tasks','Median observed CPU cores','Maximum peak RSS GB','Maximum runtime minutes','Failed tasks']),
        ('warnings.tsv', warnings, ['Process','Warning','Observed']),
    ]:
        with Path(f'{prefix}.{suffix}').open('w', encoding='utf-8', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
            writer.writeheader(); writer.writerows(rows)
    with Path(f'{prefix}.report.md').open('w', encoding='utf-8') as handle:
        handle.write('# Pipeline resource-efficiency report\n\n')
        handle.write(f'- Processes observed: {len(summary)}\n- Warnings: {len(warnings)}\n\n')
        handle.write('This report is generated from the completed Nextflow trace. Compare multiple full runs before reducing allocations further.\n')
    print(f'Wrote {prefix}.summary.tsv')
    print(f'Wrote {prefix}.warnings.tsv')
    print(f'Wrote {prefix}.report.md')


if __name__ == '__main__':
    main()
