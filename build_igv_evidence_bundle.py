#!/usr/bin/env python3
import argparse
import csv
import gzip
try:
    import pysam
except ImportError:
    pysam = None
import xml.etree.ElementTree as ET
from pathlib import Path


def open_text(path):
    return gzip.open(path, 'rt', encoding='utf-8', errors='replace') if str(path).endswith('.gz') else open(path, encoding='utf-8', errors='replace')


def sample_from(path):
    return Path(path).name.split('.')[0]



def choose_annotation(annotations, alt, ref):
    alt=(alt or '').upper(); ref=(ref or '').upper(); alleles={alt}
    if ref and alt.startswith(ref): alleles.add(alt[len(ref):] or '-')
    if ref and len(ref)>len(alt) and ref.startswith(alt): alleles.add('-')
    matched=[row for row in annotations if (row.get('Allele') or '').upper() in alleles]
    pool=matched or list(annotations)
    if not pool:return {}
    return min(pool,key=lambda row:(row.get('PICK')!='1',row.get('CANONICAL')!='YES',not bool(row.get('HGVSp') or row.get('Amino_acids'))))

def add_event(events, event_id, sample, event_class, chrom, start0, end, label, source, chrom2='', start2='', end2='', gene='', consequence='', impact='', transcript='', protein_change=''):
    events.append({
        'Event': event_id, 'Sample': sample, 'Class': event_class,
        'Chrom': chrom, 'Start0': start0, 'End': end,
        'Chrom2': chrom2, 'Start2_0': start2, 'End2': end2,
        'Label': label, 'Source': Path(source).name,
        'Gene': gene, 'Consequence': consequence, 'Impact': impact,
        'Transcript': transcript, 'ProteinChange': protein_change,
    })


def main():
    parser = argparse.ArgumentParser(description='Build IGV manifests, coordinates, sessions and compact event BAMs.')
    parser.add_argument('--genome', required=True)
    parser.add_argument('--rna-vcf', nargs='*', default=[])
    parser.add_argument('--progression-vcf', nargs='*', default=[])
    parser.add_argument('--fusion-table', nargs='*', default=[])
    parser.add_argument('--splice-table', nargs='*', default=[])
    parser.add_argument('--bam', action='append', default=[])
    parser.add_argument('--padding', type=int, default=100)
    parser.add_argument('--output-prefix', required=True)
    args = parser.parse_args()

    bams = dict(value.split('=', 1) for value in args.bam)
    if bams and pysam is None: raise RuntimeError('pysam is required when BAM inputs are supplied')
    events = []
    event_number = 0

    for event_class, paths in [('rna_variant', args.rna_vcf), ('progression_variant', args.progression_vcf)]:
        for path in paths:
            sample = sample_from(path)
            csq_fields = []
            with open_text(path) as handle:
                for line in handle:
                    if line.startswith('##INFO=<ID=CSQ'):
                        marker = 'Format: '
                        if marker in line:
                            csq_fields = line.split(marker, 1)[1].split('\">', 1)[0].rstrip('\n\r"').split('|')
                        continue
                    if line.startswith('#'):
                        continue
                    fields = line.rstrip().split('\t')
                    if len(fields) < 8:
                        continue
                    info = dict(item.split('=',1) if '=' in item else (item,'') for item in fields[7].split(';'))
                    annotations = [dict(zip(csq_fields, value.split('|'))) for value in info.get('CSQ','').split(',') if value] if csq_fields else []
                    canonical = choose_annotation(annotations, fields[4].split(',')[0], fields[3])
                    event_number += 1
                    start0 = max(0, int(fields[1]) - 1)
                    add_event(events, f'E{event_number:08d}', sample, event_class, fields[0], start0,
                              start0 + max(1, len(fields[3])), f'{fields[3]}>{fields[4]}', path,
                              gene=canonical.get('SYMBOL') or canonical.get('Gene') or '',
                              consequence=canonical.get('Consequence',''), impact=canonical.get('IMPACT',''),
                              transcript=canonical.get('Feature',''),
                              protein_change=canonical.get('HGVSp') or canonical.get('Amino_acids') or '')

    for path in args.fusion_table:
        sample = sample_from(path)
        with open(path, encoding='utf-8', errors='replace', newline='') as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                bp1, bp2 = row.get('breakpoint1', ''), row.get('breakpoint2', '')
                if ':' not in bp1 or ':' not in bp2:
                    continue
                chrom1, pos1 = bp1.rsplit(':', 1)
                chrom2, pos2 = bp2.rsplit(':', 1)
                event_number += 1
                gene1 = row.get('#gene1', row.get('gene1', ''))
                gene2 = row.get('gene2', '')
                add_event(events, f'E{event_number:08d}', sample, 'fusion', chrom1, max(0, int(pos1) - 1), int(pos1),
                          f'{gene1}--{gene2}', path, chrom2, max(0, int(pos2) - 1), int(pos2))

    for path in args.splice_table:
        sample = sample_from(path)
        with open(path, encoding='utf-8', errors='replace', newline='') as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                if row.get('Status') not in {'', 'RNA_VALIDATED'}:
                    continue
                transcript = row.get('Event') or row.get('transcript_id', '')
                junction_values = (row.get('Junctions') or '').split(';')
                for junction in filter(None, junction_values):
                    if ':' not in junction or '-' not in junction:
                        continue
                    chrom, coordinates = junction.rsplit(':', 1)
                    start, end = coordinates.split('-', 1)
                    event_number += 1
                    add_event(events, f'E{event_number:08d}', sample, 'splice_junction', chrom,
                              max(0, int(start) - 1), int(end), transcript, path)

    fields = ['Event', 'Sample', 'Class', 'Chrom', 'Start0', 'End', 'Chrom2', 'Start2_0', 'End2', 'Label', 'Source', 'Gene', 'Consequence', 'Impact', 'Transcript', 'ProteinChange']
    with open(args.output_prefix + '.events.tsv', 'w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(events)

    with open(args.output_prefix + '.events.bed', 'w', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        for event in events:
            writer.writerow([event['Chrom'], event['Start0'], event['End'], event['Event'], 0, '.'])
            if event['Chrom2']:
                writer.writerow([event['Chrom2'], event['Start2_0'], event['End2'], event['Event'] + '|B2', 0, '.'])

    with open(args.output_prefix + '.events.bedpe', 'w', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        for event in events:
            if event['Chrom2']:
                writer.writerow([event['Chrom'], event['Start0'], event['End'], event['Chrom2'], event['Start2_0'], event['End2'], event['Event'], 0, '.', '.'])

    regions = {sample: [] for sample in bams}
    for event in events:
        if event['Sample'] not in regions:
            continue
        regions[event['Sample']].append((event['Chrom'], max(0, int(event['Start0']) - args.padding), int(event['End']) + args.padding))
        if event['Chrom2']:
            regions[event['Sample']].append((event['Chrom2'], max(0, int(event['Start2_0']) - args.padding), int(event['End2']) + args.padding))

    manifest = []
    for sample, bam in sorted(bams.items()):
        output_bam = f'{args.output_prefix}.{sample}.events.bam'
        merged = []
        for chrom, start, end in sorted(set(regions.get(sample, []))):
            if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
                merged[-1] = (chrom, merged[-1][1], max(merged[-1][2], end))
            else:
                merged.append((chrom, start, end))
        region_bed = f'{args.output_prefix}.{sample}.regions.bed'
        with open(region_bed, 'w', encoding='utf-8') as handle:
            for chrom, start, end in merged:
                handle.write(f'{chrom}\t{start}\t{end}\n')
        temporary_bam=output_bam+'.unsorted.bam'
        try:
            with pysam.AlignmentFile(bam,'rb') as source:
                available=set(source.references)
                with pysam.AlignmentFile(temporary_bam,'wb',template=source) as output:
                    seen=set()
                    for chrom,start,end in merged:
                        plain=chrom[3:] if chrom.startswith('chr') else chrom
                        resolved=next((name for name in (chrom,plain,'chr'+plain) if name in available),None)
                        if resolved is None:continue
                        for read in source.fetch(resolved,start,end):
                            marker=(read.query_name,read.flag,read.reference_id,read.reference_start,read.cigarstring)
                            if marker not in seen:output.write(read);seen.add(marker)
            pysam.sort('-o',output_bam,temporary_bam);pysam.index(output_bam)
            if pysam.quickcheck(output_bam)!='':raise RuntimeError(f'HTSlib quickcheck failed: {output_bam}')
        finally:Path(temporary_bam).unlink(missing_ok=True)
        manifest.append([sample, bam, output_bam, output_bam + '.bai', len(merged)])

    with open(args.output_prefix + '.sample_manifest.tsv', 'w', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['Sample', 'Source BAM', 'Event BAM', 'Index', 'Regions'])
        writer.writerows(manifest)

    with open(args.output_prefix + '.igv.batch.txt', 'w', encoding='utf-8') as handle:
        handle.write('new\n')
        handle.write(f'genome {Path(args.genome).resolve()}\n')
        for _, _, output_bam, _, _ in manifest:
            handle.write(f'load {output_bam}\n')
        handle.write(f'load {args.output_prefix}.events.bed\n')
        handle.write(f'load {args.output_prefix}.events.bedpe\n')
        for event in events[:1000]:
            handle.write(f"goto {event['Chrom']}:{int(event['Start0']) + 1}-{event['End']}\n")
            if event['Chrom2']:
                handle.write(f"goto {event['Chrom2']}:{int(event['Start2_0']) + 1}-{event['End2']}\n")
        handle.write('exit\n')

    root = ET.Element('Session', genome=str(Path(args.genome).resolve()), version='8')
    resources = ET.SubElement(root, 'Resources')
    for _, _, output_bam, _, _ in manifest:
        ET.SubElement(resources, 'Resource', path=output_bam)
    ET.SubElement(resources, 'Resource', path=f'{args.output_prefix}.events.bed')
    ET.SubElement(resources, 'Resource', path=f'{args.output_prefix}.events.bedpe')
    ET.ElementTree(root).write(args.output_prefix + '.igv.session.xml', encoding='utf-8', xml_declaration=True)

    counts = {event_class: sum(event['Class'] == event_class for event in events)
              for event_class in ['rna_variant', 'progression_variant', 'fusion', 'splice_junction']}
    Path(args.output_prefix + '.summary.txt').write_text(
        f'Events: {len(events)}\nSamples: {len(bams)}\n'
        f'RNA variants: {counts["rna_variant"]}\nProgression variants: {counts["progression_variant"]}\n'
        f'Fusions: {counts["fusion"]}\nSplice junctions: {counts["splice_junction"]}\n', encoding='utf-8')


if __name__ == '__main__':
    main()
