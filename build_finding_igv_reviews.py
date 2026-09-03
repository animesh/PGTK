#!/usr/bin/env python3
import argparse
import bisect
import csv
import gzip
import hashlib
import re
from collections import defaultdict
from pathlib import Path
import xml.etree.ElementTree as ET

from variant_read_evidence import classify_sam_fields, evidence_status

try:
    import pysam
except ImportError:
    pysam = None

SUMMARY_FIELDS = [
    'EventID','EvidenceClasses','SourceEvents','Sources','Sample','Gene','Consequence','Impact',
    'PredictedConsequence','PredictedImpact','Transcript','ProteinChange','Chrom','Position','REF','ALT',
    'ReadValidationStatus','ValidationExplanation','CountUnit','UniqueAlignments','CallableAlignments',
    'ExactAltReads','CleanReferenceReads','ExcludedReads','AltFractionAmongClean','CallableFractionAmongExamined'
]
READ_FIELDS = ['EventID','Sample','ReadName','Class','Reason','Contig','Start','MapQ','Flag','CIGAR','Observed','Quality']
DISPLAY_FIELDS = ['EventID','Sample','Category','AlignmentKey','ReadName','Contig','Start0','Flag','CIGAR','SequenceHash']
BAM_CATEGORIES = ('exact_alt_unique','exact_alt_display','reference_display','event_display')


def parse_label(label):
    match = re.search(r'([ACGTN]+)>([ACGTN]+)', (label or '').upper())
    return match.groups() if match else ('', '')


def safe_id(event):
    ref, alt = parse_label(event.get('Label', ''))
    label = re.sub(r'[^A-Za-z0-9_.-]+', '_', event.get('Label', 'event')).strip('_')
    gene = re.sub(r'[^A-Za-z0-9_.-]+', '_', event.get('Gene', '') or 'INTERGENIC').strip('_')
    return f"{gene}_{event['Sample']}_{event['Chrom']}_{int(event['Start0']) + 1}_{ref}_{alt}_{label}"


def resolve_contig(alignment, chrom):
    names = set(alignment.references)
    text = str(chrom)
    for candidate in (text, text.removeprefix('chr'), 'chr' + text.removeprefix('chr')):
        if candidate in names:
            return candidate
    raise ValueError(f'Contig {chrom!r} is absent from BAM header')


def classify(fields, target, ref, alt, min_mapq, min_baseq):
    return classify_sam_fields(fields, target, ref, alt, min_mapq, min_baseq)


def write_table(path, fields, rows):
    with Path(path).open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def make_index(items, position_index):
    grouped = defaultdict(list)
    for item in items:
        grouped[item[0]].append(item)
    result = {}
    for chrom, values in grouped.items():
        values.sort(key=lambda value: value[position_index])
        result[chrom] = ([value[position_index] for value in values], values)
    return result


def overlapping(index, chrom, start, end, pad=0):
    data = index.get(chrom)
    if not data:
        return ()
    positions, values = data
    return values[
        bisect.bisect_left(positions, start - pad):
        bisect.bisect_right(positions, end + pad)
    ]


def browser_safe_display(read, ref='', alt=''):
    if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate or read.is_qcfail:
        return False
    if not read.cigarstring or not read.query_sequence:
        return False
    return True


def browser_display_copy(read, header):
    copied = pysam.AlignedSegment.fromstring(read.to_string(), header)
    copied.flag &= ~(1 | 2 | 4 | 8 | 32 | 64 | 128 | 256 | 512 | 1024 | 2048)
    copied.next_reference_id = -1
    copied.next_reference_start = -1
    copied.template_length = 0
    return copied


def sequence_hash(sequence):
    return hashlib.sha256((sequence or '').encode('ascii', errors='ignore')).hexdigest()[:24]


def alignment_key_values(read, reference_name=None):
    contig = reference_name if reference_name is not None else read.reference_name
    values = (
        read.query_name or '', contig or '', int(read.reference_start), int(read.flag),
        read.cigarstring or '', sequence_hash(read.query_sequence),
    )
    digest = hashlib.sha256('\x1f'.join(map(str, values)).encode('utf-8')).hexdigest()
    return digest, values


def read_observes_interval(read, start0, end0):
    if end0 <= start0:
        end0 = start0 + 1
    return any(block_start < end0 and block_end > start0 for block_start, block_end in read.get_blocks())


def finalize_bam(source, target):
    pysam.sort('-o', str(target), str(source))
    source.unlink(missing_ok=True)
    Path(str(target) + '.bai').unlink(missing_ok=True)
    target.with_suffix('.bai').unlink(missing_ok=True)
    pysam.index(str(target))


def display_row(event_id, sample, category, read, contig):
    key, values = alignment_key_values(read, contig)
    read_name, contig, start0, flag, cigar, seq_hash = values
    return {
        'EventID': event_id, 'Sample': sample, 'Category': category, 'AlignmentKey': key,
        'ReadName': read_name, 'Contig': contig, 'Start0': start0, 'Flag': flag,
        'CIGAR': cigar, 'SequenceHash': seq_hash,
    }


def main():
    parser = argparse.ArgumentParser(description='Consolidated finding review with event-exact browser identities')
    parser.add_argument('--events', required=True)
    parser.add_argument('--bam', action='append', default=[])
    parser.add_argument('--genome', required=True)
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--event-id', default='')
    parser.add_argument('--padding', type=int, default=100)
    parser.add_argument('--mapq', type=int, default=20)
    parser.add_argument('--baseq', type=int, default=20)
    parser.add_argument('--reference-display-reads', type=int, default=20)
    parser.add_argument('--alt-display-reads', type=int, default=100)
    parser.add_argument('--finding-classes', default='rna_variant,progression_variant,fusion,splice_junction')
    parser.add_argument('--primary-class-order', default='rna_variant,progression_variant')
    parser.add_argument('--priority-mode', choices=('all','filter'), default='all')
    parser.add_argument('--priority-genes', default='')
    parser.add_argument('--priority-impacts', default='')
    parser.add_argument('--priority-consequences', default='')
    parser.add_argument('--priority-limit', type=int, default=0)
    parser.add_argument('--gene-filter', default='')
    parser.add_argument('--sample-filter', default='')
    parser.add_argument('--diagnostic-read-limit', type=int, default=100000)
    parser.add_argument('--excluded-read-limit', type=int, default=10000)
    parser.add_argument('--progress-every-reads', type=int, default=1000000)
    parser.add_argument('--plan-only', action='store_true')
    args = parser.parse_args()

    if pysam is None:
        raise RuntimeError('pysam is required')
    bams = dict(item.split('=', 1) for item in args.bam)
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)

    with open(args.events, encoding='utf-8', newline='') as handle:
        events = list(csv.DictReader(handle, delimiter='\t'))

    classes = {value.strip() for value in args.finding_classes.split(',') if value.strip()}
    genes = {value.strip() for value in args.gene_filter.split(',') if value.strip()}
    samples = {value.strip() for value in args.sample_filter.split(',') if value.strip()}
    order = [value.strip() for value in args.primary_class_order.split(',') if value.strip()]
    selected = [
        (safe_id(event), event) for event in events
        if event['Class'] in classes
        and (not args.event_id or event['Event'] == args.event_id)
        and (not genes or event.get('Gene', '') in genes)
        and (not samples or event.get('Sample', '') in samples)
    ]

    grouped = defaultdict(list)
    for event_id, event in selected:
        grouped[event_id].append(event)
    consolidated = []
    for event_id in sorted(grouped):
        members = grouped[event_id]
        event = next((row for preferred in order for row in members if row['Class'] == preferred), members[0]).copy()
        event['_id'] = event_id
        event['EvidenceClasses'] = ';'.join(sorted({row['Class'] for row in members}))
        event['SourceEvents'] = ';'.join(row['Event'] for row in members)
        event['Sources'] = ';'.join(sorted({row.get('Source', '') for row in members if row.get('Source', '')}))
        consolidated.append(event)

    write_table(output / 'event_consolidation.tsv', ['EventID','EvidenceClasses','SourceEvents','Sources'], consolidated)
    if args.plan_only:
        (output / 'consolidation_summary.txt').write_text(
            f'Input selected rows: {len(selected)}\nConsolidated findings: {len(consolidated)}\nMerged duplicate rows: {len(selected)-len(consolidated)}\n'
        )
        return
    if not bams:
        raise RuntimeError('BAM inputs are required')

    by_sample = defaultdict(list)
    for event in consolidated:
        if event['Sample'] in bams:
            by_sample[event['Sample']].append(event)

    counts = {event['_id']: {'unique':0,'alt':0,'ref':0,'excluded':0} for event in consolidated}
    metadata = {}
    region_rows = []
    bed = ['track name="PGTK_consolidated_findings" itemRgb="On"']
    bam_manifest = []
    stats = []

    read_handle = (output / 'read_classification.tsv').open('w', newline='')
    excluded_handle = (output / 'excluded_reads.tsv').open('w', newline='')
    read_writer = csv.DictWriter(read_handle, fieldnames=READ_FIELDS, delimiter='\t', lineterminator='\n')
    excluded_writer = csv.DictWriter(excluded_handle, fieldnames=READ_FIELDS, delimiter='\t', lineterminator='\n')
    read_writer.writeheader()
    excluded_writer.writeheader()
    display_handle = gzip.open(output / 'display_alignment_manifest.tsv.gz', 'wt', encoding='utf-8', newline='')
    display_writer = csv.DictWriter(display_handle, fieldnames=DISPLAY_FIELDS, delimiter='\t', lineterminator='\n')
    display_writer.writeheader()
    diagnostic_rows = excluded_rows = 0

    try:
        for sample, bam_path in sorted(bams.items()):
            with pysam.AlignmentFile(bam_path, 'rb') as bam:
                variants = []
                nonvariants = []
                for event in by_sample[sample]:
                    event_id = event['_id']
                    ref, alt = parse_label(event.get('Label', ''))
                    target = int(event['Start0']) + 1
                    chrom = resolve_contig(bam, event['Chrom'])
                    start = max(1, target - args.padding)
                    end = target + max(len(ref), len(alt), 1) + args.padding
                    regions = [(chrom, start, end)]
                    if event.get('Chrom2'):
                        regions.append((
                            resolve_contig(bam, event['Chrom2']),
                            max(1, int(event['Start2_0']) + 1 - args.padding),
                            int(event['End2']) + args.padding,
                        ))
                    for region_chrom, region_start, region_end in regions:
                        region_rows.append({'EventID':event_id,'Sample':sample,'Chrom':region_chrom,'Start':region_start,'End':region_end})
                    metadata[event_id] = (event, chrom, target, ref, alt)
                    if ref and alt:
                        variants.append((chrom, target, event_id, event, start, end))
                        bed.append(f'{chrom}\t{target-1}\t{target}\t{event_id}\t1000\t.\t{target-1}\t{target}\t255,80,80')
                    else:
                        for region_chrom, region_start, region_end in regions:
                            nonvariants.append((region_chrom, region_start - 1, region_end, event_id, event))
                            bed.append(f'{region_chrom}\t{region_start-1}\t{region_end}\t{event_id} {event["Class"]}\t1000\t.\t{region_start-1}\t{region_end}\t80,120,255')

                variant_index = make_index(variants, 1)
                nonvariant_index = make_index(nonvariants, 1)
                temporary = {category: output / f'{sample}.{category}.unsorted.bam' for category in BAM_CATEGORIES}
                writers = {category: pysam.AlignmentFile(path, 'wb', template=bam) for category, path in temporary.items()}
                written_global = {category:set() for category in BAM_CATEGORIES}
                written_event = set()
                event_category_counts = defaultdict(int)
                scanned = exact_unique_count = 0
                category_counts = defaultdict(int)

                try:
                    for read in bam.fetch(until_eof=True):
                        if read.is_unmapped or read.reference_id < 0:
                            continue
                        scanned += 1
                        contig = bam.get_reference_name(read.reference_id)
                        start1 = read.reference_start + 1
                        end1 = read.reference_end or start1
                        copied = browser_display_copy(read, bam.header) if browser_safe_display(read) else None
                        copied_key = alignment_key_values(copied, contig)[0] if copied is not None else ''

                        variant_candidates = overlapping(variant_index, contig, start1, end1, 2)
                        exact_for_any = False
                        for target, event_id, event, window_start, window_end in [
                            (row[1], row[2], row[3], row[4], row[5]) for row in variant_candidates
                        ]:
                            ref, alt = parse_label(event.get('Label', ''))
                            fields = read.to_string().split('\t')
                            result_class, reason, observed, quality = classify(
                                fields, target, ref, alt, args.mapq, args.baseq
                            )
                            count = counts[event_id]
                            count['unique'] += 1
                            if result_class == 'EXACT_ALT':
                                count['alt'] += 1
                                exact_for_any = True
                            elif result_class == 'CLEAN_REFERENCE':
                                count['ref'] += 1
                            else:
                                count['excluded'] += 1

                            row = dict(zip(READ_FIELDS, [
                                event_id, sample, fields[0], result_class, reason, fields[2], fields[3],
                                fields[4], fields[1], fields[5], observed, quality,
                            ]))
                            if result_class != 'EXCLUDED' and diagnostic_rows < args.diagnostic_read_limit:
                                read_writer.writerow(row)
                                diagnostic_rows += 1
                            elif result_class == 'EXCLUDED' and excluded_rows < args.excluded_read_limit:
                                excluded_writer.writerow(row)
                                excluded_rows += 1

                            if copied is None:
                                continue
                            if result_class == 'EXACT_ALT':
                                category = 'exact_alt_display'
                                limit_ok = count['alt'] <= args.alt_display_reads
                            elif result_class == 'CLEAN_REFERENCE':
                                category = 'reference_display'
                                limit_ok = count['ref'] <= args.reference_display_reads
                            else:
                                category = None
                                limit_ok = False
                            if category and limit_ok:
                                event_key = (event_id, category, copied_key)
                                if event_key not in written_event:
                                    display_writer.writerow(display_row(event_id, sample, category, copied, contig))
                                    written_event.add(event_key)
                                    event_category_counts[(event_id, category)] += 1
                                if copied_key not in written_global[category]:
                                    writers[category].write(copied)
                                    written_global[category].add(copied_key)
                                    category_counts[category] += 1

                        if exact_for_any:
                            original_key = alignment_key_values(read, contig)[0]
                            if original_key not in written_global['exact_alt_unique']:
                                writers['exact_alt_unique'].write(read)
                                written_global['exact_alt_unique'].add(original_key)
                                exact_unique_count += 1

                        if copied is not None:
                            candidates = overlapping(
                                nonvariant_index, contig, read.reference_start,
                                read.reference_end or read.reference_start
                            )
                            for region_start, region_end, event_id, _event in [
                                (row[1], row[2], row[3], row[4]) for row in candidates
                            ]:
                                if not read_observes_interval(read, region_start, region_end):
                                    continue
                                if event_category_counts[(event_id, 'event_display')] >= args.alt_display_reads:
                                    continue
                                event_key = (event_id, 'event_display', copied_key)
                                if event_key in written_event:
                                    continue
                                display_writer.writerow(display_row(event_id, sample, 'event_display', copied, contig))
                                written_event.add(event_key)
                                event_category_counts[(event_id, 'event_display')] += 1
                                if copied_key not in written_global['event_display']:
                                    writers['event_display'].write(copied)
                                    written_global['event_display'].add(copied_key)
                                    category_counts['event_display'] += 1

                        if args.progress_every_reads and scanned % args.progress_every_reads == 0:
                            print(f'PROGRESS sample={sample} reads={scanned} exact_alt={exact_unique_count}', flush=True)
                finally:
                    for writer in writers.values():
                        writer.close()

                for category in BAM_CATEGORIES:
                    final = output / f'{sample}.{category}.bam'
                    finalize_bam(temporary[category], final)
                    count = exact_unique_count if category == 'exact_alt_unique' else category_counts[category]
                    bam_manifest.append({'Sample':sample,'Category':category,'BAM':final.name,'Index':final.name+'.bai','UniqueAlignments':count})
                stats.append((sample, scanned, len(variants), len(nonvariants), exact_unique_count))
                print(f'PROGRESS sample={sample} complete reads={scanned} variants={len(variants)} nonvariants={len(nonvariants)}', flush=True)
    finally:
        read_handle.close()
        excluded_handle.close()
        display_handle.close()

    write_table(output / 'event_regions.tsv', ['EventID','Sample','Chrom','Start','End'], region_rows)
    summaries = []
    for event in consolidated:
        event_id = event['_id']
        if event_id not in metadata:
            continue
        event, chrom, target, ref, alt = metadata[event_id]
        count = counts[event_id]
        status, explanation, callable_count, fraction = evidence_status(count['alt'], count['ref'], count['excluded'])
        summaries.append({
            'EventID':event_id,'EvidenceClasses':event['EvidenceClasses'],'SourceEvents':event['SourceEvents'],
            'Sources':event['Sources'],'Sample':event['Sample'],'Gene':event.get('Gene',''),
            'Consequence':event.get('Consequence',''),'Impact':event.get('Impact',''),
            'PredictedConsequence':event.get('Consequence',''),'PredictedImpact':event.get('Impact',''),
            'Transcript':event.get('Transcript',''),'ProteinChange':event.get('ProteinChange',''),
            'Chrom':chrom,'Position':target,'REF':ref,'ALT':alt,'ReadValidationStatus':status,
            'ValidationExplanation':explanation,'CountUnit':'primary_alignments','UniqueAlignments':count['unique'],
            'CallableAlignments':callable_count,'ExactAltReads':count['alt'],'CleanReferenceReads':count['ref'],
            'ExcludedReads':count['excluded'],'AltFractionAmongClean':f'{fraction:.6f}' if fraction is not None else 'NA',
            'CallableFractionAmongExamined':f'{callable_count/count["unique"]:.6f}' if count['unique'] else 'NA',
        })

    write_table(output / 'findings_manifest.tsv', SUMMARY_FIELDS, summaries)
    write_table(output / 'bam_manifest.tsv', ['Sample','Category','BAM','Index','UniqueAlignments'], bam_manifest)
    (output / 'support_labels.bed').write_text('\n'.join(bed) + '\n')

    priority_genes = {value.strip() for value in args.priority_genes.split(',') if value.strip()}
    priority_impacts = {value.strip() for value in args.priority_impacts.split(',') if value.strip()}
    priority_consequences = {value.strip() for value in args.priority_consequences.split(',') if value.strip()}
    priority = [
        row for row in summaries
        if args.priority_mode == 'all'
        or (priority_genes and row['Gene'] in priority_genes)
        or (priority_impacts and row['PredictedImpact'] in priority_impacts)
        or (priority_consequences and set(re.split(r'[,;&]', row['PredictedConsequence'])) & priority_consequences)
    ]
    priority = sorted(priority, key=lambda row: row['EventID'])
    if args.priority_limit:
        priority = priority[:args.priority_limit]
    write_table(output / 'priority_findings.tsv', SUMMARY_FIELDS, priority)
    priority_ids = {row['EventID'] for row in priority}
    (output / 'priority_findings.bed').write_text(
        '\n'.join([bed[0]] + [line for line in bed[1:] if line.split('\t')[3].split()[0] in priority_ids]) + '\n'
    )

    visible = [output / row['BAM'] for row in bam_manifest if row['Category'] in {'event_display','exact_alt_display','reference_display'}]
    (output / 'review.igv.batch.txt').write_text('\n'.join(
        ['new', f'genome {Path(args.genome).resolve()}', f'load {(output / "support_labels.bed").resolve()}']
        + [f'load {path.resolve()}' for path in visible]
        + [f"goto {row['Chrom']}:{row['Position']}" for row in priority[:1000]]
        + ['expand','exit']
    ) + '\n')
    session = ET.Element('Session', genome=str(Path(args.genome).resolve()), version='8')
    resources = ET.SubElement(session, 'Resources')
    ET.SubElement(resources, 'Resource', path=str((output / 'support_labels.bed').resolve()))
    for path in visible:
        ET.SubElement(resources, 'Resource', path=str(path.resolve()))
    ET.ElementTree(session).write(output / 'igv.session.xml', encoding='utf-8', xml_declaration=True)

    lines = [
        f'Input selected rows: {len(selected)}',
        f'Consolidated findings: {len(consolidated)}',
        f'Merged duplicate rows: {len(selected)-len(consolidated)}',
        f'Generated findings: {len(summaries)}',
        f'Diagnostic read rows: {diagnostic_rows}',
        f'Excluded diagnostic rows: {excluded_rows}',
    ] + [
        f'{sample} scanned reads: {scanned}; variants: {variant_count}; nonvariant regions: {nonvariant_count}; exact ALT records: {exact}'
        for sample, scanned, variant_count, nonvariant_count, exact in stats
    ]
    (output / 'consolidation_summary.txt').write_text('\n'.join(lines) + '\n')
    (output / 'README.txt').write_text(
        'Event-exact alignment identities are stored in display_alignment_manifest.tsv.gz. The sample display BAMs are storage pools only. The explorer creates temporary event-specific BAMs before invoking igv-reports.\n'
    )
    print(f'Generated consolidated IGV review data for {len(summaries)} findings under {output}')


if __name__ == '__main__':
    main()
