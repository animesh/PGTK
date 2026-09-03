#!/usr/bin/env python3
import argparse
import csv
import gzip
import hashlib
from collections import defaultdict
from pathlib import Path

try:
    import pysam
except ImportError:
    pysam = None


def sequence_hash(sequence):
    return hashlib.sha256((sequence or '').encode('ascii', errors='ignore')).hexdigest()[:24]


def alignment_key(read, reference_name=None):
    contig = reference_name if reference_name is not None else read.reference_name
    values = (
        read.query_name or '', contig or '', int(read.reference_start), int(read.flag),
        read.cigarstring or '', sequence_hash(read.query_sequence),
    )
    return hashlib.sha256('\x1f'.join(map(str, values)).encode('utf-8')).hexdigest()


def load_event_keys(path, event_id, sample):
    keys = defaultdict(set)
    with gzip.open(path, 'rt', encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        required = {'EventID','Sample','Category','AlignmentKey'}
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f'Display manifest missing columns: {sorted(missing)}')
        for row in reader:
            if row['EventID'] == event_id and row['Sample'] == sample:
                keys[row['Category']].add(row['AlignmentKey'])
    return keys


def extract(source, keys, output):
    found = set()
    with pysam.AlignmentFile(source, 'rb') as incoming:
        with pysam.AlignmentFile(output, 'wb', template=incoming) as outgoing:
            for read in incoming.fetch(until_eof=True):
                key = alignment_key(read, incoming.get_reference_name(read.reference_id))
                if key in keys and key not in found:
                    outgoing.write(read)
                    found.add(key)
    missing = keys - found
    if missing:
        raise RuntimeError(f'{source}: {len(missing)} event-specific alignment identities were not found')
    pysam.index(str(output))


def merge_bams(paths, output):
    paths = [Path(path) for path in paths]
    with pysam.AlignmentFile(paths[0], 'rb') as first:
        temporary = output.with_suffix('.unsorted.bam')
        with pysam.AlignmentFile(temporary, 'wb', template=first) as outgoing:
            seen = set()
            for path in paths:
                with pysam.AlignmentFile(path, 'rb') as incoming:
                    for read in incoming.fetch(until_eof=True):
                        key = alignment_key(read, incoming.get_reference_name(read.reference_id))
                        if key not in seen:
                            outgoing.write(read)
                            seen.add(key)
    pysam.sort('-o', str(output), str(temporary))
    temporary.unlink()
    pysam.index(str(output))


def main():
    parser = argparse.ArgumentParser(description='Build event-exact temporary BAM tracks')
    parser.add_argument('--event-id', required=True)
    parser.add_argument('--sample', required=True)
    parser.add_argument('--display-manifest', required=True)
    parser.add_argument('--bam', action='append', default=[], help='CATEGORY=PATH')
    parser.add_argument('--output-dir', required=True)
    args = parser.parse_args()

    if pysam is None:
        raise RuntimeError('pysam is required')
    sources = dict(item.split('=', 1) for item in args.bam)
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    keys = load_event_keys(args.display_manifest, args.event_id, args.sample)

    generated = {}
    for category in ('exact_alt_display','reference_display','event_display'):
        category_keys = keys.get(category, set())
        if not category_keys:
            continue
        source = sources.get(category)
        if not source:
            raise RuntimeError(f'Missing source BAM for category {category}')
        target = output / f'{category}.bam'
        extract(source, category_keys, target)
        generated[category] = target

    callable_tracks = [generated[c] for c in ('exact_alt_display','reference_display') if c in generated]
    if callable_tracks:
        merge_bams(callable_tracks, output / 'event_display.bam')
        generated['event_display'] = output / 'event_display.bam'
    elif 'event_display' not in generated:
        raise RuntimeError(f'No display alignments recorded for event {args.event_id}')

    ordered = [generated[c] for c in ('event_display','exact_alt_display','reference_display') if c in generated]
    with (output / 'tracks.tsv').open('w', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['Category','BAM','Index'])
        for path in ordered:
            writer.writerow([path.stem, str(path), str(path) + '.bai'])


if __name__ == '__main__':
    main()
