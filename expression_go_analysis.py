#!/usr/bin/env python3
import argparse
import csv
import math
import re
from collections import defaultdict
from pathlib import Path

ATTR = re.compile(r'(\S+) "([^"]*)";')
ORA_FIELDS = [
    'Analysis','Subject','Group','GO_ID','GO_Name','Namespace',
    'ForegroundGenesInTerm','ForegroundGenesTested','BackgroundGenesInTerm',
    'BackgroundGenesTested','OddsRatio','PValue','FDR','OverlapGenes'
]
RANK_FIELDS = [
    'Analysis','Subject','Group','GO_ID','GO_Name','Namespace','GenesInTerm',
    'GenesRanked','MeanScore','ZScore','PValue','FDR','LeadingGenes'
]
SUMMARY_FIELDS = [
    'Analysis','Subject','Sample','BaselineSample','Status','Message',
    'ForegroundGenes','BackgroundGenes','GOTermsTested','FDRThreshold','SignificantGOTerms',
    'RankMetric','Pseudocount','NonZeroScores','PositiveScores','NegativeScores','MinScore','MaxScore'
]


def read_tsv(path):
    with open(path, encoding='utf-8', errors='replace', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write_tsv(path, rows, fields):
    with open(path, 'w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def bh(rows, p='PValue', q='FDR'):
    if not rows:
        return
    order = sorted(range(len(rows)), key=lambda index: float(rows[index][p]))
    adjusted = [1.0] * len(rows)
    previous = 1.0
    total = len(rows)
    for rank_index in range(total - 1, -1, -1):
        index = order[rank_index]
        previous = min(previous, float(rows[index][p]) * total / (rank_index + 1))
        adjusted[index] = previous
    for row, value in zip(rows, adjusted):
        row[q] = value


def log_combination(n, k):
    if k < 0 or k > n:
        return float('-inf')
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeom_right_tail(observed, population, annotated, selected):
    lower = max(observed, 0, selected - (population - annotated))
    upper = min(selected, annotated)
    if lower > upper:
        return 0.0
    probability = math.exp(
        log_combination(annotated, lower)
        + log_combination(population - annotated, selected - lower)
        - log_combination(population, selected)
    )
    total = probability
    current = lower
    while current < upper:
        denominator = (current + 1) * (population - annotated - selected + current + 1)
        if denominator <= 0:
            break
        probability *= ((annotated - current) * (selected - current)) / denominator
        total += probability
        current += 1
    return min(1.0, max(0.0, total))


def fisher_right(a, b, c, d):
    population = a + b + c + d
    selected = a + b
    annotated = a + c
    if population <= 0 or selected < 0 or annotated < 0:
        return 1.0, 0.0
    p_value = hypergeom_right_tail(a, population, annotated, selected)
    odds = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
    return p_value, odds


def tied_ranks(values):
    order = sorted(range(len(values)), key=lambda index: (values[index], index))
    ranks = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        average_rank = (start + 1 + end) / 2.0
        for position in range(start, end):
            ranks[order[position]] = average_rank
        start = end
    return ranks


def load_go(path, namespaces):
    terms = defaultdict(set)
    names = {}
    spaces = {}
    allowed = set(namespaces.split(',')) if namespaces not in {'all', '*', ''} else None
    with open(path, encoding='utf-8', newline='') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            if allowed and row['Namespace'] not in allowed:
                continue
            go_id = row['GO_ID']
            terms[go_id].add(row['Gene'])
            names[go_id] = row['GO_Name']
            spaces[go_id] = row['Namespace']
    return terms, names, spaces


def parse_gtf(path, feature, id_attribute, symbol_attribute, biotypes):
    genes = {}
    exon_intervals = defaultdict(lambda: defaultdict(list))
    allowed = set(biotypes.split(',')) if biotypes not in {'all', '*', ''} else None
    with open(path, encoding='utf-8', errors='replace') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9:
                continue
            attrs = dict(ATTR.findall(fields[8]))
            gene_id = attrs.get(id_attribute, '')
            if not gene_id:
                continue
            biotype = attrs.get('gene_biotype') or attrs.get('gene_type') or attrs.get('transcript_biotype', '')
            if allowed and biotype not in allowed:
                continue
            if gene_id not in genes:
                genes[gene_id] = {
                    'Gene_ID': gene_id,
                    'Gene': attrs.get(symbol_attribute, gene_id),
                    'Chromosome': fields[0],
                    'Start': int(fields[3]),
                    'End': int(fields[4]),
                    'Strand': fields[6],
                    'Biotype': biotype,
                }
            else:
                genes[gene_id]['Start'] = min(genes[gene_id]['Start'], int(fields[3]))
                genes[gene_id]['End'] = max(genes[gene_id]['End'], int(fields[4]))
            if fields[2] == feature:
                exon_intervals[gene_id][fields[0]].append((int(fields[3]), int(fields[4])))
    for gene_id, record in genes.items():
        merged = []
        for values in exon_intervals.get(gene_id, {}).values():
            for start, end in sorted(values):
                if merged and start <= merged[-1][1] + 1:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], end))
                else:
                    merged.append((start, end))
        record['Length'] = sum(end - start + 1 for start, end in merged) or record['End'] - record['Start'] + 1
    return genes


def parse_featurecounts(path):
    sample = Path(path).name.split('.gene_counts.tsv')[0]
    counts = {}
    with open(path, encoding='utf-8', errors='replace') as handle:
        reader = csv.reader((line for line in handle if not line.startswith('#')), delimiter='\t')
        header = next(reader)
        gene_index = header.index('Geneid')
        count_index = len(header) - 1
        for row in reader:
            if row:
                counts[row[gene_index]] = float(row[count_index])
    return sample, counts


def cmd_merge(args):
    genes = parse_gtf(args.gtf, args.feature_type, args.id_attribute, args.symbol_attribute, args.biotypes)
    sample_counts = {}
    for path in args.counts:
        sample, counts = parse_featurecounts(path)
        if sample in sample_counts:
            raise SystemExit(f'duplicate gene count sample: {sample}')
        sample_counts[sample] = counts
    samples = sorted(sample_counts)
    totals = {sample: sum(sample_counts[sample].values()) for sample in samples}
    rows = []
    for gene_id in sorted(genes):
        gene = genes[gene_id]
        row = dict(gene)
        for sample in samples:
            count = sample_counts[sample].get(gene_id, 0.0)
            length_kb = gene['Length'] / 1000.0
            row[f'{sample}_raw_count'] = int(count) if count.is_integer() else count
            row[f'{sample}_CPM'] = count * 1e6 / totals[sample] if totals[sample] else 0.0
            row[f'{sample}_RPK'] = count / length_kb if length_kb else 0.0
        rows.append(row)
    for sample in samples:
        total_rpk = sum(float(row[f'{sample}_RPK']) for row in rows)
        for row in rows:
            row[f'{sample}_TPM'] = float(row[f'{sample}_RPK']) * 1e6 / total_rpk if total_rpk else 0.0
    fields = ['Gene_ID','Gene','Chromosome','Start','End','Strand','Biotype','Length']
    fields += [f'{sample}_{metric}' for sample in samples for metric in ('raw_count','CPM','TPM')]
    cleaned = [{key: value for key, value in row.items() if not key.endswith('_RPK')} for row in rows]
    write_tsv(args.output_prefix + '.gene_expression.tsv', cleaned, fields)
    summary = [
        {
            'Sample': sample,
            'AssignedReads': totals[sample],
            'GenesWithReads': sum(sample_counts[sample].get(gene_id, 0) > 0 for gene_id in genes),
        }
        for sample in samples
    ]
    write_tsv(args.output_prefix + '.summary.tsv', summary, ['Sample','AssignedReads','GenesWithReads'])


def matrix_columns(matrix_path):
    with open(matrix_path, encoding='utf-8', errors='replace', newline='') as handle:
        header = next(csv.reader(handle, delimiter='\t'), [])
    if not header:
        raise SystemExit(f'expression matrix is empty: {matrix_path}')
    duplicates = sorted({column for column in header if header.count(column) > 1})
    if duplicates:
        raise SystemExit(f'expression matrix has duplicate columns: {",".join(duplicates)}')
    return set(header)


def require_matrix_columns(matrix_path, columns):
    available = matrix_columns(matrix_path)
    missing = sorted(set(columns) - available)
    if missing:
        raise SystemExit(f'expression matrix lacks required columns: {",".join(missing)}')


def matrix_context(matrix_path, sample_names, background):
    matrix = read_tsv(matrix_path)
    genome = {row['Gene'] for row in matrix if row['Gene']}
    measurable = {
        row['Gene'] for row in matrix
        if row['Gene'] and any(float(row.get(f'{sample}_raw_count', 0)) > 0 for sample in sample_names)
    }
    universe = genome if background == 'genome' else measurable
    return matrix, universe


def ora(label, subject, group, foreground, universe, terms, names, spaces, min_size, max_size):
    universe = set(universe)
    foreground = set(foreground) & universe
    foreground_size = len(foreground)
    universe_size = len(universe)
    rows = []
    for go_id, genes in terms.items():
        tested = genes & universe
        term_size = len(tested)
        if not min_size <= term_size <= max_size:
            continue
        overlap = foreground & tested
        a = len(overlap)
        b = foreground_size - a
        c = term_size - a
        d = universe_size - a - b - c
        p_value, odds = fisher_right(a, b, c, d)
        rows.append({
            'Analysis': label, 'Subject': subject, 'Group': group,
            'GO_ID': go_id, 'GO_Name': names[go_id], 'Namespace': spaces[go_id],
            'ForegroundGenesInTerm': a, 'ForegroundGenesTested': foreground_size,
            'BackgroundGenesInTerm': term_size, 'BackgroundGenesTested': universe_size,
            'OddsRatio': odds, 'PValue': p_value, 'OverlapGenes': ';'.join(sorted(overlap)),
        })
    bh(rows)
    rows.sort(key=lambda row: (row['FDR'], row['PValue'], row['GO_ID']))
    return rows


def ranked(label, subject, group, scores, universe, terms, names, spaces, min_size, max_size):
    eligible = sorted((gene, float(value)) for gene, value in scores.items() if gene in universe)
    genes = [gene for gene, value in eligible]
    values = [value for gene, value in eligible]
    n = len(genes)
    gene_set = set(genes)
    ranks_array = tied_ranks(values)
    ranks = dict(zip(genes, ranks_array))
    rows = []
    for go_id, term_genes in terms.items():
        inside = term_genes & gene_set
        m = len(inside)
        if not min_size <= m <= max_size or m == n:
            continue
        rank_sum = sum(ranks[gene] for gene in inside)
        u_value = rank_sum - m * (m + 1) / 2
        mean = m * (n - m) / 2
        standard_deviation = math.sqrt(m * (n - m) * (n + 1) / 12)
        z_score = (u_value - mean) / standard_deviation if standard_deviation else 0.0
        p_value = min(1.0, math.erfc(abs(z_score) / math.sqrt(2.0)))
        mean_score = sum(scores[gene] for gene in inside) / m
        rows.append({
            'Analysis': label, 'Subject': subject, 'Group': group,
            'GO_ID': go_id, 'GO_Name': names[go_id], 'Namespace': spaces[go_id],
            'GenesInTerm': m, 'GenesRanked': n, 'MeanScore': mean_score,
            'ZScore': z_score, 'PValue': p_value,
            'LeadingGenes': ';'.join(sorted(inside, key=lambda gene: abs(scores[gene]), reverse=True)[:50]),
        })
    bh(rows)
    rows.sort(key=lambda row: (row['FDR'], row['PValue'], row['GO_ID']))
    return rows


def load_samples(path):
    with open(path, encoding='utf-8', newline='') as handle:
        rows = list(csv.DictReader(handle))
    return [
        {
            'sample': row['sample'].strip(),
            'subject': (row.get('TK') or row['sample']).strip(),
            'group': (row.get('Group') or row['sample']).strip(),
            'baseline': (row.get('baseline') or 'false').strip().lower(),
        }
        for row in rows
    ]


def cmd_sample_ora(args):
    sample_names = [item.strip() for item in args.all_samples.split(',') if item.strip()]
    matrix, universe = matrix_context(args.matrix, sample_names, args.background)
    terms, names, spaces = load_go(args.go_mapping, args.namespaces)
    expressed = {
        row['Gene'] for row in matrix
        if row['Gene']
        and float(row.get(f'{args.sample}_CPM', 0)) >= args.cpm_threshold
        and float(row.get(f'{args.sample}_TPM', 0)) >= args.tpm_threshold
    }
    label = f'{args.sample}:expressed'
    rows = ora(label, args.subject, args.group, expressed, universe, terms, names, spaces, args.go_min_size, args.go_max_size)
    summary = [{
        'Analysis': label, 'Subject': args.subject, 'Sample': args.sample, 'BaselineSample': '',
        'Status': 'OK', 'Message': '', 'ForegroundGenes': len(expressed & universe),
        'BackgroundGenes': len(universe), 'GOTermsTested': len(rows),
        'FDRThreshold': args.fdr_threshold, 'SignificantGOTerms': sum(float(row['FDR']) <= args.fdr_threshold for row in rows),
    }]
    write_tsv(args.output_prefix + '.expression_ora.tsv', rows, ORA_FIELDS)
    write_tsv(args.output_prefix + '.summary.tsv', summary, SUMMARY_FIELDS)


def cmd_ranked_go(args):
    if args.sample == args.baseline_sample:
        raise SystemExit('--sample and --baseline-sample must be different')
    if args.pseudocount <= 0:
        raise SystemExit('--pseudocount must be greater than 0')
    if args.min_nonzero_scores < 1:
        raise SystemExit('--min-nonzero-scores must be at least 1')
    sample_names = [item.strip() for item in args.all_samples.split(',') if item.strip()]
    sample_column = f'{args.sample}_TPM'
    baseline_column = f'{args.baseline_sample}_TPM'
    require_matrix_columns(args.matrix, ['Gene', sample_column, baseline_column])
    matrix, universe = matrix_context(args.matrix, sample_names, args.background)
    terms, names, spaces = load_go(args.go_mapping, args.namespaces)
    scores = {
        row['Gene']: math.log2(
            (float(row[sample_column]) + args.pseudocount)
            / (float(row[baseline_column]) + args.pseudocount)
        )
        for row in matrix if row['Gene'] in universe
    }
    nonzero_scores = sum(value != 0.0 for value in scores.values())
    if nonzero_scores < args.min_nonzero_scores:
        raise SystemExit(
            f'ranked comparison {args.sample} versus {args.baseline_sample} has '
            f'{nonzero_scores} non-zero scores; required at least {args.min_nonzero_scores}'
        )
    positive_scores = sum(value > 0.0 for value in scores.values())
    negative_scores = sum(value < 0.0 for value in scores.values())
    label = f'{args.sample}_vs_{args.baseline_sample}:{args.rank_metric}'
    rows = ranked(label, args.subject, args.group, scores, universe, terms, names, spaces, args.go_min_size, args.go_max_size)
    summary = [{
        'Analysis': label, 'Subject': args.subject, 'Sample': args.sample,
        'BaselineSample': args.baseline_sample, 'Status': 'OK', 'Message': '',
        'ForegroundGenes': len(scores), 'BackgroundGenes': len(universe),
        'GOTermsTested': len(rows), 'FDRThreshold': args.fdr_threshold,
        'SignificantGOTerms': sum(float(row['FDR']) <= args.fdr_threshold for row in rows),
        'RankMetric': args.rank_metric, 'Pseudocount': args.pseudocount,
        'NonZeroScores': nonzero_scores, 'PositiveScores': positive_scores,
        'NegativeScores': negative_scores, 'MinScore': min(scores.values()),
        'MaxScore': max(scores.values()),
    }]
    write_tsv(args.output_prefix + '.ranked_go.tsv', rows, RANK_FIELDS)
    write_tsv(args.output_prefix + '.summary.tsv', summary, SUMMARY_FIELDS)


def cmd_merge_expression_go(args):
    metadata = load_samples(args.samples)
    ora_rows = [row for path in args.ora for row in read_tsv(path)]
    ranked_rows = [row for path in args.ranked for row in read_tsv(path)]
    summary_rows = [row for path in args.summary for row in read_tsv(path)]
    by_subject = defaultdict(list)
    for item in metadata:
        by_subject[item['subject']].append(item)
    for subject, members in sorted(by_subject.items()):
        baselines = [member for member in members if member['baseline'] == 'true']
        progressions = [member for member in members if member['baseline'] == 'false']
        if len(baselines) == 1:
            continue
        reason = 'no baseline sample' if not baselines else f'multiple baseline samples: {",".join(sorted(item["sample"] for item in baselines))}'
        for member in progressions:
            summary_rows.append({
                'Analysis': f'{member["sample"]}:ranked_GO_skipped', 'Subject': subject,
                'Sample': member['sample'], 'BaselineSample': '', 'Status': 'SKIPPED',
                'Message': reason, 'ForegroundGenes': 0, 'BackgroundGenes': 0,
                'GOTermsTested': 0, 'FDRThreshold': '', 'SignificantGOTerms': 0,
            })
    ora_rows.sort(key=lambda row: (row['Analysis'], float(row['FDR']), float(row['PValue']), row['GO_ID']))
    ranked_rows.sort(key=lambda row: (row['Analysis'], float(row['FDR']), float(row['PValue']), row['GO_ID']))
    summary_rows.sort(key=lambda row: (row['Subject'], row['Sample'], row['Analysis']))
    write_tsv(args.output_prefix + '.expression_ora.tsv', ora_rows, ORA_FIELDS)
    write_tsv(args.output_prefix + '.ranked_go.tsv', ranked_rows, RANK_FIELDS)
    write_tsv(args.output_prefix + '.summary.tsv', summary_rows, SUMMARY_FIELDS)


def cmd_expression(args):
    metadata = load_samples(args.samples)
    all_samples = ','.join(item['sample'] for item in metadata)
    ora_paths = []
    ranked_paths = []
    summary_paths = []
    prefix = Path(args.output_prefix)
    for item in metadata:
        task_prefix = str(prefix) + f'.{item["sample"]}.ora'
        task_values = dict(vars(args)); task_values.update(sample=item['sample'], subject=item['subject'], group=item['group'], all_samples=all_samples, output_prefix=task_prefix); task_args = argparse.Namespace(**task_values)
        cmd_sample_ora(task_args)
        ora_paths.append(task_prefix + '.expression_ora.tsv')
        summary_paths.append(task_prefix + '.summary.tsv')
    by_subject = defaultdict(list)
    for item in metadata:
        by_subject[item['subject']].append(item)
    for subject, members in by_subject.items():
        baselines = [item for item in members if item['baseline'] == 'true']
        if len(baselines) != 1:
            continue
        baseline = baselines[0]['sample']
        for item in members:
            if item['baseline'] == 'true':
                continue
            task_prefix = str(prefix) + f'.{item["sample"]}_vs_{baseline}.ranked'
            task_values = dict(vars(args)); task_values.update(sample=item['sample'], baseline_sample=baseline, subject=subject, group=item['group'], all_samples=all_samples, output_prefix=task_prefix); task_args = argparse.Namespace(**task_values)
            cmd_ranked_go(task_args)
            ranked_paths.append(task_prefix + '.ranked_go.tsv')
            summary_paths.append(task_prefix + '.summary.tsv')
    merge_args = argparse.Namespace(samples=args.samples, ora=ora_paths, ranked=ranked_paths, summary=summary_paths, output_prefix=args.output_prefix)
    cmd_merge_expression_go(merge_args)


def cmd_variant_sets(args):
    metadata = load_samples(args.samples)
    terms, names, spaces = load_go(args.go_mapping, args.namespaces)
    genome = parse_gtf(args.gtf, 'exon', args.id_attribute, args.symbol_attribute, args.biotypes)
    universe = {record['Gene'] for record in genome.values()} & set().union(*terms.values())
    sample_genes = defaultdict(set)
    for path in args.genes:
        for row in read_tsv(path):
            sample_genes[row['Sample']].add(row['Gene'])
    all_rows = []
    summaries = []
    by_subject = defaultdict(list)
    for item in metadata:
        if item['baseline'] == 'false' and item['sample'] in sample_genes:
            by_subject[item['subject']].append(item)
    for subject, members in by_subject.items():
        sets = [sample_genes[item['sample']] for item in members]
        common = set.intersection(*sets) if sets else set()
        member_names = ','.join(sorted(item['sample'] for item in members))
        analyses = [('common_all_progression', member_names, common)]
        for item in sorted(members, key=lambda value: value['sample']):
            complete = sample_genes[item['sample']]
            others = set().union(*(sample_genes[other['sample']] for other in members if other['sample'] != item['sample'])) if len(members) > 1 else set()
            analyses.append((f'{item["sample"]}_all', item['sample'], complete))
            analyses.append((f'{item["sample"]}_exclusive', item['sample'], complete - others))
        for label, group, foreground in analyses:
            rows = ora(label, subject, group, foreground, universe, terms, names, spaces, args.go_min_size, args.go_max_size)
            all_rows.extend(rows)
            summaries.append({
                'Subject': subject, 'Analysis': label, 'Members': group,
                'ForegroundGenes': len(foreground & universe), 'BackgroundGenes': len(universe),
                'GOTermsTested': len(rows), 'FDRThreshold': args.fdr_threshold,
                'SignificantGOTerms': sum(float(row['FDR']) <= args.fdr_threshold for row in rows),
            })
    write_tsv(args.output_prefix + '.variant_set_go.tsv', all_rows, ORA_FIELDS)
    write_tsv(args.output_prefix + '.summary.tsv', summaries, ['Subject','Analysis','Members','ForegroundGenes','BackgroundGenes','GOTermsTested','FDRThreshold','SignificantGOTerms'])


def add_common_go(parser):
    parser.add_argument('--matrix', required=True)
    parser.add_argument('--go-mapping', required=True)
    parser.add_argument('--background', choices=['genome','measurable'], default='genome')
    parser.add_argument('--all-samples', required=True)
    parser.add_argument('--go-min-size', type=int, default=10)
    parser.add_argument('--go-max-size', type=int, default=500)
    parser.add_argument('--fdr-threshold', type=float, default=0.1)
    parser.add_argument('--namespaces', default='all')
    parser.add_argument('--output-prefix', required=True)


def parser():
    root = argparse.ArgumentParser()
    sub = root.add_subparsers(dest='command', required=True)

    merge = sub.add_parser('merge-counts')
    merge.add_argument('--counts', nargs='+', required=True)
    merge.add_argument('--gtf', required=True)
    merge.add_argument('--feature-type', default='exon')
    merge.add_argument('--id-attribute', default='gene_id')
    merge.add_argument('--symbol-attribute', default='gene_name')
    merge.add_argument('--biotypes', default='all')
    merge.add_argument('--output-prefix', required=True)
    merge.set_defaults(func=cmd_merge)

    sample_ora = sub.add_parser('sample-ora')
    add_common_go(sample_ora)
    sample_ora.add_argument('--sample', required=True)
    sample_ora.add_argument('--subject', required=True)
    sample_ora.add_argument('--group', required=True)
    sample_ora.add_argument('--cpm-threshold', type=float, default=1.0)
    sample_ora.add_argument('--tpm-threshold', type=float, default=0.0)
    sample_ora.set_defaults(func=cmd_sample_ora)

    ranked_go = sub.add_parser('ranked-go')
    add_common_go(ranked_go)
    ranked_go.add_argument('--sample', required=True)
    ranked_go.add_argument('--baseline-sample', required=True)
    ranked_go.add_argument('--subject', required=True)
    ranked_go.add_argument('--group', required=True)
    ranked_go.add_argument('--pseudocount', type=float, default=0.5)
    ranked_go.add_argument('--rank-metric', choices=['log2_tpm_fold_change'], default='log2_tpm_fold_change')
    ranked_go.add_argument('--min-nonzero-scores', type=int, default=1)
    ranked_go.set_defaults(func=cmd_ranked_go)

    merge_go = sub.add_parser('merge-expression-go')
    merge_go.add_argument('--samples', required=True)
    merge_go.add_argument('--ora', nargs='*', default=[])
    merge_go.add_argument('--ranked', nargs='*', default=[])
    merge_go.add_argument('--summary', nargs='*', default=[])
    merge_go.add_argument('--output-prefix', required=True)
    merge_go.set_defaults(func=cmd_merge_expression_go)

    expression = sub.add_parser('expression-go')
    expression.add_argument('--matrix', required=True)
    expression.add_argument('--samples', required=True)
    expression.add_argument('--go-mapping', required=True)
    expression.add_argument('--background', choices=['genome','measurable'], default='genome')
    expression.add_argument('--cpm-threshold', type=float, default=1.0)
    expression.add_argument('--tpm-threshold', type=float, default=0.0)
    expression.add_argument('--pseudocount', type=float, default=0.5)
    expression.add_argument('--rank-metric', choices=['log2_tpm_fold_change'], default='log2_tpm_fold_change')
    expression.add_argument('--min-nonzero-scores', type=int, default=1)
    expression.add_argument('--go-min-size', type=int, default=10)
    expression.add_argument('--go-max-size', type=int, default=500)
    expression.add_argument('--fdr-threshold', type=float, default=0.1)
    expression.add_argument('--namespaces', default='all')
    expression.add_argument('--output-prefix', required=True)
    expression.set_defaults(func=cmd_expression)

    variant = sub.add_parser('variant-sets')
    variant.add_argument('--genes', nargs='+', required=True)
    variant.add_argument('--background', choices=['genome'], default='genome')
    variant.add_argument('--samples', required=True)
    variant.add_argument('--go-mapping', required=True)
    variant.add_argument('--gtf', required=True)
    variant.add_argument('--id-attribute', default='gene_id')
    variant.add_argument('--symbol-attribute', default='gene_name')
    variant.add_argument('--biotypes', default='protein_coding')
    variant.add_argument('--go-min-size', type=int, default=10)
    variant.add_argument('--go-max-size', type=int, default=500)
    variant.add_argument('--fdr-threshold', type=float, default=0.1)
    variant.add_argument('--namespaces', default='all')
    variant.add_argument('--output-prefix', required=True)
    variant.set_defaults(func=cmd_variant_sets)
    return root


if __name__ == '__main__':
    arguments = parser().parse_args()
    if hasattr(arguments, 'fdr_threshold') and not 0.0 <= arguments.fdr_threshold <= 1.0:
        raise SystemExit('--fdr-threshold must be between 0 and 1')
    arguments.func(arguments)
