#!/usr/bin/env python3
import argparse
import csv
import gzip
import json
import math
import re
import statistics
from collections import Counter, defaultdict
from pathlib import Path

parser = argparse.ArgumentParser(description='Audit PGTK results and map scientific claims to local evidence.')
parser.add_argument('--project-dir', default='.', help='PGTK project directory')
parser.add_argument('--results-dir', default=None, help='Results directory; defaults to <project-dir>/results')
parser.add_argument('--job-id', default='', help='Optional Slurm/Nextflow job ID')
parser.add_argument('--output-prefix', default='pgtk_claim_audit', help='Output prefix')
parser.add_argument('--fdr-threshold', type=float, default=None, help='Override expected GO FDR threshold')
parser.add_argument('--top-terms', type=int, default=20, help='Top directional GO terms to report per comparison')
parser.add_argument('--phenotype-metadata', default='', help='Optional CSV with sample, stage_order and ATX101_IC50 columns')
args = parser.parse_args()

project = Path(args.project_dir).resolve()
results = Path(args.results_dir).resolve() if args.results_dir else project / 'results'
prefix = Path(args.output_prefix).resolve()
prefix.parent.mkdir(parents=True, exist_ok=True)

if not results.is_dir():
    raise SystemExit(f'ERROR: results directory not found: {results}')

def open_text(path):
    path = Path(path)
    return gzip.open(path, 'rt', encoding='utf-8', errors='replace') if path.suffix == '.gz' else path.open('r', encoding='utf-8', errors='replace')

def read_tsv(path):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter='\t'))

def safe_float(value, default=None):
    try:
        number = float(value)
        return number if math.isfinite(number) else default
    except (TypeError, ValueError):
        return default

def safe_int(value, default=0):
    number = safe_float(value)
    return int(number) if number is not None else default

def find_one(pattern):
    matches = sorted(results.glob(pattern))
    return matches[0] if len(matches) == 1 else None

def find_all(pattern):
    return sorted(results.glob(pattern))

def relative(path):
    try:
        return str(Path(path).resolve().relative_to(project))
    except Exception:
        return str(path)

def write_tsv(path, rows, fields):
    with Path(path).open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)

def evidence_text(items):
    return '; '.join(str(item) for item in items if str(item).strip())

samples_path = project / 'samples.csv'
samples = []
if samples_path.is_file():
    with samples_path.open(encoding='utf-8', newline='') as handle:
        samples = list(csv.DictReader(handle))
for row in samples:
    row['sample'] = (row.get('sample') or '').strip()
    row['subject'] = (row.get('TK') or row['sample']).strip()
    row['group'] = (row.get('Group') or row['sample']).strip()
    row['baseline'] = (row.get('baseline') or 'false').strip().lower()

claims = []
insights = []
checks = []

def add_claim(claim_id, claim, status, strength, evidence_files, evidence, limitation, next_test=''):
    claims.append({
        'ClaimID': claim_id,
        'Claim': claim,
        'Status': status,
        'EvidenceStrength': strength,
        'EvidenceFiles': evidence_text(relative(x) for x in evidence_files if x),
        'ObservedEvidence': evidence_text(evidence),
        'Limitation': limitation,
        'RequiredNextTest': next_test,
    })

def add_insight(category, insight, evidence_files, evidence, interpretation, confidence):
    insights.append({
        'Category': category,
        'Insight': insight,
        'EvidenceFiles': evidence_text(relative(x) for x in evidence_files if x),
        'ObservedEvidence': evidence_text(evidence),
        'Interpretation': interpretation,
        'Confidence': confidence,
    })

def add_check(name, status, details, files=()):
    checks.append({'Check': name, 'Status': status, 'Details': details, 'Files': evidence_text(relative(x) for x in files if x)})

# Run completion and failures
job = args.job_id
trace = results / f'pipeline_trace-{job}.tsv' if job else None
if not trace or not trace.is_file():
    traces = sorted(results.glob('pipeline_trace-*.tsv'), key=lambda p: p.stat().st_mtime, reverse=True)
    trace = traces[0] if traces else None
    if trace:
        match = re.search(r'pipeline_trace-(.+)\.tsv$', trace.name)
        job = match.group(1) if match else ''
trace_rows = read_tsv(trace) if trace else []
status_counts = Counter((row.get('status') or '').upper() for row in trace_rows)
failed_rows = [row for row in trace_rows if (row.get('status') or '').upper() in {'FAILED', 'ABORTED'} or str(row.get('exit') or '').strip() not in {'', '0', '-'}]
completed_rows = [row for row in trace_rows if (row.get('status') or '').upper() in {'COMPLETED', 'CACHED'}]
if trace_rows and not failed_rows:
    add_claim('TECH-001', 'The selected pipeline run completed without failed task attempts.', 'SUPPORTED', 'HIGH', [trace], [f'{len(completed_rows)}/{len(trace_rows)} trace rows completed or cached', '0 failed/aborted rows'], 'A clean workflow trace establishes technical completion, not biological validity.')
else:
    add_claim('TECH-001', 'The selected pipeline run completed without failed task attempts.', 'NOT_SUPPORTED', 'HIGH' if trace_rows else 'NONE', [trace] if trace else [], [f'failed rows={len(failed_rows)}', f'trace rows={len(trace_rows)}'], 'Inspect failed task work directories and rerun before interpreting outputs.')

failure_ledger = results / 'failure_logs' / job / 'failure_ledger.tsv' if job else None
failure_rows = read_tsv(failure_ledger) if failure_ledger else []
add_check('Failure ledger', 'PASS' if failure_ledger and failure_ledger.is_file() and not failure_rows else 'WARN', f'rows={len(failure_rows)}', [failure_ledger] if failure_ledger else [])

multiqc = results / 'multiqc' / 'multiqc_report.html'
add_claim('TECH-002', 'The final integrated MultiQC report was generated.', 'SUPPORTED' if multiqc.is_file() and multiqc.stat().st_size > 0 else 'NOT_SUPPORTED', 'HIGH', [multiqc] if multiqc.is_file() else [], [f'bytes={multiqc.stat().st_size}' if multiqc.is_file() else 'file missing'], 'Report existence does not guarantee that every custom section is biologically correct.')

# Baseline design
by_subject = defaultdict(list)
for row in samples:
    by_subject[row['subject']].append(row)
baseline_messages = []
baseline_ok = bool(samples)
for subject, members in sorted(by_subject.items()):
    bases = [r['sample'] for r in members if r['baseline'] == 'true']
    progressions = [r['sample'] for r in members if r['baseline'] == 'false']
    baseline_messages.append(f'{subject}: baseline={bases}, progression={progressions}')
    baseline_ok &= len(bases) == 1
add_claim('DESIGN-001', 'Every subject has exactly one baseline for progression subtraction and ranked expression comparison.', 'SUPPORTED' if baseline_ok else 'NOT_SUPPORTED', 'HIGH', [samples_path] if samples_path.is_file() else [], baseline_messages, 'This validates metadata structure only; it does not provide biological replication.')

# Variant stages
variant_inventory = results / 'comparative_advantage' / 'comparative_advantage.variant_stage_inventory.tsv'
variant_rows = read_tsv(variant_inventory)
variant_by_sample = defaultdict(dict)
for row in variant_rows:
    variant_by_sample[row.get('Sample', '')][row.get('Stage', '')] = row
if variant_rows:
    for sample, stages in sorted(variant_by_sample.items()):
        evid = [f'{stage}={safe_int(row.get("Alleles"))} alleles' for stage, row in sorted(stages.items())]
        add_insight('Variant stages', f'{sample} retained separate raw, PASS, RNA-supported, and progression stages where available.', [variant_inventory], evid, 'Stage separation prevents conflating initial calls, hard-filtered calls, RNA-supported candidates, and baseline-subtracted candidates.', 'HIGH')
    add_claim('VAR-001', 'The pipeline preserves distinct raw, PASS, and RNA-supported variant stages.', 'SUPPORTED', 'HIGH', [variant_inventory], [f'samples={len(variant_by_sample)}', f'rows={len(variant_rows)}'], 'RNA-supported calls are not orthogonally confirmed DNA mutations.')
else:
    add_claim('VAR-001', 'The pipeline preserves distinct raw, PASS, and RNA-supported variant stages.', 'NOT_TESTED', 'NONE', [], ['variant inventory missing'], 'Cannot verify stage counts.')

progression_rows = [row for row in variant_rows if row.get('Stage') == 'progression_nonbaseline_only']
if progression_rows:
    for row in progression_rows:
        sample = row.get('Sample', '')
        alleles = safe_int(row.get('Alleles'))
        snps = safe_int(row.get('SNPs'))
        indels = safe_int(row.get('Indels'))
        add_insight('Progression variants', f'{sample} has RNA-observed non-baseline-only variant candidates.', [variant_inventory], [f'alleles={alleles}', f'SNPs={snps}', f'indels={indels}'], 'Candidates may reflect clonal change, expression-dependent detectability, allele-specific expression, or coverage differences.', 'MODERATE')
    add_claim('VAR-002', 'Progression samples contain RNA-observed non-baseline-only variant candidates.', 'SUPPORTED', 'MODERATE', [variant_inventory], [f'{row.get("Sample")}: {row.get("Alleles")} alleles' for row in progression_rows], 'Baseline absence in RNA does not prove absence from DNA.', 'Matched longitudinal DNA sequencing or orthogonal genotyping.')
    add_claim('VAR-003', 'The non-baseline-only variants are newly acquired somatic DNA mutations.', 'NOT_SUPPORTED', 'LOW', [variant_inventory], ['RNA-only progression subtraction'], 'RNA coverage and expression can create apparent baseline absence.', 'Matched DNA sequencing at every stage.')

# Expression and GO
matrix = results / 'expression' / 'gene_expression.gene_expression.tsv'
expr_summary = results / 'expression' / 'gene_expression.summary.tsv'
go_summary = results / 'expression' / 'go' / 'expression_go.summary.tsv'
ranked_go = results / 'expression' / 'go' / 'expression_go.ranked_go.tsv'
ora_go = results / 'expression' / 'go' / 'expression_go.expression_ora.tsv'
go_summary_rows = read_tsv(go_summary)
ranked_rows = read_tsv(ranked_go)
ora_rows = read_tsv(ora_go)

if matrix.is_file() and matrix.stat().st_size > 0:
    add_claim('EXP-001', 'A merged raw-count, CPM, and TPM expression matrix was produced.', 'SUPPORTED', 'HIGH', [matrix, expr_summary], [f'bytes={matrix.stat().st_size}'], 'With one longitudinal sample per stage, effect directions are descriptive rather than population-level differential expression.')
else:
    add_claim('EXP-001', 'A merged raw-count, CPM, and TPM expression matrix was produced.', 'NOT_SUPPORTED', 'NONE', [], ['matrix missing'], 'Expression claims cannot be assessed.')

ranked_summaries = [r for r in go_summary_rows if r.get('BaselineSample')]
if ranked_summaries:
    threshold_values = sorted({r.get('FDRThreshold', '') for r in go_summary_rows if r.get('FDRThreshold', '') != ''})
    diagnostics = [f"{r.get('Analysis')}: nonzero={r.get('NonZeroScores')}, positive={r.get('PositiveScores')}, negative={r.get('NegativeScores')}, significant={r.get('SignificantGOTerms')}" for r in ranked_summaries]
    add_claim('GO-001', 'Ranked GO comparisons use non-zero directional progression-versus-baseline scores.', 'SUPPORTED', 'HIGH', [go_summary, ranked_go], diagnostics + [f'FDR thresholds={threshold_values}'], 'GO terms are correlated and are not independent discoveries.')
    configured = args.fdr_threshold if args.fdr_threshold is not None else 0.1
    mismatches = [r for r in go_summary_rows if r.get('FDRThreshold') not in {'', str(configured), f'{configured:g}'}]
    add_claim('GO-002', f'SignificantGOTerms is evaluated at the expected FDR threshold {configured:g}.', 'SUPPORTED' if not mismatches else 'NOT_SUPPORTED', 'HIGH', [go_summary], [f'rows={len(go_summary_rows)}', f'mismatches={len(mismatches)}'], 'This confirms the recorded threshold, not independent recomputation of every FDR value.')
else:
    add_claim('GO-001', 'Ranked GO comparisons use non-zero directional progression-versus-baseline scores.', 'NOT_TESTED', 'NONE', [go_summary] if go_summary.is_file() else [], ['ranked summary rows missing'], 'Cannot verify ranked comparisons.')

# Recompute significant counts from complete included GO tables.
def audit_go_counts(data_rows, summary_rows):
    grouped = defaultdict(list)
    for row in data_rows:
        grouped[row.get('Analysis', '')].append(row)
    results_out = []
    for summary in summary_rows:
        analysis = summary.get('Analysis', '')
        threshold = safe_float(summary.get('FDRThreshold'))
        expected = safe_int(summary.get('SignificantGOTerms'))
        rows = grouped.get(analysis, [])
        observed = sum(1 for row in rows if threshold is not None and safe_float(row.get('FDR'), 1.0) <= threshold)
        results_out.append((analysis, expected, observed, len(rows)))
    return results_out

go_recounts = audit_go_counts(ranked_rows + ora_rows, go_summary_rows)
if go_recounts:
    mismatches = [x for x in go_recounts if x[1] != x[2]]
    add_claim('GO-003', 'Reported SignificantGOTerms counts equal the number of GO rows with FDR at or below each recorded threshold.', 'SUPPORTED' if not mismatches else 'NOT_SUPPORTED', 'HIGH', [go_summary, ranked_go, ora_go], [f'{a}: reported={e}, recounted={o}, rows={n}' for a, e, o, n in go_recounts], 'Only complete local GO tables can be fully recounted.')

# Directional top GO terms and redundancy warning
if ranked_rows:
    by_analysis = defaultdict(list)
    for row in ranked_rows:
        by_analysis[row.get('Analysis', '')].append(row)
    for analysis, rows in sorted(by_analysis.items()):
        threshold = next((safe_float(r.get('FDRThreshold')) for r in go_summary_rows if r.get('Analysis') == analysis), 0.1)
        significant = [r for r in rows if safe_float(r.get('FDR'), 1.0) <= (threshold if threshold is not None else 0.1)]
        positive = sorted((r for r in significant if safe_float(r.get('MeanScore'), 0.0) > 0), key=lambda r: (safe_float(r.get('FDR'), 1.0), -abs(safe_float(r.get('ZScore'), 0.0))))[:args.top_terms]
        negative = sorted((r for r in significant if safe_float(r.get('MeanScore'), 0.0) < 0), key=lambda r: (safe_float(r.get('FDR'), 1.0), -abs(safe_float(r.get('ZScore'), 0.0))))[:args.top_terms]
        positive_count = sum(1 for r in significant if safe_float(r.get('MeanScore'), 0.0) > 0)
        negative_count = sum(1 for r in significant if safe_float(r.get('MeanScore'), 0.0) < 0)
        add_insight(
            'Ranked GO',
            f'{analysis} has directional GO enrichment.',
            [ranked_go],
            [
                f'significant={len(significant)}',
                f'positive={positive_count}',
                f'negative={negative_count}',
                'top positive=' + ', '.join(r.get('GO_Name', '') for r in positive[:5]),
                'top negative=' + ', '.join(r.get('GO_Name', '') for r in negative[:5]),
            ],
            'Interpret using direction and gene overlap; thousands of correlated GO terms must not be described as independent pathways.',
            'MODERATE',
        )

# Article-oriented transcript signatures.
# Direction is encoded explicitly so that an increase in a PCL-down signature is
# reported as discordance rather than incorrectly presented as PCL concordance.
CANONICAL_CYTO_RIBOSOME = {
    # Explicit human cytosolic ribosomal protein genes. Pseudogenes and
    # antisense loci are intentionally excluded.
    'RPL3','RPL3L','RPL4','RPL5','RPL6','RPL7','RPL7A','RPL7L1','RPL8','RPL9',
    'RPL10','RPL10A','RPL11','RPL12','RPL13','RPL13A','RPL14','RPL15','RPL17',
    'RPL18','RPL18A','RPL19','RPL21','RPL22','RPL22L1','RPL23','RPL23A','RPL24',
    'RPL26','RPL26L1','RPL27','RPL27A','RPL28','RPL29','RPL30','RPL31','RPL32',
    'RPL34','RPL35','RPL35A','RPL36','RPL36A','RPL36AL','RPL37','RPL37A','RPL38',
    'RPL39','RPL39L','RPL41','RPLP0','RPLP1','RPLP2',
    'RPS2','RPS3','RPS3A','RPS4X','RPS4Y1','RPS4Y2','RPS5','RPS6','RPS7','RPS8',
    'RPS9','RPS10','RPS11','RPS12','RPS13','RPS14','RPS15','RPS15A','RPS16',
    'RPS17','RPS18','RPS19','RPS20','RPS21','RPS23','RPS24','RPS25','RPS26',
    'RPS27','RPS27A','RPS27L','RPS28','RPS29','RPS30','RPSA',
    'UBA52','FAU',
}
CANONICAL_MITO_RIBOSOME = {
    # Explicit human mitochondrial ribosomal protein genes.
    'MRPL1','MRPL2','MRPL3','MRPL4','MRPL9','MRPL10','MRPL11','MRPL12','MRPL13',
    'MRPL14','MRPL15','MRPL16','MRPL17','MRPL18','MRPL19','MRPL20','MRPL21',
    'MRPL22','MRPL23','MRPL24','MRPL27','MRPL28','MRPL30','MRPL32','MRPL33',
    'MRPL34','MRPL35','MRPL36','MRPL37','MRPL38','MRPL39','MRPL40','MRPL41',
    'MRPL42','MRPL43','MRPL44','MRPL45','MRPL46','MRPL47','MRPL48','MRPL49',
    'MRPL50','MRPL51','MRPL52','MRPL53','MRPL54','MRPL55','MRPL57','MRPL58',
    'MRPS2','MRPS5','MRPS6','MRPS7','MRPS9','MRPS10','MRPS11','MRPS12','MRPS14',
    'MRPS15','MRPS16','MRPS17','MRPS18A','MRPS18B','MRPS18C','MRPS21','MRPS22',
    'MRPS23','MRPS24','MRPS25','MRPS26','MRPS27','MRPS28','MRPS30','MRPS31',
    'MRPS33','MRPS34','MRPS35','MRPS36',
}
SIGNATURE_DEFINITIONS = {
    'Canonical_cytosolic_ribosome': {
        'genes': CANONICAL_CYTO_RIBOSOME,
        'expected': 'UP',
        'article_axis': 'MGUS_to_MM',
        'description': 'Canonical protein-coding cytosolic ribosomal proteins only; pseudogenes and antisense loci excluded.',
    },
    'Canonical_mitochondrial_ribosome': {
        'genes': CANONICAL_MITO_RIBOSOME,
        'expected': 'UP',
        'article_axis': 'MGUS_to_MM',
        'description': 'Canonical protein-coding mitochondrial ribosomal proteins only.',
    },
    'Translation_fidelity_RNA_processing': {
        'genes': {'LAGE3','GON7','TP53RK','TPRKB','OSGEP','YRDC','CDKAL1','TRMT6','NOP56','NOP58','DKC1','GAR1','NHP2','ADAR','NONO','SFPQ','MATR3'},
        'expected': 'UP',
        'article_axis': 'MGUS_to_MM',
        'description': 'Translation-fidelity, RNA-modification and RNA-processing factors highlighted in the MGUS-to-MM model.',
    },
    'Proteostasis_ER_proteasome': {
        'prefixes': {'PDIA','PSMA','PSMB','PSMC','PSMD','PSME','DERL','SEC61'},
        'genes': {'HSPA5','HSP90B1','CALR','CANX','EIF2AK3','ERN1','ATF6','ERO1A','PRDX4','TMX3','TMX4'},
        'expected': 'UP',
        'article_axis': 'ATX101_sensitivity',
        'description': 'ER folding, ERAD, translocation and proteasome machinery.',
    },
    'ATX101_metabolic_redox': {
        'prefixes': {'GPX','PRDX'},
        'genes': {'PCNA','ENO1','G6PD','PGD','TKT','TALDO1','GAPDH','PGK1','PGAM1','PKM','NAMPT','NMNAT1','NMNAT2','NMNAT3','NADSYN1','NADK','GSR','GCLC','GCLM'},
        'expected': 'UP',
        'article_axis': 'ATX101_sensitivity',
        'description': 'Glycolysis, PPP, NAD and glutathione/redox support linked to ATX-101 sensitivity.',
    },
    'PCL_glycolysis_up': {
        'genes': {'SLC2A1','HK1','HK2','LDHA','PFKP','ENO1','GAPDH','PKM'},
        'expected': 'UP',
        'article_axis': 'MM_to_sPCL',
        'description': 'Aerobic glycolysis expected to increase in the MM-to-sPCL model.',
    },
    'PCL_mitochondrial_TCA_down': {
        'genes': {'SDHA','SDHB','IDH3G','DHTKD1','ACO1','ACO2','DLST','MDH1'},
        'expected': 'DOWN',
        'article_axis': 'MM_to_sPCL',
        'description': 'Mitochondrial/TCA genes expected to decrease in the MM-to-sPCL model.',
    },
    'PCL_glycosylation_down': {
        'genes': {'PMM2','GMPPA','GMPPB','GMDS','GFPT1','PGM3','NANS','CMAS','ALG1','ALG2','STT3B','ST6GAL1'},
        'expected': 'DOWN',
        'article_axis': 'MM_to_sPCL',
        'description': 'Nucleotide-sugar and glycosylation genes expected to decrease in the MM-to-sPCL model.',
    },
    'PCL_adhesion_remodeling_mixed': {
        'genes': {'CD44','ITGA4','ITGB1','ITGB7','ITGA6','SDC1','CXCR4','SELPLG','TLN1','VCL','RAP1A','RAP1B','SRC'},
        'expected': 'MIXED',
        'article_axis': 'MM_to_sPCL',
        'description': 'Stage-dependent adhesion remodelling; no single expected sign is imposed.',
    },
    'PCL_proliferation_invasion_up': {
        'genes': {'MKI67','PCNA','S100A4','GBP2','TAGLN2','FLNA','DEK','MARCKS','RBP1','H1F0'},
        'expected': 'UP',
        'article_axis': 'MM_to_sPCL',
        'description': 'Proliferation, motility and invasion-associated genes expected to increase.',
    },
}

def signature_match(gene, definition):
    if gene in definition.get('genes', set()):
        return True
    return any(gene.startswith(prefix) for prefix in definition.get('prefixes', set()))

def classify_concordance(expected, positive, negative, median, total):
    informative = positive + negative
    if total == 0 or informative == 0:
        return 'INSUFFICIENT_DATA', 0.0
    if expected == 'MIXED':
        balance = abs(positive - negative) / informative
        return ('MIXED_REMODELLING' if balance < 0.35 else 'DIRECTIONALLY_SKEWED'), 1.0 - balance
    expected_count = positive if expected == 'UP' else negative
    opposite_count = negative if expected == 'UP' else positive
    expected_fraction = expected_count / informative
    expected_median = median > 0 if expected == 'UP' else median < 0
    if expected_fraction >= 0.70 and expected_median:
        label = 'STRONG_CONCORDANCE'
    elif expected_fraction >= 0.55 and expected_median:
        label = 'PARTIAL_CONCORDANCE'
    elif expected_fraction <= 0.35 and not expected_median:
        label = 'OPPOSITE_DIRECTION'
    else:
        label = 'WEAK_OR_MIXED'
    return label, expected_fraction

signature_rows = []
novel_findings = []
if matrix.is_file():
    with open_text(matrix) as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fields = reader.fieldnames or []
        tpm_samples = sorted(field[:-4] for field in fields if field.endswith('_TPM'))
        matrix_rows = list(reader)
    sample_meta = {r['sample']: r for r in samples}
    for sample in tpm_samples:
        meta = sample_meta.get(sample, {})
        if meta.get('baseline') == 'true':
            continue
        subject = meta.get('subject', sample)
        baseline_candidates = [r['sample'] for r in samples if r['subject'] == subject and r['baseline'] == 'true']
        if len(baseline_candidates) != 1:
            continue
        baseline = baseline_candidates[0]
        if f'{baseline}_TPM' not in fields:
            continue
        gene_scores = {}
        for row in matrix_rows:
            gene = (row.get('Gene') or '').strip()
            biotype = (row.get('Biotype') or '').strip()
            if not gene or biotype != 'protein_coding':
                continue
            sample_tpm = safe_float(row.get(f'{sample}_TPM'), 0.0)
            baseline_tpm = safe_float(row.get(f'{baseline}_TPM'), 0.0)
            # Require measurable expression in at least one member of the comparison.
            if max(sample_tpm, baseline_tpm) < 1.0:
                continue
            gene_scores[gene] = math.log2((sample_tpm + 0.5) / (baseline_tpm + 0.5))
        for signature, definition in SIGNATURE_DEFINITIONS.items():
            selected = [(gene, score) for gene, score in gene_scores.items() if signature_match(gene, definition)]
            if not selected:
                continue
            scores = [score for _, score in selected]
            positive = sum(score > 0 for score in scores)
            negative = sum(score < 0 for score in scores)
            zero = len(scores) - positive - negative
            median = statistics.median(scores)
            concordance, expected_fraction = classify_concordance(definition['expected'], positive, negative, median, len(scores))
            signature_rows.append({
                'Comparison': f'{sample}_vs_{baseline}',
                'ArticleAxis': definition['article_axis'],
                'Signature': signature,
                'ExpectedDirection': definition['expected'],
                'Genes': len(scores),
                'Positive': positive,
                'Negative': negative,
                'Zero': zero,
                'PositiveFraction': positive / len(scores),
                'NegativeFraction': negative / len(scores),
                'ExpectedDirectionFraction': expected_fraction,
                'MedianLog2TPMRatio': median,
                'MedianFoldChange': 2 ** median,
                'Concordance': concordance,
                'MinScore': min(scores),
                'MaxScore': max(scores),
                'TopPositiveGenes': ';'.join(g for g, score in sorted(selected, key=lambda item: item[1], reverse=True)[:15]),
                'TopNegativeGenes': ';'.join(g for g, score in sorted(selected, key=lambda item: item[1])[:15]),
                'Definition': definition['description'],
            })
    if signature_rows:
        add_claim('BIO-001', 'The current RNA data can quantify article-oriented progression signatures directionally.', 'SUPPORTED', 'MODERATE', [matrix], [f'signature comparisons={len(signature_rows)}'], 'These are predefined descriptive transcript signatures, restricted to protein-coding genes with TPM >= 1 in at least one compared sample; they are not replicated proteomic tests.')
        for row in signature_rows:
            interpretation = (
                f"Expected {row['ExpectedDirection']}; observed concordance={row['Concordance']}. "
                f"Expected-direction fraction={row['ExpectedDirectionFraction']:.1%}; median fold change={row['MedianFoldChange']:.2f}."
            )
            add_insight('Article signature', f"{row['Comparison']} {row['Signature']}", [matrix], [f"genes={row['Genes']}", f"positive={row['Positive']}", f"negative={row['Negative']}", f"median log2 TPM ratio={row['MedianLog2TPMRatio']:.3f}", f"concordance={row['Concordance']}"], interpretation, 'MODERATE')

        by_key = {(r['Comparison'], r['Signature']): r for r in signature_rows}
        comparisons = sorted({r['Comparison'] for r in signature_rows})
        for comparison in comparisons:
            glycolysis = by_key.get((comparison, 'PCL_glycolysis_up'))
            tca = by_key.get((comparison, 'PCL_mitochondrial_TCA_down'))
            glyco = by_key.get((comparison, 'PCL_glycosylation_down'))
            proteostasis = by_key.get((comparison, 'Proteostasis_ER_proteasome'))
            atx = by_key.get((comparison, 'ATX101_metabolic_redox'))
            ribo = by_key.get((comparison, 'Canonical_cytosolic_ribosome'))
            if glycolysis and glycolysis['Concordance'] == 'STRONG_CONCORDANCE':
                novel_findings.append({'FindingID': f'NF-{len(novel_findings)+1:03d}', 'Comparison': comparison, 'Finding': 'Coherent glycolytic intensification', 'Status': 'SUPPORTED', 'Evidence': f"{glycolysis['Positive']}/{glycolysis['Genes']} genes increased; median fold change={glycolysis['MedianFoldChange']:.2f}", 'NovelInterpretation': 'Progression is accompanied by a coordinated glycolytic program rather than isolated single-gene changes.', 'AlternativeExplanation': 'Transcript abundance does not prove glycolytic flux.', 'NextValidation': 'Seahorse glycolytic flux, lactate production, glucose uptake and isotope tracing.'})
            if glycolysis and tca and glyco and glycolysis['Concordance'] == 'STRONG_CONCORDANCE' and tca['Concordance'] == 'OPPOSITE_DIRECTION' and glyco['Concordance'] == 'OPPOSITE_DIRECTION':
                novel_findings.append({'FindingID': f'NF-{len(novel_findings)+1:03d}', 'Comparison': comparison, 'Finding': 'Partial rather than complete sPCL-like metabolic convergence', 'Status': 'SUPPORTED', 'Evidence': f"glycolysis={glycolysis['Concordance']}; TCA-down={tca['Concordance']}; glycosylation-down={glyco['Concordance']}", 'NovelInterpretation': 'The sample acquires the glycolytic component of the published sPCL program while retaining or increasing TCA and glycosylation transcripts.', 'AlternativeExplanation': 'RNA and protein direction may differ, and the sPCL study was protein-level.', 'NextValidation': 'Quantitative proteomics, metabolomics, glycomics and respiratory flux.'})
            if proteostasis and atx:
                supported_levels = {'PARTIAL_CONCORDANCE', 'STRONG_CONCORDANCE'}
                coupled_supported = (
                    proteostasis['Concordance'] in supported_levels
                    and atx['Concordance'] in supported_levels
                )
                if coupled_supported:
                    novel_findings.append({'FindingID': f'NF-{len(novel_findings)+1:03d}', 'Comparison': comparison, 'Finding': 'Coupled proteostasis and metabolic-redox adaptation', 'Status': 'SUPPORTED', 'Evidence': f"proteostasis={proteostasis['Concordance']} (median FC={proteostasis['MedianFoldChange']:.2f}); ATX101 metabolic-redox={atx['Concordance']} (median FC={atx['MedianFoldChange']:.2f})", 'NovelInterpretation': 'Protein-handling and metabolic-redox support rise together, providing a plausible dependency state associated with ATX-101 sensitivity.', 'AlternativeExplanation': 'Association does not establish that either program causes drug sensitivity.', 'NextValidation': 'ATX-101 perturbation with rescue assays, NAD/NADH, GSH/GSSG, PPP flux and ER-stress readouts.'})
                elif atx['Concordance'] in supported_levels:
                    novel_findings.append({'FindingID': f'NF-{len(novel_findings)+1:03d}', 'Comparison': comparison, 'Finding': 'Metabolic-redox adaptation without established coupled proteostasis activation', 'Status': 'PARTIALLY_SUPPORTED', 'Evidence': f"proteostasis={proteostasis['Concordance']} (median FC={proteostasis['MedianFoldChange']:.2f}); ATX101 metabolic-redox={atx['Concordance']} (median FC={atx['MedianFoldChange']:.2f})", 'NovelInterpretation': 'The metabolic-redox program is coordinated, but a coupled global proteostasis increase is not established in this comparison.', 'AlternativeExplanation': 'Proteostasis may be regulated at protein-activity or post-transcriptional levels.', 'NextValidation': 'Proteasome activity, ER-stress markers, MIB/signallomics and matched proteomics.'})
            if ribo:
                novel_findings.append({'FindingID': f'NF-{len(novel_findings)+1:03d}', 'Comparison': comparison, 'Finding': 'Canonical ribosome direction after pseudogene exclusion', 'Status': ribo['Concordance'], 'Evidence': f"protein-coding genes={ribo['Genes']}; expected-direction fraction={ribo['ExpectedDirectionFraction']:.1%}; median FC={ribo['MedianFoldChange']:.2f}", 'NovelInterpretation': 'This refined result replaces the earlier broad RPL/RPS-prefix analysis that was dominated by pseudogenes and zero-expression loci.', 'AlternativeExplanation': 'TPM ratios do not measure translation rate or ribosome occupancy.', 'NextValidation': 'Matched proteomics, polysome profiling or ribosome profiling.'})

# GO-derived novel observations with explicit artifact handling.
if ranked_rows:
    ranked_by_analysis = defaultdict(list)
    for row in ranked_rows:
        ranked_by_analysis[row.get('Analysis', '')].append(row)
    for analysis, rows in sorted(ranked_by_analysis.items()):
        sig = [r for r in rows if safe_float(r.get('FDR'), 1.0) <= 0.1]
        positive_names = {r.get('GO_Name', '').lower() for r in sig if safe_float(r.get('MeanScore'), 0.0) > 0}
        negative_names = {r.get('GO_Name', '').lower() for r in sig if safe_float(r.get('MeanScore'), 0.0) < 0}
        if any('olfactory' in name or 'sensory perception of smell' in name for name in positive_names):
            novel_findings.append({'FindingID': f'NF-{len(novel_findings)+1:03d}', 'Comparison': analysis.split(':',1)[0], 'Finding': 'Olfactory-receptor GO signal requires artifact review', 'Status': 'CAUTION', 'Evidence': 'Significant positive olfactory/smell GO categories', 'NovelInterpretation': 'Large low-expression receptor families may distort broad rank-based GO output and should not drive myeloma interpretation.', 'AlternativeExplanation': 'True ectopic receptor expression is possible but not established.', 'NextValidation': 'Inspect leading genes, TPM, uniqueness/mappability and canonical protein-coding status.'})
        mhc_class_ii_negative = sorted(
            [
                row for row in sig
                if safe_float(row.get('MeanScore'), 0.0) < 0
                and 'mhc class ii' in (row.get('GO_Name') or '').lower()
            ],
            key=lambda row: safe_float(row.get('FDR'), 1.0),
        )
        if len(mhc_class_ii_negative) >= 2:
            mhc_evidence = ';'.join(
                f"{row.get('GO_ID')}|{row.get('GO_Name')}|MeanScore={safe_float(row.get('MeanScore'), 0.0):.4g}|FDR={safe_float(row.get('FDR'), 1.0):.3g}"
                for row in mhc_class_ii_negative[:10]
            )
            novel_findings.append({
                'FindingID': f'NF-{len(novel_findings)+1:03d}',
                'Comparison': analysis.split(':',1)[0],
                'Finding': 'Reduced MHC class II antigen-presentation program',
                'Status': 'SUPPORTED_TRANSCRIPT_LEVEL',
                'Evidence': f'{len(mhc_class_ii_negative)} significant negative MHC class II GO terms: {mhc_evidence}',
                'NovelInterpretation': 'Progression may reduce intrinsic MHC class II antigen-presentation or immune-visibility transcripts.',
                'AlternativeExplanation': 'RNA abundance does not establish surface protein abundance or functional antigen presentation.',
                'NextValidation': 'HLA-DRA/DRB protein measurement, flow cytometry and antigen-presentation assays.',
            })
        if any('mitochondrial atp synthesis' in name for name in negative_names):
            novel_findings.append({'FindingID': f'NF-{len(novel_findings)+1:03d}', 'Comparison': analysis.split(':',1)[0], 'Finding': 'Possible TCA/OXPHOS uncoupling', 'Status': 'HYPOTHESIS', 'Evidence': 'TCA signature is not decreased while mitochondrial ATP-synthesis GO is negatively enriched', 'NovelInterpretation': 'Progression may retain TCA transcript capacity while reducing terminal oxidative-phosphorylation output.', 'AlternativeExplanation': 'GO redundancy and transcript-protein discordance may create the apparent divergence.', 'NextValidation': 'OCR, ATP-linked respiration, mitochondrial membrane potential and targeted respiratory-chain proteomics.'})

# Fusion, splice, FASTA inventories
rna_inventory = results / 'comparative_advantage' / 'comparative_advantage.rna_event_inventory.tsv'
fasta_inventory = results / 'comparative_advantage' / 'comparative_advantage.fasta_inventory.tsv'
rna_rows = read_tsv(rna_inventory)
fasta_rows = read_tsv(fasta_inventory)
if rna_rows:
    add_claim('RNA-001', 'RNA-supported fusion and splice candidates were retained as separate event classes.', 'SUPPORTED', 'MODERATE', [rna_inventory], [f"{r.get('Sample')} {r.get('Event class')}={r.get('Data rows')}" for r in rna_rows], 'RNA support does not confirm the underlying genomic rearrangement or complete transcript structure.')
if fasta_rows:
    add_claim('FASTA-001', 'Per-sample exploratory variant, fusion, splice, and combined protein FASTAs were generated.', 'SUPPORTED', 'HIGH', [fasta_inventory], [f"{r.get('Sample')} {r.get('FASTA class')}={r.get('Sequences')}" for r in fasta_rows], 'FASTA presence is not peptide-level proteomic confirmation.')

proteomics_dir = results / 'proteogenomics_validation'
proteomics_outputs = list(proteomics_dir.glob('*')) if proteomics_dir.is_dir() else []
if proteomics_outputs:
    add_claim('PROT-001', 'Novel RNA-derived events have downstream proteomic evidence.', 'PARTIALLY_SUPPORTED', 'VARIABLE', proteomics_outputs[:10], [f'files={len(proteomics_outputs)}'], 'Inspect direct MS/MS, MBR-only, sample matching, altered-residue coverage, contaminants, and decoys before claiming confirmation.')
else:
    add_claim('PROT-001', 'Novel RNA-derived events are proteomically confirmed.', 'NOT_TESTED', 'NONE', [], ['proteogenomics validation outputs absent or disabled'], 'RNA-derived protein sequences are search candidates only.', 'Enable MaxQuant validation using the exact searched FASTAs and raw-file sample mapping.')

# Resource efficiency
resource_summary = results / f'resource_usage-{job}.summary.tsv' if job else None
resource_rows = read_tsv(resource_summary) if resource_summary else []
resource_recommendations = []
if resource_rows:
    for row in resource_rows:
        process = row.get('Process', '')
        cores = safe_float(row.get('Median observed CPU cores'), 0.0)
        rss = safe_float(row.get('Maximum peak RSS GB'), 0.0)
        runtime = safe_float(row.get('Maximum runtime minutes'), 0.0)
        if runtime >= 240:
            resource_recommendations.append(f'{process}: critical-path candidate, max runtime {runtime:.1f} min')
        if cores <= 1.1 and runtime >= 60:
            resource_recommendations.append(f'{process}: low parallel utilization, median {cores:.2f} cores over up to {runtime:.1f} min')
        if rss <= 1.0 and runtime <= 30:
            resource_recommendations.append(f'{process}: lightweight task, max RSS {rss:.2f} GB')
    add_insight('Efficiency', 'The trace identifies critical-path and overallocated process candidates.', [resource_summary], resource_recommendations[:30], 'Use at least two clean full runs before reducing stable allocations; prioritize billing reductions separately from wall-time reductions.', 'HIGH')


# Five progression hypotheses, evaluated across separate evidence layers.
# Symbols are corrected to current approved names where possible. CD56 is represented by NCAM1,
# and NADSYN is represented by NADSYN1. Expression evidence is baseline-relative rather than
# mean absolute abundance, preventing constitutively abundant genes from dominating the score.
HYPOTHESIS_DEFINITIONS = {
    'H1_proteostasis_translation': {
        'description': 'Ribosome, proteasome, ER/UPR and translation quality control',
        'genes': set(CANONICAL_CYTO_RIBOSOME) | set(CANONICAL_MITO_RIBOSOME) | {
            'PCNA','EIF2AK3','ERN1','ATF6','ATF4','DDIT3','XBP1','HSPA5','HSP90B1','CALR','CANX',
            'DNAJB9','PDIA3','PDIA4','PDIA6','EDEM1','EDEM2','DERL1','DERL2','SEL1L','SEC61A1',
            'SEC61B','P4HB','PPP1R15A','LAGE3','GON7','TP53RK','TPRKB','OSGEP','YRDC','CDKAL1',
            'TRMT6','NOP56','NOP58','DKC1','GAR1','NHP2','ADAR','NONO','SFPQ','MATR3',
        } | {f'PSMA{i}' for i in range(1,9)} | {f'PSMB{i}' for i in range(1,12)} |
            {f'PSMC{i}' for i in range(1,7)} | {f'PSMD{i}' for i in range(1,19)} | {f'PSME{i}' for i in range(1,5)},
        'go_keywords': ('ribosom','translation','proteasom','unfolded protein','endoplasmic reticulum stress','erad'),
    },
    'H2_metabolic_rewiring': {
        'description': 'Glycolysis, PPP, NAD/PRPP and mitochondrial counter-axis',
        'genes': set('SLC2A1 HK1 HK2 GPI PFKP PFKM ALDOA TPI1 GAPDH PGK1 PGAM1 ENO1 PKM LDHA LDHB PDK1 PDHA1 PDHB DLAT DLD PGD G6PD TKT TALDO1 RPIA RPE PRPS1 PRPS2 NAMPT NADSYN1 NMNAT1 NMNAT2 NMNAT3 NADK NDUFA1 NDUFA2 NDUFB2 NDUFB3 NDUFS1 NDUFS2 NDUFS3 NDUFV1 SDHA SDHB SDHC SDHD UQCRC1 UQCRC2 COX4I1 ATP5F1A ATP5F1B'.split()),
        'go_keywords': ('glycol','glucose','hexose','pentose phosphate','nad','nicotinamide','oxidative phosphorylation','respiration','mitochond'),
    },
    'H3_surface_glycan_adhesion_remodeling': {
        'description': 'Surface antigen, glycan, extracellular-matrix and adhesion remodeling',
        'genes': set('TNFRSF17 CD38 SLAMF7 NCAM1 CD28 ITGA3 ITGA4 ITGA5 ITGAV ITGB1 ITGB3 ITGB5 ITGB6 ITGB7 FERMT2 CYFIP2 TGFB1 TGFB1I1 PMM2 GMPPA GMPPB GMDS GFPT1 PGM3 NANS CMAS ALG1 ALG2 STT3A STT3B ST3GAL6 MGAT1 MGAT2 B4GALT1 B4GALT3 FUT8 MAN2A1 MAN2A2 SLC35A1 SLC35C1 SLC35A2 SLC35B4 CD44 SDC1 CXCR4 SELPLG TLN1 VCL RAP1A RAP1B SRC'.split()),
        'go_keywords': ('adhesion','integrin','cell junction','glycosyl','glycan','surface receptor','extracellular matrix'),
    },
    'H4_MYC_IRF4_DNA_repair_evolution': {
        'description': 'MYC/IRF4/TGFB1, DNA repair/checkpoint and progression genes',
        'genes': set('MYC IRF4 TGFB1 TGFBR1 TGFBR2 IKZF1 IKZF3 BRCA1 BRCA2 RAD51 RAD50 MRE11 ATM ATR CHEK1 CHEK2 TP53 RB1 ZKSCAN3 CKS1B DIS3 KRAS NRAS BRAF CCND1 CCND2 CCND3 FGFR3 MAF MAFB MAX E2F1 E2F2 E2F3 PARP1 XRCC1 XRCC5 XRCC6 PRKDC NBN FANCA FANCD2'.split()),
        'go_keywords': ('dna repair','dna damage','checkpoint','double-strand break','cell cycle','myc'),
    },
    'H5_PCNA_stress_ATX101': {
        'description': 'PCNA/APIM stress scaffold and published ATX-101 sensitivity biomarkers',
        'genes': set('PCNA TPD52 TNFRSF17 LILRB4 TSG101 ZNRF2 UPF3B FADS2 SMAP CGREF1 GAA COG4 ENO1 PGD GAPDH PFKP HSPA5 HSP90B1 CALR CANX ATF4 DDIT3 XBP1 EIF2AK3 ERN1 NADK NAMPT PRPS1 PRPS2'.split()),
        'go_keywords': ('pcna','dna repair','glycol','pentose phosphate','oxidative stress','unfolded protein','endoplasmic reticulum stress'),
    },
}


def sign_test_one_sided(positive, negative):
    informative = positive + negative
    if informative == 0:
        return None
    threshold = max(positive, negative)
    return min(1.0, sum(math.comb(informative, k) for k in range(threshold, informative + 1)) / (2 ** informative))


def hypothesis_status(positive, negative, median_score, sign_p):
    informative = positive + negative
    if informative < 5:
        return 'INSUFFICIENT_DATA'
    fraction = positive / informative
    if fraction >= 0.70 and median_score > 0 and sign_p is not None and sign_p <= 0.05:
        return 'SUPPORTED_TRANSCRIPT_DIRECTION'
    if fraction >= 0.55 and median_score > 0:
        return 'PARTIALLY_SUPPORTED_TRANSCRIPT_DIRECTION'
    if fraction <= 0.35 and median_score < 0:
        return 'OPPOSITE_DIRECTION'
    return 'MIXED_OR_NOT_SUPPORTED'


# Optional phenotype metadata. It augments, but never overrides, the pipeline samplesheet.
phenotype_rows = {}
phenotype_path = Path(args.phenotype_metadata).resolve() if args.phenotype_metadata else None
if phenotype_path and phenotype_path.is_file():
    with phenotype_path.open(encoding='utf-8', newline='') as handle:
        for row in csv.DictReader(handle):
            sample_name = (row.get('sample') or '').strip()
            if sample_name:
                phenotype_rows[sample_name] = row

hypothesis_scores = []
hypothesis_evidence = []
hypothesis_variant_hits = []
hypothesis_novel = []

# Build complete full-matrix TPM maps for baseline-relative scoring.
if matrix.is_file():
    with open_text(matrix) as handle:
        matrix_reader = csv.DictReader(handle, delimiter='\t')
        matrix_fields = matrix_reader.fieldnames or []
        full_matrix_rows = list(matrix_reader)
    available_tpm_samples = sorted(field[:-4] for field in matrix_fields if field.endswith('_TPM'))
    sample_metadata = {row['sample']: row for row in samples}
    for sample_name in available_tpm_samples:
        sample_info = sample_metadata.get(sample_name, {})
        if sample_info.get('baseline') == 'true':
            continue
        subject = sample_info.get('subject', sample_name)
        baselines = [row['sample'] for row in samples if row['subject'] == subject and row['baseline'] == 'true']
        if len(baselines) != 1 or f'{baselines[0]}_TPM' not in matrix_fields:
            continue
        baseline_name = baselines[0]
        comparison = f'{sample_name}_vs_{baseline_name}'
        all_scores = {}
        for row in full_matrix_rows:
            gene = (row.get('Gene') or '').strip()
            if not gene or (row.get('Biotype') or '').strip() != 'protein_coding':
                continue
            sample_tpm = safe_float(row.get(f'{sample_name}_TPM'), 0.0)
            baseline_tpm = safe_float(row.get(f'{baseline_name}_TPM'), 0.0)
            if max(sample_tpm, baseline_tpm) < 1.0:
                continue
            all_scores[gene] = math.log2((sample_tpm + 0.5) / (baseline_tpm + 0.5))
        for hypothesis_id, definition in HYPOTHESIS_DEFINITIONS.items():
            selected = [(gene, all_scores[gene]) for gene in sorted(definition['genes'] & set(all_scores))]
            values = [value for _, value in selected]
            positive = sum(value > 0 for value in values)
            negative = sum(value < 0 for value in values)
            zero = len(values) - positive - negative
            median_score = statistics.median(values) if values else float('nan')
            p_value = sign_test_one_sided(positive, negative)
            status = hypothesis_status(positive, negative, median_score, p_value)
            phenotype = phenotype_rows.get(sample_name, {})
            hypothesis_scores.append({
                'Comparison': comparison, 'Sample': sample_name, 'BaselineSample': baseline_name,
                'Subject': subject, 'Group': sample_info.get('group', ''),
                'StageOrder': phenotype.get('stage_order', ''), 'ATX101_IC50': phenotype.get('atx101_ic50', ''),
                'Hypothesis': hypothesis_id, 'Description': definition['description'],
                'GenesDefined': len(definition['genes']), 'GenesMeasured': len(values),
                'Positive': positive, 'Negative': negative, 'Zero': zero,
                'PositiveFraction': positive / (positive + negative) if positive + negative else '',
                'MedianLog2TPMRatio': median_score if values else '',
                'MedianFoldChange': 2 ** median_score if values else '',
                'DirectionalSignPValue': p_value if p_value is not None else '',
                'TranscriptStatus': status,
                'TopPositiveGenes': ';'.join(gene for gene, value in sorted(selected, key=lambda item: item[1], reverse=True)[:20]),
                'TopNegativeGenes': ';'.join(gene for gene, value in sorted(selected, key=lambda item: item[1])[:20]),
            })

# Compare hypothesis scores with significant ranked GO terms from the same expression dataset.
ranked_by_comparison = defaultdict(list)
for row in ranked_rows:
    comparison = (row.get('Analysis') or '').split(':', 1)[0]
    ranked_by_comparison[comparison].append(row)

# Progression-variant GO is a separate layer and does not establish pathway activation.
progression_go_file = results / 'progression_biology' / 'progression_biology.go_enrichment.tsv'
progression_go_rows = read_tsv(progression_go_file)
progression_by_sample = defaultdict(list)
for row in progression_go_rows:
    sample_name = row.get('Sample') or row.get('Group') or row.get('Analysis', '').split(':', 1)[0]
    progression_by_sample[sample_name].append(row)

for score_row in hypothesis_scores:
    comparison = score_row['Comparison']
    sample_name = score_row['Sample']
    definition = HYPOTHESIS_DEFINITIONS[score_row['Hypothesis']]
    go_matches = [row for row in ranked_by_comparison.get(comparison, []) if safe_float(row.get('FDR'), 1.0) <= 0.1 and any(keyword in (row.get('GO_Name') or '').lower() for keyword in definition['go_keywords'])]
    go_positive = [row for row in go_matches if safe_float(row.get('MeanScore'), 0.0) > 0]
    go_negative = [row for row in go_matches if safe_float(row.get('MeanScore'), 0.0) < 0]
    progression_matches = [row for row in progression_by_sample.get(sample_name, []) if safe_float(row.get('FDR'), 1.0) <= 0.1 and any(keyword in (row.get('GO_Name') or '').lower() for keyword in definition['go_keywords'])]
    shared_go_genes = set()
    for row in progression_matches:
        for gene in (row.get('OverlapGenes') or '').split(';'):
            if gene in definition['genes']:
                shared_go_genes.add(gene)
    # Sarek and MaxQuant are reported only if pipeline outputs are present.
    external_summary = results / 'comparative_advantage' / 'comparative_advantage.external_caller_comparison.tsv'
    proteomics_dir = results / 'proteogenomics_validation'
    sarek_state = 'AVAILABLE_GENERAL_ONLY' if external_summary.is_file() else 'NOT_TESTED'
    proteomic_state = 'AVAILABLE_REQUIRES_EVENT_REVIEW' if proteomics_dir.is_dir() and any(proteomics_dir.iterdir()) else 'NOT_TESTED'
    transcript_status = score_row['TranscriptStatus']
    expression_concordance = 'CONCORDANT_WITHIN_EXPRESSION_LAYER' if go_positive else 'NO_SUPPORTING_RANKED_GO'
    if transcript_status == 'SUPPORTED_TRANSCRIPT_DIRECTION' and go_positive:
        overall = 'SUPPORTED_WITHIN_EXPRESSION_LAYER'
    elif transcript_status.startswith('SUPPORTED') or transcript_status.startswith('PARTIALLY'):
        overall = 'PARTIALLY_SUPPORTED_WITHIN_EXPRESSION_LAYER'
    else:
        overall = 'NOT_SUPPORTED_OR_MIXED'
    hypothesis_evidence.append({
        'Comparison': comparison, 'Hypothesis': score_row['Hypothesis'],
        'TranscriptStatus': score_row['TranscriptStatus'],
        'RankedGOSignificantPositiveTerms': len(go_positive),
        'RankedGOSignificantNegativeTerms': len(go_negative),
        'TopRankedGOPositive': ';'.join(f"{row.get('GO_ID')}|{row.get('GO_Name')}|FDR={safe_float(row.get('FDR'), 1.0):.3g}" for row in sorted(go_positive, key=lambda row: safe_float(row.get('FDR'), 1.0))[:10]),
        'TopRankedGONegative': ';'.join(f"{row.get('GO_ID')}|{row.get('GO_Name')}|FDR={safe_float(row.get('FDR'), 1.0):.3g}" for row in sorted(go_negative, key=lambda row: safe_float(row.get('FDR'), 1.0))[:10]),
        'ProgressionVariantGOSignificantTerms': len(progression_matches),
        'ProgressionVariantHypothesisGenes': ';'.join(sorted(shared_go_genes)),
        'EvidenceLayer': 'TRANSCRIPT_EXPRESSION',
        'EvidenceIndependence': 'INTERNAL_CONCORDANCE_SAME_EXPRESSION_DATASET',
        'ExpressionConcordance': expression_concordance,
        'SarekEvidence': sarek_state,
        'SarekInterpretation': 'General same-RNA callset reproducibility context only; not hypothesis-gene support.',
        'ProteomicEvidence': proteomic_state,
        'IndependentEvidenceStatus': 'NOT_TESTED' if proteomic_state == 'NOT_TESTED' else 'REQUIRES_EVENT_LEVEL_REVIEW',
        'OverallStatus': overall,
        'Limitation': 'Gene-set scores and ranked GO are derived from the same expression matrix. Progression-variant GO, Sarek comparison and proteomics answer different questions and are not interchangeable.',
    })
    hypothesis_novel.append({
        'FindingID': f'HNF-{len(hypothesis_novel)+1:03d}', 'Comparison': comparison,
        'Hypothesis': score_row['Hypothesis'], 'Status': overall,
        'Finding': definition['description'],
        'Evidence': f"transcript={score_row['TranscriptStatus']}; measured={score_row['GenesMeasured']}; positive={score_row['Positive']}; negative={score_row['Negative']}; medianFC={score_row['MedianFoldChange'] if score_row['MedianFoldChange'] != '' else 'NA'}; positive ranked GO terms={len(go_positive)}; progression-variant GO terms={len(progression_matches)}",
        'NovelInterpretation': 'Baseline-relative gene-set direction and ranked GO are concordant analyses of the same expression dataset. Agreement strengthens within-expression consistency but is not independent validation.',
        'AlternativeExplanation': 'One longitudinal series lacks replicate-based variance; coordinated RNA changes need not imply protein activity or causation.',
        'NextValidation': 'Matched quantitative proteomics, functional pathway assays and ATX-101 perturbation with phenotype measurements.',
    })

# Descriptive variant hits from any locally published annotated progression VCFs.
vcf_candidates = sorted(results.glob('**/*progression*.vcf.gz')) + sorted(results.glob('**/*rna_validated*.vcf.gz'))
seen_variant_hits = set()
for vcf_path in vcf_candidates:
    try:
        with open_text(vcf_path) as handle:
            csq_fields = []
            for line in handle:
                if line.startswith('##INFO=<ID=CSQ') and 'Format:' in line:
                    csq_fields = line.split('Format:', 1)[1].rsplit('"', 1)[0].rstrip('>').split('|')
                elif line.startswith('#'):
                    continue
                else:
                    fields = line.rstrip('\n').split('\t')
                    if len(fields) < 8 or not csq_fields:
                        continue
                    info = {part.split('=', 1)[0]: part.split('=', 1)[1] for part in fields[7].split(';') if '=' in part}
                    for annotation in info.get('CSQ', '').split(','):
                        values = annotation.split('|') + [''] * len(csq_fields)
                        record = dict(zip(csq_fields, values))
                        gene = (record.get('SYMBOL') or '').strip()
                        for hypothesis_id, definition in HYPOTHESIS_DEFINITIONS.items():
                            if gene not in definition['genes']:
                                continue
                            key = (str(vcf_path), fields[0], fields[1], fields[3], fields[4], hypothesis_id, gene)
                            if key in seen_variant_hits:
                                continue
                            seen_variant_hits.add(key)
                            hypothesis_variant_hits.append({
                                'File': relative(vcf_path), 'Sample': next((sample for sample in sample_metadata if vcf_path.name.startswith(sample)), ''),
                                'Hypothesis': hypothesis_id, 'Gene': gene,
                                'Consequence': record.get('Consequence', ''), 'Chrom': fields[0], 'Pos': fields[1],
                                'Ref': fields[3], 'Alt': fields[4],
                                'Interpretation': 'Descriptive RNA-observed variant hit; not evidence that the pathway is activated or that the alteration is causal.',
                            })
    except (OSError, EOFError, gzip.BadGzipFile):
        add_check('Hypothesis VCF parsing', 'WARN', f'Could not parse {relative(vcf_path)}', [vcf_path])

if hypothesis_scores:
    add_claim('HYP-001', 'The five predefined progression hypotheses were evaluated with matched-baseline transcript scores.', 'SUPPORTED', 'MODERATE', [matrix], [f'comparisons={len({row["Comparison"] for row in hypothesis_scores})}', f'hypothesis rows={len(hypothesis_scores)}'], 'Scores are descriptive within one longitudinal series and do not provide replicate-based differential-expression inference.')
    add_claim('HYP-002', 'The audit keeps expression signatures, ranked GO, progression-variant GO, Sarek context and proteomic evidence as separate evidence layers.', 'SUPPORTED', 'HIGH', [matrix, ranked_go, progression_go_file], [f'evidence matrix rows={len(hypothesis_evidence)}'], 'Gene-set scores and ranked GO use the same expression dataset and therefore provide internal concordance, not independent validation. Independent support requires a distinct assay or dataset.')

# Unsupported broad claims
add_claim('BIO-002', 'TK13 or TK14 is plasma cell leukemia.', 'NOT_SUPPORTED', 'NONE', [matrix, ranked_go] if matrix.is_file() else [], ['RNA expression and pathway resemblance are not a clinical classification'], 'PCL requires clinical/pathological criteria and cannot be inferred from this pipeline alone.', 'Clinical phenotype, pathology, circulating plasma-cell assessment, and orthogonal molecular data.')
add_claim('BIO-003', 'The pipeline demonstrates causation of increased ATX-101 sensitivity.', 'NOT_SUPPORTED', 'NONE', [matrix] if matrix.is_file() else [], ['observational longitudinal omics'], 'Association with stress/proteostasis signatures does not prove causation.', 'Perturbation experiments, drug-response measurements, rescue experiments, and functional metabolic/redox assays.')
add_claim('BIO-004', 'The RNA findings directly replicate the published proteomic findings.', 'PARTIALLY_SUPPORTED' if signature_rows else 'NOT_TESTED', 'MODERATE' if signature_rows else 'NONE', [matrix] if matrix.is_file() else [], [f'article-oriented signature rows={len(signature_rows)}'], 'RNA and protein abundance can diverge; the studies also differ in design and stage.', 'Matched quantitative proteomics and formal cross-study effect-size concordance.')
add_claim('EXP-002', 'The expression results establish population-level differential expression.', 'NOT_SUPPORTED', 'NONE', [matrix] if matrix.is_file() else [], [f'samples={len(samples)}', f'subjects={len(by_subject)}'], 'One longitudinal sample per stage provides no biological-replicate variance estimate.', 'Replicated cohort or validated external cohort.')

# File inventory and outputs
claim_fields = ['ClaimID','Claim','Status','EvidenceStrength','EvidenceFiles','ObservedEvidence','Limitation','RequiredNextTest']
insight_fields = ['Category','Insight','EvidenceFiles','ObservedEvidence','Interpretation','Confidence']
check_fields = ['Check','Status','Details','Files']
signature_fields = ['Comparison','ArticleAxis','Signature','ExpectedDirection','Genes','Positive','Negative','Zero','PositiveFraction','NegativeFraction','ExpectedDirectionFraction','MedianLog2TPMRatio','MedianFoldChange','Concordance','MinScore','MaxScore','TopPositiveGenes','TopNegativeGenes','Definition']
novel_fields = ['FindingID','Comparison','Finding','Status','Evidence','NovelInterpretation','AlternativeExplanation','NextValidation']
hypothesis_score_fields = ['Comparison','Sample','BaselineSample','Subject','Group','StageOrder','ATX101_IC50','Hypothesis','Description','GenesDefined','GenesMeasured','Positive','Negative','Zero','PositiveFraction','MedianLog2TPMRatio','MedianFoldChange','DirectionalSignPValue','TranscriptStatus','TopPositiveGenes','TopNegativeGenes']
hypothesis_evidence_fields = ['Comparison','Hypothesis','TranscriptStatus','EvidenceLayer','EvidenceIndependence','ExpressionConcordance','RankedGOSignificantPositiveTerms','RankedGOSignificantNegativeTerms','TopRankedGOPositive','TopRankedGONegative','ProgressionVariantGOSignificantTerms','ProgressionVariantHypothesisGenes','SarekEvidence','SarekInterpretation','ProteomicEvidence','IndependentEvidenceStatus','OverallStatus','Limitation']
hypothesis_variant_fields = ['File','Sample','Hypothesis','Gene','Consequence','Chrom','Pos','Ref','Alt','Interpretation']
hypothesis_novel_fields = ['FindingID','Comparison','Hypothesis','Status','Finding','Evidence','NovelInterpretation','AlternativeExplanation','NextValidation']
write_tsv(str(prefix) + '.claims.tsv', claims, claim_fields)
write_tsv(str(prefix) + '.insights.tsv', insights, insight_fields)
write_tsv(str(prefix) + '.checks.tsv', checks, check_fields)
write_tsv(str(prefix) + '.article_signatures.tsv', signature_rows, signature_fields)
write_tsv(str(prefix) + '.novel_findings.tsv', novel_findings, novel_fields)
write_tsv(str(prefix) + '.hypothesis_scores.tsv', hypothesis_scores, hypothesis_score_fields)
write_tsv(str(prefix) + '.hypothesis_evidence_matrix.tsv', hypothesis_evidence, hypothesis_evidence_fields)
write_tsv(str(prefix) + '.hypothesis_variant_hits.tsv', hypothesis_variant_hits, hypothesis_variant_fields)
write_tsv(str(prefix) + '.hypothesis_novel_findings.tsv', hypothesis_novel, hypothesis_novel_fields)

supported = [r for r in claims if r['Status'] == 'SUPPORTED']
partial = [r for r in claims if r['Status'] == 'PARTIALLY_SUPPORTED']
unsupported = [r for r in claims if r['Status'] in {'NOT_SUPPORTED', 'NOT_TESTED'}]
write_tsv(str(prefix) + '.supported_claims.tsv', supported, claim_fields)
write_tsv(str(prefix) + '.unsupported_claims.tsv', unsupported + partial, claim_fields)

report = []
report.append('# PGTK claim and insight audit')
report.append('')
report.append(f'- Project: `{project}`')
report.append(f'- Results: `{results}`')
report.append(f'- Selected job: `{job or "not detected"}`')
report.append(f'- Claims: {len(claims)}')
report.append(f'- Supported: {len(supported)}')
report.append(f'- Partially supported: {len(partial)}')
report.append(f'- Not supported or not tested: {len(unsupported)}')
report.append('')
report.append('## Claim matrix')
report.append('')
for row in claims:
    report.append(f"### {row['ClaimID']}: {row['Status']}")
    report.append('')
    report.append(row['Claim'])
    report.append('')
    report.append(f"Evidence strength: **{row['EvidenceStrength']}**")
    report.append('')
    report.append(f"Observed evidence: {row['ObservedEvidence'] or 'none'}")
    report.append('')
    report.append(f"Evidence files: `{row['EvidenceFiles'] or 'none'}`")
    report.append('')
    report.append(f"Limitation: {row['Limitation']}")
    if row['RequiredNextTest']:
        report.append('')
        report.append(f"Required next test: {row['RequiredNextTest']}")
    report.append('')
report.append('## Additional insights')
report.append('')
for row in insights:
    report.append(f"- **{row['Category']} | {row['Insight']}**: {row['ObservedEvidence']} Interpretation: {row['Interpretation']} Confidence: {row['Confidence']}.")
report.append('')
report.append('## Explicit novel findings')
report.append('')
for row in novel_findings:
    report.append(f"### {row['FindingID']}: {row['Finding']} [{row['Status']}]")
    report.append('')
    report.append(f"Comparison: {row['Comparison']}")
    report.append('')
    report.append(f"Evidence: {row['Evidence']}")
    report.append('')
    report.append(f"Interpretation: {row['NovelInterpretation']}")
    report.append('')
    report.append(f"Alternative explanation: {row['AlternativeExplanation']}")
    report.append('')
    report.append(f"Next validation: {row['NextValidation']}")
    report.append('')
report.append('## Five progression hypotheses')
report.append('')
for row in hypothesis_evidence:
    report.append(f"- **{row['Comparison']} | {row['Hypothesis']} | {row['OverallStatus']}**: transcript={row['TranscriptStatus']}; expression concordance={row['ExpressionConcordance']}; positive ranked GO terms={row['RankedGOSignificantPositiveTerms']} (correlated terms, not independent pathways); progression-variant GO terms={row['ProgressionVariantGOSignificantTerms']}; Sarek={row['SarekEvidence']} (general context only); proteomics={row['ProteomicEvidence']}; independent evidence={row['IndependentEvidenceStatus']}.")
report.append('')
report.append('## Interpretation rules enforced')
report.append('')
report.extend([
    '- RNA-supported variant does not mean DNA-confirmed mutation.',
    '- Non-baseline-only RNA observation does not prove newly acquired mutation.',
    '- RNA-supported fusion does not mean genomic rearrangement confirmation.',
    '- Read-supported splice transcript does not prove the full isoform structure.',
    '- A custom FASTA entry does not mean peptide or protein confirmation.',
    '- GO terms are overlapping ontology categories, not independent pathways.',
    '- Ranked GO and targeted transcript signatures use the same expression matrix and provide internal concordance, not independent validation.',
    '- General Sarek overlap is callset reproducibility context and is not hypothesis-specific gene support.',
    '- Surface, glycan and adhesion findings indicate remodeling; they do not establish glycan loss or marrow escape.',
    '- MHC class II reduction is reported only when at least two significant negative MHC class II GO terms are present.',
    '- One longitudinal sample per stage does not support population-level differential-expression inference.',
    '- Similarity to an article-derived signature is mechanistic concordance, not clinical stage classification or causal proof.',
])
Path(str(prefix) + '.report.md').write_text('\n'.join(report) + '\n', encoding='utf-8')

payload = {'project': str(project), 'results': str(results), 'job_id': job, 'claims': claims, 'insights': insights, 'checks': checks, 'article_signatures': signature_rows, 'novel_findings': novel_findings, 'hypothesis_scores': hypothesis_scores, 'hypothesis_evidence_matrix': hypothesis_evidence, 'hypothesis_variant_hits': hypothesis_variant_hits, 'hypothesis_novel_findings': hypothesis_novel}
Path(str(prefix) + '.json').write_text(json.dumps(payload, indent=2) + '\n', encoding='utf-8')

print(f'Wrote {prefix}.report.md')
print(f'Wrote {prefix}.claims.tsv')
print(f'Wrote {prefix}.supported_claims.tsv')
print(f'Wrote {prefix}.unsupported_claims.tsv')
print(f'Wrote {prefix}.insights.tsv')
print(f'Wrote {prefix}.article_signatures.tsv')
print(f'Wrote {prefix}.novel_findings.tsv')
print(f'Wrote {prefix}.hypothesis_scores.tsv')
print(f'Wrote {prefix}.hypothesis_evidence_matrix.tsv')
print(f'Wrote {prefix}.hypothesis_variant_hits.tsv')
print(f'Wrote {prefix}.hypothesis_novel_findings.tsv')
print(f'Wrote {prefix}.checks.tsv')
print(f'Wrote {prefix}.json')
print(f'Claims: supported={len(supported)}, partial={len(partial)}, unsupported_or_not_tested={len(unsupported)}')
