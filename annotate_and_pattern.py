#!/usr/bin/env python3
"""
annotate_and_pattern.py

Annotate all GroupA/GroupB exclusive genes with:
  - Chromosome, start, end, strand, biotype, cytoband (arm + band)
Then run pattern analysis:
  1. Chromosomal distribution — are exclusive genes clustered on specific chromosomes?
  2. Chromosome arm bias (p vs q arm)
  3. Cytoband hotspots
  4. Biotype breakdown (protein_coding vs lncRNA vs pseudogene etc.)
  5. SV type breakdown (DEL / DUP / INV / BND / INS)
  6. Recurrence spectrum (how many samples share each exclusive gene)
"""

import csv, re, math
from pathlib import Path
from collections import defaultdict, Counter

OUT_DIR  = Path("results/analysis")
ANN_FILE = OUT_DIR / "gene_annotation_table.tsv"

# ── Load annotation table ────────────────────────────────────────────────────
ann = {}   # gene_id -> dict
with open(ANN_FILE) as fh:
    for row in csv.DictReader(fh, delimiter='\t'):
        ann[row['gene_id']] = row
# also index by gene_name for fallback
ann_by_name = {v['gene_name']: v for v in ann.values()}

# GRCh38 chromosome arm boundaries (centromere positions, Mb, approximate)
# Source: UCSC hg38 cytoBand
CENTROMERE = {
    '1': 123400000,  '2': 93900000,   '3': 90900000,
    '4': 50000000,   '5': 48800000,   '6': 59800000,
    '7': 60100000,   '8': 45200000,   '9': 43000000,
    '10': 39800000,  '11': 53400000,  '12': 35500000,
    '13': 17700000,  '14': 17200000,  '15': 19000000,
    '16': 36800000,  '17': 25100000,  '18': 18500000,
    '19': 26200000,  '20': 28100000,  '21': 12000000,
    '22': 15000000,  'X':  61000000,  'Y':  10400000,
}

def get_arm(chrom, pos):
    cen = CENTROMERE.get(chrom)
    if cen is None: return '?'
    return 'p' if pos < cen else 'q'

def get_cytoband_region(chrom, pos):
    """Rough cytoband: chrArm + 10Mb band number (e.g. 1p13)"""
    arm  = get_arm(chrom, pos)
    cen  = CENTROMERE.get(chrom, 0)
    if arm == 'p':
        dist = cen - pos
    else:
        dist = pos - cen
    band_num = int(dist / 10_000_000) + 1
    return f"{chrom}{arm}{band_num}"

def lookup_gene(gene_id, gene_name):
    a = ann.get(gene_id) or ann_by_name.get(gene_name)
    return a

# ── Load and annotate all exclusive files ────────────────────────────────────
FILES = {
    ('SV',  'GroupA'): OUT_DIR / 'cohort_structural_variants_(sv)_groupA_exclusive.tsv',
    ('SV',  'GroupB'): OUT_DIR / 'cohort_structural_variants_(sv)_groupB_exclusive.tsv',
    ('SNV', 'GroupA'): OUT_DIR / 'cohort_snv___indels_groupA_exclusive.tsv',
    ('SNV', 'GroupB'): OUT_DIR / 'cohort_snv___indels_groupB_exclusive.tsv',
}

all_annotated = []

for (vtype, group), fpath in FILES.items():
    rows = list(csv.DictReader(open(fpath), delimiter='\t'))
    for r in rows:
        a = lookup_gene(r['gene_id'], r['gene_name'])
        chrom = a['chrom']     if a else 'unknown'
        start = int(a['start']) if a else 0
        end   = int(a['end'])   if a else 0
        biotype = a['biotype'] if a else 'unknown'
        arm     = get_arm(chrom, start) if a else '?'
        cytoband = get_cytoband_region(chrom, start) if a else 'unknown'
        n_hit = int(r['n_groupA_samples']) if group == 'GroupA' else int(r['n_groupB_samples'])
        vtypes = r['groupA_vtypes'] if group == 'GroupA' else r['groupB_vtypes']
        all_annotated.append({
            'gene_id':   r['gene_id'],
            'gene_name': r['gene_name'],
            'vtype':     vtype,
            'group':     group,
            'chrom':     chrom,
            'start':     start,
            'end':       end,
            'arm':       arm,
            'cytoband':  cytoband,
            'biotype':   biotype,
            'n_samples': n_hit,
            'vtypes':    vtypes,
        })

# Write full annotated table
ann_out = OUT_DIR / 'exclusive_genes_annotated.tsv'
with open(ann_out, 'w', newline='') as fh:
    fields = ['gene_id','gene_name','vtype','group','chrom','start','end',
              'arm','cytoband','biotype','n_samples','vtypes']
    w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
    w.writeheader()
    w.writerows(all_annotated)
print(f"Annotated table written: {ann_out}  ({len(all_annotated)} rows)")

# ── Helper: pretty table ─────────────────────────────────────────────────────
def print_table(title, headers, rows, max_rows=30):
    print(f"\n{'─'*70}")
    print(f"  {title}")
    print(f"{'─'*70}")
    widths = [max(len(str(h)), max((len(str(r[i])) for r in rows[:max_rows]), default=0))
              for i, h in enumerate(headers)]
    fmt = '  ' + '  '.join(f'{{:<{w}}}' for w in widths)
    print(fmt.format(*headers))
    print('  ' + '  '.join('-'*w for w in widths))
    for row in rows[:max_rows]:
        print(fmt.format(*[str(x) for x in row]))
    if len(rows) > max_rows:
        print(f"  ... ({len(rows)-max_rows} more rows)")

# ── Analysis ─────────────────────────────────────────────────────────────────

def analyse(subset_label, subset):
    """Run all pattern analyses on a subset of annotated rows."""
    print(f"\n\n{'='*70}")
    print(f"  PATTERN ANALYSIS: {subset_label}")
    print(f"  n={len(subset)} gene-hits")
    print(f"{'='*70}")

    # filter out unknown chrom
    known = [r for r in subset if r['chrom'] not in ('unknown','')]

    # ── 1. Chromosomal distribution ──────────────────────────────────────────
    chr_order = [str(i) for i in range(1,23)] + ['X','Y']
    chrom_count = Counter(r['chrom'] for r in known)
    # normalise by expected gene count per chrom
    chrom_genes_expected = Counter(a['chrom'] for a in ann.values()
                                   if a['chrom'] in chr_order)
    chrom_rows = []
    for c in chr_order:
        n = chrom_count.get(c, 0)
        exp = chrom_genes_expected.get(c, 1)
        enrichment = round((n / exp) / (len(known) / len(ann)) , 2) if exp else 0
        chrom_rows.append((c, n, exp, enrichment))
    chrom_rows_sorted = sorted(chrom_rows, key=lambda x: -x[3])
    print_table("1. Chromosomal enrichment (observed / expected by gene density)",
                ['Chr','Exclusive_genes','Genes_on_chr','Enrichment_score'],
                chrom_rows_sorted, max_rows=24)

    # ── 2. Arm bias ──────────────────────────────────────────────────────────
    arm_count = Counter(r['arm'] for r in known)
    total_arm = sum(arm_count.values())
    arm_rows = [(arm, arm_count.get(arm,0),
                 f"{100*arm_count.get(arm,0)/total_arm:.1f}%")
                for arm in ['p','q','?']]
    print_table("2. Chromosome arm bias (p vs q)",
                ['Arm','Count','Fraction'], arm_rows)

    # ── 3. Cytoband hotspots (top 20) ────────────────────────────────────────
    band_count = Counter(r['cytoband'] for r in known)
    band_rows  = [(b, n) for b, n in band_count.most_common(20)]
    print_table("3. Cytoband hotspots (top 20)",
                ['Cytoband','Exclusive_genes'], band_rows)

    # ── 4. Biotype breakdown ─────────────────────────────────────────────────
    bio_count = Counter(r['biotype'] for r in known)
    bio_total = sum(bio_count.values())
    bio_rows  = [(b, n, f"{100*n/bio_total:.1f}%")
                 for b, n in bio_count.most_common(15)]
    print_table("4. Gene biotype breakdown",
                ['Biotype','Count','Fraction'], bio_rows)

    # ── 5. SV/SNV type breakdown ─────────────────────────────────────────────
    vtype_flat = []
    for r in subset:
        vtype_flat.extend(r['vtypes'].split(','))
    vtype_count = Counter(vtype_flat)
    vtype_rows  = [(t, n) for t, n in vtype_count.most_common()]
    print_table("5. Variant type breakdown",
                ['Variant_type','Count'], vtype_rows)

    # ── 6. Recurrence spectrum ───────────────────────────────────────────────
    rec_count = Counter(r['n_samples'] for r in subset)
    rec_rows  = [(n, rec_count.get(n,0)) for n in sorted(rec_count.keys(), reverse=True)]
    print_table("6. Recurrence (how many group-samples carry each exclusive gene)",
                ['N_samples','N_genes'], rec_rows)

    return chrom_rows_sorted, band_count

# ─────────────────────────────────────────────────────────────────────────────
# Run per vtype x group and then compare directly
for vtype in ('SV', 'SNV'):
    for group in ('GroupA', 'GroupB'):
        sub = [r for r in all_annotated if r['vtype']==vtype and r['group']==group]
        analyse(f"{vtype} — {group} exclusive", sub)

# ── Cross-group chromosome comparison ────────────────────────────────────────
print(f"\n\n{'='*70}")
print("  CROSS-GROUP CHROMOSOME COMPARISON (enrichment ratio GroupA vs GroupB)")
print("  >1 = more enriched in GroupA; <1 = more enriched in GroupB")
print(f"{'='*70}")

chr_order = [str(i) for i in range(1,23)] + ['X','Y']
for vtype in ('SV', 'SNV'):
    print(f"\n  {vtype}:")
    subA = [r for r in all_annotated if r['vtype']==vtype and r['group']=='GroupA']
    subB = [r for r in all_annotated if r['vtype']==vtype and r['group']=='GroupB']
    cA = Counter(r['chrom'] for r in subA if r['chrom'] in chr_order)
    cB = Counter(r['chrom'] for r in subB if r['chrom'] in chr_order)
    nA, nB = max(len(subA),1), max(len(subB),1)

    rows = []
    for c in chr_order:
        fa = cA.get(c,0) / nA
        fb = cB.get(c,0) / nB
        ratio = round(fa/fb, 2) if fb > 0 else ('A_only' if fa > 0 else '-')
        if cA.get(c,0) + cB.get(c,0) > 0:
            rows.append((c, cA.get(c,0), cB.get(c,0), ratio))
    rows_sorted = sorted(rows, key=lambda x: (
        -(x[3] if isinstance(x[3], float) else (999 if x[3]=='A_only' else 0))))
    print_table(f"  {vtype}: Chromosome enrichment ratio (GroupA/GroupB)",
                ['Chr','GroupA_n','GroupB_n','Ratio_A_vs_B'],
                rows_sorted, max_rows=24)

# ── Hotspot overlap between SV and SNV ───────────────────────────────────────
print(f"\n\n{'='*70}")
print("  CHROMOSOMAL HOTSPOTS: same chromosome exclusive in BOTH SV and SNV?")
print(f"{'='*70}")

for group in ('GroupA', 'GroupB'):
    sv_chroms  = Counter(r['chrom'] for r in all_annotated
                         if r['vtype']=='SV' and r['group']==group)
    snv_chroms = Counter(r['chrom'] for r in all_annotated
                         if r['vtype']=='SNV' and r['group']==group)
    overlap = set(sv_chroms) & set(snv_chroms) & set(chr_order)
    rows = []
    for c in chr_order:
        s = sv_chroms.get(c,0)
        n = snv_chroms.get(c,0)
        if s > 0 and n > 0:
            rows.append((c, s, n, s+n))
    rows_sorted = sorted(rows, key=lambda x: -x[3])
    print_table(f"\n  {group} — chromosomes with exclusive genes in BOTH SV and SNV",
                ['Chr','SV_exclusive','SNV_exclusive','Total'], rows_sorted, max_rows=24)

# ── 22q11 / 14q32 (IgH) / immunoglobulin locus check ────────────────────────
print(f"\n\n{'='*70}")
print("  SPECIAL LOCI: Immunoglobulin / known deletion hotspots")
print(f"{'='*70}")

ig_genes = [r for r in all_annotated
            if any(x in r['gene_name'] for x in ['IGL','IGH','IGK','IGLV','IGHV'])]
print(f"\n  Immunoglobulin locus genes in exclusive sets:")
ig_by_group = Counter((r['group'], r['vtype']) for r in ig_genes)
for (grp, vt), n in sorted(ig_by_group.items()):
    genes = [r['gene_name'] for r in ig_genes if r['group']==grp and r['vtype']==vt][:10]
    print(f"    {grp}/{vt}: n={n}  e.g. {', '.join(genes)}")

print(f"\n  22q11 (DiGeorge) region genes:")
q22 = [r for r in all_annotated if r['chrom']=='22' and r['start'] < 24_000_000]
for r in q22:
    print(f"    {r['gene_name']:<20} {r['group']}/{r['vtype']}  cytoband={r['cytoband']}")

print("\nDone. Full annotated table: results/analysis/exclusive_genes_annotated.tsv")
