#!/usr/bin/env bash
# =============================================================================
# group_exclusive_mutations.sh
# Find genes hit by SVs or SNVs/indels in ONLY ONE patient group.
#
# Requirements: bcftools, python3  (no bedtools needed)
# =============================================================================
set -euo pipefail

RESULTS_DIR="results"
SV_DIR="${RESULTS_DIR}/delly_final"
SNV_DIR="${RESULTS_DIR}/mutect2_merged"
GTF="Homo_sapiens.GRCh38.110.gtf"
OUT_DIR="${RESULTS_DIR}/analysis"
TMPDIR_WORK=$(mktemp -d)
trap "rm -rf ${TMPDIR_WORK}" EXIT

mkdir -p "${OUT_DIR}"

GROUPS=($(ls "${SV_DIR}"/TK*.final.sv.vcf.gz \
  | xargs -I{} basename {} .final.sv.vcf.gz | sort))
echo "Groups found: ${GROUPS[*]}"

echo ""
echo "[1/5] Building gene interval index from GTF..."
python3 << PYEOF
import re, sys

gtf_file  = "${GTF}"
out_file  = "${TMPDIR_WORK}/genes.tsv"
valid_chr = set([str(i) for i in range(1,23)] + ['X','Y'])
genes = []
with open(gtf_file) as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 9 or parts[2] != 'gene':
            continue
        chrom = parts[0].lstrip('chr')
        if chrom not in valid_chr:
            continue
        start = int(parts[3]) - 1
        end   = int(parts[4])
        attrs = parts[8]
        gid   = re.search(r'gene_id "([^"]+)"', attrs)
        gname = re.search(r'gene_name "([^"]+)"', attrs)
        gid   = gid.group(1)   if gid   else ''
        gname = gname.group(1) if gname else gid
        genes.append((chrom, start, end, gid, gname))
genes.sort(key=lambda x: (x[0], x[1]))
with open(out_file, 'w') as fh:
    for g in genes:
        fh.write('\t'.join([g[0], str(g[1]), str(g[2]), g[3], g[4]]) + '\n')
print(f"  {len(genes)} gene intervals written")
PYEOF

echo ""
echo "[2/5] Extracting PASS variant intervals from VCFs..."

for group in "${GROUPS[@]}"; do
  vcf="${SV_DIR}/${group}.final.sv.vcf.gz"
  bcftools view -f PASS "${vcf}" \
    | bcftools query -f '%CHROM\t%POS\t%INFO/SVTYPE\t%INFO/END\n' \
    | awk '{chr=$1; gsub(/^chr/,"",chr); pos=$2+0; svtype=$3; endpos=$4+0;
            if(svtype=="BND"||endpos<pos) endpos=pos;
            print chr"\t"(pos-1)"\t"endpos"\t"svtype}' \
    | sort -k1,1 -k2,2n > "${TMPDIR_WORK}/sv_${group}.bed"
  echo "  SV  ${group}: $(wc -l < "${TMPDIR_WORK}/sv_${group}.bed") PASS variants"
done

for group in "${GROUPS[@]}"; do
  vcf="${SNV_DIR}/${group}.merged.snv.vcf.gz"
  [[ -f "${vcf}" ]] || { echo "  SNV ${group}: VCF not found, skipping"; continue; }
  bcftools view -f PASS "${vcf}" \
    | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' \
    | awk '{chr=$1; gsub(/^chr/,"",chr); pos=$2+0; ref=$3; alt=$4;
            endpos=pos+length(ref)-1;
            vtype=(length(ref)==1&&length(alt)==1)?"SNV":"INDEL";
            print chr"\t"(pos-1)"\t"endpos"\t"vtype}' \
    | sort -k1,1 -k2,2n > "${TMPDIR_WORK}/snv_${group}.bed"
  echo "  SNV ${group}: $(wc -l < "${TMPDIR_WORK}/snv_${group}.bed") PASS variants"
done

echo ""
echo "[3/5] Intersecting variants with gene intervals (pure Python)..."

GROUPS_STR="${GROUPS[*]}"
python3 << PYEOF
import os, bisect, csv
from collections import defaultdict

tmpdir   = "${TMPDIR_WORK}"
groups   = "${GROUPS_STR}".split()
out_dir  = "${OUT_DIR}"

def load_genes(genes_file):
    chrom_genes = defaultdict(list)
    with open(genes_file) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            c, s, e, gid, gname = parts[0], int(parts[1]), int(parts[2]), parts[3], parts[4]
            chrom_genes[c].append((s, e, gid, gname))
    for c in chrom_genes:
        chrom_genes[c].sort()
    return chrom_genes

def intersect(bed_file, chrom_genes):
    hits = set()
    if not os.path.exists(bed_file):
        return hits
    with open(bed_file) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4:
                continue
            chrom, vstart, vend, vtype = parts[0], int(parts[1]), int(parts[2]), parts[3]
            if vend == vstart:
                vend = vstart + 1
            if chrom not in chrom_genes:
                continue
            genes = chrom_genes[chrom]
            starts = [g[0] for g in genes]
            j = max(0, bisect.bisect_left(starts, vstart) - 1)
            while j < len(genes):
                gstart, gend, gid, gname = genes[j]
                if gstart >= vend:
                    break
                if gend > vstart:
                    hits.add((gid, gname, vtype))
                j += 1
    return hits

chrom_genes = load_genes(os.path.join(tmpdir, 'genes.tsv'))
print(f"  Loaded gene intervals for {len(chrom_genes)} chromosomes")

for vtype_label in ('sv', 'snv'):
    print(f"\n  -- {vtype_label.upper()} intersections --")
    gene_hits  = defaultdict(lambda: defaultdict(set))
    gene_names = {}

    for group in groups:
        bed_file = os.path.join(tmpdir, f'{vtype_label}_{group}.bed')
        hits = intersect(bed_file, chrom_genes)
        for gid, gname, vt in hits:
            gene_hits[gid][group].add(vt)
            gene_names[gid] = gname
        unique_genes = len({h[0] for h in hits})
        print(f"    {group}: {len(hits)} overlaps, {unique_genes} unique genes hit")

    matrix_file = os.path.join(out_dir, f'gene_{vtype_label}_matrix.tsv')
    with open(matrix_file, 'w', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow(['gene_id','gene_name'] + groups)
        for gid in sorted(gene_hits.keys()):
            row = [gid, gene_names[gid]]
            for g in groups:
                vts = gene_hits[gid].get(g, set())
                row.append(','.join(sorted(vts)) if vts else '0')
            writer.writerow(row)
    print(f"    Matrix: {len(gene_hits)} genes written to {matrix_file}")

    exclusive = []
    for gid in sorted(gene_hits.keys()):
        hit_groups = [g for g in groups if gene_hits[gid].get(g)]
        if len(hit_groups) == 1:
            g = hit_groups[0]
            exclusive.append({
                'gene_id':         gid,
                'gene_name':       gene_names[gid],
                'exclusive_group': g,
                'variant_types':   ','.join(sorted(gene_hits[gid][g])),
                'mutation_class':  vtype_label.upper()
            })

    excl_file = os.path.join(out_dir, f'{vtype_label}_exclusive_genes.tsv')
    with open(excl_file, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=['gene_id','gene_name','exclusive_group','variant_types','mutation_class'], delimiter='\t')
        writer.writeheader()
        writer.writerows(exclusive)
    print(f"    Exclusive genes ({vtype_label.upper()}): {len(exclusive)}")

    from collections import Counter
    c = Counter(r['exclusive_group'] for r in exclusive)
    for grp, cnt in sorted(c.items()):
        print(f"      {grp}: {cnt} exclusive genes")
PYEOF

echo ""
echo "[4/5] Writing combined summary..."
python3 << PYEOF
import csv
from collections import defaultdict

sv_file  = "${OUT_DIR}/sv_exclusive_genes.tsv"
snv_file = "${OUT_DIR}/snv_exclusive_genes.tsv"
out_file = "${OUT_DIR}/exclusive_genes_summary.tsv"

gene_data = defaultdict(lambda: {'gene_id':'','sv_groups':[],'sv_types':[],'snv_groups':[],'snv_types':[]})
with open(sv_file) as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        d = gene_data[r['gene_name']]
        d['gene_id'] = r['gene_id']
        d['sv_groups'].append(r['exclusive_group'])
        d['sv_types'].append(r['variant_types'])
with open(snv_file) as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        d = gene_data[r['gene_name']]
        if not d['gene_id']: d['gene_id'] = r['gene_id']
        d['snv_groups'].append(r['exclusive_group'])
        d['snv_types'].append(r['variant_types'])
with open(out_file, 'w', newline='') as fh:
    writer = csv.writer(fh, delimiter='\t')
    writer.writerow(['gene_name','gene_id','sv_exclusive_group','sv_types','snv_exclusive_group','snv_types'])
    for gname in sorted(gene_data.keys()):
        d = gene_data[gname]
        writer.writerow([gname, d['gene_id'],
            ';'.join(d['sv_groups']),  ';'.join(d['sv_types']),
            ';'.join(d['snv_groups']), ';'.join(d['snv_types'])])
print(f"  Combined summary: {len(gene_data)} unique genes")
PYEOF

echo ""
echo "[5/5] Preview — exclusive SV genes per group (top 15):"
python3 << PYEOF
import csv
from collections import defaultdict
rows = list(csv.DictReader(open("${OUT_DIR}/sv_exclusive_genes.tsv"), delimiter='\t'))
by_group = defaultdict(list)
for r in rows:
    by_group[r['exclusive_group']].append(r['gene_name'])
for grp in sorted(by_group.keys()):
    print(f"  {grp} ({len(by_group[grp])} genes): {', '.join(by_group[grp][:15])}")
PYEOF

echo ""
echo "Preview — exclusive SNV/INDEL genes per group (top 15):"
python3 << PYEOF
import csv
from collections import defaultdict
rows = list(csv.DictReader(open("${OUT_DIR}/snv_exclusive_genes.tsv"), delimiter='\t'))
by_group = defaultdict(list)
for r in rows:
    by_group[r['exclusive_group']].append(r['gene_name'])
for grp in sorted(by_group.keys()):
    print(f"  {grp} ({len(by_group[grp])} genes): {', '.join(by_group[grp][:15])}")
PYEOF

echo ""
echo "============================================================"
echo " DONE. Output files in ${OUT_DIR}/:"
echo "   gene_sv_matrix.tsv           SV gene x group presence matrix"
echo "   gene_snv_matrix.tsv          SNV gene x group presence matrix"
echo "   sv_exclusive_genes.tsv       Genes with SV in exactly 1 group"
echo "   snv_exclusive_genes.tsv      Genes with SNV/indel in exactly 1 group"
echo "   exclusive_genes_summary.tsv  Combined view"
echo "============================================================"
