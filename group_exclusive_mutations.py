#!/usr/bin/env python3
"""
group_exclusive_mutations.py
Find genes hit by SVs or SNVs/indels in ONLY ONE patient group.

Requirements: bcftools (on PATH), python3

Usage: python3 group_exclusive_mutations.py
Run from: /mnt/z/Download/TK/
"""

import os
import re
import csv
import bisect
import subprocess
import tempfile
import sys
from collections import defaultdict
from pathlib import Path

# ── Paths ────────────────────────────────────────────────────────────────────
BASE_DIR  = Path(".")
SV_DIR    = BASE_DIR / "results" / "delly_final"
SNV_DIR   = BASE_DIR / "results" / "mutect2_merged"
GTF       = BASE_DIR / "Homo_sapiens.GRCh38.110.gtf"
OUT_DIR   = BASE_DIR / "results" / "analysis"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Discover groups ──────────────────────────────────────────────────────────
GROUPS = sorted([
    p.name.replace(".final.sv.vcf.gz", "")
    for p in SV_DIR.glob("TK*.final.sv.vcf.gz")
])
print(f"Groups found: {GROUPS}")

# ── Step 1: Build gene interval index ────────────────────────────────────────
print("\n[1/5] Building gene interval index from GTF...")
valid_chr = set([str(i) for i in range(1, 23)] + ["X", "Y"])
chrom_genes: dict[str, list] = defaultdict(list)  # chrom -> [(start,end,gid,gname)]
gene_starts: dict[str, list] = {}                  # chrom -> [starts] for bisect

with open(GTF) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9 or parts[2] != "gene":
            continue
        chrom = parts[0].lstrip("chr")
        if chrom not in valid_chr:
            continue
        start = int(parts[3]) - 1   # 0-based
        end   = int(parts[4])
        attrs = parts[8]
        m_id   = re.search(r'gene_id "([^"]+)"', attrs)
        m_name = re.search(r'gene_name "([^"]+)"', attrs)
        gid    = m_id.group(1)   if m_id   else ""
        gname  = m_name.group(1) if m_name else gid
        chrom_genes[chrom].append((start, end, gid, gname))

for chrom in chrom_genes:
    chrom_genes[chrom].sort()
    gene_starts[chrom] = [g[0] for g in chrom_genes[chrom]]

total_genes = sum(len(v) for v in chrom_genes.values())
print(f"  {total_genes} gene intervals across {len(chrom_genes)} chromosomes")


# ── Step 2: Extract PASS variants via bcftools ───────────────────────────────
def extract_sv_bed(vcf_path: Path, out_path: Path):
    """Extract PASS SVs as BED-like file: chrom start end svtype"""
    result = subprocess.run(
        ["bcftools", "view", "-f", "PASS", str(vcf_path)],
        capture_output=True, text=True, check=True
    )
    rows = []
    for line in result.stdout.splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        chrom   = f[0].lstrip("chr")
        pos     = int(f[1])
        info    = f[7]
        svtype  = re.search(r'SVTYPE=([^;]+)', info)
        svend   = re.search(r'(?:^|;)END=([^;]+)', info)
        svtype  = svtype.group(1) if svtype else "UNK"
        endpos  = int(svend.group(1)) if svend else pos
        if svtype == "BND" or endpos < pos:
            endpos = pos + 1
        rows.append(f"{chrom}\t{pos-1}\t{endpos}\t{svtype}\n")
    with open(out_path, "w") as fh:
        fh.writelines(rows)
    return len(rows)


def extract_snv_bed(vcf_path: Path, out_path: Path):
    """Extract PASS SNVs/indels as BED-like file: chrom start end vtype"""
    result = subprocess.run(
        ["bcftools", "view", "-f", "PASS", str(vcf_path)],
        capture_output=True, text=True, check=True
    )
    rows = []
    for line in result.stdout.splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        chrom  = f[0].lstrip("chr")
        pos    = int(f[1])
        ref    = f[3]
        alt    = f[4]
        endpos = pos + len(ref) - 1
        vtype  = "SNV" if len(ref) == 1 and len(alt) == 1 else "INDEL"
        rows.append(f"{chrom}\t{pos-1}\t{endpos}\t{vtype}\n")
    with open(out_path, "w") as fh:
        fh.writelines(rows)
    return len(rows)


# ── Step 3: Intersect variants with gene intervals ───────────────────────────
def intersect_with_genes(bed_path: Path, chrom_genes, gene_starts):
    """Return set of (gene_id, gene_name, vtype) that overlap any variant."""
    hits: set = set()
    if not bed_path.exists():
        return hits
    with open(bed_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom, vstart, vend, vtype = parts[0], int(parts[1]), int(parts[2]), parts[3]
            if vend <= vstart:
                vend = vstart + 1
            if chrom not in chrom_genes:
                continue
            genes  = chrom_genes[chrom]
            starts = gene_starts[chrom]
            # Start from the first gene that might overlap (start < vend)
            j = max(0, bisect.bisect_left(starts, vstart) - 1)
            while j < len(genes):
                gstart, gend, gid, gname = genes[j]
                if gstart >= vend:
                    break
                if gend > vstart:
                    hits.add((gid, gname, vtype))
                j += 1
    return hits


print("\n[2/5] Extracting PASS variant intervals from VCFs...")
with tempfile.TemporaryDirectory() as tmpdir:
    tmpdir = Path(tmpdir)

    for group in GROUPS:
        sv_vcf  = SV_DIR  / f"{group}.final.sv.vcf.gz"
        snv_vcf = SNV_DIR / f"{group}.merged.snv.vcf.gz"
        n_sv  = extract_sv_bed(sv_vcf,   tmpdir / f"sv_{group}.bed")
        print(f"  SV  {group}: {n_sv} PASS SVs")
        if snv_vcf.exists():
            n_snv = extract_snv_bed(snv_vcf, tmpdir / f"snv_{group}.bed")
            print(f"  SNV {group}: {n_snv} PASS SNVs/indels")
        else:
            print(f"  SNV {group}: VCF not found, skipping")

    print("\n[3/5] Intersecting variants with gene intervals...")

    for vtype_label in ("sv", "snv"):
        print(f"\n  -- {vtype_label.upper()} --")
        gene_hits:  dict = defaultdict(lambda: defaultdict(set))  # gid -> group -> {vtypes}
        gene_names: dict = {}

        for group in GROUPS:
            bed_path = tmpdir / f"{vtype_label}_{group}.bed"
            hits = intersect_with_genes(bed_path, chrom_genes, gene_starts)
            for gid, gname, vt in hits:
                gene_hits[gid][group].add(vt)
                gene_names[gid] = gname
            unique = len({h[0] for h in hits})
            print(f"    {group}: {len(hits)} overlaps, {unique} unique genes hit")

        # Write gene x group matrix
        matrix_path = OUT_DIR / f"gene_{vtype_label}_matrix.tsv"
        with open(matrix_path, "w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(["gene_id", "gene_name"] + GROUPS)
            for gid in sorted(gene_hits.keys()):
                row = [gid, gene_names[gid]]
                for g in GROUPS:
                    vts = gene_hits[gid].get(g, set())
                    row.append(",".join(sorted(vts)) if vts else "0")
                writer.writerow(row)
        print(f"    Matrix ({len(gene_hits)} genes) -> {matrix_path}")

        # Find exclusive genes
        exclusive = []
        for gid in sorted(gene_hits.keys()):
            hit_groups = [g for g in GROUPS if gene_hits[gid].get(g)]
            if len(hit_groups) == 1:
                g = hit_groups[0]
                exclusive.append({
                    "gene_id":         gid,
                    "gene_name":       gene_names[gid],
                    "exclusive_group": g,
                    "variant_types":   ",".join(sorted(gene_hits[gid][g])),
                    "mutation_class":  vtype_label.upper(),
                })

        excl_path = OUT_DIR / f"{vtype_label}_exclusive_genes.tsv"
        with open(excl_path, "w", newline="") as fh:
            writer = csv.DictWriter(fh,
                fieldnames=["gene_id","gene_name","exclusive_group","variant_types","mutation_class"],
                delimiter="\t")
            writer.writeheader()
            writer.writerows(exclusive)
        print(f"    Exclusive genes: {len(exclusive)} -> {excl_path}")

        from collections import Counter
        c = Counter(r["exclusive_group"] for r in exclusive)
        for grp, cnt in sorted(c.items()):
            print(f"      {grp}: {cnt}")

# ── Step 4: Combined summary ─────────────────────────────────────────────────
print("\n[4/5] Writing combined summary...")

sv_rows  = list(csv.DictReader(open(OUT_DIR / "sv_exclusive_genes.tsv"),  delimiter="\t"))
snv_rows = list(csv.DictReader(open(OUT_DIR / "snv_exclusive_genes.tsv"), delimiter="\t"))

combined: dict = defaultdict(lambda: {"gene_id":"","sv_groups":[],"sv_types":[],"snv_groups":[],"snv_types":[]})
for r in sv_rows:
    d = combined[r["gene_name"]]
    d["gene_id"] = r["gene_id"]
    d["sv_groups"].append(r["exclusive_group"])
    d["sv_types"].append(r["variant_types"])
for r in snv_rows:
    d = combined[r["gene_name"]]
    if not d["gene_id"]: d["gene_id"] = r["gene_id"]
    d["snv_groups"].append(r["exclusive_group"])
    d["snv_types"].append(r["variant_types"])

summary_path = OUT_DIR / "exclusive_genes_summary.tsv"
with open(summary_path, "w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")
    writer.writerow(["gene_name","gene_id","sv_exclusive_group","sv_types","snv_exclusive_group","snv_types"])
    for gname in sorted(combined.keys()):
        d = combined[gname]
        writer.writerow([gname, d["gene_id"],
            ";".join(d["sv_groups"]),  ";".join(d["sv_types"]),
            ";".join(d["snv_groups"]), ";".join(d["snv_types"])])
print(f"  Combined summary: {len(combined)} unique genes -> {summary_path}")

# ── Step 5: Print preview ────────────────────────────────────────────────────
print("\n[5/5] Preview — exclusive SV genes per group (first 15):")
by_group: dict = defaultdict(list)
for r in sv_rows:
    by_group[r["exclusive_group"]].append(r["gene_name"])
for grp in sorted(by_group.keys()):
    print(f"  {grp} ({len(by_group[grp])} genes): {', '.join(by_group[grp][:15])}")

print("\nPreview — exclusive SNV/INDEL genes per group (first 15):")
by_group2: dict = defaultdict(list)
for r in snv_rows:
    by_group2[r["exclusive_group"]].append(r["gene_name"])
for grp in sorted(by_group2.keys()):
    print(f"  {grp} ({len(by_group2[grp])} genes): {', '.join(by_group2[grp][:15])}")

print("\n" + "="*60)
print(f"Output files in {OUT_DIR}/")
print("  gene_sv_matrix.tsv           SV presence matrix (genes x groups)")
print("  gene_snv_matrix.tsv          SNV presence matrix (genes x groups)")
print("  sv_exclusive_genes.tsv       Genes with SV in exactly 1 group")
print("  snv_exclusive_genes.tsv      Genes with SNV/indel in exactly 1 group")
print("  exclusive_genes_summary.tsv  Combined view (SV + SNV)")
print("="*60)
