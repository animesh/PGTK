#!/usr/bin/env python3
"""
group_compare_two_cohorts.py

Compare two meta-groups aligned to ATX-101 sensitivity (PMID 39682151):
  GroupA = Sensitive  = TK13, TK14, TK16   (ATX-101 hypersensitive)
  GroupB = Resistant  = TK10, TK12, TK18, TK91, TK92, TK93  (less sensitive)

For each variant type (SV, SNV/indel), find genes that are:
  - Exclusive to GroupA (hit in >=1 GroupA sample, 0 GroupB samples)
  - Exclusive to GroupB (hit in >=1 GroupB sample, 0 GroupA samples)
  - Shared (hit in both groups)

Uses the pre-built gene x group matrices from the previous analysis.
"""

import csv, sys
from collections import defaultdict
from pathlib import Path

OUT_DIR  = Path("results/analysis")
MATRIX_SV  = OUT_DIR / "gene_sv_matrix.tsv"
MATRIX_SNV = OUT_DIR / "gene_snv_matrix.tsv"

# Sensitive (ATX-101 hypersensitive) vs Resistant — per PMID 39682151
GROUP_A = {"TK13", "TK14", "TK16"}                        # Sensitive
GROUP_B = {"TK10", "TK12", "TK18", "TK91", "TK92", "TK93"}  # Resistant

def compare_groups(matrix_path: Path, label: str):
    print(f"\n{'='*60}")
    print(f"  {label} — GroupA vs GroupB comparison")
    print(f"  GroupA: {sorted(GROUP_A)}")
    print(f"  GroupB: {sorted(GROUP_B)}")
    print(f"{'='*60}")

    results = []
    with open(matrix_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        all_cols = [c for c in reader.fieldnames if c not in ("gene_id", "gene_name")]
        for row in reader:
            gid   = row["gene_id"]
            gname = row["gene_name"]
            hit_A = [g for g in GROUP_A if g in all_cols and row.get(g, "0") != "0"]
            hit_B = [g for g in GROUP_B if g in all_cols and row.get(g, "0") != "0"]

            vt_A = set()
            for g in hit_A:
                vt_A.update(row[g].split(","))
            vt_B = set()
            for g in hit_B:
                vt_B.update(row[g].split(","))

            if hit_A and not hit_B:
                category = "GroupA_exclusive"
            elif hit_B and not hit_A:
                category = "GroupB_exclusive"
            elif hit_A and hit_B:
                category = "Shared"
            else:
                continue

            results.append({
                "gene_id":        gid,
                "gene_name":      gname,
                "category":       category,
                "groupA_samples": ",".join(sorted(hit_A)),
                "groupA_vtypes":  ",".join(sorted(vt_A)),
                "groupB_samples": ",".join(sorted(hit_B)),
                "groupB_vtypes":  ",".join(sorted(vt_B)),
                "n_groupA_samples": len(hit_A),
                "n_groupB_samples": len(hit_B),
            })

    excl_A = [r for r in results if r["category"] == "GroupA_exclusive"]
    excl_B = [r for r in results if r["category"] == "GroupB_exclusive"]
    shared = [r for r in results if r["category"] == "Shared"]

    print(f"\n  Total genes with any variant: {len(results)}")
    print(f"  GroupA exclusive: {len(excl_A)}")
    print(f"  GroupB exclusive: {len(excl_B)}")
    print(f"  Shared (both):    {len(shared)}")

    tag = label.lower().replace("/", "_").replace(" ", "_")

    # Write exclusive GroupA (Sensitive)
    out_A = OUT_DIR / f"cohort_{tag}_sensitive_exclusive.tsv"
    with open(out_A, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(results[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(sorted(excl_A, key=lambda r: (-r["n_groupA_samples"], r["gene_name"])))

    # Write exclusive GroupB (Resistant)
    out_B = OUT_DIR / f"cohort_{tag}_resistant_exclusive.tsv"
    with open(out_B, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(results[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(sorted(excl_B, key=lambda r: (-r["n_groupB_samples"], r["gene_name"])))

    # Write shared
    out_S = OUT_DIR / f"cohort_{tag}_shared.tsv"
    with open(out_S, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(results[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(sorted(shared, key=lambda r: (-(r["n_groupA_samples"]+r["n_groupB_samples"]), r["gene_name"])))

    print(f"\n  Output:")
    print(f"    Sensitive exclusive -> {out_A}")
    print(f"    Resistant exclusive -> {out_B}")
    print(f"    Shared              -> {out_S}")

    # ── Preview: top hits in each category ──────────────────────────────────
    print(f"\n  Top Sensitive-exclusive genes (sorted by # Sensitive samples hit):")
    print(f"  {'Gene':<20} {'Sensitive samples':<35} {'SV/SNV types'}")
    print(f"  {'-'*75}")
    for r in sorted(excl_A, key=lambda r: -r["n_groupA_samples"])[:25]:
        print(f"  {r['gene_name']:<20} {r['groupA_samples']:<35} {r['groupA_vtypes']}")

    print(f"\n  Top Resistant-exclusive genes (sorted by # Resistant samples hit):")
    print(f"  {'Gene':<20} {'Resistant samples':<35} {'SV/SNV types'}")
    print(f"  {'-'*75}")
    for r in sorted(excl_B, key=lambda r: -r["n_groupB_samples"])[:25]:
        print(f"  {r['gene_name']:<20} {r['groupB_samples']:<35} {r['groupB_vtypes']}")

    # ── Recurrence: genes hit in ALL samples of one group ───────────────────
    n_A = len(GROUP_A)
    n_B = len(GROUP_B)
    excl_A_all = [r for r in excl_A if r["n_groupA_samples"] == n_A]
    excl_B_all = [r for r in excl_B if r["n_groupB_samples"] == n_B]

    if excl_A_all:
        print(f"\n  ★ Sensitive-exclusive genes hit in ALL {n_A} Sensitive samples:")
        for r in excl_A_all:
            print(f"    {r['gene_name']} ({r['groupA_vtypes']})")
    if excl_B_all:
        print(f"\n  ★ Resistant-exclusive genes hit in ALL {n_B} Resistant samples:")
        for r in excl_B_all:
            print(f"    {r['gene_name']} ({r['groupB_vtypes']})")

    # ── Genes hit in majority of one group ───────────────────────────────────
    maj_A = [r for r in excl_A if r["n_groupA_samples"] >= max(2, n_A // 2)]
    maj_B = [r for r in excl_B if r["n_groupB_samples"] >= max(2, n_B // 2)]

    print(f"\n  ★ Sensitive-exclusive genes hit in ≥{max(2,n_A//2)}/{n_A} Sensitive samples ({len(maj_A)} genes):")
    for r in sorted(maj_A, key=lambda r: -r["n_groupA_samples"])[:30]:
        print(f"    {r['gene_name']:<22} n={r['n_groupA_samples']}  samples={r['groupA_samples']}  types={r['groupA_vtypes']}")

    print(f"\n  ★ Resistant-exclusive genes hit in ≥{max(2,n_B//2)}/{n_B} Resistant samples ({len(maj_B)} genes):")
    for r in sorted(maj_B, key=lambda r: -r["n_groupB_samples"])[:30]:
        print(f"    {r['gene_name']:<22} n={r['n_groupB_samples']}  samples={r['groupB_samples']}  types={r['groupB_vtypes']}")

    return results

# Run for both variant types
compare_groups(MATRIX_SV,  "Structural Variants (SV)")
compare_groups(MATRIX_SNV, "SNV / Indels")

print("\n" + "="*60)
print("All done. Files written to results/analysis/")
print("  cohort_sv_sensitive_exclusive.tsv    (TK13, TK14, TK16)")
print("  cohort_sv_resistant_exclusive.tsv    (TK10, TK12, TK18, TK91, TK92, TK93)")
print("  cohort_sv_shared.tsv")
print("  cohort_snv___indels_sensitive_exclusive.tsv")
print("  cohort_snv___indels_resistant_exclusive.tsv")
print("  cohort_snv___indels_shared.tsv")
print("="*60)
