#!/usr/bin/env python3
"""
plot_coverage_ideogram.py

Chromosome ideogram with:
  - One line per replicate BAM coloured by sample group
  - Normalised read depth (reads per chrom / total mapped reads * 1e6 / chrom_length_Mb)
  - Cytoband ideogram beneath each chromosome panel
  - Exclusive gene density (GroupA / GroupB / Shared) shown as bar track
  - Per-sample lines (thin) + per-group mean (thick)
"""

import csv, re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from pathlib import Path
from collections import defaultdict

# ── Config ─────────────────────────────────────────────────────────────────
IDXSTATS   = Path("results/coverage/all_idxstats.tsv")
ANN_FILE   = Path("results/analysis/exclusive_genes_annotated.tsv")
OUT_PNG    = Path("results/analysis/coverage_ideogram.png")
OUT_PNG.parent.mkdir(parents=True, exist_ok=True)

# BAM prefix  →  (sample_id, group)
SAMPLE_MAP = {
    "TK1049": ("TK10","GroupA"), "TK1050": ("TK10","GroupA"), "TK1051": ("TK10","GroupA"),
    "TK12R2": ("TK12","GroupB"), "TK12R3": ("TK12","GroupB"),
    "TK131":  ("TK13","GroupB"), "TK132":  ("TK13","GroupB"), "TK133":  ("TK13","GroupB"),
    "TK141":  ("TK14","GroupB"), "TK142":  ("TK14","GroupB"), "TK143":  ("TK14","GroupB"),
    "TK16R1": ("TK16","GroupA"), "TK16R2": ("TK16","GroupA"), "TK16R3": ("TK16","GroupA"),
    "TK18R1": ("TK18","GroupA"), "TK18R2": ("TK18","GroupA"), "TK18R3": ("TK18","GroupA"),
    "TK91L003":("TK91","GroupA"),"TK91L005":("TK91","GroupA"),"TK91L006":("TK91","GroupA"),
    "TK92L002":("TK92","GroupA"),"TK92L003":("TK92","GroupA"),
    "TK92L005":("TK92","GroupA"),"TK92L006":("TK92","GroupA"),
    "TK93L002":("TK93","GroupA"),"TK93L003":("TK93","GroupA"),
    "TK93L005":("TK93","GroupA"),"TK93L006":("TK93","GroupA"),
}

# Chromosome order and GRCh38 lengths (bp)
CHR_ORDER = [str(i) for i in range(1, 23)] + ["X", "Y"]
CHR_LEN = {
    "1":248956422,"2":242193529,"3":198295559,"4":190214555,
    "5":181538259,"6":170805979,"7":159345973,"8":145138636,
    "9":138394717,"10":133797422,"11":135086622,"12":133275309,
    "13":114364328,"14":107043718,"15":101991189,"16":90338345,
    "17":83257441, "18":80373285, "19":58617616, "20":64444167,
    "21":46709983, "22":50818468, "X":156040895, "Y":57227415,
}

# GRCh38 centromere positions (approximate bp)
CENTROMERE = {
    "1":123400000,"2":93900000,"3":90900000,"4":50000000,
    "5":48800000,"6":59800000,"7":60100000,"8":45200000,
    "9":43000000,"10":39800000,"11":53400000,"12":35500000,
    "13":17700000,"14":17200000,"15":19000000,"16":36800000,
    "17":25100000,"18":18500000,"19":26200000,"20":28100000,
    "21":12000000,"22":15000000,"X":61000000,"Y":10400000,
}

# Colours per sample (distinct but grouped by A/B)
GROUP_PALETTE = {
    "GroupA": "#E63946",  # red family
    "GroupB": "#1D70B8",  # blue family
}

# Per-sample tints
SAMPLE_COLORS = {
    "TK10": "#FF6B6B","TK16": "#E63946","TK18": "#C1121F",
    "TK91": "#FFAA5A","TK92": "#FB8500","TK93": "#E07000",
    "TK12": "#4895EF","TK13": "#1D70B8","TK14": "#023E8A",
}

# ── Load idxstats ──────────────────────────────────────────────────────────
rows = []
with open(IDXSTATS) as fh:
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 4: continue
        bam_prefix, chrom, chrom_len, mapped = parts[0], parts[1], int(parts[2]), int(parts[3])
        if chrom not in CHR_ORDER: continue
        info = SAMPLE_MAP.get(bam_prefix)
        if info is None: continue
        sample_id, group = info
        rows.append({"bam": bam_prefix, "sample": sample_id, "group": group,
                     "chrom": chrom, "chrom_len": chrom_len, "mapped": mapped})

df = pd.DataFrame(rows)

# Total mapped per BAM (across all autosomes + X + Y)
total_mapped = df.groupby("bam")["mapped"].sum().rename("total_mapped")
df = df.join(total_mapped, on="bam")

# Normalised: reads per Mb per million total reads (RPKM-style per chrom)
df["chrom_len_Mb"] = df["chrom_len"] / 1e6
df["rpm"]  = df["mapped"] / (df["total_mapped"] / 1e6)   # reads per million total
df["rpkm"] = df["rpm"]    / df["chrom_len_Mb"]            # per Mb

# Per-sample mean (average over replicates for each sample)
sample_mean = (df.groupby(["sample","group","chrom"])["rpkm"]
                 .mean().reset_index().rename(columns={"rpkm":"rpkm_mean"}))

# ── Load exclusive gene annotation ────────────────────────────────────────
ann = pd.read_csv(ANN_FILE, sep="\t")
ann = ann[ann["chrom"].isin(CHR_ORDER)]

def gene_density(group_label, vtype_label, chrom):
    sub = ann[(ann["group"]==group_label) & (ann["vtype"]==vtype_label) & (ann["chrom"]==chrom)]
    return len(sub)

# ── Plot ───────────────────────────────────────────────────────────────────
N_CHR = len(CHR_ORDER)

# Layout: for each chromosome a column, 4 row tracks:
#   0: coverage lines (tall)
#   1: cytoband ideogram (thin)
#   2: SV exclusive bar (thin)
#   3: SNV exclusive bar (thin)

fig = plt.figure(figsize=(28, 10), facecolor="#0d1117")
fig.subplots_adjust(left=0.04, right=0.98, top=0.88, bottom=0.06,
                    wspace=0.08, hspace=0.0)

# Use gridspec: 24 columns, 4 rows
outer_gs = gridspec.GridSpec(4, N_CHR, figure=fig,
                              height_ratios=[6, 0.5, 0.8, 0.8],
                              hspace=0.05, wspace=0.08,
                              left=0.04, right=0.98,
                              top=0.88, bottom=0.06)

# Pre-compute global max rpkm for consistent y-axis (exclude extreme outliers)
global_max = np.percentile(df["rpkm"].values, 98)

# Coverage axes list
cov_axes = []

for ci, chrom in enumerate(CHR_ORDER):
    ax_cov  = fig.add_subplot(outer_gs[0, ci])
    ax_ideo = fig.add_subplot(outer_gs[1, ci])
    ax_sv   = fig.add_subplot(outer_gs[2, ci])
    ax_snv  = fig.add_subplot(outer_gs[3, ci])
    cov_axes.append(ax_cov)

    chrom_df = df[df["chrom"] == chrom]
    chrom_sm = sample_mean[sample_mean["chrom"] == chrom]

    # ── Coverage: individual replicate lines ──────────────────────────────
    for bam_id, bdf in chrom_df.groupby("bam"):
        sample = bdf["sample"].iloc[0]
        group  = bdf["group"].iloc[0]
        color  = SAMPLE_COLORS.get(sample, GROUP_PALETTE[group])
        val    = bdf["rpkm"].iloc[0]
        ax_cov.axhline(val, color=color, alpha=0.35, linewidth=0.9, zorder=2)

    # Per-sample mean (thick line)
    for _, row in chrom_sm.iterrows():
        color = SAMPLE_COLORS.get(row["sample"], GROUP_PALETTE[row["group"]])
        ax_cov.axhline(row["rpkm_mean"], color=color, alpha=0.9,
                       linewidth=2.0, zorder=3)

    # Group mean (dashed white)
    for grp, gdf in chrom_sm.groupby("group"):
        gmean = gdf["rpkm_mean"].mean()
        ax_cov.axhline(gmean, color=GROUP_PALETTE[grp], alpha=1.0,
                       linewidth=2.5, linestyle="--", zorder=4)

    ax_cov.set_xlim(0, 1)
    ax_cov.set_ylim(0, global_max * 1.05)
    ax_cov.set_xticks([])
    ax_cov.set_facecolor("#0d1117")
    for sp in ax_cov.spines.values():
        sp.set_visible(False)

    if ci == 0:
        ax_cov.set_ylabel("Reads/Mb\n(norm)", color="#aaaaaa", fontsize=6, labelpad=2)
        ax_cov.tick_params(axis="y", colors="#aaaaaa", labelsize=5, length=2, pad=1)
    else:
        ax_cov.set_yticks([])

    ax_cov.set_title(f"chr{chrom}", color="#cccccc", fontsize=6.5,
                     pad=2, fontweight="bold")

    # ── Cytoband ideogram ─────────────────────────────────────────────────
    cen = CENTROMERE.get(chrom, CHR_LEN.get(chrom, 1) // 2)
    clen = CHR_LEN.get(chrom, 1)
    # p arm
    ax_ideo.barh(0, cen/clen, left=0,   height=0.6, color="#888888", align="center")
    # q arm
    ax_ideo.barh(0, 1-cen/clen, left=cen/clen, height=0.6, color="#aaaaaa", align="center")
    # centromere marker
    ax_ideo.plot([cen/clen, cen/clen], [-0.3, 0.3], color="#ff4444", lw=1.5, zorder=5)
    ax_ideo.set_xlim(0, 1); ax_ideo.set_ylim(-0.5, 0.5)
    ax_ideo.set_xticks([]); ax_ideo.set_yticks([])
    ax_ideo.set_facecolor("#0d1117")
    for sp in ax_ideo.spines.values(): sp.set_visible(False)

    # ── SV exclusive gene density bars ────────────────────────────────────
    sv_a = gene_density("GroupA","SV",chrom)
    sv_b = gene_density("GroupB","SV",chrom)
    sv_max = max(1, sv_a + sv_b)
    ax_sv.barh(0.3, sv_a/sv_max, left=0,        height=0.35,
               color=GROUP_PALETTE["GroupA"], alpha=0.85)
    ax_sv.barh(0.3, sv_b/sv_max, left=sv_a/sv_max, height=0.35,
               color=GROUP_PALETTE["GroupB"], alpha=0.85)
    ax_sv.set_xlim(0, 1); ax_sv.set_ylim(0, 0.65)
    ax_sv.set_xticks([]); ax_sv.set_yticks([])
    ax_sv.set_facecolor("#0d1117")
    for sp in ax_sv.spines.values(): sp.set_visible(False)
    if ci == 0:
        ax_sv.set_ylabel("SV excl", color="#aaaaaa", fontsize=5, labelpad=2)
        ax_sv.text(-0.18, 0.3, f"{sv_a}A\n{sv_b}B", transform=ax_sv.transAxes,
                   color="#888888", fontsize=4, va="center")

    # ── SNV exclusive gene density bars ───────────────────────────────────
    snv_a = gene_density("GroupA","SNV",chrom)
    snv_b = gene_density("GroupB","SNV",chrom)
    snv_max = max(1, snv_a + snv_b)
    ax_snv.barh(0.3, snv_a/snv_max, left=0,         height=0.35,
                color=GROUP_PALETTE["GroupA"], alpha=0.85)
    ax_snv.barh(0.3, snv_b/snv_max, left=snv_a/snv_max, height=0.35,
                color=GROUP_PALETTE["GroupB"], alpha=0.85)
    ax_snv.set_xlim(0, 1); ax_snv.set_ylim(0, 0.65)
    ax_snv.set_xticks([]); ax_snv.set_yticks([])
    ax_snv.set_facecolor("#0d1117")
    for sp in ax_snv.spines.values(): sp.set_visible(False)
    if ci == 0:
        ax_snv.set_ylabel("SNV excl", color="#aaaaaa", fontsize=5, labelpad=2)

# ── Global title and legend ────────────────────────────────────────────────
fig.text(0.5, 0.955, "Per-Sample Read Coverage across Chromosomes  |  GroupA (TK10/16/18/91/92/93) vs GroupB (TK12/13/14)",
         ha="center", va="center", color="white", fontsize=11, fontweight="bold")
fig.text(0.5, 0.925, "Lines = normalised read depth per replicate (thin) and per sample (thick)  |  Dashed = group mean  |  Bars = exclusive gene load (A=red, B=blue)",
         ha="center", va="center", color="#aaaaaa", fontsize=8)

# Legend
legend_elements = []
for sid, col in sorted(SAMPLE_COLORS.items()):
    grp = "GroupA" if sid in {"TK10","TK16","TK18","TK91","TK92","TK93"} else "GroupB"
    legend_elements.append(Line2D([0],[0], color=col, lw=2,
                                  label=f"{sid} ({grp})"))
legend_elements.append(Line2D([0],[0], color="white", lw=2, linestyle="--",
                               label="Group mean (dashed)"))
legend_elements.append(mpatches.Patch(color="#E63946", alpha=0.85, label="SV/SNV exclusive – GroupA"))
legend_elements.append(mpatches.Patch(color="#1D70B8", alpha=0.85, label="SV/SNV exclusive – GroupB"))

fig.legend(handles=legend_elements, loc="upper right",
           bbox_to_anchor=(0.99, 0.99),
           ncol=2, fontsize=6.5, framealpha=0.15,
           facecolor="#1a1a2e", edgecolor="#444444",
           labelcolor="white")

plt.savefig(OUT_PNG, dpi=180, bbox_inches="tight", facecolor="#0d1117")
print(f"Saved: {OUT_PNG}")
plt.close()
