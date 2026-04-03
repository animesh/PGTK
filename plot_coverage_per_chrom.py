#!/usr/bin/env python3
"""
plot_coverage_per_chrom.py

4-row × 6-column layout (chr 1-22 + X + Y = 24 panels).
Each panel:
  - Main: per-replicate thin lines + per-sample thick lines + group dashed mean
  - Bottom strip: cytoband (p/centromere/q) + SV/SNV exclusive gene density
  - Annotations for notable hotspots on select chromosomes
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from pathlib import Path

# ── Config ─────────────────────────────────────────────────────────────────
IDXSTATS = Path("results/coverage/all_idxstats.tsv")
ANN_FILE = Path("results/analysis/exclusive_genes_annotated.tsv")
OUT      = Path("results/analysis/coverage_per_chrom.png")

CHR_ORDER = [str(i) for i in range(1, 23)] + ["X", "Y"]
CHR_LEN = {
    "1":248956422,"2":242193529,"3":198295559,"4":190214555,
    "5":181538259,"6":170805979,"7":159345973,"8":145138636,
    "9":138394717,"10":133797422,"11":135086622,"12":133275309,
    "13":114364328,"14":107043718,"15":101991189,"16":90338345,
    "17":83257441,"18":80373285,"19":58617616,"20":64444167,
    "21":46709983,"22":50818468,"X":156040895,"Y":57227415,
}
CENTROMERE = {
    "1":123400000,"2":93900000,"3":90900000,"4":50000000,
    "5":48800000,"6":59800000,"7":60100000,"8":45200000,
    "9":43000000,"10":39800000,"11":53400000,"12":35500000,
    "13":17700000,"14":17200000,"15":19000000,"16":36800000,
    "17":25100000,"18":18500000,"19":26200000,"20":28100000,
    "21":12000000,"22":15000000,"X":61000000,"Y":10400000,
}

SAMPLE_MAP = {
    "TK1049":("TK10","A"),"TK1050":("TK10","A"),"TK1051":("TK10","A"),
    "TK12R2":("TK12","B"),"TK12R3":("TK12","B"),
    "TK131":("TK13","B"),"TK132":("TK13","B"),"TK133":("TK13","B"),
    "TK141":("TK14","B"),"TK142":("TK14","B"),"TK143":("TK14","B"),
    "TK16R1":("TK16","A"),"TK16R2":("TK16","A"),"TK16R3":("TK16","A"),
    "TK18R1":("TK18","A"),"TK18R2":("TK18","A"),"TK18R3":("TK18","A"),
    "TK91L003":("TK91","A"),"TK91L005":("TK91","A"),"TK91L006":("TK91","A"),
    "TK92L002":("TK92","A"),"TK92L003":("TK92","A"),
    "TK92L005":("TK92","A"),"TK92L006":("TK92","A"),
    "TK93L002":("TK93","A"),"TK93L003":("TK93","A"),
    "TK93L005":("TK93","A"),"TK93L006":("TK93","A"),
}

# per-sample colours (warm=GroupA, cool=GroupB)
SCOL = {
    "TK10":"#FF6B6B","TK16":"#E63946","TK18":"#C1121F",
    "TK91":"#FFAA5A","TK92":"#FB8500","TK93":"#CC6600",
    "TK12":"#4895EF","TK13":"#1D70B8","TK14":"#023E8A",
}
GCOL = {"A":"#E63946","B":"#1D70B8"}

# Notable genes to annotate on specific chromosomes
HOTSPOT_LABELS = {
    "22": [("IgL locus\n(22q11)", 0.07, "GroupA")],
    "17": [("TP53\n(17p13)", 0.05, "GroupB"), ("BRCA1\n(17q21)", 0.55, "GroupB")],
    "14": [("IgH locus\n(14q32)", 0.88, "GroupB")],
    "11": [("ATM\n(11q22)", 0.72, "GroupB")],
    "15": [("15q11\nimprinting", 0.08, "GroupA")],
    "19": [("19p13\ngene-dense", 0.06, "GroupA")],
}

# ── Load data ──────────────────────────────────────────────────────────────
rows = []
with open(IDXSTATS) as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) < 4: continue
        bam, chrom, clen, mapped = p[0], p[1], int(p[2]), int(p[3])
        if chrom not in CHR_ORDER: continue
        info = SAMPLE_MAP.get(bam)
        if info is None: continue
        sid, grp = info
        rows.append({"bam":bam,"sample":sid,"group":grp,
                     "chrom":chrom,"chrom_len":clen,"mapped":mapped})

df = pd.DataFrame(rows)
tot = df.groupby("bam")["mapped"].sum().rename("total")
df  = df.join(tot, on="bam")
df["rpkm"] = (df["mapped"] / (df["total"]/1e6)) / (df["chrom_len"]/1e6)

smean = (df.groupby(["sample","group","chrom"])["rpkm"]
           .mean().reset_index().rename(columns={"rpkm":"mean_rpkm"}))

ann = pd.read_csv(ANN_FILE, sep="\t")
ann = ann[ann["chrom"].isin(CHR_ORDER)]

def get_counts(chrom):
    sub = ann[ann["chrom"]==chrom]
    return {
        "sv_a":  len(sub[(sub.group=="GroupA")&(sub.vtype=="SV")]),
        "sv_b":  len(sub[(sub.group=="GroupB")&(sub.vtype=="SV")]),
        "snv_a": len(sub[(sub.group=="GroupA")&(sub.vtype=="SNV")]),
        "snv_b": len(sub[(sub.group=="GroupB")&(sub.vtype=="SNV")]),
    }

global_max = np.percentile(df["rpkm"], 98)

# ── Figure ─────────────────────────────────────────────────────────────────
NCOLS, NROWS = 6, 4
fig = plt.figure(figsize=(24, 20), facecolor="#0d1117")

# Title
fig.text(0.5, 0.985,
         "Read Coverage per Chromosome — GroupA (TK10·16·18·91·92·93) vs GroupB (TK12·13·14)",
         ha="center", va="top", color="white", fontsize=14, fontweight="bold")
fig.text(0.5, 0.972,
         "Thin lines = individual replicates  ·  Thick lines = per-sample mean  ·  "
         "Dashed = group mean  ·  Colour bars = exclusive gene load (A=red · B=blue · top=SV · bottom=SNV)",
         ha="center", va="top", color="#aaaaaa", fontsize=8.5)

outer = gridspec.GridSpec(NROWS, NCOLS, figure=fig,
                          hspace=0.42, wspace=0.25,
                          left=0.05, right=0.99,
                          top=0.96, bottom=0.06)

for idx, chrom in enumerate(CHR_ORDER):
    row_i = idx // NCOLS
    col_i = idx  % NCOLS

    # Each cell: 3 rows (coverage, ideogram, gene bars)
    inner = gridspec.GridSpecFromSubplotSpec(
        3, 1,
        subplot_spec=outer[row_i, col_i],
        height_ratios=[5, 0.6, 1.0],
        hspace=0.0
    )

    ax_c = fig.add_subplot(inner[0])   # coverage
    ax_i = fig.add_subplot(inner[1])   # ideogram
    ax_g = fig.add_subplot(inner[2])   # gene bars

    clen = CHR_LEN.get(chrom, 1)
    cen  = CENTROMERE.get(chrom, clen//2)

    # ── Coverage ──────────────────────────────────────────────────────────
    chrom_df = df[df["chrom"]==chrom]
    chrom_sm = smean[smean["chrom"]==chrom]

    # thin replicate lines
    for _, r in chrom_df.iterrows():
        col = SCOL.get(r["sample"], GCOL[r["group"]])
        ax_c.axhline(r["rpkm"], color=col, alpha=0.28, linewidth=0.8, zorder=2)

    # thick per-sample mean
    for _, r in chrom_sm.iterrows():
        col = SCOL.get(r["sample"], GCOL[r["group"]])
        ax_c.axhline(r["mean_rpkm"], color=col, alpha=0.85, linewidth=2.2, zorder=3)

    # dashed group mean
    for grp, gdf in chrom_sm.groupby("group"):
        gm = gdf["mean_rpkm"].mean()
        ax_c.axhline(gm, color=GCOL[grp], alpha=1.0,
                     linewidth=2.5, linestyle="--", zorder=4)

    ax_c.set_xlim(0,1)
    ax_c.set_ylim(0, global_max*1.08)
    ax_c.set_xticks([])
    ax_c.set_facecolor("#0d1117")
    ax_c.tick_params(axis="y", colors="#888888", labelsize=6, length=2, pad=1)
    ax_c.yaxis.set_major_locator(plt.MaxNLocator(3, integer=True))
    for sp in ax_c.spines.values():
        sp.set_color("#333333")
        sp.set_linewidth(0.5)
    ax_c.spines["top"].set_visible(False)
    ax_c.spines["right"].set_visible(False)
    ax_c.spines["bottom"].set_visible(False)

    # Chromosome title
    ax_c.set_title(f"chr{chrom}", color="#dddddd", fontsize=9,
                   pad=3, fontweight="bold")

    # Hotspot annotations
    if chrom in HOTSPOT_LABELS:
        for label, xpos, grp in HOTSPOT_LABELS[chrom]:
            col = GCOL.get(grp[5:7] if "Group" in grp else grp, "#ffffff")
            ax_c.axvline(xpos, color=col, alpha=0.4, lw=1, linestyle=":", zorder=1)
            ax_c.text(xpos+0.02, global_max*0.92, label,
                      color=col, fontsize=5.5, alpha=0.9,
                      va="top", ha="left", linespacing=1.2)

    # ── Ideogram ──────────────────────────────────────────────────────────
    p_frac   = cen / clen
    ax_i.barh(0, p_frac,      left=0,      height=0.55, color="#666666", align="center")
    ax_i.barh(0, 1-p_frac,    left=p_frac, height=0.55, color="#999999", align="center")
    ax_i.plot([p_frac,p_frac],[-0.35,0.35], color="#ff5555", lw=2, zorder=5)
    ax_i.set_xlim(0,1); ax_i.set_ylim(-0.5,0.5)
    ax_i.set_xticks([]); ax_i.set_yticks([])
    ax_i.set_facecolor("#0d1117")
    for sp in ax_i.spines.values(): sp.set_visible(False)
    ax_i.text(0.02, 0, "p", color="#aaaaaa", fontsize=5.5, va="center")
    ax_i.text(p_frac+0.03, 0, "q", color="#cccccc", fontsize=5.5, va="center")

    # ── Gene bars (SV top half, SNV bottom half) ──────────────────────────
    cnt = get_counts(chrom)
    sv_tot  = max(cnt["sv_a"]  + cnt["sv_b"],  1)
    snv_tot = max(cnt["snv_a"] + cnt["snv_b"], 1)

    # SV row (y=0.65)
    ax_g.barh(0.70, cnt["sv_a"]/sv_tot,  left=0,                   height=0.3,
              color=GCOL["A"], alpha=0.85)
    ax_g.barh(0.70, cnt["sv_b"]/sv_tot,  left=cnt["sv_a"]/sv_tot,  height=0.3,
              color=GCOL["B"], alpha=0.85)
    # SNV row (y=0.25)
    ax_g.barh(0.25, cnt["snv_a"]/snv_tot, left=0,                    height=0.3,
              color=GCOL["A"], alpha=0.75)
    ax_g.barh(0.25, cnt["snv_b"]/snv_tot, left=cnt["snv_a"]/snv_tot, height=0.3,
              color=GCOL["B"], alpha=0.75)

    # Labels
    ax_g.text(-0.02, 0.70, "SV",  ha="right", va="center",
              color="#aaaaaa", fontsize=5, transform=ax_g.transAxes if False else
              ax_g.transData.__class__(ax_g) if False else ax_g.transData)
    ax_g.text(1.01,  0.70, f"{cnt['sv_a']}A/{cnt['sv_b']}B",
              ha="left", va="center", color="#888888", fontsize=4.5,
              transform=ax_g.transAxes)
    ax_g.text(1.01,  0.25, f"{cnt['snv_a']}A/{cnt['snv_b']}B",
              ha="left", va="center", color="#888888", fontsize=4.5,
              transform=ax_g.transAxes)

    ax_g.set_xlim(0,1); ax_g.set_ylim(0,1)
    ax_g.set_xticks([]); ax_g.set_yticks([])
    ax_g.set_facecolor("#0d1117")
    for sp in ax_g.spines.values(): sp.set_visible(False)
    ax_g.text(-0.01, 0.70, "SV",  ha="right", va="center",
              color="#777777", fontsize=5, transform=ax_g.transAxes)
    ax_g.text(-0.01, 0.25, "SNV", ha="right", va="center",
              color="#777777", fontsize=5, transform=ax_g.transAxes)

# ── Legend ─────────────────────────────────────────────────────────────────
handles = []
for sid in ["TK10","TK16","TK18","TK91","TK92","TK93"]:
    handles.append(Line2D([0],[0],color=SCOL[sid],lw=2.5,label=f"{sid} (GroupA)"))
for sid in ["TK12","TK13","TK14"]:
    handles.append(Line2D([0],[0],color=SCOL[sid],lw=2.5,label=f"{sid} (GroupB)"))
handles.append(Line2D([0],[0],color="#E63946",lw=2,linestyle="--",label="GroupA mean"))
handles.append(Line2D([0],[0],color="#1D70B8",lw=2,linestyle="--",label="GroupB mean"))
handles.append(mpatches.Patch(color="#E63946",alpha=0.85,label="Excl. genes – GroupA"))
handles.append(mpatches.Patch(color="#1D70B8",alpha=0.85,label="Excl. genes – GroupB"))

fig.legend(handles=handles, loc="lower center",
           bbox_to_anchor=(0.5, 0.002),
           ncol=8, fontsize=7.5, framealpha=0.2,
           facecolor="#1a1a2e", edgecolor="#444444",
           labelcolor="white", handlelength=1.5)

plt.savefig(OUT, dpi=180, bbox_inches="tight", facecolor="#0d1117")
print(f"Saved → {OUT}")
plt.close()
