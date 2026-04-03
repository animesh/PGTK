#!/usr/bin/env python3
"""
mannwhitney_perchrom.py

Runs Mann-Whitney U test on per-chromosome normalised read depth (RPKM)
using existing samtools idxstats data (already collected).

No new BAM scanning needed — runs in seconds.

Output:
  results/analysis/perchrom_mannwhitney_idxstats.tsv
  results/analysis/perchrom_mannwhitney_plot.png
"""

import csv
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

# ── Config ─────────────────────────────────────────────────────────────────
IDXSTATS = Path("results/coverage/all_idxstats.tsv")
ANN_FILE = Path("results/analysis/exclusive_genes_annotated.tsv")
OUT_TSV  = Path("results/analysis/perchrom_mannwhitney_idxstats.tsv")
OUT_PNG  = Path("results/analysis/perchrom_mannwhitney_plot.png")

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
SCOL = {
    "TK10":"#FF6B6B","TK16":"#E63946","TK18":"#C1121F",
    "TK91":"#FFAA5A","TK92":"#FB8500","TK93":"#CC6600",
    "TK12":"#4895EF","TK13":"#1D70B8","TK14":"#023E8A",
}
GCOL = {"A":"#E63946","B":"#1D70B8"}

# ── Load idxstats ──────────────────────────────────────────────────────────
rows = []
with open(IDXSTATS) as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) < 4: continue
        bam, chrom, clen, mapped = p[0], p[1], int(p[2]), int(p[3])
        if chrom not in CHR_ORDER: continue
        info = SAMPLE_MAP.get(bam)
        if not info: continue
        sid, grp = info
        rows.append({"bam":bam,"sample":sid,"group":grp,
                     "chrom":chrom,"chrom_len":clen,"mapped":mapped})

df = pd.DataFrame(rows)
tot = df.groupby("bam")["mapped"].sum().rename("total")
df  = df.join(tot, on="bam")
# RPKM: reads per million total, per Mb of chrom
df["rpkm"] = (df["mapped"] / (df["total"]/1e6)) / (df["chrom_len"]/1e6)
# log2 ratio relative to per-sample genome-wide mean
bam_mean = df.groupby("bam")["rpkm"].mean().rename("bam_mean_rpkm")
df = df.join(bam_mean, on="bam")
df["log2r"] = np.log2((df["rpkm"] / df["bam_mean_rpkm"]).replace(0, np.nan))

# Per-sample mean (average replicates)
smean = (df.groupby(["sample","group","chrom"])["log2r"]
           .mean().reset_index().rename(columns={"log2r":"log2r_mean"}))

# ── Mann-Whitney U per chromosome ─────────────────────────────────────────
mw_rows = []
for chrom in CHR_ORDER:
    sub = smean[smean["chrom"]==chrom]
    A = sub[sub["group"]=="A"]["log2r_mean"].dropna().values
    B = sub[sub["group"]=="B"]["log2r_mean"].dropna().values
    if len(A) < 2 or len(B) < 2:
        mw_rows.append({"chrom":chrom,"n_A":len(A),"n_B":len(B),
                        "mean_A":np.nan,"mean_B":np.nan,"log2fc":np.nan,
                        "U":np.nan,"pval":np.nan})
        continue
    U, pval = mannwhitneyu(A, B, alternative="two-sided")
    mw_rows.append({"chrom":chrom,"n_A":len(A),"n_B":len(B),
                    "mean_A":round(float(np.mean(A)),4),
                    "mean_B":round(float(np.mean(B)),4),
                    "log2fc":round(float(np.mean(A)-np.mean(B)),4),
                    "U":round(float(U),1),"pval":pval})

mw = pd.DataFrame(mw_rows)
valid = mw["pval"].notna()
if valid.sum() > 0:
    _, padj, _, _ = multipletests(mw.loc[valid,"pval"], method="fdr_bh")
    mw.loc[valid,"padj"] = padj
mw["-log10p"]    = -np.log10(mw["pval"].clip(lower=1e-300))
mw["-log10padj"] = -np.log10(mw["padj"].clip(lower=1e-300))
mw.to_csv(OUT_TSV, sep="\t", index=False, float_format="%.4g")

# Print table
print(f"\n{'='*72}")
print("  Per-Chromosome Mann-Whitney U  (GroupA vs GroupB  |  log2 RPKM ratio)")
print(f"{'='*72}")
print(f"  {'Chr':<4}  {'n_A':>4}  {'n_B':>4}  {'mean_A':>7}  {'mean_B':>7}  "
      f"{'log2FC':>7}  {'p-val':>9}  {'FDR':>9}  sig")
print("  " + "─"*70)
for _, r in mw.iterrows():
    sig = ("***" if r.get("padj",1) < 0.001 else
           "**"  if r.get("padj",1) < 0.01  else
           "*"   if r.get("padj",1) < 0.05  else "")
    fc_arrow = "▲A" if r["log2fc"] > 0 else ("▼B" if r["log2fc"] < 0 else "  ")
    print(f"  {r['chrom']:<4}  {int(r['n_A']):>4}  {int(r['n_B']):>4}  "
          f"{r['mean_A']:>7.3f}  {r['mean_B']:>7.3f}  "
          f"{r['log2fc']:>7.3f} {fc_arrow}  {r['pval']:>9.2e}  "
          f"{r.get('padj',float('nan')):>9.2e}  {sig}")

print(f"\n  Saved → {OUT_TSV}")

# ── Load annotation ────────────────────────────────────────────────────────
ann = pd.read_csv(ANN_FILE, sep="\t")
ann = ann[ann["chrom"].isin(CHR_ORDER)]

def get_counts(c):
    s = ann[ann["chrom"]==c]
    return (len(s[(s.group=="GroupA")&(s.vtype=="SV")]),
            len(s[(s.group=="GroupB")&(s.vtype=="SV")]),
            len(s[(s.group=="GroupA")&(s.vtype=="SNV")]),
            len(s[(s.group=="GroupB")&(s.vtype=="SNV")]))

# ── Figure ─────────────────────────────────────────────────────────────────
N   = len(CHR_ORDER)
fig = plt.figure(figsize=(26, 18), facecolor="#0d1117")
gs  = gridspec.GridSpec(4, 1, figure=fig,
                        height_ratios=[5, 3, 2, 1],
                        hspace=0.08,
                        left=0.06, right=0.98,
                        top=0.93, bottom=0.07)

ax_dot  = fig.add_subplot(gs[0])   # dot plot: per-sample log2r per chrom
ax_bar  = fig.add_subplot(gs[1])   # -log10(padj) bar
ax_gene = fig.add_subplot(gs[2])   # gene count bars
ax_lbl  = fig.add_subplot(gs[3])   # chr labels

BG="#0d1117"
for ax in [ax_dot, ax_bar, ax_gene, ax_lbl]:
    ax.set_facecolor(BG)
    for sp in ax.spines.values():
        sp.set_color("#333333"); sp.set_linewidth(0.5)

x_pos = {c: i for i, c in enumerate(CHR_ORDER)}

# ── Panel 1: Dot/strip plot of log2r per sample ───────────────────────────
for _, row in smean.iterrows():
    xi = x_pos[row["chrom"]]
    col = SCOL.get(row["sample"], GCOL[row["group"]])
    jitter = np.random.uniform(-0.18, 0.18)
    ax_dot.scatter(xi + jitter, row["log2r_mean"],
                   color=col, s=28, alpha=0.8, zorder=3,
                   edgecolors="none")

# Group mean ± box
for grp, gcol in GCOL.items():
    for chrom in CHR_ORDER:
        xi = x_pos[chrom]
        vals = smean[(smean.group==grp)&(smean.chrom==chrom)]["log2r_mean"].dropna().values
        if len(vals) == 0: continue
        gm = np.mean(vals)
        gs_ = np.std(vals, ddof=1) if len(vals)>1 else 0
        ax_dot.plot([xi-0.35, xi+0.35], [gm, gm],
                    color=gcol, lw=2.5, alpha=0.9, zorder=4,
                    solid_capstyle="round")
        if gs_ > 0:
            ax_dot.fill_between([xi-0.35, xi+0.35],
                                 [gm-gs_, gm-gs_],
                                 [gm+gs_, gm+gs_],
                                 color=gcol, alpha=0.12, zorder=2)

ax_dot.axhline(0, color="#ffffff", lw=0.8, alpha=0.3, linestyle=":")
ax_dot.set_xlim(-0.7, N-0.3)
ax_dot.set_ylim(-2.0, 2.0)
ax_dot.set_xticks([])
ax_dot.set_ylabel("log₂(RPKM / genome mean)\nper sample", color="#aaaaaa", fontsize=9)
ax_dot.tick_params(axis="y", colors="#888888", labelsize=8)
ax_dot.yaxis.set_major_locator(plt.MultipleLocator(0.5))
ax_dot.spines["top"].set_visible(False)
ax_dot.spines["bottom"].set_visible(False)
ax_dot.set_title(
    "Per-Chromosome Read Depth (log₂ normalised)  ·  Mann-Whitney U  ·  "
    "GroupA (TK10/16/18/91/92/93) vs GroupB (TK12/13/14)",
    color="white", fontsize=11, fontweight="bold", pad=8)

# ── Panel 2: -log10(FDR) bar ──────────────────────────────────────────────
sig_y = -np.log10(0.05)
for _, r in mw.iterrows():
    xi = x_pos[r["chrom"]]
    h  = r["-log10padj"] if not np.isnan(r["-log10padj"]) else 0
    fc = r["log2fc"] if not np.isnan(r["log2fc"]) else 0
    col = GCOL["A"] if fc > 0 else GCOL["B"]
    ax_bar.bar(xi, h, width=0.65, color=col, alpha=0.8, zorder=2)
    # Annotate log2FC
    if h > 0.2:
        arrow = "▲" if fc > 0 else "▼"
        ax_bar.text(xi, h + 0.05, f"{arrow}{abs(fc):.2f}",
                    ha="center", va="bottom", color="#dddddd",
                    fontsize=6, zorder=5)

ax_bar.axhline(sig_y, color="#ffcc00", lw=1.5, linestyle="--",
               alpha=0.85, zorder=3, label="FDR 5% (p_adj=0.05)")
ax_bar.axhline(-np.log10(0.01), color="#ff8800", lw=1.0, linestyle=":",
               alpha=0.7, zorder=3, label="FDR 1%")

ymax = max(mw["-log10padj"].max() * 1.2, 2)
ax_bar.set_xlim(-0.7, N-0.3)
ax_bar.set_ylim(0, ymax)
ax_bar.set_xticks([])
ax_bar.set_ylabel("−log₁₀(FDR p-adj)\n▲A higher  ▼B higher", color="#aaaaaa", fontsize=8)
ax_bar.tick_params(axis="y", colors="#888888", labelsize=7)
ax_bar.spines["top"].set_visible(False)
ax_bar.spines["bottom"].set_visible(False)
ax_bar.legend(loc="upper right", fontsize=7, facecolor="#1a1a2e",
              edgecolor="#444444", labelcolor="white", framealpha=0.4)

# ── Panel 3: SV + SNV exclusive gene counts ───────────────────────────────
for chrom in CHR_ORDER:
    xi = x_pos[chrom]
    sv_a, sv_b, snv_a, snv_b = get_counts(chrom)
    total = max(sv_a+sv_b+snv_a+snv_b, 1)
    # Stacked: SV_A | SV_B | SNV_A | SNV_B
    lefts = [0]
    for n, col in [(sv_a,GCOL["A"]),(sv_b,GCOL["B"]),
                   (snv_a,"#FF9999"),(snv_b,"#99BBFF")]:
        ax_gene.bar(xi, n/total, bottom=lefts[-1], width=0.7,
                    color=col, alpha=0.85)
        lefts.append(lefts[-1] + n/total)

ax_gene.set_xlim(-0.7, N-0.3)
ax_gene.set_ylim(0, 1)
ax_gene.set_xticks([])
ax_gene.set_ylabel("Excl. genes\n(prop.)", color="#aaaaaa", fontsize=8)
ax_gene.tick_params(axis="y", colors="#888888", labelsize=7)
ax_gene.spines["top"].set_visible(False)
ax_gene.spines["bottom"].set_visible(False)

# ── Panel 4: Chr labels ────────────────────────────────────────────────────
for chrom in CHR_ORDER:
    xi = x_pos[chrom]
    ax_lbl.text(xi, 0.5, f"chr{chrom}",
                ha="center", va="center", color="#cccccc",
                fontsize=8, fontweight="bold",
                rotation=40 if len(chrom)>1 else 0)
ax_lbl.set_xlim(-0.7, N-0.3)
ax_lbl.set_ylim(0, 1)
ax_lbl.set_xticks([]); ax_lbl.set_yticks([])
for sp in ax_lbl.spines.values(): sp.set_visible(False)

# ── Legend ─────────────────────────────────────────────────────────────────
handles = []
for sid in ["TK10","TK16","TK18","TK91","TK92","TK93"]:
    handles.append(mpatches.Patch(color=SCOL[sid], label=f"{sid} (GroupA)"))
for sid in ["TK12","TK13","TK14"]:
    handles.append(mpatches.Patch(color=SCOL[sid], label=f"{sid} (GroupB)"))
handles += [
    mpatches.Patch(color=GCOL["A"], label="SV excl. GroupA"),
    mpatches.Patch(color=GCOL["B"], label="SV excl. GroupB"),
    mpatches.Patch(color="#FF9999", label="SNV excl. GroupA"),
    mpatches.Patch(color="#99BBFF", label="SNV excl. GroupB"),
]
fig.legend(handles=handles, loc="lower center",
           bbox_to_anchor=(0.5, 0.001),
           ncol=8, fontsize=7.5, framealpha=0.2,
           facecolor="#1a1a2e", edgecolor="#444444",
           labelcolor="white", handlelength=1.2)

plt.savefig(OUT_PNG, dpi=180, bbox_inches="tight", facecolor=BG)
print(f"Saved → {OUT_PNG}")
plt.close()
