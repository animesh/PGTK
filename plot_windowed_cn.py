#!/usr/bin/env python3
"""
plot_windowed_cn.py

Reads samtools bedcov output (one file per BAM) and produces:
  1. 1 Mb windowed copy-number profile — one line per sample coloured by group
  2. Per-chromosome Mann-Whitney U test (GroupA vs GroupB RPKM)
  3. Final combined figure:
       Top panel  : CN profile across genome (all samples as lines)
       Middle panel: -log10(p) per chromosome + effect size (fold change)
       Bottom strip: cytoband + SV/SNV exclusive gene bar

Run AFTER run_bedcov_parallel.sh completes (checks for DONE sentinel file).
"""

import sys, csv
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

# ── Paths ──────────────────────────────────────────────────────────────────
BASE      = Path(__file__).parent
BEDCOV_DIR = BASE / "results/coverage/bedcov"
DONE_FILE  = BEDCOV_DIR / "DONE"
ANN_FILE   = BASE / "results/analysis/exclusive_genes_annotated.tsv"
OUT_CN     = BASE / "results/analysis/windowed_cn_profile.png"
OUT_STATS  = BASE / "results/analysis/perchrom_mannwhitney.tsv"

# ── Check bedcov is complete ───────────────────────────────────────────────
if not DONE_FILE.exists():
    done = list(BEDCOV_DIR.glob("*.bedcov.tsv"))
    print(f"bedcov not fully complete yet ({len(done)}/28 BAMs done).")
    print(f"Missing: {DONE_FILE}")
    print("Re-run after run_bedcov_parallel.sh finishes.")
    sys.exit(0)

# ── Config ─────────────────────────────────────────────────────────────────
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

# ── Load bedcov files ──────────────────────────────────────────────────────
print("Loading bedcov files...")
records = []
for bam_prefix, (sample, group) in SAMPLE_MAP.items():
    f = BEDCOV_DIR / f"{bam_prefix}.bedcov.tsv"
    if not f.exists():
        print(f"  WARNING: missing {f}")
        continue
    with open(f) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4: continue
            chrom, start, end, bases = p[0], int(p[1]), int(p[2]), int(p[3])
            if chrom not in CHR_ORDER: continue
            win_len = end - start
            depth   = bases / win_len if win_len > 0 else 0
            records.append({
                "bam": bam_prefix, "sample": sample, "group": group,
                "chrom": chrom, "start": start, "end": end,
                "bases": bases, "depth": depth,
            })

df = pd.DataFrame(records)
print(f"  {len(df):,} windows loaded across {df['bam'].nunique()} BAMs")

# Normalise depth: divide each window by sample median depth (genome-wide)
# to get relative copy number (median-normalised ~ ploidy 2 = diploid)
sample_medians = (df.groupby("bam")["depth"]
                    .median()
                    .rename("median_depth"))
df = df.join(sample_medians, on="bam")
df["rel_depth"] = df["depth"] / df["median_depth"].replace(0, np.nan)
# log2 ratio (CN-like): log2(rel_depth), clipped to [-3, 3]
df["log2r"] = np.log2(df["rel_depth"].replace(0, np.nan)).clip(-3, 3)

# Per-sample mean log2r per window (average replicates → one value per sample per window)
win_sample = (df.groupby(["sample","group","chrom","start","end"])["log2r"]
                .mean().reset_index().rename(columns={"log2r":"log2r_mean"}))

# ── Genome-wide x-axis: assign cumulative position ────────────────────────
cum_offset = {}
offset = 0
for c in CHR_ORDER:
    cum_offset[c] = offset
    offset += CHR_LEN.get(c, 0) + 5_000_000  # 5Mb gap between chromosomes

total_genome = offset
df["genome_pos"] = df["start"] + df["chrom"].map(cum_offset)
win_sample["genome_pos"] = win_sample["start"] + win_sample["chrom"].map(cum_offset)

# ── Mann-Whitney U per chromosome ─────────────────────────────────────────
print("Running Mann-Whitney tests per chromosome...")
mw_rows = []
for chrom in CHR_ORDER:
    sub = win_sample[win_sample["chrom"] == chrom]
    vals_a = sub[sub["group"]=="A"]["log2r_mean"].dropna().values
    vals_b = sub[sub["group"]=="B"]["log2r_mean"].dropna().values
    if len(vals_a) < 3 or len(vals_b) < 3:
        mw_rows.append({"chrom": chrom, "n_windows": len(sub),
                        "mean_A": np.nan, "mean_B": np.nan,
                        "log2fc": np.nan, "U": np.nan, "pval": np.nan})
        continue
    U, pval = mannwhitneyu(vals_a, vals_b, alternative="two-sided")
    mw_rows.append({
        "chrom": chrom,
        "n_windows": len(sub) // win_sample["sample"].nunique(),
        "mean_A": round(float(np.mean(vals_a)), 4),
        "mean_B": round(float(np.mean(vals_b)), 4),
        "log2fc": round(float(np.mean(vals_a) - np.mean(vals_b)), 4),
        "U": round(float(U), 1),
        "pval": pval,
    })

mw_df = pd.DataFrame(mw_rows)
# BH correction
valid = mw_df["pval"].notna()
if valid.sum() > 0:
    _, padj, _, _ = multipletests(mw_df.loc[valid, "pval"], method="fdr_bh")
    mw_df.loc[valid, "padj"] = padj
else:
    mw_df["padj"] = np.nan

mw_df["-log10p"]    = -np.log10(mw_df["pval"].clip(lower=1e-300))
mw_df["-log10padj"] = -np.log10(mw_df["padj"].clip(lower=1e-300))
mw_df.to_csv(OUT_STATS, sep="\t", index=False)
print(f"  Saved stats → {OUT_STATS}")

# Print summary
print("\n  Per-chromosome Mann-Whitney (GroupA vs GroupB log2 CN ratio):")
print(f"  {'Chr':<4}  {'Mean_A':>7}  {'Mean_B':>7}  {'log2FC':>7}  {'p-val':>10}  {'padj':>10}  {'sig'}")
print("  " + "-"*65)
for _, r in mw_df.iterrows():
    sig = "***" if r.get("padj",1) < 0.001 else ("**" if r.get("padj",1) < 0.01 else
          ("*" if r.get("padj",1) < 0.05 else ""))
    print(f"  {r['chrom']:<4}  {r['mean_A']:>7.3f}  {r['mean_B']:>7.3f}  "
          f"{r['log2fc']:>7.3f}  {r['pval']:>10.2e}  "
          f"{r.get('padj',float('nan')):>10.2e}  {sig}")

# ── Load exclusive gene annotation ────────────────────────────────────────
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

# ── Figure ─────────────────────────────────────────────────────────────────
print("\nBuilding figure...")
fig = plt.figure(figsize=(32, 18), facecolor="#0d1117")
gs = gridspec.GridSpec(4, 1, figure=fig,
                       height_ratios=[7, 2.5, 1.2, 0.8],
                       hspace=0.08,
                       left=0.055, right=0.98,
                       top=0.93, bottom=0.06)

ax_cn    = fig.add_subplot(gs[0])   # CN profile
ax_stat  = fig.add_subplot(gs[1])   # stats
ax_ideo  = fig.add_subplot(gs[2])   # ideogram + gene bars
ax_label = fig.add_subplot(gs[3])   # chr labels

BG = "#0d1117"
for ax in [ax_cn, ax_stat, ax_ideo, ax_label]:
    ax.set_facecolor(BG)
    for sp in ax.spines.values():
        sp.set_color("#333333")
        sp.set_linewidth(0.5)

# ── Panel 1: CN profile ────────────────────────────────────────────────────
# Thin lines per replicate, thick per sample, dashed group mean per window

# Group mean log2r per window for ribbon
for grp, color in GCOL.items():
    sub = win_sample[win_sample["group"]==grp].copy()
    grp_mean = sub.groupby(["chrom","start","end","genome_pos"])["log2r_mean"].mean().reset_index()
    grp_mean = grp_mean.sort_values("genome_pos")
    ax_cn.fill_between(grp_mean["genome_pos"],
                       grp_mean["log2r_mean"] - 0.05,
                       grp_mean["log2r_mean"] + 0.05,
                       color=color, alpha=0.12, zorder=1)

# Individual replicate lines (thin, low alpha)
for bam, bdf in df.groupby("bam"):
    sample = SAMPLE_MAP[bam][0]
    group  = SAMPLE_MAP[bam][1]
    col    = SCOL.get(sample, GCOL[group])
    bdf_s  = bdf[bdf["chrom"].isin(CHR_ORDER)].sort_values("genome_pos")
    # Plot per chromosome segment (avoid lines crossing chr gaps)
    for chrom, cdf in bdf_s.groupby("chrom"):
        ax_cn.plot(cdf["genome_pos"], cdf["log2r"],
                   color=col, alpha=0.18, linewidth=0.5, zorder=2, rasterized=True)

# Per-sample mean lines (thicker)
for sample, sdf in win_sample.groupby("sample"):
    group = sdf["group"].iloc[0]
    col   = SCOL.get(sample, GCOL[group])
    for chrom, cdf in sdf.groupby("chrom"):
        cdf_s = cdf.sort_values("genome_pos")
        ax_cn.plot(cdf_s["genome_pos"], cdf_s["log2r_mean"],
                   color=col, alpha=0.75, linewidth=1.4, zorder=3, rasterized=True)

# Group mean line (dashed, bold)
for grp, color in GCOL.items():
    sub = win_sample[win_sample["group"]==grp]
    gm  = sub.groupby(["chrom","genome_pos"])["log2r_mean"].mean().reset_index()
    for chrom, cdf in gm.groupby("chrom"):
        cdf_s = cdf.sort_values("genome_pos")
        ax_cn.plot(cdf_s["genome_pos"], cdf_s["log2r_mean"],
                   color=color, alpha=1.0, linewidth=2.5, linestyle="--",
                   zorder=4, rasterized=True)

# Diploid reference line
ax_cn.axhline(0, color="#ffffff", linewidth=0.8, alpha=0.3, linestyle=":", zorder=5)

# Chromosome boundary shading
for ci, chrom in enumerate(CHR_ORDER):
    x0 = cum_offset[chrom]
    x1 = x0 + CHR_LEN.get(chrom, 0)
    if ci % 2 == 0:
        ax_cn.axvspan(x0, x1, color="#ffffff", alpha=0.02, zorder=0)
    # Centromere marker
    cen_x = x0 + CENTROMERE.get(chrom, CHR_LEN.get(chrom, 0)//2)
    ax_cn.axvline(cen_x, color="#ff4444", alpha=0.2, linewidth=0.6, zorder=1)

ax_cn.set_xlim(0, total_genome)
ax_cn.set_ylim(-2.5, 2.5)
ax_cn.set_ylabel("log₂ Relative Depth\n(median-normalised CN ratio)", color="#aaaaaa", fontsize=9)
ax_cn.set_xticks([])
ax_cn.tick_params(axis="y", colors="#888888", labelsize=8)
ax_cn.yaxis.set_major_locator(plt.MultipleLocator(1))
ax_cn.spines["bottom"].set_visible(False)
ax_cn.spines["top"].set_visible(False)

# ── Panel 2: Stats ─────────────────────────────────────────────────────────
# Bar per chromosome: -log10(padj), coloured by direction (A>B or B>A)
mw_df["chr_mid"] = mw_df["chrom"].apply(
    lambda c: cum_offset.get(c, 0) + CHR_LEN.get(c, 0) / 2)

sig_line = -np.log10(0.05)

for _, r in mw_df.iterrows():
    if np.isnan(r["-log10padj"]): continue
    x   = r["chr_mid"]
    h   = r["-log10padj"]
    fc  = r["log2fc"]
    col = GCOL["A"] if fc > 0 else GCOL["B"]
    ax_stat.bar(x, h,
                width=CHR_LEN.get(r["chrom"], 0) * 0.6,
                color=col, alpha=0.75, zorder=2)
    # log2FC text inside bar
    if h > 0.3:
        ax_stat.text(x, h + 0.05, f"{fc:+.2f}", ha="center", va="bottom",
                     color="#cccccc", fontsize=6.5, zorder=5)

ax_stat.axhline(sig_line, color="#ffcc00", linewidth=1.2, linestyle="--",
                alpha=0.8, zorder=3, label="FDR 5%")
ax_stat.axhline(-np.log10(0.01), color="#ff8800", linewidth=0.8, linestyle=":",
                alpha=0.6, zorder=3, label="FDR 1%")

ax_stat.set_xlim(0, total_genome)
ax_stat.set_ylim(0, max(mw_df["-log10padj"].max() * 1.15, 2))
ax_stat.set_xticks([])
ax_stat.set_ylabel("−log₁₀(p_adj)\nGroupA(+) vs GroupB(−)", color="#aaaaaa", fontsize=8)
ax_stat.tick_params(axis="y", colors="#888888", labelsize=7)
ax_stat.spines["bottom"].set_visible(False)
ax_stat.spines["top"].set_visible(False)
ax_stat.legend(loc="upper right", fontsize=7, facecolor="#1a1a2e",
               edgecolor="#444444", labelcolor="white", framealpha=0.4)

# ── Panel 3: Ideogram + gene bars ─────────────────────────────────────────
for chrom in CHR_ORDER:
    x0   = cum_offset[chrom]
    clen = CHR_LEN.get(chrom, 1)
    cen  = CENTROMERE.get(chrom, clen//2)
    x1   = x0 + clen
    xc   = x0 + cen

    # p arm
    ax_ideo.barh(0.7, cen,      left=x0, height=0.35, color="#666666", align="center")
    # q arm
    ax_ideo.barh(0.7, clen-cen, left=xc, height=0.35, color="#999999", align="center")
    # centromere
    ax_ideo.plot([xc, xc], [0.5, 0.9], color="#ff5555", lw=1.5, zorder=5)

    # SV gene bar
    cnt = get_counts(chrom)
    sv_tot  = max(cnt["sv_a"] + cnt["sv_b"], 1)
    snv_tot = max(cnt["snv_a"] + cnt["snv_b"], 1)
    bar_w   = clen * 0.9

    ax_ideo.barh(0.25, bar_w * cnt["sv_a"]/sv_tot,  left=x0, height=0.18,
                 color=GCOL["A"], alpha=0.85)
    ax_ideo.barh(0.25, bar_w * cnt["sv_b"]/sv_tot,
                 left=x0 + bar_w*cnt["sv_a"]/sv_tot, height=0.18,
                 color=GCOL["B"], alpha=0.85)

ax_ideo.set_xlim(0, total_genome)
ax_ideo.set_ylim(0, 1.1)
ax_ideo.set_xticks([])
ax_ideo.set_yticks([0.25, 0.7])
ax_ideo.set_yticklabels(["SV excl.", "Ideogram"], color="#888888", fontsize=7)
ax_ideo.spines["bottom"].set_visible(False)
ax_ideo.spines["top"].set_visible(False)

# ── Panel 4: Chromosome labels ────────────────────────────────────────────
for chrom in CHR_ORDER:
    x_mid = cum_offset[chrom] + CHR_LEN.get(chrom, 0) / 2
    ax_label.text(x_mid, 0.5, f"chr{chrom}",
                  ha="center", va="center", color="#cccccc",
                  fontsize=7.5, fontweight="bold",
                  rotation=45 if len(chrom) > 1 else 0)
ax_label.set_xlim(0, total_genome)
ax_label.set_ylim(0, 1)
ax_label.set_xticks([])
ax_label.set_yticks([])
for sp in ax_label.spines.values(): sp.set_visible(False)

# ── Title & Legend ─────────────────────────────────────────────────────────
fig.text(0.5, 0.965,
         "Windowed Copy-Number Profile (1 Mb bins)  |  GroupA (TK10·16·18·91·92·93) vs GroupB (TK12·13·14)",
         ha="center", color="white", fontsize=13, fontweight="bold")
fig.text(0.5, 0.948,
         "log₂ relative depth (median-normalised)  ·  Thin = replicate  ·  Thick = sample mean  ·  "
         "Dashed = group mean  ·  Bar chart = Mann-Whitney −log₁₀(FDR-adjusted p) per chromosome "
         "(red=GroupA higher · blue=GroupB higher)  ·  Number = log₂FC",
         ha="center", color="#aaaaaa", fontsize=8)

legend_handles = []
for sid in ["TK10","TK16","TK18","TK91","TK92","TK93"]:
    legend_handles.append(Line2D([0],[0],color=SCOL[sid],lw=2,label=f"{sid} (GroupA)"))
for sid in ["TK12","TK13","TK14"]:
    legend_handles.append(Line2D([0],[0],color=SCOL[sid],lw=2,label=f"{sid} (GroupB)"))
legend_handles += [
    Line2D([0],[0],color=GCOL["A"],lw=2.5,linestyle="--",label="GroupA mean"),
    Line2D([0],[0],color=GCOL["B"],lw=2.5,linestyle="--",label="GroupB mean"),
    Line2D([0],[0],color="#ffffff",lw=1,linestyle=":",label="Diploid (log2=0)"),
    mpatches.Patch(color=GCOL["A"],alpha=0.8,label="Excl. SV genes – GroupA"),
    mpatches.Patch(color=GCOL["B"],alpha=0.8,label="Excl. SV genes – GroupB"),
]
fig.legend(handles=legend_handles, loc="lower center",
           bbox_to_anchor=(0.5, 0.001),
           ncol=8, fontsize=7.5, framealpha=0.2,
           facecolor="#1a1a2e", edgecolor="#444444",
           labelcolor="white", handlelength=1.5)

plt.savefig(OUT_CN, dpi=180, bbox_inches="tight", facecolor=BG)
print(f"\nSaved → {OUT_CN}")
plt.close()
print("All done.")
