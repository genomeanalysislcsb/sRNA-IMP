#!/usr/bin/env python3
"""fig45_hc_figures.py — Expobiome HC demonstration figures (sRNA-Int, BHC01-05) from the tables of fig45_hc_tables.py.
Usage: fig45_hc_figures.py TABLEDIR OUTPREFIX [--type sRNA-Int]
Writes OUTPREFIX.fig4.{png,pdf} (read fate, placement, sRNA vs metaG / metaT per MAG) and
OUTPREFIX.fig5.{png,pdf} (genus recurrence, tRNA fragment classes, new layers)."""
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

ap = argparse.ArgumentParser()
ap.add_argument("tables"); ap.add_argument("out"); ap.add_argument("--type", default="sRNA-Int")
a = ap.parse_args()
T = a.type
rd = lambda n: pd.read_csv(f"{a.tables}/hc_{n}.tsv", sep="\t")

# reference palette (validated; light surface)
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK, INK2, MUTED, GRID, AXIS, SURF = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
NEUTRAL = "#c3c2b7"
DIV = LinearSegmentedColormap.from_list("div", ["#184f95", "#3987e5", "#9ec5f4", "#f0efec", "#f4a3a2", "#e34948", "#a32a2a"])
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.edgecolor": AXIS, "axes.labelcolor": INK2,
                     "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK, "axes.titlesize": 9,
                     "axes.titleweight": "bold", "axes.titlelocation": "left", "figure.facecolor": "white",
                     "axes.facecolor": "white", "axes.spines.top": False, "axes.spines.right": False,
                     "legend.frameon": False, "legend.fontsize": 7})


def tag(ax, t):
    ax.text(-0.13, 1.04, t, transform=ax.transAxes, fontsize=11, fontweight="bold", ha="right", va="bottom")


def short(s):
    return s.replace("BHC0", "P")


fate, place, mags, gen = rd("fate"), rd("placement"), rd("mags"), rd("genus")
trf, lay, rho = rd("trna_fragments"), rd("layers"), rd("spearman")
S = list(fate["sample"])

# ---------------- Fig 4 ----------------
fig = plt.figure(figsize=(7.2, 6.4))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.15], hspace=0.6, wspace=0.38, left=0.1, right=0.98, top=0.9, bottom=0.08)

# A read fate (share of trimmed reads)
ax = fig.add_subplot(gs[0, 0]); tag(ax, "a")
cats = ["tRNA", "rRNA", "CDS (sense)", "antisense", "other annotated", "intergenic", "host", "unmapped"]
cols = [CAT[0], CAT[1], CAT[2], CAT[3], CAT[4], CAT[6], CAT[7], NEUTRAL]
sh = fate[cats].div(fate[cats].sum(axis=1), axis=0) * 100
y = np.arange(len(S))[::-1]
left = np.zeros(len(S))
for c, col in zip(cats, cols):
    ax.barh(y, sh[c], left=left, color=col, height=0.7, edgecolor="white", linewidth=1)
    for yi, l, v in zip(y, left, sh[c]):
        if v >= 7:
            ax.text(l + v / 2, yi, f"{v:.0f}", ha="center", va="center", fontsize=6.5,
                    color="#0b0b0b" if col in (CAT[3], NEUTRAL) else "white")
    left += sh[c].values
ax.set_yticks(y, [short(s) for s in S]); ax.set_xlim(0, 100); ax.set_xlabel("% of trimmed reads")
ax.set_title("Read fate per participant")
ax.tick_params(axis="y", length=0)
h = [plt.Rectangle((0, 0), 1, 1, color=c) for c in cols]
ax.legend(h, cats, ncol=4, loc="upper left", bbox_to_anchor=(-0.12, -0.2), fontsize=6.5, handlelength=1, columnspacing=0.8)

# B placement in MAGs
ax = fig.add_subplot(gs[0, 1]); tag(ax, "b")
om = [T, "metaT", "metaG"]
w = 0.26
for i, o in enumerate(om):
    v = place[place.omics == o].set_index("sample").loc[S, "in_MAGs"] * 100
    ax.bar(np.arange(len(S)) + (i - 1) * w, v, width=w - 0.03, color=CAT[i], label=o)
ax.set_ylim(0, 100)
ax.set_ylabel("% of mapped reads in MAGs"); ax.set_title("Reads on binned contigs")
ax.yaxis.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True); ax.tick_params(axis="x", length=0)
ax.legend(ncol=3, loc="upper left", bbox_to_anchor=(0, 1.0))
nm = place.drop_duplicates("sample").set_index("sample").loc[S, "n_MAGs"]
ax.set_xticks(np.arange(len(S)), [f"{short(s)}\n{n} MAGs" for s, n in zip(S, nm)], fontsize=7)

# C, D sRNA vs metaG / metaT per MAG (pooled; per-participant rho)
for k, (o, t) in enumerate([("metaG", "c"), ("metaT", "d")]):
    ax = fig.add_subplot(gs[1, k]); tag(ax, t)
    x, yv = mags[f"{o}_CPM"] + 1, mags[f"{T}_CPM"] + 1
    ax.scatter(x, yv, s=10, color=CAT[0], alpha=0.55, edgecolor="white", linewidth=0.3)
    ax.set_xscale("log"); ax.set_yscale("log")
    lo, hi = 30, max(x.max(), yv.max()) * 1.5
    ax.plot([lo, hi], [lo, hi], color=MUTED, lw=0.8, ls="--")
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_xlabel(f"{o} CPM per MAG"); ax.set_ylabel(f"{T} CPM per MAG")
    r = rho[(rho.vs == o) & (rho.level == "MAG (no rRNA contigs)")].set_index("sample").loc[S, "spearman"]
    ax.set_title(f"{T} vs {o}, per MAG")
    txt = "Spearman ρ per participant\n" + "  ".join(f"{short(s)} {v:.2f}" for s, v in r.items())
    ax.text(0.02, 0.98, txt, transform=ax.transAxes, va="top", fontsize=6.3, color=INK2, linespacing=1.4)
    ax.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
fig.suptitle("Expobiome healthy participants (n = 5), intracellular small RNA", y=0.985, x=0.01, ha="left",
             fontsize=10, fontweight="bold")
for e in ("png", "pdf"):
    fig.savefig(f"{a.out}.fig4.{e}", dpi=300)
plt.close(fig)

# ---------------- Fig 5 ----------------
fig = plt.figure(figsize=(7.2, 7.6))
gs = fig.add_gridspec(2, 3, height_ratios=[1.35, 1], width_ratios=[1, 1, 0.9], hspace=0.3, wspace=0.12,
                      left=0.27, right=0.97, top=0.92, bottom=0.08)
# E/F genus heatmaps (genera in >= 4 of 5 participants)
gen["log2_sRNA_vs_metaG"] = np.log2((gen[f"{T}_CPM"] + 1) / (gen["metaG_CPM"] + 1))
n = gen.groupby("genus")["sample"].nunique()
rec = n[n >= 4].index
g = gen[gen.genus.isin(rec)]
phy = mags.dropna(subset=["genus"]).drop_duplicates("genus").set_index("genus")["classification"].str.extract(r"p__([^;]*)")[0]
pG = g.pivot(index="genus", columns="sample", values="log2_sRNA_vs_metaG").reindex(columns=S)
pT = g.pivot(index="genus", columns="sample", values="log2_sRNA_vs_metaT").reindex(columns=S)
order = pG.median(axis=1).sort_values().index
lim = 6
norm = TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim)
for k, (p, t, ttl) in enumerate([(pG, "a", "vs metaG"), (pT, "b", "vs metaT")]):
    ax = fig.add_subplot(gs[0, k]); tag(ax, t)
    m = p.loc[order].clip(-lim, lim)
    im = ax.imshow(np.ma.masked_invalid(m.values), cmap=DIV, norm=norm, aspect="auto")
    ax.set_xticks(range(len(S)), [short(s) for s in S], fontsize=6.5)
    if k == 0:
        ax.set_yticks(range(len(order)), [f"{gname}  ({phy.get(gname, '')})" for gname in order], fontsize=6.3)
    else:
        ax.set_yticks(range(len(order)), [""] * len(order))
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_xticks(np.arange(-.5, len(S)), minor=True); ax.set_yticks(np.arange(-.5, len(order)), minor=True)
    ax.grid(which="minor", color="white", lw=1.5); ax.tick_params(which="minor", length=0)
    ax.set_title(f"{T} / {ttl.split()[-1]}", fontsize=8.5)
cax = fig.add_axes([0.80, 0.62, 0.015, 0.22])
cb = fig.colorbar(im, cax=cax); cb.set_label(f"log2 CPM ratio (genus, MAGs summed)\nclipped at ±{lim}; blank = genus not binned", fontsize=6.5)
cb.outline.set_visible(False); cb.ax.tick_params(labelsize=6.5, length=0)

# G tRNA fragment classes
ax = fig.add_subplot(gs[1, 0]); tag(ax, "c"); ax.set_position([0.08, 0.08, 0.27, 0.32])
cls = ["5p_tRF", "5p_half", "5p_long", "i_tRF", "3p_half", "3p_tRF", "3p_long"]
lab = {"5p_tRF": "5′ tRF", "5p_half": "5′ half", "5p_long": "5′ long", "i_tRF": "i-tRF", "3p_half": "3′ half",
       "3p_tRF": "3′ tRF", "3p_long": "3′ long"}
p = trf.pivot(index="sample", columns="class", values="reads").reindex(index=S)[cls]
p = p.div(p.sum(axis=1), axis=0) * 100
bottom = np.zeros(len(S))
for c, col in zip(cls, CAT):
    ax.bar(range(len(S)), p[c], bottom=bottom, color=col, width=0.7, edgecolor="white", linewidth=1, label=lab[c])
    for i, (b, v) in enumerate(zip(bottom, p[c])):
        if v >= 9:
            ax.text(i, b + v / 2, f"{v:.0f}", ha="center", va="center", fontsize=6,
                    color="#0b0b0b" if col == CAT[3] else "white")
    bottom += p[c].values
ax.set_xticks(range(len(S)), [short(s) for s in S]); ax.set_ylim(0, 100); ax.set_ylabel("% of microbial tRNA reads")
ax.set_title("tRNA fragment classes"); ax.tick_params(axis="x", length=0)
ax.legend(loc="upper left", bbox_to_anchor=(-0.05, -0.1), ncol=4, fontsize=6, handlelength=1, columnspacing=0.6)

# H new layers: small multiples
mets = [("ANTI", "antisense on CDS: metaT strand"),
        ("CRISPR arrays >=10 reads", "CRISPR arrays ≥10 reads"),
        ("spacers with reads", "crRNA spacers with reads"),
        ("novel intergenic loci", "novel intergenic loci")]
sub = gs[1, 1:].subgridspec(2, 2, hspace=0.75, wspace=0.4)
for i, (c, ttl) in enumerate(mets):
    ax = fig.add_subplot(sub[i // 2, i % 2])
    if i == 0:
        tag(ax, "d")
    if c == "ANTI":
        an = pd.read_csv(f"{a.tables}/hc_antisense_cds.tsv", sep="\t")
        if "metaT_strand_class" in an:
            col, cl = "metaT_strand_class", ["sRNA strand only", "both strands", "gene strand only", "low metaT"]
            ramp = [CAT[0], CAT[1], CAT[2], NEUTRAL]
            lab = {"sRNA strand only": "metaT only on the sRNA strand", "both strands": "metaT on both strands",
                   "gene strand only": "metaT only on the gene strand", "low metaT": "<10 metaT reads"}
        else:
            col, cl = "metaT_class", ["no metaT reads", "1-9 metaT reads", ">=10 metaT reads"]
            ramp = ["#c3c2b7", "#86b6ef", "#1c5cab"]; lab = {k: k.replace(">=", "≥") for k in cl}
        p = an.groupby(["sample", col]).antisense_reads.sum().unstack().reindex(index=S, columns=cl).fillna(0)
        p = p.div(p.sum(axis=1), axis=0) * 100
        bot = np.zeros(len(S))
        for k, cc in zip(cl, ramp):
            ax.bar(range(len(S)), p[k], bottom=bot, color=cc, width=0.65, edgecolor="white", linewidth=0.8, label=lab[k])
            for i2, (b2, v2) in enumerate(zip(bot, p[k])):
                if v2 >= 12:
                    ax.text(i2, b2 + v2 / 2, f"{v2:.0f}", ha="center", va="center", fontsize=5.5,
                            color="#0b0b0b" if cc == NEUTRAL else "white")
            bot += p[k].values
        ax.set_ylim(0, 100); ax.set_ylabel("% of antisense reads", fontsize=6.5)
        fig.legend(*ax.get_legend_handles_labels(), loc="lower left", bbox_to_anchor=(0.43, 0.0), ncol=2, fontsize=6,
                   handlelength=1, columnspacing=0.8, title="antisense sRNA on CDS (≥10 reads/CDS):", title_fontsize=6)
    else:
        v = lay.set_index("sample").loc[S, c]
        ax.bar(range(len(S)), v, color=CAT[0], width=0.65)
    ax.set_xticks(range(len(S)), [short(s) for s in S], fontsize=6.3); ax.tick_params(axis="x", length=0)
    ax.set_title(ttl, fontsize=7.5, fontweight="normal", color=INK)
    ax.yaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True); ax.tick_params(axis="y", labelsize=6.3)
fig.suptitle("Expobiome healthy participants (n = 5): recurrent patterns and additional layers", y=0.985, x=0.01, ha="left",
             fontsize=10, fontweight="bold")
for e in ("png", "pdf"):
    fig.savefig(f"{a.out}.fig5.{e}", dpi=300)
plt.close(fig)
print("written", a.out + ".fig4/.fig5")
