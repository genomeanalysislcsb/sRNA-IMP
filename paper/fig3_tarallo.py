#!/usr/bin/env python3
"""fig3_tarallo.py - Fig. 3 of the sRNA-IMP paper: reanalysis of Tarallo et al. 2019 (stool small RNA, HC / adenoma / CRC).

Usage: fig3_tarallo.py OUTDIR OUTPREFIX [--type sRNA] [--cohort-dir DIR] [--compare-outdir OUTDIR2] [--compare-label "earlier strict"]

OUTDIR is an sRNA-IMP output folder after --cohort (docs/tutorial_tarallo.md). Read from it:
  Cohort/tables/cohort.sample_stats.tsv          groups, QC status, human_fraction (panel a)
  Cohort/tables/cohort.<type>.mirna.counts.tsv   human miRNAs (panels d, e)
  Cohort/tables/cohort.<type>.genus.counts.tsv   rRNA genera, Bracken (panel c)
  <ID>/Analysis/<type>/nonhuman/trna/<ID>.trna_species.counts.tsv   species from microbial tRNA reads (panel b)
Libraries excluded by the cohort QC are left out, as in the cohort statistics. Groups: HC, Adenoma, CRC.
--cohort-dir: another cohort folder of the same run (default OUTDIR/Cohort), e.g. one written with
COHORT_QC_DROP=none and COHORT_DIR=OUTDIR/Cohort_all80 to show all libraries.

Panels: a host load; b tRNA species with the change reported by Tarallo et al.; c Agathobacter among rRNA reads;
d human miRNA reads per library (with --compare-outdir: this run vs another run of the same libraries, e.g.
different host settings); e the published miRNA markers.
Tests (the cohort module's own implementations): Mann-Whitney vs HC, Spearman trend HC < adenoma < CRC, Kruskal-Wallis.
Writes OUTPREFIX.png, OUTPREFIX.pdf and OUTPREFIX.tsv (every number and p-value shown). Needs numpy and matplotlib
(the sRNA-IMP environment or envs/paper.yml)."""
import argparse
import csv
import os
import re
import sys
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
from srna_cohort import mannwhitney, spearman, kruskal  # noqa: E402

GO = ["HC", "Adenoma", "CRC"]


def read_rows(p):
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def matrix(p, samples):
    """feature -> np.array over samples (missing columns = 0)"""
    out = {}
    with open(p) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        idx = [head.index(s) if s in head else -1 for s in samples]
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out[f[0]] = np.array([num(f[i]) if i >= 0 else 0.0 for i in idx])
    return out


def colsum(m, n):
    return np.sum(list(m.values()), axis=0) if m else np.zeros(n)


def load_run(outdir, typ, cohort_dir=""):
    tab = os.path.join(cohort_dir or os.path.join(outdir, "Cohort"), "tables")
    st = read_rows(os.path.join(tab, "cohort.sample_stats.tsv"))
    types = sorted({r["type"] for r in st})
    if typ is None:
        if len(types) != 1:
            sys.exit("several library types in %s (%s): give --type" % (outdir, ", ".join(types)))
        typ = types[0]
    st = [r for r in st if r["type"] == typ and r.get("in_statistics", "yes") == "yes" and r["group"] in GO]
    samples = [r["sample"] for r in st]
    return {"type": typ, "samples": samples, "group": np.array([r["group"] for r in st]),
            "human_pct": np.array([num(r.get("human_fraction")) * 100 for r in st]),
            "mirna": matrix(os.path.join(tab, "cohort.%s.mirna.counts.tsv" % typ), samples),
            "genus": matrix(os.path.join(tab, "cohort.%s.genus.counts.tsv" % typ), samples)}


def trna_species(outdir, typ, samples):
    """species (first two words of the tRNA reference organism) -> counts over samples"""
    out = defaultdict(lambda: np.zeros(len(samples)))
    for j, s in enumerate(samples):
        p = os.path.join(outdir, s, "Analysis", typ, "nonhuman", "trna", "%s.trna_species.counts.tsv" % s)
        if not os.path.exists(p):
            print("missing", p, file=sys.stderr)
            continue
        for r in read_rows(p):
            name = re.sub(r"\s*\[.*\]$", "", r["species"])
            out[" ".join(name.split()[:2])][j] += num(r["count"])
    return dict(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir"); ap.add_argument("out")
    ap.add_argument("--type", default=None, help="library type (column 'type' of the sample sheet); default: the only one")
    ap.add_argument("--cohort-dir", default="", help="cohort folder of OUTDIR to read (default OUTDIR/Cohort)")
    ap.add_argument("--compare-outdir", default="", help="second run of the same libraries for panel d")
    ap.add_argument("--compare-label", default="earlier strict", help="x label of the second run in panel d")
    ap.add_argument("--label", default="final defaults", help="x label of this run in panel d (with --compare-outdir)")
    a = ap.parse_args()

    CAT = ["#2a78d6", "#eb6834", "#1baf7a"]
    INK, INK2, MUTED, GRID, AXIS = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.edgecolor": AXIS, "axes.labelcolor": INK2,
                         "xtick.color": MUTED, "ytick.color": MUTED, "axes.titlesize": 8.5, "axes.titleweight": "bold",
                         "axes.titlelocation": "left", "axes.spines.top": False, "axes.spines.right": False,
                         "legend.frameon": False, "legend.fontsize": 7})
    GC = {"HC": CAT[0], "Adenoma": CAT[2], "CRC": CAT[1]}

    R = load_run(a.outdir, a.type, a.cohort_dir)
    S, G = R["samples"], R["group"]
    order = np.array([GO.index(g) for g in G], float)
    n = len(S)
    rows = []

    def tag(ax, t):
        ax.text(-0.2, 1.06, t, transform=ax.transAxes, fontsize=11, fontweight="bold", ha="right", va="bottom")

    def pfmt(p):
        return "p = n/a" if p != p else ("p = %.1g" % p if p >= 0.001 else "p = %.0e" % p)

    def by(v, g):
        x = v[G == g]
        return x[~np.isnan(x)]

    def strip(ax, v, ylab, name, log=False):
        rng = np.random.default_rng(1)
        for i, g in enumerate(GO):
            x = by(v, g)
            ax.scatter(i + rng.uniform(-0.18, 0.18, len(x)), x, s=12, color=GC[g], alpha=0.75, edgecolor="white",
                       linewidth=0.3, zorder=3)
            ax.plot([i - 0.28, i + 0.28], [np.median(x)] * 2, color=INK, lw=1.6, zorder=4)
        ax.set_xticks(range(3), ["%s\n(n=%d)" % (g, len(by(v, g))) for g in GO], fontsize=7)
        ax.tick_params(axis="x", length=0)
        if log:
            ax.set_yscale("log")
            lo, hi = ax.get_ylim()
            ax.set_ylim(lo, hi * (hi / lo) ** 0.15)          # room for the p-values above the points
        ax.set_ylabel(ylab); ax.yaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
        pk = kruskal([by(v, g) for g in GO])[1]
        pt = spearman(order, v)[1]
        pc = mannwhitney(by(v, "CRC"), by(v, "HC"))[1]
        ax.text(0.02, 0.98, "CRC vs HC %s\ntrend %s" % (pfmt(pc), pfmt(pt)), transform=ax.transAxes, va="top",
                fontsize=6.5, color=INK2)
        rows.append((name, "median HC/Adenoma/CRC", "/".join("%.3g" % np.median(by(v, g)) for g in GO), pc, pt, pk))

    def foldchanges(ax, vals, labels, name, pseudo, pcol_x=None):
        """vals: list of arrays (one per feature, top to bottom); log2 mean ratio vs HC per group"""
        yy = np.arange(len(vals))[::-1]
        lfs, pcol = [], []
        for k, x in enumerate(vals):
            if x is None:
                continue
            for g, mk in (("Adenoma", "o"), ("CRC", "s")):
                lf = np.log2((by(x, g).mean() + pseudo) / (by(x, "HC").mean() + pseudo))
                pv = mannwhitney(by(x, g), by(x, "HC"))[1]
                sig = pv < 0.05 and pcol_x is None
                ax.scatter(lf, yy[k] + (0.15 if g == "Adenoma" else -0.15), s=38, marker=mk, color=GC[g],
                           edgecolor=INK if sig else "white", linewidth=0.9 if sig else 0.4, zorder=3,
                           label=g + " vs HC" if k == 0 else None)
                lfs.append(lf)
                if pcol_x is not None and g == "CRC":
                    pcol.append((yy[k] - 0.15, pfmt(pv)))
                rows.append((name + labels[k][0], g + "_vs_HC_log2FC", "%.2f" % lf, pv, spearman(order, x)[1],
                             float("nan")))
        ax.axvline(0, color=AXIS, lw=0.8)
        ax.set_yticks(yy, [lab for _, lab in labels], fontsize=6.8)
        if pcol_x is not None:                      # x range from the data, p values in a column right of it
            lo, hi = min([-2.0] + lfs) - 0.5, max([2.0] + lfs) + 0.5
            px = hi + 0.1 * (hi - lo)
            ax.set_xlim(lo, px + 0.32 * (hi - lo))
            step = 1 if hi - lo <= 8 else 2
            ax.set_xticks([t for t in range(int(np.ceil(lo)), int(np.floor(hi)) + 1) if t % step == 0])
            for y, t in pcol:
                ax.text(px, y, t, fontsize=6.3, va="center", color=INK2)
            ax.text(px, len(vals) - 1 + 0.32, "CRC vs HC", fontsize=6.3, color=INK2)
            ax.axvline(hi + 0.05 * (hi - lo), color=GRID, lw=0.5)
        ax.xaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True); ax.tick_params(axis="y", length=0)

    fig = plt.figure(figsize=(7.2, 6.8))
    gs = fig.add_gridspec(2, 3, hspace=0.62, wspace=0.62, left=0.09, right=0.98, top=0.88, bottom=0.1,
                          width_ratios=[1, 1, 1.25])
    gsb = gs[0, 1:].subgridspec(1, 2, width_ratios=[0.28, 1])
    gsl = gs[1, :].subgridspec(1, 3, width_ratios=[1, 0.9, 1.2], wspace=0.85)

    # a host load
    ax = fig.add_subplot(gs[0, 0]); tag(ax, "a")
    strip(ax, R["human_pct"], "host reads (% of trimmed)", "host_pct", log=True)
    ax.set_title("Host load")

    # b species from microbial tRNA reads (published microbial findings)
    ax = fig.add_subplot(gsb[0, 1]); tag(ax, "b")
    T = trna_species(a.outdir, R["type"], S)
    tot = colsum(T, n)
    sp = [("Escherichia coli", "↑"), ("Klebsiella pneumoniae", "↑"), ("Fusobacterium nucleatum", "↑ n.s."),
          ("Bacteroides ovatus", "↓"), ("Agathobacter rectalis", "↓")]
    vals = [T[s] / np.where(tot > 0, tot, np.nan) if s in T else None for s, _ in sp]
    foldchanges(ax, vals, [(s, "%s. %s (%s)" % (s[0], s.split()[1], p)) for s, p in sp], "tRNA_", 1e-9, pcol_x=True)
    ax.set_xlabel("log2 fold change vs HC (share of microbial tRNA reads)")
    ax.set_title("Species from microbial tRNA reads")
    ax.legend(loc="upper left", fontsize=6.5, handletextpad=0.2)
    ax.text(-0.45, -0.3, "b, e: in brackets the change reported by Tarallo et al.; e: outlined p < 0.05",
            transform=ax.transAxes, fontsize=6, color=INK2)

    # c Agathobacter among rRNA reads (Bracken genus)
    ax = fig.add_subplot(gsl[0, 0]); tag(ax, "c")
    gt = colsum(R["genus"], n)
    ag = R["genus"].get("Agathobacter", np.zeros(n)) / np.where(gt > 0, gt, np.nan) * 100
    strip(ax, ag + 0.01, "Agathobacter (% of rRNA reads)\nBracken; 0 shown at 0.01", "Agathobacter_rRNA_pct",
          log=True)             # libraries without classified rRNA reads are left out (n below the group)
    ax.set_title("Agathobacter (rRNA)")

    # d miRNA depth
    ax = fig.add_subplot(gsl[0, 1]); tag(ax, "d")
    Mn = R["mirna"]
    sn = colsum(Mn, n)
    if a.compare_outdir:
        C = load_run(a.compare_outdir, R["type"])
        so_all = dict(zip(C["samples"], colsum(C["mirna"], len(C["samples"]))))
        keep = np.array([s in so_all for s in S])
        so = np.array([so_all.get(s, np.nan) for s in S])
        for j in np.where(keep)[0]:
            ax.plot([0, 1], [so[j] + 1, sn[j] + 1], color=GC[G[j]], lw=0.6, alpha=0.6)
        ax.scatter(np.zeros(keep.sum()), so[keep] + 1, s=8, color=[GC[g] for g in G[keep]], zorder=3)
        ax.scatter(np.ones(keep.sum()), sn[keep] + 1, s=8, color=[GC[g] for g in G[keep]], zorder=3)
        for x0, v in ((0, so[keep]), (1, sn[keep])):
            ax.plot([x0 - 0.18, x0 + 0.18], [np.median(v) + 1] * 2, color=INK, lw=1.6)
            ax.text(x0, v.max() * 2.2, "median\n%.0f" % np.median(v), ha="center", fontsize=6.3, color=INK2)
        lab = lambda s: s.replace(" ", "\n", 1)  # noqa: E731
        ax.set_xticks([0, 1], [lab(a.compare_label), lab(a.label)], fontsize=7); ax.set_xlim(-0.5, 1.5)
        ax.tick_params(axis="x", length=0)
        ax.set_yscale("log"); ax.set_ylabel("human miRNA reads per library")
        ax.yaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
        ax.set_ylim(top=ax.get_ylim()[1] * 4)
        rows.append(("mirna_reads_per_library", "median %s/%s" % (a.compare_label, a.label),
                     "%.0f/%.0f" % (np.median(so[keep]), np.median(sn[keep])),
                     float("nan"), float("nan"), float("nan")))
    else:
        strip(ax, sn + 1, "human miRNA reads per library", "mirna_reads_per_library", log=True)
    ax.set_title("miRNA depth")

    # e published miRNA markers (CPM; miR-1246, the dominant miRNA, is left out of the CPM denominator of the others)
    ax = fig.add_subplot(gsl[0, 2]); tag(ax, "e")
    mk = [("hsa-miR-200b-3p", "↑"), ("hsa-miR-30c-5p", "adenoma"), ("hsa-miR-21-5p", "changed"),
          ("hsa-miR-378a-3p", "changed"), ("hsa-miR-215-5p", "changed"), ("hsa-miR-1246", "↑"), ("hsa-miR-1290", "↑")]
    den = sn - Mn.get("hsa-miR-1246", np.zeros(n))
    vals = []
    for f, _ in mk:
        if f not in Mn:
            vals.append(None)
        else:
            d = sn if f == "hsa-miR-1246" else den
            vals.append(Mn[f] / np.where(d > 0, d, np.nan) * 1e6)
    foldchanges(ax, vals, [(f, "%s (%s)" % (f.replace("hsa-", ""), p)) for f, p in mk], "miRNA_", 1.0)
    ax.set_xlabel("log2 fold change vs HC (CPM)"); ax.set_xlim(-2.2, 2.2)
    ax.set_title("Published miRNA markers", loc="right")
    ax.legend(loc="upper right", bbox_to_anchor=(1.02, -0.2), ncol=2, fontsize=6.3, handletextpad=0.2)

    fig.suptitle("Reanalysis of Tarallo et al. 2019: %d stool small RNA libraries" % n, x=0.01, ha="left",
                 fontsize=9.5, fontweight="bold")
    fig.legend([plt.Line2D([], [], ls="", marker="o", color=GC[g], ms=5) for g in GO], GO, loc="upper right",
               bbox_to_anchor=(0.98, 0.985), ncol=3, fontsize=7, handletextpad=0.1)
    for ext in ("png", "pdf"):
        fig.savefig("%s.%s" % (a.out, ext), dpi=300)
    with open(a.out + ".tsv", "w") as out:
        out.write("feature\tstatistic\tvalue\tp_vs_HC\tp_trend\tp_kruskal\n")
        for r in rows:
            out.write("\t".join(r[:3]) + "\t" + "\t".join("%.3g" % x for x in r[3:]) + "\n")
    print("written %s.{png,pdf,tsv} (%d libraries: %s)" % (a.out, n, ", ".join("%s %d" % (g, (G == g).sum()) for g in GO)))


if __name__ == "__main__":
    main()
