#!/usr/bin/env python3
"""figS2_crispr_locus.py - one transcribed CRISPR-Cas locus with crRNA reads per spacer (Supplementary Fig. S2).

Usage: figS2_crispr_locus.py OUTDIR OUTPREFIX --samples BHC01,BHC02,... [--type sRNA-Int] [--array ID]
                             [--labels BHC01=P1,BHC02=P2,...]   (names shown in the figure; default the IDs)

OUTDIR is an sRNA-IMP output folder; read from <ID>/Analysis/<type>/crispr/ (srna_crispr.py):
  crispr_arrays.tsv, crispr_units.tsv, crispr_genes.tsv, crispr_coverage.tsv, crispr_ends.tsv, crispr_tracr.tsv
The array is --array, else the binned array with Cas genes and the most reads in the first sample.
Panels: a locus (Cas genes, tracrRNA candidate, array) with small RNA coverage per strand; b reads per spacer in
transcription order; c 5' and 3' ends of spacer reads relative to the spacer; d spacers of arrays of the same species
with the same repeat (either orientation) in the other samples, and how many are shared.
Writes OUTPREFIX.png / .pdf (170 mm wide) and OUTPREFIX.tsv (the plotted numbers)."""
import argparse
import csv
import os
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, Rectangle

ACC, ORANGE, GREEN = "#2a78d6", "#eb6834", "#1baf7a"
SAMPLE_COL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6"]
INK, INK2, MUTED, GRID, AXIS = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.5, "axes.edgecolor": AXIS, "axes.labelcolor": INK2,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.titlesize": 8, "axes.titleweight": "bold",
                     "axes.titlelocation": "left", "axes.spines.top": False, "axes.spines.right": False,
                     "legend.frameon": False, "legend.fontsize": 6.5})


def rows(outdir, s, typ, kind):
    p = os.path.join(outdir, s, "Analysis", typ, "crispr", "%s.crispr_%s.tsv" % (s, kind))
    if not os.path.exists(p):
        return []
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return 0.0


def rc(x):
    return x[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


def same_repeat(a, b, min_id=0.97):
    """same repeat in either orientation: >= min_id identity (at most 1 mismatch in a 36-nt repeat) over the shorter one, ends allowed to differ by 1 nt"""
    for q in (b, rc(b)):
        for sh in (-1, 0, 1):
            x, y = (a[sh:], q) if sh >= 0 else (a, q[-sh:])
            n = min(len(x), len(y))
            if n >= 20 and sum(c1 == c2 for c1, c2 in zip(x[:n], y[:n])) / n >= min_id:
                return True
    return False


def tag(ax, t, x=-0.08):
    ax.text(x, 1.04, t, transform=ax.transAxes, fontsize=11, fontweight="bold", ha="right", va="bottom")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir"); ap.add_argument("out")
    ap.add_argument("--samples", required=True); ap.add_argument("--type", default="sRNA-Int")
    ap.add_argument("--array", default="")
    ap.add_argument("--labels", default="", help="ID=label,... shown instead of the sample IDs")
    a = ap.parse_args()
    LAB = dict(x.split("=", 1) for x in a.labels.split(",") if "=" in x)
    S = a.samples.split(",")
    s0 = S[0]
    arrays = rows(a.outdir, s0, a.type, "arrays")
    if a.array:
        ar = next(r for r in arrays if r["array_id"] == a.array)
    else:
        cand = [r for r in arrays if r["cas_type"] not in ("", "no cas genes") and r["bin"] not in ("", "unbinned")]
        ar = max(cand, key=lambda r: num(r["reads"]))
    aid, contig = ar["array_id"], ar["contig"]
    taxon = ar["classification"].split(";")[-1].replace("s__", "").replace("_", " ") or ar["bin"]
    units = [u for u in rows(a.outdir, s0, a.type, "units") if u["array_id"] == aid]
    genes = [g for g in rows(a.outdir, s0, a.type, "genes") if g["array_id"] == aid]
    cov = [c for c in rows(a.outdir, s0, a.type, "coverage") if c["array_id"] == aid]
    ends = [e for e in rows(a.outdir, s0, a.type, "ends") if e["array_id"] == aid]
    tracr = [t for t in rows(a.outdir, s0, a.type, "tracr") if t["array_id"] == aid and t.get("tracr_candidate") == "yes"]
    out_rows = []

    fig = plt.figure(figsize=(170 / 25.4, 110 / 25.4))
    gs = fig.add_gridspec(2, 3, height_ratios=[0.9, 1], hspace=0.45, wspace=0.6, left=0.08, right=0.97,
                          top=0.92, bottom=0.1, width_ratios=[1.2, 1, 1])

    # a locus
    ax = fig.add_subplot(gs[0, :])
    tag(ax, "a", -0.03)
    pos = np.array([int(c["pos"]) for c in cov]); plus = np.array([num(c["plus"]) for c in cov])
    minus = np.array([num(c["minus"]) for c in cov])
    lo = min([pos.min()] + [int(g["start"]) for g in genes]); hi = max([pos.max(), int(ar["end"])])
    x = np.arange(lo, hi + 1); P = np.zeros(len(x)); M = np.zeros(len(x))
    P[pos - lo] = plus; M[pos - lo] = minus
    f = lambda v: np.log10(v + 1)  # noqa: E731
    ax.fill_between(x, 0, f(P), color=ACC, lw=0, step="mid", label="reads on + strand")
    ax.fill_between(x, 0, -f(M), color=ORANGE, lw=0, step="mid", label="reads on − strand")
    top, bot = max(f(P).max(), 1), max(f(M).max(), 1)
    yp = list(range(2, int(np.ceil(top)) + 1, 2)); ym = list(range(2, int(np.ceil(bot)) + 1, 2))
    ax.set_yticks([-t for t in ym[::-1]] + [0] + yp,
                  ["10$^{%d}$" % t for t in ym[::-1]] + ["0"] + ["10$^{%d}$" % t for t in yp], fontsize=6.3)
    ax.set_ylim(-bot - 2.6, top + 0.9)
    ax.axhline(0, color=AXIS, lw=0.6)
    gy = -bot - 1.4
    for g in genes:
        s_, e_ = int(g["start"]), int(g["end"]); st = 1 if g["strand"] == "+" else -1
        ax.add_patch(FancyArrow(s_ if st > 0 else e_, gy, (e_ - s_) * st, 0, width=0.7, head_width=0.95,
                                head_length=min(250, (e_ - s_) * 0.3), length_includes_head=True, fc="white",
                                ec=INK2, lw=0.7, zorder=3))
        ax.text((s_ + e_) / 2 - 60 * st, gy, g["name"], ha="center", va="center", fontsize=5.8, style="italic",
                zorder=4)
    for t in tracr:
        s_, e_ = int(t["start"]), int(t["end"])
        ax.add_patch(Rectangle((s_, gy - 0.35), e_ - s_, 0.7, fc=GREEN, ec="none", zorder=3))
        ax.text((s_ + e_) / 2, gy + 0.65, "tracrRNA", fontsize=6.3, color=INK2, ha="center", va="bottom")
    s_, e_ = int(ar["start"]), int(ar["end"])
    for u in units:
        us, ue = int(u["start"]), int(u["end"])
        ax.add_patch(Rectangle((us, gy - 0.35), ue - us, 0.7, fc=INK2 if u["kind"] == "repeat" else "#f0efec",
                               ec="none", zorder=3))
    ax.text((s_ + e_) / 2, gy - 0.65, "CRISPR array (%s repeats)" % ar["repeats"], ha="center", va="top",
            fontsize=6.3, color=INK2)
    ax.set_xlim(lo - 200, hi + 100)
    ax.set_xticks([])
    ax.plot([lo, lo + 1000], [gy - 1.0, gy - 1.0], color=INK2, lw=1)
    ax.text(lo + 500, gy - 1.2, "1 kb", ha="center", va="top", fontsize=6.3, color=INK2)
    ax.set_ylabel("reads", fontsize=6.5)
    ax.yaxis.set_label_coords(-0.05, 0.68)
    ax.spines["bottom"].set_visible(False); ax.tick_params(axis="x", length=0)
    ax.set_title("Type %s CRISPR-Cas locus of %s (%s, %s)" % (ar["cas_type"], taxon, LAB.get(s0, s0), a.type))
    ax.legend(loc="upper left", ncol=2, bbox_to_anchor=(0.02, 1.02))
    out_rows.append(("locus", "reads_plus_minus", "%.0f/%.0f" % (P.sum(), M.sum())))

    # b reads per spacer
    ax = fig.add_subplot(gs[1, 0]); tag(ax, "b", -0.17)
    sp = [u for u in units if u["kind"] == "spacer"]
    v = np.array([num(u["reads"]) for u in sp])
    ax.bar(np.arange(1, len(v) + 1), v + 1, color=ACC, width=0.75)
    ax.set_yscale("log"); ax.set_xlabel("spacer (transcription order)"); ax.set_ylabel("reads + 1")
    ax.set_ylim(0.8, (v.max() + 1) * 30)
    ax.yaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
    rep = sum(num(u["reads"]) for u in units if u["kind"] == "repeat")
    ax.set_title("crRNA reads per spacer")
    ax.text(0.02, 0.98, "%d of %d spacers with reads" % ((v >= 1).sum(), len(v)), transform=ax.transAxes,
            ha="left", va="top", fontsize=6.3, color=INK2)
    out_rows.append(("array", "reads_spacers/repeat_only", "%.0f/%.0f" % (v.sum(), rep)))
    for i, (u, n) in enumerate(zip(sp, v), 1):
        out_rows.append(("spacer_%d" % i, u["sequence"], "%.1f" % n))

    # c read ends
    ax = fig.add_subplot(gs[1, 1]); tag(ax, "c", -0.17)
    sl = int(ar["spacer_len"]); rl = int(ar["repeat_len"])
    e5 = defaultdict(float); e3 = defaultdict(float)
    for e in ends:
        (e5 if e["end"] == "5p" else e3)[int(e["offset"])] += num(e["reads"])
    t5, t3 = sum(e5.values()) or 1, sum(e3.values()) or 1
    xs5 = sorted(k for k in e5 if -10 <= k <= sl); xs3 = sorted(k for k in e3 if -sl <= k <= rl)
    ax.bar(xs5, [e5[k] / t5 * 100 for k in xs5], color=ACC, width=0.85, label="5′ ends")
    ax.bar([sl + k for k in xs3], [e3[k] / t3 * 100 for k in xs3], color=ORANGE, width=0.85, label="3′ ends")
    ymax = max([e5[k] / t5 * 100 for k in xs5] + [e3[k] / t3 * 100 for k in xs3]) * 1.45
    ax.set_ylim(0, ymax)
    tr = ax.get_xaxis_transform()
    ax.add_patch(Rectangle((0, -0.09), sl, 0.05, fc="#f0efec", ec=INK2, lw=0.5, transform=tr, clip_on=False))
    ax.add_patch(Rectangle((sl, -0.09), rl, 0.05, fc=INK2, ec=INK2, lw=0.5, transform=tr, clip_on=False))
    ax.text(sl / 2, -0.12, "spacer", ha="center", va="top", fontsize=6, color=INK2, transform=tr)
    ax.text(sl + rl / 2, -0.12, "repeat", ha="center", va="top", fontsize=6, color=INK2, transform=tr)
    ax.set_xlim(-5, sl + rl + 2); ax.set_xticks([]); ax.spines["bottom"].set_visible(False)
    ax.set_ylabel("% of spacer reads"); ax.set_title("crRNA ends")
    m5 = max(e5, key=e5.get); m3 = max(e3, key=e3.get)
    ax.legend(loc="upper center", handlelength=1, bbox_to_anchor=(0.52, 0.8))
    ax.text(0.02, 0.98, "5′ end: spacer nt %d\n3′ end: repeat nt %d" % (m5 + 1, m3),
            transform=ax.transAxes, va="top", fontsize=6, color=INK2)
    out_rows += [("ends", "main_5p_offset_in_spacer", str(m5)), ("ends", "main_3p_offset_in_repeat", str(m3))]

    # d same repeat in other samples
    ax = fig.add_subplot(gs[1, 2]); tag(ax, "d", -0.1)
    R = ar["repeat_consensus"]
    keyf = lambda q: min(q, rc(q))  # noqa: E731
    sets = []
    for s in S:
        same = [r["array_id"] for r in rows(a.outdir, s, a.type, "arrays")
                if r["repeat_consensus"] and same_repeat(R, r["repeat_consensus"])
                and r["classification"].split(";")[-1] == ar["classification"].split(";")[-1]]
        seqs = set()
        for u in rows(a.outdir, s, a.type, "units"):
            if u["array_id"] in same and u["kind"] == "spacer" and u["sequence"]:
                seqs.add(keyf(u["sequence"]))
        if seqs:
            sets.append((s, seqs, len(same)))
    allsp = sorted(set().union(*[x[1] for x in sets]))
    for i, (s, seqs, na) in enumerate(sets):
        ci = SAMPLE_COL[S.index(s) % len(SAMPLE_COL)]
        idx = [j for j, q in enumerate(allsp) if q in seqs]
        shared = [j for j in idx if sum(allsp[j] in x[1] for x in sets) > 1]
        ax.scatter(idx, [i] * len(idx), marker="|", s=40, color=ci, lw=1.2)
        if shared:
            ax.scatter(shared, [i] * len(shared), marker="|", s=40, color=INK, lw=1.2)
        ax.text(len(allsp) + 2, i, "%d" % len(seqs), va="center", fontsize=6.3, color=INK2)
        out_rows.append(("same_repeat_" + s, "arrays/spacers", "%d/%d" % (na, len(seqs))))
    n_shared = sum(1 for q in allsp if sum(q in x[1] for x in sets) > 1)
    ax.set_yticks(range(len(sets)), [LAB.get(x[0], x[0]) for x in sets]); ax.set_ylim(len(sets) - 0.4, -0.6)
    ax.set_xticks([]); ax.spines["bottom"].set_visible(False); ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(-1, len(allsp) + 8)
    ax.set_xlabel("%d distinct spacers, %s shared" % (len(allsp), n_shared or "none"))
    ax.set_title("Spacers per participant")
    out_rows.append(("same_repeat", "spacers_shared", str(n_shared)))

    for ext in ("png", "pdf"):
        fig.savefig("%s.%s" % (a.out, ext), dpi=300)
    with open(a.out + ".tsv", "w") as fh:
        fh.write("item\tkey\tvalue\n")
        for r in out_rows:
            fh.write("\t".join(r) + "\n")
    print("written", a.out, aid)


if __name__ == "__main__":
    main()
