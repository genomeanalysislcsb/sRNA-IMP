#!/usr/bin/env python3
"""fig1_workflow.py — Fig. 1: the sRNA-IMP workflow (modules, inputs, host and microbial routes, outputs).
Usage: fig1_workflow.py OUTPREFIX [--panel-b | --supp]   writes OUTPREFIX.pdf and OUTPREFIX.png (170 mm wide)
  --panel-b   add panel b: three read-level sketches (host split with contig rescue, placement on the assembly,
              strand check with stranded metaT)
  --supp      only these sketches, as Supplementary Fig. S1 (170 x 66 mm)"""
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

args = [x for x in sys.argv[1:] if not x.startswith("--")]
out = args[0] if args else "fig1_workflow"
SUPP = "--supp" in sys.argv          # the read-level sketches alone (Supplementary Fig. S1)
PANEL_B = "--panel-b" in sys.argv or SUPP
BH = 66 if PANEL_B else 0          # height of panel b (mm), drawn below panel a at negative y
AH = 0 if SUPP else 150
ACC, INK, INK2, EDGE, TINT = "#2a78d6", "#0b0b0b", "#52514e", "#898781", "#f0efec"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.5})

W, H = 170 / 25.4, (AH + BH) / 25.4
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 170); ax.set_ylim(-BH, AH); ax.axis("off")
ORANGE = "#eb6834"


def box(x, y, w, h, title, lines=(), kind="plain", size=7.5):
    """x, y = centre (mm); kind: plain | accent | input | note"""
    fc, ec, lw, ls = "white", EDGE, 0.8, "-"
    if kind == "accent":
        fc, ec, lw = "#e3eefb", ACC, 1.4
    elif kind == "input":
        fc = TINT
    elif kind == "note":
        ls = (0, (3, 2))
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0,rounding_size=2",
                                fc=fc, ec=ec, lw=lw, ls=ls, zorder=2))
    n = len(lines)
    top = y + (n * 3.4) / 2
    ax.text(x, top, title, ha="center", va="center", fontsize=size, fontweight="bold", color=INK, zorder=3)
    for i, ln in enumerate(lines):
        ax.text(x, top - 3.6 * (i + 1), ln, ha="center", va="center", fontsize=size - 1, color=INK2, zorder=3)


def arrow(pts, label="", lpos=None, ha="left", color=EDGE):
    """polyline through pts (mm), arrow head at the end"""
    for (x0, y0), (x1, y1) in zip(pts[:-2], pts[1:-1]):
        ax.plot([x0, x1], [y0, y1], color=color, lw=0.8, zorder=1, solid_capstyle="round")
    ax.add_patch(FancyArrowPatch(pts[-2], pts[-1], arrowstyle="-|>", mutation_scale=7, color=color, lw=0.8,
                                 shrinkA=0, shrinkB=0, zorder=1))
    if label:
        ax.text(*lpos, label, ha=ha, va="center", fontsize=6.3, color=INK2, style="italic", zorder=3)


def panel_a():
    # inputs
    box(30, 138, 52, 12, "Small RNA reads", ["one library per sample and type"], "input")
    box(92, 138, 56, 12, "Sample assembly", ["contigs, annotation, bins, taxonomy"], "input")
    box(146, 138, 44, 12, "metaG / metaT", ["counts; stranded metaT BAM"], "input")

    # per library: QC and host split
    box(30, 118, 52, 12, "qc, contam", ["trimming, contaminants, QC verdict"])
    box(30, 96, 52, 20, "split: host or microbial", ["exact matches, ≥ 20 nt: genome,", "transcripts, mature tRNA + CCA;",
                                                     "mature miRNAs with 3′ tails"], "accent")
    arrow([(30, 132), (30, 124)])
    arrow([(30, 112), (30, 106)])
    arrow([(72, 132), (72, 101), (56, 101)], "contig rescue", (73.5, 116))

    # host and microbial classification
    box(17, 70, 30, 16, "human", ["miRNA, tRNA, rRNA,", "other ncRNA"])
    box(50, 70, 32, 16, "nonhuman", ["RNA classes, rRNA", "taxonomy, miRNA-like"])
    box(33, 50, 48, 11, "mirna_origin", ["host or microbial miRNA-like reads"])
    arrow([(22, 86), (22, 81), (17, 81), (17, 78)], "host reads", (24, 82.5))
    arrow([(56, 91), (64, 91), (64, 81), (50, 81), (50, 78)])
    arrow([(17, 62), (17, 56), (22, 56), (22, 55.5)])
    arrow([(50, 62), (50, 56), (44, 56), (44, 55.5)])

    # assembly route
    box(118, 96, 64, 16, "assembly", ["reads per feature, sense or antisense,", "contig, MAG and functional term"], "accent")
    arrow([(56, 91), (86, 91)], "microbial reads", (66, 88.5))
    arrow([(112, 132), (112, 104)])
    box(118, 70, 64, 16, "novel, fragments, bins, crispr", ["novel loci (Rfam, RNAfold), tRNA / rRNA fragments,",
                                                             "reads per MAG, crRNAs per CRISPR spacer"])
    arrow([(118, 88), (118, 78)])
    box(118, 46, 64, 16, "--compare", ["sRNA vs metaG and metaT per feature,", "contig and MAG; metaT strand of antisense loci"])
    arrow([(118, 62), (118, 54)])
    arrow([(160, 132), (160, 46), (150, 46)])

    # host-only libraries
    box(66, 30, 46, 15, "host-only libraries", ["(e.g. blood, no assembly):", "qc, contam, split, human, nonhuman"], "note",
        size=7)

    # reports and cohort
    box(85, 9, 158, 12, "--report  and  --cohort", ["per-sample report with QC verdict and generated Methods; cohort: QC exclusion, "
                                                    "group tests, models adjusted for host load, MultiQC"])
    arrow([(8, 62), (8, 15)])
    arrow([(33, 44.5), (33, 15)])
    arrow([(118, 38), (118, 15)])


if not SUPP:
    panel_a()

if PANEL_B:
    from matplotlib.patches import FancyArrow, Rectangle
    if not SUPP:
        ax.text(1, 147, "a", fontsize=11, fontweight="bold", va="center")
        ax.text(1, -4, "b", fontsize=11, fontweight="bold", va="center")

    def line(x0, x1, y, label="", c=EDGE, lw=1.6):
        ax.plot([x0, x1], [y, y], color=c, lw=lw, solid_capstyle="butt", zorder=2)
        if label:
            ax.text(x0, y + 2.3, label, fontsize=6.3, color=INK2, va="center")

    def read(x0, x1, y, c):
        ax.add_patch(Rectangle((x0, y - 0.8), x1 - x0, 1.6, fc=c, ec="none", zorder=4))

    def gene(x0, x1, y, strand, label):
        L = x1 - x0
        ax.add_patch(FancyArrow(x0 if strand > 0 else x1, y, (L - 0.01) * strand, 0, width=3.2, head_width=4.6,
                                head_length=3.2, length_includes_head=True, fc="white", ec=INK2, lw=0.7, zorder=3))
        ax.text((x0 + x1) / 2, y, label, fontsize=6, ha="center", va="center", color=INK, zorder=5)

    def title(x, t):
        ax.text(x, -6.5, t, fontsize=7.2, fontweight="bold", color=INK, va="center")

    def note(x, y, t, c=INK2, w="normal"):
        ax.text(x, y, t, fontsize=6.3, color=c, va="center", fontweight=w)

    # b1 host or microbial
    x0 = 5
    title(x0, "1  Host or microbial?")
    line(x0, x0 + 50, -17, "human genome")
    line(x0, x0 + 50, -30, "contig of the same sample")
    read(x0 + 34, x0 + 44, -17, INK); read(x0 + 34, x0 + 44, -30, INK)
    ax.plot([x0 + 39, x0 + 39], [-18, -29], color=INK2, lw=0.6, ls=(0, (2, 1.5)), zorder=1)
    note(x0, -23.5, "a 24-nt read, exact on both", INK2)
    note(x0, -39, "genome first (sMETASeq):  human")
    note(x0, -44.5, "sRNA-IMP:  microbial (contig rescue)", ACC, "bold")
    note(x0, -52, "host only with an exact host match,\n\u2265 20 nt and no equally good contig hit")

    # b2 placement
    x0 = 60
    title(x0, "2  Placement on the assembly")
    line(x0, x0 + 48, -25)
    gene(x0 + 3, x0 + 24, -25, +1, "gene A")
    gene(x0 + 32, x0 + 47, -25, -1, "gene B")
    read(x0 + 8, x0 + 16, -19.5, ACC); note(x0 + 8, -16.8, "sense", ACC, "bold")
    read(x0 + 9, x0 + 17, -30.5, ORANGE); note(x0 + 9, -33.4, "antisense", ORANGE, "bold")
    read(x0 + 25.5, x0 + 30.5, -19.5, "#898781"); note(x0 + 23, -16.8, "intergenic", INK2, "bold")
    note(x0, -39.5, "contig \u2192 MAG (GTDB taxonomy)", INK, "bold")
    note(x0, -45, "reads per feature, strand, contig, MAG")
    note(x0, -50, "and functional term; novel loci,")
    note(x0, -55, "tRNA / rRNA fragments, crRNAs")

    # b3 strand check
    x0 = 113
    title(x0, "3  Strand check with metaT")
    rows = [(-18, "sRNA strand only", (False, True)), (-32, "both strands", (True, True)),
            (-46, "gene strand only", (True, False))]
    for n_, (y, lab, (up, dn)) in enumerate(rows):
        ax.plot([x0, x0 + 28], [y, y], color=EDGE, lw=0.8, zorder=1)
        ax.add_patch(FancyArrow(x0 + 3, y, 20, 0, width=2.4, head_width=3.6, head_length=2.6,
                                length_includes_head=True, fc="white", ec=INK2, lw=0.6, zorder=3))
        read(x0 + 7, x0 + 15, y - 3.2, ORANGE)                         # the antisense small RNA
        for k in range(3):
            if up:
                read(x0 + 5 + 6 * k, x0 + 9 + 6 * k, y + 3.2, ACC)         # metaT reads, gene strand
            if dn:
                read(x0 + 17 + 3.2 * k, x0 + 19.4 + 3.2 * k, y - 3.2, ORANGE)  # metaT reads, sRNA strand
        note(x0 + 30, y, lab, INK, "bold")
    note(x0, -53.5, "long orange bar: the antisense small RNA;")
    note(x0, -57.5, "short bars: metaT reads on the sRNA strand")
    note(x0, -61.5, "(orange) or the gene strand (blue)")

for ext in ("pdf", "png"):
    fig.savefig("%s.%s" % (out, ext), dpi=300)
print("written", out)
