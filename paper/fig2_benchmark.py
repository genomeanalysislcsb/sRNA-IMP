#!/usr/bin/env python3
"""fig2_benchmark.py — Fig. 2: simulation benchmark, sRNA-IMP (final defaults) vs sMETASeq.
Usage: fig2_benchmark.py OUTPREFIX --summary final='bench/sim/*[0-9].eval.summary.tsv' \
                                     --summary smetaseq='bench/sim/*.smetaseq.eval.summary.tsv'
       fig2_benchmark.py OUTPREFIX --table panel_eval.tsv
--summary RUN=GLOB: srna_sim_eval.py summaries (docs/benchmark.md) of one run; the library is the file name up to its
first dot and must contain _h<host per mille>_ and end in _s<seed> (e.g. SIM_h020_s1.eval.summary.tsv).
--table: the same rows gathered in one file: library, run, section, metric, value (no header).
--run / --cmp: the run names of sRNA-IMP and sMETASeq [final / smetaseq].
Writes OUTPREFIX.{png,pdf} and OUTPREFIX.tsv (the plotted numbers)."""
import argparse
import glob
import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("out")
ap.add_argument("--summary", action="append", default=[], metavar="RUN=GLOB")
ap.add_argument("--table", default="")
ap.add_argument("--run", default="final"); ap.add_argument("--cmp", default="smetaseq")
a = ap.parse_args()

CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK2, MUTED, GRID, AXIS = "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.edgecolor": AXIS, "axes.labelcolor": INK2,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.titlesize": 8.5, "axes.titleweight": "bold",
                     "axes.titlelocation": "left", "axes.spines.top": False, "axes.spines.right": False,
                     "legend.frameon": False, "legend.fontsize": 7})
TOOLS = {a.run: ("sRNA-IMP", CAT[0], "o"), a.cmp: ("sMETASeq", CAT[1], "s")}

COLS = ["lib", "run", "sec", "metric", "value"]
parts = []
if a.table:
    parts.append(pd.read_csv(a.table, sep="\t", header=None, names=COLS, dtype=str))
for spec in a.summary:
    run, pat = spec.split("=", 1)
    files = sorted(glob.glob(pat))
    if not files:
        raise SystemExit("no files for --summary %s" % spec)
    for f in files:
        d = pd.read_csv(f, sep="\t", dtype=str).rename(columns={"section": "sec"})
        parts.append(d.assign(lib=os.path.basename(f).split(".")[0], run=run)[COLS])
if not parts:
    raise SystemExit("give --summary RUN=GLOB (twice) or --table")
e = pd.concat(parts, ignore_index=True)
e = e[e.run.isin(TOOLS)]
e["host"] = e.lib.str.extract(r"_h(\d+)_")[0].astype(float) / 10
e["seed"] = e.lib.str.extract(r"_s(\d+)$")[0]
num = lambda s: float(re.match(r"\s*([-\d.eE+]+)", s).group(1))


def get(metric, run, f=num):
    x = e[(e.metric == metric) & (e.run == run)]
    return x.assign(v=x.value.map(f))[["lib", "host", "seed", "v"]]


def tag(ax, t):
    ax.text(-0.18, 1.06, t, transform=ax.transAxes, fontsize=11, fontweight="bold", ha="right", va="bottom")


def host_axis(ax):
    ax.set_xscale("log"); ax.set_xticks([0.5, 2, 5, 20], ["0.5", "2", "5", "20"]); ax.minorticks_off()
    ax.set_xlabel("host share of simulated reads (%)"); ax.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)


def lines(ax, metric, f=num, scale=1.0, ylab=""):
    rows = []
    for run, (name, col, mk) in TOOLS.items():
        d = get(metric, run, f)
        if d.empty:
            continue
        d = d.assign(v=d.v * scale)
        m = d.groupby("host").v.mean()
        ax.plot(m.index, m.values, color=col, lw=2, marker=mk, ms=5, label=name, zorder=3)
        ax.scatter(d.host, d.v, color=col, s=10, alpha=0.5, edgecolor="none", zorder=2)
        rows += [(metric, name, h, v) for h, v in m.items()]
    host_axis(ax); ax.set_ylabel(ylab)
    return rows


fig = plt.figure(figsize=(7.2, 6.8))
gs = fig.add_gridspec(2, 3, hspace=0.75, wspace=0.6, left=0.09, right=0.97, top=0.9, bottom=0.13)
out = []

# a human precision
ax = fig.add_subplot(gs[0, 0]); tag(ax, "a")
out += lines(ax, "human_precision", scale=100, ylab="precision of host calls (%)")
ax.set_ylim(0, 102); ax.set_title("Host calls that are host")
ax.legend(loc="center right")

# b microbial reads called host per million microbial reads
ax = fig.add_subplot(gs[0, 1]); tag(ax, "b")
for run, (name, col, mk) in TOOLS.items():
    c = get("microbial_called_human", run)
    t1 = get("reads_microbial_assembled", run).set_index("lib").v
    t2 = get("reads_microbial_unassembled", run).set_index("lib").v
    c["v"] = c.v / (t1 + t2).reindex(c.lib).values * 1e6
    m = c.groupby("host").v.mean()
    ax.plot(m.index, m.values, color=col, lw=2, marker=mk, ms=5, label=name, zorder=3)
    ax.scatter(c.host, c.v, color=col, s=10, alpha=0.5, edgecolor="none", zorder=2)
    out += [("microbial_called_host_per_M", name, h, v) for h, v in m.items()]
ax.set_yscale("log"); host_axis(ax); ax.set_ylim(10, 1e6)
ax.set_ylabel("microbial reads called host\nper million microbial reads"); ax.set_title("Microbial reads lost to host")

# c host recall by RNA class (all lengths, pooled over libraries)
ax = fig.add_subplot(gs[0, 2]); tag(ax, "c")
cls = ["miRNA", "tRNA", "rRNA", "other"]
w = 0.38
for i, (run, (name, col, mk)) in enumerate(TOOLS.items()):
    x = e[(e.run == run) & e.metric.str.startswith("human_recall_")].copy()
    x["cls"] = x.metric.str.extract(r"human_recall_([A-Za-z]+)_")[0]
    mm = x.value.str.extract(r"\((\d+) of (\d+) missed\)").astype(float)
    x["miss"], x["tot"] = mm[0], mm[1]
    g = x.groupby("cls")[["miss", "tot"]].sum()
    r = (1 - g.miss / g.tot).reindex(cls) * 100
    ax.bar(np.arange(len(cls)) + (i - 0.5) * w, r.values, width=w - 0.04, color=col, label=name)
    out += [("host_recall_" + k, name, np.nan, v) for k, v in r.items()]
ax.set_xticks(range(len(cls)), ["miRNA", "tRNA", "rRNA", "other\nncRNA"]); ax.set_ylim(0, 100)
ax.set_ylabel("host reads recovered (%)"); ax.set_title("Host recall by RNA class")
ax.yaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True); ax.tick_params(axis="x", length=0)
ax.text(0.0, -0.24, "sRNA-IMP calls no host reads < 20 nt\nexcept mature miRNAs (≥ 17 nt)", transform=ax.transAxes,
        fontsize=6, color=INK2, va="top")

# d miRNA Spearman
ax = fig.add_subplot(gs[1, 0]); tag(ax, "d")
out += lines(ax, "mirna_spearman", ylab="Spearman ρ, miRNA counts vs truth")
ax.set_ylim(0.5, 1); ax.set_title("miRNA quantification")

# e miRNAs detected
ax = fig.add_subplot(gs[1, 1]); tag(ax, "e")
frac = lambda s: (lambda m: float(m.group(1)) / float(m.group(2)))(re.match(r"(\d+) of (\d+)", s))
out += lines(ax, "mirnas_true_detected", f=frac, scale=100, ylab="simulated miRNAs detected (%)")
ax.set_ylim(70, 100); ax.set_title("miRNA detection")

# f novel loci: hidden loci with >= 10 reads recovered, by class and length (sRNA-IMP, pooled)
ax = fig.add_subplot(gs[1, 2]); tag(ax, "f")
x = e[(e.run == a.run) & e.metric.str.startswith("recovered_ge10_")].copy()
p = x.metric.str.extract(r"recovered_ge10_([A-Za-z]+)_(.+)$")
x["cls"], x["len"] = p[0], p[1]
mm = x.value.str.extract(r"(\d+) of (\d+)").astype(float)
x["rec"], x["tot"] = mm[0], mm[1]
g = x.groupby(["cls", "len"])[["rec", "tot"]].sum()
lens = ["0-100nt", "100-250nt", ">=250nt"]
for i, (c, col) in enumerate([("ncRNA", CAT[2]), ("regulatory", CAT[6])]):
    r = [(g.loc[(c, L), "rec"] / g.loc[(c, L), "tot"] * 100) if (c, L) in g.index else np.nan for L in lens]
    ax.bar(np.arange(3) + (i - 0.5) * w, r, width=w - 0.04, color=col, label=c.replace("regulatory", "regulatory element"))
    out += [("hidden_recovered_%s_%s" % (c, L), "sRNA-IMP", np.nan, v) for L, v in zip(lens, r)]
ax.set_xticks(range(3), ["< 100", "100–250", "≥ 250"]); ax.set_xlabel("locus length (nt)"); ax.set_ylim(0, 100)
ax.set_ylabel("hidden loci recovered (%)"); ax.set_title("Unannotated loci")
ax.yaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True); ax.tick_params(axis="x", length=0)
ax.set_ylim(0, 125); ax.set_yticks([0, 25, 50, 75, 100]); ax.legend(loc="upper left", fontsize=6.5, ncol=1)
nl = e[(e.run == a.run) & (e.metric == "hidden_loci_with_ge10_reads_recovered")].value.str.extract(r"(\d+) of (\d+)").astype(float).sum()
ax.text(0.0, -0.3, "sRNA-IMP only; %.0f %% of hidden\nloci with ≥ 10 reads found" % (
    nl[0] / nl[1] * 100), transform=ax.transAxes, fontsize=6, color=INK2, va="top")

nlib = e.lib.nunique()
fig.suptitle("Simulation benchmark: %d stool-like libraries, sRNA-IMP vs sMETASeq" % nlib,
             x=0.01, ha="left", fontsize=9.5, fontweight="bold")
for ext in ("png", "pdf"):
    fig.savefig("%s.%s" % (a.out, ext), dpi=300)
pd.DataFrame(out, columns=["metric", "tool", "host_pct", "value"]).to_csv(a.out + ".tsv", sep="\t", index=False)
print("written", a.out)
