#!/usr/bin/env python3
"""srna_marker_check.py — follow one taxon / one gene marker across sRNA-Ext, sRNA-Int, metaT and metaG and relate
it to a sample variable (e.g. stool consistency).

Two marker types per sample, from the per-sample comparison outputs (--compare) and the annotation:
  taxon   share of all reads (CPM) on bins whose GTDB classification matches --taxon (regex), per source
          (<s>/Analysis/comparison/<s>.all_sources.bins.tsv)
  genes   CPM on CDS whose kegg_ko is in --ko (e.g. mcrA/B/G), per source
          (<s>/Analysis/comparison/<s>.all_sources.features.tsv + <s>/Analysis/annotation/<s>.gff)

Statistics: Spearman correlation of each marker with --var over all samples and within each group, BH-FDR; the
correlation between sources (does the sRNA signal follow metaG / metaT?).

Example (methanogens vs Bristol stool scale):
  srna_marker_check.py --outdir sRNA_v2 --metadata metadata.tsv --var bristol \\
      --taxon 'o__Methanobacteriales|o__Methanomassiliicoccales' --ko K00399,K00401,K00402 --out sRNA_v2/Cohort/marker_methanogens
"""
import argparse
import collections
import csv
import glob
import math
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import srna_cohort as SC  # noqa: E402  (spearman, bh, read_metadata)

SOURCES = ["sRNA-Ext", "sRNA-Int", "metaT", "metaG"]


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def gff_ko_ids(path, kos):
    ids = set()
    if not os.path.exists(path):
        return ids
    pat = re.compile(r"kegg_ko=([^;\t]+)")
    with open(path, errors="replace") as fh:
        for ln in fh:
            if "kegg_ko=" not in ln:
                continue
            m = pat.search(ln)
            if m and kos & set(m.group(1).split(",")):
                i = re.search(r"ID=([^;\t]+)", ln)
                if i:
                    ids.add(i.group(1))
    return ids


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True, help="pipeline output directory (one folder per sample)")
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--id-col", default="")
    ap.add_argument("--group-col", default="group")
    ap.add_argument("--var", required=True, help="numeric metadata column(s), comma-separated (e.g. bristol)")
    ap.add_argument("--taxon", default="", help="regex on the GTDB classification of the bins")
    ap.add_argument("--ko", default="", help="comma-separated KEGG KOs (gene marker)")
    ap.add_argument("--out", required=True, help="output prefix")
    a = ap.parse_args()

    hdr, meta = SC.read_metadata(a.metadata)
    idc = a.id_col or hdr[0]
    md = {r[idc].split(".")[0]: r for r in meta}
    kos = {k.strip() for k in a.ko.split(",") if k.strip()}
    tax = re.compile(a.taxon) if a.taxon else None
    rows = []
    for f in sorted(glob.glob(os.path.join(a.outdir, "*", "Analysis", "comparison", "*.all_sources.bins.tsv"))):
        s = os.path.basename(f).split(".all_sources")[0]
        if s not in md:
            continue
        r = {"sample": s, "group": md[s].get(a.group_col, "NA")}
        for v in a.var.split(","):
            r[v] = num(md[s].get(v))
        if tax:
            b = read_tsv(f)
            for src in SOURCES:
                col = src + "_CPM"
                if b and col in b[0]:
                    r["taxon_" + src] = sum(num(x[col]) for x in b if tax.search(x.get("classification", "")) and
                                            not math.isnan(num(x[col])))
        if kos:
            ff = f.replace(".all_sources.bins.tsv", ".all_sources.features.tsv")
            gff = os.path.join(a.outdir, s, "Analysis", "annotation", s + ".gff")
            ids = gff_ko_ids(gff, kos)
            r["genes_n"] = len(ids)
            if os.path.exists(ff) and ids:
                with open(ff) as fh:
                    rd = csv.DictReader(fh, delimiter="\t")
                    tot = collections.defaultdict(float)
                    for x in rd:
                        if x.get("feature") in ids:
                            for src in SOURCES:
                                v = num(x.get(src + "_CPM"))
                                if not math.isnan(v):
                                    tot[src] += v
                for src in SOURCES:
                    if (src + "_CPM") in rd.fieldnames:
                        r["genes_" + src] = tot[src]
        rows.append(r)
        if len(rows) % 25 == 0:
            sys.stderr.write("[marker] %d samples read\n" % len(rows))
    if not rows:
        sys.exit("no comparison outputs (*.all_sources.bins.tsv) under %s" % a.outdir)
    cols = ["sample", "group"] + a.var.split(",") + [k for k in ("genes_n",) if any(k in r for r in rows)] + \
        ["%s_%s" % (m, s) for m in ("taxon", "genes") for s in SOURCES if any("%s_%s" % (m, s) in r for r in rows)]
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out + ".samples.tsv", "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in rows:
            o.write("\t".join(("%.4g" % r[c]) if isinstance(r.get(c), float) else str(r.get(c, "")) for c in cols) + "\n")
    # correlations with the variable(s), overall and per group
    groups = sorted({r["group"] for r in rows})
    stats = []
    markers = [c for c in cols if c.startswith(("taxon_", "genes_")) and c != "genes_n"]
    for v in a.var.split(","):
        for m in markers:
            for g in ["all"] + groups:
                sel = [r for r in rows if (g == "all" or r["group"] == g) and not math.isnan(r.get(v, float("nan")))
                       and isinstance(r.get(m), float)]
                if len(sel) < 8:
                    continue
                x = np.array([r[v] for r in sel])
                y = np.array([r[m] for r in sel])
                rho, p, _n = SC.spearman(x, y)
                stats.append([v, m, g, len(sel), int((y > 0).sum()), rho, p])
    q = SC.bh([s_[6] for s_ in stats])
    # agreement between sources
    agree = []
    for mk in ("taxon", "genes"):
        srcs = [s for s in SOURCES if "%s_%s" % (mk, s) in cols]
        for i, s1 in enumerate(srcs):
            for s2 in srcs[i + 1:]:
                sel = [r for r in rows if isinstance(r.get("%s_%s" % (mk, s1)), float) and
                       isinstance(r.get("%s_%s" % (mk, s2)), float)]
                if len(sel) >= 8:
                    rho, p, _n = SC.spearman(np.array([r["%s_%s" % (mk, s1)] for r in sel]),
                                         np.array([r["%s_%s" % (mk, s2)] for r in sel]))
                    agree.append([mk, s1, s2, len(sel), rho, p])
    with open(a.out + ".stats.tsv", "w") as o:
        o.write("variable\tmarker\tgroup\tn\tdetected\trho\tp\tq\n")
        for s_, qq in zip(stats, q):
            o.write("\t".join(str(x) for x in s_[:5]) + "\t%.3f\t%.3g\t%.3g\n" % (s_[5], s_[6], qq))
        o.write("\n# agreement between sources (Spearman over samples)\nmarker\tsource_a\tsource_b\tn\trho\tp\n")
        for r in agree:
            o.write("\t".join(str(x) for x in r[:4]) + "\t%.3f\t%.3g\n" % (r[4], r[5]))
    # figure: marker (log10 CPM + 1) vs variable, one panel per marker and source
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        v = a.var.split(",")[0]
        ncol = len(SOURCES)
        mks = [mk for mk in ("taxon", "genes") if any(c.startswith(mk + "_") for c in markers)]
        fig, axes = plt.subplots(len(mks), ncol, figsize=(3.2 * ncol, 2.9 * len(mks)), squeeze=False)
        colors = {g: c for g, c in zip(groups, ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9"])}
        rng = np.random.default_rng(1)
        st = {(s_[0], s_[1], s_[2]): (s_[5], qq) for s_, qq in zip(stats, q)}
        for i, mk in enumerate(mks):
            for j, src in enumerate(SOURCES):
                ax = axes[i][j]
                m = "%s_%s" % (mk, src)
                for g in groups:
                    sel = [r for r in rows if r["group"] == g and isinstance(r.get(m), float) and not math.isnan(r[v])]
                    ax.scatter([r[v] + rng.uniform(-0.15, 0.15) for r in sel], [math.log10(r[m] + 1) for r in sel],
                               s=10, alpha=0.7, color=colors[g], label=g, lw=0)
                rho_q = st.get((v, m, "all"))
                ax.set_title("%s · %s%s" % ({"taxon": "taxon bins", "genes": "marker genes"}[mk], src,
                                            "\nρ %.2f, q %.2g" % rho_q if rho_q else ""), fontsize=8, loc="left")
                ax.set_xlabel(v, fontsize=7)
                ax.set_ylabel("log10(CPM + 1)", fontsize=7)
                ax.tick_params(labelsize=6)
        axes[0][0].legend(frameon=False, fontsize=6)
        fig.suptitle("%s%s vs %s" % ("taxon: %s" % a.taxon if a.taxon else "", "  genes: %s" % a.ko if a.ko else "", v),
                     fontsize=9, x=0.01, ha="left")
        fig.tight_layout()
        fig.savefig(a.out + ".pdf")
        fig.savefig(a.out + ".png", dpi=130)
    except Exception as e:                                      # the tables are the main output
        sys.stderr.write("[marker] no figure (%s)\n" % e)
    sys.stderr.write("[marker] %d samples -> %s.samples.tsv, %s.stats.tsv, %s.pdf\n" % (len(rows), a.out, a.out, a.out))


if __name__ == "__main__":
    main()
