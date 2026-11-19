#!/usr/bin/env python3
"""compare_reference_mode.py — one library run twice: on the sample's own assembly and on reference genomes
(srna_prepare_inputs.py --genomes, e.g. MGnify genomes from srna_fetch_mgnify.py). Reports what the reference set
recovers.

  python3 paper/compare_reference_mode.py OUT_ASSEMBLY OUT_REFERENCE SAMPLE TYPE OUTPREFIX [--labels own,reference]

Reads per run: Stats/<type>/<s>.<type>.read_stats.tsv (host split, contig rescue), assembly/<s>.mapping_summary.tsv
(mapped share), assembly/<s>.category_counts.tsv (sense categories), novel/<s>.assembly.summary.tsv (novel loci),
crispr/<s>.crispr_arrays.tsv (arrays with >= 10 reads) and bins/<s>.taxon_counts.tsv (reads per genus).
Writes OUTPREFIX.metrics.tsv (one row per metric, one column per run) and OUTPREFIX.genera.tsv (reads per genus in
both runs; GTDB suffixes such as _A are removed for matching, since the two taxonomies can come from different GTDB
releases), and prints the Spearman correlation of genus counts over genera seen in both.
"""
import argparse
import csv
import math
import os
import re


def rows(path):
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def kv(path, key, val):
    return {r[key]: r[val] for r in rows(path)}


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def rank(v):
    o = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(o):
        j = i
        while j + 1 < len(o) and v[o[j + 1]] == v[o[i]]:
            j += 1
        for k in range(i, j + 1):
            r[o[k]] = (i + j) / 2 + 1
        i = j + 1
    return r


def spearman(a, b):
    if len(a) < 3:
        return float("nan")
    ra, rb = rank(a), rank(b)
    ma, mb = sum(ra) / len(ra), sum(rb) / len(rb)
    num_ = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
    den = math.sqrt(sum((x - ma) ** 2 for x in ra) * sum((y - mb) ** 2 for y in rb))
    return num_ / den if den else float("nan")


def metrics(out, s, t):
    a = os.path.join(out, s, "Analysis", t)
    rs = kv(os.path.join(out, s, "Stats", t, "%s.%s.read_stats.tsv" % (s, t)), "step", "reads")
    ms = kv(os.path.join(a, "assembly", "%s.mapping_summary.tsv" % s), "category", "value")
    cat = {r["category"]: num(r["reads"]) for r in rows(os.path.join(a, "assembly", "%s.category_counts.tsv" % s))}
    nov = kv(os.path.join(a, "novel", "%s.assembly.summary.tsv" % s), "category", "reads")
    arr = rows(os.path.join(a, "crispr", "%s.crispr_arrays.tsv" % s))
    tax = [r for r in rows(os.path.join(a, "bins", "%s.taxon_counts.tsv" % s)) if r["rank"] == "g"]
    mapped = num(ms.get("mapped_reads"))
    m = {
        "trimmed reads": num(rs.get("trimmed")),
        "host reads": num(rs.get("human_total")),
        "host hits returned by contig rescue": num(rs.get("host_rescued_by_assembly")),
        "non-host reads": num(rs.get("nonhuman")),
        "non-host reads mapped": mapped,
        "share of non-host reads mapped": mapped / num(ms.get("input_reads")) if ms else float("nan"),
    }
    tot = sum(cat.values()) or float("nan")
    for c in ("tRNA", "rRNA", "CDS", "antisense_CDS", "ncRNA", "tmRNA", "CRISPR", "intergenic"):
        m["share of mapped reads: %s" % c] = cat.get(c, 0.0) / tot
    m["novel loci (intergenic)"] = num(nov.get("novel_intergenic_candidates"))
    m["novel loci (antisense to CDS)"] = num(nov.get("novel_antisense_candidates"))
    m["CRISPR arrays with >= 10 reads"] = sum(1 for r in arr if num(r.get("reads")) >= 10)
    gen = {}
    for r in tax:
        g = r["taxon"]
        if g in ("unbinned",) or g.endswith("unclassified"):
            m.setdefault("share of mapped reads not placed on a genus", 0.0)
            m["share of mapped reads not placed on a genus"] += num(r["reads"]) / mapped if mapped else 0
            continue
        k = re.sub(r"_[A-Z]+$", "", g.replace("g__", ""))
        gen[k] = gen.get(k, 0.0) + num(r["reads"])
    m["genera with reads"] = len(gen)
    return m, gen


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out_a"); ap.add_argument("out_b"); ap.add_argument("sample"); ap.add_argument("type")
    ap.add_argument("outprefix"); ap.add_argument("--labels", default="own_assembly,reference")
    a = ap.parse_args()
    la, lb = a.labels.split(",")
    ma, ga = metrics(a.out_a, a.sample, a.type)
    mb, gb = metrics(a.out_b, a.sample, a.type)
    with open(a.outprefix + ".metrics.tsv", "w") as fo:
        fo.write("metric\t%s\t%s\n" % (la, lb))
        for k in list(ma) + [k for k in mb if k not in ma]:
            fo.write("%s\t%s\t%s\n" % (k, ("%.4g" % ma.get(k, float("nan"))), ("%.4g" % mb.get(k, float("nan")))))
    both = sorted(set(ga) & set(gb))
    with open(a.outprefix + ".genera.tsv", "w") as fo:
        fo.write("genus\t%s\t%s\n" % (la, lb))
        for g in sorted(set(ga) | set(gb), key=lambda g: -(ga.get(g, 0) + gb.get(g, 0))):
            fo.write("%s\t%.1f\t%.1f\n" % (g, ga.get(g, 0), gb.get(g, 0)))
    rho = spearman([ga[g] for g in both], [gb[g] for g in both])
    sa, sb = sum(ga.values()), sum(gb.values())
    print("genera: %d %s, %d %s, %d in both; Spearman of genus counts over shared genera %.2f; shared genera carry "
          "%.0f %% of the %s and %.0f %% of the %s genus reads" % (
              len(ga), la, len(gb), lb, len(both), rho, 100 * sum(ga[g] for g in both) / sa if sa else float("nan"), la,
              100 * sum(gb[g] for g in both) / sb if sb else float("nan"), lb))
    for k in ("share of non-host reads mapped", "host hits returned by contig rescue", "novel loci (intergenic)",
              "CRISPR arrays with >= 10 reads"):
        print("%-40s %10.4g %10.4g" % (k, ma.get(k, float("nan")), mb.get(k, float("nan"))))


if __name__ == "__main__":
    main()
