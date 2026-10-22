#!/usr/bin/env python3
"""
srna_integrity.py — RNA integrity / degradation indicators of an sRNA library

Classic 5'->3' coverage skew (RIN-like gene-body metrics) assumes reads that
tile long transcripts. sRNA libraries are size-selected, so instead this
script looks at how *defined* the fragment ends on abundant RNAs are:
processed RNAs (mature tRNAs, tRNA halves, crRNAs) have sharp ends, random
(chemical / nuclease) degradation produces diffuse ends and an even spread
of fragments along long RNAs.

Inputs are outputs of the fragments and assembly modules:
  <frag>.fragments.features.tsv, .classes.tsv, [.ends.tsv]
  <asm>.length_dist.tsv, <asm>.mapping_summary.tsv   (optional)

Outputs (<out-prefix>.*):
  integrity.tsv           metric, value, features, direction, description
  integrity.profile.tsv   rRNA fragment midpoints in 10 bins along complete
                          16S / 23S genes (fraction of reads, strand-aware)
  integrity.features.tsv  per feature with enough reads: 5'/3'-end entropy
                          and share of the top-3 end positions

Indicators are relative: compare samples of the same fraction (cell-free
sRNA-Ext is expected to be more fragmented than sRNA-Int for biological
reasons) and, if available, against RIN/DV200 of the same extractions.
"""
import argparse
import collections
import math
import os
import sys

TRNA_MATURE_END = ["full_length", "5p_half", "3p_half", "5p_tRF", "3p_tRF", "5p_long", "3p_long"]


def read_tsv(path):
    if not path or not os.path.exists(path):
        return None
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(hdr, l.rstrip("\n").split("\t"))) for l in fh if l.strip()]


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return 0.0


def complete_rrna(product):
    p = product.lower()
    return "partial" not in p and "truncated" not in p


def rrna_kind(product):
    p = product.lower()
    for k in ("16s", "23s", "5s"):
        if k in p:
            return k.upper()
    return "other"


def ent_top(counts, L):
    """normalised Shannon entropy of end positions and share of the top-3 positions"""
    tot = sum(counts)
    if tot <= 0:
        return float("nan"), float("nan")
    h = -sum((c / tot) * math.log(c / tot) for c in counts if c > 0)
    hmax = math.log(max(2.0, min(float(L), tot)))      # few reads cannot reach log(L)
    top3 = sum(sorted(counts, reverse=True)[:3]) / tot
    return max(0.0, min(1.0, h / hmax)), top3


def wmean(pairs):
    s = sum(w for v, w in pairs if not math.isnan(v))
    return (sum(v * w for v, w in pairs if not math.isnan(v)) / s) if s else float("nan")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--frag-prefix", required=True, help="prefix of the fragments module outputs")
    ap.add_argument("--asm-prefix", default="", help="prefix of the assembly module outputs")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--sample", default="")
    ap.add_argument("--min-trna-reads", type=float, default=50)
    ap.add_argument("--min-rrna-reads", type=float, default=200)
    a = ap.parse_args()

    feats = read_tsv(a.frag_prefix + ".fragments.features.tsv")
    if feats is None:
        sys.exit("missing %s.fragments.features.tsv" % a.frag_prefix)
    sample = a.sample or os.path.basename(a.frag_prefix)
    M = []   # (metric, value, n_features, direction, description)

    # ---------------- length profile / mapping ----------------
    ld = read_tsv(a.asm_prefix + ".length_dist.tsv") if a.asm_prefix else None
    if ld:
        tot = collections.Counter()
        for r in ld:
            tot[int(r["length"])] += fnum(r.get("mapped_weighted")) + fnum(r.get("unmapped"))
        n = sum(tot.values())
        if n:
            short = sum(v for l, v in tot.items() if l <= 20) / n
            long_ = sum(v for l, v in tot.items() if l >= 30) / n
            acc, med = 0.0, 0
            for l in sorted(tot):
                acc += tot[l]
                if acc >= n / 2:
                    med = l
                    break
            M += [("reads_le20nt_fraction", short, "", "higher = more fragmented",
                   "fraction of reads of <=20 nt (all reads entering the assembly mapping)"),
                  ("reads_ge30nt_fraction", long_, "", "higher = longer fragments",
                   "fraction of reads of >=30 nt"),
                  ("median_read_length", med, "", "", "median read length (nt)")]
    ms = read_tsv(a.asm_prefix + ".mapping_summary.tsv") if a.asm_prefix else None
    if ms:
        d = {r["category"]: fnum(r["value"]) for r in ms}
        if d.get("input_reads"):
            M.append(("assembly_mapped_fraction", d.get("mapped_reads", 0) / d["input_reads"], "",
                      "", "fraction of reads mapping to the sample's contigs"))

    # ---------------- tRNA ends ----------------
    trna = [r for r in feats if r["type"] == "tRNA"]
    t_tot = sum(fnum(r["total"]) for r in trna)
    if t_tot:
        mature = sum(fnum(r["tRNA:" + c]) for r in trna for c in TRNA_MATURE_END)
        halves = sum(fnum(r["tRNA:5p_half"]) + fnum(r["tRNA:3p_half"]) for r in trna)
        cca_w = sum(fnum(r["tRNA_3p_with_CCA"]) for r in trna)
        cca_wo = sum(fnum(r["tRNA_3p_without_CCA"]) for r in trna)
        M += [("tRNA_mature_end_fraction", mature / t_tot, len(trna), "higher = more intact",
               "tRNA reads with at least one end at the mature 5' or 3' end"),
              ("tRNA_internal_fraction", 1 - mature / t_tot, len(trna), "higher = more degraded",
               "tRNA reads with both ends inside the tRNA (i-tRFs)"),
              ("tRNA_halves_fraction", halves / t_tot, len(trna), "biology (stress), not degradation",
               "5'/3' tRNA halves (anticodon cleavage)")]
        if cca_w + cca_wo:
            M.append(("tRNA_3p_CCA_fraction", cca_w / (cca_w + cca_wo), len(trna), "higher = more intact",
                      "3'-anchored tRNA reads that carry the terminal CCA"))

    # ---------------- rRNA spread along complete genes ----------------
    prof = collections.defaultdict(lambda: [0.0] * 10)
    even = []
    for r in feats:
        if r["type"] != "rRNA" or not complete_rrna(r["product"]):
            continue
        k = rrna_kind(r["product"])
        if k not in ("16S", "23S"):
            continue
        b = [fnum(r["midpoint_bin%d" % i]) for i in range(1, 11)]
        s = sum(b)
        if s < a.min_rrna_reads:
            continue
        for i in range(10):
            prof[k][i] += b[i]
        h = -sum((x / s) * math.log(x / s) for x in b if x > 0) / math.log(10)
        even.append((h, s))
    if even:
        M.append(("rRNA_bin_evenness", wmean(even), len(even), "higher = more random fragmentation",
                  "evenness (0-1) of fragment midpoints over 10 bins along complete 16S/23S genes"))
    for k, b in sorted(prof.items()):
        s = sum(b)
        r53 = (sum(b[:3]) + 1) / (sum(b[7:]) + 1)
        M.append(("rRNA_%s_5p_3p_log2ratio" % k, math.log2(r53), "", "",
                  "log2 of reads in the 5' vs the 3' 30%% of complete %s genes" % k))

    # ---------------- end entropy (per-position table) ----------------
    fe_rows = []
    ends = a.frag_prefix + ".fragments.ends.tsv"
    if os.path.exists(ends):
        meta = {r["feature_id"]: r for r in feats}
        cur, e5, e3, L, typ = None, [], [], 0, ""

        def flush():
            if cur is None:
                return
            n = sum(e5)
            r = meta.get(cur, {})
            if typ == "tRNA" and n >= a.min_trna_reads or \
               typ == "rRNA" and n >= a.min_rrna_reads and complete_rrna(r.get("product", "")):
                h5, t5 = ent_top(e5, L)
                h3, t3 = ent_top(e3, L)
                fe_rows.append((cur, typ, r.get("product", ""), L, n, h5, t5, h3, t3))
        with open(ends) as fh:
            fh.readline()
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[0] != cur:
                    flush()
                    cur, e5, e3, L, typ = f[0], [], [], int(f[3]), f[1]
                e5.append(fnum(f[5]))
                e3.append(fnum(f[6]))
            flush()
        for t in ("tRNA", "rRNA"):
            rows = [x for x in fe_rows if x[1] == t]
            if not rows:
                continue
            M += [("%s_5p_end_entropy" % t, wmean([(x[5], x[4]) for x in rows]), len(rows),
                   "higher = more diffuse ends",
                   "normalised entropy (0-1) of 5' end positions, read-weighted mean over %ss" % t),
                  ("%s_5p_top3_share" % t, wmean([(x[6], x[4]) for x in rows]), len(rows),
                   "higher = sharper ends", "share of reads whose 5' end is at one of the 3 main positions"),
                  ("%s_3p_top3_share" % t, wmean([(x[8], x[4]) for x in rows]), len(rows),
                   "higher = sharper ends", "share of reads whose 3' end is at one of the 3 main positions")]
    else:
        sys.stderr.write("[integrity] no %s — end entropy skipped\n" % ends)

    # ---------------- write ----------------
    p = a.out_prefix
    os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)
    with open(p + ".integrity.tsv", "w") as out:
        out.write("sample\tmetric\tvalue\tfeatures\tdirection\tdescription\n")
        for m, v, n, d, desc in M:
            out.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (sample, m, ("%.4f" % v) if isinstance(v, float) else v,
                                                    n, d, desc))
    with open(p + ".integrity.profile.tsv", "w") as out:
        out.write("rRNA\tbin\tfraction\n")
        for k, b in sorted(prof.items()):
            s = sum(b) or 1
            for i, x in enumerate(b):
                out.write("%s\t%d\t%.4f\n" % (k, i + 1, x / s))
    with open(p + ".integrity.features.tsv", "w") as out:
        out.write("feature_id\ttype\tproduct\tlength\treads\tentropy_5p\ttop3_5p\tentropy_3p\ttop3_3p\n")
        for x in sorted(fe_rows, key=lambda x: -x[4]):
            out.write("%s\t%s\t%s\t%d\t%.1f\t%.4f\t%.4f\t%.4f\t%.4f\n" % x)
    sys.stderr.write("[integrity] %d metrics, %d features with end statistics\n" % (len(M), len(fe_rows)))


if __name__ == "__main__":
    main()
