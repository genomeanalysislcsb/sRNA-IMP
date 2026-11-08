#!/usr/bin/env python3
"""srna_sim_eval.py — score a pipeline run on a simulated library (srna_simulate.py) against the truth.

Usage
  srna_sim_eval.py --truth sim/SIM01.truth.tsv.gz --hidden sim/SIM01.hidden_loci.tsv \\
      --outdir <pipeline outdir> --sample SIM01 [--type sRNA] [--min-len 18] --out sim/SIM01.eval

Writes <out>.summary.tsv (metric, value) and <out>.mirna.tsv / <out>.features.tsv / <out>.novel.tsv, and prints the
summary. Modules that were not run are skipped.

Metrics
  qc        dimers removed; reads kept vs truth inserts >= --min-len
  split     host vs microbial assignment of the trimmed reads: precision / recall of "human", microbial reads called
            human per class (the conserved-sequence errors the contig rescue is meant to remove)
  human     miRNA counts (genome path + rescued miRNA-like reads) vs truth: Spearman, detected, false miRNAs
  assembly  feature counts (sense) vs truth reads per feature (assembled, not hidden features): Spearman
  novel     hidden ncRNA / regulatory features recovered by the novel candidates (same contig and strand, overlap)
"""
import argparse
import collections
import glob
import gzip
import os
import sys


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def fastq_ids(path):
    ids = set()
    if not os.path.exists(path):
        return None
    with opener(path) as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                ids.add(line[1:].split()[0])
    return ids


def ranks(x):
    order = sorted(range(len(x)), key=lambda i: x[i])
    r = [0.0] * len(x)
    i = 0
    while i < len(x):
        j = i
        while j + 1 < len(x) and x[order[j + 1]] == x[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = (i + j) / 2.0 + 1
        i = j + 1
    return r


def spearman(a, b):
    if len(a) < 3:
        return float("nan")
    ra, rb = ranks(a), ranks(b)
    ma, mb = sum(ra) / len(ra), sum(rb) / len(rb)
    num = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
    den = (sum((x - ma) ** 2 for x in ra) * sum((y - mb) ** 2 for y in rb)) ** 0.5
    return num / den if den else float("nan")


def read_counts(path, key="feature", val="count"):
    out = collections.Counter()
    for p in glob.glob(path):
        with opener(p) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            ki, vi = hdr.index(key), hdr.index(val)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                try:
                    out[f[ki]] += float(f[vi])
                except (ValueError, IndexError):
                    pass
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--hidden", default="")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--type", default="sRNA")
    ap.add_argument("--min-len", type=int, default=18, help="MINLEN of the run [18]")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    S, T = a.sample, a.type
    pre = os.path.join(a.outdir, S, "Preprocessing", T)
    ana = os.path.join(a.outdir, S, "Analysis", T)
    M = []                                               # (section, metric, value)

    truth = {}
    with opener(a.truth) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = dict(zip(hdr, line.rstrip("\n").split("\t")))
            truth[f["read_id"]] = f
    M.append(("truth", "reads", len(truth)))
    for o, n in sorted(collections.Counter(t["origin"] for t in truth.values()).items()):
        M.append(("truth", "reads_" + o, n))

    # ---- qc
    trimmed = fastq_ids(os.path.join(pre, "%s.trimmed.fastq.gz" % S))
    if trimmed is not None:
        dim = [r for r, t in truth.items() if t["origin"] == "dimer"]
        M.append(("qc", "dimers_removed", "%.4f" % (1 - sum(r in trimmed for r in dim) / max(len(dim), 1))))
        ok = [r for r, t in truth.items() if t["origin"] != "dimer" and int(t["insert_len"]) >= a.min_len]
        M.append(("qc", "inserts_ge_minlen_kept", "%.4f" % (sum(r in trimmed for r in ok) / max(len(ok), 1))))
        M.append(("qc", "reads_kept", len(trimmed)))

    # ---- split: human vs non-human
    hum = fastq_ids(os.path.join(pre, "%s.human.fastq.gz" % S))
    if hum is not None and trimmed is not None:
        tp = fp = fn = 0
        fp_cls = collections.Counter()
        n_micro = collections.Counter()
        for r in trimmed:
            t = truth.get(r)
            if not t or t["origin"] == "dimer":
                continue
            is_h, call_h = t["origin"] == "human", r in hum
            if not is_h:
                n_micro[t["class"]] += 1
            if is_h and call_h:
                tp += 1
            elif call_h:
                fp += 1
                fp_cls[t["class"]] += 1
            elif is_h:
                fn += 1
        M.append(("split", "human_precision", "%.4f" % (tp / max(tp + fp, 1))))
        M.append(("split", "human_recall", "%.4f" % (tp / max(tp + fn, 1))))
        M.append(("split", "microbial_called_human", fp))
        for c, n in fp_cls.most_common():
            M.append(("split", "microbial_called_human_%s" % c, "%d (%.4f%% of %s)" % (n, 100.0 * n / max(n_micro[c], 1), c)))

    # ---- human miRNAs
    called = read_counts(os.path.join(ana, "human", "mirna", "%s.human.mirna.counts.tsv" % S))
    resc = read_counts(os.path.join(ana, "nonhuman", "mirna_like", "%s.mirna_like.human.counts.tsv" % S))
    if called or resc:
        tot = called + resc
        tm = collections.Counter(t["feature"] for t in truth.values() if t["origin"] == "human" and t["class"] == "miRNA")
        names = sorted(set(tm) | set(tot))
        sp = spearman([tm[n] for n in names], [tot[n] for n in names])
        M.append(("human", "mirna_truth_reads", sum(tm.values())))
        M.append(("human", "mirna_called_reads", "%.1f" % sum(tot.values())))
        M.append(("human", "mirna_called_reads_rescued", "%.1f" % sum(resc.values())))
        M.append(("human", "mirna_spearman", "%.3f" % sp))
        M.append(("human", "mirnas_true_detected", "%d of %d" % (sum(1 for n in tm if tot[n] > 0), len(tm))))
        false = [n for n in tot if tot[n] > 0 and tm[n] == 0]
        M.append(("human", "mirnas_false", "%d (%.1f reads)" % (len(false), sum(tot[n] for n in false))))
        with open(a.out + ".mirna.tsv", "w") as fh:
            fh.write("mirna\ttruth\tcalled\trescued\n")
            for n in names:
                fh.write("%s\t%d\t%.2f\t%.2f\n" % (n, tm[n], tot[n], resc[n]))

    # ---- assembly: feature counts vs truth
    fc = os.path.join(ana, "assembly", "%s.feature_counts.tsv" % S)
    if os.path.exists(fc):
        sense = read_counts(fc, "feature_id", "sense_reads")
        anti = read_counts(fc, "feature_id", "antisense_reads")
        ts, ta = collections.Counter(), collections.Counter()
        for t in truth.values():
            if t["origin"] == "microbial_assembled" and t["hidden"] == "no" and t["feature"] != "-":
                (ta if t["class"] == "CDS_antisense" else ts)[t["feature"]] += 1
        names = sorted(set(ts) | set(sense))
        M.append(("assembly", "feature_sense_spearman", "%.3f" % spearman([ts[n] for n in names], [sense[n] for n in names])))
        M.append(("assembly", "feature_sense_reads_truth_vs_called", "%d / %.0f" % (sum(ts.values()), sum(sense.values()))))
        cds = [n for n in set(ta) | set(anti) if n in ta or anti[n] > 0]
        M.append(("assembly", "cds_antisense_reads_truth_vs_called", "%d / %.0f" % (sum(ta.values()), sum(anti[n] for n in cds))))
        with open(a.out + ".features.tsv", "w") as fh:
            fh.write("feature\ttruth_sense\tcalled_sense\ttruth_antisense\tcalled_antisense\n")
            for n in sorted(set(names) | set(cds)):
                fh.write("%s\t%d\t%.2f\t%d\t%.2f\n" % (n, ts[n], sense[n], ta[n], anti[n]))

    # ---- novel: hidden loci recovered
    if a.hidden and os.path.exists(a.hidden):
        cands = []
        for p in glob.glob(os.path.join(ana, "novel", "%s.candidates.*.bed" % S)):
            if p.endswith(".full.bed"):
                continue
            with open(p) as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) >= 6:
                        cands.append((f[0], int(f[1]), int(f[2]), f[5], os.path.basename(p).split(".")[-2]))
        if cands:
            hid = []
            with open(a.hidden) as fh:
                fh.readline()
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    hid.append((f[0], f[1], int(f[2]), int(f[3]), f[4], f[5]))
            reads = collections.Counter(t["feature"] for t in truth.values() if t["hidden"] == "yes")
            by = collections.defaultdict(list)
            for c in cands:
                by[c[0]].append(c)
            found, rows = 0, []
            for fid, ctg, s, e, st, cls in hid:
                hit = [c for c in by[ctg] if c[1] < e and s < c[2] and c[3] == st]
                found += bool(hit)
                rows.append((fid, ctg, s, e, st, cls, reads[fid], len(hit)))
            supported = [r for r in rows if r[6] >= 10]
            M.append(("novel", "candidates", len(cands)))
            M.append(("novel", "hidden_loci", len(hid)))
            M.append(("novel", "hidden_loci_recovered", "%d (%.3f)" % (found, found / max(len(hid), 1))))
            M.append(("novel", "hidden_loci_with_ge10_reads_recovered", "%d of %d" % (
                sum(1 for r in supported if r[7]), len(supported))))
            on_hidden = sum(1 for c in cands if any(h[1] == c[0] and h[2] < c[2] and c[1] < h[3] for h in hid))
            M.append(("novel", "candidates_on_hidden_loci", on_hidden))
            with open(a.out + ".novel.tsv", "w") as fh:
                fh.write("feature\tcontig\tstart\tend\tstrand\tclass\ttruth_reads\tcandidates_overlapping\n")
                for r in rows:
                    fh.write("\t".join(map(str, r)) + "\n")

    with open(a.out + ".summary.tsv", "w") as fh:
        fh.write("section\tmetric\tvalue\n")
        for s, m, v in M:
            fh.write("%s\t%s\t%s\n" % (s, m, v))
            sys.stdout.write("%-9s %-42s %s\n" % (s, m, v))


if __name__ == "__main__":
    main()
