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
import shutil
import subprocess
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


def sam_mapped(path):
    """read names with at least one alignment in a SAM file (sMETASeq: host = mapped to GRCh38)"""
    ids = set()
    if not os.path.exists(path):
        return None
    with opener(path) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            f = line.split("\t", 2)
            if not int(f[1]) & 4:
                ids.add(f[0])
    return ids


def htseq_counts(path):
    """htseq-count output (name<TAB>count, '__' summary lines skipped)"""
    out = collections.Counter()
    if os.path.exists(path):
        with open(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) == 2 and not f[0].startswith("__"):
                    out[f[0]] += float(f[1])
    return out


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
    ap.add_argument("--smetaseq", default="", help="score an sMETASeq run instead (smetaseq_reimpl.sh -o PREFIX): "
                    "trimmed reads, host = reads mapped to GRCh38, miRNA counts from PREFIX.miRBase.csv")
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
    SM = a.smetaseq
    trimmed = fastq_ids(SM + ".cut.fastq" if SM else os.path.join(pre, "%s.trimmed.fastq.gz" % S))
    if trimmed is not None:
        dim = [r for r, t in truth.items() if t["origin"] == "dimer"]
        M.append(("qc", "dimers_removed", "%.4f" % (1 - sum(r in trimmed for r in dim) / max(len(dim), 1))))
        ok = [r for r, t in truth.items() if t["origin"] != "dimer" and int(t["insert_len"]) >= a.min_len]
        M.append(("qc", "inserts_ge_minlen_kept", "%.4f" % (sum(r in trimmed for r in ok) / max(len(ok), 1))))
        M.append(("qc", "reads_kept", len(trimmed)))

    # ---- split: human vs non-human
    hum = sam_mapped(SM + ".mapped.cut.sam") if SM else fastq_ids(os.path.join(pre, "%s.human.fastq.gz" % S))
    if hum is not None and trimmed is not None:
        tp = fp = fn = 0
        fp_cls = collections.Counter()
        h_all, h_miss = collections.Counter(), collections.Counter()   # human reads by (class, length bin)
        mir_in_host = collections.Counter()                            # simulated miRNA reads that reached the host fraction
        n_micro = collections.Counter()
        for r in trimmed:
            t = truth.get(r)
            if not t or t["origin"] == "dimer":
                continue
            is_h, call_h = t["origin"] == "human", r in hum
            if not is_h:
                n_micro[t["class"]] += 1
            if is_h and t["class"] == "miRNA" and call_h:
                mir_in_host[t["feature"]] += 1
            if is_h:
                L = int(t["insert_len"])
                lb = "<20" if L < 20 else "20-24" if L <= 24 else "25-34" if L <= 34 else ">=35"
                h_all[(t["class"], lb)] += 1
                if not call_h:
                    h_miss[(t["class"], lb)] += 1
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
        for k in sorted(h_all, key=lambda k: (k[0], ["<20", "20-24", "25-34", ">=35"].index(k[1]))):
            M.append(("split", "human_recall_%s_%snt" % k, "%.3f (%d of %d missed)" % (
                1 - h_miss[k] / h_all[k], h_miss[k], h_all[k])))

    # ---- human miRNAs
    if SM:
        called, resc = htseq_counts(SM + ".miRBase.csv"), collections.Counter()
    else:
        called = read_counts(os.path.join(ana, "human", "mirna", "%s.human.mirna.counts.tsv" % S))
        resc = read_counts(os.path.join(ana, "nonhuman", "mirna_like", "%s.mirna_like.human.counts.tsv" % S))
    try:
        mir_in_host
    except NameError:
        mir_in_host = None
    if called or resc:
        tot = called + resc
        tm = collections.Counter(t["feature"] for t in truth.values() if t["origin"] == "human" and t["class"] == "miRNA")
        names = sorted(set(tm) | set(tot))
        sp = spearman([tm[n] for n in names], [tot[n] for n in names])
        M.append(("human", "mirna_truth_reads", sum(tm.values())))
        M.append(("human", "mirna_called_reads", "%.1f" % sum(tot.values())))
        M.append(("human", "mirna_called_reads_rescued", "%.1f" % sum(resc.values())))
        if mir_in_host is not None:
            M.append(("human", "mirna_truth_reads_in_host_fraction", "%d (%.3f)" % (
                sum(mir_in_host.values()), sum(mir_in_host.values()) / max(sum(tm.values()), 1))))
            lost = [n for n in tm if tm[n] >= 10 and mir_in_host[n] == 0]
            M.append(("human", "mirnas_ge10_reads_none_in_host", "%d (%s)" % (len(lost), ",".join(sorted(lost)[:12]))))
        M.append(("human", "mirna_spearman", "%.3f" % sp))
        M.append(("human", "mirnas_true_detected", "%d of %d" % (sum(1 for n in tm if tot[n] > 0), len(tm))))
        false = [n for n in tot if tot[n] > 0 and tm[n] == 0]
        M.append(("human", "mirnas_false", "%d (%.1f reads)" % (len(false), sum(tot[n] for n in false))))
        with open(a.out + ".mirna.tsv", "w") as fh:
            fh.write("mirna\ttruth\tcalled\trescued\ttruth_in_host_fraction\n")
            for n in names:
                fh.write("%s\t%d\t%.2f\t%.2f\t%s\n" % (n, tm[n], tot[n], resc[n],
                                                      mir_in_host[n] if mir_in_host is not None else "-"))

    # ---- assembly: feature counts vs truth
    fc = os.path.join(ana, "assembly", "%s.feature_counts.tsv" % S)
    if os.path.exists(fc) and not SM:
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
    if a.hidden and os.path.exists(a.hidden) and not SM:
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
            for cls in sorted({r[5] for r in supported}):
                for lo, hi in ((0, 100), (100, 250), (250, 10 ** 9)):
                    sub = [r for r in supported if r[5] == cls and lo <= r[3] - r[2] < hi]
                    if sub:
                        M.append(("novel", "recovered_ge10_%s_%s" % (cls, "%d-%dnt" % (lo, hi) if hi < 10 ** 9 else
                                                                     ">=%dnt" % lo),
                                  "%d of %d" % (sum(1 for r in sub if r[7]), len(sub))))
            M.append(("novel", "candidates", len(cands)))
            M.append(("novel", "hidden_loci", len(hid)))
            M.append(("novel", "hidden_loci_recovered", "%d (%.3f)" % (found, found / max(len(hid), 1))))
            M.append(("novel", "hidden_loci_with_ge10_reads_recovered", "%d of %d" % (
                sum(1 for r in supported if r[7]), len(supported))))
            on_hidden = sum(1 for c in cands if any(h[1] == c[0] and h[2] < c[2] and c[1] < h[3] for h in hid))
            M.append(("novel", "candidates_on_hidden_loci", on_hidden))
            # what the other candidates are: the dominant truth class of the reads on them (same strand)
            hid_by = collections.defaultdict(list)
            for h in hid:
                hid_by[h[1]].append(h)
            need = {c[0] for c in cands}
            rd = collections.defaultdict(list)
            for t in truth.values():
                if t["origin"] == "microbial_assembled" and t["contig"] in need:
                    rd[t["contig"]].append((int(t["start"]), int(t["end"]), t["strand"], t["class"]))
            # Rfam hits of the candidates (cmscan --tblout; query = contig:start-end(strand) from bedtools getfasta)
            rfam = {}
            tbl = os.path.join(ana, "novel", "%s.candidates.rfam.tbl" % S)
            if os.path.exists(tbl):
                with open(tbl) as fh:
                    for line in fh:
                        if line.startswith("#") or not line.strip():
                            continue
                        f = line.split()
                        rfam.setdefault(f[2], f[0])          # best (first) hit per candidate
            def rfam_kind(name):
                n = name.lower()
                if "rrna" in n:
                    return "rRNA"
                if n.startswith("trna") or n == "tmrna" or n.startswith("tmrna"):
                    return "tRNA/tmRNA"
                return "other_Rfam"
            dom, lens, nohit = collections.Counter(), collections.Counter(), collections.Counter()
            notruth = []
            for c in cands:
                if any(h[2] < c[2] and c[1] < h[3] for h in hid_by[c[0]]):
                    continue
                k = collections.Counter(r[3] for r in rd[c[0]] if r[0] < c[2] and c[1] < r[1] and r[2] == c[3])
                d = k.most_common(1)[0][0] if k else "no_truth_reads"
                dom[d] += 1
                if d == "no_truth_reads":
                    notruth.append(c)
                if rfam or os.path.exists(tbl):
                    hit = rfam.get("%s:%d-%d(%s)" % (c[0], c[1], c[2], c[3]))
                    nohit["%s|%s" % (d, rfam_kind(hit) if hit else "no_Rfam_hit")] += 1
                L = c[2] - c[1]
                lens["<100" if L < 100 else "100-250" if L < 250 else ">=250"] += 1
            for k, n in dom.most_common():
                M.append(("novel", "other_candidates_mostly_%s" % k, n))
            for k, n in sorted(nohit.items(), key=lambda kv: -kv[1]):
                d, r = k.split("|")
                M.append(("novel", "other_candidates_%s_rfam_%s" % (d, r), n))
            for k in ("<100", "100-250", ">=250"):
                if lens[k]:
                    M.append(("novel", "other_candidates_len_%snt" % k, lens[k]))
            # where the reads on candidates without simulated reads come from: read names in the pipeline BAM
            bams = [x for x in (os.path.join(ana, "novel", "%s.sorted.bam" % S), os.path.join(ana, "assembly", "%s.bam" % S),
                                os.path.join(ana, "%s.bam" % S)) if os.path.exists(x)] + \
                sorted(glob.glob(os.path.join(ana, "*", "%s*.bam" % S)))
            bam = bams[0] if bams else ""
            if notruth and not bam:
                M.append(("novel", "no_truth_candidates_reads_from", "skipped: no BAM under %s" % ana))
            elif notruth and not shutil.which("samtools"):
                M.append(("novel", "no_truth_candidates_reads_from", "skipped: samtools not in PATH"))
            elif notruth and not os.path.exists(bam + ".bai") and not os.path.exists(bam[:-4] + ".bai"):
                subprocess.run(["samtools", "index", bam], check=False)
            if notruth and bam and shutil.which("samtools"):
                M.append(("novel", "no_truth_candidates_bam", bam))
                bed = a.out + ".notruth.bed"
                with open(bed, "w") as fh:
                    for c in notruth:
                        fh.write("%s\t%d\t%d\n" % (c[0], c[1], c[2]))
                p = subprocess.run(["samtools", "view", "-F", "4", "-L", bed, bam], capture_output=True, text=True)
                byc = collections.defaultdict(list)
                for line in p.stdout.splitlines():
                    f = line.split("\t", 10)
                    pos, L = int(f[3]) - 1, len(f[9])
                    byc[f[2]].append((pos, pos + L, f[0]))
                origin = collections.Counter()
                for c in notruth:
                    o = collections.Counter()
                    for s0, e0, rid in byc[c[0]]:
                        if s0 < c[2] and c[1] < e0:
                            t = truth.get(rid)
                            if not t:
                                o["unknown_read"] += 1
                            elif t["origin"] == "microbial_assembled":
                                o["assembled_elsewhere_" + t["class"]] += 1
                            else:
                                o["%s_%s" % (t["origin"], t["class"])] += 1
                    origin[o.most_common(1)[0][0] if o else "no_reads_in_bam"] += 1
                for k, n in origin.most_common(12):
                    M.append(("novel", "no_truth_candidates_reads_from_%s" % k, n))
                os.remove(bed)
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
