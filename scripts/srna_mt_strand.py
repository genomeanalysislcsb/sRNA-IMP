#!/usr/bin/env python3
"""srna_mt_strand.py — strand-resolved metaT reads per feature, for the antisense small RNA reads (--compare).

The metaT feature tables of IMP count only the gene strand. With a stranded metaT library the metaT BAM shows, for
every feature that carries antisense small RNA reads, whether the gene strand, the opposite strand (= the strand of
the antisense sRNA) or both are transcribed. featureCounts is run twice on the metaT BAM (-s 1 and -s 2, paired and
unpaired reads, multi-overlapping and multi-mapping reads as fractions); which of the two is the gene strand is taken
from the data (--strandedness auto: per CDS share of the reads, median over 2,000 random CDS) or given. Only the alignments over
these features (with their mates) are fetched from the indexed BAM before counting, so a sample takes minutes; without
an index, or with --full, the whole BAM is counted.

Classes per feature with >= --min-antisense antisense sRNA reads (metaT >= --min-mt reads on a strand counts):
  sRNA strand only   metaT only on the antisense (sRNA) strand: the locus is transcribed in the direction of the
                     sRNA, not of the annotated gene (often a short hypothetical ORF over an unannotated transcript)
  both strands       sense and antisense transcription (antisense to an expressed gene)
  gene strand only   metaT only on the gene strand
  low metaT          fewer than --min-mt metaT reads on both strands

Usage: srna_mt_strand.py --sample S --bam mt.reads.sorted.bam --srna TYPE=<feature_counts.tsv> [--srna ...]
                         --out-prefix OUT/S [--strandedness auto|reverse|forward] [--threads N]
Feature tables: sRNA-IMP <S>.feature_counts.tsv (feature_id contig start end strand type sense_reads antisense_reads
product ...; start 0-based). Outputs: OUT/S.<TYPE>.metaT_strand.tsv (features with >= --min-antisense antisense
reads), OUT/S.metaT_strand.features.tsv (metaT per feature and strand, all features), OUT/S.metaT_strand.summary.tsv.
Unstranded metaT (gene-strand share below --min-stranded) is reported in the summary and gets no classes."""
import argparse
import os
import random
import subprocess
import sys
import tempfile

CLASSES = ["sRNA strand only", "both strands", "gene strand only", "low metaT"]


def read_features(path):
    rows = {}
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        ix = {k: i for i, k in enumerate(head)}
        for k in ("feature_id", "contig", "start", "end", "strand", "type", "sense_reads", "antisense_reads"):
            if k not in ix:
                sys.exit("Error: %s has no column %s" % (path, k))
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(head):
                f += [""] * (len(head) - len(f))
            rows[f[ix["feature_id"]]] = {
                "contig": f[ix["contig"]], "start": int(f[ix["start"]]), "end": int(f[ix["end"]]),
                "strand": f[ix["strand"]], "type": f[ix["type"]], "sense": float(f[ix["sense_reads"]] or 0),
                "anti": float(f[ix["antisense_reads"]] or 0), "product": f[ix["product"]] if "product" in ix else ""}
    return rows


def featurecounts(fc, bam, saf, s, threads, tmp, paired):
    out = os.path.join(tmp, "fc.s%d.txt" % s)
    cmd = [fc, "-T", str(threads), "-O", "-M", "--fraction", "-s", str(s), "-F", "SAF", "-a", saf, "-o", out, bam]
    if paired:
        cmd[1:1] = ["-p"]
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        sys.exit("Error: featureCounts -s %d failed:\n%s" % (s, r.stderr[-2000:]))
    counts = {}
    with open(out) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("Geneid\t"):
                continue
            f = line.rstrip("\n").split("\t")
            v = float(f[-1])
            if v > 0:
                counts[f[0]] = v
    return counts


def indexed(bam):
    b = bam[:-4] if bam.endswith(".bam") else bam
    return any(os.path.exists(x) for x in (bam + ".bai", bam + ".csi", b + ".bai", b + ".csi"))


def is_paired(bam):
    """True if the first records of the BAM include paired reads (IMP: pairs + singletons in one BAM)."""
    try:
        p = subprocess.Popen(["samtools", "view", "-f", "1", bam], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                             text=True)
        line = p.stdout.readline()
        p.kill()
        return bool(line)
    except OSError:
        return True


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--srna", action="append", required=True, help="TYPE=<feature_counts.tsv>")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--strandedness", default="auto", choices=["auto", "reverse", "forward"],
                    help="reverse = dUTP / TruSeq Stranded (read 1 antisense, -s 2 = gene strand)")
    ap.add_argument("--featurecounts", default="featureCounts")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--min-antisense", type=float, default=10)
    ap.add_argument("--min-mt", type=float, default=10)
    ap.add_argument("--min-stranded", type=float, default=0.75,
                    help="minimum median share of a CDS's metaT reads on the gene strand to treat the library as stranded")
    ap.add_argument("--strand-cds", type=int, default=2000,
                    help="random CDS (seed 1) counted to infer the strandedness [2000]")
    ap.add_argument("--full", action="store_true",
                    help="count on the whole BAM (default: only the alignments over the counted features, "
                         "fetched with mates from the indexed BAM)")
    a = ap.parse_args()

    tabs = {}
    for spec in a.srna:
        t, p = spec.split("=", 1)
        tabs[t] = read_features(p)
    feats = {}
    for t in tabs:
        for k, v in tabs[t].items():
            feats.setdefault(k, v)
    if not feats:
        sys.exit("Error: no features")
    os.makedirs(os.path.dirname(os.path.abspath(a.out_prefix)), exist_ok=True)
    # counted: features with antisense sRNA reads in any library type + random CDS for the strandedness
    sel = {k for t in tabs for k, v in tabs[t].items() if v["anti"] >= a.min_antisense}
    allcds = sorted(k for k, v in feats.items() if v["type"] == "CDS" and v["strand"] in "+-" and v["end"] > v["start"])
    allcds = [k for k in allcds if k not in sel]
    cds = random.Random(1).sample(allcds, min(a.strand_cds, len(allcds)))
    sel |= set(cds)
    count = [k for k in feats if k in sel and feats[k]["strand"] in "+-" and feats[k]["end"] > feats[k]["start"]]
    with tempfile.TemporaryDirectory(dir=os.path.dirname(os.path.abspath(a.out_prefix))) as tmp:
        saf, bed = os.path.join(tmp, "features.saf"), os.path.join(tmp, "regions.bed")
        with open(saf, "w") as fh, open(bed, "w") as fb:
            fh.write("GeneID\tChr\tStart\tEnd\tStrand\n")
            for k in count:
                v = feats[k]
                fh.write("%s\t%s\t%d\t%d\t%s\n" % (k, v["contig"], v["start"] + 1, v["end"], v["strand"]))
                fb.write("%s\t%d\t%d\n" % (v["contig"], v["start"], v["end"]))
        bam, how = a.bam, "whole BAM"
        if not a.full and indexed(a.bam):
            sub = os.path.join(tmp, "regions.bam")
            r = subprocess.run(["samtools", "view", "-@", str(a.threads), "-b", "-M", "-P", "-L", bed, "-o", sub, a.bam],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if r.returncode == 0 and os.path.getsize(sub) > 0:
                bam, how = sub, "alignments over %d features (with mates)" % len(count)
            else:
                print("samtools view -L failed, counting on the whole BAM:\n%s" % r.stderr[-500:], file=sys.stderr)
        elif not a.full:
            print("no BAM index (.bai/.csi) — counting on the whole BAM (samtools index %s makes this fast)" % a.bam,
                  file=sys.stderr)
        paired = is_paired(bam)
        c1 = featurecounts(a.featurecounts, bam, saf, 1, a.threads, tmp, paired)
        c2 = featurecounts(a.featurecounts, bam, saf, 2, a.threads, tmp, paired)
    print("%s: counted on the %s" % (a.sample, how), file=sys.stderr)
    # strandedness: median share of the s2 reads per random CDS with >= --min-mt reads (robust to the few genes with
    # strong antisense transcription that dominate a sum); the summed shares are reported as well
    t1, t2 = sum(c1.get(k, 0) for k in cds), sum(c2.get(k, 0) for k in cds)
    sh = sorted(c2.get(k, 0) / (c1.get(k, 0) + c2.get(k, 0)) for k in cds if c1.get(k, 0) + c2.get(k, 0) >= a.min_mt)
    med2 = (sh[len(sh) // 2] if len(sh) % 2 else (sh[len(sh) // 2 - 1] + sh[len(sh) // 2]) / 2) if sh else 0.5
    if a.strandedness == "auto":
        mode = "reverse" if med2 >= 0.5 else "forward"
    else:
        mode = a.strandedness
    own, opp = (c2, c1) if mode == "reverse" else (c1, c2)
    gene_share = med2 if mode == "reverse" else 1 - med2
    tot = t1 + t2
    sum_share = ((t2 if mode == "reverse" else t1) / tot) if tot else 0
    stranded = gene_share >= a.min_stranded
    S = a.sample
    with open(a.out_prefix + ".metaT_strand.features.tsv", "w") as fh:
        fh.write("feature_id\ttype\tmetaT_gene_strand\tmetaT_opposite_strand\n")
        for k in count:
            v = feats[k]
            o, p = own.get(k, 0), opp.get(k, 0)
            if o or p:
                fh.write("%s\t%s\t%.3f\t%.3f\n" % (k, v["type"], o, p))
    summ = [("all", "strandedness", mode if stranded else "unstranded"), ("all", "paired", int(paired)),
            ("all", "cds_metaT_gene_strand_fraction_median", "%.4f" % gene_share),
            ("all", "cds_metaT_gene_strand_fraction_sum", "%.4f" % sum_share), ("all", "cds_counted", len(sh)),
            ("all", "cds_metaT_reads_s1", "%.1f" % t1), ("all", "cds_metaT_reads_s2", "%.1f" % t2)]
    for t, rows in tabs.items():
        sel = {k: v for k, v in rows.items() if v["anti"] >= a.min_antisense}
        n = {c: 0 for c in CLASSES}
        r = {c: 0.0 for c in CLASSES}
        with open("%s.%s.metaT_strand.tsv" % (a.out_prefix, t), "w") as fh:
            fh.write("feature_id\tcontig\tstart\tend\tstrand\ttype\tproduct\tsense_reads\tantisense_reads\t"
                     "metaT_gene_strand\tmetaT_opposite_strand\tmetaT_strand_class\n")
            for k, v in sorted(sel.items(), key=lambda x: -x[1]["anti"]):
                o, p = own.get(k, 0), opp.get(k, 0)
                if not stranded:
                    c = "unstranded metaT"
                elif o >= a.min_mt and p >= a.min_mt:
                    c = "both strands"
                elif p >= a.min_mt:
                    c = "sRNA strand only"
                elif o >= a.min_mt:
                    c = "gene strand only"
                else:
                    c = "low metaT"
                if c in n:
                    n[c] += 1
                    r[c] += v["anti"]
                fh.write("%s\t%s\t%d\t%d\t%s\t%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%s\n" % (
                    k, v["contig"], v["start"], v["end"], v["strand"], v["type"], v["product"], v["sense"], v["anti"],
                    o, p, c))
        summ.append((t, "features_with_antisense", len(sel)))
        summ.append((t, "antisense_reads", "%.1f" % sum(v["anti"] for v in sel.values())))
        if stranded:
            for c in CLASSES:
                summ.append((t, "features:" + c, n[c]))
                summ.append((t, "antisense_reads:" + c, "%.1f" % r[c]))
            ra = sum(r.values())
            supp = (r["sRNA strand only"] + r["both strands"]) / ra if ra else 0
            summ.append((t, "antisense_reads_with_metaT_on_sRNA_strand_fraction", "%.4f" % supp))
            print("%s %s: %d features with >= %g antisense reads; metaT supports the sRNA strand for %.1f %% of "
                  "their antisense reads" % (S, t, len(sel), a.min_antisense, 100 * supp), file=sys.stderr)
    with open(a.out_prefix + ".metaT_strand.summary.tsv", "w") as fh:
        fh.write("sample\tlibrary_type\tmetric\tvalue\n")
        for t, m, v in summ:
            fh.write("%s\t%s\t%s\t%s\n" % (S, t, m, v))
    print("%s: metaT %s (gene-strand share of CDS reads %.3f, %s)" % (
        S, mode if stranded else "unstranded", gene_share, "paired" if paired else "single"), file=sys.stderr)


if __name__ == "__main__":
    main()
