#!/usr/bin/env python3
"""
srna_fragments.py — positional analysis of tRNA / rRNA (/tmRNA) derived small RNAs

Reads all alignments of a BAM (samtools must be in PATH), assigns every
alignment to an overlapping annotated tRNA / rRNA / tmRNA feature on the
same strand and classifies the fragment by its position on the parent RNA.
Multi-mapping reads are weighted 1/n (n = number of reported alignments).

tRNA classes (positions relative to the mature tRNA gene, tolerance --tol nt):
  full_length   starts at 5' end and ends at 3' end
  5p_half       starts at 5' end, ends in the anticodon region
  3p_half       starts in the anticodon region, ends at 3' end
  5p_tRF        starts at 5' end, ends before the anticodon region
  3p_tRF        starts after the anticodon region, ends at 3' end
  i_tRF         internal
  The anticodon region is approximated as --ac-from..--ac-to (fraction of
  the tRNA length; default 0.40-0.55, i.e. ~nt 30-42 of a 76-nt tRNA).
  For fragments ending at the 3' end the presence of a terminal CCA is recorded.

rRNA / tmRNA classes: 5p_end, 3p_end, full_length, internal, plus the
position of the fragment midpoint in 10 bins along the gene.

Outputs (<prefix>.*):
  fragments.features.tsv   per feature: class counts (weighted reads)
  fragments.classes.tsv    per feature type and class: totals and fractions
  fragments.lengths.tsv    per feature type, class and read length
  fragments.ends.tsv       per position of features with >= --min-reads:
                           5' ends, 3' ends and coverage (cleavage sites)
  fragments.summary.tsv    overall numbers
"""
import argparse
import gzip
import bisect
import collections
import os
import re
import subprocess
import sys


def open_in(path):
    """open a text input, transparently gunzipping *.gz"""
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def ref_len(cigar):
    return sum(int(n) for n, op in CIGAR_RE.findall(cigar) if op in "MDN=X")


def revcomp(s):
    return s.translate(COMP)[::-1]


def parse_attrs(a):
    d = {}
    for kv in a.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k.strip()] = v.strip()
    return d


def load_features(gbed, types, contig2bin):
    """gff2bed BED: chrom start end id score strand source type phase attrs"""
    feats = []
    by_contig = collections.defaultdict(list)
    with open_in(gbed) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            ftype = f[7]
            if ftype not in types:
                continue
            chrom, start, end, fid, strand = f[0], int(f[1]), int(f[2]), f[3], f[5]
            attrs = parse_attrs(f[9]) if len(f) > 9 else {}
            if fid in (".", ""):
                fid = attrs.get("ID", "%s:%d-%d(%s)" % (chrom, start, end, strand))
            product = attrs.get("product", attrs.get("Name", attrs.get("gene", ".")))
            idx = len(feats)
            feats.append({
                "id": fid, "contig": chrom, "start": start, "end": end,
                "strand": strand, "type": ftype, "product": product,
                "bin": contig2bin.get(chrom, "unbinned"),
            })
            by_contig[chrom].append((start, end, idx))
    index = {}
    for c, lst in by_contig.items():
        lst.sort()
        starts = [s for s, _, _ in lst]
        maxlen = max(e - s for s, e, _ in lst)
        index[c] = (starts, lst, maxlen)
    return feats, index


def candidates(index, contig, rs, re_):
    if contig not in index:
        return []
    starts, lst, maxlen = index[contig]
    i = bisect.bisect_right(starts, re_ - 1)
    out = []
    j = i - 1
    while j >= 0 and starts[j] > rs - maxlen - 1:
        s, e, idx = lst[j]
        if s < re_ and e > rs:
            out.append(idx)
        j -= 1
    return out


def samtools_lines(bam, threads):
    p = subprocess.Popen(["samtools", "view", "-@", str(threads), "-F", "4", bam],
                         stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    for line in p.stdout:
        yield line
    if p.wait() != 0:
        sys.exit("samtools view failed on %s" % bam)


def classify_trna(s, e, L, tol, ac_from, ac_to):
    a5 = s <= tol
    a3 = e >= L - tol
    acl, ach = ac_from * L, ac_to * L
    if a5 and a3:
        return "full_length"
    if a5:
        if e < acl:
            return "5p_tRF"
        if e <= ach + tol:
            return "5p_half"
        return "5p_long"          # longer than a half, not reaching 3' end
    if a3:
        if s > ach:
            return "3p_tRF"
        if s >= acl - tol:
            return "3p_half"
        return "3p_long"
    return "i_tRF"


def classify_other(s, e, L, tol):
    a5 = s <= tol
    a3 = e >= L - tol
    if a5 and a3:
        return "full_length"
    if a5:
        return "5p_end"
    if a3:
        return "3p_end"
    return "internal"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--gbed", required=True, help="annotation BED from gff2bed")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--types", default="tRNA,rRNA,tmRNA")
    ap.add_argument("--contig2bin", default="", help="optional contig<TAB>bin table")
    ap.add_argument("--min-overlap", type=float, default=0.9, help="fraction of the read inside the feature")
    ap.add_argument("--tol", type=int, default=3, help="tolerance (nt) for 5'/3' ends")
    ap.add_argument("--ac-from", type=float, default=0.40)
    ap.add_argument("--ac-to", type=float, default=0.55)
    ap.add_argument("--min-reads", type=float, default=10, help="min weighted reads for per-position output")
    ap.add_argument("--threads", type=int, default=2)
    a = ap.parse_args()

    types = set(t for t in a.types.split(",") if t)
    contig2bin = read_contig2bin(a.contig2bin) if a.contig2bin else {}
    feats, index = load_features(a.gbed, types, contig2bin)
    sys.stderr.write("[fragments] %d features of types %s\n" % (len(feats), ",".join(sorted(types))))

    # alignments per read: NH tag (assembly module) or a first pass over the BAM
    nh_re = re.compile(r"\tNH:i:(\d+)")
    use_nh = False
    for line in samtools_lines(a.bam, a.threads):
        use_nh = nh_re.search(line) is not None
        break
    nali = collections.Counter()
    if not use_nh:
        for line in samtools_lines(a.bam, a.threads):
            nali[line.split("\t", 1)[0]] += 1
    total_reads = len(nali)

    trna_classes = ["full_length", "5p_half", "3p_half", "5p_tRF", "3p_tRF", "5p_long", "3p_long", "i_tRF"]
    other_classes = ["full_length", "5p_end", "3p_end", "internal"]

    fcls = collections.defaultdict(lambda: collections.defaultdict(float))   # feat -> class -> w
    fcca = collections.defaultdict(lambda: [0.0, 0.0])                       # feat -> [3'end with CCA, without]
    fbins = collections.defaultdict(lambda: [0.0] * 10)
    lens = collections.defaultdict(float)                                    # (type, class, len) -> w
    ends5 = collections.defaultdict(lambda: collections.defaultdict(float))  # feat -> pos -> w
    ends3 = collections.defaultdict(lambda: collections.defaultdict(float))
    assigned_w = 0.0
    total_w = 0.0

    # pass 2: assignment
    for line in samtools_lines(a.bam, a.threads):
        f = line.split("\t", 10)
        name, flag, contig, pos, cigar, seq = f[0], int(f[1]), f[2], int(f[3]) - 1, f[5], f[9]
        if use_nh:
            m = nh_re.search(line)
            w = 1.0 / int(m.group(1)) if m else 1.0
        else:
            w = 1.0 / nali[name]
        total_w += w
        rlen = ref_len(cigar)
        rs, re_ = pos, pos + rlen
        strand = "-" if flag & 16 else "+"
        best, best_ov = None, 0
        for idx in candidates(index, contig, rs, re_):
            ft = feats[idx]
            if ft["strand"] not in (strand, "."):
                continue
            ov = min(re_, ft["end"]) - max(rs, ft["start"])
            if ov > best_ov:
                best, best_ov = idx, ov
        if best is None or best_ov < a.min_overlap * rlen:
            continue
        ft = feats[best]
        L = ft["end"] - ft["start"]
        if ft["strand"] == "-":
            s, e = ft["end"] - re_, ft["end"] - rs
            rseq = revcomp(seq)
        else:
            s, e = rs - ft["start"], re_ - ft["start"]
            rseq = seq
        s = max(s, 0)
        e = min(e, L)
        if ft["type"] == "tRNA":
            cls = classify_trna(s, e, L, a.tol, a.ac_from, a.ac_to)
            if e >= L - a.tol:
                fcca[best][0 if rseq.upper().endswith("CCA") else 1] += w
        else:
            cls = classify_other(s, e, L, a.tol)
        mid = (s + e) / 2.0
        fbins[best][min(9, int(10 * mid / max(L, 1)))] += w
        fcls[best][cls] += w
        lens[(ft["type"], cls, len(seq))] += w
        ends5[best][s] += w
        ends3[best][max(e - 1, 0)] += w
        assigned_w += w

    p = a.out_prefix
    os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)
    # per feature
    with open(p + ".fragments.features.tsv", "w") as out:
        out.write("\t".join(["feature_id", "contig", "start", "end", "strand", "type", "product", "bin",
                             "length", "total"] + ["tRNA:" + c for c in trna_classes] +
                            ["other:" + c for c in other_classes] +
                            ["tRNA_3p_with_CCA", "tRNA_3p_without_CCA"] +
                            ["midpoint_bin%d" % i for i in range(1, 11)]) + "\n")
        for idx in sorted(fcls, key=lambda i: -sum(fcls[i].values())):
            ft = feats[idx]
            c = fcls[idx]
            tot = sum(c.values())
            row = [ft["id"], ft["contig"], ft["start"], ft["end"], ft["strand"], ft["type"], ft["product"],
                   ft["bin"], ft["end"] - ft["start"], "%.3f" % tot]
            row += ["%.3f" % (c.get(k, 0.0) if ft["type"] == "tRNA" else 0.0) for k in trna_classes]
            row += ["%.3f" % (c.get(k, 0.0) if ft["type"] != "tRNA" else 0.0) for k in other_classes]
            row += ["%.3f" % fcca[idx][0], "%.3f" % fcca[idx][1]]
            row += ["%.3f" % x for x in fbins[idx]]
            out.write("\t".join(str(x) for x in row) + "\n")

    # per type / class
    tc = collections.defaultdict(float)
    for idx, c in fcls.items():
        for k, v in c.items():
            tc[(feats[idx]["type"], k)] += v
    ttot = collections.defaultdict(float)
    for (t, k), v in tc.items():
        ttot[t] += v
    with open(p + ".fragments.classes.tsv", "w") as out:
        out.write("type\tclass\tweighted_reads\tfraction_of_type\n")
        for (t, k) in sorted(tc):
            out.write("%s\t%s\t%.3f\t%.4f\n" % (t, k, tc[(t, k)], tc[(t, k)] / ttot[t] if ttot[t] else 0))

    with open(p + ".fragments.lengths.tsv", "w") as out:
        out.write("type\tclass\tread_length\tweighted_reads\n")
        for (t, k, l) in sorted(lens):
            out.write("%s\t%s\t%d\t%.3f\n" % (t, k, l, lens[(t, k, l)]))

    with open(p + ".fragments.ends.tsv", "w") as out:
        out.write("feature_id\ttype\tproduct\tlength\tpos\tn_5p_ends\tn_3p_ends\tcoverage\n")
        for idx in sorted(fcls):
            if sum(fcls[idx].values()) < a.min_reads:
                continue
            ft = feats[idx]
            L = ft["end"] - ft["start"]
            diff = [0.0] * (L + 1)
            for pos5, v in ends5[idx].items():
                diff[pos5] += v
            for pos3, v in ends3[idx].items():
                diff[min(pos3 + 1, L)] -= v
            cov = 0.0
            for i in range(L):
                cov += diff[i]
                n5 = ends5[idx].get(i, 0.0)
                n3 = ends3[idx].get(i, 0.0)
                out.write("%s\t%s\t%s\t%d\t%d\t%.3f\t%.3f\t%.3f\n" % (ft["id"], ft["type"], ft["product"], L,
                                                                    i + 1, n5, n3, cov))

    with open(p + ".fragments.summary.tsv", "w") as out:
        out.write("metric\tvalue\n")
        out.write("mapped_reads\t%d\n" % (total_reads or int(round(total_w))))
        out.write("mapped_weighted\t%.3f\n" % total_w)
        out.write("assigned_to_%s_weighted\t%.3f\n" % ("_".join(sorted(types)), assigned_w))
        for t in sorted(ttot):
            out.write("%s_weighted\t%.3f\n" % (t, ttot[t]))
        out.write("features_with_reads\t%d\n" % len(fcls))


def read_contig2bin(path):
    d = {}
    with open_in(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and not line.startswith("#"):
                d[f[0]] = f[1]
    return d


if __name__ == "__main__":
    main()
