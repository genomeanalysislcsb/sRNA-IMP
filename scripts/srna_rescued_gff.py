#!/usr/bin/env python3
"""srna_rescued_gff.py — the contig annotation (GFF3) plus the rRNA genes rescued by cmscan (srna_rrna_rescue.py),
for genome browsers or for re-use upstream (e.g. IMP). The original GFF is not changed.

  * rescued rRNA genes of all sRNA types (rescue BEDs) are merged (identical contig / strand and > 50 % overlap:
    the highest-scoring hit is kept) and added as type rRNA, source cmscan;
  * CDS on the opposite strand overlapping a rescued gene by >= 50 % of the shorter one get
    note=opposite_rescued_rRNA (likely spurious ORFs predicted on the other strand of an unannotated rRNA).

  srna_rescued_gff.py --gff <s>.gff --rescue <type1>/<s>.rRNA_rescue.bed,<type2>/... --out <s>.rRNA_rescued.gff
"""
import argparse
import collections
import os
import re
import sys


def score(attrs):
    m = re.search(r"score ([0-9.]+)", attrs)
    return float(m.group(1)) if m else 0.0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--rescue", required=True, help="comma-separated rescue BEDs")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    hits = []
    for f in a.rescue.split(","):
        if f and os.path.exists(f):
            with open(f) as fh:
                for ln in fh:
                    c = ln.rstrip("\n").split("\t")
                    if len(c) >= 10:
                        hits.append((c[0], int(c[1]), int(c[2]), c[5], c[9]))
    hits.sort(key=lambda h: -score(h[4]))
    keep = []
    for h in hits:
        if not any(h[0] == k[0] and h[3] == k[3] and min(h[2], k[2]) - max(h[1], k[1]) > 0.5 * min(h[2] - h[1], k[2] - k[1])
                   for k in keep):
            keep.append(h)
    by = collections.defaultdict(list)
    for h in keep:
        by[h[0]].append(h)
    nflag = 0
    tmp = a.out + ".part"
    with open(a.gff) as fh, open(tmp, "w") as o:
        for ln in fh:
            if ln.startswith("##FASTA"):
                break
            c = ln.rstrip("\n").split("\t")
            if len(c) == 9 and c[2] == "CDS" and c[0] in by:
                s, e = int(c[3]) - 1, int(c[4])
                for h in by[c[0]]:
                    ov = min(e, h[2]) - max(s, h[1])
                    if h[3] != c[6] and ov > 0 and ov >= 0.5 * min(e - s, h[2] - h[1]):
                        c[8] = c[8].rstrip(";") + ";note=opposite_rescued_rRNA"
                        nflag += 1
                        ln = "\t".join(c) + "\n"
                        break
            o.write(ln)
        for i, h in enumerate(sorted(keep, key=lambda h: (h[0], h[1])), 1):
            attrs = h[4].replace(" (", "; ").replace(")", "")
            m = re.search(r"ID=([^;]+)", attrs)
            rid = m.group(1) if m else "rRNA_rescue_%d" % i
            prod = re.search(r"product=([^;]+)", attrs)
            o.write("\t".join((h[0], "cmscan", "rRNA", str(h[1] + 1), str(h[2]), "%.1f" % score(h[4]), h[3], ".",
                               "ID=%s;product=%s;note=rescued by cmscan (sRNA-IMP pipeline)" % (
                                   rid, prod.group(1).split(";")[0] if prod else "rRNA"))) + "\n")
    os.replace(tmp, a.out)
    sys.stderr.write("[rescued_gff] %d rescued rRNA genes added, %d CDS flagged -> %s\n" % (len(keep), nflag, a.out))


if __name__ == "__main__":
    main()
