#!/usr/bin/env python3
"""srna_rrna_seqs.py — read sequences with an Rfam rRNA hit (module nonhuman, step 7), for leaving rRNA reads out
of the antisense / novel analyses.

RiboDetector misses many short rRNA fragments (in stool sRNA-seq often ~40 % of the reads on rRNA genes); the
Rfam screen of the remaining reads catches most of them, but works on collapsed reads (no read names). This
script writes the sequences of those hits in BOTH orientations (one per line, gzip), so a SAM SEQ field can be
looked up directly whatever strand the read mapped to.

  srna_rrna_seqs.py -i <prefix>.rfam.hits.txt -o <prefix>.rfam_rRNA.seqs.gz [--classes rRNA]
"""
import argparse
import gzip
import sys

RC = str.maketrans("ACGTUNacgtun", "TGCAANtgcaan")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--hits", required=True, help="bowtie hits of the Rfam class screen (ref: class|RF…|region)")
    ap.add_argument("-o", "--out", required=True, help="output .gz: one sequence per line (both orientations)")
    ap.add_argument("--classes", default="rRNA", help="comma-separated Rfam classes to write [rRNA]")
    a = ap.parse_args()
    want = set(a.classes.split(","))
    seqs = set()
    reads = 0
    with open(a.hits) as fh:
        for ln in fh:
            c = ln.rstrip("\n").split("\t")
            if len(c) < 5 or c[2].split("|", 1)[0] not in want:
                continue
            s = c[4].upper().replace("U", "T")
            seqs.add(s)
            seqs.add(s.translate(RC)[::-1])
            try:
                reads += int(c[0].rsplit("-", 1)[1])
            except (IndexError, ValueError):
                reads += 1
    with gzip.open(a.out, "wt", compresslevel=3) as out:
        for s in seqs:
            out.write(s + "\n")
    sys.stderr.write("[rrna_seqs] %d reads (%d sequences, both strands) with an Rfam %s hit -> %s\n" % (
        reads, len(seqs), "/".join(sorted(want)), a.out))


if __name__ == "__main__":
    main()
