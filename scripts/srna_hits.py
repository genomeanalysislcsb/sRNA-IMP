#!/usr/bin/env python3
"""srna_hits.py — shared reader for bowtie (default output) hits of collapsed reads (module nonhuman).

Collapsed read names end in -<count>. With all equally good hits (bowtie -a --best --strata --reorder) a read has
n lines; iter_reads() yields (count, [hit columns, ...]) per read, so callers can give each hit count/n. With
-k 1 every read has one hit and the weights stay whole. The hits of one read must be adjacent (bowtie --reorder
or -p 1); a read that appears twice is counted twice.

first_out: write one line per read, the hit with the smallest reference name (deterministic, unlike bowtie's
choice among ties), in bowtie's format: a small hits file for later steps that only need one hit per read.
"""
import sys


def count_from_name(name):
    try:
        return int(name.rsplit("-", 1)[1])
    except (IndexError, ValueError):
        return 1


def iter_reads(path, first_out=None):
    fh = sys.stdin if path == "-" else open(path)
    fo = open(first_out, "w") if first_out else None
    cur, hits = None, []

    def flush():
        if hits:
            if fo:
                fo.write("\t".join(min(hits, key=lambda h: h[2])) + "\n")
            return count_from_name(cur), hits
        return None

    try:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            if cols[0] != cur:
                g = flush()
                if g:
                    yield g
                cur, hits = cols[0], []
            hits.append(cols)
        g = flush()
        if g:
            yield g
    finally:
        if fh is not sys.stdin:
            fh.close()
        if fo:
            fo.close()


def fmt(x):
    """counts: whole numbers without decimals, fractional ones with 3"""
    return str(int(round(x))) if abs(x - round(x)) < 1e-9 else "%.3f" % x
