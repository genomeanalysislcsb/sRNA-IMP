#!/usr/bin/env python3
"""srna_hits.py — shared reader for bowtie (default output) hits of collapsed reads (module nonhuman).

Collapsed read names end in -<count>. With all equally good hits (bowtie -a --best --strata) a read has
n lines; iter_reads() yields (count, [hit columns, ...]) per read, so callers can give each hit count/n. With
-k 1 every read has one hit and the weights stay whole. bowtie 1 writes all hits of a read together, also with
-p (its --reorder works only with SAM output); a read whose hits are not adjacent stops with an error instead of
being counted twice.

first_out: write one line per read, the hit with the smallest reference name (deterministic, unlike bowtie's
choice among ties), in bowtie's format: a small hits file for later steps that only need one hit per read.

Collapsed references (srna_collapse_refs.py): the index holds each distinct sequence once (u<N>) and a members
table lists the original references of each. load_members() reads it; expand() spreads the count of each
distinct reference equally over its members, once at the end (not per read). In first_out the distinct id is
replaced by its first member, so later steps see the original reference names.
"""
import sys


def count_from_name(name):
    try:
        return int(name.rsplit("-", 1)[1])
    except (IndexError, ValueError):
        return 1


def load_members(path):
    """unique id -> [original reference, ...]; empty dict without a members table (references are their own)"""
    m = {}
    if path:
        with open(path) as fh:
            next(fh, None)
            for line in fh:
                u, ref = line.rstrip("\n").split("\t", 1)
                m.setdefault(u, []).append(ref)
    return m


def expand(counts, members):
    """{reference or unique id: count} -> [(original reference, count share), ...]"""
    for u, v in counts.items():
        refs = members.get(u, (u,))
        w = v / len(refs)
        for r in refs:
            yield r, w


def iter_reads(path, first_out=None, members=None):
    fh = sys.stdin if path == "-" else open(path)
    fo = open(first_out, "w") if first_out else None
    cur, hits = None, []
    seen = set()

    def flush():
        if hits:
            if fo:
                h = list(min(hits, key=lambda h: h[2]))
                if members:
                    h[2] = members.get(h[2], [h[2]])[0]
                fo.write("\t".join(h) + "\n")
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
                if cols[0] in seen:
                    sys.exit("srna_hits: hits of read %s are not adjacent in %s (run bowtie with -p 1)" % (cols[0], path))
                seen.add(cols[0])
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
