#!/usr/bin/env python3
"""srna_collapse_refs.py — collapse identical sequences of a class reference (tRNA, Rfam) for the nonhuman screens.

The same tRNA (or Rfam region) is often identical in hundreds of genomes. With all equally good hits
(NONHUMAN_TIES=fractional) bowtie then reports every copy (one stool library: 1.1e9 tRNA alignments for 0.5 M
reads). Here every distinct sequence becomes one reference (u<N>), and a members table lists the original
headers of each; the counting scripts give each copy of a read's equally good hits the same share (a distinct
hit carried by m genomes counts m times), which is the same result as all hits on the original reference, at a
fraction of the alignments.

  srna_collapse_refs.py IN.fa OUT_PREFIX            ->  OUT_PREFIX.fa, OUT_PREFIX.members.tsv
  bowtie-build --threads 16 OUT_PREFIX.fa OUT_PREFIX
  then in the config: IDX_TRNA_ALL=OUT_PREFIX  TRNA_MEMBERS=OUT_PREFIX.members.tsv   (Rfam: IDX_RFAM_CLASSES, RFAM_MEMBERS)

The original FASTA can be recovered from an existing index: bowtie-inspect <index> > IN.fa
(bowtie-inspect writes the sequences with their full reference names; multi-line FASTA is fine).
Sequences are compared case-insensitively with U read as T.
"""
import argparse
import gzip
import sys


def read_fasta(path):
    op = gzip.open if path.endswith(".gz") else open
    name, seq = None, []
    with op(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:], []
            elif line:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fasta"); ap.add_argument("out_prefix")
    a = ap.parse_args()
    uid = {}                       # sequence -> unique id
    members = []                   # (unique id, header)
    n = 0
    for name, seq in read_fasta(a.fasta):
        n += 1
        s = seq.upper().replace("U", "T")
        if s not in uid:
            uid[s] = "u%d" % (len(uid) + 1)
        members.append((uid[s], name))
    with open(a.out_prefix + ".fa", "w") as out:
        for s, u in uid.items():
            out.write(">%s\n%s\n" % (u, s))
    members.sort(key=lambda x: (int(x[0][1:]), x[1]))
    with open(a.out_prefix + ".members.tsv", "w") as out:
        out.write("unique_id\treference\n")
        for u, name in members:
            out.write("%s\t%s\n" % (u, name))
    sys.stderr.write("[collapse_refs] %d sequences -> %d distinct (%.1f-fold) -> %s.fa, %s.members.tsv\n" % (
        n, len(uid), n / max(len(uid), 1), a.out_prefix, a.out_prefix))


if __name__ == "__main__":
    main()
