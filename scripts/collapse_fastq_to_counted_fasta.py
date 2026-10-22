#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
from collections import OrderedDict
from pathlib import Path


def open_maybe_gzip(path: str, mode: str="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True, help="Input FASTQ(.gz)")
    ap.add_argument("-o", "--output", required=True, help="Output FASTA")
    args = ap.parse_args()

    counts = OrderedDict()

    with open_maybe_gzip(args.input, "rt") as fin:
        while True:
            h = fin.readline()
            if not h:
                break
            seq = fin.readline().strip()
            plus = fin.readline()
            qual = fin.readline()
            if not qual:
                break
            counts[seq] = counts.get(seq, 0) + 1

    with open(args.output, "w") as out:
        i = 1
        for seq, c in counts.items():
            out.write(f">seq{i}-{c}\n{seq}\n")
            i += 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())