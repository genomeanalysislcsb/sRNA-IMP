#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from srna_hits import iter_reads, fmt  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True, help="Bowtie text output")
    ap.add_argument("-o", "--output", required=True, help="TSV output")
    args = ap.parse_args()

    counts = defaultdict(float)
    for c, hits in iter_reads(args.input):          # all equally good hits of a read share its count
        for cols in hits:
            counts[cols[2]] += c / len(hits)

    with open(args.output, "w") as out:
        out.write("feature\tcount\n")
        for ref in sorted(counts):
            out.write(f"{ref}\t{fmt(counts[ref])}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
    