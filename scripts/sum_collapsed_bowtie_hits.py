#!/usr/bin/env python3
from __future__ import annotations

import argparse
from collections import defaultdict


def count_from_name(name: str) -> int:
    try:
        return int(name.rsplit("-", 1)[1])
    except Exception:
        return 1


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True, help="Bowtie text output")
    ap.add_argument("-o", "--output", required=True, help="TSV output")
    args = ap.parse_args()

    counts = defaultdict(int)

    with open(args.input) as f:
        for line in f:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            qname = cols[0]
            ref = cols[2]
            counts[ref] += count_from_name(qname)

    with open(args.output, "w") as out:
        out.write("feature\tcount\n")
        for ref in sorted(counts):
            out.write(f"{ref}\t{counts[ref]}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
    