#!/usr/bin/env python3

# for human

# python3 sam_fractional_counts.py -i results/human/BXXX1/01_mirna/BXXX1.mirna.sam -o results/human/BXXX1/BXXX1.mirna.counts.tsv
# python3 sam_fractional_counts.py -i results/human/BXXX1/02_trna/BXXX1.trna.sam   -o results/human/BXXX1/BXXX1.trna.counts.tsv
# python3 sam_fractional_counts.py -i results/human/BXXX1/03_rrna/BXXX1.rrna.sam   -o results/human/BXXX1/BXXX1.rrna.counts.tsv

# for non-human

# python3 sam_fractional_counts.py -i results/nonhuman/BXXX1/03_trna/BXXX1.trna.sam        -o results/nonhuman/BXXX1/BXXX1.trna.counts.tsv
# python3 sam_fractional_counts.py -i results/nonhuman/BXXX1/04_mirna_like/BXXX1.mirna_like.sam -o results/nonhuman/BXXX1/BXXX1.mirna_like.counts.tsv
# python3 sam_fractional_counts.py -i results/nonhuman/BXXX1/05_rfam/BXXX1.rfam.sam        -o results/nonhuman/BXXX1/BXXX1.rfam.counts.tsv

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple


def count_fractional_from_sam(sam_path: Path) -> Dict[str, float]:
    """
    Fractional counting per read (QNAME):
    - Skip headers.
    - For each read, collect all mapped RNAMEs in this SAM.
    - If a read maps to N targets, each target gets +1/N.
    """
    # qname -> list of rname hits
    hits: Dict[str, List[str]] = defaultdict(list)

    with sam_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line or line.startswith("@"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            qname = cols[0]
            rname = cols[2]
            if rname == "*" or rname == "":
                continue
            hits[qname].append(rname)

    counts: Dict[str, float] = defaultdict(float)
    for qname, rnames in hits.items():
        n = len(rnames)
        if n == 0:
            continue
        w = 1.0 / n
        for r in rnames:
            counts[r] += w

    return counts


def write_counts(counts: Dict[str, float], out_path: Path, decimals: int = 6) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="\n") as out:
        out.write("feature\tcount\n")
        for feat in sorted(counts.keys()):
            out.write(f"{feat}\t{counts[feat]:.{decimals}f}\n")


def main() -> int:
    ap = argparse.ArgumentParser(description="Fractionally count SAM alignments per feature (RNAME) for uncollapsed reads.")
    ap.add_argument("-i", "--input", required=True, help="Input SAM")
    ap.add_argument("-o", "--output", required=True, help="Output TSV")
    ap.add_argument("--decimals", type=int, default=6, help="Decimals for counts (default 6)")
    args = ap.parse_args()

    counts = count_fractional_from_sam(Path(args.input))
    write_counts(counts, Path(args.output), decimals=args.decimals)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
