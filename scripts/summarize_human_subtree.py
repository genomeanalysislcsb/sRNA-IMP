#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path
from typing import Optional


def read_count_table_sum(path: Path) -> float:
    """
    Sum the 'count' column from a TSV such as:
      feature\tcount
    """
    total = 0.0
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        first = True
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 2:
                continue

            # skip header if present
            if first and cols[0].lower() in {
                "feature", "mirna", "trna", "rrna", "otherrna", "id"
            }:
                first = False
                continue
            first = False

            try:
                total += float(cols[1])
            except ValueError:
                continue
    return total


def count_fastq_reads(path: Path) -> int:
    """
    Count reads in FASTQ / FASTQ.GZ as lines / 4.
    """
    opener = gzip.open if path.suffix == ".gz" else open
    n_lines = 0
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for _ in fh:
            n_lines += 1
    return n_lines // 4


def find_one(root: Path, pattern: str) -> Optional[Path]:
    matches = sorted(root.glob(pattern))
    return matches[0] if matches else None


def require_file(path: Optional[Path], label: str) -> Path:
    if path is None:
        raise FileNotFoundError(f"[ERROR] Missing required file: {label}")
    if not path.exists():
        raise FileNotFoundError(f"[ERROR] File not found: {path} ({label})")
    if path.stat().st_size == 0:
        raise ValueError(f"[ERROR] File is empty: {path} ({label})")
    return path


def optional_file(path: Optional[Path], label: str) -> Optional[Path]:
    if path is None:
        print(f"[WARN] Missing optional file: {label}")
        return None
    if not path.exists():
        print(f"[WARN] File not found: {path} ({label})")
        return None
    if path.stat().st_size == 0:
        print(f"[WARN] File is empty: {path} ({label})")
        return None
    return path


def infer_sample_name(*paths: Optional[Path]) -> str:
    for p in paths:
        if p is None:
            continue
        m = re.match(r"(.+?)\.human\.", p.name)
        if m:
            return m.group(1)
    return "sample"


def main() -> int:
    ap = argparse.ArgumentParser(description="Summarize one human smallRNA subtree.")
    ap.add_argument(
        "human_dir",
        help="Path to human subtree, e.g. BXXX4/Analysis/sRNA-Ext/human",
    )
    ap.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output summary TSV (default: <human_dir>/<sample>.human.summary.tsv)",
    )
    args = ap.parse_args()

    root = Path(args.human_dir)
    if not root.exists():
        raise FileNotFoundError(f"[ERROR] Directory does not exist: {root}")
    if not root.is_dir():
        raise NotADirectoryError(f"[ERROR] Not a directory: {root}")

    # Required files
    mirna_counts = require_file(
        find_one(root / "mirna", "*.human.mirna.counts.tsv"),
        "miRNA counts",
    )
    trna_counts = require_file(
        find_one(root / "trna", "*.human.trna.counts.tsv"),
        "tRNA counts",
    )
    rrna_counts = require_file(
        find_one(root / "rrna", "*.human.rrna.counts.tsv"),
        "rRNA counts",
    )
    oth_counts = require_file(
        find_one(root / "otherrna", "*.human.otherrna.counts.tsv"),
        "other RNA counts",
    )
    unknown_fastq = require_file(
        find_one(root / "otherrna", "*.human.unknown.fastq.gz"),
        "final unknown FASTQ",
    )

    # Optional intermediate FASTQs
    non_mirna_fastq = optional_file(
        find_one(root / "mirna", "*.human.non_mirna.fastq.gz"),
        "remaining after miRNA",
    )
    non_trna_fastq = optional_file(
        find_one(root / "trna", "*.human.non_mirna.non_trna.fastq.gz"),
        "remaining after tRNA",
    )
    non_rrna_fastq = optional_file(
        find_one(root / "rrna", "*.human.non_mirna.non_trna.non_rrna.fastq.gz"),
        "remaining after rRNA",
    )

    sample = infer_sample_name(
        mirna_counts,
        trna_counts,
        rrna_counts,
        oth_counts,
        unknown_fastq,
        non_mirna_fastq,
        non_trna_fastq,
        non_rrna_fastq,
    )

    out_path = Path(args.output) if args.output else root / f"{sample}.human.summary.tsv"

    mirna_sum = read_count_table_sum(mirna_counts)
    trna_sum = read_count_table_sum(trna_counts)
    rrna_sum = read_count_table_sum(rrna_counts)
    oth_sum = read_count_table_sum(oth_counts)

    non_mirna_reads = count_fastq_reads(non_mirna_fastq) if non_mirna_fastq else 0
    non_trna_reads = count_fastq_reads(non_trna_fastq) if non_trna_fastq else 0
    non_rrna_reads = count_fastq_reads(non_rrna_fastq) if non_rrna_fastq else 0
    unknown_reads = count_fastq_reads(unknown_fastq)

    # reconstructed estimate
    trimmed_est = mirna_sum + trna_sum + rrna_sum + oth_sum + unknown_reads

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="\n") as out:
        out.write("sample\tcategory\treads\n")
        out.write(f"{sample}\thuman_reads\t{trimmed_est:.2f}\n")
        out.write(f"{sample}\thuman_miRNA\t{mirna_sum:.2f}\n")
        out.write(f"{sample}\thuman_tRNA\t{trna_sum:.2f}\n")
        out.write(f"{sample}\thuman_rRNA\t{rrna_sum:.2f}\n")
        out.write(f"{sample}\thuman_otherRNA\t{oth_sum:.2f}\n")
        out.write(f"{sample}\tremaining_after_miRNA\t{non_mirna_reads:.2f}\n")
        out.write(f"{sample}\tremaining_after_tRNA\t{non_trna_reads:.2f}\n")
        out.write(f"{sample}\tremaining_after_rRNA\t{non_rrna_reads:.2f}\n")
        out.write(f"{sample}\tunknown_after_all\t{unknown_reads:.2f}\n")

    print(f"Wrote summary: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())