#!/usr/bin/env python3
"""Compare two sRNA-IMP output folders for the same samples, e.g. a run with the environments and a run in the
container, or before / after a change of the pipeline.

  python3 scripts/srna_compare_runs.py REF_OUTDIR NEW_OUTDIR [-S SAMPLE ...] [--tol 1e-6] [--out diff.tsv]

For every table (*.tsv, *.tsv.gz, Kraken reports) below <outdir>/<sample>/{Stats,Analysis,Report} that exists in
both runs: rows are matched on the text columns, numbers compared. A table is 'identical', 'same' (numbers equal
within --tol, e.g. rounding), or 'differs' (rows only in one run, largest relative difference, sum of the numbers,
Spearman of the first numeric column). Read numbers per step (read_stats.tsv) are printed side by side. Tables
only in one run are listed. Exit code 1 when tables differ or are missing in the new run."""
import argparse
import glob
import gzip
import os
import sys

import pandas as pd

SKIP_PARTS = ("/tmp/", "/fastqc/", "/multiqc/", "/logs/")
MAX_BYTES = 300 * 1024 * 1024


def tables(base):
    out = {}
    for sub in ("Stats", "Analysis", "Report"):
        root = os.path.join(base, sub)
        for p in glob.glob(os.path.join(root, "**", "*"), recursive=True):
            if not os.path.isfile(p) or any(s in p + "/" for s in SKIP_PARTS):
                continue
            if p.endswith((".tsv", ".tsv.gz", ".report")) and os.path.getsize(p) <= MAX_BYTES:
                out[os.path.relpath(p, base)] = p
    return out


def read(p):
    opener = gzip.open if p.endswith(".gz") else open
    with opener(p, "rt", errors="replace") as fh:
        first = fh.readline()
    header = None if p.endswith(".report") else "infer"
    return pd.read_csv(p, sep="\t", header=header, dtype=str, keep_default_na=False, comment=None,
                       skiprows=1 if first.startswith("#") and not first.startswith("#\t") else 0)


def numeric_cols(d):
    num = []
    for c in d.columns:
        v = pd.to_numeric(d[c].replace("", None), errors="coerce")
        if len(d) and v.notna().sum() == (d[c] != "").sum() and v.notna().any():
            num.append(c)
    return num


def compare(pr, pn, tol):
    a, b = read(pr), read(pn)
    if a.equals(b):
        return "identical", ""
    if list(a.columns) != list(b.columns):
        return "differs", "columns differ: %s / %s" % (len(a.columns), len(b.columns))
    num = [c for c in numeric_cols(a) if c in numeric_cols(b)]
    keys = [c for c in a.columns if c not in num]
    if not num:
        return "differs", "text differs (%d / %d rows)" % (len(a), len(b))
    if keys:
        A = a.assign(**{c: pd.to_numeric(a[c]) for c in num}).groupby(keys, dropna=False)[num].sum()
        B = b.assign(**{c: pd.to_numeric(b[c]) for c in num}).groupby(keys, dropna=False)[num].sum()
    else:
        if len(a) != len(b):
            return "differs", "rows %d / %d" % (len(a), len(b))
        A = a[num].apply(pd.to_numeric)
        B = b[num].apply(pd.to_numeric)
    only_a, only_b = A.index.difference(B.index), B.index.difference(A.index)
    both = A.index.intersection(B.index)
    x, y = A.loc[both, num].astype(float), B.loc[both, num].astype(float)
    rel = ((x - y).abs() / x.abs().where(x.abs() > 0, 1)).max().max() if len(both) else 0.0
    if not len(only_a) and not len(only_b) and (rel <= tol or pd.isna(rel)):
        return "same", "numbers equal within %g" % tol
    c0 = num[0]
    sa, sb = A[c0].sum(), B[c0].sum()
    rho = x[c0].rank().corr(y[c0].rank()) if len(both) > 2 else float("nan")
    info = "rows %d / %d (only ref %d, only new %d), max rel. diff %.3g, sum %s %.6g / %.6g (%+.2f %%), Spearman %.4f" % (
        len(A), len(B), len(only_a), len(only_b), rel, c0, sa, sb, 100.0 * (sb - sa) / sa if sa else 0.0, rho)
    return "differs", info


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ref")
    ap.add_argument("new")
    ap.add_argument("-S", "--sample", action="append", default=[], help="sample(s)  [all samples of the new run]")
    ap.add_argument("--tol", type=float, default=1e-6, help="relative tolerance for 'same'  [1e-6]")
    ap.add_argument("--out", default="", help="write the comparison per table as TSV")
    a = ap.parse_args(argv)
    samples = a.sample or sorted(d for d in os.listdir(a.new) if os.path.isdir(os.path.join(a.new, d, "Stats"))
                                 or os.path.isdir(os.path.join(a.new, d, "Analysis")))
    rows, bad = [], 0
    for s in samples:
        R, N = tables(os.path.join(a.ref, s)), tables(os.path.join(a.new, s))
        print("== %s: %d tables in ref, %d in new" % (s, len(R), len(N)))
        for rs in sorted(glob.glob(os.path.join(a.new, s, "Stats", "*", "*.read_stats.tsv"))):
            rel = os.path.relpath(rs, os.path.join(a.new, s))
            if rel in R:
                x, y = read(R[rel]), read(rs)
                print("   read numbers %s" % rel)
                k = list(x.columns[:-1])
                m = x.merge(y, on=k, how="left", sort=False, suffixes=(" ref", " new"))      # order of the steps
                extra = y[~y.set_index(k).index.isin(x.set_index(k).index)]
                m = pd.concat([m, extra.rename(columns={y.columns[-1]: y.columns[-1] + " new"})], ignore_index=True)
                print("\n".join("     " + ln for ln in m.to_string(index=False).splitlines()))
        for rel in sorted(set(R) | set(N)):
            if rel not in N:
                st, info = "missing in new", ""
            elif rel not in R:
                st, info = "only in new", ""
            else:
                try:
                    st, info = compare(R[rel], N[rel], a.tol)
                except Exception as e:                       # one unreadable table does not stop the comparison
                    st, info = "error", "%s: %s" % (type(e).__name__, e)
            rows.append((s, rel, st, info))
            if st in ("differs", "missing in new", "error"):
                bad += 1
        counts = pd.Series([r[2] for r in rows if r[0] == s]).value_counts()
        print("   " + ", ".join("%s %d" % (k, v) for k, v in counts.items()))
        for _, rel, st, info in [r for r in rows if r[0] == s and r[2] not in ("identical", "same")]:
            print("   %-15s %s%s" % (st, rel, ("  " + info) if info else ""))
    if a.out:
        pd.DataFrame(rows, columns=["sample", "table", "status", "details"]).to_csv(a.out, sep="\t", index=False)
    print("\n%d tables compared, %d differ or are missing in the new run" % (len(rows), bad))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
