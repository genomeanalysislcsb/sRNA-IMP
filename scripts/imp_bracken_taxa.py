#!/usr/bin/env python3
"""imp_bracken_taxa.py — read-based abundance of taxa in the IMP3 metaG / metaT Bracken profiles of many samples.

  python3 imp_bracken_taxa.py --ids expobiome.ids --regex 'butyricigenes|Mediterraneibacter' \\
      --imp '/work/projects/expobiome/IMP/intra_MGMT/<id>' --level S --out cohort/imp_bracken_mediterraneibacter.tsv

<id> / <ID> in --imp are replaced by the lower- / original-case sample ID. Output: one row per sample and matching
taxon with Bracken reads and fraction for metaG (mg) and metaT (mt); samples without the taxon get 0 (taxon
listed in at least one sample).
"""
import argparse, collections, glob, gzip, os, re, sys


def read_bracken(p):
    out = {}
    if not os.path.exists(p):
        return None
    with (gzip.open(p, "rt") if p.endswith(".gz") else open(p)) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(h, line.rstrip("\n").split("\t")))
            out[r["name"].strip()] = (float(r.get("new_est_reads", 0) or 0), float(r.get("fraction_total_reads", 0) or 0))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ids", required=True, help="file with one sample ID per line")
    ap.add_argument("--imp", required=True, help="IMP3 sample directory pattern with <id> or <ID>")
    ap.add_argument("--regex", required=True, help="regex on the Bracken taxon name (case-insensitive)")
    ap.add_argument("--level", default="S", help="Bracken level file: S, G, F, ... [S]")
    ap.add_argument("--bracken-subdir", default="bracken", help="folder under Analysis/taxonomy [bracken], e.g. "
                    "bracken_k2_gtdb_genome_reps_20250609")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rx = re.compile(a.regex, re.I)
    ids = [l.split()[0] for l in open(a.ids) if l.strip() and not l.startswith("#")]
    data, names = {}, set()
    for s in ids:
        d = a.imp.replace("<ID>", s).replace("<id>", s.lower())
        for src in ("mg", "mt"):
            bd = os.path.join(d, "Analysis", "taxonomy", a.bracken_subdir)
            cand = [os.path.join(bd, src, "%s.%s.bracken%s" % (src, a.level, x)) for x in (".gz", "")]
            cand += sorted(glob.glob(os.path.join(bd, "**", "%s*.%s.bracken*" % (src, a.level)), recursive=True))
            cand += sorted(glob.glob(os.path.join(bd, "**", "%s*_%s*bracken*" % (src, a.level)), recursive=True))
            b = None
            for c in cand:
                b = read_bracken(c)
                if b is not None:
                    break
            if b is None:
                sys.stderr.write("missing %s Bracken for %s (%s)\n" % (src, s, d))
                continue
            hit = {n: v for n, v in b.items() if rx.search(n)}
            data[(s, src)] = hit
            names |= set(hit)
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w") as fh:
        fh.write("sample\ttaxon\tmg_reads\tmg_fraction\tmt_reads\tmt_fraction\n")
        for s in ids:
            for n in sorted(names):
                mg = data.get((s, "mg"), {}).get(n, (0.0, 0.0))
                mt = data.get((s, "mt"), {}).get(n, (0.0, 0.0))
                fh.write("%s\t%s\t%.0f\t%.6g\t%.0f\t%.6g\n" % (s, n, mg[0], mg[1], mt[0], mt[1]))
    sys.stderr.write("%d samples, %d matching taxa: %s\n" % (len(ids), len(names), ", ".join(sorted(names))[:500]))


if __name__ == "__main__":
    main()
