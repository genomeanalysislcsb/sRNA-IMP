#!/usr/bin/env python3
"""srna_genome_marker.py — follow one reference genome (e.g. a species missing from the read classifiers) across
all samples via the sample assemblies, for sRNA-Ext, sRNA-Int, metaT and metaG at once.

For every sample, the contigs of its assembly are aligned to the reference genome(s) (minimap2 -x asm10). Contigs
with >= --min-identity identity over >= --min-cov of their length (>= --min-len nt) count as the species. Their
reads are summed from the comparison step (<s>/Analysis/comparison/<s>.all_sources.contigs.tsv: counts and CPM per
source). The share of the genome covered by the assembly is reported as an assembly-based presence measure.

  python3 srna_genome_marker.py --outdir sRNA_v2 --ids expobiome.ids --genome M_butyricigenes.fna.gz \\
      --tag butyricigenes --out sRNA_v2/Cohort/genome_marker_butyricigenes.tsv --threads 16

Per-sample results are cached in <s>/Analysis/genome_marker/<tag>.tsv (re-runs only redo missing samples).
"""
import argparse, collections, os, subprocess, sys


def intervals_len(iv):
    tot, cur = 0, None
    for s, e in sorted(iv):
        if cur is None or s > cur[1]:
            if cur:
                tot += cur[1] - cur[0]
            cur = [s, e]
        else:
            cur[1] = max(cur[1], e)
    if cur:
        tot += cur[1] - cur[0]
    return tot


def run_sample(a, s, genome_len):
    cache = os.path.join(a.outdir, s, "Analysis", "genome_marker", "%s.tsv" % a.tag)
    if os.path.exists(cache) and not a.force:
        with open(cache) as fh:
            h = fh.readline().rstrip("\n").split("\t")
            return dict(zip(h, fh.readline().rstrip("\n").split("\t")))
    asm = os.path.join(a.outdir, s, "Assembly", "%s.mgmt.assembly.merged.fa.gz" % s)
    if not os.path.exists(asm):
        sys.stderr.write("%s: no assembly (%s)\n" % (s, asm))
        return None
    cmd = [a.minimap2, "-x", "asm10", "-c", "--secondary=no", "-t", str(a.threads), a.genome, asm]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    qlen, aligned, tiv = {}, collections.Counter(), collections.defaultdict(list)
    for line in p.stdout:
        f = line.split("\t")
        q, ql, qs, qe, t, ts, te, nm, al = f[0], int(f[1]), int(f[2]), int(f[3]), f[5], int(f[7]), int(f[8]), int(f[9]), int(f[10])
        if al == 0 or nm / al < a.min_identity:
            continue
        qlen[q] = ql
        aligned[q] += qe - qs
        tiv[q].append((t, ts, te))
    p.wait()
    if p.returncode:
        sys.stderr.write("%s: minimap2 failed\n" % s)
        return None
    hits = {q for q in aligned if qlen[q] >= a.min_len and aligned[q] / qlen[q] >= a.min_cov}
    gcov = collections.defaultdict(list)
    for q in hits:
        for t, ts, te in tiv[q]:
            gcov[t].append((ts, te))
    gfrac = sum(intervals_len(v) for v in gcov.values()) / genome_len
    res = collections.OrderedDict(sample=s, contigs=len(hits), contig_bp=sum(qlen[q] for q in hits),
                                  genome_fraction="%.4f" % gfrac)
    tab = os.path.join(a.outdir, s, "Analysis", "comparison", "%s.all_sources.contigs.tsv" % s)
    sums = collections.Counter()
    if os.path.exists(tab):
        with open(tab) as fh:
            h = fh.readline().rstrip("\n").split("\t")
            cols = [c for c in h if c.endswith(("_count", "_CPM"))]
            for line in fh:
                r = dict(zip(h, line.rstrip("\n").split("\t")))
                if r["contig"] in hits:
                    for c in cols:
                        try:
                            sums[c] += float(r[c])
                        except ValueError:
                            pass
            for c in cols:
                res[c] = "%.3f" % sums[c]
    else:
        sys.stderr.write("%s: no comparison table (%s) — only assembly-based values\n" % (s, tab))
    os.makedirs(os.path.dirname(cache), exist_ok=True)
    with open(cache, "w") as fh:
        fh.write("\t".join(res) + "\n" + "\t".join(str(v) for v in res.values()) + "\n")
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True, help="sRNA-IMP output directory")
    ap.add_argument("--ids", required=True, help="file with one sample ID per line")
    ap.add_argument("--genome", required=True, help="reference genome(s) FASTA(.gz), several concatenated is fine")
    ap.add_argument("--tag", required=True, help="short name for the cache files")
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--min-identity", type=float, default=0.95, help="alignment identity (species level) [0.95]")
    ap.add_argument("--min-cov", type=float, default=0.5, help="share of the contig aligned [0.5]")
    ap.add_argument("--min-len", type=int, default=500, help="minimum contig length [500]")
    ap.add_argument("--force", action="store_true", help="ignore cached per-sample results")
    a = ap.parse_args()
    import gzip
    glen = 0
    with (gzip.open(a.genome, "rt") if a.genome.endswith(".gz") else open(a.genome)) as fh:
        for line in fh:
            if not line.startswith(">"):
                glen += len(line.strip())
    ids = [l.split()[0] for l in open(a.ids) if l.strip() and not l.startswith("#")]
    rows = []
    for i, s in enumerate(ids, 1):
        r = run_sample(a, s, glen)
        if r:
            rows.append(r)
            sys.stderr.write("[%d/%d] %s: %s contigs, genome fraction %s\n" % (i, len(ids), s, r["contigs"],
                                                                               r["genome_fraction"]))
    cols = []
    for r in rows:
        cols += [c for c in r if c not in cols]
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    sys.stderr.write("written %s (%d samples, genome %d nt)\n" % (a.out, len(rows), glen))


if __name__ == "__main__":
    main()
