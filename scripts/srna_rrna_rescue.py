#!/usr/bin/env python3
"""srna_rrna_rescue.py — find rRNA genes missing from the contig annotation where they create antisense artefacts.

Reads that are sense to an unannotated (often partial) rRNA gene are antisense to ORFs that gene prediction placed
on the opposite strand ("hypothetical proteins" with many antisense and no sense reads). Read-based rRNA filters
(RiboDetector, exact Rfam sequence matches) miss many of them, because short fragments and sample-specific rRNA
variants are not recognised. This script

  1. selects CDS with >= --min-antisense antisense reads and sense <= --max-sense-ratio x antisense
     (assembly feature_counts.tsv; the *_non_rRNA columns if present), plus contigs rescued before;
  2. extracts those contigs and searches them with Infernal cmscan against the Rfam rRNA models only
     (5S, 5.8S, SSU, LSU; built once from Rfam.cm and cmpress-ed: --cm-cache);
  3. writes the hits as rRNA features in the annotation BED format (gff2bed, 10 columns) to --out.

  srna_rrna_rescue.py --feature-counts <asm>/<s>.feature_counts.tsv --fasta contigs.fa(.gz) --rfam-cm Rfam.cm \
      [--family-txt family.txt.gz] --cm-cache <dir>/Rfam.rRNA.cm --out <asm>/<s>.rRNA_rescue.bed [--cpu N]
"""
import argparse
import csv
import gzip
import os
import subprocess
import sys
import tempfile
import time

RRNA_ACC = ["RF00001", "RF00002", "RF00177", "RF01959", "RF01960", "RF02540", "RF02541", "RF02542", "RF02543",
            "RF02546"]


def op(p, mode="rt"):
    if p.endswith(".gz"):
        return gzip.open(p, mode, encoding="latin-1") if "t" in mode else gzip.open(p, mode)
    return open(p, mode, encoding="latin-1") if "b" not in mode else open(p, mode)


def rrna_accessions(family_txt):
    if family_txt and os.path.exists(family_txt):
        acc = []
        with op(family_txt) as fh:
            for ln in fh:
                c = ln.rstrip("\n").split("\t")
                if len(c) > 18 and "rRNA" in [x.strip() for x in c[18].split(";")]:
                    acc.append(c[0])
        if acc:
            return sorted(set(acc))
    return RRNA_ACC


def build_cm(rfam_cm, accs, out):
    """copy the models (CM block + HMM filter block) of the given accessions out of Rfam.cm, then cmpress"""
    if os.path.exists(out + ".i1m"):
        return
    lock = out + ".lock"
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    while True:
        try:
            os.mkdir(lock)
            break
        except FileExistsError:
            if os.path.exists(out + ".i1m"):
                return
            time.sleep(10)
    try:
        if os.path.exists(out + ".i1m"):
            return
        want = set(accs)
        tmp = out + ".tmp"
        n = 0
        with op(rfam_cm) as fh, open(tmp, "w") as o:
            block, keep = [], False
            for ln in fh:
                block.append(ln)
                if ln.startswith("ACC "):
                    keep = ln.split()[1] in want
                if ln.startswith("//"):
                    if keep:
                        o.writelines(block)
                        if block[0].startswith("INFERNAL"):
                            n += 1
                    block, keep = [], False
        if not n:
            raise SystemExit("no rRNA models found in %s" % rfam_cm)
        for ext in ("", ".i1m", ".i1i", ".i1f", ".i1p"):
            if os.path.exists(out + ext):
                os.remove(out + ext)
        os.rename(tmp, out)
        subprocess.run(["cmpress", out], check=True, stdout=subprocess.DEVNULL)
        sys.stderr.write("[rescue] %d rRNA models -> %s\n" % (n, out))
    finally:
        os.rmdir(lock)


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return 0.0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--feature-counts", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--rfam-cm", required=True)
    ap.add_argument("--family-txt", default="")
    ap.add_argument("--cm-cache", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--annotation", default="", help="annotation BED (gff2bed): hits overlapping an annotated rRNA on "
                    "the same strand by > 50 %% are not added")
    ap.add_argument("--min-antisense", type=float, default=20)
    ap.add_argument("--max-sense-ratio", type=float, default=0.1)
    ap.add_argument("--evalue", type=float, default=1e-5)
    ap.add_argument("--cpu", type=int, default=4)
    a = ap.parse_args()

    contigs, nanti = set(), 0.0
    with open(a.feature_counts) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        nr = "antisense_reads_non_rRNA" in (rd.fieldnames or [])
        for r in rd:
            if r.get("type") != "CDS":
                continue
            an = num(r["antisense_reads_non_rRNA" if nr else "antisense_reads"])
            se = num(r["sense_reads_non_rRNA" if nr else "sense_reads"])
            if an >= a.min_antisense and se <= a.max_sense_ratio * an:
                contigs.add(r["contig"])
                nanti += an
    old = []
    if os.path.exists(a.out):                     # keep earlier rescues (idempotent re-runs)
        with open(a.out) as fh:
            old = [ln for ln in fh if ln.strip()]
        contigs |= {ln.split("\t", 1)[0] for ln in old}
    sys.stderr.write("[rescue] %d contigs with antisense-only CDS (%.0f antisense reads)\n" % (len(contigs), nanti))
    if not contigs:
        open(a.out, "w").close()
        return
    build_cm(a.rfam_cm, rrna_accessions(a.family_txt), a.cm_cache)
    with tempfile.TemporaryDirectory() as td:
        fa = os.path.join(td, "contigs.fa")
        n = 0
        with op(a.fasta) as fh, open(fa, "w") as o:
            keep = False
            for ln in fh:
                if ln.startswith(">"):
                    keep = ln[1:].split()[0] in contigs
                    n += keep
                if keep:
                    o.write(ln)
        tbl = os.path.join(td, "hits.tbl")
        subprocess.run(["cmscan", "--cpu", str(a.cpu), "--cut_ga", "--noali", "--tblout", tbl, a.cm_cache, fa],
                       check=True, stdout=subprocess.DEVNULL)
        hits = []
        with open(tbl) as fh:
            for ln in fh:
                if ln.startswith("#"):
                    continue
                c = ln.split()
                # target name, accession, query name, acc, mdl, mdl from, to, seq from, to, strand, trunc, pass,
                # gc, bias, score, E-value, inc, description
                if len(c) < 17 or c[16] != "!" or num(c[15]) > a.evalue:
                    continue
                s, e, strand = int(c[7]), int(c[8]), c[9]
                if s > e:
                    s, e = e, s
                hits.append((c[2], s - 1, e, strand, c[0], c[1], c[14], c[15], c[10]))
    known = []
    if a.annotation and os.path.exists(a.annotation):
        with op(a.annotation) as fh:
            for ln in fh:
                c = ln.split("\t")
                if len(c) > 7 and c[0] in contigs and c[7] == "rRNA" and c[6] != "cmscan":
                    known.append((c[0], int(c[1]), int(c[2]), c[5]))

    def is_known(h):
        return any(h[0] == k[0] and h[3] == k[3] and min(h[2], k[2]) - max(h[1], k[1]) > 0.5 * (h[2] - h[1])
                   for k in known)
    hits = [h for h in hits if not is_known(h)]
    # best (highest score) non-overlapping hit per region and strand
    hits.sort(key=lambda h: -num(h[6]))
    chosen = []
    for h in hits:
        if not any(h[0] == g[0] and h[3] == g[3] and min(h[2], g[2]) - max(h[1], g[1]) > 0 for g in chosen):
            chosen.append(h)
    chosen.sort(key=lambda h: (h[0], h[1]))
    with open(a.out, "w") as o:
        for i, h in enumerate(chosen, 1):
            fid = "rRNA_rescue_%s_%d_%d_%s" % (h[0], h[1] + 1, h[2], "plus" if h[3] == "+" else "minus")
            attrs = "ID=%s;product=%s (%s; rescued by cmscan, score %s, E %s%s)" % (
                fid, h[4], h[5], h[6], h[7], "" if h[8] == "no" else ", truncated " + h[8])
            o.write("\t".join(map(str, (h[0], h[1], h[2], fid, ".", h[3], "cmscan", "rRNA", ".", attrs))) + "\n")
    sys.stderr.write("[rescue] %d contigs searched, %d rRNA genes added (%s)\n" % (n, len(chosen), a.out))


if __name__ == "__main__":
    main()
