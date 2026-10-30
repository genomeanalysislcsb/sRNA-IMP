#!/usr/bin/env python3
"""srna_mirna_like.py — where do the miRNA-like reads come from: bacteria, human, or neither?

Input: the miRNA-like screen of module nonhuman (<s>.mirna_like.hits.txt: non-human, non-tRNA reads of
MIRNA_LIKE_MINLEN-MAXLEN nt that match a mature miRNA of miRBase exactly; collapsed reads named <id>-<count>).
Every distinct sequence is then mapped (bowtie 1) against
  - all miRBase mature miRNAs (exact, all hits): every species / miRNA it is identical to
  - the human mature miRNAs (<= 1 mismatch): the nearest human miRNA
  - the human genome (<= --mismatches, both strands): could it still be human?
  - the sample's own assembly (<= --mismatches, all hits): is it encoded in the microbial contigs?
and classified:
  contig_only     on the contigs, not in the human genome  -> support for a microbial origin
  human_only      in the human genome, not on the contigs  -> probably human (or food) after all
  both            in both                                  -> ambiguous (often short, low-complexity or conserved)
  neither         in neither                                -> diet / other eukaryotes / sequence variants
Reads shorter than the host filter's minimum (HOST_MIN_LEN, 20 nt) were never tested against the human genome by
module split; this script does test them, so human_only also catches such short human miRNAs (e.g. miR-1246, 19 nt).
Matching is exact by default (--mismatches 0): even then an 18-19 nt sequence has a ~5-10% chance of a random exact
hit in the human genome, so 'both' / 'human_only' calls of single short sequences are weak evidence on their own.

Outputs (<out-prefix>.): mirna_like.sequences.tsv (one row per sequence), mirna_like.families.tsv (one row per
miRNA family), mirna_like.summary.tsv (reads and sequences per class)

Human rescue (--human-rescue, default on): sequences that are clearly human — identical to a human mature miRNA AND
exact in the human genome AND not on the contigs (human_only with an exact human miRNA) — are moved to the human
statistics by the readers (srna_report.py TypeData: report, QC, cohort):
  mirna_like.human.counts.tsv   feature (human miRNA, as in human/mirna/<s>.human.mirna.counts.tsv), count
                                (reads split equally over identical human miRNAs, e.g. let-7a / let-7c)
  mirna_like.human_rescue.tsv   category, reads: rescued_reads, rescued_sequences, rescued_below_host_min_len,
                                length_<n> (for the human length distribution)
Without --human-rescue both files are removed, so the readers use the original numbers.
"""
import argparse
import collections
import os
import re
import subprocess
import sys
import tempfile

CLASSES = ("contig_only", "human_only", "both", "neither")


def count_of(name):
    try:
        return int(name.rsplit("-", 1)[1])
    except (IndexError, ValueError):
        return 1


def n_mm(field):
    return len([x for x in field.split(",") if x.strip()]) if field.strip() else 0


def family(name):
    """'dno-miR-192-5p MIMAT…' -> ('miR-192', 'animal'); 'gma-miR166c-3p' -> ('miR166', 'plant')"""
    nm = name.split()[0]
    rest = nm.split("-", 1)[1] if "-" in nm else nm
    rest = re.sub(r"-(5p|3p)$", "", rest)
    m = re.match(r"(miR|MIR)(\d+)", rest)
    if m:
        return "miR%s" % m.group(2), "plant"
    m = re.match(r"(miR|let|lin|bantam|iab)-?(\d+)?", rest)
    if m:
        base = rest.split("-")[0] if rest.startswith(("let", "lin")) else "-".join(rest.split("-")[:2])
        base = re.sub(r"[a-z]$", "", base) if re.match(r"miR-\d+[a-z]$", base) else base
        return base, "animal"
    return rest, "other"


def bowtie(idx, fa, args, threads):
    cmd = ["bowtie", "-f", "-p", str(threads)] + args + [idx, fa]
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if r.returncode:
        sys.exit("bowtie failed on %s: %s" % (idx, r.stderr[-800:]))
    hits = collections.defaultdict(list)
    for line in r.stdout.splitlines():
        f = line.split("\t")
        if len(f) >= 8:
            hits[f[0]].append((f[2], int(f[3]), f[1], n_mm(f[7])))
    return hits


def write_rescue(a, rows):
    """clearly human sequences (identical to a human miRNA, exact in the human genome, not on the contigs)"""
    cf, sf = a.out_prefix + ".mirna_like.human.counts.tsv", a.out_prefix + ".mirna_like.human_rescue.tsv"
    if not a.human_rescue:
        for f in (cf, sf):
            if os.path.exists(f):
                os.remove(f)
        return
    cnt, lens = collections.Counter(), collections.Counter()
    n = nseq = short = 0
    for r in rows:
        if r["human_rescued"] != "yes":
            continue
        for h in r["_hsa"]:
            cnt[h] += r["reads"] / len(r["_hsa"])
        lens[r["length"]] += r["reads"]
        n += r["reads"]
        nseq += 1
        short += r["reads"] if r["length"] < a.host_min_len else 0
    with open(cf, "w") as fh:
        fh.write("feature\tcount\n")
        for k, v in sorted(cnt.items(), key=lambda kv: -kv[1]):
            fh.write("%s\t%.2f\n" % (k, v))
    with open(sf, "w") as fh:
        fh.write("category\treads\n")
        fh.write("rescued_reads\t%d\nrescued_sequences\t%d\nrescued_below_host_min_len\t%d\n" % (n, nseq, short))
        for l in sorted(lens):
            fh.write("length_%d\t%d\n" % (l, lens[l]))


def read_contig2bin(path):
    d = {}
    if path and os.path.exists(path):
        import gzip
        with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 2 and not line.startswith("#"):
                    d[f[0]] = re.sub(r"\.(fa|fna|fasta)(\.gz)?$", "", f[1])
    return d


def read_gtdb(paths):
    t = {}
    import gzip
    for p in [x for x in (paths or "").split(",") if x and os.path.exists(x)]:
        with (gzip.open(p, "rt") if p.endswith(".gz") else open(p)) as fh:
            h = None
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if h is None or f[0] in ("user_genome", "Bin", "bin"):
                    h = f
                    continue
                r = dict(zip(h, f))
                key = r.get("user_genome") or r.get("Bin") or r.get("bin") or f[0]
                cl = [x for x in r.get("classification", "").split(";") if len(x) > 3]
                t[re.sub(r"\.(fa|fna|fasta)(\.gz)?$", "", key)] = cl[-1][3:] if cl else ""
    return t


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hits", required=True, help="<s>.mirna_like.hits.txt of module nonhuman")
    ap.add_argument("--mirbase-idx", required=True, help="bowtie index of all mature miRNAs (IDX_MIRNA_EUK)")
    ap.add_argument("--hsa-mirna-idx", required=True, help="bowtie index of human mature miRNAs (IDX_HSA_MIRNA)")
    ap.add_argument("--host-idx", required=True, help="bowtie index of the human genome")
    ap.add_argument("--assembly-idx", required=True, help="bowtie index of the sample's assembly")
    ap.add_argument("--contig2bin")
    ap.add_argument("--gtdbtk")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--mismatches", type=int, default=0,
                    help="mismatches for human genome and contigs [0]; with 1 mismatch an 18-nt sequence has ~5 "
                         "random hits in a 3-Gb genome, exact ~0.1")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--host-min-len", type=int, default=20, help="HOST_MIN_LEN of module split [20]")
    ap.add_argument("--human-rescue", type=int, default=1, choices=(0, 1),
                    help="1 = write the human rescue tables (clearly human sequences -> human statistics) [1]")
    a = ap.parse_args()

    seqs = collections.OrderedDict()            # sequence -> [reads, assigned miRBase name]
    with open(a.hits) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            s = f[4].upper().replace("U", "T")
            d = seqs.setdefault(s, [0, f[2]])
            d[0] += count_of(f[0])
    os.makedirs(os.path.dirname(a.out_prefix) or ".", exist_ok=True)
    if not seqs:
        for sfx, hdr in ((".mirna_like.sequences.tsv", "sequence\treads\n"), (".mirna_like.families.tsv", "family\treads\n"),
                         (".mirna_like.summary.tsv", "class\treads\tsequences\n")):
            open(a.out_prefix + sfx, "w").write(hdr)
        write_rescue(a, [])
        return
    tmp = tempfile.mkdtemp(prefix="mirna_like.", dir=os.path.dirname(a.out_prefix) or ".")
    fa = os.path.join(tmp, "seqs.fa")
    keys = list(seqs)
    with open(fa, "w") as fh:
        for i, s in enumerate(keys):
            fh.write(">s%d\n%s\n" % (i, s))
    mb = bowtie(a.mirbase_idx, fa, ["-v", "0", "-a", "--norc", "--fullref"], a.threads)
    hs = bowtie(a.hsa_mirna_idx, fa, ["-v", "1", "-a", "--best", "--strata", "--norc", "--fullref"], a.threads)
    hg = bowtie(a.host_idx, fa, ["-v", str(a.mismatches), "-k", "2", "--best"], a.threads)
    asm = bowtie(a.assembly_idx, fa, ["-v", str(a.mismatches), "-k", "10", "--best"], a.threads)
    for x in os.listdir(tmp):
        os.remove(os.path.join(tmp, x))
    os.rmdir(tmp)
    c2b, b2t = read_contig2bin(a.contig2bin), read_gtdb(a.gtdbtk)

    rows = []
    for i, s in enumerate(keys):
        sid = "s%d" % i
        reads, assigned = seqs[s]
        names = sorted({h[0].split()[0] for h in mb.get(sid, [])}) or [assigned.split()[0]]
        fam, style = family(assigned)
        hsa_exact = [n for n in names if n.startswith("hsa-")]
        near = hs.get(sid, [])
        h_mm = min((h[3] for h in hg.get(sid, [])), default=None)
        c_hits = asm.get(sid, [])
        c_mm = min((h[3] for h in c_hits), default=None)
        best = [h for h in c_hits if h[3] == c_mm]
        contig = best[0][0].split()[0] if best else ""
        bins = sorted({c2b.get(h[0].split()[0], "unbinned") for h in best}) if best else []
        taxa = sorted({b2t.get(b, "") for b in bins if b != "unbinned"} - {""})
        if c_mm is not None and h_mm is None:
            cls = "contig_only"
        elif c_mm is None and h_mm is not None:
            cls = "human_only"
        elif c_mm is not None:
            cls = "both"
        else:
            cls = "neither"
        rows.append(dict(sequence=s, reads=reads, length=len(s), mirbase=assigned.split()[0], family=fam, style=style,
                         identical_mirnas=";".join(names[:8]) + (";…" if len(names) > 8 else ""),
                         identical_species=len({n.split("-")[0] for n in names}),
                         human_mirna_exact=";".join(hsa_exact) or "-",
                         nearest_human_mirna="%s (%d mm)" % (near[0][0].split()[0], near[0][3]) if near else "-",
                         human_genome_mm="-" if h_mm is None else h_mm, contig_hits=len(c_hits),
                         contig_mm="-" if c_mm is None else c_mm, contig=contig or "-",
                         bins=";".join(bins) or "-", taxa=";".join(taxa) or "-", origin=cls,
                         human_rescued="yes" if cls == "human_only" and near and near[0][3] == 0 else "no",
                         _hsa=sorted({h[0].split()[0] for h in near if h[3] == 0})))
    rows.sort(key=lambda r: -r["reads"])
    write_rescue(a, rows)
    cols = [c for c in rows[0] if not c.startswith("_")]
    with open(a.out_prefix + ".mirna_like.sequences.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    fam = collections.OrderedDict()
    for r in rows:
        g = fam.setdefault((r["family"], r["style"]), dict(reads=0, seqs=0, cls=collections.Counter(),
                                                            hsa=set(), taxa=collections.Counter()))
        g["reads"] += r["reads"]
        g["seqs"] += 1
        g["cls"][r["origin"]] += r["reads"]
        if r["human_mirna_exact"] != "-":
            g["hsa"].update(r["human_mirna_exact"].split(";"))
        for t in r["taxa"].split(";"):
            if t != "-":
                g["taxa"][t] += r["reads"]
    with open(a.out_prefix + ".mirna_like.families.tsv", "w") as fh:
        fh.write("family\tstyle\treads\tsequences\t" + "\t".join("reads_%s" % c for c in CLASSES) +
                 "\tshare_contig_only\thuman_mirnas_identical\ttop_taxa\n")
        for (f, st), g in sorted(fam.items(), key=lambda kv: -kv[1]["reads"]):
            fh.write("%s\t%s\t%d\t%d\t%s\t%.3f\t%s\t%s\n" % (
                f, st, g["reads"], g["seqs"], "\t".join(str(g["cls"][c]) for c in CLASSES),
                g["cls"]["contig_only"] / g["reads"], ";".join(sorted(g["hsa"])) or "-",
                ";".join(t for t, _ in g["taxa"].most_common(3)) or "-"))
    tot = collections.Counter()
    nseq = collections.Counter()
    for r in rows:
        tot[r["origin"]] += r["reads"]
        nseq[r["origin"]] += 1
    with open(a.out_prefix + ".mirna_like.summary.tsv", "w") as fh:
        fh.write("class\treads\tsequences\tshare_reads\n")
        allr = sum(tot.values())
        for c in CLASSES:
            fh.write("%s\t%d\t%d\t%.4f\n" % (c, tot[c], nseq[c], tot[c] / allr if allr else 0))
    resc = sum(r["reads"] for r in rows if r["human_rescued"] == "yes")
    sys.stderr.write("[mirna_like] human rescue: %s\n" % (
        "%d reads (%d sequences) -> human statistics" % (resc, sum(r["human_rescued"] == "yes" for r in rows))
        if a.human_rescue else "off"))
    sys.stderr.write("[mirna_like] %d sequences, %d reads: %s\n" % (len(rows), sum(tot.values()), ", ".join(
        "%s %d" % (c, tot[c]) for c in CLASSES)))


if __name__ == "__main__":
    main()
