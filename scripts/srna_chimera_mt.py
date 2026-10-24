#!/usr/bin/env python3
"""srna_chimera_mt.py — chimeras in metatranscriptomic (total RNA-seq) reads mapped by IMP3 (bwa-mem2) to the
sample assembly. Companion of srna_chimera.py (same annotation, statistics, organism and complementarity tests,
same output tables), for the longer metaT fragments:

  split reads        one read aligned in two parts (primary + supplementary / secondary alignment, -M):
                     both parts >= --min-arm nt, together covering the read (<= --max-gap nt unassigned,
                     <= --max-overlap nt shared); classes inter / intra / colinear as in srna_chimera.py, plus
                     contig_edge: the junction sits at the ends of two contigs (assembly break, not used)
  discordant pairs   R1 and R2 on different contigs (both primary, MAPQ >= --min-mapq); pairs whose unseen
                     part of the fragment could run off the contig end (mate pointing to an end < --edge nt
                     away) are assembly breaks and not used; class pair_inter

RNA orientation: from the library strandedness (--strandedness auto: R1 sense / antisense to CDS in a sample of
R1 alignments; dUTP libraries are R1-antisense). Unstranded libraries: arms are annotated on either strand.
IMP3 removes rRNA reads (SortMeRNA) before mapping, so rRNA chimeras are largely absent here.

Outputs (<out-prefix>.): chimera_summary.tsv, chimera_reads.tsv.gz, chimera_pairs.tsv, chimera_typepairs.tsv,
chimera_energy.tsv, chimera_bins.tsv (as srna_chimera.py)
"""
import argparse
import collections
import os
import random
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import srna_chimera as SC  # noqa: E402

CIG = re.compile(r"(\d+)([MIDNSHP=X])")
CLASS_RANK = ["rRNA", "tRNA", "ncRNA", "CRISPR", "CDS", "novel_sRNA", "5UTR", "3UTR", "rRNA_as", "tRNA_as",
              "ncRNA_as", "CRISPR_as", "CDS_as", "intergenic"]
SUMMARY_KEYS = ("records_scanned", "split_candidates", "split_without_partner", "split_primary_clip_at_contig_edge", "chimeric_reads_inter",
                "chimeric_reads_intra", "chimeric_reads_colinear", "chimeric_reads_contig_edge",
                "discordant_pairs", "discordant_pairs_low_mapq", "discordant_pairs_contig_edge",
                "chimeric_reads_pair_inter")


class All:
    def __contains__(self, x):
        return True


def log(m):
    SC.log("metaT: " + m)


def cigar_parts(cigar):
    """-> (left clip incl. H, query aligned length, right clip, reference length)"""
    ops = [(int(n), o) for n, o in CIG.findall(cigar)]
    i = 0
    left = 0
    while i < len(ops) and ops[i][1] in "SH":
        left += ops[i][0]
        i += 1
    j = len(ops)
    right = 0
    while j > i and ops[j - 1][1] in "SH":
        right += ops[j - 1][0]
        j -= 1
    q = sum(n for n, o in ops[i:j] if o in "MI=X")
    r = sum(n for n, o in ops[i:j] if o in "MDN=X")
    return left, q, right, r


def contig_lengths(bam):
    L = {}
    out = subprocess.run(["samtools", "view", "-H", bam], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                         text=True).stdout
    for line in out.splitlines():
        if line.startswith("@SQ"):
            d = dict(x.split(":", 1) for x in line.split("\t")[1:] if ":" in x)
            L[d["SN"]] = int(d["LN"])
    return L


def detect_strandedness(bam, A, n=300000):
    """R1 primary alignments (MAPQ >= 20) overlapping one CDS: share on the CDS strand"""
    p = subprocess.Popen(["samtools", "view", "-f", "64", "-F", "2316", "-q", "20", bam], stdout=subprocess.PIPE,
                         stderr=subprocess.DEVNULL, text=True)
    sense = anti = seen = 0
    for line in p.stdout:
        f = line.split("\t", 6)
        left, q, right, r = cigar_parts(f[5])
        s = int(f[3]) - 1
        e = s + r
        st = "-" if int(f[1]) & 16 else "+"
        cds = [x for x in A.near(f[2], s, e) if x[3] == "CDS" and min(e, x[1]) - max(s, x[0]) >= 0.5 * r]
        if len(cds) == 1:
            if cds[0][2] == st:
                sense += 1
            else:
                anti += 1
        seen += 1
        if seen >= n:
            break
    p.kill()
    tot = sense + anti
    frac = sense / tot if tot else float("nan")
    if tot < 1000:
        mode = "unstranded"
    elif frac >= 0.8:
        mode = "R1-sense"
    elif frac <= 0.2:
        mode = "R1-antisense"
    else:
        mode = "unstranded"
    log("strandedness: %d R1 on single CDS, %.3f sense -> %s" % (tot, frac, mode))
    return mode, frac, tot


def flip(arm):
    n, ref, s, e, st = arm
    return (n, ref, s, e, "-" if st == "+" else "+")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True, help="IMP3 metaT BAM (Assembly/mt.reads.sorted.bam)")
    ap.add_argument("--assembly", required=True, help="assembly FASTA(.gz) the BAM was mapped to")
    ap.add_argument("--gbed", action="append", default=[], help="annotation as gff2bed (repeatable)")
    ap.add_argument("--novel-bed", help="novel sRNA loci (BED6), e.g. from an sRNA library of the sample")
    ap.add_argument("--contig2bin")
    ap.add_argument("--gtdbtk")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--strandedness", choices=("auto", "R1-sense", "R1-antisense", "unstranded"), default="auto")
    ap.add_argument("--min-arm", type=int, default=25, help="minimum aligned length of each part of a split read [25]")
    ap.add_argument("--max-gap", type=int, default=10, help="read nt covered by neither part [10]")
    ap.add_argument("--max-overlap", type=int, default=10, help="read nt covered by both parts [10]")
    ap.add_argument("--min-mapq", type=int, default=20, help="MAPQ of primary alignments / discordant mates [20]")
    ap.add_argument("--min-mapq-split", type=int, default=1, help="MAPQ of the partner part of a split read [1]")
    ap.add_argument("--edge", type=int, default=500, help="discordant mates pointing to a contig end closer than "
                    "this: assembly break [500]")
    ap.add_argument("--edge-tol", type=int, default=10, help="split junction within this of two contig ends: "
                    "assembly break [10]")
    ap.add_argument("--colinear-max", type=int, default=1000)
    ap.add_argument("--utr5", type=int, default=150)
    ap.add_argument("--utr3", type=int, default=100)
    ap.add_argument("--flank", type=int, default=20, help="contig nt around split-read arms for the duplex [20]")
    ap.add_argument("--flank-pair", type=int, default=0, help="contig nt around discordant mates for the duplex [0]")
    ap.add_argument("--null", type=int, default=20)
    ap.add_argument("--max-pairs", type=int, default=5000)
    ap.add_argument("--rnaduplex", default="RNAduplex")
    ap.add_argument("--rnaduplex-noconv", action="store_true", help=argparse.SUPPRESS)
    ap.add_argument("--max-records", type=int, default=0, help="stop after N candidate records (tests)")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    random.seed(a.seed)
    os.makedirs(os.path.dirname(a.out_prefix) or ".", exist_ok=True)
    stats = collections.Counter()

    clen = contig_lengths(a.bam)
    if not clen:
        sys.exit("no @SQ lines in %s" % a.bam)
    A = SC.Annotation(a.gbed, a.novel_bed, All(), a.utr5, a.utr3)
    if a.strandedness == "auto":
        mode, frac, ntest = detect_strandedness(a.bam, A)
        stats["strand_test_reads"] = ntest
    else:
        mode, frac = a.strandedness, float("nan")
    a.strandedness = mode

    # ---- stream candidate records
    expr = ('flag.secondary || flag.supplementary || [SA] || cigar =~ "[1-9][0-9]+S" || '
            '(rnext != "=" && rnext != "*" && !flag.munmap)')
    cmd = ["samtools", "view", "-F", "4", "-@", str(max(1, a.threads - 1)), "-e", expr, a.bam]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    split = collections.defaultdict(list)     # (qname, mate) -> alignments
    mates = collections.defaultdict(dict)     # qname -> {mate: (ref, start, end, strand, mapq)}
    n = 0
    for line in p.stdout:
        f = line.split("\t", 11)
        flag = int(f[1])
        ref, pos0, mapq, cigar = f[2], int(f[3]) - 1, int(f[4]), f[5]
        mate = 1 if flag & 64 else 2 if flag & 128 else 0
        primary = not flag & 0x900
        rev = bool(flag & 16)
        left, q, right, r = cigar_parts(cigar)
        n += 1
        L = clen.get(ref, 0)
        keep = not primary
        if primary and max(left, right) >= a.min_arm:
            # a clip that only overhangs the contig end is an assembly edge, not a chimera
            over_l = left >= a.min_arm and pos0 < left - a.edge_tol
            over_r = right >= a.min_arm and pos0 + r + right > L + a.edge_tol
            keep = (left >= a.min_arm and not over_l) or (right >= a.min_arm and not over_r)
        has_sa = len(f) > 11 and "SA:Z:" in f[11]
        if keep or has_sa:
            split[(f[0], mate)].append((ref, pos0, r, left, q, right, rev, mapq, primary,
                                        f[9] if primary else ""))
        if primary and mate and f[6] not in ("=", "*") and not flag & 8:
            mates[f[0]][mate] = (ref, pos0, pos0 + r, "-" if rev else "+", mapq, q)
        if a.max_records and n >= a.max_records:
            p.kill()
            break
    p.wait()
    if p.returncode not in (0, -9) and not a.max_records:
        sys.exit("samtools view failed: " + p.stderr.read()[-1000:])
    stats["records_scanned"] = n
    log("%d candidate records; %d reads with split candidates, %d discordant read pairs" % (
        n, len(split), sum(1 for v in mates.values() if len(v) == 2)))

    def rna_order(first, second, mate):
        """arms in read order -> RNA 5'->3' order (5' arm, 3' arm)"""
        read_is_rna = mode == "unstranded" or mate == 0 or (mode == "R1-sense") == (mate == 1)
        return (first, second) if read_is_rna else (flip(second), flip(first))

    def annotate(arm):
        _, ref, s, e, st = arm
        if mode != "unstranded":
            return A.annotate(ref, s, e, st)
        x, y = A.annotate(ref, s, e, "+"), A.annotate(ref, s, e, "-")
        return min((x, y), key=lambda t: CLASS_RANK.index(t[0]) if t[0] in CLASS_RANK else 99)

    def at_edge(a5, a3):
        """junction at the ends of two contigs: 5' arm ends at its contig end, 3' arm starts at its contig start
        (in RNA orientation)"""
        def end_of(arm, three_prime):
            _, ref, s, e, st = arm
            L = clen.get(ref, 0)
            near_end, near_start = e >= L - a.edge_tol, s <= a.edge_tol
            if st == "+":
                return near_end if three_prime else near_start
            return near_start if three_prime else near_end
        return end_of(a5, True) and end_of(a3, False)

    agg = collections.OrderedDict()            # (class, arm5, arm3) -> [reads, representative sequence, mapq5, mapq3, ov]

    # ---- split reads
    for (qn, mate), al in split.items():
        prim = [x for x in al if x[8]]
        if not prim:
            stats["split_primary_clip_at_contig_edge"] += 1   # primary only clipped over a contig end
            continue
        if len(al) < 2:
            continue
        stats["split_candidates"] += 1
        P = prim[0]
        rlen = P[3] + P[4] + P[5]

        def qint(x):
            return (x[5], x[5] + x[4]) if x[6] else (x[3], x[3] + x[4])
        pi = qint(P)
        best = None
        for O in al:
            if O is P or O[7] < a.min_mapq_split or P[7] < a.min_mapq:
                continue
            oi = qint(O)
            if O[4] < a.min_arm or P[4] < a.min_arm:
                continue
            ov = min(pi[1], oi[1]) - max(pi[0], oi[0])
            cover = max(pi[1], oi[1]) - min(pi[0], oi[0])
            if ov > a.max_overlap or rlen - cover > a.max_gap:
                continue
            if best is None or cover > best[0]:
                best = (cover, O, oi, max(0, ov))
        if best is None:
            stats["split_without_partner"] += 1
            continue
        _, O, oi, ov = best
        arms = sorted([(pi, P), (oi, O)], key=lambda t: t[0][0])
        armt = [(x[4], x[0], x[1], x[1] + x[2], "-" if x[6] else "+") for _, x in arms]
        a5, a3 = rna_order(armt[0], armt[1], mate)
        if a5[1] != a3[1] and at_edge(a5, a3):
            cl = "contig_edge"
        else:
            cl = SC.classify(a5, a3, rlen, a.colinear_max)
        stats["chimeric_reads_" + cl] += 1
        key = (cl, a5, a3)
        if key not in agg:
            agg[key] = [0, P[9] or qn, arms[0][1][7], arms[1][1][7], ov]
        agg[key][0] += 1

    # ---- discordant pairs
    for qn, m in mates.items():
        if len(m) < 2:
            continue
        stats["discordant_pairs"] += 1
        r1, r2 = m[1], m[2]
        if r1[0] == r2[0]:
            continue
        if min(r1[4], r2[4]) < a.min_mapq:
            stats["discordant_pairs_low_mapq"] += 1
            continue

        def runs_off(x):
            ref, s, e, st = x[:4]
            L = clen.get(ref, 0)
            return (L - s < a.edge) if st == "+" else (e < a.edge)
        if runs_off(r1) or runs_off(r2):
            stats["discordant_pairs_contig_edge"] += 1
            continue
        m1 = (r1[5], r1[0], r1[1], r1[2], r1[3])
        m2 = (r2[5], r2[0], r2[1], r2[2], r2[3])
        if mode == "R1-antisense":
            a5, a3 = m2, flip(m1)
        else:
            a5, a3 = m1, flip(m2)
        stats["chimeric_reads_pair_inter"] += 1
        key = ("pair_inter", a5, a3)
        if key not in agg:
            agg[key] = [0, "pair:%s:%d|%s:%d" % (a5[1], a5[2], a3[1], a3[2]), r1[4], r2[4], 0]
        agg[key][0] += 1
    split.clear()
    mates.clear()

    rows = []
    for (cl, a5, a3), (c, rep, q5, q3, ov) in agg.items():
        stats["chimeric_sequences_" + cl] += 1
        rows.append((rep, c, cl, a5, a3, annotate(a5), annotate(a3), q5, q3, ov))
    refs = {r[3][1] for r in rows if r[2] in SC.USED} | {r[4][1] for r in rows if r[2] in SC.USED}
    C = SC.Contigs(SC.read_fasta(a.assembly, refs))
    c2b = SC.read_contig2bin(a.contig2bin)
    b2t = SC.read_gtdbtk(a.gtdbtk)

    def org(ref):
        b = c2b.get(ref, "unbinned")
        return b, b2t.get(b, "")
    log("split reads: %d inter, %d intra, %d colinear, %d contig edge; discordant pairs: %d used, %d contig edge" % (
        stats["chimeric_reads_inter"], stats["chimeric_reads_intra"], stats["chimeric_reads_colinear"],
        stats["chimeric_reads_contig_edge"], stats["chimeric_reads_pair_inter"],
        stats["discordant_pairs_contig_edge"]))
    SC.write_outputs(rows, C, org, a, stats, SUMMARY_KEYS)
    with open(a.out_prefix + ".chimera_summary.tsv", "a") as fh:
        fh.write("strand_test_fraction_sense\t%s\n" % SC.fmt(frac))


if __name__ == "__main__":
    main()
