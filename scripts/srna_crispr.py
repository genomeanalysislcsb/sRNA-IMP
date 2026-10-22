#!/usr/bin/env python3
"""
srna_crispr.py — small RNAs on CRISPR arrays (crRNAs) from the assembly BAM

For every CRISPR array of the annotation (Bakta: CRISPR, crispr-repeat, crispr-spacer; arrays have no
strand) the sRNA reads of the array are counted (weight 1/NH), the transcribed strand is inferred from the
read strands, and the reads are assigned to the spacers / repeats of the array. Mature crRNAs consist of a
spacer part and a repeat part; a read counts for the spacer it overlaps by >= MIN_SPACER_OVERLAP nt
(largest overlap), otherwise for the repeat it overlaps most. 5' and 3' read ends relative to the spacer
(in transcription orientation) show the processing pattern (e.g. type II: 5' end inside the spacer, 3'
end ~20 nt into the following repeat). Cas genes within CAS_WINDOW nt on the same contig give the
CRISPR-Cas type (cas9 -> II, cas3 -> I, cas10/csm/cmr -> III, cas12 -> V, cas13 -> VI).

Inputs : --bam (assembly module, NH tag), --gbed (gff2bed BED of the annotation with attributes),
         optional --contig2bin, --gtdbtk (bin taxonomy), --feature-counts (assembly feature counts:
         reads on the cas genes)
Outputs: <prefix>.crispr_arrays.tsv   one row per array (reads, strand, type, cas genes, bin, taxon)
         <prefix>.crispr_units.tsv    one row per spacer / repeat in transcription order (reads, sequence)
         <prefix>.crispr_ends.tsv     5'/3' read-end offsets relative to the spacer (array, end, offset, reads)
         <prefix>.crispr_lengths.tsv  read length of reads on arrays (with spacer part / repeat only)
         <prefix>.crispr_genes.tsv    Cas genes near the arrays (coordinates, strand, name, sense reads)
         <prefix>.crispr_coverage.tsv per-base read coverage (+ / - strand, weighted) of the --plot-top arrays
                                      with their Cas genes and --flank nt on each side (locus plots)
         <prefix>.crispr_tracr.tsv    tracrRNA candidates of type II / V arrays (needs --fasta): sRNA read
                                      clusters in the CRISPR-Cas locus outside the array, tested for an
                                      anti-repeat (local alignment with the crRNA repeat, G·U allowed; the reads
                                      + TRACR_EXT5 nt upstream; empirical p from shuffled sequences)
"""
import argparse
import bisect
import collections
import gzip
import os
import re
import subprocess
import sys
import urllib.parse

MIN_SPACER_OVERLAP = 10
TRACR_FLANK = 1500          # tracrRNA search: CRISPR-Cas locus (array + Cas genes) +/- this many nt
TRACR_MIN_READS = 10        # weighted reads per read cluster
TRACR_LEN = (30, 300)       # cluster length (nt)
TRACR_SHUFFLES = 500
TRACR_EXT5 = 30             # the anti-repeat is cut inside the repeat:anti-repeat duplex (RNase III) -> the mature
                            # tracrRNA often starts within it: search also this many nt upstream of the reads
CAS_WINDOW = 10000
NH_RE = re.compile(r"\tNH:i:(\d+)")
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def open_in(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def ref_len(cigar):
    return sum(int(n) for n, op in CIGAR_RE.findall(cigar) if op in "MDN=X")


def attrs_of(s):
    out = {}
    for x in s.split(";"):
        if "=" in x:
            k, v = x.split("=", 1)
            out[k.strip()] = urllib.parse.unquote(v.strip())
    return out


def cas_type(genes):
    g = " ".join(genes).lower()
    t = []
    if "cas9" in g or "csn2" in g:
        t.append("II")
    type3 = "cas10" in g or re.search(r"\bcs[mr]\d|cmr\d|csm\d", g)
    # type I: cas3 (signature), cas8 (large subunit) or Cascade cas5/cas7 without cas10 (I-A..I-G, e.g. cas5c/cas7c/cas8c = I-C)
    if "cas3" in g or "cas8" in g or re.search(r"\bcs[aey]\d|cse\d|csy\d|csd\d|csa\d", g) \
            or (re.search(r"cas[57]", g) and not type3):
        sub = re.search(r"cas8([a-g])\d?|cas[57]([a-g])\b", g)
        t.append("I" + ("-" + (sub.group(1) or sub.group(2)).upper() if sub else ""))
    if type3:
        t.append("III")
    if "cas12" in g or "cpf1" in g:
        t.append("V")
    if "cas13" in g or "c2c2" in g:
        t.append("VI")
    if not t and re.search(r"cas[12]\b|cas[12] ", g):
        t.append("adaptation only")
    return "/".join(t) if t else ("no cas genes" if not genes else "unclassified")


def load(gbed):
    arrays, units, cds = {}, collections.defaultdict(list), collections.defaultdict(list)
    with open_in(gbed) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) >= 10:                       # gff2bed: chrom start end id score strand src type phase attr
                c, s, e, strand, ftype, at = f[0], int(f[1]), int(f[2]), f[5], f[7], f[9]
            elif len(f) >= 9:                      # GFF (1-based)
                c, s, e, strand, ftype, at = f[0], int(f[3]) - 1, int(f[4]), f[6], f[2], f[8]
            else:
                continue
            if ftype not in ("CRISPR", "crispr-repeat", "crispr-spacer", "CDS"):
                continue
            a = attrs_of(at)
            if ftype == "CRISPR":
                name = a.get("Name", a.get("product", ""))
                m = re.search(r"(\d+) repeats of length (\d+).*consensus sequence ([ACGTUN]+).*spacer length (\d+)", name)
                arrays[a.get("ID", "%s:%d-%d" % (c, s, e))] = {
                    "contig": c, "start": s, "end": e, "n_repeats": m.group(1) if m else "",
                    "repeat_len": m.group(2) if m else "", "consensus": m.group(3) if m else "",
                    "spacer_len": m.group(4) if m else ""}
            elif ftype == "CDS":
                gname = a.get("gene", "")
                prod = a.get("product", "")
                if re.search(r"\bcas\d|csn2|cs[aemryd]\d|cmr\d|cse\d|csy\d|cpf1|crispr", (gname + " " + prod).lower()):
                    cds[c].append((s, e, strand, a.get("ID", ""), gname or prod))
            else:
                units[a.get("Parent", "")].append({"id": a.get("ID", ""), "kind": ftype.split("-")[1], "start": s,
                                                   "end": e, "seq": a.get("sequence", "")})
    for aid, us in units.items():
        us.sort(key=lambda u: u["start"])
    return arrays, units, cds


def samtools(args):
    p = subprocess.Popen(["samtools"] + args, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    for line in p.stdout:
        yield line
    if p.wait() != 0:
        sys.exit("samtools %s failed" % " ".join(args))


def read_c2b(path):
    d = {}
    if path:
        with open_in(path) as fh:
            for ln in fh:
                f = ln.rstrip("\n").split("\t")
                if len(f) >= 2 and not ln.startswith("#"):
                    d[f[0]] = re.sub(r"\.(fa|fna|fasta)(\.gz)?$", "", f[1])
    return d


def read_tax(paths):
    tax = {}
    for pth in [p for p in (paths or "").split(",") if p]:
        with open_in(pth) as fh:
            header = None
            for ln in fh:
                f = ln.rstrip("\n").split("\t")
                if header is None or f[0] in ("user_genome", "Bin", "bin"):
                    header = [h.lower() for h in f]
                    continue
                rec = dict(zip(header, f))
                b = re.sub(r"\.(fa|fna|fasta)(\.gz)?$", "", rec.get("user_genome", rec.get("bin", f[0])))
                if rec.get("classification"):
                    tax[b] = rec["classification"]
    return tax


# ------------------------------------------------------------------ tracrRNA candidates (type II / V)
COMP = str.maketrans("ACGTUN", "TGCAAN")


def rc(x):
    return x.upper().translate(COMP)[::-1]


def read_contigs(path, wanted):
    seqs, cur, buf = {}, None, []
    with open_in(path) as fh:
        for ln in fh:
            if ln.startswith(">"):
                if cur is not None:
                    seqs[cur] = "".join(buf)
                name = ln[1:].split()[0]
                cur, buf = (name if name in wanted else None), []
            elif cur is not None:
                buf.append(ln.strip().upper())
    if cur is not None:
        seqs[cur] = "".join(buf)
    return seqs


def pair_score(a, b):
    """a: candidate base, b: base of the reverse complement of the crRNA repeat (Watson-Crick = match);
    G·U wobble: a=G vs b=A (repeat U), a=T vs b=C (repeat G)"""
    if a == b:
        return 2
    if (a, b) in (("G", "A"), ("T", "C")):
        return 1
    return -3


def sw(x, y, gap=-5):
    """Smith-Waterman local alignment: best score, x range, y range, pairs, wobbles, alignment strings"""
    n, m = len(x), len(y)
    H = [[0] * (m + 1) for _ in range(n + 1)]
    best, bi, bj = 0, 0, 0
    for i in range(1, n + 1):
        xi, Hi, Hp = x[i - 1], H[i], H[i - 1]
        for j in range(1, m + 1):
            v = max(0, Hp[j - 1] + pair_score(xi, y[j - 1]), Hp[j] + gap, Hi[j - 1] + gap)
            Hi[j] = v
            if v > best:
                best, bi, bj = v, i, j
    i, j, ax, ay, wc, wb = bi, bj, [], [], 0, 0
    while i > 0 and j > 0 and H[i][j] > 0:
        v = H[i][j]
        sc = pair_score(x[i - 1], y[j - 1])
        if v == H[i - 1][j - 1] + sc:
            ax.append(x[i - 1]); ay.append(y[j - 1])
            wc += sc == 2
            wb += sc == 1
            i, j = i - 1, j - 1
        elif v == H[i - 1][j] + gap:
            ax.append(x[i - 1]); ay.append("-"); i -= 1
        else:
            ax.append("-"); ay.append(y[j - 1]); j -= 1
    return best, (i, bi), (j, bj), wc, wb, "".join(reversed(ax)), "".join(reversed(ay))


def shuffle_p(x, y, obs, n, rnd):
    if obs <= 0:
        return 1.0
    xs, ge = list(x), 0
    for _ in range(n):
        rnd.shuffle(xs)
        if sw("".join(xs), y)[0] >= obs:
            ge += 1
    return (ge + 1) / (n + 1)


def clusters(reads, gap=10):
    """reads: [(s, e, w)] of one strand -> [(start, end, weighted reads, [reads])]"""
    out = []
    for r in sorted(reads):
        if out and r[0] <= out[-1][1] + gap:
            c = out[-1]
            out[-1] = (c[0], max(c[1], r[1]), c[2] + r[2], c[3] + [r])
        else:
            out.append((r[0], r[1], r[2], [r]))
    return out


def tracr_candidates(arrays, units, cds, arr_rows, weighted, seqs, shuffles=TRACR_SHUFFLES):
    import random
    rnd = random.Random(1)
    out = []
    for row in arr_rows:
        aid, ctype, tstrand = row[0], row[9], row[17]
        if not ("II" in ctype.split("/") or "V" in ctype.split("/")):
            continue
        ar = arrays[aid]
        c = ar["contig"]
        if c not in seqs:
            continue
        seq = seqs[c]
        reps = [u for u in units.get(aid, []) if u["kind"] == "repeat"]
        rep = seq[reps[0]["start"]:reps[0]["end"]] if reps else ar["consensus"].replace("U", "T")
        if not rep:
            continue
        crrep = rep if tstrand == "+" else rc(rep)            # repeat as part of the crRNA (RNA orientation)
        target = rc(crrep)                                     # an anti-repeat looks like this
        near = [(s_, e_) for (s_, e_, st_, gid, g) in cds.get(c, [])
                if s_ < ar["end"] + CAS_WINDOW and e_ > ar["start"] - CAS_WINDOW]
        lo = max(0, min([ar["start"]] + [x[0] for x in near]) - TRACR_FLANK)
        hi = min(len(seq), max([ar["end"]] + [x[1] for x in near]) + TRACR_FLANK)
        for strand in "+-":
            rr = [(s_, e_, w) for (cc, s_, e_, st, w) in weighted
                  if cc == c and st == strand and s_ >= lo and e_ <= hi
                  and not (s_ < ar["end"] + 30 and e_ > ar["start"] - 30)]
            for cs, ce, cw, members in clusters(rr):
                if cw < TRACR_MIN_READS or not (TRACR_LEN[0] <= ce - cs <= TRACR_LEN[1]):
                    continue
                xs_, xe_ = (max(0, cs - TRACR_EXT5), ce) if strand == "+" else (cs, min(len(seq), ce + TRACR_EXT5))
                x = seq[xs_:xe_] if strand == "+" else rc(seq[xs_:xe_])
                ext = (cs - xs_) if strand == "+" else (xe_ - ce)  # nt of x upstream of the reads
                score, (xa, xb), (ya, yb), wc, wb, alx, aly = sw(x, target)
                ar5 = (xs_ + xa + 1) if strand == "+" else (xe_ - xa)          # anti-repeat 5' end (1-based)
                p = shuffle_p(x, target, score, shuffles, rnd)
                e5 = collections.Counter()
                e3 = collections.Counter()
                for s_, e_, w in members:
                    e5[s_ if strand == "+" else e_] += w
                    e3[e_ if strand == "+" else s_] += w
                m5, w5 = e5.most_common(1)[0]
                m3, w3 = e3.most_common(1)[0]
                alen = len(alx)
                cand = alen >= 15 and wc + wb >= 0.75 * alen and p <= 0.01
                # position relative to the locus: which side of the array / cas9
                cas9 = [(s_, e_, st_) for (s_, e_, st_, gid, g) in cds.get(c, []) if "cas9" in g.lower()
                        or "csn1" in g.lower() or "cas12" in g.lower()]
                rel = ""
                if cas9:
                    s9, e9, st9 = cas9[0]
                    up = (ce <= s9) if st9 == "+" else (cs >= e9)
                    rel = ("overlaps " if cs < e9 and ce > s9 else ("upstream of " if up else "downstream of ")) + \
                        "cas9" + (" (same strand)" if strand == st9 else " (opposite strand)")
                out.append([aid, c, cs, ce, strand, ce - cs, "%.1f" % cw,
                            m5 + (1 if strand == "+" else 0), "%.2f" % (w5 / cw), m3 + (0 if strand == "+" else 1),
                            "%.2f" % (w3 / cw), rel, score, alen, wc, wb, "%.4f" % p, "yes" if cand else "no",
                            x, crrep, alx, aly, ar5, max(0, ext - xa)])
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--gbed", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--contig2bin", default="")
    ap.add_argument("--gtdbtk", default="")
    ap.add_argument("--feature-counts", default="", help="assembly <s>.feature_counts.tsv (reads on cas genes)")
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--plot-top", type=int, default=5, help="coverage for the N arrays with most reads")
    ap.add_argument("--flank", type=int, default=500)
    ap.add_argument("--fasta", default="", help="contigs (FASTA, .gz ok): tracrRNA candidates of type II / V arrays")
    ap.add_argument("--tracr-shuffles", type=int, default=TRACR_SHUFFLES)
    a = ap.parse_args()

    arrays, units, cds = load(a.gbed)
    c2b, tax = read_c2b(a.contig2bin), read_tax(a.gtdbtk)
    fc = {}
    if a.feature_counts and os.path.exists(a.feature_counts):
        with open(a.feature_counts) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            for ln in fh:
                r = dict(zip(hdr, ln.rstrip("\n").split("\t")))
                fc[r["feature_id"]] = float(r.get("sense_reads", 0) or 0)
    # index arrays per contig
    by_contig = collections.defaultdict(list)
    for aid, ar in arrays.items():
        by_contig[ar["contig"]].append((ar["start"], ar["end"], aid))
    for c in by_contig:
        by_contig[c].sort()
    starts = {c: [x[0] for x in v] for c, v in by_contig.items()}

    # pass over the BAM: only contigs carrying an array
    reads = collections.defaultdict(list)    # aid -> [(strand, s, e, w, len)]
    use_nh = None
    nali = collections.Counter()
    contigs = sorted(by_contig)
    if not contigs:
        sys.stderr.write("[crispr] no CRISPR arrays in the annotation\n")
    recs = []
    for c in contigs:
        for ln in samtools(["view", "-F", "4", a.bam, c]):
            f = ln.split("\t", 10)
            if f[2] != c:
                continue
            if use_nh is None:
                use_nh = NH_RE.search(ln) is not None
            m = NH_RE.search(ln)
            nh = int(m.group(1)) if m else 0
            pos = int(f[3]) - 1
            e = pos + ref_len(f[5])
            strand = "-" if int(f[1]) & 16 else "+"
            recs.append((c, pos, e, strand, nh, f[0], len(f[9])))
            if not m:
                nali[f[0]] += 1
    weighted = []
    for c, s, e, strand, nh, name, L in recs:
        w = 1.0 / (nh if nh else max(nali[name], 1))
        weighted.append((c, s, e, strand, w))
        lst = by_contig[c]
        i = bisect.bisect_right(starts[c], e) - 1
        while i >= 0 and lst[i][1] > s - 20000:
            st, en, aid = lst[i]
            ov = min(e, en) - max(s, st)
            if ov >= 0.5 * (e - s):
                reads[aid].append((strand, s, e, w, L))
                break
            i -= 1

    p = a.out_prefix
    os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)
    arr_out, unit_out, end_out = [], [], collections.defaultdict(float)
    len_out = collections.defaultdict(float)
    for aid, ar in sorted(arrays.items(), key=lambda kv: -sum(r[3] for r in reads.get(kv[0], []))):
        rr = reads.get(aid, [])
        tot = sum(r[3] for r in rr)
        plus = sum(r[3] for r in rr if r[0] == "+")
        tstrand = "+" if plus >= tot - plus else "-"
        us = units.get(aid, [])
        order = us if tstrand == "+" else list(reversed(us))
        ucount = collections.defaultdict(float)
        n_sp_reads = 0.0
        for strand, s, e, w, L in rr:
            if strand != tstrand:
                continue
            best, bov, bkind = None, 0, ""
            for u in us:
                ov = min(e, u["end"]) - max(s, u["start"])
                if ov <= 0:
                    continue
                if u["kind"] == "spacer" and ov >= MIN_SPACER_OVERLAP and (bkind != "spacer" or ov > bov):
                    best, bov, bkind = u, ov, "spacer"
                elif bkind != "spacer" and ov > bov:
                    best, bov, bkind = u, ov, u["kind"]
            if best is None:
                continue
            ucount[best["id"]] += w
            len_out[(L, "with_spacer" if bkind == "spacer" else "repeat_only")] += w
            if bkind == "spacer":
                n_sp_reads += w
                if tstrand == "+":
                    o5, o3 = s - best["start"], e - best["end"]
                else:
                    o5, o3 = best["end"] - e, best["start"] - s
                end_out[(aid, "5p", o5)] += w
                end_out[(aid, "3p", o3)] += w
        genes = [g for (s_, e_, st_, gid, g) in cds.get(ar["contig"], [])
                 if s_ < ar["end"] + CAS_WINDOW and e_ > ar["start"] - CAS_WINDOW]
        cas_reads = sum(fc.get(gid, 0.0) for (s_, e_, st_, gid, g) in cds.get(ar["contig"], [])
                        if s_ < ar["end"] + CAS_WINDOW and e_ > ar["start"] - CAS_WINDOW)
        cas_strand = collections.Counter(st_ for (s_, e_, st_, gid, g) in cds.get(ar["contig"], [])
                                         if s_ < ar["end"] + CAS_WINDOW and e_ > ar["start"] - CAS_WINDOW)
        b = c2b.get(ar["contig"], "unbinned")
        sp = [u for u in us if u["kind"] == "spacer"]
        arr_out.append([aid, ar["contig"], ar["start"], ar["end"], ar["end"] - ar["start"], ar["n_repeats"],
                        ar["repeat_len"], ar["spacer_len"], ar["consensus"], cas_type(genes),
                        ",".join(sorted(set(genes))), cas_strand.most_common(1)[0][0] if cas_strand else "",
                        "%.3f" % cas_reads, b, tax.get(b, ""), "%.3f" % tot,
                        "%.3f" % (1e3 * tot / max(ar["end"] - ar["start"], 1)), tstrand,
                        "%.3f" % (max(plus, tot - plus) / tot if tot else 0), "%.3f" % n_sp_reads,
                        str(len(sp)), str(sum(1 for u in sp if ucount.get(u["id"], 0) >= 1))])
        for k, u in enumerate(order, 1):
            unit_out.append([aid, str(k), u["id"], u["kind"], ar["contig"], str(u["start"]), str(u["end"]),
                             tstrand, "%.3f" % ucount.get(u["id"], 0.0), u["seq"]])

    # Cas genes and per-base coverage of the top arrays (locus plots)
    gene_rows, cov_rows = [], []
    ranked = [r[0] for r in arr_out]
    for aid in ranked:
        ar = arrays[aid]
        for (s_, e_, st_, gid, g) in cds.get(ar["contig"], []):
            if s_ < ar["end"] + CAS_WINDOW and e_ > ar["start"] - CAS_WINDOW:
                gene_rows.append([aid, gid, ar["contig"], str(s_), str(e_), st_, g, "%.3f" % fc.get(gid, 0.0)])
    for aid in [x for x in ranked if sum(r[3] for r in reads.get(x, [])) > 0][:a.plot_top]:
        ar = arrays[aid]
        near = [(s_, e_) for (s_, e_, st_, gid, g) in cds.get(ar["contig"], [])
                if s_ < ar["end"] + CAS_WINDOW and e_ > ar["start"] - CAS_WINDOW]
        lo = max(0, min([ar["start"]] + [x[0] for x in near]) - a.flank)
        hi = max([ar["end"]] + [x[1] for x in near]) + a.flank
        n = hi - lo + 1
        dp, dm = [0.0] * (n + 1), [0.0] * (n + 1)
        for c, s_, e_, strand, w in weighted:
            if c != ar["contig"] or e_ <= lo or s_ >= hi:
                continue
            d = dp if strand == "+" else dm
            d[max(s_, lo) - lo] += w
            d[min(e_, hi) - lo] -= w
        cp = cm = 0.0
        last = None
        for i in range(n):
            cp += dp[i]
            cm += dm[i]
            cur = (round(cp, 3), round(cm, 3))
            if cur != last or i == n - 1:          # run-length: write only where coverage changes
                cov_rows.append("%s\t%s\t%d\t%.3f\t%.3f\n" % (aid, ar["contig"], lo + i, cur[0], cur[1]))
                last = cur
    with open(p + ".crispr_genes.tsv", "w") as out:
        out.write("array_id\tgene_id\tcontig\tstart\tend\tstrand\tname\tsense_reads\n")
        for r in gene_rows:
            out.write("\t".join(r) + "\n")
    with open(p + ".crispr_coverage.tsv", "w") as out:
        out.write("array_id\tcontig\tpos\tplus\tminus\n")   # coverage from pos until the next row
        out.writelines(cov_rows)
    with open(p + ".crispr_arrays.tsv", "w") as out:
        out.write("\t".join(["array_id", "contig", "start", "end", "length", "repeats", "repeat_len", "spacer_len",
                             "repeat_consensus", "cas_type", "cas_genes", "cas_strand", "cas_reads", "bin",
                             "classification", "reads", "reads_per_kb", "transcribed_strand", "strand_fraction",
                             "reads_with_spacer", "spacers", "spacers_with_reads"]) + "\n")
        for r in arr_out:
            out.write("\t".join(str(x) for x in r) + "\n")
    with open(p + ".crispr_units.tsv", "w") as out:
        out.write("array_id\tindex\tunit_id\tkind\tcontig\tstart\tend\ttranscribed_strand\treads\tsequence\n")
        for r in unit_out:
            out.write("\t".join(r) + "\n")
    with open(p + ".crispr_ends.tsv", "w") as out:
        out.write("array_id\tend\toffset\treads\n")
        for (aid, end, off), w in sorted(end_out.items()):
            out.write("%s\t%s\t%d\t%.3f\n" % (aid, end, off, w))
    with open(p + ".crispr_lengths.tsv", "w") as out:
        out.write("length\tclass\treads\n")
        for (L, cl), w in sorted(len_out.items()):
            out.write("%d\t%s\t%.3f\n" % (L, cl, w))
    tr = []
    if a.fasta:
        want = {arrays[r[0]]["contig"] for r in arr_out if ("II" in r[9].split("/") or "V" in r[9].split("/"))
                and float(r[15]) > 0}
        if want:
            tr = tracr_candidates(arrays, units, cds, arr_out, weighted, read_contigs(a.fasta, want), a.tracr_shuffles)
    with open(p + ".crispr_tracr.tsv", "w") as out:
        out.write("\t".join(["array_id", "contig", "start", "end", "strand", "length", "reads", "main_5p_end",
                             "share_5p", "main_3p_end", "share_3p", "position", "antirepeat_score", "aligned_nt",
                             "watson_crick", "wobble", "p_shuffle", "tracr_candidate", "sequence", "crRNA_repeat",
                             "align_candidate", "align_repeat_rc", "antirepeat_5p", "antirepeat_nt_before_reads"]) + "\n")
        for r in tr:
            out.write("\t".join(str(x) for x in r) + "\n")
    if a.fasta:
        sys.stderr.write("[crispr] tracrRNA: %d read clusters in type II/V loci, %d with an anti-repeat\n" % (
            len(tr), sum(1 for r in tr if r[17] == "yes")))
    sys.stderr.write("[crispr] %d arrays, %.1f reads on arrays, %d arrays with >= 10 reads\n"
                     % (len(arrays), sum(sum(r[3] for r in v) for v in reads.values()),
                        sum(1 for v in reads.values() if sum(r[3] for r in v) >= 10)))


if __name__ == "__main__":
    main()
