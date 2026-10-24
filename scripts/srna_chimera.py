#!/usr/bin/env python3
"""srna_chimera.py — chimeric reads (candidate RNA–RNA ligation products) in sRNA libraries.

Module `chimera` of sRNA-IMP. Exploratory: standard sRNA libraries are made without crosslinking and proximity
ligation (CLASH, RIL-seq), so most chimeric reads are library artefacts (template switching, intermolecular
ligation of abundant RNAs). The module therefore reports chimeras together with a test that real RNA–RNA
duplexes would pass and random joins would not: the two arms of real interaction partners should be more
complementary than arms paired at random.

Steps
  1. candidates: reads that did not map end-to-end to the sample assembly (unmapped records of the assembly BAM),
     >= 2 x --min-arm nt, not low-complexity; identical sequences collapsed
  2. seeds: the first and the last --min-arm nt of every candidate, mapped exactly (bowtie -v 0, <= --max-hits)
  3. arms: each seed extended along the contig (<= --arm-mismatches); a chimera = the two arms together cover the
     read (<= --max-gap nt unassigned; overlapping arms are allowed, the junction is then ambiguous)
  4. class: inter (different contigs), intra (same contig, not colinear), colinear (same contig and strand, 3' arm
     just downstream of the 5' arm: deletions / local template switches, reported but not used as interactions)
  5. annotation of each arm (assembly GFF as gff2bed + cmscan-rescued rRNA + novel sRNA clusters):
     rRNA, tRNA, ncRNA, CRISPR, CDS, 5UTR, 3UTR, novel_sRNA, intergenic; *_as = antisense to the feature
  6. pairs of features: reads, unique sequences, enrichment over the arms' chimera frequencies (one-sided
     hypergeometric / Fisher, as in RIL-seq), BH-FDR
  6b. organism: bin (contig2bin) and GTDB taxon of each arm; inter-contig chimeras whose arms come from the same
     bin, observed vs expected from the arms' bin frequencies (real sRNA–target pairs are mostly within one
     organism, random joins of abundant RNAs mostly between organisms)
  7. complementarity: hybridisation energy (ViennaRNA RNAduplex) of the two arms with --flank nt of contig sequence
     on both sides, against the same 5' arm paired with --null random 3' arms of other chimeras (empirical p), and
     per class of pair (rRNA–rRNA, mRNA–sRNA-like, ...) observed vs re-paired (Mann–Whitney)

Outputs (<out-prefix>.)
  chimera_summary.tsv     counts per step and class
  chimera_reads.tsv.gz    one row per chimeric sequence
  chimera_pairs.tsv       one row per pair of features (non-colinear chimeras)
  chimera_typepairs.tsv   reads per pair of arm classes
  chimera_energy.tsv      complementarity test per class of pair
"""
import argparse
import bisect
import collections
import gzip
import math
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile

COMP = str.maketrans("ACGTUNacgtun", "TGCAANtgcaan")
PRIORITY = ("rRNA", "tRNA", "ncRNA", "CRISPR", "CDS")
NCRNA_TYPES = {"ncRNA", "tmRNA", "regulatory_region", "riboswitch", "misc_RNA", "RNase_P_RNA", "SRP_RNA",
               "antisense_RNA", "sRNA"}
CRISPR_TYPES = {"CRISPR", "crispr-repeat", "crispr-spacer", "repeat_region"}
SRNA_LIKE = {"ncRNA", "novel_sRNA", "intergenic"}
MRNA_LIKE = {"5UTR", "CDS", "3UTR"}


def rc(s):
    return s.translate(COMP)[::-1]


def open_in(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def log(msg):
    print("[chimera] " + msg, file=sys.stderr, flush=True)


# ---------------------------------------------------------------- candidates
def low_complexity(s):
    c = collections.Counter(s)
    if max(c.values()) > 0.6 * len(s):
        return True
    return len({s[i:i + 3] for i in range(len(s) - 2)}) < 0.3 * (len(s) - 2)


def read_candidates(a):
    """unique sequences of unmapped reads (BAM via samtools, or a FASTQ for tests) -> Counter"""
    cnt = collections.Counter()
    stats = collections.Counter()
    if a.bam:
        p = subprocess.Popen(["samtools", "view", "-f", "4", "-@", str(max(1, a.threads - 1)), a.bam],
                             stdout=subprocess.PIPE, text=True)
        it = (line.split("\t", 10)[9] for line in p.stdout)
    else:
        def fq():
            with open_in(a.unmapped_fastq) as fh:
                for i, line in enumerate(fh):
                    if i % 4 == 1:
                        yield line.strip()
        p, it = None, fq()
    for s in it:
        stats["unmapped_reads"] += 1
        s = s.upper().replace("U", "T")
        if len(s) < 2 * a.min_arm:
            continue
        stats["unmapped_reads_long"] += 1
        if "N" in s or low_complexity(s):
            stats["low_complexity_or_N"] += 1
            continue
        cnt[s] += 1
    if p:
        p.wait()
        if p.returncode:
            sys.exit("samtools view failed")
    stats["candidate_reads"] = sum(cnt.values())
    stats["candidate_sequences"] = len(cnt)
    if a.max_sequences and len(cnt) > a.max_sequences:
        keep = dict(cnt.most_common(a.max_sequences))
        stats["candidate_sequences_dropped_cap"] = len(cnt) - len(keep)
        cnt = collections.Counter(keep)
    return cnt, stats


# ---------------------------------------------------------------- assembly
def read_fasta(path, keep=None):
    seqs, name, buf = {}, None, []
    with open_in(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None and (keep is None or name in keep):
                    seqs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None and (keep is None or name in keep):
        seqs[name] = "".join(buf).upper()
    return seqs


class Contigs:
    def __init__(self, seqs):
        self.f = seqs
        self._rc = {}

    def strand_seq(self, ref, strand):
        if strand == "+":
            return self.f[ref]
        if ref not in self._rc:
            self._rc[ref] = rc(self.f[ref])
        return self._rc[ref]

    def region(self, ref, start, end, strand, flank=0):
        L = len(self.f[ref])
        s, e = max(0, start - flank), min(L, end + flank)
        x = self.f[ref][s:e]
        return x if strand == "+" else rc(x)


# ---------------------------------------------------------------- seed mapping
def map_seeds_bowtie(seeds, a, tmp):
    """exact hits of the seed sequences: {seed: [(ref, strand, forward offset)]} (None = too many hits)"""
    fa = os.path.join(tmp, "seeds.fa")
    order = list(seeds)
    with open(fa, "w") as fh:
        for i, s in enumerate(order):
            fh.write(">%d\n%s\n" % (i, s))
    out = os.path.join(tmp, "seeds.hits")
    cmd = [a.bowtie, "-f", "-v", "0", "-k", str(a.max_hits + 1), "-p", str(a.threads), "--quiet", a.index, fa, out]
    r = subprocess.run(cmd, stderr=subprocess.PIPE, text=True)
    if r.returncode:
        sys.exit("bowtie failed: %s" % r.stderr[-2000:])
    hits = collections.defaultdict(list)
    with open(out) as fh:
        for line in fh:
            f = line.split("\t")
            hits[order[int(f[0])]].append((f[2], f[1], int(f[3])))
    return {s: (h if len(h) <= a.max_hits else None) for s, h in hits.items()}


def map_seeds_python(seeds, k, contigs, max_hits):
    """exact k-mer lookup in memory (tests / small assemblies only)"""
    want = set(seeds)
    hits = collections.defaultdict(list)
    for ref, s in contigs.f.items():
        for i in range(len(s) - k + 1):
            w = s[i:i + k]
            if w in want:
                hits[w].append((ref, "+", i))
            r = rc(w)
            if r in want:
                hits[r].append((ref, "-", i))
    return {s: (h if len(h) <= max_hits else None) for s, h in hits.items()}


# ---------------------------------------------------------------- arms
def extend5(read, k, hit, C, mm_max):
    """5' arm: seed = read[:k] at hit; extend 3'-ward. -> (len, ref, fstart, fend, strand)"""
    ref, strand, off = hit
    T = C.strand_seq(ref, strand)
    L = len(T)
    q = off if strand == "+" else L - off - k
    mm, last, i = 0, k - 1, k
    while i < len(read) and q + i < L:
        if read[i] != T[q + i]:
            mm += 1
            if mm > mm_max:
                break
        else:
            last = i
        i += 1
    n = last + 1
    s, e = q, q + n
    return (n, ref, s, e, strand) if strand == "+" else (n, ref, L - e, L - s, strand)


def extend3(read, k, hit, C, mm_max):
    """3' arm: seed = read[-k:] at hit; extend 5'-ward. -> (len, ref, fstart, fend, strand)"""
    ref, strand, off = hit
    T = C.strand_seq(ref, strand)
    L = len(T)
    r = off if strand == "+" else L - off - k
    n0 = len(read) - k
    mm, last, j = 0, 0, 1
    while n0 - j >= 0 and r - j >= 0:
        if read[n0 - j] != T[r - j]:
            mm += 1
            if mm > mm_max:
                break
        else:
            last = j
        j += 1
    n = k + last
    s, e = r - last, r + k
    return (n, ref, s, e, strand) if strand == "+" else (n, ref, L - e, L - s, strand)


def window_identity(read, T, start):
    """share of identical bases of read vs T[start:start+len(read)] (positions outside T count as mismatches)"""
    same = sum(1 for i, b in enumerate(read) if 0 <= start + i < len(T) and T[start + i] == b)
    return same / len(read)


def artefact_class(read, k, a5, a3, C, edge_tol, min_identity):
    """'contig_edge': the junction sits at the end of one contig and the start of the next (assembly break);
    'divergent': the whole read matches the contig of either arm with >= min_identity (a variant of the same
    sequence, e.g. a strain / redundant contig, not two RNAs). None otherwise."""
    _, r5, s5, e5, st5 = a5
    _, r3, s3, e3, st3 = a3
    L5, L3 = len(C.f[r5]), len(C.f[r3])
    end5 = (e5 >= L5 - edge_tol) if st5 == "+" else (s5 <= edge_tol)     # 3' end of the 5' arm at a contig end
    start3 = (s3 <= edge_tol) if st3 == "+" else (e3 >= L3 - edge_tol)   # 5' end of the 3' arm at a contig end
    if end5 and start3:
        return "contig_edge"
    n = len(read)
    T5 = C.strand_seq(r5, st5)
    q5 = s5 if st5 == "+" else L5 - e5                                   # read pos 0 on the 5' arm's strand
    T3 = C.strand_seq(r3, st3)
    q3 = (s3 if st3 == "+" else L3 - e3) - (n - a3[0])                   # read pos 0 if the 3' arm continued
    if window_identity(read, T5, q5) >= min_identity or window_identity(read, T3, q3) >= min_identity:
        return "divergent"
    return None


def classify(a5, a3, readlen, colinear_max):
    _, r5, s5, e5, st5 = a5
    _, r3, s3, e3, st3 = a3
    if r5 != r3:
        return "inter"
    if st5 == st3:
        d = s3 - e5 if st5 == "+" else s5 - e3      # distance along the transcript, 5' arm end -> 3' arm start
        if -readlen <= d <= colinear_max:
            return "colinear"
    return "intra"


# ---------------------------------------------------------------- annotation
def attr(a, keys):
    for k in keys:
        m = re.search(r"(?:^|;)%s=([^;]+)" % k, a)
        if m:
            return m.group(1)
    return ""


class Annotation:
    def __init__(self, gbeds, novel_bed, refs, utr5, utr3):
        self.utr5, self.utr3 = utr5, utr3
        feats = collections.defaultdict(list)
        for path in gbeds:
            if not path or not os.path.exists(path):
                continue
            with open_in(path) as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) < 8 or f[0] not in refs:
                        continue
                    t = f[7]
                    if t in ("rRNA",):
                        c = "rRNA"
                    elif t == "tRNA":
                        c = "tRNA"
                    elif t in NCRNA_TYPES:
                        c = "ncRNA"
                    elif t in CRISPR_TYPES:
                        c = "CRISPR"
                    elif t == "CDS":
                        c = "CDS"
                    else:
                        continue
                    at = f[9] if len(f) > 9 else ""
                    name = attr(at, ("product", "Name", "gene")) or t
                    ko = attr(at, ("kegg_ko", "KEGG", "ko"))
                    feats[f[0]].append((int(f[1]), int(f[2]), f[5], c, f[3], name, ko))
        if novel_bed and os.path.exists(novel_bed):
            with open_in(novel_bed) as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) >= 6 and f[0] in refs:
                        feats[f[0]].append((int(f[1]), int(f[2]), f[5], "novel_sRNA", f[3], "novel sRNA locus", ""))
        self.f, self.starts, self.maxlen = {}, {}, {}
        for ref, L in feats.items():
            L.sort()
            self.f[ref] = L
            self.starts[ref] = [x[0] for x in L]
            self.maxlen[ref] = max(x[1] - x[0] for x in L)

    def near(self, ref, s, e):
        if ref not in self.f:
            return []
        st, L = self.starts[ref], self.f[ref]
        j = bisect.bisect_left(st, e + max(self.utr5, self.utr3))
        lo = s - self.maxlen[ref] - max(self.utr5, self.utr3)
        out = []
        while j > 0:
            j -= 1
            if st[j] < lo:
                break
            out.append(L[j])
        return out

    def annotate(self, ref, s, e, strand):
        """-> (class, feature id, name, ko)"""
        alen = e - s
        cand = self.near(ref, s, e)
        best = None
        for (fs, fe, fst, c, fid, name, ko) in cand:
            ov = min(e, fe) - max(s, fs)
            if ov < 0.5 * alen or c == "novel_sRNA":
                continue
            sense = fst == strand or fst == "."
            rank = PRIORITY.index(c) * 2 + (0 if sense else 1)
            if best is None or rank < best[0]:
                best = (rank, c if sense else c + "_as", fid, name, ko)
        if best and not best[1].startswith("CDS_as"):
            return best[1:]
        for (fs, fe, fst, c, fid, name, ko) in cand:          # novel sRNA locus on the same strand
            if c == "novel_sRNA" and fst == strand and min(e, fe) - max(s, fs) >= 0.5 * alen:
                return ("novel_sRNA", fid, name, ko)
        if best:
            return best[1:]
        mid = (s + e) // 2
        for (fs, fe, fst, c, fid, name, ko) in cand:          # UTRs of CDS on the same strand
            if c != "CDS" or fst != strand:
                continue
            if strand == "+":
                if fs - self.utr5 <= mid < fs:
                    return ("5UTR", "5UTR|" + fid, name, ko)
                if fe <= mid < fe + self.utr3:
                    return ("3UTR", "3UTR|" + fid, name, ko)
            else:
                if fe <= mid < fe + self.utr5:
                    return ("5UTR", "5UTR|" + fid, name, ko)
                if fs - self.utr3 <= mid < fs:
                    return ("3UTR", "3UTR|" + fid, name, ko)
        return ("intergenic", "ig|%s|%s|%d" % (ref, strand, mid // 200 * 200), "", "")


# ---------------------------------------------------------------- bins / taxonomy
def strip_ext(name):
    return re.sub(r"\.(fa|fna|fasta|fsa)(\.gz)?$", "", name)


def read_contig2bin(path):
    d = {}
    if path and os.path.exists(path):
        with open_in(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 2 and not line.startswith("#"):
                    d[f[0]] = strip_ext(f[1])
    return d


def tax_short(cl):
    """deepest named GTDB rank, e.g. 's__Faecalibacterium prausnitzii' -> 'Faecalibacterium prausnitzii'"""
    for r in reversed([x for x in cl.split(";") if x]):
        if len(r) > 3:
            return r[3:]
    return ""


def read_gtdbtk(paths):
    tax = {}
    for pth in [p for p in (paths or "").split(",") if p and os.path.exists(p)]:
        with open_in(pth) as fh:
            header = None
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if header is None or f[0] in ("user_genome", "Bin", "bin"):
                    header = f
                    continue
                rec = dict(zip(header, f))
                key = rec.get("user_genome") or rec.get("Bin") or rec.get("bin") or f[0]
                tax[strip_ext(key)] = tax_short(rec.get("classification", ""))
    return tax


# ---------------------------------------------------------------- energies
def duplex_energies(pairs, a):
    """minimum free energy of the intermolecular duplex for each (seq1, seq2); RNA python module or RNAduplex"""
    if not pairs:
        return []
    try:
        import RNA
        return [RNA.duplexfold(x, y).energy for x, y in pairs]
    except ImportError:
        pass
    exe = a.rnaduplex
    if not exe or not shutil.which(exe):
        return None
    n = max(1, a.threads)
    chunks = [pairs[i::n] for i in range(n)]
    procs = []
    for ch in chunks:
        p = subprocess.Popen([exe, "--noconv"] if a.rnaduplex_noconv else [exe], stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        procs.append((p, ch))
    res = [None] * n
    rx = re.compile(r"\(\s*([-+]?\d+(?:\.\d+)?)\)\s*$")
    for i, (p, ch) in enumerate(procs):
        out, err = p.communicate("".join("%s\n%s\n" % (x, y) for x, y in ch) + "@\n")
        if p.returncode:
            sys.exit("RNAduplex failed: " + err[-1000:])
        e = [float(m.group(1)) for m in (rx.search(l) for l in out.splitlines()) if m]
        if len(e) != len(ch):
            sys.exit("RNAduplex: %d energies for %d pairs" % (len(e), len(ch)))
        res[i] = [min(0.0, x) for x in e]
    out = [None] * len(pairs)
    for i in range(n):
        out[i::n] = res[i]
    return out


def dinuc_shuffle(s, rng):
    """random sequence with exactly the dinucleotide counts of s (Altschul & Erickson 1985, Euler path)"""
    if len(s) < 4:
        return s
    edges = collections.defaultdict(list)
    for x, y in zip(s, s[1:]):
        edges[x].append(y)
    last = s[-1]
    for _ in range(1000):
        lastedge = {v: rng.choice(L) for v, L in edges.items() if v != last}
        ok = True
        for v in lastedge:                     # every vertex must reach the last character via its last edge
            seen, u = set(), v
            while u != last:
                if u in seen or u not in lastedge:
                    ok = False
                    break
                seen.add(u)
                u = lastedge[u]
            if not ok:
                break
        if ok:
            break
    else:
        return s
    lists = {}
    for v, L in edges.items():
        L = L[:]
        if v in lastedge:
            L.remove(lastedge[v])
            rng.shuffle(L)
            L.append(lastedge[v])
        else:
            rng.shuffle(L)
        lists[v] = L
    out, pos, u = [s[0]], collections.Counter(), s[0]
    for _ in range(len(s) - 1):
        u2 = lists[u][pos[u]]
        pos[u] += 1
        out.append(u2)
        u = u2
    return "".join(out)


# ---------------------------------------------------------------- statistics
def hypergeom_sf(a, nA, nB, N):
    """P(X >= a), X ~ hypergeometric: N reads, nA with feature A, nB drawn (reads with B)"""
    if a <= 0:
        return 1.0
    lo, hi = max(0, nA + nB - N), min(nA, nB)
    if a > hi:
        return 0.0

    def lp(x):
        return (math.lgamma(nA + 1) - math.lgamma(x + 1) - math.lgamma(nA - x + 1) + math.lgamma(N - nA + 1)
                - math.lgamma(nB - x + 1) - math.lgamma(N - nA - nB + x + 1) - math.lgamma(N + 1)
                + math.lgamma(nB + 1) + math.lgamma(N - nB + 1))
    tot, x = 0.0, max(a, lo)
    while x <= hi:
        t = math.exp(lp(x))
        tot += t
        if t < 1e-14 * max(tot, 1e-300):
            break
        x += 1
    return min(1.0, tot)


def bh(ps):
    idx = [i for i, p in enumerate(ps) if p is not None]
    q = [None] * len(ps)
    m = len(idx)
    prev = 1.0
    for rank, i in sorted(enumerate(sorted(idx, key=lambda i: ps[i]), 1), key=lambda x: -x[0]):
        prev = min(prev, ps[i] * m / rank)
        q[i] = prev
    return q


def mannwhitney_less(x, y):
    """one-sided p that x tends to be smaller than y (normal approximation, ties averaged)"""
    if not x or not y:
        return float("nan")
    allv = sorted([(v, 0) for v in x] + [(v, 1) for v in y])
    ranks, i = {}, 0
    rx = 0.0
    while i < len(allv):
        j = i
        while j < len(allv) and allv[j][0] == allv[i][0]:
            j += 1
        r = (i + j + 1) / 2.0
        rx += r * sum(1 for t in allv[i:j] if t[1] == 0)
        i = j
    n1, n2 = len(x), len(y)
    u = rx - n1 * (n1 + 1) / 2.0
    mu, sd = n1 * n2 / 2.0, math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12.0)
    z = (u - mu) / sd if sd else 0.0
    return 0.5 * math.erfc(-z / math.sqrt(2))


def median(v):
    v = sorted(v)
    if not v:
        return float("nan")
    m = len(v) // 2
    return v[m] if len(v) % 2 else 0.5 * (v[m - 1] + v[m])


def pair_group(ta, tb):
    t = {ta.replace("_as", ""), tb.replace("_as", "")}
    if t == {"rRNA"}:
        return "rRNA–rRNA"
    if "rRNA" in t:
        return "rRNA–other"
    if "tRNA" in t:
        return "tRNA–other"
    if "CRISPR" in t:
        return "CRISPR–other"
    if (ta in SRNA_LIKE and tb in MRNA_LIKE) or (tb in SRNA_LIKE and ta in MRNA_LIKE):
        return "sRNA-like–mRNA"
    if ta in MRNA_LIKE and tb in MRNA_LIKE:
        return "mRNA–mRNA"
    return "other"


def same_bin(x, y):
    return "NA" if "unbinned" in (x, y) else ("yes" if x == y else "no")


def fmt(x):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "NA"
    if isinstance(x, float):
        return "%.4g" % x
    return str(x)


USED = {"inter", "intra", "pair_inter"}       # classes used as interactions (pairs, statistics, energies)
INTER = {"inter", "pair_inter"}               # different contigs: same-bin test
SUMMARY_KEYS = ("unmapped_reads", "unmapped_reads_long", "low_complexity_or_N", "candidate_reads",
                "candidate_sequences", "candidate_sequences_dropped_cap", "reads_seed_unmapped_or_repetitive",
                "reads_arms_do_not_cover", "chimeric_reads_inter", "chimeric_reads_intra", "chimeric_reads_colinear",
                "chimeric_sequences_inter", "chimeric_sequences_intra", "chimeric_sequences_colinear",
                "chimeric_reads_contig_edge", "chimeric_reads_divergent")


def write_outputs(rows, C, org, a, stats, summary_keys):
    """rows: (key / sequence, reads, class, arm5, arm3, annotation5, annotation3, hits5, hits3, overlap);
    arm = (length, contig, start, end, RNA strand) in RNA 5'->3' order. Writes all <out-prefix>.chimera_* tables."""
    nch = sum(r[1] for r in rows if r[2] in USED)

    with gzip.open(a.out_prefix + ".chimera_reads.tsv.gz", "wt") as fh:
        fh.write("sequence\treads\tlength\tclass\tcontig_5\tstart_5\tend_5\tstrand_5\tlen_5\ttype_5\tfeature_5\tname_5"
                 "\tko_5\tcontig_3\tstart_3\tend_3\tstrand_3\tlen_3\ttype_3\tfeature_3\tname_3\tko_3\thits_5\thits_3"
                 "\tjunction_overlap\tbin_5\ttaxon_5\tbin_3\ttaxon_3\n")
        for s, c, cl, a5, a3, an5, an3, n5, n3, ov in sorted(rows, key=lambda r: -r[1]):
            fh.write("\t".join(map(str, [s, c, len(s), cl, a5[1], a5[2], a5[3], a5[4], a5[0], an5[0], an5[1],
                                         an5[2], an5[3], a3[1], a3[2], a3[3], a3[4], a3[0], an3[0], an3[1], an3[2],
                                         an3[3], n5, n3, ov] + list(org(a5[1])) + list(org(a3[1])))) + "\n")

    # type pairs
    tp = collections.Counter()
    tpu = collections.Counter()
    for s, c, cl, a5, a3, an5, an3, *_ in rows:
        key = tuple(sorted((an5[0], an3[0]))) + (cl,)
        tp[key] += c
        tpu[key] += 1
    with open(a.out_prefix + ".chimera_typepairs.tsv", "w") as fh:
        fh.write("type_a\ttype_b\tclass\treads\tsequences\tgroup\n")
        for key, v in tp.most_common():
            fh.write("%s\t%s\t%s\t%d\t%d\t%s\n" % (key[0], key[1], key[2], v, tpu[key], pair_group(key[0], key[1])))

    # 6. feature pairs (classes in USED)
    pr = {}
    per_feat = collections.Counter()
    for s, c, cl, a5, a3, an5, an3, *_ in rows:
        if cl not in USED:
            continue
        fa, fb = an5[1], an3[1]
        key = tuple(sorted((fa, fb)))
        d = pr.get(key)
        if d is None:
            ia, ib = (an5, an3) if fa <= fb else (an3, an5)
            ra, rb = (a5[1], a3[1]) if fa <= fb else (a3[1], a5[1])
            d = pr[key] = dict(a=ia, b=ib, cls=collections.Counter(), reads=0, seqs=0, a_first=0, rep=None,
                               org_a=org(ra), org_b=org(rb))
        d["reads"] += c
        d["seqs"] += 1
        d["cls"][cl] += c
        if fa == key[0]:
            d["a_first"] += c
        # representative for the duplex test: a read with both arms (split / ligated read) before a
        # discordant pair (mates are not the pairing regions), then the most reads
        if d["rep"] is None or (cl.startswith("pair"), -c) < (d["rep"][4].startswith("pair"), -d["rep"][1]):
            d["rep"] = (s, c, a5, a3, cl)
        per_feat[fa] += c
        if fb != fa:
            per_feat[fb] += c
    keys = list(pr)
    pv = []
    for key in keys:
        d = pr[key]
        if key[0] == key[1] or d["reads"] < 2:
            pv.append(None)
        else:
            pv.append(hypergeom_sf(d["reads"], per_feat[key[0]], per_feat[key[1]], nch))
    qv = bh(pv)

    # 7. complementarity of the arms
    # two nulls per pair: the 3' arm dinucleotide-shuffled (same composition, so GC content / organism do not
    # count as complementarity; main test) and the 5' arm with 3' arms of other chimeras (re-paired)
    test = sorted([kk for kk in keys if not pr[kk]["rep"][4].startswith("pair")],
                  key=lambda kk: -pr[kk]["reads"])[:a.max_pairs]
    arms5, arms3 = [], []
    for kk in test:
        s, c, a5, a3, cl = pr[kk]["rep"]
        arms5.append(C.region(a5[1], a5[2], a5[3], a5[4], a.flank))
        arms3.append(C.region(a3[1], a3[2], a3[3], a3[4], a.flank))
    pairs = list(zip(arms5, arms3))
    rng = random.Random(a.seed)
    nshuf, nrep = [], []
    for i in range(len(test)):
        for _ in range(a.null):
            nshuf.append((arms5[i], dinuc_shuffle(arms3[i], rng)))
        if len(test) > 1:
            for _ in range(a.null):
                j = rng.randrange(len(test) - 1)
                j += j >= i
                nrep.append((arms5[i], arms3[j]))
    energies = duplex_energies(pairs + nshuf + nrep, a)
    eobs, enull, erep = {}, {}, {}
    if energies:
        obs = energies[:len(pairs)]
        sh = energies[len(pairs):len(pairs) + len(nshuf)]
        rp = energies[len(pairs) + len(nshuf):]
        for i, kk in enumerate(test):
            eobs[kk] = obs[i]
            enull[kk] = sh[i * a.null:(i + 1) * a.null]
            erep[kk] = rp[i * a.null:(i + 1) * a.null] if rp else []
    else:
        log("no duplex energies (RNAduplex not found and no RNA python module)")

    with open(a.out_prefix + ".chimera_pairs.tsv", "w") as fh:
        fh.write("feature_a\tfeature_b\ttype_a\ttype_b\tgroup\tname_a\tname_b\tko_a\tko_b\tclass\treads\tsequences"
                 "\treads_a_5prime\treads_a_all_chimeras\treads_b_all_chimeras\tp\tq\tdG\tdG_null_median"
                 "\tp_complementarity\tdG_repaired_median\tp_vs_repaired\tbin_a\ttaxon_a\tbin_b\ttaxon_b\tsame_bin"
                 "\trep_sequence\n")
        order = sorted(range(len(keys)), key=lambda i: (pv[i] is None, pv[i] if pv[i] is not None else 0,
                                                        -pr[keys[i]]["reads"]))
        for i in order:
            kk, d = keys[i], pr[keys[i]]
            e = eobs.get(kk)
            nl = enull.get(kk, [])
            pc = (1 + sum(1 for x in nl if x <= e)) / (1 + len(nl)) if e is not None and nl else None
            nr = erep.get(kk, [])
            pr_ = (1 + sum(1 for x in nr if x <= e)) / (1 + len(nr)) if e is not None and nr else None
            fh.write("\t".join(fmt(x) for x in [
                kk[0], kk[1], d["a"][0], d["b"][0], pair_group(d["a"][0], d["b"][0]), d["a"][2], d["b"][2],
                d["a"][3], d["b"][3], "/".join(sorted(d["cls"])), d["reads"], d["seqs"], d["a_first"],
                per_feat[kk[0]], per_feat[kk[1]], pv[i], qv[i], e, median(nl) if nl else None, pc,
                median(nr) if nr else None, pr_,
                d["org_a"][0], d["org_a"][1], d["org_b"][0], d["org_b"][1], same_bin(d["org_a"][0], d["org_b"][0]),
                d["rep"][0]]) + "\n")

    # energy summary per group of pairs
    with open(a.out_prefix + ".chimera_energy.tsv", "w") as fh:
        fh.write("group\tpairs_tested\treads\tdG_median\tdG_null_median\tfraction_p_complementarity_lt_0.05"
                 "\tp_observed_lower_than_null\tdG_repaired_median\tp_observed_lower_than_repaired\n")
        groups = collections.defaultdict(list)
        for kk in eobs:
            d = pr[kk]
            groups[pair_group(d["a"][0], d["b"][0])].append(kk)
        groups["all"] = list(eobs)
        for g in sorted(groups, key=lambda x: (x == "all", x)):
            ks = groups[g]
            ob = [eobs[x] for x in ks]
            nl = [y for x in ks for y in enull[x]]
            pcs = [(1 + sum(1 for y in enull[x] if y <= eobs[x])) / (1 + len(enull[x])) for x in ks if enull[x]]
            nr = [y for x in ks for y in erep.get(x, [])]
            fh.write("\t".join(fmt(x) for x in [
                g, len(ks), sum(pr[x]["reads"] for x in ks), median(ob), median(nl),
                sum(1 for p in pcs if p < 0.05) / len(pcs) if pcs else None, mannwhitney_less(ob, nl),
                median(nr) if nr else None, mannwhitney_less(ob, nr) if nr else None]) + "\n")

    # 6b. same organism: inter-contig chimeras with both arms binned, observed vs expected same-bin fraction
    b5, b3, same = collections.Counter(), collections.Counter(), 0
    grp_same = collections.defaultdict(lambda: [0, 0])
    for s, c, cl, a5, a3, an5, an3, *_ in rows:
        if cl not in INTER:
            continue
        x, y = org(a5[1])[0], org(a3[1])[0]
        if x == "unbinned" or y == "unbinned":
            continue
        b5[x] += c
        b3[y] += c
        same += c * (x == y)
        g = grp_same[pair_group(an5[0], an3[0])]
        g[0] += c
        g[1] += c * (x == y)
    nb = sum(b5.values())
    exp = sum(b5[b] * b3[b] for b in b5) / nb ** 2 if nb else None
    stats["inter_reads_both_binned"] = nb
    stats["inter_reads_same_bin"] = same
    with open(a.out_prefix + ".chimera_bins.tsv", "w") as fh:
        fh.write("group\tinter_reads_both_binned\tsame_bin_reads\tsame_bin_fraction\texpected_fraction\n")
        for g in sorted(grp_same, key=lambda x: -grp_same[x][0]) + ["all"]:
            n_, k_ = (nb, same) if g == "all" else grp_same[g]
            fh.write("\t".join(fmt(x) for x in [g, n_, k_, k_ / n_ if n_ else None, exp]) + "\n")

    stats["chimeric_reads_used"] = nch
    stats["feature_pairs"] = len(keys)
    stats["feature_pairs_q_lt_0.05"] = sum(1 for q in qv if q is not None and q < 0.05)
    stats["feature_pairs_tested_complementarity"] = len(eobs)
    with open(a.out_prefix + ".chimera_summary.tsv", "w") as fh:
        fh.write("metric\tvalue\n")
        for kk in summary_keys + ("chimeric_reads_used", "inter_reads_both_binned", "inter_reads_same_bin",
                                  "feature_pairs", "feature_pairs_q_lt_0.05", "feature_pairs_tested_complementarity"):
            fh.write("%s\t%d\n" % (kk, stats.get(kk, 0)))
        for kk in ("min_arm", "max_hits", "arm_mismatches", "max_gap", "flank", "flank_pair", "null", "strandedness"):
            if hasattr(a, kk):
                fh.write("param_%s\t%s\n" % (kk, getattr(a, kk)))
    log("chimeras: %d reads used (%s) in %d feature pairs, %d with q < 0.05; complementarity tested for %d pairs" % (
        nch, ", ".join("%s %d" % (x, stats.get("chimeric_reads_" + x, 0)) for x in sorted(USED)
                       if stats.get("chimeric_reads_" + x)), len(keys), stats["feature_pairs_q_lt_0.05"], len(eobs)))



# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--bam", help="assembly BAM of the library (its unmapped records are the candidates)")
    src.add_argument("--unmapped-fastq", help="unmapped reads as FASTQ (tests)")
    ap.add_argument("--assembly", required=True, help="assembly FASTA(.gz)")
    ap.add_argument("--index", help="bowtie (1) index prefix of the assembly")
    ap.add_argument("--seed-mapper", choices=("bowtie", "python"), default="bowtie",
                    help="python = in-memory exact k-mer lookup (tests / small assemblies)")
    ap.add_argument("--gbed", action="append", default=[], help="annotation as gff2bed (repeatable)")
    ap.add_argument("--novel-bed", help="novel sRNA loci (BED6) of the same library")
    ap.add_argument("--contig2bin", help="contig<TAB>bin table")
    ap.add_argument("--gtdbtk", help="GTDB-Tk summary / bin stats with 'classification' (comma-separated)")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--min-arm", type=int, default=18, help="seed = minimum arm length (nt) [18]")
    ap.add_argument("--max-hits", type=int, default=5, help="seeds with more exact hits are not used [5]")
    ap.add_argument("--arm-mismatches", type=int, default=1, help="mismatches allowed while extending an arm [1]")
    ap.add_argument("--max-gap", type=int, default=2, help="read nt covered by neither arm [2]")
    ap.add_argument("--edge-tol", type=int, default=3, help="junction within this of two contig ends: contig_edge [3]")
    ap.add_argument("--min-identity", type=float, default=0.85, help="whole read >= this identical to the contig "
                    "of one arm: divergent (variant, not chimera) [0.85]")
    ap.add_argument("--colinear-max", type=int, default=1000,
                    help="same contig and strand, 3' arm <= this many nt downstream: colinear [1000]")
    ap.add_argument("--utr5", type=int, default=150, help="5' UTR window upstream of a CDS [150]")
    ap.add_argument("--utr3", type=int, default=100, help="3' UTR window downstream of a CDS [100]")
    ap.add_argument("--flank", type=int, default=20, help="contig nt added on both sides of an arm for the "
                    "duplex energy [20]")
    ap.add_argument("--null", type=int, default=20, help="random re-pairings per tested pair [20]")
    ap.add_argument("--max-pairs", type=int, default=5000, help="feature pairs tested for complementarity "
                    "(most reads first) [5000]")
    ap.add_argument("--max-sequences", type=int, default=3000000, help="cap of unique candidate sequences")
    ap.add_argument("--bowtie", default="bowtie")
    ap.add_argument("--rnaduplex", default="RNAduplex", help="RNAduplex executable (used if the RNA python "
                    "module is not available); '' = no energies")
    ap.add_argument("--rnaduplex-noconv", action="store_true", help=argparse.SUPPRESS)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    if a.seed_mapper == "bowtie" and not a.index:
        ap.error("--index is required with --seed-mapper bowtie")
    random.seed(a.seed)
    k = a.min_arm
    os.makedirs(os.path.dirname(a.out_prefix) or ".", exist_ok=True)

    # 1. candidates
    cand, stats = read_candidates(a)
    log("%d unmapped reads, %d >= %d nt, %d candidate reads (%d sequences)" % (
        stats["unmapped_reads"], stats["unmapped_reads_long"], 2 * k, stats["candidate_reads"], len(cand)))

    # 2. seeds
    seeds = set()
    for s in cand:
        seeds.add(s[:k])
        seeds.add(s[-k:])
    tmp = tempfile.mkdtemp(prefix="chimera.", dir=os.path.dirname(a.out_prefix) or ".")
    try:
        if a.seed_mapper == "bowtie":
            hits = map_seeds_bowtie(seeds, a, tmp)
            refs = {h[0] for v in hits.values() if v for h in v}
            C = Contigs(read_fasta(a.assembly, refs))
        else:
            C = Contigs(read_fasta(a.assembly))
            hits = map_seeds_python(seeds, k, C, a.max_hits)
            refs = set(C.f)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    log("%d seeds, %d with 1..%d exact hits" % (len(seeds), sum(1 for v in hits.values() if v), a.max_hits))

    # 3./4. arms and classes
    A = Annotation(a.gbed, a.novel_bed, refs, a.utr5, a.utr3)
    c2b = read_contig2bin(a.contig2bin)
    b2t = read_gtdbtk(a.gtdbtk)

    def org(ref):
        b = c2b.get(ref, "unbinned")
        return b, b2t.get(b, "")
    rows = []
    for s, c in cand.items():
        h5, h3 = hits.get(s[:k]), hits.get(s[-k:])
        if h5 is None or h3 is None:
            stats["reads_seed_unmapped_or_repetitive"] += c
            continue
        n = len(s)
        best = None
        for x in h5:
            a5 = extend5(s, k, x, C, a.arm_mismatches)
            for y in h3:
                a3 = extend3(s, k, y, C, a.arm_mismatches)
                if n - (a5[0] + a3[0]) > a.max_gap:
                    continue
                cl = classify(a5, a3, n, a.colinear_max)
                key = (cl != "colinear", -(a5[0] + a3[0]))   # a colinear explanation wins (conservative), then coverage
                if best is None or key < best[0]:
                    best = (key, a5, a3, cl)
        if best is None:
            stats["reads_arms_do_not_cover"] += c
            continue
        _, a5, a3, cl = best
        if cl != "colinear":
            cl = artefact_class(s, k, a5, a3, C, a.edge_tol, a.min_identity) or cl
        stats["chimeric_reads_" + cl] += c
        stats["chimeric_sequences_" + cl] += 1
        an5 = A.annotate(a5[1], a5[2], a5[3], a5[4])
        an3 = A.annotate(a3[1], a3[2], a3[3], a3[4])
        rows.append((s, c, cl, a5, a3, an5, an3, len(h5), len(h3), max(0, a5[0] + a3[0] - n)))
    log("chimeric reads: %d inter, %d intra, %d colinear" % (stats["chimeric_reads_inter"],
                                                              stats["chimeric_reads_intra"],
                                                              stats["chimeric_reads_colinear"]))
    write_outputs(rows, C, org, a, stats, SUMMARY_KEYS)


if __name__ == "__main__":
    main()
