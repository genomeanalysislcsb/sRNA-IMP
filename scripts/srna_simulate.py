#!/usr/bin/env python3
"""srna_simulate.py — simulate a stool-like small RNA library with known origin of every read (benchmark).

The microbial reads come from a real sample assembly (FASTA + GFF annotation, optionally bins and metaG contig
counts for abundances); human reads from human small RNA reference FASTAs. The pipeline is then run on the
simulated library with a reduced assembly / annotation, and srna_sim_eval.py scores its outputs against the truth.

What is simulated
  microbial, assembled     fragments of annotated features (tRNA, rRNA, tmRNA, ncRNA, regulatory elements, CRISPR,
                           CDS sense and antisense) and of intergenic stretches, from contigs kept in the assembly
  microbial, unassembled   the same from contigs held out of the assembly given to the pipeline (--holdout-frac):
                           reads with no contig, as in real data
  hidden loci              ncRNA / regulatory features removed from the GFF given to the pipeline (--hide-frac):
                           ground truth for the novel-locus module
  human                    mature miRNAs with 3' isomiRs, tRNA halves and fragments, rRNA and other ncRNA
                           fragments (snRNA, Y RNA, ...), at --host-frac of the reads
  adapter dimers           0-5 nt inserts (--dimer-frac)
  reads                    insert + 3' adapter + random fill to --read-len, position-dependent substitution errors

Inputs
  --assembly    contigs (FASTA, .gz ok), e.g. IMP Assembly/mgmt.assembly.merged.fa.gz
  --gff         annotation of the contigs (GFF3, .gz ok), e.g. Analysis/annotation/annotation_CDS_RNA_hmms.gff.gz
  --contig2bin  optional contig<TAB>bin (.gz ok): contigs of one bin share an abundance
  --contig-counts optional IMP Stats/mg/*.contig_read_count.txt.gz (contig start end reads): abundance = reads/length
  --human-mirna / --human-trna / --human-rrna / --human-other   human reference FASTAs (mature miRNAs, tRNAs,
                rRNAs, other ncRNAs). From the pipeline's bowtie indices: bowtie-inspect -a 0 $IDX_HSA_MIRNA > hsa_mirna.fa

Outputs (--out PREFIX)
  PREFIX.fastq.gz          the simulated library
  PREFIX.truth.tsv.gz      one row per read: origin, class, feature, contig, bin, strand, coordinates, fragment type
  PREFIX.assembly.fa.gz    the assembly given to the pipeline (held-out contigs removed)
  PREFIX.gff.gz            the annotation given to the pipeline (held-out contigs and hidden features removed)
  PREFIX.contig2bin.tsv    bins of the kept contigs (if --contig2bin)
  PREFIX.hidden_loci.tsv   features removed from the GFF (truth for the novel module)
  PREFIX.params.tsv        all parameters and the seed
  PREFIX.sample_row.tsv    a sample-sheet row for srna_pipeline.sh (absolute paths)

Example
  srna_simulate.py --assembly mgmt.assembly.merged.fa.gz --gff annotation_CDS_RNA_hmms.gff.gz \\
      --contig2bin bins.contigs.gz --contig-counts mg.assembly.contig_read_count.txt.gz \\
      --human-mirna hsa_mirna.fa --human-trna hsa_trna.fa --human-rrna hsa_rrna.fa --human-other hsa_other.fa \\
      --n-reads 5000000 --host-frac 0.02 --seed 1 --sample SIM01 --out sim/SIM01
"""
import argparse
import collections
import gzip
import os
import random
import sys

ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

# share of the microbial reads per class (normalised over the classes present in the annotation)
MICROBIAL_MIX = {"rRNA": 0.42, "tRNA": 0.26, "tmRNA": 0.02, "ncRNA": 0.05, "regulatory": 0.03, "CRISPR": 0.01,
                 "CDS": 0.11, "CDS_antisense": 0.02, "intergenic": 0.08}
HUMAN_MIX = {"miRNA": 0.35, "tRNA": 0.20, "rRNA": 0.15, "other": 0.30}
GFF_CLASS = {"rRNA": "rRNA", "tRNA": "tRNA", "tmRNA": "tmRNA", "ncRNA": "ncRNA", "regulatory_region": "regulatory",
             "riboswitch": "regulatory", "CRISPR": "CRISPR", "crispr-array": "CRISPR", "repeat_region": None,
             "CDS": "CDS"}
HIDEABLE = {"ncRNA", "regulatory"}


def rc(s):
    return s.translate(COMP)[::-1]


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with opener(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line)
    if name is not None:
        seqs[name] = "".join(buf).upper()
    return seqs


def read_gff(path, contigs):
    feats = []
    with opener(path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[0] not in contigs:
                continue
            cls = GFF_CLASS.get(f[2])
            if not cls:
                continue
            attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            fid = attrs.get("ID") or attrs.get("locus_tag") or "%s:%s-%s" % (f[0], f[3], f[4])
            s, e = int(f[3]) - 1, int(f[4])
            if e - s < 15:
                continue
            feats.append(dict(id=fid, contig=f[0], start=s, end=e, strand=f[6] if f[6] in "+-" else "+",
                              cls=cls, type=f[2], line=line))
    return feats


def read_table(path, cols=2):
    out = {}
    if path:
        with opener(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= cols and not line.startswith("#"):
                    out[f[0]] = f[1:]
    return out


def abundances(contigs, c2b, counts, rng):
    """relative abundance per contig: metaG reads per bp if given, else log-normal per bin (shared within a bin)"""
    w = {}
    if counts:
        for c, seq in contigs.items():
            r = counts.get(c)
            try:
                w[c] = float(r[-1]) / max(len(seq), 1) if r else 0.0
            except ValueError:
                w[c] = 0.0
        if sum(w.values()) > 0:
            return w
    binw = {}
    for c in contigs:
        b = c2b.get(c, [None])[0]
        if b and b != "unbinned":
            binw.setdefault(b, rng.lognormvariate(0, 1.5))
            w[c] = binw[b]
        else:
            w[c] = rng.lognormvariate(-0.5, 1.5)
    return w


class Picker:
    """weighted random choice over a fixed list (cumulative weights + bisection)"""

    def __init__(self, items, weights, rng):
        import bisect
        self.items, self.rng, self.bisect = items, rng, bisect.bisect
        tot, self.cum = 0.0, []
        for x in weights:
            tot += max(x, 0.0)
            self.cum.append(tot)
        self.tot = tot

    def __bool__(self):
        return bool(self.items) and self.tot > 0

    def pick(self):
        return self.items[min(self.bisect(self.cum, self.rng.random() * self.tot), len(self.items) - 1)]


def clamp(x, lo, hi):
    return max(lo, min(hi, x))


def frag_len(rng, mean=28, sd=8, lo=16, hi=60):
    return int(clamp(round(rng.gauss(mean, sd)), lo, hi))


def trna_fragment(rna, rng):
    """tRNA fragment types with their usual sizes; the mature 3' end carries CCA"""
    if not rna.endswith("CCA"):
        rna += "CCA"
    n = len(rna)
    kind = rng.choices(["5half", "3half", "3tRF", "5tRF", "internal", "mature"], [0.30, 0.20, 0.15, 0.10, 0.15, 0.10])[0]
    if kind == "5half":
        L = clamp(int(rng.gauss(33, 2)), 28, 38); s = 0
    elif kind == "3half":
        L = clamp(int(rng.gauss(38, 2)), 32, 44); s = max(0, n - L)
    elif kind == "3tRF":
        L = rng.randint(18, 22); s = max(0, n - L)
    elif kind == "5tRF":
        L = rng.randint(15, 30); s = 0
    elif kind == "internal":
        L = rng.randint(18, 40); s = rng.randint(0, max(0, n - L))
    else:
        L = n; s = 0
    L = min(L, n - s)
    return rna[s:s + L], s, s + L, kind


def generic_fragment(rna, rng, five_prime=0.0, mean=28, sd=8):
    n = len(rna)
    L = min(frag_len(rng, mean, sd), n)
    if five_prime and rng.random() < five_prime:
        return rna[:L], 0, L, "5end"
    s = rng.randint(0, n - L)
    return rna[s:s + L], s, s + L, "internal"


def feature_rna(f, contigs):
    seq = contigs[f["contig"]][f["start"]:f["end"]]
    return rc(seq) if f["strand"] == "-" else seq


def add_errors(read, rng, base=0.001, slope=0.004):
    n, out = len(read), list(read)
    for i in range(n):
        if rng.random() < base + slope * i / max(n - 1, 1):
            out[i] = rng.choice([b for b in "ACGT" if b != out[i]])
    return "".join(out)


def make_read(insert, rng, read_len, err):
    if len(insert) >= read_len:
        r = insert[:read_len]
    else:
        tail = ADAPTER + "".join(rng.choice("ACGT") for _ in range(read_len))
        r = (insert + tail)[:read_len]
    if err:
        r = add_errors(r, rng)
    return r


def quality(n):
    # Illumina-like: high at the start, slowly decreasing (Phred 37 -> 30)
    return "".join(chr(33 + int(37 - 7 * i / max(n - 1, 1))) for i in range(n))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--assembly", required=True)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--contig2bin", default="")
    ap.add_argument("--contig-counts", default="")
    ap.add_argument("--human-mirna", required=True)
    ap.add_argument("--human-trna", required=True)
    ap.add_argument("--human-rrna", required=True)
    ap.add_argument("--human-other", required=True)
    ap.add_argument("--n-reads", type=int, default=2000000)
    ap.add_argument("--read-len", type=int, default=51)
    ap.add_argument("--host-frac", type=float, default=0.02, help="share of human reads [0.02]")
    ap.add_argument("--dimer-frac", type=float, default=0.03, help="share of adapter dimers [0.03]")
    ap.add_argument("--holdout-frac", type=float, default=0.10,
                    help="share of contigs (by abundance-weighted length) held out of the assembly [0.10]")
    ap.add_argument("--hide-frac", type=float, default=0.30,
                    help="share of ncRNA / regulatory features removed from the GFF [0.30]")
    ap.add_argument("--min-contig", type=int, default=500, help="contigs shorter than this are not used [500]")
    ap.add_argument("--no-errors", action="store_true", help="no sequencing errors")
    ap.add_argument("--human-skew", type=float, default=1.5,
                    help="sigma of the log-normal abundance of the human references (miRNAs etc.); 0 = uniform [1.5]")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--sample", default="SIM01")
    ap.add_argument("--type", default="sRNA")
    ap.add_argument("--out", required=True, help="output prefix (directory/name)")
    a = ap.parse_args()
    rng = random.Random(a.seed)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    log = lambda m: sys.stderr.write("[simulate] %s\n" % m)

    contigs = {c: s for c, s in read_fasta(a.assembly).items() if len(s) >= a.min_contig}
    log("%d contigs >= %d bp" % (len(contigs), a.min_contig))
    feats = read_gff(a.gff, contigs)
    log("%d annotated features (%s)" % (len(feats), ", ".join(
        "%s %d" % kv for kv in sorted(collections.Counter(f["cls"] for f in feats).items()))))
    c2b = read_table(a.contig2bin)
    w = abundances(contigs, c2b, read_table(a.contig_counts, 4) if a.contig_counts else {}, rng)

    # hold out contigs: random contigs until --holdout-frac of the abundance-weighted sequence is reached
    names = sorted(contigs)
    rng.shuffle(names)
    total = sum(w[c] * len(contigs[c]) for c in names) or 1.0
    held, acc = set(), 0.0
    for c in names:
        if acc / total >= a.holdout_frac:
            break
        held.add(c)
        acc += w[c] * len(contigs[c])
    log("%d contigs held out of the assembly (%.1f%% of abundance-weighted sequence)" % (len(held), 100 * acc / total))

    # hide ncRNA / regulatory features of kept contigs from the annotation
    hidden = set()
    for f in feats:
        if f["cls"] in HIDEABLE and f["contig"] not in held and rng.random() < a.hide_frac:
            hidden.add(id(f))
    log("%d ncRNA/regulatory features hidden from the GFF" % len(hidden))

    # pickers per microbial class: features weighted by contig abundance (x length for CDS and intergenic)
    by_cls = collections.defaultdict(list)
    for f in feats:
        by_cls[f["cls"]].append(f)
    pick = {}
    for cls, fl in by_cls.items():
        L = (lambda f: f["end"] - f["start"]) if cls == "CDS" else (lambda f: 1.0)
        pick[cls] = Picker(fl, [w[f["contig"]] * L(f) for f in fl], rng)
    if "CDS" in pick:
        pick["CDS_antisense"] = pick["CDS"]
    clist = [c for c in names if w[c] > 0]
    pick["intergenic"] = Picker(clist, [w[c] * len(contigs[c]) for c in clist], rng)
    occupied = collections.defaultdict(list)
    for f in feats:
        occupied[f["contig"]].append((f["start"], f["end"]))
    mix = {k: v for k, v in MICROBIAL_MIX.items() if pick.get(k)}
    s = sum(mix.values())
    mix = {k: v / s for k, v in mix.items()}
    log("microbial class mix: " + ", ".join("%s %.2f" % kv for kv in mix.items()))
    mclasses, mweights = list(mix), [mix[k] for k in mix]

    human = {"miRNA": read_fasta(a.human_mirna), "tRNA": read_fasta(a.human_trna),
             "rRNA": read_fasta(a.human_rrna), "other": read_fasta(a.human_other)}
    # a few dominant human miRNAs / ncRNAs and a long tail, as in real libraries (log-normal weights)
    hpick = {k: Picker(list(v.items()), [rng.lognormvariate(0, a.human_skew) if a.human_skew > 0 else 1.0
                                         for _ in v], rng) for k, v in human.items() if v}
    hmix = {k: v for k, v in HUMAN_MIX.items() if k in hpick}
    hclasses, hweights = list(hmix), list(hmix.values())

    fq = gzip.open(a.out + ".fastq.gz", "wt", compresslevel=4)
    tr = gzip.open(a.out + ".truth.tsv.gz", "wt", compresslevel=4)
    tr.write("read_id\torigin\tclass\tfeature\tcontig\tbin\tstrand\tstart\tend\tinsert_len\thidden\tfragment\n")
    qual = quality(a.read_len)
    n_cls = collections.Counter()
    for i in range(a.n_reads):
        rid = "sim%09d" % (i + 1)
        u = rng.random()
        feature, contig, binname, strand, s0, e0, hid, kind = "-", "-", "-", "+", 0, 0, "no", "-"
        if u < a.dimer_frac:
            origin, cls = "dimer", "dimer"
            insert = "".join(rng.choice("ACGT") for _ in range(rng.randint(0, 5)))
        elif u < a.dimer_frac + a.host_frac:
            origin = "human"
            cls = rng.choices(hclasses, hweights)[0]
            name, rna = hpick[cls].pick()
            feature = name
            if cls == "miRNA":                          # mature miRNA with a 3' isomiR (-2..+2 nt)
                d = rng.choices([-2, -1, 0, 1, 2], [0.05, 0.2, 0.5, 0.2, 0.05])[0]
                insert = rna[:max(16, len(rna) + d)] if d <= 0 else rna + "".join(rng.choice("ACGT") for _ in range(d))
                s0, e0, kind = 0, len(insert), "isomiR%+d" % d if d else "mature"
            elif cls == "tRNA":
                insert, s0, e0, kind = trna_fragment(rna, rng)
            else:
                insert, s0, e0, kind = generic_fragment(rna, rng, five_prime=0.2)
        else:
            cls = rng.choices(mclasses, mweights)[0]
            if cls == "intergenic":
                contig = pick["intergenic"].pick()
                seq = contigs[contig]
                for _ in range(20):                     # a stretch that overlaps no annotated feature
                    L = rng.randint(18, 40)
                    s0 = rng.randint(0, max(0, len(seq) - L)); e0 = s0 + L
                    if not any(s0 < fe and fs < e0 for fs, fe in occupied[contig]):
                        break
                strand = rng.choice("+-")
                insert = seq[s0:e0] if strand == "+" else rc(seq[s0:e0])
            else:
                f = pick[cls].pick()
                contig, feature = f["contig"], f["id"]
                rna = feature_rna(f, contigs)
                if cls == "tRNA":
                    insert, a0, b0, kind = trna_fragment(rna, rng)
                elif cls == "rRNA":
                    insert, a0, b0, kind = generic_fragment(rna, rng, five_prime=0.15)
                else:
                    insert, a0, b0, kind = generic_fragment(rna, rng, mean=30, sd=9)
                b0 = min(b0, f["end"] - f["start"])          # CCA added to tRNAs is not on the contig
                # transcript coordinates -> contig coordinates
                if f["strand"] == "+":
                    s0, e0 = f["start"] + a0, f["start"] + b0
                else:
                    s0, e0 = f["end"] - b0, f["end"] - a0
                strand = f["strand"]
                if cls == "CDS_antisense":
                    insert = rc(insert)
                    strand = "-" if strand == "+" else "+"
                hid = "yes" if id(f) in hidden else "no"
            origin = "microbial_unassembled" if contig in held else "microbial_assembled"
            binname = c2b.get(contig, ["unbinned"])[0]
        n_cls[(origin, cls)] += 1
        read = make_read(insert, rng, a.read_len, not a.no_errors)
        fq.write("@%s\n%s\n+\n%s\n" % (rid, read, qual))
        tr.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\n" % (
            rid, origin, cls, feature, contig, binname, strand, s0, e0, len(insert), hid, kind))
        if (i + 1) % 1000000 == 0:
            log("%d reads" % (i + 1))
    fq.close()
    tr.close()

    # the reduced assembly / annotation given to the pipeline
    with gzip.open(a.out + ".assembly.fa.gz", "wt") as fh:
        for c in sorted(contigs):
            if c not in held:
                fh.write(">%s\n" % c)
                for k in range(0, len(contigs[c]), 80):
                    fh.write(contigs[c][k:k + 80] + "\n")
    keep_ids = {id(f) for f in feats if f["contig"] not in held and id(f) not in hidden}
    with gzip.open(a.out + ".gff.gz", "wt") as fh, open(a.out + ".hidden_loci.tsv", "w") as hl:
        fh.write("##gff-version 3\n")
        hl.write("feature\tcontig\tstart\tend\tstrand\tclass\ttype\n")
        for f in feats:
            if id(f) in keep_ids:
                fh.write(f["line"])
            elif id(f) in hidden:
                hl.write("%s\t%s\t%d\t%d\t%s\t%s\t%s\n" % (f["id"], f["contig"], f["start"], f["end"], f["strand"],
                                                         f["cls"], f["type"]))
    c2b_out = ""
    if c2b:
        c2b_out = a.out + ".contig2bin.tsv"
        with open(c2b_out, "w") as fh:
            for c in sorted(contigs):
                if c not in held and c in c2b:
                    fh.write("%s\t%s\n" % (c, c2b[c][0]))
    with open(a.out + ".params.tsv", "w") as fh:
        for k, v in sorted(vars(a).items()):
            fh.write("%s\t%s\n" % (k, v))
        fh.write("held_out_contigs\t%d\nhidden_features\t%d\n" % (len(held), len(hidden)))
        for (o, c), n in sorted(n_cls.items()):
            fh.write("reads_%s_%s\t%d\n" % (o, c, n))
    P = os.path.abspath
    with open(a.out + ".sample_row.tsv", "w") as fh:
        fh.write("\t".join([a.sample, a.type, P(a.out + ".fastq.gz"), P(a.out + ".assembly.fa.gz"), P(a.out + ".gff.gz"),
                            P(c2b_out) if c2b_out else "-", "-", "-", "-", "-", "-"]) + "\n")
    log("done: %s.{fastq.gz,truth.tsv.gz,assembly.fa.gz,gff.gz,hidden_loci.tsv,sample_row.tsv}" % a.out)
    for (o, c), n in sorted(n_cls.items()):
        log("  %-22s %-14s %9d (%.2f%%)" % (o, c, n, 100.0 * n / a.n_reads))


if __name__ == "__main__":
    main()
