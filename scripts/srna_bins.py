#!/usr/bin/env python3
"""
srna_bins.py — sRNA read counts per contig, annotated feature, bin (MAG) and taxon

Reads all alignments of a BAM (samtools in PATH; multi-mappers weighted 1/n, n from the NH tag
written by the assembly module, otherwise counted in a first pass),
assigns each alignment to the annotated feature (gff2bed BED) it overlaps
best: same strand -> "<type>", opposite strand only -> "antisense_<type>",
nothing -> "intergenic". Features of type gene/region/source/exon are ignored
when a more specific feature overlaps.

Outputs (<prefix>.*):
  contig_counts.tsv     contig, length, bin, weighted reads, reads per kb
  feature_counts.tsv    feature_id, contig, start, end, strand, type, bin, sense, antisense, gene, product
  functions.tsv         (--functions) reads on CDS per functional term: db, term, features, sense, antisense
                        (GFF attributes cog, kegg_ko, kegg_pathway, kegg_module, kegg_brite, enzyme_ec, go,
                        pfam, tigrfam, eggnog; comma-separated values each get the full read weight)
  length_dist.tsv       (--lengths) read length: weighted mapped reads, unmapped reads
  category_counts.tsv   assignment category totals (CDS, tRNA, antisense_CDS, intergenic, ...);
                        *_non_rRNA columns (feature / category / functions) without the --exclude-reads
  (--no-bins skips the two bin tables)
  bin_counts.tsv        bin, taxonomy, contigs, length, reads, fraction, per-category reads
  taxon_counts.tsv      rank, taxon, reads, fraction (GTDB ranks d p c o f g s)
"""
import argparse
import gzip
import bisect
import collections
import os
import re
import subprocess
import sys


def open_in(path):
    """open a text input, transparently gunzipping *.gz"""
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
GENERIC = {"gene", "region", "source", "exon", "biological_region", "databank_entry"}
RANKS = ["d", "p", "c", "o", "f", "g", "s"]


def ref_len(cigar):
    return sum(int(n) for n, op in CIGAR_RE.findall(cigar) if op in "MDN=X")


def strip_ext(name):
    return re.sub(r"\.(fa|fna|fasta|fsa)(\.gz)?$", "", name)


def read_contig2bin(path):
    d = {}
    if not path:
        return d
    with open_in(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                d[f[0]] = strip_ext(f[1])
    return d


def read_gtdbtk(path):
    """GTDB-Tk summary (bac120 / ar53, or a concatenation): user_genome, classification"""
    tax = {}
    if not path:
        return tax
    paths = path.split(",")
    for pth in paths:
        with open_in(pth) as fh:
            header = None
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if header is None or f[0] == "user_genome":
                    header = f
                    continue
                rec = dict(zip(header, f))
                g = strip_ext(rec.get("user_genome", f[0]))
                tax[g] = rec.get("classification", "")
    return tax


def read_binstats(path):
    """bin quality from the same table(s) as the taxonomy (CheckM/CheckM2 columns
    Completeness, Contamination, case-insensitive): bin -> (completeness, contamination, quality)
    quality (MIMAG-like): HQ >= 90 % complete and < 5 % contaminated, MQ >= 50 % and < 10 %, else LQ"""
    q = {}
    if not path:
        return q
    for pth in [p for p in path.split(",") if p]:
        with open_in(pth) as fh:
            header = None
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if header is None or f[0] == "user_genome":
                    header = [h.lower() for h in f]
                    continue
                rec = dict(zip(header, f))
                try:
                    comp = float(rec.get("completeness", ""))
                    cont = float(rec.get("contamination", ""))
                except ValueError:
                    continue
                cls = "HQ" if comp >= 90 and cont < 5 else ("MQ" if comp >= 50 and cont < 10 else "LQ")
                q[strip_ext(rec.get("user_genome", rec.get("bin", f[0])))] = ("%.2f" % comp, "%.2f" % cont, cls)
    return q


def tax_ranks(classification):
    d = {r: "" for r in RANKS}
    for part in classification.split(";"):
        part = part.strip()
        if len(part) > 3 and part[1:3] == "__" and part[0] in d:
            d[part[0]] = part
    return d


FUNC_TAGS = ["cog", "kegg_ko", "kegg_pathway", "kegg_module", "kegg_brite", "enzyme_ec", "go", "pfam",
             "tigrfam", "eggnog"]


def parse_attr(s):
    out = {}
    for x in s.split(";"):
        if "=" in x:
            k, v = x.split("=", 1)
            out[k.strip()] = v.strip()
    return out


RESCUED = []          # indices of rRNA genes rescued by cmscan (srna_rrna_rescue.py, source column 'cmscan')


def flag_features(feats):
    """CDS overlapping a rescued rRNA gene on the opposite strand by >= 50 % of the shorter of the two: ORFs
    predicted on the other strand of an rRNA gene missing from the annotation -> {idx: flag}"""
    by = collections.defaultdict(list)
    for i in RESCUED:
        by[feats[i][1]].append(feats[i])
    out = {}
    if not by:
        return out
    for idx, (fid, c, s, e, st, ty) in enumerate(feats):
        if ty != "CDS" or c not in by:
            continue
        for _, _, rs, re_, rst, _ in by[c]:
            ov = min(e, re_) - max(s, rs)
            if rst != st and ov > 0 and ov >= 0.5 * min(e - s, re_ - rs):
                out[idx] = "opposite_rescued_rRNA"
                break
    return out


def load_features(gbed, want_attrs=False):
    """gff2bed BED -> feats [(id, contig, start, end, strand, type)], index, attrs {idx: (gene, product, funcs)}"""
    attrs = {}
    feats = []
    by_contig = collections.defaultdict(list)
    with open_in(gbed) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            fid = f[3]
            if fid in (".", ""):
                m = re.search(r"(?:^|;)ID=([^;]+)", f[9]) if len(f) > 9 else None
                fid = m.group(1) if m else "%s:%s-%s(%s)" % (f[0], f[1], f[2], f[5])
            idx = len(feats)
            feats.append((fid, f[0], int(f[1]), int(f[2]), f[5], f[7]))
            if f[7] == "rRNA" and f[6] == "cmscan":
                RESCUED.append(idx)
            if want_attrs and len(f) > 9 and f[7] not in GENERIC:
                at = parse_attr(f[9])
                funcs = tuple((t, v) for t in FUNC_TAGS for v in at.get(t, "").split(",") if v) \
                    if f[7] == "CDS" else ()
                attrs[idx] = (at.get("gene", ""), at.get("product", at.get("Name", "")).replace("\t", " "), funcs)
            by_contig[f[0]].append((int(f[1]), int(f[2]), idx))
    index = {}
    for c, lst in by_contig.items():
        lst.sort()
        index[c] = ([s for s, _, _ in lst], lst, max(e - s for s, e, _ in lst))
    return feats, index, attrs


def candidates(index, contig, rs, re_):
    if contig not in index:
        return []
    starts, lst, maxlen = index[contig]
    i = bisect.bisect_right(starts, re_ - 1)
    out = []
    j = i - 1
    while j >= 0 and starts[j] > rs - maxlen - 1:
        s, e, idx = lst[j]
        if s < re_ and e > rs:
            out.append(idx)
        j -= 1
    return out


def samtools(args):
    p = subprocess.Popen(["samtools"] + args, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    for line in p.stdout:
        yield line
    if p.wait() != 0:
        sys.exit("samtools %s failed" % " ".join(args))


def read_names(path):
    """read names (first word of the header, without @ and /1) of a FASTQ(.gz)"""
    op = gzip.open if path.endswith(".gz") else open
    names = set()
    with op(path, "rt") as fh:
        for i, ln in enumerate(fh):
            if i % 4 == 0:
                n = ln[1:].split(None, 1)[0]
                names.add(n[:-2] if n.endswith(("/1", "/2")) else n)
    return names


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--gbed", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--contig2bin", default="")
    ap.add_argument("--gtdbtk", default="", help="GTDB-Tk summary tsv (comma-separated for bac120,ar53)")
    ap.add_argument("--min-overlap", type=float, default=0.5)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--functions", action="store_true", help="write <prefix>.functions.tsv (CDS functional terms)")
    ap.add_argument("--lengths", action="store_true", help="write <prefix>.length_dist.tsv")
    ap.add_argument("--no-bins", action="store_true", help="skip bin_counts.tsv / taxon_counts.tsv")
    ap.add_argument("--exclude-reads", default="",
                    help="FASTQ(.gz) of reads to leave out of the *_non_rRNA columns (e.g. the RiboDetector rRNA reads): "
                         "reads sense to unannotated rRNA would otherwise count as antisense to spurious ORFs")
    ap.add_argument("--exclude-seqs", default="",
                    help="gz/text file, one read sequence per line (both orientations; srna_rrna_seqs.py): alignments "
                         "with such a SEQ are also left out of the *_non_rRNA columns (Rfam rRNA hits missed by RiboDetector)")
    a = ap.parse_args()

    excl = read_names(a.exclude_reads) if a.exclude_reads else set()
    if a.exclude_reads:
        sys.stderr.write("[bins] %d reads excluded from the *_non_rRNA counts (%s)\n" % (len(excl), a.exclude_reads))
    xseq = set()
    if a.exclude_seqs:
        op = gzip.open if a.exclude_seqs.endswith(".gz") else open
        with op(a.exclude_seqs, "rt") as fh:
            xseq = {ln.strip() for ln in fh if ln.strip()}
        sys.stderr.write("[bins] + %d sequences (both strands) excluded from the *_non_rRNA counts (%s)\n" % (
            len(xseq), a.exclude_seqs))
    excl_any = bool(excl or xseq)
    n_x = n_xs = 0.0
    c2b = read_contig2bin(a.contig2bin)
    tax = read_gtdbtk(a.gtdbtk)
    qual = read_binstats(a.gtdbtk)
    feats, index, attrs = load_features(a.gbed, want_attrs=True)

    clen = {}
    for line in samtools(["view", "-H", a.bam]):
        if line.startswith("@SQ"):
            d = dict(x.split(":", 1) for x in line.rstrip("\n").split("\t")[1:] if ":" in x)
            clen[d.get("SN")] = int(d.get("LN", 0))

    # hits per read: NH tag (assembly module) or a first pass over the BAM
    use_nh = False
    for line in samtools(["view", "-F", "4", a.bam]):
        use_nh = "\tNH:i:" in line
        break
    nali = collections.Counter()
    if not use_nh:
        for line in samtools(["view", "-@", str(a.threads), "-F", "4", a.bam]):
            nali[line.split("\t", 1)[0]] += 1
    nh_re = re.compile(r"\tNH:i:(\d+)")
    lenw = collections.defaultdict(float)
    nreads = 0.0

    ccount = collections.defaultdict(float)
    sense = collections.defaultdict(float)
    anti = collections.defaultdict(float)
    sense_nr = collections.defaultdict(float)      # without the excluded (rRNA) reads
    anti_nr = collections.defaultdict(float)
    cat = collections.defaultdict(float)
    cat_nr = collections.defaultdict(float)
    bincat = collections.defaultdict(lambda: collections.defaultdict(float))

    for line in samtools(["view", "-@", str(a.threads), "-F", "4", a.bam]):
        f = line.split("\t", 10)
        name, flag, contig, pos, cigar = f[0], int(f[1]), f[2], int(f[3]) - 1, f[5]
        if use_nh:
            m = nh_re.search(line)
            w = 1.0 / int(m.group(1)) if m else 1.0
        else:
            w = 1.0 / nali[name]
        nreads += w
        if a.lengths:
            lenw[len(f[9])] += w
        rlen = ref_len(cigar)
        rs, re_ = pos, pos + rlen
        strand = "-" if flag & 16 else "+"
        ccount[contig] += w
        best_s = best_a = None
        ov_s = ov_a = 0
        spec_s = spec_a = False
        for idx in candidates(index, contig, rs, re_):
            fid, _, s, e, fstr, ftype = feats[idx]
            ov = min(re_, e) - max(rs, s)
            if ov < a.min_overlap * rlen:
                continue
            specific = ftype not in GENERIC
            if fstr == strand or fstr == ".":
                if (specific, ov) > (spec_s, ov_s):
                    best_s, ov_s, spec_s = idx, ov, specific
            else:
                if (specific, ov) > (spec_a, ov_a):
                    best_a, ov_a, spec_a = idx, ov, specific
        keep = name not in excl
        if not keep:
            n_x += w
        elif xseq and f[9] in xseq:
            keep = False
            n_xs += w
        if best_s is not None:
            sense[best_s] += w
            if keep:
                sense_nr[best_s] += w
            k = feats[best_s][5]
        elif best_a is not None:
            anti[best_a] += w
            if keep:
                anti_nr[best_a] += w
            k = "antisense_" + feats[best_a][5]
        else:
            k = "intergenic"
        cat[k] += w
        if keep:
            cat_nr[k] += w
        bincat[c2b.get(contig, "unbinned")][k] += w

    if excl_any:
        sys.stderr.write("[bins] left out of *_non_rRNA: %.0f weighted reads by name (RiboDetector), %.0f more by "
                         "sequence (Rfam rRNA), of %.0f mapped\n" % (n_x, n_xs, nreads))
    flags = flag_features(feats)
    if RESCUED:
        sys.stderr.write("[bins] %d rescued rRNA genes; %d CDS on their opposite strand flagged (left out of "
                         "functions.tsv)\n" % (len(RESCUED), len(flags)))
    p = a.out_prefix
    os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)
    total = sum(ccount.values())
    with open(p + ".contig_counts.tsv", "w") as out:
        out.write("contig\tlength\tbin\treads\treads_per_kb\n")
        for c in sorted(ccount, key=lambda x: -ccount[x]):
            L = clen.get(c, 0)
            out.write("%s\t%d\t%s\t%.3f\t%.3f\n" % (c, L, c2b.get(c, "unbinned"), ccount[c],
                                                   ccount[c] / (L / 1000.0) if L else 0))
    with open(p + ".feature_counts.tsv", "w") as out:
        # *_non_rRNA columns only when reads were excluded, so their presence means "rRNA reads removed"
        out.write("feature_id\tcontig\tstart\tend\tstrand\ttype\tbin\tsense_reads\tantisense_reads\tgene\tproduct" +
                  ("\tsense_reads_non_rRNA\tantisense_reads_non_rRNA" if excl_any else "") + "\tflag\n")
        for idx in sorted(set(sense) | set(anti)):
            fid, c, s, e, fstr, ftype = feats[idx]
            g, prod, _f = attrs.get(idx, ("", "", ()))
            out.write("%s\t%s\t%d\t%d\t%s\t%s\t%s\t%.3f\t%.3f\t%s\t%s" % (
                fid, c, s, e, fstr, ftype, c2b.get(c, "unbinned"), sense[idx], anti[idx], g, prod) +
                ("\t%.3f\t%.3f" % (sense_nr[idx], anti_nr[idx]) if excl_any else "") + "\t%s\n" % flags.get(idx, ""))
    if a.functions:
        fs = collections.defaultdict(lambda: [0, 0.0, 0.0, 0.0, 0.0])
        snr, anr = (sense_nr, anti_nr) if excl_any else (sense, anti)
        for idx in set(sense) | set(anti):
            if idx in flags:                       # spurious ORFs on the other strand of rescued rRNA genes
                continue
            for t, v in attrs.get(idx, ("", "", ()))[2]:
                x = fs[(t, v)]
                x[0] += 1
                x[1] += sense[idx]
                x[2] += anti[idx]
                x[3] += snr[idx]
                x[4] += anr[idx]
        with open(p + ".functions.tsv", "w") as out:
            out.write("db\tterm\tfeatures\tsense_reads\tantisense_reads" +
                      ("\tsense_reads_non_rRNA\tantisense_reads_non_rRNA" if excl_any else "") + "\n")
            for (t, v), x in sorted(fs.items(), key=lambda kv: (kv[0][0], -kv[1][1])):
                out.write("%s\t%s\t%d\t%.3f\t%.3f" % (t, v, x[0], x[1], x[2]) +
                          ("\t%.3f\t%.3f" % (x[3], x[4]) if excl_any else "") + "\n")
    if a.lengths:
        unl = collections.Counter()
        for line in samtools(["view", "-@", str(a.threads), "-f", "4", a.bam]):
            unl[len(line.split("\t", 10)[9])] += 1
        with open(p + ".length_dist.tsv", "w") as out:
            out.write("length\tmapped_weighted\tunmapped\n")
            for L in sorted(set(lenw) | set(unl)):
                out.write("%d\t%.3f\t%d\n" % (L, lenw.get(L, 0.0), unl.get(L, 0)))
    with open(p + ".category_counts.tsv", "w") as out:
        out.write("category\treads\tfraction" + ("\treads_non_rRNA" if excl_any else "") + "\n")
        for k in sorted(cat, key=lambda x: -cat[x]):
            out.write("%s\t%.3f\t%.4f" % (k, cat[k], cat[k] / total if total else 0) +
                      ("\t%.3f" % cat_nr[k] if excl_any else "") + "\n")

    if a.no_bins:
        sys.stderr.write("[bins] %.1f weighted reads, %d contigs with reads (bins skipped)\n" % (total, len(ccount)))
        return
    # bins
    bin_len = collections.defaultdict(int)
    bin_ncont = collections.defaultdict(int)
    for c, L in clen.items():
        b = c2b.get(c, "unbinned")
        bin_len[b] += L
        bin_ncont[b] += 1
    bin_reads = collections.defaultdict(float)
    for c, v in ccount.items():
        bin_reads[c2b.get(c, "unbinned")] += v
    cats = sorted(cat, key=lambda x: -cat[x])
    with open(p + ".bin_counts.tsv", "w") as out:
        out.write("\t".join(["bin", "classification", "completeness", "contamination", "quality",
                             "contigs", "length", "reads", "fraction", "reads_per_Mb"] + cats) + "\n")
        for b in sorted(bin_reads, key=lambda x: -bin_reads[x]):
            L = bin_len.get(b, 0)
            out.write("\t".join([b, tax.get(b, "unclassified" if b != "unbinned" else "")] +
                                list(qual.get(b, ("", "", "") if b == "unbinned" else ("NA", "NA", "NA"))) +
                                [str(bin_ncont.get(b, 0)), str(L), "%.3f" % bin_reads[b],
                                 "%.4f" % (bin_reads[b] / total if total else 0),
                                 "%.3f" % (bin_reads[b] / (L / 1e6) if L else 0)] +
                                ["%.3f" % bincat[b].get(k, 0.0) for k in cats]) + "\n")
    with open(p + ".taxon_counts.tsv", "w") as out:
        out.write("rank\ttaxon\treads\tfraction\n")
        for r in RANKS:
            agg = collections.defaultdict(float)
            for b, v in bin_reads.items():
                t = tax_ranks(tax.get(b, "")).get(r, "") if b in tax else ""
                if not t:
                    t = "unbinned" if b == "unbinned" else ("%s__unclassified" % r)
                agg[t] += v
            for t in sorted(agg, key=lambda x: -agg[x]):
                out.write("%s\t%s\t%.3f\t%.4f\n" % (r, t, agg[t], agg[t] / total if total else 0))
    sys.stderr.write("[bins] %d reads, %.1f weighted, %d contigs with reads, %d bins\n"
                     % (len(nali) if nali else round(nreads), total, len(ccount), len([b for b in bin_reads if b != "unbinned"])))


if __name__ == "__main__":
    main()
