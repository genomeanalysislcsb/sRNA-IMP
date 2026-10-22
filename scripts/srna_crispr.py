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
    sys.stderr.write("[crispr] %d arrays, %.1f reads on arrays, %d arrays with >= 10 reads\n"
                     % (len(arrays), sum(sum(r[3] for r in v) for v in reads.values()),
                        sum(1 for v in reads.values() if sum(r[3] for r in v) >= 10)))


if __name__ == "__main__":
    main()
