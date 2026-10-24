#!/usr/bin/env python3
"""srna_highlights.py — cross-sample highlights for the sRNA-IMP cohort report (used by srna_cohort.py).

  crispr_systems   CRISPR arrays of all libraries grouped by their repeat (same repeat = same CRISPR system, in
                   either orientation): samples / libraries with an expressed array, Cas type, taxon, expression,
                   transcribed strand, samples with a tracrRNA candidate
  novel_clusters   novel sRNA loci without Rfam hit (module novel, RNAfold set) of all libraries, clustered across
                   samples by shared sequence (exact 21-mers, content-sampled): recurrent novel sRNAs with samples,
                   libraries, expression, length, locus class, MFE, taxon and a representative sequence
  novel_rfam       Rfam families found among the novel loci (rRNA / tRNA families left out): samples, libraries
  chimeras         candidate RNA–RNA pairs (module chimera) with q < 0.05, pooled over samples by partner names,
                   and per-library complementarity / same-organism tests of the sRNA-like–mRNA chimeras
"""
import collections
import glob
import gzip
import math
import os
import re
import zlib

COMP = str.maketrans("ACGTUNacgtun", "TGCAANtgcaan")
RFAM_SKIP = re.compile(r"rRNA|tRNA|ribosomal", re.I)


def rc(s):
    return s.translate(COMP)[::-1]


def open_in(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def read_tsv(path):
    if not path or not os.path.exists(path):
        return []
    with open_in(path) as fh:
        lines = [ln.rstrip("\n").split("\t") for ln in fh if ln.strip() and not ln.startswith("#")]
    if not lines:
        return []
    h = lines[0]
    return [dict(zip(h, x)) for x in lines[1:]]


def num(x, default=0.0):
    try:
        v = float(x)
        return default if math.isnan(v) else v
    except (TypeError, ValueError):
        return default


def median(v):
    v = sorted(v)
    if not v:
        return float("nan")
    m = len(v) // 2
    return v[m] if len(v) % 2 else 0.5 * (v[m - 1] + v[m])


def tax_short(cl):
    for r in reversed([x for x in (cl or "").split(";") if x]):
        if len(r) > 3:
            return r[3:]
    return ""


def top_label(counter, n=1):
    items = [(k, v) for k, v in counter.most_common() if k]
    if not items:
        return ""
    lab = " / ".join(k for k, _ in items[:n])
    return lab + (" (+%d)" % (len(items) - n) if len(items) > n else "")


def bin_taxonomy(outdir, s, t):
    """contig -> taxon (GTDB, deepest rank) from the bins module of the library"""
    d = os.path.join(outdir, s, "Analysis", t, "bins")
    b2t = {r.get("bin", ""): tax_short(r.get("classification", "")) for r in read_tsv(os.path.join(d, "%s.bin_counts.tsv" % s))}
    c2t = {}
    p = os.path.join(d, "%s.contig_counts.tsv" % s)
    if os.path.exists(p):
        with open_in(p) as fh:
            next(fh, None)
            for ln in fh:
                f = ln.split("\t", 3)
                if len(f) >= 3 and f[2] != "unbinned":
                    c2t[f[0]] = b2t.get(f[2], "")
    return c2t


# ------------------------------------------------------------------ CRISPR
def crispr_systems(outdir, pairs, group_of, min_reads=20):
    agg = {}
    for s, t in pairs:
        d = os.path.join(outdir, s, "Analysis", t, "crispr")
        arrays = read_tsv(os.path.join(d, "%s.crispr_arrays.tsv" % s))
        if not arrays:
            continue
        tracr = {r["array_id"] for r in read_tsv(os.path.join(d, "%s.crispr_tracr.tsv" % s))
                 if r.get("tracr_candidate") == "yes"}
        for r in arrays:
            reads = num(r.get("reads"))
            if reads < min_reads:
                continue
            rep = (r.get("repeat_consensus") or "").upper().replace("U", "T")
            key = min(rep, rc(rep)) if rep else "array:%s:%s" % (s, r.get("array_id"))
            g = agg.setdefault(key, dict(repeat=rep, samples=set(), libs=collections.Counter(), rpk=[], reads=0.0,
                                         cas=collections.Counter(), tax=collections.Counter(), strand=collections.Counter(),
                                         tracr=set(), groups={}, best=None, spacers=[]))
            g["samples"].add(s)
            g["groups"][s] = group_of(s)
            g["libs"][t] += 1
            g["rpk"].append(num(r.get("reads_per_kb")))
            g["reads"] += reads
            g["cas"][r.get("cas_type") or "–"] += 1
            g["tax"][tax_short(r.get("classification", ""))] += 1
            g["strand"]["array strand" if r.get("transcribed_strand") in ("+", "-") else "?"] += 1
            g["spacers"].append(num(r.get("spacers")))
            if r.get("array_id") in tracr:
                g["tracr"].add(s)
            if g["best"] is None or num(r.get("reads_per_kb")) > g["best"][3]:
                g["best"] = (s, t, r.get("array_id"), num(r.get("reads_per_kb")))
    rows = []
    for key, g in agg.items():
        grp = collections.Counter(g["groups"].values())
        rows.append(dict(repeat=g["repeat"], repeat_len=len(g["repeat"]), cas_type=top_label(g["cas"]),
                         taxon=top_label(g["tax"]), samples=len(g["samples"]),
                         groups=" / ".join("%s %d" % (k, v) for k, v in sorted(grp.items())),
                         libraries=", ".join("%s %d" % (k.replace("sRNA-", ""), v) for k, v in sorted(g["libs"].items())),
                         median_reads_per_kb=median(g["rpk"]), max_reads_per_kb=max(g["rpk"]), reads=g["reads"],
                         median_spacers=median(g["spacers"]), tracr_samples=len(g["tracr"]),
                         example="%s %s %s" % g["best"][:3], _grp=dict(grp),
                         _samples=set(g["samples"])))
    rows.sort(key=lambda r: (-r["samples"], -r["median_reads_per_kb"]))
    return rows


def _lchoose(n, k):
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def fisher_greater(a, b, c, d):
    """one-sided Fisher exact test (2x2 [[a, b], [c, d]]): P(X >= a)"""
    n1, n2, k = a + b, c + d, a + c
    lo, hi = max(0, k - n2), min(n1, k)
    den = _lchoose(n1 + n2, k)
    return min(1.0, sum(math.exp(_lchoose(n1, x) + _lchoose(n2, k - x) - den) for x in range(max(a, lo), hi + 1)))


def bh(ps):
    m = len(ps)
    order = sorted(range(m), key=lambda i: ps[i])
    q, prev = [1.0] * m, 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        prev = min(prev, ps[i] * m / rank)
        q[i] = prev
    return q


def group_presence(items, gb, min_samples=3, max_p=0.05):
    """features present mainly in one group. items: dicts with cls (feature class), feature, info, samples (with the
    feature) and universe (samples that could have it). gb: sample -> (group, block). For each feature and group g
    with >= min_samples samples: samples of g vs samples of the other groups in the same blocks as g (e.g. the same
    site, so that a group recruited at one site is compared with the controls of that site only); one-sided Fisher
    test, BH-FDR within each class; returns tests with p < max_p, sorted by q, p."""
    tests = []
    for it in items:
        uni = [s for s in it["universe"] if gb.get(s, ("NA", ""))[0] != "NA"]
        groups = sorted({gb[s][0] for s in uni})
        if len(groups) < 2:
            continue
        has = it["samples"]
        for g in groups:
            G = [s for s in uni if gb[s][0] == g]
            blocks = {gb[s][1] for s in G}
            O = [s for s in uni if gb[s][0] != g and gb[s][1] in blocks]
            a = sum(1 for s in G if s in has)
            if a < min_samples or not O:
                continue
            c = sum(1 for s in O if s in has)
            p = fisher_greater(a, len(G) - a, c, len(O) - c)
            ref = "/".join(sorted({gb[s][0] for s in O}))
            if len(blocks) == 1 and list(blocks)[0]:
                ref += " (%s)" % list(blocks)[0]
            tests.append((it, g, a, c, len(G), len(O), ref, p))
    out = []
    for cls in sorted({t[0]["cls"] for t in tests}):
        sub = [t for t in tests if t[0]["cls"] == cls]
        for (it, g, a, c, ng, no, ref, p), q in zip(sub, bh([t[7] for t in sub])):
            if p < max_p:
                out.append(dict(cls=cls, feature=it["feature"], info=it.get("info", ""), group=g,
                                exclusive="yes" if c == 0 else "no", samples_in_group=a, group_size=ng,
                                share_in_group=a / ng, compared_with=ref, samples_other=c, other_size=no,
                                share_other=c / no, p=p, q=q, tests_in_class=len(sub)))
    out.sort(key=lambda x: (x["q"], x["p"], -x["samples_in_group"]))
    return out


def crispr_group_specific(rows, gb, universe, min_samples=3, max_p=0.05):
    """CRISPR systems expressed (>= 20 reads) mainly in one group (see group_presence)"""
    items = [dict(cls="CRISPR system", feature=r["taxon"], info=r["cas_type"], samples=r["_samples"],
                  universe=universe, _row=r) for r in rows]
    out = group_presence(items, gb, min_samples, max_p)
    byf = {}
    for it in items:
        byf.setdefault(it["feature"], it)
    res = []
    for o in out:
        r = next(it["_row"] for it in items if it["feature"] == o["feature"] and it["info"] == o["info"])
        o.update(taxon=r["taxon"], cas_type=r["cas_type"], median_reads_per_kb=r["median_reads_per_kb"],
                 tracr_samples=r["tracr_samples"], groups=r["groups"], repeat=r["repeat"])
        res.append(o)
    return res


# ------------------------------------------------------------------ novel sRNAs
HDR = re.compile(r"^>(.+):(\d+)-(\d+)\(([+-])\)")
MFE = re.compile(r"\(\s*(-?\d+(?:\.\d+)?)\)\s*$")


def novel_candidates(outdir, s, t, top):
    """novel loci without Rfam hit of one library: [(contig, start, end, strand, seq, mfe, reads, locus class)]"""
    d = os.path.join(outdir, s, "Analysis", t, "novel")
    expr = {}
    for cls, f in (("intergenic", "candidates.intergenic.bed"), ("antisense", "candidates.antisense_to_CDS.bed")):
        p = os.path.join(d, "%s.%s" % (s, f))
        if os.path.exists(p):
            with open(p) as fh:
                for ln in fh:
                    x = ln.rstrip("\n").split("\t")
                    if len(x) >= 7:
                        expr[(x[0], int(x[1]), int(x[2]), x[5])] = (num(x[6]), cls)
    out = []
    p = os.path.join(d, "%s.candidates.no_rfam.RNAfold.txt" % s)
    if not os.path.exists(p):
        return out
    with open(p) as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
    i = 0
    while i < len(lines):
        m = HDR.match(lines[i])
        if not m:
            i += 1
            continue
        seq = lines[i + 1].strip().upper().replace("U", "T") if i + 1 < len(lines) else ""
        mm = MFE.search(lines[i + 2]) if i + 2 < len(lines) else None
        key = (m.group(1), int(m.group(2)), int(m.group(3)), m.group(4))
        reads, cls = expr.get(key, (0.0, "?"))
        out.append(key + (seq, float(mm.group(1)) if mm else float("nan"), reads, cls))
        i += 3
    out.sort(key=lambda r: -r[6])
    return out[:top] if top else out


def sampled_kmers(seq, k=21, mod=6):
    """content-defined subset of the k-mers (same subset in every sample): crc32 % mod == 0"""
    return {zlib.crc32(seq[i:i + k].encode()) for i in range(len(seq) - k + 1)
            if zlib.crc32(seq[i:i + k].encode()) % mod == 0}


def novel_clusters(outdir, pairs, group_of, top_per_lib=300, max_kmer_libs=None):
    cands = []                                     # (s, t, contig, start, end, strand, seq, mfe, reads, class)
    for s, t in pairs:
        for c in novel_candidates(outdir, s, t, top_per_lib):
            if len(c[4]) >= 21:
                cands.append((s, t) + c)
    if not cands:
        return [], []
    parent = list(range(len(cands)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    owner = {}
    km_libs = collections.defaultdict(collections.Counter)
    kms = []
    for i, c in enumerate(cands):
        ks = sampled_kmers(c[6])
        kms.append(ks)
        for h in ks:
            km_libs[h][(c[0], c[1])] += 1
    # k-mers that occur in several loci of the same library are repeats (IS elements, multi-copy genes, motifs) and
    # would chain unrelated loci: they do not link
    rep = {h for h, cnt in km_libs.items() if max(cnt.values()) > 1}
    cap = max_kmer_libs or max(20, 2 * len(pairs))
    for i, ks in enumerate(kms):
        for h in ks:
            if h in rep or len(km_libs[h]) > cap:
                continue
            j = owner.setdefault(h, i)
            if j != i:
                a, b = find(i), find(j)
                if a != b:
                    parent[a] = b
    groups = collections.defaultdict(list)
    for i in range(len(cands)):
        groups[find(i)].append(i)
    taxcache = {}

    def taxon(s, t, contig):
        if (s, t) not in taxcache:
            taxcache[(s, t)] = bin_taxonomy(outdir, s, t)
        return taxcache[(s, t)].get(contig, "")
    rows = []
    for gi, idx in groups.items():
        mem = [cands[i] for i in idx]
        samples = {m[0] for m in mem}
        per_sample = collections.defaultdict(float)
        for m in mem:
            per_sample[m[0]] += m[8]
        rep = max(mem, key=lambda m: m[8])
        libs = collections.Counter(m[1] for m in mem)
        grp = collections.Counter(group_of(s) for s in samples)
        mf = [m[7] for m in mem if not math.isnan(m[7])]
        L = [len(m[6]) for m in mem]
        rows.append(dict(samples=len(samples), groups=" / ".join("%s %d" % (k, v) for k, v in sorted(grp.items())),
                         libraries=", ".join("%s %d" % (k.replace("sRNA-", ""), v) for k, v in sorted(libs.items())),
                         loci=len(mem), median_reads_per_sample=median(list(per_sample.values())),
                         reads=sum(per_sample.values()), median_length=median(L),
                         locus_class=top_label(collections.Counter(m[9] for m in mem)),
                         median_mfe=median(mf), median_mfe_per_nt=median([m[7] / len(m[6]) for m in mem
                                                                           if not math.isnan(m[7]) and m[6]]),
                         taxon="", rep_sample=rep[0], rep_type=rep[1],
                         rep_locus="%s:%d-%d(%s)" % (rep[2], rep[3], rep[4], rep[5]), rep_sequence=rep[6],
                         _rep=rep, _grp=dict(grp), _samples=set(samples)))
    rows.sort(key=lambda r: (-r["samples"], -r["reads"]))
    for i, r in enumerate(rows):
        r["cluster"] = "nsRNA_%04d" % (i + 1)
        if i < 200:                                # taxon of the representative locus (top clusters only)
            rep = r["_rep"]
            r["taxon"] = taxon(rep[0], rep[1], rep[2])
        del r["_rep"]
    return rows, cands


def novel_rfam(outdir, pairs, group_of):
    fam = {}
    for s, t in pairs:
        p = os.path.join(outdir, s, "Analysis", t, "novel", "%s.candidates.rfam.tbl" % s)
        if not os.path.exists(p):
            continue
        best = {}
        with open(p) as fh:
            for ln in fh:
                if ln.startswith("#"):
                    continue
                x = ln.split()
                if len(x) < 18:
                    continue
                q, sc = x[2], num(x[14])
                if q not in best or sc > best[q][1]:
                    best[q] = ((x[0], x[1], " ".join(x[17:])), sc)
        for q, ((name, acc, desc), sc) in best.items():
            if RFAM_SKIP.search(name) or RFAM_SKIP.search(desc):
                continue
            g = fam.setdefault((name, acc), dict(desc=desc, samples=set(), libs=collections.Counter(), loci=0))
            g["samples"].add(s)
            g["libs"][t] += 1
            g["loci"] += 1
    rows = []
    for (name, acc), g in fam.items():
        grp = collections.Counter(group_of(s) for s in g["samples"])
        rows.append(dict(family=name, accession=acc, description=g["desc"], samples=len(g["samples"]),
                         groups=" / ".join("%s %d" % (k, v) for k, v in sorted(grp.items())),
                         libraries=", ".join("%s %d" % (k.replace("sRNA-", ""), v) for k, v in sorted(g["libs"].items())),
                         loci=g["loci"], _grp=dict(grp), _samples=set(g["samples"])))
    rows.sort(key=lambda r: (-r["samples"], -r["loci"]))
    return rows


# ------------------------------------------------------------------ chimeras
def chimeras(outdir, samples, types):
    """(per-library tests, pooled candidate pairs) from module chimera, incl. Analysis/metaT/chimera"""
    lib, pool = [], {}
    for s in samples:
        for t in list(types) + ["metaT"]:
            d = os.path.join(outdir, s, "Analysis", t, "chimera")
            summ = {r["metric"]: r["value"] for r in read_tsv(os.path.join(d, "%s.chimera_summary.tsv" % s))}
            if not summ:
                continue
            en = {r["group"]: r for r in read_tsv(os.path.join(d, "%s.chimera_energy.tsv" % s))}
            bi = {r["group"]: r for r in read_tsv(os.path.join(d, "%s.chimera_bins.tsv" % s))}
            e, b = en.get("sRNA-like–mRNA", {}), bi.get("sRNA-like–mRNA", {})
            ea = en.get("all", {})
            lib.append(dict(sample=s, type=t, reads_used=num(summ.get("chimeric_reads_used")),
                            pairs_q005=num(summ.get("feature_pairs_q_lt_0.05")),
                            srna_mrna_pairs=num(e.get("pairs_tested")), dG=num(e.get("dG_median"), float("nan")),
                            dG_null=num(e.get("dG_null_median"), float("nan")),
                            p_complementarity=num(e.get("p_observed_lower_than_null"), float("nan")),
                            p_all=num(ea.get("p_observed_lower_than_null"), float("nan")),
                            same_bin=num(b.get("same_bin_fraction"), float("nan")),
                            same_bin_expected=num(b.get("expected_fraction"), float("nan"))))
            for r in read_tsv(os.path.join(d, "%s.chimera_pairs.tsv" % s)):
                q = num(r.get("q"), 1.0)
                if q >= 0.05:
                    continue
                na, nb = r.get("name_a") or r.get("type_a"), r.get("name_b") or r.get("type_b")
                key = tuple(sorted(((r.get("type_a"), na), (r.get("type_b"), nb))))
                g = pool.setdefault(key, dict(samples=set(), libs=collections.Counter(), reads=0.0, dG=[],
                                              pc=[], same=collections.Counter(), group=r.get("group"),
                                              taxa=collections.Counter()))
                g["samples"].add(s)
                g["libs"][t] += 1
                g["reads"] += num(r.get("reads"))
                if r.get("dG") not in (None, "", "NA"):
                    g["dG"].append(num(r.get("dG")))
                if r.get("p_complementarity") not in (None, "", "NA"):
                    g["pc"].append(num(r.get("p_complementarity")))
                g["same"][r.get("same_bin", "NA")] += 1
                g["taxa"][r.get("taxon_a") or ""] += 1
    rows = []
    for ((ta, na), (tb, nb)), g in pool.items():
        rows.append(dict(type_a=ta, name_a=na, type_b=tb, name_b=nb, group=g["group"], samples=len(g["samples"]),
                         libraries=", ".join("%s %d" % (k.replace("sRNA-", ""), v) for k, v in sorted(g["libs"].items())),
                         reads=g["reads"], median_dG=median(g["dG"]), median_p_complementarity=median(g["pc"]),
                         same_bin=top_label(g["same"]), taxon=top_label(g["taxa"])))
    rows.sort(key=lambda r: (r["group"] != "sRNA-like–mRNA", -r["samples"], -r["reads"]))
    return lib, rows
