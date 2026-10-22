#!/usr/bin/env python3
"""
srna_compare.py — per-sample comparison of sRNA types (e.g. sRNA-Ext, sRNA-Int),
metatranscriptomics (metaT) and metagenomics (metaG) on features, contigs and bins,
one folder per pair.

Every available pair of sources is compared on its own (e.g. sRNA-Ext_vs_sRNA-Int,
sRNA-Int_vs_metaT, sRNA-Int_vs_metaG, metaT_vs_metaG). Pair names have a fixed order:
sRNA-Ext, sRNA-Int, other sRNA types alphabetically, metaT, metaG; log2 ratios are
always log2(first / second), i.e. "per transcript" against metaT and "per DNA copy"
against metaG.

Per pair and level (feature, contig, bin):
  * kept: entries with >= --min-count raw counts in at least one of the TWO sources
  * CPM within each source (all entries of that level), log2((CPM_a+pc)/(CPM_b+pc))
  * features additionally: CPM within the feature type (CDS, tRNA, ...) and the
    log2 ratio of those — the meaningful ratio for pairs with metaT/metaG, because
    the read types are distributed very differently over the feature classes
  * summary: totals, shared entries (>= min count in both), Spearman of CPM;
    for features also per feature type, for bins also over HQ+MQ bins

Inputs
  --srna TYPE=<prefix>  srna_bins.py outputs <prefix>.feature_counts.tsv / .contig_counts.tsv
  --metat-features      featureCounts table(s) per gff feature (-g ID): file(s), IMP
                        Analysis/annotation directory or tar archive (mt.ID.<type>_counts...)
  --metag-features      same for metaG (mg.ID.<type>_counts...); default: the mg.ID tables
                        found in --metat-features
  --metat-contigs       metaT reads per contig: IMP Stats/mt/mgmt.assembly.contig_read_count.txt[.gz]
                        (contig start end count), a featureCounts table per contig, a 2-column
                        contig<TAB>count table, the IMP Stats directory or a tar of these
                        (default: summed from the feature tables)
  --metag-contigs       same for metaG; default: mg files found in --metat-contigs, else summed
                        from the metaG feature tables
  --gbed                gff2bed BED of the annotation (feature id -> contig, type)
  --omics-bin-basis     comma-separated bin count bases; the first one is the main bin level
                        (bins.tsv), each further one adds bins_<basis>.tsv tables + summary rows:
                          contig   all reads per contig, summed per bin (all sources)
                          cds      metaT/metaG: CDS feature counts only (mRNA / gene copies, no
                                   rRNA/tRNA reads); sRNA sources: all reads as in 'contig'
                          no_rrna  all sources: contigs carrying an rRNA feature left out
  --contig2bin, --gtdbtk  bins, their GTDB taxonomy and CheckM quality

Outputs (<outdir>/)
  <A>_vs_<B>/<sample>.<A>_vs_<B>.{features,contigs,bins,summary}.tsv
  <sample>.all_sources.{features,contigs,bins}.tsv   counts + CPM of all sources
  <sample>.compare.summary.tsv                        all metrics, all pairs
"""
import argparse
import glob
import io
import tarfile
import gzip
import collections
import itertools
import math
import os
import re
import sys


def open_in(path):
    """open a text input, transparently gunzipping *.gz"""
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)

FT_NONE = "unknown_type"


def strip_ext(name):
    return re.sub(r"\.(fa|fna|fasta|fsa)(\.gz)?$", "", name)


def source_order(s):
    if s == "sRNA-Ext":
        return (0, s)
    if s == "sRNA-Int":
        return (1, s)
    if s == "metaT":
        return (3, s)
    if s == "metaG":
        return (4, s)
    return (2, s)


# final IMP tables only (mt.ID.<type>_counts.tsv[.gz]); intermediate *.paired.tsv / *.single.tsv
# files, if still present, would count reads twice and are therefore not matched
IMP_RE = re.compile(r"(?:^|/)(mt|mg)\.ID\.([^/]+?)(?:_counts|_per_contig\.counts)\.tsv(?:\.gz)?$")
# feature types whose metaT counts are summed per contig when no metaT contig table is given
# (nested / overlapping types such as crispr-repeat inside CRISPR are left out)
CONTIG_TYPES = {"CDS", "tRNA", "rRNA", "ncRNA", "tmRNA", "CRISPR"}


def parse_featurecounts(fh, counts, meta, ftype=""):
    """featureCounts output: '#' comment, header Geneid Chr Start End Strand Length <count columns...>;
    all count columns (e.g. IMP's paired + single BAM) are summed"""
    for line in fh:
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if f[0] in ("Geneid", "gene_id", "ID") or len(f) < 7:
            continue
        try:
            c = sum(float(x) for x in f[6:])
        except ValueError:
            continue
        counts[f[0]] = counts.get(f[0], 0.0) + c
        meta.setdefault(f[0], (f[1].split(";")[0], f[5], ftype))


def read_featurecounts(spec, omic="mt"):
    """counts of one omic (IMP prefix 'mt' = metaT, 'mg' = metaG) from featureCounts tables.
    spec = comma-separated list of
      - featureCounts tables (.tsv / .gz)
      - directories (IMP Analysis/annotation or its per_contig_counts/): <omic>.ID.<type>_counts.tsv[.gz]
        or <omic>.ID.<type>_per_contig.counts.tsv[.gz]
      - tar archives (.tar, .tar.gz / gzipped tar) with such members
    The feature type comes from the IMP file name; a type found in several files (both naming schemes)
    is read only once. A plain table that does not follow the IMP naming is read for omic 'mt' only."""
    counts, meta, used, seen = {}, {}, [], set()

    def take(name, opener):
        m = IMP_RE.search(name)
        if m:
            if m.group(1) != omic or m.group(2) in seen:
                return
            seen.add(m.group(2))
            ftype = m.group(2)
        else:
            if omic != "mt":
                return
            ftype = ""
        with opener() as fh:
            parse_featurecounts(fh, counts, meta, ftype)
        used.append(name)

    for item in [x for x in spec.split(",") if x]:
        if os.path.isdir(item):
            for f in sorted(glob.glob(os.path.join(item, "%s.ID.*.tsv*" % omic))):
                if IMP_RE.search(f):
                    take(f, lambda f=f: open_in(f))
        elif tarfile.is_tarfile(item):
            with tarfile.open(item, "r:*") as tar:
                for mem in sorted(tar.getmembers(), key=lambda x: x.name):
                    if not (mem.isfile() and IMP_RE.search(mem.name)):
                        continue

                    def opener(mem=mem):
                        raw = tar.extractfile(mem).read()
                        if mem.name.endswith(".gz"):
                            raw = gzip.decompress(raw)
                        return io.StringIO(raw.decode())
                    take(item + ":" + mem.name, opener)
        else:
            take(item, lambda item=item: open_in(item))
    return counts, meta, used


CONTIG_FILE_RE = re.compile(r"contig_read_count[^/]*\.txt(?:\.gz)?$|contig_counts?[^/]*\.tsv(?:\.gz)?$")
CONTIG_OMIC_RE = re.compile(r"(?:^|/)(mt|mg)(?:/|\.)[^/]*$")


def parse_contig_counts(fh, counts, lengths):
    """per-contig counts, format detected per line:
    featureCounts (Geneid Chr Start End Strand Length counts...), IMP contig_read_count
    (contig start end count) or contig<TAB>count"""
    for line in fh:
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        try:
            if len(f) >= 7:
                c, L = sum(float(x) for x in f[6:]), f[5]
            elif len(f) == 4:
                c, L = float(f[3]), str(int(f[2]) - int(f[1]))
            elif len(f) == 2:
                c, L = float(f[1]), ""
            else:
                continue
        except ValueError:        # header line
            continue
        counts[f[0]] = counts.get(f[0], 0.0) + c
        if L:
            lengths.setdefault(f[0], L)


def read_contig_counts(spec, omic="mt", tagged_only=False):
    """contig counts of one omic. spec = comma-separated files / directories / tar archives.
    A file whose path says mt/ or mg/ (IMP Stats/<omic>/..., or <omic>.<name>) is used only for that
    omic; untagged files are used unless tagged_only."""
    counts, lengths, used = {}, {}, []

    def ok(name):
        m = CONTIG_OMIC_RE.search(name)
        return (m.group(1) == omic) if m else not tagged_only

    for item in [x for x in spec.split(",") if x]:
        if os.path.isdir(item):
            for f in sorted(glob.glob(os.path.join(item, "**", "*"), recursive=True)):
                if os.path.isfile(f) and CONTIG_FILE_RE.search(f) and ok(os.path.relpath(f, item)):
                    with open_in(f) as fh:
                        parse_contig_counts(fh, counts, lengths)
                    used.append(f)
        elif tarfile.is_tarfile(item):
            with tarfile.open(item, "r:*") as tar:
                for mem in sorted(tar.getmembers(), key=lambda x: x.name):
                    if not (mem.isfile() and CONTIG_FILE_RE.search(mem.name) and ok(mem.name)):
                        continue
                    raw = tar.extractfile(mem).read()
                    if mem.name.endswith(".gz"):
                        raw = gzip.decompress(raw)
                    parse_contig_counts(io.StringIO(raw.decode()), counts, lengths)
                    used.append(item + ":" + mem.name)
        elif ok(item):
            with open_in(item) as fh:
                parse_contig_counts(fh, counts, lengths)
            used.append(item)
    return counts, lengths, used


def read_table(path, key, val, extra=()):
    d, m = {}, {}
    with open_in(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ki, vi = header.index(key), [header.index(v) for v in val]
        ei = [header.index(e) for e in extra]
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header):
                continue
            d[f[ki]] = sum(float(f[i]) for i in vi)
            m[f[ki]] = [f[i] for i in ei]
    return d, m


def read_gbed(path):
    """gff2bed: chrom start end id score strand source type phase attributes"""
    meta = {}
    if path and os.path.exists(path):
        with open_in(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 8:
                    meta.setdefault(f[3], (f[0], f[7]))
    return meta


def read_contig2bin(path):
    d = {}
    if path and os.path.exists(path):
        with open_in(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 2 and not line.startswith("#"):
                    d[f[0]] = strip_ext(f[1])
    return d


def read_gtdbtk(paths):
    tax = {}
    for pth in [p for p in paths.split(",") if p]:
        with open_in(pth) as fh:
            header = None
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if header is None or f[0] == "user_genome":
                    header = f
                    continue
                rec = dict(zip(header, f))
                tax[strip_ext(rec.get("user_genome", f[0]))] = rec.get("classification", "")
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


def cpm(d, keys=None):
    keys = d.keys() if keys is None else keys
    t = sum(d.get(k, 0.0) for k in keys)
    return {k: (d.get(k, 0.0) * 1e6 / t if t else 0.0) for k in keys}, t


def spearman(x, y):
    n = len(x)
    if n < 3:
        return float("nan")

    def ranks(v):
        order = sorted(range(n), key=lambda i: v[i])
        r = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and v[order[j + 1]] == v[order[i]]:
                j += 1
            for k in range(i, j + 1):
                r[order[k]] = (i + j) / 2.0 + 1
            i = j + 1
        return r
    rx, ry = ranks(x), ranks(y)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((p - mx) * (q - my) for p, q in zip(rx, ry))
    den = math.sqrt(sum((p - mx) ** 2 for p in rx) * sum((q - my) ** 2 for q in ry))
    return num / den if den else float("nan")


def fmt(v, nd=3):
    return ("%%.%df" % nd) % v


def write_all_sources(path, level, sources, data, norm, meta_cols, meta, min_count):
    keys = set().union(*[set(data[s]) for s in sources])
    keys = [k for k in keys if max(data[s].get(k, 0.0) for s in sources) >= min_count]
    keys.sort(key=lambda k: -max(norm[s].get(k, 0.0) for s in sources))
    with open(path, "w") as out:
        out.write("\t".join([level] + meta_cols + ["%s_count" % s for s in sources] +
                            ["%s_CPM" % s for s in sources]) + "\n")
        for k in keys:
            row = [k] + list(meta.get(k, [""] * len(meta_cols)))
            row += [fmt(data[s].get(k, 0.0)) for s in sources]
            row += [fmt(norm[s].get(k, 0.0)) for s in sources]
            out.write("\t".join(row) + "\n")


def write_pair(path, level, a, b, data, norm, meta_cols, meta, pc, min_count, summary, pair,
               ftype=None, norm_t=None):
    keys = set(data[a]) | set(data[b])
    keys = [k for k in keys if max(data[a].get(k, 0.0), data[b].get(k, 0.0)) >= min_count]
    keys.sort(key=lambda k: -max(norm[a].get(k, 0.0), norm[b].get(k, 0.0)))
    cols = [level] + meta_cols + ["%s_count" % a, "%s_count" % b, "%s_CPM" % a, "%s_CPM" % b,
                                  "log2_%s_vs_%s" % (a, b)]
    if ftype is not None:
        cols += ["%s_CPM_within_type" % a, "%s_CPM_within_type" % b, "log2_within_type_%s_vs_%s" % (a, b)]
    with open(path, "w") as out:
        out.write("\t".join(cols) + "\n")
        for k in keys:
            row = [k] + list(meta.get(k, [""] * len(meta_cols)))
            ca, cb = norm[a].get(k, 0.0), norm[b].get(k, 0.0)
            row += [fmt(data[a].get(k, 0.0)), fmt(data[b].get(k, 0.0)), fmt(ca), fmt(cb),
                    fmt(math.log2((ca + pc) / (cb + pc)), 4)]
            if ftype is not None:
                ta, tb = norm_t[a].get(k, 0.0), norm_t[b].get(k, 0.0)
                row += [fmt(ta), fmt(tb), fmt(math.log2((ta + pc) / (tb + pc)), 4)]
            out.write("\t".join(row) + "\n")

    summary.append((pair, level, "all", "entries_kept", str(len(keys))))
    shared = [k for k in keys if data[a].get(k, 0.0) >= min_count and data[b].get(k, 0.0) >= min_count]
    summary.append((pair, level, "all", "shared_n", str(len(shared))))
    summary.append((pair, level, "all", "spearman_cpm",
                    fmt(spearman([norm[a][k] for k in shared], [norm[b][k] for k in shared]), 4)))
    if ftype is not None:
        by_type = collections.defaultdict(list)
        for k in shared:
            by_type[ftype.get(k, FT_NONE)].append(k)
        tot = collections.defaultdict(lambda: [0.0, 0.0])
        for k in set(data[a]) | set(data[b]):
            t = ftype.get(k, FT_NONE)
            tot[t][0] += data[a].get(k, 0.0)
            tot[t][1] += data[b].get(k, 0.0)
        for t in sorted(tot):
            summary.append((pair, level, t, "%s_count" % a, fmt(tot[t][0], 1)))
            summary.append((pair, level, t, "%s_count" % b, fmt(tot[t][1], 1)))
            ks = by_type.get(t, [])
            summary.append((pair, level, t, "shared_n", str(len(ks))))
            summary.append((pair, level, t, "spearman_cpm_within_type",
                            fmt(spearman([norm_t[a][k] for k in ks], [norm_t[b][k] for k in ks]), 4)))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--srna", action="append", default=[], help="TYPE=<prefix> of srna_bins.py outputs")
    ap.add_argument("--metat-features", default="")
    ap.add_argument("--metat-contigs", default="")
    ap.add_argument("--metag-contigs", default="",
                    help="metaG contig counts (default: mg files found in --metat-contigs)")
    ap.add_argument("--metag-features", default="",
                    help="metaG feature counts (default: mg.ID.* tables found in --metat-features)")
    ap.add_argument("--gbed", default="", help="gff2bed BED of the annotation (feature types)")
    ap.add_argument("--omics-bin-basis", default="contig,cds,no_rrna",
                    help="bin count bases (contig, cds, no_rrna); first = main bin level")
    ap.add_argument("--contig2bin", default="")
    ap.add_argument("--gtdbtk", default="")
    ap.add_argument("--pseudocount", type=float, default=1.0, help="added to CPM before log2 ratios")
    ap.add_argument("--min-count", type=float, default=5.0,
                    help="min raw count in at least one of the two sources of a pair")
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    feat, cont, fmeta, cmeta = {}, {}, {}, {}
    gbed = read_gbed(a.gbed)
    c2b = read_contig2bin(a.contig2bin)
    for spec in a.srna:
        t, pref = spec.split("=", 1)
        feat[t], m = read_table(pref + ".feature_counts.tsv", "feature_id", ["sense_reads"],
                                ["contig", "type", "bin"])
        fmeta.update(m)
        cont[t], m = read_table(pref + ".contig_counts.tsv", "contig", ["reads"], ["length", "bin"])
        cmeta.update(m)
    metat_notes = []
    ctab = {}
    for src, omic, spec, tagged in (("metaT", "mt", a.metat_contigs, False),
                                    ("metaG", "mg", a.metag_contigs, False),
                                    ("metaG", "mg", a.metat_contigs, True)):
        if not spec or src in ctab:
            continue
        cnt, lens, used = read_contig_counts(spec, omic, tagged)
        if not used:
            continue
        ctab[src] = cnt
        cont[src] = cnt
        metat_notes.append(("%s_contigs" % src, "table:" + ",".join("/".join(u.replace(":", "/").split("/")[-2:]) for u in used)))
        sys.stderr.write("[compare] %s: %s contig table(s): %s\n" % (a.sample, src, ", ".join(used)))
        for k in cnt:
            cmeta.setdefault(k, [lens.get(k, ""), c2b.get(k, "unbinned")])
    omics = []
    if a.metat_features:
        omics.append(("metaT", "mt", a.metat_features))
    if a.metag_features or a.metat_features:
        omics.append(("metaG", "mg", a.metag_features or a.metat_features))
    for src, omic, spec in omics:
        cnt, mm, used = read_featurecounts(spec, omic)
        if not used:
            continue
        feat[src] = cnt
        metat_notes.append(("%s_feature_tables" % src, str(len(used))))
        sys.stderr.write("[compare] %s: %s feature tables: %s\n" % (a.sample, src, ", ".join(used)))
        for k, (chrom, _L, t) in mm.items():
            fmeta.setdefault(k, [chrom, t or gbed.get(k, ("", ""))[1], c2b.get(chrom, "unbinned")])
        if src in ctab:
            continue
        # no contig table: feature counts summed per contig (sRNA contig counts also contain
        # intergenic reads, so compare contig levels with care)
        agg = collections.defaultdict(float)
        for k, v in cnt.items():
            chrom, _L, t = mm[k]
            if (t or gbed.get(k, ("", ""))[1]) in CONTIG_TYPES:
                agg[chrom] += v
        cont[src] = dict(agg)
        metat_notes.append(("%s_contigs" % src, "derived_from_features(%s)" % ",".join(sorted(CONTIG_TYPES))))
        for k in agg:
            cmeta.setdefault(k, ["", c2b.get(k, "unbinned")])
    # CDS flagged by srna_bins.py (opposite strand of an rRNA gene rescued by cmscan): spurious ORFs whose sRNA,
    # metaT and metaG reads are rRNA reads — left out of all feature comparisons
    flagged = set()
    for spec in a.srna:
        fp = spec.split("=", 1)[1] + ".feature_counts.tsv"
        with open_in(fp) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            if "flag" in hdr:
                fi, ki = hdr.index("flag"), hdr.index("feature_id")
                for ln in fh:
                    f = ln.rstrip("\n").split("\t")
                    if len(f) > fi and f[fi]:
                        flagged.add(f[ki])
    if flagged:
        for src in feat:
            for k in flagged:
                feat[src].pop(k, None)
        metat_notes.append(("flagged_features_excluded", str(len(flagged))))
        sys.stderr.write("[compare] %s: %d CDS opposite rescued rRNA genes left out\n" % (a.sample, len(flagged)))
    # complete feature meta from the annotation
    for k, (chrom, t) in gbed.items():
        m = fmeta.get(k)
        if m is None:
            fmeta[k] = [chrom, t, c2b.get(chrom, "unbinned")]
        elif not m[1]:
            m[1] = t
    ftype = {k: (m[1] or FT_NONE) for k, m in fmeta.items()}

    fsrc = sorted(feat, key=source_order)
    csrc = sorted(cont, key=source_order)
    if len(fsrc) < 2 and len(csrc) < 2:
        sys.stderr.write("[compare] %s: fewer than two sources — nothing to compare\n" % a.sample)
        sys.exit(0)

    summary = [("all", "omics", "all", k, v) for k, v in metat_notes]
    srna_ids = set().union(*[set(feat[s]) for s in fsrc if s not in ("metaT", "metaG")]) if fsrc else set()
    for om in ("metaT", "metaG"):
        if om in feat and srna_ids:
            hit = len(srna_ids & set(feat[om]))
            summary.append(("all", "feature", "all", "srna_feature_ids_found_in_%s" % om, "%d/%d" % (hit, len(srna_ids))))
            if hit == 0:
                sys.stderr.write("[compare] WARNING: no sRNA feature ID found in the %s feature tables — "
                                 "check that featureCounts used the same ID attribute (-g ID)\n" % om)

    # normalisation within each source: over all entries, and within feature type
    fnorm, fnorm_t = {}, {}
    for s in fsrc:
        fnorm[s], tot = cpm(feat[s])
        summary.append(("all", "feature", "all", "%s_total_counts" % s, fmt(tot, 1)))
        groups = collections.defaultdict(list)
        for k in feat[s]:
            groups[ftype.get(k, FT_NONE)].append(k)
        fnorm_t[s] = {}
        for t, ks in groups.items():
            fnorm_t[s].update(cpm(feat[s], ks)[0])
    cnorm = {}
    for s in csrc:
        cnorm[s], tot = cpm(cont[s])
        summary.append(("all", "contig", "all", "%s_total_counts" % s, fmt(tot, 1)))
    tax = read_gtdbtk(a.gtdbtk) if a.gtdbtk else {}
    bases = [x for x in a.omics_bin_basis.split(",") if x]
    for x in bases:
        if x not in ("contig", "cds", "no_rrna"):
            sys.exit("[compare] unknown --omics-bin-basis '%s' (contig, cds, no_rrna)" % x)
    if not bases:
        bases = ["contig"]

    def bin_of(c):
        return c2b.get(c) or (cmeta.get(c, ["", "unbinned"])[1] or "unbinned")

    rrna_contigs = {m[0] for k, m in fmeta.items() if m[1] == "rRNA"}
    rrna_contigs |= {c for k, (c, t) in gbed.items() if t == "rRNA"}

    def bin_counts(basis):
        bc, bn, note = {}, {}, {}
        for s in csrc:
            agg = collections.defaultdict(float)
            if basis == "cds" and s in ("metaT", "metaG") and s in feat:
                for k, v in feat[s].items():
                    if ftype.get(k) == "CDS":
                        agg[bin_of(fmeta[k][0])] += v
                note[s] = "CDS_features"
            else:
                for c, v in cont[s].items():
                    if basis == "no_rrna" and c in rrna_contigs:
                        continue
                    agg[bin_of(c)] += v
                note[s] = "contigs" + ("_without_rRNA" if basis == "no_rrna" else "")
            bc[s] = dict(agg)
            bn[s] = cpm(bc[s])[0]
        return bc, bn, note

    bins_by_basis = {}
    for x in bases:
        bins_by_basis[x] = bin_counts(x)
        for s, n in bins_by_basis[x][2].items():
            summary.append(("all", "bin" if x == bases[0] else "bin_" + x, "all", "%s_basis" % s, n))
    summary.append(("all", "omics", "all", "rRNA_contigs", str(len(rrna_contigs))))
    binc, bnorm, _ = bins_by_basis[bases[0]]
    qual = read_binstats(a.gtdbtk) if a.gtdbtk else {}
    bmeta = {b: [tax.get(b, "")] + list(qual.get(b, ("", "", "") if b == "unbinned" else ("NA", "NA", "NA")))
             for x in bases for s in csrc for b in bins_by_basis[x][0][s]}
    BCOLS = ["classification", "completeness", "contamination", "quality"]

    p = os.path.join(a.outdir, a.sample)
    if len(fsrc) >= 2:
        write_all_sources(p + ".all_sources.features.tsv", "feature", fsrc, feat, fnorm,
                          ["contig", "type", "bin"], fmeta, a.min_count)
    if len(csrc) >= 2:
        write_all_sources(p + ".all_sources.contigs.tsv", "contig", csrc, cont, cnorm,
                          ["length", "bin"], cmeta, a.min_count)
        for x in bases:
            sfx = "" if x == bases[0] else "_" + x
            write_all_sources(p + ".all_sources.bins%s.tsv" % sfx, "bin", csrc, bins_by_basis[x][0],
                              bins_by_basis[x][1], BCOLS, bmeta, a.min_count)

    pairs = []
    for x, y in itertools.combinations(sorted(set(fsrc) | set(csrc), key=source_order), 2):
        pair = "%s_vs_%s" % (x, y)
        d = os.path.join(a.outdir, pair)
        done = False
        if x in feat and y in feat:
            os.makedirs(d, exist_ok=True)
            write_pair(os.path.join(d, "%s.%s.features.tsv" % (a.sample, pair)), "feature", x, y, feat, fnorm,
                       ["contig", "type", "bin"], fmeta, a.pseudocount, a.min_count, summary, pair,
                       ftype=ftype, norm_t=fnorm_t)
            done = True
        if x in cont and y in cont:
            os.makedirs(d, exist_ok=True)
            write_pair(os.path.join(d, "%s.%s.contigs.tsv" % (a.sample, pair)), "contig", x, y, cont, cnorm,
                       ["length", "bin"], cmeta, a.pseudocount, a.min_count, summary, pair)
            for basis in bases:
                bc, bn, _ = bins_by_basis[basis]
                omic_pair = x in ("metaT", "metaG") or y in ("metaT", "metaG")
                if basis == "cds" and not omic_pair:
                    continue          # identical to 'contig' for sRNA-only pairs
                lvl = "bin" if basis == bases[0] else "bin_" + basis
                sfx = "" if basis == bases[0] else "_" + basis
                write_pair(os.path.join(d, "%s.%s.bins%s.tsv" % (a.sample, pair, sfx)), lvl, x, y, bc, bn,
                           BCOLS, bmeta, a.pseudocount, a.min_count, summary, pair)
                good = [b for b in set(bc[x]) | set(bc[y])
                        if bmeta.get(b, ["", "", "", ""])[3] in ("HQ", "MQ")
                        and bc[x].get(b, 0.0) >= a.min_count and bc[y].get(b, 0.0) >= a.min_count]
                if good:
                    summary.append((pair, lvl, "HQ+MQ", "shared_n", str(len(good))))
                    summary.append((pair, lvl, "HQ+MQ", "spearman_cpm",
                                    fmt(spearman([bn[x][b] for b in good], [bn[y][b] for b in good]), 4)))
            done = True
        if done:
            pairs.append(pair)
            with open(os.path.join(d, "%s.%s.summary.tsv" % (a.sample, pair)), "w") as out:
                out.write("pair\tlevel\tfeature_type\tmetric\tvalue\n")
                for row in summary:
                    if row[0] == pair:
                        out.write("\t".join(row) + "\n")
        else:
            sys.stderr.write("[compare] %s: %s skipped — the two sources share no level "
                             "(features vs contigs only)\n" % (a.sample, pair))

    with open(p + ".compare.summary.tsv", "w") as out:
        out.write("pair\tlevel\tfeature_type\tmetric\tvalue\n")
        for row in summary:
            out.write("\t".join(row) + "\n")
    sys.stderr.write("[compare] %s: sources %s; pairs %s\n"
                     % (a.sample, ",".join(sorted(set(fsrc) | set(csrc), key=source_order)), ",".join(pairs)))


if __name__ == "__main__":
    main()
