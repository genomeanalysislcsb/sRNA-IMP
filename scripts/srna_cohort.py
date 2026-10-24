#!/usr/bin/env python3
"""
srna_cohort.py — cross-sample (cohort) summary of the sRNA-IMP pipeline

Collects the per-sample summary tables of all samples under <outdir>, joins them with a metadata
table (groups, e.g. HC / PD / RA or HC / T1D, plus covariates like sex, age, constipation) and writes

  <cohort-dir>/tables/
    cohort.sample_stats.tsv          one row per sample and sRNA type: key numbers, integrity,
                                     wet-lab QC, CRISPR, sRNA-vs-omics correlations + metadata
    cohort.<type>.<matrix>.counts.tsv feature x sample count matrices (miRNA, genus, family,
                                     species, Rfam type, tRNA amino acid) for downstream
                                     differential analysis (also used by the adjusted models, --model)
    cohort.metadata.tsv              metadata of the matched samples (+ unmatched IDs listed)
  <cohort-dir>/stats/
    cohort.<type>.metrics_by_group.tsv   Kruskal-Wallis + pairwise Mann-Whitney, BH-FDR
    cohort.<type>.<matrix>.by_group.tsv  per feature (prevalence filtered), BH-FDR
    cohort.<type>.<matrix>.permanova.tsv Bray-Curtis PERMANOVA (group and each categorical covariate)
    cohort.paired_Int_Ext.tsv            paired Wilcoxon of the metrics, sRNA-Int vs sRNA-Ext
    cohort.covariates_by_group.tsv       balance of covariates between groups (confounding check)
    cohort.integrity_vs_qc.tsv           Spearman of integrity indicators vs wet-lab QC (RIN, ...)
  <cohort-dir>/cohort.report.pdf (+ .html)

Group tests are exploratory and unadjusted; --model adds covariate-adjusted linear models (OLS on CLR
abundances, HC3 robust SE). Only python3 + numpy + matplotlib are needed (no R).

Usage
  srna_cohort.py --outdir imp --metadata metadata.tsv [--id-col sample] [--group-col group]
                 [--covariates sex,age,constipation] [--samples samples.tsv] [--rna-qc rna_qc.tsv]
                 [--cohort-dir imp/Cohort] [--html]
"""
import argparse
import collections
import csv
import glob
import math
import os
import gzip
import random
import re
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import srna_report as R            # noqa: E402  (page layout, tables, per-sample loader)
import srna_qc as Q                # noqa: E402  (library metrics, QC rules)
import srna_highlights as HL       # noqa: E402  (cross-sample highlights: CRISPR, novel sRNAs, chimeras)
import srna_heterogeneity as HT    # noqa: E402  (dispersion, patients outside the control range)

MISSING = {"", "na", "n/a", "nan", "none", "null", "-", "?"}
# functional categories of the sRNA reads on CDS (assembly module, functions.tsv): db -> (term filter, label)
FUNC_DBS = {
    "kegg_pathway": (lambda x: re.match(r"ko\d{5}$", x) and not re.match(r"ko01[12]\d\d$", x), "KEGG pathway"),
    "kegg_module": (lambda x: re.match(r"M\d{5}$", x), "KEGG module"),
    "kegg_ko": (lambda x: re.match(r"K\d{5}$", x), "KEGG orthology (KO)"),
    "pfam": (lambda x: re.match(r"PF\d{5}$", x), "Pfam domain"),
    "cog": (lambda x: re.match(r"COG\d{4}$", x), "COG"),
}
FUNC_MX = tuple("fn_" + k for k in FUNC_DBS)              # sense reads (mRNA-derived fragments)
FUNC_AS = tuple("fnas_" + k for k in FUNC_DBS)            # antisense reads (exported, not tested)
MATRICES = ("mirna", "genus", "family", "species", "rfam_type", "trna_aa", "rna_class", "rfam_family") + FUNC_MX + FUNC_AS
# RNA classes of the non-human reads: (class, Rfam type prefix); first match wins
RNA_CLASSES = ["rRNA", "tRNA", "miRNA-like", "snoRNA", "snRNA", "sRNA", "riboswitch", "leader", "other cis-reg",
               "ribozyme", "intron", "antisense", "lncRNA", "CRISPR", "antitoxin", "tmRNA", "SRP RNA", "6S RNA",
               "RNase P", "other ncRNA gene", "unknown"]
# house-keeping ncRNAs of Rfam type 'Gene', recognised by their Rfam id
RFAM_ID_CLASS = [(r"^tmRNA", "tmRNA"), (r"SRP", "SRP RNA"), (r"^6S$", "6S RNA"), (r"RNaseP|RNase_P", "RNase P")]
RFAM_TYPE_CLASS = [("rRNA", "rRNA"), ("tRNA", "tRNA"), ("miRNA", "miRNA-like"), ("snRNA;snoRNA", "snoRNA"),
                   ("snRNA", "snRNA"), ("sRNA", "sRNA"), ("Cis-reg;riboswitch", "riboswitch"),
                   ("Cis-reg;leader", "leader"), ("Cis-reg", "other cis-reg"), ("ribozyme", "ribozyme"),
                   ("Intron", "intron"), ("antisense", "antisense"), ("lncRNA", "lncRNA"), ("CRISPR", "CRISPR"),
                   ("antitoxin", "antitoxin"), ("Gene", "other ncRNA gene")]
CLASS_COLORS = dict(zip(RNA_CLASSES, ["#8c564b", "#1f77b4", "#e377c2", "#2ca02c", "#98df8a", "#ff7f0e", "#d62728",
                                      "#ff9896", "#c49c94", "#9467bd", "#c5b0d5", "#17becf", "#9edae5", "#bcbd22",
                                      "#dbdb8d", "#393b79", "#6b6ecf", "#9c9ede", "#637939", "#7f7f7f", "#dddddd"]))


def rna_classes(d):
    """non-human reads per RNA class: rRNA / tRNA / miRNA-like / unknown from the classification, the Rfam hits
    split by their Rfam type (rRNA and tRNA hits added to those classes)"""
    out = collections.OrderedDict((c, 0.0) for c in RNA_CLASSES)
    nh = d.nonhuman
    out["rRNA"] += R.num(nh.get("rRNA_reads"))
    out["tRNA"] += R.num(nh.get("tRNA"))
    out["miRNA-like"] += R.num(nh.get("miRNA_like"))
    out["unknown"] += R.num(nh.get("unknown"))
    fam = getattr(d, "_rfam_family", None)
    if fam:
        for r in fam:
            ty = (r.get("type") or "").strip()
            rid = (r.get("rfam_id") or "").strip()
            cl = next((c for pat, c in RFAM_ID_CLASS if ty.startswith("Gene") and re.search(pat, rid)), None) or \
                next((c for pre, c in RFAM_TYPE_CLASS if ty.startswith(pre)), "other ncRNA gene")
            out[cl] += R.num(r.get("count"))
    else:                                     # no family table: Rfam reads stay one block
        out["other ncRNA gene"] += R.num(nh.get("Rfam_classified"))
    return out
# metrics checked for robust outliers within each fraction (WARN; --outlier-z, 0 = off): metric -> minimum
# absolute distance from the median that is relevant (a robust z alone flags harmless differences when the
# spread is tiny, e.g. 0.3 % vs 0.06 % adapter dimers)
OUTLIER_METRICS = {"trimmed_fraction": 0.10, "adapter_dimer_fraction": 0.02, "human_fraction": 0.002,
                   "assembly_mapped_fraction": 0.05, "nonhuman_rRNA_fraction": 0.10, "nonhuman_tRNA_fraction": 0.10,
                   "nonhuman_unknown_fraction": 0.10, "integrity_tRNA_internal_fraction": 0.10,
                   "integrity_tRNA_3p_CCA_fraction": 0.10, "rho_vs_metaG_bin": 0.15}
# bookkeeping / annotation counts that are not tested
NOT_TESTED = ("missing_modules", "missing_compare", "crispr_arrays", "antisense_rRNA_filtered")
TERM_NAMES = {}       # term -> name, from --term-names


def tname(term, n=60):
    nm = TERM_NAMES.get(term) or TERM_NAMES.get(term.replace("ko", "map", 1))
    return R.short("%s %s" % (term, nm), n) if nm else term


def load_term_names(paths):
    """term -> name from comma-separated files; formats detected per file:
       - KEGG flat files (ko, pathway, module: ENTRY / NAME / DEFINITION blocks; KO: DEFINITION, else NAME)
       - TSV term<TAB>name (e.g. KEGG REST list/ko, list/pathway, list/module; prefixes ko:, path:, md: removed)
       - Pfam-A.clans.tsv(.gz): accession, clan, clan name, short name, description -> description"""
    n0 = len(TERM_NAMES)
    for p in [x.strip() for x in (paths or "").split(",") if x.strip()]:
        if not os.path.exists(p):
            sys.stderr.write("[cohort] term names %s not found\n" % p)
            continue
        op = gzip.open if p.endswith(".gz") else open
        n1 = len(TERM_NAMES)
        with op(p, "rt", errors="replace") as fh:
            first = fh.readline()
            fh.seek(0)
            if first.startswith("ENTRY"):
                entry, name, defin = None, "", ""
                for ln in fh:
                    key = ln[:12].strip()
                    val = ln[12:].strip()
                    if key == "ENTRY":
                        entry, name, defin = val.split()[0], "", ""
                    elif key == "NAME" and entry and not name:
                        name = val
                    elif key == "DEFINITION" and entry and not defin:
                        defin = val
                    elif ln.startswith("///") and entry:
                        label = defin if entry.startswith("K") and defin else name or defin
                        if label:
                            TERM_NAMES.setdefault(entry, re.sub(r"\s+", " ", label))
                        entry = None
            else:
                for ln in fh:
                    if ln.startswith("#"):
                        continue
                    f = ln.rstrip("\n").split("\t")
                    if len(f) < 2:                    # not tab-separated: id, whitespace, name
                        f = ln.strip().split(None, 1)
                    if len(f) < 2 or not f[0]:
                        continue
                    term = re.sub(r"^(ko|path|md|map):", "", f[0].strip())
                    name = f[4] if len(f) >= 5 and term.startswith("PF") else f[1]
                    TERM_NAMES.setdefault(term, name.strip())
        sys.stderr.write("[cohort] %d term names from %s\n" % (len(TERM_NAMES) - n1, p))
    return len(TERM_NAMES) - n0


# ----------------------------------------------------------------------------- statistics
def rankdata(x):
    x = np.asarray(x, float)
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(len(x))
    xs = x[order]
    i = 0
    while i < len(x):
        j = i
        while j + 1 < len(x) and xs[j + 1] == xs[i]:
            j += 1
        ranks[order[i:j + 1]] = (i + j) / 2.0 + 1
        i = j + 1
    return ranks


def _tie_term(r):
    _, c = np.unique(r, return_counts=True)
    return float(np.sum(c ** 3 - c))


def _gammq(a, x):
    """regularised upper incomplete gamma Q(a, x) (Numerical Recipes)"""
    if x <= 0:
        return 1.0
    gln = math.lgamma(a)
    if x < a + 1:
        ap, s, d = a, 1.0 / a, 1.0 / a
        for _ in range(500):
            ap += 1
            d *= x / ap
            s += d
            if abs(d) < abs(s) * 1e-12:
                break
        return max(0.0, 1.0 - s * math.exp(-x + a * math.log(x) - gln))
    b = x + 1 - a
    c, d = 1e300, 1.0 / b
    h = d
    for i in range(1, 500):
        an = -i * (i - a)
        b += 2
        d = an * d + b
        d = 1e-300 if abs(d) < 1e-300 else d
        c = b + an / c
        c = 1e-300 if abs(c) < 1e-300 else c
        d = 1.0 / d
        de = d * c
        h *= de
        if abs(de - 1) < 1e-12:
            break
    return min(1.0, math.exp(-x + a * math.log(x) - gln) * h)


def kruskal(groups):
    groups = [np.asarray(g, float) for g in groups if len(g) > 0]
    if len(groups) < 2:
        return float("nan"), float("nan")
    allv = np.concatenate(groups)
    n = len(allv)
    if n < 3 or np.all(allv == allv[0]):
        return float("nan"), float("nan")
    r = rankdata(allv)
    h, i = 0.0, 0
    for g in groups:
        rg = r[i:i + len(g)]
        h += len(g) * (rg.mean() - (n + 1) / 2.0) ** 2
        i += len(g)
    h *= 12.0 / (n * (n + 1))
    c = 1 - _tie_term(r) / (n ** 3 - n)
    h = h / c if c > 0 else float("nan")
    return h, _gammq((len(groups) - 1) / 2.0, h / 2.0)


def _mw_exact_p(r, n1, stat):
    """exact two-sided p of the rank sum of n1 of the (mid)ranks r (DP over doubled integer ranks)"""
    rr = [int(round(2 * x)) for x in r]
    tot = sum(rr)
    # dp[k][s] = number of subsets of size k with doubled rank sum s
    dp = [collections.Counter() for _ in range(n1 + 1)]
    dp[0][0] = 1
    for x in rr:
        for k in range(min(n1, len(rr)) - 1, -1, -1):
            for sm, c in dp[k].items():
                dp[k + 1][sm + x] += c
    dist = dp[n1]
    allc = sum(dist.values())
    mean = n1 * tot / len(rr)
    obs = abs(2 * stat - mean)
    return min(1.0, sum(c for sm, c in dist.items() if abs(sm - mean) >= obs - 1e-9) / allc)


def mannwhitney(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    n1, n2 = len(a), len(b)
    if n1 < 2 or n2 < 2:
        return float("nan"), float("nan")
    r = rankdata(np.concatenate([a, b]))
    u = r[:n1].sum() - n1 * (n1 + 1) / 2.0
    n = n1 + n2
    if n <= 20:                                   # exact p for small samples (normal approx. too liberal)
        if np.all(r == r[0]):
            return float("nan"), float("nan")
        return u / (n1 * n2), _mw_exact_p(r, n1, r[:n1].sum())
    var = n1 * n2 / 12.0 * ((n + 1) - _tie_term(r) / (n * (n - 1)))
    if var <= 0:
        return float("nan"), float("nan")
    z = (u - n1 * n2 / 2.0) / math.sqrt(var)
    return u / (n1 * n2), math.erfc(abs(z) / math.sqrt(2))   # effect: P(a > b) (AUC)


def van_elteren(a, ba, b, bb):
    """stratified Wilcoxon (van Elteren) test of a vs b within blocks: values a with block labels ba, values b with
    block labels bb; only blocks with both a and b count. Returns (block-weighted P(a > b), two-sided p, blocks).
    With one informative block it is the Mann-Whitney test (exact for small n)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    ba, bb = np.asarray(ba, object), np.asarray(bb, object)
    blocks = sorted(set(ba) & set(bb))
    if not blocks:
        return float("nan"), float("nan"), 0
    if len(blocks) == 1:
        auc, p = mannwhitney(a[ba == blocks[0]], b[bb == blocks[0]])
        return auc, p, 1
    num = var = 0.0
    aucs, ws = [], []
    for bl in blocks:
        x, y = a[ba == bl], b[bb == bl]
        n1, n2 = len(x), len(y)
        n = n1 + n2
        r = rankdata(np.concatenate([x, y]))
        w = 1.0 / (n + 1)
        v = n1 * n2 / 12.0 * ((n + 1) - _tie_term(r) / (n * (n - 1)))
        num += w * (r[:n1].sum() - n1 * (n + 1) / 2.0)
        var += w * w * v
        aucs.append((r[:n1].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n2))
        ws.append(n1 * n2 / (n + 1.0))
    if var <= 0:
        return float("nan"), float("nan"), len(blocks)
    z = num / math.sqrt(var)
    return float(np.average(aucs, weights=ws)), math.erfc(abs(z) / math.sqrt(2)), len(blocks)


def wilcoxon_paired(a, b):
    d = np.asarray(a, float) - np.asarray(b, float)
    d = d[d != 0]
    n = len(d)
    if n < 5:
        return float("nan"), float("nan"), n
    r = rankdata(np.abs(d))
    w = r[d > 0].sum()
    var = n * (n + 1) * (2 * n + 1) / 24.0 - _tie_term(r) / 48.0
    z = (w - n * (n + 1) / 4.0) / math.sqrt(var)
    return float(np.median(d)), math.erfc(abs(z) / math.sqrt(2)), n


def spearman(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = ~(np.isnan(x) | np.isnan(y))
    x, y = x[ok], y[ok]
    n = len(x)
    if n < 4 or np.all(x == x[0]) or np.all(y == y[0]):
        return float("nan"), float("nan"), n
    rho = float(np.corrcoef(rankdata(x), rankdata(y))[0, 1])
    if abs(rho) >= 1:
        return rho, 0.0, n
    t = rho * math.sqrt((n - 2) / (1 - rho ** 2))         # t approximation, normal tail for df >= 30
    df = n - 2
    p = _t_sf2(abs(t), df)
    return rho, p, n


def _betacf(a, b, x):
    qab, qap, qam = a + b, a + 1, a - 1
    c, d = 1.0, 1 - qab * x / qap
    d = 1.0 / (d if abs(d) > 1e-300 else 1e-300)
    h = d
    for m in range(1, 300):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1 + aa * d
        d = 1.0 / (d if abs(d) > 1e-300 else 1e-300)
        c = 1 + aa / c
        c = c if abs(c) > 1e-300 else 1e-300
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1 + aa * d
        d = 1.0 / (d if abs(d) > 1e-300 else 1e-300)
        c = 1 + aa / c
        c = c if abs(c) > 1e-300 else 1e-300
        de = d * c
        h *= de
        if abs(de - 1) < 1e-12:
            break
    return h


def _betai(a, b, x):
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    bt = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log(1 - x))
    if x < (a + 1) / (a + b + 2):
        return bt * _betacf(a, b, x) / a
    return 1 - bt * _betacf(b, a, 1 - x) / b


def _t_sf2(t, df):
    """two-sided p of Student's t"""
    return _betai(df / 2.0, 0.5, df / (df + t * t))


def bh(pvals):
    p = np.asarray(pvals, float)
    q = np.full(len(p), np.nan)
    ok = ~np.isnan(p)
    if ok.sum() == 0:
        return q
    pv = p[ok]
    order = np.argsort(pv)
    m = len(pv)
    adj = pv[order] * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(adj, 1)
    q[ok] = out
    return q


def bray_curtis(X):
    """X: samples x features (relative abundances) -> distance matrix"""
    n = X.shape[0]
    D = np.zeros((n, n))
    for i in range(n):
        num = np.abs(X[i] - X).sum(axis=1)
        den = (X[i] + X).sum(axis=1)
        D[i] = np.where(den > 0, num / np.where(den > 0, den, 1), 0)
    return (D + D.T) / 2


def pcoa(D, k=2):
    n = D.shape[0]
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ (D ** 2) @ J
    w, v = np.linalg.eigh(B)
    idx = np.argsort(w)[::-1]
    w, v = w[idx], v[:, idx]
    pos = w[w > 0].sum() or 1
    coords = v[:, :k] * np.sqrt(np.maximum(w[:k], 0))
    return coords, w[:k] / pos


# ----------------------------------------------------------------------------- adjusted linear models
def clr(C, pseudo=0.5):
    """centred log-ratio of a count matrix (samples x features) with a pseudocount"""
    L = np.log(C + pseudo)
    return L - L.mean(axis=1, keepdims=True)


def design_matrix(rows, terms, kinds, levels_first=None):
    """rows: list of dicts; terms: variable names; kinds: name -> numeric / categorical.
    Returns (X, column names, column -> term, info) with an intercept; categorical terms as treatment
    dummies (reference = levels_first[term] or the most frequent level); numeric terms scaled to SD units.
    Linearly dependent columns (e.g. a covariate fully confounded with the group) are dropped and reported."""
    cols, names, owner, info = [np.ones(len(rows))], ["(intercept)"], ["(intercept)"], []
    for t in terms:
        vals = [r.get(t, "NA") for r in rows]
        if kinds.get(t) == "numeric":
            x = np.array([float(v) for v in vals])
            sd = x.std(ddof=1) or 1.0
            cols.append((x - x.mean()) / sd)
            names.append("%s (per SD %s)" % (t, fmt(sd, 3)))
            owner.append(t)
        else:
            lv = sorted(set(vals), key=lambda v: (-vals.count(v), v))
            ref = (levels_first or {}).get(t)
            if ref not in lv:                       # natural reference for yes/no or normal/abnormal variables
                ref = next((v for v in lv if str(v).strip().lower() in ("no", "none", "normal", "no problems", "0",
                                                                         "false", "never")), lv[0])
            for l in [x for x in lv if x != ref]:
                cols.append(np.array([1.0 if v == l else 0.0 for v in vals]))
                names.append("%s: %s vs %s" % (t, l, ref))
                owner.append(t)
    X = np.column_stack(cols)
    keep = [0]
    for j in range(1, X.shape[1]):
        if np.linalg.matrix_rank(X[:, keep + [j]]) > len(keep):
            keep.append(j)
        else:
            info.append("dropped (collinear): %s" % names[j])
    return X[:, keep], [names[j] for j in keep], [owner[j] for j in keep], info


MODEL_MIN_N = 20          # minimum libraries per model fit (and >= 3 x terms)
MODEL_MIN_N_ASIDX = 40    # antisense index: only libraries with >= 10 reads on the term count -> require more
MODEL_MIN_LEVEL = 3       # minimum libraries in and outside every 0/1 design column (e.g. a group or site level)
MODEL_MAX_LEVERAGE = 0.95 # fits with a library of leverage above this are skipped (HC3 SE undefined / collapses)


def design_ok(X, off=0):
    """enough libraries for a stable fit: n, libraries per categorical level, no library that alone carries a term"""
    n, p = X.shape
    if n < max(MODEL_MIN_N, 3 * (p - off)):
        return False
    for c in range(1 if off == 0 else off, p):
        col = X[:, c]
        u = np.unique(col[~np.isnan(col)])
        if len(u) <= 2:                                  # dummy (0/1) or a binary numeric term
            k = int(np.sum(col == u.max())) if len(u) else 0
            if min(k, n - k) < MODEL_MIN_LEVEL:
                return False
    try:
        Xc = X[:, off:]
        h = np.einsum("ij,ji->i", Xc, np.linalg.pinv(Xc.T @ Xc) @ Xc.T)
        if h.max() > MODEL_MAX_LEVERAGE:
            return False
    except np.linalg.LinAlgError:
        return False
    return True


def ols_hc3(X, Y):
    """OLS of every column of Y (n x m) on X (n x p): estimates, HC3 robust SE, t, two-sided p (t, df n-p)"""
    n, p = X.shape
    XtXi = np.linalg.pinv(X.T @ X)
    A = XtXi @ X.T                                   # p x n
    B = A @ Y                                        # p x m
    E = Y - X @ B
    h = np.einsum("ij,ji->i", X, A)                  # leverages
    w = 1.0 / np.maximum(1 - h, 1e-8) ** 2
    V = (A ** 2) @ (E ** 2 * w[:, None])             # p x m: HC3 variances
    SE = np.sqrt(np.maximum(V, 0))
    with np.errstate(divide="ignore", invalid="ignore"):
        T = np.where(SE > 0, B / SE, np.nan)
    df = n - p
    P = np.vectorize(lambda t: _t_sf2(abs(t), df) if np.isfinite(t) else np.nan)(T) if df > 0 else T * np.nan
    return B, SE, T, P, df


def ols_cluster(X, Y, clusters):
    """OLS with cluster-robust (CR1) standard errors: residuals of the same cluster (e.g. family) may be
    correlated; t-test with G - 1 degrees of freedom"""
    n, p = X.shape
    cl = np.asarray(clusters, object)
    levels = list(set(cl))
    G = len(levels)
    XtXi = np.linalg.pinv(X.T @ X)
    B = XtXi @ X.T @ Y
    E = Y - X @ B
    meat = np.zeros((p, Y.shape[1]))
    A = XtXi @ X.T                                      # p x n
    for g in levels:
        idx = np.where(cl == g)[0]
        S = A[:, idx] @ E[idx]                          # p x m: (X'X)^-1 X_g' e_g
        meat += S ** 2
    c = G / max(G - 1, 1) * (n - 1) / max(n - p, 1)
    SE = np.sqrt(np.maximum(meat * c, 0))
    with np.errstate(divide="ignore", invalid="ignore"):
        T = np.where(SE > 0, B / SE, np.nan)
    df = G - 1
    P = np.vectorize(lambda t: _t_sf2(abs(t), df) if np.isfinite(t) else np.nan)(T) if df > 0 else T * np.nan
    return B, SE, T, P, df


def demean_within(M, blocks):
    """subtract the block mean from every column (fixed effects by the within transformation)"""
    M = np.array(M, float)
    bl = np.asarray(blocks, object)
    for g in set(bl):
        idx = np.where(bl == g)[0]
        M[idx] -= M[idx].mean(axis=0)
    return M


def _pava(y):
    """isotonic (non-decreasing) regression, pool-adjacent-violators"""
    y = list(map(float, y))
    vals, wts, lens = [], [], []
    for v in y:
        vals.append(v); wts.append(1.0); lens.append(1)
        while len(vals) > 1 and vals[-2] > vals[-1]:
            w = wts[-2] + wts[-1]
            v2 = (vals[-2] * wts[-2] + vals[-1] * wts[-1]) / w
            l2 = lens[-2] + lens[-1]
            vals[-2:], wts[-2:], lens[-2:] = [v2], [w], [l2]
    return np.repeat(vals, lens)


def _stress1(D, X):
    """Kruskal stress-1 of configuration X for dissimilarities D (monotone fit, secondary ties)"""
    iu = np.triu_indices(D.shape[0], 1)
    d = D[iu]
    dv = np.sqrt(((X[:, None, :] - X[None, :, :]) ** 2).sum(axis=2))[iu]
    order = np.lexsort((dv, d))
    dh = np.empty_like(dv)
    dh[order] = _pava(dv[order])
    return float(np.sqrt(((dv - dh) ** 2).sum() / max((dv ** 2).sum(), 1e-12)))


def nmds(D, k=2, starts=10, max_iter=300, eps=1e-7, seed=1):
    """non-metric MDS (Kruskal) by SMACOF with monotone regression on the dissimilarities D (n x n).
    Ties in D: secondary approach (tied dissimilarities keep the order of the fitted distances).
    Returns (coordinates n x k, Kruskal stress-1); best of a PCoA start and `starts` random starts."""
    n = D.shape[0]
    if n < k + 2:
        return None, float("nan")
    iu = np.triu_indices(n, 1)
    d = D[iu]
    rng = np.random.default_rng(seed)
    inits = [pcoa(D, k)[0]] + [rng.normal(size=(n, k)) for _ in range(starts)]
    best = (None, float("inf"))
    for X in inits:
        X = X - X.mean(axis=0)
        sx = np.sqrt((X ** 2).sum() / n) or 1.0
        X = X / sx
        prev = float("inf")
        for _ in range(max_iter):
            diff = X[:, None, :] - X[None, :, :]
            dist = np.sqrt((diff ** 2).sum(axis=2))
            dv = dist[iu]
            order = np.lexsort((dv, d))                 # by dissimilarity, ties by current distance
            dh = np.empty_like(dv)
            dh[order] = _pava(dv[order])
            dh *= np.sqrt(len(dh) / max((dh ** 2).sum(), 1e-12))   # normalise disparities
            stress = np.sqrt(((dv - dh) ** 2).sum() / max((dv ** 2).sum(), 1e-12))
            if prev - stress < eps:
                break
            prev = stress
            Dh = np.zeros((n, n))
            Dh[iu] = dh
            Dh = Dh + Dh.T
            with np.errstate(divide="ignore", invalid="ignore"):
                Bm = -np.where(dist > 0, Dh / dist, 0.0)
            Bm[np.diag_indices(n)] = 0
            Bm[np.diag_indices(n)] = -Bm.sum(axis=1)
            X = Bm @ X / n                               # Guttman transform
        stress = _stress1(D, X)                          # stress of the configuration actually returned
        if stress < best[1]:
            best = (X.copy(), stress)
    X, st = best
    # orient: first axis = direction of largest variance
    w, v = np.linalg.eigh(np.cov(X.T))
    X = X @ v[:, np.argsort(w)[::-1]]
    return X, float(st)


def isolated_samples(D, max_frac=0.05):
    """indices of samples that are further from every other sample than any two remaining samples are from each
    other (min distance to the rest > max distance within the rest). Such a sample makes NMDS degenerate: the rest
    collapse to one point and stress goes to 0. Removed one at a time, at most max_frac of the samples."""
    keep = list(range(D.shape[0]))
    out = []
    while len(keep) > 6 and len(out) < max(1, int(max_frac * D.shape[0])):
        sub = D[np.ix_(keep, keep)]
        off = sub + np.diag([np.inf] * len(keep))
        nearest = off.min(axis=1)
        i = int(np.argmax(nearest))
        rest = [j for j in range(len(keep)) if j != i]
        if nearest[i] <= sub[np.ix_(rest, rest)].max():
            break
        out.append(keep.pop(i))
    return out


def permanova(D, labels, perms=999, seed=1, strata=None):
    labels = np.asarray(labels)
    lv = sorted(set(labels))
    n, a = len(labels), len(lv)
    if a < 2 or n <= a or min(collections.Counter(labels).values()) < 2:
        return float("nan"), float("nan"), float("nan")
    D2 = D ** 2
    sst = D2[np.triu_indices(n, 1)].sum() / n

    def ssw(lab):
        s = 0.0
        for g in lv:
            idx = np.where(lab == g)[0]
            if len(idx) > 1:
                sub = D2[np.ix_(idx, idx)]
                s += sub[np.triu_indices(len(idx), 1)].sum() / len(idx)
        return s

    w = ssw(labels)
    f = ((sst - w) / (a - 1)) / (w / (n - a)) if w > 0 else float("nan")
    rng = np.random.default_rng(seed)
    groups_ = None
    if strata is not None:                       # permutations only within strata (blocks)
        strata = np.asarray(strata, object)
        groups_ = [np.where(strata == s_)[0] for s_ in set(strata)]

    def shuffle():
        if groups_ is None:
            return rng.permutation(labels)
        out = labels.copy()
        for idx in groups_:
            out[idx] = labels[rng.permutation(idx)]
        return out
    ge = 0
    for _ in range(perms):
        wp = ssw(shuffle())
        fp = ((sst - wp) / (a - 1)) / (wp / (n - a)) if wp > 0 else 0
        ge += fp >= f
    return f, (sst - w) / sst if sst else float("nan"), (ge + 1) / (perms + 1)


def permanova_partial(D, labels, strata, perms=999, seed=1):
    """PERMANOVA of `labels` after `strata` (e.g. block after group): partial R² = SS(labels | strata) / SS total;
    labels permuted within strata. Returns (F, partial R², p)."""
    labels, strata = np.asarray(labels, object), np.asarray(strata, object)
    n = len(labels)
    D2 = D ** 2
    sst = D2[np.triu_indices(n, 1)].sum() / n

    def ssw(lab):
        s = 0.0
        for g in set(lab):
            idx = np.where(lab == g)[0]
            if len(idx) > 1:
                sub = D2[np.ix_(idx, idx)]
                s += sub[np.triu_indices(len(idx), 1)].sum() / len(idx)
        return s
    cell = lambda lab: np.array(["%s|%s" % (a_, b_) for a_, b_ in zip(strata, lab)], object)
    w_s = ssw(strata)
    a_s, a_c = len(set(strata)), len(set(cell(labels)))
    df1, df2 = a_c - a_s, n - a_c
    if df1 < 1 or df2 < 1 or not sst:
        return float("nan"), float("nan"), float("nan")
    w_c = ssw(cell(labels))
    stat = lambda wc: ((w_s - wc) / df1) / (wc / df2) if wc > 0 else float("inf")
    f = stat(w_c)
    rng = np.random.default_rng(seed)
    groups_ = [np.where(strata == s_)[0] for s_ in set(strata)]
    ge = 0
    for _ in range(perms):
        out = labels.copy()
        for idx in groups_:
            out[idx] = labels[rng.permutation(idx)]
        ge += stat(ssw(cell(out))) >= f
    return f, (w_s - w_c) / sst, (ge + 1) / (perms + 1)


# ----------------------------------------------------------------------------- metadata
def is_missing(v):
    return v is None or str(v).strip().lower() in MISSING


def read_metadata(path):
    """TSV / CSV / XLSX (first sheet; needs openpyxl) -> list of dicts (all values as strings)"""
    ext = os.path.splitext(path)[1].lower()
    if ext in (".xlsx", ".xlsm"):
        try:
            import openpyxl
        except ImportError:
            sys.exit("[cohort] reading %s needs openpyxl — export the sheet as TSV/CSV or install openpyxl" % path)
        ws = openpyxl.load_workbook(path, data_only=True, read_only=True).worksheets[0]
        rows = [["" if c is None else str(c).strip() for c in r] for r in ws.iter_rows(values_only=True)]
        rows = [r for r in rows if any(r)]
    else:
        with open(path, newline="", encoding="utf-8-sig") as fh:
            text = fh.read()
        delim = "\t" if text.count("\t") >= text.count(",") else ","
        rows = [[c.strip() for c in r] for r in csv.reader(text.splitlines(), delimiter=delim) if any(r)]
        rows = [r for r in rows if not r[0].startswith("#")]
    hdr = rows[0]
    return hdr, [dict(zip(hdr, r + [""] * (len(hdr) - len(r)))) for r in rows[1:]]


def column_kind(values):
    vals = [v for v in values if not is_missing(v)]
    if not vals:
        return "empty"
    try:
        [float(v) for v in vals]
        return "numeric" if len(set(vals)) > 2 else "categorical"
    except ValueError:
        return "categorical"


def match_ids(meta_rows, id_col, samples):
    """metadata ID -> pipeline sample: exact, case-insensitive, or before the first '.' (BXXX1.t0F)"""
    norm = {}
    for s in samples:
        norm.setdefault(s.upper(), s)
        norm.setdefault(s.split(".")[0].upper(), s)
    out, unmatched = {}, []
    for r in meta_rows:
        mid = r.get(id_col, "").strip()
        s = norm.get(mid.upper()) or norm.get(mid.split(".")[0].upper())
        if s:
            out[s] = r
        elif mid:
            unmatched.append(mid)
    return out, unmatched


# ----------------------------------------------------------------------------- collection
def find_samples(outdir, samplesheet=""):
    """(sample, type) pairs with an Analysis/<type> directory"""
    pairs = []
    wanted = None
    if samplesheet:
        wanted = set()
        with open(samplesheet) as fh:
            for ln in fh:
                f = ln.split()
                if f and not f[0].startswith("#") and f[0] != "sample_id" and len(f) > 1:
                    wanted.add((f[0], f[1]))
    for d in sorted(glob.glob(os.path.join(outdir, "*", "Analysis"))):
        s = os.path.basename(os.path.dirname(d))
        for td in sorted(glob.glob(os.path.join(d, "*"))):
            t = os.path.basename(td)
            if t in ("annotation", "comparison") or not os.path.isdir(td):
                continue
            if wanted is None or (s, t) in wanted:
                pairs.append((s, t))
    return pairs


def frac(a, b):
    return a / b if b else float("nan")


def sample_metrics(outdir, s, t):
    """one dict of key numbers per sample/type + the count vectors"""
    d = R.TypeData(os.path.join(outdir, s), s, t)
    d._rfam_family = R.read_tsv(os.path.join(outdir, s, "Analysis", t, "nonhuman", "rfam", "%s.rfam_family.counts.tsv" % s))
    m = Q.library_metrics(d, outdir, s, t)
    # count vectors
    vec = {
        "mirna": {r["feature"]: R.num(r["count"]) for r in d.mirna},
        "genus": {r["name"]: R.num(r["new_est_reads"]) for r in d.bracken["G"]},
        "family": {r["name"]: R.num(r["new_est_reads"]) for r in d.bracken["F"]},
        "species": {r["name"]: R.num(r["new_est_reads"]) for r in d.bracken["S"]},
        "rfam_type": dict(d.rfam_type),
        "trna_aa": dict(d.trna_aa),
    }
    # functional counts without the RiboDetector rRNA reads when available (assembly module with --exclude-reads):
    # reads sense to unannotated rRNA otherwise appear as antisense to spurious ORFs
    nr = bool(d.asm_functions) and "antisense_reads_non_rRNA" in d.asm_functions[0]
    m["antisense_rRNA_filtered"] = 1 if nr else 0
    sk, ak = ("sense_reads_non_rRNA", "antisense_reads_non_rRNA") if nr else ("sense_reads", "antisense_reads")
    for db, (keep, _) in FUNC_DBS.items():
        vec["fn_" + db], vec["fnas_" + db] = {}, {}
        for r in d.asm_functions:
            if r.get("db") == db and keep(r.get("term", "")):
                vec["fn_" + db][r["term"]] = R.num(r.get(sk))
                vec["fnas_" + db][r["term"]] = R.num(r.get(ak))
    vec["rna_class"] = rna_classes(d)
    # Rfam families other than rRNA / tRNA (their shares would be dominated by rRNA otherwise)
    vec["rfam_family"] = {(r.get("rfam_id") or r.get("rfam_family")): R.num(r.get("count")) for r in d._rfam_family
                          if not (r.get("type") or "").startswith(("rRNA", "tRNA"))}
    tot = sum(vec["rna_class"].values())
    for c, v in vec["rna_class"].items():
        m["class_%s" % c.replace(" ", "_").replace("-", "_")] = v / tot if tot else float("nan")
    # antisense on genes; strand error from reads antisense to rRNA / tRNA genes (should be ~0)
    cat = d.asm_category
    catf = R.read_tsv(os.path.join(outdir, s, "Analysis", t, "assembly", "%s.category_counts.tsv" % s))
    if catf and "reads_non_rRNA" in catf[0]:
        m["antisense_CDS_fraction_incl_rRNA_reads"] = (lambda a, c: a / (a + c) if a + c else float("nan"))(
            R.num(cat.get("antisense_CDS")), R.num(cat.get("CDS")))
        cat = {r["category"]: R.num(r["reads_non_rRNA"]) for r in catf}
    cds, as_cds = R.num(cat.get("CDS")), R.num(cat.get("antisense_CDS"))
    rt = R.num(cat.get("rRNA")) + R.num(cat.get("tRNA"))
    as_rt = R.num(cat.get("antisense_rRNA")) + R.num(cat.get("antisense_tRNA"))
    m["antisense_CDS_fraction"] = as_cds / (cds + as_cds) if cds + as_cds > 0 else float("nan")
    m["antisense_strand_error"] = as_rt / (rt + as_rt) if rt + as_rt > 0 else float("nan")
    if isnum(m["antisense_CDS_fraction"]) and isnum(m["antisense_strand_error"]):
        m["antisense_CDS_enrichment_log2"] = math.log2((m["antisense_CDS_fraction"] + 1e-6) /
                                                       (m["antisense_strand_error"] + 1e-6))
    m["cds_reads_with_KO"] = sum(vec["fn_kegg_ko"].values()) if vec["fn_kegg_ko"] else float("nan")
    # marker-group shares of human miRNA reads
    tot = sum(vec["mirna"].values())
    for g, _c, mirs, _r in R.MIRNA_ORIGIN:
        ms = set(mirs)
        v = sum(c for k, c in vec["mirna"].items() if k.replace("hsa-", "") in ms)
        m["miRNA_marker_%s" % "_".join(g.replace("/", " ").split())] = frac(v, tot) if tot >= 20 else float("nan")
    return m, vec


def compare_metrics(outdir, s):
    """sample-level Spearman of sRNA vs metaT / metaG (comparison step)"""
    out = {}
    p = os.path.join(outdir, s, "Analysis", "comparison", "%s.compare.summary.tsv" % s)
    for r in R.read_tsv(p):
        if r.get("metric") == "spearman_cpm" and r.get("feature_type") == "all" and \
                r.get("level") in ("feature", "contig", "bin"):
            pair = r["pair"]
            if "_vs_" not in pair:
                continue
            a, b = pair.split("_vs_", 1)
            out[(a, "rho_vs_%s_%s" % (b, r["level"]))] = R.num(r["value"], float("nan"))
            out[(b, "rho_vs_%s_%s" % (a, r["level"]))] = R.num(r["value"], float("nan"))
    return out


# ----------------------------------------------------------------------------- analysis
def fmt(v, digits=4):
    if isinstance(v, float):
        if math.isnan(v):
            return "NA"
        return ("%.*g" % (digits, v))
    return str(v)


def write_tsv(path, header, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as out:
        out.write("\t".join(header) + "\n")
        for r in rows:
            out.write("\t".join(fmt(x) for x in r) + "\n")


class Cohort:
    def __init__(self, a):
        self.a = a
        self.pairs = find_samples(a.outdir, a.samples)
        if not self.pairs:
            sys.exit("[cohort] no samples with Analysis/<type> under %s" % a.outdir)
        self.samples = sorted({s for s, _ in self.pairs})
        self.types = sorted({t for _, t in self.pairs}, key=R._src_order)
        # metadata
        self.meta, self.unmatched, self.meta_cols = {}, [], []
        self.group_col = a.group_col
        if a.metadata:
            hdr, rows = read_metadata(a.metadata)
            id_col = a.id_col or hdr[0]
            if id_col not in hdr:
                sys.exit("[cohort] id column '%s' not in metadata (columns: %s)" % (id_col, ", ".join(hdr)))
            if self.group_col and self.group_col not in hdr:
                sys.exit("[cohort] group column '%s' not in metadata (columns: %s)" % (self.group_col, ", ".join(hdr)))
            self.meta, self.unmatched = match_ids(rows, id_col, self.samples)
            covs = [c for c in (a.covariates.split(",") if a.covariates else hdr) if c and c not in (id_col,)]
            missing = [c for c in covs if c not in hdr]
            if missing:
                sys.exit("[cohort] covariate column(s) not in metadata: %s" % ", ".join(missing))
            self.meta_cols = [c for c in covs if c != self.group_col]
            self.id_col = id_col
        self.kinds = {c: column_kind([self.meta[s].get(c) for s in self.meta]) for c in self.meta_cols}
        order = [g.strip() for g in a.group_order.split(",")] if a.group_order else []
        found = sorted({self.group(s) for s in self.samples} - {"NA"})
        self.groups = [g for g in order if g in found] + [g for g in found if g not in order]
        self.gcolor = {g: R.TYPE_COLORS[i % len(R.TYPE_COLORS)] for i, g in enumerate(self.groups)}
        self.gcolor["NA"] = "#999999"
        # data
        self.m, self.vec, self.qc = {}, {}, {}
        rules, excl = Q.load_rules(a.qc_rules), Q.load_exclude(a.qc_exclude)
        drop = {"FAIL": ("FAIL",), "WARN": ("WARN", "FAIL"), "none": ()}[a.qc_drop]
        cmp_ = {s: compare_metrics(a.outdir, s) for s in self.samples}
        for i, (s, t) in enumerate(self.pairs):
            try:
                m, v = sample_metrics(a.outdir, s, t)
            except Exception as e:      # one broken sample must not stop the cohort
                sys.stderr.write("[cohort] %s/%s skipped (%s: %s)\n" % (s, t, type(e).__name__, e))
                continue
            for (tt, k), val in cmp_[s].items():
                if tt == t:
                    m[k] = val
            self.m[(s, t)], self.vec[(s, t)] = m, v
            self.qc[(s, t)] = Q.evaluate(m, s, t, rules, excl)
            if (i + 1) % 50 == 0:
                sys.stderr.write("[cohort] %d / %d sample-types read\n" % (i + 1, len(self.pairs)))
        self.flag_outliers()
        self.excluded = {k for k, (lv, _) in self.qc.items() if lv in drop}
        self.metrics = []
        for m in self.m.values():
            for k, v in m.items():
                if k not in self.metrics and not k.startswith(("_", "fastqc_")) and k not in NOT_TESTED \
                        and isinstance(v, (int, float)):
                    self.metrics.append(k)

    def flag_outliers(self):
        """robust outliers within a fraction (|x - median| / (1.4826 MAD) > --outlier-z) -> WARN 'outlier'"""
        z0 = self.a.outlier_z
        if not z0 or z0 <= 0:
            return
        for t in self.types:
            ks = [k for k in self.m if k[1] == t and not any(h[0] == "incomplete" for h in self.qc[k][1])]
            for mk, min_abs in OUTLIER_METRICS.items():
                v = [(k, self.m[k].get(mk, float("nan"))) for k in ks]
                v = [(k, x) for k, x in v if isinstance(x, (int, float)) and not math.isnan(x)]
                if len(v) < 10:
                    continue
                med = float(np.median([x for _, x in v]))
                mad = 1.4826 * float(np.median([abs(x - med) for _, x in v]))
                if mad <= 0:
                    continue
                for k, x in v:
                    z = (x - med) / mad
                    if abs(z) > z0 and abs(x - med) >= min_abs:
                        lv, hits = self.qc[k]
                        hits.append(("outlier", "WARN", x, "robust z", "%.1f" % z,
                                     "%s unusual for %s (median %s)" % (mk, t, fmt(med, 3))))
                        if lv == "PASS":
                            lv = "WARN"
                        self.qc[k] = (lv, hits)

    def group(self, s):
        if not self.group_col:
            return "all"
        v = self.meta.get(s, {}).get(self.group_col, "")
        return "NA" if is_missing(v) else v.strip()

    def mirna_types(self, min_libs=6):
        """fractions with enough human miRNA data for miRNA analyses (>= min_libs libraries with
        >= --min-mirna-reads reads); others (typically sRNA-Ext) are left out of miRNA plots and tests"""
        if not hasattr(self, "_mtypes"):
            self._mtypes = [t for t in self.types
                            if sum(1 for k in self.keys(t) if sum(self.vec[k]["mirna"].values()) >= self.a.min_mirna_reads)
                            >= min_libs]
        return self._mtypes

    def keys(self, t, all_=False):
        """libraries of a type that enter the statistics (QC-excluded ones only with all_=True)"""
        return [(s, tt) for (s, tt) in self.m if tt == t and (all_ or (s, tt) not in self.excluded)]

    # ---------------- tables
    def write_tables(self, d):
        hdr = ["sample", "type", "group", "qc", "in_statistics"] + self.meta_cols + self.metrics
        rows = []
        for (s, t), m in sorted(self.m.items()):
            md = self.meta.get(s, {})
            rows.append([s, t, self.group(s), self.qc[(s, t)][0], "no" if (s, t) in self.excluded else "yes"] +
                        [md.get(c, "") or "NA" for c in self.meta_cols] +
                        [m.get(k, float("nan")) for k in self.metrics])
        write_tsv(os.path.join(d, "tables", "cohort.qc_flags.tsv"),
                  ["sample", "type", "group", "qc", "in_statistics", "rules", "details"],
                  [[s, t, self.group(s), lv, "no" if (s, t) in self.excluded else "yes",
                    ",".join(h[0] for h in hits) or "-",
                    "; ".join("%s %s (%s %s %s)" % (h[1], h[5], R.short(fmt(h[2]), 12), h[3], h[4]) if h[3] else
                              "%s %s" % (h[1], h[5]) for h in hits) or "-"]
                   for (s, t), (lv, hits) in sorted(self.qc.items())])
        write_tsv(os.path.join(d, "tables", "cohort.sample_stats.tsv"), hdr, rows)
        for t in self.types:
            ks = sorted(self.keys(t))            # QC-excluded libraries are left out (see cohort.qc_flags.tsv)
            for mx in MATRICES:
                feats = sorted({f for k in ks for f in self.vec[k][mx]})
                if not feats:
                    continue
                write_tsv(os.path.join(d, "tables", "cohort.%s.%s.counts.tsv" % (t, mx)), ["feature"] + [s for s, _ in ks],
                          [[f] + [self.vec[k][mx].get(f, 0.0) for k in ks] for f in feats])
        if self.meta:
            cols = [self.id_col, "sample"] + ([self.group_col] if self.group_col else []) + self.meta_cols
            qcc = ["qc_" + t for t in self.types]
            write_tsv(os.path.join(d, "tables", "cohort.metadata.tsv"), ["sample"] + cols[2:] + qcc,
                      [[s] + [self.meta[s].get(c, "") or "NA" for c in cols[2:]] +
                       [("excluded_" if (s, t) in self.excluded else "") + self.qc[(s, t)][0] if (s, t) in self.qc
                        else "missing" for t in self.types] for s in sorted(self.meta)])
        with open(os.path.join(d, "tables", "cohort.unmatched.txt"), "w") as out:
            out.write("# metadata IDs without pipeline output\n" + "\n".join(self.unmatched) + "\n")
            out.write("# pipeline samples without metadata\n" + "\n".join(s for s in self.samples if s not in self.meta) + "\n")

    # ---------------- statistics
    def metric_pq(self, t, metric):
        if not hasattr(self, "_mpq"):
            self._mpq = {}
        if t not in self._mpq:
            rows, _ = self.metric_tests(t)
            self._mpq[t] = {r[0]: (r[3], r[4]) for r in rows}
        return self._mpq[t].get(metric, (float("nan"), float("nan")))

    def metric_tests(self, t):
        ks = self.keys(t)
        rows = []
        for k in self.metrics:
            by = collections.OrderedDict((g, []) for g in self.groups)
            for key in ks:
                g = self.group(key[0])
                v = self.m[key].get(k, float("nan"))
                if g in by and not math.isnan(v):
                    by[g].append(v)
            h, p = kruskal(list(by.values()))
            med = [float(np.median(v)) if v else float("nan") for v in by.values()]
            rows.append([k, sum(len(v) for v in by.values()), h, p] + med)
        q = bh([r[3] for r in rows])
        for r, qq in zip(rows, q):
            r.insert(4, qq)
        # pairwise for metrics
        pw = []
        for k in self.metrics:
            vals = {g: [self.m[key].get(k, float("nan")) for key in ks if self.group(key[0]) == g] for g in self.groups}
            vals = {g: [x for x in v if not math.isnan(x)] for g, v in vals.items()}
            for i, g1 in enumerate(self.groups):
                for g2 in self.groups[i + 1:]:
                    auc, p = mannwhitney(vals[g1], vals[g2])
                    pw.append([k, g1, g2, len(vals[g1]), len(vals[g2]), auc, p])
        qp = bh([r[6] for r in pw])
        for r, qq in zip(pw, qp):
            r.append(qq)
        return rows, pw

    def rel_matrix(self, t, mx, min_prev=0.1, min_reads=10):
        floor = self.a.min_mirna_reads if mx == "mirna" else 1   # miRNA profiles need enough reads
        ks = [k for k in sorted(self.keys(t)) if sum(self.vec[k][mx].values()) >= max(floor, 1e-9)]
        feats = sorted({f for k in ks for f in self.vec[k][mx]})
        if not ks or not feats:
            return ks, [], np.zeros((0, 0)), np.zeros((0, 0))
        C = np.array([[self.vec[k][mx].get(f, 0.0) for f in feats] for k in ks])
        Rl = C / C.sum(axis=1, keepdims=True)
        keep = ((C >= min_reads).mean(axis=0) >= min_prev)
        return ks, [f for f, kk in zip(feats, keep) if kk], C[:, keep], Rl[:, keep]

    def feature_tests(self, t, mx):
        ks, feats, C, Rl = self.rel_matrix(t, mx)
        rows = []
        if not feats or len(self.groups) < 2:
            return rows
        labels = [self.group(k[0]) for k in ks]
        for j, f in enumerate(feats):
            by = [[Rl[i, j] for i in range(len(ks)) if labels[i] == g] for g in self.groups]
            h, p = kruskal(by)
            prev = float((C[:, j] > 0).mean())
            rows.append([f, prev, h, p] + [float(np.mean(b)) if b else float("nan") for b in by])
        q = bh([r[3] for r in rows])
        for r, qq in zip(rows, q):
            r.insert(4, qq)
        rows.sort(key=lambda r: (math.isnan(r[3]), r[3]))
        return rows

    def beta(self, t, mx):
        ks, feats, C, Rl = self.rel_matrix(t, mx, min_prev=0.0, min_reads=1)
        if len(ks) < 4:
            return None
        D = bray_curtis(Rl)
        coords, var = pcoa(D)
        res = []
        tests = ([self.group_col] if self.group_col else []) + [c for c in self.meta_cols if self.kinds.get(c) == "categorical"]
        for c in tests:
            if c == self.group_col:
                lab = [self.group(k[0]) for k in ks]
            else:
                lab = [self.meta.get(k[0], {}).get(c, "") for k in ks]
            ok = [i for i, x in enumerate(lab) if not is_missing(x) and x != "NA"]
            if len(ok) < 6:
                continue
            f, r2, p = permanova(D[np.ix_(ok, ok)], [lab[i] for i in ok], perms=self.a.permutations)
            res.append([c, len(ok), len(set(lab[i] for i in ok)), f, r2, p])
        return ks, coords, var, res

    def paired(self):
        if not {"sRNA-Int", "sRNA-Ext"} <= set(self.types):
            return []
        inc = set(self.keys("sRNA-Int")) | set(self.keys("sRNA-Ext"))
        both = [s for s in self.samples if (s, "sRNA-Int") in inc and (s, "sRNA-Ext") in inc]
        rows = []
        for k in self.metrics:
            a = [self.m[(s, "sRNA-Int")].get(k, float("nan")) for s in both]
            b = [self.m[(s, "sRNA-Ext")].get(k, float("nan")) for s in both]
            ok = [i for i in range(len(both)) if not (math.isnan(a[i]) or math.isnan(b[i]))]
            md, p, n = wilcoxon_paired([a[i] for i in ok], [b[i] for i in ok])
            rows.append([k, n, float(np.median([a[i] for i in ok])) if ok else float("nan"),
                         float(np.median([b[i] for i in ok])) if ok else float("nan"), md, p])
        q = bh([r[5] for r in rows])
        for r, qq in zip(rows, q):
            r.append(qq)
        return rows

    def covariates(self):
        rows = []
        for c in self.meta_cols:
            by = collections.OrderedDict((g, []) for g in self.groups)
            for s in self.meta:
                g = self.group(s)
                v = self.meta[s].get(c, "")
                if g in by and not is_missing(v):
                    by[g].append(v)
            if self.kinds.get(c) == "numeric":
                num_ = [[float(x) for x in v] for v in by.values()]
                h, p = kruskal(num_)
                desc = ["%.1f (%.1f–%.1f), n=%d" % (np.median(v), np.percentile(v, 25), np.percentile(v, 75), len(v))
                        if v else "–" for v in num_]
                rows.append([c, "numeric", "Kruskal-Wallis", p] + desc)
            elif self.kinds.get(c) == "categorical":
                levels = sorted({x for v in by.values() for x in v})
                if len(levels) > 12:
                    continue
                tab = np.array([[sum(1 for x in v if x == l) for l in levels] for v in by.values()], float)
                p = chi2_test(tab)
                desc = [", ".join("%s %d" % (l, int(n)) for l, n in zip(levels, row) if n) or "–" for row in tab]
                rows.append([c, "categorical", "chi-square", p] + desc)
        return rows

    # ---------------- matched comparisons (e.g. site)
    def meta_val(self, sample, col):
        v = self.meta.get(sample, {}).get(col, "")
        return "NA" if is_missing(v) else v.strip()

    def block_of(self, k):
        return self.meta_val(k[0], self.a.match_col) if self.a.match_col else "NA"

    def block_design(self, t):
        """per non-reference group: blocks with the group, blocks shared with the reference, design type
        (nested: <= 2 shared blocks, e.g. site; crossed: >= 3, e.g. families)"""
        col = self.a.match_col
        if not col or len(self.groups) < 2:
            return []
        ref = self.a.reference or self.groups[0]
        ks = [k for k in self.keys(t) if self.block_of(k) != "NA"]
        rb = {self.block_of(k) for k in ks if self.group(k[0]) == ref}
        out = []
        for g in self.groups:
            if g == ref:
                continue
            gb = {self.block_of(k) for k in ks if self.group(k[0]) == g}
            sh = sorted(gb & rb)
            na = sum(1 for k in ks if self.group(k[0]) == g and self.block_of(k) in sh)
            nb = sum(1 for k in ks if self.group(k[0]) == ref and self.block_of(k) in sh)
            out.append({"group": g, "ref": ref, "blocks": len(gb), "shared": sh, "n_a": na, "n_b": nb,
                        "n_a_all": sum(1 for k in self.keys(t) if self.group(k[0]) == g),
                        "n_b_all": sum(1 for k in self.keys(t) if self.group(k[0]) == ref),
                        "design": "none" if not sh else "nested" if len(sh) <= 2 else "crossed"})
        return out

    def crossed(self, t=None):
        types = [t] if t else self.types
        return any(d["design"] == "crossed" for tt in types for d in self.block_design(tt))

    def matched_sets(self, t):
        """[(label, keys_a, keys_b, kind)]: per group 'pooled' (all vs all reference) and 'within' (only blocks that
        contain both, tested within blocks); reference-only block contrasts ('block') when the reference has <= 6
        block levels (site effect), else one 'block_all' entry (heterogeneity across many blocks, e.g. families)"""
        col = self.a.match_col
        if not col or len(self.groups) < 2:
            return []
        ref = self.a.reference or self.groups[0]
        ks = self.keys(t)
        out = []
        for d in self.block_design(t):
            g = d["group"]
            out.append(("%s vs %s (pooled)" % (g, ref), [k for k in ks if self.group(k[0]) == g],
                        [k for k in ks if self.group(k[0]) == ref], "pooled"))
            if d["shared"]:
                sh = set(d["shared"])
                where = ("%s: %s" % (col, ",".join(d["shared"]))) if len(sh) <= 3 else "%d %s levels" % (len(sh), col)
                out.append(("%s vs %s (within %s)" % (g, ref, where),
                            [k for k in ks if self.group(k[0]) == g and self.block_of(k) in sh],
                            [k for k in ks if self.group(k[0]) == ref and self.block_of(k) in sh], "within"))
        lv_ref = sorted({self.block_of(k) for k in ks if self.group(k[0]) == ref} - {"NA"})
        lv_ref = [l for l in lv_ref if sum(1 for k in ks if self.group(k[0]) == ref and self.block_of(k) == l) >= 2]
        if 2 <= len(lv_ref) <= 6:
            for i, a_ in enumerate(lv_ref):
                for b_ in lv_ref[i + 1:]:
                    ka = [k for k in ks if self.group(k[0]) == ref and self.block_of(k) == a_]
                    kb = [k for k in ks if self.group(k[0]) == ref and self.block_of(k) == b_]
                    out.append(("%s %s vs %s %s (%s effect within %s)" % (ref, a_, ref, b_, col, ref), ka, kb, "block"))
        return out

    def block_effect(self, t, mx="genus"):
        """heterogeneity of the reference group across block levels (levels with >= 2 libraries): PERMANOVA of
        Bray-Curtis (R², p) and the number of metrics differing (Kruskal-Wallis, q < SIG_Q) — the poolability check"""
        ref = self.a.reference or (self.groups[0] if self.groups else None)
        ks = [k for k in self.keys(t) if self.group(k[0]) == ref and self.block_of(k) != "NA"]
        cnt = collections.Counter(self.block_of(k) for k in ks)
        ks = [k for k in ks if cnt[self.block_of(k)] >= 2]
        scope = "reference"
        if len({self.block_of(k) for k in ks}) < 2 or len(ks) < 6:
            # e.g. one control per family: test the block on all libraries, permuted within groups
            ks = [k for k in self.keys(t) if self.block_of(k) != "NA" and self.group(k[0]) != "NA"]
            cnt = collections.Counter(self.block_of(k) for k in ks)
            ks = [k for k in ks if cnt[self.block_of(k)] >= 2]
            scope = "all"
            if len({self.block_of(k) for k in ks}) < 2 or len(ks) < 6:
                return None
        ks2, feats, Cm, Rl = self.rel_matrix(t, mx, min_prev=0.0, min_reads=1)
        idx = [i for i, k in enumerate(ks2) if k in set(ks)]
        pm = (float("nan"),) * 3
        if len(idx) >= 6:
            D = bray_curtis(Rl[idx])
            lab = [self.block_of(ks2[i]) for i in idx]
            grp = [self.group(ks2[i][0]) for i in idx]
            pm = (permanova_partial(D, lab, grp, perms=self.a.permutations)
                  if scope == "all" else permanova(D, lab, perms=self.a.permutations))
            if scope == "all" and not isnum(pm[2]):
                # one library per group and block (e.g. one patient + one control per family): block after group
                # has no residual df -> marginal block test, permuted within groups (balanced: block ⟂ group)
                pm = permanova(D, lab, perms=self.a.permutations, strata=grp)
                scope = "all_marginal"
        if scope.startswith("all"):
            return {"libraries": len(ks), "levels": len({self.block_of(k) for k in ks}), "R2": pm[1], "p": pm[2],
                    "metrics_q": 0, "metrics_tested": 0, "scope": scope}
        ps = []
        for m in self.metrics:
            by = collections.defaultdict(list)
            for k in ks:
                v = self.m[k].get(m, float("nan"))
                if isnum(v):
                    by[self.block_of(k)].append(v)
            ps.append(kruskal(list(by.values()))[1])
        q = bh(ps)
        return {"libraries": len(ks), "levels": len({self.block_of(k) for k in ks}), "R2": pm[1], "p": pm[2],
                "metrics_q": int(np.nansum(q < SIG_Q)), "metrics_tested": int(np.sum(~np.isnan(q))), "scope": "reference"}

    def two_set_tests(self, t, ka, kb, blocked=False):
        """metrics and genus / miRNA tests between two sets of libraries; Bray-Curtis PERMANOVA.
        blocked: tests within blocks (van Elteren; Mann-Whitney with one block; PERMANOVA permuted within blocks)"""
        if not hasattr(self, "_tst"):
            self._tst = {}
        ck = (t, tuple(ka), tuple(kb), blocked)
        if ck in self._tst:
            return self._tst[ck]
        bl = {k: self.block_of(k) for k in list(ka) + list(kb)} if blocked else {}

        def test(xa, ka_, xb, kb_):
            if blocked:
                auc, p, _ = van_elteren(xa, [bl[k] for k in ka_], xb, [bl[k] for k in kb_])
                return auc, p
            return mannwhitney(xa, xb)
        res = {"metrics": [], "genus": [], "mirna": [], "permanova": []}
        res.update({mx: [] for mx in FUNC_MX + ("rfam_family",)})
        for k in self.metrics:
            pa = [(x, self.m[x].get(k, float("nan"))) for x in ka]
            pb = [(x, self.m[x].get(k, float("nan"))) for x in kb]
            pa, pb = [(x, v) for x, v in pa if not math.isnan(v)], [(x, v) for x, v in pb if not math.isnan(v)]
            a, b = [v for _, v in pa], [v for _, v in pb]
            auc, p = test(a, [x for x, _ in pa], b, [x for x, _ in pb]) if a and b else (float("nan"), float("nan"))
            res["metrics"].append([k, len(a), len(b), float(np.median(a)) if a else float("nan"),
                                   float(np.median(b)) if b else float("nan"), auc, p])
        for mx in ("genus", "mirna", "rfam_family") + FUNC_MX:
            floor = self.a.min_mirna_reads if mx == "mirna" else 1e-9
            keys = [x for x in ka + kb if sum(self.vec[x][mx].values()) >= floor]
            feats = sorted({f for x in keys for f in self.vec[x][mx]})
            if len(keys) < 4 or not feats:
                continue
            C = np.array([[self.vec[x][mx].get(f, 0.0) for f in feats] for x in keys])
            Rl = C / C.sum(axis=1, keepdims=True)
            keep = (C >= 10).mean(axis=0) >= 0.1
            ina = np.array([x in set(ka) for x in keys])
            ka_ = [x for x in keys if x in set(ka)]
            kb_ = [x for x in keys if x not in set(ka)]
            for j in np.where(keep)[0]:
                auc, p = test(Rl[ina, j], ka_, Rl[~ina, j], kb_)
                res[mx].append([feats[j], float((C[:, j] > 0).mean()), float(Rl[ina, j].mean()),
                                float(Rl[~ina, j].mean()), auc, p])
            if mx == "genus" and ina.sum() >= 3 and (~ina).sum() >= 3:
                strata = [bl[x] for x in keys] if blocked and len({bl[x] for x in keys}) > 1 else None
                f, r2, p = permanova(bray_curtis(Rl), ["a" if x else "b" for x in ina], perms=self.a.permutations,
                                     strata=strata)
                res["permanova"] = [int(ina.sum()), int((~ina).sum()), f, r2, p]
        for key, pi in (("metrics", 6), ("genus", 5), ("mirna", 5), ("rfam_family", 5)) + tuple((mx, 5) for mx in FUNC_MX):
            q = bh([r[pi] for r in res[key]])
            for r, qq in zip(res[key], q):
                r.append(qq)
            res[key].sort(key=lambda r: (math.isnan(r[pi]), r[pi]))
        self._tst[ck] = res
        return res

    # ---------------- antisense
    def asidx(self, t, db, keys=None, min_reads=10):
        """antisense index AS / (S + AS) per library and term (NaN below min_reads); library CDS average"""
        ks = keys if keys is not None else sorted(self.keys(t))
        ks = [k for k in ks if self.vec[k].get("fn_" + db)]
        terms = sorted({f for k in ks for f in self.vec[k]["fn_" + db]} | {f for k in ks for f in self.vec[k]["fnas_" + db]})
        if not ks or not terms:
            return ks, [], np.zeros((0, 0)), np.zeros(0)
        S = np.array([[self.vec[k]["fn_" + db].get(f, 0.0) for f in terms] for k in ks])
        A = np.array([[self.vec[k]["fnas_" + db].get(f, 0.0) for f in terms] for k in ks])
        T = S + A
        X = np.where(T >= min_reads, A / np.where(T > 0, T, 1), np.nan)
        lib = A.sum(axis=1) / np.maximum(T.sum(axis=1), 1e-9)
        return ks, terms, X, lib

    def antisense_hotspots(self, t, db, min_libs=6):
        """terms with consistently more (or less) antisense than the CDS average of their library:
        median log2(index_term / index_library) over libraries, one-sample Wilcoxon vs 0, BH-FDR"""
        ks, terms, X, lib = self.asidx(t, db)
        rows = []
        for j, f in enumerate(terms):
            ok = ~np.isnan(X[:, j])
            if ok.sum() < min_libs:
                continue
            lr = np.log2((X[ok, j] + 0.01) / (lib[ok] + 0.01))
            md, p, n = wilcoxon_paired(lr, np.zeros(len(lr)))
            rows.append([f, TERM_NAMES.get(f, ""), int(ok.sum()), float(np.median(X[ok, j])), float(np.median(lib[ok])),
                         float(np.median(lr)), float((lr > 0).mean()), p])
        q = bh([r[7] for r in rows])
        for r, qq in zip(rows, q):
            r.append(qq)
        rows.sort(key=lambda r: (math.isnan(r[7]), r[7]))
        return rows

    def asidx_tests(self, t, db, ka=None, kb=None, min_frac=0.5):
        """group tests (Kruskal-Wallis) of the antisense index per term; with ka/kb: Mann-Whitney between two sets"""
        two = ka is not None
        keys = (list(ka) + list(kb)) if two else None
        ks, terms, X, lib = self.asidx(t, db, keys=keys)
        rows = []
        if not terms:
            return rows
        lab = [("a" if k in set(ka) else "b") if two else self.group(k[0]) for k in ks]
        levels = ["a", "b"] if two else self.groups
        for j, f in enumerate(terms):
            ok = ~np.isnan(X[:, j])
            if ok.mean() < min_frac or ok.sum() < 6:
                continue
            by = [[X[i, j] for i in range(len(ks)) if ok[i] and lab[i] == g] for g in levels]
            if two:
                auc, p = mannwhitney(by[0], by[1])
                rows.append([f, TERM_NAMES.get(f, ""), int(ok.sum())] +
                            [float(np.median(b)) if b else float("nan") for b in by] + [auc, p])
            else:
                h, p = kruskal(by)
                rows.append([f, TERM_NAMES.get(f, ""), int(ok.sum()), h, p] +
                            [float(np.median(b)) if b else float("nan") for b in by])
        pi = len(rows[0]) - 1 if (rows and two) else 4
        q = bh([r[pi] for r in rows])
        for r, qq in zip(rows, q):
            r.append(qq)
        rows.sort(key=lambda r: (math.isnan(r[pi]), r[pi]))
        return rows

    def nmds_result(self, t, mx):
        """(keys, NMDS coordinates, stress, PERMANOVA of the group) on Bray-Curtis of the profile mx (cached)"""
        if not hasattr(self, "_nmds"):
            self._nmds, self._nmds_out = {}, {}
        if (t, mx) not in self._nmds:
            res = None
            if not (mx == "mirna" and t not in self.mirna_types()):
                ks, feats, Cm, Rl = self.rel_matrix(t, mx, min_prev=0.0, min_reads=1)
                if len(ks) >= 6 and feats:
                    D = bray_curtis(Rl)
                    lab = [self.group(k[0]) for k in ks]
                    ok = [i for i, g in enumerate(lab) if g != "NA"]
                    pm = permanova(D[np.ix_(ok, ok)], [lab[i] for i in ok], perms=self.a.permutations) \
                        if len(set(lab[i] for i in ok)) > 1 else (float("nan"),) * 3
                    # samples with a completely separate composition make the NMDS degenerate: ordinate without
                    # them (PERMANOVA above keeps all samples) and name them on the plot
                    iso = isolated_samples(D)
                    if iso:
                        self._nmds_out[(t, mx)] = [ks[i][0] for i in iso]
                        sys.stderr.write("[cohort] NMDS %s %s: left out (isolated composition): %s\n" % (
                            t, mx, ", ".join(self._nmds_out[(t, mx)])))
                        kk = [i for i in range(len(ks)) if i not in iso]
                        ks, D = [ks[i] for i in kk], D[np.ix_(kk, kk)]
                    X, st = nmds(D, k=2, starts=self.a.nmds_starts)
                    res = (ks, X, st, pm)
            self._nmds[(t, mx)] = res
        return self._nmds[(t, mx)]

    # ---------------- covariate-adjusted models
    def model_specs(self):
        """[(label, formula, variant)]: comma-separated models, '+'-separated terms. With a block column in a
        crossed design (e.g. families) every model is fitted twice: 'cluster' = pooled with block-clustered SE,
        'within' = within-block fixed effects (block-centred outcome and terms); otherwise HC3 ('')."""
        if not hasattr(self, "_specs"):
            forms = [m.strip() for m in (self.a.model or "").split(",") if m.strip()]
            out = []
            for i, f in enumerate(forms):
                if self.a.match_col and self.crossed():
                    out += [("M%d pooled" % (i + 1), f, "cluster"), ("M%d within-%s" % (i + 1, self.a.match_col), f, "within")]
                else:
                    out.append(("M%d" % (i + 1), f, ""))
            self._specs = out
        return self._specs

    def models(self):
        return [f for _, f, _ in self.model_specs()]

    def model_label(self, mi):
        return self.model_specs()[mi][0]

    def model_variant(self, mi):
        return self.model_specs()[mi][2]

    def model_terms(self, mi=0):
        ms = self.models()
        if not ms or mi >= len(ms):
            return []
        return [t.strip() for t in ms[mi].split("|", 1)[0].split("+") if t.strip()]

    def model_subset(self, mi=0):
        """'terms | column=value1/value2' fits the model on the libraries of these samples only
        (e.g. 'levodopa+bristol+sex+age | group=PD') -> (column, {values}) or None"""
        ms = self.models()
        if not ms or mi >= len(ms) or "|" not in ms[mi]:
            return None
        cond = ms[mi].split("|", 1)[1].strip()
        if "=" not in cond:
            return None
        col, vals = cond.split("=", 1)
        return col.strip(), {v.strip() for v in vals.split("/") if v.strip()}

    def model_focus(self, mi=0):
        """term of interest: the first term of the formula (group in 'group+site+…', stool consistency in
        'bristol+laxative+group+…', a medication in 'levodopa+… | group=PD')"""
        tm = self.model_terms(mi)
        return tm[0] if tm else None

    def model_data(self, t, mi=0):
        """(libraries, X, column names, column owners, info) for the model in fraction t: complete cases of the
        model terms (metadata columns or library metrics such as qc_lRNA_RIN)"""
        if not hasattr(self, "_mdata"):
            self._mdata = {}
        if (t, mi) in self._mdata:
            return self._mdata[(t, mi)]
        terms = self.model_terms(mi)
        res = None
        if terms:
            kinds, rows, keys = {}, [], []
            for tm in terms:
                if tm == self.group_col:
                    kinds[tm] = "categorical"
                elif tm in self.meta_cols:
                    kinds[tm] = self.kinds.get(tm, "categorical")
                elif tm in self.metrics:
                    kinds[tm] = "numeric"
                else:
                    sys.stderr.write("[cohort] model term '%s' is neither a metadata column nor a metric — ignored\n" % tm)
            use = [tm for tm in terms if tm in kinds]
            var = self.model_variant(mi)
            subset = self.model_subset(mi)
            for k in sorted(self.keys(t)):
                if var and self.block_of(k) == "NA":
                    continue
                if subset:
                    sv = self.group(k[0]) if subset[0] == self.group_col else self.meta_val(k[0], subset[0])
                    if sv not in subset[1]:
                        continue
                r = {}
                for tm in use:
                    if tm == self.group_col:
                        v = self.group(k[0])
                    elif tm in self.meta_cols:
                        v = self.meta_val(k[0], tm)
                    else:
                        x = self.m[k].get(tm, float("nan"))
                        v = "NA" if not isnum(x) else repr(float(x))
                    r[tm] = v
                if all(v not in ("NA", "") for v in r.values()):
                    rows.append(r)
                    keys.append(k)
            if len(keys) >= max(8, len(use) + 4):
                ref = {self.group_col: self.a.reference or (self.groups[0] if self.groups else None)}
                X, names, owner, info = design_matrix(rows, use, kinds, ref)
                blocks = [self.block_of(k) for k in keys]
                if var == "within":                     # fixed block effects by the within transformation
                    Xw = demean_within(X[:, 1:], blocks)
                    keepc = [j for j in range(Xw.shape[1]) if np.abs(Xw[:, j]).max() > 1e-9]
                    for j in range(Xw.shape[1]):
                        if j not in keepc:
                            info.append("no within-%s variation: %s" % (self.a.match_col, names[j + 1]))
                    X = np.column_stack([np.zeros(len(keys))] + [Xw[:, j] for j in keepc])   # col 0 placeholder
                    names = ["(block means)"] + [names[j + 1] for j in keepc]
                    owner = ["(block means)"] + [owner[j + 1] for j in keepc]
                    info.append("%d %s levels, %d with more than one group" % (
                        len(set(blocks)), self.a.match_col,
                        sum(1 for b_ in set(blocks) if len({self.group(k[0]) for k, bb in zip(keys, blocks)
                                                            if bb == b_}) > 1)))
                if len(keys) - X.shape[1] >= 3 and X.shape[1] > 1:
                    res = (keys, X, names, owner, info, blocks, var)
        self._mdata[(t, mi)] = res
        return res

    def model_fit(self, t, what, mi=0):
        """fit the model for one outcome set: a count matrix (CLR), 'metrics' or 'asidx_<db>';
        returns rows [feature, name, n, term, estimate, se, t, p] + q (BH per term within the set)"""
        if not hasattr(self, "_mfit"):
            self._mfit = {}
        if (t, what, mi) in self._mfit:
            return self._mfit[(t, what, mi)]
        out = []
        md = self.model_data(t, mi)
        if md:
            keys, X, names, owner, info, blocks, var = md
            pos = {k: i for i, k in enumerate(keys)}

            def fit(Xs, Ys, rows_):
                if var == "within":
                    Ys = demean_within(Ys, [blocks[i] for i in rows_])
                    return ols_cluster(Xs[:, 1:], Ys, [blocks[i] for i in rows_])
                if var == "cluster":
                    return ols_cluster(Xs, Ys, [blocks[i] for i in rows_])
                return ols_hc3(Xs, Ys)
            off = 1 if var == "within" else 0              # 'within' fits without column 0
            if what == "metrics":
                skip = set(self.model_terms(mi)) | set(NOT_TESTED)
                feats = [m for m in self.metrics if m not in skip and not m.startswith("qc_")]
                Y = np.array([[self.m[k].get(f, float("nan")) for f in feats] for k in keys], float)
                outsets = [(feats, Y, None)]
            elif what.startswith("asidx_"):
                ks2, terms, Xa, lib = self.asidx(t, what[6:], keys=keys)
                outsets = [(terms, Xa, ks2)]
                feats = terms
            else:
                if what == "mirna" and t not in self.mirna_types():
                    outsets = []
                else:
                    floor = self.a.min_mirna_reads if what == "mirna" else 1e-9
                    kk = [k for k in keys if sum(self.vec[k][what].values()) >= floor]
                    feats = sorted({f for k in kk for f in self.vec[k][what]})
                    if kk and feats:
                        Cm = np.array([[self.vec[k][what].get(f, 0.0) for f in feats] for k in kk])
                        keep = (Cm >= 10).mean(axis=0) >= 0.1
                        feats = [f for f, ok in zip(feats, keep) if ok]
                        outsets = [(feats, clr(Cm[:, keep]), kk)] if feats and len(kk) >= X.shape[1] + 3 else []
                    else:
                        outsets = []
            for feats, Y, kk in outsets:
                rowsel = [pos[k] for k in kk] if kk is not None else list(range(len(keys)))
                Xs = X[rowsel]
                if Y.size == 0:
                    continue
                if var == "within" and kk is not None:      # Xs rows are already block-centred over all keys
                    Xs = np.column_stack([np.zeros(len(rowsel)), demean_within(Xs[:, 1:], [blocks[i] for i in rowsel])])
                if not np.isnan(Y).any():
                    if not design_ok(Xs, off):
                        continue
                    B, SE, T, P, df = fit(Xs, Y, rowsel)
                    for j, f in enumerate(feats):
                        for c in range(1, Xs.shape[1]):
                            out.append([f, TERM_NAMES.get(f, ""), len(rowsel), names[c], B[c - off, j], SE[c - off, j],
                                        T[c - off, j], P[c - off, j]])
                else:                                       # per feature on its complete cases
                    for j, f in enumerate(feats):
                        ok = ~np.isnan(Y[:, j])
                        if ok.sum() < Xs.shape[1] + 3 or np.nanstd(Y[ok, j]) == 0 or \
                                len(np.unique(Y[ok, j])) < 5 or (what.startswith("asidx_") and ok.sum() < MODEL_MIN_N_ASIDX):
                            continue
                        sel = [rowsel[i] for i in np.where(ok)[0]]
                        Xo = Xs[ok]
                        if var == "within":
                            Xo = np.column_stack([np.zeros(len(sel)), demean_within(Xo[:, 1:], [blocks[i] for i in sel])])
                        if np.linalg.matrix_rank(Xo[:, off:]) < Xo.shape[1] - off or not design_ok(Xo, off):
                            continue
                        B, SE, T, P, df = fit(Xo, Y[ok, j][:, None], sel)
                        for c in range(1, Xo.shape[1]):
                            out.append([f, TERM_NAMES.get(f, ""), int(ok.sum()), names[c], B[c - off, 0],
                                        SE[c - off, 0], T[c - off, 0], P[c - off, 0]])
            for term in sorted({r[3] for r in out}):
                sub = [r for r in out if r[3] == term]
                for r, q in zip(sub, bh([r[7] for r in sub])):
                    r.append(q)
            out.sort(key=lambda r: (r[3], math.isnan(r[7]), r[7]))
        self._mfit[(t, what, mi)] = out
        return out

    def group_terms(self, t, mi=0):
        md = self.model_data(t, mi)
        if not md:
            return []
        focus = self.model_focus(mi)
        return [n for n, o in zip(md[2], md[3]) if o == focus]

    def qc_summary(self):
        """per type and group: libraries, PASS / WARN / FAIL, excluded; chi-square of exclusion vs group"""
        rows, tests = [], []
        for t in self.types:
            tab = []
            for g in self.groups + (["NA"] if any(self.group(s) == "NA" for s in self.samples) else []):
                ks = [k for k in self.keys(t, all_=True) if self.group(k[0]) == g]
                c = collections.Counter(self.qc[k][0] for k in ks)
                ex = sum(1 for k in ks if k in self.excluded)
                rows.append([t, g, len(ks), c["PASS"], c["WARN"], c["FAIL"], ex])
                if g != "NA":
                    tab.append([len(ks) - ex, ex])
            p = chi2_test(np.array(tab, float)) if len(tab) > 1 else float("nan")
            tests.append([t, sum(r[6] for r in rows if r[0] == t), p])
        rule_counts = collections.Counter()
        for k, (lv, hits) in self.qc.items():
            for h in hits:
                rule_counts[(k[1], h[0], h[1])] += 1
        return rows, tests, rule_counts

    def integrity_vs_qc(self):
        rows = []
        for t in self.types:
            ks = self.keys(t)
            qcs = [k for k in self.metrics if k.startswith("qc_")]
            ints = [k for k in self.metrics if k.startswith("integrity_")]
            for qk in qcs:
                for ik in ints:
                    x = [self.m[k].get(qk, float("nan")) for k in ks]
                    y = [self.m[k].get(ik, float("nan")) for k in ks]
                    rho, p, n = spearman(x, y)
                    rows.append([t, qk[3:], ik[10:], n, rho, p])
        # FDR within each wet-lab QC measure (RIN, Qubit, ...): they answer different questions
        for qc in sorted({r[1] for r in rows}):
            sub = [r for r in rows if r[1] == qc]
            for r, qq in zip(sub, bh([r[5] for r in sub])):
                r.append(qq)
        return rows


def chi2_test(tab):
    tab = tab[tab.sum(axis=1) > 0][:, tab.sum(axis=0) > 0]
    if tab.shape[0] < 2 or tab.shape[1] < 2:
        return float("nan")
    e = tab.sum(axis=1, keepdims=True) * tab.sum(axis=0, keepdims=True) / tab.sum()
    x2 = float(((tab - e) ** 2 / np.where(e > 0, e, 1)).sum())
    return _gammq((tab.shape[0] - 1) * (tab.shape[1] - 1) / 2.0, x2 / 2.0)


# ----------------------------------------------------------------------------- report
plt = R.plt
EXPLAIN = {
    "overview":
        "Samples with pipeline output, matched to the metadata by ID. Groups are compared on all covariates: a "
        "covariate that differs between groups (e.g. constipation, medication, age) can explain apparent group "
        "differences in sRNA and must be included in adjusted models.",
    "qcflags":
        "Each library is checked against the QC rules (thresholds on depth, trimming, adapter dimers, mapping, "
        "FastQC, RIN, completeness; manual decisions from QC_EXCLUDE). Libraries with a FAIL are left out of all "
        "statistics and exported matrices of their fraction; a sample can stay in with the other fraction. If "
        "exclusion differs between groups (χ² p small), group comparisons can be biased.",
    "findings":
        "Automatic list of results that pass the thresholds, by area. It is a guide to the report pages, not a result "
        "list: tests are unadjusted for covariates and many are run; group differences that also appear in the "
        "site contrast within the reference group, or in a covariate that differs between groups, are likely confounded.",
    "integrity_sum":
        "One row per integrity indicator: typical values per fraction, whether the fractions differ (paired test), "
        "whether groups differ within a fraction, and how the indicator relates to the wet-lab RIN of the long-RNA "
        "fraction. The bullet points below the table spell out the significant and borderline results.",
    "mirna_top":
        "The most abundant human miRNAs, as share of all human miRNA reads of a library, per group. Only libraries "
        "with enough miRNA reads are used. Stars give the Kruskal-Wallis p (* < 0.05, ** < 0.01, *** < 0.001); red "
        "bold titles pass the FDR threshold across all tested miRNAs of the fraction.",
    "qc":
        "Key library metrics per group and fraction: read-fate shares (what happens to the reads, see Read "
        "definitions), depth, mapping to the own contigs, novel loci and CRISPR reads. Differences in human or rRNA shares can reflect biology (shedding of "
        "epithelial cells, microbial load) but also sampling and extraction; check them before interpreting "
        "composition. Boxes: median and IQR, points: samples.",
    "integrity":
        "Integrity indicators per group and their relation to the wet-lab RIN of the long-RNA fraction of the same "
        "extraction. A group difference in integrity is a technical covariate to include in models; a correlation "
        "with RIN supports the indicator, a missing correlation in sRNA-Ext is expected (cell-free RNA is fragmented).",
    "taxa":
        "Mean relative abundance of rRNA-derived genera (Bracken, GTDB) per group and fraction. sRNA taxonomy follows "
        "genome abundance (metaG) more than transcription, so it is a proxy for which taxa release small RNAs.",
    "nmds":
        "NMDS places samples so that the rank order of their Bray-Curtis dissimilarities is kept as well as possible in "
        "two dimensions; unlike PCoA it makes no assumption about the scale of the dissimilarities. Stress (Kruskal "
        "stress-1) measures the misfit: < 0.1 good, 0.1–0.2 fair, > 0.2 poor (then the 2-D picture is unreliable). "
        "Axes have no units and no variance share; separation of the group ellipses is tested by the PERMANOVA.",
    "beta":
        "Principal coordinates of Bray-Curtis dissimilarities of genus profiles, coloured by group. PERMANOVA R² is the "
        "share of between-sample variation explained by a variable (999 permutations; marginal, not adjusted). "
        "Different group dispersions can also give a significant PERMANOVA.",
    "mirna":
        "Human miRNAs per group: share of reads and share of miRNA reads from the origin marker sets (samples with "
        "≥ 20 miRNA reads). More epithelial markers: more shed colonocytes; immune/blood markers: inflammation, bleeding.",
    "matched":
        "A block is a second grouping of the samples (clinical site, family, batch, …). Nested design (few shared "
        "levels, e.g. all PD from one site): the pooled comparison mixes group and block, the within comparison "
        "uses only controls of the same level(s). Crossed design (many levels with both groups, e.g. families): the "
        "within comparison contrasts patients and controls of the same family and removes family effects. The "
        "poolability check tests the reference group across block levels: without a block effect, pooled results "
        "are trustworthy; with one, prefer the within comparisons.",
    "classes":
        "What the non-human small RNAs are, library by library: rRNA and tRNA fragments dominate; snoRNAs, sRNAs, "
        "riboswitches / leaders (5′ UTR elements), ribozymes, introns and other Rfam classes are minor but specific. "
        "Classes from the classification (rRNA by RiboDetector, tRNA, miRNA-like, unknown) plus the Rfam hits split by "
        "Rfam type. Within-sample differences between the fractions are on the Int vs Ext page.",
    "antisense":
        "Reads antisense to protein-coding genes can be antisense RNAs (regulatory, base-pairing with the mRNA), "
        "transcripts of overlapping genes on the other strand, or library strand errors. The strand error is "
        "measured on rRNA/tRNA genes (antisense reads there are artefacts) and is usually far below the CDS "
        "antisense share. Hotspots are functions whose genes are more antisense-covered than the average gene of "
        "the same library, consistently across libraries (one-sample Wilcoxon).",
    "functions":
        "sRNA reads on protein-coding genes are mostly fragments of mRNAs (sense) — their functional profile mirrors "
        "which genes are transcribed and degraded; antisense reads (exported, not tested) can be antisense RNAs. "
        "Shares are within each database (KEGG pathway, module, KO, Pfam, COG; global KEGG maps left out). Terms "
        "with ≥ 10 reads in ≥ 10 % of the libraries are tested. Term names need --term-names (COHORT_TERM_NAMES).",
    "models":
        "Linear models per feature with the group and the covariates of COHORT_MODEL (e.g. site, sex, age, "
        "constipation). Abundances are centred-log-ratio (CLR) transformed, so estimates are differences in log "
        "abundance relative to the sample's average feature, adjusted for all other terms (the model type of "
        "MaAsLin2, here in Python). Robust (HC3) standard errors; BH-FDR within each outcome set and term. Terms "
        "that are fully confounded with others are dropped and named in red. Related samples (families, repeated "
        "visits) are not modelled as random effects.",
    "models_overview":
        "All adjusted models at a glance. For each model the term of interest is shown — the first term of the formula "
        "(the group contrasts in 'group+…', otherwise e.g. stool consistency, laxatives, a medication within one group) "
        "— with the number of features that change at q below the threshold per outcome set and the strongest "
        "results. Compare rows: a group difference that disappears when stool consistency or constipation is added "
        "is likely explained by them. Models restricted to a subgroup are marked [group]. Small levels "
        "(few samples) give large, fragile estimates.",
    "paired":
        "sRNA-Int (cell fraction) vs sRNA-Ext (cell-free) of the same samples, paired Wilcoxon signed-rank test with "
        "BH-FDR. Consistent differences describe what is released into the cell-free fraction.",
    "omics":
        "Per-sample Spearman correlation of sRNA with metaG and metaT (comparison step) per group. A group with lower "
        "sRNA–metaG correlation releases sRNA less in proportion to taxon abundance.",
    "features":
        "Exploratory group tests (Kruskal-Wallis, BH-FDR) for metrics and for genera / miRNAs present in ≥ 10 % of "
        "samples with ≥ 10 reads. Not adjusted for covariates or compositionality — candidates for adjusted models "
        "(see the adjusted models page, --model), not results.",
}
EXPLAIN.update({
    "highlights_crispr":
        "CRISPR systems expressed across the cohort. Arrays of different samples (each with its own assembly) are "
        "joined when they share the repeat sequence, so one row is one CRISPR system, usually one species or strain "
        "group. Expressed arrays in many samples point to CRISPR–Cas systems that are active in the gut; a tracrRNA "
        "candidate (type II / V) is an sRNA cluster in the locus pairing with the repeat (per-sample report).",
    "contents":
        "The report is organised in sections: first the cohort and the key findings, then data quality (read before "
        "interpreting differences), what the libraries contain, unadjusted group comparisons, adjusted models, "
        "heterogeneity (subgroups of patients), notable single findings across samples, and the methods. The section "
        "of each page is shown in its upper right corner.",
    "models_intro":
        "Each model asks one question, given by its first term (term of interest): the group contrasts, or e.g. stool "
        "consistency, a medication or a clinical score. The other terms are adjustments. A subset restricts the model "
        "to some samples (e.g. PD only). Numeric terms are per standard deviation. Models need ≥ 20 complete "
        "libraries and ≥ 3 per level; 'design notes' list levels or terms that were dropped. Add a question text per "
        "model with --model-notes / COHORT_MODEL_NOTES ('M9=PD motor severity; …').",
    "highlights_crispr_group":
        "CRISPR systems that are expressed in one group only, or much more often than in the other groups. Expression "
        "needs the host taxon: a group-specific system usually reflects a strain or species that is more common in that "
        "group (compare its abundance in metaG), rarely a change of CRISPR activity itself. With ~60 samples per group, "
        "'only' in 3–4 samples is not significant on its own; trust rows with q < 0.1.",
    "highlights_group":
        "Presence / absence per group for every feature class at once: CRISPR systems, novel sRNA clusters, Rfam "
        "families of novel loci (per sample), and genera, species, Rfam families, human miRNAs and KOs per fraction "
        "(present = ≥ 10 reads). A feature found in many samples of one group and rarely in the others is a "
        "candidate marker; q < 0.1 survives the multiple testing within its class. Presence depends on sequencing "
        "depth and on taxon abundance — check depth by group (QC pages) and abundance in metaG; for quantitative "
        "differences use the adjusted models.",
    "heterogeneity":
        "Complex diseases may consist of subtypes: then only some patients differ from the controls, each in its own "
        "way, and group means hardly move (adjusted models find little). A larger spread of the patients around their "
        "centroid, or a higher outlier load (share of features outside the normal range of the controls), points to "
        "such heterogeneity ('Anna Karenina' pattern of dysbiosis). Depth and RNA quality also raise the spread — "
        "compare with the QC pages.",
    "heterogeneity_features":
        "Features that leave the control range in a subset of the patients (e.g. 20 % of PD above the 95th percentile "
        "of the controls, 5 % expected). Such features can define candidate subgroups: compare the listed patients "
        "with the clinical data (severity, duration, treatment, constipation). Full lists incl. patient IDs: "
        "stats/cohort.heterogeneity.outlier_features.tsv.",
    "highlights_novel":
        "Novel sRNA loci (no Rfam hit, module novel) of all libraries, clustered across samples by shared sequence. "
        "A cluster found in many samples is a conserved, reproducibly expressed sRNA of a common gut taxon and the "
        "best candidate for a new regulatory RNA; check the representative sequence (FASTA in stats/) with Rfam, "
        "BLAST or the structure. Only the top loci per library enter (--highlights-top); locus = intergenic or "
        "antisense to a CDS.",
    "highlights_rfam":
        "Known ncRNA families (Rfam) among the novel loci, without rRNA and tRNA: e.g. riboswitches, sRNAs, "
        "leaders and ribozymes that the genome annotation did not include.",
    "highlights_chimera":
        "Exploratory (module chimera): chimeric reads may join two interacting RNAs, but without crosslinking most are "
        "library artefacts. Evidence for real duplexes is a group of pairs whose arms are more complementary than "
        "re-paired arms, and more inter-contig chimeras within one organism (bin) than expected. Pairs are "
        "candidates, not a measured interactome.",
})
for k, v in EXPLAIN.items():
    R.EXPLAIN["cohort_" + k] = v

TYPE_MARK = {"sRNA-Int": "o", "sRNA-Ext": "s"}


SIG_Q = 0.1        # q below this: bold in tables / plots (set from --sig-q)


def stars(p):
    if p is None or (isinstance(p, float) and math.isnan(p)):
        return ""
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""


def bold_rows(ax, flags):
    """bold the table rows (1-based after the header) whose flag is True"""
    if not ax.tables:
        return
    for (r, c), cell in ax.tables[-1].get_celld().items():
        if r > 0 and r - 1 < len(flags) and flags[r - 1]:
            cell.set_text_props(fontweight="bold")
            cell.set_facecolor("#fff4d6")


def isnum(x):
    return isinstance(x, (int, float)) and not math.isnan(x)


COHORT_NAME = "cohort"       # report title (--cohort-name / COHORT_NAME)


def cpage(pdf, title, subtitle, key):
    return R.page(pdf, COHORT_NAME, title, subtitle, explain="cohort_" + key)


def boxes(ax, C, metric, title, log=False, pct=False, types=None):
    """boxes per group (x) and type (side by side)"""
    types = types or C.types
    nt = len(types)
    w = 0.8 / max(nt, 1)
    any_ = False
    for gi, g in enumerate(C.groups):
        for ti_, t in enumerate(types):
            ti = C.types.index(t)                     # keep the colour of the fraction
            v = [C.m[k].get(metric, float("nan")) for k in C.keys(t) if C.group(k[0]) == g]
            v = [x * (100 if pct else 1) for x in v if not math.isnan(x) and (not log or x > 0)]
            if not v:
                continue
            any_ = True
            pos = gi + (ti_ - (nt - 1) / 2) * w
            ax.boxplot([v], positions=[pos], widths=w * 0.8, showfliers=False, patch_artist=True,
                       boxprops=dict(facecolor=R.TYPE_COLORS[ti] + "40", edgecolor=R.TYPE_COLORS[ti]),
                       medianprops=dict(color=R.TYPE_COLORS[ti]), whiskerprops=dict(color=R.TYPE_COLORS[ti]),
                       capprops=dict(color=R.TYPE_COLORS[ti]))
            jit = random.Random(hash((g, t)) & 0xffff)
            ax.scatter([pos + (jit.random() - 0.5) * w * 0.5 for _ in v], v, s=6, color=R.TYPE_COLORS[ti],
                       alpha=0.6, lw=0, zorder=3)
    ax.set_xticks(range(len(C.groups)))
    ax.set_xticklabels(C.groups, fontsize=7)
    ax.set_title(title, loc="left", fontsize=8)
    if log:
        ax.set_yscale("log")
    if not any_:
        R.note(ax, "no data")
        return
    # group test per fraction: Kruskal-Wallis p (stars) and BH q of the metrics table; bold if q < SIG_Q
    if len(C.groups) > 1:
        parts, sig = [], False
        for t in types:
            p, q = C.metric_pq(t, metric)
            if isnum(p):
                parts.append("%s p=%s%s" % (t.replace("sRNA-", ""), fmt(p, 2), stars(p)))
                sig = sig or (isnum(q) and q < SIG_Q)
        if parts:
            ax.text(0.99, 0.99, "  ".join(parts) + ("  (q<%g)" % SIG_Q if sig else ""), transform=ax.transAxes,
                    ha="right", va="top", fontsize=5.8, fontweight="bold" if sig else "normal",
                    color="#B00000" if sig else "#555")


def type_legend(fig, C, x=0.97, y=0.927, types=None):
    for i, t in enumerate(reversed(types or C.types)):
        fig.text(x - i * 0.075, y, "■ " + t, color=R.TYPE_COLORS[C.types.index(t)], fontsize=7.5, ha="right",
                 va="top")


def page_overview(pdf, C):
    fig = cpage(pdf, "cohort overview", "%d samples, %d sample-types · metadata: %s" % (
        len(C.samples), len(C.m), C.a.metadata or "none"), "overview")
    rows = []
    for g in C.groups + (["NA"] if any(C.group(s) == "NA" for s in C.samples) else []):
        rows.append([g] + ["%d / %d" % (sum(1 for k in C.keys(t) if C.group(k[0]) == g),
                                        sum(1 for k in C.keys(t, all_=True) if C.group(k[0]) == g)) for t in C.types] +
                    [sum(1 for s in C.samples if C.group(s) == g)])
    ax = fig.add_axes([0.03, 0.62, 0.40, 0.24])
    R.table(ax, rows, ["group"] + C.types + ["samples"], rowscale=1.3,
            title="libraries per group and fraction (in statistics / all)")
    # why libraries are missing from the statistics: main (FAIL) reason per excluded library
    why = collections.defaultdict(collections.Counter)
    for k in C.excluded:
        lv, hits = C.qc[k]
        fails = [h for h in hits if h[1] == lv] or hits
        reason = fails[0][0] if fails else "?"
        label = {"incomplete": "pipeline not finished", "manual": "manual exclusion (QC_EXCLUDE)",
                 "outlier": "outlier"}.get(reason, fails[0][5] if fails else reason)
        why[(k[1], C.group(k[0]))][label] += 1
    ax = fig.add_axes([0.47, 0.62, 0.50, 0.24])
    if why:
        wrows = [[t, g, sum(c.values()), "; ".join("%s: %d" % (r, n) for r, n in c.most_common())]
                 for (t, g), c in sorted(why.items(), key=lambda kv: (C.types.index(kv[0][0]), kv[0][1]))]
        R.table(ax, wrows, ["type", "group", "excluded", "reason (main QC failure)"], rowscale=1.3, fontsize=6.8,
                left=(0, 1, 3), col_widths=[0.15, 0.10, 0.11, 0.64],
                title="why libraries are not in the statistics (details: flagged libraries page)")
    else:
        R.note(ax, "all libraries are in the statistics")
    um = [s for s in C.samples if s not in C.meta]
    fig.text(0.03, 0.58, "pipeline samples without metadata: %d%s · metadata IDs without output: %d" % (
        len(um), (" (" + ", ".join(um[:8]) + (" …" if len(um) > 8 else "") + ")") if um else "", len(C.unmatched)),
        fontsize=7, color="#555")
    cov = C.covariates()
    per_first, per_page = 14, 30
    chunks = [cov[:per_first]] + [cov[i:i + per_page] for i in range(per_first, len(cov), per_page)]

    def cov_table(ax, rows, title):
        R.table(ax, [[r[0], r[2], fmt(r[3], 2)] + [R.short(x, 60) for x in r[4:]] for r in rows],
                ["covariate", "test", "p"] + C.groups, title=title,
                col_widths=[0.13, 0.09, 0.05] + [0.73 / max(len(C.groups), 1)] * len(C.groups), fontsize=6.5,
                left=tuple(range(3 + len(C.groups))), rowscale=1.15)
        bold_rows(ax, [isnum(r[3]) and r[3] < 0.05 for r in rows])

    ax = fig.add_axes([0.03, 0.06, 0.94, 0.48])
    if cov:
        cov_table(ax, chunks[0], "covariates by group (median (IQR) or counts)%s" % (
            " — %d of %d, continued on the next page" % (len(chunks[0]), len(cov)) if len(chunks) > 1 else ""))
    else:
        R.note(ax, "no metadata covariates (--metadata / --covariates)")
    pdf.savefig(fig)
    plt.close(fig)
    for i, ch in enumerate(chunks[1:], 2):
        fig = R.page(pdf, COHORT_NAME, "covariates by group (%d/%d)" % (i, len(chunks)))
        ax = fig.add_axes([0.03, 0.04, 0.94, 0.86])
        cov_table(ax, ch, "")
        pdf.savefig(fig)
        plt.close(fig)


def page_qc_flags(pdf, C):
    rows, tests, rc = C.qc_summary()
    fig = cpage(pdf, "library QC and exclusions", "rules: %s · excluded from statistics: %s" % (
        C.a.qc_rules or "scripts/qc_rules.tsv", {"FAIL": "FAIL", "WARN": "WARN + FAIL", "none": "none"}[C.a.qc_drop]),
        "qcflags")
    ax = fig.add_axes([0.03, 0.50, 0.45, 0.36])
    R.table(ax, rows, ["type", "group", "libraries", "PASS", "WARN", "FAIL", "excluded"], rowscale=1.2,
            title="QC verdict per fraction and group", left=(0, 1))
    ax = fig.add_axes([0.03, 0.34, 0.45, 0.12])
    R.table(ax, [[t, n, fmt(p, 2)] for t, n, p in tests], ["type", "excluded", "p (exclusion vs group, χ²)"],
            title="is exclusion balanced between groups?", left=(0,))
    ax = fig.add_axes([0.58, 0.34, 0.39, 0.52])
    items = sorted(rc.items(), key=lambda kv: -kv[1])[:20]
    if items:
        labels = ["%s · %s" % (t.replace("sRNA-", ""), r) for (t, r, lv), _ in items]
        ax.barh(range(len(items)), [n for _, n in items], color=[Q.LEVEL_COLOR[lv] for (t, r, lv), _ in items])
        ax.set_yticks(range(len(items)))
        ax.set_yticklabels(labels, fontsize=6.5)
        ax.invert_yaxis()
        ax.set_xlabel("libraries")
        ax.set_title("triggered rules (orange WARN, red FAIL)", loc="left", fontsize=8)
    else:
        R.note(ax, "no rule triggered")
    nflag = sum(1 for k, (lv, _) in C.qc.items() if lv != "PASS")
    fig.text(0.03, 0.29, "%d libraries excluded, %d flagged in total (WARN or FAIL) — all listed on the next page(s) "
             "and in tables/cohort.qc_flags.tsv" % (len(C.excluded), nflag), fontsize=7.5, color="#333")
    pdf.savefig(fig)
    plt.close(fig)
    page_flagged(pdf, C)


def page_flagged(pdf, C):
    """all libraries with WARN or FAIL and their reasons (paginated)"""
    lib = sorted([k for k, (lv, _) in C.qc.items() if lv != "PASS"],
                 key=lambda k: (-Q.LEVELS[C.qc[k][0]], k[0], k[1]))
    if not lib:
        return
    per = 34
    for i in range(0, len(lib), per):
        fig = R.page(pdf, COHORT_NAME, "flagged libraries (%d–%d of %d)" % (i + 1, min(i + per, len(lib)), len(lib)),
                     "FAIL = left out of the statistics (with COHORT_QC_DROP=FAIL); WARN = kept, flagged")
        rows = []
        for s_, t in lib[i:i + per]:
            lv, hits = C.qc[(s_, t)]
            why = "; ".join(("%s (z %s)" % (h[5], h[4]) if h[0] == "outlier" else h[5]) for h in hits)
            rows.append([s_, t, C.group(s_), C.meta_val(s_, C.a.match_col) if C.a.match_col else "", lv,
                         "no" if (s_, t) in C.excluded else "yes", R.short(why, 150)])
        ax = fig.add_axes([0.02, 0.04, 0.96, 0.86])
        R.table(ax, rows, ["sample", "type", "group", C.a.match_col or "", "qc", "in stats", "reasons"], fontsize=6.3,
                left=tuple(range(7)), col_widths=[0.07, 0.07, 0.05, 0.06, 0.045, 0.05, 0.655], rowscale=1.12)
        if ax.tables:
            for (r, c), cell in ax.tables[-1].get_celld().items():
                if r > 0 and c == 4:
                    cell.set_facecolor(Q.LEVEL_COLOR.get(cell.get_text().get_text(), "#fff") + "55")
        pdf.savefig(fig)
        plt.close(fig)


def _unused_excluded_table(fig, C):
    ex = sorted(C.excluded)
    ax = fig.add_axes([0.03, 0.04, 0.94, 0.25])
    R.table(ax, [[s, t, C.group(s), C.qc[(s, t)][0], R.short("; ".join(h[5] for h in C.qc[(s, t)][1]), 110)]
                 for s, t in ex[:12]], ["sample", "type", "group", "qc", "reasons"], fontsize=6.5,
            left=(0, 1, 2, 3, 4), col_widths=[0.09, 0.08, 0.08, 0.05, 0.70],
            title="excluded libraries%s" % (" (first 12 of %d, all in tables/cohort.qc_flags.tsv)" % len(ex)
                                            if len(ex) > 12 else ""))


def page_grid(pdf, C, key, title, subtitle, metrics, types=None):
    metrics = [m for m in metrics if m[0] in C.metrics]
    if not metrics or (types is not None and not types):
        return False
    fig = cpage(pdf, title, subtitle, key)
    type_legend(fig, C, types=types)
    nc = 4
    for i, (mk, lab, kw) in enumerate(metrics[:12]):
        r, c = divmod(i, nc)
        ax = fig.add_axes([0.05 + c * 0.24, 0.62 - r * 0.28, 0.19, 0.21])
        boxes(ax, C, mk, lab, types=types, **kw)
    pdf.savefig(fig)
    plt.close(fig)


def page_integrity_qc(pdf, C):
    qk = [k for k in C.metrics if k.startswith("qc_") and "rin" in k.lower()]
    if not qk:
        return False
    fig = cpage(pdf, "integrity vs wet-lab RIN", "Spearman per fraction; RIN of the long-RNA fraction", "integrity")
    type_legend(fig, C)
    pairs = [("integrity_tRNA_internal_fraction", "tRNA internal fragments"),
             ("integrity_tRNA_3p_CCA_fraction", "3′-CCA share"),
             ("integrity_tRNA_5p_end_entropy", "tRNA 5′-end entropy"),
             ("integrity_rRNA_5p_end_entropy", "rRNA 5′-end entropy")]
    for i, (ik, lab) in enumerate(pairs):
        ax = fig.add_axes([0.05 + i * 0.24, 0.45, 0.19, 0.38])
        for ti, t in enumerate(C.types):
            ks = C.keys(t)
            x = [C.m[k].get(qk[0], float("nan")) for k in ks]
            y = [C.m[k].get(ik, float("nan")) for k in ks]
            ax.scatter(x, y, s=8, color=R.TYPE_COLORS[ti], alpha=0.7, lw=0)
            rho, p, n = spearman(x, y)
            ax.text(0.02, 0.98 - ti * 0.07, "%s ρ=%s p=%s%s n=%d" % (t, fmt(rho, 2), fmt(p, 2), stars(p), n),
                    transform=ax.transAxes, fontsize=6, color=R.TYPE_COLORS[ti], va="top",
                    fontweight="bold" if isnum(p) and p < 0.05 else "normal")
        ax.set_xlabel("RIN", fontsize=7)
        ax.set_title(lab, loc="left", fontsize=8)
    rows = [r for r in C.integrity_vs_qc() if not math.isnan(r[5])]
    rows.sort(key=lambda r: r[5])
    ax = fig.add_axes([0.05, 0.04, 0.9, 0.30])
    R.table(ax, [[r[0], r[1], r[2], r[3], fmt(r[4], 2), fmt(r[5], 2), fmt(r[6], 2)] for r in rows[:10]],
            ["type", "wet-lab QC", "indicator", "n", "ρ", "p", "q"], fontsize=6.5, left=(0, 1, 2),
            title="strongest associations (all QC × indicator pairs; bold: q < %g)" % SIG_Q)
    bold_rows(ax, [isnum(r[6]) and r[6] < SIG_Q for r in rows[:10]])
    pdf.savefig(fig)
    plt.close(fig)


def page_taxa(pdf, C, mx="genus", top=15):
    fig = cpage(pdf, "rRNA %s by group" % mx, "mean relative abundance (Bracken); top %d" % top, "taxa")
    bars, labels = [], []
    for t in C.types:
        for g in C.groups:
            ks = [k for k in C.keys(t) if C.group(k[0]) == g and sum(C.vec[k][mx].values()) > 0]
            if not ks:
                continue
            tot = collections.Counter()
            for k in ks:
                s = sum(C.vec[k][mx].values())
                for f, v in C.vec[k][mx].items():
                    tot[f] += v / s / len(ks)
            bars.append(tot)
            labels.append("%s\n%s (n=%d)" % (g, t.replace("sRNA-", ""), len(ks)))
    if not bars:
        plt.close(fig)
        return False
    allt = collections.Counter()
    for b in bars:
        allt.update(b)
    feats = [f for f, _ in allt.most_common(top)]
    cmap = plt.get_cmap("tab20")
    ax = fig.add_axes([0.06, 0.08, 0.62, 0.76])
    bottom = np.zeros(len(bars))
    for i, f in enumerate(feats):
        v = np.array([100 * b.get(f, 0) for b in bars])
        ax.bar(range(len(bars)), v, bottom=bottom, color=cmap(i % 20), label=f, width=0.8)
        bottom += v
    ax.bar(range(len(bars)), 100 - bottom, bottom=bottom, color="#dddddd", label="other", width=0.8)
    ax.set_xticks(range(len(bars)))
    ax.set_xticklabels(labels, fontsize=6.5)
    ax.set_ylabel("% of rRNA reads")
    ax.set_ylim(0, 100)
    ax.legend(frameon=False, fontsize=6.5, loc="upper left", bbox_to_anchor=(1.01, 1))
    pdf.savefig(fig)
    plt.close(fig)


NMDS_LABELS = {"genus": "rRNA genera", "family": "rRNA families", "mirna": "human miRNAs", "rfam_family":
               "Rfam families (non-rRNA/tRNA)", "rna_class": "RNA classes", "fn_kegg_pathway": "KEGG pathways (CDS reads)",
               "fn_kegg_ko": "KOs (CDS reads)", "fn_pfam": "Pfam domains (CDS reads)"}


def _ellipse(ax, P, color):
    """95 % normal-theory ellipse of the points P (n x 2)"""
    if len(P) < 3:
        return
    cov = np.cov(P.T)
    w, v = np.linalg.eigh(cov)
    if np.any(w <= 0):
        return
    tt = np.linspace(0, 2 * np.pi, 100)
    circ = np.stack([np.cos(tt), np.sin(tt)])
    E = (v @ np.diag(np.sqrt(w * 5.991)) @ circ).T + P.mean(axis=0)
    ax.plot(E[:, 0], E[:, 1], color=color, lw=0.9, alpha=0.8)
    ax.scatter(*P.mean(axis=0), marker="+", s=60, color=color, lw=1.5, zorder=4)


def page_nmds(pdf, C):
    mxs = [m.strip() for m in C.a.nmds.split(",") if m.strip()]
    done = False
    marks = ["o", "^", "s", "D", "v", "P"]
    for mx in mxs:
        res = [(t, C.nmds_result(t, mx)) for t in C.types]
        res = [(t, r) for t, r in res if r]
        if not res:
            continue
        fig = cpage(pdf, "NMDS · %s" % NMDS_LABELS.get(mx, mx), "non-metric multidimensional scaling of Bray-Curtis "
                    "dissimilarities (2 dimensions, best of %d starts); ellipses: 95 %% per group; + = group centroid; "
                    "shapes: %s" % (C.a.nmds_starts + 1, C.a.match_col or "–"), "nmds")
        for i, (t, (ks, X, st, pm)) in enumerate(res[:2]):
            ax = fig.add_axes([0.06 + i * 0.48, 0.10, 0.40, 0.70])
            lv = sorted({C.meta_val(k[0], C.a.match_col) for k in ks}) if C.a.match_col else [""]
            for g in C.groups + ["NA"]:
                idx = [j for j, k in enumerate(ks) if C.group(k[0]) == g]
                if not idx:
                    continue
                for li, l in enumerate(lv):
                    sub = [j for j in idx if not C.a.match_col or C.meta_val(ks[j][0], C.a.match_col) == l]
                    if sub:
                        ax.scatter(X[sub, 0], X[sub, 1], s=18, color=C.gcolor[g], marker=marks[li % len(marks)],
                                   alpha=0.85, lw=0, label="%s%s (%d)" % (g, " · " + l if C.a.match_col else "", len(sub)))
                if g != "NA":
                    _ellipse(ax, X[idx], C.gcolor[g])
            quality = "good" if st < 0.1 else "fair" if st < 0.2 else "poor — interpret with care"
            sig = isnum(pm[2]) and pm[2] < 0.05
            ax.set_title("%s   stress %.3f (%s)" % (t, st, quality), loc="left", fontsize=9)
            ax.text(0.01, 0.01, "PERMANOVA group: R² %s, p %s" % (fmt(pm[1], 2), fmt(pm[2], 2)), transform=ax.transAxes,
                    fontsize=7, fontweight="bold" if sig else "normal", color="#B00000" if sig else "#333")
            out = C._nmds_out.get((t, mx))
            if out:
                ax.text(0.01, 0.045, "not shown (composition unlike all others): %s" % R.short(", ".join(out), 60),
                        transform=ax.transAxes, fontsize=6.5, color="#B00000")
            ax.set_xlabel("NMDS1")
            ax.set_ylabel("NMDS2")
            ax.set_aspect("equal", adjustable="datalim")
            ax.legend(frameon=False, fontsize=6, loc="upper right")
        pdf.savefig(fig)
        plt.close(fig)
        done = True
    return done


def page_beta(pdf, C, mx="genus"):
    res = [(t, C.beta(t, mx)) for t in C.types]
    res = [(t, r) for t, r in res if r]
    if not res:
        return False
    fig = cpage(pdf, "beta diversity (%s)" % mx, "Bray-Curtis PCoA per fraction; PERMANOVA marginal R²", "beta")
    for i, (t, (ks, coords, var, perm)) in enumerate(res[:2]):
        ax = fig.add_axes([0.06 + i * 0.48, 0.36, 0.40, 0.48])
        for g in C.groups + ["NA"]:
            idx = [j for j, k in enumerate(ks) if C.group(k[0]) == g]
            if idx:
                ax.scatter(coords[idx, 0], coords[idx, 1], s=12, color=C.gcolor[g], label="%s (%d)" % (g, len(idx)),
                           alpha=0.8, lw=0)
        ax.set_xlabel("PCo1 (%.0f %%)" % (100 * var[0]))
        ax.set_ylabel("PCo2 (%.0f %%)" % (100 * var[1]) if len(var) > 1 else "PCo2")
        ax.set_title(t, loc="left")
        ax.legend(frameon=False, fontsize=6.5)
        tab = fig.add_axes([0.06 + i * 0.48, 0.04, 0.40, 0.24])
        R.table(tab, [[r[0], r[1], r[2], fmt(r[3], 3), fmt(r[4], 3), fmt(r[5], 3)] for r in perm],
                ["variable", "n", "levels", "pseudo-F", "R²", "p"], fontsize=6.5,
                title="PERMANOVA %s (bold: p < 0.05)" % t)
        bold_rows(tab, [isnum(r[5]) and r[5] < 0.05 for r in perm])
    pdf.savefig(fig)
    plt.close(fig)


def poolability(C, t):
    """(text, block effect detected?) for the reference group across blocks"""
    be = C.block_effect(t)
    ref = C.a.reference or C.groups[0]
    col = C.a.match_col
    if not be:
        return "too few %s libraries per %s level to check the %s effect" % (ref, col, col), None
    det = (isnum(be["p"]) and be["p"] < 0.05) or be["metrics_q"] > 0
    if str(be.get("scope", "")).startswith("all"):
        how = ("%s after group" % col if be["scope"] == "all" else "%s alone, one library per group and level" % col)
        txt = ("All libraries across %d %s levels (%d libraries; %s; permuted within groups): "
               "%sgenus R² %s, p %s. " % (be["levels"], col, be["libraries"], how,
                                          "partial " if be["scope"] == "all" else "", fmt(be["R2"], 2),
                                          fmt(be["p"], 2)))
        if not isnum(be["p"]):
            return txt + "Not testable with this design.", None
        txt += ("Members of the same %s are more alike than expected: pooled p-values are too optimistic; use the "
                "within-%s comparisons and the models with %s-clustered SE." % (col, col, col) if det else
                "No clear %s effect: pooled and within-%s results should agree." % (col, col))
        return txt, det
    txt = ("%s across %d %s levels (%d libraries): genus PERMANOVA R² %s, p %s; %d of %d metrics differ (q < %g). " % (
        ref, be["levels"], col, be["libraries"], fmt(be["R2"], 2), fmt(be["p"], 2), be["metrics_q"],
        be["metrics_tested"], SIG_Q))
    txt += ("A %s effect is detectable: prefer the within-%s comparisons and include %s in the models." % (col, col, col)
            if det else "No %s effect detectable within %s: pooled comparisons are reasonable (absence of evidence — "
                        "small groups have little power)." % (col, ref))
    return txt, det


def page_matched(pdf, C):
    if not C.a.match_col:
        return False
    done = False
    col = C.a.match_col
    for t in C.types:
        sets = C.matched_sets(t)
        des = C.block_design(t)
        if not sets:
            continue
        fig = cpage(pdf, "block design and comparisons (%s) · %s" % (col, t),
                    "block = %s; pooled: all libraries of a group vs all %s; within: only %s levels containing both, "
                    "tested within levels (van Elteren / Mann-Whitney; PERMANOVA permuted within levels)" % (
                        col, C.a.reference or C.groups[0], col), "matched")
        rows = [[d["group"], d["blocks"], len(d["shared"]), d["design"], "%d / %d" % (d["n_a"], d["n_b"]),
                 "%d / %d" % (d["n_a_all"], d["n_b_all"])] for d in des]
        ax = fig.add_axes([0.03, 0.70, 0.50, 0.14])
        R.table(ax, rows, ["group", "%s levels" % col, "shared with %s" % (C.a.reference or C.groups[0]), "design",
                           "within: a / b", "pooled: a / b"], fontsize=6.5, rowscale=1.25, left=(0, 3),
                title="design (nested: few shared levels, e.g. site · crossed: many, e.g. families)")
        txt, det = poolability(C, t)
        import textwrap
        fig.text(0.56, 0.83, "poolability of the reference group", fontsize=8, fontweight="bold", va="top")
        fig.text(0.56, 0.805, "\n".join(textwrap.wrap(txt, 78)), fontsize=7, va="top", linespacing=1.35,
                 color="#B00000" if det else "#1a6e1a" if det is False else "#555")
        rows, flags = [], []
        for label, ka, kb, kind in sets:
            r = C.two_set_tests(t, ka, kb, blocked=(kind == "within"))
            pm = r["permanova"]
            nsig = sum(1 for x in r["metrics"] if not math.isnan(x[7]) and x[7] < SIG_Q)
            gsig = sum(1 for x in r["genus"] if not math.isnan(x[6]) and x[6] < SIG_Q)
            nbl = len({C.block_of(k) for k in ka} & {C.block_of(k) for k in kb})
            test = "van Elteren" if kind == "within" and nbl > 1 else "Mann-Whitney"
            rows.append([label, kind, test, len(ka), len(kb), fmt(pm[3], 3) if pm else "–",
                         fmt(pm[4], 3) if pm else "–", nsig, gsig])
            flags.append((bool(pm) and isnum(pm[4]) and pm[4] < 0.05) or nsig > 0 or gsig > 0)
        ax = fig.add_axes([0.03, 0.40, 0.94, 0.24])
        R.table(ax, rows, ["comparison", "kind", "test", "n a", "n b", "genus R²", "PERMANOVA p",
                           "metrics q<%g" % SIG_Q, "genera q<%g" % SIG_Q],
                col_widths=[0.36, 0.07, 0.10, 0.05, 0.05, 0.08, 0.10, 0.095, 0.095], left=(0, 1, 2), fontsize=6.5,
                rowscale=1.2, title="comparisons (BH-FDR within each comparison; bold: PERMANOVA p < 0.05 or any q < %g)"
                % SIG_Q)
        bold_rows(ax, flags)
        within = [x for x in sets if x[3] in ("within", "block")][:3]
        for i, (label, ka, kb, kind) in enumerate(within):
            r = C.two_set_tests(t, ka, kb, blocked=(kind == "within"))
            ax = fig.add_axes([0.03 + i * 0.325, 0.04, 0.30, 0.27])
            R.table(ax, [[R.short(x[0].replace("integrity_", "").replace("nonhuman_", "nh_"), 26), fmt(x[3], 3),
                          fmt(x[4], 3), fmt(x[7], 2)] for x in r["metrics"][:7]],
                    ["metric", "median a", "median b", "q"], fontsize=5.8, rowscale=1.1,
                    title=R.short(label, 52), col_widths=[0.49, 0.18, 0.18, 0.15])
            bold_rows(ax, [isnum(x[7]) and x[7] < SIG_Q for x in r["metrics"][:7]])
        pdf.savefig(fig)
        plt.close(fig)
        done = True
    return done


def page_rna_classes(pdf, C):
    """composition of every library by RNA class (within samples) and group tests per class / Rfam family"""
    done = False
    for t in C.types:
        ks = sorted(C.keys(t), key=lambda k: (C.groups.index(C.group(k[0])) if C.group(k[0]) in C.groups else 99, k[0]))
        ks = [k for k in ks if sum(C.vec[k]["rna_class"].values()) > 0]
        if not ks:
            continue
        fig = cpage(pdf, "RNA classes · %s" % t, "share of non-human reads per RNA class in every library (sorted by "
                    "group) · group tests per class and per Rfam family (non-rRNA/tRNA); bold: q < %g" % SIG_Q, "classes")
        major = ("rRNA", "tRNA", "unknown")
        for pi, (classes, ylab, y0, hh) in enumerate((
                (RNA_CLASSES, "% of non-human reads", 0.50, 0.34),
                ([c for c in RNA_CLASSES if c not in major], "% of the other classes", 0.08, 0.34))):
            ax = fig.add_axes([0.05, y0, 0.68, hh])
            tot = [sum(C.vec[k]["rna_class"].get(c, 0.0) for c in classes) for k in ks]
            bottom = np.zeros(len(ks))
            for c in classes:
                v = np.array([C.vec[k]["rna_class"].get(c, 0.0) / (tt or 1) for k, tt in zip(ks, tot)])
                if v.sum() == 0:
                    continue
                ax.bar(range(len(ks)), 100 * v, bottom=100 * bottom, color=CLASS_COLORS[c], width=0.9,
                       label=c if pi == 0 or c not in major else None)
                bottom += v
            ax.set_xlim(-0.5, len(ks) - 0.5)
            ax.set_ylim(0, 100)
            ax.set_ylabel(ylab, fontsize=7)
            ax.set_xticks(range(len(ks)))
            ax.set_xticklabels([k[0] for k in ks] if pi == 1 else [], rotation=90, fontsize=4.5 if len(ks) > 60 else 5.5)
            pos = 0
            for g in C.groups + ["NA"]:
                n = sum(1 for k in ks if C.group(k[0]) == g)
                if n:
                    if pi == 0:
                        ax.text(pos + n / 2 - 0.5, 101, "%s (%d)" % (g, n), ha="center", va="bottom", fontsize=7,
                                color=C.gcolor.get(g, "#333"), fontweight="bold")
                    if pos:
                        ax.axvline(pos - 0.5, color="black", lw=0.8)
                    pos += n
            if pi == 0:
                ax.legend(frameon=False, fontsize=6, loc="upper left", bbox_to_anchor=(1.01, 1.0), ncol=1)
            else:
                ax.set_title("without rRNA, tRNA and unknown", loc="left", fontsize=7)
        pdf.savefig(fig)
        plt.close(fig)
        fig = R.page(pdf, COHORT_NAME, "RNA classes and Rfam families between groups · %s" % t,
                     "Kruskal-Wallis across %s, BH-FDR within each table; bold: q < %g; Int vs Ext of the same samples: "
                     "paired page (class_* metrics)" % (", ".join(C.groups), SIG_Q))
        # class tests (metrics class_*) and Rfam families
        keys = [k_ for k_ in C.metrics if k_.startswith("class_")]
        rows, flags = [], []
        for k_ in keys:
            p, q = C.metric_pq(t, k_)
            med = [np.median([C.m[x][k_] for x in C.keys(t) if C.group(x[0]) == g and isnum(C.m[x].get(k_))] or
                             [float("nan")]) for g in C.groups]
            if all(not isnum(x) or x == 0 for x in med):
                continue
            rows.append([k_[6:].replace("_", " "), fmt(p, 2) + stars(p), fmt(q, 2)] + [fmt(100 * x, 3) for x in med])
            flags.append(isnum(q) and q < SIG_Q)
        axc = fig.add_axes([0.03, 0.05, 0.46, 0.80])
        R.table(axc, rows, ["class", "p", "q"] + ["% " + g for g in C.groups], fontsize=6.5, left=(0,), rowscale=1.2,
                title="RNA classes between groups (median % of non-human reads)")
        bold_rows(axc, flags)
        ft = C.feature_tests(t, "rfam_family") if len(C.groups) > 1 else []
        axf = fig.add_axes([0.52, 0.05, 0.45, 0.80])
        R.table(axf, [[R.short(r[0], 26), fmt(r[3], 2) + stars(r[3]), fmt(r[4], 2)] + [fmt(100 * x, 3) for x in r[5:]]
                      for r in ft[:25]], ["Rfam family", "p", "q"] + ["% " + g for g in C.groups], fontsize=6.5,
                left=(0,), rowscale=1.2, title="Rfam families (%d tested; mean %% of other Rfam reads)" % len(ft))
        bold_rows(axf, [isnum(r[4]) and r[4] < SIG_Q for r in ft[:25]])
        pdf.savefig(fig)
        plt.close(fig)
        done = True
    return done


def page_functions(pdf, C, top=15):
    """functional categories of the sRNA reads on CDS: top KEGG pathways per group + tests per database"""
    done = False
    for t in C.types:
        ks, feats, Cm, Rl = C.rel_matrix(t, "fn_kegg_pathway")
        if len(ks) < 4 or not feats:
            continue
        fig = cpage(pdf, "functional categories · %s" % t, "sRNA reads on CDS (sense) per KEGG pathway / module / KO "
                    "and Pfam domain; share within each database; Kruskal-Wallis, BH-FDR; bold: q < %g" % SIG_Q,
                    "functions")
        labels = [C.group(k[0]) for k in ks]
        order = list(np.argsort(-Rl.mean(axis=0))[:top])
        ax = fig.add_axes([0.30, 0.42, 0.30, 0.44])
        h = 0.8 / max(len(C.groups), 1)
        for gi, g in enumerate(C.groups):
            idx = [i for i, l in enumerate(labels) if l == g]
            if not idx:
                continue
            ax.barh([y + (gi - (len(C.groups) - 1) / 2) * h for y in range(len(order))],
                    [100 * Rl[idx, j].mean() for j in order], height=h, color=C.gcolor[g], label="%s (%d)" % (g, len(idx)))
        tests = {r[0]: r for r in C.feature_tests(t, "fn_kegg_pathway")}
        ax.set_yticks(range(len(order)))
        ax.set_yticklabels(["%s%s" % (tname(feats[j], 38), stars(tests.get(feats[j], [0, 0, 0, float("nan")])[3]))
                            for j in order], fontsize=6.3)
        for lab_, j in zip(ax.get_yticklabels(), order):
            tr = tests.get(feats[j])
            if tr and isnum(tr[4]) and tr[4] < SIG_Q:
                lab_.set_fontweight("bold")
                lab_.set_color("#B00000")
        ax.invert_yaxis()
        ax.set_xlabel("% of pathway-assigned reads", fontsize=7)
        ax.set_title("top KEGG pathways (mean per group; * p<0.05)", loc="left", fontsize=8)
        ax.legend(frameon=False, fontsize=6.5, loc="lower right")
        for i, (mx, lab) in enumerate((("fn_kegg_module", "KEGG modules"), ("fn_kegg_ko", "KOs"), ("fn_pfam", "Pfam"),
                                       ("fn_cog", "COGs"))):
            ft = C.feature_tests(t, mx)
            nq = sum(1 for r in ft if isnum(r[4]) and r[4] < SIG_Q)
            axt = fig.add_axes([0.64 + (i % 2) * 0.18, 0.64 - (i // 2) * 0.27, 0.16, 0.19])
            rows = [[tname(r[0], 24), fmt(r[3], 2), fmt(r[4], 2)] for r in ft[:6]]
            R.table(axt, rows, [lab, "p", "q"], fontsize=5.8, left=(0,), col_widths=[0.62, 0.19, 0.19], rowscale=1.05,
                    title="%s\n%d tested, %d q<%g" % (lab, len(ft), nq, SIG_Q))
            bold_rows(axt, [isnum(r[4]) and r[4] < SIG_Q for r in ft[:6]])
        # matched comparisons on pathways
        rows, flags = [], []
        for label, ka, kb, kind in C.matched_sets(t):
            r = C.two_set_tests(t, ka, kb, blocked=(kind == "within"))
            cnt = [sum(1 for x in r.get(mx, []) if isnum(x[6]) and x[6] < SIG_Q) for mx in
                   ("fn_kegg_pathway", "fn_kegg_module", "fn_kegg_ko", "fn_pfam")]
            best = r.get("fn_kegg_pathway", [])[:1]
            rows.append([label, len(ka), len(kb)] + cnt + [tname(best[0][0], 40) + " (q %s)" % fmt(best[0][6], 2)
                                                         if best else "–"])
            flags.append(any(cnt))
        if rows:
            axm = fig.add_axes([0.03, 0.06, 0.94, 0.26])
            R.table(axm, rows, ["matched comparison", "n a", "n b", "pathways q<%g" % SIG_Q, "modules", "KOs", "Pfam",
                                "top pathway"], fontsize=6.3, left=(0, 7), rowscale=1.25,
                    col_widths=[0.30, 0.05, 0.05, 0.10, 0.07, 0.06, 0.06, 0.31],
                    title="block comparisons of functional profiles")
            bold_rows(axm, flags)
        pdf.savefig(fig)
        plt.close(fig)
        done = True
    return done


def page_antisense(pdf, C):
    done = False
    for t in C.types:
        if "antisense_CDS_fraction" not in C.metrics or not any(isnum(C.m[k].get("antisense_CDS_fraction"))
                                                                 for k in C.keys(t)):
            continue
        nfilt = sum(1 for k in C.keys(t) if C.m[k].get("antisense_rRNA_filtered") == 1)
        fig = cpage(pdf, "antisense transcription · %s" % t, "sRNA reads antisense to protein-coding genes; antisense "
                    "index = AS / (sense + AS) per gene function (≥ 10 reads); bold: q < %g · rRNA reads removed in %d of "
                    "%d libraries" % (SIG_Q, nfilt, len(C.keys(t))), "antisense")
        if nfilt < len(C.keys(t)):
            fig.text(0.03, 0.505, "WARNING: %d libraries include rRNA reads in the antisense counts (reads sense to "
                     "unannotated rRNA look antisense to spurious ORFs) — rerun: srna_pipeline.sh -m assembly -f" % (
                         len(C.keys(t)) - nfilt), fontsize=6.8, color="#B00000", fontweight="bold")
        for i, (mk, lab, kw) in enumerate((("antisense_CDS_fraction", "antisense on CDS (%)", dict(pct=True)),
                                           ("antisense_strand_error", "strand error (%)", dict(pct=True)),
                                           ("antisense_CDS_enrichment_log2", "log2 CDS AS / error", {}))):
            ax = fig.add_axes([0.05 + i * 0.165, 0.60, 0.13, 0.25])
            boxes(ax, C, mk, lab, types=[t], **kw)
        # hotspots: KEGG pathways and KOs, Pfam
        for i, (db, lab) in enumerate((("kegg_pathway", "KEGG pathways"), ("kegg_ko", "KOs"), ("pfam", "Pfam"))):
            hs = C.antisense_hotspots(t, db)
            up = sorted([r for r in hs if r[5] > 0], key=lambda r: (not (isnum(r[8]) and r[8] < SIG_Q), -r[5]))[:8]
            ax = fig.add_axes([0.03 + i * 0.32, 0.05, 0.30, 0.42])
            nsig = sum(1 for r in hs if isnum(r[8]) and r[8] < SIG_Q and r[5] > 0)
            R.table(ax, [[tname(r[0], 30), fmt(r[3], 2), fmt(r[4], 2), fmt(r[6], 2), fmt(r[8], 2)] for r in up],
                    [lab, "AS idx", "library", "> lib", "q"], fontsize=5.8, left=(0,), rowscale=1.15,
                    col_widths=[0.48, 0.13, 0.13, 0.12, 0.14],
                    title="%s more antisense than the library\n(%d tested, %d q<%g)" % (lab, len(hs), nsig, SIG_Q))
            bold_rows(ax, [isnum(r[8]) and r[8] < SIG_Q for r in up])
        # group / matched tests on the antisense index
        rows, flags = [], []
        if len(C.groups) > 1:
            for db, lab in (("kegg_pathway", "pathways"), ("kegg_ko", "KOs"), ("pfam", "Pfam")):
                at = C.asidx_tests(t, db)
                nq = sum(1 for r in at if isnum(r[-1]) and r[-1] < SIG_Q)
                best = at[0] if at else None
                rows.append(["groups", lab, len(at), nq, (tname(best[0], 26) + " p=%s" % fmt(best[4], 2)) if best else "–"])
                flags.append(nq > 0)
            for label, ka, kb, kind in C.matched_sets(t):
                at = C.asidx_tests(t, "kegg_ko", ka, kb)
                nq = sum(1 for r in at if isnum(r[-1]) and r[-1] < SIG_Q)
                rows.append([label.split(" (")[0], "KOs", len(at), nq,
                             (tname(at[0][0], 26) + " p=%s" % fmt(at[0][6], 2)) if at else "–"])
                flags.append(nq > 0)
        if rows:
            ax = fig.add_axes([0.55, 0.58, 0.42, 0.27])
            R.table(ax, rows, ["comparison", "level", "tested", "q<%g" % SIG_Q, "top term"], fontsize=6, left=(0, 1, 4),
                    rowscale=1.2, col_widths=[0.28, 0.12, 0.10, 0.10, 0.40],
                    title="antisense index between groups / block comparisons")
            bold_rows(ax, flags)
        pdf.savefig(fig)
        plt.close(fig)
        done = True
    return done


MODEL_SETS = ("metrics", "genus", "family", "mirna", "rfam_family", "rna_class", "fn_kegg_pathway", "fn_kegg_module",
              "fn_kegg_ko", "fn_pfam", "fn_cog", "asidx_kegg_pathway", "asidx_kegg_ko")
MODEL_LABELS = {"metrics": "library metrics", "genus": "genera (CLR)", "family": "families (CLR)",
                "mirna": "miRNAs (CLR)", "rfam_family": "Rfam families (CLR)", "rna_class": "RNA classes (CLR)",
                "fn_kegg_pathway": "KEGG pathways (CLR)", "fn_kegg_module": "KEGG modules (CLR)",
                "fn_kegg_ko": "KOs (CLR)", "fn_pfam": "Pfam (CLR)", "fn_cog": "COGs (CLR)",
                "asidx_kegg_pathway": "antisense index: pathways", "asidx_kegg_ko": "antisense index: KOs"}



OVERVIEW_SETS = (("genus", "genera"), ("family", "families"), ("rna_class", "RNA cl."),
                 ("fn_kegg_module", "KEGG mod."), ("fn_kegg_ko", "KOs"), ("fn_pfam", "Pfam"),
                 ("metrics", "metrics"), ("asidx_kegg_ko", "AS-idx KO"))


def page_models_overview(pdf, C):
    """one table over all models: for the term of interest of every model (group contrasts, or the first term)
    the number of features with q < SIG_Q per outcome set, and the strongest results"""
    if not C.models():
        return False
    rows, flags = [], []
    for mi in range(len(C.models())):
        for t in C.types:
            md = C.model_data(t, mi)
            if not md:
                continue
            for gt in C.group_terms(t, mi):
                cells, best = [], []
                for what, _ in OVERVIEW_SETS:
                    res = [r for r in C.model_fit(t, what, mi) if r[3] == gt]
                    sig = [r for r in res if isnum(r[8]) and r[8] < SIG_Q]
                    cells.append("%d / %d" % (len(sig), len(res)) if res else "–")
                    best += [(r[8], what, r) for r in sig]
                best.sort(key=lambda x: x[0])
                top = "; ".join("%s %s (q %s)" % (R.short(tname(r[0], 18), 20), "↑" if r[4] > 0 else "↓",
                                                  fmt(q, 2)) for q, w, r in best[:2])
                subset = C.model_subset(mi)
                rows.append(["%s%s" % (C.model_label(mi), " [%s]" % "/".join(sorted(subset[1])) if subset else ""),
                             t.replace("sRNA-", ""), R.short(gt.replace("%s: " % C.group_col, ""), 34),
                             str(len(md[0]))] + cells + [top or "–"])
                flags.append(bool(best))
    if not rows:
        return False
    per = 30
    for start in range(0, len(rows), per):
        fig = cpage(pdf, "adjusted models — overview%s" % ("" if start == 0 else " (continued)"),
                    "term of interest per model (the first term of the formula) · features with q < %g / "
                    "tested per outcome set · ↑/↓: higher / lower with the term (adjusted)" % SIG_Q, "models_overview")
        ax = fig.add_axes([0.03, 0.04, 0.92, 0.80])
        R.table(ax, rows[start:start + per], ["model", "fraction", "term", "n"] + [lab for _, lab in OVERVIEW_SETS] +
                ["strongest results"], fontsize=5.8, left=(0, 1, 2, 12), rowscale=1.1,
                col_widths=[0.05, 0.04, 0.12, 0.03] + [0.052] * len(OVERVIEW_SETS) + [0.344])
        bold_rows(ax, flags[start:start + per])
        pdf.savefig(fig)
        plt.close(fig)
    return True


def page_models(pdf, C):
    done = False
    for mi, formula in enumerate(C.models()):
      for t in C.types:
        md = C.model_data(t, mi)
        if not md:
            continue
        keys, X, names, owner, info = md[:5]
        gts = C.group_terms(t, mi)
        mlab = C.model_label(mi)
        se_txt = {"": "HC3 robust SE", "cluster": "SE clustered by %s" % C.a.match_col,
                  "within": "within-%s fixed effects, SE clustered by %s" % (C.a.match_col, C.a.match_col)}[C.model_variant(mi)]
        fig = cpage(pdf, "adjusted model %s · %s" % (mlab, t), "%s: %s · %d libraries (complete cases) · OLS on CLR "
                    "abundances / metrics, %s, BH-FDR per term; bold: q < %g" % (
                        mlab, formula, len(keys), se_txt, SIG_Q), "models")
        if info:
            fig.text(0.03, 0.855, "; ".join(info), fontsize=7, color="#B00000")
        rows, flags = [], []
        for what in MODEL_SETS:
            res = C.model_fit(t, what, mi)
            if not res:
                continue
            for gt in gts:
                sub = [r for r in res if r[3] == gt]
                sig = [r for r in sub if isnum(r[8]) and r[8] < SIG_Q]
                top = sorted(sub, key=lambda r: (not isnum(r[7]), r[7]))[:2]
                rows.append([MODEL_LABELS.get(what, what), gt.replace("%s: " % C.group_col, ""), len(sub), len(sig),
                             "; ".join("%s %s%s (q %s)" % (R.short(tname(r[0], 26), 30), "+" if r[4] > 0 else "−",
                                                           fmt(abs(r[4]), 2), fmt(r[8], 2)) for r in top)])
                flags.append(bool(sig))
        ax = fig.add_axes([0.02, 0.04, 0.96, 0.79])
        R.table(ax, rows[:34], ["outcomes", "contrast", "tested", "q<%g" % SIG_Q, "strongest (estimate: CLR units / "
                                "metric units, adjusted)"], fontsize=6, left=(0, 1, 4), rowscale=1.1,
                col_widths=[0.19, 0.11, 0.06, 0.06, 0.58])
        bold_rows(ax, flags)
        pdf.savefig(fig)
        plt.close(fig)
        # volcano plots: group contrasts for genera, KEGG pathways and the KO antisense index
        fig = R.page(pdf, COHORT_NAME, "adjusted model %s · volcano plots · %s" % (mlab, t), "%s: %s · red: q < %g; "
                     "labels: strongest significant features" % (mlab, formula, SIG_Q))
        for vi, (what, lab, xl) in enumerate((("genus", "genera", "adjusted CLR difference"),
                                              ("fn_kegg_pathway", "KEGG pathways", "adjusted CLR difference"),
                                              ("asidx_kegg_ko", "antisense index KOs", "adjusted index difference"))):
          res = C.model_fit(t, what, mi)
          for i, gt in enumerate(gts[:3]):
            sub = [r for r in res if r[3] == gt and isnum(r[7])]
            if not sub:
                continue
            axv = fig.add_axes([0.06 + i * 0.32, 0.66 - vi * 0.30, 0.25, 0.21])
            est = np.array([r[4] for r in sub])
            lp = -np.log10(np.maximum([r[7] for r in sub], 1e-300))
            sig = np.array([isnum(r[8]) and r[8] < SIG_Q for r in sub])
            axv.scatter(est[~sig], lp[~sig], s=5, color="#999", lw=0)
            axv.scatter(est[sig], lp[sig], s=9, color="#B00000", lw=0)
            for r, e, l, sg in sorted(zip(sub, est, lp, sig), key=lambda z: -z[2])[:4]:
                if sg:
                    axv.text(e, l, " " + R.short(tname(r[0], 22), 22), fontsize=5.5)
            axv.axvline(0, color="#ccc", lw=0.8)
            axv.set_xlabel(xl, fontsize=6.5)
            axv.set_ylabel("−log10 p", fontsize=6.5)
            axv.set_title("%s · %s (%d sig.)" % (lab, gt.replace("%s: " % C.group_col, ""), int(sig.sum())),
                          loc="left", fontsize=7.5)
        pdf.savefig(fig)
        plt.close(fig)
        done = True
    return done


def page_paired(pdf, C):
    rows = [r for r in C.paired() if not math.isnan(r[5])]
    if not rows:
        return False
    rows.sort(key=lambda r: r[5])
    fig = cpage(pdf, "sRNA-Int vs sRNA-Ext (paired)", "samples with both fractions; Wilcoxon signed-rank, BH-FDR",
                "paired")
    ax = fig.add_axes([0.03, 0.04, 0.94, 0.82])
    R.table(ax, [[r[0], r[1], fmt(r[2], 3), fmt(r[3], 3), fmt(r[4], 3), fmt(r[5], 2), fmt(r[6], 2)] for r in rows[:30]],
            ["metric", "pairs", "median Int", "median Ext", "median (Int−Ext)", "p", "q"],
            col_widths=[0.40, 0.07, 0.11, 0.11, 0.13, 0.08, 0.08], fontsize=6.5, rowscale=1.15,
            title="bold: q < %g" % SIG_Q)
    bold_rows(ax, [isnum(r[6]) and r[6] < SIG_Q for r in rows[:30]])
    pdf.savefig(fig)
    plt.close(fig)


def page_features(pdf, C):
    if len(C.groups) < 2:
        return False
    fig = cpage(pdf, "group differences (exploratory)", "Kruskal-Wallis across %s; BH-FDR within each table" %
                ", ".join(C.groups), "features")
    blocks = []
    for t in C.types:
        mt, _ = C.metric_tests(t)
        mt = sorted([r for r in mt if not math.isnan(r[3])], key=lambda r: r[3])
        blocks.append(("metrics · %s" % t, [[r[0], r[1], fmt(r[3], 2), fmt(r[4], 2)] + [fmt(x, 3) for x in r[5:]]
                                            for r in mt[:8]], "median", [r[4] for r in mt[:8]]))
        for mx in ("genus", "mirna"):
            if mx == "mirna" and t not in C.mirna_types():
                continue
            ft = C.feature_tests(t, mx)
            blocks.append(("%s · %s" % (mx, t), [[r[0], fmt(r[1], 2), fmt(r[3], 2), fmt(r[4], 2)] +
                                                 [fmt(100 * x, 3) for x in r[5:]] for r in ft[:8]], "mean %",
                           [r[4] for r in ft[:8]]))
    for i, (title, rows, what, qs) in enumerate(blocks[:6]):
        r, c = divmod(i, 2)
        ax = fig.add_axes([0.03 + c * 0.49, 0.60 - r * 0.27, 0.46, 0.22])
        ng = len(C.groups)
        R.table(ax, [[R.short(x[0].replace("integrity_", "").replace("nonhuman_", "nh_"), 36)] + x[1:] for x in rows],
                ["feature", "n/prev", "p", "q"] + [R.short(g, 10) for g in C.groups], fontsize=6,
                col_widths=[0.40, 0.09, 0.08, 0.08] + [0.35 / max(ng, 1)] * ng,
                title="%s (%s per group)" % (title, what), rowscale=1.1)
        bold_rows(ax, [isnum(q) and q < SIG_Q for q in qs])
    pdf.savefig(fig)
    plt.close(fig)


def _iqr(v):
    v = [x for x in v if isnum(x)]
    if not v:
        return "–"
    return "%s (%s–%s)" % (fmt(float(np.median(v)), 3), fmt(float(np.percentile(v, 25)), 3),
                           fmt(float(np.percentile(v, 75)), 3))


INTEGRITY_LABELS = [
    ("qc_lRNA_RIN", "wet-lab RIN (long RNA)", "+"),
    ("integrity_tRNA_mature_end_fraction", "tRNA reads with a mature end", "+"),
    ("integrity_tRNA_internal_fraction", "tRNA internal fragments", "-"),
    ("integrity_tRNA_3p_CCA_fraction", "3′-CCA share", "+"),
    ("integrity_tRNA_5p_end_entropy", "tRNA 5′-end entropy", "-"),
    ("integrity_tRNA_5p_top3_share", "tRNA 5′ ends at main positions", "+"),
    ("integrity_tRNA_3p_top3_share", "tRNA 3′ ends at main positions", "+"),
    ("integrity_rRNA_5p_end_entropy", "rRNA 5′-end entropy", "-"),
    ("integrity_rRNA_bin_evenness", "rRNA spread along 16S/23S", "-"),
    ("integrity_tRNA_halves_fraction", "tRNA halves (stress)", ""),
    ("integrity_reads_le20nt_fraction", "reads ≤ 20 nt", ""),
    ("integrity_median_read_length", "median read length", ""),
]


def page_integrity_summary(pdf, C):
    keys = [(k, lab, d) for k, lab, d in INTEGRITY_LABELS if k in C.metrics]
    if not keys:
        return False
    fig = cpage(pdf, "RNA integrity — summary", "medians (IQR) per fraction; paired Int vs Ext; group test per "
                "fraction; Spearman ρ with RIN; bold: q < %g" % SIG_Q, "integrity_sum")
    paired = {r[0]: r for r in C.paired()}
    ivq = {}
    for r in C.integrity_vs_qc():
        if "rin" in r[1].lower():
            ivq[(r[0], "integrity_" + r[2])] = r
    rows, flags, notes = [], [], []
    for k, lab, d in keys:
        row = [lab]
        sig = False
        for t in C.types:
            row.append(_iqr([C.m[x].get(k, float("nan")) for x in C.keys(t)]))
        pr = paired.get(k)
        if pr and isnum(pr[5]):
            row.append("%s%s" % (fmt(pr[5], 2), stars(pr[5])))
            if isnum(pr[6]) and pr[6] < SIG_Q:
                sig = True
                notes.append("%s differs between Int and Ext (median %s vs %s, q=%s)." % (
                    lab, fmt(pr[2], 3), fmt(pr[3], 3), fmt(pr[6], 2)))
        else:
            row.append("–")
        for t in C.types:
            p, q = C.metric_pq(t, k)
            row.append("%s%s" % (fmt(p, 2), stars(p)) if isnum(p) else "–")
            if isnum(q) and q < SIG_Q:
                sig = True
                notes.append("%s differs between groups in %s (p=%s, q=%s)." % (lab, t, fmt(p, 2), fmt(q, 2)))
        for t in C.types:
            r = ivq.get((t, k))
            if r and isnum(r[4]):
                row.append("%s%s" % (fmt(r[4], 2), stars(r[5])))
                if isnum(r[6]) and r[6] < SIG_Q:
                    sig = True
                    notes.append("%s correlates with RIN in %s (ρ=%s, q=%s)." % (lab, t, fmt(r[4], 2), fmt(r[6], 2)))
                elif isnum(r[5]) and r[5] < 0.05:
                    notes.append("%s tends to correlate with RIN in %s (ρ=%s, p=%s, not FDR-significant)." % (
                        lab, t, fmt(r[4], 2), fmt(r[5], 2)))
            else:
                row.append("–")
        rows.append(row)
        flags.append(sig)
    ax = fig.add_axes([0.02, 0.36, 0.96, 0.50])
    hdr = ["indicator"] + ["median %s" % t.replace("sRNA-", "") for t in C.types] + ["Int vs Ext p"] + \
          ["groups p %s" % t.replace("sRNA-", "") for t in C.types] + ["ρ RIN %s" % t.replace("sRNA-", "") for t in C.types]
    nt = len(C.types)
    R.table(ax, rows, hdr, fontsize=6.4, left=(0,), rowscale=1.25,
            col_widths=[0.20] + [0.15] * nt + [0.08] + [0.08] * nt + [0.08] * nt)
    bold_rows(ax, flags)
    txt = notes[:14] or ["No integrity indicator differs significantly (q < %g) between fractions or groups, and none "
                         "is FDR-significantly correlated with RIN." % SIG_Q]
    import textwrap
    fig.text(0.03, 0.31, "\n".join("• " + "\n   ".join(textwrap.wrap(n, 165)) for n in txt), fontsize=6.8,
             va="top", linespacing=1.35)
    pdf.savefig(fig)
    plt.close(fig)


def page_mirna_top(pdf, C, top=12):
    done = False
    for t in C.mirna_types():
        ks, feats, Cm, Rl = C.rel_matrix(t, "mirna")
        nall = len(C.keys(t))
        if not feats:
            continue
        tests = {r[0]: r for r in C.feature_tests(t, "mirna")}
        order = np.argsort(-Rl.mean(axis=0))[:top]
        fig = cpage(pdf, "top human miRNAs by group · %s" % t,
                    "%% of human miRNA reads; libraries with ≥ %g miRNA reads (%d of %d); Kruskal-Wallis p (stars), "
                    "bold red: q < %g" % (C.a.min_mirna_reads, len(ks), nall, SIG_Q), "mirna_top")
        labels = [C.group(k[0]) for k in ks]
        for i, j in enumerate(order):
            r_, c_ = divmod(i, 4)
            ax = fig.add_axes([0.05 + c_ * 0.24, 0.62 - r_ * 0.28, 0.19, 0.20])
            for gi, g in enumerate(C.groups):
                v = [100 * Rl[x, j] for x in range(len(ks)) if labels[x] == g]
                if not v:
                    continue
                col = C.gcolor[g]
                ax.boxplot([v], positions=[gi], widths=0.6, showfliers=False, patch_artist=True,
                           boxprops=dict(facecolor=col + "40", edgecolor=col), medianprops=dict(color=col),
                           whiskerprops=dict(color=col), capprops=dict(color=col))
                jit = random.Random(int(gi * 7 + int(j)))
                ax.scatter([gi + (jit.random() - 0.5) * 0.3 for _ in v], v, s=6, color=col, alpha=0.6, lw=0, zorder=3)
            ax.set_xticks(range(len(C.groups)))
            ax.set_xticklabels(C.groups, fontsize=7)
            name = feats[j].replace("hsa-", "")
            tr = tests.get(feats[j])
            p, q = (tr[3], tr[4]) if tr else (float("nan"), float("nan"))
            sig = isnum(q) and q < SIG_Q
            ax.set_title("%s  p=%s%s" % (name, fmt(p, 2), stars(p)), loc="left", fontsize=7.5,
                         fontweight="bold" if sig else "normal", color="#B00000" if sig else "black")
            ax.set_ylabel("% miRNA reads", fontsize=6.5)
        pdf.savefig(fig)
        plt.close(fig)
        done = True
    return done


def page_findings(pdf, C):
    """automatic summary: everything below the significance thresholds, per area"""
    items = []      # (area, finding, statistic)
    for t in C.types:
        mt, pw = C.metric_tests(t)
        for r in mt:
            if isnum(r[4]) and r[4] < SIG_Q:
                items.append(("groups · %s" % t, "%s differs between groups (medians %s)" % (
                    r[0], ", ".join("%s %s" % (g, fmt(x, 3)) for g, x in zip(C.groups, r[5:]))),
                    "p=%s q=%s" % (fmt(r[3], 2), fmt(r[4], 2))))
        for mx in ("genus", "mirna", "rfam_family") + FUNC_MX:
            if mx == "mirna" and t not in C.mirna_types():
                continue
            for r in C.feature_tests(t, mx):
                if isnum(r[4]) and r[4] < SIG_Q:
                    items.append(("groups · %s · %s" % (t, mx.replace("fn_", "")), "%s (mean %% %s)" % (
                        tname(r[0], 50), ", ".join("%s %s" % (g, fmt(100 * x, 3)) for g, x in zip(C.groups, r[5:]))),
                        "p=%s q=%s" % (fmt(r[3], 2), fmt(r[4], 2))))
            b = C.beta(t, mx) if mx == "genus" else None
            if b:
                for r in b[3]:
                    if isnum(r[5]) and r[5] < 0.05:
                        items.append(("composition · %s" % t, "genus profiles differ by %s (R² %s)" % (r[0], fmt(r[4], 2)),
                                      "PERMANOVA p=%s" % fmt(r[5], 2)))
        for label, ka, kb, kind in C.matched_sets(t):
            r = C.two_set_tests(t, ka, kb, blocked=(kind == "within"))
            if r["permanova"] and isnum(r["permanova"][4]) and r["permanova"][4] < 0.05:
                items.append(("matched · %s" % t, "%s: genus profiles differ (R² %s)" % (label, fmt(r["permanova"][3], 2)),
                              "PERMANOVA p=%s" % fmt(r["permanova"][4], 2)))
            for x in r["metrics"]:
                if isnum(x[7]) and x[7] < SIG_Q:
                    items.append(("matched · %s" % t, "%s: %s (median %s vs %s)" % (label.split(" (")[0], x[0],
                                  fmt(x[3], 3), fmt(x[4], 3)), "q=%s" % fmt(x[7], 2)))
            for mx, what in (("genus", "genera"), ("fn_kegg_pathway", "KEGG pathways"), ("fn_kegg_module", "KEGG modules"),
                             ("fn_kegg_ko", "KOs"), ("fn_pfam", "Pfam domains"), ("fn_cog", "COGs")):
                ng = sum(1 for x in r.get(mx, []) if isnum(x[6]) and x[6] < SIG_Q)
                if ng:
                    items.append(("matched · %s" % t, "%s: %d %s differ" % (label.split(" (")[0], ng, what),
                                  "q<%g" % SIG_Q))
    for t in C.types:
        for db, lab in (("kegg_pathway", "KEGG pathway"), ("kegg_ko", "KO"), ("pfam", "Pfam")):
            hs = [r for r in C.antisense_hotspots(t, db) if isnum(r[8]) and r[8] < SIG_Q and r[5] > 0]
            if hs:
                items.append(("antisense · %s" % t, "%d %ss more antisense than their library, e.g. %s" % (
                    len(hs), lab, "; ".join(tname(r[0], 30) for r in hs[:3])), "q<%g" % SIG_Q))
            if len(C.groups) > 1:
                for r in C.asidx_tests(t, db):
                    if isnum(r[-1]) and r[-1] < SIG_Q:
                        items.append(("antisense · groups · %s" % t, "%s %s antisense index differs (medians %s)" % (
                            lab, tname(r[0], 40), ", ".join("%s %s" % (g, fmt(x, 2)) for g, x in zip(C.groups, r[5:-1]))),
                            "p=%s q=%s" % (fmt(r[4], 2), fmt(r[-1], 2))))
    for mi, formula in enumerate(C.models()):
      for t in C.types:
        for what in MODEL_SETS:
            for gt in C.group_terms(t, mi):
                sig = [r for r in C.model_fit(t, what, mi) if r[3] == gt and isnum(r[8]) and r[8] < SIG_Q]
                if sig:
                    items.append(("adjusted %s · %s" % (C.model_label(mi), t), "%s, %s: %d significant, e.g. %s" % (
                        MODEL_LABELS.get(what, what), gt.replace("%s: " % C.group_col, ""), len(sig),
                        "; ".join("%s %s%s" % (tname(r[0], 30), "+" if r[4] > 0 else "−", fmt(abs(r[4]), 2))
                                  for r in sig[:3])), "q<%g (%s)" % (SIG_Q, C.model_label(mi))))
    pr = [r for r in C.paired() if isnum(r[6]) and r[6] < SIG_Q]
    if pr:
        top = sorted(pr, key=lambda r: r[5])[:6]
        items.append(("Int vs Ext", "%d metrics differ between the fractions, e.g. %s" % (
            len(pr), "; ".join("%s %s vs %s" % (r[0], fmt(r[2], 3), fmt(r[3], 3)) for r in top[:4])), "q<%g" % SIG_Q))
    for r in C.integrity_vs_qc():
        if isnum(r[6]) and r[6] < SIG_Q:
            items.append(("integrity vs QC", "%s: %s ~ %s (ρ %s)" % (r[0], r[2], r[1], fmt(r[4], 2)),
                          "p=%s q=%s" % (fmt(r[5], 2), fmt(r[6], 2))))
    for r in C.covariates():
        if isnum(r[3]) and r[3] < 0.05:
            items.append(("confounding", "%s differs between groups — include as covariate" % r[0], "p=%s" % fmt(r[3], 2)))
    _, tests, _ = C.qc_summary()
    for t, n, p in tests:
        if isnum(p) and p < 0.05:
            items.append(("QC", "QC exclusion differs between groups in %s (%d excluded)" % (t, n), "p=%s" % fmt(p, 2)))
    fig = cpage(pdf, "key findings (automatic)", "all results with q < %g (tests within tables) or p < 0.05 "
                "(PERMANOVA, covariates, exclusion); exploratory; 'adjusted' rows from --model" % SIG_Q, "findings")
    ax = fig.add_axes([0.02, 0.04, 0.96, 0.82])
    if items:
        R.table(ax, [[a, R.short(f, 140), st_] for a, f, st_ in items[:40]], ["area", "finding", "statistic"],
                fontsize=6.3, left=(0, 1, 2), col_widths=[0.16, 0.70, 0.14], rowscale=1.12,
                title="%d findings%s" % (len(items), " (first 40)" if len(items) > 40 else ""))
    else:
        R.note(ax, "no result below the thresholds")
    pdf.savefig(fig)
    plt.close(fig)


def _short(x, n):
    return R.short(str(x), n)


def page_highlights_crispr(pdf, C, rows, top=30):
    if not rows:
        return False
    fig = cpage(pdf, "highlights — CRISPR systems", "expressed CRISPR arrays (≥ 20 reads) of all libraries grouped by "
                "repeat (same repeat in either orientation = one system); sorted by number of samples; "
                "RPK = reads per kb of array", "highlights_crispr")
    n2 = sum(1 for r in rows if r["samples"] >= 2)
    fig.text(0.03, 0.87, "%d systems, %d in ≥ 2 samples, %d with a tracrRNA candidate" % (
        len(rows), n2, sum(1 for r in rows if r["tracr_samples"])), fontsize=8, fontweight="bold")
    ax = fig.add_axes([0.02, 0.04, 0.95, 0.80])
    R.table(ax, [[str(i + 1), _short(r["taxon"], 34), _short(r["cas_type"].replace("no cas genes", "none").replace("no Cas genes", "none"), 12), r["samples"], r["groups"],
                  r["libraries"], fmt(r["median_reads_per_kb"], 3), fmt(r["max_reads_per_kb"], 3),
                  fmt(r["median_spacers"], 3), r["tracr_samples"] or "–", _short(r["repeat"], 38)]
                 for i, r in enumerate(rows[:top])],
            ["#", "taxon (main)", "Cas", "samples", "samples per group", "libraries", "median RPK", "max RPK",
             "spacers", "tracrRNA", "repeat"], fontsize=5.8, left=(1, 2, 4, 5, 10), rowscale=1.1,
            col_widths=[0.025, 0.19, 0.05, 0.05, 0.12, 0.08, 0.065, 0.065, 0.05, 0.05, 0.255])
    pdf.savefig(fig)
    plt.close(fig)
    return True


def page_highlights_crispr_group(pdf, C, rows, top=30):
    if rows is None:
        return False
    fig = cpage(pdf, "highlights — CRISPR systems of one group", "systems expressed (≥ 20 reads) in ≥ %d samples of "
                "a group and more often there than in the other groups of the same site (one-sided Fisher test, "
                "p < 0.05; q = BH-FDR); 'only' = in none of the compared samples; sorted by q, p"
                % C.a.highlights_min_samples, "highlights_crispr_group")
    if not rows:
        R.note(fig.add_axes([0.03, 0.4, 0.9, 0.3]), "no CRISPR system more often expressed in one group (p < 0.05)")
        pdf.savefig(fig)
        plt.close(fig)
        return True
    n_ex = sum(1 for r in rows if r["exclusive"] == "yes")
    nq = sum(1 for r in rows if r["q"] < 0.1)
    fig.text(0.03, 0.87, "%d systems with p < 0.05 (%d only in one group); %s" % (
        len(rows), n_ex, "%d with q < 0.1" % nq if nq else "none significant after FDR correction (q ≥ 0.1): "
        "candidates only"), fontsize=8, fontweight="bold", color="#333" if nq else "#B00000")
    ax = fig.add_axes([0.02, 0.04, 0.95, 0.80])
    R.table(ax, [[r["group"], "only" if r["exclusive"] == "yes" else "enriched",
                  "%d / %d (%.0f %%)" % (r["samples_in_group"], r["group_size"], 100 * r["share_in_group"]),
                  "%s: %d / %d" % (_short(r["compared_with"], 14), r["samples_other"], r["other_size"]),
                  fmt(r["p"], 2), fmt(r["q"], 2), _short(r["taxon"], 34),
                  _short(r["cas_type"].replace("no cas genes", "none").replace("no Cas genes", "none"), 12),
                  fmt(r["median_reads_per_kb"], 3), r["tracr_samples"] or "–", _short(r["repeat"], 36)]
                 for r in rows[:top]],
            ["group", "", "samples in group", "compared with", "p", "q", "taxon (main)", "Cas", "median RPK",
             "tracrRNA", "repeat"], fontsize=5.8, left=(0, 1, 6, 7, 10), rowscale=1.1,
            col_widths=[0.04, 0.045, 0.09, 0.12, 0.045, 0.045, 0.18, 0.06, 0.055, 0.045, 0.255])
    pdf.savefig(fig)
    plt.close(fig)
    return True


PRESENCE_SETS = (("genus", "genus", 10), ("species", "species", 10), ("rfam_family", "Rfam family", 10),
                 ("mirna", "miRNA", 10), ("fn_kegg_ko", "KEGG KO", 10))


def group_specific_features(C, H, gb, uni, min_samples):
    """presence / absence per group for all feature classes: highlights (CRISPR systems, novel sRNA clusters, Rfam
    families of novel loci; per sample) and the count tables (genera, species, Rfam families, miRNAs, KOs; per
    fraction, present = >= 10 reads); with a block column (site) each group is compared within its blocks"""
    items = []
    for r in H.get("crispr") or []:
        items.append(dict(cls="CRISPR system", feature="%s · %s" % (r["taxon"] or "?", r["cas_type"]),
                          info=r["repeat"], samples=r["_samples"], universe=uni))
    for r in H.get("novel") or []:
        if r["samples"] >= min_samples:
            items.append(dict(cls="novel sRNA cluster", feature=r["cluster"], info="%s · %s nt · %s" % (
                r["taxon"] or "?", fmt(r["median_length"], 3), r["locus_class"]), samples=r["_samples"],
                universe=uni))
    for r in H.get("rfam") or []:
        items.append(dict(cls="Rfam family (novel loci)", feature=r["family"], info=r["description"],
                          samples=r["_samples"], universe=uni))
    for t in C.types:
        ks = C.keys(t)
        tu = sorted({s for s, _ in ks})
        for mx, lab, minr in PRESENCE_SETS:
            pres = collections.defaultdict(set)
            for k in ks:
                for f, v in (C.vec.get(k, {}).get(mx) or {}).items():
                    if f and v >= minr:
                        pres[f].add(k[0])
            for f, ss in pres.items():
                if len(ss) >= min_samples:
                    items.append(dict(cls="%s · %s" % (lab, t.replace("sRNA-", "")), feature=f,
                                      info=tname(f, 60).replace(f, "").strip() if mx.startswith("fn_") else "",
                                      samples=ss, universe=tu))
    return HL.group_presence(items, gb, min_samples=min_samples)


def page_highlights_group(pdf, C, rows, top=32):
    if rows is None:
        return False
    fig = cpage(pdf, "highlights — features of one group (all classes)", "presence in ≥ %d samples of a group and more "
                "often than in the other groups of the same site (Fisher test on presence, p < 0.05; q = BH-FDR per "
                "class); present = expressed / ≥ 10 reads"
                % C.a.highlights_min_samples, "highlights_group")
    if not rows:
        R.note(fig.add_axes([0.03, 0.4, 0.9, 0.3]), "no feature more often present in one group (p < 0.05)")
        pdf.savefig(fig)
        plt.close(fig)
        return True
    cls = collections.OrderedDict()
    for r in rows:
        c = cls.setdefault(r["cls"], [0, 0, r["tests_in_class"]])
        c[0] += 1
        c[1] += r["q"] < 0.1
    fig.text(0.03, 0.875, "per class — p < 0.05 / q < 0.1 / tested:  " + ";  ".join(
        "%s %d / %d / %d" % (k, v[0], v[1], v[2]) for k, v in sorted(cls.items(), key=lambda x: -x[1][1])),
        fontsize=6.3, wrap=True)
    ax = fig.add_axes([0.02, 0.04, 0.95, 0.78])
    R.table(ax, [[_short(r["cls"], 24), _short(r["feature"], 38), _short(r["info"], 40), r["group"],
                  "only" if r["exclusive"] == "yes" else "",
                  "%d / %d (%.0f %%)" % (r["samples_in_group"], r["group_size"], 100 * r["share_in_group"]),
                  _short(r["compared_with"], 16),
                  "%d / %d (%.0f %%)" % (r["samples_other"], r["other_size"], 100 * r["share_other"]), fmt(r["p"], 2),
                  fmt(r["q"], 2)] for r in rows[:top]],
            ["class", "feature", "", "group", "", "samples in group", "compared with", "samples there", "p", "q"],
            fontsize=5.8, left=(0, 1, 2, 3, 4, 6), rowscale=1.1,
            col_widths=[0.12, 0.2, 0.2, 0.04, 0.035, 0.09, 0.09, 0.09, 0.05, 0.05])
    pdf.savefig(fig)
    plt.close(fig)
    return True


def page_heterogeneity(pdf, C, HZ, top=32):
    disp, load, feat = HZ.get("disp"), HZ.get("load"), HZ.get("feat")
    if not disp and not load:
        return False
    fig = cpage(pdf, "heterogeneity — spread of the patients and outlier load", "disease group vs controls of the "
                "same block; spread: distance to the group centroid (Bray-Curtis PCoA, permutations); outlier load: "
                "share of features outside the 5–95 % control range", "heterogeneity")
    if disp:
        ax = fig.add_axes([0.03, 0.60, 0.94, 0.24])
        R.table(ax, [[r["type"].replace("sRNA-", ""), NMDS_LABELS.get(r["set"], r["set"]), r["comparison"],
                      "%d / %d" % (r["n_cases"], r["n_controls"]), fmt(r["median_dist_cases"], 3),
                      fmt(r["median_dist_controls"], 3), fmt(r["ratio"], 3), fmt(r["p"], 2), fmt(r["q"], 2)]
                     for r in disp],
                ["fraction", "profile", "comparison", "n patients / controls", "spread patients (median)",
                 "spread controls", "ratio (means)", "p", "q"], fontsize=5.6, left=(0, 1, 2),
                col_widths=[0.07, 0.2, 0.17, 0.11, 0.12, 0.1, 0.08, 0.07, 0.07], title="spread (dispersion)",
                rowscale=1.0)
        bold_rows(ax, [isnum(r["q"]) and r["q"] < SIG_Q for r in disp])
    if load:
        ax = fig.add_axes([0.03, 0.03, 0.94, 0.50])
        rows = sorted(load, key=lambda r: (r["type"], r["comparison"], r["set"]))
        R.table(ax, [[r["type"].replace("sRNA-", ""), r["set"], r["comparison"], r["features"],
                      "%.1f %%" % (100 * r["median_load_cases"]), "%.1f %%" % (100 * r["median_load_controls"]),
                      "%d / %d" % (r["cases_above_control_p90"], r["n_cases"]), fmt(r["p"], 2), fmt(r["q"], 2)]
                     for r in rows[:40]],
                ["fraction", "features", "comparison", "tested", "outlier load patients (median)", "controls",
                 "patients > 90th pct of controls", "p", "q"], fontsize=5.2, left=(0, 1, 2), rowscale=0.92,
                col_widths=[0.07, 0.12, 0.17, 0.07, 0.16, 0.09, 0.16, 0.07, 0.07], title="outlier load")
        bold_rows(ax, [isnum(r["q"]) and r["q"] < SIG_Q for r in rows[:40]])
    pdf.savefig(fig)
    plt.close(fig)
    if feat:
        fig = cpage(pdf, "heterogeneity — features changed in a subset of patients", "patients above the 95th (below "
                    "the 5th) percentile of the controls of the same block vs controls there (Fisher test; q = BH-FDR "
                    "per fraction / set / comparison / direction)", "heterogeneity_features")
        nq = sum(1 for r in feat if r["q"] < 0.1)
        fig.text(0.03, 0.875, "%d feature tests with p < 0.05, %d with q < 0.1" % (len(feat), nq), fontsize=8,
                 fontweight="bold", color="#333" if nq else "#B00000")
        ax = fig.add_axes([0.02, 0.04, 0.95, 0.80])
        R.table(ax, [[r["type"].replace("sRNA-", ""), _short(r["set"], 14), r["comparison"],
                      _short(("%s %s" % (r["feature"], r["name"])).strip(), 40), "↑" if r["direction"] == "up" else "↓",
                      "%d / %d (%.0f %%)" % (r["cases_outside"], r["n_cases"], 100 * r["share_cases"]),
                      "%d / %d" % (r["controls_outside"], r["n_controls"]), fmt(r["p"], 2), fmt(r["q"], 2),
                      _short(r["patients"], 46)] for r in feat[:top]],
                ["fraction", "set", "comparison", "feature", "", "patients outside", "controls", "p", "q",
                 "patients"], fontsize=5.6, left=(0, 1, 2, 3, 9), rowscale=1.1,
                col_widths=[0.05, 0.08, 0.13, 0.22, 0.02, 0.1, 0.06, 0.05, 0.05, 0.24])
        bold_rows(ax, [r["q"] < 0.1 for r in feat[:top]])
        pdf.savefig(fig)
        plt.close(fig)
    return True


# report structure: (key, title, what the section answers, [(page key, short description)])
REPORT_SECTIONS = [
    ("summary", "Summary", "the cohort, the libraries in the statistics and the automatic key findings",
     [("overview", "groups, libraries per fraction, exclusions, covariates by group"),
      ("findings", "key findings of all sections (automatic)")]),
    ("quality", "Data quality", "can the libraries be compared? QC verdicts, exclusions, depth, read fate and RNA "
     "integrity by group", [("qc_flags", "QC verdicts and flagged libraries"), ("qc", "library metrics by group"),
                            ("integrity", "RNA integrity indicators by group and vs wet-lab RIN")]),
    ("composition", "Composition", "what the sRNA libraries contain: taxa (rRNA), beta diversity, RNA classes, "
     "human miRNAs, functions (KEGG / Pfam / COG of CDS reads) and antisense transcription",
     [("taxa", "rRNA genera by group"), ("beta", "Bray-Curtis PCoA and PERMANOVA"), ("nmds", "NMDS per profile"),
      ("classes", "RNA classes and Rfam families"), ("mirna", "human miRNAs, origin markers"),
      ("functions", "functional categories"), ("antisense", "antisense index, hotspots")]),
    ("comparisons", "Group comparisons (unadjusted)", "exploratory differences between the groups, within blocks "
     "(site), between fractions and against metaG / metaT",
     [("matched", "block design, pooled and within-block comparisons"), ("features", "Kruskal-Wallis per feature"),
      ("paired", "sRNA-Int vs sRNA-Ext"), ("omics", "agreement with metaG / metaT")]),
    ("models", "Adjusted models", "group and covariate effects adjusted for confounders (OLS on CLR, robust SE); "
     "the models, their term of interest and all results",
     [("models_intro", "the models: formula, term of interest, adjustment, subset, libraries"),
      ("models_overview", "all models at a glance: significant features per outcome set"),
      ("models", "per model and fraction: tables and volcano plots")]),
    ("heterogeneity", "Heterogeneity and subgroups", "do only some patients differ (disease subtypes)? spread, "
     "outlier load and features changed in a subset of patients",
     [("heterogeneity", "dispersion, outlier load, features outside the control range with patient IDs")]),
    ("highlights", "Highlights across samples", "the most notable single findings: CRISPR systems, recurrent novel "
     "sRNAs, ncRNA families, features of one group, chimeras (experimental)",
     [("highlights_crispr", "CRISPR systems, also by group"), ("highlights_group", "features of one group, all classes"),
      ("highlights_novel", "recurrent novel sRNAs, Rfam families among novel loci"),
      ("highlights_chimera", "candidate RNA–RNA pairs (module chimera, experimental)")]),
    ("methods", "Methods", "statistics, parameters and the automated Methods text",
     [("methods", "methods and caveats"), ("methods_text", "Methods (journal style)")]),
]


def section_label(key):
    for i, (k, title, _, _) in enumerate(REPORT_SECTIONS, 1):
        if k == key:
            return "%d · %s" % (i, title)
    return key


def page_contents(pdf, C, sections):
    fig = cpage(pdf, "contents", "%d samples · %s · groups: %s%s" % (
        len(C.samples), ", ".join(C.types), ", ".join(C.groups) or "–",
        " · block: %s" % C.a.match_col if C.a.match_col else ""), "contents")
    y = 0.86
    for i, (k, title, what, pages) in enumerate(sections, 1):
        fig.text(0.04, y, "%d  %s" % (i, title), fontsize=11, fontweight="bold", color="#0b3d66", va="top")
        fig.text(0.30, y - 0.002, what, fontsize=7.5, color="#444", va="top", wrap=True)
        y -= 0.032
        for pk, desc in pages:
            fig.text(0.06, y, "• " + desc, fontsize=7, color="#222", va="top")
            y -= 0.019
        y -= 0.012
    pdf.savefig(fig)
    plt.close(fig)
    return True


def model_description(C, mi):
    """(term of interest, kind, adjusted for, subset, description) of a model"""
    terms = C.model_terms(mi)
    if not terms:
        return "", "", "", "", ""
    focus, others = terms[0], terms[1:]
    sub = C.model_subset(mi)
    subset = "%s = %s" % (sub[0], "/".join(sorted(sub[1]))) if sub else "all samples"
    ref = C.a.reference or (C.groups[0] if C.groups else "")
    if focus == C.group_col:
        kind = "groups"
        desc = "group differences: %s vs %s" % (", ".join(g for g in C.groups if g != ref), ref)
    else:
        k = C.kinds.get(focus, "")
        if k == "numeric":
            kind = "numeric (per SD)"
            desc = "change per SD of %s" % focus
        else:
            lv = sorted({C.meta_val(s, focus) for s in C.samples} - {"NA"})
            kind = "categorical (%s)" % "/".join(lv[:4]) if lv else "categorical"
            desc = "%s (levels vs reference)" % focus
    if sub:
        desc += " (%s only)" % "/".join(sorted(sub[1]))
    return focus, kind, ", ".join(others) or "–", subset, desc


def page_models_intro(pdf, C):
    if not C.models():
        return False
    notes = [x.strip() for x in (C.a.model_notes or "").split(";") if x.strip()]
    fig = cpage(pdf, "adjusted models — the models", "OLS on CLR abundances / metric values per feature, robust "
                "(HC3) standard errors, BH-FDR per term and outcome set; term of interest = the first term of the "
                "formula; numeric terms per SD", "models_intro")
    rows = []
    for mi, formula in enumerate(C.models()):
        focus, kind, adj, subset, desc = model_description(C, mi)
        lab = C.model_label(mi)
        note = ""
        for n in notes:
            if n.split("=", 1)[0].strip() in (lab, lab.split()[0]) and "=" in n:
                note = n.split("=", 1)[1].strip()
        if not note and mi < len(notes) and "=" not in notes[mi]:
            note = notes[mi]
        libs = []
        for t in C.types:
            md = C.model_data(t, mi)
            libs.append(str(len(md[0])) if md else "–")
        info = "; ".join(sorted({x for t in C.types for x in ((C.model_data(t, mi) or [None] * 5)[4] or [])}))
        rows.append([lab, R.short(note or desc, 70), focus, R.short(kind, 22), R.short(adj, 40), subset] + libs +
                    [R.short(info or "–", 40)])
    ax = fig.add_axes([0.02, 0.06, 0.95, 0.80])
    R.table(ax, rows, ["model", "question", "term of interest", "kind", "adjusted for", "samples"] +
            ["n %s" % t.replace("sRNA-", "") for t in C.types] + ["design notes"], fontsize=6.2,
            left=(0, 1, 2, 3, 4, 5, 6 + len(C.types)), rowscale=1.25,
            col_widths=[0.05, 0.25, 0.1, 0.1, 0.17, 0.08] + [0.04] * len(C.types) + [0.17 - 0.04 * (len(C.types) - 2)])
    pdf.savefig(fig)
    plt.close(fig)
    return True


def page_highlights_novel(pdf, C, rows, rfam, top=28):
    if not rows and not rfam:
        return False
    if rows:
        multi = [r for r in rows if r["samples"] >= 2]
        fig = cpage(pdf, "highlights — recurrent novel sRNAs", "novel sRNA loci without Rfam hit (module novel) of all "
                    "libraries, clustered across samples by shared sequence (exact 21-mers); sorted by samples; "
                    "MFE/nt: RNAfold energy per nt (lower = more structured)", "highlights_novel")
        fig.text(0.03, 0.87, "%d clusters, %d in ≥ 2 samples, %d in ≥ 10 samples" % (
            len(rows), len(multi), sum(1 for r in rows if r["samples"] >= 10)), fontsize=8, fontweight="bold")
        ax = fig.add_axes([0.02, 0.04, 0.95, 0.80])
        R.table(ax, [[r["cluster"], r["samples"], r["groups"], r["libraries"], _short(r["locus_class"], 18),
                      fmt(r["median_length"], 3), fmt(r["median_reads_per_sample"], 3), fmt(r["median_mfe_per_nt"], 2),
                      _short(r["taxon"] or "–", 30), _short(r["rep_sequence"], 44)] for r in rows[:top]],
                ["cluster", "samples", "samples per group", "libraries", "locus", "nt", "reads/sample",
                 "MFE/nt", "taxon (representative)", "representative sequence (5′→3′)"], fontsize=5.6,
                left=(0, 2, 3, 4, 8, 9), rowscale=1.1,
                col_widths=[0.07, 0.045, 0.11, 0.075, 0.085, 0.045, 0.06, 0.045, 0.16, 0.305])
        pdf.savefig(fig)
        plt.close(fig)
    if rfam:
        fig = cpage(pdf, "highlights — known ncRNA families among the novel loci", "best Rfam hit (cmscan) of the "
                    "novel loci of all libraries, rRNA / tRNA families left out; sorted by number of samples",
                    "highlights_rfam")
        ax = fig.add_axes([0.02, 0.04, 0.95, 0.82])
        R.table(ax, [[r["family"], r["accession"], _short(r["description"], 60), r["samples"], r["groups"],
                      r["libraries"], r["loci"]] for r in rfam[:32]],
                ["family", "Rfam", "description", "samples", "samples per group", "libraries", "loci"],
                fontsize=6, left=(0, 1, 2, 4, 5), rowscale=1.1, col_widths=[0.15, 0.07, 0.38, 0.06, 0.15, 0.11, 0.06])
        pdf.savefig(fig)
        plt.close(fig)
    return True


def page_highlights_chimera(pdf, C, lib, rows, top=24):
    if not lib:
        return False
    fig = cpage(pdf, "highlights — chimeric reads (candidate RNA–RNA pairs)", "module chimera; sRNA-like–mRNA = "
                "ncRNA / novel / intergenic arm with a 5′ UTR / CDS arm; complementarity: arms vs re-paired arms "
                "(Mann–Whitney); same bin: inter-contig chimeras with both arms in one bin vs expected",
                "highlights_chimera")
    summ = []
    for t in sorted({x["type"] for x in lib}, key=lambda t: (t == "metaT", t)):
        L = [x for x in lib if x["type"] == t]
        pc = [x["p_complementarity"] for x in L if not math.isnan(x["p_complementarity"])]
        sb = [(x["same_bin"], x["same_bin_expected"]) for x in L if not math.isnan(x["same_bin"])]
        summ.append([t, len(L), fmt(float(np.median([x["reads_used"] for x in L])), 3),
                     fmt(float(np.median([x["pairs_q005"] for x in L])), 3), len(pc),
                     sum(1 for p in pc if p < 0.05),
                     fmt(float(np.median([x[0] for x in sb])), 2) if sb else "NA",
                     fmt(float(np.median([x[1] for x in sb])), 2) if sb else "NA"])
    ax = fig.add_axes([0.03, 0.70, 0.94, 0.17])
    R.table(ax, summ, ["data", "libraries", "median chimeric reads", "median pairs q<0.05",
                       "libraries tested (sRNA–mRNA)", "complementarity p<0.05", "median same-bin share",
                       "median expected"], fontsize=6.5, rowscale=1.2,
            col_widths=[0.10, 0.08, 0.13, 0.13, 0.15, 0.13, 0.14, 0.12])
    if rows:
        ax = fig.add_axes([0.02, 0.04, 0.95, 0.60])
        R.table(ax, [[_short(r["group"], 16), _short("%s (%s)" % (r["name_a"], r["type_a"]), 40),
                      _short("%s (%s)" % (r["name_b"], r["type_b"]), 40), r["samples"], r["libraries"],
                      fmt(r["reads"], 3), fmt(r["median_dG"], 3), fmt(r["median_p_complementarity"], 2),
                      r["same_bin"], _short(r["taxon"], 26)] for r in rows[:top]],
                ["group", "partner a", "partner b", "samples", "libraries", "reads", "ΔG", "p compl.", "same bin",
                 "taxon"], fontsize=5.8, left=(0, 1, 2, 4, 8, 9), rowscale=1.1,
                col_widths=[0.09, 0.22, 0.22, 0.05, 0.08, 0.05, 0.05, 0.06, 0.05, 0.13],
                title="candidate pairs with q < 0.05, pooled over libraries by partner (sRNA-like–mRNA first)")
    pdf.savefig(fig)
    plt.close(fig)
    return True


def page_methods(pdf, C):
    fig = R.page(pdf, COHORT_NAME, "methods and caveats")
    txt = [
        "Inputs: per-sample pipeline outputs (Stats, Analysis/<type>/*, comparison) of %d samples; metadata %s "
        "(ID column '%s', group column '%s')." % (len(C.samples), C.a.metadata or "none", getattr(C, "id_col", "-"),
                                                  C.group_col or "-"),
        "Group tests: Kruskal-Wallis with tie correction, pairwise Mann-Whitney (effect: P(group1 > group2)), "
        "Benjamini-Hochberg FDR within each table. Paired Int vs Ext: Wilcoxon signed-rank (normal approximation). "
        "Correlations: Spearman (integrity vs wet-lab QC: BH-FDR within each QC measure). Covariate balance: Kruskal-Wallis (numeric) or chi-square (categorical).",
        "QC: every library is checked against the QC rules (%s) and manual decisions (%s); libraries with "
        "%s are left out of all statistics and exported matrices of their fraction (%d of %d libraries). Paired "
        "Int vs Ext tests use samples with both fractions retained." % (
            C.a.qc_rules or "scripts/qc_rules.tsv", C.a.qc_exclude or "none",
            {"FAIL": "a FAIL", "WARN": "a WARN or FAIL", "none": "(none)"}[C.a.qc_drop], len(C.excluded), len(C.m)),
        "Block design (--block %s): per group the number of block levels shared with the reference decides the "
        "design — nested (≤ 2, e.g. site) or crossed (≥ 3, e.g. family). Comparisons: pooled (all libraries, "
        "Mann-Whitney, exact for n ≤ 20; PERMANOVA unrestricted) and within (only levels with both groups; van "
        "Elteren stratified Wilcoxon, weights 1/(n+1), normal approximation; PERMANOVA permuted within levels); "
        "nested designs also contrast the reference subgroups (block effect alone). Poolability: block PERMANOVA "
        "within the reference group, or — with too few reference libraries per level — on all libraries, block after "
        "group (partial R²), permuted within groups. BH-FDR within each comparison. Crossed designs fit each model "
        "twice: pooled with block-clustered SE (CR1, t with G − 1 df) and within-block (fixed effects by demeaning, "
        "clustered SE)." % (C.a.match_col or "not used"),
        "miRNAs: analysed only in fractions with ≥ 6 libraries of ≥ %g human miRNA reads (%s); other fractions "
        "(typically cell-free sRNA-Ext, with almost no human miRNAs) are left out of miRNA plots and tests." % (
            C.a.min_mirna_reads, ", ".join(C.mirna_types()) or "none"),
        "Adjusted models (M1, M2, … = --model %s): per feature OLS of the CLR-transformed abundance (pseudocount 0.5; features with "
        "≥ 10 reads in ≥ 10 %% of the libraries) or of the metric on the model terms; categorical terms as treatment "
        "contrasts (group reference %s), numeric terms per SD; HC3 robust standard errors, t-test, BH-FDR per term "
        "within each outcome set; complete cases only; collinear terms dropped." % (
            C.a.model or "not used", C.a.reference or (C.groups[0] if C.groups else "-")),
        "Composition: relative abundances of Bracken genus estimates; features kept if ≥ 10 reads in ≥ 10 %% of the "
        "samples of a fraction. Beta diversity: Bray-Curtis, PCoA, PERMANOVA with %d permutations (marginal per "
        "variable)." % C.a.permutations,
        "Caveats: apart from the adjusted models (--model), the tests are exploratory and not adjusted for covariates "
        "(sex, age, constipation, medication, integrity/RIN, sequencing depth). Human miRNA counts are low (often < 1000 "
        "reads per sample); marker shares are only computed with ≥ 20 miRNA reads.",
        "Disease-associated factors (e.g. constipation and transit time, medication, diet) change the stool "
        "microbiome and potentially sRNA integrity by themselves; related or repeated samples (families, time "
        "points) are not independent.",
    ]
    y = 0.88
    for t in txt:
        import textwrap
        w = "\n".join(textwrap.wrap(t, 170))
        fig.text(0.04, y, w, fontsize=8, va="top", linespacing=1.4)
        y -= 0.05 + 0.025 * w.count("\n")
    pdf.savefig(fig)
    plt.close(fig)


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True, help="pipeline output directory (-o of srna_pipeline.sh)")
    ap.add_argument("--metadata", default="", help="metadata TSV/CSV/XLSX (one row per sample)")
    ap.add_argument("--id-col", default="", help="sample ID column (default: first column)")
    ap.add_argument("--group-col", default="group", help="group column ('' = no groups)  [group]")
    ap.add_argument("--group-order", default="", help="comma-separated group order (e.g. HC,PD,RA or HC,T1D)")
    ap.add_argument("--covariates", default="", help="comma-separated covariate columns (default: all columns)")
    ap.add_argument("--samples", default="", help="restrict to the samples/types of this sample sheet")
    ap.add_argument("--rna-qc", default="", help="wet-lab RNA QC TSV (as for the per-sample report)")
    ap.add_argument("--mirna-markers", default="", help="miRNA origin marker TSV (default: next to the script)")
    ap.add_argument("--cohort-dir", default="", help="output directory  [<outdir>/Cohort]")
    ap.add_argument("--permutations", type=int, default=999)
    ap.add_argument("--block", "--match-col", dest="match_col", default="",
                    help="block column: a second grouping such as site, family or batch (nested or crossed with the "
                         "groups; pooled and within-block comparisons, poolability check, block-aware models)")
    ap.add_argument("--reference", default="", help="reference group for matched comparisons (default: first group)")
    ap.add_argument("--qc-rules", default="", help="QC rules TSV (default scripts/qc_rules.tsv)")
    ap.add_argument("--qc-exclude", default="", help="manual QC decisions: sample, type (or *), level, reason")
    ap.add_argument("--model", default="", help="covariate-adjusted models, comma-separated, terms joined by '+', "
                                                  "e.g. 'group+site, group+site+sex+age' (metadata columns or library "
                                                  "metrics such as qc_lRNA_RIN); '' = none")
    ap.add_argument("--cohort-name", default="cohort", help="name of the cohort, used as report title [cohort]")
    ap.add_argument("--model-notes", default="", help="short question per model, ';'-separated, 'M9=PD motor "
                    "severity; M13=RA immunosuppressants' or in model order (contents of the models page)")
    ap.add_argument("--no-heterogeneity", action="store_true",
                    help="skip the heterogeneity pages (dispersion, patients outside the control range)")
    ap.add_argument("--highlights-min-samples", type=int, default=3,
                    help="group-specific CRISPR systems: minimum samples of the group with the system [3]")
    ap.add_argument("--highlights-top", type=int, default=300,
                    help="novel loci per library (most reads) used for the cross-sample novel sRNA clusters [300]; 0 = all")
    ap.add_argument("--nmds", default="genus,fn_kegg_pathway,mirna,rfam_family",
                    help="profiles for NMDS pages (genus, family, mirna, rfam_family, rna_class, fn_kegg_pathway, "
                         "fn_kegg_ko, fn_pfam; '' = none)")
    ap.add_argument("--nmds-starts", type=int, default=10, help="random starts for NMDS (plus one PCoA start)")
    ap.add_argument("--term-names", default="", help="comma-separated TSVs term<TAB>name for KEGG / Pfam labels")
    ap.add_argument("--min-mirna-reads", type=float, default=50,
                    help="libraries with fewer human miRNA reads are left out of the miRNA tests and PCoA")
    ap.add_argument("--sig-q", type=float, default=0.1, help="FDR q below which results are bold / highlighted")
    ap.add_argument("--outlier-z", type=float, default=4.0,
                    help="flag libraries with |robust z| above this in key metrics of their fraction as WARN (0 = off)")
    ap.add_argument("--qc-drop", default="FAIL", choices=("FAIL", "WARN", "none"),
                    help="libraries left out of the statistics: FAIL (default), WARN = WARN and FAIL, none")
    ap.add_argument("--html", action="store_true", help="also write cohort.report.html")
    ap.add_argument("--methods-style", default="", help="add the automated Methods (srna_methods.py) in this style "
                    "(star / nature) as last pages; built after the cohort tables")
    ap.add_argument("--methods-params", default="", help="NAME<TAB>VALUE configuration for the Methods")
    ap.add_argument("--methods-envs", default="", help="environments for the tool versions (conda-meta)")
    ap.add_argument("--methods-resources", default="", help="resources TSV (default scripts/methods_resources.tsv)")
    a = ap.parse_args()
    global COHORT_NAME
    COHORT_NAME = a.cohort_name or "cohort"
    if not a.metadata:
        a.group_col = ""

    R.load_markers(a.mirna_markers or os.path.join(os.path.dirname(os.path.abspath(__file__)), "mirna_origin_markers.tsv"))
    R.load_rna_qc(a.rna_qc)
    global SIG_Q
    SIG_Q = a.sig_q
    load_term_names(a.term_names)
    C = Cohort(a)
    d = a.cohort_dir or os.path.join(a.outdir, "Cohort")
    os.makedirs(os.path.join(d, "stats"), exist_ok=True)
    sys.stderr.write("[cohort] %d samples, %d sample-types (%d excluded by QC), groups: %s\n" % (
        len(C.samples), len(C.m), len(C.excluded), ", ".join(C.groups)))
    C.write_tables(d)

    # statistics tables
    for t in C.types:
        rows, pw = C.metric_tests(t)
        write_tsv(os.path.join(d, "stats", "cohort.%s.metrics_by_group.tsv" % t),
                  ["metric", "n", "H", "p", "q"] + ["median_%s" % g for g in C.groups], rows)
        write_tsv(os.path.join(d, "stats", "cohort.%s.metrics_pairwise.tsv" % t),
                  ["metric", "group1", "group2", "n1", "n2", "P(g1>g2)", "p", "q"], pw)
        for mx in MATRICES:
            if mx == "mirna" and t not in C.mirna_types():
                continue                       # too few miRNA reads in this fraction (counts are still exported)
            if mx in FUNC_AS:
                continue                       # antisense function counts: exported only
            ft = C.feature_tests(t, mx)
            if ft:
                if mx in FUNC_MX:
                    ft = [[r[0], TERM_NAMES.get(r[0], "")] + r[1:] for r in ft]
                write_tsv(os.path.join(d, "stats", "cohort.%s.%s.by_group.tsv" % (t, mx)),
                          ["feature"] + (["name"] if mx in FUNC_MX else []) + ["prevalence", "H", "p", "q"] +
                          ["mean_rel_%s" % g for g in C.groups], ft)
            if mx in ("genus", "mirna", "family", "fn_kegg_pathway", "fn_kegg_ko"):
                b = C.beta(t, mx)
                if b and b[3]:
                    write_tsv(os.path.join(d, "stats", "cohort.%s.%s.permanova.tsv" % (t, mx)),
                              ["variable", "n", "levels", "pseudo_F", "R2", "p"], b[3])
    for t in C.types:
        for label, ka, kb, kind in C.matched_sets(t):
            r = C.two_set_tests(t, ka, kb, blocked=(kind == "within"))
            tag = re.sub(r"[^A-Za-z0-9]+", "_", label.split(" (")[0]).strip("_")
            write_tsv(os.path.join(d, "stats", "cohort.%s.matched.%s.metrics.tsv" % (t, tag)),
                      ["metric", "n_a", "n_b", "median_a", "median_b", "P(a>b)", "p", "q"], r["metrics"])
            for mx in ("genus", "mirna", "rfam_family") + FUNC_MX:
                if r[mx]:
                    write_tsv(os.path.join(d, "stats", "cohort.%s.matched.%s.%s.tsv" % (t, tag, mx)),
                              ["feature", "prevalence", "mean_rel_a", "mean_rel_b", "P(a>b)", "p", "q"], r[mx])
            if r["permanova"]:
                write_tsv(os.path.join(d, "stats", "cohort.%s.matched.%s.permanova.tsv" % (t, tag)),
                          ["comparison", "n_a", "n_b", "pseudo_F", "R2", "p"], [[label] + r["permanova"]])
    AS_DBS = ("kegg_pathway", "kegg_module", "kegg_ko", "pfam", "cog")
    for t in C.types:
        for db in AS_DBS:
            hs = C.antisense_hotspots(t, db)
            if hs:
                write_tsv(os.path.join(d, "stats", "cohort.%s.antisense_hotspots.%s.tsv" % (t, db)),
                          ["term", "name", "libraries", "median_antisense_index", "median_library_index",
                           "median_log2_ratio", "fraction_libraries_above", "p", "q"], hs)
            if len(C.groups) > 1:
                at = C.asidx_tests(t, db)
                if at:
                    write_tsv(os.path.join(d, "stats", "cohort.%s.antisense_index.%s.by_group.tsv" % (t, db)),
                              ["term", "name", "libraries", "H", "p"] + ["median_%s" % g for g in C.groups] + ["q"], at)
            for label, ka, kb, kind in C.matched_sets(t):
                at = C.asidx_tests(t, db, ka, kb)
                if at:
                    tag = re.sub(r"[^A-Za-z0-9]+", "_", label.split(" (")[0]).strip("_")
                    write_tsv(os.path.join(d, "stats", "cohort.%s.matched.%s.antisense_index.%s.tsv" % (t, tag, db)),
                              ["term", "name", "libraries", "median_a", "median_b", "P(a>b)", "p", "q"], at)
            ks, terms, X, lib = C.asidx(t, db)
            if terms:
                write_tsv(os.path.join(d, "tables", "cohort.%s.antisense_index.%s.tsv" % (t, db)),
                          ["term"] + [k[0] for k in ks], [[f] + list(X[:, j]) for j, f in enumerate(terms)])
    for t in C.types:
        for mx in [m.strip() for m in a.nmds.split(",") if m.strip()]:
            r = C.nmds_result(t, mx)
            if r:
                ks, X, st, pm = r
                write_tsv(os.path.join(d, "stats", "cohort.%s.%s.nmds.tsv" % (t, mx)),
                          ["sample", "group", "NMDS1", "NMDS2", "stress", "permanova_R2", "permanova_p"],
                          [[k[0], C.group(k[0]), X[j, 0], X[j, 1], st, pm[1], pm[2]] for j, k in enumerate(ks)])
    if C.models():
        idx_rows = []
        for mi, formula in enumerate(C.models()):
            mlab = C.model_label(mi).replace(" ", "_")
            for t in C.types:
                md = C.model_data(t, mi)
                if not md:
                    sys.stderr.write("[cohort] %s: too few complete cases for model %s (%s)\n" % (t, mlab, formula))
                    idx_rows.append([mlab, formula, t, 0, "too few complete cases"])
                    continue
                idx_rows.append([mlab, formula, t, len(md[0]), "; ".join(md[4]) or "-"])
                write_tsv(os.path.join(d, "stats", "cohort.%s.model_%s.design.tsv" % (t, mlab)), ["sample"] + md[2],
                          [[k[0]] + list(md[1][i]) for i, k in enumerate(md[0])])
                for what in MODEL_SETS:
                    res = C.model_fit(t, what, mi)
                    if res:
                        write_tsv(os.path.join(d, "stats", "cohort.%s.model_%s.%s.tsv" % (t, mlab, what)),
                                  ["feature", "name", "n", "term", "estimate", "se", "t", "p", "q"], res)
                if md[4]:
                    sys.stderr.write("[cohort] %s model %s: %s\n" % (t, mlab, "; ".join(md[4])))
        write_tsv(os.path.join(d, "stats", "cohort.models.tsv"), ["model", "formula", "type", "libraries", "notes"],
                  idx_rows)
    pr = C.paired()
    if pr:
        write_tsv(os.path.join(d, "stats", "cohort.paired_Int_Ext.tsv"),
                  ["metric", "pairs", "median_Int", "median_Ext", "median_diff", "p", "q"], pr)
    cv = C.covariates()
    if cv:
        write_tsv(os.path.join(d, "stats", "cohort.covariates_by_group.tsv"),
                  ["covariate", "kind", "test", "p"] + C.groups, cv)
    iq = C.integrity_vs_qc()
    if iq:
        write_tsv(os.path.join(d, "stats", "cohort.integrity_vs_qc.tsv"),
                  ["type", "qc", "indicator", "n", "rho", "p", "q"], iq)

    # cross-sample highlights: CRISPR systems, recurrent novel sRNAs, chimeras
    H = {}
    try:
        keep = [k for k in C.pairs if k not in C.excluded]
        H["crispr"] = HL.crispr_systems(a.outdir, keep, C.group)
        uni = sorted({s for s, _ in keep})
        gb = {s: (C.group(s), C.meta_val(s, a.match_col) if a.match_col else "") for s in C.samples}
        H["crispr_group"] = HL.crispr_group_specific(H["crispr"], gb, uni, min_samples=a.highlights_min_samples)
        H["group_specific"] = group_specific_features(C, H, gb, uni, a.highlights_min_samples)
        H["novel"], _ = HL.novel_clusters(a.outdir, keep, C.group, top_per_lib=a.highlights_top)
        H["rfam"] = HL.novel_rfam(a.outdir, keep, C.group)
        H["chim_lib"], H["chim"] = HL.chimeras(a.outdir, sorted({s for s, _ in keep}), C.types)
        for key, name in (("crispr", "crispr_systems"), ("crispr_group", "crispr_group_specific"),
                          ("group_specific", "group_specific_features"), ("novel", "novel_srna_clusters"), ("rfam", "novel_rfam_families"),
                          ("chim_lib", "chimera_libraries"), ("chim", "chimera_pairs")):
            rows = H[key]
            if rows:
                cols = [c for c in rows[0] if not c.startswith("_")]
                write_tsv(os.path.join(d, "stats", "cohort.highlights.%s.tsv" % name), cols,
                          [[r[c] for c in cols] for r in rows])
        if H["novel"]:
            with open(os.path.join(d, "stats", "cohort.highlights.novel_srna_clusters.fa"), "w") as fh:
                for r in H["novel"]:
                    fh.write(">%s samples=%d locus=%s sample=%s/%s\n%s\n" % (
                        r["cluster"], r["samples"], r["rep_locus"], r["rep_sample"], r["rep_type"], r["rep_sequence"]))
        sys.stderr.write("[cohort] highlights: %d CRISPR systems, %d novel sRNA clusters (%d in >= 2 samples), "
                         "%d Rfam families, %d chimera libraries\n" % (
                             len(H["crispr"]), len(H["novel"]), sum(1 for r in H["novel"] if r["samples"] >= 2),
                             len(H["rfam"]), len(H["chim_lib"])))
    except Exception as e:                       # highlights must not stop the cohort report
        sys.stderr.write("[cohort] highlights skipped (%s: %s)\n" % (type(e).__name__, e))

    # heterogeneity of the disease groups (subtypes): dispersion and patients outside the control range
    HZ = {}
    if len(C.groups) > 1 and not a.no_heterogeneity:
        try:
            HZ["disp"], disp_d = HT.dispersion(C, sys.modules[__name__], perms=a.permutations)
            HZ["feat"], load_rows, HZ["load"] = HT.outliers(C, sys.modules[__name__], HL,
                                                             name_of=lambda f: tname(f, 60).replace(f, "").strip())
            for key, rows_, name in (("disp", HZ["disp"], "dispersion"), ("disp_d", disp_d, "dispersion_distances"),
                                     ("feat", HZ["feat"], "outlier_features"), ("load", HZ["load"], "outlier_load"),
                                     ("load_rows", load_rows, "outlier_load_per_library")):
                if rows_:
                    cols = list(rows_[0])
                    write_tsv(os.path.join(d, "stats", "cohort.heterogeneity.%s.tsv" % name), cols,
                              [[r[c] for c in cols] for r in rows_])
            sys.stderr.write("[cohort] heterogeneity: %d dispersion tests, %d features outside the control range "
                             "(p < 0.05, %d with q < 0.1)\n" % (len(HZ["disp"]), len(HZ["feat"]),
                                                                sum(1 for r in HZ["feat"] if r["q"] < 0.1)))
        except Exception as e:
            sys.stderr.write("[cohort] heterogeneity skipped (%s: %s)\n" % (type(e).__name__, e))

    # automated Methods (needs the tables above)
    mfile = ""
    if a.methods_style and a.methods_style != "none":
        mfile = os.path.join(d, "cohort.methods.%s.md" % a.methods_style)
        cmd = [sys.executable, os.path.join(os.path.dirname(os.path.abspath(__file__)), "srna_methods.py"), "--mode",
               "all", "--style", a.methods_style, "--outdir", a.outdir, "--cohort-dir", d, "--out", mfile]
        for opt, v in (("--params", a.methods_params), ("--envs", a.methods_envs),
                       ("--resources", a.methods_resources)):
            if v:
                cmd += [opt, v]
        r = subprocess.run(cmd, capture_output=True, text=True)
        sys.stderr.write(r.stderr)
        if r.returncode != 0:
            mfile = ""

    # report
    pdf_path = os.path.join(d, "cohort.report.pdf")
    with R.PdfPages(pdf_path) as pf:
        pdf = R.Sink(pf, html=a.html)
        pct = dict(pct=True)
        jobs = [
            ("contents", lambda: page_contents(pdf, C, REPORT_SECTIONS)),
            ("summary|overview", lambda: page_overview(pdf, C)),
            ("summary|findings", lambda: page_findings(pdf, C)),
            ("quality|qc_flags", lambda: page_qc_flags(pdf, C)),
            ("quality|qc", lambda: page_grid(pdf, C, "qc", "library metrics by group", "read fate shares, depth, mapping, "
                                     "novel loci, CRISPR; boxes per fraction", [
                ("trimmed_fraction", "trimmed (% of raw)", pct), ("adapter_dimer_fraction", "adapter dimers (%)", pct),
                ("human_fraction", "human (% of trimmed)", pct), ("assembly_mapped_fraction", "mapped to contigs (%)", pct),
                ("nonhuman_rRNA_fraction", "rRNA (% non-human)", pct), ("nonhuman_tRNA_fraction", "tRNA (% non-human)", pct),
                ("nonhuman_Rfam_fraction", "other Rfam (% non-human)", pct), ("nonhuman_unknown_fraction", "unknown (% non-human)", pct),
                ("raw_reads", "raw reads", dict(log=True)), ("novel_loci", "novel candidate loci", {}),
                ("crispr_reads_per_million_mapped", "CRISPR reads per M mapped", dict(log=True)),
                ("crispr_arrays_with_10_reads", "CRISPR arrays ≥ 10 reads", {})])),
            ("quality|integrity", lambda: page_grid(pdf, C, "integrity", "RNA integrity by group", "indicators of module fragments", [
                ("qc_lRNA_RIN", "wet-lab RIN (long RNA)", {}),
                ("integrity_tRNA_internal_fraction", "tRNA internal fragments", {}),
                ("integrity_tRNA_3p_CCA_fraction", "3′-CCA share", {}),
                ("integrity_tRNA_5p_end_entropy", "tRNA 5′-end entropy", {}),
                ("integrity_tRNA_5p_top3_share", "tRNA 5′ ends at main positions", {}),
                ("integrity_rRNA_5p_end_entropy", "rRNA 5′-end entropy", {}),
                ("integrity_rRNA_bin_evenness", "rRNA spread along 16S/23S", {}),
                ("integrity_median_read_length", "median read length", {}),
                ("integrity_tRNA_halves_fraction", "tRNA halves (stress)", {})])),
            ("quality|integrity_qc", lambda: page_integrity_qc(pdf, C)),
            ("quality|integrity_summary", lambda: page_integrity_summary(pdf, C)),
            ("composition|taxa", lambda: page_taxa(pdf, C)),
            ("composition|beta", lambda: page_beta(pdf, C)),
            ("composition|nmds", lambda: page_nmds(pdf, C)),
            ("composition|classes", lambda: page_rna_classes(pdf, C)),
            ("composition|mirna", lambda: page_grid(pdf, C, "mirna", "human miRNAs by group",
                                        "share of reads; marker sets; fractions with enough miRNA reads: %s" % (
                                            ", ".join(C.mirna_types()) or "none"), [
                ("human_miRNA_fraction_of_trimmed", "human miRNA (% trimmed)", dict(pct=True)),
                ("human_miRNA_reads", "human miRNA reads", dict(log=True))] +
                [(k, k.replace("miRNA_marker_", "marker: ").replace("_", " ") + " (% miRNA)", pct)
                 for k in C.metrics if k.startswith("miRNA_marker_")], types=C.mirna_types())),
            ("composition|mirna_top", lambda: page_mirna_top(pdf, C)),
            ("composition|functions", lambda: page_functions(pdf, C)),
            ("composition|antisense", lambda: page_antisense(pdf, C)),
            ("comparisons|matched", lambda: page_matched(pdf, C)),
            ("comparisons|features", lambda: page_features(pdf, C)),
            ("comparisons|paired", lambda: page_paired(pdf, C)),
            ("comparisons|omics", lambda: page_grid(pdf, C, "omics", "sRNA vs metaG / metaT", "per-sample Spearman ρ", [
                (k, k.replace("rho_vs_", "ρ vs ").replace("_", " "), {}) for k in C.metrics
                if k.startswith("rho_vs_meta") or k == "rho_vs_sRNA-Int_bin"])),
            ("models|models_intro", lambda: page_models_intro(pdf, C)),
            ("models|models_overview", lambda: page_models_overview(pdf, C)),
            ("models|models", lambda: page_models(pdf, C)),
            ("heterogeneity|heterogeneity", lambda: page_heterogeneity(pdf, C, HZ)),
            ("highlights|highlights_crispr", lambda: page_highlights_crispr(pdf, C, H.get("crispr"))),
            ("highlights|highlights_crispr_group", lambda: page_highlights_crispr_group(pdf, C, H.get("crispr_group"))),
            ("highlights|highlights_group", lambda: page_highlights_group(pdf, C, H.get("group_specific"))),
            ("highlights|highlights_novel", lambda: page_highlights_novel(pdf, C, H.get("novel"), H.get("rfam"))),
            ("highlights|highlights_chimera", lambda: page_highlights_chimera(pdf, C, H.get("chim_lib"), H.get("chim"))),
            ("methods|methods", lambda: page_methods(pdf, C)),
            ("methods|methods_text", lambda: R.pages_markdown(pdf, COHORT_NAME, mfile, "methods (%s)" % a.methods_style)
             if mfile else False),
        ]
        meta = pf.infodict()
        meta["Title"] = "%s cohort report — %s" % (R.PIPELINE_NAME, COHORT_NAME)
        done = []
        for name, fn in jobs:
            sk, _, name = name.rpartition("|")
            pdf.section = section_label(sk) if sk else ""
            try:
                if fn() is not False:
                    done.append(name)
            except Exception as e:
                import traceback
                traceback.print_exc()
                sys.stderr.write("[cohort] page %s skipped (%s: %s)\n" % (name, type(e).__name__, e))
                plt.close("all")
    if a.html:
        R.write_html(os.path.join(d, "cohort.report.html"), COHORT_NAME, C.types, pdf.pages, os.path.basename(pdf_path),
                     sections=pdf.sections)
    sys.stderr.write("[cohort] %s (pages: %s)\n" % (pdf_path, ", ".join(done)))


if __name__ == "__main__":
    main()
