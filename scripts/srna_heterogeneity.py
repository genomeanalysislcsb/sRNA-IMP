#!/usr/bin/env python3
"""srna_heterogeneity.py — heterogeneity of the disease groups for the sRNA-IMP cohort report (used by srna_cohort.py).

Complex diseases may contain subtypes: then only some patients differ from the controls, in different directions,
and group means hardly move. Two analyses look for that, each disease group against the reference group (HC) of the
same block (site), per fraction:

  dispersion      PERMDISP-like: distance of every library to its group centroid in the Bray-Curtis PCoA space
                  (positive axes); larger spread of the patients = more heterogeneous profiles (permutation test of
                  the difference in mean distance)
  outlier range   per feature, the share of patients above the 95th (below the 5th) percentile of the controls,
                  against the share of controls there (one-sided Fisher test, BH-FDR per set / comparison /
                  direction) -> features changed in a subset of patients, with the patients; and per library the
                  share of features outside the control range ("outlier load"; controls leave-one-out), compared
                  between patients and controls (Mann-Whitney)
"""
import collections
import math
import random

import numpy as np

DISP_SETS = ("genus", "fn_kegg_pathway", "rfam_family")
OUT_SETS = (("genus", "genus"), ("species", "species"), ("rfam_family", "Rfam family"), ("rna_class", "RNA class"),
            ("mirna", "miRNA"), ("fn_kegg_pathway", "KEGG pathway"), ("fn_kegg_ko", "KEGG KO"))
MIN_N = 8


def comparisons(C, t):
    """[(label, disease group, block, case keys, control keys)] of a fraction; within blocks if a block column"""
    ref = C.a.reference or (C.groups[0] if C.groups else None)
    ks = C.keys(t)
    out = []
    for g in C.groups:
        if g == ref:
            continue
        cases = [k for k in ks if C.group(k[0]) == g]
        if C.a.match_col:
            for b in sorted({C.block_of(k) for k in cases} - {"NA"}):
                ca = [k for k in cases if C.block_of(k) == b]
                co = [k for k in ks if C.group(k[0]) == ref and C.block_of(k) == b]
                if len(ca) >= MIN_N and len(co) >= MIN_N:
                    out.append(("%s vs %s (%s)" % (g, ref, b), g, b, ca, co))
        else:
            co = [k for k in ks if C.group(k[0]) == ref]
            if len(cases) >= MIN_N and len(co) >= MIN_N:
                out.append(("%s vs %s" % (g, ref), g, "", cases, co))
    return out


def _pcoa_full(D):
    n = D.shape[0]
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ (D ** 2) @ J
    w, v = np.linalg.eigh(B)
    keep = w > 1e-10
    return v[:, keep] * np.sqrt(w[keep])


def _centroid_dist(X, lab):
    d = np.zeros(len(lab))
    for g in set(lab):
        idx = [i for i, x in enumerate(lab) if x == g]
        c = X[idx].mean(axis=0)
        d[idx] = np.sqrt(((X[idx] - c) ** 2).sum(axis=1))
    return d


def dispersion(C, SC, perms=999, seed=1):
    rows, dist = [], []
    rng = random.Random(seed)
    for t in C.types:
        for mx in DISP_SETS:
            ks, feats, Cm, Rl = C.rel_matrix(t, mx, min_prev=0.0, min_reads=1)
            if not feats:
                continue
            pos = {k: i for i, k in enumerate(ks)}
            for label, g, b, ca, co in comparisons(C, t):
                idx = [pos[k] for k in ca + co if k in pos]
                lab = [1 if ks[i] in set(ca) else 0 for i in idx]
                if sum(lab) < MIN_N or len(lab) - sum(lab) < MIN_N:
                    continue
                R = Rl[idx]
                R = R / np.where(R.sum(axis=1, keepdims=True) > 0, R.sum(axis=1, keepdims=True), 1)
                X = _pcoa_full(SC.bray_curtis(R))
                d = _centroid_dist(X, lab)
                la = np.array(lab)
                stat = d[la == 1].mean() - d[la == 0].mean()
                ge = 0
                perm = list(lab)
                for _ in range(perms):
                    rng.shuffle(perm)
                    dp = _centroid_dist(X, perm)
                    pa = np.array(perm)
                    if abs(dp[pa == 1].mean() - dp[pa == 0].mean()) >= abs(stat) - 1e-12:
                        ge += 1
                rows.append(dict(type=t, set=mx, comparison=label, group=g, block=b, n_cases=int(la.sum()),
                                 n_controls=int((la == 0).sum()), median_dist_cases=float(np.median(d[la == 1])),
                                 median_dist_controls=float(np.median(d[la == 0])),
                                 ratio=float(d[la == 1].mean() / max(d[la == 0].mean(), 1e-12)),
                                 p=(ge + 1) / (perms + 1)))
                for i, v in zip(idx, d):
                    dist.append(dict(type=t, set=mx, comparison=label, sample=ks[i][0],
                                     role="case" if ks[i] in set(ca) else "control", distance=float(v)))
    q = SC.bh([r["p"] for r in rows]) if rows else []
    for r, qq in zip(rows, q):
        r["q"] = float(qq)
    return rows, dist


def _matrix(C, keys, mx):
    feats = sorted({f for k in keys for f in (C.vec.get(k, {}).get(mx) or {}) if f})
    if not feats:
        return [], None, None
    M = np.array([[(C.vec[k].get(mx) or {}).get(f, 0.0) for f in feats] for k in keys], float)
    tot = M.sum(axis=1, keepdims=True)
    return feats, M, M / np.where(tot > 0, tot, 1)


def outliers(C, SC, HL, name_of=lambda f: "", min_present=3, min_reads=10, max_p=0.05, max_ids=12):
    feat_rows, load_rows, load_sum = [], [], []
    mtypes = set(C.mirna_types()) if hasattr(C, "mirna_types") else set(C.types)
    for t in C.types:
        for label, g, b, ca, co in comparisons(C, t):
            for mx, lab in OUT_SETS:
                if mx == "mirna" and t not in mtypes:
                    continue
                keys = ca + co
                feats, M, Rl = _matrix(C, keys, mx)
                if not feats:
                    continue
                pres = (M >= min_reads).sum(axis=0) >= min_present
                if pres.sum() == 0:
                    continue
                feats = [f for f, k in zip(feats, pres) if k]
                Rl = Rl[:, pres]
                nca = len(ca)
                A, H = Rl[:nca], Rl[nca:]
                hi = np.percentile(H, 95, axis=0)
                lo = np.percentile(H, 5, axis=0)
                up_a, up_h = A > hi, H > hi
                dn_ok = lo > 0
                dn_a, dn_h = (A < lo) & dn_ok, (H < lo) & dn_ok
                for direction, xa, xh, ok in (("up", up_a, up_h, np.ones(len(feats), bool)), ("down", dn_a, dn_h, dn_ok)):
                    tests = []
                    for j in np.where(ok)[0]:
                        ka, kh = int(xa[:, j].sum()), int(xh[:, j].sum())
                        if ka < 2:
                            continue
                        p = HL.fisher_greater(ka, nca - ka, kh, len(co) - kh)
                        tests.append((j, ka, kh, p))
                    if not tests:
                        continue
                    qs = SC.bh([x[3] for x in tests])
                    for (j, ka, kh, p), q in zip(tests, qs):
                        if p < max_p:
                            ids = [ca[i][0] for i in np.where(xa[:, j])[0]]
                            feat_rows.append(dict(type=t, set=lab, comparison=label, group=g, feature=feats[j],
                                                  name=name_of(feats[j]) if mx.startswith("fn_") else "",
                                                  direction=direction, cases_outside=ka, n_cases=nca,
                                                  share_cases=ka / nca, controls_outside=kh, n_controls=len(co),
                                                  share_controls=kh / len(co), p=p, q=float(q),
                                                  tests=len(tests), patients=",".join(ids[:max_ids]) +
                                                  (",…" if len(ids) > max_ids else "")))
                # outlier load per library (controls leave-one-out)
                load_a = ((A > hi) | ((A < lo) & dn_ok)).mean(axis=1)
                load_h = np.zeros(len(co))
                for i in range(len(co)):
                    Ho = np.delete(H, i, axis=0)
                    h1, l1 = np.percentile(Ho, 95, axis=0), np.percentile(Ho, 5, axis=0)
                    load_h[i] = ((H[i] > h1) | ((H[i] < l1) & (l1 > 0))).mean()
                for k, v in zip(ca, load_a):
                    load_rows.append(dict(type=t, set=lab, comparison=label, sample=k[0], role="case", load=float(v)))
                for k, v in zip(co, load_h):
                    load_rows.append(dict(type=t, set=lab, comparison=label, sample=k[0], role="control",
                                          load=float(v)))
                u, p = SC.mannwhitney(load_a, load_h)
                load_sum.append(dict(type=t, set=lab, comparison=label, group=g, features=len(feats),
                                     median_load_cases=float(np.median(load_a)),
                                     median_load_controls=float(np.median(load_h)),
                                     cases_above_control_p90=int((load_a > np.percentile(load_h, 90)).sum()),
                                     n_cases=nca, p=float(p)))
    q = SC.bh([r["p"] for r in load_sum]) if load_sum else []
    for r, qq in zip(load_sum, q):
        r["q"] = float(qq)
    feat_rows.sort(key=lambda r: (r["q"], r["p"]))
    return feat_rows, load_rows, load_sum
