#!/usr/bin/env python3
"""
srna_qc.py — per-library QC metrics, FastQC module status and rule-based QC flags (PASS / WARN / FAIL)

Used by srna_report.py (per sample) and srna_cohort.py (all samples). Rules come from a TSV
(scripts/qc_rules.tsv or QC_RULES), manual decisions from an exclusion TSV (QC_EXCLUDE:
sample, type, level, reason; level FAIL = exclude, WARN = flag, PASS = keep despite failed rules).

Standalone:  srna_qc.py --outdir imp [--rules qc_rules.tsv] [--exclude qc_exclude.tsv] > qc_flags.tsv
"""
import collections
import glob
import math
import os
import re
import sys
import zipfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import srna_report as R            # noqa: E402

LEVELS = {"PASS": 0, "WARN": 1, "FAIL": 2}
OPS = {"<": lambda a, b: a < b, "<=": lambda a, b: a <= b, ">": lambda a, b: a > b, ">=": lambda a, b: a >= b,
       "==": lambda a, b: a == b, "!=": lambda a, b: a != b}
MODULES_DEFAULT = ["qc", "contam", "prep", "split", "human", "nonhuman", "assembly", "novel", "fragments", "bins",
                   "crispr"]


def frac(a, b):
    return a / b if b else float("nan")


def fastqc_status(stats_dir, sample):
    """{'raw_per_base_sequence_quality': 'PASS', ...} from <sample>.<raw|trimmed>_fastqc.zip (summary.txt)"""
    out = {}
    for zp in glob.glob(os.path.join(stats_dir, "fastqc", "*_fastqc.zip")):
        stage = "raw" if ".raw_fastqc" in zp else "trimmed" if ".trimmed_fastqc" in zp else None
        if not stage:
            continue
        try:
            with zipfile.ZipFile(zp) as z:
                name = [n for n in z.namelist() if n.endswith("/summary.txt")]
                if not name:
                    continue
                for ln in z.read(name[0]).decode("utf-8", "replace").splitlines():
                    f = ln.split("\t")
                    if len(f) >= 2:
                        key = re.sub(r"[^a-z0-9]+", "_", f[1].lower()).strip("_")
                        out["%s_%s" % (stage, key)] = f[0].strip()
        except (zipfile.BadZipFile, OSError):
            continue
    return out


_PACKED = {}


def _packed_done(outdir):
    if outdir not in _PACKED:
        d = collections.defaultdict(set)
        p = os.path.join(outdir, "module_status.txt")
        if os.path.exists(p):
            with open(p) as fh:
                for ln in fh:
                    m = re.match(r"DONE\s+(\S+)/log_cluster_jobs/([^/]+)/\.([^.]+)\.done", ln.strip())
                    if m:
                        d[(m.group(1).split("/")[-1], m.group(2))].add(m.group(3))
        _PACKED[outdir] = d
    return _PACKED[outdir]


def library_metrics(d, outdir, s, t, modules=None):
    """key numbers of one library (TypeData d) — shared by report and cohort"""
    m = collections.OrderedDict()
    raw, trimmed = d.rs("raw"), d.rs("trimmed")
    nh = R.num(d.nonhuman.get("input_reads")) or d.rs("nonhuman")
    hum = d.rs("human_total")
    m["raw_reads"] = raw or float("nan")
    m["trimmed_reads"] = trimmed or float("nan")
    m["trimmed_fraction"] = frac(trimmed, raw)
    m["adapter_dimer_fraction"] = frac(R.num(d.qc.get("adapter_dimers_raw")), raw)
    m["distinct_seq_fraction"] = frac(R.num(d.qc.get("distinct_sequences")), trimmed)
    m["human_fraction"] = frac(hum, trimmed)
    m["nonhuman_reads"] = nh or float("nan")
    m["human_miRNA_fraction_of_trimmed"] = frac(R.num(d.human.get("human_miRNA")), trimmed)
    m["human_miRNA_reads"] = R.num(d.human.get("human_miRNA"))
    for k, lab in (("rRNA_reads", "rRNA"), ("tRNA", "tRNA"), ("Rfam_classified", "Rfam"),
                   ("miRNA_like", "miRNA_like"), ("unknown", "unknown")):
        m["nonhuman_%s_fraction" % lab] = frac(R.num(d.nonhuman.get(k)), nh)
    am = d.asm_summary
    m["assembly_mapped_fraction"] = frac(R.num(am.get("mapped_reads")), R.num(am.get("input_reads")))
    nl = len(d.cand["intergenic"]) + len(d.cand["antisense_to_CDS"])
    m["novel_loci"] = nl if nl else float("nan")
    for r in d.integrity:
        m["integrity_" + r["metric"]] = R.num(r["value"], float("nan"))
    if d.cr_arrays:
        rd = [R.num(r.get("reads")) for r in d.cr_arrays]
        m["crispr_arrays"] = len(rd)
        m["crispr_arrays_with_10_reads"] = sum(1 for x in rd if x >= 10)
        m["crispr_reads_per_million_mapped"] = 1e6 * sum(rd) / (R.num(am.get("mapped_reads")) or float("nan"))
    for k, v in R.rna_qc_for(s, t).items():
        if k not in R.RNA_QC_KEYS and "note" not in k.lower():
            m["qc_" + k] = R.num(v, float("nan"))
    # completeness
    base = os.path.join(outdir, s)
    mods = modules or MODULES_DEFAULT
    done = {x for x in mods if os.path.exists(os.path.join(base, "log_cluster_jobs", t, ".%s.done" % x))}
    if not done:                                   # packed summary (srna_pack_summary.sh): module_status.txt
        done = _packed_done(outdir).get((s, t), set())
    missing = [x for x in mods if x not in done]
    m["missing_modules"] = len(missing)
    m["_missing_modules"] = ",".join(missing)
    m["missing_compare"] = 0 if os.path.exists(os.path.join(base, "Analysis", "comparison", "%s.compare.summary.tsv" % s)) else 1
    for k, v in fastqc_status(os.path.join(base, "Stats", t), s).items():
        m["fastqc_" + k] = v
    return m


def load_rules(path=""):
    path = path or os.path.join(os.path.dirname(os.path.abspath(__file__)), "qc_rules.tsv")
    rules = []
    if not os.path.exists(path):
        sys.stderr.write("[qc] rules %s not found — no QC flags\n" % path)
        return rules
    with open(path) as fh:
        lines = [l for l in fh if l.strip() and not l.startswith("#")]
    hdr = lines[0].rstrip("\n").split("\t")
    for l in lines[1:]:
        r = dict(zip(hdr, l.rstrip("\n").split("\t")))
        if r.get("op") not in OPS or r.get("level", "").upper() not in ("WARN", "FAIL"):
            sys.stderr.write("[qc] rule ignored (op / level): %s" % l)
            continue
        r["level"] = r["level"].upper()
        r["types"] = [x.strip() for x in r.get("types", "*").split(",") if x.strip()] or ["*"]
        rules.append(r)
    return rules


def load_exclude(path=""):
    """(sample, type) -> (level, reason); type * = all types of the sample"""
    out = {}
    if path and os.path.exists(path):
        with open(path) as fh:
            for l in fh:
                if not l.strip() or l.startswith("#") or l.lower().startswith("sample\t"):
                    continue
                f = l.rstrip("\n").split("\t")
                if len(f) >= 2:
                    lv = (f[2].strip().upper() if len(f) > 2 and f[2].strip() else "FAIL")
                    out[(f[0].strip(), f[1].strip() or "*")] = (lv if lv in LEVELS else "FAIL",
                                                                f[3].strip() if len(f) > 3 else "manual")
    return out


def evaluate(m, s, t, rules, exclude=None):
    """-> (overall level, [(rule, level, value, op, threshold, description)])"""
    hits = []
    for r in rules:
        if "*" not in r["types"] and t not in r["types"]:
            continue
        v = m.get(r["metric"])
        if v is None or (isinstance(v, float) and math.isnan(v)):
            continue
        thr = r["threshold"]
        try:
            ok = OPS[r["op"]](float(v), float(thr))
        except (TypeError, ValueError):
            ok = OPS[r["op"]](str(v).upper(), str(thr).upper())
        if ok:
            desc = r.get("description", "")
            if r["metric"] == "missing_modules" and m.get("_missing_modules"):
                desc += " (%s)" % m["_missing_modules"]
            hits.append((r["rule"], r["level"], v, r["op"], thr, desc))
    level = max([h[1] for h in hits], key=lambda x: LEVELS[x], default="PASS")
    ex = (exclude or {}).get((s, t)) or (exclude or {}).get((s, "*"))
    if ex:
        hits.append(("manual", ex[0], "", "", "", ex[1]))
        level = ex[0]                              # a manual decision overrides the rules
    return level, hits


LEVEL_COLOR = {"PASS": "#009E73", "WARN": "#E69F00", "FAIL": "#D55E00", "": "#bbbbbb"}


def main():
    import argparse
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--rules", default="")
    ap.add_argument("--exclude", default="")
    ap.add_argument("--rna-qc", default="")
    a = ap.parse_args()
    R.load_rna_qc(a.rna_qc)
    rules, ex = load_rules(a.rules), load_exclude(a.exclude)
    print("sample\ttype\tqc\treasons")
    for sd in sorted(glob.glob(os.path.join(a.outdir, "*", "Analysis"))):
        s = os.path.basename(os.path.dirname(sd))
        for td in sorted(glob.glob(os.path.join(sd, "*"))):
            t = os.path.basename(td)
            if t in ("annotation", "comparison") or not os.path.isdir(td):
                continue
            d = R.TypeData(os.path.join(a.outdir, s), s, t)
            lv, hits = evaluate(library_metrics(d, a.outdir, s, t), s, t, rules, ex)
            print("%s\t%s\t%s\t%s" % (s, t, lv, "; ".join("%s:%s" % (h[0], h[1]) for h in hits)))


if __name__ == "__main__":
    main()
