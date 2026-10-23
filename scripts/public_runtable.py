#!/usr/bin/env python3
"""Fetch ENA run tables for public studies and pair sRNA-seq with metaG runs by sample.

Example (Tarallo et al. 2019, mSystems; sRNA-seq GSE132236, shotgun metaG SRP136711):
    python3 public_runtable.py --srna GSE132236 --metag SRP136711 --out tarallo2019

Writes <out>.srna.tsv, <out>.metag.tsv (all ENA fields), <out>.pairs.tsv (proposed pairing,
check it!) and <out>.download.sh (wget of the FASTQ files from ENA).
Run on a node with internet access (the HPC access node).
"""
import argparse, csv, io, json, os, re, sys, urllib.parse, urllib.request

ENA = "https://www.ebi.ac.uk/ena/portal/api"
FIELDS = ("run_accession,sample_accession,secondary_sample_accession,experiment_accession,study_accession,"
          "secondary_study_accession,sample_alias,sample_title,experiment_title,library_strategy,library_source,"
          "library_layout,instrument_model,read_count,base_count,fastq_ftp,fastq_md5")


def get(url):
    with urllib.request.urlopen(url, timeout=120) as r:
        return r.read().decode()


def study_of(acc):
    """GEO series -> SRA/ENA study (SRP / PRJNA): the 'Series_relation' lines of the GEO record,
    falling back to an E-utilities search"""
    if not acc.upper().startswith("GSE"):
        return acc
    try:
        soft = get("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s&targ=self&form=text&view=brief" % acc)
        m = re.search(r"\b(SRP\d+)", soft) or re.search(r"\b(PRJNA\d+)", soft)
        if m:
            print("%s -> %s (GEO record)" % (acc, m.group(1)), file=sys.stderr)
            return m.group(1)
    except Exception as e:
        print("GEO lookup failed (%s), trying E-utilities" % e, file=sys.stderr)
    for db, term in (("bioproject", acc), ("sra", acc), ("gds", acc + "[ACCN]")):
        q = urllib.parse.urlencode({"db": db, "term": term, "retmax": 5, "retmode": "json"})
        ids = json.loads(get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + q))["esearchresult"]["idlist"]
        if not ids:
            continue
        txt = get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=%s&id=%s" % (db, ",".join(ids)))
        m = re.search(r"\b(SRP\d+)", txt) or re.search(r"\b(PRJNA\d+)", txt)
        if m:
            print("%s -> %s (%s)" % (acc, m.group(1), db), file=sys.stderr)
            return m.group(1)
    sys.exit("could not resolve the SRA study of %s; pass the SRP/PRJNA accession from the GEO page "
             "(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s, 'Relations') with --srna" % (acc, acc))


def runs(acc):
    q = urllib.parse.urlencode({"accession": study_of(acc), "result": "read_run", "fields": FIELDS,
                                "format": "tsv", "download": "false", "limit": 0})
    rows = list(csv.DictReader(io.StringIO(get(ENA + "/filereport?" + q)), delimiter="\t"))
    if not rows:
        sys.exit("ENA returned no runs for %s" % acc)
    return rows


ACC = re.compile(r"^(gsm|srs|srr|srx|sam[nea]|ers|err|erx|drs|drr)\d+$")


def key(r, regex=None):
    """subject key from sample title / alias (accessions such as GSM/SRS/SAMN are ignored):
    --pair-regex if given (first group), otherwise the last digit-bearing token, e.g.
    'sRNA-Seq of stool sample from the healthy subject H4 - VF10' -> 'vf10'"""
    for f in ("sample_title", "experiment_title", "sample_alias"):
        v = (r.get(f) or "").lower()
        if regex:
            m = re.search(regex, v, re.I)
            if m:
                return re.sub(r"[^a-z0-9]", "", m.group(1).lower())
            continue
        toks = [re.sub(r"[^a-z0-9]", "", t) for t in re.split(r"[\s;,:/|_-]+", v)]
        toks = [t for t in toks if re.search(r"\d", t) and not ACC.match(t)]
        if toks:
            return toks[-1]
    return ""


def write(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS.split(","), delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--srna", required=True, help="sRNA-seq study (GSE..., SRP..., PRJNA...)")
    ap.add_argument("--metag", required=True, help="metagenomic study (SRP..., PRJNA...)")
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--pair-regex", default=None, help="regex with one group for the subject ID in the sample "
                    "title / alias of both studies, e.g. '([vg]f\\w+)'")
    ap.add_argument("--from-tables", action="store_true", help="re-pair from existing <out>.srna.tsv and "
                    "<out>.metag.tsv without downloading")
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    if a.from_tables:
        s, g = [list(csv.DictReader(open(a.out + x), delimiter="\t")) for x in (".srna.tsv", ".metag.tsv")]
    else:
        s, g = runs(a.srna), runs(a.metag)
        write(a.out + ".srna.tsv", s)
        write(a.out + ".metag.tsv", g)
    gk = {}
    for r in g:
        gk.setdefault(key(r, a.pair_regex), []).append(r)
    n = 0
    with open(a.out + ".pairs.tsv", "w") as fh:
        fh.write("key\tsrna_run\tsrna_alias\tsrna_title\tmetag_runs\tmetag_alias\tmetag_title\n")
        for r in s:
            k = key(r, a.pair_regex)
            m = gk.get(k, [])
            n += bool(m)
            fh.write("\t".join([k, r["run_accession"], r["sample_alias"], r["sample_title"],
                                ",".join(x["run_accession"] for x in m), m[0]["sample_alias"] if m else "",
                                m[0]["sample_title"] if m else ""]) + "\n")
    with open(a.out + ".download.sh", "w") as fh:
        fh.write("#!/bin/bash\n# FASTQ downloads from ENA (resumable: wget -c)\nset -e\nmkdir -p srna metag\n")
        for sub, rows in (("srna", s), ("metag", g)):
            for r in rows:
                for u in filter(None, (r.get("fastq_ftp") or "").split(";")):
                    fh.write("wget -c -q -P %s https://%s\n" % (sub, u))
    print("sRNA runs: %d · metaG runs: %d (%d studies' samples) · sRNA runs paired by key: %d\n"
          "check %s.pairs.tsv before use; the key is a heuristic on sample alias/title" % (
              len(s), len(g), len({r['sample_accession'] for r in g}), n, a.out), file=sys.stderr)


if __name__ == "__main__":
    main()
