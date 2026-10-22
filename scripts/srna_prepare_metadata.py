#!/usr/bin/env python3
"""
srna_prepare_metadata.py — turn the full study export (REDCap CSV) into the cohort metadata table

Keeps a small set of covariates that matter for sRNA / microbiome comparisons, derives a few
(Bristol score, disease duration, age at onset), normalises values and writes one row per sample:

  sample_id group site sex age bmi bristol stool laxative laxative_type constipation_status probiotics
  diet gluten_free disease_duration_y age_at_onset age_at_diagnosis longterm_treatment levodopa
  immunosuppressant tnf_inhibitor jak_inhibitor antirheumatic ppi statin antidepressant hypoglycemic
  mds_updrs3 cdai haq bact_count spike_mg spike_mt

Empty medication entries count as "not taken" (--meds-missing no, default; NA keeps them unknown).

Not exported on purpose: birth date and other identifying fields, lab barcodes/numbers, free-text
fields (supplements, medication names), single questionnaire items. Add further raw columns with
--extra (name or new=raw). A data dictionary (<out>.dictionary.tsv) documents every column.

Usage: srna_prepare_metadata.py data/2026-06-04_metadata.csv -o metadata.tsv [--extra nmss_total,who5=who5_total_score_points]
"""
import argparse
import csv
import datetime
import re
import sys
from collections import Counter

NA = "NA"


def clean(v):
    v = (v or "").strip()
    return NA if v == "" or v.upper() in ("NA", "N/A", "NA.NO") else v


def yesno(v):
    v = clean(v).lower()
    return "yes" if v in ("yes", "ja", "y", "1", "true") else "no" if v in ("no", "nein", "n", "0", "false") else NA


def num(v, digits=2):
    try:
        x = float(clean(v))
        return ("%.*f" % (digits, x)).rstrip("0").rstrip(".") if digits else str(int(round(x)))
    except ValueError:
        return NA


def date(v):
    try:                                   # strptime: also on Python < 3.7 (no date.fromisoformat)
        return datetime.datetime.strptime(clean(v)[:10], "%Y-%m-%d").date()
    except ValueError:
        return None


def years(d0, d1):
    return (d1 - d0).days / 365.25 if d0 and d1 else None


def bristol(v):
    m = re.search(r"(?:Typ|Type)\s*([1-7])", v or "")
    return m.group(1) if m else NA


# (output column, description, function(row) -> value)
def build_columns(args):
    med = (lambda v: "no" if yesno(v) == NA else yesno(v)) if args.meds_missing == "no" else yesno
    mnote = " (empty entry = not taken)" if args.meds_missing == "no" else ""

    def dur(r, key):
        y = years(date(r.get(key)), date(r.get("lab_date")))
        return "%.1f" % y if y is not None and y >= 0 else NA

    def disease_date(r, pd_key, ra_key):
        return pd_key if r.get("diagnosis") == "PD" else ra_key if r.get("diagnosis") == "RA" else None

    def duration(r):
        k = disease_date(r, "pd_diagnosis", "ra_diagnosis")
        return dur(r, k) if k else NA

    def age_minus(r, k):
        if not k:
            return NA
        y = years(date(r.get(k)), date(r.get("lab_date")))
        try:
            return "%.1f" % (float(r["age"]) - y) if y is not None and y >= 0 else NA
        except (ValueError, KeyError):
            return NA

    def levodopa(r):
        v = yesno(r.get("Levodopa") or r.get("levodopa"))
        if v == NA and (r.get("diagnosis") != "PD" or args.meds_missing == "no"):
            return "no"
        return v

    cols = [
        ("sample_id", "record_id", lambda r: clean(r.get("record_id"))),
        ("group", "diagnosis (HC / PD / RA)", lambda r: clean(r.get("diagnosis"))),
        ("site", "clinical_site", lambda r: clean(r.get("clinical_site") or r.get("site"))),
        ("sex", "sex", lambda r: clean(r.get("sex")).lower()),
        ("age", "age at the visit (years)", lambda r: num(r.get("age"), 0)),
        ("bmi", "body-mass index", lambda r: num(r.get("bmi"), 1)),
        ("bristol", "Bristol stool type 1-7 (from bristol_stool_type)", lambda r: bristol(r.get("bristol_stool_type"))),
        ("stool", "bss: Normal / Constipation / Diarrhea", lambda r: clean(r.get("bss"))),
        ("laxative", "laxative use", lambda r: yesno(r.get("laxative"))),
        ("laxative_type", "laxative.simplified (No / Flax seed / Psyllium seed / PEG)", lambda r: clean(r.get("laxative.simplified"))),
        ("constipation_status", "laxative.bss (e.g. untreated / resolved / non-resolved constipation)", lambda r: clean(r.get("laxative.bss"))),
        ("probiotics", "probiotics", lambda r: yesno(r.get("probiotics"))),
        ("diet", "Omnivore / Vegetarian / Pescetarian / Vegan", lambda r: clean(r.get("diet"))),
        ("gluten_free", "gluten-free diet", lambda r: yesno(r.get("diet_glutenfree"))),
        ("disease_duration_y", "years from diagnosis (pd_diagnosis / ra_diagnosis) to the visit (lab_date)", duration),
        ("age_at_onset", "age minus years since symptom start (pd_start_date / ra_start)",
         lambda r: age_minus(r, disease_date(r, "pd_start_date", "ra_start"))),
        ("age_at_diagnosis", "age minus years since diagnosis", lambda r: age_minus(r, disease_date(r, "pd_diagnosis", "ra_diagnosis"))),
        ("longterm_treatment", "No longterm / PD_treatment / RA_treatment", lambda r: clean(r.get("longterm_treatment"))),
        ("levodopa", "levodopa (empty: no, except PD with --meds-missing NA)", levodopa),
        ("immunosuppressant", "medication class" + mnote, lambda r: med(r.get("Immunosuppressant"))),
        ("tnf_inhibitor", "medication class" + mnote, lambda r: med(r.get("TNF Inhibitor"))),
        ("jak_inhibitor", "medication class" + mnote, lambda r: med(r.get("JAK Inhibitor"))),
        ("antirheumatic", "medication class" + mnote, lambda r: med(r.get("Antirheumatic"))),
        ("ppi", "proton pump inhibitor" + mnote, lambda r: med(r.get("Proton Pump Inhibitor"))),
        ("statin", "medication class" + mnote, lambda r: med(r.get("Statin"))),
        ("antidepressant", "medication class" + mnote, lambda r: med(r.get("Antidepressant"))),
        ("hypoglycemic", "medication class" + mnote, lambda r: med(r.get("Hypoglycemic"))),
        ("mds_updrs3", "MDS-UPDRS part III total (PD motor)", lambda r: num(r.get("mds_3_total_score"), 0)),
        ("cdai", "clinical disease activity index (RA)", lambda r: num(r.get("js_cdai"), 1)),
        ("haq", "HAQ disability score", lambda r: num(r.get("haq_score"), 3)),
        ("bact_count", "bact_count (technical: microbial load)", lambda r: num(r.get("bact_count"), 3)),
        ("spike_mg", "spike_mg (technical)", lambda r: num(r.get("spike_mg"), 3)),
        ("spike_mt", "spike_mt (technical)", lambda r: num(r.get("spike_mt"), 3)),
    ]
    for e in [x.strip() for x in (args.extra or "").split(",") if x.strip()]:
        new, raw = (e.split("=", 1) + [None])[:2] if "=" in e else (e, e)
        cols.append((new, "raw column %s" % raw, lambda r, raw=raw: clean(r.get(raw))))
    return cols


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("export", help="study export (CSV or TSV)")
    ap.add_argument("-o", "--out", required=True, help="output metadata TSV")
    ap.add_argument("--extra", default="", help="additional raw columns: name or new_name=raw_name, comma-separated")
    ap.add_argument("--meds-missing", default="no", choices=("no", "NA"),
                    help="empty medication entries: no = not taken (default), NA = unknown")
    ap.add_argument("--event", default="", help="keep only rows whose redcap_event_name contains this (e.g. t0)")
    a = ap.parse_args()

    with open(a.export, encoding="utf-8-sig", errors="replace", newline="") as fh:
        text = fh.read()
    delim = "\t" if text.count("\t") > text.count(",") else ";" if text.count(";") > text.count(",") else ","
    rows = list(csv.DictReader(text.splitlines(), delimiter=delim))
    if a.event:
        rows = [r for r in rows if a.event in r.get("redcap_event_name", "")]
    hdr = rows[0].keys() if rows else []
    missing = [c for c in ("record_id", "diagnosis") if c not in hdr]
    if missing:
        sys.exit("[metadata] required column(s) missing: %s" % ", ".join(missing))
    for e in [x.split("=")[-1].strip() for x in a.extra.split(",") if x.strip()]:
        if e not in hdr:
            sys.exit("[metadata] --extra column not in the export: %s" % e)

    cols = build_columns(a)
    out_rows, seen = [], Counter()
    for r in rows:
        sid = clean(r.get("record_id"))
        if sid == NA:
            continue
        seen[sid] += 1
        out_rows.append([fn(r) for _, _, fn in cols])
    dups = [k for k, v in seen.items() if v > 1]
    if dups:
        sys.stderr.write("[metadata] WARNING: %d sample IDs occur more than once (events / repeats?): %s — "
                         "use --event\n" % (len(dups), ", ".join(dups[:10])))
    with open(a.out, "w") as out:
        out.write("\t".join(c for c, _, _ in cols) + "\n")
        for r in out_rows:
            out.write("\t".join(r) + "\n")
    dpath = re.sub(r"\.(tsv|txt|csv)$", "", a.out) + ".dictionary.tsv"
    with open(dpath, "w") as out:
        out.write("column\tsource / meaning\tnon_missing\n")
        for i, (c, desc, _) in enumerate(cols):
            out.write("%s\t%s\t%d\n" % (c, desc, sum(1 for r in out_rows if r[i] != NA)))
    # summary without individual data
    gi = 1
    g = Counter(r[gi] for r in out_rows)
    sys.stderr.write("[metadata] %d samples -> %s (%s); dictionary: %s\n" % (
        len(out_rows), a.out, ", ".join("%s %d" % kv for kv in sorted(g.items())), dpath))
    site = Counter((r[1], r[2]) for r in out_rows)
    sys.stderr.write("[metadata] group x site: %s\n" % ", ".join("%s/%s %d" % (k[0], k[1], v) for k, v in sorted(site.items())))


if __name__ == "__main__":
    main()
