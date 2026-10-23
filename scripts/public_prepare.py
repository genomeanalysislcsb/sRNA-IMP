#!/usr/bin/env python3
"""Prepare a public sRNA-seq + metaG dataset for IMP3 (metaG-only) and sRNA-IMP.

Input: the <prefix>.srna.tsv, <prefix>.metag.tsv and <prefix>.pairs.tsv written by public_runtable.py.
Output in --outdir:
  ids.txt              one sample ID per line (the subject code, e.g. VF10)
  metadata.tsv         sample_id, group, subject, sRNA / metaG run, read counts (for --cohort)
  download.sh          ENA downloads of the paired runs only, named by sample ID:
                         srna/<ID>.fastq.gz, metag/<ID>_R1.fastq.gz, metag/<ID>_R2.fastq.gz (md5-checked)
  imp3_samples.tsv     sample ID, metaG R1, metaG R2 (to build the IMP3 metaG-only configs)
  sample.template.tsv  sRNA-IMP sample template; expand with
                         srna_make_samples.sh -c -o samples.tsv sample.template.tsv ids.txt

Example (Tarallo et al. 2019):
  python3 public_prepare.py --prefix tarallo2019/raw --outdir tarallo2019 \\
      --base /work/projects/expobiome/public/tarallo2019 --type sRNA
"""
import argparse, csv, os, re, sys

GROUPS = ((r"healthy", "HC"), (r"adenoma", "Adenoma"), (r"colorectal cancer|carcinoma|\bcrc\b", "CRC"))


def read(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def group_of(title):
    for rx, g in GROUPS:
        if re.search(rx, title, re.I):
            return g
    return "NA"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--prefix", required=True, help="prefix used with public_runtable.py --out")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--base", required=True, help="absolute data directory on the HPC (fastq/ and IMP/ go below it)")
    ap.add_argument("--type", default="sRNA", help="sRNA library type label in the sample sheet (default sRNA)")
    ap.add_argument("--imp-assembly", default="Assembly/mg.assembly.merged.fa.gz",
                    help="assembly path inside an IMP3 sample dir (check after the IMP3 run)")
    ap.add_argument("--imp-gff", default="Analysis/annotation/annotation_CDS_RNA_hmms.gff.gz")
    ap.add_argument("--imp-bins", default="Binning/MetaWrap/bins.contigs.gz")
    ap.add_argument("--imp-gtdb", default="Stats/all_bin_stats.tsv.gz")
    ap.add_argument("--imp-mg-contigs", default="Stats/mg/mg.assembly.contig_read_count.txt.gz")
    a = ap.parse_args()

    srna = {r["run_accession"]: r for r in read(a.prefix + ".srna.tsv")}
    metag = {r["run_accession"]: r for r in read(a.prefix + ".metag.tsv")}
    pairs = [p for p in read(a.prefix + ".pairs.tsv") if p["metag_runs"]]
    os.makedirs(a.outdir, exist_ok=True)

    rows, seen, problems = [], set(), []
    for p in pairs:
        s = srna[p["srna_run"]]
        mruns = p["metag_runs"].split(",")
        sid = (p["metag_alias"] or p["key"]).strip()
        sid = re.sub(r"[^A-Za-z0-9_.-]", "_", sid)
        if sid in seen:
            problems.append("duplicate sample ID %s" % sid)
            continue
        if len(mruns) != 1:
            problems.append("%s: %d metaG runs (%s), only the first is used" % (sid, len(mruns), p["metag_runs"]))
        m = metag[mruns[0]]
        if len([u for u in m["fastq_ftp"].split(";") if u]) != 2:
            problems.append("%s: metaG run %s is not paired-end (%s)" % (sid, m["run_accession"], m["fastq_ftp"]))
        subj = re.search(r"subject\s+(\S+)|patient\s+(\S+)", s["sample_title"])
        rows.append(dict(sample_id=sid, group=group_of(s["sample_title"]),
                         subject=next(g for g in subj.groups() if g) if subj else "",
                         srna_run=s["run_accession"], srna_reads=s["read_count"],
                         metag_run=m["run_accession"], metag_reads=m["read_count"], srna_title=s["sample_title"],
                         _srna=s, _metag=m))
        seen.add(sid)

    with open(os.path.join(a.outdir, "ids.txt"), "w") as fh:
        fh.writelines(r["sample_id"] + "\n" for r in rows)
    cols = ["sample_id", "group", "subject", "srna_run", "srna_reads", "metag_run", "metag_reads", "srna_title"]
    with open(os.path.join(a.outdir, "metadata.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        fh.writelines("\t".join(str(r[c]) for c in cols) + "\n" for r in rows)

    fq = os.path.join(a.base, "fastq")
    with open(os.path.join(a.outdir, "download.sh"), "w") as fh:
        fh.write("#!/bin/bash\n# ENA downloads of the paired runs (resumable; md5 checked)\nset -euo pipefail\n"
                 "mkdir -p %s/srna %s/metag\n"
                 "get() { wget -c -q -O \"$2\" \"https://$1\"; echo \"$3  $2\" | md5sum -c --quiet - || "
                 "{ echo \"md5 mismatch: $2\" >&2; exit 1; }; }\n" % (fq, fq))
        for r in rows:
            s, m = r["_srna"], r["_metag"]
            su, smd = s["fastq_ftp"].split(";")[0], s["fastq_md5"].split(";")[0]
            fh.write("get %s %s/srna/%s.fastq.gz %s\n" % (su, fq, r["sample_id"], smd))
            for i, (u, md) in enumerate(zip(m["fastq_ftp"].split(";"), m["fastq_md5"].split(";")), 1):
                fh.write("get %s %s/metag/%s_R%d.fastq.gz %s\n" % (u, fq, r["sample_id"], i, md))
    os.chmod(os.path.join(a.outdir, "download.sh"), 0o755)

    with open(os.path.join(a.outdir, "imp3_samples.tsv"), "w") as fh:
        fh.write("sample_id\tmetaG_R1\tmetaG_R2\toutdir\n")
        for r in rows:
            fh.write("%s\t%s/metag/%s_R1.fastq.gz\t%s/metag/%s_R2.fastq.gz\t%s/IMP/%s\n" % (
                r["sample_id"], fq, r["sample_id"], fq, r["sample_id"], a.base, r["sample_id"]))

    imp = os.path.join(a.base, "IMP", "<ID>")
    j = lambda x: os.path.join(imp, x)
    with open(os.path.join(a.outdir, "sample.template.tsv"), "w") as fh:
        fh.write("# sRNA-IMP sample template for a public dataset (one sRNA library per stool sample, metaG-only IMP3)\n"
                 "# expand: srna_make_samples.sh -c -o samples.tsv sample.template.tsv ids.txt\n"
                 "# check the IMP3 file names (--imp-*) against one finished IMP3 sample first\n")
        fh.write("\t".join(["sample_id", "type", "raw_fastq", "assembly", "gff", "contig2bin", "gtdbtk",
                            "metaT_features", "metaT_contigs", "metaG_features", "metaG_contigs"]) + "\n")
        fh.write("\t".join(["<ID>", a.type, os.path.join(fq, "srna", "<ID>.fastq.gz"), j(a.imp_assembly),
                            j(a.imp_gff), j(a.imp_bins), j(a.imp_gtdb), "-", "-", j("Analysis/annotation"),
                            j(a.imp_mg_contigs)]) + "\n")

    n = {}
    for r in rows:
        n[r["group"]] = n.get(r["group"], 0) + 1
    print("%d samples (%s) -> %s" % (len(rows), ", ".join("%s %d" % kv for kv in sorted(n.items())), a.outdir),
          file=sys.stderr)
    for p in problems:
        print("warning: " + p, file=sys.stderr)


if __name__ == "__main__":
    main()
