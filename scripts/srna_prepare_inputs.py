#!/usr/bin/env python3
"""srna_prepare_inputs.py — turn assembly, annotation, binning and read-mapping outputs of any workflow into the
per-sample inputs of sRNA-IMP (sample sheet columns 4-11). IMP3 outputs need no conversion; use this for
MEGAHIT/metaSPAdes + Prokka/Bakta + MetaBAT2/MetaWRAP + CheckM2 + GTDB-Tk, or similar.

Usage:
  srna_prepare_inputs.py --sample S --outdir prep --contigs contigs.fa[.gz] --gff annotation.gff3[.gz]
      [--bins-dir DIR | --contig2bin TSV] [--gtdbtk bac120.summary.tsv[,ar53.summary.tsv]] [--checkm quality_report.tsv]
      [--metag-bam mg.bam] [--metat-bam mt.bam --metat-strand reverse|forward|unstranded] [--threads 8]
      [--srna-fastq S.fastq.gz --type sRNA]

Writes prep/S/:
  assembly.fa.gz       the contigs (names unchanged; must match the GFF and the BAMs)
  annotation.gff.gz    GFF3 without a ##FASTA section, every feature with an ID; Prokka CRISPR repeat_region -> CRISPR
  contig2bin.tsv.gz    contig<TAB>bin (from --bins-dir: one FASTA per bin, bin = file name without extension)
  bins.tsv             user_genome, classification (GTDB-Tk), Completeness, Contamination (CheckM/CheckM2)
  counts/<mt|mg>.ID.<type>_counts.tsv   featureCounts per feature type, keyed by the GFF ID
  counts/<mt|mg>/contig_counts.tsv      featureCounts per contig (all reads, incl. intergenic)
  sample_row.tsv       the sample-sheet row (with --srna-fastq; else with <raw_fastq> to fill in)
Counting as in IMP3: -M --fraction (multi-mappers split), paired and unpaired reads in one BAM, metaT with the
library's strandedness (dUTP / TruSeq Stranded = reverse = -s 2), metaG unstranded. Needs featureCounts (subread)
for the counts; the rest is pure Python. For the metaT strand check of --compare also set COMPARE_MT_BAM to the
metaT BAM (sorted and indexed)."""
import argparse
import glob
import gzip
import os
import re
import shutil
import subprocess
import sys
import tempfile

TYPES = ["CDS", "tRNA", "rRNA", "ncRNA", "tmRNA", "CRISPR", "crispr-repeat", "crispr-spacer", "regulatory_region",
         "oriC"]
FASTA_EXT = re.compile(r"\.(fa|fna|fasta|fsa)(\.gz)?$")


def op(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def contigs_and_lengths(path, out):
    lens, name, n = {}, None, 0
    with op(path) as fh, gzip.open(out, "wt") as fo:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                lens[name] = 0
                fo.write(">%s\n" % name)
            else:
                s = line.strip()
                lens[name] += len(s)
                fo.write(s + "\n")
    return lens


def clean_gff(path, out, contigs):
    """drop ##FASTA, give every feature an ID, map Prokka CRISPR repeat_region -> CRISPR, check contig names"""
    n_type, missing, new_ids = {}, set(), 0
    with op(path) as fh, gzip.open(out, "wt") as fo:
        fo.write("##gff-version 3\n")
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9:
                continue
            if f[0] not in contigs:
                missing.add(f[0])
                continue
            if f[2] == "repeat_region" and re.search(r"rpt_family=CRISPR", f[8], re.I):
                f[2] = "CRISPR"
            if f[6] not in "+-.":
                f[6] = "."
            if not re.search(r"(?:^|;)ID=", f[8]):
                new_ids += 1
                f[8] = "ID=feature_%d;%s" % (new_ids, f[8].strip(";")) if f[8] not in ("", ".") else "ID=feature_%d" % new_ids
            n_type[f[2]] = n_type.get(f[2], 0) + 1
            fo.write("\t".join(f) + "\n")
    if missing:
        print("warning: %d GFF contig names not in the assembly (e.g. %s) — features dropped" % (
            len(missing), ", ".join(sorted(missing)[:3])), file=sys.stderr)
    return n_type, new_ids


def features_by_type(gff):
    by = {}
    with op(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) == 9 and f[2] in TYPES:
                fid = re.search(r"(?:^|;)ID=([^;]+)", f[8]).group(1)
                by.setdefault(f[2], []).append((fid, f[0], f[3], f[4], f[6] if f[6] in "+-" else "+"))
    return by


def contig2bin(bins_dir, table, out):
    rows = []
    if bins_dir:
        for fa in sorted(glob.glob(os.path.join(bins_dir, "*"))):
            if not FASTA_EXT.search(fa):
                continue
            b = FASTA_EXT.sub("", os.path.basename(fa))
            with op(fa) as fh:
                rows += [(line[1:].split()[0], b) for line in fh if line.startswith(">")]
    elif table:
        with op(table) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 2 and not line.startswith("#"):
                    rows.append((f[0], FASTA_EXT.sub("", f[1])))
    if rows:
        with gzip.open(out, "wt") as fo:
            fo.writelines("%s\t%s\n" % r for r in rows)
    return rows


def bin_table(gtdbtk, checkm, out):
    tax, qual = {}, {}
    for p in [x for x in (gtdbtk or "").split(",") if x]:
        with op(p) as fh:
            head = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                r = dict(zip(head, line.rstrip("\n").split("\t")))
                tax[FASTA_EXT.sub("", r.get("user_genome", ""))] = r.get("classification", "")
    if checkm:
        with op(checkm) as fh:
            head = [h.strip().lower() for h in fh.readline().rstrip("\n").split("\t")]
            for line in fh:
                f = line.rstrip("\n").split("\t")
                r = dict(zip(head, f))
                b = FASTA_EXT.sub("", r.get("name", r.get("bin id", f[0])))
                qual[b] = (r.get("completeness", ""), r.get("contamination", ""))
    bins = sorted(set(tax) | set(qual))
    if bins:
        with open(out, "w") as fo:
            fo.write("user_genome\tclassification\tCompleteness\tContamination\n")
            for b in bins:
                c, k = qual.get(b, ("", ""))
                fo.write("%s\t%s\t%s\t%s\n" % (b, tax.get(b, ""), c, k))
    return bins


def featurecounts(fc, bam, saf, strand, threads, out):
    paired = False
    try:
        p = subprocess.Popen(["samtools", "view", "-f", "1", bam], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                             text=True)
        paired = bool(p.stdout.readline())
        p.kill()
    except OSError:
        paired = True
    cmd = [fc, "-T", str(threads), "-M", "--fraction", "-s", str(strand), "-F", "SAF", "-a", saf, "-o", out, bam]
    if paired:
        cmd[1:1] = ["-p", "--countReadPairs"]
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0 and paired:          # featureCounts < 2.0.2 has no --countReadPairs
        cmd.remove("--countReadPairs")
        r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        sys.exit("Error: featureCounts failed on %s:\n%s" % (bam, r.stderr[-1500:]))
    if os.path.exists(out + ".summary"):
        os.remove(out + ".summary")


def count(omic, bam, strand, by_type, lens, cdir, fc, threads):
    os.makedirs(os.path.join(cdir, omic), exist_ok=True)
    with tempfile.TemporaryDirectory(dir=cdir) as tmp:
        for t, feats in by_type.items():
            saf = os.path.join(tmp, t + ".saf")
            with open(saf, "w") as fo:
                fo.write("GeneID\tChr\tStart\tEnd\tStrand\n")
                fo.writelines("%s\t%s\t%s\t%s\t%s\n" % x for x in feats)
            featurecounts(fc, bam, saf, strand, threads, os.path.join(cdir, "%s.ID.%s_counts.tsv" % (omic, t)))
        saf = os.path.join(tmp, "contigs.saf")
        with open(saf, "w") as fo:
            fo.write("GeneID\tChr\tStart\tEnd\tStrand\n")
            fo.writelines("%s\t%s\t1\t%d\t+\n" % (c, c, L) for c, L in lens.items() if L > 0)
        featurecounts(fc, bam, saf, 0, threads, os.path.join(cdir, omic, "contig_counts.tsv"))
    print("%s: %d feature types and contig counts from %s" % (omic, len(by_type), bam), file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True); ap.add_argument("--outdir", required=True)
    ap.add_argument("--contigs", required=True); ap.add_argument("--gff", required=True)
    ap.add_argument("--bins-dir", default=""); ap.add_argument("--contig2bin", default="")
    ap.add_argument("--gtdbtk", default=""); ap.add_argument("--checkm", default="")
    ap.add_argument("--metag-bam", default=""); ap.add_argument("--metat-bam", default="")
    ap.add_argument("--metat-strand", default="reverse", choices=["reverse", "forward", "unstranded"])
    ap.add_argument("--featurecounts", default="featureCounts"); ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--srna-fastq", default=""); ap.add_argument("--type", default="sRNA")
    a = ap.parse_args()

    d = os.path.abspath(os.path.join(a.outdir, a.sample))
    os.makedirs(d, exist_ok=True)
    asm, gff = os.path.join(d, "assembly.fa.gz"), os.path.join(d, "annotation.gff.gz")
    lens = contigs_and_lengths(a.contigs, asm)
    n_type, new_ids = clean_gff(a.gff, gff, lens)
    print("%s: %d contigs; GFF types %s%s" % (a.sample, len(lens), ", ".join("%s %d" % kv for kv in sorted(n_type.items())),
          "; %d IDs added" % new_ids if new_ids else ""), file=sys.stderr)
    for t, why in (("CDS", "no gene-level counts or antisense analysis"),
                   ("tRNA", "tRNAs and tmRNAs are then predicted with ARAGORN in the novel module"),
                   ("rRNA", "rRNA reads are then recognised by RiboDetector / Rfam only, and missing rRNA genes "
                            "are recovered by the assembly module")):
        if t not in n_type:
            print("warning: no %s features in the GFF: %s" % (t, why), file=sys.stderr)
    c2b = contig2bin(a.bins_dir, a.contig2bin, os.path.join(d, "contig2bin.tsv.gz"))
    bins = bin_table(a.gtdbtk, a.checkm, os.path.join(d, "bins.tsv"))
    by_type = features_by_type(gff)
    cdir = os.path.join(d, "counts")
    if a.metat_bam or a.metag_bam:
        if not shutil.which(a.featurecounts):
            sys.exit("Error: %s not found (conda install -c bioconda subread)" % a.featurecounts)
        os.makedirs(cdir, exist_ok=True)
    if a.metat_bam:
        count("mt", a.metat_bam, {"reverse": 2, "forward": 1, "unstranded": 0}[a.metat_strand], by_type, lens, cdir,
              a.featurecounts, a.threads)
    if a.metag_bam:
        count("mg", a.metag_bam, 0, by_type, lens, cdir, a.featurecounts, a.threads)
    has = os.path.isdir(cdir)
    row = [a.sample, a.type, os.path.abspath(a.srna_fastq) if a.srna_fastq else "<raw_fastq>", asm, gff,
           os.path.join(d, "contig2bin.tsv.gz") if c2b else "-", os.path.join(d, "bins.tsv") if bins else "-",
           cdir if has else "-", cdir if has else "-", "-", "-"]
    with open(os.path.join(d, "sample_row.tsv"), "w") as fo:
        fo.write("\t".join(row) + "\n")
    print("%s: %d binned contigs in %d bins; sample-sheet row in %s" % (
        a.sample, len(c2b), len(set(b for _, b in c2b)), os.path.join(d, "sample_row.tsv")), file=sys.stderr)


if __name__ == "__main__":
    main()
