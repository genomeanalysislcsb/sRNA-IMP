#!/usr/bin/env python3
"""srna_fetch_mgnify.py — pick species-representative genomes from an MGnify genome catalogue and download them, as
reference genomes for sRNA-IMP (srna_prepare_inputs.py --genomes) when a sample has no metagenome or
metatranscriptome.

MGnify catalogues (https://www.ebi.ac.uk/metagenomics/browse/genomes) cover host-associated biomes such as the human
gut, oral cavity and vagina, mouse gut, pig and chicken gut, cow rumen and others. Each has a metadata table
(genomes-all_metadata.tsv) with one row per genome: quality, species representative, GTDB lineage and a download URL
of the genome's GFF (Prokka annotation with the sequence in a ##FASTA section).

Usage:
  srna_fetch_mgnify.py --catalogue human-gut --version v2.0.2 --outdir genomes  SELECTION
  srna_fetch_mgnify.py --metadata genomes-all_metadata.tsv --outdir genomes      SELECTION

SELECTION (combinable; species representatives only unless --all-genomes):
  --genomes MGYG000000001,...       accessions (or --genome-list FILE, one per line)
  --species FILE                    species names, one per line ("Agathobacter rectalis" or "s__Agathobacter rectalis")
  --bracken FILE[,FILE...]          Bracken species tables of sRNA-IMP (nonhuman/rrna_tax/<s>.rRNA.bracken.S.tsv):
                                    species with >= --min-reads estimated rRNA reads in any of them
  --taxa g__Bacteroides,f__...      every representative whose lineage contains one of these ranks

Writes to --outdir: <accession>.gff[.gz] per genome (skipped if present), genome_taxonomy.tsv (accession, lineage),
genome_quality.tsv (CheckM-style Name, Completeness, Contamination), selection.tsv (why each genome was chosen) and
unmatched_species.txt (requested species without a representative; names can differ between GTDB releases).
Then:
  srna_prepare_inputs.py --sample REF --outdir prep --genomes genomes \\
      --genome-taxonomy genomes/genome_taxonomy.tsv --checkm genomes/genome_quality.tsv
and ASSEMBLY_SOURCE=reference in the configuration.
"""
import argparse
import csv
import gzip
import io
import os
import re
import shutil
import sys
import time
import urllib.request

BASE = "https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes"


def norm(name):
    """species name without rank prefix, single spaces, lower case"""
    name = re.sub(r"^s__", "", name.strip())
    return re.sub(r"\s+", " ", name).lower()


def species_of(lineage):
    last = lineage.split(";")[-1].strip()
    return norm(last) if last.startswith("s__") and len(last) > 3 else ""


def fetch(url, path, tries=3):
    tmp = path + ".part"
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r, open(tmp, "wb") as fo:
                shutil.copyfileobj(r, fo)
            os.replace(tmp, path)
            return True
        except Exception as e:                 # noqa: BLE001 — report and retry any network error
            err = e
            time.sleep(2 * (i + 1))
    sys.stderr.write("warning: could not download %s (%s)\n" % (url, err))
    if os.path.exists(tmp):
        os.remove(tmp)
    return False


def read_metadata(src):
    if re.match(r"https?://|ftp://", src):
        with urllib.request.urlopen(src, timeout=300) as r:
            text = r.read().decode()
    else:
        with (gzip.open(src, "rt") if src.endswith(".gz") else open(src)) as fh:
            text = fh.read()
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    if not rows:
        sys.exit("Error: empty metadata table %s" % src)
    cols = {c.lower(): c for c in rows[0]}

    def col(*names):
        for n in names:
            if n.lower() in cols:
                return cols[n.lower()]
        return None
    need = {"genome": col("Genome", "Genome_accession"), "rep": col("Species_rep"), "lineage": col("Lineage"),
            "url": col("FTP_download", "Download", "GFF_download")}
    miss = [k for k in ("genome", "lineage") if not need[k]]
    if miss:
        sys.exit("Error: metadata table %s lacks columns for %s (has: %s)" % (src, ", ".join(miss), ", ".join(rows[0])))
    comp, cont = col("Completeness"), col("Contamination")
    out = []
    for r in rows:
        out.append({"genome": r[need["genome"]], "rep": r.get(need["rep"], "") if need["rep"] else r[need["genome"]],
                    "lineage": r[need["lineage"]], "url": r.get(need["url"], "") if need["url"] else "",
                    "completeness": r.get(comp, "") if comp else "", "contamination": r.get(cont, "") if cont else ""})
    return out


def bracken_species(paths, min_reads):
    sp = {}
    for p in paths:
        with open(p) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                try:
                    n = float(r.get("new_est_reads", r.get("kraken_assigned_reads", 0)) or 0)
                except ValueError:
                    continue
                if r.get("taxonomy_lvl", "S") == "S" and n >= min_reads:
                    k = norm(r["name"])
                    sp[k] = max(sp.get(k, 0), n)
    return sp


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--catalogue", default="", help="MGnify biome, e.g. human-gut, human-oral, human-vaginal, mouse-gut")
    ap.add_argument("--version", default="", help="catalogue version, e.g. v2.0.2 (see the MGnify genomes FTP)")
    ap.add_argument("--metadata", default="", help="local genomes-all_metadata.tsv (instead of --catalogue/--version)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--genomes", default=""); ap.add_argument("--genome-list", default="")
    ap.add_argument("--species", default=""); ap.add_argument("--taxa", default="")
    ap.add_argument("--bracken", default=""); ap.add_argument("--min-reads", type=float, default=10)
    ap.add_argument("--all-genomes", action="store_true", help="allow genomes that are not species representatives")
    ap.add_argument("--max-genomes", type=int, default=500, help="stop if the selection is larger  [500]")
    ap.add_argument("--dry-run", action="store_true", help="write the tables, download nothing")
    a = ap.parse_args()

    if a.metadata:
        src = a.metadata
    elif a.catalogue and a.version:
        src = "%s/%s/%s/genomes-all_metadata.tsv" % (BASE, a.catalogue, a.version)
    else:
        ap.error("give --metadata, or --catalogue and --version")
    meta = read_metadata(src)
    pool = meta if a.all_genomes else [m for m in meta if m["genome"] == m["rep"]]
    by_acc = {m["genome"]: m for m in meta}
    by_sp = {}
    for m in pool:
        s = species_of(m["lineage"])
        if s:
            by_sp.setdefault(s, m)

    chosen, why, unmatched = {}, {}, []
    accs = [x for x in a.genomes.split(",") if x]
    if a.genome_list:
        accs += [l.split()[0] for l in open(a.genome_list) if l.strip() and not l.startswith("#")]
    for acc in accs:
        if acc in by_acc:
            chosen[acc] = by_acc[acc]; why.setdefault(acc, []).append("accession")
        else:
            unmatched.append(acc)
    want = {}
    if a.species:
        want.update({norm(l): None for l in open(a.species) if l.strip() and not l.startswith("#")})
    if a.bracken:
        want.update(bracken_species([p for p in a.bracken.split(",") if p], a.min_reads))
    for s, n in want.items():
        m = by_sp.get(s)
        if m:
            chosen[m["genome"]] = m
            why.setdefault(m["genome"], []).append("species %s%s" % (s, "" if n is None else " (%g rRNA reads)" % n))
        else:
            unmatched.append(s)
    for t in [x.strip() for x in a.taxa.split(",") if x.strip()]:
        for m in pool:
            if t in [x.strip() for x in m["lineage"].split(";")]:
                chosen[m["genome"]] = m; why.setdefault(m["genome"], []).append("taxon %s" % t)
    if not chosen:
        sys.exit("Error: no genomes selected (%d requested names without a match)" % len(unmatched))
    if len(chosen) > a.max_genomes:
        sys.exit("Error: %d genomes selected (> --max-genomes %d); narrow the selection" % (len(chosen), a.max_genomes))

    os.makedirs(a.outdir, exist_ok=True)
    with open(os.path.join(a.outdir, "genome_taxonomy.tsv"), "w") as fo:
        fo.write("genome\tclassification\n")
        fo.writelines("%s\t%s\n" % (g, m["lineage"]) for g, m in sorted(chosen.items()))
    with open(os.path.join(a.outdir, "genome_quality.tsv"), "w") as fo:
        fo.write("Name\tCompleteness\tContamination\n")
        fo.writelines("%s\t%s\t%s\n" % (g, m["completeness"], m["contamination"]) for g, m in sorted(chosen.items())
                      if m["completeness"] and m["contamination"])
    with open(os.path.join(a.outdir, "selection.tsv"), "w") as fo:
        fo.write("genome\tspecies_rep\tlineage\treason\n")
        fo.writelines("%s\t%s\t%s\t%s\n" % (g, m["rep"], m["lineage"], "; ".join(why[g])) for g, m in sorted(chosen.items()))
    with open(os.path.join(a.outdir, "unmatched_species.txt"), "w") as fo:
        fo.writelines(x + "\n" for x in sorted(set(unmatched)))
    print("[mgnify] %d genomes selected from %s; %d requested names without a match (unmatched_species.txt)" % (
        len(chosen), src, len(set(unmatched))), file=sys.stderr)
    if a.dry_run:
        return
    ok = 0
    for g, m in sorted(chosen.items()):
        url = m["url"]
        if not url:
            sys.exit("Error: the metadata table has no download URL column (FTP_download); download the GFFs by hand")
        name = os.path.basename(url.split("?")[0]) or g + ".gff"
        path = os.path.join(a.outdir, name if name.startswith(g) else g + (".gff.gz" if name.endswith(".gz") else ".gff"))
        if os.path.exists(path) and os.path.getsize(path) > 0:
            ok += 1
            continue
        ok += fetch(url.replace("ftp://", "https://", 1) if url.startswith("ftp://ftp.ebi.ac.uk") else url, path)
    print("[mgnify] %d of %d genomes in %s" % (ok, len(chosen), a.outdir), file=sys.stderr)
    if ok < len(chosen):
        sys.exit(1)


if __name__ == "__main__":
    main()
