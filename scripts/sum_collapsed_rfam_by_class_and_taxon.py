#!/usr/bin/env python3

import argparse
from collections import defaultdict


def count_from_name(name):
    """
    Count is expected as the suffix after the last '-'
    Example: seq123-17 -> 17
    """
    try:
        return int(name.rsplit("-", 1)[1])
    except Exception:
        return 1

def read_rfamseq_reduced(path):
    acc2desc = {}
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")

            if len(cols) == 2:
                acc, desc = cols
            elif len(cols) >= 3:
                acc, desc = cols[0], cols[2]  # <-- use NAME
            else:
                continue

            acc2desc[acc] = desc

    return acc2desc

def read_rfamseq(path):
    """
    Build:
      accession.version -> description/source text

    Expected rfamseq.txt columns include:
      col1 = accession.version
      col7 = description/source text
    """
    acc2desc = {}
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 7:
                continue
            acc = cols[0].strip()
            desc = cols[6].strip() if cols[6].strip() else "unknown"
            acc2desc[acc] = desc
    return acc2desc


def read_family_txt(path):
    """Rfam family.txt(.gz): col1 rfam_acc, col2 rfam_id, col4 description, col19 type
    (e.g. 'Gene; snRNA; snoRNA; CD-box;'). Returns acc -> (rfam_id, description, type)."""
    import gzip
    fam = {}
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 19 or not c[0].startswith("RF"):
                continue
            fam[c[0]] = (c[1], c[3], clean_type(c[18]))
    return fam


def clean_type(t):
    """'Gene; snRNA; snoRNA; CD-box;' -> 'snRNA;snoRNA;CD-box', 'Cis-reg; riboswitch;' ->
    'Cis-reg;riboswitch', 'Gene;' -> 'Gene'"""
    parts = [x.strip() for x in t.split(";") if x.strip()]
    if len(parts) > 1 and parts[0] == "Gene":
        parts = parts[1:]
    return ";".join(parts) or "unknown"


def parse_ref(ref):
    """
    Expected relabeled reference format:
      class|RFxxxxx|REGION_ID

    Example:
      antisense|RF00059|AMSQ01000006.1/153499-153391
    """
    parts = ref.split("|", 2)

    rfam_class = parts[0] if len(parts) > 0 else "unknown"
    rfam_family = parts[1] if len(parts) > 1 else "unknown"
    region_id = parts[2] if len(parts) > 2 else "unknown"

    seq_acc = region_id.split("/", 1)[0] if "/" in region_id else region_id
    return rfam_class, rfam_family, region_id, seq_acc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True, help="Bowtie hits text")
    ap.add_argument("-s", "--rfamseq", required=True, help="rfamseq.txt")
    ap.add_argument("--class-out", required=True, help="TSV: counts by Rfam class")
    ap.add_argument("--taxon-out", required=True, help="TSV: counts by organism/source text")
    ap.add_argument("--class-taxon-out", required=True, help="TSV: counts by class and organism/source text")
    ap.add_argument("--family-out", default=None, help="Optional TSV: counts by RF family")
    ap.add_argument("--accession-out", default=None, help="Optional TSV: counts by sequence accession")
    ap.add_argument("--family-txt", default=None, help="Optional Rfam family.txt(.gz): adds rfam_id, "
                    "description and type to the family table and enables --type-out")
    ap.add_argument("--type-out", default=None, help="Optional TSV: counts by Rfam type (needs --family-txt)")
    ap.add_argument("--type-taxon-out", default=None, help="Optional TSV: counts by Rfam type and taxon")
    args = ap.parse_args()

    acc2desc = read_rfamseq_reduced(args.rfamseq)
    fam = read_family_txt(args.family_txt) if args.family_txt else {}
    type_counts = defaultdict(int)
    type_taxon_counts = defaultdict(int)

    class_counts = defaultdict(int)
    taxon_counts = defaultdict(int)
    class_taxon_counts = defaultdict(int)
    family_counts = defaultdict(int)
    accession_counts = defaultdict(int)

    with open(args.input) as f:
        for line in f:
            if not line.strip():
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue

            qname = cols[0]
            ref = cols[2]
            c = count_from_name(qname)

            rfam_class, rfam_family, region_id, seq_acc = parse_ref(ref)
            desc = acc2desc.get(seq_acc, "unknown")

            class_counts[rfam_class] += c
            taxon_counts[desc] += c
            class_taxon_counts[(rfam_class, desc)] += c
            family_counts[rfam_family] += c
            accession_counts[seq_acc] += c
            if fam:
                rt = fam.get(rfam_family, ("", "", "unknown"))[2]
                type_counts[rt] += c
                type_taxon_counts[(rt, desc)] += c

    with open(args.class_out, "w") as out:
        out.write("rfam_class\tcount\n")
        for k in sorted(class_counts):
            out.write(f"{k}\t{class_counts[k]}\n")

    with open(args.taxon_out, "w") as out:
        out.write("taxon\tcount\n")
        for k in sorted(taxon_counts):
            out.write(f"{k}\t{taxon_counts[k]}\n")

    with open(args.class_taxon_out, "w") as out:
        out.write("rfam_class\ttaxon\tcount\n")
        for cls, taxon in sorted(class_taxon_counts):
            out.write(f"{cls}\t{taxon}\t{class_taxon_counts[(cls, taxon)]}\n")

    if args.family_out:
        with open(args.family_out, "w") as out:
            if fam:
                out.write("rfam_family\trfam_id\tdescription\ttype\tcount\n")
                for k in sorted(family_counts):
                    fid, fdesc, ft = fam.get(k, ("", "", "unknown"))
                    out.write(f"{k}\t{fid}\t{fdesc}\t{ft}\t{family_counts[k]}\n")
            else:
                out.write("rfam_family\tcount\n")
                for k in sorted(family_counts):
                    out.write(f"{k}\t{family_counts[k]}\n")

    if fam and args.type_out:
        with open(args.type_out, "w") as out:
            out.write("rfam_type\tcount\n")
            for k in sorted(type_counts, key=lambda x: -type_counts[x]):
                out.write(f"{k}\t{type_counts[k]}\n")
    if fam and args.type_taxon_out:
        with open(args.type_taxon_out, "w") as out:
            out.write("rfam_type\ttaxon\tcount\n")
            for t, taxon in sorted(type_taxon_counts):
                out.write(f"{t}\t{taxon}\t{type_taxon_counts[(t, taxon)]}\n")

    if args.accession_out:
        with open(args.accession_out, "w") as out:
            out.write("accession\tcount\n")
            for k in sorted(accession_counts):
                out.write(f"{k}\t{accession_counts[k]}\n")


if __name__ == "__main__":
    main()