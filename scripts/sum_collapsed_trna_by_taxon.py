#!/usr/bin/env python3

import argparse
from collections import defaultdict


def count_from_name(name):
    """
    Collapsed read headers look like:
      seq123-17
    so count = suffix after last '-'
    """
    try:
        return int(name.rsplit("-", 1)[1])
    except Exception:
        return 1


def parse_trna_ref(ref):
    """
    Expected tRNA reference header format:
    SequenceID|GenomeID|Phylum/Class|Species|start|end|AA|Anticodon|...

    Returns a dict with the main fields.
    """
    parts = ref.split("|")

    data = {
        "sequence_id": parts[0] if len(parts) > 0 else "unknown",
        "genome_id": parts[1] if len(parts) > 1 else "unknown",
        "phylum_class": parts[2] if len(parts) > 2 else "unknown",
        "species": parts[3] if len(parts) > 3 else "unknown",
        "start": parts[4] if len(parts) > 4 else "unknown",
        "end": parts[5] if len(parts) > 5 else "unknown",
        "aa": parts[6] if len(parts) > 6 else "unknown",
        "anticodon": parts[7] if len(parts) > 7 else "unknown",
    }
    return data


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True, help="Bowtie hits text")
    ap.add_argument("--feature-out", required=True, help="Counts by full tRNA feature")
    ap.add_argument("--phylum-out", required=True, help="Counts by Phylum/Class")
    ap.add_argument("--species-out", required=True, help="Counts by Species")
    ap.add_argument("--aa-out", required=True, help="Counts by amino acid")
    ap.add_argument("--anticodon-out", required=True, help="Counts by anticodon")
    ap.add_argument("--phylum-species-out", required=True, help="Counts by Phylum/Class and Species")
    args = ap.parse_args()

    feature_counts = defaultdict(int)
    phylum_counts = defaultdict(int)
    species_counts = defaultdict(int)
    aa_counts = defaultdict(int)
    anticodon_counts = defaultdict(int)
    phylum_species_counts = defaultdict(int)

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

            meta = parse_trna_ref(ref)

            feature_counts[ref] += c
            phylum_counts[meta["phylum_class"]] += c
            species_counts[meta["species"]] += c
            aa_counts[meta["aa"]] += c
            anticodon_counts[meta["anticodon"]] += c
            phylum_species_counts[(meta["phylum_class"], meta["species"])] += c

    with open(args.feature_out, "w") as out:
        out.write("feature\tcount\n")
        for k in sorted(feature_counts):
            out.write(f"{k}\t{feature_counts[k]}\n")

    with open(args.phylum_out, "w") as out:
        out.write("phylum_class\tcount\n")
        for k in sorted(phylum_counts):
            out.write(f"{k}\t{phylum_counts[k]}\n")

    with open(args.species_out, "w") as out:
        out.write("species\tcount\n")
        for k in sorted(species_counts):
            out.write(f"{k}\t{species_counts[k]}\n")

    with open(args.aa_out, "w") as out:
        out.write("aa\tcount\n")
        for k in sorted(aa_counts):
            out.write(f"{k}\t{aa_counts[k]}\n")

    with open(args.anticodon_out, "w") as out:
        out.write("anticodon\tcount\n")
        for k in sorted(anticodon_counts):
            out.write(f"{k}\t{anticodon_counts[k]}\n")

    with open(args.phylum_species_out, "w") as out:
        out.write("phylum_class\tspecies\tcount\n")
        for (phylum, species) in sorted(phylum_species_counts):
            out.write(f"{phylum}\t{species}\t{phylum_species_counts[(phylum, species)]}\n")


if __name__ == "__main__":
    main()