# Database download and setup guide for sRNA-IMP

This bundle documents and automates setup of the major databases used by sRNA-IMP.

## Covered resources

### Human arm
- Human miRNA reference
- Human tRNA reference
- Human rRNA reference
- Human other ncRNA reference

### Non-human arm
- RiboDetector model setup (tool installation, no large external DB required)
- Kraken2 standard taxonomy database
- Bracken k-mer distribution files
- Rfam covariance models
- Optional non-human miRNA references
- Optional combined tRNA references

### Novel ncRNA discovery
- Rfam covariance models
- ARAGORN / tRNAscan-SE installation via conda environment

## Recommended top-level directory layout

```bash
databases/
├── human/
│   ├── mirna/
│   ├── trna/
│   ├── rrna/
│   └── otherrna/
├── nonhuman/
│   ├── kraken2/
│   ├── bracken/
│   ├── mirna/
│   ├── trna/
│   └── rfam/
└── shared/
    └── rfam/
```

## Quick-start

1. Create the conda / micromamba environment that contains the required tools.
2. Edit `config_examples/database_paths.example.sh`
3. Run:
   ```bash
   bash scripts/download_and_setup_databases.sh
   ```

## Notes on data sources

### Rfam
The official Rfam genome-annotation guide recommends downloading `Rfam.cm.gz` and `Rfam.clanin`, then indexing `Rfam.cm` with `cmpress`. Use `cmscan --rfam --cut_ga` for annotation against these models.

### Kraken2
The official Kraken2 manual documents `kraken2-build --standard --db <DB>` for building the standard database.

### Bracken
Bracken requires a built Kraken database plus generation of the `database*kmer_distrib` files with `bracken-build`.

### tRNAscan-SE
tRNAscan-SE 2.0 is available from the UCSC Lowe Lab GitHub and website. It depends on Infernal.

### ARAGORN
ARAGORN is also available through Bioconda, which is the simplest reproducible installation route for pipeline use.

### RiboDetector
RiboDetector is distributed via Bioconda and GitHub. It is a tool rather than a large downloadable reference database.

### miRNA references
For miRNA references, use either:
- **miRBase** (`mature.fa`, `hairpin.fa`) when you want broad species coverage
- **MirGeneDB** when you want a more conservative, curated animal miRNA set

## Human arm reference recommendations

### Human miRNA
Use one of:
- `miRBase mature.fa` filtered to `hsa-*`
- MirGeneDB mature sequences for human

Build a Bowtie index:
```bash
bowtie-build hsa_mature.fa hsa_mature
```

### Human tRNA
Preferred sources:
- GtRNAdb exports for human nuclear tRNAs
- human mitochondrial tRNAs from a curated FASTA

Build:
```bash
cat human_tRNA.fa human_mt_tRNA.fa > hsa_tRNA_all.fa
bowtie-build hsa_tRNA_all.fa hsa_tRNA_all
```

### Human rRNA
Use curated human cytosolic + mitochondrial rRNA FASTA, then build:
```bash
bowtie-build hsa_rRNA.fa hsa_rRNA
```

### Human other ncRNA
A practical source is Ensembl ncRNA FASTA:
```bash
bowtie-build Homo_sapiens.GRCh38.ncrna.fa hsa_other_ncRNA
```

## Non-human arm reference recommendations

### Kraken2 + Bracken
- Build Kraken2 standard DB
- Build Bracken files for your read length (50 bp typical for your pipeline)

### Rfam
Use the full `Rfam.cm`, then optionally extract a smaller subset later if runtime becomes an issue.

### Optional non-human miRNA reference
Use `miRBase mature.fa`, exclude `hsa`, and optionally subset to:
- plants
- fungi
- protozoa
- viruses

### Optional tRNA reference
For broader non-human tRNA screening, create a combined FASTA containing:
- bacterial tRNAs
- fungal / eukaryotic tRNAs
- plant chloroplast tRNAs
- plant mitochondrial tRNAs

## Runtime / storage expectations

| Resource | Approximate size | Comment |
|---|---:|---|
| Kraken2 standard DB | ~100 GB during build | largest resource |
| Rfam.cm + indices | several GB | shared across workflows |
| Bracken files | depends on Kraken DB | build once per DB |
| Bowtie miRNA/tRNA/rRNA refs | MB to low GB | lightweight |

## Verification checklist

After setup, verify:
- `cmpress` generated `Rfam.cm.i1f`, `i1i`, `i1m`, `i1p`
- Kraken2 DB contains `hash.k2d`, `opts.k2d`, `taxo.k2d`
- Bracken DB contains `database50mers.kmer_distrib` (or matching read length)
- each Bowtie index has `.1.ebwt` ... `.4.ebwt` files

## Suggested provenance file

Create a plain text manifest such as:
```text
Rfam=15.1
Kraken2_standard_built=2026-03-27
Bracken_read_length=50
miRBase=22.1
MirGeneDB_downloaded=2026-03-27
```