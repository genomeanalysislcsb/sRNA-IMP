# Databases

Typical database inputs include:
- host Bowtie indexes or host reference FASTAs
- host small RNA reference indexes
- Kraken2 database
- Bracken-compatible Kraken database
- non-host tRNA reference index
- non-host miRNA-like reference index
- Rfam class index
- `rfamseq.reduced.tsv`
- optional `Rfam.cm` covariance model file

All database locations are configured in YAML, not hard-coded per run.
