# Databases

`sRNA-IMP` expects the database layer to be prepared ahead of pipeline execution. The pipeline config should point to the prepared reference files or indexes, but the documentation here explains the structure and build logic without relying on internal absolute paths.

## Overview

The database setup used for this project is organized into the following groups:
- `human/`
- `nonhuman/miRNAs/`
- `nonhuman/tRNA/`
- `rfam/`
- `silva/`
- `mirbase/`

## Human References

The human setup combines miRBase-derived miRNA references with Ensembl ncRNA references.

### miRNA references

The `mirbase/` setup downloads:
- `hairpin.fa`
- `mature.fa`

Human entries are extracted into:
- `hsa_hairpin_miRNA.fa`
- `hsa_mature_miRNA.fa`

and Bowtie indexes are built for both.

### Ensembl ncRNA references

The human README shows the use of the Ensembl human ncRNA FASTA to split records by `gene_biotype`, including:
- `lncRNA`
- `miRNA`
- `misc_RNA`
- `Mt_rRNA`
- `Mt_tRNA`
- `ribozyme`
- `rRNA`
- `scaRNA`
- `snoRNA`
- `snRNA`
- `sRNA`
- `vault_RNA`

From these subtype FASTAs, the setup builds consolidated references used by the pipeline:
- `hsa_tRNA_all`
- `hsa_rRNA_all`
- `hsa_otherRNA`
- `hsa_mature_miRNA`

These are the main human indexes expected by the host-specific workflow.

## Non-human References

### Non-human miRNA-like references

The `nonhuman/miRNAs/` setup contains:
- `nonhuman_miRNA.fa`
- the Bowtie index `mirna_euk_all`

This is used in the non-host branch for miRNA-like sequence assignment.

### Non-human tRNA references

The `nonhuman/tRNA/` setup contains multiple source FASTA files from different origins, such as:
- archaeal
- bacterial
- fungal
- plant
- viral
- phage
- plasmid
- environmental
- SRA-derived
- chloroplast-related

These are merged into:
- `trna_all.fa`
- Bowtie index `trna_all`

This merged index is the main non-host tRNA reference used by the pipeline.

## Rfam Setup

The `rfam/` setup is more elaborate and supports both class-level lightweight classification and covariance-model-based ncRNA discovery support.

### Main downloaded resources

The README indicates use of:
- `Rfam.fa.gz`
- `Rfam.cm.gz`
- `family.txt.gz`
- `rfamseq.txt`
- `full_region.txt.gz`

### Class-level Rfam reference

The class-level Rfam setup normalizes family annotations into practical groups such as:
- `rRNA`
- `tRNA`
- `miRNA`
- `snRNA`
- `snoRNA`
- `antisense`
- `ribozyme`
- `leader`
- `CRISPR`
- `riboswitch`
- `other`

The documented build flow is:
1. derive family-to-class mappings from `family.txt.gz`
2. normalize those classes into a cleaned table
3. build a combined class FASTA with `build_rfam_class_fasta.py`
4. extract a reduced accession-to-taxid table for downstream summarization

Important derived files include:
- `rfam_classes.fa`
- `rfam_classes.accessions.txt`
- `rfamseq.reduced.tsv`
- `rfamseq.reduced.named.tsv`

The pipeline uses the Bowtie index built from `rfam_classes.fa` for fast non-host Rfam class assignment.

### Full covariance models

For novel ncRNA discovery, the setup also uses the complete `Rfam.cm` covariance model database. The README notes compression of the covariance model file with:
- `cmpress`

This full CM file is appropriate for `cmscan`-based candidate filtering.

## Rfam Taxonomy Repository

The database tree also includes the `rfam-taxonomy/` repository, which provides domain-specific Rfam family subsets and clan files.

Important files there include:
- `rfam-taxonomy.py`
- `scripts/rfam_db.py`
- `domains/all-domains.csv`
- `domains/bacteria.csv`
- `domains/eukaryota.csv`
- `domains/viruses.csv`
- matching `.clanin` files per domain

The bundled README explains that these outputs can be used to:
- create domain-specific CM subsets with `cmfetch`
- run `cmscan` with matching `--clanin` files
- support stricter domain-focused annotation strategies

This is especially useful when you want domain-specific Rfam covariance model subsets rather than the full `Rfam.cm` file.

## SILVA

The `silva/` setup includes download of SILVA SSU reference FASTA data. This supports broader rRNA-oriented reference work, even though the current native pipeline primarily uses `ribodetector_cpu` plus Kraken2 for the non-host rRNA branch.

## Helper Scripts Used During Database Preparation and Analysis

The setup around these databases uses several local Python helper scripts that are worth documenting explicitly:
- `build_rfam_class_fasta.py`
- `taxid_to_name.py`
- `collapse_fastq_to_counted_fasta.py`
- `sum_collapsed_bowtie_hits.py`
- `sum_collapsed_trna_by_taxon.py`
- `sum_collapsed_rfam_by_class_and_taxon.py`
- `sam_fractional_counts.py`
- `summarize_human_subtree.py`

These scripts are part of the broader project tooling and support both reference preparation and downstream summarization.

## Tools Mentioned in the Database Build Notes

The README files and companion repository indicate use of tools such as:
- `wget`
- `gunzip` / `zcat`
- `grep`, `cut`, `sort`, `uniq`, `awk`, `sed`
- `seqkit`
- `bowtie-build`
- `cmpress`
- `cmfetch`
- `cmscan`
- `curl`
- Python helper scripts

## Practical Recommendation for Public Documentation

For public docs, it is best to describe databases by logical role and filename pattern rather than by internal infrastructure paths. A good pattern is to document:
- the source database or release family
- the derived FASTA or CM artifact names
- the indexing tool used
- which workflow consumes the result

That keeps the setup reproducible without exposing site-specific storage layout.
