# sRNA-IMP

**sRNA-IMP** is a modular pipeline for small-RNA sequencing (sRNA-seq) of microbiome samples, integrated with the
metagenomic and metatranscriptomic analysis of the same samples by [IMP3](https://github.com/genomeanalysislcsb/IMP3).
It processes intracellular (sRNA-Int) and extracellular (sRNA-Ext) libraries, separates human from microbial reads,
classifies both against curated small-RNA databases, places the microbial small RNAs on the sample's own contigs,
genes, bins and functions, compares them with metaT and metaG, and summarises a whole cohort with group statistics.

**Documentation:** <https://genomeanalysislcsb.github.io/sRNA-IMP/>

## What it does

| level | steps |
|---|---|
| per library (one sample-sheet row) | trimming + FastQC → library QC → human / non-human split → human classes (miRNA, tRNA, rRNA, other) → microbial classes (RiboDetector rRNA + Kraken 2/Bracken taxonomy, tRNA, miRNA-like, Rfam) → mapping to the sample's co-assembly (sense / antisense per gene and function) → novel ncRNA loci (incl. rescue of unannotated rRNA genes) → tRNA/rRNA fragments and RNA integrity → bins → CRISPR arrays |
| per sample | sRNA-Ext vs sRNA-Int vs metaT vs metaG comparisons; PDF/HTML report with QC verdict, MultiQC tables and automated Methods (STAR Protocols or Nature style) |
| cohort | metadata, QC-based exclusion, group tests, Bray–Curtis / PERMANOVA / NMDS, miRNAs, RNA classes, KEGG/Pfam functions, antisense index, block designs (site, family), covariate-adjusted models — all in Python (no R) |

## Quick start

```bash
# environments (micromamba or conda)
micromamba create -p envs/gut_srna_pipeline -f envs/gut_srna_pipeline.yml
micromamba create -p envs/ncrna_discovery   -f envs/ncrna_discovery.yml
micromamba create -p envs/ribodetector      -f envs/ribodetector.yml
micromamba create -p envs/multiqc           -f envs/multiqc.yml      # optional

# set database / environment paths in srna_pipeline.config.sh, then
./srna_pipeline.sh -s samples.tsv -o out --check                 # environment, databases, inputs
./submit_srna_pipeline.sh -s samples.tsv -o out -n 16 --mem 64G  # SLURM: all rows, then --compare --report --cohort
./srna_status.sh -s samples.tsv -o out                           # done / running / failed per row
```

See [Installation](https://genomeanalysislcsb.github.io/sRNA-IMP/installation/),
[Running](https://genomeanalysislcsb.github.io/sRNA-IMP/running/) and the
[configuration reference](https://genomeanalysislcsb.github.io/sRNA-IMP/configuration/).

## Repository layout

| path | content |
|---|---|
| `srna_pipeline.sh` | the pipeline (modules per row, `--compare`, `--report`, `--cohort`) |
| `srna_pipeline.config.sh` | databases, environments and parameters (adapt the paths to your system) |
| `submit_srna_pipeline.sh`, `srna_status.sh` | SLURM submission and run status / resubmission |
| `srna_make_samples.sh`, `sample.template.tsv`, `samples.example.tsv` | sample sheet from a template and an ID list |
| `metadata.template.tsv`, `scripts/srna_prepare_metadata.py` | cohort metadata |
| `srna_pack_summary.sh` | pack the summary tables and reports of a run |
| `scripts/` | Python helpers (counting, fragments, integrity, CRISPR, comparisons, reports, cohort statistics, Methods) |
| `envs/` | conda / micromamba environments |
| `docs/`, `mkdocs.yml` | documentation (MkDocs Material, deployed to GitHub Pages) |

## Citation and license

LCSB BioCore, University of Luxembourg. Please cite the repository and the tools listed in the automated Methods
of your run. License: see [LICENSE](LICENSE).

## Related tools

- [IMP3](https://github.com/genomeanalysislcsb/IMP3) ([docs](https://genomeanalysislcsb.github.io/IMP3/)): the metagenomic and metatranscriptomic pipeline whose assembly, annotation and bins sRNA-IMP uses
- [IMP explorer](https://github.com/genomeanalysislcsb/imp-explorer) ([docs](https://genomeanalysislcsb.github.io/imp-explorer/)): interactive study-level exploration of IMP3 results
