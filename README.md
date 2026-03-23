# sRNA-IMP

`sRNA-IMP` is a native, modular small-RNA pipeline for microbiome-associated smallRNA-seq data. It integrates the optimized logic from the existing preprocessing, host, non-host, assembly-based, and novel ncRNA discovery workflows into one config-driven entrypoint.

## What It Does

For each sample, `sRNA-IMP` can:
- trim and QC raw single-end smallRNA FASTQ input
- split reads into host and non-host fractions using required host genome or transcriptome references
- classify host small RNAs
- classify non-host small RNAs
- quantify reads against an optional assembly and annotation
- predict novel ncRNA candidates from non-host reads and an optional assembly
- aggregate summaries and build a per-sample MultiQC report

## Input Model

Required per sample:
- one raw, non-filtered single-end smallRNA FASTQ file
- one or more host reference FASTA files, or prebuilt host Bowtie indexes, for host/non-host splitting

Optional per sample:
- assembly FASTA
- assembly GFF

The host does not have to be human. It can be mouse, yeast, tomato, or another relevant host genome or transcriptome.

## Main Features

- native modular workflows instead of wrapping the legacy bash drivers
- single-sample and sequential multi-sample execution
- sample names defined in config and used consistently in the output tree
- automatic preprocessing prerequisite for downstream workflows
- resume support through per-workflow completion markers
- shared host Bowtie indexes built next to the original reference FASTA files and reused across samples
- automatic adapter trimming through `fastp` by default
- per-sample MultiQC reporting

## Workflows

Available workflows:
- `preprocessing`
- `host_specific`
- `nonhost_specific`
- `assembly_based`
- `novel_ncrna`
- `summary`
- `all`

## Config Styles

### Single sample

Use [config.example.yaml](/scratch/users/pmay/sRNA/sRNA-IMP/config/config.example.yaml) as the template for one sample per config file.

### Multi-sample

Use [multisample.example.yaml](/scratch/users/pmay/sRNA/sRNA-IMP/config/multisample.example.yaml) for a shared project-level config with a `samples:` list. Shared workflow, tool, and preprocessing settings stay at the top level, while each sample provides its own `name`, `fastq`, and optional assembly inputs.

### TSV-Driven Multi-sample

If you prefer to keep sample-specific inputs in a table, use [multisample.tsv_config.example.yaml](/scratch/users/pmay/sRNA/sRNA-IMP/config/multisample.tsv_config.example.yaml) together with [samples.example.tsv](/scratch/users/pmay/sRNA/sRNA-IMP/config/samples.example.tsv). In that mode, the YAML file keeps workflow parameters, database paths, and tool settings, while the TSV supplies `name`, `type`, `fastq`, `assembly_fasta`, and `assembly_gff` per sample.

## Run

Single sample:

```bash
python scripts/srna_imp.py --config config/BHC01.yaml --workflow all
```

Multi-sample sequential run:

```bash
python scripts/srna_imp.py --config config/multisample.example.yaml --workflow all
```

TSV-driven multi-sample run:

```bash
python scripts/srna_imp.py --config config/multisample.tsv_config.example.yaml --workflow all
```

Run only selected samples from a multi-sample config:

```bash
python scripts/srna_imp.py --config config/multisample.example.yaml --workflow all --sample BHC01 --sample BHC02
```

Keep going when one sample fails:

```bash
python scripts/srna_imp.py --config config/multisample.example.yaml --workflow all --keep-going
```

## Output Layout

```text
results/<TYPE>/<sample_name>/
├── Analysis/
│   ├── annotation/
│   ├── taxonomy/kraken/
│   └── <TYPE>/
│       ├── human/
│       ├── nonhuman/
│       └── novel/
├── Assembly/
├── Preprocessing/<TYPE>/
├── Stats/<TYPE>/
├── log_cluster_jobs/<TYPE>/
└── logs_cluster_jobs/
```

## GitHub Readiness

The project now includes:
- a repo-oriented `.gitignore`
- a multi-sample example config
- a basic GitHub Actions syntax-check workflow at `.github/workflows/python-syntax.yml`

A license and remote GitHub repository are still a project decision, so they should be chosen explicitly before publishing.

## Notes

- downstream workflows automatically run preprocessing first when needed
- assembly-based and novel workflows require `input.assembly_fasta` and `input.assembly_gff`
- helper Python scripts from `/scratch/users/pmay/sRNAs/scripts` are still reused where that is efficient and stable
- `ribodetector_cpu` is expected from the Conda environment and is not configured per sample
