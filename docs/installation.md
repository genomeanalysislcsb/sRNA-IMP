# Installation

## Recommended Environment

The project ships with a micromamba/conda environment definition at `envs/srna-imp.yaml`.

```bash
micromamba create -p /scratch/users/pmay/sRNA/sRNA-IMP/.micromamba/envs/srna-imp -f envs/srna-imp.yaml
micromamba activate /scratch/users/pmay/sRNA/sRNA-IMP/.micromamba/envs/srna-imp
```

## Core Runtime Tools

Important tools include:
- `fastp`
- `fastqc`
- `bowtie` and `bowtie-build`
- `samtools`
- `bedtools`
- `bwa`
- `srnaMapper`
- `kraken2`
- `ribodetector_cpu`
- `cmscan`
- `RNAfold`
- `multiqc`

## Source Code

The main entrypoint is `scripts/srna_imp.py`.
