# Installation & databases

The pipeline is a set of Bash scripts plus Python helpers; there is nothing to compile. It needs three
conda/micromamba environments (main, ncRNA discovery, RiboDetector), a set of reference databases, and the helper scripts in `SCRIPTS_DIR`.

## Environments

Every module runs in its own environment, which the pipeline activates automatically
(see [Automatic activation](#automatic-activation)).

| environment | config variable | modules | main tools |
| --- | --- | --- | --- |
| `gut_srna_pipeline` | `ENV_DEFAULT` (`ENV_QC`, `ENV_CONTAM`, `ENV_PREP`, `ENV_SPLIT`, `ENV_HUMAN`, `ENV_NONHUMAN`, `ENV_ASSEMBLY`, `ENV_COMPARE`, `ENV_REPORT`) | qc, contam, prep, split, human, nonhuman, assembly, compare, report | cutadapt ≥ 4.4, FastQC (+ Java), bowtie 1, samtools, seqkit, pigz, kraken2, bracken, BEDOPS (gff2bed), python3, matplotlib |
| `ribodetector` | `RIBO_DETECTOR_CPU` (full path, not activated) | nonhuman (rRNA filter) | RiboDetector (CPU, PyTorch) |
| `ncrna_discovery` | `ENV_NOVEL` (`ENV_FRAGMENTS`, `ENV_BINS`) | novel, fragments, bins | samtools, bedtools, seqkit, ARAGORN, Infernal (`cmscan`), ViennaRNA (`RNAfold`), optionally tRNAscan-SE |

RiboDetector is called by its full path (`RIBO_DETECTOR_CPU`) and may live in any environment.

### Creating the environments

Environment files are in `envs/`: only the tools the pipeline calls, pinned to the versions of the production
environments (2026-10), so a rebuild reproduces the results:

```bash
micromamba create -p /path/to/envs/gut_srna_pipeline -f envs/gut_srna_pipeline.yml
micromamba create -p /path/to/envs/ncrna_discovery   -f envs/ncrna_discovery.yml
micromamba create -p /path/to/envs/ribodetector      -f envs/ribodetector.yml
micromamba create -p /path/to/envs/multiqc           -f envs/multiqc.yml        # optional, --cohort
# conda: conda env create -p <path> -f <file>
```

Then point `ENV_DEFAULT` and `ENV_NOVEL` in the config to these paths (or use `-n <name>` and the
names), and set `RIBO_DETECTOR_CPU=/path/to/envs/ribodetector/bin/ribodetector_cpu`. RiboDetector has its own
environment because of PyTorch (large, with strict Python/numpy constraints); the pipeline only calls it by path.

| file | tools | used by |
| --- | --- | --- |
| `gut_srna_pipeline.yml` | python 3.11, cutadapt ≥ 4.4, FastQC, bowtie 1.3, samtools, seqkit, kraken2, bracken, matplotlib, numpy, openpyxl, pigz, BEDOPS, gawk, pandoc (optional, `.docx` Methods) | qc, contam, prep, split, human, nonhuman, assembly, compare, report |
| `ncrna_discovery.yml` | python 3.11, samtools ≥ 1.18, bedtools ≥ 2.31, seqkit, ARAGORN, Infernal, ViennaRNA, tRNAscan-SE (optional), gawk, coreutils | novel, fragments, bins |
| `ribodetector.yml` | RiboDetector (CPU) | nonhuman (rRNA filter) |
| `multiqc.yml` | MultiQC (optional) | --cohort: FastQC / cutadapt overview of all samples |

### Automatic activation

For each module the environment `ENV_<MODULE>` (default `ENV_DEFAULT`) is activated, in this order of
attempts:

1. already active → nothing to do
2. `ENV_MODULES` set → `module load $ENV_MODULES` (Lmod / environment modules)
3. micromamba: `MAMBA_BIN`, `PATH`, `$MAMBA_EXE`, `~/.local/bin`, `~/bin`, `~/micromamba/bin`
4. conda: `CONDA_SH`, `$CONDA_EXE`, `~/miniforge3`, `~/mambaforge`, `~/miniconda3`, `~/anaconda3`
5. last resort: an environment given as a path is used by putting its `bin/` first in `PATH` (with a warning)

If activation fails the run stops with a message naming `MAMBA_BIN`, `CONDA_SH` and `ENV_MODULES`.
`-e/--env` overrides `ENV_DEFAULT` for a single run.

## Databases and references

| purpose | variable | content |
| --- | --- | --- |
| human genome | `HOST_GENOME_FA` | T2T-CHM13v2 genome FASTA (incl. chrM) |
| human transcriptome | `HOST_TX_FA` | T2T-CHM13v2.0 transcripts; `''` skips the transcriptome step |
| human ncRNA classes | `IDX_HSA_MIRNA`, `IDX_HSA_TRNA`, `IDX_HSA_RRNA`, `IDX_HSA_NCRNA_OTHER` | bowtie indices: mature miRNA, nuclear + mt tRNA, rRNA + mt rRNA, sn/sno/Y/vault/… RNA |
| non-human tRNA | `IDX_TRNA_ALL` | bowtie index of microbial tRNAs (with taxonomy in the headers) |
| non-human miRNA-like | `IDX_MIRNA_EUK` | bowtie index of eukaryotic mature miRNAs |
| Rfam classes | `IDX_RFAM_CLASSES`, `RFAMSEQ_TXT`, `RFAM_FAMILY_TXT` | bowtie index of Rfam sequences (headers `description\|RFxxxxx\|accession/region`) + accession → taxon table + Rfam `family.txt` for family types |
| rRNA taxonomy | `KRAKEN_DB`, `BRACKEN_DB` | Kraken2 database (GTDB genome representatives) with `database50mers.kmer_distrib` |
| Rfam models | `RFAM_CM` | `Rfam.cm`, pressed with `cmpress` (`Rfam.cm.i1m` etc.) |
| contaminants | `CONTAM_FA` | optional, space-separated FASTAs (UniVec, PhiX, kit contaminants) |

### Bowtie indices built on demand

Indices for `HOST_GENOME_FA`, `HOST_TX_FA` and `CONTAM_FA` are searched at `HOST_*_IDX` or next to the FASTA
(`<name>.bowtie.*`). If missing they are built **once**: next to the FASTA when that directory passes a
real write test, otherwise in `HOST_INDEX_TMPDIR` (default `<outdir>/tmp/host_index`). Concurrent jobs wait
on a lock directory instead of building twice. Building the human genome index takes a few hours.

The bowtie index of each sample assembly (`Assembly/<sample>.mgmt.assembly.merged.fa.*.ebwt`) is built once
per sample and shared by the `split` and `novel` modules.

## Checking an installation

```bash
./srna_pipeline.sh -s samples.tsv -o <outdir> --check-env
```

lists every tool, helper script and database of the selected modules as `OK` or `MISSING`, each checked
inside the module's own environment. Every normal run performs the same check first and stops if anything
is missing (`--skip-env-check` bypasses it).
