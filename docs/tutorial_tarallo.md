# Tutorial: reanalysing Tarallo et al. 2019

This tutorial runs the whole way from public data to a cohort report, on the stool small RNA and shotgun
metagenome data of [Tarallo et al. 2019](https://doi.org/10.1128/mSystems.00289-19) (mSystems 4:e00289-19):
80 participants (24 healthy, 27 adenoma, 29 colorectal cancer) with one small RNA library (GEO
[GSE132236](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132236)) and one paired-end metagenome
(SRA [SRP136711](https://www.ncbi.nlm.nih.gov/sra/?term=SRP136711)) each. It is the dataset of the
validation in the sRNA-IMP paper.

| step | what | tool |
| --- | --- | --- |
| 1 | find the runs and pair each small RNA library with its metagenome | `scripts/public_runtable.py` |
| 2 | sample IDs, group metadata, download script, IMP3 input list, sample template | `scripts/public_prepare.py` |
| 3 | download the reads (md5-checked) | `download.sh` |
| 4 | assemble, annotate and bin each metagenome | IMP3, metagenome-only |
| 5 | sample sheet and configuration | `srna_make_samples.sh`, `srna_pipeline.config.sh` |
| 6 | run all libraries, then comparison, reports and cohort | `submit_srna_pipeline.sh` |
| 7 | follow, resubmit, collect the results | `srna_status.sh`, `srna_pack_summary.sh` |

You need sRNA-IMP with its environments and databases ([Installation](installation.md)), IMP3 with its databases
([IMP3 installation](https://genomeanalysislcsb.github.io/IMP3/)), a SLURM cluster, and internet access on the node
where you run steps 1–3 (on most clusters the access node).

## 0. Folders

```bash
P=/path/to/sRNA-IMP                  # sRNA-IMP checkout
I=/path/to/IMP3                      # IMP3 checkout
B=/path/to/tarallo2019               # data and results (fastq/, IMP/, sRNA-IMP output); a few hundred GB
mkdir -p $B && cd $B
```

## 1. Runs and pairing

```bash
python3 $P/scripts/public_runtable.py --srna GSE132236 --metag SRP136711 --out tarallo2019
```

The script resolves the GEO series to its SRA study, fetches both ENA run tables (`tarallo2019.srna.tsv`,
`tarallo2019.metag.tsv`) and proposes a pairing by sample (`tarallo2019.pairs.tsv`). **Check the pairing**: one
line per participant, one small RNA run and one metagenome run each.

## 2. IDs, metadata, download script, templates

```bash
python3 $P/scripts/public_prepare.py --prefix tarallo2019 --outdir . --base $B --type sRNA
```

It writes, in the current folder:

| file | content |
| --- | --- |
| `ids.txt` | one sample ID per line (the subject code, e.g. `VF10`) |
| `metadata.tsv` | `sample_id`, `group` (HC / Adenoma / CRC from the GEO titles), subject, runs, read counts: the cohort metadata |
| `download.sh` | ENA downloads named by sample ID: `fastq/srna/<ID>.fastq.gz`, `fastq/metag/<ID>_R1/_R2.fastq.gz` |
| `imp3_samples.tsv` | sample ID, metagenome R1 and R2, IMP3 output folder |
| `sample.template.tsv` | the sRNA-IMP sample-sheet template, with the IMP3 output paths |

It reports the number of samples per group (24 HC, 27 Adenoma, 29 CRC) and warns about runs that are not paired-end
or samples with several metagenome runs.

## 3. Download

```bash
bash download.sh            # resumable (wget -c); every file is md5-checked against ENA
```

Check the adapter of the small RNA libraries before running anything:

```bash
zcat fastq/srna/$(head -1 ids.txt).fastq.gz | head -400000 | awk 'NR % 4 == 2' | grep -c AGATCGGAAGAGCACACGTCT
```

About 97 % of the 100,000 reads contain the TruSeq read-1 adapter, which is the default `ADAPTER_3P`.

## 4. IMP3, metagenome only

Run IMP3 once per sample with an empty metatranscriptome. A config template (keys as in the
[IMP3 configuration](https://genomeanalysislcsb.github.io/IMP3/); check them against your IMP3 version):

```yaml
# imp3.template.yaml
sample: SAMPLE
outputdir: "BASE/IMP/SAMPLE"
db_path: "/path/to/IMP3_databases"
tmp_dir: "BASE/tmp"
steps: "preprocessing assembly analysis binning taxonomy summary"
raws:
  Metagenomics: "BASE/fastq/metag/SAMPLE_R1.fastq.gz BASE/fastq/metag/SAMPLE_R2.fastq.gz"
  Metatranscriptomics: ""
assembly:
  assembler: megahit
```

One config per sample and one SLURM submission each, with the IMP3 cluster profile:

```bash
mkdir -p imp3_configs
while read id; do
  sed -e "s|SAMPLE|$id|g" -e "s|BASE|$B|g" imp3.template.yaml > imp3_configs/$id.yaml
  snakemake -s $I/Snakefile --profile $I/cluster_profile/ --configfile imp3_configs/$id.yaml
done < ids.txt
```

(or IMP3's `runscripts/runIMP3_submit.sh` for many samples). When the first sample has finished, check that the
paths in `sample.template.tsv` exist in its output folder; a metagenome-only IMP3 run names them

| sample-sheet column | IMP3 output |
| --- | --- |
| `assembly` | `Assembly/mg.assembly.merged.fa.gz` |
| `gff` | `Analysis/annotation/annotation_CDS_RNA_hmms.gff.gz` |
| `contig2bin` | `Binning/MetaWrap/bins.contigs.gz` |
| `gtdbtk` | `Stats/all_bin_stats.tsv.gz` |
| `metaG_features` | `Analysis/annotation` (the `mg.ID.*_counts.tsv` tables) |
| `metaG_contigs` | `Stats/mg/mg.assembly.contig_read_count.txt.gz` |

If a name differs, rerun step 2 with the matching `--imp-assembly`, `--imp-gff`, `--imp-bins`, `--imp-gtdb` or
`--imp-mg-contigs`. Without IMP3, assemble, annotate, bin and map with other tools and convert the outputs with
`scripts/srna_prepare_inputs.py` ([inputs](inputs.md#where-the-assembly-inputs-come-from)).

## 5. Sample sheet and configuration

```bash
$P/srna_make_samples.sh -c -o samples.tsv sample.template.tsv ids.txt    # -c: reports missing files
cp $P/srna_pipeline.config.sh tarallo2019.config.sh
```

Set the database and environment paths as for any run ([Configuration](configuration.md)), and for this cohort:

```bash
SCRIPTS_DIR=$P/scripts
SAMPLE_TYPE=stool
COHORT_METADATA=$B/metadata.tsv
COHORT_GROUP_COL=group
COHORT_GROUP_ORDER=HC,Adenoma,CRC
COHORT_REFERENCE=HC
COHORT_BLOCK=""                  # one site, no blocks
COHORT_MODEL=""                  # no age or sex in the public metadata
```

All other settings stay at their defaults (exact host matches, 20 nt host minimum, assembly rescue, the tRNA CCA
and miRNA host steps). Check environments, databases and inputs:

```bash
$P/srna_pipeline.sh -s samples.tsv -o $B/sRNA-IMP -c tarallo2019.config.sh --check
```

## 6. Run

From the access node (the script submits the jobs itself):

```bash
$P/submit_srna_pipeline.sh -s samples.tsv -o $B/sRNA-IMP -c tarallo2019.config.sh \
    -n 16 -T 12:00:00 --mem 64G -J tarallo
```

This submits one array task per library (all modules from `qc` to `crispr`) and one final job that starts when the
array is done and runs `--compare --report --cohort`. There is no metatranscriptome, so the comparison is small RNA
vs metagenome only, and the reports and Methods say so.

## 7. Follow, resubmit, collect

```bash
$P/srna_status.sh -s samples.tsv -o $B/sRNA-IMP -J tarallo              # done / running / failed per library
$P/srna_status.sh -s samples.tsv -o $B/sRNA-IMP -J tarallo --resubmit   # failed libraries only
$P/srna_pack_summary.sh -o $B/sRNA-IMP                                  # small result files in one tarball
```

The results:

| where | what |
| --- | --- |
| `sRNA-IMP/<ID>/Report/<ID>.report.html` | per-sample report with QC verdict and generated Methods |
| `sRNA-IMP/Cohort/cohort.report.html` | cohort report: QC, host load, RNA classes, taxa, miRNAs, functions, group tests and models |
| `sRNA-IMP/Cohort/*.tsv` | the tables behind every cohort figure |

With the default settings the cohort report should show what the paper reports: host reads a median 0.75 % of
trimmed reads in healthy participants and 3.4 % in cancer; *Escherichia coli* and *Klebsiella pneumoniae* up in
cancer (log2 about +2 among microbial tRNA reads) and *Agathobacter rectalis* down; miR-200b-3p and miR-30c-5p up;
miR-21-5p, miR-378a-3p and miR-215-5p unchanged.
