# Tutorial 1: reanalysing Tarallo et al. 2019

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
| 4 | assemble, annotate and bin each metagenome | IMP3 v3.0.1 (Apptainer image), metagenome-only |
| 5 | sample sheet and configuration | `srna_make_samples.sh`, `srna_pipeline.config.sh` |
| 6 | run all libraries, then comparison, reports and cohort | `submit_srna_pipeline.sh` |
| 7 | follow, resubmit, collect the results | `srna_status.sh`, `srna_pack_summary.sh` |
| 8 | Fig. 3 of the paper | `paper/fig3_tarallo.py` |

You need sRNA-IMP with its environments and databases ([Installation](installation.md)), the IMP3 databases
([IMP3 databases](https://genomeanalysislcsb.github.io/IMP3/databases/)), Apptainer on the compute nodes, a SLURM cluster, and internet access on the node
where you run steps 1–3 (on most clusters the access node).

## 0. Folders

```bash
P=/path/to/sRNA-IMP                  # sRNA-IMP checkout
B=/path/to/tarallo2019               # data and results (fastq/, imp3_v3.0.0/, sRNA-IMP/); a few hundred GB
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
python3 $P/scripts/public_prepare.py --prefix tarallo2019 --outdir . --base $B --type sRNA \
    --imp-root $B/imp3_v3.0.0/imp
```

`--imp-root` is where the IMP3 outputs will be: `imp/<ID>` inside the IMP3 project of step 4. Reads already downloaded elsewhere: `--fastq-dir <dir>` (with `srna/` and `metag/` below it).

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

The paper used IMP3 v3.0.1 (the v3.0.0 workflow with launcher fixes and Mantis working outside `/dev/shm`) and the v3.0.0 Apptainer image from Zenodo
([doi:10.5281/zenodo.24003091](https://doi.org/10.5281/zenodo.24003091)). Snakemake runs on the host, every rule in
the container; the databases stay outside ([IMP3 databases](https://genomeanalysislcsb.github.io/IMP3/databases/)).

```bash
mkdir -p $B/imp3 && cd $B/imp3
git clone https://github.com/genomeanalysislcsb/IMP3 && git -C IMP3 checkout v3.0.1
R=https://zenodo.org/records/24003091/files
wget $R/imp3_3.0.0.sif $R/imp3_submodules_3.0.0.tar.gz $R/SHA256SUMS
sha256sum -c SHA256SUMS                                    # both files: OK
tar -xzf imp3_submodules_3.0.0.tar.gz -C IMP3              # submodules/ next to the Snakefile
I=$B/imp3/IMP3
SIF=$B/imp3/imp3_3.0.0.sif
DB=/path/to/imp3_databases
```

One IMP3 project for all samples, with the metagenomes from step 2 as its sample sheet:

```bash
$I/imp3 init $B/imp3_v3.0.0 --executor slurm --site ulhpc --container $SIF --db $DB    # other clusters: --site generic
awk -F'\t' 'BEGIN{OFS="\t"; print "sample","mg_r1","mg_r2","mt_r1","mt_r2"} NR>1{print $1,$2,$3,"",""}' \
    $B/imp3_samples.tsv > $B/imp3_v3.0.0/samples.tsv
cd $B/imp3_v3.0.0
sed -i "s|^sample_sheet:.*|sample_sheet: $B/imp3_v3.0.0/samples.tsv|" imp3.yaml
grep -n "^imp3_dir\|^container\|^db_path\|^sample_sheet\|^tmp_root\|^pre_commands\|^parallel_samples\|^submit_master" imp3.yaml
```

`imp3 init` takes the container mount points from the image, puts the temporary files on `$SCRATCH` (or into the
project), and on ULHPC loads the Apptainer module in the jobs (`pre_commands`). Check that `submit_master` uses a
QoS and time your account allows. Mantis (functional annotation) needs more memory than a standard node gives
for some of these samples; in `profile/cluster.yaml`, uncomment the `mantis_annot` block (ULHPC: iris `bigmem`, see
the [IMP3 FAQ](https://genomeanalysislcsb.github.io/IMP3/faq/)).

Run three samples first, one per group, and check their outputs:

```bash
$I/imp3 run --samples VF00TL,GF12,VF114 --dry-run
$I/imp3 run --samples VF00TL,GF12,VF114 --submit
$I/imp3 status
```

Then all samples. `--submit` sends one Snakemake master per sample as a Slurm job, at most `-j` at a time (the others
wait in chains), so the masters leave room in the account's job limit for the rule jobs they submit:

```bash
$I/imp3 run --submit -j 15
$I/imp3 status | tail -1                                   # e.g. done 80
```

A sample whose master stops (time limit, a failed rule) is continued from its finished steps by
`$I/imp3 run --samples <ID> --submit`. A metagenome-only run writes the files `sample.template.tsv` expects:

| sample-sheet column | IMP3 output (`imp3_v3.0.0/imp/<ID>/`) |
| --- | --- |
| `assembly` | `Assembly/mg.assembly.merged.fa.gz` |
| `gff` | `Analysis/annotation/annotation_CDS_RNA_hmms.gff.gz` |
| `contig2bin` | `Binning/MetaWrap/bins.contigs.gz` |
| `gtdbtk` | `Stats/all_bin_stats.tsv.gz` |
| `metaG_features` | `Analysis/annotation` (the `mg.ID.*_counts.tsv.gz` tables) |
| `metaG_contigs` | `Stats/mg/mg.assembly.contig_read_count.txt.gz` |

Assembly and read mapping are reproducible. In 30 samples assembled twice (an earlier IMP3 run and this one), the
assembly size was identical in all 30, N50 in 29, the contig number differed by at most 5 (of 6,000–186,000), and the
metagenome reads on contigs were identical in all 30. Binning is not fully deterministic: the number of bins was the
same in 17 of 30 samples (± 2 otherwise), single bins moved between quality classes, and the GTDB genera of the
high- and medium-quality bins overlapped with a median Jaccard index of 1.00 (0.87–1.00).

Without IMP3, assemble, annotate, bin and map with other tools and convert the outputs with
`scripts/srna_prepare_inputs.py` ([inputs](inputs.md#where-the-assembly-inputs-come-from)).

## 5. Sample sheet and configuration

```bash
cd $B
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

`--check` notes when the collapsed tRNA and Rfam references (`TRNA_MEMBERS`, `RFAM_MEMBERS`) are not set. The first
`nonhuman` job then builds them next to `IDX_TRNA_ALL` and `IDX_RFAM_CLASSES` (a few minutes; other jobs wait on its
lock) and every later run reuses them. To build them before submitting, or once for a shared database, see
[Installation](installation.md#collapsed-trna-and-rfam-references). Without them, a read's equally good hits on the
hundreds of identical tRNA copies make the tRNA screen of a deep library take hours.

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
| `sRNA-IMP/Cohort/tables/`, `Cohort/stats/` | count tables per library, QC verdicts (`cohort.qc_flags.tsv`), group tests and model coefficients |

Eleven libraries with fewer than 1 M raw reads get a FAIL verdict and are left out of the statistics (69 of 80:
21 healthy, 23 adenoma, 25 cancer; list with `awk -F'\t' '$5=="no"' $B/sRNA-IMP/Cohort/tables/cohort.qc_flags.tsv`).
With the default settings the cohort report should show what the paper reports: host reads a median 0.76 % of
trimmed reads in healthy participants and 3.5 % in cancer; *Escherichia coli* and *Klebsiella pneumoniae* up in
cancer (log2 about +4 among microbial tRNA reads) and *Agathobacter rectalis* down; miR-30c-5p up in adenoma and
cancer and miR-200b-3p in adenoma; miR-21-5p, miR-378a-3p and miR-215-5p unchanged.

To repeat the cohort statistics without some libraries (here the five with more than 20 % host reads), list them as
manual FAIL decisions and write the cohort to a separate folder:

```bash
printf '%s\t*\tFAIL\thost > 20 %%\n' VF107 VF166 VF208 VF65 VF97 > hihost_exclude.tsv
cp tarallo2019.config.sh hihost.config.sh
printf 'QC_EXCLUDE=%s\nCOHORT_DIR=%s\n' "$PWD/hihost_exclude.tsv" "$B/sRNA-IMP/Cohort_nohighhuman" >> hihost.config.sh
$P/srna_pipeline.sh -s samples.tsv -o $B/sRNA-IMP -c hihost.config.sh --cohort --skip-env-check
```

## 8. Fig. 3 of the paper

```bash
python3 $P/paper/fig3_tarallo.py $B/sRNA-IMP $B/fig3_tarallo
```

The script reads the cohort tables (`Cohort/tables/cohort.sample_stats.tsv`, `cohort.sRNA.mirna.counts.tsv`,
`cohort.sRNA.genus.counts.tsv`) and the per-sample tRNA species counts
(`<ID>/Analysis/sRNA/nonhuman/trna/<ID>.trna_species.counts.tsv`), and writes `fig3_tarallo.png`, `.pdf` and
`.tsv` (every value and p-value in the figure). Libraries excluded by the cohort QC are left out, and libraries
without classified rRNA reads are left out of panel c. The tests are those of the cohort module (Mann-Whitney vs
HC, Spearman trend HC < adenoma < CRC, Kruskal-Wallis).

The original study reports no library exclusion. To show the same panels with all 80 libraries, write a second
cohort folder without the QC exclusion and point the script at it:

```bash
cp tarallo2019.config.sh all80.config.sh
printf 'COHORT_QC_DROP=none\nCOHORT_DIR=%s\n' "$B/sRNA-IMP/Cohort_all80" >> all80.config.sh
$P/srna_pipeline.sh -s samples.tsv -o $B/sRNA-IMP -c all80.config.sh --cohort --skip-env-check
python3 $P/paper/fig3_tarallo.py $B/sRNA-IMP $B/fig3_tarallo_all80 --cohort-dir $B/sRNA-IMP/Cohort_all80
```

Panel d shows human miRNA reads per library by group. To compare two runs of the same libraries instead (the paper
compares the earlier strict host settings with the defaults), give the second output folder:

```bash
python3 $P/paper/fig3_tarallo.py $B/sRNA-IMP $B/fig3_tarallo \
    --compare-outdir $B/sRNA-IMP_strict --compare-label "earlier strict" --label "final defaults"
```
