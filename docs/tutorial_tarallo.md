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
| 4 | assemble, annotate and bin each metagenome | IMP3 v3.0.1 (Apptainer image), metagenome-only |
| 5 | sample sheet and configuration | `srna_make_samples.sh`, `srna_pipeline.config.sh` |
| 6 | run all libraries, then comparison, reports and cohort | `submit_srna_pipeline.sh` |
| 7 | follow, resubmit, collect the results | `srna_status.sh`, `srna_pack_summary.sh` |

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

`--imp-root` is where the IMP3 outputs will be: `imp/<ID>` inside the IMP3 project of step 4.

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

The paper used IMP3 v3.0.1 (the v3.0.0 workflow with launcher fixes) and the v3.0.0 Apptainer image from Zenodo
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

Assembly and read mapping are reproducible: a rerun of the same sample gives the same contigs (± 2 of ~50,000) and
the same reads on contigs. Binning is not fully deterministic; a rerun can move single bins between quality classes.

Without IMP3, assemble, annotate, bin and map with other tools and convert the outputs with
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
