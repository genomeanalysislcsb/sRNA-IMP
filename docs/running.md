# Running the pipeline

## Command line

### Building the samples table from a template

For a whole study, `srna_make_samples.sh` writes the samples table from a template (header, comments and one
row per sRNA type, with `<ID>` = ID as listed, `<id>` = lower-case ID) and an ID list (one per line):

```bash
./srna_make_samples.sh -c -o samples.tsv sample.template.tsv expobiome.ids
```

`-c` checks that every path exists (missing files are listed, exit status 2; the table is still written).


```text
srna_pipeline.sh -s <samples.tsv> [options]

  -s, --samples FILE   sample sheet (required)
  -o, --outdir DIR     base output directory            [$PWD/imp]
  -c, --config FILE    config file                      [srna_pipeline.config.sh next to the script]
  -e, --env ENV        environment (name or path)       [ENV_DEFAULT]
  -m, --modules LIST   comma-separated modules          [all]
  -t, --threads N      threads                          [$SLURM_CPUS_PER_TASK or 16]
  -r, --row N          only data row N (used by SLURM arrays)
  -S, --sample ID      only rows of this sample
  -T, --type TYPE      only rows of this type
  -f, --force          re-run finished modules
      --check          check environment, tools, databases and all inputs
      --check-env      check environment, tools and databases only
      --skip-env-check skip the environment check
      --compare        run the per-sample comparison
      --report         write the per-sample PDF report (after --compare when both are given)
      --count-rows     print the number of data rows
```

Always pass `-c` explicitly on the HPC: without it the config **next to the script** is used, which may
be a different copy. The second line of every log (`# samples=… config=…`) shows the config used.

## Interactive

```bash
./srna_pipeline.sh -s samples.tsv -o /mnt/aiongpfs/projects/expobiome/sRNAs/imp \
  -c srna_pipeline.config.sh --check
./srna_pipeline.sh -s samples.tsv -o /mnt/aiongpfs/projects/expobiome/sRNAs/imp \
  -c srna_pipeline.config.sh -r 1 -t 16
```

## SLURM

```bash
./submit_srna_pipeline.sh -s samples.tsv -o <outdir> -c srna_pipeline.config.sh \
  -n 16 -T 2-00:00:00 --mem 64G -p batch -q normal
```

* one array task per sample-sheet row (`--rows 1-10,15`, `--max-parallel 20`),
* a comparison job after the array (`--dependency afterany`), or `--compare-only`,
* inputs and environment are checked before anything is queued; on failures you are asked
  whether to submit anyway,
* every job activates the environment itself and repeats the check on its node,
* logs in `<outdir>/slurm_logs/` (`%x_%A_%a.out/.err`), `--dry-run` prints the job script.

## Re-running

Finished modules are skipped (`log_cluster_jobs/<type>/.<module>.done`). After changing parameters, re-run
the affected module **and everything after it** with `-f`:

| change | re-run |
| --- | --- |
| adapter, quality, `MINLEN`/`MAXLEN` | `-m qc,contam,split,human,nonhuman,assembly,novel,fragments,bins -f` (all except prep) |
| host rules (`HOST_*`) | `-m split,human,nonhuman,assembly,novel,fragments,bins -f` |
| Kraken/Bracken settings | `-m nonhuman -f` |
| fragment/bin parameters | `-m fragments,bins -f` |

`prep` depends only on contigs and annotation and never needs to be repeated for read-level changes.

After an update of the counting (not of the mapping), the existing BAM can be reused:

```bash
ASSEMBLY_RECOUNT_ONLY=1 ./srna_pipeline.sh -s samples.tsv -o out -S BXXX1 -m assembly,novel -f
```

The log then says `ASSEMBLY_RECOUNT_ONLY=1: keeping the existing mapping`; if the BAM or
`log_cluster_jobs/<type>/assembly/bowtie.log` is missing, it warns and maps again. Several parameters go in
front of the command separated by spaces, see
[Setting parameters for a single run](configuration.md#setting-parameters-for-a-single-run).

## Failures

A failed module prints `FAILED (exit n)` and the last 25 lines of its most recent tool log; the remaining
modules of that row are not run, other rows continue. Shared resources (contigs, indices) are protected by
lock directories (`*.lock`); a lock left behind by a killed job can be removed if no job is running.

## Per-sample report

```bash
./srna_pipeline.sh -s samples.tsv -o /path/out --report              # all samples
./srna_pipeline.sh -s samples.tsv -o /path/out --report -S BXXX1     # one sample
./srna_pipeline.sh -s samples.tsv -o /path/out --compare --report    # comparison, then report
```

`scripts/srna_report.py` reads only the summary tables of the modules and writes
`<sample>/Report/<sample>.report.pdf` (A4 landscape, all sRNA types of the sample side by side) and, with
`REPORT_HTML=1` (default), the same pages as a self-contained `<sample>.report.html` (images embedded,
page menu; ~1.5 MB) for quick checks in a browser:

| page | content |
| --- | --- |
| overview | key numbers raw → trimmed → human / non-human → classes → assembly; [read fate](definitions.md#read-fate); inputs |
| QC | trimming, adapter dimers, contaminants, most frequent sequences |
| human | human share, assembly rescue, length distribution, classes, top miRNAs, miRNA origin from marker miRNAs (intestinal/epithelial, blood, immune, muscle) |
| non-human | Bracken genus and family, Rfam types, tRNA by amino acid |
| species | deepest Kraken2 rank per rRNA read, top Bracken species with directly assigned reads and distinct minimizers, weak calls flagged |
| assembly | reads per feature category, tRNA / rRNA fragment classes and lengths |
| QC verdict | PASS / WARN / FAIL per library with the triggered QC rules, FastQC module grid (raw / trimmed) |
| RNA integrity | wet-lab QC of the extraction (`RNA_QC`: RIN, Qubit, A260 ratios), integrity indicators per fraction (tRNA end definition and CCA, end entropy, rRNA spread along 16S/23S, read lengths) and 16S/23S fragment profiles. RIN is measured on the long-RNA fraction of the same extraction and is of limited meaning for cell-free RNA |
| assembly mapping | mapping rate, reads per feature category (sense / antisense), mapped vs unmapped read lengths, top protein-coding genes, top KEGG pathways |
| ncRNAs and regulatory elements | Bakta ncRNA / regulatory_region / tmRNA features by Rfam type (riboswitch, leader, sRNA, ribozyme, …) and top elements, sense / antisense reads |
| CRISPR | top CRISPR arrays (taxon, Cas type, reads and transcribed strand), reads per spacer of the top array, crRNA read ends, read lengths |
| CRISPR loci | locus view of the two top arrays: Cas genes, repeats / spacers and stranded sRNA coverage per fraction (full locus and array only) |
| novel | loci by Rfam class (rRNA / tRNA / tmRNA / other Rfam / no hit), read support, other Rfam families, top candidates other than rRNA/tRNA/tmRNA with RNAfold MFE |
| bins / omics | Spearman table (all pairs, bin bases, features within type), unbinned share, scatter plots, top bins |

Every page starts with a short biological explanation of what it shows; a last page lists the miRNA marker sets and their references. Pages whose inputs are missing (module not run) are skipped. The report also writes MultiQC
custom-content tables to `<sample>/Report/multiqc/` (`*_mqc.tsv`: general statistics, read fate,
Bracken genera), so that `multiqc <outdir>` gives an overview over all samples. The SLURM wrapper runs
`--compare --report` in its final job. Needs `matplotlib` in the environment (`ENV_REPORT`, default
`ENV_DEFAULT`).

### Human miRNA origin markers

The marker groups of the human page come from `scripts/mirna_origin_markers.tsv` (or the file in
`REPORT_MIRNA_MARKERS`), one tab-separated row per group:

```text
group                    color    mirnas                                   references
intestinal / epithelial  #009E73  miR-192-5p,miR-194-5p,miR-215-5p,...     Landgraf P et al. (2007) ... | Ludwig N et al. (2016) ...
blood                    #D55E00  miR-451a,miR-486-5p,miR-144-3p,...       Pritchard CC et al. (2012) ...
_fecal miRNA origin                                                        Liu S et al. (2016) ...
```

* `mirnas`: mature names without `hsa-` (miRBase naming); exact match, one group per miRNA;
* `references`: full citations separated by ` | `; they are numbered on the last report page and cited as
  `[n]` in the legend of the human page;
* rows whose group starts with `_` carry only references (shown on the reference page);
* add groups (e.g. liver) or change lists without touching the code; colors are optional.

## Cohort summary (all samples)

`./srna_pipeline.sh -s samples.tsv -o imp --cohort` — see [Cohort analysis](cohort.md).

## Automated Methods

`--report` writes `Report/<sample>.methods.<style>.md` (sample processing) and `--cohort` writes
`Cohort/cohort.methods.<style>.md` (processing + cohort analysis), in the styles of `METHODS_STYLE`:

| style | layout |
|---|---|
| `nature` | Methods with run-in subsection headings, numbered references (superscript, order of first citation), *Data availability* and *Code availability* |
| `star` | STAR Protocols (Cell Press): summary, before you begin, **key resources table**, step-by-step method details (one step per module with its command), expected outcomes, quantification and statistical analysis, limitations, troubleshooting, resource availability; numbered Cell Press references |

The Methods in `METHODS_REPORT_STYLE` (default **STAR Protocols**) are also added as the last pages of
`<sample>.report.pdf` and `cohort.report.pdf` (two-column text pages; selectable text in the HTML reports).

The text is built from the run itself:

- **parameters** — the resolved configuration of the run (adapter, mismatches, thresholds, QC rules,
  cohort settings such as block, models, NMDS, permutations); only steps that were enabled are described;
- **tool versions** — read from `conda-meta/` of the environments (installed versions), else from the pinned
  versions in `envs/*.yml`;
- **databases and upstream workflow** (genome, miRBase, GtRNAdb, Rfam, Kraken 2 database, IMP3, annotation,
  binning) — from `scripts/methods_resources.tsv`: set names, versions and identifiers there **once** for
  your installation;
- **cohort numbers** — samples, libraries per fraction and group, excluded libraries;
- **references** — Nature or Cell Press format with DOIs; the repository and documentation URLs are filled in.

Everything that could not be determined is written as `[CHECK: …]` — search for it before submission.
To regenerate without rerunning anything:

```bash
python3 scripts/srna_methods.py --mode all --style nature --config srna_pipeline.config.sh --outdir <outdir> \
  --cohort-dir <outdir>/Cohort --envs <env1>,<env2> --out methods.md [--docx]
```

(`--config` reads the parameters from the config file, with environment overrides; `--params` takes a `NAME<TAB>VALUE` file instead.)

## Packing results for a check

```bash
./srna_pack_summary.sh -o /path/out                    # all samples
./srna_pack_summary.sh -o /path/out BXXX1 BXXX2        # selected samples
./srna_pack_summary.sh -o /path/out --clean-legacy     # also list old srnaMapper/CMC files (dry run)
./srna_pack_summary.sh -o /path/out --clean-legacy --yes   # ... and delete them
```

Packs the summary tables of all modules (Stats, human, nonhuman, assembly, novel, fragments, bins, crispr,
comparison), the reports, `inputs.tsv`, the logs with numbers (bowtie, Kraken/Bracken, compare, report) and a
`module_status.txt` (`.done` stamps, warnings / errors). BAM, FASTQ and SAM files are never packed; tables
larger than `--max-size` (5M) are listed in `skipped_large_files.tsv`.
