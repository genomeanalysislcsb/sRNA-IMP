# Configuration parameters

All parameters are defined in `srna_pipeline.config.sh`, which is sourced by the pipeline. Variables written as `${VAR:-default}` can also be overridden from the environment.

### Setting parameters for a single run

Instead of editing the config, put one or more `NAME=value` pairs **in front of the command, separated by
spaces** (no semicolons). They apply to this run only:

```bash
ASSEMBLY_RECOUNT_ONLY=1 RRNA_EXCLUDE_RFAM=1 ./srna_pipeline.sh -s samples.tsv -o out -S BXXX1 -m assembly,novel -f
```

- `A=1; B=2; ./srna_pipeline.sh …` does **not** work: with semicolons the values are plain shell variables
  that the pipeline does not see.
- `export A=1 B=2` on a separate line works too, but the values stay set for every later run in that shell
  (`unset A B` afterwards).
- For SLURM, the variables are passed on to the jobs when set in front of `submit_srna_pipeline.sh` or
  exported before it (`sbatch` exports the environment by default).
- This works for all variables written as `${VAR-default}` / `${VAR:-default}` in the config; a value set
  without a default (`VAR=value`) in the config always wins.

## General

| parameter | default | meaning |
| --- | --- | --- |
| `SCRIPTS_DIR` | `/scratch/users/pmay/sRNAs/scripts` | Python helper scripts |
| `TMP_ROOT` | `$SLURM_TMPDIR` / `$TMPDIR` / `/tmp` | temporary files |
| `STAGE_MODE` | `copy` | `copy` or `symlink` raw reads |
| `MINLEN` | 18 | minimum read length after trimming |
| `MIN_LIBRARY_READS` | 10000 | libraries with fewer trimmed reads (failed libraries) are not analysed beyond `split`: later modules are skipped and marked done; the QC rules mark them FAIL |
| `MAXLEN` | 50 | maximum read length after trimming (use the read length, e.g. 100) |
| `ASSEMBLY_SUBDIR` | `assembly` | sub-folder of the assembly results (`""` = `Analysis/<type>/`) |
| `KEEP_INTERMEDIATE` | 0 | 1 = keep SAMs, hit tables etc. |

## Environments

| parameter | default |
| --- | --- |
| `ENV_DEFAULT` | `…/envs/gut_srna_pipeline` |
| `ENV_QC`, `ENV_CONTAM`, `ENV_PREP`, `ENV_SPLIT`, `ENV_HUMAN`, `ENV_NONHUMAN`, `ENV_ASSEMBLY`, `ENV_COMPARE`, `ENV_REPORT` | `ENV_DEFAULT` |
| `ENV_NOVEL`, `ENV_FRAGMENTS`, `ENV_BINS` | `…/envs/ncrna_discovery` |
| `RIBO_DETECTOR_CPU` | full path of `ribodetector_cpu` (own env, `envs/ribodetector.yml`) |
| `ENV_MODULES`, `MAMBA_BIN`, `CONDA_SH` | empty — how to find micromamba/conda |

## qc

| parameter | default | meaning |
| --- | --- | --- |
| `ADAPTER_3P` | `AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC` | 3′ adapter (Illumina TruSeq/NEBNext read 1) |
| `ADAPTER_3P_TRUNC` | 10 | also trim adapter copies missing up to N leading nt (`-a ${ADAPTER_3P:N}`); catches adapter/index-primer dimers that start inside the adapter. 0 = off |
| `TRIM_QUALITY` | 15 | cutadapt `-q` |
| `CUTADAPT_EXTRA` | `""` | extra options, e.g. `--poly-a --trim-n` |

## split (human vs non-human)

| parameter | default | meaning |
| --- | --- | --- |
| `HOST_GENOME_FA` | T2T-CHM13v2 genome | human genome |
| `HOST_TX_FA` | T2T-CHM13v2.0 transcripts | human transcriptome (`''` = skip) |
| `HOST_GENOME_IDX`, `HOST_TX_IDX` | `""` | existing bowtie index prefixes |
| `HOST_INDEX_TMPDIR` | `<outdir>/tmp/host_index` | where indices are built if the FASTA directory is read-only |
| `HOST_MIN_LEN` | 20 | reads shorter are never human |
| `HOST_SHORT_LEN` | 22 | reads shorter must match … |
| `HOST_SHORT_MISMATCHES` | 0 | … with this many mismatches |
| `HOST_MISMATCHES` | 0 | mismatches for reads ≥ `HOST_SHORT_LEN` (1 = more tolerant, more chance hits) |
| `HOST_ASSEMBLY_RESCUE` | 1 | human hits matching the contigs as well → non-human |

## human

| parameter | meaning |
| --- | --- |
| `IDX_HSA_MIRNA`, `IDX_HSA_TRNA`, `IDX_HSA_RRNA`, `IDX_HSA_NCRNA_OTHER` | bowtie indices of the human classes |

## nonhuman

| parameter | default | meaning |
| --- | --- | --- |
| `RIBO_DETECTOR_CPU`, `RIBO_CHUNK_SIZE` | path, 256 | RiboDetector |
| `IDX_TRNA_ALL`, `IDX_MIRNA_EUK`, `IDX_RFAM_CLASSES`, `RFAMSEQ_TXT` | | class references (DNA alphabet: indexes without T are reported by the check) |
| `RFAM_FAMILY_TXT` | `…/rfam/family.txt.gz` | Rfam `family.txt` of the same release; adds `rfam_type.counts.tsv` (rRNA, sRNA, snoRNA, riboswitch, …) and rfam_id/description/type to the family table; `""` = off |
| `MIRNA_LIKE_MINLEN`, `MIRNA_LIKE_MAXLEN` | 18, 26 | length window of the miRNA-like search |
| `KRAKEN_DB`, `BRACKEN_DB` | GTDB genome reps 2025-06-09 | rRNA taxonomy |
| `KRAKEN_CONF` | 0.5 | `--confidence` |
| `KRAKEN_MINHITGROUPS` | 3 | `--minimum-hit-groups` |
| `KRAKEN_MINIMIZER_DATA`, `KRAKEN_MEMORY_MAPPING` | 1, 1 | report minimizer columns; memory-map the DB |
| `KRAKEN_REPORT_MINREADS` | 10 | report lines with fewer clade reads are dropped before Bracken |
| `KRAKEN_MIN_DISTINCT_MINIMIZERS` | 30 | species with fewer distinct minimizers are dropped too (needs `KRAKEN_MINIMIZER_DATA=1`); removes hits that rest on one short conserved stretch; 0 = off |
| `KRAKEN_SOFT_FAIL` | 0 | 0 = a failed kraken2 run (e.g. out of memory) fails the nonhuman module, which is then not marked done; 1 = warn and continue without rRNA taxonomy |
| `BRACKEN_READ_LENGTH` | 50 | Bracken `-r` (needs `database50mers.kmer_distrib`) |
| `BRACKEN_LEVELS` | `D P C O F G S` | ranks (D = domain; `K` is mapped to `D`, Kraken2 has no kingdom rank for prokaryotes) |
| `BRACKEN_MIN_READS` | 10 | Bracken `-t` |

## assembly

| parameter | default | meaning |
| --- | --- | --- |
| `ASSEMBLY_MISMATCHES` | 1 | bowtie `-v` for the mapping on the contigs |
| `ASSEMBLY_MAX_HITS` | 50 | at most N equally good hits per read (`-k N --best --strata`); 0 = all hits (`-a`) |
| `ASSEMBLY_RECOUNT_ONLY` | 0 | 1 = keep the existing BAM (and `bowtie.log`) and only recount, e.g. after a counting update: `ASSEMBLY_RECOUNT_ONLY=1 … -m assembly,novel -f` |
| `RRNA_EXCLUDE_RFAM` | 1 | leave reads with an Rfam rRNA hit (by sequence, `rfam/<s>.rfam_rRNA.seqs.gz`) out of the `*_non_rRNA` counts and the novel loci, in addition to the RiboDetector rRNA reads |
| `DB_NORC` | 1 | stranded reads: match only the forward strand of the RNA databases (`bowtie --norc`; human miRNA/tRNA/rRNA/other, microbial tRNA, miRNA, Rfam); `""` = both strands |

## novel

| parameter | default | meaning |
| --- | --- | --- |
| `NOVEL_EXCLUDE_RRNA` | 1 | leave rRNA reads (RiboDetector, and Rfam rRNA hits if `RRNA_EXCLUDE_RFAM=1`) out of the clustering |
| `RRNA_RESCUE` | 1 | cmscan with the Rfam rRNA models on contigs whose CDS have only antisense reads; hits become rRNA genes (`assembly/<s>.rRNA_rescue.bed`) and the assembly counts are redone |
| `RRNA_RESCUE_MIN_ANTISENSE`, `RRNA_RESCUE_MAX_SENSE_RATIO` | 20, 0.1 | which CDS are suspicious |
| `RRNA_RESCUE_CM` | `""` | rRNA subset of `RFAM_CM`, built once; `""` = `<outdir>/.shared/Rfam.rRNA.cm` |
| `NOVEL_MERGEDIST` | 10 | max. gap within a cluster |
| `NOVEL_MINLEN`, `NOVEL_MAXLEN` | 18, 250 | cluster length window |
| `NOVEL_MINREADS` | 10 | min. supporting reads (multi-mappers weighted 1/hits) |
| `CMSCAN_OPTS` | `--rfam --cut_ga` | cmscan options for the candidates |
| `RUN_ARAGORN`, `RUN_TRNASCAN`, `RUN_CMSCAN`, `RUN_RNAFOLD` | 1, 0, 1, 1 | optional tools |
| `RFAM_CM` | `Rfam.cm` | Rfam covariance models |

## contam, fragments, bins, compare

| parameter | default | meaning |
| --- | --- | --- |
| `CONTAM_FA` | empty | contaminant FASTAs (space-separated) |
| `CONTAM_MISMATCHES` | 1 | bowtie `-v` |
| `QC_TOP_SEQS`, `QC_SORT_MEM` | 50, 2G | overrepresented sequences; sort buffer |
| `FRAG_TYPES` | `tRNA,rRNA,tmRNA` | analysed feature types |
| `FRAG_TOL` | 3 | nt tolerance at the ends |
| `FRAG_AC_FROM`, `FRAG_AC_TO` | 0.40, 0.55 | anticodon region (fraction of tRNA length) |
| `FRAG_MIN_OVERLAP` | 0.9 | fraction of the read inside the feature |
| `FRAG_MIN_READS` | 10 | features with per-position profiles |
| `BINS_MIN_OVERLAP` | 0.5 | fraction of the read overlapping a feature |
| `COMPARE_PSEUDOCOUNT` | 1 | added to CPM before log2 ratios |
| `COMPARE_MIN_COUNT` | 5 | min. raw count in at least one of the two sources of a pair |
| `REPORT_HTML` | 1 | `--report` also writes `<sample>.report.html`; `""` = PDF only |
| `METHODS_STYLE` | `nature,star` | automated Methods written by `--report` and `--cohort`: `nature`, `star` (STAR Protocols) or both; `none` = off |
| `METHODS_REPORT_STYLE` | `star` | Methods added as the last pages of the sample and cohort reports (PDF; real text in the HTML): `star`, `nature` or `none` |
| `METHODS_RESOURCES` | `""` | databases and upstream tools (name, version, identifier, references) for the Methods; `""` = `scripts/methods_resources.tsv` |
| `METHODS_DOCX` | `""` | `1` = also `.docx` (needs pandoc) |
| `REPORT_MIRNA_MARKERS` | `""` | TSV of human miRNA origin marker sets with references; `""` = `scripts/mirna_origin_markers.tsv` |
| `RNA_QC` | `""` | wet-lab RNA QC table (columns `sample`, `sample_id`, `type` = sRNA-Int / sRNA-Ext, then any QC columns such as `lRNA_RIN`, DV200, Qubit, A260/280); shown on the RNA integrity page and in the MultiQC tables |
| `QC_RULES` | `""` | QC rules TSV (PASS / WARN / FAIL per library); `""` = `scripts/qc_rules.tsv` |
| `QC_EXCLUDE` | `""` | manual QC decisions (sample, type or `*`, level, reason) |
| `COHORT_BLOCK` / `COHORT_REFERENCE` | `""` | block column: `site` (nested, Expobiome) or `family` (crossed, T1D) — pooled + within-block comparisons, poolability check, block-aware models (old name `COHORT_MATCH_COL`); reference group |
| `COHORT_OUTLIER_Z` | `4` | outlier flag (WARN) for libraries with \|robust z\| above this in key metrics; `0` = off |
| `COHORT_QC_DROP` | `FAIL` | libraries left out of the cohort statistics: `FAIL`, `WARN` (WARN + FAIL) or `none` |
| `MULTIQC` | `multiqc` | MultiQC command run by `--cohort`; `""` = skip |
| `COHORT_METADATA` | `""` | metadata table for `--cohort` (TSV/CSV/XLSX, one row per sample) |
| `COHORT_ID_COL` / `COHORT_GROUP_COL` | `""` / `group` | ID column (default first) and group column |
| `COHORT_GROUP_ORDER` | `""` | group order, e.g. `HC,PD,RA` or `HC,T1D` |
| `COHORT_COVARIATES` | `""` | covariate columns (default: all other columns) |
| `COHORT_PERMUTATIONS` | `999` | PERMANOVA permutations |
| `COMPARE_BIN_BASIS` | `contig,cds,no_rrna` | how reads are summed per bin; the first is the main bin level (`bins.tsv`), the others add `bins_<basis>.tsv` and summary rows `bin_<basis>` (see workflow §11) |

## chimera (optional module, `-m chimera`)

| variable | default | meaning |
|---|---|---|
| `CHIMERA_MIN_ARM` | 18 | minimum arm = exact seed length; reads ≥ 2 × this are tested |
| `CHIMERA_MAX_HITS` | 5 | seeds with more exact hits on the assembly are not used |
| `CHIMERA_ARM_MISMATCHES` | 1 | mismatches allowed while extending an arm |
| `CHIMERA_MAX_GAP` | 2 | read nt covered by neither arm |
| `CHIMERA_COLINEAR_MAX` | 1000 | same contig and strand, 3′ arm ≤ N nt downstream → colinear (not an interaction) |
| `CHIMERA_UTR5` / `CHIMERA_UTR3` | 150 / 100 | UTR windows around CDS |
| `CHIMERA_FLANK` | 20 | contig nt added on both sides of an arm for the duplex energy |
| `CHIMERA_NULL` | 20 | random re-pairings per pair |
| `CHIMERA_MAX_PAIRS` | 5000 | feature pairs tested for complementarity |
| `RNADUPLEX` | – | RNAduplex executable; default: PATH, then `$ENV_NOVEL/bin/RNAduplex` |
| `CHIMERA_MT` | 1 | also analyse the IMP3 metaT BAM (once per sample, `Analysis/metaT/chimera/`) |
| `CHIMERA_MT_BAM` | – | metaT BAM (`<ID>`/`<id>` placeholders); default from the `metaT_features` IMP directory |
| `CHIMERA_MT_STRAND` | auto | `R1-sense`, `R1-antisense` (dUTP), `unstranded` |
| `CHIMERA_MT_MIN_ARM` / `_MIN_MAPQ` / `_EDGE` | 25 / 20 / 500 | split-read part length, MAPQ, contig-end distance of discordant mates |
| `ENV_CHIMERA` | `ENV_DEFAULT` | environment (bowtie, samtools) |

Run after `assembly` (and `novel`, for novel sRNA loci): `./srna_pipeline.sh -s samples.tsv -o out -S BXXX1 -m chimera`.

