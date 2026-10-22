# Read definitions

This page defines the read categories used in all output tables. Numbers in `read_stats.tsv`, the human and non-human summaries and all count tables refer to these definitions.

## Raw read

A read of the input FASTQ as delivered by the sequencing facility (e.g. 52 nt for the expobiome runs).

## Adapter dimer

A raw read that **starts** with the first 12 nt of the 3′ adapter (`ADAPTER_3P`), i.e. a ligation product
without insert. Counted by `contam`; removed by trimming because it becomes shorter than `MINLEN`.

## Trimmed read

A raw read after cutadapt that

1. had the 3′ adapter removed if present (≥ 3 nt overlap, up to 10 % errors),
2. had low-quality 3′ bases removed (`-q TRIM_QUALITY`, default 15),
3. optionally had poly-A tails and terminal Ns removed (`CUTADAPT_EXTRA="--poly-a --trim-n"`),
4. is **between `MINLEN` and `MAXLEN` nt** long afterwards (default 18–50; 18–100 for the 52-nt runs, so
   reads without adapter are kept).

All later steps start from the trimmed reads.

## Human read

A trimmed read is **human** if *all* of the following hold:

| # | rule | parameter (default) |
| --- | --- | --- |
| 1 | it is at least `HOST_MIN_LEN` nt long | `HOST_MIN_LEN=20` |
| 2 | it aligns **end-to-end** (no soft-clipping) on either strand to the human genome T2T-CHM13v2 (incl. chrM), or — if not — to the human transcriptome | `HOST_GENOME_FA`, `HOST_TX_FA` |
| 3 | with **0 mismatches** (default at every length; reads ≥ `HOST_SHORT_LEN` may be given `HOST_MISMATCHES` > 0) | `HOST_SHORT_LEN=22`, `HOST_SHORT_MISMATCHES=0`, `HOST_MISMATCHES=0` |
| 4 | it does **not** align to the sample's own assembly with the same or fewer mismatches than to the human reference | `HOST_ASSEMBLY_RESCUE=1` |

Bowtie 1 reports one best alignment per read (`-k 1 --best`); the mismatch count of that alignment
(`NM`) is compared with the best alignment against the contigs (`-v HOST_MISMATCHES -k 1 --best`).

`*.human.fastq.gz` contains all human reads; `*.human.rmDup.fastq.gz` the same reads with duplicate read
names removed (input of the `human` module).

### Why the length rules?

A read that is *random* with respect to the human genome still matches it somewhere by chance if it is
short. The expected number of chance hits of a random read of length *L* with up to *v* mismatches in a
genome of size *G* (both strands) is

$$
E(L, v) = 2G \cdot \frac{\sum_{k=0}^{v} \binom{L}{k}\, 3^k}{4^L}
$$

which for the 3.1 Gb human genome gives

| read length | E, 0 mismatches | E, 1 mismatch |
| --- | --- | --- |
| 18 nt | 0.09 | 5 |
| 19 nt | 0.02 | 1.3 |
| 20 nt | 0.006 | 0.34 |
| 22 nt | 0.0004 | 0.02 |
| 25 nt | 0.000006 | 0.0004 |

With one mismatch, almost every 18–19-nt microbial read would map to the human genome. Exact matching
and no host assignment below 20 nt keep the expected chance hits per read at or below 0.006 (20 nt) and
below 0.0004 from 22 nt. Allowing one mismatch from 22 nt (`HOST_MISMATCHES=1`) raises this to 0.02 and,
in BXXX2 sRNA-Int, roughly doubled the number of 22–35-nt "human" reads. In the BXXX2 sRNA-Int
library, 86 % of the reads counted as human with the original rule (`-v 1` from 18 nt) were 18–21 nt long
and almost none were classified as human small RNA (see [QC](qc.md#false-human-reads)).

### Why the assembly check?

Bacterial rRNA and tRNA fragments are conserved and can be identical to human — in particular
mitochondrial — rRNA/tRNA, and bacterial 4.5S SRP RNA is related to human 7SL RNA. A read that matches the
sample's own contigs at least as well as the human reference is therefore considered microbial
(`host_rescued_by_assembly` in `read_stats.tsv`). The check relies on the assembly: a microbial sequence
missing from the contigs cannot be rescued.

## Non-human read

Every trimmed read that is not human:

* reads shorter than `HOST_MIN_LEN` (`below_host_min_len_20nt`),
* reads without an alignment to the human genome or transcriptome under the rules above,
* human hits rescued by the assembly check.

`*.nonhuman.fastq.gz` is the input of `nonhuman`, `assembly` and `novel`. "Non-human" means *not
attributable to the human host*; it is not a positive identification as microbial (see the classes in
[Classification](classification.md#non-human-reads)).

## Unknown

| where | meaning |
| --- | --- |
| human `unknown_after_all` | human read (genome/transcriptome hit) that matches none of the human miRNA, tRNA, rRNA or other-ncRNA references — e.g. mRNA/intron/repeat fragments, or residual chance hits |
| non-human `unknown` | non-human non-rRNA sequence that matches none of tRNA, miRNA-like (18–26 nt) or Rfam classes |
| assembly `unmapped` | non-human read that bowtie (`-v ASSEMBLY_MISMATCHES`) cannot place on the contigs |

## Read fate

The **read fate** of a library is what happens to each of its raw reads: every raw read ends up in exactly
one category, so the categories add up to the raw reads (100 %). It is the stacked bar on the overview page of
the per-sample report and the `srna_read_fate` table for MultiQC.

| category | reads | source |
| --- | --- | --- |
| removed in trimming | raw − trimmed: adapter dimers, reads shorter than `MINLEN` or longer than `MAXLEN` after trimming | cutadapt |
| human | [human reads](#human-read) | split |
| rRNA | non-human rRNA (RiboDetector) | nonhuman |
| tRNA | non-human, non-rRNA reads matching tRNAs | nonhuman |
| Rfam | non-human, non-rRNA reads matching other Rfam families | nonhuman |
| miRNA-like | non-human, non-rRNA 18–26 nt reads matching miRNAs | nonhuman |
| unknown | non-human reads in none of these classes | nonhuman |

Reads shorter than `HOST_MIN_LEN` are not tested against the human references and therefore count as non-human
(they appear in the non-human classes, not as a category of their own). Fractions are given relative to the raw
reads in the read fate, but relative to the trimmed or non-human reads in the class tables — the reference is
always stated. The *cutadapt* outcome (reads with adapter, too short, too long, passing) is a separate QC
table, not the read fate.

## Count types

| count | meaning | used in |
| --- | --- | --- |
| raw count | number of reads (collapsed sequences of the non-human classification expanded by their count) | all summaries |
| fractional (weighted) count | a read with *n* equally good hits adds 1/*n* to each hit (`NH` tag of the assembly BAM) | human classes, assembly, novel, fragments, bins |
| CPM | counts per million within one source | compare |
