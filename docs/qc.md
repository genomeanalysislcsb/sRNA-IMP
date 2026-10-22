# Quality control & troubleshooting

## What to check for every library

| check | file | expected |
| --- | --- | --- |
| reads with adapter | `Stats/<type>/<s>.<type>.cutadapt.log` → "Reads with adapters" | well above 50 % for short-insert libraries |
| trimmed vs raw | `read_stats.tsv` | trimmed clearly **below** raw (adapter dimers, inserts < 18 nt removed) |
| loss by trimming | cutadapt log → "too short" / "too long" | mostly *too short*; *too long* ≈ 0 if `MAXLEN` ≥ read length |
| adapter dimers, complexity | `qc_summary.tsv` | low dimer fraction; top-10 sequences not dominating |
| human length distribution | `<s>.<type>.human_length_dist.tsv` | no steep decline from 18–20 nt; peaks at ~22 nt (miRNA) and 30–35 nt (tRNA halves) if real human RNA is present |
| human classes | human summary | miRNA / tRNA / rRNA / other present, `unknown_after_all` not dominant |
| rRNA share, taxonomy | `nonhuman.summary.tsv`, Bracken tables | plausible gut taxa |

## Wrong adapter

**Symptom:** almost no reads have an adapter removed (BXXX2: 1.2 % and 0.2 %), "trimmed" equals "raw",
very few reads map anywhere — or, with `MAXLEN` below the read length, almost all reads are discarded as
*too long* (BXXX2: 96 % lost with `MAXLEN=50` and 52-nt reads).

**Check** which adapter the raw reads contain:

```bash
R=raw.fastq.gz
for a in TGGAATTCTCGG AGATCGGAAGAGC AACTGTAGGCAC CTGTAGGCACCA; do
  echo -n "$a  "; zcat $R | head -400000 | awk 'NR%4==2' | grep -c $a
done
```

| sequence | kit |
| --- | --- |
| `TGGAATTCTCGG` | Illumina TruSeq Small RNA |
| `AGATCGGAAGAGC` | Illumina TruSeq / NEBNext (expobiome: 85 % and 37 % of reads) |
| `AACTGTAGGCAC` | QIAseq miRNA (UMIs!) |
| `CTGTAGGCACCA` | Illumina v1.5 small RNA |

A count of `AAAAAAAAAAGATCGGAAGAGC` close to the `AGATCGGAAGAGC` count would indicate poly-A-tailing
chemistry (e.g. SMARTer smRNA-seq: also `-u 3`); for expobiome it is 0.

**Fix:** set `ADAPTER_3P`, then re-run everything except `prep` with `-f`. Untrimmed adapters make reads
unmappable end-to-end, so *all* downstream results of such a run are invalid.

## False human reads

**Symptom:** many human reads, almost none of them human miRNA/tRNA/rRNA, `unknown_after_all` dominant,
human length distribution falling steeply from 18 nt.

**Cause:** short microbial reads match the 3 Gb human genome by chance (an 18-mer with one mismatch is
expected to match about 5 times at random, see [Read definitions](definitions.md#why-the-length-rules)), and
conserved bacterial rRNA/tRNA/SRP fragments resemble human (mitochondrial) sequences.

BXXX2 sRNA-Int, human reads under the successive rules:

| host rule | human reads | share of trimmed |
| --- | --- | --- |
| `-v 1` from 18 nt (original) | 385,565 | 5.6 % |
| exact below 22 nt + assembly check | 68,833 | 1.0 % |
| + no host assignment below 20 nt (`HOST_MIN_LEN=20`) | ≈ 54,000 | ≈ 0.8 % |
| + exact at every length (`HOST_MISMATCHES=0`, **current default**) | **30,560** | **0.44 %** |

The intermediate rows were computed while the BXXX2 reads were accidentally run against the BXXX1
assembly; only the assembly check depends on that, so they illustrate the trend, not exact values.

With the first rule, 86 % of the human reads were 18–21 nt, and the human classes were 1 miRNA, 80 tRNA,
678 rRNA, 59,683 other ncRNA and 325,123 unknown reads. With the second rule the length distribution still
fell from 9,376 (18 nt) to 3,384 (21 nt) reads and jumped to 5,743 at 22 nt, where one mismatch is allowed
again — hence `HOST_MIN_LEN=20` and `HOST_MISMATCHES=0`. With both, the distribution falls smoothly from
4,008 reads (20 nt) through 2,610 (22 nt) to 1,140 (30 nt), without the break at 22 nt and without a miRNA
or tRNA-half peak: the intracellular library contains essentially no human small RNA.

Final `read_stats.tsv` of BXXX2 sRNA-Int:

| step | reads |
| --- | --- |
| raw | 8,918,925 |
| trimmed | 6,909,507 |
| below_host_min_len_20nt | 502,410 |
| host_genome_mapped | 32,802 |
| host_tx_mapped | 14 |
| host_rescued_by_assembly | 2,256 |
| human_total | 30,560 |
| nonhuman | 6,878,947 |

The 502,410 reads below 20 nt (7.3 % of trimmed) are always non-human: at that length a genuine human read
cannot be told apart from a chance match. With exact matching the assembly check moves fewer reads (7 % of
host hits instead of 25 %) because most conserved bacterial fragments differed from human by one mismatch;
it still catches identical sequences.

### Intracellular vs extracellular (BXXX2, default settings)

| step | sRNA-Int | sRNA-Ext |
| --- | --- | --- |
| raw | 8,918,925 | 10,604,499 |
| trimmed | 6,909,507 (77 %) | 9,500,172 (90 %) |
| below_host_min_len_20nt | 502,410 | 251,476 |
| host_genome_mapped | 32,802 | 4,333 |
| host_tx_mapped | 14 | 2 |
| host_rescued_by_assembly | 2,256 (7 % of host hits) | 1,668 (38 % of host hits) |
| **human_total** | **30,560 (0.44 %)** | **2,667 (0.03 %)** |
| nonhuman | 6,878,947 | 9,497,505 |

* Both fractions are almost entirely microbial; the extracellular fraction contains even less human RNA
  than the intracellular one.
* Neither human length distribution shows a small-RNA signature: sRNA-Ext falls steadily from 339 reads
  at 20 nt to about 40 nt, with no peak at 22 nt (miRNA) or 30–35 nt (tRNA halves). The 228 reads at 51 nt
  are reads without adapter (insert ≥ read length).
* In sRNA-Ext a third of the host hits also match the sample's contigs exactly — conserved rRNA/tRNA
  fragments, which dominate extracellular RNA. The assembly check is therefore essential for this fraction.
* `HOST_MIN_LEN=20` loses no human miRNAs: the `miRNA_like` class of the non-human reads (which compares
  against eukaryotic miRNAs incl. human) contains **no** `hsa-` miRNA in BXXX2 sRNA-Ext.

!!! tip "Check for every new data set"
    ```bash
    grep -i "hsa-" <outdir>/<sample>/Analysis/<type>/nonhuman/mirna_like/<sample>.mirna_like.counts.tsv
    ```
    Clear `hsa-` counts there would mean human miRNA isoforms < 20 nt are being counted as non-human; then
    lower `HOST_MIN_LEN` to 19 (an exact 19-nt match happens by chance for about 2 % of random reads).

## Environment problems

`--check-env` names every missing tool per module and environment.

| message | cause / fix |
| --- | --- |
| `MISSING fastqc` although installed | another copy of the environment is active — check `ENV_DEFAULT` against `micromamba env list` |
| `MISSING env <name> … does it exist?` | wrong path in `ENV_<MODULE>` |
| `cannot activate … no micromamba/conda found` | set `MAMBA_BIN` (output of `which micromamba` on a compute node), `CONDA_SH` or `ENV_MODULES` |
| FastQC fails after 0 s | no Java in the environment, or too little memory for `-t` threads |

## A module is skipped although parameters changed

Finished modules leave `.<module>.done`; the log shows `already done … skipping`. Re-run with `-f`
(see [Re-running](running.md#re-running)). Also check that the job used the intended config
(`config=` in the second log line).
