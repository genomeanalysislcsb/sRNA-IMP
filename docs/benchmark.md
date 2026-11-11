# Simulation benchmark

`scripts/srna_simulate.py` builds a stool-like small RNA library in which the origin of every read is known, from a
real sample assembly and human small RNA references. The pipeline runs on it like on a real library, and
`scripts/srna_sim_eval.py` scores its outputs against the truth.

## What is simulated

| source | how | truth for |
| --- | --- | --- |
| microbial, assembled | fragments of annotated tRNA, rRNA, tmRNA, ncRNA, regulatory elements, CRISPR, CDS (sense and antisense) and intergenic stretches; contig abundance from metaG contig counts, else log-normal per bin | classification, feature counts, antisense |
| microbial, unassembled | the same from contigs held out of the assembly given to the pipeline (`--holdout-frac`, default 10 %) | reads without a contig |
| hidden loci | ncRNA and regulatory features removed from the annotation given to the pipeline (`--hide-frac`, default 30 %) | novel-locus recall |
| human | mature miRNAs with 3′ isomiRs, tRNA halves and fragments, rRNA and other ncRNA fragments (`--host-frac`, default 2 %) | host vs microbial split, miRNA counts |
| adapter dimers | 0–5 nt inserts (`--dimer-frac`, default 3 %) | trimming |

tRNA fragments follow the usual types and sizes (5′ halves, 3′ halves with CCA, 3′ and 5′ tRFs, internal
fragments, mature tRNAs); reads are insert + 3′ adapter + random fill to `--read-len` (default 51) with
position-dependent substitution errors (`--no-errors` to switch off). Every read's origin, class, feature, contig,
bin, strand, coordinates and fragment type go to `<out>.truth.tsv.gz`.

## Run

Human reference FASTAs from the pipeline's bowtie indices (once):

```bash
for x in MIRNA TRNA RRNA NCRNA_OTHER; do
  idx=$(bash -c 'source srna_pipeline.config.sh >/dev/null 2>&1; v=IDX_HSA_'$x'; echo ${!v}')
  bowtie-inspect -a 0 "$idx" > bench/hsa_${x,,}.fa
done
```

Simulate from one sample (here an IMP3 metaG + metaT sample):

```bash
I=/path/to/imp/bxxx1
python3 scripts/srna_simulate.py \
  --assembly $I/Assembly/mgmt.assembly.merged.fa.gz \
  --gff $I/Analysis/annotation/annotation_CDS_RNA_hmms.gff.gz \
  --contig2bin $I/Binning/MetaWrap/bins.contigs.gz \
  --contig-counts $I/Stats/mg/mgmt.assembly.contig_read_count.txt.gz \
  --human-mirna bench/hsa_mirna.fa --human-trna bench/hsa_trna.fa \
  --human-rrna bench/hsa_rrna.fa --human-other bench/hsa_ncrna_other.fa \
  --n-reads 5000000 --host-frac 0.02 --seed 1 --sample SIM01 --out bench/sim/SIM01
```

`<out>.sample_row.tsv` is a sample-sheet row that points at the simulated reads and the reduced assembly and
annotation. Run the pipeline on it as usual (same config as for real data), then score:

```bash
(printf 'sample_id\ttype\traw_fastq\tassembly\tgff\tcontig2bin\tgtdbtk\tmetaT_features\tmetaT_contigs\tmetaG_features\tmetaG_contigs\n';
 cat bench/sim/SIM*.sample_row.tsv) > bench/samples.sim.tsv
./srna_pipeline.sh -s bench/samples.sim.tsv -o bench/run -c srna_pipeline.config.sh -r 1

python3 scripts/srna_sim_eval.py --truth bench/sim/SIM01.truth.tsv.gz --hidden bench/sim/SIM01.hidden_loci.tsv \
  --outdir bench/run --sample SIM01 --out bench/sim/SIM01.eval
```

## Metrics (`<out>.summary.tsv`)

| section | metric |
| --- | --- |
| qc | share of adapter dimers removed; share of inserts ≥ `MINLEN` kept |
| split | precision and recall of the human calls; microbial reads called human, per RNA class |
| human | human miRNA counts (genome path + rescued miRNA-like reads) vs truth: Spearman, true miRNAs detected, false miRNAs |
| assembly | sense reads per feature vs truth (Spearman, totals); antisense reads on CDS |
| novel | hidden loci recovered by the novel candidates (same contig and strand, overlapping), also for loci with ≥ 10 reads |

Per-miRNA, per-feature and per-locus tables: `<out>.mirna.tsv`, `<out>.features.tsv`, `<out>.novel.tsv`.

For a parameter sweep, simulate once and run the pipeline into separate output folders with different configs
(e.g. `HOST_MIN_LEN`, `HOST_MISMATCHES`, `HOST_ASSEMBLY_RESCUE`), then score each run against the same truth.

## Comparison with sMETASeq

`scripts/smetaseq_reimpl.sh` re-implements [sMETASeq](https://github.com/MjelleLab/sMETASeq) (Mjelle et al. 2020,
iScience 23:101131): cutadapt (adapter only, `-m 1`), bowtie2 `-k10` to GRCh38, htseq-count on miRBase and
RNAcentral GFFs, kraken2 on the unmapped reads. The commands and options are those of the published scripts; the
only changes are the input file name of bowtie2, which in the published scripts does not match the cutadapt output, and
htseq-count reading only the aligned records (`samtools view -b -F 4`), because current HTSeq fails on the
unaligned records of the 1-nt reads that `cutadapt -m 1` keeps; the feature counts are unchanged.
Databases as in the sMETASeq README (GRCh38 no_alt bowtie2 index, miRBase `hsa.gff3`, RNAcentral GRCh38 GFF3 with
`chr` prefix, a kraken2 standard database); tools from bioconda (`cutadapt bowtie2 samtools htseq kraken2`).

```bash
scripts/smetaseq_reimpl.sh -i bench/sim/SIM01.fastq.gz -o bench/smetaseq/SIM01 -t 16 \
  -x GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -m hsa.gff3 \
  -r homo_sapiens.GRCh38.chr.gff3 -k k2_standard
python3 scripts/srna_sim_eval.py --truth bench/sim/SIM01.truth.tsv.gz --outdir . --sample SIM01 \
  --smetaseq bench/smetaseq/SIM01 --out bench/sim/SIM01.smetaseq.eval
```

With `--smetaseq` the evaluation scores the trimmed reads, the host calls (reads mapped to GRCh38) and the miRNA
counts (`PREFIX.miRBase.csv`) against the same truth; sMETASeq has no assembly-based outputs.
