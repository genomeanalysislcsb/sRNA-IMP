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
