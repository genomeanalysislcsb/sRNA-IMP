#!/usr/bin/env bash
# smetaseq_reimpl.sh — faithful re-implementation of sMETASeq (Mjelle et al. 2020, iScience 23:101131;
# https://github.com/MjelleLab/sMETASeq, scripts as of commit 8a95005, Feb 2023) for the benchmark.
#
# The published scripts (sMETASeq_NEBNext.sh / _TruSeq.sh / _NEXTFlex.sh) run, per sample:
#   cutadapt -a <adapter> -m 1                      -> FILE.cut.fastq
#   bowtie2 -p10 -k10 -x <GRCh38 no_alt> -q FILE.cut_trimmed.fastq --un FILE.Unmapped.cut.fastq -S FILE.mapped.cut.sam
#   htseq-count -i Name -t transcript SAM <RNAcentral GRCh38 gff3, chr-prefixed>  -> FILE.RNACentral.csv
#   htseq-count -i Name -t miRNA      SAM <miRBase hsa.gff3>                      -> FILE.miRBase.csv
#   unmapped FASTQ -> FASTA; kraken2 --db <standard> --report FILE.kraken_report.csv --use-mpa-style
# They do not run as published: cutadapt writes FILE.cut.fastq, bowtie2 reads FILE.cut_trimmed.fastq. This script
# runs the same commands with that file name fixed (deviation 1); options and databases are unchanged.
#
# Usage: smetaseq_reimpl.sh -i reads.fastq.gz -o out/PREFIX -x <bowtie2 index> -m hsa.gff3 -r rnacentral.chr.gff3
#                           -k <kraken2 db> [-a ADAPTER] [-t threads]
#   -a  3' adapter [AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC = sMETASeq_NEBNext.sh; TruSeq small RNA: TGGAATTCTCGGGTGCCAAGG]
# Needs cutadapt, bowtie2, samtools, htseq-count, kraken2 in PATH (bioconda: cutadapt bowtie2 samtools htseq kraken2).
# Deviation 2: htseq-count reads only the aligned records (samtools view -b -F 4), because current HTSeq fails on
# the unaligned records of 1-nt reads; feature counts are unchanged.
# Databases as in the sMETASeq README: GRCh38 no_alt bowtie2 index (NCBI), miRBase hsa.gff3, RNAcentral
# homo_sapiens.GRCh38.gff3 with "chr" prefix added, a kraken2 standard database.
set -euo pipefail
usage() { sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }
IN="" OUT="" IDX="" MIRBASE="" RNAC="" KDB="" ADAPTER=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC T=10
while getopts "i:o:x:m:r:k:a:t:h" o; do
  case $o in
    i) IN=$OPTARG ;; o) OUT=$OPTARG ;; x) IDX=$OPTARG ;; m) MIRBASE=$OPTARG ;; r) RNAC=$OPTARG ;;
    k) KDB=$OPTARG ;; a) ADAPTER=$OPTARG ;; t) T=$OPTARG ;; h) usage 0 ;; *) usage 1 ;;
  esac
done
[[ -n "$IN" && -n "$OUT" && -n "$IDX" && -n "$MIRBASE" && -n "$RNAC" && -n "$KDB" ]] || usage 1
for t in cutadapt bowtie2 samtools htseq-count kraken2; do command -v "$t" >/dev/null || { echo "Error: $t not in PATH" >&2; exit 1; }; done
mkdir -p "$(dirname "$OUT")"
F=$OUT
echo "[smetaseq] $(date '+%F %T') $IN -> $F.*" >&2

cutadapt -a "$ADAPTER" -m 1 -o "$F.cut.fastq" "$IN" > "$F.cutadapt.log"
# (published scripts read $FILE.cut_trimmed.fastq here, which is never written)
bowtie2 -p "$T" -k10 -x "$IDX" -q "$F.cut.fastq" --un "$F.Unmapped.cut.fastq" -S "$F.mapped.cut.sam" 2> "$F.bowtie2.log"
# HTSeq 2.1.2 / pysam 0.24.1 stop with a UnicodeDecodeError on the unaligned records of 1-nt reads that
# cutadapt -m 1 keeps (bowtie2 writes them with YF:Z:LN). htseq-count counts only aligned reads, so it gets the
# aligned records only (samtools view -F 4): the feature counts are unchanged, only __not_aligned is not filled.
samtools view -@ "$T" -b -F 4 -o "$F.mapped.cut.bam" "$F.mapped.cut.sam"
htseq-count -f bam -i Name -t transcript "$F.mapped.cut.bam" "$RNAC" > "$F.RNACentral.csv" 2> "$F.htseq_rnacentral.log"
htseq-count -f bam -i Name -t miRNA "$F.mapped.cut.bam" "$MIRBASE" > "$F.miRBase.csv" 2> "$F.htseq_mirbase.log"
sed -n '1~4s/^@/>/p;2~4p' "$F.Unmapped.cut.fastq" > "$F.Unmapped.cut.fasta"
kraken2 --db "$KDB" --threads "$T" --report "$F.kraken_report.csv" --use-mpa-style --output "$F.kraken_mpa.csv" \
  "$F.Unmapped.cut.fasta" 2> "$F.kraken2.log"
echo "[smetaseq] $(date '+%F %T') done" >&2
