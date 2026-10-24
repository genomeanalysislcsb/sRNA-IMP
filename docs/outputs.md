# Output files

## Directory layout

One tree per sample, one sub-tree per sRNA type; `Assembly/` and `Analysis/annotation/` are shared by all
types of a sample.

```text
<outdir>/
├── <sample>/
│   ├── Preprocessing/<type>/     raw, trimmed, human, human.rmDup, nonhuman FASTQs
│   ├── Stats/<type>/             FastQC, cutadapt log, read_stats, QC summary, seqkit stats
│   ├── Assembly/                 contigs (+ bwa and bowtie indices)
│   ├── Analysis/
│   │   ├── annotation/           <sample>.gff, <sample>.bed
│   │   ├── <type>/
│   │   │   ├── human/            mirna/ trna/ rrna/ otherrna/ unknown/
│   │   │   ├── nonhuman/         rrna/ rrna_tax/ collapsed/ trna/ mirna_like/ rfam/ unknown/
│   │   │   ├── crispr/           crRNAs per CRISPR array / spacer
│   │   │   ├── chimera/          chimeric reads (optional module)
│   │   │   ├── assembly/         contig mapping (BAM), feature / category / functional counts
│   │   │   ├── novel/            BAM, clusters, candidates, Rfam, RNAfold
│   │   │   ├── fragments/        tRNA / rRNA fragment classes and profiles, RNA integrity
│   │   │   └── bins/             contig / feature / category / bin / taxon counts
│   │   └── comparison/           --compare: <A>_vs_<B>/ per pair + all_sources tables
│   ├── Report/                   --report: <sample>.report.pdf / .html, multiqc/*_mqc.tsv
│   └── log_cluster_jobs/<type>/  inputs.tsv, .<module>.done, <module>/ tool logs
├── Cohort/                       --cohort: tables/, stats/, cohort.report.pdf (see Cohort analysis)
├── log_cohort/                   cohort.log
├── slurm_logs/                   SLURM stdout/stderr, input checks, job scripts
└── tmp/host_index/               human/contaminant bowtie indices if not writable next to the FASTA
```

## Key files per module

| module | file (under `<outdir>/<sample>/`) | content |
| --- | --- | --- |
| qc | `Preprocessing/<type>/<s>.raw.fastq.gz`, `.trimmed.fastq.gz` | staged raw reads, trimmed reads |
| qc | `Stats/<type>/fastqc/`, `<s>.<type>.cutadapt.log` | FastQC reports, cutadapt summary |
| contam | `Stats/<type>/<s>.<type>.qc_summary.tsv` | adapter dimers, cutadapt outcome, complexity, contaminants |
| contam | `Stats/<type>/<s>.<type>.overrepresented.tsv` | most frequent sequences |
| split | `Preprocessing/<type>/<s>.human.fastq.gz`, `.human.rmDup.fastq.gz`, `.nonhuman.fastq.gz` | read sets by origin |
| split | `Stats/<type>/<s>.<type>.read_stats.tsv` | read numbers per step (below) |
| split | `Stats/<type>/<s>.<type>.human_length_dist.tsv` | length distribution of the human reads |
| human | `Analysis/<type>/human/*/<s>.human.*.counts.tsv` | fractional counts per human reference |
| human | `Analysis/<type>/human/` summary | reads per human class |
| nonhuman | `Analysis/<type>/nonhuman/<s>.nonhuman.summary.tsv` | reads per non-human class |
| nonhuman | `nonhuman/rfam/<s>.rfam_{class,family,type,type_taxon,taxon,accession}.counts.tsv` | Rfam hits per description, family (with id/type), Rfam type, taxon |
| nonhuman | `nonhuman/trna/<s>.trna_{aa,anticodon,species,phylum_class}.counts.tsv` | microbial tRNA hits |
| nonhuman | `nonhuman/rrna_tax/<s>.rRNA.kraken.report`, `.kraken.min10.report`, `.bracken.<rank>.tsv` | rRNA taxonomy |
| nonhuman | `nonhuman/trna/*.counts.tsv`, `rfam/*.counts.tsv`, `mirna_like/*.counts.tsv` | class counts |
| assembly | `Analysis/<type>/assembly/<s>.bam` (+ `.bai`) | the one mapping of the non-human reads on the contigs (NH tag), used by novel, fragments, bins |
| assembly | `assembly/<s>.{feature,category,contig}_counts.tsv` | reads per feature (sense / antisense, gene, product), category, contig |
| assembly | `assembly/<s>.functions.tsv`, `.length_dist.tsv`, `.mapping_summary.tsv` | reads per functional term of the CDS; read lengths mapped / unmapped; mapping numbers |
| novel | `Analysis/<type>/novel/<s>.sorted.bam` | all alignments on contigs (used by fragments, bins) |
| novel | `novel/<s>.candidates.{intergenic,antisense_to_CDS}.bed`, `.candidates.no_rfam.RNAfold.txt`, `.assembly.summary.tsv` | novel ncRNA candidates |
| fragments | `Analysis/<type>/fragments/<s>.fragments.{features,classes,lengths,ends,summary}.tsv` | fragment classes and profiles |
| fragments | `Analysis/<type>/fragments/<s>.integrity{,.profile,.features}.tsv` | RNA integrity indicators |
| bins | `Analysis/<type>/bins/<s>.{contig,feature,category,bin,taxon}_counts.tsv` | counts by contig, feature, bin, taxon; `bin_counts.tsv` also has GTDB classification, completeness, contamination and quality (HQ ≥ 90 % / < 5 %, MQ ≥ 50 % / < 10 %, else LQ) |
| crispr | `Analysis/<type>/crispr/<s>.crispr_{arrays,units,ends,lengths}.tsv` | reads per CRISPR array (strand, Cas type, bin, taxon), per spacer / repeat in transcription order, read-end offsets, read lengths |
| chimera | `Analysis/<type>/chimera/<s>.chimera_{summary,pairs,typepairs,energy,bins}.tsv`, `.chimera_reads.tsv.gz` | optional module: chimeric reads, feature pairs with enrichment (p, q) and arm complementarity (duplex ΔG vs re-paired arms), bin / taxon of both arms and same-organism share vs expected |
| compare | `Analysis/comparison/<A>_vs_<B>/<s>.<A>_vs_<B>.{features,contigs,bins,bins_cds,bins_no_rrna,summary}.tsv` | one folder per pair (`sRNA-Ext_vs_sRNA-Int`, `sRNA-Ext_vs_metaT`, `sRNA-Int_vs_metaT`, `…_vs_metaG`, `metaT_vs_metaG`): counts, CPM, log2 ratio (first / second); features also CPM and log2 ratio within feature type; summary with Spearman per level and per feature type |
| compare | `Analysis/comparison/<s>.all_sources.{features,contigs,bins}.tsv`, `<s>.compare.summary.tsv` | counts + CPM of all sources side by side; all metrics of all pairs |
| report | `Report/<s>.report.{pdf,html}`, `Report/multiqc/*_mqc.tsv` | per-sample PDF report; MultiQC custom-content tables (see Running → Per-sample report) |

## read_stats.tsv

| step | meaning |
| --- | --- |
| `raw` | raw reads |
| `trimmed` | reads after cutadapt |
| `below_host_min_len_20nt` | trimmed reads shorter than `HOST_MIN_LEN`, never human |
| `host_genome_mapped` | reads aligned to the human genome under the length/mismatch rules |
| `host_tx_mapped` | genome-unmapped reads aligned to the human transcriptome |
| `host_rescued_by_assembly` | host hits that match the sample's contigs at least as well → non-human |
| `human_total` | final human reads = genome + transcriptome − rescued |
| `nonhuman` | trimmed − human_total |

Example, BXXX2 sRNA-Int with the default settings (`HOST_MIN_LEN=20`, exact host matches):

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

## rRNA rescue (module novel)

| file | content |
|---|---|
| `Analysis/<type>/assembly/<s>.rRNA_rescue.bed` | rRNA genes missing from the annotation, found by cmscan on contigs with antisense-only CDS (gff2bed format, source `cmscan`) |
| `Analysis/annotation/<s>.rRNA_rescued.gff` | original GFF3 + rescued rRNA genes of all types; flagged CDS carry `note=opposite_rescued_rRNA` |
| `feature_counts.tsv`, column `flag` | `opposite_rescued_rRNA` for CDS on the other strand of a rescued gene (left out of functions and comparisons) |

