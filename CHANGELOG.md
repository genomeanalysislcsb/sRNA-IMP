# Changelog

## v1.0.0 (2026-11-19)

First release, used in the sRNA-IMP paper (Microbiome, submitted).

**Host and microbial reads**
- Host split with exact matches, a 20-nt minimum, host transcripts, mature tRNAs with their 3′ CCA and short mature
  miRNAs with up to two non-templated 3′ nucleotides; host hits that match the sample's contigs equally well return
  to the microbial fraction (contig rescue).
- Human small RNA classes (miRNAs with 3′ isomiRs, tRNAs, rRNAs, other ncRNAs); microbial RNA classes (RiboDetector,
  Kraken 2 / Bracken rRNA taxonomy, tRNA, miRNA-like, Rfam classes) with equally good hits sharing a read
  (`NONHUMAN_TIES=fractional`) on collapsed tRNA and Rfam references built on first use; origin of miRNA-like reads.

**Assembly-based layers**
- Mapping to an assembly of the same sample from any workflow (IMP3 directly; others via
  `srna_prepare_inputs.py`) or, without metagenome and metatranscriptome, to annotated reference genomes
  (`srna_prepare_inputs.py --genomes`, MGnify catalogue genomes chosen from the library's rRNA species with
  `srna_fetch_mgnify.py`; `ASSEMBLY_SOURCE=reference`).
- Counts per feature and strand, contig, MAG, taxon and functional term; recovery of unannotated rRNA genes;
  novel small RNA loci (Rfam, RNAfold); tRNA / rRNA fragment classes and RNA integrity; crRNAs per CRISPR spacer with
  Cas type and tracrRNA candidates.
- Comparison with metagenomic and metatranscriptomic coverage; strand check of antisense loci with a stranded
  metatranscriptome.

**Running, QC and reporting**
- Per-library modules with completion stamps; SLURM array submission, status and resubmission
  (`submit_srna_pipeline.sh`, `srna_status.sh`); conda environments (`envs/`) and an Apptainer image (`container/`).
- Rule-based library QC (PASS / WARN / FAIL); per-sample PDF / HTML reports; cohort analysis with group tests,
  ordinations and covariate-adjusted models; Methods text generated from the configuration, environments and outputs
  (`srna_methods.py`, Nature and STAR Protocols styles).
- Simulator and evaluation (`srna_simulate.py`, `srna_sim_eval.py`), sMETASeq re-implementation for the benchmark,
  figure scripts of the paper (`paper/`), and two tutorials on the public Tarallo et al. 2019 cohort.
