# =====================================================================
# srna_pipeline.config.sh — configuration for srna_pipeline.sh
# ---------------------------------------------------------------------
# This file is *sourced* by srna_pipeline.sh. Every variable can be
# overridden by passing another config with -c, or by exporting it in
# the environment before the run (env values win over this file only
# where the ${VAR:-default} form is used).
# =====================================================================

# ------------------------- helper scripts ----------------------------
SCRIPTS_DIR=${SCRIPTS_DIR:-/scratch/users/pmay/sRNAs/scripts}

# ------------------------- environments ------------------------------
# micromamba/conda env (name or prefix path). ENV_DEFAULT
# (gut_srna_pipeline) is activated at start of srna_pipeline.sh and in every
# SLURM job; a module whose ENV_<MODULE> differs runs in (and is checked in)
# its own env.
# Leave empty to use whatever is active when the script starts.
# How the env is activated automatically (all optional; tried in this order):
#   ENV_MODULES  modules to 'module load' first, e.g. "tools/micromamba"
#   MAMBA_BIN    path to the micromamba executable
#   CONDA_SH     path to <conda base>/etc/profile.d/conda.sh
# Without these, micromamba/conda are searched in PATH and the usual home
# locations; as last resort an env given as a path is used via its bin/.
ENV_MODULES=${ENV_MODULES:-}
MAMBA_BIN=${MAMBA_BIN:-}
CONDA_SH=${CONDA_SH:-}
ENVS_DIR=${ENVS_DIR:-/work/projects/expobiome/tools/envs}
ENV_DEFAULT=${ENV_DEFAULT:-$ENVS_DIR/gut_srna_pipeline_v2}
ENV_QC=${ENV_QC:-$ENV_DEFAULT}
ENV_SPLIT=${ENV_SPLIT:-$ENV_DEFAULT}
ENV_HUMAN=${ENV_HUMAN:-$ENV_DEFAULT}
ENV_NONHUMAN=${ENV_NONHUMAN:-$ENV_DEFAULT}
# own env: novel module (ARAGORN, cmscan, RNAfold, samtools, bedtools)
ENV_NOVEL=${ENV_NOVEL:-$ENVS_DIR/ncrna_discovery_v2}
# prep (gff2bed from bedops): main env
ENV_PREP=${ENV_PREP:-$ENV_DEFAULT}
# assembly module (bowtie + samtools + python3): main env
ENV_ASSEMBLY=${ENV_ASSEMBLY:-$ENV_DEFAULT}
# new modules: QC in the main env; fragments/bins need samtools + python3
# (taken from the novel env); compare only needs python3
ENV_CONTAM=${ENV_CONTAM:-$ENV_DEFAULT}
ENV_FRAGMENTS=${ENV_FRAGMENTS:-$ENV_NOVEL}
ENV_BINS=${ENV_BINS:-$ENV_NOVEL}
ENV_CRISPR=${ENV_CRISPR:-$ENV_DEFAULT}
ENV_CHIMERA=${ENV_CHIMERA:-$ENV_DEFAULT}   # bowtie + samtools; RNAduplex is taken from ENV_NOVEL if not in this env
ENV_COMPARE=${ENV_COMPARE:-$ENV_DEFAULT}

# ------------------------- temp space --------------------------------
TMP_ROOT=${TMP_ROOT:-${SLURM_TMPDIR:-${TMPDIR:-/tmp}}}

# ------------------------- input staging -----------------------------
# copy   : copy the raw FASTQ into <sample>/Preprocessing/<type>
# symlink: symlink gzipped inputs instead of copying (saves space)
STAGE_MODE=${STAGE_MODE:-copy}


# ------------------------- general -----------------------------------
MINLEN=${MINLEN:-18}          # minimum read length (cutadapt -m, assembly filter)
MAXLEN=${MAXLEN:-100}          # maximum read length (cutadapt -M, length distributions)
MIN_LIBRARY_READS="${MIN_LIBRARY_READS-10000}"   # libraries with fewer trimmed reads are not analysed beyond split (failed library; QC FAIL)

# ------------------------- qc module (cutadapt) ----------------------
ADAPTER_3P=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC   # Illumina TruSeq/NEBNext read-1 adapter, found in the raw reads (same for sRNA-Ext and sRNA-Int, no UMIs)
ADAPTER_3P_TRUNC=10                             # also trim adapter copies lacking up to N leading nt (dimers starting inside the adapter); 0 = off
TRIM_QUALITY=15                                 # cutadapt -q (3' quality trimming); 0 = off
CUTADAPT_EXTRA="--poly-a --trim-n"              # e.g. "--discard-untrimmed --trim-n"

# ------------------------- split module (host = human) ---------------
# Reference FASTAs (plain or .gz). Bowtie indices are looked for at
# HOST_*_IDX (if set) or next to the FASTA as <name>.bowtie.*; if missing
# they are built there when the directory is writable, otherwise in
# HOST_INDEX_TMPDIR (default <outdir>/tmp/host_index, shared by all jobs,
# built only once). Set HOST_TX_FA='' to skip the transcriptome step.
HOST_GENOME_FA=/work/projects/imp/db/imp3/filtering/T2TCHM13v2genomic.fa
HOST_TX_FA=/work/projects/imp/db/imp3/filtering/T2T-CHM13v2.0_transcripts.fa
HOST_GENOME_IDX=""            # optional existing bowtie index prefix
HOST_TX_IDX=""                # optional existing bowtie index prefix
HOST_INDEX_TMPDIR=${HOST_INDEX_TMPDIR:-}
HOST_MISMATCHES=0             # bowtie -v for reads >= HOST_SHORT_LEN (0 = exact; 1 lets ~half of the
                              # 22-35 nt "human" reads in by chance, see docs/qc.md)
# Short reads hit the 3 Gb human genome by chance (an 18-mer with 1
# mismatch matches ~5 times at random), so they must match exactly:
HOST_MIN_LEN=20               # reads shorter than this are never assigned to the host
HOST_SHORT_LEN=22             # reads shorter than this ...
HOST_SHORT_MISMATCHES=0       # ... are mapped to the host with bowtie -v 0
# Host-mapped reads that also map to the sample's own assembly with
# <= mismatches are kept as non-human (needs module prep before split).
HOST_ASSEMBLY_RESCUE=1

# Sub-directory of Analysis/<type>/ for the assembly-mapping results.
# Set to "" to write them directly into Analysis/<type>/ (old layout).
ASSEMBLY_SUBDIR=${ASSEMBLY_SUBDIR-assembly}

# 1 = keep large intermediate files (SAMs, hit tables, ...)
KEEP_INTERMEDIATE=${KEEP_INTERMEDIATE:-0}

# ------------------------- human module ------------------------------
IDX_HSA_MIRNA=/work/projects/imp/db/smallRNA/human/hsa_mature_miRNA
IDX_HSA_TRNA=/work/projects/imp/db/smallRNA/human/hsa_tRNA_all
IDX_HSA_RRNA=/work/projects/imp/db/smallRNA/human/hsa_rRNA_all
IDX_HSA_NCRNA_OTHER=/work/projects/imp/db/smallRNA/human/hsa_otherRNA

# ------------------------- non-human module --------------------------
RIBO_DETECTOR_CPU=${RIBO_DETECTOR_CPU:-$ENVS_DIR/ribodetector/bin/ribodetector_cpu}
RIBO_CHUNK_SIZE=256

IDX_RFAM_CLASSES=/work/projects/imp/db/smallRNA/rfam/rfam_classes
IDX_TRNA_ALL=/work/projects/imp/db/smallRNA/nonhuman/tRNA/trna_all
IDX_MIRNA_EUK=/work/projects/imp/db/smallRNA/nonhuman/miRNAs/mirna_euk_all
RFAMSEQ_TXT=/work/projects/imp/db/smallRNA/rfam/rfamseq.reduced.tsv
RFAM_FAMILY_TXT=/work/projects/imp/db/smallRNA/rfam/family.txt.gz   # Rfam family.txt (same release): adds rfam_type counts (rRNA, sRNA, snoRNA, riboswitch, ...); "" = off

MIRNA_LIKE_MINLEN=18
MIRNA_LIKE_MAXLEN=26

# Kraken2 / Bracken on rRNA reads (settings from runKraken.sh)
#KRAKEN_DB=/work/projects/imp/db/imp3/k2_core_nt_20250609
KRAKEN_DB=${KRAKEN_DB:-/mnt/isilon/projects/ecosystem_biology/databases/imp3/k2_gtdb_genome_reps_20250609}
BRACKEN_DB=${BRACKEN_DB:-$KRAKEN_DB}
KRAKEN_CONF=0.5
KRAKEN_MINHITGROUPS=3
KRAKEN_MINIMIZER_DATA=1       # 1 = --report-minimizer-data (as in runKraken.sh)
KRAKEN_MEMORY_MAPPING=1       # 1 = --memory-mapping
BRACKEN_READ_LENGTH=50        # needs database${BRACKEN_READ_LENGTH}mers.kmer_distrib
KRAKEN_REPORT_MINREADS=10     # drop report lines with < N clade reads (Bracken uses the filtered report)
KRAKEN_MIN_DISTINCT_MINIMIZERS=30   # needs KRAKEN_MINIMIZER_DATA=1: drop species with fewer distinct minimizers (single-locus false positives); 0 = off
KRAKEN_SOFT_FAIL=0            # 0 = a failed kraken2 run fails the nonhuman module (not marked done); 1 = warn and continue without rRNA taxonomy
BRACKEN_LEVELS="D P C O F G S"   # D = domain (Kraken2 rank code); K is mapped to D
BRACKEN_MIN_READS=10

# ------------------------- novel ncRNA module ------------------------
RFAM_CM=/work/projects/imp/db/smallRNA/rfam/Rfam.cm
NOVEL_EXCLUDE_RRNA="${NOVEL_EXCLUDE_RRNA-1}"   # 1 = leave the RiboDetector rRNA reads out of the novel loci (fewer rRNA candidates, faster cmscan)
RRNA_EXCLUDE_RFAM="${RRNA_EXCLUDE_RFAM-1}"   # 1 = also leave out reads with an Rfam rRNA hit (by sequence) from the *_non_rRNA counts and novel loci:
                                            # RiboDetector misses many short rRNA fragments (test library sRNA-Int: 2.4 M of 5.9 M rRNA reads)
RRNA_RESCUE="${RRNA_RESCUE-1}"               # module novel: cmscan (Rfam rRNA models) on contigs whose CDS have only antisense reads;
                                            # hits are added as rRNA genes (<s>.rRNA_rescue.bed) and the assembly counts are redone
RRNA_RESCUE_MIN_ANTISENSE="${RRNA_RESCUE_MIN_ANTISENSE-20}"     # CDS with >= N antisense reads ...
RRNA_RESCUE_MAX_SENSE_RATIO="${RRNA_RESCUE_MAX_SENSE_RATIO-0.1}" # ... and sense <= ratio x antisense are checked
RRNA_RESCUE_CM="${RRNA_RESCUE_CM-}"         # rRNA subset of RFAM_CM (built once, cmpress); "" = <outdir>/.shared/Rfam.rRNA.cm
NOVEL_MINLEN=18
NOVEL_MAXLEN=250
NOVEL_MINREADS=10            # min. read support per locus; multi-mapping reads count 1/hits per alignment
NOVEL_MERGEDIST=10
RUN_ARAGORN=1
RUN_TRNASCAN=0
RUN_CMSCAN=1
CMSCAN_OPTS="--rfam --cut_ga --toponly"   # candidates are extracted in read (RNA) orientation: search that strand only; truncated hits allowed
RUN_RNAFOLD=1

# ------------------------- assembly module ---------------------------
# one bowtie mapping of the uncollapsed non-human reads to the contigs, used by assembly, novel,
# fragments, bins and compare (BAM: Analysis/<type>/assembly/<sample>.bam, NH tag = hits per read)
ASSEMBLY_MISMATCHES=1         # bowtie -v
ASSEMBLY_MAX_HITS=50          # bowtie -k N --best --strata (0 = all hits, -a)
ASSEMBLY_RECOUNT_ONLY="${ASSEMBLY_RECOUNT_ONLY-0}"  # 1 = module assembly reuses an existing BAM and only recounts (e.g. -m assembly -f after a counting update)

# ------------------------- strandedness -------------------------------
# sRNA libraries are stranded (read = RNA strand). With DB_NORC=1 the reads are matched only to the
# forward strand of the RNA databases (human miRNA/tRNA/rRNA/other ncRNA, microbial tRNA, miRNA,
# Rfam), i.e. reverse-complement (antisense) hits are not counted as that RNA. "" = both strands.
DB_NORC=1

# ------------------------- contam module (library QC) ----------------
# Space-separated list of contaminant FASTAs (plain or .gz), each screened
# separately; bowtie indices are found/built like the host indices.
# e.g. CONTAM_FA="/path/UniVec_Core.fa /path/phiX174.fa /path/kit_contaminants.fa"
CONTAM_FA="/work/projects/expobiome/sRNAs/PhiX.fa /work/projects/expobiome/sRNAs/UniVec_Core.fa" #${CONTAM_FA:-}
CONTAM_MISMATCHES=1
QC_TOP_SEQS=50                # overrepresented sequences reported
QC_SORT_MEM=2G                # sort buffer for the complexity count

# ------------------------- fragments module --------------------------
FRAG_TYPES="tRNA,rRNA,tmRNA"  # gff feature types analysed
FRAG_TOL=3                    # nt tolerance for "starts/ends at the 5'/3' end"
FRAG_AC_FROM=0.40             # anticodon region (fraction of tRNA length)
FRAG_AC_TO=0.55
FRAG_MIN_OVERLAP=0.9          # fraction of the read inside the feature
FRAG_MIN_READS=10             # features with >= N reads get per-position profiles

# ------------------------- bins module -------------------------------
BINS_MIN_OVERLAP=0.5          # fraction of the read overlapping a feature

# ------------------------- chimera module (optional: -m chimera) -------
# chimeric reads = candidate RNA-RNA ligation products (exploratory without crosslinking / proximity ligation)
CHIMERA_MIN_ARM=${CHIMERA_MIN_ARM:-18}            # minimum arm (= exact seed) length; reads >= 2x this are tested
CHIMERA_MAX_HITS=${CHIMERA_MAX_HITS:-5}           # seeds with more exact hits on the assembly are not used
CHIMERA_ARM_MISMATCHES=${CHIMERA_ARM_MISMATCHES:-1}  # mismatches allowed while extending an arm
CHIMERA_MAX_GAP=${CHIMERA_MAX_GAP:-2}             # read nt covered by neither arm
CHIMERA_COLINEAR_MAX=${CHIMERA_COLINEAR_MAX:-1000}  # same contig/strand, 3' arm <= N nt downstream: colinear (not used)
CHIMERA_UTR5=${CHIMERA_UTR5:-150}                 # 5' UTR window upstream of a CDS (sRNA targets)
CHIMERA_UTR3=${CHIMERA_UTR3:-100}                 # 3' UTR window downstream of a CDS
CHIMERA_FLANK=${CHIMERA_FLANK:-20}                # contig nt on both sides of an arm for the duplex energy
CHIMERA_NULL=${CHIMERA_NULL:-20}                  # random re-pairings per pair (complementarity null)
CHIMERA_MAX_PAIRS=${CHIMERA_MAX_PAIRS:-5000}      # feature pairs tested for complementarity (most reads first)
RNADUPLEX=${RNADUPLEX:-}                          # RNAduplex executable; "" = PATH, then $ENV_NOVEL/bin/RNAduplex
# metaT chimeras (same module, once per sample): split reads and discordant pairs in the IMP3 metaT BAM
CHIMERA_MT=${CHIMERA_MT:-1}                       # 1 = also analyse the metaT BAM
CHIMERA_MT_BAM=${CHIMERA_MT_BAM:-}                # "" = <IMP dir of metaT_features>/Assembly/mt.reads.sorted.bam; <ID>/<id> allowed
CHIMERA_MT_STRAND=${CHIMERA_MT_STRAND:-auto}      # auto | R1-sense | R1-antisense (dUTP) | unstranded
CHIMERA_MT_MIN_ARM=${CHIMERA_MT_MIN_ARM:-25}      # minimum aligned length of each part of a split read
CHIMERA_MT_MIN_MAPQ=${CHIMERA_MT_MIN_MAPQ:-20}    # MAPQ of primary alignments and discordant mates
CHIMERA_MT_EDGE=${CHIMERA_MT_EDGE:-500}           # discordant mate pointing to a contig end closer than this: assembly break

# ------------------------- compare (--compare) -----------------------
COMPARE_PSEUDOCOUNT=1         # added to CPM before log2 ratios
COMPARE_MIN_COUNT=5           # min raw count in at least one of the two sources of a pair
COMPARE_BIN_BASIS=contig,cds,no_rrna   # bin count bases; first = bins.tsv, others = bins_<basis>.tsv (cds: metaT/metaG CDS reads only; no_rrna: rRNA contigs left out)

# ------------------------- report (--report) ------------------------
REPORT_HTML=1                 # also write <sample>.report.html (same pages, self-contained); "" = PDF only
REPORT_MIRNA_MARKERS=""          # TSV of human miRNA origin markers + references; "" = scripts/mirna_origin_markers.tsv
RNA_QC=""                     # wet-lab RNA QC TSV (sample, sample_id, type, lRNA_RIN, ... ; type = sRNA-Int/sRNA-Ext), e.g. rna_qc.tsv; "" = none
QC_RULES="${QC_RULES-}"        # QC rules TSV (thresholds -> PASS / WARN / FAIL per library); "" = scripts/qc_rules.tsv
QC_EXCLUDE="${QC_EXCLUDE-}"    # manual QC decisions: sample, type (or *), level (FAIL/WARN/PASS), reason; "" = none

# ------------------------- cohort (--cohort) ------------------------
COHORT_METADATA="${COHORT_METADATA-}"            # metadata TSV/CSV/XLSX, one row per sample (xlsx needs openpyxl); "" = no groups
COHORT_ID_COL="${COHORT_ID_COL-}"              # sample ID column; "" = first column (IDs like BXXX1 or BXXX1.t0F both match)
COHORT_GROUP_COL="${COHORT_GROUP_COL-group}"      # group column (e.g. HC / PD / RA, or HC / T1D)
COHORT_GROUP_ORDER="${COHORT_GROUP_ORDER-}"         # e.g. "HC,PD,RA" or "HC,T1D" (order in tables and plots)
COHORT_COVARIATES="${COHORT_COVARIATES-}"          # comma-separated covariate columns; "" = all other metadata columns
COHORT_PERMUTATIONS="${COHORT_PERMUTATIONS-999}"       # PERMANOVA permutations
COHORT_DIR="${COHORT_DIR-}"                 # "" = <outdir>/Cohort
COHORT_BLOCK="${COHORT_BLOCK-${COHORT_MATCH_COL-}}"  # block column: site (Expobiome: nested, few levels) or family (T1D: crossed, many levels);
                                            # adds pooled + within-block comparisons, the poolability check of the reference group,
                                            # and (crossed) block-clustered / within-block model variants. "" = off
COHORT_MATCH_COL="$COHORT_BLOCK"            # old name, kept for compatibility
COHORT_REFERENCE="${COHORT_REFERENCE-}"     # reference group for the matched comparisons; "" = first group of COHORT_GROUP_ORDER
COHORT_OUTLIER_Z="${COHORT_OUTLIER_Z-4}"     # libraries with |robust z| > this in key metrics of their fraction get WARN "outlier"; 0 = off
COHORT_MIN_MIRNA_READS="${COHORT_MIN_MIRNA_READS-50}"  # libraries with fewer human miRNA reads are left out of the miRNA tests
COHORT_SIG_Q="${COHORT_SIG_Q-0.1}"           # FDR q below which results are bold / listed as key findings
COHORT_TERM_NAMES="${COHORT_TERM_NAMES-}"   # optional TSVs term<TAB>name for readable KEGG / Pfam labels, comma-separated
                               # (e.g. KEGG REST list/ko, list/pathway, list/module; Pfam-A.clans.tsv)
COHORT_NMDS="${COHORT_NMDS-genus,fn_kegg_pathway,mirna,rfam_family}"  # NMDS (Bray-Curtis) pages; also family, rna_class, fn_kegg_ko, fn_pfam; "" = none
COHORT_NMDS_STARTS="${COHORT_NMDS_STARTS-10}"   # random starts per NMDS (+ 1 PCoA start); ~30 s per NMDS for 180 libraries
COHORT_NAME="${COHORT_NAME-cohort}"   # name of the cohort: title of the cohort report
COHORT_MODEL="${COHORT_MODEL-}"   # covariate-adjusted models (OLS on CLR, HC3); several models comma-separated, terms joined by +,
                                  # e.g. "group+site, group+site+sex+age, group+site+sex+age+constipation_status" (= M1, M2, M3);
                                  # terms: metadata columns or library metrics (e.g. qc_lRNA_RIN); "" = none
COHORT_MODEL_NOTES="${COHORT_MODEL_NOTES-}"   # optional question per model for the models page, e.g. "M9=PD motor severity; M13=RA immunosuppressants"
COHORT_QC_DROP="${COHORT_QC_DROP-FAIL}"   # libraries left out of the cohort statistics: FAIL | WARN (= WARN+FAIL) | none
MULTIQC="${MULTIQC-multiqc}"   # MultiQC command (or full path, e.g. $ENVS_DIR/multiqc/bin/multiqc) run by --cohort over all
METHODS_STYLE="${METHODS_STYLE-nature,star}"   # automated Methods with --report / --cohort: nature, star (STAR Protocols) or both; none = off
METHODS_REPORT_STYLE="${METHODS_REPORT_STYLE-star}"   # Methods added as last pages of the sample and cohort reports: star, nature or none
METHODS_RESOURCES="${METHODS_RESOURCES-}"      # databases / versions for the Methods; "" = scripts/methods_resources.tsv (edit once)
METHODS_DOCX="${METHODS_DOCX-}"                # 1 = also .docx (needs pandoc)
                               # FastQC / cutadapt outputs and the report tables; "" = skip
