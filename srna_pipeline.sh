#!/usr/bin/env bash
# =====================================================================
# srna_pipeline.sh — small RNA analysis pipeline (expobiome)
# ---------------------------------------------------------------------
# Per sample-sheet row (= one sample + one sRNA type), default order:
#   qc        FastQC raw -> cutadapt (adapter/quality/length) -> FastQC trimmed
#   contam    library QC: adapter dimers, complexity, overrepresented
#             sequences, contaminant screen (UniVec/PhiX/...)
#   prep      assembly + annotation staging (shared by all types of a sample)
#   split     host (human) genome -> host transcriptome mapping (bowtie),
#             human / non-human FASTQs; reads < HOST_MIN_LEN nt are never
#             human, reads < HOST_SHORT_LEN nt must match
#             exactly, host hits that also match the sample assembly stay
#             non-human (HOST_ASSEMBLY_RESCUE)
#   human     human ncRNA classification (miRNA, tRNA, rRNA, other)
#   nonhuman  rRNA filter, Kraken2/Bracken, tRNA / miRNA-like / Rfam classes
#   assembly  mapping on the sample's contigs, CMC weighting, feature counts
#   novel     novel ncRNA loci on the contigs
#   fragments tRNA / rRNA fragment classes, cleavage sites, coverage profiles
#   bins      counts per contig, feature, bin (MAG) and GTDB taxon
#   crispr    sRNAs on CRISPR arrays (crRNAs): per array / spacer, strand, read ends, Cas type
# Per sample (all types together, run with --compare):
#   compare   sRNA-Ext vs sRNA-Int vs metaT vs metaG on features, contigs and bins
#             (any combination of >= 2 available sources)
#
# Directory layout (one tree per sample, one sub-tree per sRNA type):
#
#   <outdir>/<sample>/
#     Preprocessing/<type>/   <sample>.raw.fastq.gz, .trimmed.fastq.gz,
#                             .human.fastq.gz, .human.rmDup.fastq.gz, .nonhuman.fastq.gz
#     Stats/<type>/           fastqc/, cutadapt log, read_stats.tsv, seqkit statistics
#     Assembly/               <sample>.mgmt.assembly.merged.fa.gz (+ bwa/bowtie idx)
#     Analysis/annotation/    <sample>.gff, <sample>.bed      (shared by all types)
#     Analysis/<type>/
#        human/               mirna trna rrna otherrna unknown
#        nonhuman/            rrna rrna_tax collapsed trna mirna_like rfam unknown
#        assembly/            mapping on contigs + feature/db count matrices
#        novel/               novel ncRNA discovery
#        fragments/           tRNA/rRNA fragment analysis
#        bins/                per contig / feature / bin / taxon counts
#     Analysis/comparison/    Ext / Int / metaT / metaG comparison (--compare)
#     log_cluster_jobs/<type>/<module>/
#
# Sample sheet (tab- or space-separated, '#' comments, optional header):
#   sample_id type raw_fastq assembly gff [contig2bin gtdbtk metaT_features metaT_contigs metaG_features metaG_contigs]
# raw_fastq, assembly and gff are required (full paths). The optional
# columns 6-11 (full paths or "-") are per sample and may be given on any
# row of that sample:
#   contig2bin      contig<TAB>bin table
#   gtdbtk          GTDB-Tk summary tsv (comma-separated for bac120,ar53)
#   metaT_features  featureCounts table of metaT reads per gff feature
#   metaT_contigs   metaT reads per contig: IMP Stats/mt/mgmt.assembly.contig_read_count.txt.gz,
#                   the IMP Stats dir, a tar of it, or a featureCounts table per contig
#   metaG_features  IMP mg.ID.* featureCounts tables (dir/tar/list); default:
#                   mg tables found in metaT_features (IMP annotation dir or tar)
#   metaG_contigs   same as metaT_contigs for metaG; default: mg files in metaT_contigs
# The same sample can appear with several types; Assembly/ and
# Analysis/annotation/ are shared.
# =====================================================================

set -euo pipefail

SCRIPT_PATH=$(readlink -f "${BASH_SOURCE[0]}")
PIPELINE_DIR=$(dirname "$SCRIPT_PATH")

ALL_MODULES="qc,contam,prep,split,human,nonhuman,assembly,novel,fragments,bins,crispr"

usage() {
cat <<EOF
Usage: $(basename "$0") -s <samples.tsv> [options]

Required:
  -s, --samples FILE     sample sheet (see header of this script)

Options:
  -o, --outdir DIR       base output directory               [\$PWD/imp]
  -c, --config FILE      config file                          [$PIPELINE_DIR/srna_pipeline.config.sh]
  -e, --env ENV          conda/micromamba env (name or path)  [ENV_DEFAULT from config]
  -m, --modules LIST     comma-separated modules              [$ALL_MODULES]
  -t, --threads N        threads                              [\$SLURM_CPUS_PER_TASK or 16]
  -r, --row N            process only data row N (1-based; used by SLURM arrays)
  -S, --sample ID        process only rows with this sample id
  -T, --type TYPE        process only rows with this type
  -f, --force            re-run modules that already finished (.done stamps)
      --compare          run the per-sample comparison (Ext / Int / metaT / metaG) for all
                         samples (or --sample) instead of the per-row modules
      --report           write the per-sample PDF report <outdir>/<sample>/Report/<sample>.report.pdf
                         (+ MultiQC tables) for all samples (or --sample); with --compare: after it
      --cohort           cross-sample summary of all samples of the sheet with the metadata
                         (COHORT_METADATA): tables, group statistics, <outdir>/Cohort/cohort.report.pdf;
                         runs after --compare / --report when combined
      --check            only validate environment, tools, databases and inputs
      --check-env        only check environment, tools and databases, then exit
      --skip-env-check   do not run the environment/tool check before starting
      --count-rows       print the number of data rows and exit
  -h, --help             this help
EOF
}

# ------------------------------------------------------------------
# argument parsing
# ------------------------------------------------------------------
SAMPLESHEET=""
OUTDIR="$PWD/imp"
CONFIG="$PIPELINE_DIR/srna_pipeline.config.sh"
MODULES="$ALL_MODULES"
THREADS="${SLURM_CPUS_PER_TASK:-16}"
ROW=""
SAMPLE_FILTER=""
TYPE_FILTER=""
FORCE=0
CHECK_ONLY=0
COUNT_ONLY=0
ENV_CHECK_ONLY=0
SKIP_ENV_CHECK=0
COMPARE=0
REPORT=0
COHORT=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    -s|--samples)  SAMPLESHEET=$2; shift 2 ;;
    -o|--outdir)   OUTDIR=$2; shift 2 ;;
    -c|--config)   CONFIG=$2; shift 2 ;;
    -e|--env)      export ENV_DEFAULT=$2; shift 2 ;;
    -m|--modules)  MODULES=$2; shift 2 ;;
    -t|--threads)  THREADS=$2; shift 2 ;;
    -r|--row)      ROW=$2; shift 2 ;;
    -S|--sample)   SAMPLE_FILTER=$2; shift 2 ;;
    -T|--type)     TYPE_FILTER=$2; shift 2 ;;
    -f|--force)    FORCE=1; shift ;;
    --check)       CHECK_ONLY=1; shift ;;
    --check-env)   ENV_CHECK_ONLY=1; shift ;;
    --skip-env-check) SKIP_ENV_CHECK=1; shift ;;
    --compare)     COMPARE=1; shift ;;
    --report)      REPORT=1; shift ;;
    --cohort)      COHORT=1; shift ;;
    --count-rows)  COUNT_ONLY=1; shift ;;
    -h|--help)     usage; exit 0 ;;
    *) echo "Error: unknown option '$1'" >&2; usage >&2; exit 1 ;;
  esac
done

[[ -n "$SAMPLESHEET" ]] || { echo "Error: -s/--samples is required" >&2; usage >&2; exit 1; }
[[ -s "$SAMPLESHEET" ]] || { echo "Error: sample sheet '$SAMPLESHEET' not found or empty" >&2; exit 1; }
[[ -s "$CONFIG" ]]      || { echo "Error: config '$CONFIG' not found" >&2; exit 1; }

SAMPLESHEET=$(readlink -f "$SAMPLESHEET")
CONFIG=$(readlink -f "$CONFIG")

# shellcheck source=srna_pipeline.config.sh
source "$CONFIG"

for m in ${MODULES//,/ }; do
  # validated against ALL_MODULES (single source of the module list)
  [[ ",$ALL_MODULES," == *",$m,"* ]] || { echo "Error: unknown module '$m' (allowed: $ALL_MODULES)" >&2; exit 1; }
done
# run the selected modules in pipeline order, whatever order -m lists them in
_ord=""
for _m in ${ALL_MODULES//,/ }; do [[ ",$MODULES," == *",$_m,"* ]] && _ord+="${_ord:+,}$_m"; done
[[ ",$MODULES," == *",compare,"* ]] && _ord+="${_ord:+,}compare"
MODULES=$_ord
# --compare / --report: only these per-sample steps (and their environment check)
if [[ $COMPARE -eq 1 || $REPORT -eq 1 || $COHORT -eq 1 ]]; then
  MODULES=""
  [[ $COMPARE -eq 1 ]] && MODULES=compare
  [[ $REPORT -eq 1 ]] && MODULES+="${MODULES:+,}report"
  [[ $COHORT -eq 1 ]] && MODULES+="${MODULES:+,}cohort"
fi

# ------------------------------------------------------------------
# generic helpers
# ------------------------------------------------------------------
log()  { echo "[$(date '+%F %T')][${PREFIX:-}][${TYPE:-}][${MODULE:-main}] $*"; }
warn() { log "WARNING: $*" >&2; }
die()  { log "ERROR: $*" >&2; exit 1; }

# read all data rows of the sample sheet into ROWS[]
ROWS=()
while IFS= read -r line || [[ -n "$line" ]]; do
  line=${line%$'\r'}
  [[ -z "${line//[[:space:]]/}" ]] && continue
  [[ "$line" =~ ^[[:space:]]*# ]] && continue
  first=$(echo "$line" | awk '{print tolower($1)}')
  [[ "$first" == "sample_id" || "$first" == "sample" || "$first" == "sampleid" ]] && continue
  ROWS+=("$line")
done < "$SAMPLESHEET"

if [[ $COUNT_ONLY -eq 1 ]]; then echo "${#ROWS[@]}"; exit 0; fi
[[ ${#ROWS[@]} -gt 0 ]] || { echo "Error: no data rows in $SAMPLESHEET" >&2; exit 1; }

is_gz() { [[ "$1" == *.gz ]]; }

# stage a FASTQ/FASTA so that the destination is always gzipped
stage_gz() {        # src dst(.gz)
  local src=$1 dst=$2
  if [[ -s "$dst" && $FORCE -eq 0 ]]; then log "  exists, keeping: $dst"; return 0; fi
  [[ -s "$src" ]] || die "input not found or empty: $src"
  rm -f "$dst"
  if is_gz "$src"; then
    if [[ "$STAGE_MODE" == "symlink" ]]; then ln -s "$(readlink -f "$src")" "$dst"
    else cp -p "$src" "$dst.part" && mv "$dst.part" "$dst"; fi
  else
    gzip -c "$src" > "$dst.part" && mv "$dst.part" "$dst"
  fi
  log "  staged: $src -> $dst"
}

# stage a text file uncompressed
stage_plain() {     # src dst
  local src=$1 dst=$2
  [[ -s "$src" ]] || die "input not found or empty: $src"
  if is_gz "$src"; then zcat "$src" > "$dst.part"; else cat "$src" > "$dst.part"; fi
  mv "$dst.part" "$dst"
}

# --- portable lock (mkdir is atomic on NFS/GPFS/Lustre) -------------
HELD_LOCKS=()
acquire_lock() {    # lockdir [max_wait_seconds]
  local l=$1 max=${2:-86400} waited=0
  until mkdir "$l" 2>/dev/null; do
    if [[ $waited -eq 0 ]]; then log "  waiting for lock $l (held by: $(cat "$l/owner" 2>/dev/null || echo '?'))"; fi
    sleep 20; waited=$((waited + 20))
    [[ $waited -ge $max ]] && die "timeout waiting for lock $l — remove it if no other job is running"
  done
  echo "${SLURM_JOB_ID:-local}:$(hostname):$$" > "$l/owner"
  HELD_LOCKS+=("$l")
}
release_lock() { rm -rf "$1"; }

cleanup_module() {
  local l
  for l in "${HELD_LOCKS[@]+"${HELD_LOCKS[@]}"}"; do rm -rf "$l"; done
  [[ -n "${TEMPDIR:-}" && -d "${TEMPDIR:-}" ]] && rm -rf "$TEMPDIR"
  return 0
}

make_tempdir() {
  mkdir -p "$TMP_ROOT"
  TEMPDIR=$(mktemp -d -p "$TMP_ROOT" "sRNA.${PREFIX}.${TYPE}.${MODULE}.XXXXXX") \
    || die "could not create temporary directory in $TMP_ROOT"
}

# activation failure: fatal, unless ACT_SOFT=1 (used by the env check)
_act_fail() { if [[ "${ACT_SOFT:-0}" -eq 1 ]]; then warn "$*"; else die "$*"; fi; }

activate_env() {    # env name or prefix
  local env=$1
  [[ -z "$env" ]] && return 0
  # already active?
  if [[ -n "${CONDA_PREFIX:-}" ]] && \
     [[ "$CONDA_PREFIX" == "${env%/}" || "$(basename "$CONDA_PREFIX")" == "$env" ]]; then
    return 0
  fi
  set +u
  # 1) optional environment modules (e.g. ENV_MODULES="tools/micromamba")
  if [[ -n "${ENV_MODULES:-}" ]]; then
    if ! type module >/dev/null 2>&1; then
      for f in /etc/profile.d/lmod.sh /etc/profile.d/modules.sh /usr/share/lmod/lmod/init/bash; do
        [[ -r "$f" ]] && { source "$f"; break; }
      done
    fi
    type module >/dev/null 2>&1 && module load $ENV_MODULES \
      || warn "could not 'module load $ENV_MODULES'"
  fi
  # 2) locate micromamba / conda even in non-interactive (batch) shells
  local mm="" cs="" c
  for c in "${MAMBA_BIN:-}" "$(command -v micromamba 2>/dev/null)" "${MAMBA_EXE:-}" \
           "$HOME/.local/bin/micromamba" "$HOME/bin/micromamba" "$HOME/micromamba/bin/micromamba"; do
    [[ -n "$c" && -x "$c" ]] && { mm=$c; break; }
  done
  if [[ -z "$mm" ]]; then
    for c in "${CONDA_SH:-}" "${CONDA_EXE:+$(dirname "$(dirname "$CONDA_EXE")")/etc/profile.d/conda.sh}" \
             "$HOME/miniforge3/etc/profile.d/conda.sh" "$HOME/mambaforge/etc/profile.d/conda.sh" \
             "$HOME/miniconda3/etc/profile.d/conda.sh" "$HOME/anaconda3/etc/profile.d/conda.sh"; do
      [[ -n "$c" && -r "$c" ]] && { cs=$c; break; }
    done
  fi
  # 3) activate
  if [[ "$(type -t micromamba 2>/dev/null)" == function ]]; then
    # shell hook already loaded (e.g. by a previous activation)
    micromamba activate "$env" || { set -u; _act_fail "cannot activate env '$env' with micromamba"; return 1; }
  elif [[ -n "$mm" ]]; then
    eval "$("$mm" shell hook --shell bash)"
    micromamba activate "$env" || { set -u; _act_fail "cannot activate env '$env' with $mm"; return 1; }
  elif [[ "$(type -t conda 2>/dev/null)" == function ]]; then
    conda activate "$env" || { set -u; _act_fail "cannot activate conda env '$env'"; return 1; }
  elif [[ -n "$cs" ]] || command -v conda >/dev/null 2>&1; then
    if [[ -n "$cs" ]]; then source "$cs"; else eval "$(conda shell.bash hook)"; fi
    conda activate "$env" || { set -u; _act_fail "cannot activate conda env '$env'"; return 1; }
  elif [[ -d "$env/bin" ]]; then
    # 4) fallback: env given as a path — put its bin/ first in PATH
    warn "neither micromamba nor conda found — using $env/bin via PATH"
    export PATH="$env/bin:$PATH" CONDA_PREFIX="${env%/}"
  else
    set -u
    _act_fail "cannot activate '$env': no micromamba/conda found (set MAMBA_BIN, CONDA_SH or ENV_MODULES in the config)"; return 1
  fi
  set -u
  log "activated environment: ${CONDA_PREFIX:-?}"
}

need() { local e; for e in "$@"; do command -v "$e" >/dev/null 2>&1 || die "required executable '$e' not found in PATH"; done; }

count_counted_fasta_reads() {   # sum of '-<count>' suffixes in collapsed fasta
  local FA=$1
  [[ -s "$FA" ]] || { echo 0; return; }
  awk '/^>/{n=split(substr($0,2),a,"-"); s+=a[n]} END{print s+0}' "$FA"
}
gz_cmd() { if command -v pigz >/dev/null 2>&1; then echo "pigz -p $THREADS"; else echo "gzip"; fi; }

# --- bowtie index for the host reference ---------------------------
bowtie_index_ok() {   # prefix
  local p=$1 ext
  for ext in ebwt ebwtl; do
    [[ -s "$p.1.$ext" && -s "$p.2.$ext" && -s "$p.3.$ext" && -s "$p.4.$ext" \
       && -s "$p.rev.1.$ext" && -s "$p.rev.2.$ext" ]] && return 0
  done
  return 1
}
dir_writable() {      # dir — real write test (ACLs / read-only mounts)
  local d=$1 t
  [[ -d "$d" && -w "$d" ]] || return 1
  t="$d/.srna_write_test.$$.$RANDOM"
  ( : > "$t" ) 2>/dev/null || return 1
  rm -f "$t" 2>/dev/null || true
  return 0
}
# ensure_host_index <fasta> <configured prefix or ""> <label>  -> sets HOST_IDX
#   1. use <prefix> (or <fasta dir>/<fasta name>.bowtie) if the index exists
#   2. else use an index already built in HOST_INDEX_TMPDIR
#   3. else build it next to the fasta if that directory is writable,
#      otherwise in HOST_INDEX_TMPDIR (lock-protected, built only once)
ensure_host_index() {
  local fa=$1 prefix=$2 label=$3 base fb target
  [[ -s "$fa" ]] || die "$label fasta not found: $fa"
  if [[ -z "$prefix" ]]; then
    base=$(basename "$fa"); base=${base%.gz}; base=${base%.fa}; base=${base%.fasta}; base=${base%.fna}
    prefix=$(dirname "$(readlink -f "$fa")")/$base.bowtie
  fi
  fb=$HOST_INDEX_TMPDIR/$(basename "$prefix")
  if   bowtie_index_ok "$prefix"; then HOST_IDX=$prefix; log "  $label index: $HOST_IDX"; return 0
  elif bowtie_index_ok "$fb";     then HOST_IDX=$fb;     log "  $label index: $HOST_IDX"; return 0
  fi
  if dir_writable "$(dirname "$prefix")"; then
    target=$prefix
  else
    mkdir -p "$HOST_INDEX_TMPDIR"
    target=$fb
    log "  $(dirname "$prefix") not writable — building $label index in $HOST_INDEX_TMPDIR"
  fi
  acquire_lock "$target.build.lock" 172800
  if ! bowtie_index_ok "$target"; then
    log "  building $label bowtie index $target (this takes a while) ..."
    bowtie-build --threads "$THREADS" "$fa" "$target" > "$target.build.log" 2>&1 \
      || { release_lock "$target.build.lock"; die "bowtie-build failed for $fa (see $target.build.log)"; }
  fi
  release_lock "$target.build.lock"
  bowtie_index_ok "$target" || die "$label index incomplete after build: $target"
  HOST_IDX=$target
  log "  $label index: $HOST_IDX"
}

# bowtie index of the sample assembly (shared by split + novel, built once)
ensure_assembly_bowtie_index() {
  local LOCK=$ASSEMBLY_DIR/.bowtie.lock
  acquire_lock "$LOCK"
  if bowtie_index_ok "$ASSEMBLY"; then
    log "  assembly bowtie index present"
  else
    log "  building bowtie index of the assembly ..."
    bowtie-build --threads "$THREADS" "$ASSEMBLY.gz" "$ASSEMBLY" > "$ASSEMBLY_DIR/$PREFIX.bowtie-build.log" 2>&1 \
      || { release_lock "$LOCK"; die "bowtie-build of the assembly failed (see $ASSEMBLY_DIR/$PREFIX.bowtie-build.log)"; }
  fi
  release_lock "$LOCK"
}

fastq_reads() { [[ -s "$1" ]] && zcat -f "$1" | awk 'END{print int(NR/4)}' || echo 0; }

# ------------------------------------------------------------------
# per-row path setup
# ------------------------------------------------------------------
set_paths() {
  SAMPLE_DIR=$OUTDIR/$PREFIX
  PREPROC_DIR=$SAMPLE_DIR/Preprocessing/$TYPE
  STATS_DIR=$SAMPLE_DIR/Stats/$TYPE
  ASSEMBLY_DIR=$SAMPLE_DIR/Assembly
  ANNO_DIR=$SAMPLE_DIR/Analysis/annotation
  ANALYSIS_DIR=$SAMPLE_DIR/Analysis/$TYPE
  LOGDIR=$SAMPLE_DIR/log_cluster_jobs/$TYPE

  HUMANDIR=$ANALYSIS_DIR/human
  NONHUMANDIR=$ANALYSIS_DIR/nonhuman
  NOVELDIR=$ANALYSIS_DIR/novel
  FRAGDIR=$ANALYSIS_DIR/fragments
  BINSDIR=$ANALYSIS_DIR/bins
  CRISPRDIR=$ANALYSIS_DIR/crispr
  COMPDIR=$SAMPLE_DIR/Analysis/comparison
  if [[ -n "$ASSEMBLY_SUBDIR" ]]; then ASMDIR=$ANALYSIS_DIR/$ASSEMBLY_SUBDIR; else ASMDIR=$ANALYSIS_DIR; fi

  RAWREADS=$PREPROC_DIR/$PREFIX.raw.fastq.gz
  TRIMMEDREADS=$PREPROC_DIR/$PREFIX.trimmed.fastq.gz
  HUMANREADS_RAW=$PREPROC_DIR/$PREFIX.human.fastq.gz
  HUMANREADS=$PREPROC_DIR/$PREFIX.human.rmDup          # .fastq.gz
  NONHUMANREADS=$PREPROC_DIR/$PREFIX.nonhuman          # .fastq.gz

  ASSEMBLY=$ASSEMBLY_DIR/$PREFIX.mgmt.assembly.merged.fa   # stored as $ASSEMBLY.gz
  GFF=$ANNO_DIR/$PREFIX.gff
  GBED=$ANNO_DIR/$PREFIX.bed
}

make_dirs() {
  mkdir -p "$PREPROC_DIR" "$STATS_DIR" "$ASSEMBLY_DIR" "$ANNO_DIR" "$ANALYSIS_DIR" "$LOGDIR"
}

# ==================================================================
# MODULE: qc — FastQC on raw reads, cutadapt trimming, FastQC trimmed
# ==================================================================
module_qc() {
  need cutadapt fastqc
  local L=$LOGDIR/qc QCDIR=$STATS_DIR/fastqc
  mkdir -p "$L" "$QCDIR"

  log "staging raw reads ..."
  stage_gz "$IN_RAW" "$RAWREADS"

  log "FastQC on raw reads ..."
  fastqc -t "$THREADS" --outdir "$QCDIR" "$RAWREADS" > "$L/$PREFIX.fastqc.raw.log" 2>&1

  log "cutadapt: adapter $ADAPTER_3P, quality >= $TRIM_QUALITY, length $MINLEN-$MAXLEN nt ..."
  local part=$PREPROC_DIR/$PREFIX.trimmed.part.fastq.gz
  # shellcheck disable=SC2086
  # adapter copies missing their first <= ADAPTER_3P_TRUNC nt (adapter/primer dimers that start
  # inside the adapter, e.g. GGAAGAGCACACG...) are not found by a regular 3' adapter: add the
  # adapter suffix as a second adapter; such reads are trimmed to <= ADAPTER_3P_TRUNC nt and dropped
  # (-m). Normal reads still match the full adapter (more matching bases wins).
  local adapt=(-a "$ADAPTER_3P")
  [[ "${ADAPTER_3P_TRUNC:-0}" -gt 0 && ${#ADAPTER_3P} -gt $((ADAPTER_3P_TRUNC + 15)) ]] \
    && adapt+=(-a "${ADAPTER_3P:$ADAPTER_3P_TRUNC}")
  cutadapt "${adapt[@]}" -q "$TRIM_QUALITY" -m "$MINLEN" -M "$MAXLEN" -j "$THREADS" \
    $CUTADAPT_EXTRA -o "$part" "$RAWREADS" \
    > "$STATS_DIR/$PREFIX.$TYPE.cutadapt.log" 2> "$L/$PREFIX.cutadapt.err"
  mv "$part" "$TRIMMEDREADS"

  log "FastQC on trimmed reads ..."
  fastqc -t "$THREADS" --outdir "$QCDIR" "$TRIMMEDREADS" > "$L/$PREFIX.fastqc.trimmed.log" 2>&1
}

# ==================================================================
# MODULE: split — host (human) vs non-host reads
#   trimmed -> bowtie host genome -> unmapped -> bowtie host transcriptome
#   human    = genome-mapped + transcriptome-mapped
#   nonhuman = unmapped to both
# ==================================================================
module_split() {
  need bowtie bowtie-build samtools seqkit
  local L=$LOGDIR/split
  mkdir -p "$L"
  [[ -s "$TRIMMEDREADS" ]] || die "missing $TRIMMEDREADS — run module qc first"
  local SHORT_LEN=${HOST_SHORT_LEN:-22} SHORT_MM=${HOST_SHORT_MISMATCHES:-0}
  local HMIN=${HOST_MIN_LEN:-20}
  (( SHORT_LEN < HMIN )) && SHORT_LEN=$HMIN
  local N_TOOSHORT=0
  local RESCUE=${HOST_ASSEMBLY_RESCUE:-1}
  if [[ "$RESCUE" -eq 1 ]]; then
    [[ -s "$ASSEMBLY.gz" ]] || die "missing $ASSEMBLY.gz — run module prep before split (needed for HOST_ASSEMBLY_RESCUE=1)"
  fi
  make_tempdir
  local Z; Z=$(gz_cmd)

  log "checking host indices ..."
  ensure_host_index "$HOST_GENOME_FA" "${HOST_GENOME_IDX:-}" "host genome"
  local GIDX=$HOST_IDX TIDX=""
  if [[ -n "${HOST_TX_FA:-}" ]]; then
    ensure_host_index "$HOST_TX_FA" "${HOST_TX_IDX:-}" "host transcriptome"
    TIDX=$HOST_IDX
  fi

  local G_UN=$TEMPDIR/genome_unmapped.fastq HG=$TEMPDIR/host_genome.fastq
  local HT=$TEMPDIR/host_tx.fastq NONHOST=$TEMPDIR/nonhost.fastq
  local NM_HOST=$TEMPDIR/host.nm.tsv
  : > "$NM_HOST"

  # host_map <index> <in.fastq[.gz]> <mapped.fastq> <unmapped.fastq> <label>
  #   reads <  HMIN nt     : never host (too short to be assigned reliably)
  #   reads <  SHORT_LEN nt: bowtie -v SHORT_MM   (short reads hit the human
  #   reads >= SHORT_LEN nt: bowtie -v HOST_MISMATCHES   genome by chance)
  #   appends "<read>\t<mismatches>" of mapped reads to NM_HOST
  host_map() {
    local idx=$1 in=$2 out_m=$3 out_u=$4 lab=$5 t v f
    : > "$out_m"; : > "$out_u"
    : > "$TEMPDIR/$lab.tooshort.fq"; : > "$TEMPDIR/$lab.short.fq"
    if (( HMIN > 1 )); then
      seqkit seq -j "$THREADS" -M $((HMIN - 1)) "$in" > "$TEMPDIR/$lab.tooshort.fq" 2>/dev/null
    fi
    if (( SHORT_LEN > HMIN )); then
      seqkit seq -j "$THREADS" -m "$HMIN" -M $((SHORT_LEN - 1)) "$in" > "$TEMPDIR/$lab.short.fq" 2>/dev/null
    fi
    seqkit seq -j "$THREADS" -m "$SHORT_LEN" "$in" > "$TEMPDIR/$lab.long.fq" 2>/dev/null
    [[ "$lab" == host_genome ]] && N_TOOSHORT=$(fastq_reads "$TEMPDIR/$lab.tooshort.fq")
    cat "$TEMPDIR/$lab.tooshort.fq" >> "$out_u"; rm -f "$TEMPDIR/$lab.tooshort.fq"
    for t in short long; do
      f=$TEMPDIR/$lab.$t.fq
      [[ "$t" == short ]] && v=$SHORT_MM || v=$HOST_MISMATCHES
      [[ -s "$f" ]] || continue
      bowtie -p "$THREADS" -v "$v" -k 1 --best --sam --un "$f.un" \
        "$idx" "$f" > "$f.sam" 2> "$L/$PREFIX.bowtie.$lab.$t.log"
      samtools fastq -F 4 "$f.sam" >> "$out_m" 2>> "$L/$PREFIX.samtools_fastq.$lab.log"
      awk -F '\t' '!/^@/ && int($2/4)%2==0 { nm=0; if (match($0, /NM:i:[0-9]+/)) nm=substr($0, RSTART+5, RLENGTH-5); print $1"\t"nm }' \
        "$f.sam" >> "$NM_HOST"
      [[ -e "$f.un" ]] && cat "$f.un" >> "$out_u"
      rm -f "$f" "$f.sam" "$f.un"
    done
  }

  log "mapping to host genome (< $HMIN nt: never host, $HMIN-$((SHORT_LEN-1)) nt: -v $SHORT_MM, >= $SHORT_LEN nt: -v $HOST_MISMATCHES; -k 1 --best) ..."
  host_map "$GIDX" "$TRIMMEDREADS" "$HG" "$G_UN" host_genome

  if [[ -n "$TIDX" ]]; then
    log "mapping genome-unmapped reads to host transcriptome ..."
    host_map "$TIDX" "$G_UN" "$HT" "$NONHOST" host_tx
  else
    log "no HOST_TX_FA configured — transcriptome step skipped"
    mv "$G_UN" "$NONHOST"; : > "$HT"
  fi

  local HUMAN_ALL=$TEMPDIR/human_all.fastq RESCUED=$TEMPDIR/rescued.fastq
  cat "$HG" "$HT" > "$HUMAN_ALL"; : > "$RESCUED"
  local N_RESCUED=0

  # ---- assembly rescue: host-mapped reads that map to the sample's own
  #      contigs at least as well (<= mismatches) are non-human ----
  if [[ "$RESCUE" -eq 1 && -s "$HUMAN_ALL" ]]; then
    ensure_assembly_bowtie_index
    log "checking host-mapped reads against the sample assembly ..."
    bowtie -p "$THREADS" -v "$HOST_MISMATCHES" -k 1 --best --sam \
      "$ASSEMBLY" "$HUMAN_ALL" > "$TEMPDIR/asm.sam" 2> "$L/$PREFIX.bowtie.assembly_rescue.log"
    awk -F '\t' '!/^@/ && int($2/4)%2==0 { nm=0; if (match($0, /NM:i:[0-9]+/)) nm=substr($0, RSTART+5, RLENGTH-5); print $1"\t"nm }' \
      "$TEMPDIR/asm.sam" > "$TEMPDIR/asm.nm.tsv"
    rm -f "$TEMPDIR/asm.sam"
    # rescued = mismatches(assembly) <= mismatches(host)
    awk -F '\t' 'NR==FNR { h[$1]=$2; next } ($1 in h) && $2 <= h[$1] { print $1 }' \
      "$NM_HOST" "$TEMPDIR/asm.nm.tsv" | sort -u > "$TEMPDIR/rescued.ids"
    if [[ -s "$TEMPDIR/rescued.ids" ]]; then
      seqkit grep -j "$THREADS" -f "$TEMPDIR/rescued.ids" "$HUMAN_ALL" > "$RESCUED"
      seqkit grep -j "$THREADS" -v -f "$TEMPDIR/rescued.ids" "$HUMAN_ALL" > "$HUMAN_ALL.keep"
      mv "$HUMAN_ALL.keep" "$HUMAN_ALL"
      cat "$RESCUED" >> "$NONHOST"
    fi
    N_RESCUED=$(fastq_reads "$RESCUED")
    log "  $N_RESCUED host-mapped reads also match the assembly -> non-human"
    cp "$TEMPDIR/rescued.ids" "$L/$PREFIX.assembly_rescued.ids" 2>/dev/null || true
  elif [[ "$RESCUE" -ne 1 ]]; then
    log "assembly rescue disabled (HOST_ASSEMBLY_RESCUE=0)"
  fi

  log "writing human / non-human FASTQs ..."
  $Z -c "$HUMAN_ALL" > "$HUMANREADS_RAW.part" && mv "$HUMANREADS_RAW.part" "$HUMANREADS_RAW"
  $Z -c "$NONHOST" > "$NONHUMANREADS.fastq.gz.part" && mv "$NONHUMANREADS.fastq.gz.part" "$NONHUMANREADS.fastq.gz"

  log "removing duplicate read names from human reads ..."
  seqkit rmdup --by-name -j "$THREADS" "$HUMANREADS_RAW" -o "$HUMANREADS.fastq.gz" \
    > "$L/$PREFIX.seqkit.rmdup.log" 2>&1

  log "read statistics ..."
  seqkit stats -j "$THREADS" "$RAWREADS" "$TRIMMEDREADS" "$HUMANREADS_RAW" "$HUMANREADS.fastq.gz" \
    "$NONHUMANREADS.fastq.gz" > "$STATS_DIR/$PREFIX.$TYPE.seqkit.txt"
  {
    echo -e "sample\tstep\treads"
    echo -e "$PREFIX\traw\t$(fastq_reads "$RAWREADS")"
    echo -e "$PREFIX\ttrimmed\t$(fastq_reads "$TRIMMEDREADS")"
    echo -e "$PREFIX\tbelow_host_min_len_${HMIN}nt\t$N_TOOSHORT"
    echo -e "$PREFIX\thost_genome_mapped\t$(fastq_reads "$HG")"
    echo -e "$PREFIX\thost_tx_mapped\t$(fastq_reads "$HT")"
    echo -e "$PREFIX\thost_rescued_by_assembly\t$N_RESCUED"
    echo -e "$PREFIX\thuman_total\t$(fastq_reads "$HUMANREADS_RAW")"
    echo -e "$PREFIX\tnonhuman\t$(fastq_reads "$NONHUMANREADS.fastq.gz")"
  } > "$STATS_DIR/$PREFIX.$TYPE.read_stats.tsv"
  {
    echo -e "length\thuman_reads"
    zcat "$HUMANREADS_RAW" | awk 'NR%4==2{c[length($0)]++} END{for(l in c) print l"\t"c[l]}' | sort -n
  } > "$STATS_DIR/$PREFIX.$TYPE.human_length_dist.tsv"
  log "read stats: $STATS_DIR/$PREFIX.$TYPE.read_stats.tsv"
}

# ==================================================================
# MODULE: prep — assembly + annotation (shared between types)
# ==================================================================
module_prep() {
  need gff2bed
  local L=$LOGDIR/prep; mkdir -p "$L"

  # ---- shared per-sample resources (several types may run at once) ----
  local LOCK=$ASSEMBLY_DIR/.prep.lock
  acquire_lock "$LOCK"
  if [[ ! -s "$ASSEMBLY.gz" || $FORCE -eq 1 ]]; then
    log "staging assembly ..."
    stage_gz "$IN_ASSEMBLY" "$ASSEMBLY.gz"
  else
    log "assembly already present: $ASSEMBLY.gz"
  fi
  if [[ ! -s "$GFF" || ! -s "$GBED" || $FORCE -eq 1 ]]; then
    log "staging annotation (gff -> bed) ..."
    stage_plain "$IN_GFF" "$GFF.tmp"
    # change "?" strand annotations to "."
    sed 's/\t?\t/\t.\t/g' "$GFF.tmp" > "$GFF" && rm -f "$GFF.tmp"
    gff2bed < "$GFF" > "$GBED.part" && mv "$GBED.part" "$GBED"
  else
    log "annotation already present: $GFF"
  fi
  release_lock "$LOCK"
}

# ==================================================================
# MODULE: human — sequential bowtie screen miRNA > tRNA > rRNA > other
# ==================================================================
module_human() {
  need bowtie python3
  local L=$LOGDIR/human
  mkdir -p "$L" "$HUMANDIR"/{mirna,trna,rrna,otherrna,unknown}
  [[ -s "$HUMANREADS.fastq.gz" ]] || die "missing $HUMANREADS.fastq.gz — run module split first"
  make_tempdir

  # step <name> <index> <input.fastq.gz> <unmapped basename> <outdir>
  human_step() {
    local name=$1 idx=$2 in=$3 unname=$4 outsub=$5
    log "$name searching ..."
    bowtie -p "$THREADS" ${DB_NORC:+--norc} -v 1 -a --best --sam --un "$TEMPDIR/$unname.fastq" \
      "$idx" "$in" > "$TEMPDIR/$PREFIX.human.$name.sam" 2> "$L/$PREFIX.bowtie.$name.log"
    gzip -f "$TEMPDIR/$unname.fastq"
    mv "$TEMPDIR/$unname.fastq.gz" "$HUMANDIR/$outsub/$unname.fastq.gz"
    log "$name searching done"
  }

  human_step mirna    "$IDX_HSA_MIRNA"       "$HUMANREADS.fastq.gz" \
             "$PREFIX.human.non_mirna" mirna
  human_step trna     "$IDX_HSA_TRNA"        "$HUMANDIR/mirna/$PREFIX.human.non_mirna.fastq.gz" \
             "$PREFIX.human.non_mirna.non_trna" trna
  human_step rrna     "$IDX_HSA_RRNA"        "$HUMANDIR/trna/$PREFIX.human.non_mirna.non_trna.fastq.gz" \
             "$PREFIX.human.non_mirna.non_trna.non_rrna" rrna
  human_step otherrna "$IDX_HSA_NCRNA_OTHER" "$HUMANDIR/rrna/$PREFIX.human.non_mirna.non_trna.non_rrna.fastq.gz" \
             "$PREFIX.human.unknown" otherrna

  log "feature counting ..."
  local c
  for c in mirna trna rrna otherrna; do
    python3 "$SCRIPTS_DIR/sam_fractional_counts.py" \
      -i "$TEMPDIR/$PREFIX.human.$c.sam" -o "$TEMPDIR/$PREFIX.human.$c.counts.tsv"
    mv "$TEMPDIR/$PREFIX.human.$c.counts.tsv" "$HUMANDIR/$c/"
  done

  log "summarizing human non-coding RNA counts ..."
  python3 "$SCRIPTS_DIR/summarize_human_subtree.py" "$HUMANDIR"
  log "human non-coding RNA search done"
}

# ==================================================================
# MODULE: nonhuman — RiboDetector, Kraken2/Bracken (runKraken.sh),
#                    tRNA, miRNA-like, Rfam screens, summary
# ==================================================================
module_nonhuman() {
  need seqkit bowtie python3 kraken2 bracken
  [[ -x "$RIBO_DETECTOR_CPU" ]] || need "$RIBO_DETECTOR_CPU"
  local L=$LOGDIR/nonhuman
  mkdir -p "$L" "$NONHUMANDIR"/{rrna,rrna_tax,collapsed,trna,mirna_like,rfam,unknown}
  [[ -s "$NONHUMANREADS.fastq.gz" ]] || die "missing $NONHUMANREADS.fastq.gz — run module split first"
  make_tempdir

  local COLLAPSE_PY="$SCRIPTS_DIR/collapse_fastq_to_counted_fasta.py"
  local SUM_BOWTIE_PY="$SCRIPTS_DIR/sum_collapsed_bowtie_hits.py"
  local SUM_TRNA_TAX_PY="$SCRIPTS_DIR/sum_collapsed_trna_by_taxon.py"
  local SUM_RFAM_TAX_PY="$SCRIPTS_DIR/sum_collapsed_rfam_by_class_and_taxon.py"

  # ---------------- 1) rRNA filtering with RiboDetector ----------------
  log "calculating average read length ..."
  local R_MEAN_L
  R_MEAN_L=$(seqkit stats -T "$NONHUMANREADS.fastq.gz" | cut -f 7 | tail -n +2 | awk '{printf("%d\n", $1 + 0.5)}')
  log "ribodetector rRNA filtering (mean read length $R_MEAN_L) ..."

  local RRNA_FASTQ=$NONHUMANDIR/rrna/$PREFIX.rRNA.fastq
  local NONRRNA=$NONHUMANDIR/rrna/$PREFIX.nonrRNA.fastq
  "$RIBO_DETECTOR_CPU" \
    --threads "$THREADS" \
    --len "$R_MEAN_L" \
    --input "$NONHUMANREADS.fastq.gz" \
    --ensure none \
    --chunk_size "$RIBO_CHUNK_SIZE" \
    --output "$NONRRNA" \
    --rrna "$RRNA_FASTQ" > "$L/$PREFIX.ribodetector.log" 2>&1
  gzip -f "$NONRRNA" "$RRNA_FASTQ"
  local RRNA_FASTQ_GZ=$RRNA_FASTQ.gz NONRRNA_FASTQ_GZ=$NONRRNA.gz
  log "ribodetector rRNA filtering done"

  # ---------------- 2) rRNA taxonomy: Kraken2 + Bracken ----------------
  # (latest version, from runKraken.sh)
  log "Kraken2 taxonomy on rRNA reads ..."
  local RRNA_TAX_DIR=$NONHUMANDIR/rrna_tax
  local KRAKEN_REPORT=$RRNA_TAX_DIR/$PREFIX.rRNA.kraken.report
  local KRAKEN_OUT=$RRNA_TAX_DIR/$PREFIX.rRNA.kraken.out
  local KRAKEN_LOG=$L/$PREFIX.kraken.log
  local BRACKEN_LOG=$L/$PREFIX.bracken.log
  local kopts=()
  [[ "$KRAKEN_MEMORY_MAPPING" -eq 1 ]] && kopts+=(--memory-mapping)
  [[ "$KRAKEN_MINIMIZER_DATA" -eq 1 ]] && kopts+=(--report-minimizer-data)

  # results of an earlier run must not survive a failed re-run
  rm -f "$KRAKEN_REPORT" "$KRAKEN_OUT" "$KRAKEN_OUT.gz" "$RRNA_TAX_DIR/$PREFIX.rRNA.kraken.min"*.report \
        "$RRNA_TAX_DIR/$PREFIX.rRNA.bracken."*
  : > "$BRACKEN_LOG"
  if [[ -s "$RRNA_FASTQ_GZ" ]]; then
    if kraken2 --db "$KRAKEN_DB" \
        --confidence "$KRAKEN_CONF" \
        --threads "$THREADS" \
        --use-names \
        --minimum-hit-groups "$KRAKEN_MINHITGROUPS" \
        --gzip-compressed \
        ${kopts[@]+"${kopts[@]}"} \
        --report "$KRAKEN_REPORT" \
        --output "$KRAKEN_OUT" \
        "$RRNA_FASTQ_GZ" > "$KRAKEN_LOG" 2>&1; then
      # filter the report: keep taxa with >= KRAKEN_REPORT_MINREADS clade reads
      # (column 2; same position with/without --report-minimizer-data).
      # The unclassified line (taxid 0) is always kept. Clade counts never
      # increase towards the leaves, so the filtered tree stays consistent.
      log "filtering Kraken2 report (clade reads >= $KRAKEN_REPORT_MINREADS) ..."
      local KRAKEN_REPORT_FILT=$RRNA_TAX_DIR/$PREFIX.rRNA.kraken.min${KRAKEN_REPORT_MINREADS}.report
      # with --report-minimizer-data, species (rank S*) whose reads cover fewer than
      # KRAKEN_MIN_DISTINCT_MINIMIZERS distinct minimizers are dropped as well: such hits come from
      # one short conserved stretch (seen: 21,635 reads on 19 minimizers), not from the taxon's rRNA
      awk -F '\t' -v min="$KRAKEN_REPORT_MINREADS" -v mind="${KRAKEN_MIN_DISTINCT_MINIMIZERS:-0}" '
        { taxid = (NF >= 8) ? $7 : $5 }
        taxid == 0 { print; next }
        $2 < min { next }
        NF >= 8 && $6 ~ /^S/ && $5 < mind { dropped++; next }
        { print }
        END { if (dropped) printf "%d species below %d distinct minimizers dropped\n", dropped, mind > "/dev/stderr" }' \
        "$KRAKEN_REPORT" > "$KRAKEN_REPORT_FILT" 2>> "$KRAKEN_LOG"
      log "  kept $(wc -l < "$KRAKEN_REPORT_FILT") of $(wc -l < "$KRAKEN_REPORT") report lines"

      : > "$BRACKEN_LOG"
      local lvl
      for lvl in $BRACKEN_LEVELS; do
        # Kraken2 reports the (super)kingdom of Bacteria/Archaea as rank D (domain); Bracken's K
        # (kingdom) is empty for prokaryotes and only produces "Error: no reads found"
        [[ "$lvl" == K ]] && { log "  Bracken level K -> D (domain)"; lvl=D; }
        bracken -d "$BRACKEN_DB" -i "$KRAKEN_REPORT_FILT" \
          -o "$RRNA_TAX_DIR/$PREFIX.rRNA.bracken.$lvl.tsv" \
          -w "$RRNA_TAX_DIR/$PREFIX.rRNA.bracken.$lvl.report" \
          -r "$BRACKEN_READ_LENGTH" -l "$lvl" -t "$BRACKEN_MIN_READS" \
          >> "$BRACKEN_LOG" 2>&1 || warn "bracken level $lvl failed (see $BRACKEN_LOG)"
        [[ -s "$RRNA_TAX_DIR/$PREFIX.rRNA.bracken.$lvl.tsv" ]] \
          || warn "bracken level $lvl: no output (no taxa with >= $BRACKEN_MIN_READS reads at this rank?)"
      done
      [[ "$KEEP_INTERMEDIATE" -eq 1 ]] || gzip -f "$KRAKEN_OUT"
    else
      # e.g. killed (memory / time) right after "Loading database information"
      if [[ "${KRAKEN_SOFT_FAIL:-0}" -eq 1 ]]; then
        warn "kraken2 failed (see $KRAKEN_LOG) — continuing without rRNA taxonomy (KRAKEN_SOFT_FAIL=1)"
      else
        die "kraken2 failed (see $KRAKEN_LOG; often memory: the DB must fit in RAM or page cache) — module not marked done"
      fi
    fi
  else
    warn "no rRNA reads — Kraken2/Bracken skipped"
  fi
  log "Kraken2 taxonomy on rRNA done"

  # ---------------- 3) collapse non-rRNA FASTQ -> counted FASTA --------
  log "collapsing non-rRNA reads ..."
  local COLLAPSED_FA=$NONHUMANDIR/collapsed/$PREFIX.non_rRNA.collapsed.fa
  python3 "$COLLAPSE_PY" -i "$NONRRNA_FASTQ_GZ" -o "$COLLAPSED_FA"

  # ---------------- 4) tRNA screen ------------------------------------
  log "tRNA screen ..."
  local TRNA_DIR=$NONHUMANDIR/trna
  local TRNA_HITS=$TRNA_DIR/$PREFIX.trna.hits.txt
  local TRNA_UN_FA=$TRNA_DIR/$PREFIX.non_tRNA.fa
  local TRNA_COUNTS=$TRNA_DIR/$PREFIX.trna.counts.tsv
  local TRNA_PHYLUM_COUNTS=$TRNA_DIR/$PREFIX.trna_phylum_class.counts.tsv
  local TRNA_SPECIES_COUNTS=$TRNA_DIR/$PREFIX.trna_species.counts.tsv
  local TRNA_AA_COUNTS=$TRNA_DIR/$PREFIX.trna_aa.counts.tsv
  local TRNA_ANTICODON_COUNTS=$TRNA_DIR/$PREFIX.trna_anticodon.counts.tsv
  local TRNA_PHYLUM_SPECIES_COUNTS=$TRNA_DIR/$PREFIX.trna_phylum_species.counts.tsv

  # --fullref: the tRNA / Rfam / miRNA headers contain spaces (species names, Rfam descriptions);
  # without it bowtie cuts the reference name at the first space and aa/anticodon/class/family are lost
  bowtie -p "$THREADS" ${DB_NORC:+--norc} -v 1 -k 1 --best --fullref -f --un "$TRNA_UN_FA" \
    "$IDX_TRNA_ALL" "$COLLAPSED_FA" > "$TRNA_HITS" 2> "$L/$PREFIX.bowtie.trna.log"
  [[ -e "$TRNA_UN_FA" ]] || : > "$TRNA_UN_FA"

  if [[ -s "$TRNA_HITS" ]]; then
    python3 "$SUM_TRNA_TAX_PY" -i "$TRNA_HITS" \
      --feature-out "$TRNA_COUNTS" \
      --phylum-out "$TRNA_PHYLUM_COUNTS" \
      --species-out "$TRNA_SPECIES_COUNTS" \
      --aa-out "$TRNA_AA_COUNTS" \
      --anticodon-out "$TRNA_ANTICODON_COUNTS" \
      --phylum-species-out "$TRNA_PHYLUM_SPECIES_COUNTS"
  else
    printf "feature\tcount\n"              > "$TRNA_COUNTS"
    printf "phylum_class\tcount\n"         > "$TRNA_PHYLUM_COUNTS"
    printf "species\tcount\n"              > "$TRNA_SPECIES_COUNTS"
    printf "aa\tcount\n"                   > "$TRNA_AA_COUNTS"
    printf "anticodon\tcount\n"            > "$TRNA_ANTICODON_COUNTS"
    printf "phylum_class\tspecies\tcount\n" > "$TRNA_PHYLUM_SPECIES_COUNTS"
  fi

  # ---------------- 5) split tRNA-unmatched by length -----------------
  log "splitting tRNA-unmatched reads by length ..."
  local MIRNA_DIR=$NONHUMANDIR/mirna_like
  local MIRNA_IN_FA=$MIRNA_DIR/$PREFIX.non_tRNA.L${MIRNA_LIKE_MINLEN}_${MIRNA_LIKE_MAXLEN}.fa
  local LONGER_IN_FA=$TEMPDIR/$PREFIX.non_tRNA.GT${MIRNA_LIKE_MAXLEN}.fa
  seqkit seq -g -m "$MIRNA_LIKE_MINLEN" -M "$MIRNA_LIKE_MAXLEN" "$TRNA_UN_FA" > "$MIRNA_IN_FA"
  seqkit seq -g -m $((MIRNA_LIKE_MAXLEN + 1)) "$TRNA_UN_FA" > "$LONGER_IN_FA"

  # ---------------- 6) miRNA-like screen ------------------------------
  log "miRNA-like screen ..."
  local MIRNA_HITS=$MIRNA_DIR/$PREFIX.mirna_like.hits.txt
  local MIRNA_UN_FA=$MIRNA_DIR/$PREFIX.non_tRNA_non_miRNA.L${MIRNA_LIKE_MINLEN}_${MIRNA_LIKE_MAXLEN}.fa
  local MIRNA_COUNTS=$MIRNA_DIR/$PREFIX.mirna_like.counts.tsv
  if [[ -s "$MIRNA_IN_FA" ]]; then
    bowtie -p "$THREADS" ${DB_NORC:+--norc} -v 0 -k 1 --best --fullref -f --un "$MIRNA_UN_FA" \
      "$IDX_MIRNA_EUK" "$MIRNA_IN_FA" > "$MIRNA_HITS" 2> "$L/$PREFIX.bowtie.mirna_like.log"
    [[ -e "$MIRNA_UN_FA" ]] || : > "$MIRNA_UN_FA"
    python3 "$SUM_BOWTIE_PY" -i "$MIRNA_HITS" -o "$MIRNA_COUNTS"
  else
    : > "$MIRNA_HITS"; : > "$MIRNA_UN_FA"
    printf "feature\tcount\n" > "$MIRNA_COUNTS"
  fi

  # ---------------- 7) Rfam class screen on remainder -----------------
  log "Rfam class screen ..."
  local RFAM_DIR=$NONHUMANDIR/rfam
  local RFAM_INPUT_FA=$RFAM_DIR/$PREFIX.rfam_input.fa
  local RFAM_HITS=$RFAM_DIR/$PREFIX.rfam.hits.txt
  local RFAM_UN_FA=$NONHUMANDIR/unknown/$PREFIX.unknown.fa
  local RFAM_CLASS_COUNTS=$RFAM_DIR/$PREFIX.rfam_class.counts.tsv
  local RFAM_TAXON_COUNTS=$RFAM_DIR/$PREFIX.rfam_taxon.counts.tsv
  local RFAM_CLASS_TAXON_COUNTS=$RFAM_DIR/$PREFIX.rfam_class_taxon.counts.tsv
  local RFAM_FAMILY_COUNTS=$RFAM_DIR/$PREFIX.rfam_family.counts.tsv
  local RFAM_ACCESSION_COUNTS=$RFAM_DIR/$PREFIX.rfam_accession.counts.tsv

  cat "$MIRNA_UN_FA" "$LONGER_IN_FA" > "$RFAM_INPUT_FA"
  if [[ -s "$RFAM_INPUT_FA" ]]; then
    bowtie -p "$THREADS" ${DB_NORC:+--norc} -v 1 -k 1 --best --fullref -f --un "$RFAM_UN_FA" \
      "$IDX_RFAM_CLASSES" "$RFAM_INPUT_FA" > "$RFAM_HITS" 2> "$L/$PREFIX.bowtie.rfam.log"
    [[ -e "$RFAM_UN_FA" ]] || : > "$RFAM_UN_FA"
    python3 "$SUM_RFAM_TAX_PY" -i "$RFAM_HITS" -s "$RFAMSEQ_TXT" \
      --class-out "$RFAM_CLASS_COUNTS" \
      --taxon-out "$RFAM_TAXON_COUNTS" \
      --class-taxon-out "$RFAM_CLASS_TAXON_COUNTS" \
      --family-out "$RFAM_FAMILY_COUNTS" \
      --accession-out "$RFAM_ACCESSION_COUNTS" \
      ${RFAM_FAMILY_TXT:+--family-txt "$RFAM_FAMILY_TXT" --type-out "$RFAM_DIR/$PREFIX.rfam_type.counts.tsv" \
        --type-taxon-out "$RFAM_DIR/$PREFIX.rfam_type_taxon.counts.tsv"}
    # reads with an Rfam rRNA hit that RiboDetector missed (mostly short LSU/SSU fragments): their sequences,
    # used with the RiboDetector reads to leave rRNA out of the antisense counts and the novel loci
    python3 "$SCRIPTS_DIR/srna_rrna_seqs.py" -i "$RFAM_HITS" -o "$RFAM_DIR/$PREFIX.rfam_rRNA.seqs.gz" \
      2>> "$L/$PREFIX.bowtie.rfam.log" || warn "srna_rrna_seqs.py failed — Rfam rRNA reads will not be excluded"
  else
    : > "$RFAM_HITS"; : > "$RFAM_UN_FA"
    printf "rfam_class\tcount\n"        > "$RFAM_CLASS_COUNTS"
    printf "taxon\tcount\n"             > "$RFAM_TAXON_COUNTS"
    printf "rfam_class\ttaxon\tcount\n" > "$RFAM_CLASS_TAXON_COUNTS"
    printf "rfam_family\tcount\n"       > "$RFAM_FAMILY_COUNTS"
    printf "accession\tcount\n"         > "$RFAM_ACCESSION_COUNTS"
    [[ -n "${RFAM_FAMILY_TXT:-}" ]] && printf "rfam_type\tcount\n" > "$RFAM_DIR/$PREFIX.rfam_type.counts.tsv"
  fi

  # ---------------- 8) summary -----------------------------------------
  log "writing summary ..."
  local SUMMARY=$NONHUMANDIR/$PREFIX.nonhuman.summary.tsv
  {
    echo -e "sample\tcategory\tcount"
    echo -e "$PREFIX\tinput_reads\t$(fastq_reads "$NONHUMANREADS.fastq.gz")"
    echo -e "$PREFIX\trRNA_reads\t$(fastq_reads "$RRNA_FASTQ_GZ")"
    echo -e "$PREFIX\tnon_rRNA_reads\t$(fastq_reads "$NONRRNA_FASTQ_GZ")"
    echo -e "$PREFIX\tnon_rRNA_collapsed_reads\t$(count_counted_fasta_reads "$COLLAPSED_FA")"
    echo -e "$PREFIX\ttRNA\t$(awk -F'\t' 'NR>1{s+=$NF} END{print s+0}' "$TRNA_COUNTS")"
    echo -e "$PREFIX\tmiRNA_like\t$(awk -F'\t' 'NR>1{s+=$NF} END{print s+0}' "$MIRNA_COUNTS")"
    echo -e "$PREFIX\tRfam_classified\t$(awk -F'\t' 'NR>1{s+=$NF} END{print s+0}' "$RFAM_CLASS_COUNTS")"
    echo -e "$PREFIX\tunknown\t$(count_counted_fasta_reads "$RFAM_UN_FA")"
  } > "$SUMMARY"

  # ---------------- 9) clean up ----------------------------------------
  rm -f "$MIRNA_UN_FA" "$RFAM_INPUT_FA"
  log "summary: $SUMMARY"
}

# ==================================================================
# MODULE: assembly — map non-human reads on contigs, weight multi-
#                    mappers (CMC), intersect with annotation, counts
# ==================================================================
# sequences of reads with an Rfam rRNA hit (module nonhuman); built from the hits file of older runs if missing.
# Prints the path (empty if not available or RRNA_EXCLUDE_RFAM=0).
rfam_rrna_seqs() {
  [[ "${RRNA_EXCLUDE_RFAM:-1}" -eq 1 ]] || return 0
  local F=$NONHUMANDIR/rfam/$PREFIX.rfam_rRNA.seqs.gz H=$NONHUMANDIR/rfam/$PREFIX.rfam.hits.txt
  if [[ ! -s "$F" && -s "$H" ]]; then
    python3 "$SCRIPTS_DIR/srna_rrna_seqs.py" -i "$H" -o "$F" >&2 || { rm -f "$F"; return 0; }
  fi
  [[ -s "$F" ]] && echo "$F"
  return 0
}

# annotation_bed: $GBED, or $GBED + the rRNA genes rescued by cmscan for this type (module novel) in $TEMPDIR
annotation_bed() {
  local RESC=$ASMDIR/$PREFIX.rRNA_rescue.bed gb=$TEMPDIR/annotation.rescued.bed
  if [[ -s "$RESC" ]]; then
    [[ -s "$gb" ]] || cat "$GBED" "$RESC" | sort -k1,1 -k2,2n -T "$TEMPDIR" > "$gb"
    echo "$gb"
  else
    echo "$GBED"
  fi
}

# assembly_counts LOG: feature, category, contig and functional counts of the assembly BAM ($ASMDIR/$PREFIX.*).
# Annotation = $GBED + rRNA genes rescued by cmscan (module novel, <s>.rRNA_rescue.bed) if present.
# Reads classified as rRNA by RiboDetector, and reads with an Rfam rRNA hit, are left out of the *_non_rRNA columns:
# reads sense to unannotated (partial) rRNA would otherwise count as antisense to spurious ORFs on the other strand.
assembly_counts() {
  local log=$1 P=$ASMDIR/$PREFIX BAM=$ASMDIR/$PREFIX.bam gb RESC=$ASMDIR/$PREFIX.rRNA_rescue.bed
  rm -f "$TEMPDIR/annotation.rescued.bed"
  gb=$(annotation_bed)
  [[ -s "$RESC" ]] && log "  + $(wc -l < "$RESC") rRNA genes rescued by cmscan ($RESC)"
  local excl=() RR=$NONHUMANDIR/rrna/$PREFIX.rRNA.fastq.gz
  [[ -s "$RR" ]] && excl=(--exclude-reads "$RR") || warn "no RiboDetector rRNA reads ($RR) — *_non_rRNA columns = all reads"
  local RS; RS=$(rfam_rrna_seqs)
  [[ -n "$RS" ]] && excl+=(--exclude-seqs "$RS") || warn "no Rfam rRNA sequence list — only RiboDetector rRNA reads are excluded"
  python3 "$SCRIPTS_DIR/srna_bins.py" --bam "$BAM" --gbed "$gb" --out-prefix "$P" \
    --functions --lengths --no-bins --threads "$THREADS" ${excl[@]+"${excl[@]}"} > "$log" 2>&1 \
    || die "srna_bins.py failed (see $log)"
}

module_assembly() {
  need bowtie bowtie-build samtools python3
  local L=$LOGDIR/assembly
  mkdir -p "$L" "$ASMDIR"
  local READS=$NONHUMANREADS.fastq.gz
  [[ -s "$READS" ]] || die "missing $READS — run module split first"
  [[ -s "$ASSEMBLY.gz" && -s "$GBED" ]] || die "missing assembly/annotation — run module prep first"
  make_tempdir
  local P=$ASMDIR/$PREFIX BAM=$ASMDIR/$PREFIX.bam

  if [[ "${ASSEMBLY_RECOUNT_ONLY:-0}" == 1 && -s "$BAM" && -s "$L/bowtie.log" ]]; then
    log "ASSEMBLY_RECOUNT_ONLY=1: keeping the existing mapping ($BAM), recounting only"
  else
    if [[ "${ASSEMBLY_RECOUNT_ONLY:-0}" == 1 ]]; then
      warn "ASSEMBLY_RECOUNT_ONLY=1 but $( [[ -s "$BAM" ]] && echo "no $L/bowtie.log" || echo "no BAM $BAM" ) — mapping again"
    else
      log "ASSEMBLY_RECOUNT_ONLY=${ASSEMBLY_RECOUNT_ONLY:-0} — full mapping"
    fi
    log "bowtie index of the assembly ..."
    ensure_assembly_bowtie_index

    # One mapping for all assembly-based modules (novel, fragments, bins, compare):
    # uncollapsed non-human reads, bowtie -v <mm> --best --strata and either all hits (-a) or at most
    # ASSEMBLY_MAX_HITS hits per read. All hits of a read are consecutive in bowtie's SAM output, so the
    # number of hits is added as NH:i:<n> on the fly; downstream tools weight each alignment 1/NH.
    local mm=${ASSEMBLY_MISMATCHES:-${NOVEL_MISMATCHES:-1}} kopt=(-a)
    [[ "${ASSEMBLY_MAX_HITS:-0}" -gt 0 ]] && kopt=(-k "$ASSEMBLY_MAX_HITS")
    log "mapping non-human reads (bowtie -v $mm ${kopt[*]} --best --strata) ..."
    bowtie -v "$mm" "${kopt[@]}" --best --strata -S -p "$THREADS" "$ASSEMBLY" "$READS" 2> "$L/bowtie.log" \
    | awk -v mx="${ASSEMBLY_MAX_HITS:-0}" -v cnt="$L/reads_at_max_hits.txt" 'BEGIN{FS=OFS="\t"}
        function flush(   i) { if (mx > 0 && n >= mx) atmax++
                               for (i = 1; i <= n; i++) print buf[i] "\tNH:i:" n; n = 0 }
        /^@/ { print; next }
        int($2/4)%2 == 1 { flush(); q = ""; print; next }
        $1 != q { flush(); q = $1 }
        { buf[++n] = $0 }
        END { flush(); print atmax + 0 > cnt }' \
    | samtools view -@ "$THREADS" -b - \
    | samtools sort -@ "$THREADS" -T "$TEMPDIR/sort" -o "$BAM" 2> "$L/samtools_sort.log"
    samtools index "$BAM"
    grep -q 'reads processed' "$L/bowtie.log" || die "bowtie failed (see $L/bowtie.log)"
  fi

  log "feature, category and functional counts ..."
  assembly_counts "$L/$PREFIX.counts.log"

  {
    echo -e "sample\tcategory\tvalue"
    echo -e "$PREFIX\tinput_reads\t$(awk -F': ' '/reads processed/{print $2+0}' "$L/bowtie.log")"
    echo -e "$PREFIX\tmapped_reads\t$(awk -F': ' '/reads with at least one/{split($2,x," "); print x[1]+0}' "$L/bowtie.log")"
    echo -e "$PREFIX\tunmapped_reads\t$(awk -F': ' '/reads that failed to align/{split($2,x," "); print x[1]+0}' "$L/bowtie.log")"
    echo -e "$PREFIX\treads_hitting_max_hits\t$(cat "$L/reads_at_max_hits.txt" 2>/dev/null || echo NA)"
    echo -e "$PREFIX\talignments\t$(samtools view -c -F 4 "$BAM")"
    echo -e "$PREFIX\tmismatches\t$mm"
    echo -e "$PREFIX\tmax_hits\t${ASSEMBLY_MAX_HITS:-all}"
  } > "$P.mapping_summary.tsv"
  log "mapping summary: $P.mapping_summary.tsv"
}

# ==================================================================
# MODULE: novel — discovery of novel ncRNA loci on the assembly
# ==================================================================
module_novel() {
  need samtools bedtools seqkit sort cut uniq
  [[ "$RUN_ARAGORN"  -eq 1 ]] && need aragorn
  [[ "$RUN_TRNASCAN" -eq 1 ]] && need tRNAscan-SE
  [[ "$RUN_CMSCAN"   -eq 1 ]] && { need cmscan; [[ -s "$RFAM_CM" ]] || die "RFAM_CM not found: $RFAM_CM"; }
  [[ "$RUN_RNAFOLD"  -eq 1 ]] && need RNAfold

  local L=$LOGDIR/novel OUT=$NOVELDIR
  mkdir -p "$L" "$OUT"
  local BAM=$ASMDIR/$PREFIX.bam
  [[ -s "$BAM" ]] || die "missing $BAM — run module assembly first"
  [[ -s "$ASSEMBLY.gz" && -s "$GFF" ]] || die "missing assembly/annotation — run module prep first"
  make_tempdir
  local P=$OUT/$PREFIX

  # uncompressed reference for ARAGORN / tRNAscan / getfasta (private copy)
  local REF=$TEMPDIR/$PREFIX.assembly.fa
  zcat "$ASSEMBLY.gz" > "$REF"

  # rRNA genes missing from the annotation (antisense artefacts): cmscan with the Rfam rRNA models on contigs whose
  # CDS have (almost) only antisense reads; hits are added as rRNA features and the assembly counts are redone
  local RESC=$ASMDIR/$PREFIX.rRNA_rescue.bed FC=$ASMDIR/$PREFIX.feature_counts.tsv
  if [[ "${RRNA_RESCUE:-1}" -eq 1 && -s "$FC" && -s "$RFAM_CM" ]] && command -v cmscan >/dev/null && command -v cmpress >/dev/null; then
    log "[0/14] rRNA rescue: cmscan (Rfam rRNA models) on contigs with antisense-only CDS ..."
    local before; before=$(awk -F'\t' 'NR>1 && $1=="antisense_CDS"{print $NF}' "$ASMDIR/$PREFIX.category_counts.tsv" 2>/dev/null)
    python3 "$SCRIPTS_DIR/srna_rrna_rescue.py" --feature-counts "$FC" --fasta "$REF" --rfam-cm "$RFAM_CM" \
      ${RFAM_FAMILY_TXT:+--family-txt "$RFAM_FAMILY_TXT"} --annotation "$GBED" --out "$RESC" --cpu "$THREADS" \
      --cm-cache "${RRNA_RESCUE_CM:-$OUTDIR/.shared/Rfam.rRNA.cm}" \
      --min-antisense "${RRNA_RESCUE_MIN_ANTISENSE:-20}" --max-sense-ratio "${RRNA_RESCUE_MAX_SENSE_RATIO:-0.1}" \
      > "$L/$PREFIX.rrna_rescue.log" 2>&1 || die "srna_rrna_rescue.py failed (see $L/$PREFIX.rrna_rescue.log)"
    sed 's/^/    | /' "$L/$PREFIX.rrna_rescue.log"
    if [[ -s "$RESC" ]]; then
      log "[0/14] recounting the assembly mapping with the rescued rRNA genes ..."
      assembly_counts "$LOGDIR/assembly/$PREFIX.counts.log"
      log "[0/14] antisense_CDS (without rRNA reads): ${before:-?} -> $(awk -F'\t' 'NR>1 && $1=="antisense_CDS"{print $NF}' "$ASMDIR/$PREFIX.category_counts.tsv")"
    fi
    # annotation with the rescued genes of all types (GFF3 next to the original; the original is not changed)
    local rl=() f
    for f in "$SAMPLE_DIR"/Analysis/*/"${ASSEMBLY_SUBDIR:+$ASSEMBLY_SUBDIR/}$PREFIX.rRNA_rescue.bed"; do if [[ -s "$f" ]]; then rl+=("$f"); fi; done
    if [[ ${#rl[@]} -gt 0 ]]; then
      acquire_lock "$GFF.rescued.lock" 600
      python3 "$SCRIPTS_DIR/srna_rescued_gff.py" --gff "$GFF" --rescue "$(IFS=,; echo "${rl[*]}")" \
        --out "${GFF%.gff}.rRNA_rescued.gff" 2>> "$L/$PREFIX.rrna_rescue.log" \
        && log "[0/14] annotation with rescued rRNA genes: ${GFF%.gff}.rRNA_rescued.gff" \
        || warn "srna_rescued_gff.py failed (see $L/$PREFIX.rrna_rescue.log)"
      release_lock "$GFF.rescued.lock"
    fi
  elif [[ "${RRNA_RESCUE:-1}" -eq 1 ]]; then
    warn "rRNA rescue skipped (needs cmscan/cmpress, RFAM_CM and $FC)"
  fi

  # reads come from the assembly module's BAM (one mapping for all modules); every alignment
  # gets the weight 1/NH, so a multi-mapping read adds up to one read and repeated loci
  # (rRNA operons, multi-copy tRNAs) are not inflated
  # rRNA reads (RiboDetector, module nonhuman) are left out (NOVEL_EXCLUDE_RRNA=1): they otherwise form
  # candidate loci on unannotated / partial rRNA genes (about a third of all candidates), cost cmscan time
  # and are not novel. tRNA / tmRNA loci are already excluded by the annotation + ARAGORN.
  local RRIDS=$TEMPDIR/rrna.ids RR=$NONHUMANDIR/rrna/$PREFIX.rRNA.fastq.gz
  : > "$RRIDS"
  if [[ "${NOVEL_EXCLUDE_RRNA:-1}" -eq 1 ]]; then
    if [[ -s "$RR" ]]; then
      zcat "$RR" | awk 'NR % 4 == 1 { n = substr($1, 2); sub(/\/[12]$/, "", n); print n }' > "$RRIDS"
      log "[1-4/14] leaving out $(wc -l < "$RRIDS") RiboDetector rRNA reads"
    else
      warn "no RiboDetector rRNA reads ($RR) — rRNA reads are not excluded from the novel loci"
    fi
    local RS; RS=$(rfam_rrna_seqs)
    if [[ -n "$RS" ]]; then
      zcat "$RS" > "$TEMPDIR/rrna.seqs"
      log "[1-4/14] + reads with an Rfam rRNA hit ($(( $(wc -l < "$TEMPDIR/rrna.seqs") / 2 )) sequences)"
    fi
  fi
  [[ -e "$TEMPDIR/rrna.seqs" ]] || : > "$TEMPDIR/rrna.seqs"
  log "[1-4/14] alignments from $BAM (weight 1/NH) ..."
  samtools view -@ "$THREADS" -F 4 "$BAM" \
  | awk -v plus="$TEMPDIR/plus.bed" -v minus="$TEMPDIR/minus.bed" -v ids="$RRIDS" -v seqs="$TEMPDIR/rrna.seqs" 'BEGIN{OFS="\t"
        while ((getline n < ids) > 0) rr[n] = 1
        while ((getline n < seqs) > 0) rs[n] = 1 }
      ($1 in rr) || ($10 in rs) { next }
      { nh = 1; for (i = 12; i <= NF; i++) if ($i ~ /^NH:i:/) { nh = substr($i, 6) + 0; break }
        sd = (int($2/16)%2 == 1) ? "-" : "+"
        print $3, $4 - 1, $4 - 1 + length($10), $1, 1 / nh, sd > (sd == "-" ? minus : plus) }'
  touch "$TEMPDIR/plus.bed" "$TEMPDIR/minus.bed"
  sort -k1,1 -k2,2n -T "$TEMPDIR" "$TEMPDIR/plus.bed"  -o "$P.plus.bed"
  sort -k1,1 -k2,2n -T "$TEMPDIR" "$TEMPDIR/minus.bed" -o "$P.minus.bed"

  log "[5/14] clustering reads into candidate loci ..."
  local s
  for s in plus minus; do
    if [[ -s "$P.$s.bed" ]]; then
      bedtools merge -s -d "$NOVEL_MERGEDIST" -i "$P.$s.bed" -c 6 -o distinct \
        | awk -v s="$s" 'BEGIN{OFS="\t"} {print $1,$2,$3,"cluster_"s"_"NR,0,$4}' > "$P.$s.clusters.bed"
    else
      : > "$P.$s.clusters.bed"
    fi
  done
  cat "$P.plus.clusters.bed" "$P.minus.clusters.bed" | sort -k1,1 -k2,2n > "$P.clusters.bed"

  log "[6/14] filtering loci by length ($NOVEL_MINLEN-$NOVEL_MAXLEN nt) ..."
  awk -v minlen="$NOVEL_MINLEN" -v maxlen="$NOVEL_MAXLEN" '($3-$2)>=minlen && ($3-$2)<=maxlen' \
    "$P.clusters.bed" > "$P.clusters.lenfiltered.bed"
  awk '$6=="+"' "$P.clusters.lenfiltered.bed" > "$P.clusters.lenfiltered.plus.bed"
  awk '$6=="-"' "$P.clusters.lenfiltered.bed" > "$P.clusters.lenfiltered.minus.bed"
  for s in plus minus; do
    if [[ -s "$P.clusters.lenfiltered.$s.bed" ]]; then
      # col 7 = weighted reads (sum of 1/hits), col 8 = alignments
      bedtools map -a "$P.clusters.lenfiltered.$s.bed" -b "$P.$s.bed" -c 5,5 -o sum,count \
        | awk 'BEGIN{FS=OFS="\t"} {if ($7==".") $7=0; printf "%s\t%s\t%s\t%s\t%s\t%s\t%.2f\t%d\n", $1,$2,$3,$4,$5,$6,$7,$8}' \
        > "$P.clusters.counts.$s.bed"
    else
      : > "$P.clusters.counts.$s.bed"
    fi
  done
  cat "$P.clusters.counts.plus.bed" "$P.clusters.counts.minus.bed" | sort -k1,1 -k2,2n > "$P.clusters.counts.bed"

  log "[7/14] filtering loci by read support (>= $NOVEL_MINREADS) ..."
  awk -v minreads="$NOVEL_MINREADS" '$7>=minreads' "$P.clusters.counts.bed" > "$P.clusters.supported.bed"

  log "[8/14] GFF -> BED annotation ..."
  awk 'BEGIN{FS=OFS="\t"} $0 !~ /^#/ && NF>=8 && $3!="region" && $3!="source" {print $1,$4-1,$5,$3,".",$7}' "$GFF" \
  | awk 'BEGIN{OFS="\t"} NF>=6 { if ($2>$3){t=$2;$2=$3;$3=t} if ($2<0) $2=0; print }' \
  | sort -k1,1 -k2,2n > "$P.annot.bed"
  if [[ -s "$RESC" ]]; then
    awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"rRNA",".",$6}' "$RESC" | cat "$P.annot.bed" - | sort -k1,1 -k2,2n > "$TEMPDIR/annot.bed" \
      && mv "$TEMPDIR/annot.bed" "$P.annot.bed"
  fi

  if [[ "$RUN_ARAGORN" -eq 1 ]]; then
    log "[9/14] running ARAGORN ..."
    aragorn -gcstd -w -o "$P.aragorn.txt" "$REF"
    awk 'BEGIN{OFS="\t"}
      /^>/ { seq=$1; sub(/^>/,"",seq) }
      /tRNA|tmRNA/ {
        strand="+"; if ($0 ~ /c\[/) strand="-"
        if (match($0, /\[[0-9]+,[0-9]+\]/)) {
          split(substr($0, RSTART+1, RLENGTH-2), a, ",")
          s=a[1]+0; e=a[2]+0; if (s>e){t=s;s=e;e=t}
          type="ARAGORN_RNA"; if ($0 ~ /tmRNA/) type="tmRNA"; else if ($0 ~ /tRNA/) type="tRNA"
          print seq, s-1, e, type, ".", strand
        }
      }' "$P.aragorn.txt" | sort -k1,1 -k2,2n > "$P.aragorn.bed"
    if [[ -s "$P.aragorn.bed" ]]; then
      cat "$P.annot.bed" "$P.aragorn.bed" | sort -k1,1 -k2,2n > "$TEMPDIR/annot.bed" && mv "$TEMPDIR/annot.bed" "$P.annot.bed"
    fi
  else
    log "[9/14] ARAGORN skipped"
  fi

  if [[ "$RUN_TRNASCAN" -eq 1 ]]; then
    log "[10/14] running tRNAscan-SE ..."
    tRNAscan-SE --thread "$THREADS" -o "$P.trnascan.out" "$REF" > "$L/trnascan.log" 2>&1
    awk 'BEGIN{FS=OFS="\t"} NR>3 && $1 !~ /^Sequence/ && NF>=9 {
           seq=$1; s=$3; e=$4; strand="+"; if (s>e){t=s;s=e;e=t;strand="-"}
           gsub(/ /,"",seq); print seq, s-1, e, "tRNA", ".", strand }' "$P.trnascan.out" \
      | sort -k1,1 -k2,2n > "$P.trnascan.bed" || true
    if [[ -s "$P.trnascan.bed" ]]; then
      cat "$P.annot.bed" "$P.trnascan.bed" | sort -k1,1 -k2,2n > "$TEMPDIR/annot.bed" && mv "$TEMPDIR/annot.bed" "$P.annot.bed"
    fi
  else
    log "[10/14] tRNAscan-SE skipped"
  fi

  log "[11/14] identifying candidate novel loci ..."
  awk '$4=="CDS"' "$P.annot.bed" > "$P.cds.bed"
  # intergenic: no same-strand overlap with any annotated feature
  # stranded reads: a cluster is a candidate only where no annotated feature lies on its own strand;
  # "intergenic" = no feature on either strand, "antisense to CDS" = opposite a CDS (own strand free)
  bedtools intersect -s -v -a "$P.clusters.supported.bed" -b "$P.annot.bed" > "$TEMPDIR/nosense.bed"
  bedtools intersect -v -a "$TEMPDIR/nosense.bed" -b "$P.annot.bed" > "$P.candidates.intergenic.bed"
  # antisense to CDS: overlap a CDS on the opposite strand
  bedtools intersect -S -wa -wb -a "$TEMPDIR/nosense.bed" -b "$P.cds.bed" > "$P.candidates.antisense_to_CDS.full.bed"
  cut -f1-8 "$P.candidates.antisense_to_CDS.full.bed" | sort -k1,1 -k2,2n -k3,3n -k6,6 | uniq \
    > "$P.candidates.antisense_to_CDS.bed"

  log "[12/14] extracting candidate sequences ..."
  : > "$P.candidates.intergenic.fa"; : > "$P.candidates.antisense_to_CDS.fa"
  [[ -s "$P.candidates.intergenic.bed" ]] && \
    bedtools getfasta -s -fi "$REF" -bed "$P.candidates.intergenic.bed" -fo "$P.candidates.intergenic.fa"
  [[ -s "$P.candidates.antisense_to_CDS.bed" ]] && \
    bedtools getfasta -s -fi "$REF" -bed "$P.candidates.antisense_to_CDS.bed" -fo "$P.candidates.antisense_to_CDS.fa"
  cat "$P.candidates.intergenic.fa" "$P.candidates.antisense_to_CDS.fa" > "$P.candidates.all.fa"

  if [[ -s "$P.candidates.all.fa" ]]; then
    if [[ "$RUN_CMSCAN" -eq 1 ]]; then
      log "[13/14] cmscan (Rfam) on candidates ..."
      # shellcheck disable=SC2086
      cmscan --cpu "$THREADS" ${CMSCAN_OPTS:---rfam --cut_ga} --tblout "$P.candidates.rfam.tbl" \
        "$RFAM_CM" "$P.candidates.all.fa" > "$P.candidates.cmscan.out"
      { awk '$1 !~ /^#/ {print $3}' "$P.candidates.rfam.tbl" | sort -u > "$P.rfam_hit_ids.txt"; } || true
      if [[ -s "$P.rfam_hit_ids.txt" ]]; then
        seqkit grep -v -f "$P.rfam_hit_ids.txt" "$P.candidates.all.fa" > "$P.candidates.no_rfam.fa"
      else
        cp "$P.candidates.all.fa" "$P.candidates.no_rfam.fa"
      fi
    else
      log "[13/14] cmscan skipped"
      cp "$P.candidates.all.fa" "$P.candidates.no_rfam.fa"
    fi

    if [[ -s "$P.candidates.no_rfam.fa" && "$RUN_RNAFOLD" -eq 1 ]]; then
      log "[14/14] RNAfold on candidates without Rfam hit ..."
      RNAfold --noPS --jobs="$THREADS" < "$P.candidates.no_rfam.fa" > "$P.candidates.no_rfam.RNAfold.txt"
    else
      log "[14/14] RNAfold skipped (disabled or no candidates)"
    fi
  else
    log "[13+14/14] no candidate sequences — cmscan/RNAfold skipped"
  fi

  {
    echo -e "sample\tcategory\treads"
    # bowtie -a --best --strata reports every equally good hit: alignments >> reads
    echo -e "$PREFIX\tinput_reads\t$(awk -F': ' '/reads processed/{print $2+0}' "$LOGDIR/assembly/bowtie.log")"
    echo -e "$PREFIX\tmapped_reads\t$(awk -F': ' '/reads with at least one/{split($2,x," "); print x[1]+0}' "$LOGDIR/assembly/bowtie.log")"
    echo -e "$PREFIX\tmapped_alignments\t$(samtools view -c -F 4 "$BAM")"
    echo -e "$PREFIX\tnovel_intergenic_candidates\t$(grep -cve '^\s*$' "$P.candidates.intergenic.bed" || true)"
    echo -e "$PREFIX\tnovel_antisense_candidates\t$(grep -cve '^\s*$' "$P.candidates.antisense_to_CDS.bed" || true)"
  } > "$P.assembly.summary.tsv"

  if [[ "$KEEP_INTERMEDIATE" -eq 0 ]]; then
    rm -f "$P.plus.bed" "$P.minus.bed" "$P.plus.clusters.bed" "$P.minus.clusters.bed" \
          "$P.clusters.bed" "$P.clusters.lenfiltered.bed"
  fi
  log "summary: $P.assembly.summary.tsv"
}

# ==================================================================
# driver
# ==================================================================
# ==================================================================
# MODULE: contam — library QC on raw/trimmed reads
#   adapter dimers, cutadapt fate, complexity, overrepresented sequences,
#   reads hitting contaminant references (CONTAM_FA)
# ==================================================================
module_contam() {
  need bowtie bowtie-build samtools
  local L=$LOGDIR/contam
  mkdir -p "$L"
  [[ -s "$TRIMMEDREADS" && -s "$RAWREADS" ]] || die "missing raw/trimmed reads — run module qc first"
  make_tempdir
  local QC=$STATS_DIR/$PREFIX.$TYPE.qc_summary.tsv
  local OVR=$STATS_DIR/$PREFIX.$TYPE.overrepresented.tsv
  local CUTLOG=$STATS_DIR/$PREFIX.$TYPE.cutadapt.log

  # dimer = the read starts with the adapter, or with the adapter missing up to ADAPTER_3P_TRUNC nt
  local tr=${ADAPTER_3P_TRUNC:-0}
  log "adapter dimers in raw reads (adapter core ${ADAPTER_3P:$tr:12} within the first $((tr + 1)) nt) ..."
  local raw_n dimer_n
  read -r raw_n dimer_n < <(zcat "$RAWREADS" | awk -v a="${ADAPTER_3P:$tr:12}" -v m=$((tr + 1)) \
    'NR%4==2{n++; p=index($0,a); if (p>=1 && p<=m) d++} END{print n+0, d+0}')

  cut_val() { { grep -m1 "$1" "$CUTLOG" 2>/dev/null | awk -F: '{gsub(/[ ,]/,"",$2); sub(/\(.*/,"",$2); print $2+0}'; } || true; }
  local c_tot c_ad c_short c_long c_pass
  c_tot=$(cut_val "Total reads processed"); c_ad=$(cut_val "Reads with adapters")
  c_short=$(cut_val "Reads that were too short"); c_long=$(cut_val "Reads that were too long")
  c_pass=$(cut_val "Reads written (passing filters)")

  log "library complexity / overrepresented sequences ..."
  zcat "$TRIMMEDREADS" | awk 'NR%4==2' | sort -S "${QC_SORT_MEM:-2G}" -T "$TEMPDIR" | uniq -c \
    | sort -k1,1nr -S "${QC_SORT_MEM:-2G}" -T "$TEMPDIR" > "$TEMPDIR/seqcounts"
  local trim_n distinct singletons top10
  read -r trim_n distinct singletons top10 < <(awk '{t+=$1; d++; if($1==1) s++; if(NR<=10) top+=$1}
      END{print t+0, d+0, s+0, top+0}' "$TEMPDIR/seqcounts")
  awk -v t="$trim_n" -v n="$QC_TOP_SEQS" 'BEGIN{OFS="\t"; print "rank","count","fraction","length","sequence"}
      NR<=n{print NR, $1, (t? $1/t : 0), length($2), $2}' "$TEMPDIR/seqcounts" > "$OVR"

  local rows=()
  local fa idx db hits
  for fa in ${CONTAM_FA:-}; do
    db=$(basename "$fa"); db=${db%.gz}; db=${db%.fa}; db=${db%.fasta}; db=${db%.fna}
    log "contaminant screen: $db ..."
    ensure_host_index "$fa" "" "contaminant $db"
    idx=$HOST_IDX
    bowtie -p "$THREADS" -v "$CONTAM_MISMATCHES" -k 1 --best --sam "$idx" "$TRIMMEDREADS" \
      2> "$L/$PREFIX.bowtie.$db.log" | samtools view -F 4 - | cut -f3 | sort -T "$TEMPDIR" | uniq -c \
      | sort -k1,1nr > "$TEMPDIR/$db.hits" || true
    hits=$(awk '{s+=$1} END{print s+0}' "$TEMPDIR/$db.hits")
    awk -v t="$trim_n" 'BEGIN{OFS="\t"; print "target","reads","fraction_of_trimmed"} {print $2,$1,(t?$1/t:0)}' \
      "$TEMPDIR/$db.hits" > "$STATS_DIR/$PREFIX.$TYPE.contaminants.$db.tsv"
    rows+=("contaminant_$db" "$hits")
  done

  frac() { awk -v a="$1" -v b="$2" 'BEGIN{printf "%.4f", (b>0 ? a/b : 0)}'; }
  {
    echo -e "metric\tvalue\tfraction"
    echo -e "raw_reads\t$raw_n\t1"
    echo -e "adapter_dimers_raw\t$dimer_n\t$(frac "$dimer_n" "$raw_n")"
    echo -e "cutadapt_reads_with_adapter\t${c_ad:-NA}\t$(frac "${c_ad:-0}" "${c_tot:-0}")"
    echo -e "cutadapt_too_short\t${c_short:-NA}\t$(frac "${c_short:-0}" "${c_tot:-0}")"
    echo -e "cutadapt_too_long\t${c_long:-NA}\t$(frac "${c_long:-0}" "${c_tot:-0}")"
    echo -e "cutadapt_passing\t${c_pass:-NA}\t$(frac "${c_pass:-0}" "${c_tot:-0}")"
    echo -e "trimmed_reads\t$trim_n\t$(frac "$trim_n" "$raw_n")"
    echo -e "distinct_sequences\t$distinct\t$(frac "$distinct" "$trim_n")"
    echo -e "singleton_sequences\t$singletons\t$(frac "$singletons" "$distinct")"
    echo -e "top10_sequences_reads\t$top10\t$(frac "$top10" "$trim_n")"
    local i
    for ((i=0; i<${#rows[@]}; i+=2)); do
      echo -e "${rows[i]}\t${rows[i+1]}\t$(frac "${rows[i+1]}" "$trim_n")"
    done
  } > "$QC"
  [[ -n "${CONTAM_FA:-}" ]] || log "CONTAM_FA empty — contaminant screen skipped"
  log "QC summary: $QC"
}

# ==================================================================
# MODULE: fragments — tRNA / rRNA / tmRNA fragment analysis
#   uses the bowtie alignments of the novel module (all non-human reads)
# ==================================================================
module_fragments() {
  need samtools python3
  local L=$LOGDIR/fragments BAM=$ASMDIR/$PREFIX.bam
  mkdir -p "$L" "$FRAGDIR"
  [[ -s "$BAM" ]] || die "missing $BAM — run module novel first"
  [[ -s "$GBED" ]] || die "missing $GBED — run module prep first"
  make_tempdir
  local FBED; FBED=$(annotation_bed)            # + rRNA genes rescued by cmscan (module novel)
  # no tRNAs in the annotation? use the ARAGORN tRNA/tmRNA predictions of the novel module
  if ! awk -F'\t' '$8=="tRNA"{f=1; exit} END{exit !f}' "$GBED" && [[ -s "$NOVELDIR/$PREFIX.aragorn.bed" ]]; then
    log "no tRNA features in the annotation — adding ARAGORN predictions"
    local base=$FBED; FBED=$TEMPDIR/features.bed
    { cat "$base"
      awk 'BEGIN{OFS="\t"} {id=sprintf("aragorn_%s_%d_%d_%s",$1,$2,$3,$6);
            print $1,$2,$3,id,".",$6,"ARAGORN",$4,".","ID="id";product="$4}' "$NOVELDIR/$PREFIX.aragorn.bed"
    } > "$FBED"
  fi
  local binopt=()
  [[ -n "$IN_BINS" ]] && binopt=(--contig2bin "$IN_BINS")
  log "tRNA/rRNA fragment analysis ..."
  python3 "$SCRIPTS_DIR/srna_fragments.py" --bam "$BAM" --gbed "$FBED" --out-prefix "$FRAGDIR/$PREFIX" \
    --types "$FRAG_TYPES" --tol "$FRAG_TOL" --ac-from "$FRAG_AC_FROM" --ac-to "$FRAG_AC_TO" \
    --min-overlap "$FRAG_MIN_OVERLAP" --min-reads "$FRAG_MIN_READS" --threads "$THREADS" \
    ${binopt[@]+"${binopt[@]}"} 2> "$L/$PREFIX.srna_fragments.log"
  log "fragment classes: $FRAGDIR/$PREFIX.fragments.classes.tsv"
  log "RNA integrity indicators ..."
  python3 "$SCRIPTS_DIR/srna_integrity.py" --frag-prefix "$FRAGDIR/$PREFIX" --asm-prefix "$ASMDIR/$PREFIX" \
    --out-prefix "$FRAGDIR/$PREFIX" --sample "$PREFIX" 2>> "$L/$PREFIX.srna_fragments.log" \
    || warn "integrity indicators failed (see $L/$PREFIX.srna_fragments.log)"
  log "integrity: $FRAGDIR/$PREFIX.integrity.tsv"
}

# ==================================================================
# MODULE: bins — counts per contig, feature, bin (MAG) and taxon
# ==================================================================
module_bins() {
  need samtools python3
  local L=$LOGDIR/bins BAM=$ASMDIR/$PREFIX.bam
  mkdir -p "$L" "$BINSDIR"
  [[ -s "$BAM" ]] || die "missing $BAM — run module assembly first"
  [[ -s "$GBED" ]] || die "missing $GBED — run module prep first"
  local opt=()
  [[ -n "$IN_BINS" ]] && opt+=(--contig2bin "$IN_BINS") || warn "no contig2bin given — all contigs 'unbinned'"
  [[ -n "$IN_GTDB" ]] && opt+=(--gtdbtk "$IN_GTDB")
  make_tempdir
  local gb; gb=$(annotation_bed)                # + rRNA genes rescued by cmscan (module novel)
  log "counting reads per contig / feature / bin / taxon ..."
  python3 "$SCRIPTS_DIR/srna_bins.py" --bam "$BAM" --gbed "$gb" --out-prefix "$BINSDIR/$PREFIX" \
    --min-overlap "$BINS_MIN_OVERLAP" --threads "$THREADS" ${opt[@]+"${opt[@]}"} 2> "$L/$PREFIX.srna_bins.log"
  log "bin counts: $BINSDIR/$PREFIX.bin_counts.tsv"
}

# ==================================================================
# MODULE: crispr — sRNAs on CRISPR arrays (crRNAs)
# ==================================================================
module_crispr() {
  need samtools python3
  local L=$LOGDIR/crispr BAM=$ASMDIR/$PREFIX.bam
  mkdir -p "$L" "$CRISPRDIR"
  [[ -s "$BAM" ]] || die "missing $BAM — run module assembly first"
  [[ -s "$GBED" ]] || die "missing $GBED — run module prep first"
  local opt=()
  [[ -n "$IN_BINS" ]] && opt+=(--contig2bin "$IN_BINS")
  [[ -n "$IN_GTDB" ]] && opt+=(--gtdbtk "$IN_GTDB")
  [[ -s "$ASMDIR/$PREFIX.feature_counts.tsv" ]] && opt+=(--feature-counts "$ASMDIR/$PREFIX.feature_counts.tsv")
  log "reads on CRISPR arrays ..."
  python3 "$SCRIPTS_DIR/srna_crispr.py" --bam "$BAM" --gbed "$GBED" --out-prefix "$CRISPRDIR/$PREFIX" \
    --threads "$THREADS" ${opt[@]+"${opt[@]}"} 2> "$L/$PREFIX.srna_crispr.log" \
    || die "srna_crispr.py failed (see $L/$PREFIX.srna_crispr.log)"
  tail -n 1 "$L/$PREFIX.srna_crispr.log" | sed 's/^/    /'
  log "CRISPR arrays: $CRISPRDIR/$PREFIX.crispr_arrays.tsv"
}

# ==================================================================
# COMPARE (per sample, --compare): sRNA types vs metaT / metaG
# ==================================================================
run_compare() {   # sample ; uses CMP_* arrays filled by the driver
  local S=$1
  local SDIR=$OUTDIR/$S CDIR=$OUTDIR/$S/Analysis/comparison L=$OUTDIR/$S/log_cluster_jobs/compare
  mkdir -p "$CDIR" "$L"
  local args=() t n=0
  for t in ${CMP_TYPES[$S]}; do
    if [[ -s "$SDIR/Analysis/$t/bins/$S.feature_counts.tsv" ]]; then
      args+=(--srna "$t=$SDIR/Analysis/$t/bins/$S"); n=$((n+1))
    else
      warn "$S/$t: no bins output (run module bins) — type left out of the comparison"
    fi
  done
  [[ -n "${CMP_MTF[$S]:-}" ]] && { args+=(--metat-features "${CMP_MTF[$S]}"); n=$((n+1)); }
  [[ -n "${CMP_MTC[$S]:-}" ]] && args+=(--metat-contigs "${CMP_MTC[$S]}")
  [[ -z "${CMP_MTF[$S]:-}" && -n "${CMP_MTC[$S]:-}" ]] && n=$((n+1))
  [[ -n "${CMP_MGC[$S]:-}" ]] && args+=(--metag-contigs "${CMP_MGC[$S]}")
  [[ -z "${CMP_MGF[$S]:-}${CMP_MTF[$S]:-}" && -n "${CMP_MGC[$S]:-}${CMP_MTC[$S]:-}" ]] && n=$((n+1))
  # metaG: explicit column 10, or mg.ID tables auto-detected in metaT_features by srna_compare.py
  if [[ -n "${CMP_MGF[$S]:-}" ]]; then
    args+=(--metag-features "${CMP_MGF[$S]}"); n=$((n+1))
  elif [[ -n "${CMP_MTF[$S]:-}" && -z "${CMP_MGC[$S]:-}" ]]; then
    n=$((n+1))   # possibly; the python script skips the sample itself if no mg tables are found
  fi
  [[ -s "$SDIR/Analysis/annotation/$S.bed" ]] && args+=(--gbed "$SDIR/Analysis/annotation/$S.bed")
  [[ -n "${CMP_BIN[$S]:-}" ]] && args+=(--contig2bin "${CMP_BIN[$S]}")
  [[ -n "${CMP_GTDB[$S]:-}" ]] && args+=(--gtdbtk "${CMP_GTDB[$S]}")
  if [[ $n -lt 2 ]]; then
    log "$S: fewer than two sources available (sRNA types: ${CMP_TYPES[$S]}; metaT: ${CMP_MTF[$S]:-${CMP_MTC[$S]:-none}}; metaG: ${CMP_MGF[$S]:-none}) — skipped"
    return 0
  fi
  local mg=""; [[ -n "${CMP_MGF[$S]:-}${CMP_MGC[$S]:-}" ]] && mg=" metaG"; [[ -z "$mg" && -n "${CMP_MTF[$S]:-}" ]] && mg=" metaG(auto-detected from metaT_features)"
  log "$S: comparing ${CMP_TYPES[$S]}${CMP_MTF[$S]:+ metaT}$mg ..."
  set +e
  (
    set -euo pipefail
    activate_env "${ENV_COMPARE:-${ENV_DEFAULT:-}}"
    python3 "$SCRIPTS_DIR/srna_compare.py" --sample "$S" --outdir "$CDIR" \
      --pseudocount "$COMPARE_PSEUDOCOUNT" --min-count "$COMPARE_MIN_COUNT" \
      --omics-bin-basis "${COMPARE_BIN_BASIS:-contig,cds,no_rrna}" "${args[@]}" \
      > "$L/$S.compare.log" 2>&1
  )
  local rc=$?
  set -e
  if [[ $rc -eq 0 ]]; then
    log "$S: comparison written to $CDIR"
  else
    log "$S: comparison FAILED — see $L/$S.compare.log"; tail -n 15 "$L/$S.compare.log" | sed 's/^/    | /'
    CMP_FAILED=1
  fi
}

# ==================================================================
# REPORT (per sample, --report): PDF with the key numbers and plots of all types
# ==================================================================
# ---------------------------------------------------------------- automated Methods (srna_methods.py)
# write_params FILE: the resolved configuration (every NAME= of the config file) as NAME<TAB>VALUE
write_params() {
  local f=$1 v
  : > "$f"
  for v in $(grep -oE '^[A-Z_][A-Z0-9_]*=' "$CONFIG" 2>/dev/null | tr -d = | sort -u); do
    printf '%s\t%s\n' "$v" "${!v-}" >> "$f"
  done
}
# methods_envs: environment prefixes whose conda-meta gives the tool versions
methods_envs() {
  local e out=() rd
  for e in "${ENV_DEFAULT:-}" "${ENV_NOVEL:-}" "${ENV_QC:-}" "${ENV_NONHUMAN:-}" "${ENV_ASSEMBLY:-}"; do
    [[ -n "$e" && -d "$e/conda-meta" ]] && out+=("$e")
  done
  rd=$(dirname "$(dirname "${RIBO_DETECTOR_CPU:-/x/x}")"); [[ -d "$rd/conda-meta" ]] && out+=("$rd")
  if [[ "${MULTIQC:-}" == */* ]]; then rd=$(dirname "$(dirname "$MULTIQC")"); [[ -d "$rd/conda-meta" ]] && out+=("$rd"); fi
  [[ -n "${CONDA_PREFIX:-}" && -d "$CONDA_PREFIX/conda-meta" ]] && out+=("$CONDA_PREFIX")
  local IFS=,; echo "${out[*]}"
}
# write_methods MODE OUT_PREFIX [extra args]: one Markdown file per style in METHODS_STYLE (nature,star; none = off)
write_methods() {
  local mode=$1 pre=$2; shift 2
  [[ ( "${METHODS_STYLE:-nature,star}" == none || -z "${METHODS_STYLE-x}" ) && \
     ( -z "${METHODS_REPORT_STYLE:-}" || "${METHODS_REPORT_STYLE}" == none ) ]] && return 0
  local pf=$pre.params.tsv st styles=${METHODS_STYLE:-nature,star}
  [[ -n "${METHODS_REPORT_STYLE:-}" && "${METHODS_REPORT_STYLE}" != none && ",$styles," != *",${METHODS_REPORT_STYLE},"* ]] \
    && styles+=",$METHODS_REPORT_STYLE"
  write_params "$pf"
  for st in ${styles//,/ }; do
    python3 "$SCRIPTS_DIR/srna_methods.py" --mode "$mode" --style "$st" --params "$pf" --envs "$(methods_envs)" \
      ${METHODS_RESOURCES:+--resources "$METHODS_RESOURCES"} --outdir "$OUTDIR" --out "$pre.methods.$st.md" \
      ${METHODS_DOCX:+--docx} "$@" \
      || warn "Methods ($st) failed"
  done
  rm -f "$pf"
}

run_report() {   # sample
  local S=$1
  local RDIR=$OUTDIR/$S/Report L=$OUTDIR/$S/log_cluster_jobs/report
  mkdir -p "$RDIR" "$L"
  log "$S: writing report ..."
  local types=() d mfile=""
  for d in "$OUTDIR/$S"/Analysis/*/; do d=$(basename "$d"); [[ $d == sRNA* ]] && types+=("$d"); done
  ( activate_env "${ENV_REPORT:-${ENV_DEFAULT:-}}"
    write_methods sample "$RDIR/$S" --samples "$S" --types "$(IFS=,; echo "${types[*]}")" ) \
    > "$L/$S.methods.log" 2>&1 || warn "$S: Methods failed (see $L/$S.methods.log)"
  [[ -n "${METHODS_REPORT_STYLE:-}" && "$METHODS_REPORT_STYLE" != none && -s "$RDIR/$S.methods.$METHODS_REPORT_STYLE.md" ]] \
    && mfile=$RDIR/$S.methods.$METHODS_REPORT_STYLE.md
  set +e
  (
    set -euo pipefail
    activate_env "${ENV_REPORT:-${ENV_DEFAULT:-}}"
    python3 "$SCRIPTS_DIR/srna_report.py" --outdir "$OUTDIR" --sample "$S" \
      --pdf "$RDIR/$S.report.pdf" ${REPORT_HTML:+--html "$RDIR/$S.report.html"} \
      ${REPORT_MIRNA_MARKERS:+--mirna-markers "$REPORT_MIRNA_MARKERS"} \
      ${RFAM_FAMILY_TXT:+--rfam-family "$RFAM_FAMILY_TXT"} \
      ${RNA_QC:+--rna-qc "$RNA_QC"} \
      ${QC_RULES:+--qc-rules "$QC_RULES"} ${QC_EXCLUDE:+--qc-exclude "$QC_EXCLUDE"} \
      ${mfile:+--methods "$mfile"} --mqc-dir "$RDIR/multiqc" > "$L/$S.report.log" 2>&1
  )
  local rc=$?
  set -e
  if [[ $rc -eq 0 ]]; then
    log "$S: report written to $RDIR/$S.report.pdf${REPORT_HTML:+ and .html}${mfile:+ (with Methods, $METHODS_REPORT_STYLE)}"
    grep -h 'skipped' "$L/$S.report.log" | sed 's/^/    | /' || true
  else
    log "$S: report FAILED — see $L/$S.report.log"; tail -n 15 "$L/$S.report.log" | sed 's/^/    | /'
    RPT_FAILED=1
  fi
}

# NB: run_module must NOT be called inside if/||/&& — bash would then
# silently disable 'set -e' inside the module. The result is returned
# in MODULE_RC instead.
run_module() {
  local m=$1
  local stamp=$LOGDIR/.$m.done
  MODULE=$m
  MODULE_RC=0
  if [[ -e "$stamp" && $FORCE -eq 0 ]]; then
    log "already done ($stamp) — skipping (use -f to re-run)"
    return 0
  fi
  local envvar="ENV_${m^^}"
  local env="${!envvar:-${ENV_DEFAULT:-}}"
  log "start (env: ${env:-current})"
  local t0=$SECONDS
  local running=$LOGDIR/.$m.running failed=$LOGDIR/.$m.failed
  rm -f "$stamp" "$failed"
  # marker for srna_status.sh: host, pid, start, SLURM job (removed when the module ends;
  # a marker left behind means the job was killed)
  mkdir -p "$LOGDIR"
  printf '%s\t%s\t%s\t%s\n' "$(hostname -s 2>/dev/null || hostname)" "$$" "$(date '+%F %T')" \
    "${SLURM_JOB_ID:+$SLURM_JOB_ID${SLURM_ARRAY_TASK_ID:+_$SLURM_ARRAY_TASK_ID}}" > "$running"
  set +e
  (
    set -euo pipefail
    trap cleanup_module EXIT
    activate_env "$env"
    "module_$m"
  )
  MODULE_RC=$?
  set -e
  rm -f "$running"
  if [[ $MODULE_RC -eq 0 ]]; then
    printf '%s\t%s\n' "$(date '+%F %T')" "$(( SECONDS - t0 ))" > "$stamp"   # date, run time (s; srna_methods.py)
    log "finished in $(( SECONDS - t0 )) s"
  else
    log "FAILED (exit $MODULE_RC) after $(( SECONDS - t0 )) s"
    printf '%s\texit %s\n' "$(date '+%F %T')" "$MODULE_RC" > "$failed"
    # show the tail of the most recent tool log of this module
    local last
    last=$(ls -t "$LOGDIR/$m"/* 2>/dev/null | head -n 1 || true)
    if [[ -n "$last" && -s "$last" ]]; then
      log "last tool log: $last"
      tail -n 25 "$last" | sed 's/^/    | /'
    fi
  fi
  return 0
}

check_row() {
  local ok=0 f
  log "inputs:"
  for f in "raw_fastq:$IN_RAW" "assembly:$IN_ASSEMBLY" "gff:$IN_GFF" \
           "host_genome_fa:$HOST_GENOME_FA" "host_tx_fa:${HOST_TX_FA:-}"; do
    local k=${f%%:*} v=${f#*:}
    [[ -z "$v" ]] && { [[ "$k" == host_tx_fa ]] && log "  $k: (none)" || { log "  $k: NOT SET"; ok=1; }; continue; }
    if [[ -s "$v" ]]; then log "  $k: OK  $v"
    else log "  $k: MISSING $v"; ok=1; fi
  done
  for f in "contig2bin:$IN_BINS" "gtdbtk:$IN_GTDB" "metaT_features:$IN_MTF" "metaT_contigs:$IN_MTC" "metaG_features:$IN_MGF" "metaG_contigs:$IN_MGC"; do
    local k=${f%%:*} v=${f#*:} x
    if [[ -z "$v" ]]; then log "  $k: (not given)"; continue; fi
    for x in ${v//,/ }; do
      if [[ -s "$x" ]]; then log "  $k: OK  $x"; else log "  $k: MISSING $x"; ok=1; fi
    done
  done
  log "outputs: $ANALYSIS_DIR/{human,nonhuman,${ASSEMBLY_SUBDIR:-.},novel,fragments,bins}"
  local m
  for m in ${MODULES//,/ }; do
    [[ -e "$LOGDIR/.$m.done" ]] && log "  module $m: done" || log "  module $m: pending"
  done
  return $ok
}

N_ROWS=${#ROWS[@]}
SEL=()
if [[ -n "$ROW" ]]; then
  [[ "$ROW" =~ ^[0-9]+$ && "$ROW" -ge 1 && "$ROW" -le $N_ROWS ]] || { echo "Error: row $ROW out of range 1..$N_ROWS" >&2; exit 1; }
  SEL=("${ROWS[$((ROW-1))]}")
else
  SEL=("${ROWS[@]}")
fi

mkdir -p "$OUTDIR"
OUTDIR=$(readlink -f "$OUTDIR")
HOST_INDEX_TMPDIR=${HOST_INDEX_TMPDIR:-$OUTDIR/tmp/host_index}
echo "# srna_pipeline.sh  $(date '+%F %T')  host=$(hostname)  job=${SLURM_JOB_ID:-local}${SLURM_ARRAY_TASK_ID:+_$SLURM_ARRAY_TASK_ID}"
echo "# samples=$SAMPLESHEET outdir=$OUTDIR config=$CONFIG modules=$MODULES threads=$THREADS"

# activate the main environment (gut_srna_pipeline) once for the whole run
MODULE=main
if [[ -n "${ENV_DEFAULT:-}" ]]; then
  activate_env "$ENV_DEFAULT"
  echo "# environment: ${CONDA_PREFIX:-not activated}"
fi

# ------------------------------------------------------------------
# environment / tool / database check (for the selected modules)
# ------------------------------------------------------------------
ENVCHK_MISSING=0
_ok()   { printf '#   %-8s %-34s %s\n' "OK" "$1" "${2:-}"; }
_miss() { printf '#   %-8s %-34s %s\n' "MISSING" "$1" "${2:-}"; ENVCHK_MISSING=$((ENVCHK_MISSING+1)); }
chk_tool() { local t p; for t in "$@"; do p=$(command -v "$t" 2>/dev/null) && _ok "$t" "$p" || _miss "$t" "not in PATH"; done; }
chk_file() { [[ -s "$2" ]] && _ok "$1" "$2" || _miss "$1" "$2"; }
chk_dir()  { [[ -d "$2" ]] && _ok "$1" "$2" || _miss "$1" "$2"; }
chk_bt()   {
  if ! bowtie_index_ok "$2"; then _miss "$1" "$2.*.ebwt"; return 0; fi
  # an index built from RNA sequences (U instead of T, e.g. miRBase mature.fa) silently turns every
  # U into N in bowtie-build, and such references never match a read.
  # 1) FASTA next to the index (cheap, also finds a few U records in an otherwise DNA file)
  local fa="" ext nt nu
  for ext in .fa .fasta .fna .fa.gz .fasta.gz; do [[ -s "$2$ext" ]] && { fa=$2$ext; break; }; done
  if [[ -n "$fa" ]]; then
    read -r nt nu < <( (set +o pipefail; { [[ $fa == *.gz ]] && zcat "$fa" || cat "$fa"; } 2>/dev/null \
                        | awk '!/^>/{t += gsub(/[Tt]/, "&"); if (/[Uu]/) u++} END{print t+0, u+0}') )
    if [[ "${nt:-0}" -eq 0 ]]; then
      _miss "$1" "$fa: no T — RNA (U) sequences; convert (sed '/^>/!y/Uu/Tt/') and rebuild the index"
      return 0
    fi
    if [[ "${nu:-0}" -gt 0 ]]; then
      echo "#   WARN     $(printf '%-34s' "$1") $fa: $nu sequence line(s) with U — these references are N in the index; convert U->T and rebuild"
    fi
    _ok "$1" "$2"
    return 0
  fi
  # 2) no FASTA: look at the index itself; an empty bowtie-inspect output (killed, out of memory,
  #    old bowtie) is reported as a warning, not as a missing index
  if command -v bowtie-inspect >/dev/null 2>&1; then
    local out
    out=$( (set +o pipefail; timeout 120 bowtie-inspect "$2" 2>/dev/null | grep -v '^>' | head -c 200000) )
    if [[ -z "$out" ]]; then
      echo "#   WARN     $(printf '%-34s' "$1") $2: bowtie-inspect gave no output — T/U check skipped"
    elif [[ "$(tr -cd 'Tt' <<< "$out" | wc -c)" -eq 0 ]]; then
      _miss "$1" "$2: no T in the index — built from RNA (U) sequences? rebuild after U->T (sed '/^>/!y/Uu/Tt/')"
      return 0
    fi
  fi
  _ok "$1" "$2"
}
chk_py()   { local m; for m in "$@"; do python3 -c "import $m" 2>/dev/null && _ok "python: $m" || _miss "python: $m" "import failed"; done; }
chk_script() { local f; for f in "$@"; do chk_file "script $f" "$SCRIPTS_DIR/$f"; done; }
has_module() { [[ ",$MODULES," == *",$1,"* ]]; }

# --- per-module tool/database checks (used by check_environment) ---
chk_mod_qc() {
  chk_tool cutadapt fastqc
  return 0
}
chk_mod_split() {
  chk_tool bowtie bowtie-build samtools seqkit
  chk_file "HOST_GENOME_FA" "$HOST_GENOME_FA"
  [[ -n "${HOST_TX_FA:-}" ]] && chk_file "HOST_TX_FA" "$HOST_TX_FA"
  [[ -n "${HOST_GENOME_IDX:-}" ]] && chk_bt "HOST_GENOME_IDX" "$HOST_GENOME_IDX"
  [[ -n "${HOST_TX_IDX:-}" ]] && chk_bt "HOST_TX_IDX" "$HOST_TX_IDX"
  return 0
}
chk_mod_prep() {
  chk_tool gff2bed
  return 0
}
chk_mod_human() {
  chk_tool bowtie
  chk_bt IDX_HSA_MIRNA "$IDX_HSA_MIRNA"
  chk_bt IDX_HSA_TRNA "$IDX_HSA_TRNA"
  chk_bt IDX_HSA_RRNA "$IDX_HSA_RRNA"
  chk_bt IDX_HSA_NCRNA_OTHER "$IDX_HSA_NCRNA_OTHER"
  chk_script sam_fractional_counts.py summarize_human_subtree.py
  return 0
}
chk_mod_nonhuman() {
  chk_tool seqkit bowtie kraken2 bracken
  if [[ -x "$RIBO_DETECTOR_CPU" ]] || command -v "$RIBO_DETECTOR_CPU" >/dev/null 2>&1; then
    _ok ribodetector_cpu "$(command -v "$RIBO_DETECTOR_CPU")"
  else
    _miss ribodetector_cpu "$RIBO_DETECTOR_CPU"
  fi
  chk_bt IDX_TRNA_ALL "$IDX_TRNA_ALL"
  chk_bt IDX_MIRNA_EUK "$IDX_MIRNA_EUK"
  chk_bt IDX_RFAM_CLASSES "$IDX_RFAM_CLASSES"
  chk_file RFAMSEQ_TXT "$RFAMSEQ_TXT"
  [[ -n "${RFAM_FAMILY_TXT:-}" ]] && chk_file RFAM_FAMILY_TXT "$RFAM_FAMILY_TXT"
  chk_file "KRAKEN_DB hash.k2d" "$KRAKEN_DB/hash.k2d"
  chk_file "BRACKEN ${BRACKEN_READ_LENGTH}mers" "$BRACKEN_DB/database${BRACKEN_READ_LENGTH}mers.kmer_distrib"
  chk_script collapse_fastq_to_counted_fasta.py sum_collapsed_bowtie_hits.py \
             sum_collapsed_trna_by_taxon.py sum_collapsed_rfam_by_class_and_taxon.py
  return 0
}
chk_mod_assembly() {
  chk_tool bowtie bowtie-build samtools python3
  chk_script srna_bins.py
  return 0
}
chk_mod_novel() {
  chk_tool samtools bedtools seqkit
  [[ "$RUN_ARAGORN"  -eq 1 ]] && chk_tool aragorn
  [[ "$RUN_TRNASCAN" -eq 1 ]] && chk_tool tRNAscan-SE
  if [[ "$RUN_CMSCAN" -eq 1 ]]; then
    chk_tool cmscan
    chk_file RFAM_CM "$RFAM_CM"
    chk_file "RFAM_CM pressed (.i1m)" "$RFAM_CM.i1m"
  fi
  [[ "$RUN_RNAFOLD"  -eq 1 ]] && chk_tool RNAfold
  return 0
}

chk_mod_contam() {
  chk_tool bowtie bowtie-build samtools
  local fa
  for fa in ${CONTAM_FA:-}; do chk_file "CONTAM_FA" "$fa"; done
  [[ -n "${CONTAM_FA:-}" ]] || echo "#   (CONTAM_FA empty — only QC metrics, no contaminant screen)"
  return 0
}
chk_mod_fragments() {
  chk_tool samtools python3
  chk_script srna_fragments.py
  chk_script srna_integrity.py
  return 0
}
chk_mod_bins() {
  chk_tool samtools python3
  chk_script srna_bins.py
  return 0
}
chk_mod_crispr() {
  chk_tool samtools python3
  chk_script srna_crispr.py
  return 0
}
chk_mod_compare() {
  chk_tool python3
  chk_script srna_compare.py
  return 0
}
chk_mod_cohort() {
  chk_tool python3
  chk_py matplotlib
  chk_py numpy
  chk_script srna_cohort.py
  chk_script srna_report.py
  chk_script srna_qc.py
  if [[ -n "${QC_RULES:-}" ]]; then chk_file QC_RULES "$QC_RULES"; else chk_script qc_rules.tsv; fi
  [[ -n "${QC_EXCLUDE:-}" ]] && chk_file QC_EXCLUDE "$QC_EXCLUDE"
  if [[ -n "${MULTIQC:-}" ]]; then
    command -v "$MULTIQC" >/dev/null 2>&1 && _ok MultiQC "$(command -v "$MULTIQC")" \
      || echo "#   (MultiQC '$MULTIQC' not found — the cohort step skips it)"
  fi
  if [[ -n "${COHORT_METADATA:-}" ]]; then chk_file COHORT_METADATA "$COHORT_METADATA"
    case "$COHORT_METADATA" in *.xlsx|*.xlsm) chk_py openpyxl ;; esac
  fi
  return 0
}
chk_mod_report() {
  chk_tool python3
  chk_py matplotlib
  chk_script srna_report.py
  chk_script srna_qc.py
  if [[ -n "${REPORT_MIRNA_MARKERS:-}" ]]; then chk_file REPORT_MIRNA_MARKERS "$REPORT_MIRNA_MARKERS"
  else chk_script mirna_origin_markers.tsv; fi
  [[ -n "${RNA_QC:-}" ]] && chk_file RNA_QC "$RNA_QC"
  return 0
}

check_environment() {
  echo "# ---------------- environment / tool check ----------------"
  echo "# conda env:"
  if [[ -n "${ENV_DEFAULT:-}" ]]; then
    if [[ -n "${CONDA_PREFIX:-}" ]] && \
       [[ "$CONDA_PREFIX" == "${ENV_DEFAULT%/}" || "$(basename "$CONDA_PREFIX")" == "$(basename "${ENV_DEFAULT%/}")" ]]; then
      _ok "env $(basename "$ENV_DEFAULT")" "$CONDA_PREFIX"
    else
      _miss "env $(basename "$ENV_DEFAULT")" "active: ${CONDA_PREFIX:-none}"
    fi
  else
    echo "#   (ENV_DEFAULT empty — using current PATH)"
  fi

  echo "# general tools:"
  chk_tool python3 awk sed sort gzip zcat
  command -v pigz >/dev/null 2>&1 && _ok pigz "$(command -v pigz)" || echo "#   (pigz not found — gzip is used)"

  # per-module checks; modules with their own env (ENV_<MODULE>) are
  # checked inside that env, in a subshell
  local m envvar menv rc
  for m in ${MODULES//,/ }; do
    envvar="ENV_${m^^}"; menv="${!envvar:-${ENV_DEFAULT:-}}"
    if [[ -n "$menv" && "$menv" != "${ENV_DEFAULT:-}" ]]; then
      echo "# $m:  (own env: $menv)"
      set +e
      (
        ENVCHK_MISSING=0
        ACT_SOFT=1 activate_env "$menv" >/dev/null || true
        if [[ -n "${CONDA_PREFIX:-}" ]] && \
           [[ "$CONDA_PREFIX" == "${menv%/}" || "$(basename "$CONDA_PREFIX")" == "$(basename "${menv%/}")" ]]; then
          _ok "env $(basename "$menv")" "$CONDA_PREFIX"
          "chk_mod_$m"
        else
          _miss "env $(basename "$menv")" "could not activate $menv — does it exist? (micromamba env list)"
          echo "#            tools of module $m not checked"
        fi
        exit "$ENVCHK_MISSING"
      )
      rc=$?
      set -e
      ENVCHK_MISSING=$((ENVCHK_MISSING + rc))
    else
      echo "# $m:"
      "chk_mod_$m"
    fi
  done
  if [[ $ENVCHK_MISSING -eq 0 ]]; then
    echo "# environment/tool check: all OK"
  else
    echo "# environment/tool check: $ENVCHK_MISSING item(s) MISSING"
  fi
  echo "# ----------------------------------------------------------"
  return 0
}

if [[ $SKIP_ENV_CHECK -eq 0 || $ENV_CHECK_ONLY -eq 1 ]]; then
  check_environment
  if [[ $ENV_CHECK_ONLY -eq 1 ]]; then
    [[ $ENVCHK_MISSING -eq 0 ]] && exit 0 || exit 1
  fi
  if [[ $ENVCHK_MISSING -gt 0 && $CHECK_ONLY -eq 0 ]]; then
    echo "Error: environment/tool check failed — fix the MISSING items above (or use --skip-env-check)" >&2
    exit 1
  fi
fi

# ------------------------------------------------------------------
# COHORT (--cohort): cross-sample summary with the metadata
# ------------------------------------------------------------------
run_cohort() {
  local CDIR=${COHORT_DIR:-$OUTDIR/Cohort} L=$OUTDIR/log_cohort
  mkdir -p "$CDIR" "$L"
  PREFIX=cohort; TYPE=all; MODULE=cohort
  [[ -n "${COHORT_METADATA:-}" ]] || warn "COHORT_METADATA not set — cohort summary without groups"
  log "cross-sample summary -> $CDIR ..."
  local mstyle=""
  if [[ -n "${METHODS_REPORT_STYLE:-}" && "$METHODS_REPORT_STYLE" != none ]]; then
    mstyle=$METHODS_REPORT_STYLE; write_params "$CDIR/cohort.params.tsv"
  fi
  set +e
  (
    set -euo pipefail
    activate_env "${ENV_REPORT:-${ENV_DEFAULT:-}}"
    python3 "$SCRIPTS_DIR/srna_cohort.py" --outdir "$OUTDIR" --samples "$SAMPLESHEET" --cohort-dir "$CDIR" \
      ${COHORT_METADATA:+--metadata "$COHORT_METADATA"} \
      ${COHORT_ID_COL:+--id-col "$COHORT_ID_COL"} \
      --group-col "${COHORT_GROUP_COL-group}" \
      ${COHORT_GROUP_ORDER:+--group-order "$COHORT_GROUP_ORDER"} \
      ${COHORT_COVARIATES:+--covariates "$COHORT_COVARIATES"} \
      ${RNA_QC:+--rna-qc "$RNA_QC"} \
      ${QC_RULES:+--qc-rules "$QC_RULES"} ${QC_EXCLUDE:+--qc-exclude "$QC_EXCLUDE"} \
      --qc-drop "${COHORT_QC_DROP:-FAIL}" --outlier-z "${COHORT_OUTLIER_Z:-4}" \
      --min-mirna-reads "${COHORT_MIN_MIRNA_READS:-50}" --sig-q "${COHORT_SIG_Q:-0.1}" \
      ${COHORT_TERM_NAMES:+--term-names "$COHORT_TERM_NAMES"} \
      --nmds "${COHORT_NMDS-genus,fn_kegg_pathway,mirna,rfam_family}" --nmds-starts "${COHORT_NMDS_STARTS:-10}" \
      ${COHORT_MODEL:+--model "$COHORT_MODEL"} \
      ${COHORT_BLOCK:+--block "$COHORT_BLOCK"} ${COHORT_REFERENCE:+--reference "$COHORT_REFERENCE"} \
      ${REPORT_MIRNA_MARKERS:+--mirna-markers "$REPORT_MIRNA_MARKERS"} \
      --permutations "${COHORT_PERMUTATIONS:-999}" ${REPORT_HTML:+--html} \
      ${mstyle:+--methods-style "$mstyle" --methods-params "$CDIR/cohort.params.tsv" --methods-envs "$(methods_envs)"} \
      ${mstyle:+${METHODS_RESOURCES:+--methods-resources "$METHODS_RESOURCES"}} > "$L/cohort.log" 2>&1
  )
  local rc=$?
  set -e
  if [[ $rc -ne 0 ]]; then warn "cohort summary failed (see $L/cohort.log)"; return 1; fi
  log "cohort report: $CDIR/cohort.report.pdf"
  rm -f "$CDIR/cohort.params.tsv"
  ( activate_env "${ENV_REPORT:-${ENV_DEFAULT:-}}"
    write_methods all "$CDIR/cohort" --cohort-dir "$CDIR" ) >> "$L/cohort.log" 2>&1 \
    && log "cohort Methods: $CDIR/cohort.methods.{${METHODS_STYLE:-nature,star}}.md"
  # MultiQC over all FastQC / cutadapt outputs + the per-sample report tables (*_mqc.tsv)
  if [[ -n "${MULTIQC:-}" ]] && command -v "$MULTIQC" >/dev/null 2>&1; then
    log "MultiQC -> $CDIR/multiqc ..."
    local dirs=() S
    for S in "${CMP_ORDER[@]}"; do
      [[ -d "$OUTDIR/$S/Stats" ]] && dirs+=("$OUTDIR/$S/Stats")
      [[ -d "$OUTDIR/$S/Report/multiqc" ]] && dirs+=("$OUTDIR/$S/Report/multiqc")
    done
    # --dirs: sample names carry <sample>/Stats/<type>, otherwise Int and Ext of a sample collide
    "$MULTIQC" -f -q --dirs --dirs-depth 3 -o "$CDIR/multiqc" -n multiqc_report.html "${dirs[@]}" \
      > "$L/multiqc.log" 2>&1 && log "MultiQC report: $CDIR/multiqc/multiqc_report.html" \
      || warn "MultiQC failed (see $L/multiqc.log)"
  elif [[ -n "${MULTIQC:-}" ]]; then
    log "MultiQC ($MULTIQC) not found or not executable — skipped (envs/multiqc.yml, MULTIQC=<path>)"
  else
    log "MULTIQC is empty — MultiQC skipped"
  fi
}

# per-sample optional columns (6-11) may be given on any row of a sample
declare -A SAMPLE_OPT=()
declare -A CMP_TYPES=() CMP_BIN=() CMP_GTDB=() CMP_MTF=() CMP_MTC=() CMP_MGF=() CMP_MGC=()
CMP_ORDER=()
for line in "${ROWS[@]}"; do
  read -r _s _t _r _a _g _b _x _f _c _m _n _rest <<< "$line" || true
  [[ -n "${_s:-}" ]] || continue
  for pair in "C_BIN:${_b:-}" "C_GTDB:${_x:-}" "C_MTF:${_f:-}" "C_MTC:${_c:-}" "C_MGF:${_m:-}" "C_MGC:${_n:-}"; do
    k=${pair%%:*}; v=${pair#*:}
    [[ -n "$v" && "$v" != "-" && -z "${SAMPLE_OPT["$_s:$k"]:-}" ]] && SAMPLE_OPT["$_s:$k"]=$v
  done
  [[ -n "${CMP_TYPES[$_s]+x}" ]] || CMP_ORDER+=("$_s")
  [[ " ${CMP_TYPES[$_s]:-} " == *" ${_t:-} "* ]] || CMP_TYPES[$_s]="${CMP_TYPES[$_s]:-}${CMP_TYPES[$_s]:+ }${_t:-}"
done
unset _s _t _r _a _g _b _x _f _c _m _n _rest pair k v

if [[ $COMPARE -eq 1 || $REPORT -eq 1 ]]; then
  CMP_FAILED=0; RPT_FAILED=0
  for S in "${CMP_ORDER[@]}"; do
    [[ -n "$SAMPLE_FILTER" && "$S" != "$SAMPLE_FILTER" ]] && continue
    PREFIX=$S; TYPE=all
    if [[ $COMPARE -eq 1 ]]; then
      MODULE=compare
      CMP_BIN[$S]=${SAMPLE_OPT["$S:C_BIN"]:-}; CMP_GTDB[$S]=${SAMPLE_OPT["$S:C_GTDB"]:-}
      CMP_MTF[$S]=${SAMPLE_OPT["$S:C_MTF"]:-}; CMP_MTC[$S]=${SAMPLE_OPT["$S:C_MTC"]:-}
      CMP_MGF[$S]=${SAMPLE_OPT["$S:C_MGF"]:-}; CMP_MGC[$S]=${SAMPLE_OPT["$S:C_MGC"]:-}
      run_compare "$S"
    fi
    if [[ $REPORT -eq 1 ]]; then
      MODULE=report
      run_report "$S"
    fi
  done
fi
if [[ $COHORT -eq 1 ]]; then
  COH_FAILED=0
  run_cohort || COH_FAILED=1
fi
if [[ $COMPARE -eq 1 || $REPORT -eq 1 || $COHORT -eq 1 ]]; then
  exit $(( ${CMP_FAILED:-0} | ${RPT_FAILED:-0} | ${COH_FAILED:-0} ))
fi

FAILED=0
[[ $CHECK_ONLY -eq 1 && $ENVCHK_MISSING -gt 0 ]] && FAILED=1
for line in "${SEL[@]}"; do
  read -r PREFIX TYPE C_RAW C_ASM C_GFF C_BIN C_GTDB C_MTF C_MTC C_MGF C_MGC _rest <<< "$line" || true
  PREFIX=${PREFIX:-}; TYPE=${TYPE:-}
  [[ "$PREFIX" =~ ^[A-Za-z0-9._-]+$ ]] || { echo "Error: bad sample id in row: $line" >&2; FAILED=1; continue; }
  [[ "$TYPE"   =~ ^[A-Za-z0-9._-]+$ ]] || { echo "Error: bad/missing type in row: $line" >&2; FAILED=1; continue; }
  [[ -n "$SAMPLE_FILTER" && "$PREFIX" != "$SAMPLE_FILTER" ]] && continue
  [[ -n "$TYPE_FILTER"   && "$TYPE"   != "$TYPE_FILTER"   ]] && continue

  # raw_fastq, assembly and gff: required, full (absolute) paths
  bad=""
  [[ "${C_RAW:-}" == /* ]] || bad+=" raw_fastq(col3)"
  [[ "${C_ASM:-}" == /* ]] || bad+=" assembly(col4)"
  [[ "${C_GFF:-}" == /* ]] || bad+=" gff(col5)"
  if [[ -n "$bad" ]]; then
    echo "Error: $PREFIX/$TYPE: missing or not a full path:$bad" >&2; FAILED=1; continue
  fi
  IN_RAW=$C_RAW; IN_ASSEMBLY=$C_ASM; IN_GFF=$C_GFF
  # optional per-sample columns (full path or "-"); fall back to another row of the same sample
  for v in C_BIN C_GTDB C_MTF C_MTC C_MGF C_MGC; do
    x=${!v:-}; [[ "$x" == "-" ]] && x=""
    [[ -z "$x" ]] && x=${SAMPLE_OPT["$PREFIX:$v"]:-}
    if [[ -n "$x" && "$x" != /* ]]; then
      echo "Error: $PREFIX/$TYPE: $v must be a full path or '-': $x" >&2; FAILED=1; x=""
    fi
    printf -v "$v" '%s' "$x"
  done
  IN_BINS=$C_BIN; IN_GTDB=$C_GTDB; IN_MTF=$C_MTF; IN_MTC=$C_MTC; IN_MGF=$C_MGF; IN_MGC=$C_MGC

  set_paths
  MODULE=main
  log "=============== sample $PREFIX / type $TYPE ==============="

  if [[ $CHECK_ONLY -eq 1 ]]; then
    check_row || FAILED=1
    continue
  fi

  make_dirs
  {
    echo -e "sample\t$PREFIX"; echo -e "type\t$TYPE"
    echo -e "raw_fastq\t$IN_RAW"
    echo -e "assembly\t$IN_ASSEMBLY"; echo -e "gff\t$IN_GFF"
    echo -e "host_genome_fa\t$HOST_GENOME_FA"; echo -e "host_tx_fa\t${HOST_TX_FA:-}"
    echo -e "contig2bin\t$IN_BINS"; echo -e "gtdbtk\t$IN_GTDB"
    echo -e "metaT_features\t$IN_MTF"; echo -e "metaT_contigs\t$IN_MTC"; echo -e "metaG_features\t$IN_MGF"; echo -e "metaG_contigs\t$IN_MGC"
    echo -e "config\t$CONFIG"; echo -e "date\t$(date '+%F %T')"
  } > "$LOGDIR/inputs.tsv"

  for m in ${MODULES//,/ }; do
    run_module "$m"
    if [[ $MODULE_RC -ne 0 ]]; then
      MODULE=$m; log "remaining modules for $PREFIX/$TYPE not run"
      FAILED=1
      break
    fi
  done
  MODULE=main
  log "=============== done $PREFIX / $TYPE ==============="
done

exit $FAILED
