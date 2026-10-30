#!/usr/bin/env bash
# =====================================================================
# submit_srna_pipeline.sh — submit srna_pipeline.sh as a SLURM job array
# (one array task per data row of the sample sheet = one sample/type),
# followed by one per-sample comparison job (Ext / Int / metaT) that
# starts when the whole array has finished.
# =====================================================================
set -euo pipefail

SCRIPT_DIR=$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")
RUN=$SCRIPT_DIR/srna_pipeline.sh

usage() {
cat <<EOF
Usage: $(basename "$0") -s <samples.tsv> -o <outdir> [options] [-- extra sbatch args]

Pipeline options (passed to srna_pipeline.sh):
  -s, --samples FILE     sample sheet                                  (required)
  -o, --outdir DIR       base output directory                         (required)
  -c, --config FILE      config file          [$SCRIPT_DIR/srna_pipeline.config.sh]
  -m, --modules LIST     modules              [qc,contam,prep,split,human,nonhuman,mirna_origin,assembly,novel,fragments,bins,crispr]
  -f, --force            re-run finished modules
  -e, --env ENV          conda/micromamba env activated in the job   [ENV_DEFAULT from config = gut_srna_pipeline]
      --no-compare       do not submit the comparison + report job after the array
      --compare-only     submit only the comparison + report job (no array)
      --compare-time T   wall time of the comparison + report job (all samples, sequential)  [1-00:00:00]
      --compare-mem M    memory of the comparison + report + cohort job                     [32G]

SLURM options:
  -n, --cpus N           cpus per task                                 [16]
  -T, --time T           wall time                                     [2-00:00:00]
      --mem MEM          memory per node (e.g. 64G); default: partition default
  -p, --partition P      partition                                     [batch]
  -q, --qos Q            QoS                                           [normal]
  -A, --account A        account
  -J, --job-name NAME    job name                                      [sRNA]
      --rows LIST        array rows to run, e.g. 1-10,15               [all]
      --max-parallel N   max simultaneously running tasks              [20]
      --mail-user ADDR   send mail on END/FAIL
      --dry-run          print the job script and sbatch command only
  -h, --help
EOF
}

SAMPLES=""; OUTDIR=""; CONFIG="$SCRIPT_DIR/srna_pipeline.config.sh"
MODULES="qc,contam,prep,split,human,nonhuman,mirna_origin,assembly,novel,fragments,bins,crispr"; FORCE=""
DO_COMPARE=1; COMPARE_ONLY=0
CPUS=16; TIME="2-00:00:00"; MEM=""; PARTITION="batch"; QOS="normal"; ACCOUNT=""
ENVNAME=""; JOBNAME="sRNA"; ROWS=""; MAXPAR=20; MAIL=""; DRYRUN=0
EXTRA=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    -s|--samples)    SAMPLES=$2; shift 2 ;;
    -o|--outdir)     OUTDIR=$2; shift 2 ;;
    -c|--config)     CONFIG=$2; shift 2 ;;
    -m|--modules)    MODULES=$2; shift 2 ;;
    -f|--force)      FORCE="--force"; shift ;;
    -e|--env)        ENVNAME=$2; shift 2 ;;
    --no-compare)    DO_COMPARE=0; shift ;;
    --compare-only)  COMPARE_ONLY=1; shift ;;
    --compare-time)  COMPARE_TIME=$2; shift 2 ;;
    --compare-mem)   COMPARE_MEM=$2; shift 2 ;;
    -n|--cpus)       CPUS=$2; shift 2 ;;
    -T|--time)       TIME=$2; shift 2 ;;
    --mem)           MEM=$2; shift 2 ;;
    -p|--partition)  PARTITION=$2; shift 2 ;;
    -q|--qos)        QOS=$2; shift 2 ;;
    -A|--account)    ACCOUNT=$2; shift 2 ;;
    -J|--job-name)   JOBNAME=$2; shift 2 ;;
    --rows)          ROWS=$2; shift 2 ;;
    --max-parallel)  MAXPAR=$2; shift 2 ;;
    --mail-user)     MAIL=$2; shift 2 ;;
    --dry-run)       DRYRUN=1; shift ;;
    -h|--help)       usage; exit 0 ;;
    --)              shift; EXTRA=("$@"); break ;;
    *) echo "Error: unknown option '$1'" >&2; usage >&2; exit 1 ;;
  esac
done

[[ -n "$SAMPLES" && -n "$OUTDIR" ]] || { usage >&2; exit 1; }
[[ -x "$RUN" ]] || chmod +x "$RUN" 2>/dev/null || true
SAMPLES=$(readlink -f "$SAMPLES")
CONFIG=$(readlink -f "$CONFIG")
mkdir -p "$OUTDIR/slurm_logs"
OUTDIR=$(readlink -f "$OUTDIR")

# environment to activate in the job (default: ENV_DEFAULT from the config)
if [[ -z "$ENVNAME" ]]; then
  ENVNAME=$(bash -c 'source "$1" >/dev/null 2>&1; echo "${ENV_DEFAULT:-}"' _ "$CONFIG")
fi
[[ -n "$ENVNAME" ]] || ENVNAME=gut_srna_pipeline

N=$(bash "$RUN" -s "$SAMPLES" -c "$CONFIG" --count-rows)
[[ "$N" -gt 0 ]] || { echo "Error: no data rows in $SAMPLES" >&2; exit 1; }
ARRAY="${ROWS:-1-$N}%$MAXPAR"

# quick validation of all inputs before anything is queued
echo "# validating $N row(s) ..."
CHECKLOG="$OUTDIR/slurm_logs/check.$(date +%Y%m%d_%H%M%S).log"
if ! bash "$RUN" -s "$SAMPLES" -o "$OUTDIR" -c "$CONFIG" -m "$MODULES" --check > "$CHECKLOG" 2>&1; then
  echo "WARNING: environment, tools, databases or inputs missing — see $CHECKLOG" >&2
  grep -E "MISSING|Error" "$CHECKLOG" | head -20 >&2 || true
  if [[ $DRYRUN -eq 0 && -t 0 ]]; then
    read -r -p "Submit anyway? [y/N] " ans
    [[ "$ans" == [yY]* ]] || { echo "not submitted"; exit 1; }
  fi
fi

SB=(--job-name "$JOBNAME" --array "$ARRAY" --nodes 1 --ntasks 1 --cpus-per-task "$CPUS"
    --time "$TIME" --partition "$PARTITION" --qos "$QOS"
    --output "$OUTDIR/slurm_logs/%x_%A_%a.out" --error "$OUTDIR/slurm_logs/%x_%A_%a.err")
[[ -n "$MEM" ]]     && SB+=(--mem "$MEM")
[[ -n "$ACCOUNT" ]] && SB+=(--account "$ACCOUNT")
[[ -n "$MAIL" ]]    && SB+=(--mail-user "$MAIL" --mail-type END,FAIL)
SB+=(${EXTRA[@]+"${EXTRA[@]}"})

JOBSCRIPT=$(cat <<EOF
#!/bin/bash -l
# generated by submit_srna_pipeline.sh on $(date '+%F %T')
set -euo pipefail
echo "# job \$SLURM_JOB_ID array task \$SLURM_ARRAY_TASK_ID on \$(hostname)"

# the environment ($ENVNAME) is activated automatically by srna_pipeline.sh (-e)
# node-local temp space if the cluster provides one
export TMP_ROOT="\${TMP_ROOT:-\${SLURM_TMPDIR:-\${TMPDIR:-/tmp}}}"
# one line on purpose: backslash-newline inside a here-document in \$( ) is handled differently by bash versions
bash "$RUN" -s "$SAMPLES" -o "$OUTDIR" -c "$CONFIG" -e "$ENVNAME" -m "$MODULES" -t "\$SLURM_CPUS_PER_TASK" -r "\$SLURM_ARRAY_TASK_ID" $FORCE
EOF
)

CMPSCRIPT=$(cat <<EOF
#!/bin/bash -l
# generated by submit_srna_pipeline.sh on $(date '+%F %T') — per-sample comparison
set -euo pipefail
echo "# job \$SLURM_JOB_ID (compare) on \$(hostname)"
bash "$RUN" -s "$SAMPLES" -o "$OUTDIR" -c "$CONFIG" -e "$ENVNAME" --compare --report --cohort
EOF
)
CB=(--job-name "${JOBNAME}_compare" --nodes 1 --ntasks 1 --cpus-per-task 4 --mem "${COMPARE_MEM:-32G}"
    --time "${COMPARE_TIME:-1-00:00:00}"
    --partition "$PARTITION" --qos "$QOS"
    --output "$OUTDIR/slurm_logs/%x_%j.out" --error "$OUTDIR/slurm_logs/%x_%j.err")
[[ -n "$ACCOUNT" ]] && CB+=(--account "$ACCOUNT")
[[ -n "$MAIL" ]]    && CB+=(--mail-user "$MAIL" --mail-type END,FAIL)

echo "# env: $ENVNAME"
echo "# rows: $N  array: $ARRAY  cpus: $CPUS  time: $TIME  partition: $PARTITION  qos: $QOS"
if [[ $DRYRUN -eq 1 ]]; then
  if [[ $COMPARE_ONLY -eq 0 ]]; then
    echo "sbatch ${SB[*]}"
    echo "----- job script -----"
    echo "$JOBSCRIPT"
  fi
  if [[ $DO_COMPARE -eq 1 || $COMPARE_ONLY -eq 1 ]]; then
    dep_txt=""; [[ $COMPARE_ONLY -eq 0 ]] && dep_txt=" --dependency afterany:<array job id>"
    echo "sbatch ${CB[*]}$dep_txt"
    echo "----- compare job script -----"
    echo "$CMPSCRIPT"
  fi
  exit 0
fi

DEP=()
if [[ $COMPARE_ONLY -eq 0 ]]; then
  echo "$JOBSCRIPT" > "$OUTDIR/slurm_logs/${JOBNAME}.jobscript.sh"
  JID=$(sbatch --parsable "${SB[@]}" <<< "$JOBSCRIPT")
  JID=${JID%%;*}
  echo "Submitted array job $JID"
  DEP=(--dependency "afterany:$JID")
fi
if [[ $DO_COMPARE -eq 1 || $COMPARE_ONLY -eq 1 ]]; then
  echo "$CMPSCRIPT" > "$OUTDIR/slurm_logs/${JOBNAME}_compare.jobscript.sh"
  CID=$(sbatch --parsable "${CB[@]}" ${DEP[@]+"${DEP[@]}"} <<< "$CMPSCRIPT")
  echo "Submitted compare job ${CID%%;*}${DEP[*]:+ (after array $JID)}"
fi
