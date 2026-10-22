#!/usr/bin/env bash
# srna_status.sh — check which rows of a run are finished, running or failed, and resubmit the failed ones
#
# Usage: srna_status.sh -s samples.tsv -o outdir [options] [-- extra options for submit_srna_pipeline.sh]
#   -m LIST       modules expected per row          [qc,contam,prep,split,human,nonhuman,assembly,novel,fragments,bins,crispr]
#   -J NAME       SLURM job name of the run (to detect running/pending tasks; as -J of submit)  [sRNA]
#   -v            list every row (default: only rows that are not complete)
#   --resubmit    resubmit the incomplete, not running rows (submit_srna_pipeline.sh --rows ... --no-compare)
#   --dry-run     with --resubmit: only print the submit command
# A row is: done (all .done stamps); running (a module's .running marker: interactive or batch;
# for markers from another host the SLURM job in the marker is checked when squeue exists);
# pending (in squeue); failed (a module wrote .failed, or a .running marker of a dead process/job
# = killed); not started (nothing yet). The last ERROR line of the row's newest log is shown.
set -euo pipefail
DIR=$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")
usage() { sed -n '2,13p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }

SAMPLES="" OUTDIR="" MODULES="qc,contam,prep,split,human,nonhuman,assembly,novel,fragments,bins,crispr"
JOB="sRNA" VERBOSE=0 RESUB=0 DRY=0 EXTRA=()
while [[ $# -gt 0 ]]; do
  case $1 in
    -s) SAMPLES=$2; shift 2 ;;
    -o) OUTDIR=$2; shift 2 ;;
    -m) MODULES=$2; shift 2 ;;
    -J) JOB=$2; shift 2 ;;
    -v) VERBOSE=1; shift ;;
    --resubmit) RESUB=1; shift ;;
    --dry-run) DRY=1; shift ;;
    -h|--help) usage 0 ;;
    --) shift; EXTRA=("$@"); break ;;
    *) echo "Error: unknown option $1" >&2; usage 1 ;;
  esac
done
[[ -s "$SAMPLES" && -n "$OUTDIR" ]] || usage 1

# data rows exactly as srna_pipeline.sh counts them (row number = array task id)
mapfile -t ROWS < <(awk '{sub(/\r$/,"")} /^[[:space:]]*$/ || /^[[:space:]]*#/ {next}
                         tolower($1) ~ /^(sample_id|sample|sampleid)$/ {next} {print $1"\t"$2}' "$SAMPLES")

# running / pending array tasks. A task belongs to this run if its array job name is $JOB or its
# array job id appears in the SLURM log names of this outdir (<name>_<arrayjobid>_<task>.out), so a
# different -J or a resubmission is still recognised. Pending ranges (e.g. 11-179%20) are expanded.
declare -A QUEUED=() RUNIDS=()
for f in "$OUTDIR"/slurm_logs/*_*_*.out "$OUTDIR"/slurm_logs/*_*_*.err; do
  [[ -e $f ]] || continue
  b=$(basename "$f"); b=${b%.*}; a=${b%_*}; a=${a##*_}
  [[ $a =~ ^[0-9]+$ ]] && RUNIDS[$a]=1
done
NQ=0
if command -v squeue >/dev/null 2>&1; then
  while read -r arr tid name st; do
    [[ $name == "$JOB" || -n "${RUNIDS[$arr]:-}" ]] || continue
    tid=${tid%%\%*}; tid=${tid#[}; tid=${tid%]}
    for part in ${tid//,/ }; do
      if [[ $part =~ ^([0-9]+)-([0-9]+)$ ]]; then
        for ((k = BASH_REMATCH[1]; k <= BASH_REMATCH[2]; k++)); do QUEUED[$k]=$st; NQ=$((NQ + 1)); done
      elif [[ $part =~ ^[0-9]+$ ]]; then
        QUEUED[$part]=$st; NQ=$((NQ + 1))
      fi
    done
  done < <(squeue -h -u "${USER:-$(id -un)}" -r -o "%F %K %j %T" 2>/dev/null || true)
  echo "# squeue: $NQ array task(s) of this run queued or running"
else
  echo "# squeue not available — running state from the module markers only"
fi

compress() {   # 1 2 3 5 -> 1-3,5
  awk '{for(i=1;i<=NF;i++){v=$i; if(s==""){s=v;p=v} else if(v==p+1){p=v} else {o=o (o?",":"") (s==p?s:s"-"p); s=v;p=v}}}
       END{if(s!="") o=o (o?",":"") (s==p?s:s"-"p); print o}' <<< "$*"
}

n_done=0 n_run=0 n_fail=0 n_new=0 FAILED=() NEW=()
HOSTNAME_S=$(hostname -s 2>/dev/null || hostname)
HAVE_SQ=0; command -v squeue >/dev/null 2>&1 && HAVE_SQ=1
NMOD=$(tr ',' '\n' <<< "$MODULES" | grep -c .)
printf "%-5s %-12s %-10s %-9s %s\n" row sample type status "missing modules / last error"
for i in "${!ROWS[@]}"; do
  r=$((i + 1)); IFS=$'\t' read -r s t <<< "${ROWS[$i]}"
  miss=()
  for m in ${MODULES//,/ }; do
    [[ -e "$OUTDIR/$s/log_cluster_jobs/$t/.$m.done" ]] || miss+=("$m")
  done
  if [[ ${#miss[@]} -eq 0 ]]; then
    n_done=$((n_done + 1)); [[ $VERBOSE -eq 1 ]] && printf "%-5s %-12s %-10s %-9s\n" "$r" "$s" "$t" done
    continue
  fi
  # module markers written by srna_pipeline.sh
  LD="$OUTDIR/$s/log_cluster_jobs/$t"
  run_m="" fail_m="" killed=""
  for m in "${miss[@]}"; do
    if [[ -e "$LD/.$m.running" ]]; then
      IFS=$'\t' read -r mh mp ms mj < "$LD/.$m.running" || true
      alive=unknown
      if [[ $mh == "$HOSTNAME_S" ]]; then
        kill -0 "$mp" 2>/dev/null && alive=yes || alive=no
      elif [[ -n "${mj:-}" && $HAVE_SQ -eq 1 ]]; then
        squeue -h -j "${mj%%_*}" >/dev/null 2>&1 && [[ -n "$(squeue -h -j "${mj%%_*}" 2>/dev/null)" ]] && alive=yes || alive=no
      fi
      if [[ $alive == no ]]; then
        killed="$m (started $ms on $mh)"
      else
        run_m="$m since $ms on $mh${mj:+, job $mj}"
        [[ $alive == unknown ]] && run_m+=" (other host, not verifiable)"
      fi
      break
    fi
    [[ -e "$LD/.$m.failed" ]] && { fail_m=$m; break; }
  done
  if [[ -n "$run_m" ]]; then
    n_run=$((n_run + 1))
    printf "%-5s %-12s %-10s %-9s %s\n" "$r" "$s" "$t" running "$run_m"
    continue
  fi
  if [[ -z "$fail_m$killed" && -n "${QUEUED[$r]:-}" ]]; then
    n_run=$((n_run + 1))
    printf "%-5s %-12s %-10s %-9s %s\n" "$r" "$s" "$t" "${QUEUED[$r],,}" "${miss[*]}"
    continue
  fi
  if [[ -z "$fail_m$killed" && ${#miss[@]} -eq $NMOD ]]; then
    n_new=$((n_new + 1)); NEW+=("$r")
    [[ $VERBOSE -eq 1 ]] && printf "%-5s %-12s %-10s %-9s\n" "$r" "$s" "$t" "not started"
    continue
  fi
  n_fail=$((n_fail + 1)); FAILED+=("$r")
  err=""
  log=$(ls -t "$OUTDIR"/slurm_logs/"${JOB}"_*_"$r".{out,err} 2>/dev/null | head -n 1 || true)
  if [[ -n "$log" ]]; then
    base=${log%.*}
    err=$(cat "$base.out" "$base.err" 2>/dev/null | grep -E "ERROR|FAILED|Killed|CANCELLED|TIME LIMIT|oom-kill|Out Of Memory" \
          | tail -n 1 | cut -c1-150 || true)
  fi
  if [[ -z "$err" ]]; then                       # interactive runs: newest tool log of the module
    lm=${fail_m:-${killed%% *}}
    tl=$(ls -t "$LD/${lm:-none}"/* 2>/dev/null | head -n 1 || true)
    [[ -n "$tl" ]] && err=$(grep -E "ERROR|Error|error:|Killed|MemoryError|Traceback" "$tl" 2>/dev/null | tail -n 1 | cut -c1-150 || true)
  fi
  why=${killed:+killed in $killed}; why=${why:-${fail_m:+failed in $fail_m}}
  printf "%-5s %-12s %-10s %-9s %s\n" "$r" "$s" "$t" failed "${why:+$why; }missing: ${miss[*]}${err:+ | $err}"
done
echo "# rows: ${#ROWS[@]}  done: $n_done  running/pending: $n_run  failed: $n_fail  not started: $n_new"
[[ $n_new -gt 0 ]] && echo "# not started: $(compress "${NEW[@]}")"

# per-sample steps after the array
nrep=$(compgen -G "$OUTDIR/*/Report/*.report.pdf" | wc -l || true)
ncmp=$(compgen -G "$OUTDIR/*/Analysis/comparison/*.compare.summary.tsv" | wc -l || true)
echo "# comparisons: $ncmp  reports: $nrep  cohort report: $([[ -s $OUTDIR/Cohort/cohort.report.pdf ]] && echo yes || echo no)"

if [[ ${#FAILED[@]} -gt 0 ]]; then
  LIST=$(compress "${FAILED[@]}")
  echo "# rows to rerun: $LIST"
  CMD=("$DIR/submit_srna_pipeline.sh" -s "$SAMPLES" -o "$OUTDIR" -m "$MODULES" -J "$JOB" --rows "$LIST" --no-compare
       ${EXTRA[@]+"${EXTRA[@]}"})
  if [[ $RESUB -eq 1 ]]; then
    [[ $DRY -eq 1 ]] && CMD+=(--dry-run)
    echo "# ${CMD[*]}"
    "${CMD[@]}"
  else
    echo "# resubmit with:  $0 -s $SAMPLES -o $OUTDIR -J $JOB --resubmit"
  fi
  echo "# afterwards:     $DIR/submit_srna_pipeline.sh -s $SAMPLES -o $OUTDIR --compare-only   (compare + report + cohort)"
fi
