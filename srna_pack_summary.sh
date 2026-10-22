#!/usr/bin/env bash
# =====================================================================
# srna_pack_summary.sh — pack the small summary files of srna_pipeline.sh
# runs (all modules, comparison, reports, module status) into one tarball,
# e.g. to check results elsewhere. Large per-read files (BAM, FASTQ, SAM)
# are never included; files above --max-size are listed, not packed.
#
# Optional: --clean-legacy removes leftovers of the old srnaMapper/CMC
# assembly mapping and the old novel BAM (dry run unless --yes).
# =====================================================================
set -euo pipefail

usage() {
cat <<USAGE
Usage: $(basename "$0") -o <pipeline outdir> [options] [sample ...]

  -o, --outdir DIR     output directory of srna_pipeline.sh (-o)          (required)
  -t, --tar FILE       tarball to write              [srna_summary_<samples>_<date>.tgz in outdir]
      --max-size SIZE  largest file to pack (find -size syntax)              [5M]
      --no-reports     leave out Report/*.pdf / *.html
      --clean-legacy   list old srnaMapper/CMC files and novel/*.sorted.bam (not packed)
      --yes            with --clean-legacy: really delete them
  -h, --help
  sample ...           sample directories to include             [all with Analysis/]
USAGE
}

OUT=""; TAR=""; MAXS=5M; REPORTS=1; CLEAN=0; YES=0; SAMPLES=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    -o|--outdir) OUT=$2; shift 2 ;;
    -t|--tar) TAR=$2; shift 2 ;;
    --max-size) MAXS=$2; shift 2 ;;
    --no-reports) REPORTS=0; shift ;;
    --clean-legacy) CLEAN=1; shift ;;
    --yes) YES=1; shift ;;
    -h|--help) usage; exit 0 ;;
    -*) echo "unknown option $1" >&2; usage >&2; exit 1 ;;
    *) SAMPLES+=("$1"); shift ;;
  esac
done
[[ -n "$OUT" && -d "$OUT" ]] || { usage >&2; exit 1; }
cd "$OUT"
if [[ ${#SAMPLES[@]} -eq 0 ]]; then
  for d in */; do [[ -d "$d/Analysis" ]] && SAMPLES+=("${d%/}"); done
fi
[[ ${#SAMPLES[@]} -gt 0 ]] || { echo "no sample directories with Analysis/ in $OUT" >&2; exit 1; }

# ---------------- optional clean-up of legacy files ----------------
if [[ $CLEAN -eq 1 ]]; then
  mapfile -t old < <(find "${SAMPLES[@]}" -type f \( \
      -path '*/Analysis/*/assembly/*.mapped.*' -o -path '*/Analysis/*/assembly/*.unmapped.*' \
      -o -path '*/Analysis/*/assembly/*.mapping.sam' -o -path '*/Analysis/*/assembly/*.matrix.tsv' \
      -o -path '*/Analysis/*/novel/*.sorted.bam' -o -path '*/Analysis/*/novel/*.sorted.bam.bai' \
      -o -path '*/Preprocessing/*/*.nonhuman.uniq*' \) 2>/dev/null)
  if [[ ${#old[@]} -eq 0 ]]; then
    echo "# no legacy files found"
  else
    du -ch "${old[@]}" | tail -n 1 | sed 's/total/# legacy files, total size:/'
    printf '  %s\n' "${old[@]:0:20}"
    if [[ ${#old[@]} -gt 20 ]]; then echo "  ... ($(( ${#old[@]} - 20 )) more)"; fi
    if [[ $YES -eq 1 ]]; then rm -f "${old[@]}"; echo "# deleted ${#old[@]} files"
    else echo "# dry run — add --yes to delete"; fi
  fi
fi

# ---------------- collect summary files ----------------
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
list=$tmp/files.list
{
  for s in "${SAMPLES[@]}"; do
    # QC / split statistics
    find "$s/Stats" -type f \( -name '*.tsv' -o -name '*.cutadapt.log' -o -name '*.txt' -o -name '*_fastqc.zip' \) -size -"$MAXS" 2>/dev/null || true
    # all analysis modules (human, nonhuman, assembly, novel, fragments, bins, crispr, annotation-free tables)
    find "$s/Analysis" -path "$s/Analysis/annotation" -prune -o -type f \
         \( -name '*.tsv' -o -name '*.report' -o -name '*.tbl' -o -name '*.bed' -o -name '*.RNAfold.txt' \
            -o -name '*_ids.txt' \) -size -"$MAXS" -print 2>/dev/null \
      | grep -vE '/assembly/.*\.(mapped|unmapped)\.|\.matrix\.tsv$' || true
    # inputs and the logs that carry numbers / errors
    find "$s/log_cluster_jobs" -type f \( -name 'inputs.tsv' -o -name '*kraken.log' -o -name '*bracken.log' \
         -o -name 'bowtie.log' -o -name '*.compare.log' -o -name '*.report.log' -o -name '*.srna_crispr.log' \
         -o -name '*.counts.log' \) 2>/dev/null || true
    if [[ $REPORTS -eq 1 && -d "$s/Report" ]]; then
      find "$s/Report" -type f \( -name '*.pdf' -o -name '*.html' -o -name '*_mqc.tsv' -o -name '*.methods.*.md' \)
    fi
  done
} | sort -u > "$list"

# skipped (too large) files and module status
for s in "${SAMPLES[@]}"; do
  find "$s/Analysis" "$s/Stats" -type f \( -name '*.tsv' -o -name '*.bed' -o -name '*.txt' \) -size +"$MAXS" \
       -printf '%s\t%p\n' 2>/dev/null || true
done > "$tmp/skipped_large_files.tsv"
for s in "${SAMPLES[@]}"; do
  find "$s/log_cluster_jobs" -name '.*.done' 2>/dev/null | sed 's/^/DONE /' || true
  grep -rHiE 'error|warn' "$s/log_cluster_jobs" 2>/dev/null | grep -viE 'fastqc|failed to align' | tail -n 100 || true
done > "$tmp/module_status.txt"

n=$(wc -l < "$list")
[[ -n "$TAR" ]] || TAR="srna_summary_$(IFS=_; echo "${SAMPLES[*]:0:4}")$([[ ${#SAMPLES[@]} -gt 4 ]] && echo _etc)_$(date +%Y%m%d_%H%M).tgz"
tar czf "$TAR" -T "$list" -C "$tmp" skipped_large_files.tsv module_status.txt
echo "# samples : ${SAMPLES[*]}"
echo "# files   : $n packed, $(wc -l < "$tmp/skipped_large_files.tsv") larger than $MAXS skipped (listed in skipped_large_files.tsv)"
echo "# modules : $(grep -c '^DONE' "$tmp/module_status.txt") done, $(grep -vc '^DONE' "$tmp/module_status.txt" || true) warning/error lines"
echo "# tarball : $(readlink -f "$TAR") ($(du -h "$TAR" | cut -f1))"
