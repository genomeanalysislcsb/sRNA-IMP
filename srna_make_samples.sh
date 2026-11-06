#!/usr/bin/env bash
# srna_make_samples.sh — expand a sample template into a samples TSV for srna_pipeline.sh
#
# Usage: srna_make_samples.sh [-c] [-o out.tsv] template.tsv ids_file
#   template.tsv  comment lines (#...) and the header are copied once; every other line is a
#                 template row and is written once per ID, with
#                   <ID>  replaced by the ID as given in the ids file (e.g. BXXX1)
#                   <id>  replaced by the lower-case ID (e.g. bxxx1)
#   ids_file      one ID per line (blank lines and #comments ignored, duplicates skipped)
#   -o FILE       output file (default: stdout)
#   -c            check that every path column (starting with /) exists; missing files are
#                 reported on stderr and the exit status is 2 (the file is still written)
#   -x            drop rows whose raw_fastq (column 3) does not exist, e.g. samples with only
#                 sRNA-Int or only sRNA-Ext; dropped rows are listed on stderr
set -euo pipefail

usage() { sed -n '2,14p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }
OUT="" CHECK=0 SKIP=0 ARGS=()
while [[ $# -gt 0 ]]; do                       # options may come before or after the files
  case $1 in
    -o) [[ $# -ge 2 ]] || usage 1; OUT=$2; shift 2 ;;
    -c) CHECK=1; shift ;;
    -x) SKIP=1; shift ;;
    -h|--help) usage 0 ;;
    -*) echo "Error: unknown option $1" >&2; usage 1 ;;
    *) ARGS+=("$1"); shift ;;
  esac
done
set -- ${ARGS[@]+"${ARGS[@]}"}
[[ $# -eq 2 ]] || usage 1
TPL=$1 IDS=$2
[[ -s "$TPL" ]] || { echo "Error: template $TPL not found or empty" >&2; exit 1; }
[[ -s "$IDS" ]] || { echo "Error: ids file $IDS not found or empty" >&2; exit 1; }

tmp=$(mktemp); trap 'rm -f "$tmp"' EXIT
awk '
  { sub(/\r$/, "") }
  NR == FNR {                                   # ids file
    id = $0; sub(/[ \t]*#.*/, "", id); gsub(/^[ \t]+|[ \t]+$/, "", id)
    if (id != "" && !(id in seen)) { seen[id] = 1; ids[++n] = id }
    next
  }
  /^[ \t]*#/ || /^[ \t]*$/ { if (!body) print; next }
  !hdr && $1 == "sample_id" { hdr = 1; print; next }   # header line (copied once)
  { body = 1; tpl[++m] = $0 }
  END {
    if (!n) { print "Error: no IDs in ids file" > "/dev/stderr"; exit 1 }
    if (!m) { print "Error: no template rows (with <ID>/<id>) in template" > "/dev/stderr"; exit 1 }
    for (i = 1; i <= n; i++) {
      lo = tolower(ids[i])
      for (j = 1; j <= m; j++) {
        line = tpl[j]; gsub(/<ID>/, ids[i], line); gsub(/<id>/, lo, line); print line
      }
    }
    printf "[make_samples] %d IDs x %d template rows = %d rows\n", n, m, n * m > "/dev/stderr"
  }' "$IDS" "$TPL" > "$tmp"

if [[ $SKIP -eq 1 ]]; then
  kept=$(mktemp); drop=0
  while IFS= read -r line; do
    read -r -a f <<< "$line"
    if [[ ${#f[@]} -ge 3 && ${f[0]} != \#* && ${f[0]} != sample_id && ! -e ${f[2]} ]]; then
      echo "dropped (no fastq): ${f[0]} ${f[1]}: ${f[2]}" >&2; drop=$((drop + 1)); continue
    fi
    printf '%s\n' "$line"
  done < "$tmp" > "$kept"
  mv "$kept" "$tmp"
  echo "[make_samples] $drop row(s) dropped, $(grep -vc -e '^#' -e '^sample_id' "$tmp" || true) kept" >&2
fi

rc=0
if [[ $CHECK -eq 1 ]]; then
  miss=0
  while read -r -a f; do                        # tab or space separated
    [[ ${#f[@]} -eq 0 || ${f[0]} == \#* || ${f[0]} == sample_id ]] && continue
    for p in "${f[@]:2}"; do
      [[ $p == /* && ! -e $p ]] && { echo "missing: ${f[0]} ${f[1]}: $p" >&2; miss=$((miss + 1)); }
    done
  done < "$tmp"
  if [[ $miss -gt 0 ]]; then echo "[make_samples] $miss missing file(s)" >&2; rc=2
  else echo "[make_samples] all paths exist" >&2; fi
fi

if [[ -n "$OUT" ]]; then cp "$tmp" "$OUT"; echo "[make_samples] written: $OUT" >&2; else cat "$tmp"; fi
exit $rc
