#!/usr/bin/env bash
# Read-only staleness audit for the PRODUCTION results/ tree.
# Flags products DOWNSTREAM of the clean p-value layer that were built BEFORE it
# (mtime < pvalue.done). Does NOT delete or rebuild anything.
set -u
RES="workflow/results"
REF="$RES/pvalue.done"

if [ ! -f "$REF" ]; then echo "ERROR: $REF not found (run from repo root on cluster)"; exit 1; fi
REF_EPOCH=$(stat -c %Y "$REF")
echo "Reference: pvalue.done = $(date -d @"$REF_EPOCH" '+%Y-%m-%d %H:%M:%S')"
echo "Rule: a DOWNSTREAM file older than this was built from the poisoned layer."
echo

# Classification of immediate children of results/.
# DOWNSTREAM = consumes pvalue.done / diff tables -> STALE if older than REF.
# UPSTREAM   = feeds the p-value layer -> older is EXPECTED, not stale.
# PARALLEL   = independent branch (own inputs/own pvalue) -> not main-layer stale.
declare -A ROLE
for d in regional loops figures paper_figures crosslocus; do ROLE[$d]=DOWNSTREAM; done
for f in decomposition.done; do ROLE[$f]=DOWNSTREAM; done
for d in blacklist score compare perm; do ROLE[$d]=UPSTREAM; done
for f in score.done compare.done region_qc.done pvalue.done pvalue_samples.done samples.bearing.tsv; do ROLE[$f]=UPSTREAM; done
for d in calibration benchmark consensusQ hic; do ROLE[$d]=PARALLEL; done
for f in calibration.done; do ROLE[$f]=PARALLEL; done

printf "%-34s %-10s %-8s %s\n" "CHILD" "ROLE" "OLD/TOT" "VERDICT"
printf "%-34s %-10s %-8s %s\n" "-----" "----" "-------" "-------"
for child in $(ls -1 "$RES" 2>/dev/null | sort); do
  path="$RES/$child"
  role="${ROLE[$child]:-UNCLASSIFIED}"
  if [ -d "$path" ]; then
    tot=$(find "$path" -type f 2>/dev/null | wc -l)
    old=$(find "$path" -type f ! -newer "$REF" 2>/dev/null | wc -l)
  else
    tot=1
    if [ "$(stat -c %Y "$path")" -le "$REF_EPOCH" ]; then old=1; else old=0; fi
  fi
  verdict=""
  if [ "$role" = "DOWNSTREAM" ] && [ "$old" -gt 0 ]; then verdict="<<< STALE: $old file(s) predate clean layer"; fi
  if [ "$role" = "UNCLASSIFIED" ] && [ "$old" -gt 0 ]; then verdict="??? eyeball ($old old)"; fi
  printf "%-34s %-10s %-8s %s\n" "$child" "$role" "$old/$tot" "$verdict"
done

echo
echo "=== DOWNSTREAM files older than the clean layer (the actual rebuild list) ==="
for child in regional loops figures paper_figures crosslocus; do
  p="$RES/$child"
  [ -d "$p" ] || continue
  find "$p" -type f ! -newer "$REF" -printf '%TY-%Tm-%Td %TH:%TM  %p\n' 2>/dev/null
done | sort
