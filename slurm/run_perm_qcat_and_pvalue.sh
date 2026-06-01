#!/usr/bin/env bash
set -euo pipefail

# Build null qcat files from existing perm*/ shifted BigWigs and run bearing_pvalue.py
# for each sample listed in a sheet with columns: sample, bw, out.
#
# Usage:
#   ./run_perm_qcat_and_pvalue.sh [sheet.tsv] [n_perms] [perm_prefix] [results_dir]
#
# Example:
#   ./run_perm_qcat_and_pvalue.sh samples_template.tsv 4 perm results

SHEET="${1:-samples_template.tsv}"
N_PERMS="${2:-4}"
PERM_PREFIX="${3:-perm}"
RESULTS_DIR="${4:-results}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIGWIG_TO_QCAT="$SCRIPT_DIR/bigwig_to_qcat.py"
BEARING_PVALUE="$SCRIPT_DIR/bearing_pvalue.py"

if [[ ! -f "$BIGWIG_TO_QCAT" ]]; then
  echo "ERROR: missing helper script: $BIGWIG_TO_QCAT" >&2
  exit 1
fi

if [[ ! -f "$BEARING_PVALUE" ]]; then
  echo "ERROR: missing helper script: $BEARING_PVALUE" >&2
  exit 1
fi

if [[ ! -f "$SHEET" ]]; then
  echo "ERROR: sheet not found: $SHEET" >&2
  exit 1
fi

mkdir -p "$RESULTS_DIR"

awk -F '\t' '
  BEGIN {
    OFS="\t"
    sample_i=0
    bw_i=0
    out_i=0
  }

  /^[[:space:]]*#/ || NF==0 { next }

  sample_i==0 {
    for (i=1; i<=NF; i++) {
      h=$i
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", h)
      hl=tolower(h)
      if (hl=="sample" || hl=="name") sample_i=i
      else if (hl=="bw") bw_i=i
      else if (hl=="out") out_i=i
    }
    if (sample_i==0 || bw_i==0 || out_i==0) {
      print "ERROR: sheet must contain columns: sample, bw, out" > "/dev/stderr"
      exit 2
    }
    next
  }

  {
    sample=$sample_i
    bw=$bw_i
    out=$out_i
    gsub(/^[[:space:]]+|[[:space:]]+$/, "", sample)
    gsub(/^[[:space:]]+|[[:space:]]+$/, "", bw)
    gsub(/^[[:space:]]+|[[:space:]]+$/, "", out)
    if (sample=="" || bw=="" || out=="") next
    print sample, bw, out
  }
' "$SHEET" | while IFS=$'\t' read -r sample bw_csv out; do
  if [[ -z "$sample" || -z "$bw_csv" || -z "$out" ]]; then
    echo "Skipping malformed row: sample=$sample out=$out" >&2
    continue
  fi

  out_name=${out##*/}
  base=${out_name%.qcat.bgz}
  nulls=()

  IFS=',' read -r -a local_bws <<< "$bw_csv"

  for ((p=1; p<=N_PERMS; p++)); do
    shifted=()
    for bw in "${local_bws[@]}"; do
      fn=${bw##*/}
      stem=${fn%.bw}
      shifted+=("${PERM_PREFIX}${p}/${sample}/${stem}_perm${p}.bw")
    done

    null_qcat="${PERM_PREFIX}${p}/${sample}/${base}_perm${p}.qcat.bgz"
    echo "[${sample}] building $null_qcat"
    python "$BIGWIG_TO_QCAT" --bw "${shifted[@]}" --out "$null_qcat"
    nulls+=("$null_qcat")
  done

  echo "[${sample}] running p-value with ${#nulls[@]} null qcats"
  python "$BEARING_PVALUE" \
    --qcat "$out" \
    --null-qcat "${nulls[@]}" \
    --out-prefix "$RESULTS_DIR/$sample" \
    --fdr 0.05 --score-plot
done

echo "Done. Outputs written under: $RESULTS_DIR/"
