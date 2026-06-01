#!/usr/bin/env bash
set -euo pipefail

# Automated p-value computation for all differential comparisons
# Processes each diff_*.qcat.bgz against permuted nulls

COMPARISON_DIR="${1:-comparison_vdj}"
N_PERMS="${2:-4}"
PERM_PREFIX="${3:-perm}"
RESULTS_DIR="${4:-results}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# bearing_pvalue.py lives at the repo root; this script may sit in slurm/.
# Prefer a sibling copy, else look one directory up (repo root).
if [[ -f "$SCRIPT_DIR/bearing_pvalue.py" ]]; then
  BEARING_PVALUE="$SCRIPT_DIR/bearing_pvalue.py"
elif [[ -f "$SCRIPT_DIR/../bearing_pvalue.py" ]]; then
  BEARING_PVALUE="$(cd "$SCRIPT_DIR/.." && pwd)/bearing_pvalue.py"
else
  BEARING_PVALUE="$SCRIPT_DIR/bearing_pvalue.py"
fi

if [[ ! -f "$BEARING_PVALUE" ]]; then
  echo "ERROR: missing script: $BEARING_PVALUE" >&2
  exit 1
fi

if [[ ! -d "$COMPARISON_DIR" ]]; then
  echo "ERROR: comparison directory not found: $COMPARISON_DIR" >&2
  exit 1
fi

mkdir -p "$RESULTS_DIR"

# Extract all unique comparison names from the comparison directory
comparisons=$(ls -1 "$COMPARISON_DIR"/diff_*.qcat.bgz 2>/dev/null | sed 's|.*/diff_||; s|\.qcat\.bgz||' | sort)

if [[ -z "$comparisons" ]]; then
  echo "ERROR: no diff_*.qcat.bgz files found in $COMPARISON_DIR" >&2
  exit 1
fi

echo "Found comparisons:"
echo "$comparisons"
echo ""

# Process each comparison
while IFS= read -r comp; do
  if [[ -z "$comp" ]]; then
    continue
  fi
  
  main_qcat="$COMPARISON_DIR/diff_${comp}.qcat.bgz"
  out_prefix="$RESULTS_DIR/diff_${comp}"
  
  if [[ ! -f "$main_qcat" ]]; then
    echo "WARNING: main qcat not found: $main_qcat" >&2
    continue
  fi
  
  # Build array of null qcats from all permutations
  nulls=()
  for ((p=1; p<=N_PERMS; p++)); do
    null_qcat="${PERM_PREFIX}${p}/diff_comparison/diff_${comp}.qcat.bgz"
    if [[ ! -f "$null_qcat" ]]; then
      echo "WARNING: null qcat not found for perm$p: $null_qcat" >&2
      continue
    fi
    nulls+=("$null_qcat")
  done
  
  if [[ ${#nulls[@]} -eq 0 ]]; then
    echo "ERROR: no null qcats found for comparison: $comp" >&2
    continue
  fi
  
  echo "[$(date)] Computing p-values for: $comp (using ${#nulls[@]} null qcats)"
  if python "$BEARING_PVALUE" \
    --qcat "$main_qcat" \
    --null-qcat "${nulls[@]}" \
    --diff \
    --out-prefix "$out_prefix" \
    --fdr 0.05 --score-plot; then
    echo "[$(date)] Completed: $comp"
  else
    echo "[$(date)] WARNING: p-value computation skipped for $comp (no bins with signal, or insufficient nulls); continuing." >&2
  fi
  echo ""
done <<< "$comparisons"

echo "Done. All results in: $RESULTS_DIR/"
