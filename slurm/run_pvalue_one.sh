#!/usr/bin/env bash
set -euo pipefail

# Compute differential p-values for ONE comparison against its permutation nulls.
# Companion to run_perm_diff_pvalue.sh (which loops over all comparisons serially);
# this single-comparison form lets Snakemake run the 10 comparisons concurrently
# (rule pvalue_diff_one), dropping wall time from ~37h to roughly one comparison.
#
# Usage:
#   run_pvalue_one.sh COMPARISON_DIR COMP N_PERMS PERM_PREFIX RESULTS_DIR \
#                     [SCORE_METHOD] [PV_MINSIG]
# Set FORCE_PVALUE=1 to recompute even if RESULTS_DIR/diff_COMP.stats.tsv exists.
# ASCII only.

COMPARISON_DIR="${1:?comparison dir}"
COMP="${2:?comparison name}"
N_PERMS="${3:?n_perms}"
PERM_PREFIX="${4:?perm prefix}"
RESULTS_DIR="${5:?results dir}"
SCORE_METHOD="${6:-kl}"
PV_MINSIG="${7:-}"

PV_FLAGS="--score-method ${SCORE_METHOD}"
if [[ -n "$PV_MINSIG" ]]; then
  PV_FLAGS="$PV_FLAGS --min-signal ${PV_MINSIG}"
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
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

mkdir -p "$RESULTS_DIR"
main_qcat="$COMPARISON_DIR/diff_${COMP}.qcat.bgz"
out_prefix="$RESULTS_DIR/diff_${COMP}"

if [[ ! -f "$main_qcat" ]]; then
  echo "ERROR: main qcat not found: $main_qcat" >&2
  exit 1
fi

# Resume guard: skip if this comparison already finished (non-empty stats.tsv).
if [[ "${FORCE_PVALUE:-0}" != "1" && -s "${out_prefix}.stats.tsv" ]]; then
  echo "[resume] ${COMP}: ${out_prefix}.stats.tsv exists, skipping"
  exit 0
fi

# Collect available null qcats across permutations.
nulls=()
for ((p=1; p<=N_PERMS; p++)); do
  null_qcat="${PERM_PREFIX}${p}/diff_comparison/diff_${COMP}.qcat.bgz"
  if [[ ! -f "$null_qcat" ]]; then
    echo "WARNING: null qcat not found for perm$p: $null_qcat" >&2
    continue
  fi
  nulls+=("$null_qcat")
done
if [[ ${#nulls[@]} -eq 0 ]]; then
  echo "ERROR: no null qcats found for comparison: $COMP" >&2
  exit 1
fi

echo "[$(date)] Computing p-values for: $COMP (using ${#nulls[@]} null qcats)"
python "$BEARING_PVALUE" \
  --qcat "$main_qcat" \
  --null-qcat "${nulls[@]}" \
  --diff \
  --out-prefix "$out_prefix" \
  --fdr 0.05 --score-plot $PV_FLAGS
echo "[$(date)] Completed: $COMP"
