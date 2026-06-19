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

# Collect the null qcats across permutations. A PARTIAL perm set must NOT
# silently produce p-values (that is how stale/degenerate stats were written
# before): require all N_PERMS to be present and non-empty unless the caller
# explicitly sets ALLOW_PARTIAL_NULL=1.
nulls=()
missing=()
empty=()
for ((p=1; p<=N_PERMS; p++)); do
  null_qcat="${PERM_PREFIX}${p}/diff_comparison/diff_${COMP}.qcat.bgz"
  if [[ ! -f "$null_qcat" ]]; then
    missing+=("perm$p")
    continue
  fi
  if [[ ! -s "$null_qcat" ]]; then
    empty+=("perm$p")
    continue
  fi
  nulls+=("$null_qcat")
done
if [[ ${#missing[@]} -gt 0 || ${#empty[@]} -gt 0 ]]; then
  echo "ERROR: incomplete permutation null for $COMP:" >&2
  echo "  expected $N_PERMS perms; have ${#nulls[@]};" \
       "missing ${#missing[@]} (${missing[*]:-none});" \
       "empty ${#empty[@]} (${empty[*]:-none})." >&2
  if [[ "${ALLOW_PARTIAL_NULL:-0}" != "1" ]]; then
    echo "  Refusing to compute p-values against a partial null." \
         "Rebuild the missing perms, or set ALLOW_PARTIAL_NULL=1 to override." >&2
    exit 1
  fi
  echo "  ALLOW_PARTIAL_NULL=1 set -- proceeding with ${#nulls[@]} perms (NOT for production)." >&2
fi
if [[ ${#nulls[@]} -eq 0 ]]; then
  echo "ERROR: no null qcats found for comparison: $COMP" >&2
  exit 1
fi

# Pass the expected perm count to bearing_pvalue.py as a second, independent
# gate (count must match unless partial was explicitly allowed).
EXPECT_FLAG=""
if [[ "${ALLOW_PARTIAL_NULL:-0}" != "1" ]]; then
  EXPECT_FLAG="--expect-n-perms $N_PERMS"
fi

echo "[$(date)] Computing p-values for: $COMP (using ${#nulls[@]} null qcats)"
python "$BEARING_PVALUE" \
  --qcat "$main_qcat" \
  --null-qcat "${nulls[@]}" \
  --diff \
  --out-prefix "$out_prefix" \
  --fdr 0.05 --score-plot $EXPECT_FLAG $PV_FLAGS
echo "[$(date)] Completed: $COMP"
