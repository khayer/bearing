#!/usr/bin/env bash
# ----------------------------------------------------------------------
# File     : benchmark/run_complementarity_sweep.sh
# Version  : 1.0.0
# Date     : 2026-10-11
# Authors  : Katharina E. Hayer (katharinaehayer@gmail.com) and Claude
#            (Anthropic), co-created
# ----------------------------------------------------------------------
# Two-condition complementarity benchmark (BEARING vs per-assay test) over
# depth x composition prior x seed, then one pooled table
# (complementarity_table.tsv, candidate Table S19) and a per-configuration
# table (complementarity_per_config.tsv).
#
# SIMULATED data only. CPU only, ~15 s and < 2 GB per cell at the defaults
# (20 Mb, 100 permutations); 18 cells ~ 5 min. No cluster needed.
#
# Usage (from the repository root):
#   bash benchmark/run_complementarity_sweep.sh [OUTROOT]
# Grids can be overridden: DEPTHS="2 8 32" PRIORS="0 1" SEEDS="1 2 3"
#   NPERMS=100 LENGTH=20000000
# ASCII-only.

set -euo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
PY=${PYTHON:-python3}
OUTROOT=${1:-complementarity_sweep}
DEPTHS=${DEPTHS:-"2 8 32"}
PRIORS=${PRIORS:-"0 1"}
SEEDS=${SEEDS:-"1 2 3"}
NPERMS=${NPERMS:-100}
LENGTH=${LENGTH:-20000000}

mkdir -p "$OUTROOT"
for a in $PRIORS; do
  for d in $DEPTHS; do
    for s in $SEEDS; do
      od="$OUTROOT/prior${a}_depth${d}_seed${s}"
      if [ -s "$od/per_type.tsv" ]; then echo "[skip] $od"; continue; fi
      echo ">>> prior=$a depth=$d seed=$s"
      "$PY" "$REPO/benchmark/complementarity_benchmark.py" --outdir "$od" \
          --prior-strength "$a" --depth "$d" --seed "$s" \
          --n-perms "$NPERMS" --length "$LENGTH" > "$od.log" 2>&1 \
        || { echo "FAILED: see $od.log" >&2; exit 1; }
      sed -n 3,5p "$od/summary.txt" | sed 's/^/    /'
    done
  done
done

"$PY" "$REPO/benchmark/summarize_complementarity.py" --root "$OUTROOT" \
    --out "$OUTROOT/complementarity_table.tsv" \
    --per-config-out "$OUTROOT/complementarity_per_config.tsv"
