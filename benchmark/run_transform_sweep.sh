#!/usr/bin/env bash
# ----------------------------------------------------------------------
# File     : benchmark/run_transform_sweep.sh
# Version  : 1.0.0
# Date     : 2026-11-08
# Authors  : Katharina E. Hayer (katharinaehayer@gmail.com) and Claude
#            (Anthropic), co-created
# ----------------------------------------------------------------------
# Complementarity benchmark (1.2.0) over input transform x RNA tail x depth
# x seed, then one table (transform_sweep_table.tsv) built by
# summarize_transform_sweep.py.
#
# SIMULATED data only. CPU only, ~10-35 s and < 2 GB per cell (20 Mb, 100
# permutations, prior 1, min_signal 0.1, floor 0.5, min-shift 1.5 Mb).
# 4 x 2 x 3 x 3 = 72 cells; ~15 min with JOBS=4.
#
# Usage (from the repository root):
#   bash benchmark/run_transform_sweep.sh [OUTROOT]
# Grids: TRANSFORMS="none nzq log1p_rna asinh_auto" TAILS="0 1.5"
#        DEPTHS="2 8 32" SEEDS="1 2 3" NPERMS=100 JOBS=2
# ASCII-only.

set -euo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
PY=${PYTHON:-python3}
OUTROOT=${1:-transform_sweep}
TRANSFORMS=${TRANSFORMS:-"none nzq log1p_rna asinh_auto"}
TAILS=${TAILS:-"0 1.5"}
DEPTHS=${DEPTHS:-"2 8 32"}
SEEDS=${SEEDS:-"1 2 3"}
NPERMS=${NPERMS:-100}
JOBS=${JOBS:-2}

mkdir -p "$OUTROOT"
cells="$OUTROOT/cells.txt"
: > "$cells"
for tf in $TRANSFORMS; do
  for tl in $TAILS; do
    for d in $DEPTHS; do
      for s in $SEEDS; do
        echo "$tf $tl $d $s" >> "$cells"
      done
    done
  done
done

run_one() {
  tf=$1; tl=$2; d=$3; s=$4
  od="$OUTROOT/tf_${tf}_tail${tl}_d${d}_s${s}"
  if [ -s "$od/per_type.tsv" ]; then echo "[skip] $od"; return 0; fi
  "$PY" "$REPO/benchmark/complementarity_benchmark.py" --outdir "$od" \
      --transform "$tf" --rna-tail-sigma "$tl" --depth "$d" --seed "$s" \
      --n-perms "$NPERMS" --prior-strength 1 --min-signal 0.1 \
      --pvalue-min-signal 0.5 --min-shift 1500000 > "$od.log" 2>&1 \
    || { echo "FAILED: see $od.log" >&2; return 1; }
  echo "[done] $od"
}
export -f run_one
export PY REPO OUTROOT NPERMS
xargs -P "$JOBS" -L 1 bash -c 'run_one $0 $1 $2 $3' < "$cells"

"$PY" "$REPO/benchmark/summarize_transform_sweep.py" --root "$OUTROOT" \
    --out "$OUTROOT/transform_sweep_table.tsv"
