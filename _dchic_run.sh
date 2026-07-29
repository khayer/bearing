#!/usr/bin/env bash
#
# _dchic_run.sh -- shared dcHiC driver, sourced by the run_all_*.sh scripts.
#
# Not run directly. Each caller sets:
#   INPUT_PREFIX   e.g. input_v1          (expects input_v1.txt in the res dir)
#   DIFF_PREFIX    e.g. v1                 (diffdir becomes v1_100kb, v1_50kb)
#   RESOLUTIONS    e.g. "100kb 50kb"
# then calls   dchic_run_all
# and optionally  dchic_region  <region> <label> <expA> <expB>  per resolution.
#
# KEY LESSON BAKED IN: dcHiC's cis step dies with
#   "Two levels of parallelism are used"
# if BOTH --cthread and --pthread are > 1. So cthread is pinned to 1 and only
# pthread parallelizes. Do not "optimize" this back to cthread 2.

set -euo pipefail

# ------- locate dcHiC ----------------------------------------------------
: "${DCHIC:=$HOME/tools/dcHiC/dchicf.r}"
if [ ! -f "$DCHIC" ]; then
  echo "FATAL: dchicf.r not found at $DCHIC. Set DCHIC=/path/to/dchicf.r" >&2
  exit 1
fi

# ------- require the dchic conda env -------------------------------------
if ! command -v Rscript >/dev/null 2>&1; then
  echo "FATAL: Rscript not on PATH. Run: conda activate dchic" >&2
  exit 1
fi
if ! Rscript -e 'library(functionsdchic)' >/dev/null 2>&1; then
  echo "FATAL: R cannot load functionsdchic. Run: conda activate dchic" >&2
  exit 1
fi

# cis parallelism. dcHiC's cis step raises "Two levels of parallelism are used"
# whenever an outer pool (pthread>1) is live while its inner eigendecomposition
# threads. At 100 kb the inner work is fast enough to slip through with pthread=4;
# at 50 kb (4x larger matrices) it collides and the step dies partway. The only
# combination guaranteed free of nesting is fully serial. Default to that; a
# caller whose build tolerates more can export PTHREAD=4 for the coarse runs.
PTHREAD="${PTHREAD:-1}"
GENOME="${GENOME:-mm10}"

# absolute path to the repo root (where dchic_region_result.py and paper/ live).
# The region-extract cds into dchic_in_<res>/, so a relative path would break.
REPO_ROOT="${REPO_ROOT:-$(pwd)}"
# Parent dir for the dchic_in_<res> scratch (matches run_dchic_convert.sh /
# build_dchic_inputs.sh). Defaults to the repo root; the Snakemake rule points
# it at the results outdir.
DCHIC_WORKDIR="${DCHIC_WORKDIR:-$REPO_ROOT}"
REGION_SCRIPT="${REGION_SCRIPT:-$REPO_ROOT/dchic_region_result.py}"
TABLE_OUT="${TABLE_OUT:-$REPO_ROOT/paper/table_sources}"

dchic_run_all() {
  for RES in $RESOLUTIONS; do
    local INPUT="${INPUT_PREFIX}.txt"
    local DIFF="${DIFF_PREFIX}_${RES}"
    local DIR="$DCHIC_WORKDIR/dchic_in_${RES}"

    echo ""
    echo "############################################################"
    echo "# ${DIFF_PREFIX}  @  ${RES}   (input ${DIR}/${INPUT})"
    echo "############################################################"

    if [ ! -d "$DIR" ]; then
      echo "  SKIP: $DIR does not exist (convert this resolution first)"; continue
    fi
    if [ ! -s "$DIR/$INPUT" ]; then
      echo "  SKIP: $DIR/$INPUT missing or empty"; continue
    fi

    ( cd "$DIR"

      echo "  [cis]     $(date +%H:%M:%S)  (cthread 1, pthread $PTHREAD)"
      Rscript "$DCHIC" --file "$INPUT" --pcatype cis \
          --dirovwt T --cthread 1 --pthread "$PTHREAD" > "cis_${DIFF}.log" 2>&1
      local npc
      npc=$(find . -name "*.pc.txt" | wc -l)
      echo "            .pc.txt files: $npc"
      if [ "$npc" -eq 0 ]; then
        echo "  cis produced no PC files -- see $DIR/cis_${DIFF}.log" >&2
        tail -5 "cis_${DIFF}.log" >&2; exit 1
      fi

      echo "  [select]  $(date +%H:%M:%S)"
      Rscript "$DCHIC" --file "$INPUT" --pcatype select \
          --dirovwt T --genome "$GENOME" > "select_${DIFF}.log" 2>&1

      echo "  [analyze] $(date +%H:%M:%S)"
      Rscript "$DCHIC" --file "$INPUT" --pcatype analyze \
          --dirovwt T --diffdir "$DIFF" > "analyze_${DIFF}.log" 2>&1

      echo "  [viz]     $(date +%H:%M:%S)"
      Rscript "$DCHIC" --file "$INPUT" --pcatype viz \
          --diffdir "$DIFF" --genome "$GENOME" > "viz_${DIFF}.log" 2>&1

      local BG="DifferentialResult/${DIFF}/fdr_result/differential.intra_sample_combined.pcQnm.bedGraph"
      if [ -s "$BG" ]; then
        local nsig
        nsig=$(tail -n +2 "DifferentialResult/${DIFF}/fdr_result/differential.intra_sample_combined.Filtered.pcQnm.bedGraph" 2>/dev/null | wc -l)
        echo "  DONE: $BG"
        echo "        genome-wide significant bins (Filtered): $nsig"
      else
        echo "  WARNING: expected bedGraph not found: $BG" >&2
      fi
    )
  done
}

# dchic_region <region> <label> <expA> <expB>  -- run per resolution
dchic_region() {
  local REGION="$1" LABEL="$2" EXPA="$3" EXPB="$4"
  for RES in $RESOLUTIONS; do
    local DIFF="${DIFF_PREFIX}_${RES}"
    local DIR="$DCHIC_WORKDIR/dchic_in_${RES}"
    local BG="$DIR/DifferentialResult/${DIFF}/fdr_result/differential.intra_sample_combined.pcQnm.bedGraph"
    [ -s "$BG" ] || { echo "  region: no bedGraph for $DIFF, skipping"; continue; }
    echo ""
    echo "=== ${LABEL}  ${EXPA} vs ${EXPB}  @ ${RES} ==="
    if [ ! -f "$REGION_SCRIPT" ]; then
      echo "  region: $REGION_SCRIPT not found; set REPO_ROOT or REGION_SCRIPT" >&2
      continue
    fi
    mkdir -p "$TABLE_OUT"
    python3 "$REGION_SCRIPT" \
        --bedgraph "$BG" --region "$REGION" --label "$LABEL" \
        --exp-a "$EXPA" --exp-b "$EXPB" \
        --out "${TABLE_OUT}/${LABEL}_compartment_dchic_${EXPB}_${RES}.tsv" \
      || echo "  (region extract failed for $DIFF)"
  done
}
