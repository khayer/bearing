#!/usr/bin/env bash
#
# run_dchic_convert.sh -- convert per-replicate HiCExplorer matrices to dcHiC
# input, for every condition, at 100 kb and 50 kb.
#
# Safe to leave running. It:
#   - DISCOVERS the files rather than assuming names (the exact per-replicate
#     sample names in 04mat_cor were never confirmed, so hard-coding them would
#     silently skip conditions)
#   - skips anything already converted, so it can be re-run after an interrupt
#   - keeps going if one sample fails, and lists the failures at the end
#   - logs everything to dchic_convert.log
#
# Reads:  04mat_cor/raw_corrected_KR_<SAMPLE>_bs_<RES>.h5   (per-replicate, KR-balanced)
# Writes: dchic_in_<RES-label>/<SAMPLE>.matrix + .bed, and input.txt
#
# NOTE ON CHROMOSOMES: restricted to the nine chromosomes the existing
# compartment analysis used (from the hic_pipe run_compartment_analysis rule),
# so this run is comparable to the published PC1 tracks and carries a lighter
# multiple-testing burden.
#
# NOTE ON BALANCED INPUT: these are KR-corrected values, not raw counts. That
# matches what hicPCA saw, but dcHiC's spec calls its third column "count".
# If the answer comes out marginal, re-run from the raw masters as a check.

set -u
set -o pipefail

# ---------------------------------------------------------------- configuration
MATDIR="../../projects/HiC_V31_NT_allele/data/endpoints/HiC_explorer_mm10/04mat_cor"
CONVERTER="h5_to_dchic.py"
CHROMS="chr2,chr6,chr12,chr13,chr14,chr16,chr17,chr18,chr19"
LOG="dchic_convert.log"
# Parent dir for the dchic_in_<res> scratch. Defaults to the repo root (cwd) for
# a standalone run; the Snakemake rule points it at the results outdir.
WORK="${DCHIC_WORKDIR:-$(pwd)}"
mkdir -p "$WORK"

# Which samples to include. Generous on purpose -- edit after reading the
# DISCOVERY block in the log if something is missing or something unwanted
# slipped in. Matched against the sample name, case-insensitively.
# Matched anywhere in the sample name: an earlier version anchored on (^|_) and
# silently skipped every "ArimaDN_rep1"-style name, which is most of them.
INCLUDE='(DN|DP|EBKO|EbKO|ProB|S3T3|s3T3|dV1P|dV1CTCF|V1PRO|V1SCR|TKO)'
# Non-Tcrb experiments in the same directory, plus the alternative Pro-B and DN
# preps (BAA/BDA/BA...) that are NOT the ones behind merged_corrected_KR_ProB /
# _DN. If a condition you need is missing from the DISCOVERY block, loosen this.
EXCLUDE='(degron|MNK3|4GyIR|P5424|S6312|V31NT|d5PC|EbPrss2|ArimaBAA|ArimaBDA|ArimaBAProB|ArimaBADN|ArimaBAEbDN|ArimaBA)'

declare -A RES_LABEL=( [100000]=100kb [50000]=50kb [250000]=250kb )

# ---------------------------------------------------------------- sanity checks
exec > >(tee -a "$LOG") 2>&1
echo "=============================================================="
echo "run_dchic_convert.sh   $(date)"
echo "=============================================================="

if [ ! -d "$MATDIR" ]; then
  echo "FATAL: matrix directory not found: $MATDIR"; exit 1
fi
if [ ! -f "$CONVERTER" ]; then
  echo "FATAL: $CONVERTER not found in $(pwd)"; exit 1
fi
python3 -c "import hicmatrix" 2>/dev/null || {
  echo "FATAL: python cannot import 'hicmatrix'."
  echo "       activate the environment that has HiCExplorer, then re-run."
  exit 1
}

FAILED=()
CONVERTED=0
SKIPPED=0

# ---------------------------------------------------------------- main loop
for RES in 100000 50000 250000; do
  LAB=${RES_LABEL[$RES]}
  OUTDIR="$WORK/dchic_in_${LAB}"
  mkdir -p "$OUTDIR"

  echo ""
  echo "##############################################################"
  echo "# RESOLUTION $LAB  ->  $OUTDIR"
  echo "##############################################################"

  # ---- discovery -------------------------------------------------------
  mapfile -t ALL < <(ls -1 "$MATDIR"/raw_corrected_KR_*_bs_${RES}.h5 2>/dev/null | sort)
  if [ ${#ALL[@]} -eq 0 ]; then
    echo "  no files match $MATDIR/raw_corrected_KR_*_bs_${RES}.h5 -- skipping $LAB"
    continue
  fi

  echo ""
  echo "  DISCOVERY: ${#ALL[@]} file(s) at ${RES} bp"
  KEEP=()
  for F in "${ALL[@]}"; do
    B=$(basename "$F")
    S=${B#raw_corrected_KR_}
    S=${S%_bs_${RES}.h5}
    if [[ "$S" =~ $EXCLUDE ]]; then
      printf "    %-42s  skip (excluded)\n" "$S"; continue
    fi
    if [[ ! "$S" =~ $INCLUDE ]]; then
      printf "    %-42s  skip (not in INCLUDE)\n" "$S"; continue
    fi
    printf "    %-42s  INCLUDE\n" "$S"
    KEEP+=("$S")
  done

  if [ ${#KEEP[@]} -eq 0 ]; then
    echo "  nothing selected at $LAB -- check the INCLUDE pattern against the list above"
    continue
  fi

  # ---- convert ---------------------------------------------------------
  echo ""
  echo "  CONVERTING ${#KEEP[@]} sample(s)"
  for S in "${KEEP[@]}"; do
    F="$MATDIR/raw_corrected_KR_${S}_bs_${RES}.h5"
    if [ -s "$OUTDIR/${S}.matrix" ] && [ -s "$OUTDIR/${S}.bed" ]; then
      echo "    [skip] $S (already converted)"
      SKIPPED=$((SKIPPED+1))
      continue
    fi
    echo "    [run ] $S"
    if python3 "$CONVERTER" --h5 "$F" --sample "$S" \
         --outdir "$OUTDIR" --chroms "$CHROMS"; then
      CONVERTED=$((CONVERTED+1))
    else
      echo "    [FAIL] $S"
      FAILED+=("$LAB/$S")
      rm -f "$OUTDIR/${S}.matrix" "$OUTDIR/${S}.bed"   # no half-written inputs
    fi
  done

  # ---- build the dcHiC input file --------------------------------------
  echo ""
  echo "  BUILDING $OUTDIR/input.txt"
  python3 "$CONVERTER" --write-input "$OUTDIR" \
      --out "$OUTDIR/input.txt" --res-label "$LAB" || \
      echo "    WARNING: could not build input.txt for $LAB"
done

# ---------------------------------------------------------------- summary
echo ""
echo "=============================================================="
echo "DONE  $(date)"
echo "  converted: $CONVERTED"
echo "  skipped  : $SKIPPED (already present)"
if [ ${#FAILED[@]} -gt 0 ]; then
  echo "  FAILED   : ${#FAILED[@]}"
  for f in "${FAILED[@]}"; do echo "     $f"; done
else
  echo "  failed   : 0"
fi
echo ""
echo "NEXT -- check input.txt lists the conditions you expect, with >=2"
echo "replicates each, then run dcHiC per resolution:"
echo ""
for LAB in 100kb 50kb; do
cat <<EOF
  Rscript dchicf.r --file $WORK/dchic_in_${LAB}/input.txt --pcatype cis     --dirovwt T --cthread 2 --pthread 4
  Rscript dchicf.r --file $WORK/dchic_in_${LAB}/input.txt --pcatype select  --dirovwt T --genome mm10
  Rscript dchicf.r --file $WORK/dchic_in_${LAB}/input.txt --pcatype analyze --dirovwt T --diffdir comp_${LAB}
  Rscript dchicf.r --file $WORK/dchic_in_${LAB}/input.txt --pcatype viz     --diffdir comp_${LAB} --genome mm10

EOF
done
echo "Then, for the V1 question:"
echo "  python3 dchic_region_result.py \\"
echo "     --bedgraph DifferentialResult/comp_100kb/fdr_result/differential.intra_sample_combined.bedGraph \\"
echo "     --region chr6:40872100-40908238 --label V1 --exp-a <CONTROL> --exp-b dV1P \\"
echo "     --out paper/table_sources/v1_compartment_dchic_100kb.tsv"
echo "=============================================================="
