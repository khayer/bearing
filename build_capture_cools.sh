#!/usr/bin/env bash
#
# build_capture_cools.sh -- build per-condition .cool files from the V1 capture
# HiC-Pro valid-pairs BAMs, for the BEARING 3D side.
#
# Per rep:  samtools view BAM | validpairs_to_pairs.py | sort | bgzip -> pairs.gz
#           cooler cload pairs  (at each resolution)        -> per-rep .cool
# Per cond: cooler merge reps                               -> merged .cool
#           cooler balance --cis-only  (only with --balance)
#
# Output:   <OUTDIR>/cool/merged_capture_<LABEL>_bs_<RES>.cool
#           (LABEL matches samples_v1.tsv: DN, V1P, RCTKO, V1TxS)
#
#   bash build_capture_cools.sh                 # resolutions 2000,5000; no balance
#   bash build_capture_cools.sh --balance       # add cis-only balancing weights
#   bash build_capture_cools.sh --res 1000,2000,5000
#
# Requires: samtools, cooler, bgzip on PATH, and validpairs_to_pairs.py beside
# this script. ASCII only.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"

# --- settings to confirm --------------------------------------------------
BAMDIR="/mnt/isilon/bassing_lab/projects/capture_hic/Brittney_V1/03aln"
BAMPAT="raw_ArimaHTS_S065_%s_rep%s_bs_250_master_valid_pairs.bam"   # cond, rep
CHROMSIZES="/mnt/isilon/bassing_lab/integration_paper/bearing/workflow/resources/mm10.chrom.sizes"
OUTDIR="/mnt/isilon/bassing_lab/integration_paper/capture_cools"
RES="2000,5000"          # ~ the 3 kb usable floor from the depth probe
REPS="1 2"
BALANCE=0
# genotype token in the BAM name  ->  label used in samples_v1.tsv
declare -A LABEL=( [R1KO]=DN [V1PKO]=V1P [RCTKO]=RCTKO [V1TxS]=V1TxS )

while [ $# -gt 0 ]; do
  case "$1" in
    --balance) BALANCE=1; shift ;;
    --res)     RES="$2"; shift 2 ;;
    --res=*)   RES="${1#*=}"; shift ;;
    *) echo "unknown arg: $1" >&2; exit 2 ;;
  esac
done
IFS=',' read -r -a RESARR <<< "$RES"

for t in samtools cooler bgzip; do
  command -v "$t" >/dev/null || { echo "ERROR: $t not on PATH" >&2; exit 1; }
done
mkdir -p "$OUTDIR/pairs" "$OUTDIR/cool"

echo "== capture cools: res=${RES}  balance=${BALANCE}  -> $OUTDIR/cool =="

for COND in "${!LABEL[@]}"; do
  LAB="${LABEL[$COND]}"
  REPCOOLS_BY_RES=()
  for REP in $REPS; do
    BAM="$BAMDIR/$(printf "$BAMPAT" "$COND" "$REP")"
    if [ ! -f "$BAM" ]; then
      echo "  [skip] missing $BAM" >&2
      continue
    fi
    PAIRS="$OUTDIR/pairs/${LAB}_rep${REP}.pairs.gz"
    echo "  [pairs] $COND rep$REP -> $PAIRS"
    samtools view "$BAM" \
      | python3 "$HERE/validpairs_to_pairs.py" "$CHROMSIZES" \
      | sort -k1,1n -k2,2n -k5,5n -k7,7n \
      | bgzip > "$PAIRS"
    for R in "${RESARR[@]}"; do
      RCOOL="$OUTDIR/cool/${LAB}_rep${REP}_bs_${R}.cool"
      echo "  [cool ] $LAB rep$REP @ ${R}bp"
      cooler cload pairs -c1 4 -p1 5 -c2 6 -p2 7 \
        --assembly mm10 "$CHROMSIZES:$R" "$PAIRS" "$RCOOL"
    done
  done
  # merge reps per resolution
  for R in "${RESARR[@]}"; do
    MERGED="$OUTDIR/cool/merged_capture_${LAB}_bs_${R}.cool"
    REPLIST=( "$OUTDIR/cool/${LAB}"_rep*_bs_${R}.cool )
    [ -e "${REPLIST[0]}" ] || { echo "  [skip merge] no rep cools for $LAB @ ${R}" >&2; continue; }
    echo "  [merge] $LAB @ ${R}bp -> $MERGED"
    cooler merge "$MERGED" "${REPLIST[@]}"
    if [ "$BALANCE" = "1" ]; then
      # cis-only: genome-wide balancing is invalid for capture (off-panel empty)
      echo "  [balance] $LAB @ ${R}bp (cis-only)"
      cooler balance --cis-only --force "$MERGED"
    fi
  done
done

echo "== done.  point the cool column of samples_v1.tsv at:"
echo "   ../capture_cools/cool/merged_capture_{COND}_bs_{res}.cool   (COND in DN,V1P,RCTKO,V1TxS)"
