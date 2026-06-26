#!/usr/bin/env bash
#
# build_capture_cools.sh -- per-condition .cool files from the V1 capture
# HiC-Pro valid-pairs BAMs, for visualization / the BEARING 3D side.
#
# Pass A:  samtools view BAM | validpairs_to_pairs.py | sort | bgzip -> pairs.gz
#          (and count pooled pairs per condition)
# --downsample: floor = min condition total; each condition kept to floor so all
#          matrices share the SAME pooled depth (same scheme the benchmark used,
#          keep = floor/condition_total) -- REQUIRED for a fair DN-vs-mutant figure.
# Then:    cooler cload pairs (each resolution) -> per-rep .cool
#          cooler merge reps -> merged_capture_<LABEL>_bs_<RES>.cool
#          cooler balance --cis-only  (only with --balance)
#
#   bash build_capture_cools.sh --downsample        # depth-matched (use for figures)
#   bash build_capture_cools.sh                      # native depth (single-cond views)
#   bash build_capture_cools.sh --downsample --res 2000 --balance
#
# Requires: samtools, cooler, bgzip on PATH; validpairs_to_pairs.py beside this.
# ASCII only.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"

BAMDIR="/mnt/isilon/bassing_lab/projects/capture_hic/Brittney_V1/03aln"
BAMPAT="raw_ArimaHTS_S065_%s_rep%s_bs_250_master_valid_pairs.bam"   # cond, rep
CHROMSIZES="/mnt/isilon/bassing_lab/integration_paper/bearing/workflow/resources/mm10.chrom.sizes"
OUTDIR="/mnt/isilon/bassing_lab/integration_paper/capture_cools"
RES="2000,5000"
REPS="1 2"
BALANCE=0
DOWNSAMPLE=0
SEED=42
declare -A LABEL=( [R1KO]=DN [V1PKO]=V1P [RCTKO]=RCTKO [V1TxS]=V1TxS )

while [ $# -gt 0 ]; do
  case "$1" in
    --balance)    BALANCE=1; shift ;;
    --downsample) DOWNSAMPLE=1; shift ;;
    --res)        RES="$2"; shift 2 ;;
    --res=*)      RES="${1#*=}"; shift ;;
    --seed)       SEED="$2"; shift 2 ;;
    *) echo "unknown arg: $1" >&2; exit 2 ;;
  esac
done
IFS=',' read -r -a RESARR <<< "$RES"
for t in samtools cooler bgzip; do
  command -v "$t" >/dev/null || { echo "ERROR: $t not on PATH" >&2; exit 1; }
done
mkdir -p "$OUTDIR/pairs" "$OUTDIR/cool"

# ---- PASS A: build per-rep pairs, count pooled pairs per condition -----------
declare -A CONDTOTAL
echo "== pass A: pairs + per-condition depth =="
for COND in "${!LABEL[@]}"; do
  LAB="${LABEL[$COND]}"
  tot=0
  for REP in $REPS; do
    BAM="$BAMDIR/$(printf "$BAMPAT" "$COND" "$REP")"
    [ -f "$BAM" ] || { echo "  [skip] missing $BAM" >&2; continue; }
    PAIRS="$OUTDIR/pairs/${LAB}_rep${REP}.pairs.gz"
    samtools view "$BAM" \
      | python3 "$HERE/validpairs_to_pairs.py" "$CHROMSIZES" \
      | sort -k1,1n -k2,2n -k5,5n -k7,7n | bgzip > "$PAIRS"
    n=$(zcat "$PAIRS" | wc -l)
    tot=$((tot + n))
  done
  CONDTOTAL[$LAB]=$tot
  printf "  %-6s pooled pairs: %d\n" "$LAB" "$tot"
done

FLOOR=""
for v in "${CONDTOTAL[@]}"; do
  [ "$v" -gt 0 ] || continue
  if [ -z "$FLOOR" ] || [ "$v" -lt "$FLOOR" ]; then FLOOR="$v"; fi
done
[ "$DOWNSAMPLE" = "1" ] && echo "== downsample floor = $FLOOR pooled pairs (all conditions matched) =="

# ---- PASS B: (downsample) -> cload -> merge ---------------------------------
echo "== pass B: bin + merge =="
for COND in "${!LABEL[@]}"; do
  LAB="${LABEL[$COND]}"
  [ "${CONDTOTAL[$LAB]:-0}" -gt 0 ] || continue
  KEEP=1.0
  if [ "$DOWNSAMPLE" = "1" ]; then
    KEEP=$(awk -v f="$FLOOR" -v t="${CONDTOTAL[$LAB]}" 'BEGIN{printf "%.6f", f/t}')
  fi
  printf "  %-6s keep=%.3f\n" "$LAB" "$KEEP"
  for REP in $REPS; do
    PAIRS="$OUTDIR/pairs/${LAB}_rep${REP}.pairs.gz"
    [ -f "$PAIRS" ] || continue
    USE="$PAIRS"
    if [ "$DOWNSAMPLE" = "1" ] && awk -v k="$KEEP" 'BEGIN{exit !(k<0.999999)}'; then
      DS="$OUTDIR/pairs/${LAB}_rep${REP}.ds.gz"
      zcat "$PAIRS" | awk -v k="$KEEP" -v s="$SEED" 'BEGIN{srand(s)} rand()<k' \
        | bgzip > "$DS"
      USE="$DS"
    fi
    for R in "${RESARR[@]}"; do
      cooler cload pairs -c1 4 -p1 5 -c2 6 -p2 7 --assembly mm10 \
        "$CHROMSIZES:$R" "$USE" "$OUTDIR/cool/${LAB}_rep${REP}_bs_${R}.cool"
    done
  done
  for R in "${RESARR[@]}"; do
    MERGED="$OUTDIR/cool/merged_capture_${LAB}_bs_${R}.cool"
    REPLIST=( "$OUTDIR/cool/${LAB}"_rep*_bs_${R}.cool )
    [ -e "${REPLIST[0]}" ] || continue
    cooler merge "$MERGED" "${REPLIST[@]}"
    [ "$BALANCE" = "1" ] && cooler balance --cis-only --force "$MERGED"
    echo "  [done] $MERGED"
  done
done
echo "== point samples_v1 cool column at merged_capture_{COND}_bs_{res}.cool =="
