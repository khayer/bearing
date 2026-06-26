#!/usr/bin/env bash
#
# run_mcc_grid_benchmark.sh -- answer the proposal's adaptive-vs-uniform question
# on the MCC, now that the V1 BES is done. Runs benchmark_grid_vs_standard.py
# once per GRID; each run reports adaptive(grid) / uniform(grid mean width) /
# standard(--std-binsize), with rho_contact + feat_conc, for both pairs in
# comparisons_v1.tsv (DN-vs-RCTKO, DN-vs-V1P).
#
# Grids tested:
#   native   seg_v1_chr6_real.bed       ~350 bp, feature-anchored (the loophole:
#            1D-projected contact change can use bins finer than the 2D floor)
#   floored  superbins_v1r_1000.bed     ~1 kb, depth-matched (flooring inverted
#            the anchoring -> confirmatory)
#
# Read: in the NATIVE run, does adaptive beat its own uniform? That is the
# proposal answer. If yes, feature-anchoring carries on the 1D projection
# despite the 2D depth floor; if adaptive ~ uniform, it does not.
#
#   bash run_mcc_grid_benchmark.sh
# Run from the data dir. Requires pysam/numpy/scipy and the V1 BES present.
# ASCII only.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"

BAMDIR="/mnt/isilon/bassing_lab/projects/capture_hic/Brittney_V1/03aln"
BENCH="$HERE/benchmark_grid_vs_standard.py"
ANNOT="annotations/AgRgenes_mm10_s.bed"
COMPARISONS="comparisons_v1.tsv"
OUT="mcc_grid_benchmark"
STD_BINSIZE=3000
MIN_DIST=10000
MAX_DIST=1000000
MIN_SHIFT=50000
# label (matches comparisons condA/condB)  ->  BAM genotype token
declare -A TOK=( [DN]=R1KO [V1P]=V1PKO [RCTKO]=RCTKO )
# grids to test:  name  ->  bed
declare -A GRID=( [native]="seg_v1_chr6_real.bed" [floored]="superbins_v1r_1000.bed" )

mkdir -p "$OUT"

# 1. pair counts (first-in-pair primaries = templates), labelled to match BES
PC="$OUT/pair_counts_v1.tsv"
: > "$PC"
for LAB in "${!TOK[@]}"; do
  for REP in 1 2; do
    BAM="$BAMDIR/raw_ArimaHTS_S065_${TOK[$LAB]}_rep${REP}_bs_250_master_valid_pairs.bam"
    [ -f "$BAM" ] || { echo "  [skip] missing $BAM" >&2; continue; }
    N=$(samtools view -c -f 0x40 -F 0x900 "$BAM")
    printf "%s\t%s\t%d\n" "$LAB" "$(basename "$BAM")" "$N" >> "$PC"
  done
done
echo "== pair counts =="; cat "$PC"

# 2. one benchmark run per grid (both comparisons inside)
for NAME in "${!GRID[@]}"; do
  BED="${GRID[$NAME]}"
  [ -f "$BED" ] || { echo "  [skip grid] missing $BED" >&2; continue; }
  echo "== grid: $NAME ($BED) =="
  python "$BENCH" --superbins "$BED" --counts "$PC" --bam-dir "$BAMDIR" \
    --comparisons "$COMPARISONS" --annot "$ANNOT" \
    --std-binsize "$STD_BINSIZE" --min-distance "$MIN_DIST" \
    --max-distance "$MAX_DIST" --min-shift "$MIN_SHIFT" \
    --out "$OUT/grid_benchmark_${NAME}.tsv"
done

# 3. combined readout: adaptive vs uniform per grid/pair
echo; echo "== summary: rho_contact (and feat_conc) by grid / pair / placement =="
python3 - "$OUT" << 'PY'
import csv, glob, os, sys
outdir = sys.argv[1]
rows = []
for f in sorted(glob.glob(os.path.join(outdir, "grid_benchmark_*.tsv"))):
    grid = os.path.basename(f).replace("grid_benchmark_", "").replace(".tsv", "")
    with open(f) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            r["_gridrun"] = grid
            rows.append(r)
if not rows:
    print("  (no benchmark output found)"); sys.exit()
print("  %-8s %-12s %-9s %8s %9s %9s" %
      ("gridrun", "pair", "placement", "rho", "p_emp", "featconc"))
for r in rows:
    print("  %-8s %-12s %-9s %8s %9s %9s" % (
        r.get("_gridrun", "?"), r.get("comparison", "?"), r.get("grid", "?"),
        r.get("rho_contact", "NA"), r.get("p_emp_contact", "NA"),
        r.get("feat_concentration", "NA")))
print("\n  READ: in the 'native' gridrun, compare placement=adaptive vs uniform.")
print("  adaptive rho > uniform rho (with feat_conc>1) => feature-anchoring helps the MCC.")
PY

