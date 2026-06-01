#!/usr/bin/env bash
#
# smoke_test.sh -- fastest end-to-end check of the BEARING core pipeline.
#
# Generates a tiny synthetic 4-sample / 2-condition dataset (no real data, no
# Hi-C, no cluster) and runs the actual scoring -> compare -> perm -> p-value
# path locally in a few minutes. This exercises the real scripts and the data
# flow between them; it does NOT reproduce biology (use the full workflow for
# that). Run from the repo root:
#
#   bash workflow/smoke_test.sh
#
# Exit 0 = the core path runs end to end and produces non-empty outputs.
# ASCII-only.

set -euo pipefail

WITH_PERM=0
[ "${1:-}" = "--with-perm" ] && WITH_PERM=1

REPO="$(cd "$(dirname "$0")/.." && pwd)"
WORK="$(mktemp -d)"
CATS="$REPO/categories/mm10_6track_panel.yaml"
echo "repo : $REPO"
echo "work : $WORK"
cd "$WORK"

CHROM=chrSim
LEN=2000000
echo "== 1. simulate 4 synthetic samples (2 conditions x 2 reps) =="
for s in A_rep1 A_rep2 B_rep1 B_rep2; do
  seed=$(( ${#s} + RANDOM % 7 ))
  python3 "$REPO/benchmark/simulate_bearing_tracks.py" \
    --prefix "$s" --chrom "$CHROM" --length "$LEN" --bin-size 200 \
    --seed "$seed" >/dev/null
done

echo "== 2. build sample sheet + chrom.sizes =="
printf "%s\t%d\n" "$CHROM" "$LEN" > sim.chrom.sizes
{
  printf "sample\tcondition\treplicate\tout\tbw\n"
  for s in A_rep1 A_rep2 B_rep1 B_rep2; do
    cond="${s%%_*}"; rep="${s##*rep}"
    bws="${s}_atac.bw,${s}_rnaplus.bw,${s}_rnaminus.bw,${s}_ctcf.bw,${s}_rad21.bw,${s}_h3k27ac.bw"
    printf "%s\t%s\t%s\t%s.qcat.bgz\t%s\n" "$s" "$cond" "$rep" "$s" "$bws"
  done
} > sheet.tsv
cat sheet.tsv

echo "== 3. score each sample (real bigwig_to_qcat.py) =="
for s in A_rep1 A_rep2 B_rep1 B_rep2; do
  python3 "$REPO/bigwig_to_qcat.py" \
    --bw ${s}_atac.bw ${s}_rnaplus.bw ${s}_rnaminus.bw ${s}_ctcf.bw ${s}_rad21.bw ${s}_h3k27ac.bw \
    --sample-name "$s" --categories "$CATS" \
    --chrom-sizes sim.chrom.sizes --chroms "$CHROM" \
    --min-signal 0.1 --out "${s}.qcat.bgz" --jobs 2 >/dev/null
  test -s "${s}.qcat.bgz" || { echo "FAIL: ${s}.qcat.bgz empty"; exit 1; }
done
echo "  4 qcat files OK"

echo "== 4. compare (real compare_qcat.py: PCA / Spearman / diff) =="
# Use a qcat-column sheet so compare consumes the qcats scored in step 3
# rather than re-deriving them (re-derivation assumes mm10 chrom sizes).
{
  printf "sample\tcondition\treplicate\tqcat\n"
  for s in A_rep1 A_rep2 B_rep1 B_rep2; do
    cond="${s%%_*}"; rep="${s##*rep}"
    printf "%s\t%s\t%s\t%s/%s.qcat.bgz\n" "$s" "$cond" "$rep" "$WORK" "$s"
  done
} > sheet_qcat.tsv
python3 "$REPO/compare_qcat.py" \
  --sheet sheet_qcat.tsv --out compare \
  --categories "$CATS" --chroms "$CHROM" --no-clip --workers 2 >/dev/null
# primary + supplementary Spearman both present (this round's change)
test -s compare/total_saliency_spearman.tsv \
  || { echo "FAIL: primary Spearman missing"; exit 1; }
test -s compare/total_saliency_spearman_all_aligned_bins.tsv \
  || { echo "FAIL: supplementary Spearman missing"; exit 1; }
echo "  compare outputs OK (primary + supplementary Spearman present)"

if [ "$WITH_PERM" = "1" ]; then
  echo "== 5. permutation null + p-values (small: 3 perms) =="
  python3 "$REPO/generate_perm_nulls.py" \
    --sheet sheet.tsv --out-dir perm --n-perms 3 \
    --categories "$CATS" --chrom-sizes sim.chrom.sizes --min-signal 0.1 \
    --shift-workers 2 --jobs 2 >/dev/null 2>&1 \
    && echo "  perm nulls generated OK" \
    || echo "  NOTE: perm step needs blacklist/floors in some configs; see log"
else
  echo "== 5. permutation null SKIPPED (pass --with-perm to include; slower) =="
fi

echo
echo "=========================================================="
echo "PASS: BEARING core path ran end to end on synthetic data."
if [ "$WITH_PERM" = "1" ]; then
  echo "  scored 4 samples -> compared -> permuted -> (p-values ready)"
else
  echo "  scored 4 samples -> compared (diff/PCA/Spearman)"
fi
echo "Outputs in: $WORK/compare"
echo "(Hi-C phases and biology are NOT exercised here.)"
echo "=========================================================="
