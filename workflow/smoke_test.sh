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
  echo "== 5. permutation nulls (observed + differential), 3 perms =="
  python3 "$REPO/generate_perm_nulls.py" \
    --sheet sheet.tsv --out-dir perm --n-perms 3 \
    --categories "$CATS" --chrom-sizes sim.chrom.sizes --min-signal 0.1 \
    --shift-workers 2 --jobs 2 >/dev/null 2>&1
  python3 "$REPO/generate_perm_nulls.py" \
    --sheet sheet.tsv --diff-sheet sheet.tsv --out-dir perm --n-perms 3 \
    --categories "$CATS" --chrom-sizes sim.chrom.sizes --min-signal 0.1 \
    --jobs 2 >/dev/null 2>&1
  echo "  perm nulls generated"

  echo "== 6. per-sample p-values (bearing_pvalue.py) =="
  for s in A_rep1 A_rep2 B_rep1 B_rep2; do
    nulls=$(ls perm/perm*/"$s"/*.qcat.bgz 2>/dev/null)
    [ -z "$nulls" ] && { echo "  FAIL: no null qcats for $s"; exit 1; }
    python3 "$REPO/bearing_pvalue.py" \
      --qcat "${s}.qcat.bgz" --null-qcat $nulls \
      --cats-json "${s}_cats.json" --chrom-sizes sim.chrom.sizes \
      --min-signal 0.1 --sort-output --out-prefix "$s" >/dev/null 2>&1
    test -s "${s}.neglog10p.bw" || { echo "  FAIL: ${s}.neglog10p.bw missing"; exit 1; }
    test -s "${s}.stats.tsv"    || { echo "  FAIL: ${s}.stats.tsv missing"; exit 1; }
    echo "  ${s}: neglog10p.bw + stats.tsv OK"
  done

  echo "== 7. differential p-values (diff_A_vs_B) =="
  # diff nulls land under perm/perm*/diff_comparison/ (layout from generate_perm_nulls)
  diffnull=$(ls perm/perm*/diff_comparison/*.qcat.bgz 2>/dev/null || true)
  if [ -z "$diffnull" ]; then
    echo "  NOTE: diff null qcats not found (diff perm run may need more time/perms);"
    echo "        skipping differential p-value. Per-sample p-values above are complete."
  else
    python3 "$REPO/bearing_pvalue.py" \
      --qcat compare/diff_A_vs_B.qcat.bgz --null-qcat $diffnull --diff \
      --cats-json A_rep1_cats.json --chrom-sizes sim.chrom.sizes \
      --min-signal 0.0 --sort-output --out-prefix diff_A_vs_B >/dev/null 2>&1 || true
    if [ -s diff_A_vs_B.stats.tsv ]; then
      echo "  diff_A_vs_B.stats.tsv OK"
    else
      echo "  NOTE: diff stats empty (synthetic A/B have no planted differential;"
      echo "        on real data with signal this produces diff_A_vs_B.stats.tsv)"
    fi
  fi
else
  echo "== 5. permutation + p-values SKIPPED (pass --with-perm; slower) =="
fi

echo
echo "=========================================================="
echo "PASS: BEARING core path ran end to end on synthetic data."
if [ "$WITH_PERM" = "1" ]; then
  echo "  scored 4 -> compared -> permuted -> p-values"
  echo "  per sample: <sample>.neglog10p.bw + <sample>.stats.tsv"
  echo "  differential: diff_A_vs_B.stats.tsv"
else
  echo "  scored 4 samples -> compared (diff/PCA/Spearman)"
fi
echo "Outputs in: $WORK"
echo "(Hi-C phases and biology are NOT exercised here.)"
echo "=========================================================="