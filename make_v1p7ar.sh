#!/usr/bin/env bash
# ----------------------------------------------------------------------
# File     : make_v1p7ar.sh  (run once from the bearing repo root)
# Version  : 1.0.0
# Date     : 2026-11-04
# Authors  : Katharina E. Hayer (katharinaehayer@gmail.com) and Claude
#            (Anthropic), co-created
# ----------------------------------------------------------------------
# A21: V1P 7-track l1 run with the DN Hi-C side rebuilt from Arima-only
# libraries (DNAr, depth-matched to V1P Arima), so the 4C RC channel no longer
# compares merged conventional+Arima DN against Arima-only V1P.
#
# Writes (refuses to overwrite):
#   workflow/config/samples_v1p7ar.tsv          copy of samples_v1p_7track.tsv,
#       in DN rows only: merged_corrected_KR_DN_bs_ -> merged_corrected_KR_DNAr_bs_
#       (cool, TAD score, compartments, RC O/E 4C bw) and
#       merged_DN_loops_all_resolutions -> merged_DNAr_loops_all_resolutions
#   workflow/config/config_v1p7ar_prior1_l1.yaml copy of config_v1p7_prior1_l1.yaml,
#       new samples_sheet + outdir results_v1p7ar_prior1_l1 (nothing else changes)
# results_v1p7_prior1_l1 is NOT touched (it holds the numbers in C4 / B10).
# Checks every path in the new sheet exists ({res} checked at 10000).
# ASCII-only.
set -euo pipefail
C=workflow/config
OLD_SHEET=$C/samples_v1p_7track.tsv
OLD_CFG=$C/config_v1p7_prior1_l1.yaml
NEW_SHEET=$C/samples_v1p7ar.tsv
NEW_CFG=$C/config_v1p7ar_prior1_l1.yaml
OUTDIR=results_v1p7ar_prior1_l1

for f in $NEW_SHEET $NEW_CFG; do
  [ -e "$f" ] && { echo "exists, not overwriting: $f" >&2; exit 1; }
done
for f in $OLD_SHEET $OLD_CFG; do
  [ -f "$f" ] || { echo "missing: $f" >&2; exit 1; }
done
[ -e workflow/$OUTDIR ] && { echo "exists: workflow/$OUTDIR (pick another outdir)" >&2; exit 1; }

# 1. sample sheet
sed -e 's#merged_corrected_KR_DN_bs_#merged_corrected_KR_DNAr_bs_#g' \
    -e 's#merged_DN_loops_all_resolutions#merged_DNAr_loops_all_resolutions#g' \
    $OLD_SHEET > $NEW_SHEET
sed -i '1i # A21 sheet made from samples_v1p_7track.tsv by make_v1p7ar.sh 1.0.0 (2026-11-04): DN Hi-C inputs -> Arima-only DNAr.' $NEW_SHEET

echo "== sheet changes (old < > new)"
diff $OLD_SHEET $NEW_SHEET | grep '^[<>]' | tr '\t,' '\n\n' | grep -e '^[<>]' -e DNAr -e 'KR_DN_' -e merged_DN_ || true
NCH=$(diff $OLD_SHEET $NEW_SHEET | grep -c '^>' || true)
[ "$NCH" -ge 2 ] || { echo "ERROR: fewer than 2 changed lines ($NCH): DN paths not found under the expected names - paste $OLD_SHEET to Claude" >&2; rm -f $NEW_SHEET; exit 1; }
if grep -v '^#' $NEW_SHEET | grep -q -e 'KR_DN_bs_' -e 'merged_DN_loops'; then
  echo "ERROR: merged-DN paths remain in $NEW_SHEET" >&2; rm -f $NEW_SHEET; exit 1
fi

# 2. every path in the new sheet must exist (relative to the repo root = data_dir)
echo "== path check"
MISS=0
for p in $(grep -v '^#' $NEW_SHEET | tail -n +2 | tr '\t,' '\n\n' | grep -e '/' -e '\.bw$' -e '\.cool$' -e '\.bedpe$' -e '\.bm$' | sort -u); do
  q=${p//\{res\}/10000}
  if [ -e "$q" ]; then echo "   ok   $q"; else echo "   MISS $q"; MISS=1; fi
done
[ "$MISS" = 0 ] || { echo "ERROR: missing inputs (build the DNAr loops / 4C first). $NEW_SHEET removed." >&2; rm -f $NEW_SHEET; exit 1; }

# 3. config
sed -e "s#^samples_sheet:.*#samples_sheet: \"config/samples_v1p7ar.tsv\"#" \
    -e "s#^outdir:.*#outdir: \"$OUTDIR\"#" \
    $OLD_CFG > $NEW_CFG
sed -i "1i # config_v1p7ar_prior1_l1.yaml: generated 2026-11-04 by make_v1p7ar.sh 1.0.0 (K. E. Hayer + Claude) from config_v1p7_prior1_l1.yaml; DN Hi-C inputs Arima-only (A21)." $NEW_CFG
echo "== config"
diff $OLD_CFG $NEW_CFG || true
grep -nE '^(samples_sheet|categories|outdir|diff_stat|prior_strength|n_perms|min_signal|pvalue_min_signal):' $NEW_CFG

# 4. suggested targets = the top-level sentinels of the old run, new outdir, absolute paths (A10)
echo "== suggested targets (same sentinels as results_v1p7_prior1_l1)"
for f in workflow/results_v1p7_prior1_l1/*.done; do
  [ -e "$f" ] && echo "   $PWD/workflow/$OUTDIR/$(basename "$f")"
done
echo "Done."
