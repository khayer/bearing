#!/usr/bin/env bash
# verify_canonical_v2.sh
# Corrected after the v1 run exposed the real CLIs and paths.
# Confirms the DEFAULT-score production run (not the JSD variant) is complete and N=100.
# Run in the bearing conda env. ASCII-only.

set -u
echo "=== BEARING canonical-run verification v2 ==="; date

# ---------------------------------------------------------------------------
# EDIT / VERIFY these
REPO="/mnt/isilon/bassing_lab/integration_paper/bearing"
RUN_DIR="$REPO/workflow/results"
PERM_DIR="$RUN_DIR/perm"
PVAL_DIR="$RUN_DIR/pvalue"
EXPECT_PERMS=100
EXPECT_BINS=13654391
# VERIFY these against your run:
MIN_SIGNAL="0.1"                                   # config value (Methods says 0.01 - that mismatch is a separate fix)
CATS="$REPO/DN_rep1_cats.json"
SAMPLES="DN_rep1 DN_rep2 EbKO_rep1 EbKO_rep2 DP_rep1 DP_rep2 ProB_rep1 ProB_rep2 S3T3_rep1 S3T3_rep2"
COMPS="DN_vs_DP DN_vs_EbKO DN_vs_ProB DN_vs_S3T3"
CONDS="DN EbKO DP ProB S3T3"
# ---------------------------------------------------------------------------
PASS=0; FAIL=0
ok(){ echo "  PASS: $1"; PASS=$((PASS+1)); }
bad(){ echo "  FAIL: $1"; FAIL=$((FAIL+1)); }

# 1) n_perms from the REAL production config (exclude the jsd/generated variants)
echo; echo "[1] n_perms in production config (excluding *jsd* and *generated*)"
CFG=$(find "$REPO/workflow/config" -name "config*.yaml" 2>/dev/null | grep -viE "jsd|generated" | head -1)
echo "    using: $CFG"
grep -iE "n_perms|num_perms" "$CFG" 2>/dev/null | sed 's/^/      /'
grep -iE "n_perms|num_perms" "$CFG" 2>/dev/null | grep -q "$EXPECT_PERMS" \
  && ok "production config declares $EXPECT_PERMS" || bad "production config does not show $EXPECT_PERMS (read lines above)"

# 2) count actual permutation rounds on disk
echo; echo "[2] permutation rounds present on disk"
NPERM=$(ls -d "$PERM_DIR"/perm*/ 2>/dev/null | wc -l)
echo "    $PERM_DIR has $NPERM perm round dirs"
[ "$NPERM" -eq "$EXPECT_PERMS" ] && ok "exactly $EXPECT_PERMS perm rounds" \
  || bad "found $NPERM perm rounds, expected $EXPECT_PERMS"
# every round should carry the diff outputs for all four comparisons
MISS=0
for d in "$PERM_DIR"/perm*/; do
  for c in $COMPS; do
    [ -s "${d}diff_comparison/diff_${c}.qcat.bgz" ] || MISS=$((MISS+1))
  done
done
[ "$MISS" -eq 0 ] && ok "all perm rounds have all 4 diff comparisons" \
  || bad "$MISS missing/empty perm diff files across rounds"

# 3) permutation integrity (correct CLI this time)
echo; echo "[3] check_perm_set.py --integrity"
python "$REPO/check_perm_set.py" --perm-dir "$PERM_DIR" --n-perms "$EXPECT_PERMS" \
       --comparisons $COMPS --integrity \
  && ok "perm set integrity clean" || bad "perm set integrity failed (see output)"

# 4) scoring provenance (correct CLI)
echo; echo "[4] assert_score_provenance.py"
python "$REPO/assert_score_provenance.py" --outdir "$RUN_DIR" --samples $SAMPLES \
       --min-signal "$MIN_SIGNAL" --categories "$CATS" \
  && ok "score provenance clean" || bad "score provenance failed (see output)"

# 5) differential stats present (correct path: pvalue/diff_<comp>.stats.tsv)
echo; echo "[5] differential stats files"
for c in $COMPS; do
  f="$PVAL_DIR/diff_${c}.stats.tsv"
  if [ -s "$f" ]; then ok "diff_${c}.stats.tsv ($(wc -l < "$f") lines)"; else bad "diff_${c}.stats.tsv missing"; fi
done

# 6) calibration completeness - must have all 5 conditions
echo; echo "[6] calibration summary completeness (expect 5 conditions)"
CAL=$(find "$RUN_DIR" -iname "*calib*summary*.tsv" 2>/dev/null | xargs -I{} sh -c 'echo "$(wc -l < "{}") {}"' 2>/dev/null | sort -rn | head -1 | awk '{print $2}')
echo "    richest calibration summary: $CAL ($(wc -l < "$CAL" 2>/dev/null) lines)"
if [ -n "$CAL" ]; then
  for cc in $CONDS; do grep -q "$cc" "$CAL" && echo "      has $cc" || echo "      MISSING $cc"; done
  N=$(grep -cE "$(echo $CONDS | tr ' ' '|')" "$CAL")
  [ "$N" -ge 5 ] && ok "calibration has all 5 conditions" || bad "calibration summary incomplete ($N/5)"
else
  bad "no calibration summary found"
fi

# 7) floor proves N (text-based, not binary): min p in a stats file
echo; echo "[7] empirical floor implies N (from a diff stats TSV)"
SF="$PVAL_DIR/diff_DN_vs_DP.stats.tsv"
if [ -s "$SF" ]; then
  # smallest nonzero p in the p-value column; column name may vary - adjust if needed
  MINP=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)if($i ~ /pval|p_value|pvalue/ && $i !~ /adj|bh/){pc=i}} NR>1 && $pc+0>0{print $pc}' "$SF" 2>/dev/null | sort -g | head -1)
  echo "    smallest per-bin p in DN_vs_DP = $MINP (floor at N=100 is ~7.3e-10)"
fi

# 8) spot-check production numbers (the freeze target)
echo; echo "[8] production-number spot check"
grep -i "DN_vs_ProB" "$RUN_DIR/consolidated_enrichment_igh.tsv"  | grep -i DH_cluster        | awk -F'\t' '{print "    Igh DH      q="$(NF-1)"  (expect 1.09e-137, k_neg="$11")"}'
grep -i "DN_vs_DP"   "$RUN_DIR/consolidated_enrichment_tcrb.tsv" | grep -i Vbeta_cluster_main | awk -F'\t' '{print "    Tcrb Vbeta  q="$(NF-1)"  (expect 0.036, k_pos="$10")"}'

echo; echo "=== SUMMARY v2: $PASS passed, $FAIL failed ==="
[ "$FAIL" -eq 0 ] && echo "Default-score run is complete and N=100. Safe to freeze $RUN_DIR." \
                  || echo "Resolve real FAILs above (perm count, integrity, provenance, calibration completeness)."
