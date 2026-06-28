#!/usr/bin/env bash
# verify_canonical.sh
# Confirm a BEARING results folder is complete and was run at n_perms = 100.
# Run on the cluster in the bearing conda env. ASCII-only on purpose.
# Usage: edit RUN_DIR (and REPO if scripts are elsewhere), then: bash verify_canonical.sh

set -u
echo "=== BEARING canonical-run verification ==="
date

# ---------------------------------------------------------------------------
# EDIT THESE TWO PATHS
RUN_DIR="/mnt/isilon/bassing_lab/integration_paper/bearing/RESULTS_FROZEN"   # the folder you want to freeze
REPO="/mnt/isilon/bassing_lab/integration_paper/bearing"                      # repo with the scripts
# ---------------------------------------------------------------------------

EXPECT_PERMS=100
EXPECT_BINS=13654391           # genome-wide bin count after blacklist/low-signal mask
EXPECT_FLOOR=7.32e-10          # = 1 / (100 * 13,654,391 + 1)
PASS=0; FAIL=0
ok(){ echo "  PASS: $1"; PASS=$((PASS+1)); }
bad(){ echo "  FAIL: $1"; FAIL=$((FAIL+1)); }

echo; echo "RUN_DIR = $RUN_DIR"
[ -d "$RUN_DIR" ] || { echo "RUN_DIR does not exist - fix the path."; exit 2; }

# 1) n_perms from config -----------------------------------------------------
echo; echo "[1] n_perms declared in config"
CFG=$(find "$RUN_DIR" "$REPO" -maxdepth 3 -name "config*.yaml" 2>/dev/null | head -1)
if [ -n "$CFG" ]; then
  NP=$(grep -iE "n_perms|num_perms|permutations" "$CFG" | head -1)
  echo "    $CFG -> $NP"
  echo "$NP" | grep -q "$EXPECT_PERMS" && ok "config declares $EXPECT_PERMS perms" || bad "config does not say $EXPECT_PERMS (check above)"
else
  bad "no config*.yaml found under RUN_DIR or REPO"
fi

# 2) n_perms proven from the empirical floor ---------------------------------
# floor p = 1 / (N_null + 1); N_null = N_perms * n_bins. Recover N_perms from any floor.
echo; echo "[2] n_perms proven from the permutation floor"
FLOORFILE=$(grep -rilE "floor|min.*achievable|7\.3.?e-10" "$RUN_DIR" 2>/dev/null | grep -iE "calib|summary|pvalue|diff" | head -1)
if [ -n "$FLOORFILE" ]; then
  FLOOR=$(grep -oiE "[0-9.]+e-1[0-9]" "$FLOORFILE" | head -1)
  echo "    found floor token '$FLOOR' in $FLOORFILE"
  if [ -n "$FLOOR" ]; then
    NREC=$(python3 -c "f=float('$FLOOR'); print(round((1.0/f)/$EXPECT_BINS))" 2>/dev/null)
    echo "    recovered N_perms = round( (1/$FLOOR) / $EXPECT_BINS ) = $NREC"
    [ "$NREC" = "$EXPECT_PERMS" ] && ok "floor implies N=$EXPECT_PERMS" || bad "floor implies N=$NREC, not $EXPECT_PERMS"
  fi
else
  echo "    (no floor token auto-found; check DN_calibration.pdf footer reads p=7.32e-10)"
fi

# 3) permutation set integrity (NOT just file size) --------------------------
echo; echo "[3] permutation set integrity (full-stream, not size)"
if [ -f "$REPO/check_perm_set.py" ]; then
  python3 "$REPO/check_perm_set.py" --integrity --expect-perms "$EXPECT_PERMS" "$RUN_DIR" \
    && ok "check_perm_set.py --integrity clean" \
    || bad "check_perm_set.py --integrity reported a problem (truncated/empty perm, or wrong count)"
else
  echo "    check_perm_set.py not found; manual fallback:"
  NB=$(find "$RUN_DIR" -name "*perm*.bgz" 2>/dev/null | wc -l)
  EMPTY=$(find "$RUN_DIR" -name "*perm*.bgz" -size -1k 2>/dev/null | wc -l)
  echo "    perm .bgz files: $NB ; suspiciously small (<1k): $EMPTY"
  echo "    NOTE size is not enough - an 84MB perm file can still hold zero usable bins."
fi

# 4) scoring provenance sidecars match --------------------------------------
echo; echo "[4] scoring provenance (.bgz.sig sidecars)"
if [ -f "$REPO/assert_score_provenance.py" ]; then
  python3 "$REPO/assert_score_provenance.py" "$RUN_DIR" \
    && ok "assert_score_provenance.py clean (signatures match)" \
    || bad "provenance mismatch - scores were not produced by the stamped pipeline"
else
  SIG=$(find "$RUN_DIR" -name "*.bgz.sig" 2>/dev/null | wc -l)
  BGZ=$(find "$RUN_DIR" -name "*.bgz" 2>/dev/null | wc -l)
  echo "    .bgz=$BGZ  .bgz.sig=$SIG (every scored .bgz should have a sidecar)"
fi

# 5) required deliverables present and non-empty -----------------------------
echo; echo "[5] required output files present and non-empty"
need(){ f=$(find "$RUN_DIR" -name "$1" 2>/dev/null | head -1); if [ -n "$f" ] && [ -s "$f" ]; then ok "$1 ($(wc -l < "$f") lines)"; else bad "$1 missing or empty"; fi; }
need "consolidated_enrichment_tcrb.tsv"
need "consolidated_enrichment_igh.tsv"
need "calibration_summary.tsv"
for c in DP EbKO ProB S3T3; do
  f=$(find "$RUN_DIR" -name "diff_DN_vs_${c}*.stats" 2>/dev/null | head -1)
  if [ -n "$f" ] && [ -s "$f" ]; then ok "diff_DN_vs_${c} stats"; else bad "diff_DN_vs_${c} stats missing"; fi
done

# 6) sanity-check the regional TSVs are the production numbers we are freezing
echo; echo "[6] spot-check that regional TSVs carry the production numbers"
IGH=$(find "$RUN_DIR" -name "consolidated_enrichment_igh.tsv" 2>/dev/null | head -1)
if [ -n "$IGH" ]; then
  DH=$(grep -i "DN_vs_ProB" "$IGH" | grep -i "DH_cluster" | awk -F'\t' '{print $(NF-1)}')
  echo "    Igh DN_vs_ProB DH_cluster q_combined = $DH  (expect ~1.09e-137, 187/187)"
fi
TCRB=$(find "$RUN_DIR" -name "consolidated_enrichment_tcrb.tsv" 2>/dev/null | head -1)
if [ -n "$TCRB" ]; then
  VB=$(grep -i "DN_vs_DP" "$TCRB" | grep -i "Vbeta_cluster_main" | awk -F'\t' '{print $(NF-1)}')
  echo "    Tcrb DN_vs_DP Vbeta q_combined      = $VB   (expect ~0.036, 24/24)"
fi

echo; echo "=== SUMMARY: $PASS passed, $FAIL failed ==="
[ "$FAIL" -eq 0 ] && echo "Folder looks complete and consistent with N=100. Safe to freeze." \
                  || echo "Resolve the FAIL items before freezing."
