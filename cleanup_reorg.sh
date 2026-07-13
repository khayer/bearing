#!/usr/bin/env bash
#
# BEARING repository cleanup -- Sections C, D, E, F of CLEANUP_MANIFEST.md
# -----------------------------------------------------------------------------
# Safe by default: this runs in DRY-RUN mode and only PRINTS the git commands.
# Review the output, then run for real with:
#
#     DRY_RUN=0 bash cleanup_reorg.sh
#
# It operates through `git mv` / `git rm`, so every change is staged and fully
# reversible with `git reset --hard` (or per-file `git checkout --`) until you
# commit. It creates a dedicated branch first and refuses to run on a dirty
# tree. Run it from the repository root.
#
# Prerequisites verified before this script was generated:
#   - The import guard (tests/test_import_guard.py) is committed and green.
#   - Every deletion/move below was checked against Snakefiles, shell wrappers,
#     configs, and Python imports -- not just the docs.
# -----------------------------------------------------------------------------

set -euo pipefail
DRY_RUN="${DRY_RUN:-1}"
BRANCH="repo-cleanup-reorg"

run() {
  echo "  + $*"
  if [ "$DRY_RUN" = "0" ]; then "$@"; fi
}

# Delete a tracked file; skip cleanly if it is already gone (idempotent re-runs).
gone_rm() {
  for f in "$@"; do
    if git ls-files --error-unmatch "$f" >/dev/null 2>&1; then
      run git rm -q -- "$f"
    else
      echo "  skip (not tracked, already removed?): $f"
    fi
  done
}

# Move a tracked root file into dev/ ; skip cleanly if already gone.
mv_dev() {
  for f in "$@"; do
    if git ls-files --error-unmatch "$f" >/dev/null 2>&1; then
      run git mv -- "$f" "dev/$f"
    else
      echo "  skip (not tracked, already moved?): $f"
    fi
  done
}

# --- Safety checks -----------------------------------------------------------
if [ ! -d .git ] || [ ! -f pyproject.toml ]; then
  echo "ERROR: run this from the bearing repository root (no .git/pyproject.toml here)."
  exit 1
fi
if [ "$DRY_RUN" = "0" ]; then
  if [ -n "$(git status --porcelain)" ]; then
    echo "ERROR: working tree is not clean. Commit or stash first."
    exit 1
  fi
  echo "==> creating/using branch: $BRANCH"
  git checkout -B "$BRANCH"
fi

echo
echo "###############################################################"
echo "#  DRY_RUN=$DRY_RUN   (set DRY_RUN=0 to execute)"
echo "###############################################################"

# =============================================================================
# SECTION D -- delete OBSOLETE divergent root copies (keep mcc/scripts/)
# Confirmed: mcc/Snakefile invokes these via {S}=scriptdir=mcc/scripts.
# =============================================================================
echo
echo "== Section D: divergent duplicates, delete obsolete root copy =="
gone_rm partial_corr.py
gone_rm plot_1d_vs_3d.py

# ---------------------------------------------------------------------------
# DECISION NEEDED -- build_adaptive_segmentation.py is NOT deleted here.
# Unlike the two above, workflow/Snakefile line ~284 invokes the ROOT copy via
#   {BIN}=config["bearing_dir"] (the repo root), and the two copies DIVERGE.
# So the root copy is the one the main workflow runs -- deleting it would break
# rule build_adaptive_segmentation. Decide one of:
#   (a) root copy is canonical  -> delete mcc/scripts/build_adaptive_segmentation.py instead
#   (b) mcc copy is canonical   -> point workflow/Snakefile at mcc/scripts/ then delete root
# Diff them first:  git diff --no-index build_adaptive_segmentation.py mcc/scripts/build_adaptive_segmentation.py
# Then uncomment ONE of:
#   gone_rm mcc/scripts/build_adaptive_segmentation.py      # option (a)
#   gone_rm build_adaptive_segmentation.py                  # option (b), after editing the Snakefile
# ---------------------------------------------------------------------------

# =============================================================================
# SECTION C -- delete IDENTICAL duplicate root copies (keep canonical location)
# =============================================================================
echo
echo "== Section C: identical duplicates, delete root copy =="
# standalone identical copies (canonical copy lives in the subdir shown)
gone_rm build_superbin_grid.py          # keep mcc/scripts/build_superbin_grid.py
gone_rm superbin_feature_overlap.py     # keep mcc/scripts/superbin_feature_overlap.py
gone_rm make_mcc_bearing_ini.py         # keep mcc/scripts/make_mcc_bearing_ini.py
gone_rm minsig_depth_check.py           # keep workflow/minsig_depth_check.py

echo
echo "== Section C (coupled MCC wrapper+helper pairs) =="
# Root .sh wrapper calls the sibling root .py; both are byte-identical to the
# mcc/scripts/ pair, and mcc/Snakefile uses the mcc/scripts/ copies. Delete
# each root pair TOGETHER so no wrapper is left pointing at a deleted helper.
gone_rm validpairs_to_pairs.py build_capture_cools.sh           # keep mcc/scripts/ pair
gone_rm benchmark_grid_vs_standard.py run_mcc_grid_benchmark.sh # keep mcc/scripts/ pair

# =============================================================================
# SECTION E -- deduplicate annotation BEDs
# workflow config + Snakefile read root annotations/ (joined with BIN).
# The only refs to workflow/annotations/ are inside Section-F scratch scripts
# being retired below, so the workflow/annotations/ copy is redundant.
# =============================================================================
echo
echo "== Section E: drop duplicate workflow/annotations/ (keep root annotations/) =="
if git ls-files --error-unmatch workflow/annotations >/dev/null 2>&1; then
  run git rm -q -r -- workflow/annotations
else
  echo "  skip: workflow/annotations already removed"
fi

# =============================================================================
# SECTION F -- move unreferenced dev/scratch scripts to dev/
# EXCLUDED from this list (kept at root on purpose):
#   check_perm_set.py     -- called by verify_canonical.sh / verify_canonical_v2.sh
#   fetch_geo_bigwigs.py  -- reviewer-useful; DOCUMENT instead of hide
#   convert_categories.py -- user-facing utility; DOCUMENT instead of hide
#   score_provenance.py, bearing_hic_plot_pval_overlay.py -- load-bearing (Section B)
# =============================================================================
echo
echo "== Section F: move scratch/dev scripts to dev/ =="
run mkdir -p dev
mv_dev \
  baseline_comparison.py \
  bearing_diff_concordance.py \
  bearing_specificity_profile.py \
  bearing_tad_decomposition.py \
  bin_hic_superbins.py \
  cbe_across_runs.py \
  compare_regimes.py \
  diagnose_flat_pvalue.py \
  extract_track_diff_bed.py \
  feature_changes_across_runs.py \
  hic_probe.py \
  make_bearing_bins_bed.py \
  plot_bearing_stats_summary.py \
  plot_hic_superbins.py \
  plot_sig_density.py \
  probe_capture_depth.py \
  prototype_adaptive_binning.py \
  regen_fig7_panelB.py \
  region_match_ebko.py \
  region_match_ebko2.py \
  region_match_ebko3.py \
  regional_null_calibration.py \
  regional_null_calibration_v2.py \
  replicate_stability.py \
  robustness_summary.py \
  run_bes_hic_matrix.py \
  score_autocorrelation.py \
  tail_diag.py \
  track_ablation.py

if [ "$DRY_RUN" = "0" ]; then
  cat > dev/README.md <<'EOF'
# dev/ -- development and scratch scripts

These scripts are retained for provenance but are NOT part of the documented
BEARING pipeline. They are unreferenced by the docs, workflow, paper figure
map, and tests. Some contain hardcoded paths or were one-off diagnostics.
Nothing in the supported pipeline imports them.
EOF
  git add dev/README.md
fi

# =============================================================================
# Verify
# =============================================================================
echo
if [ "$DRY_RUN" = "0" ]; then
  echo "== reinstall + run tests (import guard must stay green) =="
  pip install -e . -q
  python -m pytest -q
  echo
  echo "Done. Review with:  git status   and   git diff --staged --stat"
  echo "Commit when satisfied:  git commit -m 'Clean up repo: dedupe + retire scratch scripts'"
else
  echo "DRY RUN complete -- nothing changed. Re-run with DRY_RUN=0 to apply."
fi

# =============================================================================
# FOLLOW-UPS this script deliberately does NOT do (need a human):
#   1. build_adaptive_segmentation.py -- resolve the divergent-copy decision above.
#   2. Stale comment: regions_manuscript.tsv line ~8 documents
#      `python baseline_comparison.py ...` which now lives in dev/. Update or drop.
#   3. Section G doc fixes: broken docs/ROADMAP.md links (code_todo.md,
#      manuscript_todo.md); plotting entry-point mismatch (README vs
#      working_examples); unpinned requirements.txt; sync_manifest.sha256.
#   4. Add QUICKSTART.md to root and link it from the top of README.
#   5. DOCUMENT fetch_geo_bigwigs.py and convert_categories.py (kept at root).
# =============================================================================
