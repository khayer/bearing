#!/usr/bin/env bash
#
# reproduce_all.sh -- reproduce the BEARING analyses.
#
# READ THIS FIRST
# ---------------
# The pipeline IS the Snakemake workflow. This script does not re-implement it
# and does not re-declare its parameters. Everything that determines whether a
# rerun matches the paper lives in ONE place:
#
#     workflow/config/config.yaml
#
# The previous version of this file carried its own copy of MIN_SIGNAL,
# NORMALIZE, N_PERMS and MIN_SHIFT, with comments asking which values had
# produced the published numbers. That copy drifted: it said N_PERMS=10 while
# production ran 100, so anyone following it got different numbers. Parameters
# are not duplicated here again. If you need to know what was run, read
# config.yaml; if you need to change it, change config.yaml.
#
# What this file IS: the ordered list of commands, in three phases --
#   Phase 1  the workflow (one snakemake invocation; builds nearly everything)
#   Phase 2  post-hoc analyses that are NOT yet workflow rules
#   Phase 3  assembly of the manuscript tables (not yet scripted -- see below)
#
# Step through it phase by phase. Do not run it blindly.
# SLURM: partitions dbhiq,defq. ASCII only.

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"

CONFIG="workflow/config/config.yaml"
OUT="workflow/results"          # matches `outdir:` in the config; GITIGNORED
SHEET="workflow/config/samples.tsv"
CHROM_SIZES="workflow/resources/mm10.chrom.sizes"

# Phase 2 outputs go HERE, not into $OUT.
#
# $OUT is gitignored (it holds hundreds of GB of qcats, nulls and per-bin
# tables). The Phase 2 outputs are a few hundred KB and are the direct inputs to
# paper/build_tables.py, so they are TRACKED: that is what lets a fresh clone
# rebuild the tables workbook without re-running the pipeline. If these land in
# $OUT instead, they vanish from the clone and build_tables.py reports MISSING.
#
# This path must match --sources-dir in build_tables.py. Changing one without
# the other silently breaks the loop.
SOURCES="paper/table_sources"
mkdir -p "$SOURCES/sens"

echo "repo:   $REPO"
echo "config: $CONFIG   <- the single source of truth for all parameters"
echo

# ===========================================================================
# Phase 0: environment and inputs
# ===========================================================================
#
# conda env `bearing` (python 3.12). Inputs (BigWigs, Hi-C cools) resolve
# against `data_dir:` in the config. If they are not mounted locally, stage
# them first:
#
#   python3 workflow/stage_inputs.py --sheet $SHEET \
#       --out-sheet workflow/config/samples.local.tsv --cache-dir resources/staged
#   # then run snakemake with: --config samples_sheet=workflow/config/samples.local.tsv
#
# The GTF named by `gtf:` in the config is a processed file
# (gencode.vM23.annotation_modified_overlaps_removed_sorted.gtf). See
# resources/README for how it is derived from the GENCODE release.
#
# Sanity-check the config and inputs before submitting anything:
#
#   python3 workflow/preflight.py --configfile $CONFIG

# ===========================================================================
# Phase 1: the workflow
# ===========================================================================
#
# ONE invocation builds: blacklist -> per-sample scoring -> compare ->
# permutation nulls -> p-values -> regional enrichment -> Hi-C analyses ->
# calibration -> the automated paper figures.
#
# Dry run first. It should list only what is missing; if it wants to rebuild
# the world, something upstream changed -- find out what before you submit.
#
#   snakemake --configfile $CONFIG -n -p
#
# Local:
#   snakemake --configfile $CONFIG --cores 16
#
# Cluster (each rule is submitted to SLURM):
#   bash workflow/run_cluster.sh
#
# Node-local /tmp is only 8 GB on these nodes; point Snakemake's tmpdir at
# scratch or large rules will fail:
#   export TMPDIR=/scr1/users/$USER/snakemake
#
# Figures the workflow builds automatically (with source data alongside):
#   $OUT/paper_figures/suppS1_jsd_decomposition.pdf   (+ .png, + _per_track.tsv)
#   $OUT/paper_figures/fig6_panelB_igh_decomp.pdf
#   $OUT/paper_figures/graphical_abstract_panel3.png
#   $OUT/paper_figures/qc_crosslocus_ebko_survey_rawQ.pdf
#   $OUT/paper_figures/fig3_panelB_oe_grid.png          (Hi-C configs only)
#   $OUT/paper_figures/fig4_panelA_self_contact_grid.png
#   $OUT/paper_figures/suppS7_panelB_compartment_res.png
#
# After a rerun, check nothing downstream is stale relative to the p-value layer:
#   bash staleness_audit.sh
#   bash verify_canonical_v2.sh

# ===========================================================================
# Phase 2: post-hoc analyses (NOT yet workflow rules)
# ===========================================================================
#
# These produce numbers that appear in the supplementary tables but are run by
# hand. They are therefore NOT in the Snakemake DAG: nothing tells you when
# their outputs are stale relative to $OUT/pvalue.done. Re-run them after any
# change to the p-value layer, and check with staleness_audit.sh.
#
# TODO: promote these to rules so the DAG covers them.

# --- Table S10: baseline comparison ---------------------------------------
# python3 dev/baseline_comparison.py --sheet $SHEET --repo . \
#     --chrom-sizes $CHROM_SIZES --regions regions_manuscript.tsv \
#     --min-signal 0.1 --out $SOURCES/baseline_comparison.tsv

# --- Table S11: regional-null calibration (empirical p) --------------------
# The contrast supplied MUST itself be null: use a within-condition replicate
# differential written by bearing_calibration.py (the `calibration` rule).
# --mode background draws INDEPENDENT loci elsewhere in the genome; --mode
# within samples inside the analysis locus and collapses to ~4 effective
# placements for a large region, so background is the default and the mode used
# for the published numbers.
# NOTE: dev/regional_null_calibration.py (no --mode, no --verify, n=1000 within
# the locus) is the superseded prototype with that flaw. Use the v2 script.
#
# for LOCUS_TAG in tcrb igh; do
#   python3 dev/regional_null_calibration_v2.py \
#       --stats   $OUT/calibration/DN/DN/DN_perbin.tsv.gz \
#       --regions workflow/annotations/${LOCUS_TAG}_regions_v5.bed \
#       --locus   chr6:40790000-41690000 \
#       --chrom-sizes $CHROM_SIZES \
#       --mode background --n-random 500 --p-thresh 0.05 --seed 42 \
#       --verify 5 --null-contrast --repo . \
#       --out $SOURCES/regional_null_calibration_${LOCUS_TAG}_DNrep_bg.tsv
# done
#
# OPEN ISSUE (do not treat the current numbers as final): the background loci
# are matched on scored-bin count to `ref_bins`, which is taken from the NULL
# contrast at the analysis locus (~4500 bins for the DN replicate differential),
# not from the PRODUCTION contrast being tested (1952 bins for DN-vs-Pro-B).
# The null therefore has more bins, hence more power, hence reaches smaller p,
# hence more null regions clear the real region's p -- i.e. the reported
# empirical p is CONSERVATIVE. Matching would move every value DOWN. Two calls
# are close enough for that to matter (Igh VH-proximal at 0.012, Tcrb DN-vs-DP
# Vbeta at 0.064); the rest are floor-limited at 1/n_random or far from 0.05.
# The correct matching must be specified BEFORE it is run, not chosen after
# seeing which answer it gives.

# --- Table S12: track ablation --------------------------------------------
# python3 dev/track_ablation.py --repo . \
#     --diff-qcat $OUT/compare/diff_DN_vs_DP.qcat.bgz \
#     --perm-glob "$OUT/perm/perm*/diff_comparison" \
#     --regions regions_manuscript.tsv --locus chr6:40790000-41690000 \
#     --p-thresh 0.05 --emit-stats \
#     --out $SOURCES/track_ablation_DN_vs_DP.tsv

# --- Table S13: replicate stability ---------------------------------------
# for CMP in DP EbKO ProB S3T3; do
#   srun -p dbhiq,defq --mem=64G --time=2:00:00 \
#     python3 dev/replicate_stability.py --repo . \
#       --sheet $SHEET --chrom-sizes $CHROM_SIZES \
#       --regions regions_manuscript.tsv \
#       --cond-a "DN:DN_rep1,DN_rep2" --cond-b "${CMP}:${CMP}_rep1,${CMP}_rep2" \
#       --out $SOURCES/replicate_stability_DN_vs_${CMP}.tsv
# done

# --- Table S14 + Supplementary Figure S11: differential-floor sensitivity ---
# Reuses the production diff qcats and the production permutation null; only the
# floor varies. Do NOT subsample the null: with ~5M tests, BH significance needs
# p ~1e-8 and the smallest attainable p is 1/(n_null+1).
#
# srun -p dbhiq,defq --mem=32G --cpus-per-task=12 --time=2:00:00 \
#   python3 pvminsig_sweep.py \
#       --diff-dir $OUT/compare \
#       --perm-glob "$OUT/perm/perm*/diff_comparison" \
#       --comparisons DN_vs_DP DN_vs_EbKO DN_vs_ProB DN_vs_S3T3 \
#       --floors 0.1 0.25 0.5 1.0 1.5 --threads 12 \
#       --out $SOURCES/sens/floor_sweep.tsv --verify $OUT/pvalue
#
# --verify must report [OK] at floor 0.5 for every comparison: that arm IS the
# production floor, so it must reproduce production exactly.
#
# TODO: Supplementary Figure S11 itself has no script yet; it was built ad hoc.
# See PROMPT_S11_reproducible_pipeline.md.

# --- score autocorrelation (quoted in Methods: ~967 bp independence length) --
# python3 dev/score_autocorrelation.py \
#     --score-col bearing_score --bin-size 200 --max-lag-bp 5000 \
#     --main-chroms-only < $OUT/pvalue/diff_DN_vs_DP.stats.tsv

# ===========================================================================
# Phase 3: manuscript tables and figures
# ===========================================================================
#
# Check what resolves before building; it names every missing source:
#
#   python3 paper/build_tables.py --results-dir $OUT --sources-dir $SOURCES \
#       --curated paper/tables_curated.yaml --dry-run
#
# Then build. It refuses to write a workbook containing silently empty tables;
# --allow-missing writes an explicit MISSING-SOURCE banner sheet instead.
#
#   python3 paper/build_tables.py \
#       --results-dir $OUT --sources-dir $SOURCES \
#       --config $CONFIG --sheet $SHEET \
#       --curated paper/tables_curated.yaml --repo . \
#       --out BEARING_tables.xlsx
#
# Sheets are generated, not transcribed:
#   Table S9  from $CONFIG plus the defaults parsed out of the scripts, so the
#             table cannot state a parameter the pipeline did not use;
#   Table S4  from $SHEET;
#   S1-S14    from the Phase 1/2 TSVs;
#   Table 1 and Table S5 are curated metadata and literature -- they are not
#             computable, so they live in paper/tables_curated.yaml, version
#             controlled and reviewable in a diff.
#
# The workbook carries a Provenance sheet: every source with its SHA256, mtime
# and row count, flagged STALE if it predates $OUT/pvalue.done.
#
# STILL A GAP: the figure deck assembles panels from $OUT/paper_figures/ plus
# pyGenomeTracks output by hand, and Supplementary Figure S11 has no script at
# all (see PROMPT_S11_reproducible_pipeline.md).

echo "This script is documentation, not an entrypoint. Step through the phases."
echo "Parameters live in $CONFIG -- not here."
