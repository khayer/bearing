#!/usr/bin/env bash
# ----------------------------------------------------------------------
# File     : dev/run_contrast_matched_regional_null.sh
# Version  : 1.0.0
# Date     : 2026-10-09
# Authors  : Katharina E. Hayer (katharinaehayer@gmail.com) and Claude
#            (Anthropic), co-created
# ----------------------------------------------------------------------
# Contrast-matched empirical null for the regional test (Table S11 at
# prior_strength 1).
#
# WHY: the existing Table S11 null scores size-matched background regions with
# the DN rep1-vs-rep2 differential (calibration/DN/DN/DN_perbin.tsv.gz). That
# file has every bin (no |d| >= 0.5 testing floor) and almost no p < 0.05 bins,
# while the real calls run on the floor-filtered cross-condition tables
# (118-1002 tested bins per locus at alpha 1). Most null regions therefore get
# k = 0 and p_combined = 1, and the "empirical p" collapses to the fraction of
# null regions with any hit (e.g. DN-DP Eb/CBE3: nominal p 0.23 -> "empirical"
# 0.008). This mismatch was already flagged [TO FINALIZE] in Table S11 v25.
#
# HERE: the null is drawn from the SAME contrast's stats.tsv (same bin set,
# same testing floor, same g_dir): background loci of the same span and
# matched tested-bin count (+/-25%) elsewhere in the genome, antigen-receptor
# loci excluded, with a region of the same tested-bin size placed at random.
# The empirical p is the fraction of those regions at least as extreme as the
# real call. Background loci can carry real biology, so this is a specificity
# null ("more extreme than a random matched genomic locus in this contrast"),
# conservative relative to a pure type-I null. No --null-contrast flag.
#
# Run from workflow/ on a compute node:
#   bash ../dev/run_contrast_matched_regional_null.sh results_prior1_n100
# ASCII-only.

set -euo pipefail
OUT="${1:-results_prior1_n100}"
REPO="$(cd "$(dirname "$0")/.." && pwd)"
CS="$REPO/workflow/resources/mm10.chrom.sizes"
NRAND="${NRAND:-500}"
DEST="$OUT/table_sources/regional_null_contrast_matched"
mkdir -p "$DEST"

declare -A BED=( [tcrb]="$REPO/annotations/tcrb_regions_v5.bed"
                 [igh]="$REPO/annotations/igh_regions_v4.bed" )
declare -A WIN=( [tcrb]="chr6:40790000-41690000"
                 [igh]="chr12:113000000-116100000" )

for L in tcrb igh; do
  for B in EbKO DP ProB S3T3; do
    CMP="DN_vs_$B"
    echo "== $L $CMP"
    python "$REPO/dev/regional_null_calibration.py" \
      --stats "$OUT/pvalue/diff_${CMP}.stats.tsv" \
      --regions "${BED[$L]}" --locus "${WIN[$L]}" --chrom-sizes "$CS" \
      --mode background --n-random "$NRAND" --p-thresh 0.05 --seed 42 \
      --verify 5 --repo "$REPO" \
      --out "$DEST/${L}_${CMP}.tsv" 2>&1 | grep -v "^!\|WARNING: --null-contrast\|inside the locus, the rejection"
  done
done
echo "done -> $DEST"
