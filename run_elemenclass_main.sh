#!/usr/bin/env bash
# Element-class regional enrichment on the MAIN 6-track panel comparisons
# (DN_vs_DP / DN_vs_EbKO / DN_vs_ProB / DN_vs_3T3), then one joint-BH family.
# ASCII only. No Snakemake, no pvalue.done gate - runs regional_enrichment.py
# directly against the existing main-panel diff tables.
#
# USAGE:  bash run_elemclass_main.sh [MAIN_OUTDIR]
#   MAIN_OUTDIR defaults to workflow/results  (NOT results_v1p - that is the
#   separate 8-track V1P panel with the 4C channels).
set -euo pipefail

MAIN_OUT="${1:-workflow/results}"
ECDIR="annotations/tcrb_element_classes"
LOCUS="chr6:40500000-41800000"
PTHRESH="0.05"
DIFFDIR="$MAIN_OUT/pvalue"
OUTDIR="$MAIN_OUT/regional"
mkdir -p "$OUTDIR"

# --- auto-detect comparisons from the diff tables actually present ---
mapfile -t DIFFS < <(ls "$DIFFDIR"/diff_DN_vs_*.stats.tsv 2>/dev/null | grep -v "diff_DN_vs_V1P" || true)
if [ "${#DIFFS[@]}" -eq 0 ]; then
    echo "ERROR: no diff_DN_vs_*.stats.tsv in $DIFFDIR"
    echo "       point MAIN_OUTDIR at the main 6-track outdir (see: grep outdir workflow/config/config.yaml)"
    exit 1
fi
echo "# comparisons found in $DIFFDIR:"
for d in "${DIFFS[@]}"; do echo "   $(basename "$d")"; done

# NOTE: if 'python3 regional_enrichment.py batch --help' does NOT list
# --region-assign, delete that flag from the two calls below.
CLASSES=(CBE promoter RSS Vbeta_functional Vbeta_pseudogene lncRNA DJC_recomb_center)

TSVS=()
for d in "${DIFFS[@]}"; do
    b="$(basename "$d")"; cmp="${b#diff_}"; cmp="${cmp%.stats.tsv}"
    for cls in "${CLASSES[@]}"; do
	outf="$OUTDIR/elemclass_${cls}_${cmp}.tsv"
	echo "# enrich  $cls  x  $cmp"
	python3 regional_enrichment.py batch \
		--diff-table "$d" \
		--regions "$ECDIR/${cls}.bed" --region-assign overlap \
		--locus "$LOCUS" --p-thresh "$PTHRESH" --bh-by none \
		--out "$outf"
	TSVS+=("$outf")
    done
done

# --- one joint-BH family across ALL classes x comparisons ---
CMPS="$(for d in "${DIFFS[@]}"; do b="$(basename "$d")"; c="${b#diff_}"; echo "${c%.stats.tsv}"; done | tr '\n
