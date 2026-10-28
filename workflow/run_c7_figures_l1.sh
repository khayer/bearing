#!/usr/bin/env bash
# ----------------------------------------------------------------------
# File     : workflow/run_c7_figures_l1.sh
# Version  : 1.0.0
# Date     : 2026-10-28
# Authors  : Katharina E. Hayer (katharinaehayer@gmail.com) and Claude
#            (Anthropic), co-created
# ----------------------------------------------------------------------
#
# run_c7_figures_l1.sh -- build (or confirm) the manuscript figures that must
# come from the production l1 run, then pack the small outputs for upload.
#
#   Fig 7B   paper_figures/fig7b/Fig7B_<cmp>.{png,json}   (rule fig7b_all)
#   Fig 7E   paper_figures/fig7e/Fig7E_regional_empirical.*  (NEW rule; the
#            calibrated q_testable from table_sources/regional_empirical_*.tsv)
#   Fig S11  figures/SuppFig_S11_floor_sensitivity.{png,source.tsv}
#   Fig 5    table_sources/tcrb_domain_decomp.per_region.tsv (domain shares)
#            + figures/combined/DN_vs_DP/ (combined Hi-C + BEARING browser view)
#
# Not here: the V1P figure (waits on decision B10 / the A21 rebuild).
#
# USAGE (repository root, on the login node)
#   module purge; conda activate bearing
#   bash workflow/run_c7_figures_l1.sh            # inventory + dry run only
#   bash workflow/run_c7_figures_l1.sh --submit   # dry run, then submit to SLURM
#   bash workflow/run_c7_figures_l1.sh --pack     # tar what exists -> c7_l1_figures.tgz
#
# Options: --run <dir under workflow/, or absolute> (default results_prior1_n100_l1)
#          --config <yaml> (default workflow/config/config_<run>.generated.yaml,
#          written by run_cluster_l1.sh)
#
# SAFETY: aborts if the dry run would run any rule outside the figure rules
# listed in ALLOWED (e.g. a re-score or p-value rebuild triggered by a newer
# config or script mtime). Do not run while another job uses this checkout.
# Targets are passed as ABSOLUTE paths (relative --target paths raise
# MissingRuleException; see task A10).
#
# ASCII-only.

set -euo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
WF="$REPO/workflow"
PROFILE="$WF/profiles/slurm"
RUN="results_prior1_n100_l1"; CONFIG=""; MODE="dry"; SMK=${SNAKEMAKE:-snakemake}

if [ -n "${PYTHONPATH:-}" ]; then
  echo "NOTE: clearing PYTHONPATH (was: $PYTHONPATH)"; unset PYTHONPATH
fi
while [ $# -gt 0 ]; do
  case "$1" in
    --run) RUN="${2%/}"; shift 2 ;;
    --config) CONFIG="$2"; shift 2 ;;
    --submit) MODE="submit"; shift ;;
    --pack) MODE="pack"; shift ;;
    -h|--help) sed -n 2,40p "$0"; exit 0 ;;
    *) echo "unknown arg: $1" >&2; exit 2 ;;
  esac
done
case "$RUN" in /*) OUT="$RUN"; RUN="$(basename "$RUN")" ;; *) OUT="$WF/$RUN" ;; esac
[ -n "$CONFIG" ] || CONFIG="$WF/config/config_${RUN}.generated.yaml"
[ -d "$OUT" ] || { echo "ERROR: run folder not found: $OUT" >&2; exit 1; }
[ -f "$CONFIG" ] || { echo "ERROR: config not found: $CONFIG" >&2; exit 1; }
[ -e "$OUT/pvalue.done" ] || { echo "ERROR: $OUT/pvalue.done missing (p-values not finished)" >&2; exit 1; }

CMPS="DN_vs_DP DN_vs_EbKO DN_vs_ProB DN_vs_S3T3"
TARGETS=(
  "$OUT/paper_figures/fig7b/fig7b.done"
  "$OUT/paper_figures/fig7e/Fig7E_regional_empirical.pdf"
  "$OUT/figures/SuppFig_S11_floor_sensitivity.pdf"
  "$OUT/table_sources/tcrb_domain_decomp.per_region.tsv"
  "$OUT/figures/combined/DN_vs_DP/combined.done"
)
ALLOWED=" fig7b_null_one fig7b_all fig7e_regional_empirical regional_empirical_table \
plot_floor_sensitivity tcrb_domain_decomp bearing_hic_combined_one "

echo "== C7 inventory in $OUT"
for t in "${TARGETS[@]}"; do
  if [ -e "$t" ]; then echo "   have  ${t#$OUT/}"; else echo "   MISS  ${t#$OUT/}"; fi
done
for c in $CMPS; do
  f="$OUT/paper_figures/fig7b/Fig7B_${c}.json"
  [ -e "$f" ] && echo "   have  ${f#$OUT/}" || echo "   MISS  ${f#$OUT/}"
done

pack() {
  local tgz="$WF/c7_l1_figures.tgz"
  ( cd "$OUT" && ls -d \
      paper_figures/fig7b/Fig7B_*.png paper_figures/fig7b/Fig7B_*.json \
      paper_figures/fig7e/Fig7E_regional_empirical.* \
      figures/SuppFig_S11_floor_sensitivity.png figures/SuppFig_S11_floor_sensitivity.source.tsv \
      table_sources/tcrb_domain_decomp*.tsv \
      table_sources/regional_empirical_tcrb.tsv table_sources/regional_empirical_igh.tsv \
      figures/combined/DN_vs_DP/*.pdf 2>/dev/null > "$WF/.c7_pack.list" || true
    [ -s "$WF/.c7_pack.list" ] || { echo "nothing to pack yet" >&2; exit 1; }
    tar czhf "$tgz" -T "$WF/.c7_pack.list" )
  echo "== packed $(wc -l < "$WF/.c7_pack.list") files -> $tgz ($(du -h "$tgz" | cut -f1))"
  echo "   upload this file to the chat."
}
[ "$MODE" = "pack" ] && { pack; exit 0; }

echo; echo "== dry run"
cd "$WF"
DRYLOG="$OUT/c7_dryrun.log"
"$SMK" -s "$WF/Snakefile" --configfile "$CONFIG" --rerun-triggers mtime -n "${TARGETS[@]}" \
  > "$DRYLOG" 2>&1 || { tail -30 "$DRYLOG"; echo "ERROR: dry run failed ($DRYLOG)" >&2; exit 1; }
JOBS=$(awk '/^Job stats:/{f=1;next} f&&/^total/{f=0} f&&NF==2&&$2~/^[0-9]+$/{print $1}' "$DRYLOG" | sort -u)
awk '/^Job stats:/{f=1} f{print "   " $0} f&&/^total/{exit}' "$DRYLOG"
[ -z "$JOBS" ] && { echo "   nothing to build: every C7 target exists. Next: --pack"; exit 0; }
BAD=""
for r in $JOBS; do
  case "$ALLOWED" in *" $r "*) ;; *) BAD="$BAD $r" ;; esac
done
if [ -n "$BAD" ]; then
  echo "ERROR: the dry run would also run:$BAD" >&2
  echo "       (not figure rules - probably a newer config/script mtime). Not submitting." >&2
  echo "       Read $DRYLOG ('reason:' lines) and paste them to the chat." >&2
  exit 1
fi
echo "   OK: only figure rules would run."
[ "$MODE" = "dry" ] && { echo; echo "Dry run only. Re-run with --submit."; exit 0; }

echo; echo "== submitting C7 targets"
"$SMK" -s "$WF/Snakefile" --configfile "$CONFIG" --rerun-triggers mtime \
  --profile "$PROFILE" --keep-going "${TARGETS[@]}"
echo; echo "Done. Next: bash workflow/run_c7_figures_l1.sh --pack"
