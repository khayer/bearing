#!/usr/bin/env bash
# ----------------------------------------------------------------------
# File     : figure_render_commands_2026-11-13.sh
# Version  : 1.0.0
# Date     : 2026-11-13
# Authors  : Katharina E. Hayer (katharinaehayer@gmail.com) and Claude
#            (Anthropic), co-created
# ----------------------------------------------------------------------
# The three panels deck v26 still lacks, all from results_prior1_n100_l1:
#   S11  Supp Fig S11 floor sensitivity (rule plot_floor_sensitivity;
#        panel C on DN_vs_DP as the legend describes - decision still open)
#   S4   Supp Fig S4 per-track input tracks (rule suppS4_tracks)
#   F5   Fig 5 re-render with bearing_hic_kl_track_plot.py v11
#        (needs bearing_fig5_colorbar_2026-11-13.patch applied first)
# Nothing is submitted to SLURM; everything runs inside the interactive
# session. The script does NOT run git.
#
# USAGE (repository root):
#   squeue -u $USER                       # nothing may be using the checkout
#   git am claude/bearing_fig5_colorbar_2026-11-13.patch
#   srun -p defq,dbhiq --mem=64G -c 4 --time=6:00:00 --pty bash
#   module purge; conda activate bearing
#   bash figure_render_commands_2026-11-13.sh dry 2>&1 | tee figs_v27_dry.log
#       -> paste the log; check that only the expected rules are listed
#   bash figure_render_commands_2026-11-13.sh run 2>&1 | tee figs_v27.log
#       -> upload workflow/results_prior1_n100_l1/figs_v27.tgz + figs_v27.log
# ASCII only.
set -uo pipefail
MODE="${1:-dry}"
REPO="$(pwd)"
OUT="$REPO/workflow/results_prior1_n100_l1"
CFG="$REPO/workflow/config/config_results_prior1_n100_l1.generated.yaml"
DATA="/mnt/isilon/bassing_lab/integration_paper/bearing"
GTF="$REPO/workflow/resources/gencode.vM23.annotation_modified_overlaps_removed_sorted.gtf"
HICF="$DATA/../hic_files"
TAD="$DATA/../../projects/HiC_V31_NT_allele/data/endpoints/HiC_explorer_mm10/07tad"
FIG="$OUT/figs_v27"
S11="$OUT/figures/SuppFig_S11_floor_sensitivity.png"
S4="$OUT/paper_figures/suppS4_tracks/suppS4.done"
SNK="env -u PYTHONPATH snakemake --snakefile $REPO/workflow/Snakefile --directory $REPO/workflow --configfile $CFG --rerun-triggers mtime --cores 4"

[ -d "$OUT/pvalue" ] || { echo "ERROR: $OUT not found - run from the repo root"; exit 1; }
[ -e "$CFG" ] || { echo "ERROR: config not found: $CFG"; exit 1; }

echo "== 0. inputs"
for f in "$OUT/sens/floor_sweep.tsv" "$OUT/sens/per_bin" "$OUT/qc/Q_DN.json" \
         "$REPO/workflow/config/floor_example_bins.tsv" "$OUT/samples.qcat.tsv" \
         "$OUT/pvalue_samples.done" "$OUT/pvalue/diff_DN_vs_DP.stats.tsv" \
         "$OUT/DN_rep1_cats.json" "$HICF/merged_corrected_KR_DN_bs_10000.cool"; do
  if [ -e "$f" ]; then echo "  ok      $f"; else echo "  MISSING $f"; fi
done
echo "  plot script version line(s):"
grep -n "^Version: \|v11 changes" "$REPO/bearing_hic_kl_track_plot.py" | head -3
grep -q "v11 changes" "$REPO/bearing_hic_kl_track_plot.py" \
  || echo "  WARNING: bearing_hic_kl_track_plot.py is not v11 - git am the patch before 'run'"
grep -n "Version\|version" "$REPO/plot_floor_sensitivity.py" | head -3

if [ "$MODE" = "dry" ]; then
  echo; echo "== dry run (Job stats tables only). Expected: plot_floor_sensitivity 1, suppS4_tracks 1."
  echo "   If pvalue_sample* or other rules appear, STOP and paste this log."
  $SNK -n "$S11" "$S4" 2>&1 | sed -n '/Job stats/,/^$/p;/^total/p;/Nothing to be done/p;/Error\|Exception/p'
  $SNK -n "$S11" "$S4" 2>&1 | grep -E "^rule |reason:" | sort | uniq -c | head -40
  exit 0
fi

mkdir -p "$FIG/panels"
echo; echo "== 1. Supp Fig S11"
$SNK "$S11" || echo "  S11 FAILED"
cp -L "$OUT"/figures/SuppFig_S11_floor_sensitivity.* "$FIG/panels/" 2>/dev/null

echo; echo "== 2. Supp Fig S4"
$SNK "$S4" || echo "  S4 FAILED"
ls -la "$OUT/paper_figures/suppS4_tracks/" | head -40
for s in DN_rep1 DN_rep2 S3T3_rep1 S3T3_rep2; do
  cp -L "$OUT"/paper_figures/suppS4_tracks/*"$s"* "$FIG/panels/" 2>/dev/null
done

echo; echo "== 3. Fig 5 (v11 plot script)"
python "$REPO/bearing_hic_kl_track_plot.py" \
  --diff "$OUT/pvalue/diff_DN_vs_DP.stats.tsv" --categories "$OUT/DN_rep1_cats.json" \
  --region chr6:40400000-42400000 \
  --cool-A "$HICF/merged_corrected_KR_DN_bs_10000.cool" --cool-B "$HICF/merged_corrected_KR_DP_bs_10000.cool" \
  --insul-A "$TAD/merged_corrected_KR_DN_bs_25000_tad_score.bm" \
  --insul-B "$TAD/merged_corrected_KR_DP_bs_25000_tad_score.bm" \
  --regions-bed "$REPO/annotations/tcrb_regions_v5.bed" --gtf "$GTF" \
  --label-a DN --label-b DP --max-distance 2000000 \
  --track-pctile-low 99.0 --track-pctile-high 99.5 --threshold-scope genome \
  --out "$FIG/panels/Fig5_tcrb_domain_DN_vs_DP_l1_v11.pdf" || echo "  Fig 5 FAILED"

echo; echo "== pack"
ls -la "$FIG/panels"
( cd "$OUT" && tar czf figs_v27.tgz figs_v27 ) && ls -lh "$OUT/figs_v27.tgz"
echo "Upload: $OUT/figs_v27.tgz and figs_v27.log"
