#!/usr/bin/env bash
# ----------------------------------------------------------------------
# File     : figure_update_commands_2026-11-06.sh
# Version  : 1.0.0
# Date     : 2026-11-06
# Authors  : Katharina E. Hayer (katharinaehayer@gmail.com) and Claude
#            (Anthropic), co-created
# ----------------------------------------------------------------------
# Bring every manuscript figure panel and every remaining [L1 NEEDED] value in
# line with workflow/results_prior1_n100_l1 (prior 1, 100 perms, diff_stat l1).
#
# Almost everything already EXISTS on disk: run_cluster_l1.sh built the full
# Snakemake `all` target in results_prior1_n100_l1 (tables sources, regional,
# element classes, Hi-C overlays, combined browser figures, compare/ linked from
# results_prior1_n100). So this script mostly COLLECTS; the only new compute is
#   - Fig 5 (one bearing_hic_kl_track_plot.py call, ~10-30 min, needs ~60 GB)
#   - Fig 7E re-render with fig7e_regional_empirical.py v1.1 (seconds)
#   - Supp Fig S12 (make_figS12_elemclass.py, seconds)
#   - Supp Fig S7B 500 kb row (compartment_analysis.py, minutes)
#   - section 2 value extraction (awk/python over existing tables, minutes)
# Nothing here writes into the results tables; outputs go to $FIG.
# It does NOT run git. Do not `git pull` / `git am` while a cluster job uses
# the checkout (the v1p7ar calibration). Section 4 (S11) needs a patch and is
# therefore separate: run it only after those jobs have finished.
#
# USAGE (interactive node, from the repository root):
#   srun -p defq --mem=64G -c 4 --time=4:00:00 --pty bash
#   module purge; conda activate bearing
#   cd <bearing repo>
#   cp <download>/make_figS12_elemclass.py paper/figures/    # from the project
#   bash figure_update_commands_2026-11-06.sh 2>&1 | tee figs_v25.log
#   -> upload workflow/results_prior1_n100_l1/figs_v25.tgz and figs_v25.log
# ASCII only.
set -uo pipefail

REPO="$(pwd)"
OUT="$REPO/workflow/results_prior1_n100_l1"
SRC="$REPO/workflow/results_prior1_n100"          # alpha-1 sum run (swapnull lives here)
DATA="/mnt/isilon/bassing_lab/integration_paper/bearing"
GTF="$REPO/workflow/resources/gencode.vM23.annotation_modified_overlaps_removed_sorted.gtf"
HICF="$DATA/../hic_files"
TAD="$DATA/../../projects/HiC_V31_NT_allele/data/endpoints/HiC_explorer_mm10/07tad"
COMP="$DATA/../comp"
FIG="$OUT/figs_v25"
VAL="$FIG/values_l1.txt"
CMPS="DN_vs_DP DN_vs_EbKO DN_vs_ProB DN_vs_S3T3"
mkdir -p "$FIG"/{panels,tables}
[ -d "$OUT/pvalue" ] || { echo "ERROR: $OUT not found - run from the repo root"; exit 1; }
: > "$VAL"
have() { if ls $1 >/dev/null 2>&1; then echo "  ok      $1"; return 0; else echo "  MISSING $1"; return 1; fi; }
grab() { for f in $1; do [ -e "$f" ] && cp -L "$f" "$FIG/$2/" ; done; }

echo "== 0. what the l1 run already has"
for f in "$OUT/paper_tables/BEARING_tables.xlsx" "$OUT/calibration/calibration_summary.tsv" \
         "$OUT/table_sources/regional_empirical_tcrb.tsv" "$OUT/regional/elemclass_consolidated.tsv" \
         "$OUT/hic/crosslocus.done" "$OUT/hic/bes_hic_within.done" "$OUT/figures_combined.done" \
         "$OUT/paper_figures/collect.done" "$OUT/sens/floor_sweep.tsv"; do have "$f"; done

###############################################################################
echo; echo "== 1. TABLE SOURCES (C6 values; copied, not recomputed)"
###############################################################################
grab "$OUT/table_sources/baseline_comparison_*.tsv" tables            # Table S10
grab "$OUT/table_sources/track_ablation_*.tsv" tables                 # Table S12
grab "$OUT/table_sources/regional_null_contrast_matched/ablation_*_summary.tsv" tables
grab "$OUT/table_sources/replicate_stability_*.tsv" tables            # Table S13
grab "$OUT/table_sources/regional_empirical_*.tsv" tables             # Tables S2/S11
grab "$OUT/table_sources/focal_flank/*.tsv" tables
grab "$OUT/table_sources/tcrb_domain_decomp*.tsv" tables
grab "$OUT/regional/enrich_DN_vs_*_cbe.tsv" tables                    # Table S3
grab "$OUT/regional/enrich_DN_vs_*_tcrb.tsv" tables
grab "$OUT/regional/enrich_DN_vs_*_igh.tsv" tables
grab "$OUT/regional/elemclass_*.tsv" tables                           # Table S16
grab "$OUT/regional/significant_bins_summary.tsv" tables
grab "$OUT/hic/crosslocus_p95_dcontact*" tables                       # Fig S7C / cross-locus
grab "$OUT/hic/crosslocus_sets_v2*" tables
grab "$OUT/hic/bes_hic_within_DN_vs_DP*" tables
grab "$OUT/compare/total_saliency_spearman*.tsv" tables               # Fig S2B (statistic-independent)
grab "$OUT/compare/pca_assay_loadings*.tsv" tables
grab "$OUT/compare/q_pair_jsd*.tsv" tables
grab "$OUT/compare/regions/metrics/new_wide_*" tables                 # Fig S6
grab "$OUT/compare/regions/metrics/igh_*" tables                      # Fig S5D
grab "$OUT/calibration/calibration_summary.tsv" tables
grab "$OUT/sens/floor_sweep.tsv" tables
grab "$OUT/paper_tables/BEARING_tables.xlsx" tables
ls "$OUT/compare" > "$FIG/tables/compare_listing.txt" 2>/dev/null
ls "$OUT/compare/regions" "$OUT/compare/regions/metrics" >> "$FIG/tables/compare_listing.txt" 2>/dev/null

###############################################################################
echo; echo "== 2. VALUES for the [L1 NEEDED] markers -> $VAL"
###############################################################################
{
echo "## 2a g_dir per comparison (Results regional subsection; Table S2)"
for c in $CMPS; do
  for L in tcrb igh; do
    f="$OUT/regional/enrich_${c}_${L}.tsv"
    [ -e "$f" ] && awk -F'\t' -v c="$c" -v L="$L" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}
      {print c, L, "g_dir=" $h["g_dir"]; exit}' "$f"
  done
done

echo; echo "## 2b elevation ratio + FDR bins (Fig 5 legend, Results): per_locus_summary on tcrb_domain.bed"
python "$REPO/hic/per_locus_summary.py" --loci-bed "$REPO/annotations/tcrb_domain.bed" \
  --diffs $(for c in $CMPS; do printf "%s=%s " "$c" "$OUT/pvalue/diff_${c}.stats.tsv"; done) \
  --out "$FIG/tables/per_locus_summary_tcrb_domain.tsv" >/dev/null 2>&1 \
  && cut -f1-8 "$FIG/tables/per_locus_summary_tcrb_domain.tsv" | column -t \
  || echo "per_locus_summary failed - see python hic/per_locus_summary.py --help"

echo; echo "## 2c annotation-independent block, DN-vs-DP, chr6:41,566,000-41,602,600 (Results)"
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}
  $1=="chr6" && $2>=41566000 && $3<=41602600 {
    n++; p=$h["pval"]; q=$h["pval_adj_bh"]; s=$h["bearing_score"];
    if(p<0.01){k++; if(s<0)neg++; else pos++; if(q<mq||mq=="")mq=q}
    if(q<aq||aq=="")aq=q; tot2+= ($h["kl_2"]<0?-$h["kl_2"]:$h["kl_2"]);
    for(j=1;j<=6;j++){v=$h["kl_" j]; tot+= (v<0?-v:v)} }
  END{printf "tested %d; p<0.01: %d (DN-higher %d, DP-higher %d); min q (p<0.01 bins) %s; min q (all) %s; RNA+ share of summed |kl| %.1f%%\n",
      n,k,pos,neg,mq,aq,(tot>0?100*tot2/tot:0)}' "$OUT/pvalue/diff_DN_vs_DP.stats.tsv"

echo; echo "## 2d floor-crossing share of the l1 BH-significant bins (M.3.1)"
python - "$OUT" "$SRC" <<'PY'
import csv, gzip, sys, os
out, src = sys.argv[1:3]
for cmp_ in ("DN_vs_ProB", "DN_vs_S3T3"):
    sig = set()
    with open(os.path.join(out, "pvalue", "diff_%s.stats.tsv" % cmp_)) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if str(r["significant_fdr0.05"]).strip() in ("1", "True", "true"):
                sig.add((r["chrom"], int(r["start"])))
    f = os.path.join(src, "swapnull", "%s_genome.perbin.tsv.gz" % cmp_)
    if not os.path.exists(f):
        print(cmp_, "BH bins", len(sig), "- swapnull per-bin file missing:", f); continue
    hit = cross = 0
    with gzip.open(f, "rt") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        cols = rd.fieldnames
        kc = [c for c in cols if c.lower() in ("chrom", "chr")][0]
        ks = [c for c in cols if c.lower() in ("start", "bin_start")][0]
        for r in rd:
            if (r[kc], int(float(r[ks]))) in sig:
                hit += 1
                a = [float(r["a1"]), float(r["a2"])]; b = [float(r["b1"]), float(r["b2"])]
                if (a[0] == 0 and a[1] == 0) or (b[0] == 0 and b[1] == 0):
                    cross += 1
    print("%s: BH bins %d, matched %d, floor crossings %d (%.1f%%)"
          % (cmp_, len(sig), hit, cross, 100.0 * cross / max(hit, 1)))
PY

echo; echo "## 2e Table S3 CBE counts (DN-vs-DP, DN-vs-EbKO): CBEs, with tested bins, with k>0, q<0.05"
for c in DN_vs_DP DN_vs_EbKO; do
  f="$OUT/regional/enrich_${c}_cbe.tsv"
  [ -e "$f" ] && awk -F'\t' -v c="$c" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}
    {n++; if($h["L_region_bins"]>0)t++; if($h["k"]>0)k++; if($h["q_combined"]<0.05)q++}
    END{printf "%s: %d CBEs, %d with tested bins, %d with a below-threshold bin, %d q<0.05\n",c,n,t,k,q}' "$f"
done

echo; echo "## 2f Spearman / PCA (compare/, statistic-independent; prior 1)"
for f in "$OUT"/compare/total_saliency_spearman.tsv; do [ -e "$f" ] && head -20 "$f"; done
grep -h -i "explained\|variance" "$OUT"/compare/*pca*.tsv 2>/dev/null | head
} 2>&1 | tee -a "$VAL"

###############################################################################
echo; echo "== 3. FIGURE PANELS -> $FIG/panels"
###############################################################################
P="$FIG/panels"
# --- Fig 2B (JSD heatmap; same Q as before, re-export for consistency)
grab "$OUT/compare/q_pair_jsd_heatmap.pdf" panels
# --- F3A: TAD/insulation landscape, 6 conditions (no dV1CTCF), 25 kb
grab "$OUT/hic/tcrb_tad_25000*" panels
# --- F3B: O/E grid (already correct; re-export)
grab "$OUT/paper_figures/fig3_panelB_oe_grid.png" panels
# --- F4B / S7B: compartments without dV1CTCF; 100 kb (core + dV1P), 250 kb core, + 500 kb core
grab "$OUT/paper_figures/suppS7b_compartment_100kb_v1p.compartment_landscape.*" panels
grab "$OUT/paper_figures/suppS7b_compartment_250kb_core.compartment_landscape.*" panels
python "$REPO/hic/compartment_analysis.py" \
  --pc1 DN=$COMP/merged_corrected_KR_DN_bs_500000_pca1.bw DP=$COMP/merged_corrected_KR_DP_bs_500000_pca1.bw \
        ProB=$COMP/merged_corrected_KR_ProB_bs_500000_pca1.bw EbKO=$COMP/merged_corrected_KR_EbKO_bs_500000_pca1.bw \
        S3T3=$COMP/merged_corrected_KR_S3T3_bs_500000_pca1.bw \
  --region chr6:39000000-55000000 --features-bed "$REPO/annotations/tcrb_extension_features.bed" \
  --tad-dir "$DATA/../tads" --tad-resolution 25000 --orient-with-gtf "$GTF" \
  --highlight chr6:40872100-41558435 --out-prefix "$P/suppS7b_compartment_500kb_core" \
  || echo "  S7B 500 kb failed (check that the 500 kb pca1 bigwigs exist in $COMP)"
# --- S7A: contact isolation incl. dV1P (no dV1CTCF)
grab "$OUT/hic/tcrb_isolation_OE_perm*" panels
# --- F5: Tcrb domain browser figure, DN vs DP, from the l1 stats (no FDR bins expected)
python "$REPO/bearing_hic_kl_track_plot.py" \
  --diff "$OUT/pvalue/diff_DN_vs_DP.stats.tsv" --categories "$OUT/DN_rep1_cats.json" \
  --region chr6:40400000-42400000 \
  --cool-A "$HICF/merged_corrected_KR_DN_bs_10000.cool" --cool-B "$HICF/merged_corrected_KR_DP_bs_10000.cool" \
  --insul-A "$TAD/merged_corrected_KR_DN_bs_25000_tad_score.bm" \
  --insul-B "$TAD/merged_corrected_KR_DP_bs_25000_tad_score.bm" \
  --regions-bed "$REPO/annotations/tcrb_regions_v5.bed" --gtf "$GTF" \
  --label-a DN --label-b DP --max-distance 2000000 \
  --track-pctile-low 99.0 --track-pctile-high 99.5 --threshold-scope genome \
  --out "$P/Fig5_tcrb_domain_DN_vs_DP_l1.pdf" || echo "  Fig 5 failed"
# --- F6, S3, S5A-C, S8: combined browser views (built by rule bearing_hic_combined_one)
for c in $CMPS; do
  for reg in tcrb igh rag; do grab "$OUT/figures/combined/$c/${reg}_combined.pdf" panels; done
  for f in "$P"/tcrb_combined.pdf "$P"/igh_combined.pdf "$P"/rag_combined.pdf; do
    [ -e "$f" ] && mv "$f" "${f%.pdf}_${c}.pdf"; done
done
ls "$OUT"/figures/combined/DN_vs_DP/ > "$P/combined_DN_vs_DP_listing.txt" 2>/dev/null
# --- F7B: l1 null histograms; F7E: re-render with v1.1 (labels q 0.0075 / 0.04997 / <= 5.0e-04)
grab "$OUT/paper_figures/fig7b/Fig7B_DN_vs_*" panels
python "$REPO/paper/figures/fig7e_regional_empirical.py" \
  --tables "$OUT/table_sources/regional_empirical_tcrb.tsv" "$OUT/table_sources/regional_empirical_igh.tsv" \
  --comparisons $CMPS --out-prefix "$P/Fig7E_regional_empirical_v1.1" || echo "  Fig 7E failed"
# --- S1: per-track JSD (A) + stacked bars (B)
grab "$OUT/paper_figures/suppS1_jsd_decomposition*" panels
# --- S2: genome-wide PCA + Spearman (compare/, prior 1)
grab "$OUT/compare/*pca*.pdf" panels
grab "$OUT/compare/*spearman*.pdf" panels
# --- S4: per-track input tracks over the domain (DN_rep1/2, S3T3_rep1/2)
grab "$OUT/paper_figures/suppS4_tracks/*" panels
# --- S5D / S6 / S9: per-region compare outputs (igh, new_wide)
grab "$OUT/compare/regions/metrics/igh_*.pdf" panels
grab "$OUT/compare/regions/metrics/new_wide_*.pdf" panels
grab "$OUT/compare/regions/new_wide_compare.pdf" panels
# --- S5E: cross-locus DN-vs-EbKO survey
grab "$OUT/paper_figures/qc_crosslocus_ebko_survey_rawQ*" panels
# --- S7C: within-Tcrb BES vs Hi-C (stats.json has the top-decile OR) + cross-locus
grab "$OUT/hic/bes_hic_within_DN_vs_DP*" panels
# --- S12: element classes
python "$REPO/paper/figures/make_figS12_elemclass.py" --regional-dir "$OUT/regional" \
  --out "$P/S12_elemclass_l1" | tee -a "$VAL" || echo "  S12 failed (copy make_figS12_elemclass.py first)"
# --- S13: insulation specificity (quantile) panels + table
grab "$OUT/hic/insulation_specificity/*_quantile_hl.png" panels
grab "$OUT/hic/insulation_specificity/specificity_quantile.*" panels

echo; echo "== pack"
( cd "$OUT" && tar czf figs_v25.tgz figs_v25 ) && ls -lh "$OUT/figs_v25.tgz"
echo "Upload: $OUT/figs_v25.tgz (+ this log). Then Claude assembles deck v26 and fills the markers."

###############################################################################
# 4. Supp Fig S11 (and the C6 workbook build) - ONLY after the v1p7ar jobs end
###############################################################################
# squeue -u $USER                          # must show no job using this checkout
# git am claude_bearing_c6_c7_tables_s11_2026-10-29.patch   # plot_floor_sensitivity 1.2,
#                                                          # build_tables S7/S8 fixes
# cd workflow
# snakemake --configfile config/config_results_prior1_n100_l1.generated.yaml \
#   --rerun-triggers mtime -R plot_floor_sensitivity paper_tables \
#   results_prior1_n100_l1/figures/SuppFig_S11_floor_sensitivity.png \
#   results_prior1_n100_l1/paper_tables/BEARING_tables.xlsx --cores 4
# (Panel C still uses --mechanism-comparison DN_vs_DP: decision pending.)
#
# 5. Table S2 Panel B (regional threshold sensitivity, p_thresh 0.01 / 0.10) - optional,
#    only if you keep that sentence; otherwise delete it from M.7.
# for T in 0.01 0.10; do for L in tcrb igh; do for c in DN_vs_DP DN_vs_EbKO DN_vs_ProB DN_vs_S3T3; do
#   sbatch -p defq --mem=32G --time=10:00:00 --wrap "python $REPO/dev/regional_null_calibration.py \
#     --stats $OUT/pvalue/diff_${c}.stats.tsv --regions $REPO/annotations/<tcrb_regions_v5|igh_regions_v4>.bed \
#     --locus <chr6:40790000-41690000|chr12:113000000-116100000> --chrom-sizes $REPO/workflow/resources/mm10.chrom.sizes \
#     --mode background --n-random <2000|500> --p-thresh $T --seed 42 --verify 0 --repo $REPO \
#     --out $OUT/sens/pthresh_$T/cmnull_${L}_${c}.tsv --summary-out $OUT/sens/pthresh_$T/cmnull_${L}_${c}_summary.tsv"
# done; done; done
