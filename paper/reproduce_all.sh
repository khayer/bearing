#!/usr/bin/env bash
#
# reproduce_all.sh -- end-to-end reproduction of the BEARING analyses.
#
# This is a consolidated, ordered version of the real run commands. It is NOT
# meant to be run blindly with `bash reproduce_all.sh`; step through it phase by
# phase. SLURM submissions use the dbhiq,defq partitions.
#
# >>> CONFIRM THE FOUR VALUES IN THE CONFIG BLOCK BEFORE RUNNING <<<
# They determine whether the rerun reproduces the published numbers (see the
# notes next to each). ASCII-only.

set -euo pipefail

# ============================ CONFIG (CONFIRM) =============================
BEARING=~/data/tools/scripts/bigwig_to_qcat   # path to the cleaned repo
SHEET=samples_calib.tsv                        # the real sample sheet
REGIONS=regions_template.tsv                    # region set for compare/plots
GENES=mm10_genes.bed
GTF=../bearing_score/gencode.vM23.annotation_modified_overlaps_removed_sorted.gtf
CATS=$BEARING/categories/mm10_6track_panel.yaml
RESULTS=results_rerun                           # output dir for diff stats

# --- the four values to confirm ---
# (1) MIN_SIGNAL: your notes use 0.1 (and the dir is "..min_signal_test"), but
#     Methods M.4 states the default 0.01. Which produced the PUBLISHED numbers?
MIN_SIGNAL=0.1
# (2) NORMALIZE: most commands omit it (default OFF, matching Methods M.2 which
#     says normalization was NOT used); one variant used nonzero-quantile.
#     Set to "" for none, or "--normalize-method nonzero-quantile" to turn on.
NORMALIZE=""
# (3) N_PERMS: notes show 10 in some commands and 5 in others.
N_PERMS=10
# (4) MIN_SHIFT: notes show 1500000 (perm) and 1000000 (diff perm).
MIN_SHIFT=1500000
SEED=42
BLACKLIST=detected_blacklist_merged_with_encode_blacklist.bed
# ==========================================================================

# ---------------- Phase 0: blacklist (detect -> merge -> +ENCODE) ----------
sbatch -p dbhiq,defq -t 24:00:00 -c 4 --mem=64G --wrap="python3 $BEARING/bigwig_to_qcat.py \
    --sheet $SHEET --out placeholder.qcat.bgz \
    --detect-unmappable --unmappable-zero-frac 0.90 \
    --unmappable-out detected_blacklist.bed --jobs 4"
# after it finishes:
cat detected_blacklist.bed | sort -k1,1 -k2,2n | bedtools merge -d 3003 > detected_blacklist_merged.bed
cat detected_blacklist_merged.bed <(cut -f1-3 mm10-blacklist.v2.bed) | sort -k1,1 -k2,2n | bedtools merge > $BLACKLIST

# ---------------- Phase 1: observed scoring (+ write floors) ---------------
# Writes per-sample qcat AND floors.tsv (reused by the permutation nulls).
sbatch -p dbhiq,defq -t 24:00:00 -c 18 --mem=164G --wrap="python3 $BEARING/bigwig_to_qcat.py \
    --sheet $SHEET --out placeholder.qcat.bgz \
    --blacklist $BLACKLIST --min-signal $MIN_SIGNAL $NORMALIZE \
    --write-floors-tsv floors.tsv \
    --jobs 16 --signal-plots --stats median p90 p99 p80 --stats-plots \
    --bed AgRgenes_mm10_s.bed $BLACKLIST $GENES"

# sanity check (optional)
python3 $BEARING/batch_pygenometracks.py --sheet $SHEET --regions-file $REGIONS \
    --outdir pygenometracks --threads 4 --run

# ---------------- Phase 2: compare (CHANGED) + perm nulls (parallel) -------
# compare_qcat now emits PRIMARY Spearman (global-nonzero shared set) plus the
# supplementary spearman_all_aligned_bins.* -- diff the TSVs against the paper.
sbatch -p dbhiq,defq -t 24:00:00 -c 12 --mem=200G --wrap="python3 $BEARING/compare_qcat.py \
    --sheet $SHEET --out $RESULTS \
    --regions-file $REGIONS --genes $GENES --no-clip --workers 10"

sbatch -p dbhiq,defq -t 24:00:00 -c 12 --mem=200G --wrap="python3 $BEARING/generate_perm_nulls.py \
    --sheet $SHEET --out-dir perm --n-perms $N_PERMS \
    --blacklist $BLACKLIST --floors-tsv floors.tsv --min-signal $MIN_SIGNAL \
    --min-shift $MIN_SHIFT --seed $SEED --shift-workers 10 --jobs 10"

# DATA QUIRK: DP rep2 reuses rep1's CTCF perm track. Replicate the manual copy
# for each perm round (S25 -> S29) before p-values, or fix the sheet so DP_rep2
# CTCF points to the intended file:
# for k in $(seq 1 $N_PERMS); do cp perm/perm$k/..._S25_perm$k.bw perm/perm$k/..._S29_perm$k.bw; done

# ---------------- Phase 3: per-bin p-values --------------------------------
COMPARISON_DIR_OVERRIDE=. $BEARING/submit_perm_qcat_pvalue_slurm.sh \
    $SHEET $N_PERMS perm $RESULTS \
    4 16G 04:00:00  8 24G 04:00:00  4 24G 02:00:00  2 24G 02:00:00 \
    $MIN_SHIFT $BLACKLIST  floors.tsv $MIN_SIGNAL  $SEED 1000000 $SHEET

# ---------------- Phase 4: differential nulls + diff p-values --------------
sbatch -p dbhiq,defq -t 24:00:00 -c 12 --mem=200G --wrap="python $BEARING/generate_perm_nulls.py \
    --sheet $SHEET --diff-sheet $SHEET --n-perms $N_PERMS \
    --out-dir perm --blacklist $BLACKLIST --floors-tsv floors.tsv --jobs 10"
sbatch -p dbhiq,defq --cpus-per-task=4 --mem=256G --time=24:00:00 \
    --wrap 'export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1; \
    bash '"$BEARING"'/run_perm_diff_pvalue.sh . '"$N_PERMS"' perm/perm '"$RESULTS"

# ---------------- Phase 5: regional enrichment (all comparisons) -----------
for cmp in DN_vs_DP DN_vs_EbKO DN_vs_ProB DN_vs_S3T3; do
  sbatch -p dbhiq,defq -t 24:00:00 -c 2 --mem=100G --wrap="python3 $BEARING/regional_enrichment.py batch \
      --diff-table $RESULTS/diff_${cmp}.stats.tsv \
      --regions tcrb_regions_v5.bed --locus chr6:40790000-41690000 \
      --out enrich_${cmp}_tcrb_regions_v5.tsv"
done
python3 $BEARING/consolidate_regional_enrichment.py \
    --tsvs "enrich_DN_vs_*_tcrb_regions_v5.tsv" "enrich_DN_vs_*_igh_*.tsv" \
    --out consolidated_enrichment.tsv
python3 $BEARING/render_regional_enrichment_heatmap.py \
    --tsv consolidated_enrichment.tsv --out panel_E_heatmap.pdf \
    --comparisons DN_vs_DP DN_vs_EbKO DN_vs_Pro-B DN_vs_3T3

# ---------------- Phase 6: BEARING + Hi-C figures (batch) -------------------
sbatch -p dbhiq,defq -t 24:00:00 -c 2 --mem=100G --wrap="python3 $BEARING/batch_bearing_hic_plots.py \
    --sheet $SHEET --regions-file $REGIONS --reference-condition DN \
    --contact DN=../hic_files/merged_corrected_KR_DN_bs_10000.cool \
    --contact DP=../hic_files/merged_corrected_KR_DP_bs_10000.cool \
    --contact EbKO=../hic_files/merged_corrected_KR_EBKO_bs_10000.cool \
    --contact ProB=../hic_files/merged_corrected_KR_ProB_bs_10000.cool \
    --contact S3T3=../hic_files/merged_corrected_KR_s3T3_bs_10000.cool \
    --comparison-dir . --gtf $GTF --results-dir $RESULTS/ \
    --outdir hic_batch --triangle --run"

# ---------------- Phase 7: Hi-C analyses (Figs 3, 4, S7) -------------------
# 7-condition set incl. the two V1-perturbation Hi-C backgrounds.
CONDS="DN DP ProB EbKO S3T3 dV1P dV1CTCF"
for res in 100000 25000 250000 500000; do
  python3 $BEARING/hic/tad_extension_analysis.py --tad-dir ../tads/ --resolution $res \
      --conditions $CONDS --anchor chr6:41550000 \
      --features-bed tcrb_extension_features.bed --out-prefix tcrb_tad_${res} --plot
done
for res in 100000 250000 500000; do
  python3 $BEARING/hic/compartment_analysis.py \
      --pc1 $(for c in $CONDS; do echo $c=../comp/merged_corrected_KR_${c}_bs_${res}_pca1.bw; done | tr '\n' ' ') \
      --region chr6:39000000-55000000 --features-bed tcrb_extension_features.bed \
      --tad-dir ../tads/ --tad-resolution $res --orient-with-gtf $GTF \
      --out-prefix tcrb_compartment_${res}
done
python3 $BEARING/hic/tcrb_contact_isolation.py \
    --cool $(for c in $CONDS; do echo $c=../hic_files/merged_corrected_KR_${c}_bs_10000.cool; done | tr '\n' ' ') \
    --target chr6:40850000-41600000 --window-size 750000 --max-distance 6000000 \
    --perm-n 200 --out-prefix tcrb_isolation_OE_perm_10K
# cross-locus BES-Hi-C co-localization (pooled AR loci)
python3 $BEARING/hic/bes_hic_crosslocus.py \
    --diffs DN_vs_DP=$RESULTS/diff_DN_vs_DP.stats.tsv DN_vs_ProB=$RESULTS/diff_DN_vs_ProB.stats.tsv \
    --cool-a ../hic_files/merged_corrected_KR_DN_bs_10000.cool \
    --cool-b DN_vs_DP=../hic_files/merged_corrected_KR_DP_bs_10000.cool \
             DN_vs_ProB=../hic_files/merged_corrected_KR_ProB_bs_10000.cool \
    --targets targets.tsv --hic-bin 10000 --aggregation p95 --contact-metric delta_contact \
    --n-controls 200 --n-panels 1000 --blacklist mm10-blacklist.v2.bed \
    --gtf $GTF --match-gene-density --out-prefix crosslocus_p95_dcontact

# ---------------- Phase 8: within-condition FDR calibration (Table S6) -----
bash $BEARING/run_replicate_calibration.sh   # uses the same null machinery

# ---------------- Phase 9: cross-comparison decomposition report -----------
python3 $BEARING/bearing_decomposition_report.py \
    --diff-dir $RESULTS --categories DN_rep1_cats.json --out cross_comparison --top-n 500

# ---------------- Phase 10: synthetic benchmark (simulated data) -----------
python3 $BEARING/benchmark/simulate_bearing_tracks.py --prefix sim --seed 1
python3 $BEARING/benchmark/evaluate_bearing_recovery.py --prefix sim --out-prefix recovery
bash $BEARING/benchmark/run_bearing_recovery_sweep.sh
python3 $BEARING/benchmark/plot_bearing_sweep.py --master sweep_master.tsv --out-prefix sweep
