#!/usr/bin/env bash
# Pin the beat-2 boundary two ways:
#  (1) COARSER PASS  -- rebin contacts to 2kb and 5kb super-bins; does V1P stripe
#      loss rise above the V1TxS noise floor when per-bin SNR improves?
#  (2) REPLICATE SPLIT -- DN rep1 vs DN rep2, a pure within-condition null; does a
#      known-null contrast pile into the RC just like V1TxS and the mutants?
#
# Run from REPO ROOT with the bearing env active. V1TxS must be COMMENTED in
# comparisons_v1_4track.tsv (the 4-track run is V1P/RCTKO-only). ASCII only.
set -uo pipefail

BAMDIR=/mnt/isilon/bassing_lab/projects/capture_hic/Brittney_V1/03aln
COMMON="--bam-dir $BAMDIR --annot anchors_v1.bed --std-binsize 3000 \
  --min-distance 10000 --max-distance 1000000 --min-shift 50000"
CORE="chr6:40800000-41650000"
OUT=mcc_boundary
mkdir -p $OUT

run_bench () {  # grid counts comparisons prefix outtsv
    python mcc/scripts/benchmark_grid_vs_standard.py --superbins "$1" \
        --counts "$2" --comparisons "$3" $COMMON --out "$5" --dump-bins "$4" \
        2> "$5.log" || echo "  [warn] benchmark exited nonzero (see $5.log) -- dumps before the failing comparison are still written"
}

core_partials () {  # prefix
    for pair in DN_vs_V1P DN_vs_RCTKO DN_vs_V1TxS; do
        f="$1_${pair}.tsv"
        if [ -f "$f" ]; then
            echo "  -- $pair --"
            python mcc/scripts/partial_corr.py "$f" --region $CORE \
                | grep -E "rho_raw|rho \| depth |well-covered|premise" | sed 's/^/    /'
        fi
    done
}

echo "############################################################"
echo "# (1) COARSER PASS"
echo "############################################################"
for W in 2000 5000; do
    echo; echo "======== ${W}bp super-bin grid ========"
    GRID=$OUT/sb${W}.bed
    python mcc/scripts/build_superbin_grid.py --seg seg_v1_chr6_real.bed \
        --min-width $W --out-bed $GRID
    run_bench $GRID capture_cools/pair_counts_v1.tsv \
        mcc/config/comparisons_v1_4track.tsv $OUT/bw${W} $OUT/gb_${W}.tsv
    run_bench $GRID capture_cools/pair_counts_v1txs.tsv \
        mcc/config/comparisons_v1txs.tsv $OUT/bw${W} $OUT/gb_${W}_txs.tsv
    echo "[core partials @ ${W}bp]"; core_partials $OUT/bw${W}
    echo "[tail vs V1TxS null @ ${W}bp]"
    python3 mcc/scripts/tail_diag_generic.py --null $OUT/bw${W}_DN_vs_V1TxS.tsv \
        --test $OUT/bw${W}_DN_vs_V1P.tsv $OUT/bw${W}_DN_vs_RCTKO.tsv --region $CORE
done

echo; echo "############################################################"
echo "# (2) REPLICATE-SPLIT NULL  (DN rep1 vs DN rep2, native grid)"
echo "#     NOTE: floor = rep1 depth (2.6M) < mutant floor (8.3M),"
echo "#     so this null is NOISIER -- read its SHAPE (RC pile-up), not magnitude."
echo "############################################################"
printf "DN1\traw_ArimaHTS_S065_R1KO_rep1_bs_250_master_valid_pairs.bam\t2599078\n"  > $OUT/pc_rep.tsv
printf "DN2\traw_ArimaHTS_S065_R1KO_rep2_bs_250_master_valid_pairs.bam\t5746499\n" >> $OUT/pc_rep.tsv
printf "locus\tregion\tcondA\tcondB\tbes\n" > $OUT/comp_rep.tsv
printf "Tcrb\t%s\tDN1\tDN2\tworkflow/results_v1_4track/pvalue/diff_DN_vs_V1P.stats.tsv\n" \
    "chr6:39771104-42660754" >> $OUT/comp_rep.tsv
run_bench seg_v1_chr6_real.bed $OUT/pc_rep.tsv $OUT/comp_rep.tsv $OUT/bw_rep $OUT/gb_rep.tsv
echo "[tail: replicate-split null + V1P vs the V1TxS reference null, native grid]"
python3 mcc/scripts/tail_diag_generic.py --null mcc_results_v1txs/binwise_DN_vs_V1TxS.tsv \
    --test $OUT/bw_rep_DN1_vs_DN2.tsv mcc_results_v1_4track/binwise_DN_vs_V1P.tsv --region $CORE

echo; echo "DONE. Key reads:"
echo " * COARSER PASS: if V1P %>p95 climbs well above 5% as the grid coarsens"
echo "   (and its core rho_raw grows more negative), real stripe loss was being"
echo "   buried by 200bp noise -> a weak beat 2 survives at honest resolution."
echo "   If V1P stays ~5% at every grid, the limit is fundamental."
echo " * REPLICATE SPLIT: if DN1-vs-DN2 also concentrates in the RC like V1TxS,"
echo "   the RC pile-up is confirmed as where MCC noise lives, not biology."
