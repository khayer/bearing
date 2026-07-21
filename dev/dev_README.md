# dev/ -- analyses that are not workflow rules

This directory holds two kinds of script, and they must not be confused:

  A. LOAD-BEARING analyses. They produce numbers, tables, or figures that appear
     in the manuscript, but they are run BY HAND rather than as Snakemake rules.
     Each is mapped below to exactly what it produces. These are part of the
     paper's provenance and must remain runnable.

  B. SCRATCH and superseded scripts. One-off diagnostics, and earlier versions
     of the load-bearing scripts that have been replaced. Retained for history
     only. Nothing in the paper depends on them. They are listed at the bottom so
     no one runs a superseded version by mistake.

Because the load-bearing scripts are not in the Snakemake DAG, nothing detects
when their outputs are stale relative to workflow/results/pvalue.done. Re-run
them after any change to the p-value layer and check with staleness_audit.sh.
Their outputs belong in paper/table_sources/ (tracked in git, small), which is
where paper/build_tables.py reads them from -- NOT in workflow/results/ (which is
gitignored, so a clone could not rebuild the tables from it).

===========================================================================
A. LOAD-BEARING -- maps to a manuscript table, figure, or quoted number
===========================================================================

replicate_stability.py
    Produces: Table S13 (leave-one-replicate-out stability).
    Run once per non-DN condition:
        --cond-a "DN:DN_rep1,DN_rep2" --cond-b "DP:DP_rep1,DP_rep2"   (and ProB, S3T3)
    Output: paper/table_sources/replicate_stability_DN_vs_{DP,ProB,S3T3}.tsv
    NOTE: EbKO is NOT a contrast here. The "Ebeta/Trbv31/CBE3 (EbKO deletion)"
    rows are a REGION analysed within the DN-vs-DP contrast; there is correctly
    no replicate_stability_DN_vs_EbKO.tsv.

track_ablation.py
    Produces: Table S12 (per-track leave-one-out attribution).
    Output: paper/table_sources/track_ablation_*.tsv

baseline_comparison.py
    Produces: Table S10 (baseline method comparison).
    Output: paper/table_sources/baseline_comparison*.tsv

regional_null_calibration.py            (was regional_null_calibration_v2.py)
    Produces: Table S11 (regional enrichment vs autocorrelation-corrected null).
    Run with --mode background --n-random 500 --seed 42, per locus:
    Output: paper/table_sources/regional_null_calibration_{tcrb,igh}_DNrep_bg.tsv
    OPEN ISSUE (see paper/reproduce_all.sh): background loci are matched on
    scored-bin count to the NULL contrast, not the production contrast. The
    correct matching must be settled before the empirical p-values are final.

pvminsig_sweep.py                       (repo root, not dev/)
    Produces: Table S14 and Supplementary Figure S11 (differential-floor
    sensitivity). Output: paper/table_sources/sens/floor_sweep.tsv

region_match_ebko.py            (canonical; supersedes region_match_ebko2/3.py)
    Produces: the EbKO recombination-centre characterisation quoted in Results
    (per-track share of the joint differential, 36.3% RNA+; 100% of the 55
    edgeR-significant plus-strand RC bins DN-directed; edgeR-logFC vs BEARING
    RNA+ dKL Spearman rho=0.43, p=7e-4). Uses a 200 bp grid join that absorbs
    the edgeR 1-based / BEARING 0-based offset (the earlier ebko2 raw-start join
    matched zero bins; ebko3 fixed the join but dropped the composition block --
    both retired).
    Reads DN-vs-EbKO diff plus rna_concordance_DN_vs_EbKO_{pos,neg}_edgeR_allbins.csv
    (produced by rna_concordance_stranded.R).
    Output: paper/table_sources/region_match_ebko_rc.tsv (via --out).

===========================================================================
B. SCRATCH and SUPERSEDED -- not used by the paper; do not run for results
===========================================================================

Superseded (an A-script replaced these; kept only for history):
    (region_match_ebko2.py and region_match_ebko3.py were retired once
     region_match_ebko.py merged the correct join with the composition block.)

Diagnostic that INFORMED a methodological choice (its number is not quoted):
    score_autocorrelation.py
        Measures the spatial decorrelation length of the per-bin BEARING score
        (L(1/e) ~ 395 bp, L(0.5) ~ 237 bp on DN-vs-DP; ~2x the 200 bp grid).
        This is the quantitative basis for the Methods statement that within-
        track autocorrelation makes the binomial spatial null unreliable, so
        the regional test is calibrated empirically instead. The manuscript
        makes that argument QUALITATIVELY ("autocorrelation at the kb scale")
        and does not cite the number, so this is a diagnostic, not the source
        of a quoted value. --out / --summary-out persist it to
        paper/table_sources/ for anyone who wants to check the claim.

One-off diagnostics (never produced a manuscript number):
    tail_diag.py, diagnose_flat_pvalue.py, probe_capture_depth.py,
    feature_changes_across_runs.py, cbe_across_runs.py, compare_regimes.py,
    robustness_summary.py, hic_probe.py, bearing_specificity_profile.py,
    prototype_adaptive_binning.py, tail_diag.py, plot_sig_density.py,
    plot_bearing_stats_summary.py, bin_hic_superbins.py, plot_hic_superbins.py,
    run_bes_hic_matrix.py, bearing_tad_decomposition.py, make_bearing_bins_bed.py,
    extract_track_diff_bed.py, bearing_diff_concordance.py, regen_fig7_panelB.py

If you promote any B-script to produce a paper number, MOVE it up to section A
with its output path, or the mapping stops being trustworthy.

===========================================================================
The goal
===========================================================================
Every number in the manuscript traces to a script. For the tables that is
enforced by paper/build_tables.py (each sheet names its source file, and the
Provenance sheet records the SHA256). For the numbers quoted in the TEXT, this
README is the map. region_match_ebko.py (the canonical EbKO RC script) now
writes its numbers to a file, and score_autocorrelation.py -- a diagnostic
whose number is not quoted -- also persists its output. Every value quoted in
the text now traces to a file: "every number has a file" holds.
