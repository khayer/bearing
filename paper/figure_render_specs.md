# BEARING figure render specifications -- v3, reconciled + gaps closed
# Generated 2026-05-31. Reconciled from canonical scripts in:
#   <bearing_repo>/hic/             (Hi-C analyses)
#   <bearing_repo>/                 (BEARING / qcat scripts at root)
#   <bearing_repo>/paper/figures/   (composite grid wrappers, NEW in v3)
#
# Track palette: USE categories/mm10_6track_panel.yaml. Cohesin = #8b0000.
# Ignore prior inline --colors specs.
#
# v3 CHANGE: the three layout/stitching gaps from v2 are now CLOSED by wrapper
# scripts under paper/figures/. Two are auto-run by the Snakemake workflow (they
# are part of the default `all` target when the sheet defines Hi-C conditions);
# the third (S9) is handled by the compare_qcat region plot. No manual stitching
# required for S3B/S4A/S7B.

================================================================================
SPEC 1. Fig 3 -- Tcrb recombination-domain TAD boundary
================================================================================

Panel A -- Insulation-score landscape across 7 conditions
  SCRIPT:  hic/tad_extension_analysis.py
  COMMAND:
    python hic/tad_extension_analysis.py \
      --tad-dir ../tads \
      --resolution 25000 \
      --conditions DN,DP,EbKO,ProB,S3T3,dV1P,dV1CTCF \
      --region chr6:40400000-42400000 \
      --features-bed annotations/tcrb_extension_features.bed \
      --anchor chr6:41550000 \
      --out-prefix fig3_panelA_insulation \
      --plot
  NOTES:   Multi-condition insulation + boundary landscape; shared y is the
           default under --plot. Anchor chr6:41,550,000 = recombination-center
           boundary.
  OUTPUT:  fig3_panelA_insulation.png

Panel B -- O/E contact heatmap grid (5-condition tidy matrix)   [GAP CLOSED]
  SCRIPT:  paper/figures/oe_contact_grid.py
  AUTO:    yes -- Snakemake rule `fig3b_oe_grid`
  COMMAND:
    python paper/figures/oe_contact_grid.py \
      --cool DN=../hic_files/...DN..._25000.cool \
             DP=...DP... EbKO=...EbKO... ProB=...ProB... S3T3=...S3T3... \
      --region chr6:40400000-42400000 \
      --resolution 25000 --vmin -2 --vmax 2 --ncols 5 \
      --annotate-oe-region chr6:41000000-41600000 \
      --out results/paper_figures/fig3_panelB_oe_grid.png
  NOTES:   Reads .cool via cooler; computes log2 O/E by per-diagonal expected
           (same method as tcrb_contact_isolation.py); one shared colorbar; mean
           linear O/E annotated per panel over --annotate-oe-region.
  OUTPUT:  results/paper_figures/fig3_panelB_oe_grid.png  (1800 x ~600 px)

================================================================================
SPEC 2. Fig 4 -- Self-contact unit + A/B compartments
================================================================================

Panel A -- Per-condition O/E heatmaps + self-contact O/E annotation [GAP CLOSED]
  SCRIPT:  paper/figures/oe_contact_grid.py
  AUTO:    yes -- Snakemake rule `fig4a_self_contact_grid`
  COMMAND:
    python paper/figures/oe_contact_grid.py \
      --cool DN=... DP=... EbKO=... ProB=... S3T3=... \
      --region chr6:39900000-42900000 \
      --resolution 25000 --vmin -2 --vmax 2 --ncols 5 \
      --annotate-oe-region chr6:40400000-42400000 \
      --out results/paper_figures/fig4_panelA_self_contact_grid.png
  NOTES:   Caption O/E values (DN 1.71, DP 2.03, EbKO 1.61, Pro-B 1.99,
           S3T3 0.93) appear as the per-panel annotation (mean linear O/E over
           the Tcrb self-contact block). For the permutation p-values use
           hic/tcrb_contact_isolation.py --perm-n (same as Supp S7A).
  OUTPUT:  results/paper_figures/fig4_panelA_self_contact_grid.png  (1500 x 750 px)

Panel B -- A/B compartment eigenvector at 100 kb
  SCRIPT:  hic/compartment_analysis.py
  COMMAND:
    python hic/compartment_analysis.py \
      --pc1 DN=../comp/...DN..._100000_pca1.bw DP=... EbKO=... ProB=... S3T3=... \
      --region chr6:39900000-42900000 \
      --features-bed annotations/tcrb_extension_features.bed \
      --orient-with-gtf resources/gencode.vM23...gtf \
      --out-prefix fig4_panelB_compartments
  NOTES:   --highlight via --features-bed; PC1 sign via --orient-with-gtf
           (gene-density), matching paper convention.
  OUTPUT:  fig4_panelB_compartments.png

================================================================================
SPEC 3. Fig 6 -- Cross-locus extension to Igh, DN vs Pro-B
================================================================================
(Same panels appear in journal slide9 and portrait slide18.)

RECOMMENDED COMBINED RENDER (Panels A + C in one call):
  SCRIPT:  bearing_hic_combined_plot.py
  COMMAND:
    python bearing_hic_combined_plot.py \
      --contact-a ../hic_files/...DN..._25000.cool \
      --contact-b ../hic_files/...ProB..._25000.cool \
      --label-a DN --label-b ProB \
      --region chr12:113000000-116100000 --resolution 25000 \
      --diff-qcat diff_DN_vs_ProB.qcat.bgz --rgb-hic \
      --gtf resources/gencode.vM23...gtf \
      --out fig6_combined_AC.png
  NOTES:   Overlays differential bristle tracks AND Hi-C contact difference in
           one figure. Use for combined A+C; drop the standalone Panel A call.
  OUTPUT:  fig6_combined_AC.png

ALTERNATIVE (Panel A standalone, bristle only):
  SCRIPT:  bearing_kl_track_plot.py  (canonical bristle script)
  COMMAND:
    python bearing_kl_track_plot.py \
      --diff diff_DN_vs_ProB.qcat.bgz \
      --categories categories/mm10_6track_panel.yaml \
      --region chr12:113000000-116100000 \
      --regions-bed annotations/igh_regions_v4.bed \
      --out fig6_panelA_igh_bristle.png
  NOTES:   Tracks + colors from the YAML; do NOT pass inline --colors. Per-region
           q-values (VH_proximal q=2.2e-13 etc.) come from regional_enrichment.py
           (Phase 5); annotate via --regions-bed name column or deck-assembly.
  OUTPUT:  fig6_panelA_igh_bristle.png

Panel B -- Per-track decomposition heatmap (6 tracks x 5 regions)
  SCRIPT:  bearing_region_decomposition.py
  COMMAND:
    python bearing_region_decomposition.py \
      --diff-dir results/compare \
      --regions annotations/igh_regions_v4.bed \
      --categories DN_rep1_cats.json \
      --filename-pattern "diff_{comparison}.stats.tsv" \
      --out fig6_panelB_igh_decomp
  NOTES:   No --locus flag; region set defines the locus. Comparison selected
           via --filename-pattern + diff dir contents.
  OUTPUT:  fig6_panelB_igh_decomp.png

================================================================================
SPEC 4. Supp S6 -- Domain-scale Tcrb reproducibility
================================================================================

Panel A -- Per-region PCA at domain scale (chr6:40,400,000-42,400,000)
  SCRIPTS: rebin_qcat.py (only if changing 200 bp resolution; usually not needed)
           + compare_qcat.py (PCA global; per-region output is the panel)
  COMMAND:
    python compare_qcat.py \
      --sheet results/samples.bearing.tsv \
      --out results/compare_domain \
      --categories categories/mm10_6track_panel.yaml \
      --regions-file <(grep -E "^(name|new_wide)" workflow/config/regions.tsv) \
      --no-clip --workers 10
  NOTES:   Domain-scoped PCA = per-region output for the 'new_wide' row
           (chr6:40,400,000-42,400,000). Writes *_pca.pdf per per-region output.
  OUTPUT:  suppS6_panelA_domain_pca.png  (new_wide pca.pdf)

Panel B -- Per-track Spearman matrix at domain scale
  SCRIPT:  compare_qcat.py
  COMMAND: Same as Panel A (one invocation produces PCA and Spearman).
  NOTES:   PRIMARY Spearman uses the global-nonzero shared bin set (Methods M.9).
           USE THIS, not spearman_all_aligned_bins.*. Files: spearman.pdf +
           total_saliency_spearman.tsv.
  OUTPUT:  suppS6_panelB_domain_spearman.png  (new_wide spearman.pdf)

================================================================================
SPEC 5. Supp S7 -- Hi-C robustness controls
================================================================================

Panel A -- V1-perturbation contact robustness (dV1P, dV1CTCF)
  SCRIPT:  hic/tcrb_contact_isolation.py
  AUTO:    yes -- Snakemake rule `hic_contact_isolation` (covers the perm p-values)
  COMMAND:
    python hic/tcrb_contact_isolation.py \
      --cool dV1P=...dV1P..._25000.cool dV1CTCF=...dV1CTCF..._25000.cool \
      --target chr6:40400000-42400000 --window-size 750000 \
      --max-distance 6000000 --perm-n 200 --out-prefix suppS7_panelA_v1pert
  NOTES:   Caption values (dV1P O/E 1.93 p=0.01; dV1CTCF 1.68 p=0.02) come from
           --perm-n output (TSV).
  OUTPUT:  suppS7_panelA_v1pert.png

Panel B -- Compartment-resolution robustness (25/50/100 kb x 5 cond)  [GAP CLOSED]
  SCRIPT:  paper/figures/compartment_resolution_grid.py
  AUTO:    yes -- Snakemake rule `suppS7b_compartment_grid`
  COMMAND:
    python paper/figures/compartment_resolution_grid.py \
      --pc1 DN=../comp/...DN..._{res}_pca1.bw DP=... EbKO=... ProB=... S3T3=... \
      --resolutions 25000 50000 100000 \
      --region chr6:39000000-43500000 \
      --landmarks Cntnap2:chr6:46400000:B Gpnmb:chr6:48900000:A Tcrb:chr6:41000000:A \
      --out results/paper_figures/suppS7_panelB_compartment_res.png
  NOTES:   {res} in each --pc1 path is substituted per row. Rows = resolutions,
           cols = conditions; A=red / B=blue fill; landmarks dotted. Single
           composite figure (no stitching).
  OUTPUT:  results/paper_figures/suppS7_panelB_compartment_res.png  (1200 x 800 px)

================================================================================
SPEC 6. Supp S9 -- Per-replicate consistency at recombination domain
================================================================================

  RECOMMENDED:  compare_qcat.py region plot (renders all 10 replicates for a region)
  COMMAND:
    python compare_qcat.py \
      --sheet results/samples.bearing.tsv \
      --out results/compare_domain \
      --categories categories/mm10_6track_panel.yaml \
      --regions-file <(grep -E "^(name|new_wide)" workflow/config/regions.tsv) \
      --no-clip --workers 10
  NOTES:   The multi-sample row stack IS the compare_qcat region plot at the
           new_wide row -- same invocation as Supp S6. No separate render needed.
  FALLBACK:  paper/figures/montage_panels.py -- stitch per-replicate PNGs if you
             batch bearing_kl_track_plot.py once per replicate instead. Manual.
  OUTPUT:  suppS9_perreplicate_domain.png  (new_wide region plot)

================================================================================
SPEC 7. Supp S10 -- Within-condition replicate-differential FDR calibration
================================================================================

  SCRIPT:  bearing_calibration.py (exact match)
  AUTO:    yes (per-condition) -- Snakemake `calibration` rule wraps
           run_replicate_calibration.sh. 5-panel composite arranged in
           deck-assembly OR via paper/figures/montage_panels.py.
  COMMAND (per condition; loop the 5):
    python bearing_calibration.py \
      --diff-qcat results/compare/diff_{cond}_rep1_vs_rep2.qcat.bgz \
      --null-qcat results/perm/perm*/diff_{cond}_rep1_vs_rep2.null.qcat.bgz \
      --bearing-pvalue bearing_pvalue.py \
      --min-signal 0.1 --fdr 0.05 --label {cond} \
      --out-prefix suppS10_fdr_{cond}
  NOTES:   Behind Supp Table S6 (lambda DN 0.37, EbKO 0.41, DP 0.41, Pro-B 0.45,
           S3T3 0.35; BH q<0.05 bins 384-1759). Emits QQ/histogram + lambda +
           sig-bin counts per condition.
  OUTPUT:  suppS10_fdr_calibration.png (composite via montage_panels.py)

================================================================================
SPEC 8. Journal Fig 6 == Spec 3.  Same files; no separate render.
================================================================================

================================================================================
SUMMARY: invented name -> real script
================================================================================
plot_insulation_landscape.py        -> hic/tad_extension_analysis.py
plot_hic_oe_grid.py                 -> paper/figures/oe_contact_grid.py
plot_self_contact_panels.py         -> hic/tcrb_contact_isolation.py (+ oe_contact_grid.py for the grid)
plot_compartment_tracks.py          -> hic/compartment_analysis.py
plot_compartment_resolution_grid.py -> paper/figures/compartment_resolution_grid.py
plot_bristle.py                     -> bearing_kl_track_plot.py
plot_per_track_decomp.py            -> bearing_region_decomposition.py
plot_hic_diff_overlay.py            -> bearing_hic_combined_plot.py
rebin_and_pca.py                    -> rebin_qcat.py + compare_qcat.py
spearman_per_track.py               -> compare_qcat.py
plot_pvalue_calibration.py          -> bearing_calibration.py

================================================================================
GAP STATUS (was three open in v2; all addressed in v3)
================================================================================
1. Fig 3B / Fig 4A O/E matrix grid   -> CLOSED: paper/figures/oe_contact_grid.py
   (auto: rules fig3b_oe_grid, fig4a_self_contact_grid)
2. Supp S7B compartment grid          -> CLOSED: paper/figures/compartment_resolution_grid.py
   (auto: rule suppS7b_compartment_grid)
3. Supp S9 replicate stack            -> CLOSED: use compare_qcat region plot
   (paper/figures/montage_panels.py is the manual fallback)

Remaining non-stitching item (not a gap): Fig 6A per-region q-value annotation
comes from regional_enrichment.py and is overlaid at deck-assembly time.
