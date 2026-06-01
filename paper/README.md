# Reproducing the manuscript figures

Each manuscript figure maps to specific scripts in this repository. Inputs are
the processed BigWigs and Hi-C matrices (GEO/SRA, see Data Availability) plus
the annotation BEDs in `../annotations/`. Run the core pipeline
(`../workflow/Snakefile`) first to produce the qcat, differential, and p-value
files the figure scripts consume.

| Figure | Script(s) | Notes |
|--------|-----------|-------|
| Fig 1  | conceptual (schematic + equations) | graphical_abstract art |
| Fig 2  | `../compare_qcat.py` (Q-vector JSD) | per-sample background JSD heatmap |
| Fig 3  | `../hic/tad_extension_analysis.py`, `../bearing_tad_compare.py`, `../sanity_check_tad_boundaries.py` | insulation TAD landscape |
| Fig 4  | `../hic/tcrb_contact_isolation.py` (OE self-contact), `../hic/compartment_analysis.py` (PC1) | OE at 10/25 kb; compartments at 100 kb |
| Fig 5  | `../bearing_hic_combined_plot.py` (two-condition entry point), `../diffuse_architecture_top1pct.py`, `../render_regional_enrichment_heatmap.py` | DN vs DP; diffuse/focal split; regional heatmap |
| Fig 6  | `../bearing_hic_combined_plot.py` | DN vs Pro-B at Igh |
| Fig 7  | `../shift_bigwig.py`, `../generate_perm_nulls.py`, `../bearing_pvalue.py`, `../regional_enrichment.py`, `../consolidate_regional_enrichment.py` | permutation null + regional enrichment |
| S1     | `../compare_qcat.py --per-track-jsd` | per-track JSD decomposition |
| S2     | `../compare_qcat.py` (PCA + Spearman) | genome-wide reproducibility |
| S3     | `../bearing_kl_track_plot.py` | focal Tcrb zoom |
| S4     | `../bearing_hic_kl_track_plot.py` | per-track input deconstruction |
| S5     | `../bearing_hic_combined_plot.py` | Igh specificity controls |
| S6     | `../compare_qcat.py` (per-region PCA/Spearman) | domain-scale reproducibility |
| S7     | `../bes_hic_correlation.py` (within-Tcrb), `../hic/bes_hic_crosslocus.py` (pooled AR loci) | BES-Hi-C co-localization |
| S8     | `../bearing_hic_combined_plot.py`, `../compare_qcat.py` | Rag1/Rag2 dual control |
| S9     | `../bearing_hic_kl_track_plot.py` | per-replicate consistency |
| S10    | `../bearing_calibration.py` | within-condition replicate-differential FDR calibration |

| Table | Script |
|-------|--------|
| S2 (regional q-values) | `../regional_enrichment.py` + `../consolidate_regional_enrichment.py` |
| S6 (FDR calibration)   | `../bearing_calibration.py` |
| Hi-C data tables (TAD boundaries, anchor TAD, compartment switch, OE) | `../hic/` scripts |
| cross-locus null       | `../hic/bes_hic_crosslocus.py` |

## Synthetic benchmark

The synthetic-data recovery benchmark (simulated tracks with known ground
truth) lives in `../benchmark/` and is independent of the figure data above.
See `../benchmark/README.md`.
