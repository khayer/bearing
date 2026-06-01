# Paper figure grid wrappers

Composite-panel scripts for layouts that no single analysis script emits as one
call. The analysis (O/E values, compartment PC1, p-values) is done upstream;
these only arrange results into the published grids.

| Script | Produces | Auto-run by workflow? |
|--------|----------|-----------------------|
| `oe_contact_grid.py` | Fig 3B, Fig 4A (O/E heatmap grid, shared scale + colorbar) | Yes (rules `fig3b_oe_grid`, `fig4a_self_contact_grid`) |
| `compartment_resolution_grid.py` | Supp S7B (resolution x condition compartment grid) | Yes (rule `suppS7b_compartment_grid`) |
| `montage_panels.py` | Supp S9 / any composite (stitch rendered PNGs) | No -- manual; see note |

Outputs land in `results/paper_figures/` and are part of the default `all`
target whenever the sample sheet defines Hi-C conditions.

`montage_panels.py` is a manual fallback. For Supp S9 (10-replicate stack) the
recommended source is the `compare_qcat` region plot at the `new_wide` region,
which already renders all samples; use the montage only if you batch per-replicate
panels separately.

Dependencies: cooler (oe grid), pyBigWig (compartment grid), Pillow (montage).
