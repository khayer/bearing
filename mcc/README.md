# MCC adaptive-vs-uniform binning (V1-region Tcrb)

Applies BEARING to the Allyn et al. 2025 V1-region Tcrb mutants on capture-Hi-C
(MCC), to answer one question: **does feature-anchored (adaptive) binning recover
the 1D-3D relationship better than uniform binning?** It does, but only under
specific conditions, and the sign of the relationship depends on the window --
both documented below.

This module is self-contained. The 1D BEARING differential (`diff_*.neglog10p.bw`)
must already exist (run `../workflow/Snakefile` with `config_v1` first); everything
3D is built here, on depth-matched data.

---

## 1. The cool files

`capture_cools/cool/merged_capture_{DN,V1P,RCTKO}_bs_{res}.cool`

Per-condition contact matrices built from the HiC-Pro valid-pairs BAMs
(`build_capture_cools.sh`): `samtools view | validpairs_to_pairs.py | sort | bgzip`
then `cooler cload pairs` per replicate and `cooler merge` across replicates.

**These are DEPTH-MATCHED.** Raw pooled depth differs ~2x across conditions
(DN ~8.3M, V1P ~12.5M, RCTKO ~15.7M pairs). `--downsample` subsamples each
condition's pairs to the shallowest (DN) before binning:

| label | pooled pairs | keep fraction |
|-------|--------------|---------------|
| DN    | 8,345,577    | 1.000         |
| V1P   | 12,517,693   | 0.667         |
| RCTKO | 15,671,462   | 0.533         |

Without this, a side-by-side DN-vs-mutant matrix is confounded: the deeper
condition looks denser everywhere from sequencing, not biology. The benchmark
applies the *same* keep fractions internally, so the statistics and the figures
see the same depth-matched data. **Always build with `--downsample` for any
cross-condition comparison.**

---

## 2. Binning schemes

| scheme   | grid                                   | mean width | use |
|----------|----------------------------------------|------------|-----|
| adaptive | `seg_v1_chr6_real.bed` (native)        | ~496 bp    | feature-anchored, the method under test |
| uniform  | mean-width tiling (matched resolution) | ~495 bp    | resolution-matched control |
| standard | fixed bins                             | 3000 bp    | conventional Hi-C resolution |
| floored  | `superbins_v1r_1000.bed`               | ~976 bp    | adaptive floored to contact-pixel scale |

**Adaptive is genuinely feature-anchored.** Super-bins overlapping CTCF/CBE
anchors are finer than off-anchor bins (median 50 vs 350 bp; width ratio 0.143,
feature-shift p=0.009; `superbin_feature_overlap.py`). This is real punctate AgR
structure, not a budget artifact -- but it only appears with a realistic
segmentation budget. **Build chr6 with `--chroms chr6 --target-bins 110000`**
(its ~5.5% genome share); running a genome-wide budget on one chromosome pins
every bin to the 50 bp floor and erases the anchoring.

**Flooring inverts it.** Floor the grid to a contact-resolvable scale (~1 kb) and
on-anchor bins become *wider* than off-anchor (ratio 1.75). The fine structure
lives below ~1 kb, so the adaptive advantage is recoverable only through the
**1D projection** of the contact change (the per-bin marginal), which aggregates
contacts across the bin and tolerates bins far finer than the 2D pixel floor.

---

## 3. The benchmark

`benchmark_grid_vs_standard.py` -- per comparison, on each grid:

- bins the (downsampled) pairs, computes the band-limited (10 kb-1 Mb) projected
  contact change per bin (`delta_contact`),
- correlates it (Spearman `rho_contact`) with the per-bin BEARING differential
  (`bes`), against a circular-shift null (`p_emp`),
- reports `feat_concentration` (how concentrated the contact change is on the
  finest bins) and `zero_marginal_bins` (depth-empty bins; <2% on the native
  grid, so the fine bins are populated via the projection, not empty).

`--dump-bins PREFIX` writes per-bin `chrom start end bes delta_contact on_feature
width marginal` (sorted by BES) for the adaptive grid -- the substrate for the
confound controls and the figures.

**Headline (DN-vs-V1P, full window):** native adaptive rho **-0.169** (p=0.001,
feat_conc 7.5) vs uniform **-0.060** at matched resolution vs floored **-0.055**
(feat_conc ~1.0, the inversion) vs standard 3 kb **-0.103**. Native adaptive is
the most sensitive grid; flooring destroys it; the advantage does **not** transfer
to genome-wide Hi-C (there adaptive ~= uniform -- the basis for running genome-wide
integration on fixed grids).

---

## 4. Confound controls (`partial_corr.py`)

The raw correlation is sensitive to per-bin contact depth and bin width, so every
estimate is checked:

- `rho | depth` -- partial Spearman controlling the bin marginal. Kills a coverage
  artifact (high-BES bins sitting off-target at low depth).
- `rho | depth,width` -- also controls width (over-corrects, since width *is* the
  adaptive mechanism, but a useful skeptic's floor).
- `[premise] rho(BES, depth)` -- if strongly negative, high-BES bins are low-depth
  and the coverage confound is live. **Here it is +0.25 -- the wrong sign for a
  confound**, so the relationship is not coverage-driven.

Depth-controlled, V1P *strengthens* to -0.30 and RCTKO surfaces at -0.168 (its
raw -0.058 was depth-suppressed, not null). Under depth+width, V1P holds at -0.106
while RCTKO fades to -0.037 -- V1P robust, RCTKO real-but-width-dependent. Both
survive CBE masking (V1P -0.132; RCTKO -0.168), so neither rides on the deleted
elements.

---

## 5. The sign depends on the window (read before interpreting)

`partial_corr.py --region` and the figure's in-window partial reveal that the
**sign of the 1D-3D correlation is not uniform across the locus**:

- **core** (chr6:40.8-41.65M, the Vb + RC active blocks): **positive** (BES and
  contact change co-localize at the functional anchors),
- **full** (chr6:39.77-42.66M): **negative** -- the anti-localization the benchmark
  reports lives in the flanks and the intervening inactive chromatin.

This matches Allyn 2025: V1P removes a long-range **stripe** of contacts from Trbv1
to the DJb/RC. The 1D change is punctate (at the Trbv1 promoter and the RC anchors,
where it co-localizes with the contact change -> positive core), while the stripe's
3D change is distributed along the intervening Prss/inactive chromatin that carries
no 1D signal (-> negative full). The benchmark's full-window negative is therefore
the signature of a distributed long-range stripe, not a failure of co-localization
at the regulatory elements. **Run `SIGN_MAP.tsv` before writing any directional
claim.** See `../docs/` and the Snakefile `sign_map` rule.

---

## 6. Reproduce

```
# RUN FROM THE REPO ROOT (all paths resolve from there).
# 1D first (separate): run workflow/ with config_v1 -> diff_*.neglog10p.bw
# then, depth-matched 3D + benchmark + controls + figures:
snakemake -n -s mcc/Snakefile --configfile mcc/config.yaml     # dry-run
snakemake    -s mcc/Snakefile --configfile mcc/config.yaml --cores 8
```

Outputs: `grid_benchmark_{native,floored,native_cbemask}.tsv`, per-window
`partial_*` reports + `SIGN_MAP.tsv`, `figures/fig_1d_vs_3d_*.pdf` (displacement),
`figures/mcc_bearing_*.pdf` (pyGenomeTracks supplement on depth-matched matrices),
`feature_geometry.txt`.

Edit `config.yaml` for cluster paths (BAM dir, chrom.sizes, annotation BEDs).

---

## 7. Script inventory (`scripts/`)

- `validpairs_to_pairs.py` -- valid-pairs BAM -> cooler pairs text
- `build_capture_cools.sh` -- per-condition depth-matched cools (`--downsample`)
- `build_adaptive_segmentation.py` / `build_superbin_grid.py` -- native / floored grids
- `benchmark_grid_vs_standard.py` -- adaptive vs uniform vs standard (`--dump-bins`)
- `partial_corr.py` -- depth/width-controlled partials (`--region`)
- `plot_1d_vs_3d.py` -- 1D-vs-3D displacement figure (in-window partial)
- `make_mcc_bearing_ini.py` -- pyGenomeTracks ini (matrices + BES + 5-track panel)
- `superbin_feature_overlap.py` -- feature-anchoring geometry test
- `run_mcc_grid_benchmark.sh` -- standalone native+floored benchmark driver

---

## 8. Four-track panel: drop-CTCF robustness + V1TxS control

The five-track panel includes CTCF, which is mixed-platform across conditions
(CUT&RUN for DN/RCTKO, ChIP for V1P) and therefore not safe for cross-condition
claims. The four-track panel (RNA +/-, Cohesin, NIPBL) is platform-consistent
and serves two purposes.

**(a) Drop-CTCF robustness check** -- does the core anti-localization survive
removing the mixed-platform track? Rebuild the 1D differential without CTCF, then
re-run the benchmark + partial:

```
# 1. derive the 4-track 1D config from YOUR config_v1.yaml (preserves tuned params)
sed -e 's#samples_sheet:.*#samples_sheet: "config/samples_v1_4track.tsv"#' \
    -e 's#categories:.*#categories: "../categories/v1_4track_panel.yaml"#' \
    -e 's#outdir:.*#outdir: "results_v1_4track"#' \
    workflow/config/config_v1.yaml > workflow/config/config_v1_4track.yaml
echo 'paper_figures: false' >> workflow/config/config_v1_4track.yaml   # skip manuscript-only figures

# 2. run the 1D BEARING pipeline on the 4-track panel (DN/V1P/RCTKO)
snakemake -s workflow/Snakefile --configfile workflow/config/config_v1_4track.yaml ...

# 3. benchmark + partial on the 4-track BES (3D contact side is unchanged)
snakemake -s mcc/Snakefile --configfile mcc/config_4track.yaml --profile ... \
  mcc_results_v1_4track/SIGN_MAP.tsv
```

Compare the core partial to the five-track core (V1P -0.26/-0.33). If it holds,
the anti-localization is clean of the CTCF platform artifact. If it weakens
materially, CTCF was carrying it -- report that.

**(b) V1TxS specificity control** -- V1TxS (transcription stop) changes 1D
(transcription) while leaving 3D architecture flat (Allyn 2025), so a correct
integration must return ~0 co-localization: the "true negative" that converts
"anti-localization exists" to "anti-localization is specific". To add it:

1. Fill the V1TxS bigwig paths in `samples_v1_4track.tsv` (replace the FIXME
   placeholders with the real RNA/Rad21/NIPBL paths from your bigwigs dir / GEO
   manifest) and uncomment the two V1TxS rows.
2. Uncomment `V1TxS: V1TxS` (conditions) and `- [DN, V1TxS]` (comparisons) in
   `mcc/config_4track.yaml`, and the DN-vs-V1TxS line in
   `comparisons_v1_4track.tsv`.
3. The V1TxS capture BAM (`raw_ArimaHTS_S065_V1TxS_rep{1,2}`) is already handled
   by `build_capture_cools.sh` (LABEL map) -- regenerate pair_counts to include it.

Prediction: DN-vs-V1TxS core partial near zero (1D moves, 3D flat), in contrast
to V1P/RCTKO. That is the specificity result.

---

## 9. V1TxS specificity control -- ISOLATED run (do not merge into the 4-track)

The downsample floor is GLOBAL: `floor = min(pooled pairs over all conditions
in the counts file)`. Adding V1TxS to the shared `pair_counts_v1.tsv` would, if
V1TxS is shallower than DN (current floor 8,345,577), lower the floor and
re-downsample DN/V1P/RCTKO -- silently moving the locked V1P (-0.254) and
RCTKO (-0.196) core partials. So V1TxS runs on its OWN counts/floor/outroot.

```
# 1. confirm the V1TxS 1D BES exists (from the 4-track 1D run)
ls workflow/results_v1_4track/pvalue/diff_DN_vs_V1TxS.stats.tsv

# 2. build a DN+V1TxS counts file (V1TxS counted the SAME way: first-in-pair primaries)
BAMDIR=/mnt/isilon/bassing_lab/projects/capture_hic/Brittney_V1/03aln
grep -P '^DN\t' capture_cools/pair_counts_v1.tsv > capture_cools/pair_counts_v1txs.tsv
for REP in 1 2; do
  BAM="$BAMDIR/raw_ArimaHTS_S065_V1TxS_rep${REP}_bs_250_master_valid_pairs.bam"
  N=$(samtools view -c -f 0x40 -F 0x900 "$BAM")
  printf "V1TxS\t%s\t%d\n" "$(basename "$BAM")" "$N" >> capture_cools/pair_counts_v1txs.tsv
done
cat capture_cools/pair_counts_v1txs.tsv

# 3. run the isolated V1TxS benchmark (own outroot, own floor)
snakemake -s mcc/Snakefile --configfile mcc/config_v1txs.yaml --profile workflow/profiles/slurm
head mcc_results_v1txs/SIGN_MAP.tsv
```

Prediction: DN-vs-V1TxS core partial near zero (transcription moves in 1D,
architecture stays flat in 3D), in contrast to V1P/RCTKO -- the specificity
result. Keep DN-vs-V1TxS COMMENTED in comparisons_v1_4track.tsv and
config_4track.yaml so the 4-track run stays V1P/RCTKO-only.
