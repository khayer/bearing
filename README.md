# BEARING

![BEARING logo](bearing_logo_v2.svg)

**B**in-level **E**pigenomic **A**ttribution **R**eference-model-free
**I**ntegrated **N**ormalized **G**enomic score.

BEARING converts continuous BigWig signal tracks into per-bin clamped-KL score
vectors, attributes each per-bin score to its component assays, compares those
vectors across conditions, and renders publication figures integrating Hi-C,
per-track decompositions, p-value overlays, and gene annotations. It requires
no pre-trained chromatin-state model and no replicates (significance comes from
an empirical permutation null).

## Install

```bash
# conda / Bioconda (recommended)
conda env create -f environment.yml
conda activate bearing          # installs the package + bearing-* commands

# HPC container -- pull the prebuilt image (built in CI, hosted on GHCR)
bash workflow/get_container.sh          # writes ./bearing.sif
# requires apptainer or singularity on PATH
# or build it yourself:
apptainer build bearing.sif Apptainer.def
```

`pip install -e .` exposes console entry points: `bearing-score`,
`bearing-compare`, `bearing-pvalue`, `bearing-regional`, `bearing-perm`,
`bearing-shift`, `bearing-rebin`.

## Run the pipeline

The end-to-end pipeline is a Snakemake workflow over five steps; edit
`workflow/config/config.yaml` and:

```bash
snakemake -s workflow/Snakefile --configfile workflow/config/config.yaml --cores 16
```

See `workflow/README.md` for the step-by-step breakdown and `docs/how_I_run_it.md`
for the manual / SLURM equivalents.

## Two plotting entry points

- `bearing_hic_combined_plot.py` -- **two-condition** comparison figure
  (Hi-C + per-track differential + p-value overlay), e.g. DN vs DP. This is the
  entry point for the main per-comparison figures.
- `compare_qcat.py` -- plots **all conditions** over a given set of regions, and
  produces the differential qcat, PCA, Spearman, and JSD diagnostics.

## Repository map

```
bigwig_to_qcat.py        BigWig -> qcat scoring (per-bin clamped KL)
compare_qcat.py          multi-sample QC + differential qcat (all conditions)
bearing_pvalue.py        empirical permutation p-values (--diff for differentials)
regional_enrichment.py   pre-specified regional enrichment test
consolidate_regional_enrichment.py / render_regional_enrichment_heatmap.py
generate_perm_nulls.py   permutation-null orchestration
shift_bigwig.py          circular-shift null + replicate-shift control
rebin_qcat.py            multi-resolution rendering
bearing_calibration.py   within-condition replicate-differential FDR calibration
diffuse_architecture_top1pct.py   per-track diffuse-vs-focal quantification
bes_hic_correlation.py         within-locus BES-Hi-C co-localization

bearing_hic_combined_plot.py      two-condition figure (entry point)
bearing_hic_plot.py / *_pval_overlay.py / *_triangle.py / bearing_hic_kl_track_plot.py
bearing_kl_track_plot.py          focal-locus per-track tracks
batch_bearing_hic_plots.py / batch_pygenometracks.py

hic/                     Hi-C analysis (Figs 3, 4, S7)
  tad_extension_analysis.py, tcrb_contact_isolation.py, compartment_analysis.py,
  bes_hic_crosslocus.py, define_wide_loci.py, per_locus_summary.py, top_diff_sites.py
benchmark/               SYNTHETIC-DATA recovery benchmark (simulated, see its README)
workflow/                Snakemake pipeline + config
categories/              track-panel definitions (mm10_6track_panel.yaml is canonical)
annotations/             region BEDs (Tcrb, Igh, feature sets)
paper/                   figure -> code reproducibility map
docs/                    scoring method, roadmap, run notes
slurm/                   HPC submission scripts
tests/                   unit + integration tests
bearing/                 shared library modules
```

## Validation

BEARING is validated two ways: (1) a **synthetic benchmark** on simulated data
with known per-bin ground truth (`benchmark/`), which tests whether the method
recovers a planted composition and is stable to its free parameters; and (2)
biological recapitulation on real data, anchored to the previously characterized
DN-vs-Eb-knockout Tcrb comparison and extended to new sample types and the Igh
locus.

## Data availability

Processed BigWig/qcat files and Hi-C matrices: see manuscript Data Availability.
Annotation BEDs used in analyses are in `annotations/`.
