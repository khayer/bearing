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

> **New here?** [`QUICKSTART.md`](QUICKSTART.md) is a 5-minute, no-data check
> of a local install: run the tests, score a simulated chromosome, and compare
> your plots with the reference images there. CI runs the same path.

## Install

```bash
# conda / Bioconda (recommended)
conda env create -f environment.yml
conda activate bearing          # installs the package + bearing-* commands

# HPC container -- pull the prebuilt image (built in CI, hosted on GHCR)
# requires apptainer OR singularity on PATH
bash workflow/get_container.sh          # writes ./bearing.sif
# or build it yourself:
apptainer build bearing.sif Apptainer.def
```

`pip install -e .` exposes console entry points: `bearing-score`,
`bearing-compare`, `bearing-pvalue`, `bearing-regional`, `bearing-perm`,
`bearing-shift`, `bearing-rebin`. For tests, use `pip install -e ".[test]"`
(the conda environment already has everything) and run `pytest -q`.

## Run the pipeline

The end-to-end pipeline is a Snakemake workflow over five steps; edit
`workflow/config/config.yaml` and:

```bash
snakemake -s workflow/Snakefile --configfile workflow/config/config.yaml --cores 16
```

See `workflow/README.md` for the step-by-step breakdown and
`paper/reproduce_all.sh` for the explicit staged-SLURM equivalents.

### Two different `--min-signal` flags (do not conflate)

Two scripts expose a flag named `--min-signal`, with **different meanings and
different defaults**. They are intentionally not renamed (renaming would break
the CLI); know which is which:

| Flag | Meaning | Default | Workflow config key |
|---|---|---|---|
| `bigwig_to_qcat.py --min-signal` | per-sample **raw-signal** floor; gates scoring | `0.01` | `min_signal` (set to `0.1`) |
| `bearing_pvalue.py --min-signal` | **\|score\|** floor on the differential; lower bins are not tested | `0.5` (KL), `0.05` (JSD) | `pvalue_min_signal` |

The workflow sets `pvalue_min_signal` explicitly so results stay reproducible
even if the library default changes.

### Composition prior for weak bins

The score stage also supports an optional composition prior via
`bigwig_to_qcat.py --prior-strength` / `config prior_strength`. The default is
`0.0` (off, identical to the legacy scoring path). When `prior_strength > 0`,
each bin is shrunk toward the background composition `Q` before scoring:

```text
P = (x + alpha * Q) / (sum(x) + alpha)
```

This smooths the hard cliff created by the raw-signal floor for barely-detectable
bins while preserving bins with real signal; it is intended for sensitivity
scans rather than the default production run. The workflow records the setting
via `score_provenance` and passes it to the permutation-null scoring, so
observed and null scores stay on the same scale.

Set it in the config (`prior_strength: 1.0`) or on the command line:

```bash
python bigwig_to_qcat.py \
  --bw sample1.bw sample2.bw \
  --sample-name sample \
  --categories categories/mm10_6track_panel.yaml \
  --chrom-sizes chrom.sizes --out sample.qcat.bgz \
  --min-signal 0 --prior-strength 1.0
```

A dedicated workflow helper is included for a focused prior sweep:

```bash
bash workflow/run_cluster_prior.sh --alpha 1 --n-perms 10
```

This writes a derived config under `workflow/config/` and reruns the default
p-value path with the prior on.

## Plotting entry points

- `bearing_hic_combined_plot.py` -- **two-condition** comparison figure
  (Hi-C + per-track differential + p-value overlay), e.g. DN vs DP. This is the
  high-level entry point for the main per-comparison figures.
- `compare_qcat.py` -- plots **all conditions** over a given set of regions, and
  produces the differential qcat, PCA, Spearman, and JSD diagnostics.
- `bearing_hic_plot.py` -- the lower-level **single-region** plotter that the
  combined figure builds on. Use it directly to render one region from explicit
  `--contact`/`--qcat`/`--diff-qcat` inputs; this is the script used in the
  reviewer walkthrough in `docs/working_examples.md`.

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

bearing_hic_*_plot.py, bearing_kl_track_plot.py, batch_*.py
                         figure scripts (see "Plotting entry points" above)

convert_categories.py    convert a JSON/YAML categories file to the standard categories YAML
fetch_geo_bigwigs.py     fetch processed bigWigs from a GEO series into a BEARING-ready
                         layout (writes a manifest + a resumable download script)

hic/                     Hi-C analysis (Figs 3, 4, S7)
  tad_extension_analysis.py, tcrb_contact_isolation.py, compartment_analysis.py,
  bes_hic_crosslocus.py, define_wide_loci.py, per_locus_summary.py, top_diff_sites.py
benchmark/               synthetic-data recovery benchmark (see Validation)
workflow/                Snakemake pipeline + config
categories/              track-panel definitions (mm10_6track_panel.yaml is canonical)
annotations/             region BEDs (Tcrb, Igh, feature sets)
paper/                   figure -> code reproducibility map
docs/                    scoring method, statistics/significance, roadmap, run notes
slurm/                   HPC submission scripts
tests/                   unit + integration tests (tests/data/: small simulated fixtures)
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
