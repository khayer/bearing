# Synthetic benchmark (simulated data with known ground truth)

**Everything in this directory operates on SIMULATED data, not experimental
data.** No biological sample, BigWig, or Hi-C matrix from the manuscript is used
here. The point of the synthetic benchmark is to test BEARING against a *known*
per-bin composition that the method must recover, which biological
recapitulation alone cannot provide.

## Why this exists

Showing that BEARING reproduces known Tcrb/Igh biology does not, by itself,
rule out the objection that the method was tuned to recover what was already
known to be there. The synthetic benchmark answers that objection directly: we
plant a known per-bin probability vector

```
P = (1 - s) * Q + s * a
```

where `Q` is the background composition, `a` is a known allocation over the six
tracks (ATAC, RNA+, RNA-, CTCF, RAD21/Cohesin, H3K27ac), and `s` is a known
deviation strength. Null bins have `P = Q` and should score low. The real
pipeline is then run on the simulated BigWigs and asked to recover `a`.

The simulator models realistic per-bin coverage (heavy-tailed depth), a
signal-to-noise scalar, and negative-binomial track noise, so recovery is
tested under realistic noise rather than on a noise-free toy.

## What it validates

1. **Recovery of planted composition.** Does the recovered per-track allocation
   match the planted `a`? Because BEARING's per-track score is a *clamped* KL
   contribution (one-sided enrichment), recovery is scored on the enriched
   subset (tracks with `a_i > Q_i`); suppressed tracks are reported as
   expected-zero, not error. Compositional error is reported as a CLR
   (Aitchison) distance in addition to the Pearson/Spearman reviewers expect.
2. **Parameter stability.** A sweep over the three free parameters (low-signal
   mask, zero-clamp pseudocount, normalization quantile) shows whether a
   quantitative claim such as the cohesin/CTCF compositional split is a real
   feature or an artifact of those choices.

## Contents

- `simulate_bearing_tracks.py`  ground-truth generator (writes synthetic BigWigs + a truth-table BED)
- `evaluate_bearing_recovery.py`  runs the real pipeline and scores recovery
- `plot_bearing_sweep.py`  renders the SNR x dispersion x pseudocount sweep
- `visualize_bearing_tracks.py`  track-stack sanity view (planted blocks vs BEARING score)
- `run_bearing_recovery_sweep.sh`  driver for the full coarse-grid sweep

## Run

```bash
# 1) generate synthetic tracks with a planted truth table
python3 simulate_bearing_tracks.py --out-prefix bearing_sim

# 2) run the real pipeline on them and score recovery
python3 evaluate_bearing_recovery.py --prefix bearing_sim --out-prefix recovery

# 3) full parameter sweep + figures
bash run_bearing_recovery_sweep.sh
python3 plot_bearing_sweep.py --master sweep_master.tsv --out-prefix sweep
```
