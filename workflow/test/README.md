# Simulated end-to-end workflow test

Runs the **whole** BEARING Snakemake workflow (`rule all`) on a small
simulated dataset, so you can check that an install, or a refactor, still works
end to end without real data.

```bash
conda activate bearing
bash workflow/test/run_sim_test.sh sim_workflow_test 4    # WORKDIR, CORES
# ... ~6-10 min ...
# PASS: simulated end-to-end workflow
```

## What it does

`make_sim_inputs.sh` simulates 4 samples (conditions **A** and **B**, two
replicates each) of the 6-track panel on one 2 Mb chromosome (`chrSim`).
B's planted blocks are shifted relative to A's, so A vs B has real differences.
It also writes the sample sheet, regions, gene/CBE/element-class BEDs, and
`config.sim.yaml`, which is layered on top of `workflow/config/config.yaml`.

`run_sim_test.sh` then runs all 52 workflow jobs: scoring, provenance check,
compare, region QC (pyGenomeTracks), permutation nulls, per-sample and
differential p-values, floor sweep, calibration, regional and element-class
enrichment, CBE / feature-track changes, decomposition report, the synthetic
recovery benchmark, and the benchmark summary. It then checks the key outputs
and that A vs B has significant bins (the reference run gives 487 bins at
FDR < 0.05, all inside the planted blocks).

## Settings that differ from a real run

| Key | Test value | Why |
|---|---|---|
| `n_perms` | 3 | speed (real runs: 100) |
| `pvalue_min_signal` | 0.05 | small simulated differentials |
| `pvalue_null_min_bytes` | 0 | null qcats are ~350 KB, below the 1 MB genome-wide truncation gate |
| `paper_figures`, `paper_tables` | false | these need the manuscript conditions |
| `regional_loci`, `element_class_*`, `cbe_bed`, `agr_genes_bed` | simulated BEDs | defaults point at mm10 Tcrb/Igh annotations |

## Not covered

The Hi-C phases (TADs, compartments, dcHiC, contact isolation, Hi-C figures)
are not exercised: they only run when the sheet has `cool` paths, and they
expect manuscript inputs (TAD/insulation/PCA/loop files, Tcrb anchors). The
`.cool` loading code itself is covered by the unit tests (`pytest -q`).
