# BEARING Quickstart (5 minutes, no data required)

This runs the **real** BEARING scoring pipeline end-to-end on a small simulated
chromosome, so a reviewer can confirm the code works without downloading any
BigWig, Hi-C, or GEO files. It is the same path exercised by continuous
integration (`.github/workflows/test.yml`).

## 1. Install (pip is enough for the core scoring path)

```bash
git clone https://github.com/khayer/bearing.git
cd bearing
pip install -e . pytest        # core scientific stack + entry points
```

For the Hi-C figure scripts and the full workflow you need the conda
environment instead (adds hicexplorer, pygenometracks, snakemake, bgzip/tabix,
and pytest):

```bash
conda env create -f environment.yml
conda activate bearing        # includes pytest; `pytest -q` works immediately
```

## 2. Run the test suite (~6 s)

```bash
pytest -q
# expected: 62 passed
```

## 3. Simulate a tiny chromosome with known ground truth

```bash
python benchmark/simulate_bearing_tracks.py \
  --prefix sim --chrom chrSim --length 2000000 --bin-size 200 --seed 1
```

This writes six synthetic BigWigs (`sim_atac.bw`, `sim_rnaplus.bw`,
`sim_rnaminus.bw`, `sim_ctcf.bw`, `sim_rad21.bw`, `sim_h3k27ac.bw`), a
`sim.chrom.sizes`, and a `sim_truth.tsv` recording the planted composition.

## 4. Score the BigWigs into a qcat file

```bash
python bigwig_to_qcat.py \
  --bw sim_atac.bw sim_rnaplus.bw sim_rnaminus.bw \
       sim_ctcf.bw sim_rad21.bw sim_h3k27ac.bw \
  --sample-name sim \
  --categories categories/mm10_6track_panel.yaml \
  --chrom-sizes sim.chrom.sizes --chroms chrSim \
  --out sim.qcat.bgz --jobs 2
```

## 5. Confirm the output has real (non-zero) scores

```bash
test -s sim.qcat.bgz
zcat sim.qcat.bgz | head -2
# each line: chrSim  start  end  id:N,qcat:[[score,track], ...]
```

That is the full core loop: **simulate -> score -> inspect**. Everything else
in the repository builds on the qcat files produced in step 4.

---

## Next steps on real data

Once the synthetic run works, the same commands apply to real BigWigs via a
sample sheet. See `docs/working_examples.md` for the full templates. The short
version:

```bash
# score a whole sample sheet
python bigwig_to_qcat.py --sheet samples_template.tsv --out placeholder.qcat.bgz --jobs 8

# multi-sample QC + differential qcat (JSD, Spearman, PCA, diff files)
python compare_qcat.py --sheet samples_template.tsv --out comparison \
  --regions-file regions_template.tsv --consensus-q --workers 4

# empirical permutation p-values
python generate_perm_nulls.py --sheet samples_template.tsv \
  --diff-sheet samples_template.tsv --compare-consensus-q --n-perms 5 --out-dir perm
python bearing_pvalue.py --qcat comparison/diff_WT_vs_KO.qcat.bgz \
  --null-qcat perm/perm1/diff_comparison/diff_WT_vs_KO.qcat.bgz \
  --diff --out-prefix results/diff_WT_vs_KO --fdr 0.05 --score-plot
```

## Full pipeline (Snakemake)

Edit `workflow/config/config.yaml`, then:

```bash
snakemake -s workflow/Snakefile --configfile workflow/config/config.yaml --cores 16
```

## Reproducing the manuscript figures

`paper/README.md` maps every figure and table to the exact script(s) that
produce it. Run the core workflow first to generate the qcat / differential /
p-value inputs those scripts consume.
