# Working Examples

This file collects runnable command templates for the main workflows in the repository.

## Reviewer Smoke Test

If a reviewer wants to test the code on their own data, the simplest path is
to start from their own BigWig inputs, score them with `bigwig_to_qcat.py`,
and then run the comparison workflow on a single locus.

1. Start from `reviewer_samples.tsv`, which is a small DN vs 3T3 BigWig sheet
  with paths that should be replaced by the reviewer’s own files.
2. Run `bigwig_to_qcat.py` on that sheet to generate qcat files.
3. Reuse the same sheet with `compare_qcat.py`, which will convert the BigWigs
  on the fly if qcat files are not already present.
4. Use `reviewer_regions.template` and restrict the run to `chr2` so the
  reviewer only exercises the Rag locus and one other locus on the same
  chromosome.
5. Run `bearing_hic_plot.py` only if the reviewer also has Hi-C contacts for
  the same conditions and wants to verify the final figure path.

Example scoring command:

```bash
python bigwig_to_qcat.py \
  --sheet reviewer_samples.tsv \
  --out reviewer_placeholder.qcat.bgz \
  --jobs 2
```

Example comparison command:

```bash
python compare_qcat.py \
  --sheet reviewer_samples.tsv \
  --out reviewer_comparison \
  --regions-file reviewer_regions.template \
  --chroms chr2 \
  --consensus-q \
  --workers 2
```

Example single-region figure command:

```bash
python bearing_hic_plot.py \
  --contact-a DN.cool \
  --contact-b 3T3.cool \
  --qcat-a DN_rep1.qcat.bgz \
  --qcat-b 3T3_rep1.qcat.bgz \
  --diff-qcat reviewer_comparison/diff_DN_vs_3T3.qcat.bgz \
  --region chr2:50700000-50950000 \
  --resolution 10000 \
  --out reviewer_dn_vs_3t3_tcrb.pdf
```

Expected outputs to check:
- `reviewer_comparison/diff_DN_vs_3T3.qcat.bgz`
- `reviewer_comparison/jsd_heatmap.pdf`
- `reviewer_comparison/spearman.pdf`
- `reviewer_comparison/pca.pdf`
- `reviewer_comparison/regions/rag_compare.pdf` and `hoxd_compare.pdf` when
  region plotting is enabled

## 1. Score BigWig tracks into qcat files

Use the sheet-driven scoring path when you have one row per sample:

```bash
python bigwig_to_qcat.py \
  --sheet samples_template.tsv \
  --out placeholder.qcat.bgz \
  --jobs 8
```

If you are using a blacklist, add it to the scoring run:

```bash
python bigwig_to_qcat.py \
  --sheet samples_template.tsv \
  --out placeholder.qcat.bgz \
  --blacklist detected_blacklist_merged.bed \
  --jobs 8
```

## 2. Compare qcat profiles across samples

This produces JSD, Spearman, PCA, differential qcat files, and the comparison INI:

```bash
python compare_qcat.py \
  --sheet samples_template.tsv \
  --out comparison \
  --regions-file regions_template.tsv \
  --consensus-q \
  --workers 4
```

If you want the INI to use the original qcat values rather than plotting-only clipped copies, add `--no-clip`.

## 3. Generate permutation nulls

Run the permutation pipeline after the observed qcat files exist:

```bash
python generate_perm_nulls.py \
  --sheet samples_template.tsv \
  --diff-sheet samples_template.tsv \
  --compare-consensus-q \
  --n-perms 5 \
  --out-dir perm
```

## 4. Compute p-values from observed and null qcats

Example for a single sample:

```bash
python bearing_pvalue.py \
  --qcat WT_rep1.qcat.bgz \
  --null-qcat perm/perm1/WT_rep1/null_perm1.qcat.bgz perm/perm2/WT_rep1/null_perm2.qcat.bgz \
  --out-prefix results/WT_rep1 \
  --fdr 0.05 \
  --score-plot
```

Example for a differential qcat:

```bash
python bearing_pvalue.py \
  --qcat comparison/diff_WT_vs_KO.qcat.bgz \
  --null-qcat perm/perm1/diff_comparison/diff_WT_vs_KO.qcat.bgz perm/perm2/diff_comparison/diff_WT_vs_KO.qcat.bgz \
  --diff \
  --out-prefix results/diff_WT_vs_KO \
  --fdr 0.05 \
  --score-plot
```

## 5. Plot a single Hi-C region

`bearing_hic_plot.py` is the single-region plotter (the `bearing_hic_combined_plot.py`
two-condition entry point in the README builds on it). Use it directly when you
want to render one region from explicit inputs:

```bash
python bearing_hic_plot.py \
  --contact-a WT.cool \
  --contact-b KO.cool \
  --qcat-a WT_rep1.qcat.bgz \
  --qcat-b KO_rep1.qcat.bgz \
  --diff-qcat comparison/diff_WT_vs_KO.qcat.bgz \
  --pval-a results/WT_rep1.neglog10p.bw \
  --pval-b results/KO_rep1.neglog10p.bw \
  --region chr6:40793981-41688054 \
  --resolution 10000 \
  --out figure.pdf
```

## 6. Batch Hi-C plotting across regions

```bash
python batch_bearing_hic_plots.py \
  --sheet samples_template.tsv \
  --regions-file regions_template.tsv \
  --reference-condition WT \
  --contact WT=WT.cool --contact KO=KO.cool \
  --results-dir results \
  --comparison-dir comparison \
  --outdir hic_batch \
  --run
```

## Validation

A fast regression check for the compare pipeline is:

```bash
python -m unittest tests.test_compare_qcat
```

If you are changing the scoring path, also run the matching unit tests for the affected script before publishing.

## Utilities

Two helper scripts live at the repository root and are used outside the core
pipeline.

`fetch_geo_bigwigs.py` pulls processed bigWig tracks from a GEO series and
organizes them into a BEARING-ready layout. It reads the GEO SOFT file,
classifies each sample by condition/assay/strand from its title, and writes a
manifest plus a resumable download script. Run it where there is internet (your
cluster), inspect the manifest, then download:

```bash
python3 fetch_geo_bigwigs.py --gse GSE296315 --out geo_v1 \
    --conditions WT,Rag1,R1KO,513TKO,RCTKO,V1PKO,V1TxS \
    --assays CTCF,RAD21,NIPBL,RNA,GRO,H3K27ac,H3K4me3,Pol2
# inspect geo_v1/geo_bigwig_manifest.tsv, then:
bash geo_v1/download_bigwigs.sh
```

`convert_categories.py` converts a JSON or YAML categories file into the
standard categories YAML consumed by the scoring and plotting scripts:

```bash
python convert_categories.py input.json output.yaml
```
