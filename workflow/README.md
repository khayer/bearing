# BEARING workflow

One `snakemake` invocation runs the whole pipeline as a single DAG. Hi-C phases
run only if the sample sheet defines Hi-C conditions, so BEARING works with or
without Hi-C inputs.

## Run

```bash
# local
snakemake -s workflow/Snakefile --configfile workflow/config/config.yaml \
          --cores 16 --use-apptainer

# cluster -- Snakemake submits every rule to SLURM (no manual sbatch):
snakemake -s workflow/Snakefile --configfile workflow/config/config.yaml \
          --profile workflow/profiles/slurm
```

The profile (`profiles/slurm/config.yaml`) sets the executor, partitions
(dbhiq,defq), and default resources; each rule declares its own threads / mem /
runtime, so per-step sizing matches the old hand-written sbatch calls.

## Juicer .hic support (hic-straw)

`.cool` input needs nothing extra. For Juicer `.hic` input, `environment.yml`
installs `hic-straw` (note the hyphen; the bare name `hicstraw` does not exist
on PyPI -- that was the original install failure). It has no prebuilt wheel and
compiles a C++ extension against libcurl, so the env also pulls `cxx-compiler`,
`libcurl`, and `curl` from conda-forge to satisfy the build.

Verify after creating the env:
```bash
conda env create -f environment.yml
conda activate bearing
python -c "import hicstraw; print('hic-straw OK', hicstraw.__file__)"
```
(The import name is `hicstraw`; the install name is `hic-straw`.)

If the compile still fails on your system, install the prebuilt binary instead:
```bash
conda install -c bioconda hic-straw    # verify the exact package name for your platform
```

## Run on SLURM (Snakemake submits each rule as its own sbatch job)

Preflight + launch in one step (run from the directory your sheet's relative
paths resolve against, with the `bearing` conda env active):

```bash
bash workflow/run_cluster.sh             # full pipeline
bash workflow/run_cluster.sh --core-only # BEARING only, skip Hi-C
bash workflow/run_cluster.sh --dry       # preflight + dry run, no submission
```

This (1) runs `preflight.py` to confirm every BigWig / Hi-C / reference file in
the sheet exists -- failing in seconds if one is missing rather than hours into
the run -- then (2) submits via the slurm profile, where Snakemake issues one
`sbatch` per rule instance (partitions dbhiq,defq; per-rule mem/cpu/runtime).

Equivalent manual command:
```bash
python3 workflow/preflight.py --configfile workflow/config/config.yaml
snakemake -s workflow/Snakefile --configfile workflow/config/config.yaml \
          --profile workflow/profiles/slurm
```

To run jobs inside the Apptainer image instead of the active conda env, add
`--use-apptainer` (pull the image first with `workflow/get_container.sh`).

Preflight standalone (any time, before committing to a run):
```bash
python3 workflow/preflight.py --configfile workflow/config/config.yaml
python3 workflow/preflight.py --configfile workflow/config/config.yaml --core-only
# resolve relative paths against a specific data dir:
python3 workflow/preflight.py --configfile workflow/config/config.yaml --base-dir /path/to/run
```

## Staging inputs from S3 (no aws / mount-s3 on the cluster)

If the cluster cannot mount S3 or run `aws s3 sync`, host the inputs behind
HTTPS and pull them with `wget`. Two facts drive the approach:

- `cooler` can read `s3://` directly, but `pyBigWig` (used for every `bw` track
  and the `pca1` compartment tracks) cannot. So a uniform "URLs in the sheet"
  scheme needs HTTPS, not `s3://`.
- BEARING reads each BigWig densely and re-reads them per permutation round, so
  streaming per-bin over the network is slow and costly; stage once to local
  disk instead.

`workflow/stage_inputs.py` does this: put HTTPS URLs (public bucket URLs or
presigned S3 URLs) in the sheet, then

```bash
python3 workflow/stage_inputs.py \
    --sheet workflow/config/samples.tsv \
    --resolutions 10000 25000 100000 250000 500000 \
    --cache-dir resources/staged \
    --out-sheet workflow/config/samples.local.tsv

snakemake ... --config samples_sheet=workflow/config/samples.local.tsv
```

It downloads each remote file with `wget` ONLY if not already in the cache
(idempotent; atomic .part-then-rename so an interrupted pull re-downloads
cleanly), passes local paths through unchanged, expands `{res}` over the given
resolutions, and writes a localized sheet pointing at the cache. `s3://` URLs
are rejected with guidance -- convert to https/presigned first. Presigned URLs
expire, so regenerate them before a run.

## Testing the pipeline before a real run

Two checks, fastest first.

**1. Dry run (seconds) -- validates the DAG, no execution.**
```bash
snakemake -s workflow/Snakefile --configfile workflow/config/config.yaml -n
```
Lints every rule, resolves wildcards, and prints the job table (28 jobs with
Hi-C, fewer without). Catches config/path/rule errors before any compute. Add
`--rulegraph | dot -Tpng -o dag.png` to redraw the DAG.

**2. Smoke test (~40 s) -- runs the real core scripts on synthetic data.**
```bash
bash workflow/smoke_test.sh              # score -> compare (fast)
bash workflow/smoke_test.sh --with-perm  # also permutation null (slower)
```
Generates a tiny 4-sample / 2-condition synthetic dataset (no real data, no
Hi-C, no cluster) and runs bigwig_to_qcat -> compare_qcat end to end, asserting
non-empty qcats and that both the primary and supplementary Spearman outputs
appear. This exercises the actual scripts and the data flow between them; it
does NOT reproduce biology. Use it after any code change to confirm the core
path still runs.

Note: the dry run checks structure; the smoke test checks execution. Run both.

## Configure

Edit `config/config.yaml`. The four parameters under "to CONFIRM" determine
whether the rerun reproduces the published numbers (min_signal, normalize,
n_perms, min_shift). Inputs (`resources/`) are the data files referenced there.

## Sample sheet (`config/samples.tsv`)

Per-replicate BEARING rows carry `bw`; per-condition Hi-C inputs go in the
optional `cool`, `insul`, `pca1` columns as patterns with a `{res}` placeholder
(this absorbs the mixed EBKO/s3T3/EbKO filenames). Hi-C-only conditions
(dV1P, dV1CTCF) have no `bw`. Leave all three Hi-C columns blank to run without
Hi-C; the Hi-C rules then drop out of the DAG automatically.

## Phases (one rule group each)

blacklist (detect + merge + ENCODE) -> score (+ floors) -> compare (primary +
supplementary Spearman) -> perm nulls (observed + differential) -> differential
p-values -> regional enrichment -> [Hi-C: figures, TAD extension, compartments,
contact isolation, cross-locus] -> calibration -> decomposition -> synthetic
benchmark.

Heavy multi-file steps use `.done` sentinels so the DAG orders correctly without
predicting every filename. The differential p-value step wraps the sbatch-free
`run_perm_diff_pvalue.sh`; the old `submit_perm_qcat_pvalue_slurm.sh` is bypassed
because it submitted its own sbatch jobs (which would fight Snakemake).

## Confirm against a live run

The perm-null directory layout (`{out}/perm/perm*`) and the differential file
names (`diff_<A>_vs_<B>.qcat.bgz` / `.stats.tsv`) are taken from your run notes;
verify the first `pvalue_diff` job finds the diff files it expects.
