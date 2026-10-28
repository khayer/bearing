<!--
  TRY_BEARING.md - one-page "install and test" guide for lab testers
  Authors: Katharina E. Hayer with Claude (Anthropic), co-created
  Version: 1.1
  Date: 2026-10-28
-->

# Try BEARING on the practice data (about 10 minutes)

Thanks for testing! The goal: install BEARING on your laptop (or cluster
account), run it on the built-in **practice dataset**, and report how it
went. The practice data are simulated, so no downloads are needed. You do not
need to know anything about the biology.

## What you need

- macOS or Linux. On **Windows, use WSL2** (Ubuntu); pysam does not install
  on native Windows.
- Python 3.10 or newer (`python3 --version`), git, and about 200 MB of disk.

## 1. Install

```bash
git clone https://github.com/khayer/bearing.git
cd bearing
python3 -m venv .venv && source .venv/bin/activate
pip install -e . pytest
```

If `pip` fails while building **pyBigWig** or **pysam**, use conda instead:

```bash
conda create -n bearing -c conda-forge -c bioconda python=3.11 pybigwig pysam pytest
conda activate bearing
pip install -e .
```

## 2. Check the install (about 20 s)

```bash
pytest -q
```

Expected: every test passes (181 passed as of 2026-10-28), with no failures.

## 3. Run the practice dataset (about 1 min)

```bash
BEARING_OUT=bearing_test bash workflow/smoke_test.sh --with-perm
```

This simulates 2 conditions x 2 replicates of a 6-track panel (ATAC,
RNA +/-, CTCF, RAD21, H3K27ac), then runs the real pipeline:
**score -> compare -> permutation null -> p-values**.

Expected last lines:

```
PASS: BEARING core path ran end to end on synthetic data.
```

Outputs are written to `bearing_test/`. The most useful ones:

| File | What it is |
|------|------------|
| `A_rep1.qcat.bgz` | per-bin, per-track BEARING scores (one per sample) |
| `compare/pca.pdf`, `compare/spearman.pdf` | sample QC |
| `diff_A_vs_B.stats.tsv` | differential test: p, BH q, direction, per-track contribution |

## 4. Optional: look at a picture

```bash
cd bearing_test
python ../benchmark/visualize_bearing_tracks.py --prefix A_rep1 \
  --qcat A_rep1.qcat.bgz --out A_rep1_tracks.png
```

The grey boxes are the planted changes. The BEARING score should rise inside
them, and the per-track colours show which assay drives each one.

## 5. Report back

Either way, please report: success or failure, both are useful.
Open a GitHub issue (template **"Test report"**) at
https://github.com/khayer/bearing/issues or email the maintainer (Katharina Hayer). Include:

1. Your OS and chip (for example, macOS 14 / Apple M2, Ubuntu 22.04, WSL2).
2. The output of:
   ```bash
   python -c "import sys,platform,numpy,pyBigWig,pysam; print(platform.platform(), sys.version.split()[0], numpy.__version__, pysam.__version__)"
   ```
3. Install route (pip or conda) and roughly how long each step took.
4. The `pytest -q` summary line.
5. Whether you saw `PASS` in step 3. If not, paste the last ~30 lines of output.
6. Anything confusing in this page.

Real data come later. This round only checks that BEARING installs and runs
for people other than its author.
