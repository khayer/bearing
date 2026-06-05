#!/usr/bin/env python3
"""
shift_bigwig.py -- Generate circularly-shifted BigWig files for the Bearing
                   permutation null.

Each input BigWig is shifted by a random (or specified) number of base pairs
along each chromosome. The signal wraps around at chromosome ends so the
marginal distribution of signal values is preserved exactly -- only the
spatial arrangement changes.

Purpose: break the co-localisation of signal across tracks at biologically
meaningful loci while preserving the per-track signal distribution. The
resulting BigWigs are used as input to bigwig_to_qcat.py to generate a
permuted qcat.bgz, which is then passed to bearing_pvalue.py via
--null-qcat to compute empirical p-values.

USAGE
-----
  # Single file, auto-random shift:
  python shift_bigwig.py --bw atac.bw --out atac_perm1.bw

  # All files in one command, independent random shifts per file:
  python shift_bigwig.py \\
    --bw atac.bw ctcf.bw rad21.bw h3k27ac.bw rnapos.bw rnaneg.bw \\
    --out-dir permuted/ \\
    --suffix _perm1

  # Sheet mode -- mirrors the samples.tsv format used by bigwig_to_qcat.py:
  python shift_bigwig.py \\
    --sheet samples.tsv \\
    --out-dir permuted/ \\
    --suffix _perm1
  # This shifts all BigWigs for every sample in the sheet, writes shifted
  # files to permuted/<sample>/, and prints the bigwig_to_qcat.py commands
  # needed to generate the null qcats.

  # Multiple permutations for a denser null (sheet mode):
  for i in 1 2 3; do
    python shift_bigwig.py --sheet samples.tsv --out-dir perm${i}/ --suffix _perm${i}
    # then run the printed bigwig_to_qcat.py commands, or use --run-qcat
  done
  python bearing_pvalue.py \\
    --qcat DN_rep1.qcat.bgz \\
    --null-qcat perm1/DN_rep1/null.qcat.bgz perm2/DN_rep1/null.qcat.bgz \\
    --out-prefix results/DN_rep1 \\
    --score-plot

SAMPLE SHEET FORMAT
-------------------
TSV with header row. Matches the format accepted by bigwig_to_qcat.py:

  sample    bw                                    out
  DN_rep1   atac.bw,ctcf.bw,rad21.bw,h3k27ac.bw  DN_rep1.qcat.bgz
  DN_rep2   atac2.bw,ctcf2.bw,...                 DN_rep2.qcat.bgz

Columns:
  sample  : sample label (used for output subdirectory name)
  bw      : comma-separated BigWig paths (or bw1, bw2, ... columns)
  out     : original qcat output name (optional; used to name the null qcat)

DESIGN NOTES
------------
- Each track is shifted independently by a random amount drawn uniformly
  from [min_shift, chrom_len - min_shift]. This breaks multi-track
  co-localisation while preserving single-track autocorrelation.
- The shift amount is the same across all chromosomes within one file
  to preserve inter-chromosomal relationships. Use --per-chrom-shift
  for independent shifts per chromosome.
- Bins without signal in the original BigWig remain empty after shifting.
- The output BigWig uses the same chromosome sizes as the input.

DEPENDENCIES
------------
  pip install pyBigWig numpy
"""

import argparse
import concurrent.futures
import math
import sys
from pathlib import Path

import numpy as np
import pyBigWig

from bearing.sheet import extract_bw_values, sample_name_from_row
from bearing.tsv import read_tsv_table


def _read_binned_means(bw_in, chrom, chrom_len, n_bins, bin_size):
    """Read chromosome means as bin averages.

    Uses a single bulk values() read plus a numpy reshape-mean, which is ~10x
    faster than asking pyBigWig for nBins separate bin means via stats(): the
    stats(nBins=...) path issues one summary query per bin and dominates the
    runtime on dense genome-wide tracks. Falls back to stats(), then to a
    per-bin loop, for robustness across pyBigWig builds.
    """
    try:
        vals = bw_in.values(chrom, 0, chrom_len, numpy=True)
        vals = np.nan_to_num(vals, nan=0.0).astype(np.float32)
        pad = (-len(vals)) % bin_size
        if pad:
            vals = np.concatenate([vals, np.zeros(pad, dtype=np.float32)])
        binned = vals.reshape(-1, bin_size).mean(axis=1).astype(np.float32)
        # reshape yields ceil(len/bin_size) bins; align to n_bins exactly
        if len(binned) >= n_bins:
            return binned[:n_bins]
        out = np.zeros(n_bins, dtype=np.float32)
        out[:len(binned)] = binned
        return out
    except Exception:
        pass
    try:
        vals = bw_in.stats(chrom, 0, chrom_len, type="mean", nBins=n_bins)
        arr = np.asarray(vals, dtype=np.float32)
        arr[np.isnan(arr)] = 0.0
        return arr
    except Exception:
        # Conservative fallback keeps behavior compatible with older pyBigWig builds.
        out = np.zeros(n_bins, dtype=np.float32)
        for i in range(n_bins):
            s = i * bin_size
            e = min(s + bin_size, chrom_len)
            try:
                v = bw_in.stats(chrom, s, e, type="mean")[0]
                out[i] = 0.0 if (v is None or math.isnan(v)) else float(v)
            except Exception:
                out[i] = 0.0
        return out


def circular_shift_bigwig(in_path, out_path, shift=None, bin_size=200,
                           min_shift=1_000_000, seed=None,
                           per_chrom_shift=False):
    """
    Read a BigWig at bin_size resolution, apply a circular shift per
    chromosome, and write a new BigWig.

    Parameters
    ----------
    in_path    : str or Path -- input BigWig
    out_path   : str or Path -- output BigWig
    shift      : int or None -- fixed shift in bp; if None, random
    bin_size   : int         -- bin resolution (should match scoring pipeline)
    min_shift  : int         -- minimum random shift (avoids near-zero shifts)
    seed       : int or None -- random seed for reproducibility
    per_chrom_shift : bool   -- if True, use a different random shift per chrom
    """
    rng = np.random.default_rng(seed)
    fixed_shift_bins = None

    with pyBigWig.open(str(in_path)) as bw_in:
        chrom_sizes = bw_in.chroms()
        if not chrom_sizes:
            sys.exit(f"ERROR: no chromosomes in {in_path}")

        header = [(c, s) for c, s in chrom_sizes.items()]

        with pyBigWig.open(str(out_path), "w") as bw_out:
            bw_out.addHeader(header)

            for chrom, chrom_len in chrom_sizes.items():
                n_bins = math.ceil(chrom_len / bin_size)
                if n_bins == 0:
                    continue

                # Read all bins for the chromosome in one call for speed.
                vals = _read_binned_means(bw_in, chrom, chrom_len, n_bins, bin_size)

                # Determine shift for this chromosome
                if shift is not None:
                    chrom_shift_bins = (shift // bin_size) % n_bins
                elif per_chrom_shift:
                    max_shift = max(1, n_bins - min_shift // bin_size)
                    min_shift_bins = min_shift // bin_size
                    chrom_shift_bins = int(rng.integers(
                        min_shift_bins, max(min_shift_bins + 1, max_shift)))
                else:
                    # Reuse one random shift for this file across chromosomes.
                    if fixed_shift_bins is None:
                        max_s = max(1, n_bins - min_shift // bin_size)
                        min_s = min_shift // bin_size
                        fixed_shift_bins = int(
                            rng.integers(min_s, max(min_s + 1, max_s)))
                    chrom_shift_bins = fixed_shift_bins

                # Apply circular shift
                shifted = np.roll(vals, chrom_shift_bins)

                # Write back to BigWig (skip zero-signal bins for efficiency)
                nonzero = np.where(shifted != 0.0)[0]
                if len(nonzero) == 0:
                    continue

                starts = (nonzero * bin_size).tolist()
                ends   = [min(s + bin_size, chrom_len) for s in starts]
                chroms = [chrom] * len(starts)
                values = shifted[nonzero].tolist()

                bw_out.addEntries(chroms, starts, ends=ends, values=values)


def load_sheet(path):
    """
    Parse a samples TSV in the same format as bigwig_to_qcat.py --sheet.

    Returns a list of dicts:
      { "sample": str, "bw_paths": [Path, ...], "out": str or None }
    """
    samples = []
    fields, rows = read_tsv_table(Path(path))
    if "bw" not in fields and not any(
        h.startswith("bw") and h[2:].isdigit() for h in fields
    ):
        sys.exit(
            "ERROR: sheet must have a 'bw' column (comma-separated paths) or bw1, bw2, ... columns."
        )

    for row in rows:
        sample = sample_name_from_row(row)
        out = row.get("out") or None
        bw_list = extract_bw_values(row)

        if not bw_list:
            print(
                f"  WARNING: no BigWig files listed for sample '{sample}', skipping.",
                file=sys.stderr,
            )
            continue

        samples.append({
            "sample": sample,
            "bw_paths": [Path(p) for p in bw_list],
            "out": out,
        })

    if not samples:
        sys.exit("ERROR: no usable rows found in sheet.")
    return samples


def _null_qcat_name(out_field, sample, suffix):
    """Derive the null qcat filename from the original out field or sample name."""
    if out_field:
        base = out_field
        for ext in (".qcat.bgz", ".bgz"):
            if base.endswith(ext):
                base = base[:-len(ext)]
                break
        return Path(base).name + suffix + ".qcat.bgz"
    elif sample:
        return sample + suffix + ".qcat.bgz"
    return "null" + suffix + ".qcat.bgz"


def main():
    ap = argparse.ArgumentParser(
        description="Circularly shift BigWig files for Bearing permutation null.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # Input: explicit files OR sheet
    input_grp = ap.add_mutually_exclusive_group(required=True)
    input_grp.add_argument("--bw", nargs="+", metavar="FILE",
                           help="Input BigWig file(s).")
    input_grp.add_argument("--sheet", metavar="FILE",
                           help="Sample sheet TSV (same format as "
                                "bigwig_to_qcat.py --sheet). Each row's "
                                "BigWigs are shifted independently and "
                                "written to out-dir/<sample>/.")

    ap.add_argument("--out", metavar="FILE", default=None,
                    help="Output path for a single file. "
                         "Mutually exclusive with --out-dir and --sheet.")
    ap.add_argument("--out-dir", metavar="DIR", default=None,
                    help="Output directory. Required with --sheet. "
                         "For --bw, files are named <stem><suffix>.bw.")
    ap.add_argument("--suffix", default="_perm",
                    help="Suffix appended to output filenames and null qcat "
                         "names (default: _perm). Change per permutation "
                         "round, e.g. _perm1, _perm2.")
    ap.add_argument("--shift", type=int, default=None, metavar="BP",
                    help="Fixed circular shift in base pairs. "
                         "If omitted, a random shift >= --min-shift is used.")
    ap.add_argument("--min-shift", type=int, default=1_000_000, metavar="BP",
                    help="Minimum random shift in bp (default: 1,000,000).")
    ap.add_argument("--bin-size", type=int, default=200, metavar="BP",
                    help="Bin resolution matching the scoring pipeline "
                         "(default: 200 bp).")
    ap.add_argument("--seed", type=int, default=None,
                    help="Master random seed for reproducibility.")
    ap.add_argument("--per-chrom-shift", action="store_true",
                    help="Use an independent random shift per chromosome.")
    ap.add_argument("--workers", type=int, default=1, metavar="N",
                    help="Parallel workers for shifting BigWig files (default: 1).")
    args = ap.parse_args()

    # Validate output args
    if args.bw:
        if args.out and len(args.bw) > 1:
            sys.exit("ERROR: --out can only be used with a single --bw file.")
        if args.out and args.out_dir:
            sys.exit("ERROR: --out and --out-dir are mutually exclusive.")
        if not args.out and not args.out_dir:
            sys.exit("ERROR: specify --out (single file) or --out-dir.")
    else:  # --sheet
        if not args.out_dir:
            sys.exit("ERROR: --sheet requires --out-dir.")
        if args.out:
            sys.exit("ERROR: --out cannot be used with --sheet.")

    rng = np.random.default_rng(args.seed)

    # ── Build job list ─────────────────────────────────────────────────────
    # Each job: { "bw_paths", "out_paths", "sample", "null_qcat_name" }
    jobs = []

    if args.bw:
        out_dir = Path(args.out_dir) if args.out_dir else None
        if out_dir:
            out_dir.mkdir(parents=True, exist_ok=True)
        bw_paths = [Path(p) for p in args.bw]
        out_paths = []
        for p in bw_paths:
            if args.out:
                out_paths.append(Path(args.out))
            else:
                out_paths.append(out_dir / (p.stem + args.suffix + ".bw"))
        null_name = _null_qcat_name(None, "sample", args.suffix)
        jobs.append({
            "sample": "",
            "bw_paths": bw_paths,
            "out_paths": out_paths,
            "out_dir": out_dir or Path("."),
            "null_qcat_name": null_name,
        })

    else:  # --sheet
        sheet = load_sheet(args.sheet)
        base_out_dir = Path(args.out_dir)
        for row in sheet:
            sample_dir = base_out_dir / (row["sample"] or "sample")
            sample_dir.mkdir(parents=True, exist_ok=True)
            out_paths = [
                sample_dir / (p.stem + args.suffix + ".bw")
                for p in row["bw_paths"]
            ]
            null_name = _null_qcat_name(row["out"], row["sample"], args.suffix)
            jobs.append({
                "sample": row["sample"],
                "bw_paths": row["bw_paths"],
                "out_paths": out_paths,
                "out_dir": sample_dir,
                "null_qcat_name": null_name,
            })

    # ── Process jobs ───────────────────────────────────────────────────────
    total_files = sum(len(j["bw_paths"]) for j in jobs)
    print(f"shift_bigwig.py: {len(jobs)} sample(s), "
          f"{total_files} BigWig file(s) total.")

    next_step_cmds = []

    # Build one global task list so workers are shared across all samples.
    shifted_paths_by_job = [[None] * len(job["bw_paths"]) for job in jobs]
    all_tasks = []

    for job_idx, job in enumerate(jobs):
        sample_label = job["sample"] or "sample"
        print(f"\n[{sample_label}]")
        for file_idx, (in_path, out_path) in enumerate(zip(job["bw_paths"], job["out_paths"])):
            if not in_path.exists():
                sys.exit(f"ERROR: file not found: {in_path}")

            file_seed = int(rng.integers(0, 2**31))
            shift_desc = (f"{args.shift:,} bp (fixed)" if args.shift
                          else f"random (seed={file_seed})")
            print(f"  {in_path.name}  ->  {out_path.name}  [{shift_desc}]")
            all_tasks.append((job_idx, file_idx, in_path, out_path, file_seed))

    def _run_shift(task):
        job_idx, file_idx, in_path, out_path, file_seed = task
        circular_shift_bigwig(
            in_path, out_path,
            shift=args.shift,
            bin_size=args.bin_size,
            min_shift=args.min_shift,
            seed=file_seed,
            per_chrom_shift=args.per_chrom_shift,
        )
        return job_idx, file_idx, out_path

    workers = max(1, int(args.workers))
    if workers > 1 and len(all_tasks) > 1:
        max_workers = min(workers, len(all_tasks))
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as ex:
            futures = [ex.submit(_run_shift, t) for t in all_tasks]
            for fut in concurrent.futures.as_completed(futures):
                job_idx, file_idx, out_path = fut.result()
                shifted_paths_by_job[job_idx][file_idx] = out_path
    else:
        for task in all_tasks:
            job_idx, file_idx, out_path = _run_shift(task)
            shifted_paths_by_job[job_idx][file_idx] = out_path

    for job_idx, job in enumerate(jobs):
        sample_label = job["sample"] or "sample"
        null_qcat_path = job["out_dir"] / job["null_qcat_name"]
        shifted_paths = shifted_paths_by_job[job_idx]
        bw_args = " ".join(str(p) for p in shifted_paths)
        cmd = (f"python bigwig_to_qcat.py \\\n"
               f"  --bw {bw_args} \\\n"
               f"  --out {null_qcat_path}")
        next_step_cmds.append((sample_label, cmd, null_qcat_path))

    # ── Print next steps ───────────────────────────────────────────────────
    print(f"\n{'='*60}")
    print("Next step: generate null qcat(s) from the shifted BigWigs:")
    print()
    for sample_label, cmd, _ in next_step_cmds:
        if sample_label:
            print(f"# {sample_label}")
        print(cmd)
        print()

    if len(next_step_cmds) > 1:
        print("Then run bearing_pvalue.py for each sample, e.g.:")
        sample_label, _, null_path = next_step_cmds[0]
        orig_name = sample_label or "sample"
        print(f"python bearing_pvalue.py \\")
        print(f"  --qcat {orig_name}.qcat.bgz \\")
        print(f"  --null-qcat {null_path} \\")
        print(f"  --out-prefix results/{orig_name} \\")
        print(f"  --score-plot")
    else:
        _, _, null_path = next_step_cmds[0]
        print("Then run bearing_pvalue.py:")
        print(f"python bearing_pvalue.py \\")
        print(f"  --qcat sample.qcat.bgz \\")
        print(f"  --null-qcat {null_path} \\")
        print(f"  --out-prefix results/sample \\")
        print(f"  --score-plot")


if __name__ == "__main__":
    main()