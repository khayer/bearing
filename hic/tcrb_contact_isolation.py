#!/usr/bin/env python3
"""
tcrb_contact_isolation_v3.py

v3 changes from v2:
  - Memory-efficient chromosome distance-decay computation. v2 fetched
    the full chromosome as a dense matrix (e.g. ~1.8 GB at 10 kb on
    chr6), causing OOM kills. v3 walks the chromosome in chunks of
    --chunk-size-bins (default 2000) with overlap = max_k, processing
    each chunk's diagonals incrementally. Peak memory now ~50-100 MB
    per chunk.
  - New --perm-n flag: for each (condition, distance) cell, draws
    distance-matched random control pairs of equivalent-sized intervals
    from the same chromosome (excluding the Tcrb target itself). The
    empirical p-value compares observed Tcrb -> flank contact to this
    null. Significant points (p_emp < 0.05) are marked with stars on
    the plot. Default 200 permutations per cell; 0 disables.
  - TSV adds columns: null_median, null_p95, p_emp.

v2 changes from v1:
  - Default normalization is observed/expected (OE).

Tests how isolated the narrow Tcrb locus (chr6:40.85-41.6M) is from
flanking chromosomal regions, and whether perturbations alter this
isolation pattern.

Method:
  1. Define a target region (default: chr6:40850000-41600000, narrow Tcrb).
  2. Define a series of flanking windows of equal size at increasing
     distances upstream and downstream. For a 750 kb target with 500 kb
     step, windows at -3.5, -3.0, -2.5, ..., +0.5, +1.0, ..., +9.0 Mb.
  3. For each .cool file (one per condition), extract the balanced
     submatrix containing the target + all flanking windows.
  4. For each (target, flanking_window) pair, compute the mean balanced
     contact frequency across all bin-pairs in the rectangle.
  5. Plot mean contact frequency vs. signed distance from target center,
     one line per condition.

Outputs:
  <prefix>.contact_by_distance.tsv    long-format per-condition data
  <prefix>.contact_isolation.pdf      one-line-per-condition plot

Usage:
  python tcrb_contact_isolation.py \\
      --cool DN=../hic_files/merged_corrected_KR_DN_bs_25000.cool \\
             DP=../hic_files/merged_corrected_KR_DP_bs_25000.cool \\
             ProB=../hic_files/merged_corrected_KR_ProB_bs_25000.cool \\
             EbKO=../hic_files/merged_corrected_KR_EBKO_bs_25000.cool \\
             S3T3=../hic_files/merged_corrected_KR_s3T3_bs_25000.cool \\
             dV1P=../hic_files/merged_corrected_KR_dV1P_bs_25000.cool \\
             dV1CTCF=../hic_files/merged_corrected_KR_dV1CTCF_bs_25000.cool \\
      --target chr6:40850000-41600000 \\
      --window-size 750000 \\
      --max-distance 6000000 \\
      --out-prefix tcrb_contact_isolation
"""
import argparse
import os
import sys
import numpy as np
import pandas as pd


def fetch_cool_matrix(path, chrom, start, end, balance=True):
    """Fetch balanced contact matrix as a numpy array."""
    import cooler
    c = cooler.Cooler(path)
    region = "{}:{}-{}".format(chrom, start, end)
    mat = c.matrix(balance=balance, sparse=False).fetch(region)
    bins = c.bins().fetch(region)
    binsize = int(c.binsize)
    return mat, bins["start"].to_numpy(), binsize


def compute_chrom_distance_decay(path, chrom, balance=True,
                                    max_distance_bins=600,
                                    chunk_size_bins=2000):
    """Compute mean balanced contact frequency per genomic distance for
    one chromosome.

    Memory-efficient: chunks the chromosome into overlapping slices of
    chunk_size_bins, fetching each as a dense matrix and accumulating
    diagonal sums per distance. Overlap = max_distance_bins so we don't
    miss bin-pairs that span chunk boundaries.

    At 10 kb resolution, chr6 has ~15000 bins. A full dense matrix is
    ~1.8 GB. With chunk_size_bins=2000 and overlap of max_distance_bins
    (default 600), each fetch is at most 2600 x 2600 float64 = ~54 MB,
    and we make ~8 fetches per chromosome -- much more memory-friendly.
    """
    import cooler
    c = cooler.Cooler(path)
    bins = c.bins().fetch(chrom)
    n_bins = len(bins)
    if n_bins == 0:
        return np.full(max_distance_bins + 1, np.nan)
    binsize = int(c.binsize)
    max_k = min(max_distance_bins, n_bins - 1)
    sums = np.zeros(max_k + 1, dtype=np.float64)
    counts = np.zeros(max_k + 1, dtype=np.int64)

    # Walk along the chromosome in non-overlapping primary windows of
    # chunk_size_bins, but fetch a slightly larger window (with right
    # padding = max_k) so all bin-pairs within distance max_k from the
    # primary window are captured. Each pair is counted exactly once
    # by restricting bin1_idx to the primary window range.
    step = chunk_size_bins
    pad = max_k
    chr_start_bp = int(bins["start"].iloc[0])
    chr_end_bp = int(bins["end"].iloc[-1])

    n_chunks = (n_bins + step - 1) // step
    for chunk_i in range(n_chunks):
        i_lo = chunk_i * step
        i_hi = min(n_bins, i_lo + step)
        # Fetch primary [i_lo, i_hi) plus right padding [i_hi, i_hi+pad)
        fetch_i_hi = min(n_bins, i_hi + pad)
        fetch_start_bp = int(bins["start"].iloc[i_lo])
        fetch_end_bp = int(bins["end"].iloc[fetch_i_hi - 1])
        region = "{}:{}-{}".format(chrom, fetch_start_bp, fetch_end_bp)
        mat = c.matrix(balance=balance, sparse=False).fetch(region)
        # primary rows are [0, i_hi - i_lo)
        primary_n = i_hi - i_lo
        # Accumulate over diagonals
        for k in range(max_k + 1):
            if k >= mat.shape[1]:
                break
            # Cells (i, i+k) for i in [0, primary_n) where i+k < mat.shape[1]
            i_max = min(primary_n, mat.shape[0], mat.shape[1] - k)
            if i_max <= 0:
                continue
            ii = np.arange(i_max)
            vals = mat[ii, ii + k]
            finite = np.isfinite(vals)
            if finite.any():
                sums[k] += float(np.nansum(vals[finite]))
                counts[k] += int(finite.sum())
        del mat

    expected = np.full(max_k + 1, np.nan)
    valid = counts > 0
    expected[valid] = sums[valid] / counts[valid]
    return expected


def random_control_contacts(matrix, bin_starts, target_start, target_end,
                              distance_bp, window_size_bp,
                              n_perms=200, rng=None,
                              expected_by_dist=None, binsize=None,
                              chrom_min_bp=None, chrom_max_bp=None):
    """For a given target-to-flank distance, draw n_perms random pairs of
    (target_size, window_size) intervals separated by the same distance,
    anywhere on the matrix, and return the distribution of their mean
    contact (OE-normalized if expected_by_dist given).

    Used to ask: at this distance, is the observed Tcrb -> flanking
    contact unusually high/low relative to random pairs at the same
    distance on the same chromosome?
    """
    if rng is None:
        rng = np.random.default_rng(42)
    t_width = target_end - target_start
    bs = bin_starts
    bs_min = int(bs.min())
    bs_max = int(bs.max())
    if chrom_min_bp is None:
        chrom_min_bp = bs_min
    if chrom_max_bp is None:
        chrom_max_bp = bs_max + binsize if binsize else bs_max

    perm_vals = []
    n_tries = 0
    max_tries = n_perms * 20
    while len(perm_vals) < n_perms and n_tries < max_tries:
        n_tries += 1
        sign = 1 if rng.random() < 0.5 else -1
        if sign > 0:
            t_lo = chrom_min_bp
            t_hi = chrom_max_bp - t_width - distance_bp - window_size_bp // 2
        else:
            t_lo = chrom_min_bp + distance_bp + window_size_bp // 2
            t_hi = chrom_max_bp - t_width
        if t_hi <= t_lo:
            continue
        rt_start = int(rng.integers(t_lo, t_hi))
        rt_end = rt_start + t_width
        rt_mid = (rt_start + rt_end) // 2
        rw_center = rt_mid + sign * distance_bp
        rw_start = rw_center - window_size_bp // 2
        rw_end = rw_center + window_size_bp // 2
        if rw_start < chrom_min_bp or rw_end > chrom_max_bp:
            continue
        if not (rt_end < target_start or rt_start > target_end):
            continue
        v, n = mean_contacts_between(matrix, bin_starts,
                                       rt_start, rt_end,
                                       rw_start, rw_end,
                                       expected_by_dist=expected_by_dist,
                                       binsize=binsize)
        if np.isfinite(v):
            perm_vals.append(v)
    return np.array(perm_vals)


def mean_contacts_between(matrix, bin_starts, target_start, target_end,
                            window_start, window_end,
                            expected_by_dist=None, binsize=None):
    """Mean balanced contact frequency between two genomic intervals.

    If expected_by_dist is provided, returns observed/expected (OE) ratio:
    each (i,j) cell is divided by expected[|i-j|] before averaging.
    Cells with non-finite values are skipped.
    """
    bs = bin_starts
    target_mask = (bs >= target_start) & (bs < target_end)
    window_mask = (bs >= window_start) & (bs < window_end)
    if not target_mask.any() or not window_mask.any():
        return float("nan"), 0
    sub = matrix[np.ix_(target_mask, window_mask)]
    if expected_by_dist is None:
        finite = np.isfinite(sub)
        if not finite.any():
            return float("nan"), 0
        return float(np.nanmean(sub[finite])), int(finite.sum())

    # OE: build distance matrix in bin units
    t_idx = np.where(target_mask)[0]
    w_idx = np.where(window_mask)[0]
    # |i - j| for each cell in the submatrix
    dist_mat = np.abs(t_idx[:, None] - w_idx[None, :])
    max_d = len(expected_by_dist) - 1
    dist_mat = np.minimum(dist_mat, max_d)
    exp_sub = expected_by_dist[dist_mat]
    # avoid divide by zero / nan in expected
    with np.errstate(divide="ignore", invalid="ignore"):
        oe = sub / exp_sub
    finite = np.isfinite(oe)
    if not finite.any():
        return float("nan"), 0
    return float(np.nanmean(oe[finite])), int(finite.sum())


def main():
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cool", required=True, nargs="+",
                    help="CONDITION=PATH entries for each .cool file")
    ap.add_argument("--target", default="chr6:40850000-41600000",
                    help="target region in chrom:start-end format "
                         "(default chr6:40850000-41600000 = narrow Tcrb)")
    ap.add_argument("--window-size", type=int, default=750000,
                    help="size of flanking windows in bp (default 750000, "
                         "approx. matches default target size)")
    ap.add_argument("--max-distance", type=int, default=6000000,
                    help="how far up- and downstream to probe in bp "
                         "(default 6000000 = 6 Mb)")
    ap.add_argument("--no-balance", action="store_true",
                    help="skip cooler balance (default: balance applied)")
    ap.add_argument("--log-y", action="store_true",
                    help="use log scale on y-axis (often useful for "
                         "contact frequency)")
    ap.add_argument("--no-oe", action="store_true",
                    help="skip observed/expected normalization. By "
                         "default this script divides each contact by "
                         "the chromosome-wide mean at the same genomic "
                         "distance, giving fold-enrichment over distance "
                         "decay (OE). Use --no-oe to plot raw balanced "
                         "contact frequency.")
    ap.add_argument("--perm-n", type=int, default=200,
                    help="number of random distance-matched control pairs "
                         "per (condition, distance) cell, for the empirical "
                         "p-value. Default 200; set 0 to disable.")
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    # Parse target
    chrom, rng = args.target.split(":")
    t_start, t_end = map(int, rng.split("-"))
    t_mid = (t_start + t_end) // 2
    t_width = t_end - t_start
    print("Target: {}:{:,}-{:,} ({:.0f} kb wide, center {:,})"
          .format(chrom, t_start, t_end, t_width / 1000, t_mid), flush=True)

    # Parse cool paths
    cond_paths = []
    for spec in args.cool:
        if "=" not in spec:
            sys.exit("invalid --cool spec: " + spec)
        name, path = spec.split("=", 1)
        cond_paths.append((name.strip(), path.strip()))

    # Build flanking windows. Use windows that DO NOT overlap target.
    # Window centers stepped by --window-size starting one window-size
    # away from the target edges, on both sides.
    upstream_centers = []
    downstream_centers = []
    step = args.window_size
    # Upstream: starting one step from target start, going to -max_distance
    pos = t_start - step // 2
    while (t_mid - pos) <= args.max_distance and pos - step // 2 >= 0:
        upstream_centers.append(pos)
        pos -= step
    # Downstream: starting one step past target end
    pos = t_end + step // 2
    while (pos - t_mid) <= args.max_distance:
        downstream_centers.append(pos)
        pos += step

    # Convert to (start, end) tuples
    flanking = []
    for c in upstream_centers:
        flanking.append((c - step // 2, c + step // 2,
                          c - t_mid))  # signed distance from target center
    for c in downstream_centers:
        flanking.append((c - step // 2, c + step // 2, c - t_mid))
    flanking.sort(key=lambda x: x[2])
    print("Probing {} flanking windows ({} upstream + {} downstream)"
          .format(len(flanking), len(upstream_centers),
                  len(downstream_centers)))

    # Determine extraction region: covers target + furthest flanking window
    min_pos = min([t_start] + [f[0] for f in flanking])
    max_pos = max([t_end] + [f[1] for f in flanking])
    # Pad to bin boundary; cooler handles this
    extract_start = max(0, min_pos - 100000)
    extract_end = max_pos + 100000

    # Run per condition
    rows = []
    for cond, path in cond_paths:
        if not os.path.isfile(path):
            print("WARNING: {} not found, skipping".format(path))
            continue
        print("\nLoading", cond, "from", path)
        try:
            mat, bin_starts, binsize = fetch_cool_matrix(
                path, chrom, extract_start, extract_end,
                balance=(not args.no_balance))
        except Exception as e:
            print("  ERROR loading {}: {}".format(cond, e))
            continue
        print("  matrix shape {} (binsize {} bp, region {:,}-{:,})".format(
            mat.shape, binsize, extract_start, extract_end))

        # Per-condition expected-by-distance (chrom-wide), for OE
        expected_by_dist = None
        if not args.no_oe:
            # Max distance we'll ever query: ~max_distance + window/2
            max_d_bp = args.max_distance + args.window_size
            max_d_bins = int(np.ceil(max_d_bp / binsize)) + 10
            print("  computing chrom-wide distance decay (up to {} bins) ...".format(
                max_d_bins), flush=True)
            expected_by_dist = compute_chrom_distance_decay(
                path, chrom, balance=(not args.no_balance),
                max_distance_bins=max_d_bins)
            # Quick sanity print
            valid_exp = expected_by_dist[np.isfinite(expected_by_dist)]
            if len(valid_exp) > 0:
                print("    expected[0]={:.3e}  expected[10]={:.3e}  "
                      "expected[100]={:.3e}".format(
                          expected_by_dist[0] if np.isfinite(expected_by_dist[0]) else float("nan"),
                          expected_by_dist[min(10, len(expected_by_dist)-1)] if np.isfinite(expected_by_dist[min(10, len(expected_by_dist)-1)]) else float("nan"),
                          expected_by_dist[min(100, len(expected_by_dist)-1)] if np.isfinite(expected_by_dist[min(100, len(expected_by_dist)-1)]) else float("nan")))

        for f_start, f_end, dist in flanking:
            mean_c, n_cells = mean_contacts_between(
                mat, bin_starts, t_start, t_end, f_start, f_end,
                expected_by_dist=expected_by_dist, binsize=binsize)

            # Distance-matched random control on the same matrix
            p_emp = float("nan")
            null_median = float("nan")
            null_p95 = float("nan")
            if args.perm_n > 0 and np.isfinite(mean_c):
                rng = np.random.default_rng(seed=abs(hash(cond + str(dist))) % (2**32))
                nulls = random_control_contacts(
                    mat, bin_starts, t_start, t_end,
                    abs(dist), args.window_size,
                    n_perms=args.perm_n, rng=rng,
                    expected_by_dist=expected_by_dist, binsize=binsize,
                    chrom_min_bp=extract_start, chrom_max_bp=extract_end)
                if len(nulls) >= 20:
                    # Two-sided empirical p
                    n_extreme = int(np.sum(nulls >= mean_c))
                    n_total = len(nulls)
                    # one-sided "observed >= null" then doubled
                    p_one = (n_extreme + 1) / (n_total + 1)
                    p_emp = min(1.0, 2 * min(p_one, 1 - p_one + 1.0/(n_total+1)))
                    null_median = float(np.median(nulls))
                    null_p95 = float(np.percentile(nulls, 95))

            rows.append({
                "condition": cond,
                "chrom": chrom,
                "target_start": t_start,
                "target_end": t_end,
                "flank_start": f_start,
                "flank_end": f_end,
                "distance_to_target_center": dist,
                "abs_distance": abs(dist),
                "side": "downstream" if dist > 0 else "upstream",
                "mean_contact": mean_c,
                "null_median": null_median,
                "null_p95": null_p95,
                "p_emp": p_emp,
                "n_cells": n_cells,
            })

    df = pd.DataFrame(rows)
    if len(df) == 0:
        sys.exit("no data collected")

    out_tsv = args.out_prefix + ".contact_by_distance.tsv"
    df.to_csv(out_tsv, sep="\t", index=False, float_format="%.6e")
    print("\nWrote", out_tsv)

    # Plot
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter

    fig, (ax1, ax2) = plt.subplots(
        2, 1, sharex=True, figsize=(11, 8),
        gridspec_kw={"height_ratios": [3, 1]})

    # 7-color palette extended; adjust if needed
    palette = {
        "DN": "#e63946",
        "DP": "#1d3557",
        "ProB": "#2a9d8f",
        "EbKO": "#f4a261",
        "S3T3": "#9d4edd",
        "dV1P": "#e76f51",
        "dV1CTCF": "#264653",
    }
    fallback = ["#0072B2", "#D55E00", "#009E73", "#F0E442",
                  "#56B4E9", "#CC79A7", "#E69F00"]

    for i, cond in enumerate(df["condition"].unique()):
        sub = df[df["condition"] == cond].sort_values("distance_to_target_center")
        color = palette.get(cond, fallback[i % len(fallback)])
        ax1.plot(sub["distance_to_target_center"] / 1e6,
                  sub["mean_contact"],
                  marker="o", markersize=3, lw=1.2,
                  color=color, label=cond, alpha=0.8)
        # Star significant points (p_emp < 0.05)
        if "p_emp" in sub.columns:
            sig = sub[sub["p_emp"] < 0.05]
            if len(sig) > 0:
                ax1.scatter(sig["distance_to_target_center"] / 1e6,
                              sig["mean_contact"],
                              marker="*", s=70, color=color,
                              edgecolor="black", linewidth=0.6,
                              zorder=10)

    ax1.axvline(0, color="gray", lw=0.5, ls="--", alpha=0.5)
    ax1.axvspan(-t_width/2/1e6, t_width/2/1e6,
                  color="gold", alpha=0.15, zorder=0,
                  label="target ({:.0f} kb)".format(t_width/1000))
    if args.no_oe:
        y_label = "Mean balanced contact frequency"
    else:
        y_label = "Mean observed/expected contact"
        # Add reference line at OE=1 (no enrichment)
        ax1.axhline(1.0, color="black", lw=0.4, ls=":", alpha=0.5)
    ax1.set_ylabel(y_label)
    if args.log_y:
        ax1.set_yscale("log")
    ax1.legend(fontsize=8, loc="best", ncol=2)
    title_norm = "observed/expected" if not args.no_oe else "balanced contacts"
    ax1.set_title(
        "Contact ({}) from {}:{:,}-{:,} to flanking windows "
        "({:.0f} kb each)".format(title_norm, chrom, t_start, t_end,
                                    args.window_size/1000),
        fontsize=10)
    for sp in ("top", "right"):
        ax1.spines[sp].set_visible(False)

    # Fold-change panel: each condition normalized to DN at the same distance
    # (or first condition if DN missing)
    ref_cond = "DN" if "DN" in df["condition"].unique() else df["condition"].iloc[0]
    ref = df[df["condition"] == ref_cond].set_index(
        "distance_to_target_center")["mean_contact"]
    for i, cond in enumerate(df["condition"].unique()):
        if cond == ref_cond:
            continue
        sub = df[df["condition"] == cond].sort_values("distance_to_target_center")
        sub = sub.copy()
        sub["ref"] = sub["distance_to_target_center"].map(ref)
        with np.errstate(divide="ignore", invalid="ignore"):
            sub["log2_vs_ref"] = np.log2(sub["mean_contact"] / sub["ref"])
        color = palette.get(cond, fallback[i % len(fallback)])
        ax2.plot(sub["distance_to_target_center"] / 1e6,
                  sub["log2_vs_ref"],
                  marker="o", markersize=3, lw=1.0,
                  color=color, label=cond)
    ax2.axhline(0, color="gray", lw=0.5, ls="--", alpha=0.5)
    ax2.axvline(0, color="gray", lw=0.5, ls="--", alpha=0.5)
    ax2.axvspan(-t_width/2/1e6, t_width/2/1e6,
                  color="gold", alpha=0.15, zorder=0)
    ax2.set_ylabel("log2(condition / {})".format(ref_cond), fontsize=9)
    ax2.set_xlabel("Signed distance from {} center (Mb)".format(args.target))
    for sp in ("top", "right"):
        ax2.spines[sp].set_visible(False)

    plt.tight_layout()
    out_pdf = args.out_prefix + ".contact_isolation.pdf"
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()
    print("Wrote", out_pdf)

    # Brief stdout summary: contacts at +0.5, +1, +2, +5 Mb compared
    # across conditions
    print("\n=== Contact frequency at selected downstream distances ===")
    targets = [args.window_size, 2 * args.window_size,
                3 * args.window_size, 5 * args.window_size]
    cond_list = list(df["condition"].unique())
    header = "  " + " " * 10 + "".join("{:>10}".format(c) for c in cond_list)
    print(header)
    for tgt in targets:
        # find closest distance bin
        avail = df["distance_to_target_center"].unique()
        nearest = avail[np.argmin(np.abs(avail - tgt))]
        line = "  +{:>4.1f} Mb  ".format(nearest / 1e6)
        for c in cond_list:
            r = df[(df["condition"] == c)
                    & (df["distance_to_target_center"] == nearest)]
            if len(r):
                v = float(r["mean_contact"].iloc[0])
                line += "{:>10.2e}".format(v)
            else:
                line += "{:>10}".format("NA")
        print(line)


if __name__ == "__main__":
    main()
