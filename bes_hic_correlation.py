#!/usr/bin/env python3
"""
bes_hic_correlation_v5.py

Test BEARING manuscript claim 1: do BEARING 1D compositional shifts (BES)
co-localize with Hi-C 3D contact-pattern changes between the same two
conditions?

v5 changes from v4:
  - Per-track stratification uses |kl_i| not kl_i. The signed kl_i values
    can be negative (e.g., kl_2 around -2 at Eb-Trbv31 means RNA+ is gained
    in B vs A). Using signed p95 misses bins dominated by "negative-side"
    tracks entirely. Now per-track p95 is over abs(kl_i), so dominance
    reflects magnitude of contribution regardless of direction.

v4 changes from v3:
  - Header-aware BES reader: handles BEARING per-bin output table with
    named columns (chrom, start, end, bearing_score, ..., kl_1..kl_6)
  - Column-by-name selection (--bes-score-col, default bearing_score)
  - Type coercion to silence DtypeWarning

v3 changes from v2:
  - balance=True by default (KR weights confirmed in cooler 'weight' column)

v2 changes from v1:
  - .cool input via cooler (replaces .h5 / hicmatrix)
  - Percentile aggregation (p95 primary, p75, median); no mean
  - Multi-column .bm support: row-wise mean across all score columns

INPUT (one pairwise comparison at a time, e.g. DN vs DP):
  --bes               per-bin BES bedgraph (BEARING output)
  --insul-A, --insul-B  HiCExplorer .bm files for the two conditions
  --cool-A,  --cool-B   .cool contact matrices (optional but recommended)
  --region            chr:start-end (default: whole BES file)
  --hic-bin           target bin size for aggregation (default 10000)
  --balance           apply cooler weight balancing (default: False;
                      assumes KR is already baked into matrix values)

USAGE
  python bes_hic_correlation_v2.py \\
    --bes        results_v6/DN_vs_DP_bes.bedgraph \\
    --insul-A    .../07tad/merged_corrected_KR_DN_bs_10000_tad_score.bm \\
    --insul-B    .../07tad/merged_corrected_KR_DP_bs_10000_tad_score.bm \\
    --cool-A     .../cool/merged_corrected_KR_DN.10000.cool \\
    --cool-B     .../cool/merged_corrected_KR_DP.10000.cool \\
    --region     chr6:40000000-47400000 \\
    --hic-bin    10000 \\
    --prefix     DN_vs_DP_corr_10kb

OUTPUTS
  PREFIX.tsv            per-bin merged table
  PREFIX.stats.json     correlations, null p-values, top-decile enrichment
  PREFIX.scatter.pdf    tracks + scatter + null distributions
"""

import argparse
import json
import os
import sys
import warnings

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import cooler
    HAVE_COOLER = True
except ImportError:
    HAVE_COOLER = False


# ---------------------------------------------------------------------
# IO
# ---------------------------------------------------------------------

def parse_region(s):
    chrom, rest = s.split(":")
    a, b = rest.replace(",", "").split("-")
    return chrom, int(a), int(b)


def read_bes_bedgraph(path, region=None, score_col="bearing_score"):
    """Read BEARING per-bin output table.

    Auto-detects whether the file has a header line. If header is present,
    columns are addressed by name; otherwise falls back to positional
    (chrom, start, end, bes, kl_1...).

    Required columns (by name): chrom, start, end, and either
    'bearing_score' (default) or 'bes'. Per-track kl_* columns are
    auto-detected by name and used for stratification.

    BES is non-negative from KL; abs() is a safety no-op."""
    # Sniff first non-blank line to detect a header
    with open(path) as fh:
        first = fh.readline().strip()
        while first.startswith("#") or first == "":
            first = fh.readline().strip()
    fields = first.split("\t")
    has_header = fields[0].lower() in ("chrom", "chr", "#chrom", "seqnames")

    if has_header:
        df = pd.read_csv(path, sep="\t", header=0, comment="#",
                          low_memory=False)
        df.columns = [c.strip() for c in df.columns]
    else:
        df = pd.read_csv(path, sep="\t", header=None, comment="#",
                          low_memory=False)
        n = df.shape[1]
        df.columns = ["chrom", "start", "end", "bearing_score"] + \
                     ["kl_{}".format(i + 1) for i in range(n - 4)]

    # Rename score column to canonical "bes"
    if score_col in df.columns:
        df = df.rename(columns={score_col: "bes"})
    elif "bes" in df.columns:
        pass
    else:
        sys.exit("Could not find '{}' or 'bes' column in {}. "
                 "Got columns: {}".format(score_col, path,
                                          list(df.columns)))

    # Coerce types (handles mixed-types warning)
    for c in ("start", "end", "bes"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["bes", "start", "end"]).copy()
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)

    # Identify kl_* track columns; coerce to numeric
    track_cols = sorted(
        [c for c in df.columns if c.startswith("kl_")],
        key=lambda c: int(c.split("_")[1]) if c.split("_")[1].isdigit()
                        else 999,
    )
    for tc in track_cols:
        df[tc] = pd.to_numeric(df[tc], errors="coerce").fillna(0.0)

    df["abs_bes"] = df["bes"].abs()

    if region is not None:
        chrom, a, b = region
        df = df[(df["chrom"] == chrom) & (df["end"] > a)
                & (df["start"] < b)]

    return df.reset_index(drop=True), track_cols


def read_bm(path, region=None):
    """Read HiCExplorer .bm. Auto-detects number of score columns and
    returns the row-wise MEAN across all score columns."""
    df = pd.read_csv(path, sep="\t", header=None, comment="#")
    if df.shape[1] < 4:
        sys.exit("{}: .bm needs >=4 columns".format(path))
    out = df.iloc[:, :3].copy()
    out.columns = ["chrom", "start", "end"]
    score_cols = df.iloc[:, 3:]
    out["score"] = score_cols.mean(axis=1).values
    out["n_score_cols"] = score_cols.shape[1]
    if region is not None:
        chrom, a, b = region
        out = out[(out["chrom"] == chrom) & (out["end"] > a) & (out["start"] < b)]
    return out.reset_index(drop=True), score_cols.shape[1]


# ---------------------------------------------------------------------
# Hi-C bin grid + BES aggregation
# ---------------------------------------------------------------------

def make_hic_bins(chrom, start, end, binsize):
    bin_start = (start // binsize) * binsize
    bin_end = ((end + binsize - 1) // binsize) * binsize
    starts = np.arange(bin_start, bin_end, binsize)
    return pd.DataFrame({
        "chrom": chrom,
        "start": starts,
        "end": starts + binsize,
        "bin_idx": np.arange(len(starts)),
    })


def aggregate_to_hic_bins(bes_df, hic_bins, binsize, track_cols):
    """Assign each BEARING bin to a Hi-C bin by midpoint; aggregate using
    median, p75, p95 (no mean; max as bonus)."""
    bes = bes_df.copy()
    bes["mid"] = (bes["start"] + bes["end"]) // 2
    base_start = int(hic_bins.iloc[0]["start"])
    bes["bin_idx"] = ((bes["mid"] - base_start) // binsize).astype(int)
    bes = bes[(bes["bin_idx"] >= 0) & (bes["bin_idx"] < len(hic_bins))]

    def p75(x):
        return np.percentile(x, 75) if len(x) else np.nan

    def p95(x):
        return np.percentile(x, 95) if len(x) else np.nan

    grouped = bes.groupby("bin_idx")["abs_bes"]
    agg = pd.DataFrame({
        "med_bes": grouped.median(),
        "p75_bes": grouped.apply(p75),
        "p95_bes": grouped.apply(p95),
        "max_bes": grouped.max(),
        "n_bins": grouped.count(),
    }).reset_index()

    # per-track p95 of |contribution| (kl_i can be signed by condition;
    # we want which track has the largest magnitude contribution)
    for tc in track_cols:
        abs_col = "abs_{}".format(tc)
        bes[abs_col] = bes[tc].abs()
        g = bes.groupby("bin_idx")[abs_col].apply(p95)
        agg["{}_abs_p95".format(tc)] = agg["bin_idx"].map(g)

    out = hic_bins.merge(agg, on="bin_idx", how="left")
    for c in ["med_bes", "p75_bes", "p95_bes", "max_bes", "n_bins"]:
        out[c] = out[c].fillna(0.0)
    return out


def aggregate_bm_to_hic_bins(bm_df, hic_bins, binsize, score_name):
    """Insulation already at ~Hi-C resolution; average if multiple .bm bins
    map to one Hi-C bin."""
    bm = bm_df.copy()
    bm["mid"] = (bm["start"] + bm["end"]) // 2
    base_start = int(hic_bins.iloc[0]["start"])
    bm["bin_idx"] = ((bm["mid"] - base_start) // binsize).astype(int)
    bm = bm[(bm["bin_idx"] >= 0) & (bm["bin_idx"] < len(hic_bins))]
    g = bm.groupby("bin_idx")["score"].mean().rename(score_name)
    return hic_bins.merge(g, on="bin_idx", how="left")


# ---------------------------------------------------------------------
# Hi-C contact-difference per bin (cooler)
# ---------------------------------------------------------------------

def fetch_cool_matrix(path, chrom, start, end, balance):
    if not HAVE_COOLER:
        sys.exit("cooler not installed: pip install cooler")
    c = cooler.Cooler(path)
    binsize = c.binsize
    region = "{}:{}-{}".format(chrom, start, end)
    selector = c.matrix(balance=balance, sparse=True)
    M = selector.fetch(region).tocsr()
    bins = c.bins().fetch(region).reset_index(drop=True)
    return M, bins, binsize


def per_bin_delta_contact(M_A, M_B, binsize, min_distance, max_distance):
    """Per bin i, sum over j of |M_A[i,j] - M_B[i,j]| for j such that
    min_distance <= (j - i) * binsize <= max_distance.
    Sums to both i and j so each bin accumulates symmetrically."""
    M_A = M_A.astype(float)
    M_B = M_B.astype(float)
    diff = abs(M_A - M_B)

    diff_coo = diff.tocoo()
    rows = diff_coo.row
    cols = diff_coo.col
    vals = diff_coo.data

    # drop NaN from balancing
    finite = np.isfinite(vals)
    rows = rows[finite]
    cols = cols[finite]
    vals = vals[finite]

    min_bins = max(1, int(min_distance // binsize))
    max_bins = int(max_distance // binsize)

    dist = cols - rows
    keep = (dist >= min_bins) & (dist <= max_bins)
    rows = rows[keep]
    cols = cols[keep]
    vals = vals[keep]

    n = diff.shape[0]
    delta = (np.bincount(rows, weights=vals, minlength=n)
             + np.bincount(cols, weights=vals, minlength=n))
    return delta


# ---------------------------------------------------------------------
# Correlations + null
# ---------------------------------------------------------------------

def corr_pair(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 20:
        return {"n": int(m.sum()), "spearman_rho": None,
                "spearman_p": None, "pearson_r": None, "pearson_p": None}
    rho_s, p_s = stats.spearmanr(x[m], y[m])
    rho_p, p_p = stats.pearsonr(x[m], y[m])
    return {"n": int(m.sum()),
            "spearman_rho": float(rho_s), "spearman_p": float(p_s),
            "pearson_r": float(rho_p), "pearson_p": float(p_p)}


def circular_shift_null(x, y, n_perm, min_shift_bins, seed=42):
    rng = np.random.default_rng(seed)
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    xx = x[m]
    yy = y[m]
    n = len(xx)
    if n < 20:
        return None, None, None
    obs, _ = stats.spearmanr(xx, yy)
    null = np.empty(n_perm)
    lo = max(min_shift_bins, 1)
    hi = max(n - min_shift_bins, lo + 1)
    for k in range(n_perm):
        s = rng.integers(lo, hi)
        x_shift = np.roll(xx, s)
        r, _ = stats.spearmanr(x_shift, yy)
        null[k] = r
    p_emp = float((np.abs(null) >= abs(obs)).mean())
    p_emp = max(p_emp, 1.0 / n_perm)
    return float(obs), null, p_emp


def permutation_null(x, y, n_perm, seed=42):
    """Permutation null for a NON-CONTIGUOUS subset (e.g. a per-track stratum),
    where the circular shift would be degenerate (too few valid shifts) and
    there is no contiguity to preserve. Shuffles y, recomputes Spearman.
    Slightly less conservative than the shift null if residual autocorrelation
    remains, but properly sampled at small n."""
    rng = np.random.default_rng(seed)
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    xx, yy = x[m], y[m]
    n = len(xx)
    if n < 20:
        return None, None, None
    obs, _ = stats.spearmanr(xx, yy)
    null = np.empty(n_perm)
    for k in range(n_perm):
        null[k], _ = stats.spearmanr(xx, rng.permutation(yy))
    p_emp = max(float((np.abs(null) >= abs(obs)).mean()), 1.0 / n_perm)
    return float(obs), null, p_emp


def top_decile_enrichment(x, y, q=0.9):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    xx = x[m]
    yy = y[m]
    if len(xx) < 20:
        return None
    qx = np.quantile(xx, q)
    qy = np.quantile(yy, q)
    a = int(((xx >= qx) & (yy >= qy)).sum())
    b = int(((xx >= qx) & (yy < qy)).sum())
    c = int(((xx < qx) & (yy >= qy)).sum())
    d = int(((xx < qx) & (yy < qy)).sum())
    try:
        odds, p = stats.fisher_exact([[a, b], [c, d]])
    except Exception:
        odds, p = float("nan"), float("nan")
    return {"quantile": float(q), "n": int(len(xx)),
            "n_both_top": a, "n_x_only": b,
            "n_y_only": c, "n_neither": d,
            "odds_ratio": float(odds), "fisher_p": float(p)}


# ---------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------

def make_plots(df, region, results, nulls, out_pdf):
    fig, axes = plt.subplots(5, 2, figsize=(11, 14))

    chrom, rstart, rend = region
    sub = df[(df["chrom"] == chrom) & (df["end"] > rstart)
             & (df["start"] < rend)].copy()
    mids = (sub["start"] + sub["end"]) / 2 / 1e6

    # row 0: BES tracks (p95 and max)
    ax = axes[0, 0]
    ax.plot(mids, sub["p95_bes"], color="black", lw=0.9)
    ax.fill_between(mids, 0, sub["p95_bes"], alpha=0.2, color="black")
    ax.set_title("p95 |BES| per Hi-C bin")
    ax.set_xlim(rstart / 1e6, rend / 1e6)

    ax = axes[0, 1]
    ax.plot(mids, sub["p75_bes"], color="gray", lw=0.9)
    ax.fill_between(mids, 0, sub["p75_bes"], alpha=0.2, color="gray")
    ax.set_title("p75 |BES| per Hi-C bin")
    ax.set_xlim(rstart / 1e6, rend / 1e6)

    # row 1: Hi-C delta tracks
    ax = axes[1, 0]
    if "delta_insulation" in sub.columns:
        ax.plot(mids, sub["delta_insulation"], color="C0", lw=0.9)
        ax.fill_between(mids, 0, sub["delta_insulation"], alpha=0.2)
        ax.set_title("|delta insulation|")
        ax.set_xlim(rstart / 1e6, rend / 1e6)
    ax = axes[1, 1]
    if "delta_contact" in sub.columns and sub["delta_contact"].notna().any():
        ax.plot(mids, sub["delta_contact"], color="C3", lw=0.9)
        ax.fill_between(mids, 0, sub["delta_contact"], alpha=0.2, color="C3")
        ax.set_title("|delta contact|")
    else:
        ax.text(0.5, 0.5, "no Hi-C matrix provided",
                ha="center", va="center", transform=ax.transAxes)
    ax.set_xlim(rstart / 1e6, rend / 1e6)

    # rows 2-3: scatter (4 panels)
    pairs = [("p95_bes", "delta_insulation"),
             ("p95_bes", "delta_contact"),
             ("med_bes", "delta_insulation"),
             ("med_bes", "delta_contact")]
    for k, (xc, yc) in enumerate(pairs):
        ax = axes[2 + k // 2, k % 2]
        if yc not in sub.columns or sub[yc].isna().all():
            ax.text(0.5, 0.5, "no data", ha="center", va="center",
                    transform=ax.transAxes)
            ax.set_title("{} vs {} (n/a)".format(xc, yc))
            continue
        x = sub[xc].values
        y = sub[yc].values
        m = np.isfinite(x) & np.isfinite(y)
        ax.scatter(x[m], y[m], s=6, alpha=0.5, edgecolor="none")
        r = results.get("{}__{}".format(xc, yc), {})
        rho = r.get("spearman_rho")
        p = r.get("empirical_p")
        title = "{} vs {}".format(xc, yc)
        if rho is not None:
            title += "  rho={:.3f}".format(rho)
        if p is not None:
            title += "  p_emp={:.3g}".format(p)
        ax.set_title(title, fontsize=9)
        ax.set_xlabel(xc)
        ax.set_ylabel(yc)

    # row 4: null histograms for p95 correlations
    for k, yc in enumerate(["delta_insulation", "delta_contact"]):
        ax = axes[4, k]
        key = "p95_bes__{}".format(yc)
        null = nulls.get(key)
        if null is None:
            ax.text(0.5, 0.5, "no null", ha="center", va="center",
                    transform=ax.transAxes)
            ax.set_title("null for p95_bes vs {} (n/a)".format(yc))
            continue
        ax.hist(null, bins=40, color="gray", alpha=0.7)
        obs = results[key].get("spearman_rho")
        if obs is not None:
            ax.axvline(obs, color="red", lw=2)
        ax.set_title("null: p95_bes vs {}".format(yc), fontsize=9)
        ax.set_xlabel("Spearman rho")

    fig.suptitle("BES vs Hi-C change  {}:{:,}-{:,}".format(*region),
                 fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out_pdf)
    plt.close(fig)


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bes", required=True)
    ap.add_argument("--bes-score-col", default="bearing_score",
                    help="name of BES column in --bes file "
                         "(default: bearing_score; alt: score_normalised, "
                         "bearing_score_tested)")
    ap.add_argument("--insul-A", required=True)
    ap.add_argument("--insul-B", required=True)
    ap.add_argument("--cool-A", default=None)
    ap.add_argument("--cool-B", default=None)
    ap.add_argument("--region", default=None,
                    help="chrom:start-end; required if BES is genome-wide")
    ap.add_argument("--hic-bin", type=int, default=10000)
    ap.add_argument("--min-distance", type=int, default=50000)
    ap.add_argument("--max-distance", type=int, default=500000)
    ap.add_argument("--no-balance", action="store_true",
                    help="skip cooler balance (default: balance IS applied; "
                         "KR weights are in cooler 'weight' column)")
    ap.add_argument("--n-perm", type=int, default=1000)
    ap.add_argument("--min-shift", type=int, default=100000)
    ap.add_argument("--prefix", default="bes_hic_corr")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    balance = not args.no_balance

    # Region
    if args.region:
        region = parse_region(args.region)
    else:
        # Peek at file: skip header if present
        with open(args.bes) as fh:
            first = fh.readline().strip()
            while first.startswith("#") or first == "":
                first = fh.readline().strip()
        has_hdr = first.split("\t")[0].lower() in (
            "chrom", "chr", "#chrom", "seqnames")
        df0 = pd.read_csv(args.bes, sep="\t",
                          header=0 if has_hdr else None,
                          comment="#", usecols=[0, 1, 2],
                          low_memory=False)
        df0.columns = ["chrom", "start", "end"]
        df0["start"] = pd.to_numeric(df0["start"], errors="coerce")
        df0["end"] = pd.to_numeric(df0["end"], errors="coerce")
        df0 = df0.dropna()
        region = (df0.iloc[0]["chrom"], int(df0["start"].min()),
                  int(df0["end"].max()))
        print("Inferred region: {}:{:,}-{:,}".format(*region))

    chrom, rstart, rend = region

    # BES
    print("Loading BES ...")
    bes, track_cols = read_bes_bedgraph(args.bes, region,
                                         score_col=args.bes_score_col)
    print("  {} BEARING bins; {} per-track cols ({})"
          .format(len(bes), len(track_cols), ",".join(track_cols)))

    # Hi-C bin grid
    hic_bins = make_hic_bins(chrom, rstart, rend, args.hic_bin)
    print("Hi-C grid: {} bins of {} bp".format(len(hic_bins), args.hic_bin))

    merged = aggregate_to_hic_bins(bes, hic_bins, args.hic_bin, track_cols)

    # Insulation
    print("Loading insulation ...")
    ins_A, ncol_A = read_bm(args.insul_A, region)
    ins_B, ncol_B = read_bm(args.insul_B, region)
    print("  insul A: {} bins, {} score cols (mean used)"
          .format(len(ins_A), ncol_A))
    print("  insul B: {} bins, {} score cols (mean used)"
          .format(len(ins_B), ncol_B))
    merged = aggregate_bm_to_hic_bins(ins_A, merged, args.hic_bin, "insul_A")
    merged = aggregate_bm_to_hic_bins(ins_B, merged, args.hic_bin, "insul_B")
    merged["delta_insulation"] = (merged["insul_A"] - merged["insul_B"]).abs()

    # Hi-C delta contact
    if args.cool_A and args.cool_B:
        print("Loading .cool matrices ...")
        M_A, bins_A, bs_A = fetch_cool_matrix(args.cool_A, chrom,
                                              rstart, rend, balance)
        M_B, bins_B, bs_B = fetch_cool_matrix(args.cool_B, chrom,
                                              rstart, rend, balance)
        if bs_A != bs_B:
            sys.exit("cool A binsize ({}) != cool B binsize ({})"
                     .format(bs_A, bs_B))
        if bs_A != args.hic_bin:
            warnings.warn("cool binsize ({}) != --hic-bin ({}); "
                          "using cool binsize for delta_contact"
                          .format(bs_A, args.hic_bin))
        if M_A.shape != M_B.shape:
            sys.exit("cool matrix shapes differ: {} vs {}"
                     .format(M_A.shape, M_B.shape))
        delta = per_bin_delta_contact(M_A, M_B, bs_A,
                                       args.min_distance, args.max_distance)
        # Map cool bins to merged grid
        bins_A = bins_A.copy()
        bins_A["delta_contact"] = delta
        bins_A["mid"] = (bins_A["start"] + bins_A["end"]) // 2
        base_start = int(hic_bins.iloc[0]["start"])
        bins_A["bin_idx"] = ((bins_A["mid"] - base_start) // args.hic_bin)\
            .astype(int)
        g = bins_A.groupby("bin_idx")["delta_contact"].sum()
        merged = merged.merge(g.rename("delta_contact"),
                              on="bin_idx", how="left")
    else:
        merged["delta_contact"] = np.nan
        print("Skipping delta_contact (no --cool-A/--cool-B given)")

    # Correlations + null
    min_shift_bins = max(1, args.min_shift // args.hic_bin)
    x_cols = ["p95_bes", "p75_bes", "med_bes", "max_bes"]
    y_cols = ["delta_insulation"]
    if merged["delta_contact"].notna().any():
        y_cols.append("delta_contact")

    results = {}
    nulls = {}
    for xc in x_cols:
        for yc in y_cols:
            cp = corr_pair(merged[xc].values, merged[yc].values)
            obs, null, p_emp = circular_shift_null(
                merged[xc].values, merged[yc].values,
                args.n_perm, min_shift_bins, args.seed)
            cp["empirical_p"] = p_emp
            if null is not None:
                cp["null_mean"] = float(np.nanmean(null))
                cp["null_sd"] = float(np.nanstd(null))
                nulls["{}__{}".format(xc, yc)] = null
            cp["top_decile_enrichment"] = top_decile_enrichment(
                merged[xc].values, merged[yc].values, q=0.9)
            results["{}__{}".format(xc, yc)] = cp

    # Per-track stratification using |kl_i| (kl_i is signed; magnitude
    # is what determines dominance)
    track_results = {}
    if track_cols:
        track_p95_cols = ["{}_abs_p95".format(t) for t in track_cols]
        tot = merged[track_p95_cols].sum(axis=1).replace(0, np.nan)
        for tc, col in zip(track_cols, track_p95_cols):
            frac = merged[col] / tot
            dom = frac > 0.5
            if dom.sum() < 20:
                continue
            sub = merged.loc[dom]
            row = {"n_bins_dominated": int(dom.sum())}
            for yc in y_cols:
                cp = corr_pair(sub["p95_bes"].values, sub[yc].values)
                # permutation null per stratum: a track stratum is a
                # non-contiguous subset of bins, so the circular shift would be
                # degenerate (too few valid shifts at small n) and there is no
                # contiguity to preserve. Reported as p_emp.
                _obs, _null, p_emp = permutation_null(
                    sub["p95_bes"].values, sub[yc].values,
                    args.n_perm, args.seed)
                cp["empirical_p"] = p_emp
                row[yc] = cp
            track_results[tc] = row

    # Outputs
    tsv_path = args.prefix + ".tsv"
    merged.to_csv(tsv_path, sep="\t", index=False)
    print("Wrote per-bin TSV ->", tsv_path)

    summary = {
        "region": "{}:{}-{}".format(*region),
        "hic_bin": args.hic_bin,
        "min_distance": args.min_distance,
        "max_distance": args.max_distance,
        "balance": balance,
        "n_perm": args.n_perm,
        "min_shift_bins": min_shift_bins,
        "n_hic_bins": int(len(merged)),
        "insul_A_score_cols": ncol_A,
        "insul_B_score_cols": ncol_B,
        "correlations": results,
        "track_stratified": track_results,
    }
    with open(args.prefix + ".stats.json", "w") as fh:
        json.dump(summary, fh, indent=2, default=str)
    print("Wrote stats JSON ->", args.prefix + ".stats.json")

    pdf_path = args.prefix + ".scatter.pdf"
    make_plots(merged, region, results, nulls, pdf_path)
    print("Wrote PDF ->", pdf_path)

    # Print summary
    print()
    print("=" * 72)
    print("CORRELATION SUMMARY")
    print("=" * 72)
    for k, v in results.items():
        if v.get("spearman_rho") is None:
            continue
        print()
        print(k)
        print("  n={}  rho={:.3f}  p_param={:.3g}  p_emp={:.3g}".format(
            v["n"], v["spearman_rho"], v["spearman_p"], v["empirical_p"]))
        enr = v.get("top_decile_enrichment")
        if enr:
            print("  top-decile OR={:.2f}  fisher_p={:.3g}".format(
                enr["odds_ratio"], enr["fisher_p"]))

    if track_results:
        print()
        print("PER-TRACK STRATIFIED (bins where one track > 50% of p95 BES):")
        for tc, row in track_results.items():
            print(" ", tc, "n_dom={}".format(row["n_bins_dominated"]))
            for yc, cp in row.items():
                if yc == "n_bins_dominated":
                    continue
                if cp.get("spearman_rho") is None:
                    continue
                print("    {} vs p95_bes  rho={:.3f}  p_param={:.3g}  p_emp={:.3g}".format(
                    yc, cp["spearman_rho"], cp["spearman_p"],
                    cp.get("empirical_p") if cp.get("empirical_p") is not None else float("nan")))


if __name__ == "__main__":
    main()
