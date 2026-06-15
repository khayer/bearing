#!/usr/bin/env python3
"""
suppfig_s1_jsd_decomposition.py -- Supplementary Figure S1.

Per-track decomposition of the pairwise Q-vector Jensen-Shannon divergence (JSD).
S1 is the companion to main Figure 2B, which shows the scalar 45-pair JSD heatmap
(the pipeline's q_pair_jsd_heatmap / calibration_observed_q_pair_jsd_heatmap). S1
breaks that scalar into one heatmap per assay track by exploiting the additive
form of the JSD: for a sample pair (A, B) with track-probability vectors q_A, q_B
(each normalized to sum to 1) and the midpoint M = (q_A + q_B) / 2, the per-track
contribution is

    JSD_t(A,B) = 0.5 * q_A[t] * log2(q_A[t]/M[t])
               + 0.5 * q_B[t] * log2(q_B[t]/M[t])

with the convention 0 * log2(0/x) = 0. Each JSD_t is non-negative and, summed
over tracks, equals the scalar JSD(A,B) shown in Figure 2B. Base-2 logarithms put
the total JSD in [0, 1]. The per-track panels therefore SUM to Figure 2B; the
script asserts this reconstruction to within 1e-9 and writes the residuals.

To stay numerically and visually consistent with Figure 2B, the Q-vector loader,
scalar-JSD computation, per-track decomposition, condition colors, and track
palette are imported from compare_qcat.py (the script that produces Fig 2B)
rather than reimplemented.

Outputs (written to --outdir):
  {prefix}.pdf / {prefix}.png        6-panel per-track JSD heatmap (one per track)
  {prefix}_stacked.pdf / .png        stacked-bar per-track decomposition (45 pairs)
  suppS1_jsd_per_track.tsv           per-pair total + per-track contributions +
                                     reconstruction residual

QC thresholds: the scalar JSD warn/error levels are 0.05 / 0.15. S1's per-track
warn level is the warn level split across the tracks (0.05 / 6 = 0.0083 for the
6-track panel); it is drawn on the shared colorbar. No per-track error line is
drawn.

ASCII only. Uses numpy, pandas, matplotlib only; deterministic; no network.
"""
import argparse
import glob
import os
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

# Reuse the exact Fig 2B machinery so S1 reconstructs it by construction.
from compare_qcat import (
    load_sample_sheet,
    load_q_from_cats_json,
    compute_q_pair_jsd_matrix,
    compute_per_track_jsd_decomposition,
    condition_color_map,
    get_default_track_names,
    plot_per_track_jsd_stacked_bar,
    _track_colors_for,
)

SCALAR_WARN = 0.05   # scalar JSD warn threshold (Fig 2B traffic-light)
SCALAR_ERROR = 0.15  # scalar JSD error threshold


def _cats_path_for(base):
    """The _cats.json beside a qcat base path (mirrors load_q_from_cats_json)."""
    b = str(base).replace(".qcat.bgz", "").replace(".bgz", "")
    return b + "_cats.json"


def load_track_names_from_cats(qcat_path, n_tracks):
    """Read assay track names from a sample's _cats.json (numeric-state order).
    Falls back to the canonical default names if the file lacks usable names."""
    import json
    cats_path = _cats_path_for(qcat_path)
    if os.path.exists(cats_path):
        try:
            with open(cats_path) as fh:
                cats = json.load(fh)
            categories = cats.get("categories") if isinstance(cats, dict) else None
            if isinstance(categories, dict) and categories:
                names = []
                for k in sorted(categories.keys(), key=lambda x: int(x)):
                    entry = categories[k]
                    if isinstance(entry, dict):
                        names.append(entry.get("name"))
                    elif isinstance(entry, (list, tuple)) and entry:
                        names.append(entry[0])
                    else:
                        names.append(None)
                if names and all(nm is not None for nm in names) and len(names) == n_tracks:
                    return names
        except (ValueError, KeyError, TypeError):
            pass
    return get_default_track_names(n_tracks)


def collect_samples(args):
    """Return (sample_names, conditions, qcat_paths) from a sheet or a qdir."""
    if args.qcat_sheet:
        samples = load_sample_sheet(args.qcat_sheet)
        names, conds, qcats = [], [], []
        for s in samples:
            qcat = s.get("qcat")
            if qcat is None:
                sys.exit("[ERROR] sheet row lacks a qcat path: %s" % s.get("sample"))
            nm = s.get("sample") or Path(str(qcat)).stem
            names.append(nm)
            conds.append(s["condition"])
            qcats.append(qcat)
        return names, conds, qcats

    # --qdir: every *_cats.json in the directory is one sample.
    cats_files = sorted(glob.glob(os.path.join(args.qdir, "*_cats.json")))
    if not cats_files:
        sys.exit("[ERROR] no *_cats.json files in %s" % args.qdir)
    names, conds, qcats = [], [], []
    for cf in cats_files:
        base = cf[:-len("_cats.json")]
        nm = os.path.basename(base)
        cond = nm.rsplit("_rep", 1)[0] if "_rep" in nm else nm
        names.append(nm)
        conds.append(cond)
        qcats.append(base)  # loader re-appends _cats.json
    return names, conds, qcats


def apply_sample_order(order_str, names, conds, qcats):
    """Reorder samples to a pinned order (comma-separated names) to match Fig 2B."""
    want = [t.strip() for t in order_str.split(",") if t.strip()]
    idx = {nm: i for i, nm in enumerate(names)}
    missing = [w for w in want if w not in idx]
    if missing:
        sys.exit("[ERROR] --sample-order names not found: %s\n  available: %s"
                 % (", ".join(missing), ", ".join(names)))
    if len(want) != len(names):
        sys.exit("[ERROR] --sample-order lists %d samples but %d are present"
                 % (len(want), len(names)))
    order = [idx[w] for w in want]
    return ([names[i] for i in order], [conds[i] for i in order],
            [qcats[i] for i in order])


def per_track_matrices(per_track_records, sample_names, track_names):
    """Build a (K, n, n) symmetric array of per-track JSD from the records.
    Diagonal is NaN (masked in the heatmap)."""
    n = len(sample_names)
    K = len(track_names)
    pos = {nm: i for i, nm in enumerate(sample_names)}
    mats = np.full((K, n, n), np.nan, dtype=np.float64)
    for k in range(K):
        for i in range(n):
            mats[k, i, i] = np.nan
    for rec in per_track_records:
        i, j = pos[rec["sample_A"]], pos[rec["sample_B"]]
        for k, tname in enumerate(track_names):
            v = rec["jsd_per_track"][tname]
            mats[k, i, j] = v
            mats[k, j, i] = v
    return mats


def plot_six_panel(mats, sample_names, conditions, track_names, track_colors,
                   warn_level, out_pdf, out_png, dpi=300):
    """6-panel per-track JSD heatmap, shared sequential colorbar 0..max, with the
    per-track warn level marked on the colorbar. Diagonal masked; cells annotated
    to 3 decimals; tick labels colored by condition; panel titles = track names."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(sample_names)
    K = len(track_names)
    cmap_cond = condition_color_map(conditions)

    # Shared scale: 0 -> global max per-track contribution (off-diagonal).
    finite = mats[np.isfinite(mats)]
    vmax = float(np.nanmax(finite)) if finite.size else 1.0
    if vmax <= 0:
        vmax = 1.0
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#e6e6e6")  # masked diagonal

    ncols = 3
    nrows = int(np.ceil(K / ncols))
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(ncols * 3.6 + 0.6, nrows * 3.3 + 0.6),
                             squeeze=False)

    im = None
    for k in range(K):
        ax = axes[k // ncols][k % ncols]
        mat = mats[k]
        im = ax.imshow(mat, cmap=cmap, vmin=0.0, vmax=vmax, aspect="equal")
        ax.set_xticks(range(n))
        ax.set_yticks(range(n))
        ax.set_xticklabels(sample_names, rotation=45, ha="right", fontsize=7)
        ax.set_yticklabels(sample_names, fontsize=7)
        ax.set_title(track_names[k], fontsize=11, weight="bold",
                     color=track_colors[k])
        for i in range(n):
            for j in range(n):
                if i == j or not np.isfinite(mat[i, j]):
                    txt, tc = "-", "#666666"
                else:
                    txt = "%.3f" % mat[i, j]
                    tc = "white" if mat[i, j] < 0.6 * vmax else "black"
                ax.text(j, i, txt, ha="center", va="center", fontsize=6, color=tc)
        for tick, c in zip(ax.get_xticklabels(),
                           [cmap_cond[cc] for cc in conditions]):
            tick.set_color(c)
        for tick, c in zip(ax.get_yticklabels(),
                           [cmap_cond[cc] for cc in conditions]):
            tick.set_color(c)

    for k in range(K, nrows * ncols):
        axes[k // ncols][k % ncols].set_visible(False)

    # Shared colorbar with the per-track warn level marked.
    cbar_ax = fig.add_axes([0.92, 0.18, 0.016, 0.64])
    cbar = fig.colorbar(im, cax=cbar_ax, label="per-track JSD contribution (bits)")
    if 0.0 <= warn_level <= vmax:
        cbar.ax.axhline(warn_level, color="black", lw=1.2, ls="--")
        cbar.ax.text(1.8, warn_level, "warn %.4f" % warn_level,
                     transform=cbar.ax.get_yaxis_transform(),
                     va="center", ha="left", fontsize=7)

    fig.suptitle("Supplementary Figure S1: per-track Q-vector JSD decomposition\n"
                 "(panels sum to the scalar JSD of Figure 2B)",
                 fontsize=12, y=0.99)
    fig.text(0.5, 0.005,
             "Per-track warn level = scalar warn %.2f / %d tracks = %.4f. "
             "Diagonal masked; tick labels colored by condition."
             % (SCALAR_WARN, K, warn_level),
             ha="center", fontsize=8, style="italic")
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        fig.tight_layout(rect=[0, 0.03, 0.90, 0.95])
    Path(out_pdf).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def write_tsv(per_track_records, track_names, jsd_matrix, sample_names, out_path):
    """Write suppS1_jsd_per_track.tsv: per-pair total + per-track + residual vs the
    directly-computed scalar JSD (the Figure 2B value). Returns max abs residual."""
    pos = {nm: i for i, nm in enumerate(sample_names)}
    rows = []
    max_resid = 0.0
    for rec in sorted(per_track_records, key=lambda r: -r["jsd_total"]):
        i, j = pos[rec["sample_A"]], pos[rec["sample_B"]]
        scalar_direct = float(jsd_matrix[i, j])
        resid = rec["jsd_total"] - scalar_direct
        max_resid = max(max_resid, abs(resid))
        row = {"sample_A": rec["sample_A"], "sample_B": rec["sample_B"],
               "total_jsd": rec["jsd_total"]}
        for t in track_names:
            row["JSD_%s" % t] = rec["jsd_per_track"][t]
        row["scalar_jsd_direct"] = scalar_direct
        row["residual"] = resid
        rows.append(row)
    cols = (["sample_A", "sample_B", "total_jsd"]
            + ["JSD_%s" % t for t in track_names]
            + ["scalar_jsd_direct", "residual"])
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=cols).to_csv(out_path, sep="\t", index=False,
                                            float_format="%.10g")
    return max_resid


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--qcat-sheet", help="sample sheet TSV (sample, condition, qcat)")
    src.add_argument("--qdir", help="directory of per-sample *_cats.json files")
    ap.add_argument("--outdir", required=True, help="output directory (paper_figures)")
    ap.add_argument("--prefix", default="suppS1_jsd_decomposition",
                    help="output filename prefix (default suppS1_jsd_decomposition)")
    ap.add_argument("--sample-order", default=None,
                    help="comma-separated sample names to pin ordering to Fig 2B")
    ap.add_argument("--dpi", type=int, default=300, help="raster DPI (default 300)")
    ap.add_argument("--tol", type=float, default=1e-9,
                    help="max allowed reconstruction residual (default 1e-9)")
    ap.add_argument("--fig2b-tsv", default=None,
                    help="optional q_pair_jsd.tsv to cross-check scalar values")
    args = ap.parse_args()

    names, conds, qcats = collect_samples(args)
    if args.sample_order:
        names, conds, qcats = apply_sample_order(args.sample_order, names, conds, qcats)

    per_sample_q = [load_q_from_cats_json(q) for q in qcats]
    bad = [names[i] for i, q in enumerate(per_sample_q) if q is None]
    if bad:
        sys.exit("[ERROR] could not load Q vector for: %s" % ", ".join(bad))
    n_tracks = len(per_sample_q[0])
    if any(len(q) != n_tracks for q in per_sample_q):
        sys.exit("[ERROR] samples have differing track counts")

    track_names = load_track_names_from_cats(qcats[0], n_tracks)
    track_colors = _track_colors_for(track_names)
    warn_level = SCALAR_WARN / max(1, n_tracks)
    n_pairs = len(list(combinations(range(len(names)), 2)))

    print("S1 decomposition: %d samples, %d tracks, %d pairs"
          % (len(names), n_tracks, n_pairs))
    print("  tracks: %s" % ", ".join(track_names))

    # Scalar JSD (Fig 2B) and per-track decomposition, both from compare_qcat.
    jsd_matrix, _scalar_records = compute_q_pair_jsd_matrix(per_sample_q, names)
    per_track_records = compute_per_track_jsd_decomposition(
        per_sample_q, names, track_names=track_names)

    outdir = Path(args.outdir)
    out_pdf = outdir / ("%s.pdf" % args.prefix)
    out_png = outdir / ("%s.png" % args.prefix)
    tsv_path = outdir / "suppS1_jsd_per_track.tsv"

    # --- Validation: per-track panels must sum to the scalar Fig 2B JSD ---------
    max_resid = write_tsv(per_track_records, track_names, jsd_matrix, names, tsv_path)
    status = "PASS" if max_resid < args.tol else "FAIL"
    print("VALIDATION [%s]: max |sum_t JSD_t - scalar JSD| = %.3e (tol %.1e)"
          % (status, max_resid, args.tol))

    if args.fig2b_tsv and os.path.exists(args.fig2b_tsv):
        ext = pd.read_csv(args.fig2b_tsv, sep="\t")
        key = {frozenset((r["sample_A"], r["sample_B"])): float(r["jsd"])
               for _, r in ext.iterrows()}
        pos = {nm: i for i, nm in enumerate(names)}
        d = 0.0
        for rec in per_track_records:
            kk = frozenset((rec["sample_A"], rec["sample_B"]))
            if kk in key:
                d = max(d, abs(rec["jsd_total"] - key[kk]))
        print("CROSS-CHECK vs Fig 2B TSV: max |S1 total - Fig2B scalar| = %.3e" % d)

    if max_resid >= args.tol:
        sys.exit("[ERROR] reconstruction residual exceeds tolerance; S1 does not "
                 "sum to Fig 2B -- not writing figures")

    # --- Figures ----------------------------------------------------------------
    mats = per_track_matrices(per_track_records, names, track_names)
    plot_six_panel(mats, names, conds, track_names, track_colors, warn_level,
                   out_pdf, out_png, dpi=args.dpi)
    print("  wrote %s" % out_pdf)
    print("  wrote %s" % out_png)

    stacked_pdf = outdir / ("%s_stacked.pdf" % args.prefix)
    plot_per_track_jsd_stacked_bar(per_track_records, track_names, stacked_pdf,
                                   warn_threshold=SCALAR_WARN,
                                   track_colors=track_colors)
    print("  wrote %s" % stacked_pdf)
    print("  wrote %s" % tsv_path)


if __name__ == "__main__":
    main()
