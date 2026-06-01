#!/usr/bin/env python3
"""
plot_bearing_sweep.py
=====================
Turn the coarse-grid sweep_master.tsv (from run_bearing_recovery_sweep.sh) into
figures. The grid spans SNR x NB-dispersion x pseudocount; this script shows
how the recovery metrics move (or do not) across it.

Three things are plotted:
  1. Recovery metrics vs SNR, one line per NB-dispersion, for each metric
     (per-block slope vs expected, per-block Pearson, rank recovery). This is
     the main view: does recovery hold up as signal-to-noise and noisiness
     change.
  2. Informative-block RAD21+CTCF fraction vs SNR (with its min-max band per
     cell shown as error bars), one line per dispersion. This is the
     composition-stability view tied to the real 55-71% claim.
  3. A pseudocount-invariance check: for every (SNR, dispersion) cell, the
     spread of each metric across the pseudocount values. If pseudocount has no
     effect (as expected for the zero-clamp epsilon), every bar is ~0.

USAGE
    python3 plot_bearing_sweep.py --master sweep_master.tsv --out-prefix sweep

DEPENDENCIES
    numpy, matplotlib
"""

import argparse
import csv
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


METRICS = [
    ("perblock_slope_vs_expected", "per-block slope vs expected", (0.0, 1.1)),
    ("perblock_pearson_vs_expected", "per-block Pearson vs expected", (0.0, 1.0)),
    ("rank_recovery", "rank recovery (top-track match)", (0.0, 1.0)),
]


def load_master(path):
    rows = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                rows.append({
                    "snr": float(r["snr"]),
                    "disp": float(r["disp"]),
                    "pseudo": float(r["pseudocount"]),
                    "slope": float(r["perblock_slope_vs_expected"]),
                    "pearson": float(r["perblock_pearson_vs_expected"]),
                    "rank": float(r["rank_recovery"]),
                    "frac_min": float(r["informative_frac_min"]),
                    "frac_max": float(r["informative_frac_max"]),
                })
            except (ValueError, KeyError):
                # skip NA / malformed rows but keep going
                continue
    return rows


def metric_key(name):
    return {"perblock_slope_vs_expected": "slope",
            "perblock_pearson_vs_expected": "pearson",
            "rank_recovery": "rank"}[name]


def mean_over_pseudo(rows, snr, disp, key):
    """Mean of a metric over all pseudocount values at fixed (snr, disp)."""
    vals = [r[key] for r in rows if r["snr"] == snr and r["disp"] == disp]
    return float(np.mean(vals)) if vals else float("nan")


def spread_over_pseudo(rows, snr, disp, key):
    """Max-min of a metric across pseudocount values at fixed (snr, disp)."""
    vals = [r[key] for r in rows if r["snr"] == snr and r["disp"] == disp]
    return (max(vals) - min(vals)) if vals else float("nan")


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Plot BEARING sweep recovery metrics from sweep_master.tsv.")
    ap.add_argument("--master", required=True, help="sweep_master.tsv")
    ap.add_argument("--out-prefix", default="sweep",
                    help="output filename prefix (default: sweep)")
    args = ap.parse_args(argv)

    rows = load_master(args.master)
    if not rows:
        raise SystemExit("ERROR: no usable rows parsed from %s" % args.master)

    snrs = sorted(set(r["snr"] for r in rows))
    disps = sorted(set(r["disp"] for r in rows))
    pseudos = sorted(set(r["pseudo"] for r in rows))
    disp_colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(disps)))

    # ----- Figure 1: metrics vs SNR, one line per dispersion -----------------
    fig1, axes = plt.subplots(1, len(METRICS), figsize=(5.0 * len(METRICS), 4.4))
    if len(METRICS) == 1:
        axes = [axes]
    for ax, (mname, mlabel, ylim) in zip(axes, METRICS):
        key = metric_key(mname)
        for di, disp in enumerate(disps):
            ys = [mean_over_pseudo(rows, s, disp, key) for s in snrs]
            ax.plot(snrs, ys, "o-", color=disp_colors[di],
                    label="dispersion=%g" % disp)
        ax.set_xlabel("SNR")
        ax.set_ylabel(mlabel)
        ax.set_title(mlabel)
        ax.set_ylim(*ylim)
        ax.set_xscale("log")
        ax.set_xticks(snrs)
        ax.set_xticklabels([("%g" % s) for s in snrs])
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
    fig1.suptitle("Recovery metrics across the SNR x dispersion grid "
                  "(mean over pseudocount; ideal slope ~1, high Pearson/rank)",
                  fontsize=12)
    fig1.tight_layout(rect=[0, 0, 1, 0.95])
    out1 = "%s_metrics_vs_snr.png" % args.out_prefix
    fig1.savefig(out1, dpi=150)
    plt.close(fig1)

    # ----- Figure 2: informative-block RAD21+CTCF fraction vs SNR ------------
    fig2, ax = plt.subplots(figsize=(6.4, 4.6))
    for di, disp in enumerate(disps):
        ys, lo, hi = [], [], []
        for s in snrs:
            cells = [r for r in rows if r["snr"] == s and r["disp"] == disp]
            if not cells:
                ys.append(np.nan); lo.append(np.nan); hi.append(np.nan); continue
            fmin = np.mean([c["frac_min"] for c in cells])
            fmax = np.mean([c["frac_max"] for c in cells])
            mid = 0.5 * (fmin + fmax)
            ys.append(mid); lo.append(mid - fmin); hi.append(fmax - mid)
        ax.errorbar(snrs, ys, yerr=[lo, hi], fmt="o-", color=disp_colors[di],
                    capsize=3, label="dispersion=%g" % disp)
    ax.set_xlabel("SNR")
    ax.set_ylabel("informative-block RAD21+CTCF fraction")
    ax.set_title("Composition stability of the informative block\n"
                 "(error bars = within-cell min-max; brackets the real "
                 "55-71% split)")
    ax.set_ylim(0, 1)
    ax.set_xscale("log")
    ax.set_xticks(snrs)
    ax.set_xticklabels([("%g" % s) for s in snrs])
    ax.axhspan(0.55, 0.71, color="0.85", zorder=0, label="real 55-71% band")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig2.tight_layout()
    out2 = "%s_composition_vs_snr.png" % args.out_prefix
    fig2.savefig(out2, dpi=150)
    plt.close(fig2)

    # ----- Figure 3: pseudocount-invariance check ----------------------------
    # For each metric, the max spread across pseudocount over all (snr, disp).
    fig3, ax = plt.subplots(figsize=(6.4, 4.2))
    metric_names = [m[1] for m in METRICS]
    max_spreads = []
    for mname, mlabel, _ in METRICS:
        key = metric_key(mname)
        spreads = [spread_over_pseudo(rows, s, d, key)
                   for s in snrs for d in disps]
        spreads = [v for v in spreads if v == v]
        max_spreads.append(max(spreads) if spreads else float("nan"))
    bars = ax.bar(range(len(METRICS)), max_spreads, color="#4c78a8")
    ax.set_xticks(range(len(METRICS)))
    ax.set_xticklabels(metric_names, rotation=15, ha="right", fontsize=8)
    ax.set_ylabel("max spread across pseudocount\n(over all SNR x dispersion cells)")
    ax.set_title("Pseudocount (zero-clamp epsilon) invariance\n"
                 "(bars ~0 => recovery does not depend on pseudocount; "
                 "pseudo values tested: %s)"
                 % ", ".join("%g" % p for p in pseudos))
    ymax = max([v for v in max_spreads if v == v] + [0.01])
    ax.set_ylim(0, max(0.05, ymax * 1.3))
    for b, v in zip(bars, max_spreads):
        ax.text(b.get_x() + b.get_width() / 2, v + ymax * 0.04,
                "%.4f" % v, ha="center", fontsize=9)
    fig3.tight_layout()
    out3 = "%s_pseudocount_invariance.png" % args.out_prefix
    fig3.savefig(out3, dpi=150)
    plt.close(fig3)

    # ----- console summary ----------------------------------------------------
    print("Parsed %d cells: SNR=%s, dispersion=%s, pseudocount=%s"
          % (len(rows), snrs, disps, pseudos))
    print("Wrote:")
    for f in (out1, out2, out3):
        print("  %s" % f)
    # report the pseudocount invariance explicitly
    print("\nPseudocount invariance (max spread across eps, over all cells):")
    for (mname, mlabel, _), sp in zip(METRICS, max_spreads):
        print("  %-32s %.4f" % (mlabel + ":", sp))


if __name__ == "__main__":
    main()
