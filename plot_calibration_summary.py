#!/usr/bin/env python3
"""
plot_calibration_summary.py -- Supplementary Figure S10 composite: cross-condition
replicate-differential FDR calibration. Reads the aggregated calibration_summary.tsv
(one row per condition; the per-condition bearing_calibration.py outputs feed it) and
draws the three-panel summary the S10 legend describes:

  (A) genomic-inflation lambda per condition (dashed reference at 1.0)
  (B) genome-wide BH-significant bin count per condition
  (C) fraction of bins at the permutation floor (dashed reference at 0.5)

This is the cross-condition summary, distinct from the per-condition four-panel
diagnostic ({cond}_calibration.pdf). Handles thousands separators in the tsv.

  python plot_calibration_summary.py \
    --summary-tsv results/calibration/calibration_summary.tsv \
    --out results/calibration/suppS10_calibration_summary.pdf

ASCII only.
"""
import argparse
import csv
import sys

# canonical condition order; any conditions not listed are appended in file order
CANON = ["DN", "EbKO", "DP", "ProB", "S3T3"]


def _num(x):
    return float(x.replace(",", "").replace('"', "").strip())


def load(path):
    rows = {}
    with open(path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            c = r["condition"].strip()
            try:
                rows[c] = dict(
                    lam=_num(r["lambda"]),
                    floor=_num(r["pct_at_floor"]),
                    bh=_num(r["bh_significant"]),
                )
            except (KeyError, ValueError):
                continue
    return rows


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--summary-tsv", required=True)
    ap.add_argument("--out", required=True, help="output PDF (PNG also written)")
    args = ap.parse_args()

    rows = load(args.summary_tsv)
    if not rows:
        sys.exit("[ERROR] no rows parsed from %s" % args.summary_tsv)

    order = [c for c in CANON if c in rows] + [c for c in rows if c not in CANON]
    lam = [rows[c]["lam"] for c in order]
    bh = [rows[c]["bh"] for c in order]
    floor = [rows[c]["floor"] for c in order]

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    x = np.arange(len(order))
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))

    ax = axes[0]
    ax.bar(x, lam, color="#4c72b0")
    ax.axhline(1.0, ls="--", lw=1.0, color="black")
    ax.set_ylim(0, max(1.05, max(lam) * 1.15))
    ax.set_ylabel("genomic inflation lambda")
    ax.set_title("(A) Inflation lambda (<1 = conservative)")

    ax = axes[1]
    ax.bar(x, bh, color="#55a868")
    ax.set_ylabel("BH-significant bins (q < 0.05)")
    ax.set_title("(B) Replicate-differential BH hits")
    for xi, v in zip(x, bh):
        ax.text(xi, v, "%d" % int(v), ha="center", va="bottom", fontsize=8)

    ax = axes[2]
    ax.bar(x, floor, color="#c44e52")
    ax.axhline(0.5, ls="--", lw=1.0, color="black")
    ax.set_ylim(0, max(0.55, max(floor) * 1.3 if max(floor) > 0 else 0.55))
    ax.set_ylabel("fraction of bins at permutation floor")
    ax.set_title("(C) Permutation-floor fraction")

    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels(order, rotation=0)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)

    fig.tight_layout()
    out = args.out
    fig.savefig(out, dpi=150, bbox_inches="tight")
    png = out[:-4] + ".png" if out.lower().endswith(".pdf") else out + ".png"
    fig.savefig(png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote %s and %s" % (out, png))


if __name__ == "__main__":
    main()
