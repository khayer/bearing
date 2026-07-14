#!/usr/bin/env python3
"""
plot_bearing_stats_summary.py
=============================
Summarize per-condition differences vs a reference (DN) at a locus into one
overview figure with three rows per comparison:

  1. flagged-bin counts  - diverging stacked bars; up = enriched in the
     reference, down = enriched in the condition; light = top-1% tier,
     solid = top-0.1% tier.
  2. summed |KL| per track - total magnitude of difference attributable to
     each chromatin track over the region (magnitude, not count).
  3. median signed KL per track - the typical per-bin effect size and its
     direction (robust to outliers).

Rows 2 and 3 require the FULL diff stats TSVs (all per-bin kl_* values), so
pass them with --diff LABEL=PATH. The lighter --stats mode reads the
companion *_stats.tsv sidecars but can only draw row 1 (counts), since the
sidecars contain only the flagged extremes.

Usage (full, recommended):
  plot_bearing_stats_summary.py \\
    --diff DP=results_v6_paper/diff_DN_vs_DP.stats.tsv \\
    --diff EbKO=results_v6_paper/diff_DN_vs_EbKO.stats.tsv \\
    --diff ProB=results_v6_paper/diff_DN_vs_ProB.stats.tsv \\
    --diff 3T3=results_v6_paper/diff_DN_vs_S3T3.stats.tsv \\
    --region chr6:40400000-42550000 \\
    --categories DN_rep1_cats.json \\
    --title "TCRb locus: difference vs DN" \\
    --out bearing_stats_summary.pdf

Usage (counts-only from sidecars):
  plot_bearing_stats_summary.py \\
    --stats tcrb_wide_2Mb_DN_vs_DP_combined.stats.tsv \\
    --stats tcrb_wide_2Mb_DN_vs_EbKO_combined.stats.tsv \\
    --categories DN_rep1_cats.json --out summary_counts.pdf
"""

import argparse
import csv
import json
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

AUTOSOMES = {"chr{}".format(i) for i in range(1, 23)}
DEFAULT_ORDER = ["ATAC", "RNAseq +", "RNAseq -", "CTCF", "Cohesin", "H3K27ac"]
DEFAULT_COLORS = {
    "ATAC": "#b39ddb", "RNAseq +": "#6495ed", "RNAseq -": "#1a3a8f",
    "CTCF": "#ff2200", "Cohesin": "#8b0000", "H3K27ac": "#00c864",
}


def load_categories_json(path):
    with open(path) as f:
        d = json.load(f)
    cats = d["categories"]
    pairs = [(cats[k][0], cats[k][1]) for k in sorted(cats.keys(), key=int)]
    return [p[0] for p in pairs], {p[0]: p[1] for p in pairs}


def parse_region(s):
    m = re.match(r"^([^:]+):([0-9,]+)-([0-9,]+)$", s)
    if not m:
        raise SystemExit("Cannot parse --region: " + s)
    return m.group(1), int(m.group(2).replace(",", "")), int(m.group(3).replace(",", ""))


def _detect_kl_columns(path):
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
    kl = sorted([c for c in header if c.startswith("kl_")], key=lambda c: int(c.split("_")[1]))
    return header, kl


def metrics_from_full_diff(path, names, region, pctile_low, pctile_high,
                           scope, diff_mult):
    """Per-track region metrics from a full diff stats TSV.
    Returns dict track -> {pos_low,pos_high,neg_low,neg_high,sum_abs,median}."""
    chrom, start, end = region
    header, kl_cols = _detect_kl_columns(path)
    if not kl_cols:
        raise SystemExit("{}: no kl_* columns".format(path))
    n_use = min(len(kl_cols), len(names))
    kl_cols = kl_cols[:n_use]
    use_names = names[:n_use]
    want = [c for c in ("chrom", "start", "end") if c in header] + kl_cols
    df = pd.read_csv(path, sep="\t", usecols=want, low_memory=False)

    if scope == "genome":
        df_ref = df[df["chrom"].isin(AUTOSOMES)]
    else:
        df_ref = df[df["chrom"] == chrom]
    mask = (df["chrom"] == chrom) & (df["end"] > start) & (df["start"] < end)
    df_r = df[mask]

    out = {}
    for col, name in zip(kl_cols, use_names):
        ref_abs = df_ref[col].abs().to_numpy()
        ref_abs = ref_abs[np.isfinite(ref_abs)]
        thr_low = float(np.nanpercentile(ref_abs, pctile_low)) if len(ref_abs) else np.inf
        thr_high = float(np.nanpercentile(ref_abs, pctile_high)) if len(ref_abs) else np.inf
        vals = diff_mult * df_r[col].to_numpy(dtype=float)
        av = np.abs(vals)
        low = (av >= thr_low) & (av < thr_high)
        high = av >= thr_high
        active = np.isfinite(vals) & (av > 1e-9)   # non-masked bins
        vnz = vals[active]
        out[name] = {
            "pos_low": int(np.sum(low & (vals > 0))),
            "pos_high": int(np.sum(high & (vals > 0))),
            "neg_low": int(np.sum(low & (vals < 0))),
            "neg_high": int(np.sum(high & (vals < 0))),
            "n_active": int(active.sum()),
            "sum_abs": float(np.nansum(av)),
            # Effect-size summaries computed over ACTIVE bins so the sea of
            # masked zeros does not drag everything to ~0.
            "median_signed_nz": float(np.median(vnz)) if vnz.size else 0.0,
            "median_abs_nz": float(np.median(np.abs(vnz))) if vnz.size else 0.0,
            "mean_signed": float(np.nanmean(vals)) if len(vals) else 0.0,
        }
    return out


def counts_from_sidecar(path):
    rows = list(csv.DictReader(open(path, newline=""), delimiter="\t"))
    rows = [{k: (v.strip() if isinstance(v, str) else v) for k, v in r.items()} for r in rows]
    cond = None
    per = {}
    for r in rows:
        if r.get("feature") != "kl_track":
            continue
        if cond is None and " - " in r.get("direction", ""):
            cond = r["direction"].split(" - ", 1)[1]
        t = r["track"]
        d = per.setdefault(t, {"pos_low": 0, "pos_high": 0, "neg_low": 0,
                               "neg_high": 0, "sum_abs": float("nan"), "median": float("nan")})
        try:
            v = float(r["value"])
        except (ValueError, KeyError):
            continue
        high = "0.1%" in r.get("tier", "")
        sign = "pos" if v >= 0 else "neg"
        d["{}_{}".format(sign, "high" if high else "low")] += 1
    return cond, per


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--diff", action="append", default=[], metavar="LABEL=PATH",
                    help="Full diff stats TSV per comparison (enables all 3 rows).")
    ap.add_argument("--stats", action="append", default=[], metavar="TSV",
                    help="Sidecar *_stats.tsv (counts-only fallback).")
    ap.add_argument("--region", default=None, help="chr:start-end (required with --diff).")
    ap.add_argument("--categories", default=None, help="Cats JSON for colors/order.")
    ap.add_argument("--reference", default="DN")
    ap.add_argument("--track-pctile-low", type=float, default=99.0)
    ap.add_argument("--track-pctile-high", type=float, default=99.9)
    ap.add_argument("--threshold-scope", choices=["genome", "chrom"], default="genome")
    ap.add_argument("--diff-direction", choices=["a_minus_b", "b_minus_a"],
                    default="a_minus_b")
    ap.add_argument("--effect-metric",
                    choices=["median_signed_nz", "median_abs_nz", "mean_signed"],
                    default="median_signed_nz",
                    help="Row 3 effect-size metric, computed over ACTIVE (non-zero) "
                         "bins so masked zeros do not flatten it. "
                         "median_signed_nz (default): typical signed KL where the "
                         "track is active (direction-aware). median_abs_nz: typical "
                         "magnitude where active. mean_signed: net signed mean over "
                         "all region bins.")
    ap.add_argument("--title", default=None)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    if args.categories:
        order, colors = load_categories_json(args.categories)
    else:
        order, colors = DEFAULT_ORDER, DEFAULT_COLORS
    diff_mult = 1.0 if args.diff_direction == "a_minus_b" else -1.0

    parsed = []  # list of (cond_label, per_track_dict)
    full_mode = bool(args.diff)
    if full_mode:
        if not args.region:
            raise SystemExit("--region is required with --diff")
        region = parse_region(args.region)
        for item in args.diff:
            if "=" not in item:
                raise SystemExit("--diff expects LABEL=PATH: " + item)
            label, path = item.split("=", 1)
            per = metrics_from_full_diff(path, order, region,
                                         args.track_pctile_low, args.track_pctile_high,
                                         args.threshold_scope, diff_mult)
            parsed.append((label, per))
    else:
        if not args.stats:
            raise SystemExit("Provide --diff (full) or --stats (sidecar) inputs.")
        for path in args.stats:
            cond, per = counts_from_sidecar(path)
            parsed.append((cond or path, per))

    seen = set()
    for _, per in parsed:
        seen.update(per.keys())
    tracks = [t for t in order if t in seen] + [t for t in sorted(seen) if t not in order]
    track_colors = [colors.get(t, "#777777") for t in tracks]
    x = np.arange(len(tracks))
    n_cond = len(parsed)

    n_rows = 3 if full_mode else 1
    height = (2.9 if full_mode else 3.4) * (n_rows / 3.0 * 3) + 1.2
    fig, axes = plt.subplots(
        n_rows, n_cond, figsize=(2.8 * n_cond + 1.0, 2.7 * n_rows + 0.8),
        squeeze=False, gridspec_kw={"hspace": 0.5, "wspace": 0.28})

    def shared_lim(getter):
        m = 1e-9
        for _, per in parsed:
            for t in tracks:
                m = max(m, abs(getter(per.get(t, {}))))
        return m

    cmax = shared_lim(lambda d: max(d.get("pos_low", 0) + d.get("pos_high", 0),
                                    d.get("neg_low", 0) + d.get("neg_high", 0)))
    cmax = int(cmax * 1.15) + 1
    eff_key = args.effect_metric
    eff_signed = eff_key in ("median_signed_nz", "mean_signed")
    eff_labels = {
        "median_signed_nz": "median KL effect\n(active bins; up {} / down cond)".format(args.reference),
        "median_abs_nz": "median |KL|\n(active bins)",
        "mean_signed": "mean KL effect\n(region; up {} / down cond)".format(args.reference),
    }
    if full_mode:
        smax = shared_lim(lambda d: d.get("sum_abs", 0.0)) * 1.15
        emax = shared_lim(lambda d: abs(d.get(eff_key, 0.0))) * 1.15
        emax = max(emax, 1e-6)

    for ci, (cond, per) in enumerate(parsed):
        # Row 1: diverging flagged-bin counts
        ax = axes[0][ci]
        for xi, (t, col) in enumerate(zip(tracks, track_colors)):
            d = per.get(t, {})
            pl, ph = d.get("pos_low", 0), d.get("pos_high", 0)
            nl, nh = d.get("neg_low", 0), d.get("neg_high", 0)
            ax.bar(xi, pl, 0.74, color=col, alpha=0.45, edgecolor="none")
            ax.bar(xi, ph, 0.74, bottom=pl, color=col, alpha=1.0, edgecolor="black", linewidth=0.3)
            ax.bar(xi, -nl, 0.74, color=col, alpha=0.45, edgecolor="none")
            ax.bar(xi, -nh, 0.74, bottom=-nl, color=col, alpha=1.0, edgecolor="black", linewidth=0.3)
        ax.axhline(0, color="black", lw=0.7)
        ax.set_ylim(-cmax, cmax)
        ax.set_title("{} vs {}".format(args.reference, cond), fontsize=10, fontweight="bold")
        _style(ax, tracks, x, labels=(not full_mode))
        if ci == 0:
            ax.set_ylabel("flagged bins\n(up: {}  down: cond)".format(args.reference), fontsize=8)

        if not full_mode:
            continue

        # Row 2: summed |KL|
        ax2 = axes[1][ci]
        for xi, (t, col) in enumerate(zip(tracks, track_colors)):
            ax2.bar(xi, per.get(t, {}).get("sum_abs", 0.0), 0.74, color=col,
                    edgecolor="black", linewidth=0.3)
        ax2.set_ylim(0, smax)
        _style(ax2, tracks, x, labels=False)
        if ci == 0:
            ax2.set_ylabel("summed |KL|\n(region)", fontsize=8)

        # Row 3: effect size (computed over active bins)
        ax3 = axes[2][ci]
        for xi, (t, col) in enumerate(zip(tracks, track_colors)):
            ax3.bar(xi, per.get(t, {}).get(eff_key, 0.0), 0.74, color=col,
                    edgecolor="black", linewidth=0.3)
        if eff_signed:
            ax3.axhline(0, color="black", lw=0.7)
            ax3.set_ylim(-emax, emax)
        else:
            ax3.set_ylim(0, emax)
        _style(ax3, tracks, x, labels=True)
        if ci == 0:
            ax3.set_ylabel(eff_labels[eff_key], fontsize=8)

    handles = [
        Patch(facecolor="gray", alpha=1.0, edgecolor="black", label="top 0.1% (solid)"),
        Patch(facecolor="gray", alpha=0.45, edgecolor="none", label="top 1% (light)"),
    ]
    axes[0][-1].legend(handles=handles, loc="upper right", fontsize=6.5, framealpha=0.85)

    fig.suptitle(args.title or "Difference vs {}".format(args.reference),
                 fontsize=12, y=0.995)
    fig.savefig(args.out, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("Wrote", args.out)


def _style(ax, tracks, x, labels):
    ax.set_xticks(x)
    if labels:
        ax.set_xticklabels(tracks, rotation=45, ha="right", fontsize=7)
    else:
        ax.set_xticklabels([])
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


if __name__ == "__main__":
    main()
