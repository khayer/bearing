#!/usr/bin/env python3
"""
plot_1d_vs_3d.py -- make the anti-localization visible. Reads a --dump-bins file
(chrom start end bes delta_contact on_feature width marginal) and draws two
stacked locus tracks over a window:

  top    BES (1D differential)        -- "where 1D moved"
  bottom delta_contact (3D change)    -- "where 3D moved"

The top-K bins of each track are marked, and EACH track's hotspots are projected
as faint vertical lines onto the OTHER track -- so a reader sees the 1D peaks
land where the 3D track is low, and vice versa. Title carries the raw in-window
Spearman and (if given) the depth-controlled partial.

  python plot_1d_vs_3d.py --dump mcc_grid_benchmark/binwise_DN_vs_V1P.tsv \
    --region chr6:40800000-41650000 --features annotations/cbe_mm10.bed \
    --highlight chr6:41504000:41507000 --partial-rho -0.30 \
    --title "DN vs V1P -- Tcrb MCC" --out fig_1d_vs_3d_V1P.pdf

ASCII only. numpy + matplotlib.
"""
import argparse
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def rankdata(a):
    a = np.asarray(a, float)
    order = np.argsort(a, kind="mergesort")
    r = np.empty(len(a)); r[order] = np.arange(1, len(a) + 1)
    return r


def spearman(a, b):
    ra, rb = rankdata(a) - rankdata(a).mean(), rankdata(b) - rankdata(b).mean()
    d = np.sqrt((ra * ra).sum() * (rb * rb).sum())
    return float((ra * rb).sum() / d) if d > 0 else float("nan")


def partial_spearman(a, b, ctrl):
    """partial Spearman(a,b | ctrl), all rank-transformed."""
    ra, rb, rc = rankdata(a), rankdata(b), rankdata(ctrl)
    X = np.column_stack([np.ones(len(rc)), rc])
    def res(y):
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        return y - X @ beta
    return spearman_resid(res(ra), res(rb))


def spearman_resid(x, y):
    x = x - x.mean(); y = y - y.mean()
    d = np.sqrt((x * x).sum() * (y * y).sum())
    return float((x * y).sum() / d) if d > 0 else float("nan")


def load_region(path, chrom, a, b):
    rows = [l.rstrip("\n").split("\t") for l in open(path)]
    idx = {h: i for i, h in enumerate(rows[0])}
    has_marg = "marginal" in idx
    out = []
    n_drop = 0
    for r in rows[1:]:
        if r[idx["chrom"]] != chrom:
            continue
        s, e = int(r[idx["start"]]), int(r[idx["end"]])
        if e <= a or s >= b:
            continue
        bes = float(r[idx["bes"]])
        delta = float(r[idx["delta_contact"]])
        m = float(r[idx["marginal"]]) if has_marg else float("nan")
        # drop non-finite bins (match partial_corr: bes, delta, marginal finite)
        if not (np.isfinite(bes) and np.isfinite(delta) and np.isfinite(m)):
            n_drop += 1
            continue
        out.append(((s + e) / 2.0, bes, delta, m))
    out.sort()
    if n_drop:
        sys.stderr.write("[load] dropped %d non-finite bins\n" % n_drop)
    c = np.array([o[0] for o in out])
    return (c, np.array([o[1] for o in out]), np.array([o[2] for o in out]),
            np.array([o[3] for o in out]))


def load_bed(path, chrom, a, b):
    iv = []
    if not path:
        return iv
    for line in open(path):
        if not line.strip() or line.startswith(("#", "track")):
            continue
        f = line.split()
        if f[0] != chrom:
            continue
        s, e = int(f[1]), int(f[2])
        if e > a and s < b:
            iv.append((s, e))
    return iv


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dump", required=True)
    ap.add_argument("--region", required=True)
    ap.add_argument("--features", default=None)
    ap.add_argument("--highlight", default=None, help="chrom:start:end, e.g. V1 promoter")
    ap.add_argument("--topk", type=int, default=15)
    ap.add_argument("--partial-rho", type=float, default=None)
    ap.add_argument("--title", default="1D vs 3D")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    chrom, rng = args.region.split(":")
    a, b = (int(x) for x in rng.split("-"))
    c, bes, delta, marg = load_region(args.dump, chrom, a, b)
    if len(c) == 0:
        sys.exit("[ERROR] no bins in region")
    feats = load_bed(args.features, chrom, a, b)
    rho = spearman(bes, delta)
    pr = (partial_spearman(bes, delta, marg)
          if np.isfinite(marg).all() and np.ptp(marg) > 0 else float("nan"))

    kb = c / 1e6
    top_b = c[np.argsort(-bes)[:args.topk]] / 1e6
    top_d = c[np.argsort(-delta)[:args.topk]] / 1e6

    fig, (axb, axd) = plt.subplots(2, 1, figsize=(11, 5.2), sharex=True,
                                   gridspec_kw=dict(hspace=0.12))
    GOLD, BLUE = "#d4a017", "#1f6fb2"

    # shared decorations
    for ax in (axb, axd):
        if args.highlight:
            hc, hs, he = args.highlight.split(":")
            if hc == chrom:
                ax.axvspan(int(hs) / 1e6, int(he) / 1e6, color="#ffe08a",
                           alpha=0.5, lw=0, zorder=0)
        for fs, fe in feats:
            ax.axvspan(fs / 1e6, fe / 1e6, color="#cccccc", alpha=0.6, lw=0, zorder=0)

    # top: BES, with 3D hotspots projected as faint blue lines
    axb.fill_between(kb, bes, color="#333333", step="mid", alpha=0.85, lw=0)
    for x in top_d:
        axb.axvline(x, color=BLUE, alpha=0.35, lw=0.8, zorder=1)
    axb.plot(top_b, bes[np.argsort(-bes)[:args.topk]], "v", color=GOLD,
             ms=7, mec="k", mew=0.4, zorder=3, label="top 1D bins")
    axb.set_ylabel("BES (1D)\n-log10 p")
    axb.legend(loc="upper right", fontsize=8, framealpha=0.9)
    axb.set_title("%s   |   in-window Spearman=%.3f   partial|depth=%s"
                  % (args.title, rho,
                     ("%.3f" % pr) if np.isfinite(pr)
                     else ("%.2f*" % args.partial_rho
                           if args.partial_rho is not None else "n/a")),
                  fontsize=11)

    # bottom: delta_contact, with 1D hotspots projected as faint gold lines
    axd.fill_between(kb, delta, color=BLUE, step="mid", alpha=0.7, lw=0)
    for x in top_b:
        axd.axvline(x, color=GOLD, alpha=0.5, lw=0.8, zorder=1)
    axd.plot(top_d, delta[np.argsort(-delta)[:args.topk]], "^", color=BLUE,
             ms=7, mec="k", mew=0.4, zorder=3, label="top 3D bins")
    axd.set_ylabel("dContact (3D)\nband-limited")
    axd.set_xlabel("%s position (Mb)" % chrom)
    axd.legend(loc="upper right", fontsize=8, framealpha=0.9)

    fig.text(0.013, 0.5, "gold = 1D hotspots   blue = 3D hotspots   "
             "(note they fall on different bins)", rotation=90, va="center",
             fontsize=8, color="#555555")
    fig.savefig(args.out, bbox_inches="tight", dpi=200)
    sys.stderr.write("[fig] %s  (%d bins, in-window rho=%.3f, partial|depth=%.3f)\n"
                     % (args.out, len(c), rho, pr))


if __name__ == "__main__":
    main()
