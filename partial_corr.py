#!/usr/bin/env python3
"""
partial_corr.py -- is the BES <-> delta_contact relationship real, or a
coverage/width confound? Reads a --dump-bins file (chrom start end bes
delta_contact on_feature width marginal) and reports:

  rho_raw                Spearman(bes, delta)                      (what the benchmark reports)
  rho | depth            partial Spearman controlling marginal     (kills a coverage artifact)
  rho | depth,width      partial controlling marginal AND width    (kills width coupling too)
  rho (well-covered)     Spearman(bes, delta) on top-50% marginal  (artifact can't hide here)

If the negative survives "| depth" and the well-covered subset, it is a real
1D-3D relationship. If it collapses toward 0, it was the capture enriching
contacts on-target while high-BES bins sat off-target (low depth -> low |delta|).

  python partial_corr.py mcc_grid_benchmark/binwise_DN_vs_V1P.tsv
ASCII only. numpy only.
"""
import sys
import numpy as np


def rankdata(a):
    a = np.asarray(a, float)
    order = np.argsort(a, kind="mergesort")
    ranks = np.empty(len(a), float)
    ranks[order] = np.arange(1, len(a) + 1)
    # average ties
    sa = a[order]
    i = 0
    while i < len(sa):
        j = i
        while j + 1 < len(sa) and sa[j + 1] == sa[i]:
            j += 1
        if j > i:
            ranks[order[i:j + 1]] = (i + 1 + j + 1) / 2.0
        i = j + 1
    return ranks


def resid(y, X):
    X1 = np.column_stack([np.ones(len(y)), X])
    beta, *_ = np.linalg.lstsq(X1, y, rcond=None)
    return y - X1 @ beta


def pearson(x, y):
    x = x - x.mean(); y = y - y.mean()
    d = np.sqrt((x * x).sum() * (y * y).sum())
    return float((x * y).sum() / d) if d > 0 else float("nan")


def spearman(a, b):
    return pearson(rankdata(a), rankdata(b))


def partial_spearman(a, b, controls):
    ra, rb = rankdata(a), rankdata(b)
    C = np.column_stack([rankdata(c) for c in controls])
    return pearson(resid(ra, C), resid(rb, C))


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: partial_corr.py binwise_<A>_vs_<B>.tsv")
    path = sys.argv[1]
    rows = [l.rstrip("\n").split("\t") for l in open(path)]
    hdr = rows[0]
    idx = {h: i for i, h in enumerate(hdr)}
    for need in ("bes", "delta_contact", "marginal", "width"):
        if need not in idx:
            sys.exit("[ERROR] column '%s' missing -- rerun benchmark with the "
                     "updated --dump-bins (adds width+marginal)" % need)
    data = np.array([[float(r[idx[c]]) for c in ("bes", "delta_contact",
                     "marginal", "width", "on_feature")] for r in rows[1:]])
    bes, delta, marg, width, onf = data.T
    ok = np.isfinite(bes) & np.isfinite(delta) & np.isfinite(marg)
    bes, delta, marg, width, onf = bes[ok], delta[ok], marg[ok], width[ok], onf[ok]
    n = len(bes)

    print("file              : %s" % path)
    print("bins              : %d  (%d on-feature, %d zero-depth)"
          % (n, int(onf.sum()), int((marg == 0).sum())))
    print("rho_raw           : %+.4f   Spearman(BES, delta)" % spearman(bes, delta))
    print("rho | depth       : %+.4f   partial, control marginal"
          % partial_spearman(bes, delta, [marg]))
    print("rho | depth,width : %+.4f   partial, control marginal+width"
          % partial_spearman(bes, delta, [marg, width]))

    hi = marg >= np.median(marg)
    print("rho (well-covered): %+.4f   top-50%% by depth (n=%d)"
          % (spearman(bes[hi], delta[hi]), int(hi.sum())))
    # also: does BES itself track depth? (the confound's premise)
    print("[premise] rho(BES, depth) = %+.4f   (strong negative => high-BES bins "
          "are low-depth, i.e. confound plausible)" % spearman(bes, marg))


if __name__ == "__main__":
    main()
