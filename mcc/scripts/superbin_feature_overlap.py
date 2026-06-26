#!/usr/bin/env python3
"""
superbin_feature_overlap.py -- is the adaptive super-bin grid FINE AT FEATURES
(signal-anchored, the proposal's claim) or just fine wherever the genome is
mappable (blacklist fragmenting a repetitive locus into islands)? Both give the
same width histogram; this separates them, using only geometry (no contacts, no
BES required).

Within a window it splits super-bins into those overlapping a feature (optionally
padded) and those not, and compares widths. If the grid is feature-anchored,
feature bins are NARROWER. Significance comes from a null that circularly shifts
the feature set within the window many times and recomputes the gap, so a grid
that is uniformly fine everywhere (or fine in the wrong places) returns p ~ 1.

  width_ratio = median(width | on feature) / median(width | off feature)
                < 1  => resolution concentrated on features (feature-anchored)

INPUT
  --superbins  superbins_*.bed   (chrom start end ...)
  --features   bed of features   (chrom start end [name]); cat AgR + CBE if wanted
  --region     chrom:start-end   (the capture window)
  --pad        bp to extend each feature on both sides (default 0)
  --n-perm / --seed

ASCII only. numpy only.
"""
import argparse
import sys
import numpy as np


def load_bed_region(path, chrom, a, b):
    iv = []
    op = open
    if path.endswith(".gz"):
        import gzip
        op = gzip.open
    with op(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.split()
            if f[0] != chrom:
                continue
            s, e = int(f[1]), int(f[2])
            if e <= a or s >= b:
                continue
            iv.append((max(s, a), min(e, b)))
    return iv


def overlaps_any(s, e, feats):
    for fs, fe in feats:
        if s < fe and e > fs:
            return True
    return False


def median_gap(starts, ends, feats):
    widths = ends - starts
    on = np.array([overlaps_any(s, e, feats) for s, e in zip(starts, ends)])
    if on.sum() == 0 or (~on).sum() == 0:
        return None, None, None, on
    wf = float(np.median(widths[on]))
    wn = float(np.median(widths[~on]))
    return wf, wn, wf - wn, on


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--superbins", required=True)
    ap.add_argument("--features", required=True)
    ap.add_argument("--region", required=True)
    ap.add_argument("--pad", type=int, default=0)
    ap.add_argument("--n-perm", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    chrom, rng = args.region.split(":")
    a, b = (int(x.replace(",", "")) for x in rng.split("-"))
    W = b - a

    sb = load_bed_region(args.superbins, chrom, a, b)
    if not sb:
        sys.exit("[ERROR] no super-bins in region")
    starts = np.array([s for s, _ in sb])
    ends = np.array([e for _, e in sb])

    feats = load_bed_region(args.features, chrom, a, b)
    if args.pad:
        feats = [(max(a, fs - args.pad), min(b, fe + args.pad)) for fs, fe in feats]
    if not feats:
        sys.exit("[ERROR] no features in region")

    wf, wn, obs, on = median_gap(starts, ends, feats)
    if wf is None:
        sys.exit("[ERROR] all bins on one side of the feature split")

    # null: circularly shift the feature set within the window
    rng_ = np.random.default_rng(args.seed)
    flo = np.array([fs - a for fs, _ in feats])
    fhi = np.array([fe - a for _, fe in feats])
    n_le = 0
    for _ in range(args.n_perm):
        sh = int(rng_.integers(0, W))
        ps = (flo + sh) % W + a
        pe = ps + (fhi - flo)
        # split wrapped intervals at the window edge
        perm = []
        for s0, e0 in zip(ps, pe):
            if e0 <= b:
                perm.append((s0, e0))
            else:
                perm.append((s0, b)); perm.append((a, a + (e0 - b)))
        _, _, d, _ = median_gap(starts, ends, perm)
        if d is not None and d <= obs:   # as or more feature-concentrated
            n_le += 1
    p_emp = (n_le + 1) / (args.n_perm + 1)

    feat_bp = sum(e - s for s, e in feats)
    print("region            : %s  (%.2f Mb)" % (args.region, W / 1e6))
    print("super-bins         : %d  (%d on features, %d off)"
          % (len(sb), int(on.sum()), int((~on).sum())))
    print("feature coverage   : %.1f%% of window (pad=%d)" % (100 * feat_bp / W, args.pad))
    print("median width ON  features: %d bp" % wf)
    print("median width OFF features: %d bp" % wn)
    print("width ratio (on/off)     : %.3f   (<1 => finer at features)" % (wf / wn))
    print("feature-shift null p_emp : %.4f   (small => anchoring is not by chance)" % p_emp)
    verdict = ("FEATURE-ANCHORED: resolution concentrates on features"
               if (wf < wn and p_emp < 0.05)
               else "NOT feature-anchored: fragmentation looks mappability-driven, "
                    "not signal-driven")
    print("verdict            :", verdict)


if __name__ == "__main__":
    main()
