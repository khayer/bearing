#!/usr/bin/env python3
"""
vbeta_factor_quant.py -- is an architectural factor (Cohesin/CTCF/NIPBL) reduced
ACROSS THE Vb CLUSTER in a mutant vs DN, not just at the perturbed element?

Reads the per-track bigwigs for two conditions from the V1 sheet, computes mean
signal per Vb segment, and tests whether the mutant is systematically lower.

NORMALIZATION (critical): does NOT normalize to the window total -- that would
divide out a cluster-wide reduction. Normalizes each condition to a CONTROL
region the perturbation leaves untouched (for V1P, Allyn shows the 3'CBE/RC is
promoter-independent -> use --control-region there). Per-segment ratios are then
comparable across conditions, and a sign test asks if the cluster as a whole drops.

  python vbeta_factor_quant.py --sheet mcc/config/samples_v1.tsv \
    --cats mcc/config/v1_5track_panel.yaml \
    --conditions DN,V1P --tracks Cohesin,CTCF,NIPBL \
    --segments annotations/vbeta_segments.bed \
    --control-region chr6:41500000-41620000 \
    --out vbeta_cohesin_quant.tsv

Requires pyBigWig (cluster env). ASCII only. Stats are stdlib (sign test +
median log2 ratio); no scipy dependency.
"""
import argparse
import json
import math
import sys


def read_sheet(path):
    cond = {}
    with open(path) as fh:
        header = None
        for line in fh:
            s = line.rstrip("\n")
            if not s.strip() or s.startswith("#"):
                continue
            header = s.split("\t")
            break
        idx = {h: i for i, h in enumerate(header)}
        ci, bi = idx.get("condition", 1), idx.get("bw", 4)
        for line in fh:
            s = line.rstrip("\n")
            if not s.strip() or s.startswith("#"):
                continue
            f = s.split("\t")
            if len(f) <= bi:
                continue
            cond.setdefault(f[ci], []).append(
                [p.strip() for p in f[bi].split(",") if p.strip()])
    return cond  # {cond: [ [bw per track] per replicate ]}


def read_cats(path):
    if path.endswith(".json"):
        d = json.load(open(path))["categories"]
        return [d[k][0].strip() for k in sorted(d, key=lambda x: int(x))]
    names = []
    for line in open(path):
        t = line.strip()
        if t.startswith("- name:"):
            names.append(t.split(":", 1)[1].strip().strip('"'))
    return names


def read_bed(path):
    out = []
    for line in open(path):
        if not line.strip() or line.startswith(("#", "track")):
            continue
        f = line.split()
        name = f[3] if len(f) > 3 else "%s:%s-%s" % (f[0], f[1], f[2])
        out.append((f[0], int(f[1]), int(f[2]), name))
    return out


def sign_test(deltas):
    """two-sided sign test: are values systematically != 0?"""
    pos = sum(1 for d in deltas if d > 0)
    neg = sum(1 for d in deltas if d < 0)
    n = pos + neg
    if n == 0:
        return float("nan"), 0, 0
    k = min(pos, neg)
    # exact two-sided binomial tail at p=0.5
    from math import comb
    p = 2.0 * sum(comb(n, i) for i in range(0, k + 1)) / (2 ** n)
    return min(p, 1.0), pos, neg


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sheet", required=True)
    ap.add_argument("--cats", required=True)
    ap.add_argument("--conditions", required=True, help="REF,TEST e.g. DN,V1P")
    ap.add_argument("--tracks", required=True, help="comma list, e.g. Cohesin,CTCF,NIPBL")
    ap.add_argument("--segments", required=True, help="BED of Vb segments")
    ap.add_argument("--control-region", required=True, help="chrom:start-end, "
                    "perturbation-independent (V1P: 3'CBE/RC)")
    ap.add_argument("--pad", type=int, default=0)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    import pyBigWig
    ref, test = args.conditions.split(",")
    tracks = args.tracks.split(",")
    cats = read_cats(args.cats)
    sheet = read_sheet(args.sheet)
    segs = read_bed(args.segments)
    cc, crng = args.control_region.split(":")
    ca, cb = (int(x) for x in crng.split("-"))

    def track_idx(name):
        for i, c in enumerate(cats):
            if c.replace(" ", "").lower() == name.replace(" ", "").lower():
                return i
        sys.exit("[ERROR] track %s not in cats (%s)" % (name, cats))

    def mean_signal(cond, ti, chrom, a, b):
        """mean over replicates of the per-bp mean signal in [a,b]."""
        vals = []
        for rep_bws in sheet[cond]:
            bw = pyBigWig.open(rep_bws[ti])
            v = bw.stats(chrom, a, b, type="mean")[0]
            bw.close()
            if v is not None:
                vals.append(v)
        return sum(vals) / len(vals) if vals else float("nan")

    rows = []
    summary = []
    for tname in tracks:
        ti = track_idx(tname)
        # control normalization factor per condition
        ref_ctrl = mean_signal(ref, ti, cc, ca, cb)
        test_ctrl = mean_signal(test, ti, cc, ca, cb)
        deltas = []
        for chrom, s, e, name in segs:
            rv = mean_signal(ref, ti, chrom, max(0, s - args.pad), e + args.pad)
            tv = mean_signal(test, ti, chrom, max(0, s - args.pad), e + args.pad)
            if not (rv and tv and ref_ctrl and test_ctrl) or rv <= 0 or tv <= 0:
                continue
            # normalize each condition to its control, then ratio TEST/REF
            rn = rv / ref_ctrl
            tn = tv / test_ctrl
            l2 = math.log2(tn / rn)
            deltas.append(l2)
            rows.append([tname, name, "%.4g" % rv, "%.4g" % tv,
                         "%.4g" % ref_ctrl, "%.4g" % test_ctrl, "%.3f" % l2])
        if deltas:
            deltas_sorted = sorted(deltas)
            med = deltas_sorted[len(deltas_sorted) // 2]
            p, pos, neg = sign_test(deltas)
            summary.append((tname, len(deltas), med, neg, pos, p))

    with open(args.out, "w") as fh:
        fh.write("track\tsegment\tref_signal\ttest_signal\tref_ctrl\ttest_ctrl\t"
                 "log2_test_over_ref\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")

    print("control-normalized %s -> %s across %d Vb segments\n" % (ref, test, len(segs)))
    print("%-8s %5s %12s %14s %10s" % ("track", "n", "median_log2", "n_down/n_up", "sign_p"))
    for tname, n, med, neg, pos, p in summary:
        flag = "  <-- reduced" if (med < 0 and p < 0.05) else ""
        print("%-8s %5d %12.3f %8d/%-5d %10.2g%s"
              % (tname, n, med, neg, pos, p, flag))
    print("\n(median_log2 < 0 with small sign_p => factor systematically reduced "
          "across the Vb cluster, not just at the perturbed element)")


if __name__ == "__main__":
    main()
