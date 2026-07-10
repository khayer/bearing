#!/usr/bin/env python3
# minsig_depth_check.py
#
# Question: is the between-sample spread in the min_signal=0.1 masked fraction
# (31.4% - 51.2%, with Pro-B replicates 17.5 points apart) driven by library
# depth / coverage breadth rather than by cell state?
#
# min_signal is an ABSOLUTE cut on summed raw signal and the production run has
# normalization OFF, so a sample with less raw signal should clear the 0.1 floor
# in fewer bins. NOTE: the RNAseq+/- bigWigs are CPM-normalized (*_CPM.bw) while
# ATAC / CTCF / Cohesin / H3K27ac are not, so if this mechanism is real the
# correlation should be carried by the four un-normalized tracks.
#
# Reads only bigWig HEADERS via bigWigInfo (fast; no data scan).
#   basesCovered = non-zero genome coverage; mean = mean over covered bases
#   total signal ~ mean * basesCovered
#
# Run from the repo root:
#   python minsig_depth_check.py --sheet workflow/config/samples.tsv
#
# ASCII only. Reads real data; fabricates nothing.

import argparse, re, subprocess, sys

# active (non-zero) bins at min_signal=0.1, measured from workflow/results/*.qcat.bgz
ACTIVE_AT_0P1 = {
    "DN_rep1":   6666644,
    "DN_rep2":   7180525,
    "DP_rep2":   8730566,
    "ProB_rep1": 9360604,
    "ProB_rep2": 6975463,
    "S3T3_rep1": 8071861,
    "S3T3_rep2": 7520016,
    # fill in once measured:
    # "DP_rep1": ..., "EbKO_rep1": ..., "EbKO_rep2": ...,
}
TRACKS = ["ATAC", "RNAseq+", "RNAseq-", "CTCF", "Cohesin", "H3K27ac"]
# indices of tracks that are NOT depth-normalized (RNAseq +/- are CPM)
RAW_IDX = [0, 3, 4, 5]


def bigwig_info(path):
    """Return (basesCovered, mean) from the bigWig header, or (None, None)."""
    try:
        out = subprocess.run(["bigWigInfo", path], capture_output=True,
                             text=True, timeout=120)
    except Exception as e:
        print("  bigWigInfo failed on %s: %s" % (path, e), file=sys.stderr)
        return None, None
    if out.returncode != 0:
        print("  bigWigInfo rc=%d on %s" % (out.returncode, path), file=sys.stderr)
        return None, None
    bc = mn = None
    for line in out.stdout.splitlines():
        m = re.match(r"\s*basesCovered:\s*([\d,]+)", line)
        if m:
            bc = int(m.group(1).replace(",", ""))
        m = re.match(r"\s*mean:\s*([-\d.eE+]+)", line)
        if m:
            mn = float(m.group(1))
    return bc, mn


def pearson(xs, ys):
    n = len(xs)
    if n < 3:
        return float("nan")
    mx = sum(xs) / n
    my = sum(ys) / n
    num = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    dx = sum((a - mx) ** 2 for a in xs) ** 0.5
    dy = sum((b - my) ** 2 for b in ys) ** 0.5
    return num / (dx * dy) if dx and dy else float("nan")


def spearman(xs, ys):
    def ranks(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        for pos, i in enumerate(order):
            r[i] = pos + 1.0
        return r
    return pearson(ranks(xs), ranks(ys))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sheet", default="workflow/config/samples.tsv")
    a = ap.parse_args()

    rows = []
    with open(a.sheet) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5 or f[0].strip() in ("", "sample"):
                continue
            sample, bw = f[0].strip(), f[4].strip()
            if not bw:
                continue
            rows.append((sample, [p.strip() for p in bw.split(",")]))

    print("%-11s %13s %13s %13s  %s"
          % ("sample", "total_signal", "raw_signal", "raw_bases", "per-track basesCovered"))
    table = []
    for sample, paths in rows:
        tot = raw_sig = raw_bases = 0.0
        percov = []
        ok = True
        for i, p in enumerate(paths[:6]):
            bc, mn = bigwig_info(p)
            if bc is None or mn is None:
                ok = False
                percov.append("NA")
                continue
            sig = bc * mn
            tot += sig
            percov.append("%.0fM" % (bc / 1e6))
            if i in RAW_IDX:
                raw_sig += sig
                raw_bases += bc
        if not ok:
            print("%-11s  (missing bigWig header info -- skipped)" % sample)
            continue
        print("%-11s %13.3e %13.3e %13.0f  %s"
              % (sample, tot, raw_sig, raw_bases, " ".join(percov)))
        table.append((sample, tot, raw_sig, raw_bases))

    # correlate against active bins at 0.1
    have = [t for t in table if t[0] in ACTIVE_AT_0P1]
    if len(have) < 3:
        print("\nNeed >=3 samples with known active-bin counts to correlate.")
        print("Fill ACTIVE_AT_0P1 with the remaining samples and rerun.")
        return
    act = [ACTIVE_AT_0P1[t[0]] for t in have]
    print("\nCorrelation vs active (non-zero) bins at min_signal=0.1, n=%d samples:" % len(have))
    for label, idx in (("total signal (all 6 tracks)", 1),
                       ("raw signal (4 un-normalized)", 2),
                       ("raw basesCovered (4 un-norm)", 3)):
        v = [t[idx] for t in have]
        print("  %-30s Pearson r=%+.3f   Spearman rho=%+.3f"
              % (label, pearson(v, act), spearman(v, act)))
    print("\nInterpretation:")
    print("  strong positive r for the RAW tracks (and weaker for total) supports")
    print("  'absolute threshold on un-normalized signal' as the mechanism.")
    print("  near-zero r means the spread is NOT a simple depth artifact -- in that")
    print("  case report only the observed ranges and make no causal claim.")


if __name__ == "__main__":
    main()
