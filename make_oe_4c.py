#!/usr/bin/env python3
"""
make_oe_4c.py

Turn a Hi-C / capture-Hi-C .cool (or .mcool) into per-viewpoint,
observed/expected "virtual 4C" bigwig tracks that drop straight into the
bigwig_to_qcat panel.

Why O/E: a raw 4C profile is dominated by the near-viewpoint (proximity)
signal, which is high in every genotype and does not change. Dividing each
bin by the expected contact at its genomic distance flattens the decay curve
so a real loop (e.g. Trbv1 <-> RC) becomes a per-bin enrichment peak that sits
on a ~0 background -- a sparse, peaky track that BEARING handles well
(especially with normalize = nonzero-quantile, "zeros stay zero").

Output value modes:
  enrichment (default): max(O/E - 1, 0)  -> 0 background, positive at loops
  oe                  : O/E              -> ~1 background
  log2                : log2(O/E)        -> signed (note: qcat clips negatives)

Usage (single cool):
  python make_oe_4c.py --cool sample.cool --viewpoints vp.bed --outdir out/

Usage (iterate a BEARING samples sheet, using its `cool` column):
  python make_oe_4c.py --samples workflow/config/samples.tsv \
      --data-dir /mnt/.../bearing --viewpoints viewpoints.bed \
      --resolution 10000 --outdir ../bigwigs/OE4C

viewpoints BED (tab-separated, no header): chrom  start  end  name
  chr6  41530000  41560000  RC
  chr6  40890000  40895000  Trbv1
"""

import argparse
import os
import sys
import numpy as np

try:
    import cooler
except ImportError:
    sys.exit("ERROR: need `cooler` (pip install cooler)")
try:
    import pyBigWig
except ImportError:
    sys.exit("ERROR: need `pyBigWig` (pip install pyBigWig)")


def log(msg):
    sys.stderr.write(msg + "\n")


def read_viewpoints(path):
    vps = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 4:
                f = line.split()
            if len(f) < 4:
                raise ValueError("viewpoints BED needs 4 cols: chrom start end name")
            vps.append((f[0], int(f[1]), int(f[2]), f[3]))
    return vps


def read_samples_cools(sheet, data_dir):
    """Return list of (sample_name, cool_path) from a BEARING samples.tsv."""
    rows = []
    header = None
    with open(sheet) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            f = line.split("\t")
            if header is None:
                header = f
                continue
            rec = dict(zip(header, f))
            cool = rec.get("cool", "").strip()
            if not cool:
                continue
            name = rec.get("sample", os.path.basename(cool)).strip()
            if not os.path.isabs(cool):
                cool = os.path.join(data_dir, cool)
            rows.append((name, cool))
    return rows


def open_cool(path, resolution):
    if "::" in path:
        return cooler.Cooler(path)
    if path.endswith(".mcool"):
        if resolution is None:
            raise ValueError("%s is multi-res; pass --resolution" % path)
        return cooler.Cooler("%s::/resolutions/%d" % (path, resolution))
    clr = cooler.Cooler(path)
    if resolution is not None and clr.binsize != resolution:
        log("  WARNING: %s binsize=%s but --resolution=%s; using file binsize"
            % (path, clr.binsize, resolution))
    return clr


def expected_by_diagonal(mat):
    """Mean over each diagonal (ignoring NaNs); index = bin offset."""
    n = mat.shape[0]
    exp = np.full(n, np.nan, dtype=float)
    for k in range(n):
        d = np.diagonal(mat, k)
        if np.any(np.isfinite(d)):
            exp[k] = np.nanmean(d)
    return exp


def oe_4c_for_viewpoint(clr, chrom, vstart, vend, max_distance,
                        self_mask, balance):
    """Return (starts, values_OE) for the windowed region around a viewpoint."""
    binsize = clr.binsize
    chromsize = int(clr.chromsizes[chrom])
    vcenter = (vstart + vend) // 2
    wstart = max(0, vcenter - max_distance)
    wend = min(chromsize, vcenter + max_distance)
    # snap to bin edges
    wstart -= wstart % binsize
    if wend % binsize:
        wend += binsize - (wend % binsize)
    wend = min(wend, chromsize)
    region = "%s:%d-%d" % (chrom, wstart, wend)

    try:
        mat = clr.matrix(balance=balance).fetch(region)
    except Exception as e:
        if balance:
            log("  balance failed (%s); falling back to raw counts" % e)
            mat = clr.matrix(balance=False).fetch(region)
        else:
            raise
    mat = np.asarray(mat, dtype=float)

    binsdf = clr.bins().fetch(region)
    starts = binsdf["start"].values.astype(int)
    n = len(starts)
    if mat.shape[0] != n:
        n = min(n, mat.shape[0])
        mat = mat[:n, :n]
        starts = starts[:n]

    exp = expected_by_diagonal(mat)
    ends = starts + binsize
    vp_idx = np.where((starts < vend) & (ends > vstart))[0]
    if len(vp_idx) == 0:
        return starts, np.zeros(n)

    cols = np.arange(n)
    accum = np.zeros(n)
    cnt = np.zeros(n)
    for v in vp_idx:
        offs = np.abs(cols - v)
        e = exp[offs]
        obs = mat[v, :]
        with np.errstate(divide="ignore", invalid="ignore"):
            row = obs / e
        ok = np.isfinite(row)
        accum[ok] += row[ok]
        cnt[ok] += 1.0
    oe = np.where(cnt > 0, accum / np.maximum(cnt, 1.0), np.nan)

    # mask self-ligation / near-viewpoint diagonal
    bin_centers = starts + binsize / 2.0
    oe[np.abs(bin_centers - vcenter) < self_mask] = np.nan
    return starts, oe


def transform(oe, mode):
    v = np.array(oe, dtype=float)
    if mode == "enrichment":
        out = np.clip(v - 1.0, 0.0, None)
    elif mode == "oe":
        out = v
    elif mode == "log2":
        with np.errstate(divide="ignore", invalid="ignore"):
            out = np.log2(v)
    else:
        raise ValueError("mode must be enrichment|oe|log2")
    out[~np.isfinite(out)] = 0.0
    return out


def write_bigwig(path, chromsizes, chrom, starts, values, binsize):
    bw = pyBigWig.open(path, "w")
    header = sorted(((c, int(s)) for c, s in chromsizes.items()))
    bw.addHeader(header)
    starts = np.asarray(starts, dtype=int)
    values = np.asarray(values, dtype=float)
    ends = np.minimum(starts + binsize, int(chromsizes[chrom]))
    keep = ends > starts
    bw.addEntries([chrom] * int(keep.sum()),
                  starts[keep].tolist(),
                  ends=ends[keep].tolist(),
                  values=values[keep].astype(float).tolist())
    bw.close()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--cool", help="single .cool/.mcool file")
    src.add_argument("--samples", help="BEARING samples.tsv (uses its `cool` column)")
    ap.add_argument("--viewpoints", required=True, help="BED: chrom start end name")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--data-dir", default=".", help="resolves relative cool paths in --samples")
    ap.add_argument("--resolution", type=int, default=None, help="required for .mcool")
    ap.add_argument("--max-distance", type=int, default=6000000,
                    help="half-window around viewpoint (bp); default 6e6")
    ap.add_argument("--self-mask", type=int, default=20000,
                    help="mask +/- this many bp around the viewpoint (bp); default 2e4")
    ap.add_argument("--mode", default="enrichment", choices=["enrichment", "oe", "log2"])
    ap.add_argument("--no-balance", action="store_true", help="use raw counts, not balanced")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    vps = read_viewpoints(args.viewpoints)
    log("viewpoints: %s" % ", ".join(v[3] for v in vps))

    if args.cool:
        base = os.path.basename(args.cool).split("::")[0]
        for ext in (".mcool", ".cool"):
            if base.endswith(ext):
                base = base[: -len(ext)]
        cools = [(base, args.cool)]
    else:
        cools = read_samples_cools(args.samples, args.data_dir)
    if not cools:
        sys.exit("No cool files found (empty `cool` column?)")

    manifest = []
    for sample, coolpath in cools:
        log("[%s] %s" % (sample, coolpath))
        clr = open_cool(coolpath, args.resolution)
        chromsizes = dict(clr.chromsizes.items())
        for chrom, vstart, vend, name in vps:
            if chrom not in chromsizes:
                log("  skip %s: %s not in cooler" % (name, chrom))
                continue
            starts, oe = oe_4c_for_viewpoint(
                clr, chrom, vstart, vend,
                args.max_distance, args.self_mask,
                balance=not args.no_balance)
            vals = transform(oe, args.mode)
            out = os.path.join(args.outdir, "%s_%s_OE4C.bw" % (sample, name))
            write_bigwig(out, chromsizes, chrom, starts, vals, clr.binsize)
            log("  wrote %s  (%d bins, max=%.3f)" % (out, len(starts), float(np.max(vals)) if len(vals) else 0.0))
            manifest.append(out)

    log("\nDone. %d bigwigs:" % len(manifest))
    for m in manifest:
        print(m)


if __name__ == "__main__":
    main()
