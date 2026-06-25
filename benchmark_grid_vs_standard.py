#!/usr/bin/env python3
"""
benchmark_grid_vs_standard.py -- is the adaptive super-bin Hi-C grid actually
better than standard binning for the BEARING integration, or just different?

For each (locus x condition-pair) it bins the SAME downsampled valid-pairs onto
three grids in ONE pass, with identical width-density normalization, so the only
variable is where the bin edges fall:

  adaptive  : the BEARING super-bin grid (variable width, edges on signal
              transitions, gap-aware)
  uniform   : fixed-width bins at the adaptive grid's MEAN width over the window
              (matched effective resolution; isolates the value of signal-
              anchored placement vs uniform placement)
  standard  : fixed 10 kb bins (the traditional workflow)

Two metrics per grid:
  rho_contact / p_emp : claim-1 co-localization. p95(BES) per bin vs the band-
        limited |density_A - density_B| per bin, with a circular-shift null
        (same definitions as bes_hic_correlation.py, made variable-bin-safe by
        using bp midpoints and width-density).
  feat_conc : fraction of total |delta_contact| magnitude falling in bins that
        overlap the external AgR features, divided by the length fraction those
        bins occupy. >1 means differential contact signal is concentrated in the
        features beyond what their length alone would give. Coverage-normalized
        so a grid is not rewarded merely for putting more bins on features.

HONEST CAVEATS (print and keep in the writeup): the adaptive grid is DERIVED
from the 1D signal and the AgR features track that signal, so part of any
adaptive advantage is coherent-unit denoising of both the 1D and 3D summaries,
which is a binning benefit, not a uniquely-3D one. Single loci, small n. The
benchmark is built to report "no better than uniform" if that is the result.

INPUT
  --superbins   superbins_5kb.bed (chrom start end superbin_id width n_bearing)
  --counts      pair_counts_by_lib.tsv (cond bam pair_count); sets the downsample
  --bam-dir     dir holding the bam basenames in --counts
  --comparisons TSV (header) columns: locus region condA condB bes
                where bes = per-bin BEARING differential table for that pair
                (chrom start end bearing_score ...)
  --annot       AgR_mm10.bed (external features for feat_conc)
  --out         output TSV

ASCII only. Requires pysam, numpy, scipy.
"""
import argparse
import gzip
import os
import sys
import zlib

import numpy as np
from scipy import stats


# --------------------------------------------------------------- grid helpers
def load_superbins_window(path, chrom, a, b):
    starts, ends = [], []
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.split()
            if f[0] != chrom:
                continue
            s, e = int(f[1]), int(f[2])
            if e <= a or s >= b:
                continue
            starts.append(s); ends.append(e)
    o = np.argsort(starts)
    return np.array(starts)[o], np.array(ends)[o]


def uniform_grid(a, b, width):
    edges = np.arange(a, b, width, dtype=np.int64)
    edges = np.append(edges, b)
    return edges[:-1], edges[1:]


def assign(pos, starts, ends):
    i = int(np.searchsorted(starts, pos, side="right")) - 1
    if i < 0 or pos >= ends[i]:
        return -1
    return i


# ------------------------------------------------------------- pair binning
def keep_pair(name, keep_p, seedstr):
    return (zlib.crc32((seedstr + name).encode()) / 4294967296.0) < keep_p


def bin_condition(bam_paths, chrom, a, b, keep_p, grids, seed):
    """One pass over each BAM; assign within-window cis pairs to every grid.
    grids: list of (starts, ends). Returns list of counts dicts (one per grid)."""
    import pysam
    seedstr = "%d:" % seed
    counts = [dict() for _ in grids]
    for path in bam_paths:
        bam = pysam.AlignmentFile(path, "rb")
        so = bam.header.to_dict().get("HD", {}).get("SO", "")
        use_fetch = (so == "coordinate") and bam.has_index()
        itr = bam.fetch(chrom, a, b) if use_fetch else bam.fetch(until_eof=True)
        sys.stderr.write("    %s : %s\n" % (os.path.basename(path),
                                            "fetch" if use_fetch else "scan"))
        for r in itr:
            if (r.flag & 0xD00) or not (r.flag & 0x40):
                continue
            if not use_fetch and r.reference_name != chrom:
                continue
            if r.next_reference_id < 0 or r.next_reference_name != chrom:
                continue
            p1, p2 = r.reference_start, r.next_reference_start
            if p1 < a or p1 >= b or p2 < a or p2 >= b:
                continue
            if not keep_pair(r.query_name, keep_p, seedstr):
                continue
            for gi, (gs, ge) in enumerate(grids):
                i = assign(p1, gs, ge)
                if i < 0:
                    continue
                j = assign(p2, gs, ge)
                if j < 0:
                    continue
                if i > j:
                    i, j = j, i
                counts[gi][(i, j)] = counts[gi].get((i, j), 0) + 1
        bam.close()
    return counts


def dense(counts, n):
    M = np.zeros((n, n), dtype=np.float64)
    for (i, j), c in counts.items():
        M[i, j] += c
        if i != j:
            M[j, i] += c
    return M


# --------------------------------------------------------------- metrics
def delta_contact(matA, matB, mids, w, min_d, max_d):
    densA = matA / np.outer(w, w)
    densB = matB / np.outer(w, w)
    diff = np.abs(densA - densB)
    D = np.abs(mids[:, None] - mids[None, :])
    band = (D >= min_d) & (D <= max_d)
    u = np.triu(np.where(band, diff, 0.0), k=1)
    return u.sum(axis=1) + u.sum(axis=0)


def read_bes(path, chrom, a, b, col="bearing_score"):
    op = gzip.open if path.endswith(".gz") else open
    mids, vals, hdr = [], [], None
    with op(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if hdr is None:
                hdr = [c.strip() for c in f]
                ci = hdr.index(col) if col in hdr else 3
                continue
            if f[0] != chrom:
                continue
            s, e = int(f[1]), int(f[2])
            if e <= a or s >= b:
                continue
            try:
                v = float(f[ci])
            except (ValueError, IndexError):
                continue
            mids.append((s + e) // 2); vals.append(abs(v))
    return np.array(mids), np.array(vals)


def bes_per_bin(bes_mids, bes_vals, starts, ends):
    n = len(starts)
    out = np.full(n, np.nan)
    if len(bes_mids) == 0:
        return out
    idx = np.searchsorted(starts, bes_mids, side="right") - 1
    for b in range(n):
        sel = (idx == b) & (bes_mids < ends[b])
        if sel.any():
            out[b] = np.percentile(bes_vals[sel], 95)
    return out


def shift_null(x, y, n_perm, min_shift, seed=42):
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    n = len(x)
    if n < 20:
        return None, None, n
    obs, _ = stats.spearmanr(x, y)
    rng = np.random.default_rng(seed)
    lo = max(min_shift, 1); hi = max(n - min_shift, lo + 1)
    null = np.empty(n_perm)
    for k in range(n_perm):
        null[k], _ = stats.spearmanr(np.roll(x, rng.integers(lo, hi)), y)
    p_emp = max(float((np.abs(null) >= abs(obs)).mean()), 1.0 / n_perm)
    return float(obs), p_emp, n


def read_features(path, chrom, a, b):
    op = gzip.open if path.endswith(".gz") else open
    feats = []
    with op(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t") if "\t" in line else line.split()
            if f[0] != chrom:
                continue
            s, e = int(f[1]), int(f[2])
            if e > a and s < b:
                feats.append((s, e))
    return feats


def feat_concentration(delta, starts, ends, feats):
    w = (ends - starts).astype(float)
    in_feat = np.zeros(len(starts), dtype=bool)
    for s, e in feats:
        in_feat |= (ends > s) & (starts < e)
    sig = np.nan_to_num(delta)
    tot = sig.sum()
    if tot <= 0 or not in_feat.any():
        return float("nan")
    frac_sig = sig[in_feat].sum() / tot
    frac_len = w[in_feat].sum() / w.sum()
    return float(frac_sig / frac_len) if frac_len > 0 else float("nan")


# --------------------------------------------------------------- driver
def read_comparisons(path):
    rows, hdr = [], None
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = [x.strip() for x in line.rstrip("\n").split("\t")]
            if hdr is None:
                hdr = [c.lower() for c in f]
                continue
            rows.append(dict(zip(hdr, f)))
    return rows


def pooled_floor(counts_path):
    pooled, libs = {}, {}
    with open(counts_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.split("\t") if "\t" in line else line.split()
            if len(f) < 3:
                continue
            cond, bam, cnt = f[0], f[1], int(f[2])
            pooled[cond] = pooled.get(cond, 0) + cnt
            libs.setdefault(cond, []).append(bam)
    return pooled, libs, min(pooled.values())


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--superbins", required=True)
    ap.add_argument("--counts", required=True)
    ap.add_argument("--bam-dir", default=".")
    ap.add_argument("--comparisons", required=True)
    ap.add_argument("--annot", required=True)
    ap.add_argument("--bes-score-col", default="bearing_score")
    ap.add_argument("--std-binsize", type=int, default=10000)
    ap.add_argument("--min-distance", type=int, default=50000)
    ap.add_argument("--max-distance", type=int, default=500000)
    ap.add_argument("--n-perm", type=int, default=1000)
    ap.add_argument("--min-shift", type=int, default=100000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", required=True)
    ap.add_argument("--dump-bins", default=None,
                    help="prefix; write per-bin BES vs delta_contact (adaptive grid, "
                         "sorted by BES) to <prefix>_<A>_vs_<B>.tsv for eyeballing")
    args = ap.parse_args()

    pooled, libs, floor = pooled_floor(args.counts)
    sys.stderr.write("[floor] %d pooled pairs\n" % floor)

    def bam_paths(cond):
        return [b if os.path.isabs(b) else os.path.join(args.bam_dir, b)
                for b in libs[cond]]

    memo = {}  # (cond, region) -> (grids, counts_list, mids_list, w_list)

    rows_out = []
    for rec in read_comparisons(args.comparisons):
        locus, region = rec["locus"], rec["region"]
        A, B = rec["conda"], rec["condb"]
        bes = rec["bes"]
        chrom, rng = region.split(":")
        a, b = (int(x.replace(",", "")) for x in rng.split("-"))

        # build the three grids for this window
        ad_s, ad_e = load_superbins_window(args.superbins, chrom, a, b)
        if len(ad_s) == 0:
            sys.stderr.write("[skip] no super-bins in %s\n" % region); continue
        mean_w = int(round((ad_e - ad_s).mean()))
        un_s, un_e = uniform_grid(a, b, mean_w)
        st_s, st_e = uniform_grid(a, b, args.std_binsize)
        grids = [("adaptive", ad_s, ad_e), ("uniform", un_s, un_e),
                 ("standard", st_s, st_e)]
        grid_arrays = [(s, e) for _, s, e in grids]

        feats = read_features(args.annot, chrom, a, b)
        bes_mids, bes_vals = read_bes(bes, chrom, a, b, args.bes_score_col)

        cond_counts = {}
        for cond in (A, B):
            key = (cond, region)
            if key not in memo:
                kp = floor / pooled[cond]
                sys.stderr.write("[bin] %s @ %s (keep=%.3f)\n" % (cond, locus, kp))
                memo[key] = bin_condition(bam_paths(cond), chrom, a, b, kp,
                                          grid_arrays, args.seed)
            cond_counts[cond] = memo[key]

        for gi, (gname, gs, ge) in enumerate(grids):
            n = len(gs)
            w = (ge - gs).astype(float)
            mids = (gs + ge) / 2.0
            MA = dense(cond_counts[A][gi], n)
            MB = dense(cond_counts[B][gi], n)
            delta = delta_contact(MA, MB, mids, w, args.min_distance, args.max_distance)
            besb = bes_per_bin(bes_mids, bes_vals, gs, ge)
            rho, p_emp, nused = shift_null(
                besb, delta, args.n_perm,
                max(1, args.min_shift // max(int(w.mean()), 1)), args.seed)
            conc = feat_concentration(delta, gs, ge, feats)
            if args.dump_bins and gname == "adaptive":
                fa = np.array([any(s < fe and e > fs for fs, fe in feats)
                               for s, e in zip(gs, ge)])
                order = np.argsort(-besb)
                dp = "%s_%s_vs_%s.tsv" % (args.dump_bins, A, B)
                with open(dp, "w") as dh:
                    dh.write("chrom\tstart\tend\tbes\tdelta_contact\ton_feature\n")
                    for i in order:
                        dh.write("%s\t%d\t%d\t%.4g\t%.4g\t%d\n"
                                 % (chrom, gs[i], ge[i], besb[i], delta[i], int(fa[i])))
                sys.stderr.write("[dump] %s (%d bins, sorted by BES desc)\n"
                                 % (dp, len(gs)))
            marg = MA.sum(axis=1) + MB.sum(axis=1)
            rows_out.append([
                "%s_vs_%s" % (A, B), locus, gname, n, int(w.mean()),
                int((marg == 0).sum()),
                "NA" if rho is None else "%.4f" % rho,
                "NA" if p_emp is None else "%.4g" % p_emp,
                "%.3f" % conc if conc == conc else "NA",
                nused,
            ])
            sys.stderr.write("  [%s] %s n=%d meanw=%d rho=%s p_emp=%s conc=%s\n" % (
                gname, locus, n, int(w.mean()),
                "NA" if rho is None else "%.3f" % rho,
                "NA" if p_emp is None else "%.3g" % p_emp,
                "%.2f" % conc if conc == conc else "NA"))

    hdr = ["comparison", "locus", "grid", "n_bins", "mean_width_bp",
           "zero_marginal_bins", "rho_contact", "p_emp_contact",
           "feat_concentration", "n_used"]
    with open(args.out, "w") as fh:
        fh.write("\t".join(hdr) + "\n")
        for r in rows_out:
            fh.write("\t".join(str(x) for x in r) + "\n")
    sys.stderr.write("[out] %s (%d rows)\n" % (args.out, len(rows_out)))


if __name__ == "__main__":
    main()
