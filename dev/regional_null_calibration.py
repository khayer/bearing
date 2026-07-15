#!/usr/bin/env python3
"""
regional_null_calibration.py
============================
Empirical calibration of the regional-enrichment test.

WHY
regional_enrichment.py combines a spatial-concentration binomial with a
directional-concordance binomial via Fisher's method, assuming bins are
exchangeable within a region. Bins are spatially autocorrelated (differential
score independence length ~967 bp on a 200 bp grid), so exchangeability is
violated. The manuscript ASSERTS the test is conservative under this violation.
This script measures it.

TWO REQUIREMENTS, BOTH EASY TO GET WRONG
1. The contrast must itself be NULL. p_spatial is conditional on the locus: it
   asks whether the locus's significant bins concentrate in the sub-region. Run
   it on DN-vs-DP at Tcrb and a random window overlapping the Vbeta cluster is
   CORRECTLY significant. Use a within-condition replicate differential --
   bearing_calibration.py writes results/calibration/<C>/<C>/<C>_perbin.tsv.gz
   with exactly the columns needed.
2. The random regions must be INDEPENDENT. Sampling windows inside one 900 kb
   locus gives, for a 1250-bin region in a 4500-bin locus, only ~3.6
   non-overlapping placements: "0 of 1000 windows" is really "0 of ~4".

MODES
  --mode within      random windows inside the analysis locus (fast, but the
                     effective sample size collapses for large regions)
  --mode background  one matched window per randomly drawn background locus of
                     the same span elsewhere in the genome, excluding the
                     antigen-receptor loci. Draws are independent. DEFAULT.

FAST PATH + EQUIVALENCE
compute_regional_enrichment() rescans all ~13.6M bins on every call (once for
the locus slice, once for g_dir), which is far too slow for hundreds of
background loci. This script reimplements the identical formulas on sorted
numpy arrays and, with --verify N, checks them against the production function
on N randomly chosen (locus, region) pairs. Any mismatch is a hard error.

READING THE OUTPUT
  FPR@alpha ~ alpha   -> calibrated
  FPR@alpha << alpha  -> conservative, as the manuscript claims
  FPR@alpha >> alpha  -> anti-conservative; regional q-values must be revised
Because the component p-values are discrete binomials with atoms (a large atom
at p = 1 whenever k = 0), the rate at a given alpha can be misleading. The
column that matters is `<=real_p`: the fraction of independent null regions
whose p_combined is at least as extreme as the real region's.

USAGE
  python regional_null_calibration.py \\
      --stats   workflow/results/calibration/DN/DN/DN_perbin.tsv.gz \\
      --regions workflow/annotations/tcrb_regions_v5.bed \\
      --locus   chr6:40790000-41690000 \\
      --chrom-sizes workflow/resources/mm10.chrom.sizes \\
      --mode background --n-random 500 --p-thresh 0.05 --seed 42 \\
      --verify 5 --null-contrast \\
      --out regional_null_calibration_tcrb_DNrep.tsv

ASCII only. Reads real data; fabricates nothing.
"""

import argparse
import gzip
import math
import os
import sys

import numpy as np
import scipy.stats

# Antigen-receptor loci and other regions with expected lymphoid signal.
# Background loci must not overlap these.
DEFAULT_EXCLUDE = [
    ("chr6", 40000000, 42500000),    # Tcrb
    ("chr6", 67000000, 71000000),    # Igk
    ("chr12", 112900000, 116200000),  # Igh
    ("chr14", 52000000, 55000000),   # Tcra/Tcrd
    ("chr16", 19000000, 19400000),   # Igl
]


def parse_ucsc(s):
    chrom, rest = s.split(":")
    a, b = rest.split("-")
    return chrom, int(a), int(b)


def load_stats(path, p_thresh):
    """chrom -> (starts[int64] sorted, pvals[float], signs[int8]). Also g_dir."""
    opener = gzip.open if path.endswith(".gz") else open
    per = {}
    n = 0
    with opener(path, "rt") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        ci, si = head.index("chrom"), head.index("start")
        pi, bi = head.index("pval"), head.index("bearing_score")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            per.setdefault(f[ci], [[], [], []])
            per[f[ci]][0].append(int(f[si]))
            per[f[ci]][1].append(float(f[pi]))
            per[f[ci]][2].append(1 if float(f[bi]) > 0 else -1)
            n += 1
    out = {}
    sig_pos = sig_tot = 0
    for c, (s, p, g) in per.items():
        s = np.asarray(s, dtype=np.int64)
        p = np.asarray(p, dtype=np.float64)
        g = np.asarray(g, dtype=np.int8)
        o = np.argsort(s, kind="stable")
        s, p, g = s[o], p[o], g[o]
        out[c] = (s, p, g)
        m = p < p_thresh
        sig_tot += int(m.sum())
        sig_pos += int((g[m] > 0).sum())
    g_dir = (sig_pos / sig_tot) if sig_tot else 0.5
    g_dir = min(max(g_dir, 1e-6), 1.0 - 1e-6)
    return out, n, g_dir, sig_tot


def slice_locus(idx, chrom, lo, hi):
    if chrom not in idx:
        return None
    s, p, g = idx[chrom]
    a = np.searchsorted(s, lo, "left")
    b = np.searchsorted(s, hi, "left")
    if b - a < 1:
        return None
    return s[a:b], p[a:b], g[a:b]


def regional_test(sl, i0, i1, p_thresh, g_dir):
    """EXACTLY regional_enrichment.compute_regional_enrichment's math, for the
    window covering scored-bin indices [i0, i1) of the locus slice `sl`."""
    s, p, g = sl
    L_locus_bins = s.size
    n_locus = int((p < p_thresh).sum())
    L_region_bins = i1 - i0

    rp, rg = p[i0:i1], g[i0:i1]
    m = rp < p_thresh
    k = int(m.sum())
    k_pos = int((rg[m] > 0).sum())
    k_neg = k - k_pos

    pi = (L_region_bins / L_locus_bins) if L_locus_bins > 0 else 0.0
    if n_locus == 0:
        p_spatial = 1.0
    else:
        p_spatial = min(float(scipy.stats.binom.sf(k - 1, n=n_locus, p=pi)), 1.0)

    if k == 0:
        p_directional = 1.0
    elif k_pos >= k_neg:
        p_directional = min(float(scipy.stats.binom.sf(k_pos - 1, n=k, p=g_dir)), 1.0)
    else:
        p_directional = min(float(scipy.stats.binom.sf(k_neg - 1, n=k, p=1.0 - g_dir)), 1.0)

    if p_spatial == 0 or p_directional == 0:
        p_combined = 0.0
    else:
        T = -2.0 * (math.log(p_spatial) + math.log(p_directional))
        p_combined = float(scipy.stats.chi2.sf(T, df=4))
    return p_spatial, p_directional, p_combined, k, n_locus


def overlaps(chrom, lo, hi, blocks):
    return any(c == chrom and not (hi <= a or lo >= b) for c, a, b in blocks)


def sample_background_loci(idx, sizes, span, n_want, ref_bins, rng, exclude, tol=0.25):
    """Random loci of the same span, matched on scored-bin count within +/- tol."""
    chroms = [c for c in sizes if c in idx and "_" not in c and c not in ("chrM", "chrX", "chrY")]
    if not chroms:
        return []
    w = np.array([sizes[c] for c in chroms], dtype=float)
    w /= w.sum()
    out, tries, max_tries = [], 0, n_want * 200
    while len(out) < n_want and tries < max_tries:
        tries += 1
        c = chroms[int(rng.choice(len(chroms), p=w))]
        if sizes[c] <= span:
            continue
        lo = int(rng.integers(0, sizes[c] - span))
        hi = lo + span
        if overlaps(c, lo, hi, exclude):
            continue
        sl = slice_locus(idx, c, lo, hi)
        if sl is None:
            continue
        nb = sl[0].size
        if abs(nb - ref_bins) > tol * ref_bins:
            continue
        out.append((c, lo, hi, sl))
    return out


def verify_against_production(repo, stats_path, regions_path, locus, p_thresh, idx, g_dir, n):
    """Hard check: our fast math must equal compute_regional_enrichment()."""
    sys.path.insert(0, repo)
    from regional_enrichment import (parse_diff_qcat, parse_diff_pvals,
                                     parse_regions_bed, compute_regional_enrichment)
    print("  verifying fast path against regional_enrichment.compute_regional_enrichment ...")
    if stats_path.endswith(".gz"):
        import shutil, tempfile
        tmp = os.path.join(tempfile.mkdtemp(prefix="regcal_"), "perbin.tsv")
        with gzip.open(stats_path, "rt") as s, open(tmp, "w") as d:
            shutil.copyfileobj(s, d)
        stats_path = tmp
    qcat = parse_diff_qcat(stats_path)
    pv = parse_diff_pvals(stats_path)
    real = parse_regions_bed(regions_path)
    lc, ls, le = locus
    prod = {r["region_name"]: r for r in
            compute_regional_enrichment(qcat, pv, real, lc, ls, le, p_thresh)}
    sl = slice_locus(idx, lc, ls, le)
    worst = 0.0
    for (c, s, e, name) in real[:n]:
        i0 = int(np.searchsorted(sl[0], s, "left"))
        i1 = int(np.searchsorted(sl[0], e, "left"))
        ps, pd, pc, k, nl = regional_test(sl, i0, i1, p_thresh, g_dir)
        pr = prod[name]
        d = max(abs(ps - pr["p_spatial"]), abs(pd - pr["p_directional"]),
                abs(pc - pr["p_combined"]))
        worst = max(worst, d)
        print("    %-28s k=%-4d prod p=%-10.4g fast p=%-10.4g  |diff|=%.2g"
              % (name[:28], k, pr["p_combined"], pc, d))
    if worst > 1e-9:
        sys.exit("FAST PATH DISAGREES WITH PRODUCTION (max |diff| = %.3g). Aborting." % worst)
    print("    OK: max |diff| = %.2g\n" % worst)
    return prod


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stats", required=True)
    ap.add_argument("--regions", required=True)
    ap.add_argument("--locus", required=True)
    ap.add_argument("--chrom-sizes", required=True)
    ap.add_argument("--mode", choices=["within", "background"], default="background")
    ap.add_argument("--p-thresh", type=float, default=0.05)
    ap.add_argument("--n-random", type=int, default=500)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--repo", default=".")
    ap.add_argument("--verify", type=int, default=5,
                    help="check the fast path against production on N real regions (0 to skip)")
    ap.add_argument("--null-contrast", action="store_true")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    if not a.null_contrast:
        print("!" * 72)
        print("WARNING: --null-contrast not set. If this contrast carries real signal")
        print("inside the locus, the rejection rates below are NOT type-I errors.")
        print("!" * 72)

    lc, ls, le = parse_ucsc(a.locus)
    sizes = {}
    with open(a.chrom_sizes) as fh:
        for line in fh:
            f = line.split()
            if len(f) >= 2:
                sizes[f[0]] = int(f[1])

    print("loading %s ..." % a.stats)
    idx, nbins, g_dir, nsig = load_stats(a.stats, a.p_thresh)
    print("  %d bins, %d below p<%.3g, genome-wide g_dir = %.4f" % (nbins, nsig, a.p_thresh, g_dir))

    prod = {}
    if a.verify:
        prod = verify_against_production(os.path.abspath(a.repo), a.stats, a.regions,
                                         (lc, ls, le), a.p_thresh, idx, g_dir, a.verify)

    real = []
    with open(a.regions) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4:
                real.append((f[0], int(f[1]), int(f[2]), f[3]))

    sl = slice_locus(idx, lc, ls, le)
    if sl is None:
        sys.exit("no scored bins in the analysis locus")
    ref_bins = sl[0].size
    print("  analysis locus holds %d scored bins" % ref_bins)

    rng = np.random.default_rng(a.seed)
    span = le - ls
    bg = []
    if a.mode == "background":
        bg = sample_background_loci(idx, sizes, span, a.n_random, ref_bins, rng, DEFAULT_EXCLUDE)
        print("  drew %d independent background loci (span %d bp, scored bins within 25%%)\n"
              % (len(bg), span))
        if len(bg) < 50:
            sys.exit("too few background loci matched; loosen --tol or reduce --n-random")

    rows = []
    print("%-28s %7s %7s %10s %10s %10s %10s"
          % ("region (size class)", "n_bins", "n_null", "FPR@0.05", "FPR@0.01", "real p", "<=real_p"))
    for (c, s, e, name) in real:
        if c != lc:
            continue
        i0 = int(np.searchsorted(sl[0], s, "left"))
        i1 = int(np.searchsorted(sl[0], e, "left"))
        nb = i1 - i0
        if nb < 5:
            continue
        rp = prod[name]["p_combined"] if name in prod else \
            regional_test(sl, i0, i1, a.p_thresh, g_dir)[2]

        ps_list = []
        if a.mode == "background":
            for (bc, blo, bhi, bsl) in bg:
                L = bsl[0].size
                if L <= nb:
                    continue
                j = int(rng.integers(0, L - nb))
                _, _, pc, _, _ = regional_test(bsl, j, j + nb, a.p_thresh, g_dir)
                ps_list.append(pc)
                rows.append((name, nb, bc, int(bsl[0][j]), pc))
        else:
            L = sl[0].size
            for _ in range(a.n_random):
                j = int(rng.integers(0, L - nb))
                _, _, pc, _, _ = regional_test(sl, j, j + nb, a.p_thresh, g_dir)
                ps_list.append(pc)
                rows.append((name, nb, lc, int(sl[0][j]), pc))

        p = np.asarray(ps_list)
        if p.size < 20:
            print("%-28s %7d   (only %d null draws; skipped)" % (name[:28], nb, p.size))
            continue
        print("%-28s %7d %7d %10.4f %10.4f %10.4g %10.4f"
              % (name[:28], nb, p.size, float((p < 0.05).mean()), float((p < 0.01).mean()),
                 rp, float((p <= rp).mean())))

    with open(a.out, "w") as fh:
        fh.write("real_region\tn_bins\tnull_chrom\tnull_start\tp_combined\n")
        for r in rows:
            fh.write("%s\t%d\t%s\t%d\t%.6g\n" % r)
    print("\nwrote %d null-region tests -> %s" % (len(rows), a.out))
    print("\n'<=real_p' is the fraction of INDEPENDENT null regions of the same size whose")
    print("p_combined is at least as extreme as the real region's. This is the number to")
    print("quote: the discrete binomial atoms make FPR at a fixed alpha hard to interpret.")


if __name__ == "__main__":
    main()
