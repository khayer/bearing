#!/usr/bin/env python3
"""
track_ablation.py
=================
Does each regional call survive when the PUBLIC architectural tracks are removed?

WHY
For Pro-B and S3T3, the CTCF, RAD21 (cohesin) and H3K27ac tracks are public
ChIP-seq from other labs -- and two RAD21 tracks are from MEF, not 3T3 (Table
S4). A cross-condition differential could therefore be driven by platform/cell-
type differences in those tracks rather than by biology. The reviewer asked for
a drop-one-(track-group) ablation. This runs it.

HOW (no re-scoring needed)
The production per-bin differential score is the SUM of the six per-track KL
contributions (bearing_pvalue.py: score_total = sum(per_track.values())). Those
contributions are already stored per bin in the qcat. Dropping a track means
summing only the retained tracks -- for BOTH the observed scores AND the
permutation-null scores, so the null is re-derived for the reduced panel rather
than reused. The regional test is then applied via the SAME production code.

Default retained panel: the three IN-HOUSE tracks (ATAC, RNAseq+, RNAseq-),
i.e. drop CTCF + Cohesin + H3K27ac (indices 4,5,6 in the 6-track panel).

  conda activate bearing
  python track_ablation.py \\
      --diff-qcat  workflow/results/compare/diff_DN_vs_ProB.qcat.bgz \\
      --perm-glob 'workflow/results/perm/perm*/diff_comparison/diff_DN_vs_ProB.qcat.bgz' \\
      --regions   workflow/annotations/igh_regions_v4.bed \\
      --locus     chr12:113000000-116100000 \\
      --keep 1,2,3 --p-thresh 0.05 \\
      --out track_ablation_igh_DN_vs_ProB_inhouse.tsv

Compare the p_combined / q_combined it prints against the full-panel production
values. A call that stays significant on ATAC+RNA alone is NOT a public-ChIP-seq
artifact.

ASCII only. Reads real data; fabricates nothing.
"""

import argparse
import glob
import gzip
import math
import os
import sys

import json

import numpy as np
import scipy.stats


def parse_qcat_locus(path, keep, lc, ls, le):
    """Yield (start, sum_over_keep, sign) for bins on chrom lc within [ls, le).

    Filters on the tab-delimited chrom/start BEFORE json-parsing the payload, so
    only in-locus lines are decoded -- essential when the perm qcats are
    genome-wide (13.6M rows) but the locus is < 10k bins. Payload format matches
    bearing_pvalue.parse_qcat: a 4th column that is either a JSON object with a
    "qcat" key, or "qcat:[[score,idx],...]" optionally followed by ",raw:".
    keep : 1-based track indices to retain.
    """
    keepset = set(keep)
    lc_tab = lc + "\t"
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if not line.startswith(lc_tab):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            start = int(f[1])
            if start < ls or start >= le:
                continue
            col = f[3]
            if col.startswith("{"):
                pairs = json.loads(col).get("qcat", [])
            else:
                qs = col.find("qcat:")
                if qs == -1:
                    continue
                rs = col.find(",raw:", qs)
                payload = col[qs + 5:rs] if rs >= 0 else col[qs + 5:]
                pairs = json.loads(payload)
            s = 0.0
            for score, idx in pairs:
                if int(idx) in keepset:
                    s += float(score)
            yield start, s, (1 if s > 0 else -1)


def _verify_parser(repo, path, keep, lc, ls, le, n=200):
    """Sanity: the fast locus reader must match bearing_pvalue.parse_qcat exactly."""
    sys.path.insert(0, repo)
    from bearing_pvalue import parse_qcat
    keepset = set(keep)
    ref = {}
    for chrom, start, _e, _t, per in parse_qcat(path, min_signal=0.0):
        if chrom == lc and ls <= start < le:
            ref[start] = sum(v for i, v in per.items() if i in keepset)
            if len(ref) >= n:
                break
    checked = 0
    for start, s, _sg in parse_qcat_locus(path, keep, lc, ls, le):
        if start in ref:
            if abs(ref[start] - s) > 1e-9:
                sys.exit("locus reader disagrees with production parse_qcat at %d" % start)
            checked += 1
            if checked >= len(ref):
                break
    print("  parser check: %d bins match production parse_qcat exactly" % checked)


def parse_qcat_genome(path, keep):
    """Yield (reduced_score, sign) for EVERY bin, summing over kept tracks.
    Genome-wide pass; used only to derive g_dir exactly as production does."""
    keepset = set(keep)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            tab1 = line.find("\t")
            if tab1 < 0:
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            col = f[3]
            if col.startswith("{"):
                pairs = json.loads(col).get("qcat", [])
            else:
                qs = col.find("qcat:")
                if qs == -1:
                    continue
                rs = col.find(",raw:", qs)
                payload = col[qs + 5:rs] if rs >= 0 else col[qs + 5:]
                pairs = json.loads(payload)
            s = 0.0
            for score, idx in pairs:
                if int(idx) in keepset:
                    s += float(score)
            yield s, (1 if s > 0 else -1)


def genome_wide_g_dir(diff_path, keep, null_sorted, p_thresh):
    """Reduced-panel g_dir = P(reduced_score > 0 | bin genome-wide significant),
    matching regional_enrichment.py (which derives g_dir from ALL genome-wide
    significant bins, not locus bins). Computing it from locus bins gives a
    degenerate baseline and makes negative controls light up.

    Null equivalence: the permutation nulls are per-track CIRCULAR SHIFTS, which
    decorrelate position from signal, so the reduced-score null pooled from the
    locus is identically distributed to the genome-wide reduced null. We reuse
    the locus-pooled null (null_sorted) as the genome-wide reference, avoiding a
    second genome-wide pass over the 100 perm files."""
    pos = tot = 0
    n = null_sorted.size
    for s, sign in parse_qcat_genome(diff_path, keep):
        ge = n - np.searchsorted(null_sorted, abs(s), side="left")
        p = (1.0 + ge) / (1.0 + n)
        if p < p_thresh:
            tot += 1
            if sign > 0:
                pos += 1
    g = (pos / tot) if tot else 0.5
    return min(max(g, 1e-6), 1.0 - 1e-6), tot


def load_observed(path, keep, locus):
    lc, ls, le = locus
    starts, scores, signs = [], [], []
    for start, s, sign in parse_qcat_locus(path, keep, lc, ls, le):
        starts.append(start); scores.append(s); signs.append(sign)
    o = np.argsort(starts, kind="stable")
    return (np.asarray(starts)[o], np.asarray(scores)[o], np.asarray(signs)[o])


def load_null(perm_paths, keep, locus, max_perms):
    """Pooled null of |reduced score| within the locus, from the permutation qcats."""
    lc, ls, le = locus
    vals = []
    for i, p in enumerate(perm_paths[:max_perms]):
        for start, s, _sign in parse_qcat_locus(p, keep, lc, ls, le):
            vals.append(abs(s))
        if (i + 1) % 10 == 0:
            print("    ... %d/%d perm files, %d null bins so far"
                  % (i + 1, min(len(perm_paths), max_perms), len(vals)))
    return np.sort(np.asarray(vals))


def empirical_p(obs_abs, null_sorted):
    """(1 + #null >= obs) / (1 + N), the production convention."""
    n = null_sorted.size
    ge = n - np.searchsorted(null_sorted, obs_abs, side="left")
    return (1.0 + ge) / (1.0 + n)


def regional_test(starts, pvals, signs, region, p_thresh, g_dir):
    """regional_enrichment.compute_regional_enrichment math for one region."""
    _c, rs, re_, _name = region
    L_locus = starts.size
    n_locus = int((pvals < p_thresh).sum())
    inR = (starts >= rs) & (starts < re_)
    L_region = int(inR.sum())
    m = inR & (pvals < p_thresh)
    k = int(m.sum())
    k_pos = int((signs[m] > 0).sum())
    k_neg = k - k_pos
    pi = (L_region / L_locus) if L_locus else 0.0

    p_spatial = 1.0 if n_locus == 0 else min(
        float(scipy.stats.binom.sf(k - 1, n=n_locus, p=pi)), 1.0)
    if k == 0:
        p_dir = 1.0
    elif k_pos >= k_neg:
        p_dir = min(float(scipy.stats.binom.sf(k_pos - 1, n=k, p=g_dir)), 1.0)
    else:
        p_dir = min(float(scipy.stats.binom.sf(k_neg - 1, n=k, p=1.0 - g_dir)), 1.0)
    if p_spatial == 0 or p_dir == 0:
        p_comb = 0.0
    else:
        T = -2.0 * (math.log(p_spatial) + math.log(p_dir))
        p_comb = float(scipy.stats.chi2.sf(T, df=4))
    return dict(L_region=L_region, k=k, k_pos=k_pos, k_neg=k_neg,
                p_spatial=p_spatial, p_directional=p_dir, p_combined=p_comb)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--diff-qcat", required=True)
    ap.add_argument("--perm-glob", required=True,
                    help="glob for the DIFFERENTIAL permutation qcats of THIS contrast, e.g. "
                         "'workflow/results/perm/perm*/diff_comparison/diff_DN_vs_ProB.qcat.bgz'")
    ap.add_argument("--repo", default=".", help="repo root (imports bearing_pvalue.parse_qcat)")
    ap.add_argument("--regions", required=True)
    ap.add_argument("--locus", required=True)
    ap.add_argument("--keep", default="1,2,3",
                    help="1-based track indices to retain (default in-house: ATAC,RNA+,RNA-)")
    ap.add_argument("--p-thresh", type=float, default=0.05)
    ap.add_argument("--max-perms", type=int, default=100)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    keep = [int(x) for x in a.keep.split(",")]
    lc, rest = a.locus.split(":"); ls, le = (int(v) for v in rest.split("-"))
    locus = (lc, ls, le)
    print("retaining tracks %s; locus %s" % (keep, a.locus))

    print("loading observed reduced scores ...")
    _verify_parser(os.path.abspath(a.repo), a.diff_qcat, keep, lc, ls, le)
    starts, scores, signs = load_observed(a.diff_qcat, keep, locus)
    print("  %d bins in locus" % starts.size)

    perm_paths = sorted(glob.glob(a.perm_glob))
    if not perm_paths:
        sys.exit("no permutation qcats matched --perm-glob")
    print("loading pooled null from %d permutation files (cap %d) ..."
          % (len(perm_paths), a.max_perms))
    null = load_null(perm_paths, keep, locus, a.max_perms)
    print("  %d null bins" % null.size)

    # per-bin empirical p on the reduced score (locus bins), for the spatial test
    pvals = empirical_p(np.abs(scores), null)

    # g_dir MUST be genome-wide (production uses all genome-wide significant bins,
    # not locus bins) or the directional test is computed against a degenerate
    # baseline and negative controls light up. Stream the observed diff qcat
    # genome-wide, thresholding |reduced score| against the pooled reduced null.
    print("  computing genome-wide reduced-panel g_dir (one pass over observed diff) ...")
    g_dir, n_gw_sig = genome_wide_g_dir(a.diff_qcat, keep, null, a.p_thresh)
    print("  reduced-panel g_dir = %.4f (%d genome-wide significant bins)\n"
          % (g_dir, n_gw_sig))

    real = []
    with open(a.regions) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4 and f[0] == lc:
                real.append((f[0], int(f[1]), int(f[2]), f[3]))

    print("%-26s %6s %5s %14s %14s"
          % ("region", "k", "dir", "p_spatial", "p_combined"))
    rows = []
    for region in real:
        r = regional_test(starts, pvals, signs, region, a.p_thresh, g_dir)
        print("%-26s %6d %5s %14.3g %14.3g"
              % (region[3][:26], r["k"], "%d+/%d-" % (r["k_pos"], r["k_neg"]),
                 r["p_spatial"], r["p_combined"]))
        rows.append((region[3], r["L_region"], r["k"], r["k_pos"], r["k_neg"],
                     r["p_spatial"], r["p_directional"], r["p_combined"]))

    with open(a.out, "w") as fh:
        fh.write("region\tL_region\tk\tk_pos\tk_neg\tp_spatial\tp_directional\tp_combined\n")
        for r in rows:
            fh.write("%s\t%d\t%d\t%d\t%d\t%.6g\t%.6g\t%.6g\n" % r)
    print("\nwrote %s" % a.out)
    print("Compare p_combined against the full-panel production values. A call that")
    print("stays significant on the in-house tracks alone is not a public-ChIP artifact.")
    print("NOTE: nominal p_combined is anti-conservative (see Table S11); for a call")
    print("near threshold, re-run regional_null_calibration_v2.py on the reduced panel.")


if __name__ == "__main__":
    main()
