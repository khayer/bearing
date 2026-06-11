#!/usr/bin/env python3
"""
prototype_adaptive_binning.py -- single-locus proof-of-concept for adaptive,
equal-coverage binning as an alternative to the fixed 200 bp grid.

Goal (advisor's idea): instead of uniform 200 bp bins, derive ONE consensus
segmentation in which every bin carries comparable coverage, computed once on
POOLED data (all samples, all tracks). Wide bins in sparse regions, narrow bins
where signal is dense -> every bin carries roughly equal statistical weight.

This is a deliberately scoped prototype to show before/after on one locus
(default Tcrb). It demonstrates the three things that matter for whether the
BEARING math survives non-uniform bins:

  1. ONE shared grid for all samples (so the differential A - B stays defined).
     The segmentation is built from pooled signal, NOT per-sample.
  2. Width-weighted background Q (so Q is not biased by how regions were chopped).
  3. The KL / JSD differential computed on variable-width bins.

It also reports the key calibration number: the coefficient of variation (CV)
of coverage-per-bin, fixed vs adaptive. Lower CV = more equal statistical
weight per bin, which is the mechanism that could improve the per-sample null
calibration that normalization alone did not fully fix.

LIMITATION (stated honestly): on a single locus the background Q is computed
LOCALLY (width-weighted over the locus), not genome-wide as in production. Both
the fixed and adaptive arms use the SAME local-Q convention, so the before/after
isolates the binning. Extending to a genome-wide Q is a follow-up if this looks
promising.

Real data:
  python prototype_adaptive_binning.py \
    --sheet workflow/results/samples.bearing.tsv \
    --categories categories/mm10_6track_panel.yaml \
    --chrom-sizes workflow/resources/mm10.chrom.sizes \
    --blacklist workflow/resources/mm10-blacklist.v2.bed \
    --locus chr6:40790000-41688054 \
    --compare DN,ProB --target-bins 1500 \
    --score-method kl --out-prefix results/adaptive_proto/tcrb

Self-test (no BigWigs, runs anywhere):
  python prototype_adaptive_binning.py --demo --out-prefix /tmp/adaptive_demo

ASCII only.
"""
import argparse
import os
import sys

import numpy as np

PSEUDOCOUNT = 1e-6
FINE_BP = 200            # fine grid used to build the coverage profile
NEG_STRAND_COL = 2       # 0-indexed column for RNAseq- (state 3): abs() applied


# --------------------------------------------------------------------------
# segmentation
# --------------------------------------------------------------------------
def equal_coverage_segments(coverage_fine, fine_edges, n_bins, keep_mask=None):
    """
    Build ~equal-coverage segments over the RETAINED fine bins.

    Segments break at any gap (run of dropped bins), so a bin never spans
    blacklisted / below-floor dead space -- the dead regions become gaps, exactly
    as the production min_signal floor excludes them, instead of being swallowed
    into one giant bin.

    coverage_fine : (n_fine,) pooled signal per fine bin (>= 0)
    fine_edges    : (n_fine + 1,) genomic coordinates of fine-bin boundaries
    keep_mask     : (n_fine,) bool; None = keep all
    returns       : list of (start_bp, end_bp, lo_fine, hi_fine) with hi exclusive
    """
    cov = np.asarray(coverage_fine, dtype=np.float64)
    n = len(cov)
    keep_idx = np.flatnonzero(np.ones(n, dtype=bool) if keep_mask is None
                              else np.asarray(keep_mask, dtype=bool))
    if keep_idx.size == 0:
        return []
    r_cov = cov[keep_idx]
    total = r_cov.sum()
    if total <= 0:
        bnd = np.linspace(0, len(keep_idx), n_bins + 1).round().astype(int)
    else:
        cum = np.concatenate([[0.0], np.cumsum(r_cov)])
        targets = np.linspace(0.0, total, n_bins + 1)
        bnd = np.searchsorted(cum, targets, side="left")
        bnd[0] = 0
        bnd[-1] = len(keep_idx)
        for k in range(1, len(bnd)):          # strictly increasing
            if bnd[k] <= bnd[k - 1]:
                bnd[k] = min(bnd[k - 1] + 1, len(keep_idx))
        bnd = np.minimum(bnd, len(keep_idx))

    segs = []
    for k in range(len(bnd) - 1):
        a, b = int(bnd[k]), int(bnd[k + 1])
        if b <= a:
            continue
        block = keep_idx[a:b]                  # retained fine indices in this quota
        # split the quota into contiguous runs so a bin never spans a gap
        cuts = np.flatnonzero(np.diff(block) != 1)
        starts = [0] + (cuts + 1).tolist()
        ends = (cuts + 1).tolist() + [len(block)]
        for s, e in zip(starts, ends):
            lo = int(block[s])
            hi = int(block[e - 1]) + 1
            segs.append((int(fine_edges[lo]), int(fine_edges[hi]), lo, hi))
    return [s for s in segs if s[3] > s[2]]


def aggregate_to_segments(signal_fine, fine_widths, segments):
    """
    Width-weighted mean signal per track per segment.

    signal_fine : (n_fine, n_tracks)
    fine_widths : (n_fine,) bp width of each fine bin
    returns     : (n_seg, n_tracks), (n_seg,) segment widths in bp
    """
    n_tracks = signal_fine.shape[1]
    out = np.zeros((len(segments), n_tracks), dtype=np.float64)
    widths = np.zeros(len(segments), dtype=np.float64)
    for j, (_s, _e, lo, hi) in enumerate(segments):
        w = fine_widths[lo:hi]
        wsum = w.sum()
        widths[j] = wsum
        if wsum > 0:
            out[j] = (signal_fine[lo:hi] * w[:, None]).sum(axis=0) / wsum
    return out, widths


# --------------------------------------------------------------------------
# scoring (mirrors bigwig_to_qcat: compositional P vs Q, clamped KL or JSD)
# --------------------------------------------------------------------------
def to_prob(signal_matrix):
    m = np.clip(signal_matrix, 0.0, None) + PSEUDOCOUNT
    return m / m.sum(axis=1, keepdims=True)


def width_weighted_Q(P, widths):
    w = np.asarray(widths, dtype=np.float64)
    return (P * w[:, None]).sum(axis=0) / w.sum()


def bearing_scores(P, Q, method="kl"):
    Qb = Q[np.newaxis, :]
    if method == "jsd":
        M = 0.5 * (P + Qb)
        termP = P * np.log2((P + 1e-300) / (M + 1e-300))
        termQ = Qb * np.log2((Qb + 1e-300) / (M + 1e-300))
        s = 0.5 * (termP + termQ)
        s = np.where(P > Qb, s, 0.0)
    else:
        s = P * np.log2(P / (Qb + 1e-300) + 1e-300)
    return np.clip(s, 0.0, None)


def cv(x):
    x = np.asarray(x, dtype=np.float64)
    m = x.mean()
    return float(x.std() / m) if m > 0 else float("nan")


def spearman(a, b):
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3:
        return float("nan")
    ra = np.argsort(np.argsort(a[ok]))
    rb = np.argsort(np.argsort(b[ok]))
    ra = ra - ra.mean()
    rb = rb - rb.mean()
    denom = np.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return float((ra * rb).sum() / denom) if denom > 0 else float("nan")


# --------------------------------------------------------------------------
# signal extraction
# --------------------------------------------------------------------------
def extract_real(sheet, categories_yaml, chrom_sizes_path, blacklist_path,
                 chrom, start, end, fine_bp=FINE_BP):
    """Return {sample: (n_fine, n_tracks) signal}, fine_edges. Reuses pipeline loaders."""
    from bigwig_to_qcat import mean_signal_in_bins, NEGATIVE_STRAND_STATES
    import pyBigWig

    neg_cols = sorted(s - 1 for s in NEGATIVE_STRAND_STATES)
    # fine grid over the locus
    edges = list(range(start, end, fine_bp)) + [end]
    bins = [(edges[i], edges[i + 1]) for i in range(len(edges) - 1)]
    fine_edges = np.array(edges, dtype=np.int64)

    per_sample = {}
    for name, bw_paths in sheet:
        cols = []
        for p in bw_paths:
            bw = pyBigWig.open(p)
            v = mean_signal_in_bins(bw, chrom, bins)
            bw.close()
            cols.append(np.asarray(v, dtype=np.float64))
        mat = np.vstack(cols).T            # (n_fine, n_tracks)
        for c in neg_cols:
            if c < mat.shape[1]:
                mat[:, c] = np.abs(mat[:, c])
        per_sample[name] = mat
    return per_sample, fine_edges


def make_demo(seed=42):
    """Synthetic locus: a focal dense region + broad sparse flanks, 2 conditions."""
    rng = np.random.default_rng(seed)
    n_fine = 1000                          # 200 kb at 200 bp
    start = 40790000
    fine_edges = np.array([start + i * FINE_BP for i in range(n_fine + 1)], dtype=np.int64)
    x = np.arange(n_fine)
    # pooled coverage shape: low baseline + two sharp peaks (dense regions)
    base = 0.3 + 0.0 * x
    peak1 = 6.0 * np.exp(-((x - 350) ** 2) / (2 * 12 ** 2))
    peak2 = 4.0 * np.exp(-((x - 640) ** 2) / (2 * 20 ** 2))
    shape = base + peak1 + peak2

    def sample(track_focus, diff_boost=0.0):
        mat = np.zeros((n_fine, 6))
        for t in range(6):
            amp = shape * (0.5 + 0.5 * ((t == track_focus)))
            if t == track_focus:
                amp = amp + diff_boost * peak1   # condition-specific enrichment at peak1
            mat[:, t] = np.clip(amp + rng.normal(0, 0.05, n_fine), 0, None)
        return mat

    per_sample = {
        "A_rep1": sample(3, diff_boost=2.5),   # CTCF-focused, enriched at peak1
        "A_rep2": sample(3, diff_boost=2.3),
        "B_rep1": sample(3, diff_boost=0.0),   # same focus, no peak1 enrichment
        "B_rep2": sample(3, diff_boost=0.2),
    }
    return per_sample, fine_edges, ("A", "B")


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
def condition_of(sample_name):
    return sample_name.split("_rep")[0].split("_")[0]


def mean_over_condition(per_sample, cond):
    mats = [m for s, m in per_sample.items() if condition_of(s) == cond]
    if not mats:
        raise SystemExit("no samples for condition '%s' (have %s)"
                         % (cond, sorted({condition_of(s) for s in per_sample})))
    return np.mean(mats, axis=0)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--demo", action="store_true", help="synthetic self-test (no BigWigs)")
    ap.add_argument("--sheet")
    ap.add_argument("--categories")
    ap.add_argument("--chrom-sizes")
    ap.add_argument("--blacklist")
    ap.add_argument("--locus", help="chrom:start-end")
    ap.add_argument("--compare", default="A,B", help="two conditions, e.g. DN,ProB")
    ap.add_argument("--target-bins", type=int, default=1500,
                    help="number of adaptive bins (consensus segmentation)")
    ap.add_argument("--fine-bp", type=int, default=200,
                    help="fine-grid resolution (bp) used to build the coverage "
                         "profile; this is the FLOOR on adaptive bin width, so "
                         "lower it (e.g. 50) to allow sub-200 bins in dense "
                         "regions. Ignored in --demo mode (always 200).")
    ap.add_argument("--min-signal", type=float, default=0.0,
                    help="drop fine bins below this summed-signal floor in EVERY "
                         "sample (production excludes dead regions this way; "
                         "use your config value, e.g. 0.1). Default 0 = off.")
    ap.add_argument("--score-method", choices=["kl", "jsd"], default="kl")
    ap.add_argument("--sweep", action="store_true",
                    help="run several --target-bins values and plot coverage-CV "
                         "and differential-preservation (Spearman vs fixed) vs "
                         "bin count, in one figure (resolution-vs-calibration curve).")
    ap.add_argument("--sweep-bins", default="400,800,1500,3000,6000",
                    help="comma-separated target-bin counts for --sweep")
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    if args.demo:
        per_sample, fine_edges, (condA, condB) = make_demo()
        chrom = "chr6"
    else:
        from bigwig_to_qcat import load_categories_yaml  # noqa: F401 (validates env)
        need = [args.sheet, args.categories, args.chrom_sizes, args.locus]
        if any(v is None for v in need):
            ap.error("real mode needs --sheet --categories --chrom-sizes --locus")
        sheet_dir = os.path.dirname(os.path.abspath(args.sheet))
        sheet = []
        with open(args.sheet) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            idx = {h: i for i, h in enumerate(header)}
            si, bi = idx.get("sample", 0), idx.get("bw", 1)
            for line in fh:
                line = line.rstrip("\n")
                if not line or line.startswith("#"):
                    continue
                f = line.split("\t")
                if len(f) <= bi:
                    continue
                paths = [p if os.path.isabs(p) else os.path.join(sheet_dir, p)
                         for p in (q.strip() for q in f[bi].split(",")) if p]
                sheet.append((f[si], paths))
        chrom, rng = args.locus.split(":")
        start, end = (int(v) for v in rng.split("-"))
        per_sample, fine_edges = extract_real(sheet, args.categories, args.chrom_sizes,
                                              args.blacklist, chrom, start, end,
                                              fine_bp=args.fine_bp)
        condA, condB = args.compare.split(",")

    fine_bp = 200 if args.demo else args.fine_bp

    fine_widths = np.diff(fine_edges).astype(np.float64)
    n_fine = len(fine_widths)

    # pooled coverage across ALL samples and tracks -> consensus segmentation
    pooled = np.zeros(n_fine, dtype=np.float64)
    for mat in per_sample.values():
        pooled += mat.sum(axis=1)

    # keep_mask: drop fine bins that production would not score, so equal-coverage
    # bins never span dead/blacklisted space (the source of giant merged bins).
    #  - blacklist: artifact / unmappable regions (when --blacklist given, non-demo)
    #  - min_signal floor: a bin is dropped only if it is below floor in EVERY
    #    sample (dead everywhere); a bin active in any sample is kept, because
    #    that is exactly where a differential can live.
    keep_mask = np.ones(n_fine, dtype=bool)
    n_bl = n_floor = 0
    if (not args.demo) and args.blacklist:
        from bigwig_to_qcat import load_blacklist, bins_overlapping_blacklist
        bl = load_blacklist(args.blacklist)
        fine_bins = [(int(fine_edges[i]), int(fine_edges[i + 1])) for i in range(n_fine)]
        bl_mask = bins_overlapping_blacklist(fine_bins, bl, chrom)
        n_bl = int(bl_mask.sum())
        keep_mask &= ~bl_mask
    if args.min_signal > 0:
        per_sample_total = np.vstack([m.sum(axis=1) for m in per_sample.values()])
        alive_anywhere = per_sample_total.max(axis=0) >= args.min_signal
        n_floor = int((~alive_anywhere & keep_mask).sum())
        keep_mask &= alive_anywhere

    segments = equal_coverage_segments(pooled, fine_edges, args.target_bins, keep_mask)

    keep_idx = np.flatnonzero(keep_mask)
    # fixed grid = each RETAINED fine bin (production scores only these)
    fixed_segments = [(int(fine_edges[i]), int(fine_edges[i + 1]), int(i), int(i + 1))
                      for i in keep_idx]

    # equal-COUNT uniform reference over the retained span: same bin count as
    # adaptive, uniform width, gap-aware. Isolates strategy from bin count.
    uniform_segments = equal_coverage_segments(
        np.ones(n_fine), fine_edges, len(segments), keep_mask)

    # coverage per bin, all grids
    cov_fixed = np.array([pooled[lo:hi].sum() for (_a, _b, lo, hi) in fixed_segments])
    cov_adapt = np.array([pooled[lo:hi].sum() for (_a, _b, lo, hi) in segments])
    cov_uni = np.array([pooled[lo:hi].sum() for (_a, _b, lo, hi) in uniform_segments])

    # per-condition mean signal -> P -> width-weighted Q -> scores, both grids
    def score_grid(segs):
        widths = np.array([s[1] - s[0] for s in segs], dtype=np.float64)
        out = {}
        for cond in (condA, condB):
            mat = mean_over_condition(per_sample, cond)        # (n_fine, tracks)
            agg, _w = aggregate_to_segments(mat, fine_widths, segs)
            P = to_prob(agg)
            Q = width_weighted_Q(P, widths)
            out[cond] = bearing_scores(P, Q, args.score_method).sum(axis=1)  # per-bin total
        diff = out[condA] - out[condB]
        return widths, out, diff

    w_fix, sc_fix, diff_fix = score_grid(fixed_segments)

    # ---- sweep mode: resolution-vs-calibration curve ----------------------
    if args.sweep:
        fixed_diff_keep = diff_fix              # per retained fine bin (keep_idx order)
        targets = [int(t) for t in args.sweep_bins.split(",")]
        rows = []
        for t in targets:
            segs_t = equal_coverage_segments(pooled, fine_edges, t, keep_mask)
            wt, _sc, diff_t = score_grid(segs_t)
            # map adaptive diff back onto the fine grid, then onto retained bins
            adiff_fine = np.full(n_fine, np.nan)
            for j, (_s, _e, lo, hi) in enumerate(segs_t):
                adiff_fine[lo:hi] = diff_t[j]
            rho = spearman(adiff_fine[keep_idx], fixed_diff_keep)
            cov_t = np.array([pooled[lo:hi].sum() for (_a, _b, lo, hi) in segs_t])
            rows.append((t, len(segs_t), cv(cov_t), rho, int(np.median(wt)), int(wt.max())))
            print("target=%5d  bins=%5d  covCV=%.3f  diff_rho=%.3f  medW=%d  maxW=%d"
                  % rows[-1])
        os.makedirs(os.path.dirname(os.path.abspath(args.out_prefix)) or ".", exist_ok=True)
        stsv = args.out_prefix + "_sweep.tsv"
        with open(stsv, "w") as fh:
            fh.write("target_bins\tactual_bins\tcoverage_cv\tdiff_spearman\tmedian_width\tmax_width\n")
            for r in rows:
                fh.write("%d\t%d\t%.5f\t%.5f\t%d\t%d\n" % r)
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            nb = [r[1] for r in rows]
            cvs = [r[2] for r in rows]
            rhos = [r[3] for r in rows]
            fig, axL = plt.subplots(figsize=(8, 5))
            axR = axL.twinx()
            l1 = axL.plot(nb, cvs, "o-", color="#C2410C", label="coverage-per-bin CV")
            l2 = axR.plot(nb, rhos, "s-", color="#1D5C8A",
                          label="differential preservation (Spearman vs fixed)")
            axL.set_xlabel("number of adaptive bins")
            axL.set_ylabel("coverage CV (lower = better calibration)", color="#C2410C")
            axR.set_ylabel("diff Spearman (higher = better resolution)", color="#1D5C8A")
            axL.set_title("Resolution vs calibration tradeoff (%s, %s-%s)"
                          % (chrom, condA, condB), fontweight="bold")
            axL.grid(alpha=0.25)
            lines = l1 + l2
            axL.legend(lines, [ln.get_label() for ln in lines], fontsize=9, loc="center right")
            fig.tight_layout()
            fig.savefig(args.out_prefix + "_sweep.pdf", dpi=150)
            fig.savefig(args.out_prefix + "_sweep.png", dpi=130)
            print("wrote %s, %s_sweep.pdf/png" % (stsv, args.out_prefix))
        except Exception as exc:               # noqa: BLE001
            print("wrote %s (figure skipped: %s)" % (stsv, exc))
        return

    w_ada, sc_ada, diff_ada = score_grid(segments)

    cv_fixed, cv_adapt = cv(cov_fixed), cv(cov_adapt)
    cv_uni = cv(cov_uni)
    print("locus            : %s:%d-%d" % (chrom, fine_edges[0], fine_edges[-1]))
    print("excluded fine bins: %d blacklisted, %d below-floor (dead in all samples); "
          "%d of %d retained" % (n_bl, n_floor, int(keep_mask.sum()), n_fine))
    print("fixed bins       : %d (uniform %d bp, retained only)" % (len(fixed_segments), fine_bp))
    print("adaptive bins    : %d (median width %d bp, range %d-%d)"
          % (len(segments), int(np.median(w_ada)), int(w_ada.min()), int(w_ada.max())))
    print("coverage-per-bin CV:")
    print("  fixed %3d bp           (%5d bins) : %.3f" % (fine_bp, len(fixed_segments), cv_fixed))
    print("  uniform, matched count (%5d bins) : %.3f" % (len(uniform_segments), cv_uni))
    print("  adaptive equal-coverage(%5d bins) : %.3f  <- most equal weight" % (len(segments), cv_adapt))
    print("  adaptive vs matched-uniform CV reduction: %.1fx"
          % (cv_uni / cv_adapt if cv_adapt > 0 else float("nan")))
    print("differential (%s - %s) preserved on shared adaptive grid: "
          "max|diff| fixed=%.3f adaptive=%.3f" % (condA, condB,
          float(np.abs(diff_fix).max()), float(np.abs(diff_ada).max())))

    # write adaptive bins TSV
    os.makedirs(os.path.dirname(os.path.abspath(args.out_prefix)) or ".", exist_ok=True)
    tsv = args.out_prefix + "_adaptive_bins.tsv"
    with open(tsv, "w") as fh:
        fh.write("chrom\tstart\tend\twidth_bp\tpooled_cov\t%s_score\t%s_score\tdiff\n"
                 % (condA, condB))
        for j, (s, e, _lo, _hi) in enumerate(segments):
            fh.write("%s\t%d\t%d\t%d\t%.4f\t%.5f\t%.5f\t%.5f\n"
                     % (chrom, s, e, e - s, cov_adapt[j],
                        sc_ada[condA][j], sc_ada[condB][j], diff_ada[j]))

    # figure
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        xt = (fine_edges[:-1] + fine_edges[1:]) / 2.0 / 1e6   # Mb
        seg_mid = np.array([(s + e) / 2.0 for (s, e, _l, _h) in segments]) / 1e6
        seg_w_mb = w_ada / 1e6
        fig, ax = plt.subplots(4, 1, figsize=(11, 10), sharex=True)

        ax[0].fill_between(xt, pooled, step="mid", color="#0E7C7B", alpha=0.85, lw=0)
        ax[0].set_ylabel("pooled\ncoverage")
        ax[0].set_title("Equal-coverage adaptive binning prototype (%s:%d-%d)"
                        % (chrom, fine_edges[0], fine_edges[-1]), fontweight="bold")

        for b in fine_edges[::max(1, n_fine // 80)]:
            ax[1].axvline(b / 1e6, color="#94a3b8", lw=0.3)
        for (s, _e, _l, _h) in segments:
            ax[1].axvline(s / 1e6, color="#C2410C", lw=0.4)
        ax[1].set_ylabel("bin edges")
        ax[1].set_yticks([])
        ax[1].text(0.005, 0.78, "grey = fixed %d bp" % fine_bp, transform=ax[1].transAxes,
                   fontsize=8, color="#64748b")
        ax[1].text(0.005, 0.55, "orange = adaptive (dense where coverage high)",
                   transform=ax[1].transAxes, fontsize=8, color="#C2410C")

        # map retained-only fixed arrays back onto the full grid so they align
        # with xt and break (NaN) across dropped/dead regions
        cov_fixed_full = np.full(n_fine, np.nan)
        cov_fixed_full[keep_idx] = cov_fixed
        diff_fix_full = np.full(n_fine, np.nan)
        diff_fix_full[keep_idx] = diff_fix

        ax[2].plot(xt, cov_fixed_full, color="#94a3b8", lw=0.8,
                   label="fixed %d bp (CV=%.2f)" % (fine_bp, cv_fixed))
        ax[2].bar(seg_mid, cov_adapt, width=seg_w_mb, color="#C2410C", alpha=0.55,
                  label="adaptive (CV=%.2f)" % cv_adapt)
        ax[2].set_ylabel("coverage\nper bin")
        ax[2].legend(fontsize=8, loc="upper right")

        ax[3].step(xt, diff_fix_full, where="mid", color="#aab2bd", lw=0.7,
                   alpha=0.7, label="fixed %d bp" % fine_bp)
        ax[3].bar(seg_mid, diff_ada, width=seg_w_mb, color="#C2410C", alpha=0.8,
                  label="adaptive")
        ax[3].axhline(0, color="k", lw=0.5)
        ax[3].set_ylabel("BEARING diff\n(%s - %s)" % (condA, condB))
        ax[3].set_xlabel("position (Mb)")
        ax[3].legend(fontsize=8, loc="upper right")

        fig.tight_layout()
        pdf = args.out_prefix + "_prototype.pdf"
        fig.savefig(pdf, dpi=150)
        png = args.out_prefix + "_prototype.png"
        fig.savefig(png, dpi=130)
        print("wrote %s, %s, %s" % (tsv, pdf, png))
    except Exception as exc:               # noqa: BLE001
        print("wrote %s (figure skipped: %s)" % (tsv, exc))


if __name__ == "__main__":
    main()
