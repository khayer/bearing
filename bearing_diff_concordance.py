#!/usr/bin/env python3
"""
bearing_diff_concordance.py -- concordance between BEARING per-track differential
bins and an independent count-based differential test (edgeR/Rsubread allbins CSV).

Built to make the RNA (and CTCF/H3K27ac) triangulation a single command rather
than ad-hoc parsing. It handles, automatically:

  * the edgeR allbins CSV export quirks: a UTF-8 BOM, thousands separators and
    quotes in the coordinate columns (" 71,336,801 "), and trailing all-blank
    rows that some write.csv/Excel round-trips leave behind;
  * the 0-based-BED vs 1-based-GRanges off-by-one: BEARING BEDs are 0-based
    half-open (start = ...800) while coords that went through import(BED) in R
    come back 1-based (start = ...801). The bin END coincides in both conventions
    for a contiguous bin, so bins are matched on (chrom, end). (This assumes both
    sides are on the SAME bin grid, which they are when the edgeR bins were built
    from make_bearing_bins_bed.py.)

Reports, restricted to bins edgeR actually tested (post-filterByExpr universe):
  * recovery   : BEARING sig bins that are edgeR-significant (FDR < --fdr)
  * direction  : agreement between BEARING up/dn (BED col4 ..._up/..._dn) and the
                 sign of edgeR logFC (flags a flipped factor convention)
  * hypergeom  : enrichment of BEARING bins among edgeR-sig bins
  * top-X%     : rank-based enrichment (BEARING bins in edgeR's top --rank-pct by
                 p-value) -- the statistic to trust when edgeR's universe is
                 saturated (a large fraction FDR-significant makes hypergeom
                 uninformative).

  python bearing_diff_concordance.py \
    --bearing-bed bearing_DNvsDP_RNApos_significant_bins.bed \
    --edger-csv   rna_concordance_DNvsDP_edgeR_allbins.csv \
    --fdr 0.05

ASCII only.
"""
import argparse
import csv
import math
import sys


def _num(x):
    return x.replace(",", "").replace('"', "").strip()


def load_edger(path):
    """Parse an edgeR allbins CSV. Returns list of (chrom, start, end, logFC,
    pvalue, fdr), skipping blank/garbage rows. Column lookup is by header name."""
    rows = []
    with open(path, newline="") as fh:
        reader = csv.reader(fh)
        header = None
        for raw in reader:
            if header is None:
                header = [h.strip().lstrip("\ufeff").lower() for h in raw]
                idx = {h: i for i, h in enumerate(header)}
                need = ("logfc", "pvalue", "fdr", "chr", "start", "end")
                miss = [c for c in need if c not in idx]
                if miss:
                    sys.exit("[ERROR] edgeR CSV missing column(s): %s\n"
                             "        header was: %s" % (", ".join(miss), header))
                continue
            if len(raw) < len(header):
                continue
            chrom = _num(raw[idx["chr"]])
            if not chrom:
                continue
            try:
                s = int(_num(raw[idx["start"]]))
                e = int(_num(raw[idx["end"]]))
                lfc = float(_num(raw[idx["logfc"]]))
                pv = float(_num(raw[idx["pvalue"]]))
                fdr = float(_num(raw[idx["fdr"]]))
            except ValueError:
                continue
            rows.append((chrom, s, e, lfc, pv, fdr))
    return rows


def load_bearing(path, chrom=None):
    """Parse a directional BEARING BED6. Returns list of (chrom, start, end, up)
    where up=True for ..._up (A/first-condition enriched), False for ..._dn."""
    out = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            c = f[0]
            if chrom and c != chrom:
                continue
            try:
                s, e = int(f[1]), int(f[2])
            except ValueError:
                continue
            name = f[3]
            up = name.endswith("_up")
            dn = name.endswith("_dn")
            if not (up or dn):
                up = None  # undirected; counted but excluded from direction calc
            out.append((c, s, e, up))
    return out


def _hyperge_sf(k, N, K, n):
    """P(X >= k) for hypergeometric(N, K, n)."""
    if n == 0 or K == 0:
        return float("nan")

    def lc(a, b):
        if b < 0 or b > a:
            return -math.inf
        return math.lgamma(a + 1) - math.lgamma(b + 1) - math.lgamma(a - b + 1)
    denom = lc(N, n)
    s = 0.0
    for i in range(k, min(K, n) + 1):
        s += math.exp(lc(K, i) + lc(N - K, n - i) - denom)
    return min(1.0, s)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bearing-bed", required=True,
                    help="directional BEARING sig BED (col4 ..._up/..._dn)")
    ap.add_argument("--edger-csv", required=True,
                    help="edgeR allbins CSV (logFC,...,FDR,chr,start,end)")
    ap.add_argument("--fdr", type=float, default=0.05,
                    help="edgeR FDR cutoff for significance (default 0.05)")
    ap.add_argument("--rank-pct", type=float, default=10.0,
                    help="top%% of edgeR universe by p-value for rank enrichment "
                         "(default 10)")
    ap.add_argument("--chrom", default=None,
                    help="restrict to one chromosome (default: all shared)")
    ap.add_argument("--out", default=None, help="optional TSV summary path")
    ap.add_argument("--dump-missed", default=None,
                    help="optional TSV: BEARING-testable bins edgeR did NOT call "
                         "significant, with their edgeR logFC/pvalue/FDR and BEARING "
                         "direction. Use to check whether misses cluster just above "
                         "the FDR line (power ceiling) vs scatter to FDR~1 (no signal).")
    args = ap.parse_args()

    edger = load_edger(args.edger_csv)
    if args.chrom:
        edger = [r for r in edger if r[0] == args.chrom]
    if not edger:
        sys.exit("[ERROR] no edgeR rows parsed (after optional --chrom filter).")

    # match key: (chrom, end) -- robust to 0-based BED vs 1-based GRanges start
    tested = {}      # (chrom,end) -> (logFC, pvalue, fdr)
    widths = {}
    for c, s, e, lfc, pv, fdr in edger:
        tested[(c, e)] = (lfc, pv, fdr)
        widths[e - s] = widths.get(e - s, 0) + 1
    sig = {ke: v for ke, v in tested.items() if v[2] < args.fdr}

    # edgeR top-rank set by p-value
    n_top = max(1, int(round(len(tested) * args.rank_pct / 100.0)))
    top_keys = set(k for k, _ in sorted(
        tested.items(), key=lambda kv: kv[1][1])[:n_top])

    bear = load_bearing(args.bearing_bed, chrom=args.chrom)
    in_u = [b for b in bear if (b[0], b[2]) in tested]
    hit = [b for b in in_u if (b[0], b[2]) in sig]
    top_hit = [b for b in in_u if (b[0], b[2]) in top_keys]

    if args.dump_missed:
        missed = [b for b in in_u if (b[0], b[2]) not in sig]
        with open(args.dump_missed, "w") as fh:
            fh.write("chrom\tstart\tend\tbearing_dir\tedger_logFC\t"
                     "edger_pvalue\tedger_fdr\n")
            for c, s, e, up in sorted(missed, key=lambda b: tested[(b[0], b[2])][2]):
                lfc, pv, fdr = tested[(c, e)]
                bdir = "up" if up is True else ("dn" if up is False else "NA")
                fh.write("%s\t%d\t%d\t%s\t%.4f\t%.4g\t%.4g\n"
                         % (c, s, e, bdir, lfc, pv, fdr))
        print("wrote %d missed bins -> %s" % (len(missed), args.dump_missed))

    # direction agreement among recovered bins
    dir_tot = dir_agree = 0
    for b in hit:
        if b[3] is None:
            continue
        lfc = sig[(b[0], b[2])][0]
        dir_tot += 1
        dir_agree += int(b[3] == (lfc > 0))

    N, K, n, k = len(tested), len(sig), len(in_u), len(hit)
    p_hyper = _hyperge_sf(k, N, K, n)
    exp_rand = n * K / N if N else float("nan")
    # rank enrichment: BEARING in-universe bins landing in edgeR top-rank set
    exp_top = n * args.rank_pct / 100.0
    p_top = _hyperge_sf(len(top_hit), N, n_top, n)

    common_w = max(widths, key=widths.get) if widths else 0
    scope = args.chrom if args.chrom else "all shared chroms"

    print("=" * 76)
    print("BEARING vs edgeR differential concordance  (%s)" % scope)
    print("-" * 76)
    print("  edgeR tested universe (post-filter): %7d  bin width mode %d bp" % (N, common_w))
    print("  edgeR significant (FDR < %.3g):       %7d  (%.0f%% of universe)"
          % (args.fdr, K, 100.0 * K / N if N else float("nan")))
    print("  BEARING sig bins (this scope):        %7d" % len(bear))
    print("    in edgeR tested universe:           %7d  (%.0f%%)"
          % (n, 100.0 * n / len(bear) if bear else float("nan")))
    print("    dropped by edgeR filterByExpr:      %7d  (low-count -> not testable)"
          % (len(bear) - n))
    print("-" * 76)
    print("  RECOVERY  edgeR-significant:   %d / %d  (%.1f%%)"
          % (k, n, 100.0 * k / n if n else float("nan")))
    if dir_tot:
        frac = dir_agree / dir_tot
        flip = "  [LOOKS FLIPPED: edgeR factor convention reversed]" if frac < 0.5 else ""
        print("  DIRECTION agreement:           %d / %d  (%.0f%%)%s"
              % (dir_agree, dir_tot, 100.0 * frac, flip))
    print("  HYPERGEOM enrichment:          p = %.3e   (obs %d vs exp %.1f, %.1fx)"
          % (p_hyper, k, exp_rand, k / exp_rand if exp_rand else float("nan")))
    if K > 0.5 * N:
        print("     ^ universe is %.0f%% significant (saturated) -- trust the rank"
              " enrichment below, not this." % (100.0 * K / N))
    print("  RANK enrichment (top %.0f%%):     p = %.3e   (obs %d vs exp %.1f, %.1fx)"
          % (args.rank_pct, p_top, len(top_hit), exp_top,
             len(top_hit) / exp_top if exp_top else float("nan")))
    print("=" * 76)

    if args.out:
        with open(args.out, "w") as fh:
            fh.write("scope\tedger_tested\tedger_sig\tedger_sig_frac\tbearing_sig\t"
                     "in_universe\trecovered\trecovery_frac\tdir_agree\tdir_total\t"
                     "hyperg_p\trank_pct\trank_hit\trank_p\n")
            fh.write("%s\t%d\t%d\t%.4f\t%d\t%d\t%d\t%.4f\t%d\t%d\t%.3e\t%.1f\t%d\t%.3e\n"
                     % (scope, N, K, K / N if N else float("nan"), len(bear), n, k,
                        k / n if n else float("nan"), dir_agree, dir_tot, p_hyper,
                        args.rank_pct, len(top_hit), p_top))
        print("wrote %s" % args.out)


if __name__ == "__main__":
    main()
