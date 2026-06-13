#!/usr/bin/env python3
"""
compare_adaptive_vs_default.py -- compare adaptive (variable-width) differential
calls against the default fixed-200bp calls.

The two runs do NOT share a bin grid (a 600 bp adaptive bin spans three fixed
bins; a 50 bp one sits inside one), so per-bin set-Jaccard (compare_default_vs_
qnorm.py) is undefined here. Instead this compares on genomic COORDINATES: does a
different binning flag the same differential REGIONS, in the same direction?

Significance is a SHARED FDR cutoff (pval_adj_bh < --fdr) applied identically to
both runs -- not each run's own significant_fdr0.05 flag -- because the adaptive
grid has a different bin count and therefore a different BH burden, and a single
alpha keeps the significance rule consistent across the two. (The BH burden still
differs with bin count; that is an inherent property of the binning, reported via
the separate significant-bp totals.)

Both sides are restricted to main chromosomes (^chr([0-9]+|X|Y)$), since the
adaptive grid is main-chrom only and the fixed run is genome-wide.

Per comparison it reports:
  sig_bp_def / sig_bp_ada : significant base pairs each side
  bp_jaccard              : shared significant bp / union (any direction)
  bp_jaccard_dir          : shared significant bp with SAME direction / union
  dir_agree               : of bp significant in both, fraction same direction
  regions_def/ada         : significant regions (merged adjacent sig bins)
  recall_def->ada         : fraction of default sig regions hit by an adaptive
                            sig region, and vice versa

  python compare_adaptive_vs_default.py \
    --default-dir  workflow/results/pvalue \
    --adaptive-dir workflow/results_adaptive/pvalue \
    --fdr 0.05

ASCII only.
"""
import argparse
import glob
import os
import re
import sys

MAIN_CHROM = re.compile(r"^chr([0-9]+|X|Y)$")
# diff stats columns (0-based)
C_CHROM, C_START, C_END, C_SCORE, C_QADJ = 0, 1, 2, 3, 8


def load_sig(path, fdr):
    """Stream a diff stats TSV; return {chrom: [(start, end, dir_sign), ...]} for
    main-chrom bins with pval_adj_bh < fdr."""
    out = {}
    with open(path) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= C_QADJ:
                continue
            chrom = f[C_CHROM]
            if not MAIN_CHROM.match(chrom):
                continue
            try:
                q = float(f[C_QADJ])
                if q >= fdr:
                    continue
                s, e = int(f[C_START]), int(f[C_END])
                sc = float(f[C_SCORE])
            except ValueError:
                continue
            sign = 1 if sc > 0 else (-1 if sc < 0 else 0)
            out.setdefault(chrom, []).append((s, e, sign))
    return out


def merge(intervals):
    """Merge overlapping/adjacent (start, end) intervals."""
    if not intervals:
        return []
    iv = sorted(intervals)
    out = [list(iv[0])]
    for s, e in iv[1:]:
        if s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def overlap_len(a, b):
    """Total intersection length of two sorted, non-overlapping interval lists."""
    i = j = tot = 0
    while i < len(a) and j < len(b):
        s = max(a[i][0], b[j][0])
        e = min(a[i][1], b[j][1])
        if e > s:
            tot += e - s
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return tot


def total_len(intervals):
    return sum(e - s for s, e in intervals)


def n_regions_hit(regions, other_merged):
    """How many of `regions` overlap any interval in other_merged (>0 bp)."""
    if not regions or not other_merged:
        return 0
    hit = 0
    j = 0
    om = other_merged
    for (s, e) in regions:
        while j < len(om) and om[j][1] <= s:
            j += 1
        k = j
        ov = False
        while k < len(om) and om[k][0] < e:
            if min(e, om[k][1]) > max(s, om[k][0]):
                ov = True
                break
            k += 1
        hit += int(ov)
    return hit


def compare_one(def_sig, ada_sig):
    chroms = set(def_sig) | set(ada_sig)
    sig_bp_def = sig_bp_ada = 0
    ov_any = ov_dir = union = 0
    reg_def = reg_ada = 0
    hit_def = hit_ada = 0
    for c in chroms:
        d = def_sig.get(c, [])
        a = ada_sig.get(c, [])
        d_all = merge([(s, e) for (s, e, _g) in d])
        a_all = merge([(s, e) for (s, e, _g) in a])
        d_up = merge([(s, e) for (s, e, g) in d if g > 0])
        d_dn = merge([(s, e) for (s, e, g) in d if g < 0])
        a_up = merge([(s, e) for (s, e, g) in a if g > 0])
        a_dn = merge([(s, e) for (s, e, g) in a if g < 0])

        sig_bp_def += total_len(d_all)
        sig_bp_ada += total_len(a_all)
        oa = overlap_len(d_all, a_all)
        od = overlap_len(d_up, a_up) + overlap_len(d_dn, a_dn)
        ov_any += oa
        ov_dir += od
        union += total_len(d_all) + total_len(a_all) - oa
        reg_def += len(d_all)
        reg_ada += len(a_all)
        hit_def += n_regions_hit(d_all, a_all)
        hit_ada += n_regions_hit(a_all, d_all)

    bpj = ov_any / union if union else float("nan")
    bpj_dir = ov_dir / union if union else float("nan")
    dir_agree = ov_dir / ov_any if ov_any else float("nan")
    rec_def = hit_def / reg_def if reg_def else float("nan")
    rec_ada = hit_ada / reg_ada if reg_ada else float("nan")
    return dict(sig_bp_def=sig_bp_def, sig_bp_ada=sig_bp_ada, bpj=bpj,
                bpj_dir=bpj_dir, dir_agree=dir_agree, reg_def=reg_def,
                reg_ada=reg_ada, rec_def=rec_def, rec_ada=rec_ada)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--default-dir", required=True)
    ap.add_argument("--adaptive-dir", required=True)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--comparisons", default=None,
                    help="comma-separated comparison names (e.g. DN_vs_DP); "
                         "default = all diff_*.stats.tsv present in both dirs")
    ap.add_argument("--out", default=None, help="optional TSV path")
    args = ap.parse_args()

    def comps_in(d):
        return {os.path.basename(p)[len("diff_"):-len(".stats.tsv")]
                for p in glob.glob(os.path.join(d, "diff_*.stats.tsv"))}
    if args.comparisons:
        comps = args.comparisons.split(",")
    else:
        comps = sorted(comps_in(args.default_dir) & comps_in(args.adaptive_dir))
    if not comps:
        sys.exit("[ERROR] no shared diff_*.stats.tsv comparisons in both dirs")

    print("=" * 100)
    print("ADAPTIVE vs DEFAULT differential calls -- region overlap "
          "(shared FDR < %.3g, main chromosomes)" % args.fdr)
    print("-" * 100)
    hdr = ("comparison", "sigbp_def", "sigbp_ada", "bpJacc", "bpJacc_dir",
           "dir_agr", "reg_def", "reg_ada", "rec_d>a", "rec_a>d")
    print("%-14s %10s %10s %7s %10s %8s %8s %8s %8s %8s" % hdr)

    rows = []
    for comp in comps:
        dp = os.path.join(args.default_dir, "diff_%s.stats.tsv" % comp)
        ap_ = os.path.join(args.adaptive_dir, "diff_%s.stats.tsv" % comp)
        if not (os.path.exists(dp) and os.path.exists(ap_)):
            continue
        r = compare_one(load_sig(dp, args.fdr), load_sig(ap_, args.fdr))
        rows.append((comp, r))
        print("%-14s %10d %10d %7.3f %10.3f %8.3f %8d %8d %8.3f %8.3f"
              % (comp, r["sig_bp_def"], r["sig_bp_ada"], r["bpj"], r["bpj_dir"],
                 r["dir_agree"], r["reg_def"], r["reg_ada"], r["rec_def"], r["rec_ada"]))

    print("=" * 100)
    print("bpJacc = shared significant bp / union; bpJacc_dir requires same "
          "direction; dir_agr = same-direction fraction of shared bp; rec_d>a = "
          "fraction of default sig regions hit by an adaptive sig region.")

    if args.out:
        with open(args.out, "w") as fh:
            fh.write("comparison\tsig_bp_def\tsig_bp_ada\tbp_jaccard\tbp_jaccard_dir\t"
                     "dir_agree\tregions_def\tregions_ada\trecall_def_to_ada\t"
                     "recall_ada_to_def\n")
            for comp, r in rows:
                fh.write("%s\t%d\t%d\t%.5f\t%.5f\t%.5f\t%d\t%d\t%.5f\t%.5f\n"
                         % (comp, r["sig_bp_def"], r["sig_bp_ada"], r["bpj"],
                            r["bpj_dir"], r["dir_agree"], r["reg_def"], r["reg_ada"],
                            r["rec_def"], r["rec_ada"]))
        print("wrote %s" % args.out)


if __name__ == "__main__":
    main()
