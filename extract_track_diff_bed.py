#!/usr/bin/env python3
"""
extract_track_diff_bed.py -- pull one track's differential bins (e.g. CTCF) out
of a BEARING diff *_stats.tsv into a directional BED, for bin-to-bin concordance
(the bearing_sig_ctcf_bed input to the CTCF edgeR check).

BEARING computes a per-bin TOTAL p-value, not a per-track one, and a bin can be
significant overall while contributing zero on a given track. So a track-specific
BED is defined from the per-track signed contribution kl_<i> (i = the track's
category index): a CTCF-differential bin is one where the CTCF contribution kl_4
is non-zero and large / the bin is significant. Sign of kl_<i> gives direction
(> 0 = first-condition / A-enriched, e.g. DN; < 0 = second-condition / B-enriched,
e.g. ProB) under BEARING's A-B convention.

Two ways to select (combine freely; the file is large, so always filter):
  --significant            keep bins flagged significant_fdr0.05 == 1
  --fdr 0.05               keep bins with pval_adj_bh < 0.05
  --top-percent 5          keep the top 5% (A-enriched) and bottom 5%
                           (B-enriched) bins by the track's kl value
  --min-abs 1.0            also require |kl_track| >= 1.0
If no selection flag is given, defaults to --top-percent 5. When --top-percent is
combined with a significance filter, the percentile is taken over the significant
bins. Bins with kl_track == 0 are never emitted (no CTCF signal there).

Output BED6: chrom  start  end  <track>_up|<track>_dn  <signed kl>  .
  - column 4 (name): direction label
  - column 5 (score): the signed kl_track value (sign = direction, |.| = strength)
  - column 6 (strand): "." (direction is in name/score, not strand, to avoid
    strand-aware overlap surprises downstream)

Example (CTCF, significant bins that actually have CTCF signal):
  python extract_track_diff_bed.py \
    --diff-tsv diff_DN_vs_ProB_chr_subset_stats.tsv \
    --categories ProB_rep1_cats.json --track CTCF \
    --significant --min-abs 0.5 \
    --out bearing_DNvsProB_CTCF_significant_bins.bed

Example (top 5% up/down by CTCF, significance-agnostic):
  python extract_track_diff_bed.py \
    --diff-tsv diff_DN_vs_ProB_chr_subset_stats.tsv \
    --categories ProB_rep1_cats.json --track CTCF \
    --top-percent 5 \
    --out bearing_DNvsProB_CTCF_top5pct.bed

ASCII only.
"""
import argparse
import json
import sys

import numpy as np

# fixed column layout of the diff stats TSV
COL = {"chrom": 0, "start": 1, "end": 2, "pval_adj_bh": 8,
       "significant": 9, "direction": 10}
KL_BASE = 10   # kl_<i> is at column KL_BASE + i (kl_1 -> 11, kl_4 -> 14)


def track_index(categories_json, track):
    """Map a track name to its category index using the cats JSON."""
    with open(categories_json) as fh:
        cats = json.load(fh)["categories"]
    want = track.strip().lower().replace(" ", "")
    for idx, val in cats.items():
        name = val[0] if isinstance(val, (list, tuple)) else val
        if name.strip().lower().replace(" ", "") == want:
            return int(idx), name
    avail = ", ".join("%s=%s" % (k, (v[0] if isinstance(v, (list, tuple)) else v))
                      for k, v in cats.items())
    sys.exit("[ERROR] track '%s' not found in %s. Available: %s"
             % (track, categories_json, avail))


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--diff-tsv", required=True)
    ap.add_argument("--categories", required=True, help="*_cats.json for this run")
    ap.add_argument("--track", required=True, help="track name, e.g. CTCF")
    ap.add_argument("--significant", action="store_true",
                    help="keep significant_fdr0.05 == 1")
    ap.add_argument("--fdr", type=float, default=None,
                    help="keep pval_adj_bh < this (overrides --significant)")
    ap.add_argument("--top-percent", type=float, default=None,
                    help="keep top/bottom this %% by kl_track")
    ap.add_argument("--min-abs", type=float, default=0.0,
                    help="require |kl_track| >= this (default 0 = any non-zero)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    idx, track_name = track_index(args.categories, args.track)
    kl_col = KL_BASE + idx
    use_top = args.top_percent
    if args.fdr is None and not args.significant and use_top is None:
        use_top = 5.0   # default selection

    def sig_pass(f):
        if args.fdr is not None:
            try:
                return float(f[COL["pval_adj_bh"]]) < args.fdr
            except ValueError:
                return False
        if args.significant:
            return f[COL["significant"]] == "1"
        return True

    # pass 1 (only if a percentile cutoff is needed): collect kl over the
    # track-ACTIVE (kl != 0) bins in the significance-filtered universe, so
    # "top N%" means the strongest N% of bins that actually have this track's
    # signal -- not N% of the whole genome (mostly zero on any one track).
    hi = lo = None
    if use_top is not None:
        vals = []
        with open(args.diff_tsv) as fh:
            fh.readline()
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) <= kl_col or not sig_pass(f):
                    continue
                try:
                    v = float(f[kl_col])
                except ValueError:
                    continue
                if v != 0.0:
                    vals.append(v)
        if not vals:
            sys.exit("[ERROR] no track-active bins pass the filter for percentiles")
        vals = np.asarray(vals, dtype=np.float64)
        hi = float(np.percentile(vals, 100.0 - use_top))
        lo = float(np.percentile(vals, use_top))

    n_up = n_dn = 0
    with open(args.diff_tsv) as fh, open(args.out, "w") as out:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= kl_col or not sig_pass(f):
                continue
            try:
                kl = float(f[kl_col])
            except ValueError:
                continue
            if kl == 0.0 or abs(kl) < args.min_abs:
                continue
            if use_top is not None:
                if not (kl >= hi or kl <= lo):
                    continue
            up = kl > 0
            name = "%s_%s" % (track_name.replace(" ", ""), "up" if up else "dn")
            out.write("%s\t%s\t%s\t%s\t%.4f\t.\n"
                      % (f[COL["chrom"]], f[COL["start"]], f[COL["end"]], name, kl))
            n_up += int(up)
            n_dn += int(not up)

    msg = ("[done] %s: %d bins (%d up / A-enriched, %d dn / B-enriched) -> %s"
           % (track_name, n_up + n_dn, n_up, n_dn, args.out))
    if use_top is not None:
        msg += " | top-%.1f%% cutoffs: up>=%.3f dn<=%.3f" % (use_top, hi, lo)
    if args.fdr is not None:
        msg += " | pval_adj_bh<%.3g" % args.fdr
    elif args.significant:
        msg += " | significant_fdr0.05"
    if args.min_abs > 0:
        msg += " | |kl|>=%.3g" % args.min_abs
    sys.stderr.write(msg + "\n")


if __name__ == "__main__":
    main()
