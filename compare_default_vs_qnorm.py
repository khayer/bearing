#!/usr/bin/env python3
"""Compare the default (un-normalized) BEARING run against the cohort-wide
quantile-normalized run, to test whether results -- especially the differential
calls -- are robust to normalization.

Three reports:

  1. DIFF CONCORDANCE (per comparison): for each diff_<cmp>.stats.tsv present in
     both runs, FDR-significant bin counts in each run, the Jaccard overlap and
     directional agreement of the two FDR-significant bin SETS, and Spearman
     correlation of the per-bin BEARING scores over shared bins. High Jaccard +
     high score correlation => the differential layer is normalization-robust.

  2. SIGNIFICANT-BIN COUNT CONCORDANCE: a compact table of FDR-sig counts
     default vs qnorm per comparison, plus the ratio.

  3. PER-SAMPLE NULL CALIBRATION before/after: for the samples whose per-sample
     null was inflated under the default run (e.g. ProB_rep1, S3T3_rep*), the
     null-survival at matched score thresholds and the FDR-significant bin count
     in each run. If cohort normalization works, the inflated tails collapse and
     FDR-significant bins appear.

Reads only the *.stats.tsv tables (chrom,start,end,bearing_score,pval,
pval_adj_bh,[direction,]significant_fdr*). No third-party deps beyond numpy.

ASCII-only.
"""
import argparse
import glob
import os
import sys
import numpy as np


def _norm_chrom(c):
    return c[3:] if c.startswith("chr") else c


def parse_region(s):
    """Parse a region string into (chrom, start, end). Accepts 'chr6:40000000-
    42000000' (windowed) or 'chr6' (whole chromosome -> start/end None). Returns
    None for None/empty input (genome-wide)."""
    if not s:
        return None
    if ":" in s:
        chrom, rng = s.split(":", 1)
        a, b = rng.replace(",", "").split("-")
        return (_norm_chrom(chrom), int(a), int(b))
    return (_norm_chrom(s), None, None)


def _in_region(chrom, start, region):
    """True if a bin (chrom, start) falls in region (chrom, lo, hi)."""
    if region is None:
        return True
    if chrom != region[0]:
        return False
    if region[1] is not None and not (region[1] <= start <= region[2]):
        return False
    return True


def load_stats(path, region=None):
    """Return dict with per-bin arrays keyed by (chrom,start). Keeps score,
    pval_adj, significance, and direction sign if present."""
    keys, score, padj, sig, sign = [], [], [], [], []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}
        ci, si, ei = idx.get("chrom"), idx.get("start"), idx.get("end")
        sc_i = idx.get("bearing_score")
        padj_i = idx.get("pval_adj_bh")
        dir_i = idx.get("direction")
        sig_col = next((h for h in header if h.lower().startswith("significant_fdr")), None)
        sig_i = idx.get(sig_col) if sig_col else None
        if None in (ci, si, sc_i):
            return None
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(ci, si, sc_i):
                continue
            try:
                key = (_norm_chrom(f[ci]), int(f[si]))
                sv = float(f[sc_i])
            except (ValueError, IndexError):
                continue
            if not _in_region(key[0], key[1], region):
                continue
            keys.append(key)
            score.append(sv)
            padj.append(_safe(f, padj_i))
            if sig_i is not None and sig_i < len(f):
                sig.append(f[sig_i].strip().lower() in ("1", "true", "yes", "t"))
            else:
                sig.append(False)
            if dir_i is not None and dir_i < len(f):
                sign.append(1 if f[dir_i].strip() in ("+", "1", "A") else -1)
            else:
                sign.append(1 if sv >= 0 else -1)
    return {
        "key": keys,
        "score": np.array(score, dtype=float),
        "padj": np.array(padj, dtype=float),
        "sig": np.array(sig, dtype=bool),
        "sign": np.array(sign, dtype=int),
    }


def _safe(f, i):
    if i is None or i >= len(f):
        return np.nan
    try:
        return float(f[i])
    except ValueError:
        return np.nan


def spearman(a, b):
    """Spearman rho without scipy: Pearson on ranks."""
    if a.size < 3:
        return float("nan")
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean(); rb -= rb.mean()
    denom = np.sqrt((ra * ra).sum() * (rb * rb).sum())
    return float((ra * rb).sum() / denom) if denom > 0 else float("nan")


def sig_keyset(d, signed=False):
    if signed:
        return set((k, int(s)) for k, sg, s in zip(d["key"], d["sig"], d["sign"]) if sg)
    return set(k for k, sg in zip(d["key"], d["sig"]) if sg)


def jaccard(a, b):
    if not a and not b:
        return 1.0
    u = len(a | b)
    return (len(a & b) / u) if u else 0.0


def diff_concordance(default_dir, qnorm_dir, region=None):
    print("=" * 78)
    print("1. DIFF CONCORDANCE (per comparison): default vs cohort-qnorm")
    print("-" * 78)
    files = sorted(glob.glob(os.path.join(default_dir, "diff_*.stats.tsv")))
    cmps = [os.path.basename(p)[len("diff_"):-len(".stats.tsv")] for p in files]
    hdr = ("%-16s %8s %8s %8s %8s %8s %8s" %
           ("comparison", "sig_def", "sig_qn", "jacc", "jacc_dir", "rho_sc", "n_shared"))
    print(hdr)
    rows = []
    for cmp in cmps:
        pd = os.path.join(default_dir, "diff_%s.stats.tsv" % cmp)
        pq = os.path.join(qnorm_dir, "diff_%s.stats.tsv" % cmp)
        if not os.path.exists(pq):
            print("%-16s   (no qnorm file yet)" % cmp)
            continue
        dd, dq = load_stats(pd, region=region), load_stats(pq, region=region)
        if dd is None or dq is None:
            print("%-16s   (unreadable)" % cmp)
            continue
        sd, sq = sig_keyset(dd), sig_keyset(dq)
        sdd, sqd = sig_keyset(dd, signed=True), sig_keyset(dq, signed=True)
        j = jaccard(sd, sq)
        jd = jaccard(sdd, sqd)
        md = {k: v for k, v in zip(dd["key"], dd["score"])}
        shared = [(md[k], v) for k, v in zip(dq["key"], dq["score"]) if k in md]
        if shared:
            a = np.array([x for x, _ in shared]); b = np.array([y for _, y in shared])
            rho = spearman(a, b)
        else:
            rho = float("nan")
        n_sd, n_sq = int(dd["sig"].sum()), int(dq["sig"].sum())
        print("%-16s %8d %8d %8.3f %8.3f %8.3f %8d"
              % (cmp, n_sd, n_sq, j, jd, rho, len(shared)))
        rows.append((cmp, n_sd, n_sq, j, jd, rho))
    print("")
    print("jacc = overlap of FDR-significant bin sets; jacc_dir = overlap when")
    print("also requiring the same direction; rho_sc = Spearman of per-bin scores")
    print("over shared bins. High jacc + high rho_sc => differential calls are")
    print("normalization-robust. Low jacc with high rho_sc => calls reshuffle near")
    print("the threshold but the underlying score landscape is preserved.")
    return rows


def count_concordance(rows):
    print("=" * 78)
    print("2. SIGNIFICANT-BIN COUNT CONCORDANCE")
    print("-" * 78)
    print("%-16s %10s %10s %8s" % ("comparison", "default", "qnorm", "qn/def"))
    for cmp, n_sd, n_sq, _j, _jd, _rho in rows:
        ratio = (n_sq / n_sd) if n_sd else float("inf") if n_sq else float("nan")
        print("%-16s %10d %10d %8.2f" % (cmp, n_sd, n_sq, ratio))


def null_survival_one(d, thresholds):
    score, padj = d["score"], d["padj"]
    order = np.argsort(score)
    s_sorted, p_sorted = score[order], padj[order]
    out = []
    for t in thresholds:
        j = np.searchsorted(s_sorted, t, side="left")
        if j >= len(s_sorted):
            out.append(float("nan")); continue
        tail = p_sorted[j:]; tail = tail[np.isfinite(tail)]
        out.append(float(np.nanmax(tail)) if tail.size else float("nan"))
    return out


def persample_calibration(default_dir, qnorm_dir, samples, thresholds):
    print("=" * 78)
    print("3. PER-SAMPLE NULL CALIBRATION before/after (BH q-survival at score)")
    print("-" * 78)
    print("(uses pval_adj_bh as the tail proxy; lower-after = sharper null)")
    head = "%-14s %-7s" % ("sample", "run") + "".join("  q>=%-4.1f" % t for t in thresholds) + "   FDR_sig"
    print(head)
    for s in samples:
        for label, d in (("default", default_dir), ("qnorm", qnorm_dir)):
            p = os.path.join(d, "%s.stats.tsv" % s)
            if not os.path.exists(p):
                print("%-14s %-7s   (missing)" % (s, label)); continue
            dd = load_stats(p)
            if dd is None:
                print("%-14s %-7s   (unreadable)" % (s, label)); continue
            surv = null_survival_one(dd, thresholds)
            nsig = int(dd["sig"].sum())
            print("%-14s %-7s" % (s, label)
                  + "".join("  %8.1e" % v for v in surv)
                  + "  %8d" % nsig)
        print("")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--default-dir", required=True,
                    help="results/pvalue from the default (un-normalized) run.")
    ap.add_argument("--qnorm-dir", required=True,
                    help="results_qnorm/pvalue from the cohort-qnorm run.")
    ap.add_argument("--samples", nargs="+",
                    default=["ProB_rep1", "S3T3_rep1", "S3T3_rep2", "ProB_rep2", "DN_rep1"],
                    help="Per-sample tracks to check for null-calibration change.")
    ap.add_argument("--region", default=None,
                    help="Restrict the diff concordance to a region: 'chr6:40000000-"
                         "42000000' or a whole chromosome 'chr6'. Omit for genome-wide. "
                         "Per-sample calibration stays genome-wide.")
    args = ap.parse_args()

    region = parse_region(args.region)
    if region is not None:
        rtxt = region[0] if region[1] is None else "%s:%d-%d" % region
        print("[scope] diff concordance restricted to %s\n" % rtxt)

    if not os.path.isdir(args.qnorm_dir):
        print("[ERROR] qnorm dir not found: %s (has the qnorm run finished?)"
              % args.qnorm_dir, file=sys.stderr)
        sys.exit(1)

    rows = diff_concordance(args.default_dir, args.qnorm_dir, region=region)
    if rows:
        count_concordance(rows)
    persample_calibration(args.default_dir, args.qnorm_dir, args.samples,
                          [2.0, 2.5, 3.0, 3.5])


if __name__ == "__main__":
    main()
