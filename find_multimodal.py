#!/usr/bin/env python3
"""
find_multimodal.py

Scan a BEARING differential .stats.tsv and surface bins whose signal is spread
across SEVERAL assay channels (multi-modal) rather than dominated by one. These
are the "you couldn't see this in any single track" examples - the case that
justifies an integrated score over eyeballing one assay.

Per bin: total = sum|kl_i|; share p_i = |kl_i|/total; effective #channels =
1 / sum(p_i^2) (inverse-Simpson). ~1 = one assay dominates; >=2.5 = genuinely
multi-assay. Ranks the strongest bins by effective #channels.

Usage:
  python find_multimodal.py diff_DN_vs_DP.stats.tsv \
      --labels "ATAC,RNAseq+,RNAseq-,CTCF,Cohesin,H3K27ac" \
      --min-score 1.0 --min-eff 2.0 --top 25
  # add --sig to keep only FDR-significant bins
"""
import argparse, sys

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tsv")
    ap.add_argument("--labels", default="", help="comma-separated channel names in panel order (for kl_1..kl_n)")
    ap.add_argument("--min-score", type=float, default=1.0, help="min |bearing_score| to consider")
    ap.add_argument("--min-eff", type=float, default=2.0, help="min effective #channels")
    ap.add_argument("--sig", action="store_true", help="only significant_fdr0.05 == 1 bins")
    ap.add_argument("--padj", type=float, default=None, help="optional max pval_adj_bh")
    ap.add_argument("--top", type=int, default=25)
    args = ap.parse_args()

    labels = [s.strip() for s in args.labels.split(",")] if args.labels else None
    rows = []
    with open(args.tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        kl_idx = [i for i, name in enumerate(header) if name.lower().startswith("kl")]
        if not kl_idx:
            sys.exit("no kl_ columns found in header: %s" % header)
        need = ["chrom", "start", "end", "bearing_score", "pval_adj_bh"]
        for n in need:
            if n not in col:
                sys.exit("missing column: %s" % n)
        sig_i = col.get("significant_fdr0.05")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(kl_idx):
                continue
            try:
                score = float(f[col["bearing_score"]])
                padj = float(f[col["pval_adj_bh"]])
                kl = [float(f[i]) if f[i] not in ("", "NA") else 0.0 for i in kl_idx]
            except ValueError:
                continue
            if abs(score) < args.min_score:
                continue
            if args.sig and sig_i is not None and f[sig_i].strip() not in ("1", "True", "true"):
                continue
            if args.padj is not None and padj > args.padj:
                continue
            tot = sum(abs(x) for x in kl)
            if tot <= 0:
                continue
            shares = [abs(x) / tot for x in kl]
            eff = 1.0 / sum(s * s for s in shares)
            if eff < args.min_eff:
                continue
            rows.append((eff, abs(score), f[col["chrom"]], int(f[col["start"]]),
                         int(f[col["end"]]), score, padj, kl))

    rows.sort(key=lambda r: (r[0], r[1]), reverse=True)
    if not rows:
        print("No multi-modal bins passed the filters. Loosen --min-score / --min-eff, or drop --sig.")
        return
    n = len(kl_idx)
    names = labels if labels and len(labels) == n else ["kl_%d" % (i + 1) for i in range(n)]
    print("# top %d multi-modal bins (eff = effective #channels of %d)" % (min(args.top, len(rows)), n))
    print("# %-22s %7s %6s %10s   per-channel |kl| (sorted)" % ("locus", "score", "eff", "padj"))
    for eff, _, chrom, start, end, score, padj, kl in rows[:args.top]:
        pairs = sorted(zip(names, kl), key=lambda p: abs(p[1]), reverse=True)
        brk = "  ".join("%s=%+.2f" % (nm, v) for nm, v in pairs if abs(v) > 0.05)
        print("%-24s %7.2f %6.2f %10.1e   %s" % ("%s:%d-%d" % (chrom, start, end),
                                                 score, eff, padj, brk))

if __name__ == "__main__":
    main()
