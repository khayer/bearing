#!/usr/bin/env python3
"""
find_multimodal.py

Scan a BEARING differential .stats.tsv and describe, per bin, how the signal is
distributed across assay channels. For each bin:
  total = sum|kl_i|
  share p_i = |kl_i|/total
  eff  = 1/sum(p_i^2)          effective #channels (1 = one assay dominates;
                               ~n = spread evenly; HIGH eff often = diffuse/noise)
  coh  = |sum(kl_i)|/total     directional coherence (1 = all same sign =
                               coordinated; ~0 = channels cancel = incoherent)

A clean, biologically real coordinated change tends to be HIGH score + moderate
eff (~2-3.5, a few assays) + HIGH coh (same direction). A diffuse noise bin tends
to be HIGH eff (~n) + LOW-ish coh + low score. Default sort is by |score|.

Usage:
  python find_multimodal.py diff_DN_vs_DP.stats.tsv \
      --labels "ATAC,RNAseq+,RNAseq-,CTCF,Cohesin,H3K27ac" \
      --sort score --min-eff 2.0 --max-eff 3.5 --min-coh 0.6 --min-score 1.0 --top 25
  # add --sig to keep only FDR-significant bins
"""
import argparse, sys

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tsv")
    ap.add_argument("--labels", default="")
    ap.add_argument("--sort", choices=["score", "eff", "coh"], default="score")
    ap.add_argument("--min-score", type=float, default=1.0)
    ap.add_argument("--min-eff", type=float, default=0.0)
    ap.add_argument("--max-eff", type=float, default=1e9, help="cap eff to exclude diffuse/noise bins")
    ap.add_argument("--min-coh", type=float, default=0.0, help="min directional coherence 0..1")
    ap.add_argument("--sig", action="store_true")
    ap.add_argument("--padj", type=float, default=None)
    ap.add_argument("--top", type=int, default=25)
    args = ap.parse_args()

    labels = [s.strip() for s in args.labels.split(",")] if args.labels else None
    rows = []
    with open(args.tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        kl_idx = [i for i, name in enumerate(header) if name.lower().startswith("kl")]
        if not kl_idx:
            sys.exit("no kl_ columns in header")
        for n in ["chrom", "start", "end", "bearing_score", "pval_adj_bh"]:
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
            eff = 1.0 / sum((abs(x) / tot) ** 2 for x in kl)
            coh = abs(sum(kl)) / tot
            if eff < args.min_eff or eff > args.max_eff or coh < args.min_coh:
                continue
            rows.append((score, eff, coh, padj, f[col["chrom"]],
                         int(f[col["start"]]), int(f[col["end"]]), kl))

    keymap = {"score": lambda r: abs(r[0]), "eff": lambda r: r[1], "coh": lambda r: r[2]}
    rows.sort(key=keymap[args.sort], reverse=True)
    if not rows:
        print("No bins passed the filters. Loosen --min-score / --min-coh / widen eff, or drop --sig.")
        return
    n = len(kl_idx)
    names = labels if labels and len(labels) == n else ["kl_%d" % (i + 1) for i in range(n)]
    print("# %-22s %7s %5s %5s %10s   per-channel kl (signed, sorted by |kl|)" %
          ("locus", "score", "eff", "coh", "padj"))
    for score, eff, coh, padj, chrom, start, end, kl in rows[:args.top]:
        pairs = sorted(zip(names, kl), key=lambda p: abs(p[1]), reverse=True)
        brk = "  ".join("%s=%+.2f" % (nm, v) for nm, v in pairs if abs(v) > 0.05)
        print("%-24s %7.2f %5.2f %5.2f %10.1e   %s" %
              ("%s:%d-%d" % (chrom, start, end), score, eff, coh, padj, brk))

if __name__ == "__main__":
    main()
