#!/usr/bin/env python3
"""
region_match_ebko.py -- EbKO recombination-centre characterisation (canonical).

Produces the EbKO RC numbers quoted in Results:
  - per-track share of the joint differential across the RC (the ~36% RNA+ figure)
  - fraction of edgeR-significant RC bins where BEARING's RNA attribution is
    DN-directed
  - Spearman correlation of edgeR logFC vs BEARING RNA dKL and vs joint score

PROVENANCE / WHY THIS SUPERSEDES ebko2 AND ebko3
------------------------------------------------
This merges the two earlier prototypes and drops both:

  region_match_ebko2.py had the track-composition block (the ~36% number) but
  joined BEARING to edgeR on the RAW start coordinate. edgeR/featureCounts output
  is 1-based (start=4807801) while BEARING bins are 0-based (start=4807800), so
  `merge(on="start")` matched NOTHING: "edgeR-tested bins in RC: 0". Every
  edgeR-concordance number it produced was empty.

  region_match_ebko3.py fixed the join -- it maps both sides to a 200 bp grid via
  (edgeR_start - 1)//BIN, absorbing the 1-based offset -- but had DROPPED the
  track-composition block, so it could not produce the ~36% figure.

Neither alone is correct-and-complete. This script uses ebko3's grid join AND
ebko2's composition block. The offset was confirmed against the real files:
edgeR start 4807801 and BEARING start 4807800 are the same 200 bp bin.

The computation is otherwise unchanged from the two originals; only the join and
the CLI/output plumbing are new.

USAGE
-----
    python3 dev/region_match_ebko.py \
        --out paper/table_sources/region_match_ebko_rc.tsv

Defaults assume the repo root as CWD and the edgeR CSVs in paper/table_sources/.
ASCII only.
"""

import argparse
import os
import sys

import pandas as pd
from scipy.stats import spearmanr

BIN = 200
RNA_POS, RNA_NEG = "kl_2", "kl_3"   # from mm10_6track_panel.yaml
TRACKS = {"kl_1": "ATAC", "kl_2": "RNA+", "kl_3": "RNA-",
          "kl_4": "CTCF", "kl_5": "RAD21", "kl_6": "H3K27ac"}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--diff",
                    default="workflow/results/pvalue/diff_DN_vs_EbKO.stats.tsv")
    ap.add_argument("--edger-pos",
                    default="paper/table_sources/rna_concordance_DN_vs_EbKO_pos_edgeR_allbins.csv")
    ap.add_argument("--edger-neg",
                    default="paper/table_sources/rna_concordance_DN_vs_EbKO_neg_edgeR_allbins.csv")
    ap.add_argument("--chrom", default="chr6")
    ap.add_argument("--lo", type=int, default=41500000)
    ap.add_argument("--hi", type=int, default=41567200)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--out", default=None,
                    help="Write the characterisation to this TSV (recommended).")
    a = ap.parse_args()

    for path in (a.diff, a.edger_pos, a.edger_neg):
        if not os.path.exists(path):
            sys.exit("ERROR: input not found: %s\n"
                     "  --diff comes from the pipeline (rule pvalue).\n"
                     "  the edgeR CSVs come from rna_concordance_stranded.R.\n"
                     "  Pass explicit paths if yours are elsewhere." % path)

    # RC neighbourhood from the big BEARING differential
    keep = []
    for ch in pd.read_csv(a.diff, sep="\t", chunksize=500000, low_memory=False):
        m = (ch["chrom"] == a.chrom) & (ch["start"] >= a.lo) & (ch["end"] <= a.hi)
        if m.any():
            keep.append(ch.loc[m])
    if not keep:
        sys.exit("ERROR: no BEARING bins in %s:%d-%d" % (a.chrom, a.lo, a.hi))
    b = pd.concat(keep, ignore_index=True)
    b["start"] = b["start"].astype(int)
    b["bid"] = b["start"] // BIN
    print("RC BEARING bins:", len(b))

    out_rows = []

    # --- track composition across the RC (from ebko2) ----------------------
    comp = {TRACKS[k]: float(b[k].abs().mean()) for k in TRACKS}
    tot = sum(comp.values()) or 1.0
    print("\nMean |dKL| per track across RC (share of joint):")
    for name, v in sorted(comp.items(), key=lambda x: -x[1]):
        share = 100 * v / tot
        print("  %-7s %.4f  (%4.1f%%)" % (name, v, share))
        out_rows.append(("track_composition", "mean_abs_dKL", name, "%.6f" % v))
        out_rows.append(("track_composition", "share_pct", name, "%.2f" % share))
    out_rows.append(("track_composition", "n_rc_bins", "", str(len(b))))

    # --- edgeR concordance with the CORRECT grid join (from ebko3) ---------
    for strand, egf, rnacol in (("pos", a.edger_pos, RNA_POS),
                                ("neg", a.edger_neg, RNA_NEG)):
        tname = TRACKS[rnacol]
        print("\n=============== strand %s  (BEARING %s=%s) ==============="
              % (strand, rnacol, tname))
        eg = pd.read_csv(egf)
        eg = eg[(eg["chr"] == a.chrom) & (eg["start"] >= a.lo - BIN)
                & (eg["end"] <= a.hi + BIN)].copy()
        eg["start"] = eg["start"].astype(int)
        # map edgeR (1-based) start to the 0-based 200 bp grid
        eg["bid"] = (eg["start"] - 1) // BIN
        n_grid = int(b["bid"].isin(set(eg["bid"])).sum())
        print("BEARING bins whose bin-id matches an edgeR bin-id:", n_grid)
        m = b.merge(eg[["bid", "logFC", "FDR"]], on="bid", how="inner")
        n_sig = int((m["FDR"] <= a.fdr).sum())
        print("merged bins:", len(m), "| edgeR sig (FDR<=%.2f): %d" % (a.fdr, n_sig))
        sec = "edger_match_%s" % strand
        out_rows.append((sec, "grid_matched_bins", strand, str(n_grid)))
        out_rows.append((sec, "merged_bins", strand, str(len(m))))
        out_rows.append((sec, "edger_sig_bins", strand, str(n_sig)))
        if len(m) == 0:
            print("  EMPTY MERGE -- offset not absorbed; do not trust this strand.")
            continue

        sig = m[m["FDR"] <= a.fdr]
        if len(sig):
            pos_frac = float((sig[rnacol] > 0).mean())
            n_raw = int((sig["pval"] <= 0.05).sum())
            print("At edgeR-sig bins (n=%d):" % len(sig))
            print("  %s dKL > 0 (DN-direction): %.0f%%" % (tname, 100 * pos_frac))
            print("  median %s dKL = %.4f (range %.3f..%.3f)"
                  % (tname, sig[rnacol].median(), sig[rnacol].min(), sig[rnacol].max()))
            print("  median joint bearing_score = %.4f" % sig["bearing_score"].median())
            print("  raw joint pval<=0.05: %d/%d; min joint pval=%.2e"
                  % (n_raw, len(sig), sig["pval"].min()))
            out_rows += [
                (sec, "dn_direction_pct", strand, "%.1f" % (100 * pos_frac)),
                (sec, "median_rna_dKL", strand, "%.6f" % sig[rnacol].median()),
                (sec, "median_joint_score", strand, "%.6f" % sig["bearing_score"].median()),
                (sec, "n_raw_pval_le_0.05", strand, str(n_raw)),
                (sec, "n_sig", strand, str(len(sig))),
            ]
        rho, p = spearmanr(m["logFC"], m[rnacol])
        rho2, p2 = spearmanr(m["logFC"], m["bearing_score"])
        print("  Spearman edgeR logFC vs BEARING %s dKL: rho=%.2f (p=%.1e)" % (tname, rho, p))
        print("  Spearman edgeR logFC vs BEARING joint score:        rho=%.2f (p=%.1e)" % (rho2, p2))
        out_rows += [
            (sec, "spearman_logFC_vs_rna_rho", strand, "%.4f" % rho),
            (sec, "spearman_logFC_vs_rna_p", strand, "%.3e" % p),
            (sec, "spearman_logFC_vs_joint_rho", strand, "%.4f" % rho2),
            (sec, "spearman_logFC_vs_joint_p", strand, "%.3e" % p2),
        ]

    if a.out:
        os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
        with open(a.out, "w") as fh:
            fh.write("# region_match_ebko.py -- EbKO RC characterisation\n")
            fh.write("# diff=%s\n# edger_pos=%s\n# edger_neg=%s\n"
                     % (a.diff, a.edger_pos, a.edger_neg))
            fh.write("# window=%s:%d-%d fdr=%s bin=%d\n"
                     % (a.chrom, a.lo, a.hi, a.fdr, BIN))
            fh.write("# join: edgeR (start-1)//%d vs BEARING start//%d "
                     "(absorbs 1-based vs 0-based offset)\n" % (BIN, BIN))
            fh.write("section\tmetric\tstrand\tvalue\n")
            for r in out_rows:
                fh.write("\t".join(r) + "\n")
        print("\nwrote %s (%d rows)" % (a.out, len(out_rows)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
