#!/usr/bin/env python3
"""
region_match_ebko2.py -- EbKO recombination-centre characterisation.

Resolves dilution vs genuine insensitivity at the EbKO recombination centre, and
produces the numbers quoted in Results: the per-track share of the joint
differential across the RC (the "~36% RNA+" figure), the fraction of edgeR-
significant RC bins in which BEARING's RNA attribution moves in the DN direction,
and the Spearman correlation between edgeR logFC and BEARING RNA / joint scores.

WHAT CHANGED
------------
Previously this printed its results and read its inputs from hard-coded paths at
the repo root. That made it (a) unrunnable once the edgeR CSVs were moved into
paper/table_sources/, and (b) the least reproducible number in the paper: the
values reached the manuscript by reading stdout. It now takes paths on the CLI
(defaulting to the tracked locations) and writes a TSV with --out. The
computation is unchanged; the printed summary is unchanged.

USAGE
-----
    python3 dev/region_match_ebko2.py \
        --out paper/table_sources/region_match_ebko_rc.tsv

Defaults assume the repo root as CWD and inputs in paper/table_sources/.
ASCII only.
"""

import argparse
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

RNA_POS, RNA_NEG = "kl_2", "kl_3"   # from mm10_6track_panel.yaml
TRACKS = {"kl_1": "ATAC", "kl_2": "RNA+", "kl_3": "RNA-",
          "kl_4": "CTCF", "kl_5": "RAD21", "kl_6": "H3K27ac"}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--diff",
                    default="workflow/results/pvalue/diff_DN_vs_EbKO.stats.tsv",
                    help="DN-vs-EbKO per-bin BEARING stats table.")
    ap.add_argument("--edger-pos",
                    default="paper/table_sources/rna_concordance_DN_vs_EbKO_pos_edgeR_allbins.csv")
    ap.add_argument("--edger-neg",
                    default="paper/table_sources/rna_concordance_DN_vs_EbKO_neg_edgeR_allbins.csv")
    ap.add_argument("--chrom", default="chr6")
    ap.add_argument("--lo", type=int, default=41500000)
    ap.add_argument("--hi", type=int, default=41567200)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--out", default=None,
                    help="Write the characterisation to this TSV. If omitted, "
                         "results are printed only (legacy behaviour).")
    a = ap.parse_args()

    for path in (a.diff, a.edger_pos, a.edger_neg):
        import os
        if not os.path.exists(path):
            sys.exit("ERROR: input not found: %s\n"
                     "  --diff comes from the pipeline (rule pvalue).\n"
                     "  the edgeR CSVs come from rna_concordance_stranded.R.\n"
                     "  Pass explicit paths if yours are elsewhere." % path)

    # pull RC neighbourhood from the big BEARING file
    keep = []
    for ch in pd.read_csv(a.diff, sep="\t", chunksize=500000, low_memory=False):
        m = (ch["chrom"] == a.chrom) & (ch["start"] >= a.lo) & (ch["end"] <= a.hi)
        if m.any():
            keep.append(ch.loc[m])
    if not keep:
        sys.exit("ERROR: no BEARING bins in %s:%d-%d" % (a.chrom, a.lo, a.hi))
    b = pd.concat(keep, ignore_index=True)
    print("RC BEARING bins:", len(b))

    out_rows = []  # tidy long-form: section, metric, strand, value

    # track composition across the whole RC window (mean absolute differential KL)
    comp = {TRACKS[k]: float(b[k].abs().mean()) for k in TRACKS}
    tot = sum(comp.values()) or 1.0
    print("\nMean |dKL| per track across RC (share of joint):")
    for name, v in sorted(comp.items(), key=lambda x: -x[1]):
        share = 100 * v / tot
        print("  %-7s %.4f  (%4.1f%%)" % (name, v, share))
        out_rows.append(("track_composition", "mean_abs_dKL", name, "%.6f" % v))
        out_rows.append(("track_composition", "share_pct", name, "%.2f" % share))
    out_rows.append(("track_composition", "n_rc_bins", "", str(len(b))))

    for strand, egf, rnacol in (("pos", a.edger_pos, RNA_POS),
                                ("neg", a.edger_neg, RNA_NEG)):
        tname = TRACKS[rnacol]
        print("\n=============== strand %s  (BEARING %s=%s) ==============="
              % (strand, rnacol, tname))
        eg = pd.read_csv(egf)
        eg = eg[(eg["chr"] == a.chrom) & (eg["start"] >= a.lo) & (eg["end"] <= a.hi)]
        m = b.merge(eg[["start", "logFC", "FDR"]], on="start", how="inner")
        n_sig = int((m["FDR"] <= a.fdr).sum())
        print("edgeR-tested bins in RC:", len(m), "| edgeR sig:", n_sig)
        out_rows.append(("edger_match_%s" % strand, "edger_tested_bins", strand, str(len(m))))
        out_rows.append(("edger_match_%s" % strand, "edger_sig_bins", strand, str(n_sig)))
        if len(m) == 0:
            continue
        sig = m[m["FDR"] <= a.fdr]
        pos_frac = float((sig[rnacol] > 0).mean()) if len(sig) else float("nan")
        print("At edgeR-sig bins (n=%d):" % len(sig))
        print("  BEARING %s dKL > 0 (DN-direction): %.0f%% of bins"
              % (tname, 100 * pos_frac))
        print("  median %s dKL = %.4f (range %.3f .. %.3f)"
              % (tname, sig[rnacol].median(), sig[rnacol].min(), sig[rnacol].max()))
        print("  median joint bearing_score = %.4f" % sig["bearing_score"].median())
        n_raw = int((sig["pval"] <= 0.05).sum())
        print("  joint pval: min=%.2e median=%.2e (none cleared BH; %d/%d have raw pval<=0.05)"
              % (sig["pval"].min(), sig["pval"].median(), n_raw, len(sig)))
        rho, p = spearmanr(m["logFC"], m[rnacol])
        rho2, p2 = spearmanr(m["logFC"], m["bearing_score"])
        print("  Spearman edgeR logFC vs BEARING %s dKL: rho=%.2f (p=%.1e)" % (tname, rho, p))
        print("  Spearman edgeR logFC vs BEARING joint score:        rho=%.2f (p=%.1e)" % (rho2, p2))

        sec = "edger_match_%s" % strand
        out_rows += [
            (sec, "dn_direction_pct", strand, "%.1f" % (100 * pos_frac)),
            (sec, "median_rna_dKL", strand, "%.6f" % sig[rnacol].median()),
            (sec, "median_joint_score", strand, "%.6f" % sig["bearing_score"].median()),
            (sec, "n_raw_pval_le_0.05", strand, str(n_raw)),
            (sec, "n_sig", strand, str(len(sig))),
            (sec, "spearman_logFC_vs_rna_rho", strand, "%.4f" % rho),
            (sec, "spearman_logFC_vs_rna_p", strand, "%.3e" % p),
            (sec, "spearman_logFC_vs_joint_rho", strand, "%.4f" % rho2),
            (sec, "spearman_logFC_vs_joint_p", strand, "%.3e" % p2),
        ]

    if a.out:
        import os
        os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
        with open(a.out, "w") as fh:
            fh.write("# region_match_ebko2.py -- EbKO RC characterisation\n")
            fh.write("# diff=%s\n# edger_pos=%s\n# edger_neg=%s\n"
                     % (a.diff, a.edger_pos, a.edger_neg))
            fh.write("# window=%s:%d-%d fdr=%s\n" % (a.chrom, a.lo, a.hi, a.fdr))
            fh.write("section\tmetric\tstrand\tvalue\n")
            for r in out_rows:
                fh.write("\t".join(r) + "\n")
        print("\nwrote %s (%d rows)" % (a.out, len(out_rows)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
