#!/usr/bin/env python3
# region_match_ebko2.py  -- run from the bearing repo root
# Resolves dilution vs genuine insensitivity at the EbKO recombination center.
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

DIFF = "workflow/results/pvalue/diff_DN_vs_EbKO.stats.tsv"
EDGER_POS = "rna_concordance_DN_vs_EbKO_pos_edgeR_allbins.csv"
EDGER_NEG = "rna_concordance_DN_vs_EbKO_neg_edgeR_allbins.csv"
CHROM, LO, HI = "chr6", 41500000, 41567200
FDR = 0.05
RNA_POS, RNA_NEG = "kl_2", "kl_3"   # from mm10_6track_panel.yaml
TRACKS = {"kl_1": "ATAC", "kl_2": "RNA+", "kl_3": "RNA-",
          "kl_4": "CTCF", "kl_5": "RAD21", "kl_6": "H3K27ac"}

# pull RC neighbourhood from the big BEARING file
keep = []
for ch in pd.read_csv(DIFF, sep="\t", chunksize=500000, low_memory=False):
    m = (ch["chrom"] == CHROM) & (ch["start"] >= LO) & (ch["end"] <= HI)
    if m.any():
        keep.append(ch.loc[m])
b = pd.concat(keep, ignore_index=True)
print("RC BEARING bins:", len(b))

# track composition across the whole RC window (mean absolute differential KL)
comp = {TRACKS[k]: float(b[k].abs().mean()) for k in TRACKS}
tot = sum(comp.values()) or 1.0
print("\nMean |dKL| per track across RC (share of joint):")
for name, v in sorted(comp.items(), key=lambda x: -x[1]):
    print(f"  {name:7s} {v:.4f}  ({100*v/tot:4.1f}%)")

for strand, egf, rnacol in (("pos", EDGER_POS, RNA_POS), ("neg", EDGER_NEG, RNA_NEG)):
    print(f"\n=============== strand {strand}  (BEARING {rnacol}={TRACKS[rnacol]}) ===============")
    eg = pd.read_csv(egf)
    eg = eg[(eg["chr"] == CHROM) & (eg["start"] >= LO) & (eg["end"] <= HI)]
    m = b.merge(eg[["start", "logFC", "FDR"]], on="start", how="inner")
    print("edgeR-tested bins in RC:", len(m), "| edgeR sig:", int((m["FDR"] <= FDR).sum()))
    if len(m) == 0:
        continue
    sig = m[m["FDR"] <= FDR]
    # does BEARING's RNA attribution move in DN-direction at edgeR-sig bins?
    pos_frac = float((sig[rnacol] > 0).mean()) if len(sig) else float("nan")
    print(f"At edgeR-sig bins (n={len(sig)}):")
    print(f"  BEARING {TRACKS[rnacol]} dKL > 0 (DN-direction): {100*pos_frac:.0f}% of bins")
    print(f"  median {TRACKS[rnacol]} dKL = {sig[rnacol].median():.4f} "
          f"(range {sig[rnacol].min():.3f} .. {sig[rnacol].max():.3f})")
    print(f"  median joint bearing_score = {sig['bearing_score'].median():.4f}")
    print(f"  joint pval: min={sig['pval'].min():.2e} median={sig['pval'].median():.2e} "
          f"(none cleared BH; how far from raw 0.05: "
          f"{int((sig['pval']<=0.05).sum())}/{len(sig)} have raw pval<=0.05)")
    # do edgeR effect size and BEARING RNA attribution track each other?
    rho, p = spearmanr(m["logFC"], m[rnacol])
    rho2, p2 = spearmanr(m["logFC"], m["bearing_score"])
    print(f"  Spearman edgeR logFC vs BEARING {TRACKS[rnacol]} dKL: rho={rho:.2f} (p={p:.1e})")
    print(f"  Spearman edgeR logFC vs BEARING joint score:        rho={rho2:.2f} (p={p2:.1e})")