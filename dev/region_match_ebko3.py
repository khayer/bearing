#!/usr/bin/env python3
# region_match_ebko3.py  -- run from the bearing repo root
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

DIFF = "workflow/results/pvalue/diff_DN_vs_EbKO.stats.tsv"
EDGER_POS = "rna_concordance_DN_vs_EbKO_pos_edgeR_allbins.csv"
EDGER_NEG = "rna_concordance_DN_vs_EbKO_neg_edgeR_allbins.csv"
CHROM, LO, HI = "chr6", 41500000, 41567200
FDR = 0.05
BIN = 200
TRACKS = {"kl_1": "ATAC", "kl_2": "RNA+", "kl_3": "RNA-",
          "kl_4": "CTCF", "kl_5": "RAD21", "kl_6": "H3K27ac"}

keep = []
for ch in pd.read_csv(DIFF, sep="\t", chunksize=500000, low_memory=False):
    m = (ch["chrom"] == CHROM) & (ch["start"] >= LO) & (ch["end"] <= HI)
    if m.any():
        keep.append(ch.loc[m])
b = pd.concat(keep, ignore_index=True)
b["start"] = b["start"].astype(int)
b["bid"] = b["start"] // BIN
print("RC BEARING bins:", len(b))
print("BEARING start head:", sorted(b['start'].tolist())[:5])

for strand, egf, rnacol in (("pos", EDGER_POS, "kl_2"), ("neg", EDGER_NEG, "kl_3")):
    print(f"\n=============== strand {strand}  (BEARING {rnacol}={TRACKS[rnacol]}) ===============")
    eg = pd.read_csv(egf)
    eg = eg[(eg["chr"] == CHROM) & (eg["start"] >= LO - BIN) & (eg["end"] <= HI + BIN)].copy()
    eg["start"] = eg["start"].astype(int)
    print("edgeR start head:", sorted(eg['start'].tolist())[:5])
    # map edgeR start to the 200bp grid; (start-1)//BIN absorbs a 1-based offset
    eg["bid"] = (eg["start"] - 1) // BIN
    # sanity: which key aligns best
    n_grid = b["bid"].isin(set(eg["bid"])).sum()
    print(f"BEARING bins whose bin-id matches an edgeR bin-id: {n_grid}")
    m = b.merge(eg[["bid", "logFC", "FDR"]], on="bid", how="inner")
    print("merged bins:", len(m), "| edgeR sig (FDR<=0.05):", int((m["FDR"] <= FDR).sum()))
    if len(m) == 0:
        print("  STILL EMPTY -- paste the two head lines above so we can see the offset")
        continue
    sig = m[m["FDR"] <= FDR]
    if len(sig):
        print(f"At edgeR-sig bins (n={len(sig)}):")
        print(f"  {TRACKS[rnacol]} dKL > 0 (DN-direction): {100*(sig[rnacol]>0).mean():.0f}%")
        print(f"  median {TRACKS[rnacol]} dKL = {sig[rnacol].median():.4f} "
              f"(range {sig[rnacol].min():.3f}..{sig[rnacol].max():.3f})")
        print(f"  median joint bearing_score = {sig['bearing_score'].median():.4f}")
        print(f"  raw joint pval<=0.05: {int((sig['pval']<=0.05).sum())}/{len(sig)}; "
              f"min joint pval={sig['pval'].min():.2e}")
    rho, p = spearmanr(m["logFC"], m[rnacol])
    rho2, p2 = spearmanr(m["logFC"], m["bearing_score"])
    print(f"  Spearman edgeR logFC vs BEARING {TRACKS[rnacol]} dKL: rho={rho:.2f} (p={p:.1e})")
    print(f"  Spearman edgeR logFC vs BEARING joint score:        rho={rho2:.2f} (p={p2:.1e})")