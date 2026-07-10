#!/usr/bin/env python3
# region_match_ebko.py  -- run from the bearing repo root
import sys, re
import numpy as np
import pandas as pd

DIFF = "workflow/results/pvalue/diff_DN_vs_EbKO.stats.tsv"
EDGER = {  # genome-wide edgeR CSVs written by rna_concordance_stranded.R
    "pos": "rna_concordance_DN_vs_EbKO_pos_edgeR_allbins.csv",
    "neg": "rna_concordance_DN_vs_EbKO_neg_edgeR_allbins.csv",
}
WINDOWS = {
    "rc_edger":  (41525453, 41567200),   # regions.tsv 'rc' (matches the edgeR file)
    "DJ_RC_ms":  (41500000, 41551500),   # manuscript DJ_RC
}
CHROM = "chr6"
FDR = 0.05

def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    q = np.empty(n)
    q[order] = np.minimum.accumulate((p[order] * n / (np.arange(n) + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)

# 1. pull only the RC neighbourhood out of the big BEARING file
lo = min(w[0] for w in WINDOWS.values())
hi = max(w[1] for w in WINDOWS.values())
chunks = []
for ch in pd.read_csv(DIFF, sep="\t", chunksize=500000, low_memory=False):
    m = (ch["chrom"] == CHROM) & (ch["start"] >= lo) & (ch["end"] <= hi)
    if m.any():
        chunks.append(ch.loc[m])
b = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame()
print("BEARING bins pulled near RC:", len(b))
print("BEARING columns:", list(b.columns))

# auto-detect the RNA+ per-track p-value column, if present
def find_track_p(cols, want_pos):
    pat = re.compile(r"^pval_(perm_|self_)?(?!adj_).*rna", re.I)
    cand = [c for c in cols if pat.search(c) and "adj" not in c.lower()]
    if want_pos:
        cand = [c for c in cand if re.search(r"(\+|pos|plus|fwd)", c, re.I)]
    else:
        cand = [c for c in cand if re.search(r"(-|neg|minus|rev)", c, re.I)]
    return cand[0] if cand else None

for strand, want_pos in (("pos", True), ("neg", False)):
    print("\n================= strand:", strand, "=================")
    eg = pd.read_csv(EDGER[strand])
    eg = eg[eg["chr"] == CHROM]
    for wname, (wlo, whi) in WINDOWS.items():
        bw = b[(b["start"] >= wlo) & (b["end"] <= whi)].copy()
        egw = eg[(eg["start"] >= wlo) & (eg["end"] <= whi)].copy()
        # edgeR side (genome-wide FDR already in file)
        eg_tested = len(egw)
        eg_sig = int((egw["FDR"] <= FDR).sum())
        # (A) BEARING genome-wide BH, counted in window, DN-direction
        a_sig = int(((bw["pval_adj_bh"] <= FDR) & (bw["direction"] == "+")).sum())
        # (B) BEARING region-restricted BH on the JOINT pval
        bw["q_region"] = bh(bw["pval"].values)
        b_sig = int(((bw["q_region"] <= FDR) & (bw["direction"] == "+")).sum())
        # (B') region-restricted BH on joint pval, on edgeR's tested universe only
        merged = bw.merge(egw[["start"]], on="start", how="inner")
        if len(merged):
            merged["q_u"] = bh(merged["pval"].values)
            bp_sig = int(((merged["q_u"] <= FDR) & (merged["direction"] == "+")).sum())
        else:
            bp_sig = 0
        # (C) BEARING RNA per-track p, region-restricted BH, if the column exists
        tcol = find_track_p(b.columns, want_pos)
        if tcol:
            bw["q_track"] = bh(bw[tcol].values)
            c_sig = int((bw["q_track"] <= FDR).sum())
            ctxt = f"{c_sig} (col {tcol})"
        else:
            ctxt = "no per-track p column in file"
        print(f"[{wname} {wlo}-{whi}]  edgeR tested={eg_tested} sig={eg_sig} | "
              f"BEARING bins={len(bw)} | (A)gw-BH={a_sig} (B)region-BH-joint={b_sig} "
              f"(B')region-BH-joint/edgeR-univ={bp_sig} (C)region-BH-RNA-track={ctxt}")