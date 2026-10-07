#!/usr/bin/env python3
# ----------------------------------------------------------------------
# File     : hic/bench_normalization.py
# Version  : 2.0.0
# Date     : 2026-10-07
# Authors  : Katharina E. Hayer (katharinaehayer@gmail.com) and Claude
#            (Anthropic), co-created
# History  : v1 (2026-08-22/24) was written in a chat session and never reached
#            the repo; it is lost. This is a REBUILD from the documented method
#            and numbers (project docs normalization_benchmark.md "Benchmark 1"
#            and fig67_crosscheck_and_v1p_supplement.md sections 5, 6 and 9).
#            --expect compares every documented number with the rebuild.
# ----------------------------------------------------------------------
"""
bench_normalization.py
======================
Which insulation normalization (raw / detrend / quantile) should the
DN-versus-condition insulation comparison use? Shares the per-bin delta, the
normalizations and the top-5% autosomal threshold with insulation_specificity.py
(imported, not reimplemented), and adds three benchmark readouts:

  1. GROUND-TRUTH RANKING (Benchmark 1). V1P and EbKO are chr6 cis-deletions,
     so their target bins are real and every other autosomal bin is a
     guaranteed technical null. For each target (condition, region) the bins'
     |delta| are scored against all other autosomal bins: AUROC (Mann-Whitney)
     and average precision (AP). The aggregate is the mean over targets.
  2. REGION COUNTS. For each (condition, region): top-5% flagged bins / total
     and max |delta|  (fig67 sections 6 and 9).
  3. GENOME BLOCKS. Runs of >= --min-block consecutive flagged bins genome-wide
     per condition (fig67 section 5/9: "genome blocks >= 10"). Reported over
     autosomes AND over all chromosomes, because v1's convention is not
     documented; --expect shows which one reproduces v1.

INPUTS
  --insul COND=path.bm ...  (incl. the reference)   --reference DN
  --norms raw,detrend,quantile
  --regions TSV  name<TAB>chrom:start-end           (count table regions)
  --targets TSV  condition<TAB>region_name          (ground-truth targets)
  --expect  TSV  metric norm condition region expected tol   (published v1)

OUTPUTS (--out-dir D)
  D/bench_thresholds.tsv   norm x condition top-5% |delta| threshold
  D/bench_regions.tsv      norm x condition x region: sig, total, max|d|
  D/bench_blocks.tsv       norm x condition: blocks (autosomal / all chroms)
  D/bench_ranking.tsv      norm x target: AUROC, AP, n_target, n_null;
                           plus norm x "aggregate" (mean over targets)
  D/bench_summary.md       readable summary
  D/bench_expect_check.tsv (with --expect)

ASCII only. Reads real data; fabricates nothing.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import insulation_specificity as ins  # noqa: E402  (shared method)


def read_regions(path):
    out = []
    with open(path) as fh:
        for ln in fh:
            if not ln.strip() or ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if f[0] == "name":
                continue
            c, s, e = ins.parse_region(f[1])
            out.append({"name": f[0], "chrom": c, "start": s, "end": e})
    return out


def read_targets(path):
    out = []
    with open(path) as fh:
        for ln in fh:
            if not ln.strip() or ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if f[0] == "condition":
                continue
            out.append((f[0].strip(), f[1].strip()))
    return out


def auroc_ap(scores_pos, scores_neg):
    """AUROC (Mann-Whitney, ties = 0.5) and average precision."""
    pos = np.asarray(scores_pos, dtype=float)
    neg = np.asarray(scores_neg, dtype=float)
    allv = np.concatenate([pos, neg])
    ranks = pd.Series(allv).rank(method="average").to_numpy()
    n1, n0 = len(pos), len(neg)
    auc = (ranks[:n1].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0)
    # average precision: mean over positives of precision at their rank
    lab = np.concatenate([np.ones(n1), np.zeros(n0)])
    order = np.argsort(-allv, kind="mergesort")
    lab = lab[order]
    tp = np.cumsum(lab)
    prec = tp / np.arange(1, len(lab) + 1)
    ap = float((prec * lab).sum() / n1)
    return float(auc), ap


def count_blocks(m, sig, min_block, chroms=None):
    """Runs of consecutive flagged bins (adjacent starts on one chromosome)."""
    n = 0
    for c, idx in m.groupby("chrom", sort=False).indices.items():
        if chroms is not None and c not in chroms:
            continue
        idx = np.sort(idx)
        st = m["start"].to_numpy()[idx]
        sg = sig[idx]
        run = 0
        prev = None
        for s_, g in zip(st, sg):
            contiguous = prev is not None and s_ - prev == BIN
            if g and (run == 0 or contiguous):
                run += 1
            elif g:
                if run >= min_block:
                    n += 1
                run = 1
            else:
                if run >= min_block:
                    n += 1
                run = 0
            prev = s_
        if run >= min_block:
            n += 1
    return n


BIN = 25000


def main():
    global BIN
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--insul", nargs="+", required=True)
    ap.add_argument("--reference", default="DN")
    ap.add_argument("--conditions", default=None)
    ap.add_argument("--norms", default="raw,detrend,quantile")
    ap.add_argument("--detrend-window", type=int, default=2_000_000)
    ap.add_argument("--top-pct", type=float, default=95.0)
    ap.add_argument("--min-block", type=int, default=10)
    ap.add_argument("--regions", required=True)
    ap.add_argument("--targets", required=True)
    ap.add_argument("--expect", default=None)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()

    paths = dict(s.split("=", 1) for s in a.insul)
    if a.reference not in paths:
        sys.exit("reference %s not among --insul" % a.reference)
    conds = ([c.strip() for c in a.conditions.split(",")] if a.conditions
             else [c for c in paths if c != a.reference])
    norms = [n.strip() for n in a.norms.split(",") if n.strip()]
    regions = read_regions(a.regions)
    reg_by = {r["name"]: r for r in regions}
    targets = read_targets(a.targets)
    for c, r in targets:
        if r not in reg_by:
            sys.exit("target region %s not in --regions" % r)
    os.makedirs(a.out_dir, exist_ok=True)

    ref = ins.load_bm(paths[a.reference])
    BIN = int(np.median(ref["end"] - ref["start"]))
    dbins = max(3, int(round(a.detrend_window / BIN)))
    cond_df = {c: ins.load_bm(paths[c]) for c in conds}
    print("bin %d bp; conditions %s; norms %s" % (BIN, conds, norms), flush=True)

    thr_rows, reg_rows, blk_rows, rank_rows = [], [], [], []
    for nm in norms:
        for c in conds:
            m = ins.build_delta(ref, cond_df[c], nm, dbins)
            auto_mask = m["chrom"].isin(ins.AUTOSOMES).to_numpy()
            absd = m["absd"].to_numpy()
            thr = float(np.percentile(absd[auto_mask & np.isfinite(absd)], a.top_pct))
            sig = absd >= thr
            thr_rows.append({"norm": nm, "condition": c, "threshold": thr})
            for r in regions:
                mk = ins.core_mask(m, r)
                d = absd[mk]
                reg_rows.append({"norm": nm, "condition": c, "region": r["name"],
                                 "sig": int((d >= thr).sum()), "total": int(mk.sum()),
                                 "max_absd": float(d.max()) if d.size else np.nan})
            blk_rows.append({"norm": nm, "condition": c,
                             "blocks_autosomal": count_blocks(m, sig, a.min_block,
                                                              ins.AUTOSOMES),
                             "blocks_all": count_blocks(m, sig, a.min_block, None)})
            for tc, tr in targets:
                if tc != c:
                    continue
                tmk = ins.core_mask(m, reg_by[tr])
                pos = absd[tmk & auto_mask]
                neg = absd[(~tmk) & auto_mask & np.isfinite(absd)]
                if pos.size == 0:
                    continue
                auc, apv = auroc_ap(pos, neg)
                rank_rows.append({"norm": nm, "target": "%s:%s" % (tc, tr),
                                  "auroc": auc, "ap": apv,
                                  "n_target": int(pos.size), "n_null": int(neg.size)})
            print("  %-9s %-6s threshold %.3f" % (nm, c, thr), flush=True)
        sub = [r for r in rank_rows if r["norm"] == nm]
        if sub:
            rank_rows.append({"norm": nm, "target": "aggregate",
                              "auroc": float(np.mean([r["auroc"] for r in sub])),
                              "ap": float(np.mean([r["ap"] for r in sub])),
                              "n_target": int(sum(r["n_target"] for r in sub)),
                              "n_null": np.nan})

    T = pd.DataFrame(thr_rows); R = pd.DataFrame(reg_rows)
    B = pd.DataFrame(blk_rows); K = pd.DataFrame(rank_rows)
    D = a.out_dir
    T.to_csv(os.path.join(D, "bench_thresholds.tsv"), sep="\t", index=False, float_format="%.4g")
    R.to_csv(os.path.join(D, "bench_regions.tsv"), sep="\t", index=False, float_format="%.4g")
    B.to_csv(os.path.join(D, "bench_blocks.tsv"), sep="\t", index=False)
    K.to_csv(os.path.join(D, "bench_ranking.tsv"), sep="\t", index=False, float_format="%.4g")

    with open(os.path.join(D, "bench_summary.md"), "w") as fh:
        fh.write("# Insulation normalization benchmark (rebuild of v1)\n\n")
        fh.write("## Ground-truth ranking (|delta| of target bins vs all other autosomal bins)\n\n")
        fh.write("| target | " + " | ".join("%s AUROC / AP" % n for n in norms) + " |\n")
        fh.write("|---|" + "---|" * len(norms) + "\n")
        for t in K["target"].unique():
            cells = []
            for n in norms:
                q = K[(K["target"] == t) & (K["norm"] == n)]
                cells.append("%.3f / %.4f" % (q["auroc"].iloc[0], q["ap"].iloc[0])
                             if len(q) else "NA")
            fh.write("| %s | %s |\n" % (t, " | ".join(cells)))
        fh.write("\n## Flagged (top-%g%%) bins per region, sig/total (max |d|)\n\n"
                 % (100 - a.top_pct))
        fh.write("| norm | condition | threshold | " + " | ".join(r["name"] for r in regions)
                 + " | blocks>=%d (auto / all) |\n" % a.min_block)
        fh.write("|---|---|---|" + "---|" * len(regions) + "---|\n")
        for n in norms:
            for c in conds:
                th = T[(T["norm"] == n) & (T["condition"] == c)]["threshold"].iloc[0]
                cells = []
                for r in regions:
                    q = R[(R["norm"] == n) & (R["condition"] == c) & (R["region"] == r["name"])].iloc[0]
                    cells.append("%d/%d (%.2f)" % (q["sig"], q["total"], q["max_absd"])
                                 if q["total"] else "-")
                b = B[(B["norm"] == n) & (B["condition"] == c)].iloc[0]
                fh.write("| %s | %s | %.3f | %s | %d / %d |\n"
                         % (n, c, th, " | ".join(cells), b["blocks_autosomal"], b["blocks_all"]))
    print("wrote bench_thresholds / bench_regions / bench_blocks / bench_ranking .tsv + bench_summary.md")

    if not a.expect:
        return
    exp = pd.read_csv(a.expect, sep="\t", comment="#", dtype=str)
    out = []
    for _, e in exp.iterrows():
        metric, nm, c, reg = e["metric"], e["norm"], e["condition"], e.get("region", "")
        want = float(e["expected"]); tol = float(e["tol"])
        got = np.nan
        if metric == "threshold":
            q = T[(T["norm"] == nm) & (T["condition"] == c)]
            got = q["threshold"].iloc[0] if len(q) else np.nan
        elif metric in ("sig", "total", "max_absd"):
            q = R[(R["norm"] == nm) & (R["condition"] == c) & (R["region"] == reg)]
            got = q[metric].iloc[0] if len(q) else np.nan
        elif metric in ("blocks_autosomal", "blocks_all"):
            q = B[(B["norm"] == nm) & (B["condition"] == c)]
            got = q[metric].iloc[0] if len(q) else np.nan
        elif metric in ("auroc", "ap"):
            tgt = "aggregate" if c == "aggregate" else "%s:%s" % (c, reg)
            q = K[(K["norm"] == nm) & (K["target"] == tgt)]
            got = q[metric].iloc[0] if len(q) else np.nan
        ok = bool(np.isfinite(got) and abs(got - want) <= tol)
        out.append({"metric": metric, "norm": nm, "condition": c, "region": reg,
                    "published": want, "rebuilt": got, "tol": tol, "ok": ok,
                    "note": e.get("note", "")})
    X = pd.DataFrame(out)
    X.to_csv(os.path.join(D, "bench_expect_check.tsv"), sep="\t", index=False, float_format="%.4g")
    print("\nEXPECT CHECK vs published v1: %d / %d within tolerance" % (int(X["ok"].sum()), len(X)))
    if (~X["ok"]).any():
        print(X.loc[~X["ok"], ["metric", "norm", "condition", "region", "published",
                               "rebuilt", "tol", "note"]].to_string(index=False))


if __name__ == "__main__":
    main()
