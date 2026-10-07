#!/usr/bin/env python3
# ----------------------------------------------------------------------
# File     : hic/insulation_specificity.py
# Version  : 2.0.0
# Date     : 2026-10-07
# Authors  : Katharina E. Hayer (katharinaehayer@gmail.com) and Claude
#            (Anthropic), co-created
# History  : v1 (2026-08-22/24, with bench_specificity.py and
#            bench_normalization.py) was written in a chat session and never
#            reached the repo; it and its *_qn_hl.pdf outputs are lost. This is
#            a REBUILD from the documented method (project docs
#            normalization_benchmark.md and fig67_crosscheck_and_v1p_supplement.md,
#            sections 5-10). --expect compares the result with the published
#            specificity table so any difference from v1 is visible, not silent.
# ----------------------------------------------------------------------
"""
insulation_specificity.py
=========================
Insulation-comparison specificity controls: does a DN-versus-condition
insulation comparison report an architectural change at the right locus only in
the right condition?

METHOD (as documented for v1)
-----------------------------
Per condition C (reference R = DN):
  1. Insulation = HiCExplorer TAD-separation score (.bm, 25 kb bins); when the
     .bm has several window columns, their mean is used (the repo convention,
     see hic/bes_hic_crosslocus.py load_insulation_bm).
  2. Normalization (--norm):
       raw       no change
       quantile  R and C quantile-normalized to a common distribution (mean of
                 the two sorted genome-wide vectors) -- removes GLOBAL scale /
                 shape differences. Preferred (v1 benchmark).
       detrend   delta minus a per-chromosome rolling median of delta
                 (--detrend-window, default 2 Mb) -- removes REGIONAL offsets.
  3. delta = C - R per bin (after normalization).
  4. Genome-wide reference: |delta| over autosomal bins (chr1-chr19).
     Top-5% threshold = 95th percentile of that distribution; a bin is "sig"
     if |delta| >= threshold (a PERCENTILE FLAG, 5% by construction, NOT a
     significance test).
  5. Element core (gene-centered 125 kb = 5 x 25 kb bins): bins whose start
     lies in [start, end). Specificity metric = percentile rank of the core's
     median |delta| within the condition's own genome-wide autosomal |delta|
     (scale-free, so conditions with different noise are comparable).
  6. Per element: target condition (positives), target percentile, margin =
     target minus the best other condition, rank of the target. Negatives: the
     maximum percentile any condition reaches, and which.

OUTPUTS (--out-dir D, --norm N)
-------------------------------
  D/specificity_N.tsv      element x condition: percentile, median|d|,
                           sig_bins/total, max|d|; plus target/margin/rank
  D/specificity_N.md       the same as the published-style summary table
  D/<element>_N_hl.pdf/.png  per element: per condition, R and C normalized
                           insulation over the window + delta with the top-5%
                           threshold; grey band = the compared core
  D/expect_check_N.tsv     (with --expect) observed vs published percentile

INPUTS
------
  --insul COND=path.bm ...   one per condition, including the reference
  --reference DN
  --elements TSV             name, core (chrom:start-end), role
                             (positive/negative/excluded), target condition
                             (positives only; blank otherwise)
  --flank BP                 plot window = core +/- flank (default 1 Mb)

USAGE
  python hic/insulation_specificity.py \\
     --insul DN=..._DN_bs_25000_tad_score.bm dV1P=..._dV1P_bs_25000_tad_score.bm \\
             EbKO=... DP=... ProB=... S3T3=... \\
     --reference DN --elements annotations/insulation_controls.tsv \\
     --norm quantile --expect paper/expected/insulation_specificity_quantile.tsv \\
     --out-dir results/hic/insulation_specificity

ASCII only. Reads real data; fabricates nothing.
"""

import argparse
import os
import re
import sys

import numpy as np
import pandas as pd

AUTOSOMES = {"chr%d" % i for i in range(1, 20)}


# ----------------------------------------------------------------------
# I/O
# ----------------------------------------------------------------------
def load_bm(path):
    """HiCExplorer .bm -> DataFrame chrom,start,end,score (mean of score cols)."""
    rows = []
    with open(path) as fh:
        for ln in fh:
            if not ln.strip() or ln.startswith("#") or ln.startswith("track"):
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            try:
                s = int(parts[1]); e = int(parts[2])
            except ValueError:
                continue
            vals = []
            for v in parts[3:]:
                try:
                    fv = float(v)
                except ValueError:
                    continue
                if np.isfinite(fv):
                    vals.append(fv)
            if vals:
                rows.append((parts[0], s, e, float(np.mean(vals))))
    if not rows:
        sys.exit("no usable rows in %s" % path)
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "score"])


def parse_region(txt):
    m = re.match(r"^(\w+):([\d,]+)-([\d,]+)$", txt.strip())
    if not m:
        sys.exit("bad region %r (need chrom:start-end)" % txt)
    return m.group(1), int(m.group(2).replace(",", "")), int(m.group(3).replace(",", ""))


def read_elements(path):
    out = []
    with open(path) as fh:
        for ln in fh:
            if not ln.strip() or ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if f[0] == "name":
                continue
            while len(f) < 4:
                f.append("")
            c, s, e = parse_region(f[1])
            out.append({"name": f[0], "chrom": c, "start": s, "end": e,
                        "role": f[2].strip() or "test", "target": f[3].strip()})
    return out


# ----------------------------------------------------------------------
# Normalization and delta
# ----------------------------------------------------------------------
def quantile_normalize_pair(a, b):
    """Map a and b onto the mean of their sorted values (ties by rank order)."""
    ra = np.argsort(np.argsort(a, kind="mergesort"), kind="mergesort")
    rb = np.argsort(np.argsort(b, kind="mergesort"), kind="mergesort")
    ref = (np.sort(a) + np.sort(b)) / 2.0
    return ref[ra], ref[rb]


def rolling_median_by_chrom(df, col, window_bins):
    out = np.empty(len(df))
    for _, idx in df.groupby("chrom", sort=False).indices.items():
        v = pd.Series(df[col].to_numpy()[idx])
        out[idx] = v.rolling(window_bins, center=True, min_periods=1).median().to_numpy()
    return out


def build_delta(ref_df, cond_df, norm, detrend_bins):
    m = ref_df.merge(cond_df, on=["chrom", "start", "end"], suffixes=("_r", "_c"))
    m = m.sort_values(["chrom", "start"]).reset_index(drop=True)
    r = m["score_r"].to_numpy(dtype=float)
    c = m["score_c"].to_numpy(dtype=float)
    if norm == "quantile":
        r, c = quantile_normalize_pair(r, c)
    m["r"] = r
    m["c"] = c
    m["delta"] = c - r
    if norm == "detrend":
        m["delta"] = m["delta"] - rolling_median_by_chrom(m, "delta", detrend_bins)
    m["absd"] = np.abs(m["delta"])
    return m


def core_mask(m, el):
    return ((m["chrom"] == el["chrom"]) & (m["start"] >= el["start"])
            & (m["start"] < el["end"])).to_numpy()


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--insul", nargs="+", required=True, help="COND=path.bm")
    ap.add_argument("--reference", default="DN")
    ap.add_argument("--elements", required=True)
    ap.add_argument("--norm", default="quantile", choices=["raw", "quantile", "detrend"])
    ap.add_argument("--detrend-window", type=int, default=2_000_000)
    ap.add_argument("--top-pct", type=float, default=95.0,
                    help="percentile for the per-condition sig-bin flag (default 95)")
    ap.add_argument("--flank", type=int, default=1_000_000)
    ap.add_argument("--conditions", default=None,
                    help="comma list (order of columns/rows); default = all non-reference")
    ap.add_argument("--expect", default=None,
                    help="TSV element<TAB>condition<TAB>percentile (published v1 values)")
    ap.add_argument("--expect-tol", type=float, default=0.03)
    ap.add_argument("--no-plots", action="store_true")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()

    paths = {}
    for spec in a.insul:
        if "=" not in spec:
            sys.exit("--insul needs COND=path, got %r" % spec)
        k, v = spec.split("=", 1)
        paths[k.strip()] = v.strip()
    if a.reference not in paths:
        sys.exit("reference %s not among --insul" % a.reference)
    conds = ([c.strip() for c in a.conditions.split(",")] if a.conditions
             else [c for c in paths if c != a.reference])
    missing = [c for c in conds if c not in paths]
    if missing:
        sys.exit("conditions without --insul: %s" % missing)
    os.makedirs(a.out_dir, exist_ok=True)
    elements = read_elements(a.elements)

    print("loading reference %s" % a.reference, flush=True)
    ref = load_bm(paths[a.reference])
    binsize = int(np.median(ref["end"] - ref["start"]))
    detrend_bins = max(3, int(round(a.detrend_window / binsize)))
    print("  bin size %d bp; detrend window %d bins" % (binsize, detrend_bins))

    per = {}      # cond -> merged frame
    thr = {}
    gw = {}       # cond -> sorted autosomal |delta|
    for c in conds:
        print("condition %s" % c, flush=True)
        m = build_delta(ref, load_bm(paths[c]), a.norm, detrend_bins)
        auto = m["absd"][m["chrom"].isin(AUTOSOMES)].to_numpy()
        auto = auto[np.isfinite(auto)]
        if auto.size == 0:
            sys.exit("no autosomal bins for %s" % c)
        gw[c] = np.sort(auto)
        thr[c] = float(np.percentile(auto, a.top_pct))
        per[c] = m
        print("  %d common bins; top-%g%% threshold |d| >= %.3f"
              % (len(m), 100 - a.top_pct, thr[c]))

    # ---- per element x condition ------------------------------------------
    rows = []
    for el in elements:
        for c in conds:
            m = per[c]
            mk = core_mask(m, el)
            n = int(mk.sum())
            if n == 0:
                rows.append({"element": el["name"], "role": el["role"],
                             "target": el["target"], "condition": c,
                             "core_bins": 0, "median_absd": np.nan,
                             "percentile": np.nan, "sig_bins": 0, "max_absd": np.nan})
                continue
            d = m.loc[mk, "absd"].to_numpy()
            med = float(np.median(d))
            pct = float(np.searchsorted(gw[c], med, side="right")) / gw[c].size
            rows.append({"element": el["name"], "role": el["role"],
                         "target": el["target"], "condition": c, "core_bins": n,
                         "median_absd": med, "percentile": pct,
                         "sig_bins": int((d >= thr[c]).sum()),
                         "max_absd": float(d.max())})
    long = pd.DataFrame(rows)

    # ---- summary: target / margin / rank (positives); max (negatives) -------
    summ = []
    for el in elements:
        sub = long[long["element"] == el["name"]].set_index("condition")
        pcts = {c: sub.loc[c, "percentile"] for c in conds}
        rec = {"element": el["name"], "role": el["role"], "target": el["target"]}
        rec.update({c: pcts[c] for c in conds})
        if el["target"] and el["target"] in pcts and np.isfinite(pcts[el["target"]]):
            t = pcts[el["target"]]
            others = [v for k, v in pcts.items() if k != el["target"] and np.isfinite(v)]
            rec["margin"] = t - max(others) if others else np.nan
            rec["rank"] = 1 + sum(1 for v in others if v > t)
            rec["max_pct"] = max([t] + others)
            rec["max_cond"] = max(pcts, key=lambda k: pcts[k] if np.isfinite(pcts[k]) else -1)
        else:
            fin = {k: v for k, v in pcts.items() if np.isfinite(v)}
            rec["margin"] = np.nan
            rec["rank"] = np.nan
            rec["max_pct"] = max(fin.values()) if fin else np.nan
            rec["max_cond"] = max(fin, key=fin.get) if fin else ""
        summ.append(rec)
    summ = pd.DataFrame(summ)

    base = os.path.join(a.out_dir, "specificity_%s" % a.norm)
    long.to_csv(base + "_long.tsv", sep="\t", index=False, float_format="%.4g")
    summ.to_csv(base + ".tsv", sep="\t", index=False, float_format="%.4g")
    with open(base + ".md", "w") as fh:
        fh.write("| element | " + " | ".join(conds) + " | target | margin | rank | max (cond) |\n")
        fh.write("|" + "---|" * (len(conds) + 5) + "\n")
        for _, r in summ.iterrows():
            cells = []
            for c in conds:
                v = r[c]
                s = "%.3g" % v if np.isfinite(v) else "NA"
                cells.append("**%s**" % s if c == r["target"] else s)
            mg = "%+.2f" % r["margin"] if np.isfinite(r["margin"]) else ""
            rk = "%d" % r["rank"] if np.isfinite(r["rank"]) else ""
            mx = ("%.2f (%s)" % (r["max_pct"], r["max_cond"])
                  if np.isfinite(r["max_pct"]) else "")
            name = r["element"] if r["role"] == "positive" else "_%s (%s)_" % (r["element"], r["role"])
            fh.write("| %s | %s | %s | %s | %s | %s |\n"
                     % (name, " | ".join(cells), r["target"] or "none", mg, rk, mx))
    with open(base + "_thresholds.tsv", "w") as fh:
        fh.write("condition\ttop_pct\tthreshold_absd\tn_autosomal_bins\n")
        for c in conds:
            fh.write("%s\t%g\t%.6g\t%d\n" % (c, a.top_pct, thr[c], gw[c].size))
    print("wrote %s.tsv / .md / _long.tsv / _thresholds.tsv" % base)

    # ---- comparison with the published v1 table -----------------------------
    if a.expect:
        exp = pd.read_csv(a.expect, sep="\t", comment="#")
        need = {"element", "condition", "percentile"}
        if not need.issubset(exp.columns):
            sys.exit("--expect needs columns %s" % sorted(need))
        j = exp.merge(long[["element", "condition", "percentile"]],
                      on=["element", "condition"], how="left",
                      suffixes=("_published", "_rebuilt"))
        j["abs_diff"] = (j["percentile_rebuilt"] - j["percentile_published"]).abs()
        j["within_tol"] = j["abs_diff"] <= a.expect_tol
        j.to_csv(os.path.join(a.out_dir, "expect_check_%s.tsv" % a.norm),
                 sep="\t", index=False, float_format="%.4g")
        ok = int(j["within_tol"].sum())
        print("\nEXPECT CHECK vs published v1 (%s): %d / %d cells within +/- %.2f; "
              "max |diff| %.3f" % (a.norm, ok, len(j), a.expect_tol, j["abs_diff"].max()))
        if ok < len(j):
            print("  cells outside tolerance:")
            print(j.loc[~j["within_tol"], ["element", "condition", "percentile_published",
                                            "percentile_rebuilt", "abs_diff"]]
                  .to_string(index=False))

    if a.no_plots:
        return
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print("plots skipped: %s" % e)
        return
    for el in elements:
        w0 = max(0, el["start"] - a.flank)
        w1 = el["end"] + a.flank
        k = len(conds)
        fig, axes = plt.subplots(2 * k, 1, figsize=(7.5, 1.55 * 2 * k), sharex=True,
                                 gridspec_kw={"height_ratios": [1.6, 1] * k})
        for i, c in enumerate(conds):
            m = per[c]
            w = m[(m["chrom"] == el["chrom"]) & (m["start"] >= w0) & (m["start"] < w1)]
            x = (w["start"] + binsize / 2.0) / 1e6
            ax, axd = axes[2 * i], axes[2 * i + 1]
            for aa in (ax, axd):
                aa.axvspan(el["start"] / 1e6, el["end"] / 1e6, color="#cccccc", alpha=0.6, lw=0)
            ax.plot(x, w["r"], color="#444444", lw=1.1, label=a.reference)
            ax.plot(x, w["c"], color="#d62728", lw=1.1, label=c)
            pr = long[(long["element"] == el["name"]) & (long["condition"] == c)]
            pct = pr["percentile"].iloc[0] if len(pr) else np.nan
            tag = " (target)" if c == el["target"] else ""
            ax.set_title("%s vs %s%s   core median |d| percentile %.3f, sig %d/%d"
                         % (c, a.reference, tag, pct,
                            int(pr["sig_bins"].iloc[0]) if len(pr) else 0,
                            int(pr["core_bins"].iloc[0]) if len(pr) else 0),
                         fontsize=7.5, loc="left")
            ax.set_ylabel("insulation\n(%s)" % a.norm, fontsize=7)
            ax.legend(fontsize=6, loc="upper right", frameon=False)
            axd.bar(x, w["delta"], width=binsize / 1e6 * 0.9,
                    color=np.where(w["absd"] >= thr[c], "#d62728", "#999999"))
            axd.axhline(thr[c], color="#d62728", lw=0.6, ls="--")
            axd.axhline(-thr[c], color="#d62728", lw=0.6, ls="--")
            axd.axhline(0, color="black", lw=0.5)
            axd.set_ylabel("delta", fontsize=7)
            for aa in (ax, axd):
                aa.tick_params(labelsize=6)
                for sp in ("top", "right"):
                    aa.spines[sp].set_visible(False)
        axes[-1].set_xlabel("%s (Mb)" % el["chrom"], fontsize=7)
        fig.suptitle("%s (%s%s) - insulation, %s-normalized; grey = compared %d kb core"
                     % (el["name"], el["role"],
                        (", target " + el["target"]) if el["target"] else "",
                        a.norm, (el["end"] - el["start"]) // 1000), fontsize=8.5)
        fig.tight_layout(rect=(0, 0, 1, 0.97))
        stem = os.path.join(a.out_dir, "%s_%s_hl" % (el["name"].lower(), a.norm))
        fig.savefig(stem + ".pdf")
        fig.savefig(stem + ".png", dpi=200)
        plt.close(fig)
    print("wrote %d element figures (*_%s_hl.pdf/.png)" % (len(elements), a.norm))


if __name__ == "__main__":
    main()
