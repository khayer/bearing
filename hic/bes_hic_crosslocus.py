#!/usr/bin/env python3
# ----------------------------------------------------------------------
# File     : bes_hic_crosslocus.py
# Version  : 2.0.0
# Date     : 2026-09-26
# Authors  : Katharina E. Hayer (katharinaehayer@gmail.com) and Claude
#            (Anthropic), co-created
# Changes  : v2.0.0 adds --score-sets: pre-registered track-combination
#            scores (sum of |kl_i| over a chosen set of tracks) tested
#            on the SAME paired control panels as the summed BES, with
#            Holm correction across the pre-registered sets. Default
#            behaviour (no --score-sets) is unchanged from v1.
# ----------------------------------------------------------------------
"""
bes_hic_crosslocus.py

Cross-locus test: do AR loci (as a class) show stronger top-decile
BEARING-Hi-C co-localization than size-matched control regions?

This converts the weak single-locus correlation (one noisy locus, ~20
Hi-C bins, spatial autocorrelation) into a population-level claim with
real n and a clean empirical null.

DESIGN
------
For each TARGET (contrast, region) pair where AR biology is expected to
be active (e.g. Tcrb in DN_vs_DP, Igh in DN_vs_ProB):
  1. Bin the region into Hi-C bins of --hic-bin bp.
  2. Per Hi-C bin: aggregate constituent 200 bp |BES| (default p95),
     compute delta_contact (sum over partner bins within
     [min_distance, max_distance] of |balanced_A - balanced_B|), and
     delta_insulation (|insul_A - insul_B|, mean over score columns).
  3. Flag the WITHIN-REGION top decile of BES and of delta_contact.
  4. Build a 2x2 (both_top, bes_only, contact_only, neither).

Pool the 2x2 across all targets -> pooled OR_target.

NULL: for each target, draw --n-controls size-matched control regions
from the SAME contrast (autosomal, blacklist-excluded, not overlapping
any target, optionally gene-density-matched). A "control panel" is one
control per target; pool its 2x2 -> OR_control. Repeat to build a null
distribution of pooled OR. Empirical p = fraction of control panels
with OR >= OR_target.

QA CHECK (do this first): run with a single target = Tcrb_wide in
DN_vs_DP at 10 kb and confirm the per-region 2x2 matches your
bes_hic_correlation_v5.py output for the same region
(p95_bes x delta_contact: n_both_top=5, n_x_only=15, n_y_only=15,
n_neither=165 at 10 kb / chr6:40.4-42.4M). If it matches, the binning
is consistent and the cross-locus pooling is validated. The within-
region decile is rank-based, so sum-vs-mean contact aggregation does
not change the 2x2.

USAGE
-----
python bes_hic_crosslocus.py \\
    --diffs DN_vs_DP=results_v6/diff_DN_vs_DP.stats.tsv \\
            DN_vs_ProB=results_v6/diff_DN_vs_ProB.stats.tsv \\
    --cool-a ../hic_files/merged_corrected_KR_DN_bs_10000.cool \\
    --cool-b DN_vs_DP=../hic_files/merged_corrected_KR_DP_bs_10000.cool \\
             DN_vs_ProB=../hic_files/merged_corrected_KR_ProB_bs_10000.cool \\
    --insul-a .../merged_corrected_KR_DN_bs_10000_tad_score.bm \\
    --insul-b DN_vs_DP=.../merged_corrected_KR_DP_bs_10000_tad_score.bm \\
              DN_vs_ProB=.../merged_corrected_KR_ProB_bs_10000_tad_score.bm \\
    --targets targets.tsv \\
    --hic-bin 10000 --min-distance 50000 --max-distance 500000 \\
    --aggregation p95 \\
    --n-controls 200 --n-panels 1000 \\
    --blacklist mm10-blacklist.v2.bed \\
    --gtf gencode.vM23.annotation_modified_overlaps_removed_sorted.gtf \\
    --match-gene-density \\
    --out-prefix crosslocus_DNvsDP_DNvsProB

TRACK-COMBINATION SETS (v2, --score-sets)
-----------------------------------------
Each set is NAME=tok,tok,... . A token is either a per-track column
suffix as written by bearing_pvalue.py (e.g. CTCF -> kl_CTCF) or a
1-based track index (4 -> the 4th kl_* column in the diff TSV). The
set score per 200 bp bin is sum_i |kl_i| over the set; it is then
aggregated per Hi-C bin exactly like |BES| (--aggregation).
The resolved column names are printed so the mapping can be checked.

Pre-registered sets for the 6-track panel
(1 ATAC, 2 RNAseq+, 3 RNAseq-, 4 CTCF, 5 Cohesin, 6 H3K27ac):
    --score-sets architectural=4,5 transcriptional=2,3,6 \
                 accessible_enhancer=1,6

Summed |BES| is always run as the reference row and is NOT part of the
Holm family. Every score (BES and each set) is evaluated on identical
control panels (same control drawn per target per panel), so the sets
are directly comparable to each other and to BES.

Extra outputs with --score-sets:
    <prefix>.score_sets.tsv          per target x score 2x2 + OR
    <prefix>.score_sets_summary.tsv  pooled OR, null median/p95,
                                     empirical p, Holm-adjusted p

targets.tsv (tab-separated, header required):
    contrast    chrom   start   end     name
    DN_vs_DP    chr6    40400000  42400000  Tcrb_wide
    DN_vs_DP    chr11   96000000  97000000  Cd4
    DN_vs_ProB  chr12   113200000 116000000 Igh
    ...
(contrast must match a key in --diffs / --cool-b / --insul-b)
"""
import argparse
import os
import sys
import json
import numpy as np
import pandas as pd


# ----------------------------------------------------------------------
# Loaders
# ----------------------------------------------------------------------
def resolve_set_tokens(columns, tokens):
    """Map tokens (track-name suffix or 1-based index) to kl_* columns."""
    kl_cols = [c for c in columns if c.startswith("kl_")]
    out = []
    for tok in tokens:
        tok = tok.strip()
        if not tok:
            continue
        if ("kl_" + tok) in columns:
            out.append("kl_" + tok)
        elif tok in columns and tok.startswith("kl_"):
            out.append(tok)
        elif tok.isdigit():
            i = int(tok)
            if i < 1 or i > len(kl_cols):
                sys.exit("track index {} out of range; diff TSV has {} kl_* "
                         "columns: {}".format(i, len(kl_cols), kl_cols))
            out.append(kl_cols[i - 1])
        else:
            sys.exit("cannot resolve track token '{}'; kl_* columns are: {}"
                     .format(tok, kl_cols))
    if not out:
        sys.exit("empty score set")
    if len(set(out)) != len(out):
        sys.exit("duplicate tracks in score set: {}".format(out))
    return out


def load_bes(path, score_sets=None, resolved_log=None):
    """Load a BEARING diff TSV. Returns DataFrame with chrom,start,end,
    abs_bes, plus one column per score set (sum of |kl_i| over the set).
    Uses bearing_score column; falls back to bearing_score_tested."""
    header = pd.read_csv(path, sep="\t", nrows=0).columns.tolist()
    needed = {"chrom", "start", "end"}
    if not needed.issubset(header):
        sys.exit("BES TSV {} missing chrom/start/end".format(path))
    if "bearing_score" in header:
        bcol = "bearing_score"
    elif "bearing_score_tested" in header:
        bcol = "bearing_score_tested"
    else:
        sys.exit("BES TSV {} has no bearing_score column".format(path))
    set_cols = {}
    for name, toks in (score_sets or []):
        set_cols[name] = resolve_set_tokens(header, toks)
        if resolved_log is not None:
            resolved_log.setdefault(name, {})[path] = set_cols[name]
    use = ["chrom", "start", "end", bcol]
    for cols in set_cols.values():
        for c in cols:
            if c not in use:
                use.append(c)
    df = pd.read_csv(path, sep="\t", usecols=use)
    df["abs_bes"] = pd.to_numeric(df[bcol], errors="coerce").abs()
    keep = ["chrom", "start", "end", "abs_bes"]
    for name, cols in set_cols.items():
        tot = np.zeros(len(df), dtype=float)
        for c in cols:
            v = pd.to_numeric(df[c], errors="coerce").abs().to_numpy()
            tot += np.nan_to_num(v, nan=0.0)
        df["set__" + name] = tot
        keep.append("set__" + name)
    return df[keep].copy()


def load_insulation_bm(path):
    """Load HiCExplorer .bm insulation/TAD-score bedGraph-like file.
    Returns DataFrame chrom,start,end,score (mean over numeric score
    columns). Tolerates variable column counts."""
    rows = []
    with open(path) as fh:
        for ln in fh:
            ln = ln.rstrip("\n")
            if not ln or ln.startswith("#") or ln.startswith("track"):
                continue
            parts = ln.split("\t")
            if len(parts) < 4:
                continue
            chrom = parts[0]
            try:
                s = int(parts[1]); e = int(parts[2])
            except ValueError:
                continue
            scores = []
            for v in parts[3:]:
                try:
                    scores.append(float(v))
                except ValueError:
                    pass
            if not scores:
                continue
            rows.append({"chrom": chrom, "start": s, "end": e,
                          "score": float(np.mean(scores))})
    return pd.DataFrame(rows)


def load_bed(path, chrom_filter=None):
    rows = []
    with open(path) as fh:
        for ln in fh:
            ln = ln.rstrip("\n")
            if not ln or ln.startswith("#") or ln.startswith("track"):
                continue
            parts = ln.split("\t")
            if len(parts) < 3:
                continue
            try:
                s = int(parts[1]); e = int(parts[2])
            except ValueError:
                continue
            if chrom_filter and parts[0] != chrom_filter:
                continue
            rows.append({"chrom": parts[0], "start": s, "end": e})
    return pd.DataFrame(rows)


def load_gene_starts(gtf_path, biotype="protein_coding"):
    """Return dict chrom -> sorted np.array of gene TSS (0-based).
    Used for gene-density matching of control regions."""
    import re
    starts = {}
    keep_biotype = (biotype and biotype != "all")
    with open(gtf_path) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            if keep_biotype:
                m = re.search(r'(gene_type|gene_biotype)\s+"([^"]+)"', parts[8])
                if not m or m.group(2) != biotype:
                    continue
            chrom = parts[0]
            try:
                s = int(parts[3]) - 1
            except ValueError:
                continue
            starts.setdefault(chrom, []).append(s)
    for c in starts:
        starts[c] = np.array(sorted(starts[c]))
    return starts


def gene_density(gene_starts, chrom, start, end):
    """Genes per Mb in [chrom, start, end)."""
    arr = gene_starts.get(chrom)
    if arr is None or len(arr) == 0:
        return 0.0
    lo = np.searchsorted(arr, start, side="left")
    hi = np.searchsorted(arr, end, side="left")
    width_mb = max(1e-9, (end - start) / 1e6)
    return (hi - lo) / width_mb


# ----------------------------------------------------------------------
# Per-region computation
# ----------------------------------------------------------------------
def bin_region(bes_df, cool_a, cool_b, insul_a, insul_b,
                chrom, start, end, hic_bin,
                min_distance, max_distance, aggregation,
                contact_agg="sum", score_cols=None):
    """Return per-Hi-C-bin DataFrame with columns:
    bin_start, agg_bes, delta_contact, delta_insulation, and one
    agg_<col> per extra score column (v2 score sets).

    cool_a, cool_b are open cooler.Cooler objects (same binsize=hic_bin).
    """
    import cooler
    region = "{}:{}-{}".format(chrom, start, end)
    mat_a = cool_a.matrix(balance=True, sparse=False).fetch(region)
    mat_b = cool_b.matrix(balance=True, sparse=False).fetch(region)
    bins = cool_a.bins().fetch(region)
    bin_starts = bins["start"].to_numpy()
    n = len(bin_starts)
    if n < 5 or mat_a.shape[0] != n or mat_b.shape[0] != n:
        return None

    # delta_contact per bin: sum/mean over partners within distance range
    min_k = max(1, int(round(min_distance / hic_bin)))
    max_k = max(min_k, int(round(max_distance / hic_bin)))
    diff = np.abs(mat_a - mat_b)
    delta_contact = np.full(n, np.nan)
    for i in range(n):
        lo = max(0, i - max_k)
        hi = min(n, i + max_k + 1)
        partners = []
        for j in range(lo, hi):
            if abs(i - j) >= min_k:
                v = diff[i, j]
                if np.isfinite(v):
                    partners.append(v)
        if partners:
            delta_contact[i] = (np.sum(partners) if contact_agg == "sum"
                                  else np.mean(partners))
        else:
            delta_contact[i] = np.nan

    # delta_insulation per bin (mean score already collapsed)
    def insul_vec(insul_df):
        out = np.full(n, np.nan)
        sub = insul_df[insul_df["chrom"] == chrom]
        if len(sub) == 0:
            return out
        for idx, bs in enumerate(bin_starts):
            hit = sub[(sub["start"] <= bs) & (sub["end"] > bs)]
            if len(hit):
                out[idx] = float(hit["score"].iloc[0])
        return out
    ins_a = insul_vec(insul_a)
    ins_b = insul_vec(insul_b)
    delta_insulation = np.abs(ins_a - ins_b)

    # BES (and score-set) aggregation per Hi-C bin. Vectorised in v2;
    # same semantics as v1: 200 bp bins are assigned by start in
    # [bs, bs + hic_bin); a Hi-C bin with no 200 bp bins gets 0.0.
    sub_bes = bes_df[(bes_df["chrom"] == chrom)
                      & (bes_df["end"] > start)
                      & (bes_df["start"] < end)]
    cols = ["abs_bes"] + list(score_cols or [])
    sb = sub_bes["start"].to_numpy()
    idx = np.searchsorted(bin_starts, sb, side="right") - 1
    ok = (idx >= 0) & (idx < n)
    ok[ok] = sb[ok] < bin_starts[idx[ok]] + hic_bin
    idx = idx[ok]
    order = np.argsort(idx, kind="stable")
    idx_sorted = idx[order]
    uniq, first = np.unique(idx_sorted, return_index=True)
    bounds = list(first[1:]) + [len(idx_sorted)]
    out_cols = {}
    for col in cols:
        vals_all = sub_bes[col].to_numpy()[ok][order]
        agg = np.zeros(n, dtype=float)
        for u, a, b in zip(uniq, first, bounds):
            vals = vals_all[a:b]
            if aggregation == "p95":
                agg[u] = float(np.percentile(vals, 95))
            elif aggregation == "p75":
                agg[u] = float(np.percentile(vals, 75))
            elif aggregation == "max":
                agg[u] = float(np.max(vals))
            elif aggregation == "median":
                agg[u] = float(np.median(vals))
            else:
                agg[u] = float(np.mean(vals))
        out_cols["agg_bes" if col == "abs_bes" else "agg_" + col] = agg

    res = pd.DataFrame({
        "bin_start": bin_starts,
        "agg_bes": out_cols["agg_bes"],
        "delta_contact": delta_contact,
        "delta_insulation": delta_insulation,
    })
    for k, v in out_cols.items():
        if k != "agg_bes":
            res[k] = v
    return res


def region_2x2(binned, contact_col="delta_contact", quantile=0.9,
               x_col="agg_bes"):
    """Within-region top-decile co-localization 2x2 between x_col
    (agg_bes by default) and the chosen Hi-C metric.
    Returns (both, x_only, y_only, neither, n)."""
    if binned is None:
        return None
    df = binned.dropna(subset=[x_col, contact_col])
    n = len(df)
    if n < 10:
        return None
    bes_thr = np.quantile(df[x_col], quantile)
    hic_thr = np.quantile(df[contact_col], quantile)
    x_top = df[x_col] >= bes_thr
    y_top = df[contact_col] >= hic_thr
    both = int((x_top & y_top).sum())
    x_only = int((x_top & ~y_top).sum())
    y_only = int((~x_top & y_top).sum())
    neither = int((~x_top & ~y_top).sum())
    return (both, x_only, y_only, neither, n)


def odds_ratio(both, x_only, y_only, neither):
    """Haldane-corrected odds ratio for a 2x2."""
    a, b, c, d = both, x_only, y_only, neither
    if min(a, b, c, d) == 0:
        a += 0.5; b += 0.5; c += 0.5; d += 0.5
    return (a * d) / (b * c)


# ----------------------------------------------------------------------
# Control region generation
# ----------------------------------------------------------------------
def overlaps_any(chrom, start, end, exclude_df):
    if exclude_df is None or len(exclude_df) == 0:
        return False
    sub = exclude_df[exclude_df["chrom"] == chrom]
    if len(sub) == 0:
        return False
    return bool(((sub["start"] < end) & (sub["end"] > start)).any())


def generate_controls(width, contrast_chrom_sizes, exclude_df,
                        n_controls, rng,
                        gene_starts=None, target_density=None,
                        density_tol=0.5, autosomes_only=True):
    """Generate n_controls random regions of the given width, on
    autosomes, not overlapping exclude_df, optionally matched on gene
    density to within density_tol (relative) of target_density."""
    controls = []
    chroms = list(contrast_chrom_sizes.keys())
    if autosomes_only:
        chroms = [c for c in chroms
                   if c.replace("chr", "").isdigit()]
    if not chroms:
        return controls
    sizes = np.array([contrast_chrom_sizes[c] for c in chroms], dtype=float)
    weights = sizes / sizes.sum()
    max_tries = n_controls * 200
    tries = 0
    while len(controls) < n_controls and tries < max_tries:
        tries += 1
        c = rng.choice(chroms, p=weights)
        clen = contrast_chrom_sizes[c]
        if clen <= width:
            continue
        s = int(rng.integers(0, clen - width))
        e = s + width
        if overlaps_any(c, s, e, exclude_df):
            continue
        if gene_starts is not None and target_density is not None:
            d = gene_density(gene_starts, c, s, e)
            if target_density <= 0:
                # match low-density target: accept only low-density ctrl
                if d > 0.5:
                    continue
            else:
                rel = abs(d - target_density) / target_density
                if rel > density_tol:
                    continue
        controls.append((c, s, e))
    return controls


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def holm(pvals):
    """Holm step-down adjusted p-values (same order as input)."""
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    adj = [0.0] * m
    running = 0.0
    for rank, i in enumerate(order):
        val = min(1.0, (m - rank) * pvals[i])
        running = max(running, val)
        adj[i] = running
    return adj


def parse_score_sets(specs):
    """NAME=tok,tok,... -> list of (NAME, [tok, ...]) in given order."""
    out = []
    seen = set()
    for spec in specs or []:
        if "=" not in spec:
            sys.exit("invalid --score-sets spec (need NAME=tok,tok): " + spec)
        k, v = spec.split("=", 1)
        k = k.strip()
        if not k or k == "BES" or k in seen:
            sys.exit("score-set name must be unique, non-empty, not 'BES': "
                     + spec)
        seen.add(k)
        out.append((k, [t for t in v.split(",") if t.strip()]))
    return out


def parse_named(args_list, what):
    out = {}
    for spec in args_list:
        if "=" not in spec:
            sys.exit("invalid --{} spec (need NAME=value): {}".format(what, spec))
        k, v = spec.split("=", 1)
        out[k.strip()] = v.strip()
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--diffs", nargs="+", required=True,
                    help="CONTRAST=path BEARING diff TSVs")
    ap.add_argument("--cool-a", required=True,
                    help="shared condition-A .cool (e.g. DN)")
    ap.add_argument("--cool-b", nargs="+", required=True,
                    help="CONTRAST=path condition-B .cool per contrast")
    ap.add_argument("--insul-a", required=True,
                    help="shared condition-A insulation .bm")
    ap.add_argument("--insul-b", nargs="+", required=True,
                    help="CONTRAST=path condition-B insulation .bm per contrast")
    ap.add_argument("--targets", required=True,
                    help="TSV: contrast, chrom, start, end, name (header)")
    ap.add_argument("--hic-bin", type=int, default=10000)
    ap.add_argument("--min-distance", type=int, default=50000)
    ap.add_argument("--max-distance", type=int, default=500000)
    ap.add_argument("--aggregation", default="p95",
                    choices=["p95", "p75", "max", "median", "mean"])
    ap.add_argument("--contact-metric", default="delta_contact",
                    choices=["delta_contact", "delta_insulation"])
    ap.add_argument("--contact-agg", default="sum", choices=["sum", "mean"])
    ap.add_argument("--quantile", type=float, default=0.9,
                    help="top-quantile for co-localization (default 0.9)")
    ap.add_argument("--n-controls", type=int, default=200,
                    help="control regions generated per target")
    ap.add_argument("--n-panels", type=int, default=1000,
                    help="random control panels for the null")
    ap.add_argument("--blacklist", default=None, help="BED to exclude")
    ap.add_argument("--gtf", default=None,
                    help="GTF for gene-density matching of controls")
    ap.add_argument("--match-gene-density", action="store_true")
    ap.add_argument("--density-tol", type=float, default=0.5)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--score-sets", nargs="+", default=None,
                    help="v2: NAME=tok,tok ... pre-registered track "
                         "combinations (track name suffix or 1-based "
                         "kl_* index). Holm-corrected across sets.")
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()
    score_sets = parse_score_sets(args.score_sets)
    # score list: (label, column in bin_region output, column in bes_df)
    scores = [("BES", "agg_bes", "abs_bes")]
    for nm, _ in score_sets:
        scores.append((nm, "agg_set__" + nm, "set__" + nm))
    extra_cols = [c for _, _, c in scores[1:]]
    resolved_log = {}

    import cooler

    diffs = parse_named(args.diffs, "diffs")
    cool_b_paths = parse_named(args.cool_b, "cool-b")
    insul_b_paths = parse_named(args.insul_b, "insul-b")

    # Validate contrast keys line up
    targets = pd.read_csv(args.targets, sep="\t")
    need_cols = {"contrast", "chrom", "start", "end", "name"}
    if not need_cols.issubset(targets.columns):
        sys.exit("targets TSV needs columns: " + ",".join(sorted(need_cols)))
    contrasts_used = sorted(targets["contrast"].unique())
    for ct in contrasts_used:
        for d, label in [(diffs, "diffs"), (cool_b_paths, "cool-b"),
                          (insul_b_paths, "insul-b")]:
            if ct not in d:
                sys.exit("contrast {} in targets but missing from --{}"
                         .format(ct, label))

    rng = np.random.default_rng(args.seed)

    # Load shared A resources
    print("Loading condition-A cool + insulation ...", flush=True)
    cool_a = cooler.Cooler(args.cool_a)
    insul_a = load_insulation_bm(args.insul_a)
    chrom_sizes_a = dict(cool_a.chromsizes)

    # Per-contrast resources (lazy load + cache)
    diff_cache = {}
    coolb_cache = {}
    insb_cache = {}
    def get_contrast(ct):
        if ct not in diff_cache:
            print("Loading diff for", ct, flush=True)
            diff_cache[ct] = load_bes(diffs[ct], score_sets, resolved_log)
            for nm, _ in score_sets:
                print("  score set {:<22} -> {}".format(
                    nm, ",".join(resolved_log[nm][diffs[ct]])), flush=True)
        if ct not in coolb_cache:
            coolb_cache[ct] = cooler.Cooler(cool_b_paths[ct])
        if ct not in insb_cache:
            insb_cache[ct] = load_insulation_bm(insul_b_paths[ct])
        return diff_cache[ct], coolb_cache[ct], insb_cache[ct]

    # Blacklist + targets become the exclusion set for controls
    blacklist = load_bed(args.blacklist) if args.blacklist else pd.DataFrame(
        columns=["chrom", "start", "end"])
    target_regions_df = targets[["chrom", "start", "end"]].copy()
    exclude_df = pd.concat([blacklist, target_regions_df], ignore_index=True)

    gene_starts = None
    if args.match_gene_density:
        if not args.gtf:
            sys.exit("--match-gene-density requires --gtf")
        print("Loading GTF gene starts for density matching ...", flush=True)
        gene_starts = load_gene_starts(args.gtf)

    # ---- Compute target 2x2s (all scores) ----
    print("\n=== TARGET regions ===")
    target_rows = []                      # BES rows (v1-compatible file)
    set_rows = []                         # long format, every score
    pooled = {lab: np.zeros(4, dtype=float) for lab, _, _ in scores}
    target_controls = {}                  # name -> list of {label: 2x2}
    for _, t in targets.iterrows():
        ct = t["contrast"]; chrom = t["chrom"]
        s = int(t["start"]); e = int(t["end"]); name = t["name"]
        width = e - s
        bes_df, cool_b, insul_b = get_contrast(ct)
        binned = bin_region(bes_df, cool_a, cool_b, insul_a, insul_b,
                              chrom, s, e, args.hic_bin,
                              args.min_distance, args.max_distance,
                              args.aggregation, args.contact_agg,
                              score_cols=extra_cols)
        tabs = {lab: region_2x2(binned, args.contact_metric, args.quantile,
                                x_col=xc) for lab, xc, _ in scores}
        if any(v is None for v in tabs.values()):
            print("  SKIP {} ({}): too few valid bins".format(name, ct))
            continue
        for lab, _, _ in scores:
            both, xo, yo, ne, nb = tabs[lab]
            orr = odds_ratio(both, xo, yo, ne)
            print("  {:<14} {:<11} {:<22} n={:<4} both={:<3} OR={:.2f}".format(
                name, ct, lab, nb, both, orr))
            row = {"name": name, "contrast": ct, "chrom": chrom,
                   "start": s, "end": e, "width": width,
                   "n_bins": nb, "both_top": both,
                   "bes_only": xo, "contact_only": yo,
                   "neither": ne, "odds_ratio": orr}
            if lab == "BES":
                target_rows.append(row)
            srow = {"score_set": lab}
            srow.update(row)
            set_rows.append(srow)
            pooled[lab] += np.array([both, xo, yo, ne], dtype=float)

        # Matched controls for this target; every score is computed on
        # the same control region so the null panels are paired.
        tgt_density = None
        if gene_starts is not None:
            tgt_density = gene_density(gene_starts, chrom, s, e)
        ctrl_regions = generate_controls(
            width, chrom_sizes_a, exclude_df, args.n_controls, rng,
            gene_starts=gene_starts, target_density=tgt_density,
            density_tol=args.density_tol)
        ctrl_tabs = []
        for (cc, cs, ce) in ctrl_regions:
            cbin = bin_region(bes_df, cool_a, cool_b, insul_a, insul_b,
                               cc, cs, ce, args.hic_bin,
                               args.min_distance, args.max_distance,
                               args.aggregation, args.contact_agg,
                               score_cols=extra_cols)
            ctab = {lab: region_2x2(cbin, args.contact_metric,
                                    args.quantile, x_col=xc)
                    for lab, xc, _ in scores}
            if all(v is not None for v in ctab.values()):
                ctrl_tabs.append({k: v[:4] for k, v in ctab.items()})
        target_controls[name] = ctrl_tabs
        print("    generated {} usable controls".format(len(ctrl_tabs)))

    if not target_rows:
        sys.exit("no usable target regions")

    target_df = pd.DataFrame(target_rows)
    pooled_or = {lab: odds_ratio(*pooled[lab].astype(int))
                 for lab, _, _ in scores}
    for lab, _, _ in scores:
        print("\nPooled target 2x2 [{}]: both={:.0f} bes_only={:.0f} "
              "contact_only={:.0f} neither={:.0f}".format(lab, *pooled[lab]))
        print("Pooled target OR [{}] = {:.3f}".format(lab, pooled_or[lab]))

    # ---- Null: one control per target, pooled, paired across scores ----
    print("\n=== NULL: {} control panels ===".format(args.n_panels))
    null_ors = {lab: [] for lab, _, _ in scores}
    names_with_ctrls = [r["name"] for r in target_rows
                         if target_controls.get(r["name"])]
    if len(names_with_ctrls) < len(target_rows):
        print("  WARNING: {} target(s) have no usable controls and are "
              "left out of the null panels".format(
                  len(target_rows) - len(names_with_ctrls)))
    for _ in range(args.n_panels):
        picks = [target_controls[nm][rng.integers(0, len(target_controls[nm]))]
                 for nm in names_with_ctrls]
        if not picks:
            break
        for lab, _, _ in scores:
            pp = np.zeros(4, dtype=float)
            for pk in picks:
                pp += np.array(pk[lab], dtype=float)
            if pp.sum() > 0:
                null_ors[lab].append(odds_ratio(*pp.astype(int)))
    stats = {}
    for lab, _, _ in scores:
        arr = np.array(null_ors[lab])
        if len(arr) == 0:
            sys.exit("no null panels could be built (insufficient controls)")
        stats[lab] = {
            "pooled_OR": pooled_or[lab],
            "null_n_panels": int(len(arr)),
            "null_OR_median": float(np.median(arr)),
            "null_OR_p95": float(np.percentile(arr, 95)),
            "empirical_p": float((np.sum(arr >= pooled_or[lab]) + 1)
                                 / (len(arr) + 1)),
            "pooled_2x2": [int(x) for x in pooled[lab]],
        }
        print("[{}] null OR median={:.3f} 95th={:.3f}  empirical p={:.4f}"
              .format(lab, stats[lab]["null_OR_median"],
                      stats[lab]["null_OR_p95"], stats[lab]["empirical_p"]))
    set_labels = [lab for lab, _, _ in scores[1:]]
    if set_labels:
        adj = holm([stats[l]["empirical_p"] for l in set_labels])
        for l, a in zip(set_labels, adj):
            stats[l]["holm_p"] = a
        print("\nHolm across {} pre-registered set(s): ".format(len(set_labels))
              + ", ".join("{}={:.4f}".format(l, stats[l]["holm_p"])
                          for l in set_labels))
        print("(Minimum achievable empirical p = {:.4f} at --n-panels {}; "
              "with Holm over {} sets the smallest adjusted p is {:.4f})"
              .format(1.0 / (args.n_panels + 1), args.n_panels,
                      len(set_labels),
                      min(1.0, len(set_labels) / (args.n_panels + 1))))

    # ---- Outputs ----
    target_df.to_csv(args.out_prefix + ".target_regions.tsv",
                      sep="\t", index=False, float_format="%.4f")
    b = stats["BES"]
    summary = {
        "script_version": "2.0.0",
        "aggregation": args.aggregation,
        "contact_metric": args.contact_metric,
        "hic_bin": args.hic_bin,
        "quantile": args.quantile,
        "n_targets": len(target_rows),
        "pooled_target_2x2": {
            "both_top": b["pooled_2x2"][0],
            "bes_only": b["pooled_2x2"][1],
            "contact_only": b["pooled_2x2"][2],
            "neither": b["pooled_2x2"][3],
        },
        "pooled_target_OR": b["pooled_OR"],
        "null_n_panels": b["null_n_panels"],
        "null_OR_median": b["null_OR_median"],
        "null_OR_p95": b["null_OR_p95"],
        "empirical_p": b["empirical_p"],
        "match_gene_density": bool(args.match_gene_density),
        "n_controls_per_target": args.n_controls,
    }
    if set_labels:
        summary["score_sets"] = {
            l: dict(stats[l], tracks=resolved_log.get(l, {}))
            for l in set_labels}
    with open(args.out_prefix + ".summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    print("\nWrote", args.out_prefix + ".target_regions.tsv")
    print("Wrote", args.out_prefix + ".summary.json")
    if set_labels:
        pd.DataFrame(set_rows).to_csv(args.out_prefix + ".score_sets.tsv",
                                      sep="\t", index=False,
                                      float_format="%.4f")
        srows = []
        for lab, _, _ in scores:
            st = stats[lab]
            tracks = sorted({c for v in resolved_log.get(lab, {}).values()
                             for c in v}) if lab != "BES" else ["all"]
            srows.append({
                "score_set": lab, "tracks": ",".join(tracks),
                "in_holm_family": lab != "BES",
                "both_top": st["pooled_2x2"][0],
                "x_only": st["pooled_2x2"][1],
                "contact_only": st["pooled_2x2"][2],
                "neither": st["pooled_2x2"][3],
                "pooled_OR": st["pooled_OR"],
                "null_OR_median": st["null_OR_median"],
                "null_OR_p95": st["null_OR_p95"],
                "empirical_p": st["empirical_p"],
                "holm_p": st.get("holm_p", float("nan")),
                "null_n_panels": st["null_n_panels"],
            })
        pd.DataFrame(srows).to_csv(args.out_prefix
                                   + ".score_sets_summary.tsv",
                                   sep="\t", index=False,
                                   float_format="%.4g")
        print("Wrote", args.out_prefix + ".score_sets.tsv")
        print("Wrote", args.out_prefix + ".score_sets_summary.tsv")

    # Plot null distribution(s) with observed line
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        k = len(scores)
        fig, axes = plt.subplots(1, k, figsize=(4.2 * k, 3.8), squeeze=False)
        for ax, (lab, _, _) in zip(axes[0], scores):
            st = stats[lab]
            ax.hist(null_ors[lab], bins=40, color="#9bb", edgecolor="white")
            ax.axvline(st["pooled_OR"], color="#e63946", lw=2,
                       label="observed OR = {:.2f}".format(st["pooled_OR"]))
            ax.axvline(st["null_OR_median"], color="gray", lw=1, ls="--",
                       label="null median = {:.2f}".format(
                           st["null_OR_median"]))
            ttl = "{}\nemp. p = {:.3f}".format(lab, st["empirical_p"])
            if "holm_p" in st:
                ttl += "  Holm = {:.3f}".format(st["holm_p"])
            ax.set_title(ttl, fontsize=9)
            ax.set_xlabel("Pooled top-decile OR (x {})".format(
                args.contact_metric), fontsize=8)
            ax.legend(fontsize=7)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
        axes[0][0].set_ylabel("control panels")
        fig.suptitle("Cross-locus score-Hi-C co-localization, AR loci vs "
                     "{}matched controls".format(
                         "gene-density-" if args.match_gene_density else ""),
                     fontsize=10)
        plt.tight_layout()
        plt.savefig(args.out_prefix + ".null_distribution.pdf",
                     bbox_inches="tight")
        plt.close()
        print("Wrote", args.out_prefix + ".null_distribution.pdf")
    except Exception as e:
        print("plot skipped:", e)


if __name__ == "__main__":
    main()
