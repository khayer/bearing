#!/usr/bin/env python3
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
def load_bes(path):
    """Load a BEARING diff TSV. Returns DataFrame with chrom,start,end,
    abs_bes. Uses bearing_score column; falls back to |bearing_score|."""
    df = pd.read_csv(path, sep="\t")
    # Expected columns include chrom,start,end,bearing_score
    needed = {"chrom", "start", "end"}
    if not needed.issubset(df.columns):
        sys.exit("BES TSV {} missing chrom/start/end".format(path))
    if "bearing_score" in df.columns:
        df["abs_bes"] = df["bearing_score"].abs()
    elif "bearing_score_tested" in df.columns:
        df["abs_bes"] = df["bearing_score_tested"].abs()
    else:
        sys.exit("BES TSV {} has no bearing_score column".format(path))
    return df[["chrom", "start", "end", "abs_bes"]].copy()


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
                contact_agg="sum"):
    """Return per-Hi-C-bin DataFrame with columns:
    bin_start, agg_bes, delta_contact, delta_insulation.

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

    # BES aggregation per Hi-C bin
    sub_bes = bes_df[(bes_df["chrom"] == chrom)
                      & (bes_df["end"] > start)
                      & (bes_df["start"] < end)]
    agg_bes = np.full(n, np.nan)
    for idx, bs in enumerate(bin_starts):
        be = bs + hic_bin
        vals = sub_bes[(sub_bes["start"] >= bs)
                        & (sub_bes["start"] < be)]["abs_bes"].to_numpy()
        if len(vals) == 0:
            agg_bes[idx] = 0.0
            continue
        if aggregation == "p95":
            agg_bes[idx] = float(np.percentile(vals, 95))
        elif aggregation == "p75":
            agg_bes[idx] = float(np.percentile(vals, 75))
        elif aggregation == "max":
            agg_bes[idx] = float(np.max(vals))
        elif aggregation == "median":
            agg_bes[idx] = float(np.median(vals))
        else:
            agg_bes[idx] = float(np.mean(vals))

    return pd.DataFrame({
        "bin_start": bin_starts,
        "agg_bes": agg_bes,
        "delta_contact": delta_contact,
        "delta_insulation": delta_insulation,
    })


def region_2x2(binned, contact_col="delta_contact", quantile=0.9):
    """Within-region top-decile co-localization 2x2 between agg_bes and
    the chosen Hi-C metric. Returns (both, x_only, y_only, neither, n)."""
    if binned is None:
        return None
    df = binned.dropna(subset=["agg_bes", contact_col])
    n = len(df)
    if n < 10:
        return None
    bes_thr = np.quantile(df["agg_bes"], quantile)
    hic_thr = np.quantile(df[contact_col], quantile)
    x_top = df["agg_bes"] >= bes_thr
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
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

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
            diff_cache[ct] = load_bes(diffs[ct])
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

    # ---- Compute target 2x2s ----
    print("\n=== TARGET regions ===")
    target_rows = []
    target_2x2_sum = np.zeros(4, dtype=float)
    target_controls = {}  # name -> list of control 2x2 tuples
    for _, t in targets.iterrows():
        ct = t["contrast"]; chrom = t["chrom"]
        s = int(t["start"]); e = int(t["end"]); name = t["name"]
        width = e - s
        bes_df, cool_b, insul_b = get_contrast(ct)
        binned = bin_region(bes_df, cool_a, cool_b, insul_a, insul_b,
                              chrom, s, e, args.hic_bin,
                              args.min_distance, args.max_distance,
                              args.aggregation, args.contact_agg)
        tab = region_2x2(binned, args.contact_metric, args.quantile)
        if tab is None:
            print("  SKIP {} ({}): too few valid bins".format(name, ct))
            continue
        both, xo, yo, ne, nb = tab
        orr = odds_ratio(both, xo, yo, ne)
        print("  {:<14} {:<11} n={:<4} both={:<3} OR={:.2f}".format(
            name, ct, nb, both, orr))
        target_rows.append({"name": name, "contrast": ct, "chrom": chrom,
                             "start": s, "end": e, "width": width,
                             "n_bins": nb, "both_top": both,
                             "bes_only": xo, "contact_only": yo,
                             "neither": ne, "odds_ratio": orr})
        target_2x2_sum += np.array([both, xo, yo, ne], dtype=float)

        # Generate matched controls for this target now
        tgt_density = None
        if gene_starts is not None:
            tgt_density = gene_density(gene_starts, chrom, s, e)
        # Use the B-condition chrom sizes intersect A (same assembly)
        ctrl_regions = generate_controls(
            width, chrom_sizes_a, exclude_df, args.n_controls, rng,
            gene_starts=gene_starts, target_density=tgt_density,
            density_tol=args.density_tol)
        ctrl_tabs = []
        for (cc, cs, ce) in ctrl_regions:
            cbin = bin_region(bes_df, cool_a, cool_b, insul_a, insul_b,
                               cc, cs, ce, args.hic_bin,
                               args.min_distance, args.max_distance,
                               args.aggregation, args.contact_agg)
            ctab = region_2x2(cbin, args.contact_metric, args.quantile)
            if ctab is not None:
                ctrl_tabs.append(ctab[:4])
        target_controls[name] = ctrl_tabs
        print("    generated {} usable controls".format(len(ctrl_tabs)))

    if not target_rows:
        sys.exit("no usable target regions")

    target_df = pd.DataFrame(target_rows)
    pooled_target_or = odds_ratio(*target_2x2_sum.astype(int))
    print("\nPooled target 2x2: both={:.0f} bes_only={:.0f} "
          "contact_only={:.0f} neither={:.0f}".format(*target_2x2_sum))
    print("Pooled target OR = {:.3f}".format(pooled_target_or))

    # ---- Build null: one control per target, pool, repeat ----
    print("\n=== NULL: {} control panels ===".format(args.n_panels))
    null_ors = []
    names_with_ctrls = [r["name"] for r in target_rows
                         if target_controls.get(r["name"])]
    for _ in range(args.n_panels):
        pooled = np.zeros(4, dtype=float)
        ok = True
        for nm in names_with_ctrls:
            tabs = target_controls[nm]
            if not tabs:
                ok = False
                break
            pick = tabs[rng.integers(0, len(tabs))]
            pooled += np.array(pick, dtype=float)
        if ok and pooled.sum() > 0:
            null_ors.append(odds_ratio(*pooled.astype(int)))
    null_ors = np.array(null_ors)
    if len(null_ors) == 0:
        sys.exit("no null panels could be built (insufficient controls)")

    emp_p = float((np.sum(null_ors >= pooled_target_or) + 1)
                   / (len(null_ors) + 1))
    null_median = float(np.median(null_ors))
    null_p95 = float(np.percentile(null_ors, 95))

    print("Null OR: median={:.3f}  95th pct={:.3f}".format(
        null_median, null_p95))
    print("Empirical p (target OR >= null) = {:.4f}".format(emp_p))

    # ---- Outputs ----
    target_df.to_csv(args.out_prefix + ".target_regions.tsv",
                      sep="\t", index=False, float_format="%.4f")
    summary = {
        "aggregation": args.aggregation,
        "contact_metric": args.contact_metric,
        "hic_bin": args.hic_bin,
        "quantile": args.quantile,
        "n_targets": len(target_rows),
        "pooled_target_2x2": {
            "both_top": int(target_2x2_sum[0]),
            "bes_only": int(target_2x2_sum[1]),
            "contact_only": int(target_2x2_sum[2]),
            "neither": int(target_2x2_sum[3]),
        },
        "pooled_target_OR": pooled_target_or,
        "null_n_panels": int(len(null_ors)),
        "null_OR_median": null_median,
        "null_OR_p95": null_p95,
        "empirical_p": emp_p,
        "match_gene_density": bool(args.match_gene_density),
        "n_controls_per_target": args.n_controls,
    }
    with open(args.out_prefix + ".summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    print("\nWrote", args.out_prefix + ".target_regions.tsv")
    print("Wrote", args.out_prefix + ".summary.json")

    # Plot null distribution with observed line
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.hist(null_ors, bins=40, color="#9bb", edgecolor="white")
        ax.axvline(pooled_target_or, color="#e63946", lw=2,
                     label="AR loci pooled OR = {:.2f}".format(pooled_target_or))
        ax.axvline(null_median, color="gray", lw=1, ls="--",
                     label="null median = {:.2f}".format(null_median))
        ax.set_xlabel("Pooled top-decile OR (BES x {})".format(args.contact_metric))
        ax.set_ylabel("control panels")
        ttl = ("Cross-locus BES-Hi-C co-localization\n"
               "AR loci vs {}matched controls (empirical p = {:.3f})".format(
                   "gene-density-" if args.match_gene_density else "",
                   emp_p))
        ax.set_title(ttl, fontsize=10)
        ax.legend(fontsize=8)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        plt.tight_layout()
        plt.savefig(args.out_prefix + ".null_distribution.pdf",
                     bbox_inches="tight")
        plt.close()
        print("Wrote", args.out_prefix + ".null_distribution.pdf")
    except Exception as e:
        print("plot skipped:", e)


if __name__ == "__main__":
    main()
