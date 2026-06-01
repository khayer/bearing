#!/usr/bin/env python3
"""
compartment_analysis.py

Compare A/B compartment assignments (PC1 eigenvector tracks) across
conditions to test whether a chromosomal region undergoes compartment
switching - the megabase-scale architectural signature relevant to
the Tcrb-Cntnap2-Gpnmb extended-domain hypothesis.

INPUTS:
  - One PC1 bigwig per condition (typical naming:
    merged_<COND>_bs_100000_pca1.bw). 100 kb resolution recommended
    for compartments - finer scales are noisy.
  - Optional features BED for labels and TADs BED files for overlay.

PC1 SIGN HANDLING:
  PC1 sign is arbitrary per chromosome unless oriented by gene density
  or GC content. This script offers three modes:

    --orient-with-gtf <GTF>  (PREFERRED)
        Counts protein-coding gene starts per PC1 bin. For each
        chromosome, computes Spearman correlation between gene
        density and raw PC1. If negative, flips PC1 sign for that
        chromosome - so HIGH GENE DENSITY = POSITIVE PC1 = A.
        After orientation, sign IS interpretable as A vs B.

    --normalize-sign
        Flips per-chromosome so the chromosomal mean is positive.
        Only useful for visual consistency; does NOT establish A vs B.

    (default: no flipping)
        Use the bigwig sign as-is. Only the switching analysis is
        biologically meaningful in this mode.

  Compartment SWITCHES across conditions are biologically meaningful
  regardless of orientation method.

OUTPUTS:
  - <prefix>.pc1_per_feature.tsv : PC1 value at each feature per cond
  - <prefix>.switch_matrix.tsv   : per-bin condition x condition
                                    switching count for the region
  - <prefix>.compartment_landscape.pdf : combined figure showing
       PC1 tracks per condition + optional TAD overlay + features

Usage:
  python compartment_analysis.py \\
      --pc1 DN=../comp/merged_DN_bs_100000_pca1.bw \\
            DP=../comp/merged_DP_bs_100000_pca1.bw \\
            ProB=../comp/merged_ProB_bs_100000_pca1.bw \\
            EbKO=../comp/merged_EbKO_bs_100000_pca1.bw \\
            S3T3=../comp/merged_S3T3_bs_100000_pca1.bw \\
            dV1P=../comp/merged_dV1P_bs_100000_pca1.bw \\
            dV1CTCF=../comp/merged_dV1CTCF_bs_100000_pca1.bw \\
      --region chr6:39000000-55000000 \\
      --features-bed tcrb_extension_features.bed \\
      --tad-dir ../tads/ --tad-resolution 100000 \\
      --normalize-sign \\
      --out-prefix tcrb_compartments
"""
import argparse
import os
import sys
import numpy as np
import pandas as pd


def read_bigwig_region(path, chrom, start, end):
    """Read a bigwig in [chrom, start, end). Tries pyBigWig first,
    falls back to bigwigtobedgraph + parse if not available."""
    try:
        import pyBigWig
        bw = pyBigWig.open(path)
        # Returns one value per bp; we want per-bin values, so use
        # intervals which gives (start, end, value) tuples
        intervals = bw.intervals(chrom, start, end) or []
        bw.close()
        if not intervals:
            return pd.DataFrame(columns=["start", "end", "value"])
        df = pd.DataFrame(intervals, columns=["start", "end", "value"])
        return df
    except ImportError:
        # Fall back: try bigWigToBedGraph if on PATH
        import subprocess
        import tempfile
        with tempfile.NamedTemporaryFile(suffix=".bg", delete=False) as tf:
            tmp = tf.name
        try:
            subprocess.run(
                ["bigWigToBedGraph",
                 "-chrom={}".format(chrom),
                 "-start={}".format(start),
                 "-end={}".format(end),
                 path, tmp],
                check=True, capture_output=True)
        except (FileNotFoundError, subprocess.CalledProcessError) as e:
            sys.exit("Cannot read bigwig {}: pyBigWig not installed AND "
                     "bigWigToBedGraph not available. Install pyBigWig "
                     "with `pip install pyBigWig` or put bigWigToBedGraph "
                     "on PATH. Error: {}".format(path, e))
        df = pd.read_csv(tmp, sep="\t", header=None,
                         names=["chrom", "start", "end", "value"])
        os.unlink(tmp)
        return df[["start", "end", "value"]]


def read_bigwig_chrom(path, chrom):
    """Read entire chromosome from bigwig (for sign normalization)."""
    try:
        import pyBigWig
        bw = pyBigWig.open(path)
        chroms = bw.chroms()
        if chrom not in chroms:
            bw.close()
            return pd.DataFrame(columns=["start", "end", "value"])
        chrom_len = chroms[chrom]
        intervals = bw.intervals(chrom, 0, chrom_len) or []
        bw.close()
        if not intervals:
            return pd.DataFrame(columns=["start", "end", "value"])
        return pd.DataFrame(intervals, columns=["start", "end", "value"])
    except ImportError:
        return read_bigwig_region(path, chrom, 0, 250000000)


def load_features_bed(path, chrom_filter=None):
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
            name = parts[3] if len(parts) >= 4 else ""
            rows.append({"chrom": parts[0], "start": s, "end": e,
                          "name": name})
    return pd.DataFrame(rows)


def load_tads_bed(path, chrom_filter=None):
    """Same loader as the TAD script (kept local for portability)."""
    return load_features_bed(path, chrom_filter=chrom_filter)


def value_at_position(df_pc1, pos):
    """Return PC1 value at a bp position (interpolated from binned bw)."""
    if len(df_pc1) == 0:
        return float("nan")
    hits = df_pc1[(df_pc1["start"] <= pos) & (df_pc1["end"] > pos)]
    if len(hits) == 0:
        return float("nan")
    return float(hits["value"].iloc[0])


def load_gene_counts_per_bin(gtf_path, chrom, bin_starts, bin_ends,
                                biotype="protein_coding"):
    """Count protein-coding gene start positions falling in each PC1 bin.

    Uses gene_type or gene_biotype attribute. Handles chr vs no-chr
    prefix. Uses awk for fast region extraction; falls back to Python.
    Returns array of len(bin_starts) with gene counts.
    """
    import subprocess
    import re

    # Handle chr-prefix variants
    chrom_variants = [chrom]
    if chrom.startswith("chr"):
        chrom_variants.append(chrom[3:])
    else:
        chrom_variants.append("chr" + chrom)

    keep_biotype = (biotype and biotype != "all")
    gene_starts = []
    for chrom_try in chrom_variants:
        awk_prog = ('$1 == "' + chrom_try + '" && $3 == "gene" { print $4 "\t" $9 }')
        try:
            proc = subprocess.run(
                ["awk", "-F", "\t", awk_prog, gtf_path],
                capture_output=True, text=True, check=False)
            if proc.stdout.strip():
                for ln in proc.stdout.strip().split("\n"):
                    parts = ln.split("\t", 1)
                    if len(parts) < 2:
                        continue
                    try:
                        s = int(parts[0])
                    except ValueError:
                        continue
                    attrs = parts[1]
                    if keep_biotype:
                        m = re.search(r'(gene_type|gene_biotype)\s+"([^"]+)"',
                                       attrs)
                        if not m or m.group(2) != biotype:
                            continue
                    gene_starts.append(s - 1)  # to 0-based
                break
        except FileNotFoundError:
            # Pure Python fallback
            with open(gtf_path) as fh:
                for ln in fh:
                    if ln.startswith("#"):
                        continue
                    parts = ln.rstrip("\n").split("\t")
                    if len(parts) < 9 or parts[0] != chrom_try or parts[2] != "gene":
                        continue
                    try:
                        s = int(parts[3]) - 1
                    except ValueError:
                        continue
                    if keep_biotype:
                        m = re.search(r'(gene_type|gene_biotype)\s+"([^"]+)"',
                                       parts[8])
                        if not m or m.group(2) != biotype:
                            continue
                    gene_starts.append(s)
            if gene_starts:
                break

    if not gene_starts:
        return np.zeros(len(bin_starts), dtype=int)

    gene_starts = np.array(sorted(gene_starts))
    # For each bin, count how many gene_starts fall in [bin_start, bin_end)
    counts = np.zeros(len(bin_starts), dtype=int)
    for i, (s, e) in enumerate(zip(bin_starts, bin_ends)):
        lo = np.searchsorted(gene_starts, s, side="left")
        hi = np.searchsorted(gene_starts, e, side="left")
        counts[i] = hi - lo
    return counts


def orient_pc1_by_gene_density(df_chrom_pc1, gene_counts):
    """Decide whether to flip PC1 for this chromosome based on gene
    density correlation. Returns (flip_decision: bool, spearman_rho: float).

    Convention: positive correlation between gene density and PC1
    means PC1 already follows A=positive convention; no flip needed.
    Negative correlation means we should flip.
    """
    if len(df_chrom_pc1) == 0:
        return False, float("nan")
    # Align lengths defensively
    n = min(len(df_chrom_pc1), len(gene_counts))
    pc1_vals = df_chrom_pc1["value"].values[:n]
    gene_vals = gene_counts[:n].astype(float)
    # Drop NaNs from PC1
    mask = np.isfinite(pc1_vals)
    if mask.sum() < 10:
        return False, float("nan")
    pc1_clean = pc1_vals[mask]
    gene_clean = gene_vals[mask]
    # Prefer Spearman (rank-based, robust to outliers), fall back to Pearson
    try:
        from scipy.stats import spearmanr
        rho, _ = spearmanr(pc1_clean, gene_clean)
        method = "Spearman"
    except ImportError:
        # numpy Pearson - sufficient for sign determination
        if np.std(pc1_clean) == 0 or np.std(gene_clean) == 0:
            return False, float("nan")
        rho = float(np.corrcoef(pc1_clean, gene_clean)[0, 1])
        method = "Pearson (scipy unavailable)"
    if not np.isfinite(rho):
        return False, float("nan")
    return (rho < 0), float(rho)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pc1", required=True, nargs="+",
                    help="CONDITION=PATH entries for each PC1 bigwig")
    ap.add_argument("--region", required=True,
                    help="chrom:start-end to analyze")
    ap.add_argument("--features-bed", default=None,
                    help="optional BED for feature labels")
    ap.add_argument("--tad-dir", default=None,
                    help="optional directory with TAD *_domains.bed files "
                         "to overlay")
    ap.add_argument("--tad-resolution", type=int, default=100000,
                    help="TAD binsize for overlay (default 100000)")
    ap.add_argument("--normalize-sign", action="store_true",
                    help="flip per-chrom PC1 so the chromosomal mean is "
                         "positive (for cross-condition display only; "
                         "does not establish A vs B identity)")
    ap.add_argument("--orient-with-gtf", default=None,
                    help="GTF file with gene annotations. If provided, "
                         "each chromosome's PC1 sign is oriented so that "
                         "gene-dense regions have POSITIVE PC1 (the "
                         "standard A=positive, B=negative convention). "
                         "Uses Spearman correlation between per-bin gene "
                         "count and raw PC1 to decide flip. Overrides "
                         "--normalize-sign when both are given.")
    ap.add_argument("--orient-biotype", default="protein_coding",
                    help="biotype filter for GTF (default protein_coding). "
                         "Use 'all' to disable.")
    ap.add_argument("--switch-tolerance", type=float, default=0.0,
                    help="treat PC1 values within +/- this of 0 as "
                         "compartmentally ambiguous (no switch called "
                         "across them). Default 0.0.")
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    # Parse region
    chrom, rng = args.region.split(":")
    r_start, r_end = map(int, rng.split("-"))

    # Parse PC1 paths
    cond_paths = []
    for spec in args.pc1:
        if "=" not in spec:
            sys.exit("invalid --pc1 spec: " + spec)
        name, path = spec.split("=", 1)
        cond_paths.append((name.strip(), path.strip()))

    # Load PC1 per condition for the displayed region
    pc1_region = {}
    pc1_chrom_mean = {}
    for cond, path in cond_paths:
        if not os.path.isfile(path):
            print("WARNING: {} not found, skipping {}".format(path, cond))
            continue
        print("Loading", cond, "PC1 from", path, flush=True)
        df_region = read_bigwig_region(path, chrom, r_start, r_end)

        # Decide flip orientation
        # Priority: --orient-with-gtf > --normalize-sign > no flip
        flip = False
        if args.orient_with_gtf:
            df_chrom = read_bigwig_chrom(path, chrom)
            if len(df_chrom) == 0:
                print("  WARNING: no PC1 data on {} for orientation".format(chrom))
            else:
                df_chrom_sorted = df_chrom.sort_values(
                    "start").reset_index(drop=True)
                bin_starts = df_chrom_sorted["start"].values
                bin_ends = df_chrom_sorted["end"].values
                gene_counts = load_gene_counts_per_bin(
                    args.orient_with_gtf, chrom, bin_starts, bin_ends,
                    biotype=args.orient_biotype)
                flip, rho = orient_pc1_by_gene_density(
                    df_chrom_sorted, gene_counts)
                print("  GTF orientation: Spearman(gene density, PC1) = "
                      "{:.3f}; {}".format(rho,
                                            "FLIP" if flip else "keep sign"))
        elif args.normalize_sign:
            df_chrom = read_bigwig_chrom(path, chrom)
            if len(df_chrom) > 0:
                mu = float(df_chrom["value"].mean())
                pc1_chrom_mean[cond] = mu
                if mu < 0:
                    flip = True
                    print("  flipping sign (chrom mean = {:.4f})".format(mu))

        if flip:
            df_region = df_region.copy()
            df_region["value"] = -df_region["value"]

        pc1_region[cond] = df_region

    if not pc1_region:
        sys.exit("no PC1 data loaded")

    cond_list = list(pc1_region.keys())

    # ---- 1. PC1 at each feature per condition ----
    if args.features_bed:
        print("\n=== PC1 at each feature ===")
        features = load_features_bed(args.features_bed, chrom_filter=chrom)
        rows = []
        for _, feat in features.iterrows():
            mid = (int(feat["start"]) + int(feat["end"])) // 2
            row = {"feature": feat["name"], "chrom": chrom,
                   "feat_start": int(feat["start"]),
                   "feat_end": int(feat["end"]),
                   "feat_mid": mid}
            vals = []
            for cond in cond_list:
                v = value_at_position(pc1_region[cond], mid)
                row[cond + "_PC1"] = v
                row[cond + "_sign"] = ("+" if v > args.switch_tolerance
                                        else ("-" if v < -args.switch_tolerance
                                              else "."))
                vals.append(v)
            rows.append(row)
            # Print compact summary line
            signs = " ".join("{}={:>2}".format(c, row[c+"_sign"]) for c in cond_list)
            print("  {:<14} {}".format(feat["name"], signs))
        feat_df = pd.DataFrame(rows)
        feat_path = args.out_prefix + ".pc1_per_feature.tsv"
        feat_df.to_csv(feat_path, sep="\t", index=False, float_format="%.4f")
        print("Wrote", feat_path)

    # ---- 2. Pairwise switch counts in the region ----
    # Align all conditions onto a common bin grid by taking the union of
    # interval starts (assumes regular 100 kb binning).
    print("\n=== Compartment switches across conditions ({}:{:,}-{:,}) ===".format(
        chrom, r_start, r_end))
    # Build common bin grid from the first condition with data
    ref_cond = cond_list[0]
    ref_df = pc1_region[ref_cond]
    if len(ref_df) == 0:
        sys.exit("reference condition {} has no PC1 data in region".format(ref_cond))
    grid = ref_df[["start", "end"]].drop_duplicates().sort_values(
        "start").reset_index(drop=True)
    grid["mid"] = (grid["start"] + grid["end"]) // 2

    # PC1 per (cond, bin)
    pc1_matrix = pd.DataFrame({"mid": grid["mid"],
                                  "start": grid["start"],
                                  "end": grid["end"]})
    for cond in cond_list:
        df = pc1_region[cond]
        # Map each grid mid to nearest interval; fast via merge_asof
        df_sorted = df.sort_values("start").reset_index(drop=True)
        vals = []
        for m in grid["mid"]:
            hit = df_sorted[(df_sorted["start"] <= m) & (df_sorted["end"] > m)]
            vals.append(float(hit["value"].iloc[0]) if len(hit) else float("nan"))
        pc1_matrix[cond + "_PC1"] = vals
        pc1_matrix[cond + "_sign"] = pc1_matrix[cond + "_PC1"].apply(
            lambda v: ("+" if v > args.switch_tolerance
                         else ("-" if v < -args.switch_tolerance else ".")))

    # Pairwise switch counts
    switch_counts = {}
    print("\nPairwise switch counts (sign change at same bin):")
    print("  cond_A vs cond_B : n_switches / n_bins_compared")
    for i in range(len(cond_list)):
        for j in range(i+1, len(cond_list)):
            ca, cb = cond_list[i], cond_list[j]
            sa = pc1_matrix[ca + "_sign"]
            sb = pc1_matrix[cb + "_sign"]
            valid = (sa != ".") & (sb != ".")
            n_valid = int(valid.sum())
            n_switch = int(((sa != sb) & valid).sum())
            switch_counts[(ca, cb)] = (n_switch, n_valid)
            print("  {:<10} vs {:<10} : {:>3} / {:>3}".format(
                ca, cb, n_switch, n_valid))

    # Switch matrix TSV
    rows = []
    for (ca, cb), (n_sw, n_v) in switch_counts.items():
        rows.append({"cond_A": ca, "cond_B": cb,
                      "n_switches": n_sw, "n_bins_valid": n_v,
                      "switch_fraction": (n_sw / n_v) if n_v else float("nan")})
    sw_df = pd.DataFrame(rows)
    sw_path = args.out_prefix + ".switch_matrix.tsv"
    sw_df.to_csv(sw_path, sep="\t", index=False, float_format="%.4f")
    print("Wrote", sw_path)

    # Also write the full PC1 matrix for inspection
    matrix_path = args.out_prefix + ".pc1_matrix.tsv"
    pc1_matrix.to_csv(matrix_path, sep="\t", index=False, float_format="%.4f")
    print("Wrote PC1 per-bin matrix to", matrix_path)

    # ---- 3. Combined figure ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from matplotlib.ticker import FuncFormatter

    # Optionally load TADs
    tads_by_cond = {}
    if args.tad_dir:
        for cond in cond_list:
            candidates = [
                "merged_corrected_KR_{}_bs_{}_domains.bed".format(
                    cond, args.tad_resolution),
            ]
            for c in candidates:
                full = os.path.join(args.tad_dir, c)
                if os.path.isfile(full):
                    tads_by_cond[cond] = load_tads_bed(full, chrom_filter=chrom)
                    break

    n_cond = len(cond_list)
    has_tads = bool(tads_by_cond)
    has_feat = bool(args.features_bed)

    panel_h_pc1 = 0.6
    panel_h_tad = 0.25 if has_tads else 0
    rows_per_cond = panel_h_pc1 + panel_h_tad
    fig_h = 1.0 + n_cond * rows_per_cond + (0.5 if has_feat else 0)

    fig, axes = plt.subplots(
        n_cond * (2 if has_tads else 1) + (1 if has_feat else 0),
        1, sharex=True,
        figsize=(13, fig_h),
        gridspec_kw={"hspace": 0.05,
                       "height_ratios": (
                           ([panel_h_pc1, panel_h_tad] * n_cond
                              if has_tads else [panel_h_pc1] * n_cond)
                           + ([0.4] if has_feat else []))})
    axes = list(axes) if hasattr(axes, "__iter__") else [axes]

    # Decide a consistent y-range across PC1 panels
    all_vals = np.concatenate([
        pc1_matrix[c + "_PC1"].dropna().values for c in cond_list
    ]) if len(cond_list) else np.array([])
    if len(all_vals):
        ymin = float(np.nanpercentile(all_vals, 1))
        ymax = float(np.nanpercentile(all_vals, 99))
        ybound = max(abs(ymin), abs(ymax))
    else:
        ybound = 1.0

    ax_idx = 0
    for cond in cond_list:
        ax = axes[ax_idx]; ax_idx += 1
        df = pc1_region[cond]
        if len(df) > 0:
            mids = (df["start"] + df["end"]) / 2
            vals = df["value"].values
            pos = vals > 0
            neg = vals < 0
            ax.fill_between(mids, 0, vals, where=pos,
                              color="#d62828", alpha=0.7, linewidth=0,
                              step="mid")
            ax.fill_between(mids, 0, vals, where=neg,
                              color="#003049", alpha=0.7, linewidth=0,
                              step="mid")
        ax.axhline(0, color="black", lw=0.5, alpha=0.5)
        ax.set_ylim(-ybound * 1.1, ybound * 1.1)
        ax.set_xlim(r_start, r_end)
        ax.set_yticks([0])
        ax.set_ylabel(cond + "\nPC1", fontsize=8, rotation=0,
                        ha="right", va="center")
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)

        if has_tads:
            ax_t = axes[ax_idx]; ax_idx += 1
            ax_t.set_xlim(r_start, r_end)
            ax_t.set_yticks([])
            ax_t.set_ylim(0, 1)
            for sp in ("top", "right", "left"):
                ax_t.spines[sp].set_visible(False)
            if cond in tads_by_cond:
                tdf = tads_by_cond[cond]
                tdf = tdf[(tdf["end"] > r_start) & (tdf["start"] < r_end)]
                for i, (_, tad) in enumerate(tdf.iterrows()):
                    s = int(tad["start"]); e = int(tad["end"])
                    color = "#88aaff" if i % 2 == 0 else "#a8c4ff"
                    ax_t.add_patch(Rectangle(
                        (s, 0.2), e - s, 0.6,
                        facecolor=color, edgecolor="#1a3a8f",
                        linewidth=0.4))

    # Features panel
    if has_feat:
        ax_f = axes[ax_idx]
        ax_f.set_xlim(r_start, r_end)
        ax_f.set_yticks([])
        ax_f.set_ylim(0, 1)
        for sp in ("top", "right", "left"):
            ax_f.spines[sp].set_visible(False)
        features = load_features_bed(args.features_bed, chrom_filter=chrom)
        features = features[(features["end"] > r_start)
                              & (features["start"] < r_end)]
        for _, feat in features.iterrows():
            s = int(feat["start"]); e = int(feat["end"])
            mid = (s + e) / 2
            w = max(e - s, (r_end - r_start) * 0.002)
            ax_f.add_patch(Rectangle(
                (mid - w/2, 0.55), w, 0.3,
                facecolor="#ffb703", edgecolor="black", linewidth=0.5))
            ax_f.text(mid, 0.4, feat["name"],
                        ha="center", va="top", fontsize=7,
                        style="italic", rotation=15)
        ax_f.set_ylabel("features", fontsize=8, rotation=0,
                          ha="right", va="center")

    axes[-1].xaxis.set_major_formatter(
        FuncFormatter(lambda x, _: "{:.1f}".format(x / 1e6)))
    axes[-1].set_xlabel(chrom + " position (Mb)")

    res_lbl = args.tad_resolution // 1000
    title = "PC1 compartments across conditions ({}:{:,}-{:,})".format(
        chrom, r_start, r_end)
    if has_tads:
        title += " - TADs at {} kb".format(res_lbl)
    if args.orient_with_gtf:
        title += " - oriented (A=positive)"
    elif args.normalize_sign:
        title += " - sign normalized"
    axes[0].set_title(title, fontsize=10)

    out_path = args.out_prefix + ".compartment_landscape.pdf"
    plt.savefig(out_path, bbox_inches="tight")
    plt.close()
    print("\nWrote", out_path)


if __name__ == "__main__":
    main()
