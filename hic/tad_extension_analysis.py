#!/usr/bin/env python3
"""
tad_extension_analysis_v2.py

v2 changes from v1:
  - --region is now OPTIONAL. If omitted, the display region is
    auto-computed as the union of anchor TADs across all conditions
    plus --region-padding bp (default 2 Mb). This guarantees the
    display covers every condition's anchor TAD even if one is much
    larger than expected.
  - New --plot flag produces a PDF figure with TADs as colored bars
    per condition, anchor line, and feature annotations.
  - New --resolution-label for plot/stdout titles.

Tests the hypothesis that the Tcrb wide locus is part of an extended
architectural domain reaching into Cntnap2 (chr6:45.06-47.30 Mb) and
possibly Gpnmb (chr6:~49 Mb), and that V1 perturbations from Allyn 2024
(dV1P, dV1CTCF) collapse this organization.

Compares HiCExplorer TAD calls across conditions:
  - DN, DP, ProB, EbKO, S3T3   (cell-state contrasts)
  - dV1P, dV1CTCF              (Allyn 2024 V1 perturbations)

For each condition, identifies:
  1. The "Tcrb-anchor TAD" containing chr6:41,550,000 (RC / Eb-Trbv31)
  2. Where its 3' boundary terminates
  3. The chain of TADs from chr6:40M to chr6:50M
  4. Shared vs unique boundaries across conditions

Reports:
  - Per-condition Tcrb-anchor TAD coordinates (TSV)
  - Per-locus TAD assignment (does Cntnap2 share a TAD with Tcrb? With Gpnmb?)
  - Boundary alignment table (cross-condition presence/absence)
  - Optional ASCII boundary diagram for the manuscript

Usage:
  python tad_extension_analysis.py \\
      --tad-dir ../tads/ \\
      --resolution 25000 \\
      --conditions DN DP ProB EbKO S3T3 dV1P dV1CTCF \\
      --anchor chr6:41550000 \\
      --region chr6:39000000-50000000 \\
      --out-prefix tcrb_tad_extension
"""
import argparse
import os
import sys
import pandas as pd
import numpy as np


def load_domains_bed(path, chrom_filter=None):
    """Load HiCExplorer *_domains.bed file. Tolerates variable column counts."""
    rows = []
    with open(path) as fh:
        for ln in fh:
            ln = ln.rstrip("\n")
            if not ln or ln.startswith("#") or ln.startswith("track"):
                continue
            parts = ln.split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            try:
                s = int(parts[1])
                e = int(parts[2])
            except ValueError:
                continue
            if chrom_filter and chrom != chrom_filter:
                continue
            name = parts[3] if len(parts) >= 4 else ""
            rows.append({"chrom": chrom, "start": s, "end": e, "name": name})
    df = pd.DataFrame(rows).sort_values(["chrom", "start"]).reset_index(drop=True)
    return df


def find_containing_tad(tads_df, chrom, pos):
    """Return TAD row containing pos, or None."""
    hits = tads_df[(tads_df["chrom"] == chrom)
                    & (tads_df["start"] <= pos)
                    & (tads_df["end"] >= pos)]
    if len(hits) == 0:
        return None
    return hits.iloc[0]  # there should only be one


def tads_in_region(tads_df, chrom, start, end):
    """Return TADs overlapping the region."""
    return tads_df[(tads_df["chrom"] == chrom)
                    & (tads_df["end"] > start)
                    & (tads_df["start"] < end)].reset_index(drop=True)


def get_boundary_set(tads_df, chrom, start, end, tolerance):
    """Return sorted list of unique boundary positions in [start, end]
    (boundaries = start and end of each TAD), merged within tolerance."""
    sub = tads_in_region(tads_df, chrom, start, end)
    pts = sorted(set(sub["start"].to_list() + sub["end"].to_list()))
    pts = [p for p in pts if start <= p <= end]
    if not pts:
        return []
    merged = [pts[0]]
    for p in pts[1:]:
        if p - merged[-1] > tolerance:
            merged.append(p)
    return merged


def main():
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tad-dir", required=True,
                    help="directory with *_bs_<RES>_domains.bed files")
    ap.add_argument("--resolution", type=int, default=25000,
                    help="binsize used for TAD calls (default 25000)")
    ap.add_argument("--conditions", nargs="+", required=True,
                    help="condition names matching filename slot, e.g. "
                         "DN DP ProB EbKO S3T3 dV1P dV1CTCF")
    ap.add_argument("--anchor", default="chr6:41550000",
                    help="anchor position in chrom:pos format (default "
                         "chr6:41550000 = Tcrb RC / Eb-Trbv31 region)")
    ap.add_argument("--region", default=None,
                    help="display region (chrom:start-end). If omitted, "
                         "auto-computed as the union of all anchor TADs "
                         "across conditions, with --region-padding bp "
                         "on each side. Use this to limit the view.")
    ap.add_argument("--region-padding", type=int, default=2000000,
                    help="padding in bp around auto-detected region "
                         "(default 2000000 = 2 Mb)")
    ap.add_argument("--boundary-tolerance", type=int, default=50000,
                    help="merge boundaries within this many bp when "
                         "comparing across conditions (default 50000)")
    ap.add_argument("--features-bed", default=None,
                    help="optional BED of named features to label "
                         "(e.g. Tcrb V1, RC, Cntnap2, Gpnmb)")
    ap.add_argument("--out-prefix", required=True,
                    help="output filename prefix")
    ap.add_argument("--plot", action="store_true",
                    help="produce PDF figure showing TADs across conditions")
    ap.add_argument("--plot-out", default=None,
                    help="explicit path for plot PDF "
                         "(default: <out-prefix>.tad_landscape.pdf)")
    ap.add_argument("--resolution-label", default=None,
                    help="optional label for plot title and stdout "
                         "(default: derived from --resolution)")
    args = ap.parse_args()

    # Parse anchor; region parsed later (may be auto-detected)
    a_chrom, a_pos = args.anchor.split(":")
    a_pos = int(a_pos)
    chrom = a_chrom
    r_start = r_end = None
    if args.region:
        r_chrom, r_range = args.region.split(":")
        if r_chrom != a_chrom:
            sys.exit("anchor and region chromosomes must match")
        r_start, r_end = map(int, r_range.split("-"))

    # Load TAD calls per condition
    tads_by_cond = {}
    missing = []
    for cond in args.conditions:
        # Try several common filename patterns
        candidates = [
            "merged_corrected_KR_{}_bs_{}_domains.bed".format(cond, args.resolution),
            "merged_corrected_KR_{}_bs_{}_domains.bed_{}.bed".format(
                cond, args.resolution, chrom),
        ]
        path = None
        for c in candidates:
            full = os.path.join(args.tad_dir, c)
            if os.path.isfile(full):
                path = full
                break
        if path is None:
            missing.append((cond, candidates))
            continue
        print("Loading", cond, "from", path, flush=True)
        df = load_domains_bed(path, chrom_filter=chrom)
        if len(df) == 0:
            print("  WARNING: no", chrom, "TADs in", path)
        tads_by_cond[cond] = df

    if missing:
        print("\nWARNING: could not find TAD files for some conditions:")
        for cond, cands in missing:
            print("  " + cond + ": tried " + str(cands))

    if not tads_by_cond:
        sys.exit("no TAD data loaded")

    # ---- 1. Anchor TAD per condition ----
    print("\n=== Anchor TAD ({}, position {:,}) ===".format(chrom, a_pos))
    anchor_rows = []
    for cond, df in tads_by_cond.items():
        tad = find_containing_tad(df, chrom, a_pos)
        if tad is None:
            print("  {}: anchor not in any TAD".format(cond))
            anchor_rows.append({"condition": cond, "chrom": chrom,
                                "start": None, "end": None, "width_bp": None})
            continue
        w = int(tad["end"]) - int(tad["start"])
        print("  {:<10} TAD: {}:{:>11,}-{:>11,}  width={:>9,} bp".format(
            cond, chrom, int(tad["start"]), int(tad["end"]), w))
        anchor_rows.append({"condition": cond, "chrom": chrom,
                            "start": int(tad["start"]),
                            "end": int(tad["end"]), "width_bp": w})
    anchor_df = pd.DataFrame(anchor_rows)
    anchor_path = args.out_prefix + ".anchor_tad.tsv"
    anchor_df.to_csv(anchor_path, sep="\t", index=False)
    print("Wrote anchor TADs to", anchor_path)

    # Auto-detect display region if not provided
    if r_start is None or r_end is None:
        valid = anchor_df.dropna(subset=["start", "end"])
        if len(valid) > 0:
            r_start = int(valid["start"].min()) - args.region_padding
            r_end = int(valid["end"].max()) + args.region_padding
        else:
            # Fall back to a wide window centered on anchor
            r_start = a_pos - 5000000
            r_end = a_pos + 5000000
        # Don't go negative
        r_start = max(0, r_start)
        print("\nAuto-detected display region: {}:{:,}-{:,} ({:.1f} Mb wide)"
              .format(chrom, r_start, r_end, (r_end - r_start) / 1e6))

    # ---- 2. Per-feature TAD assignment ----
    if args.features_bed:
        print("\n=== Per-feature TAD assignment ===")
        features = load_domains_bed(args.features_bed, chrom_filter=chrom)
        feat_rows = []
        for _, feat in features.iterrows():
            mid = (int(feat["start"]) + int(feat["end"])) // 2
            row = {"feature": feat["name"], "chrom": chrom,
                   "feat_start": int(feat["start"]),
                   "feat_end": int(feat["end"]),
                   "feat_mid": mid}
            for cond, df in tads_by_cond.items():
                tad = find_containing_tad(df, chrom, mid)
                if tad is None:
                    row[cond + "_tad"] = ""
                    row[cond + "_tad_width"] = 0
                else:
                    row[cond + "_tad"] = "{}:{}-{}".format(
                        chrom, int(tad["start"]), int(tad["end"]))
                    row[cond + "_tad_width"] = int(tad["end"]) - int(tad["start"])
            feat_rows.append(row)
        feat_df = pd.DataFrame(feat_rows)
        feat_path = args.out_prefix + ".per_feature_tad.tsv"
        feat_df.to_csv(feat_path, sep="\t", index=False)
        print("Wrote per-feature TAD assignment to", feat_path)

        # Quick "does Tcrb share a TAD with X" cross-tab
        print("\n  Shared-TAD pairs (within each condition):")
        anchor_feat = None
        for i, feat in features.iterrows():
            if "Tcrb" in feat["name"] or "RC" in feat["name"]:
                anchor_feat = feat["name"]
                break
        if anchor_feat:
            anchor_idx = features[features["name"] == anchor_feat].index[0]
            for cond in tads_by_cond.keys():
                tad_anchor = feat_df.at[anchor_idx, cond + "_tad"]
                shared = []
                for i, row in feat_df.iterrows():
                    if i == anchor_idx:
                        continue
                    if row[cond + "_tad"] == tad_anchor and tad_anchor:
                        shared.append(row["feature"])
                if shared:
                    print("  {:<10} {} shares TAD with: {}".format(
                        cond, anchor_feat, ", ".join(shared)))
                else:
                    print("  {:<10} {} TAD does not contain other features".format(
                        cond, anchor_feat))

    # ---- 3. Boundary alignment across conditions ----
    print("\n=== Boundary alignment ({}:{}-{}) ===".format(
        chrom, r_start, r_end))
    # Collect every boundary from every condition in the region
    all_pts = []
    for cond, df in tads_by_cond.items():
        for p in get_boundary_set(df, chrom, r_start, r_end,
                                    args.boundary_tolerance):
            all_pts.append((p, cond))
    all_pts.sort()

    # Cluster boundaries within tolerance across conditions
    clusters = []  # list of (rep_pos, dict of cond -> nearest_pos)
    used = [False] * len(all_pts)
    for i, (p, cond) in enumerate(all_pts):
        if used[i]:
            continue
        cluster = {cond: p}
        cluster_positions = [p]
        used[i] = True
        for j in range(i + 1, len(all_pts)):
            if used[j]:
                continue
            p2, cond2 = all_pts[j]
            if p2 - p > args.boundary_tolerance * 2:
                break
            if abs(p2 - np.mean(cluster_positions)) <= args.boundary_tolerance:
                if cond2 not in cluster:
                    cluster[cond2] = p2
                    cluster_positions.append(p2)
                    used[j] = True
        rep_pos = int(np.mean(cluster_positions))
        clusters.append((rep_pos, cluster))

    # Build boundary matrix table
    cond_list = list(tads_by_cond.keys())
    rows = []
    for rep_pos, cl in clusters:
        row = {"chrom": chrom, "boundary_rep_pos": rep_pos}
        n_present = 0
        for c in cond_list:
            if c in cl:
                row[c] = cl[c]
                n_present += 1
            else:
                row[c] = ""
        row["n_conditions"] = n_present
        rows.append(row)
    bound_df = pd.DataFrame(rows).sort_values("boundary_rep_pos").reset_index(drop=True)
    bound_path = args.out_prefix + ".boundary_alignment.tsv"
    bound_df.to_csv(bound_path, sep="\t", index=False)
    print("Wrote boundary alignment to", bound_path,
          "(" + str(len(bound_df)) + " unique boundary clusters)")

    # ---- 4. ASCII diagram ----
    print("\n=== TAD landscape ({}:{:,}-{:,}) ===\n".format(
        chrom, r_start, r_end))
    WIDTH = 80
    span = r_end - r_start
    def pos_to_col(p):
        return int((p - r_start) / span * WIDTH)

    # Header tick line
    n_ticks = 6
    tick_positions = [r_start + i * span // (n_ticks - 1) for i in range(n_ticks)]
    tick_line = list(" " * WIDTH)
    label_line = list(" " * WIDTH)
    for tp in tick_positions:
        c = min(WIDTH - 1, max(0, pos_to_col(tp)))
        tick_line[c] = "|"
        lbl = "{:.1f}M".format(tp / 1e6)
        for k, ch in enumerate(lbl):
            if c + k < WIDTH:
                label_line[c + k] = ch
    print("           " + "".join(tick_line))
    print("           " + "".join(label_line))

    # Per-condition row showing TADs as [====]
    for cond in cond_list:
        df = tads_by_cond[cond]
        sub = tads_in_region(df, chrom, r_start, r_end)
        line = list("-" * WIDTH)
        for _, tad in sub.iterrows():
            s = max(0, pos_to_col(int(tad["start"])))
            e = min(WIDTH - 1, pos_to_col(int(tad["end"])))
            if e <= s:
                e = s + 1
            line[s] = "["
            if e < WIDTH:
                line[e] = "]"
            for k in range(s + 1, e):
                if line[k] == "-":
                    line[k] = "="
        # Mark anchor
        ac = pos_to_col(a_pos)
        if 0 <= ac < WIDTH:
            line[ac] = "*"
        print("  {:<9} ".format(cond) + "".join(line))

    # Mark feature positions on a bottom row if features given
    if args.features_bed:
        features = load_domains_bed(args.features_bed, chrom_filter=chrom)
        flabel = list(" " * WIDTH)
        for _, feat in features.iterrows():
            mid = (int(feat["start"]) + int(feat["end"])) // 2
            c = pos_to_col(mid)
            if 0 <= c < WIDTH:
                flabel[c] = "v"
        print("           " + "".join(flabel))
        # feature name row(s) - simple, last-name-wins on collision
        flabels2 = list(" " * WIDTH)
        for _, feat in features.iterrows():
            mid = (int(feat["start"]) + int(feat["end"])) // 2
            c = pos_to_col(mid)
            name = feat["name"][:10]
            for k, ch in enumerate(name):
                if 0 <= c + k < WIDTH:
                    flabels2[c + k] = ch
        print("           " + "".join(flabels2))

    print("\nLegend: [====] = TAD,  * = anchor ({:,}),  v = feature".format(a_pos))

    # ---- 5. Matplotlib figure ----
    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.patches import Rectangle

        res_label = args.resolution_label or "{} kb".format(args.resolution // 1000)
        n_cond = len(cond_list)
        fig_height = 0.7 + 0.4 * n_cond + (0.4 if args.features_bed else 0)
        fig, ax = plt.subplots(figsize=(13, fig_height))

        # Per-condition TAD rendering
        y_step = 1.0
        for row_idx, cond in enumerate(cond_list):
            y = n_cond - 1 - row_idx
            df = tads_by_cond[cond]
            sub = tads_in_region(df, chrom, r_start, r_end)
            for i, (_, tad) in enumerate(sub.iterrows()):
                s = int(tad["start"]); e = int(tad["end"])
                color = "#3a86ff" if i % 2 == 0 else "#83b5ff"
                ax.add_patch(Rectangle(
                    (s, y + 0.15), e - s, 0.7,
                    facecolor=color, edgecolor="#1a3a8f", linewidth=0.6))
            ax.text(r_start - (r_end - r_start) * 0.005, y + 0.5,
                      cond, ha="right", va="center", fontsize=9,
                      fontweight="bold")

        # Anchor line
        ax.axvline(a_pos, color="#e63946", linewidth=1.0,
                     linestyle="--", alpha=0.7, zorder=5,
                     label="anchor ({:,})".format(a_pos))

        # Features
        if args.features_bed:
            features = load_domains_bed(args.features_bed, chrom_filter=chrom)
            features = features[(features["end"] > r_start)
                                  & (features["start"] < r_end)]
            feat_y = -0.6
            for _, feat in features.iterrows():
                s = int(feat["start"]); e = int(feat["end"])
                mid = (s + e) / 2
                w = max(e - s, (r_end - r_start) * 0.002)
                ax.add_patch(Rectangle(
                    (mid - w/2, feat_y - 0.05), w, 0.25,
                    facecolor="#ffb703", edgecolor="black", linewidth=0.5))
                ax.text(mid, feat_y - 0.18, feat["name"],
                          ha="center", va="top", fontsize=7,
                          style="italic", rotation=15)
            ax.set_ylim(feat_y - 0.7, n_cond + 0.2)
        else:
            ax.set_ylim(-0.3, n_cond + 0.2)

        ax.set_xlim(r_start, r_end)
        ax.set_yticks([])
        ax.set_xlabel("{} position (Mb)".format(chrom))
        # Format x-axis as Mb
        from matplotlib.ticker import FuncFormatter
        ax.xaxis.set_major_formatter(
            FuncFormatter(lambda x, _: "{:.1f}".format(x / 1e6)))
        ax.set_title("TAD landscape across conditions ({}, res = {})".format(
            chrom, res_label), fontsize=10)
        ax.legend(loc="upper right", fontsize=7)

        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)

        plot_path = args.plot_out or (args.out_prefix + ".tad_landscape.pdf")
        plt.tight_layout()
        plt.savefig(plot_path, bbox_inches="tight")
        plt.close()
        print("Wrote TAD landscape plot to", plot_path)


if __name__ == "__main__":
    main()
