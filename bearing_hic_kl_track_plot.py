#!/usr/bin/env python3
"""
bearing_hic_kl_track_plot_v10.py

v10 changes from v9:
  - Gene panel now appears BELOW the Hi-C panels (was: above), so
    gene coordinates sit close to the per-track KL panels they need
    to be compared with. The region annotations BED panel moves with
    it. Hi-C panels remain at the top of the figure.
  - New --gtf-max-rows flag (default 3) caps the gene stack panel.
    Genes that would push past the cap are dropped; a small italic
    "+N more not shown" note is added to the panel and the count
    also prints to stdout. Set --gtf-max-rows 0 to disable the cap.
    Useful for gene-dense regions like Cd4/Cd8a flanks or wide views.

v9 changes from v8:
  - Legend / stdout percent labels use {:g} format so fractional
    thresholds (e.g. 0.01%) render correctly.

v7 changes from v6:
  - Legend labels reflect --threshold-scope ('genome' or chrom name).
  - Inset colorbars on each Hi-C panel (A, B, diff) with min/max ticks.

v6 changes from v5:
  - Percentile thresholds default to genome-wide autosomes via
    --threshold-scope genome (was chrom).

Extends bearing_kl_track_plot.py by adding Hi-C panels above the KL stack.

Usage:
  python bearing_hic_kl_track_plot.py \\
    --diff       results_v6/diff_DN_vs_DP.stats.tsv \\
    --categories DN_rep1_cats.json \\
    --region     chr6:40790000-41690000 \\
    --cool-A     ../hic_files/merged_corrected_KR_DN_bs_10000.cool \\
    --cool-B     ../hic_files/merged_corrected_KR_DP_bs_10000.cool \\
    --insul-A    .../merged_corrected_KR_DN_bs_10000_tad_score.bm \\
    --insul-B    .../merged_corrected_KR_DP_bs_10000_tad_score.bm \\
    --regions-bed tcrb_regions_v5.bed \\
    --label-a    DN  --label-b DP \\
    --max-distance 500000 \\
    --out        tcrb_DN_vs_DP_browser.pdf

--show-hic controls which Hi-C panels appear (any subset of A,B,diff;
default A,B,diff). Cool files are KR-balanced via cooler's weight column
unless --no-balance is given.
"""

import argparse
import json
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize, TwoSlopeNorm

try:
    import cooler
    HAVE_COOLER = True
except ImportError:
    HAVE_COOLER = False


# ---------------------------------------------------------------------
# Helpers (mirrored from bearing_kl_track_plot.py)
# ---------------------------------------------------------------------

def parse_categories(cats_json_path):
    with open(cats_json_path) as f:
        d = json.load(f)
    cats = d["categories"]
    return [(cats[k][0], cats[k][1]) for k in sorted(cats.keys(), key=int)]


def parse_region_str(s):
    m = re.match(r"^([^:]+):([0-9,]+)-([0-9,]+)$", s)
    if not m:
        raise argparse.ArgumentTypeError("Cannot parse region: " + s)
    return m.group(1), int(m.group(2).replace(",", "")), \
        int(m.group(3).replace(",", ""))


def load_regions_bed(bed_path):
    regions = []
    with open(bed_path) as f:
        for i, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            name = parts[3] if len(parts) >= 4 else "r{}".format(i)
            regions.append((parts[0], int(parts[1]), int(parts[2]), name))
    return regions


# ---------------------------------------------------------------------
# Hi-C loading + rotated triangle plot
# ---------------------------------------------------------------------

def fetch_cool_region(path, chrom, start, end, balance=True):
    if not HAVE_COOLER:
        sys.exit("cooler not installed: pip install cooler")
    c = cooler.Cooler(path)
    binsize = c.binsize
    region = "{}:{}-{}".format(chrom, start, end)
    M = c.matrix(balance=balance, sparse=True).fetch(region).toarray()
    bins = c.bins().fetch(region).reset_index(drop=True)
    return M, bins["start"].to_numpy(), binsize


def plot_hic_triangle(ax, mat, bin_starts, binsize, region_start, region_end,
                      max_distance, cmap, norm, label=None):
    """Plot mat as 45-degree rotated triangle (genomic axis horizontal,
    contact distance vertical, upper triangle only)."""
    n = mat.shape[0]
    if n == 0:
        return None
    max_d_bins = max(1, int(max_distance / binsize))
    hw = binsize / 2.0
    polys = []
    values = []
    for i in range(n):
        bi = bin_starts[i]
        jmax = min(n, i + max_d_bins + 1)
        for j in range(i, jmax):
            v = mat[i, j]
            if not np.isfinite(v) or v == 0:
                continue
            bj = bin_starts[j]
            cx = (bi + bj) / 2.0 + hw
            cy = (bj - bi) / 2.0
            polys.append(((cx - hw, cy), (cx, cy + hw),
                          (cx + hw, cy), (cx, cy - hw)))
            values.append(v)
    if not polys:
        return None
    pc = PolyCollection(polys, array=np.asarray(values),
                         cmap=cmap, norm=norm,
                         edgecolor="none", linewidth=0, antialiased=False)
    ax.add_collection(pc)
    ax.set_xlim(region_start, region_end)
    # data y extent: 0 to (max_distance / 2) + binsize/2 (top of furthest diamond)
    ax.set_ylim(0, max_distance / 2.0 + binsize)
    ax.set_yticks([])
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    if label:
        ax.set_ylabel(label, fontsize=8, rotation=0,
                       ha="right", va="center")
    return pc


# ---------------------------------------------------------------------
# Insulation
# ---------------------------------------------------------------------

def load_insul_bm(path, chrom, start, end):
    df = pd.read_csv(path, sep="\t", header=None, comment="#")
    if df.shape[1] < 4:
        sys.exit(path + ": .bm needs >=4 columns")
    out = df.iloc[:, :3].copy()
    out.columns = ["chrom", "start", "end"]
    out["score"] = df.iloc[:, 3:].mean(axis=1).values  # row-mean of depths
    out["center"] = (out["start"] + out["end"]) / 2
    out = out[(out["chrom"] == chrom) & (out["end"] > start)
              & (out["start"] < end)]
    return out.reset_index(drop=True)


def load_insul_bm_chr(path, chrom):
    """Load full chromosome of .bm data (no region subsetting) for use as
    the null distribution context."""
    df = pd.read_csv(path, sep="\t", header=None, comment="#")
    if df.shape[1] < 4:
        sys.exit(path + ": .bm needs >=4 columns")
    out = df.iloc[:, :3].copy()
    out.columns = ["chrom", "start", "end"]
    out["score"] = df.iloc[:, 3:].mean(axis=1).values
    out["center"] = (out["start"] + out["end"]) / 2
    out = out[out["chrom"] == chrom].reset_index(drop=True)
    return out


def load_insul_bm_genome(path, autosomes_only=True):
    """Load all chromosomes of .bm data for use as the genome-wide null
    distribution context. By default restricts to autosomes (chr1-chr19
    in mouse) to exclude sex-chromosome edge effects."""
    df = pd.read_csv(path, sep="\t", header=None, comment="#")
    if df.shape[1] < 4:
        sys.exit(path + ": .bm needs >=4 columns")
    out = df.iloc[:, :3].copy()
    out.columns = ["chrom", "start", "end"]
    out["score"] = df.iloc[:, 3:].mean(axis=1).values
    out["center"] = (out["start"] + out["end"]) / 2
    if autosomes_only:
        # mouse mm10 autosomes; also handles human chr1-22 trivially since
        # those aren't in mouse bm files but the filter is harmless
        keep = {"chr{}".format(i) for i in range(1, 23)}
        out = out[out["chrom"].isin(keep)]
    return out.reset_index(drop=True)


def compute_insul_pctile_threshold(ins_A, ins_B, pctile=95.0):
    """Per-bin |delta insul| significance via percentile threshold.

    Operates on whatever subset of bins is passed in - chrom-restricted
    or genome-wide. A bin is flagged 'sig' if its |insul_A - insul_B| is
    at or above the given percentile of the same statistic computed
    across the input bins.

    Returns DataFrame with columns: chrom, center, score_A, score_B, obs,
    threshold, pctile, sig.
    """
    a = ins_A[["chrom", "center", "score"]].rename(columns={"score": "score_A"})
    b = ins_B[["chrom", "center", "score"]].rename(columns={"score": "score_B"})
    merged = pd.merge(a, b, on=["chrom", "center"]).sort_values(
        ["chrom", "center"]).reset_index(drop=True)
    obs = (merged["score_A"] - merged["score_B"]).abs().values
    finite = obs[np.isfinite(obs)]
    if len(finite) == 0:
        threshold = np.nan
    else:
        threshold = float(np.nanpercentile(finite, pctile))
    merged["obs"] = obs
    merged["threshold"] = threshold
    merged["pctile"] = pctile
    merged["sig"] = np.isfinite(obs) & (obs >= threshold)
    return merged


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def _parse_gtf_attrs(attr_str):
    """Parse GTF attribute string 'key \"value\"; key \"value\"; ...' to dict.
    Also supports GFF3 'key=value;key=value;' format as a fallback.
    """
    d = {}
    if '"' in attr_str:
        # GTF
        import re
        for m in re.finditer(r'(\S+)\s+"([^"]*)"', attr_str):
            d[m.group(1)] = m.group(2)
    elif "=" in attr_str:
        # GFF3
        for pair in attr_str.strip().rstrip(";").split(";"):
            if "=" in pair:
                k, v = pair.strip().split("=", 1)
                d[k.strip()] = v.strip()
    return d


def parse_gtf_region(gtf_path, chrom, start, end,
                      biotype_filter="protein_coding"):
    """Extract gene/transcript/exon records from a GTF for a region.

    Uses awk for fast filtering on large genome GTFs; falls back to
    Python parsing if awk is unavailable. Handles 'chr' prefix mismatch
    (GTF can use either 'chr6' or '6').

    biotype_filter: keep only records with matching gene_type or
        gene_biotype attribute. Set to None or 'all' to disable.

    Returns DataFrame with columns: chrom, start (0-based), end, strand,
    feature (gene/transcript/exon), gene_id, gene_name, transcript_id,
    gene_type. Empty DataFrame if no matches.
    """
    import subprocess

    # GTF may use 'chr6' or '6'; try both
    chrom_variants = [chrom]
    if chrom.startswith("chr"):
        chrom_variants.append(chrom[3:])
    else:
        chrom_variants.append("chr" + chrom)

    lines = []
    for chrom_try in chrom_variants:
        awk_prog = (
            '$1 == "' + chrom_try + '" && $4 <= ' + str(end)
            + ' && $5 >= ' + str(start) + ' { print }'
        )
        try:
            proc = subprocess.run(
                ["awk", "-F", "\t", awk_prog, gtf_path],
                capture_output=True, text=True, check=False,
            )
            out = proc.stdout.strip()
            if out:
                lines = out.split("\n")
                break
        except FileNotFoundError:
            # awk not available; fall back to Python parse
            with open(gtf_path) as fh:
                for ln in fh:
                    if ln.startswith("#"):
                        continue
                    parts = ln.rstrip("\n").split("\t")
                    if len(parts) < 9:
                        continue
                    if parts[0] != chrom_try:
                        continue
                    try:
                        s = int(parts[3]); e = int(parts[4])
                    except ValueError:
                        continue
                    if s <= end and e >= start:
                        lines.append(ln.rstrip("\n"))
            if lines:
                break

    if not lines:
        return pd.DataFrame()

    keep_biotype = biotype_filter and biotype_filter != "all"
    records = []
    for ln in lines:
        parts = ln.split("\t")
        if len(parts) < 9:
            continue
        feature = parts[2]
        if feature not in ("gene", "transcript", "exon"):
            continue
        try:
            s = int(parts[3]) - 1  # GTF is 1-based inclusive -> 0-based half-open
            e = int(parts[4])
        except ValueError:
            continue
        attrs = _parse_gtf_attrs(parts[8])
        gt = attrs.get("gene_type", attrs.get("gene_biotype", ""))
        if keep_biotype and gt and gt != biotype_filter:
            continue
        records.append({
            "chrom": parts[0],
            "start": s,
            "end": e,
            "strand": parts[6],
            "feature": feature,
            "gene_id": attrs.get("gene_id", ""),
            "gene_name": attrs.get("gene_name",
                                     attrs.get("Name",
                                                attrs.get("gene_id", ""))),
            "transcript_id": attrs.get("transcript_id", ""),
            "gene_type": gt,
        })
    return pd.DataFrame(records)


def select_canonical_transcripts(gtf_df):
    """For each gene, pick the longest transcript by total exon coverage.

    Returns (transcripts_df, exons_df) where transcripts_df has one row
    per gene (the canonical) and exons_df has the exon rows belonging
    to those transcripts.
    """
    if len(gtf_df) == 0:
        return pd.DataFrame(), pd.DataFrame()

    tx_df = gtf_df[gtf_df["feature"] == "transcript"].copy()
    exon_df = gtf_df[gtf_df["feature"] == "exon"].copy()

    if len(tx_df) == 0:
        # Only gene records present; treat each as a single block
        gene_df = gtf_df[gtf_df["feature"] == "gene"].copy()
        return gene_df, pd.DataFrame()

    # Sum exon length per transcript
    if len(exon_df) > 0:
        exon_df["exon_len"] = exon_df["end"] - exon_df["start"]
        tx_lens = exon_df.groupby("transcript_id")["exon_len"].sum()
        tx_df["tx_len"] = tx_df["transcript_id"].map(tx_lens).fillna(
            tx_df["end"] - tx_df["start"])
    else:
        tx_df["tx_len"] = tx_df["end"] - tx_df["start"]

    # Pick longest transcript per gene
    canonical = (tx_df.sort_values("tx_len", ascending=False)
                  .drop_duplicates(subset=["gene_id"], keep="first")
                  .copy())
    canon_tx_ids = set(canonical["transcript_id"])
    canon_exons = exon_df[exon_df["transcript_id"].isin(canon_tx_ids)].copy()
    return canonical.reset_index(drop=True), canon_exons.reset_index(drop=True)


def stack_intervals(starts, ends, gap=0, max_rows=0):
    """Greedy interval packing. Returns (row_per_item, n_rows, dropped).

    Items placed in lowest row where they don't overlap (with gap).
    If max_rows > 0 and a new row would be needed beyond max_rows, the
    item is dropped (row assignment = -1).
    """
    n = len(starts)
    if n == 0:
        return [], 0, []
    order = sorted(range(n), key=lambda i: starts[i])
    row_ends = []
    assignment = [0] * n
    dropped = []
    for idx in order:
        s = starts[idx]
        e = ends[idx]
        placed = False
        for r, row_end in enumerate(row_ends):
            if s >= row_end + gap:
                row_ends[r] = e
                assignment[idx] = r
                placed = True
                break
        if not placed:
            if max_rows > 0 and len(row_ends) >= max_rows:
                # Cap exceeded: drop this item
                assignment[idx] = -1
                dropped.append(idx)
            else:
                assignment[idx] = len(row_ends)
                row_ends.append(e)
    return assignment, len(row_ends), dropped


def render_genes_panel(ax, canonical_df, exons_df,
                        region_start, region_end,
                        max_labels=40,
                        max_rows=3,
                        plus_color="#1f4e79",
                        minus_color="#a4262c"):
    """Render canonical transcripts as stacked lines with exon boxes.

    Genes packed into minimum rows up to max_rows; overflow dropped.
    Each gene drawn as a thin horizontal line spanning the transcript,
    with thicker filled boxes for exons. Strand encoded by color.
    Labels suppressed if too many genes for legibility (controlled by
    max_labels). The number of dropped genes (if any) is annotated in
    the bottom-right of the panel.
    """
    ax.set_xlim(region_start, region_end)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.set_yticks([])
    ax.set_ylabel("Genes", fontsize=8, rotation=0,
                   ha="right", va="center")

    if canonical_df is None or len(canonical_df) == 0:
        ax.text((region_start + region_end) / 2, 0.5,
                  "no genes in region",
                  ha="center", va="center", fontsize=6, color="gray")
        ax.set_ylim(0, 1)
        return

    starts = canonical_df["start"].to_list()
    ends = canonical_df["end"].to_list()
    region_w = region_end - region_start
    label_gap = region_w * 0.01  # small horizontal gap to keep labels readable
    row_assign, n_rows, dropped_idx = stack_intervals(
        starts, ends, gap=label_gap, max_rows=max_rows)

    # Don't count dropped genes when deciding whether to show labels —
    # what matters is what's actually being rendered.
    n_rendered = sum(1 for r in row_assign if r >= 0)
    show_labels = n_rendered <= max_labels
    row_height = 1.0
    exon_h_frac = 0.55
    line_h_frac = 0.06
    label_pad_frac = 0.08

    if exons_df is not None and len(exons_df) > 0:
        exons_by_tx = exons_df.groupby("transcript_id")
    else:
        exons_by_tx = None

    # Track the last label center x within each row so labels that would
    # overprint horizontally get nudged to a higher tier instead.
    last_label_x_by_row = {}
    min_label_sep = region_w * 0.05  # labels closer than this collide

    for i in range(len(canonical_df)):
        row = row_assign[i]
        if row < 0:
            continue  # dropped due to max_rows cap
        rec = canonical_df.iloc[i]
        s = rec["start"]; e = rec["end"]
        strand = rec["strand"]
        name = rec.get("gene_name", rec.get("gene_id", ""))
        tx_id = rec.get("transcript_id", "")
        y_center = (n_rows - 1 - row) * row_height + row_height / 2
        color = plus_color if strand == "+" else minus_color

        ax.add_patch(plt.Rectangle(
            (s, y_center - line_h_frac / 2),
            e - s, line_h_frac,
            facecolor=color, edgecolor="none", alpha=0.9, zorder=2))

        drew_exons = False
        if exons_by_tx is not None and tx_id and tx_id in exons_by_tx.groups:
            for _, exon in exons_by_tx.get_group(tx_id).iterrows():
                ax.add_patch(plt.Rectangle(
                    (exon["start"], y_center - exon_h_frac / 2),
                    exon["end"] - exon["start"], exon_h_frac,
                    facecolor=color, edgecolor="none", alpha=0.95, zorder=3))
                drew_exons = True
        if not drew_exons:
            ax.add_patch(plt.Rectangle(
                (s, y_center - exon_h_frac / 2),
                e - s, exon_h_frac,
                facecolor=color, edgecolor="none", alpha=0.95, zorder=3))

        mid = max(region_start, min(region_end, (s + e) / 2))
        marker = ">" if strand == "+" else "<"
        ax.plot([mid], [y_center], marker=marker, color="white",
                  markersize=4, markeredgecolor=color, markeredgewidth=0.4,
                  zorder=4)

        if show_labels:
            lx = max(region_start, min(region_end, (s + e) / 2))
            # Stagger label height if it would collide with the previous
            # label in the same row.
            prev_x = last_label_x_by_row.get(row)
            tier = 0
            if prev_x is not None and abs(lx - prev_x) < min_label_sep:
                tier = last_label_x_by_row.get((row, "tier"), 0) + 1
                tier = tier % 3
            last_label_x_by_row[row] = lx
            last_label_x_by_row[(row, "tier")] = tier
            ax.text(
                lx,
                y_center + exon_h_frac / 2 + label_pad_frac + 0.32 * tier,
                name, ha="center", va="bottom", fontsize=6,
                style="italic", clip_on=True, zorder=5,
            )

    if dropped_idx:
        ax.text(
            region_end - region_w * 0.005, 0.05,
            "+{} more not shown".format(len(dropped_idx)),
            ha="right", va="bottom", fontsize=5.5,
            color="gray", style="italic", zorder=6,
        )

    # Extra headroom above the top row so staggered labels (up to ~2 tiers
    # at +0.32 each) are not clipped.
    ax.set_ylim(0, n_rows * row_height + 0.8)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--diff", required=True)
    ap.add_argument("--categories", required=True)
    ap.add_argument("--region", required=True, type=parse_region_str)
    ap.add_argument("--cool-A", required=True)
    ap.add_argument("--cool-B", required=True)
    ap.add_argument("--insul-A", default=None)
    ap.add_argument("--insul-B", default=None)
    ap.add_argument("--regions-bed", default=None)
    ap.add_argument("--gtf", default=None,
                    help="GTF file with gene annotations (optional). If "
                         "provided, adds a gene track panel between the "
                         "regions and Hi-C panels.")
    ap.add_argument("--gtf-biotype", default="protein_coding",
                    help="Filter genes by gene_type / gene_biotype "
                         "(default: protein_coding). Set 'all' to keep all.")
    ap.add_argument("--gtf-max-labels", type=int, default=40,
                    help="Suppress gene name labels if more than N genes "
                         "in region (default: 40).")
    ap.add_argument("--gtf-max-rows", type=int, default=3,
                    help="Cap gene panel at N stack rows. Overflow genes "
                         "are dropped with a note printed to stdout. Set "
                         "to 0 to disable cap (default: 3).")
    ap.add_argument("--out", required=True)
    ap.add_argument("--tracks", nargs="+", default=None,
                    help="subset of tracks to plot (default: all)")
    ap.add_argument("--fdr-only", action="store_true")
    ap.add_argument("--label-a", default="A")
    ap.add_argument("--label-b", default="B")
    ap.add_argument("--track-pctile-low", type=float, default=99.0,
                    help="percentile threshold for light-red "
                         "per-track markers (default 99.0 = top 1%)")
    ap.add_argument("--track-pctile-high", type=float, default=99.5,
                    help="percentile threshold for dark-red "
                         "per-track markers (default 99.5 = top 0.5%)")
    ap.add_argument("--bin-width", type=int, default=200)
    ap.add_argument("--max-distance", type=int, default=500000,
                    help="max Hi-C contact distance to plot (default 500 kb)")
    ap.add_argument("--no-balance", action="store_true",
                    help="skip cooler balance (default: balance applied)")
    ap.add_argument("--hic-vmax", type=float, default=None,
                    help="color scale max for raw Hi-C (default: 99th pctile)")
    ap.add_argument("--diff-hic-vmax", type=float, default=None,
                    help="color scale max for diff Hi-C (default: 99th pctile "
                         "of |A-B|)")
    ap.add_argument("--show-hic", default="A,B,diff",
                    help="any subset of A,B,diff (default: A,B,diff)")
    ap.add_argument("--insul-pctile", type=float, default=95.0,
                    help="percentile threshold for |delta insul| "
                         "markers (default 95.0; set 0 or 100 to disable)")
    ap.add_argument("--pval-cutoff", type=float, default=0.05,
                    help="reference significance line on the differential "
                         "p-value track (default 0.05; the track plots "
                         "-log10(pval) and a dashed line at this cutoff)")
    ap.add_argument("--threshold-scope", choices=["genome", "chrom"],
                    default="genome",
                    help="reference distribution for percentile thresholds "
                         "(BES tracks and |delta insul|). 'genome' uses all "
                         "autosomes (default; better when perturbation is "
                         "cis-restricted, e.g. EbKO); 'chrom' uses only the "
                         "displayed chromosome")
    args = ap.parse_args()

    balance = not args.no_balance
    # Short scope label for legends and panel annotations
    if args.threshold_scope == "genome":
        scope_short = "genome"
        scope_long = "genome (autosomes)"
    else:
        scope_short = None  # filled in later once chrom is known
        scope_long = None
    show_hic = [s.strip() for s in args.show_hic.split(",") if s.strip()]
    for s in show_hic:
        if s not in ("A", "B", "diff"):
            sys.exit("--show-hic values must be in {A, B, diff}: " + s)

    cats = parse_categories(args.categories)
    all_track_names = [c[0] for c in cats]
    all_track_colors = [c[1] for c in cats]

    if args.tracks:
        plot_tracks = [t for t in all_track_names if t in args.tracks]
        if not plot_tracks:
            sys.exit("--tracks {} matched no known tracks: {}".format(
                args.tracks, all_track_names))
    else:
        plot_tracks = all_track_names
    track_indices = [all_track_names.index(t) for t in plot_tracks]
    track_colors = [all_track_colors[i] for i in track_indices]

    chrom, start, end = args.region
    if args.threshold_scope != "genome":
        scope_short = chrom
        scope_long = chrom

    # BEARING per-bin
    print("Loading", args.diff, "...", flush=True)
    df = pd.read_csv(args.diff, sep="\t", low_memory=False)
    kl_cols_all = sorted([c for c in df.columns if c.startswith("kl_")],
                          key=lambda c: int(c.split("_")[1]))
    if len(kl_cols_all) != len(all_track_names):
        sys.exit("Expected {} kl_* cols, got {}".format(
            len(all_track_names), len(kl_cols_all)))
    kl_cols_plot = [kl_cols_all[i] for i in track_indices]
    mask = (df["chrom"] == chrom) & (df["end"] > start) & (df["start"] < end)
    df_r = df[mask].copy().reset_index(drop=True)
    if args.fdr_only:
        df_r = df_r[df_r["significant_fdr0.05"].astype(int) == 1]
    print("  {} BEARING bins in region".format(len(df_r)))

    # Hi-C
    print("Loading cool A ...", flush=True)
    M_A, bin_starts, bs_A = fetch_cool_region(args.cool_A, chrom, start, end,
                                               balance)
    print("Loading cool B ...", flush=True)
    M_B, bin_starts_B, bs_B = fetch_cool_region(args.cool_B, chrom, start, end,
                                                  balance)
    if bs_A != bs_B:
        sys.exit("cool A binsize ({}) != cool B binsize ({})".format(
            bs_A, bs_B))
    if M_A.shape != M_B.shape:
        sys.exit("cool matrix shapes differ: {} vs {}".format(
            M_A.shape, M_B.shape))
    M_diff = M_A - M_B
    binsize = bs_A
    print("  binsize={}, {} bins".format(binsize, M_A.shape[0]))

    # Color norms
    if args.hic_vmax is None:
        pos = M_A[np.isfinite(M_A) & (M_A > 0)]
        hic_vmax = float(np.nanpercentile(pos, 99)) if len(pos) else 1.0
    else:
        hic_vmax = args.hic_vmax
    hic_norm = Normalize(vmin=0, vmax=hic_vmax)

    if args.diff_hic_vmax is None:
        d_abs = np.abs(M_diff[np.isfinite(M_diff)])
        diff_vmax = float(np.nanpercentile(d_abs, 99)) if len(d_abs) else 1.0
    else:
        diff_vmax = args.diff_hic_vmax
    if diff_vmax <= 0:
        diff_vmax = max(hic_vmax * 0.5, 1e-9)
    diff_norm = TwoSlopeNorm(vmin=-diff_vmax, vcenter=0, vmax=diff_vmax)

    # Insulation
    ins_A = ins_B = None
    insul_sig_df = None
    do_insul_sig = (0 < args.insul_pctile < 100)
    if args.insul_A and args.insul_B:
        print("Loading insulation ...", flush=True)
        ins_A = load_insul_bm(args.insul_A, chrom, start, end)
        ins_B = load_insul_bm(args.insul_B, chrom, start, end)

        if do_insul_sig:
            scope_label = ("genome (autosomes)"
                           if args.threshold_scope == "genome" else chrom)
            print("Computing {} |delta insul| {}th percentile ..."
                  .format(scope_label, args.insul_pctile), flush=True)
            if args.threshold_scope == "genome":
                ins_A_ref = load_insul_bm_genome(args.insul_A)
                ins_B_ref = load_insul_bm_genome(args.insul_B)
            else:
                ins_A_ref = load_insul_bm_chr(args.insul_A, chrom)
                ins_B_ref = load_insul_bm_chr(args.insul_B, chrom)
            insul_sig_df = compute_insul_pctile_threshold(
                ins_A_ref, ins_B_ref, pctile=args.insul_pctile)
            # For display, restrict to the displayed chrom
            insul_sig_df_chr = insul_sig_df[
                insul_sig_df["chrom"] == chrom].reset_index(drop=True)
            n_sig_all = int(insul_sig_df["sig"].sum())
            n_sig_chr = int(insul_sig_df_chr["sig"].sum())
            sub = insul_sig_df_chr[
                (insul_sig_df_chr["center"] >= start)
                & (insul_sig_df_chr["center"] <= end)
                & insul_sig_df_chr["sig"]
            ]
            print("  threshold |delta insul| = {:.3f}  ({} sig bins {}-wide; "
                  "{} on {}; {} in displayed region)".format(
                      insul_sig_df["threshold"].iloc[0],
                      n_sig_all, scope_label,
                      n_sig_chr, chrom, len(sub)))
            # Replace insul_sig_df with the chrom-restricted version for
            # downstream rendering (the panel only shows the displayed chrom)
            insul_sig_df = insul_sig_df_chr

    # Annotations
    anno_regions = []
    if args.regions_bed:
        for c2, s2, e2, n2 in load_regions_bed(args.regions_bed):
            if c2 == chrom and e2 > start and s2 < end:
                anno_regions.append((c2, s2, e2, n2))

    # GTF gene annotations
    gtf_canon = gtf_exons = None
    if args.gtf:
        print("Loading GTF {} ...".format(args.gtf), flush=True)
        biotype = (None if args.gtf_biotype in (None, "", "all")
                    else args.gtf_biotype)
        gtf_df = parse_gtf_region(args.gtf, chrom, start, end,
                                    biotype_filter=biotype)
        if len(gtf_df) > 0:
            gtf_canon, gtf_exons = select_canonical_transcripts(gtf_df)
            n_genes = len(gtf_canon)
            biotype_label = biotype if biotype else "any biotype"
            print("  {} genes ({}) in region, {} exons (canonical "
                  "transcripts)".format(n_genes, biotype_label,
                                          len(gtf_exons)))
        else:
            print("  no GTF records in region")

    # Panel layout. Hi-C panel height auto-computed so the rotated triangle
    # renders at correct 1:1 data aspect (data height = max_distance/2,
    # data width = region width; approximate panel width = 11.5 inches).
    region_w = end - start
    hic_data_h = args.max_distance / 2.0
    PANEL_WIDTH_IN = 11.5  # approx panel width inside the 13" figure
    hic_height_in = PANEL_WIDTH_IN * hic_data_h / region_w
    # Convert other panel logical sizes to inches (0.6 in per logical unit)
    PER_UNIT_IN = 0.6
    panels = []
    # Hi-C panels come first (top of figure)
    for k in show_hic:
        panels.append(("hic", hic_height_in, k))
    # BEARING score (BES) track sits directly under the Hi-C maps so the
    # domain-scale contact difference and its per-bin BEARING-score
    # quantification read together as one block.
    panels.append(("bes", 0.9 * PER_UNIT_IN, None))
    # Region annotations and gene track sit RIGHT BELOW Hi-C so coordinates
    # are visually close to the chromatin maps. They previously sat at the
    # very top of the figure (too far from the per-track panels below).
    if anno_regions:
        panels.append(("anno", 0.4 * PER_UNIT_IN, None))
    if gtf_canon is not None and len(gtf_canon) > 0:
        # Gene panel height scales with number of stack rows. Cap at
        # --gtf-max-rows (default 3). Overflow genes are dropped so the
        # panel doesn't dominate the figure when regions have hundreds
        # of genes. The actual stacking happens again in render time
        # (must use same gap value for consistency).
        max_rows = args.gtf_max_rows if args.gtf_max_rows > 0 else 0
        _, est_rows, dropped_idx = stack_intervals(
            gtf_canon["start"].to_list(),
            gtf_canon["end"].to_list(),
            gap=(end - start) * 0.01,
            max_rows=max_rows)
        if dropped_idx:
            print("  NOTE: capped gene panel at {} rows; {} genes dropped "
                  "for legibility (use --gtf-max-rows 0 to disable)".format(
                      max_rows, len(dropped_idx)))
        # Roughly 0.25 in per row, min 0.6 in, max 2.5 in
        gene_panel_h = max(0.6, min(2.5, 0.25 * est_rows + 0.3))
        panels.append(("genes", gene_panel_h, None))
    if ins_A is not None and ins_B is not None:
        panels.append(("insul", 1.0 * PER_UNIT_IN, None))
        panels.append(("delta_insul", 0.8 * PER_UNIT_IN, None))
    # Differential p-value track sits just above the per-track KL panels so
    # the significance evidence reads next to the per-assay attribution.
    if "pval" in df_r.columns:
        panels.append(("pval", 0.9 * PER_UNIT_IN, None))
    for t, col, color in zip(plot_tracks, kl_cols_plot, track_colors):
        panels.append(("track", 1.0 * PER_UNIT_IN, (t, col, color)))

    heights = [h for _, h, _ in panels]
    n_panels = len(panels)
    total_height_in = sum(heights) + 1.5  # title + bottom xlabel + margin
    fig, axes = plt.subplots(
        n_panels, 1,
        figsize=(13, total_height_in),
        sharex=True,
        gridspec_kw={"hspace": 0.12, "height_ratios": heights},
    )
    if n_panels == 1:
        axes = [axes]

    # Per-track percentile thresholds for two-tier markers. Scope set
    # by --threshold-scope: 'genome' uses all autosomes (better when the
    # perturbation is cis-restricted, since other chromosomes serve as a
    # natural negative control distribution); 'chrom' uses only the
    # displayed chromosome.
    track_thresholds = {}
    if args.threshold_scope == "genome":
        autosomes = {"chr{}".format(i) for i in range(1, 23)}
        df_ref = df[df["chrom"].isin(autosomes)]
        scope_label = "genome (autosomes)"
    else:
        df_ref = df[df["chrom"] == chrom]
        scope_label = chrom
    print("Per-track {} |kl| thresholds (top {:g}% / top {:g}%):"
          .format(scope_label,
                  100 - args.track_pctile_low,
                  100 - args.track_pctile_high))
    for col, name in zip(kl_cols_plot, plot_tracks):
        abs_vals = df_ref[col].abs().to_numpy()
        finite = abs_vals[np.isfinite(abs_vals)]
        if len(finite) == 0:
            track_thresholds[col] = (np.inf, np.inf)
            print("  {:<10} (no data)".format(name))
        else:
            thr_low = float(np.nanpercentile(finite, args.track_pctile_low))
            thr_high = float(np.nanpercentile(finite, args.track_pctile_high))
            track_thresholds[col] = (thr_low, thr_high)
            print("  {:<10} >= {:.2f} (light) / >= {:.2f} (dark)"
                  .format(name, thr_low, thr_high))

    centers_kl = (df_r["start"].to_numpy() + df_r["end"].to_numpy()) / 2
    first_track_drawn = False

    # Render
    for ax, (kind, _, payload) in zip(axes, panels):
        if kind == "anno":
            # Sort regions by position so the vertical-stagger alternation
            # follows genomic order, then give narrow/close labels stepped
            # heights to avoid horizontal overprinting when several small
            # regions cluster (e.g. the Trbv/Trbc series).
            anno_sorted = sorted(anno_regions, key=lambda r: (r[1] + r[2]) / 2)
            region_w_anno = end - start
            last_label_x = None
            stagger = 0
            for i, (_, s2, e2, n2) in enumerate(anno_sorted):
                ax.add_patch(plt.Rectangle((s2, 0.2), e2 - s2, 0.6,
                                             facecolor="#d0e8ff",
                                             edgecolor="navy", linewidth=1.0))
                w_frac = (e2 - s2) / region_w_anno
                cx = (s2 + e2) / 2
                # Decide stagger: if this label's center is within ~6% of the
                # previous label's center, bump it to the next height tier.
                if last_label_x is not None and \
                        abs(cx - last_label_x) < region_w_anno * 0.06:
                    stagger = (stagger + 1) % 3
                else:
                    stagger = 0
                last_label_x = cx
                if w_frac < 0.05:
                    ytext = 1.1 + 0.5 * stagger
                    ax.annotate(n2, xy=(cx, 0.8),
                                 xytext=(cx, ytext),
                                 ha="center", va="bottom", fontsize=7,
                                 fontweight="bold",
                                 arrowprops=dict(arrowstyle="-",
                                                  color="navy", linewidth=0.5))
                else:
                    ax.text(cx, 0.5, n2, ha="center", va="center",
                              fontsize=8, fontweight="bold")
            ax.set_ylim(0, 2.6)
            ax.set_yticks([])
            ax.set_ylabel("Regions", fontsize=8, rotation=0,
                            ha="right", va="center")
            for sp in ("top", "right", "left"):
                ax.spines[sp].set_visible(False)

        elif kind == "genes":
            render_genes_panel(
                ax, gtf_canon, gtf_exons, start, end,
                max_labels=args.gtf_max_labels,
                max_rows=args.gtf_max_rows,
            )

        elif kind == "hic":
            if payload == "A":
                pc_hic = plot_hic_triangle(ax, M_A, bin_starts, binsize, start, end,
                                    args.max_distance, "Reds", hic_norm,
                                    label="{}\nHi-C".format(args.label_a))
                cbar_label = "contact"
            elif payload == "B":
                pc_hic = plot_hic_triangle(ax, M_B, bin_starts, binsize, start, end,
                                    args.max_distance, "Reds", hic_norm,
                                    label="{}\nHi-C".format(args.label_b))
                cbar_label = "contact"
            else:  # diff
                pc_hic = plot_hic_triangle(ax, M_diff, bin_starts, binsize,
                                    start, end, args.max_distance,
                                    "RdBu_r", diff_norm,
                                    label="{} - {}\nHi-C diff".format(
                                        args.label_a, args.label_b))
                cbar_label = "{} - {}".format(args.label_a, args.label_b)
            # Small inset colorbar in upper-right of each Hi-C panel
            if pc_hic is not None:
                from mpl_toolkits.axes_grid1.inset_locator import inset_axes
                cax = inset_axes(ax, width="10%", height="5%",
                                   loc="upper right", borderpad=0.6)
                cb = plt.colorbar(pc_hic, cax=cax, orientation="horizontal")
                cb.ax.tick_params(labelsize=5, length=2, pad=1)
                cb.set_label(cbar_label, fontsize=6, labelpad=1)
                # Reduce tick count: just min and max
                vmin, vmax = pc_hic.get_clim()
                cb.set_ticks([vmin, vmax])
                cb.ax.set_xticklabels(["{:.2g}".format(vmin),
                                         "{:.2g}".format(vmax)])

        elif kind == "insul":
            ax.plot(ins_A["center"], ins_A["score"],
                     color="#1f77b4", label=args.label_a, lw=1.2)
            ax.plot(ins_B["center"], ins_B["score"],
                     color="#d62728", label=args.label_b, lw=1.2)
            ax.legend(loc="upper right", fontsize=7)
            ax.axhline(0, color="gray", lw=0.5, ls="--")
            ax.set_ylabel("Insulation", fontsize=8, rotation=0,
                            ha="right", va="center")

        elif kind == "delta_insul":
            merged = pd.merge(ins_A[["center", "score"]],
                              ins_B[["center", "score"]],
                              on="center", suffixes=("_A", "_B"))
            delta = (merged["score_A"] - merged["score_B"]).abs()
            ax.fill_between(merged["center"], 0, delta,
                              color="steelblue", alpha=0.6)
            ax.plot(merged["center"], delta, color="steelblue", lw=0.9)
            ax.set_ylabel("|delta\ninsul|", fontsize=8, rotation=0,
                            ha="right", va="center")

            # Percentile-based significance markers
            if insul_sig_df is not None and len(insul_sig_df) > 0:
                thr = float(insul_sig_df["threshold"].iloc[0])
                ax.axhline(thr, color="gold", lw=0.8, ls="--",
                            alpha=0.8, zorder=5)
                sub = insul_sig_df[
                    (insul_sig_df["center"] >= start)
                    & (insul_sig_df["center"] <= end)
                    & insul_sig_df["sig"]
                ]
                if not sub.empty:
                    y_offset = max(0.05, float(np.nanmax(delta)) * 0.08)
                    ax.scatter(
                        sub["center"], sub["obs"] + y_offset,
                        marker="*", color="gold", s=40,
                        edgecolor="black", linewidth=0.5,
                        zorder=10,
                        label="top {:g}% {}".format(
                            100 - args.insul_pctile, scope_short),
                    )
                    ax.legend(loc="upper right", fontsize=7)

        elif kind == "bes":
            bes = df_r["bearing_score"].to_numpy()
            ax.bar(centers_kl, bes, width=args.bin_width,
                     color="#404040", edgecolor="none", alpha=0.9)
            ax.axhline(0, color="black", lw=0.5)
            ax.set_ylabel("BES\n({} vs {})".format(
                args.label_a, args.label_b),
                fontsize=8, rotation=0, ha="right", va="center")
            if "significant_fdr0.05" in df_r.columns:
                sig = df_r["significant_fdr0.05"].astype(int) == 1
                if sig.any():
                    ax.scatter(centers_kl[sig], bes[sig] + 0.3,
                                marker="*", color="gold", s=40,
                                edgecolor="black", linewidth=0.5,
                                zorder=10, label="FDR<0.05")
                    ax.legend(loc="upper right", fontsize=7)

        elif kind == "pval":
            # Differential significance as -log10(pval) bars, with a dashed
            # reference line at the chosen cutoff and gold stars on bins
            # passing significant_fdr0.05. Reads the per-bin 'pval' column.
            pv = df_r["pval"].to_numpy(dtype=float)
            with np.errstate(divide="ignore"):
                neglog = -np.log10(np.clip(pv, 1e-300, 1.0))
            ax.bar(centers_kl, neglog, width=args.bin_width,
                     color="#404040", edgecolor="none", alpha=0.9)
            cutoff_y = -np.log10(args.pval_cutoff)
            ax.axhline(cutoff_y, color="gold", lw=0.8, ls="--", alpha=0.9,
                        zorder=5, label="p = {:g}".format(args.pval_cutoff))
            ax.set_ylabel("-log10\np-value", fontsize=8, rotation=0,
                            ha="right", va="center")
            if "significant_fdr0.05" in df_r.columns:
                sig = df_r["significant_fdr0.05"].astype(int).to_numpy() == 1
                if sig.any():
                    y_off = max(0.1, float(np.nanmax(neglog)) * 0.06)
                    ax.scatter(centers_kl[sig], neglog[sig] + y_off,
                                marker="*", color="gold", s=36,
                                edgecolor="black", linewidth=0.5,
                                zorder=10, label="FDR < 0.05")
            ax.legend(loc="upper right", fontsize=6.5, framealpha=0.8)

        elif kind == "track":
            t, col, color = payload
            vals = df_r[col].to_numpy()
            nonzero = vals != 0
            pos = nonzero & (vals > 0)
            neg = nonzero & (vals < 0)
            ax.bar(centers_kl[pos], vals[pos], width=args.bin_width,
                     color=color, edgecolor="none", alpha=0.9)
            ax.bar(centers_kl[neg], vals[neg], width=args.bin_width,
                     color=color, edgecolor="black", linewidth=0.3,
                     alpha=0.5)
            ax.axhline(0, color="black", lw=0.5)
            # Two-tier markers: top pctile_low (light red), top pctile_high
            # (dark red, plotted second so overlapping bins read as dark).
            thr_low, thr_high = track_thresholds.get(col, (np.inf, np.inf))
            abs_vals = np.abs(vals)
            light = (abs_vals >= thr_low) & (abs_vals < thr_high)
            dark = abs_vals >= thr_high
            if light.any():
                ax.scatter(centers_kl[light],
                             np.sign(vals[light]) *
                             (abs_vals[light] + 0.2),
                             marker="v", color="#ff9999",
                             edgecolor="#cc0000", linewidth=0.3,
                             s=18, zorder=10,
                             label=(None if first_track_drawn
                                    else "top {:g}% {}".format(
                                        100 - args.track_pctile_low,
                                        scope_short)))
            if dark.any():
                ax.scatter(centers_kl[dark],
                             np.sign(vals[dark]) *
                             (abs_vals[dark] + 0.2),
                             marker="v", color="#990000",
                             edgecolor="black", linewidth=0.4,
                             s=24, zorder=11,
                             label=(None if first_track_drawn
                                    else "top {:g}% {}".format(
                                        100 - args.track_pctile_high,
                                        scope_short)))
            if not first_track_drawn and (light.any() or dark.any()):
                ax.legend(loc="upper right", fontsize=6, ncol=2,
                            handletextpad=0.2, columnspacing=0.5,
                            framealpha=0.8)
                first_track_drawn = True
            ax.set_ylabel("{}\nkl".format(t), fontsize=8, rotation=0,
                            ha="right", va="center",
                            color=color, fontweight="bold")

    axes[-1].set_xlim(start, end)
    axes[-1].set_xlabel("{} position (bp)".format(chrom), fontsize=9)
    axes[-1].ticklabel_format(useOffset=False, style="plain", axis="x")

    fig.suptitle(
        "BEARING + Hi-C  {}:{:,}-{:,}\n{} vs {}".format(
            chrom, start, end, args.label_a, args.label_b),
        fontsize=11, y=0.995)
    plt.tight_layout()
    plt.savefig(args.out, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print("Wrote", args.out)


if __name__ == "__main__":
    main()
