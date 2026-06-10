#!/usr/bin/env python3
"""
bearing_hic_plot.py
===================
Composite publication-quality figure combining Hi-C contact maps with
Bearing epilogos tracks for two conditions at a genomic locus.

LAYOUT
------

    5 columns x 5 rows.

    All vertical panels in row 0 share the same Y orientation as the
    Hi-C square: genomic position increases top-to-bottom (region start
    at top, region end at bottom). Track labels in row 0 are horizontal.

         col 0           col 1        col 2        col 3       col 4
       +---------------+-----------+------------+-----------+---------+
    r0 | Hi-C square   | Cond B    | Diff (B-A) | Genes     | Axis    |
       | A=lower-left  | epilogos  | epilogos   | (vertical,|(vertical|
       | B=upper-right | (vertical)| (vertical) | horiz.    | coords, |
       | origin=top    | top=start | top=start  | labels)   | horiz.  |
       |               | horiz.lbl | horiz.lbl  |           | labels) |
       +---------------+-----------+------------+-----------+---------+
    r1 | Cond A epilogos (horizontal, left=start, right=end)           |
       +----------------------------------------------------------------+
    r2 | Diff (A-B) (horizontal)                                        |
       +----------------------------------------------------------------+
    r3 | Gene track (horizontal)              | legend                  |
       +--------------------------------------+                         |
    r4 | Genomic axis (horizontal)            | legend                  |
       +--------------------------------------+-------------------------+

    Vertical orientation: all vertical tracks in row 0 run top=region_start
    to bottom=region_end, matching the Hi-C matrix (imshow origin="upper").
    This means the vertical epilogos and gene tracks are directly aligned
    with the Hi-C rows without any mental flipping.

The Hi-C square shows condition A in the lower-left triangle and condition B
in the upper-right triangle, separated by the diagonal. The differential track
(--diff-qcat) appears both horizontally and vertically using the same sign
convention (A-B), making enrichment changes readable in both orientations.
Loop arcs (--loops BEDPE) are shown as crosshair markers on the Hi-C matrix.
Highlight regions (--highlights BED4) appear as consistent colored bands across
all panels. The gene track and axis run both horizontally (rows 3-4) and
vertically (row 0 cols 3-4) for direct alignment with the Hi-C square.

USAGE
-----
    # Single region mode (.cool/.mcool preferred):
    python bearing_hic_plot.py \
        --contact-a condA.mcool --contact-b condB.mcool \
        --qcat-a condA.qcat.bgz --qcat-b condB.qcat.bgz \
        --region chr6:50000000-52000000 \
        --resolution 10000 \
        --out figure.pdf

    # Batch mode (multiple regions from TSV):
    python bearing_hic_plot.py \
        --contact-a condA.mcool --contact-b condB.mcool \
        --qcat-a condA.qcat.bgz --qcat-b condB.qcat.bgz \
        --regions-file regions.tsv \
        --outdir figures/

    # .hic format (Juicer, auto-detected from extension):
    python bearing_hic_plot.py \
        --contact-a condA.hic --contact-b condB.hic \
        --qcat-a condA.qcat.bgz --qcat-b condB.qcat.bgz \
        --region chr6:50000000-52000000 \
        --resolution 10000 \
        --out figure.pdf

    # With optional annotations:
    python bearing_hic_plot.py \
        --contact-a condA.mcool --contact-b condB.mcool \
        --qcat-a condA.qcat.bgz --qcat-b condB.qcat.bgz \
        --region chr6:50000000-52000000 --resolution 10000 \
        --loops loops.bedpe --genes genes.bed \
        --highlights highlights.bed --diff-qcat diff.qcat.bgz \
        --out figure.pdf

    # Use KL + p-value overlay mode (compact twin-axis significance track):
    python bearing_hic_plot.py \
        --contact-a condA.mcool --contact-b condB.mcool \
        --qcat-a condA.qcat.bgz --qcat-b condB.qcat.bgz \
        --pval-a condA.neglog10p.bw --pval-b condB.neglog10p.bw \
        --pval-overlay --pval-cutoff 0.05 \
        --region chr6:50000000-52000000 --resolution 10000 \
        --out figure_overlay.pdf

INPUT FORMATS
-------------
--contact-a/b : Hi-C contact file. .cool/.mcool (cooler, default/preferred)
                or .hic (Juicer). Format auto-detected from extension.
                Override with --format cool|hic.
--qcat-a/b    : Bearing qcat.bgz files (tabix-indexed, from bigwig_to_qcat.py)
--diff-qcat   : Differential qcat.bgz (from compare_qcat.py, optional)
--pval-diff   : Signed -log10(p-value) track from bearing_pvalue.py --diff
                (optional). Accepts either a BigWig (*.neglog10p.bw) or a
                stats TSV (*.stats.tsv). Positive values indicate A>B
                enrichment, negative values indicate B>A enrichment.
--loops       : BEDPE file: chr1 s1 e1 chr2 s2 e2 [score]
--genes       : BED6 file: chr start end name score strand
--gtf         : GTF annotation file; labels are taken from gene_name
                (fallback: gene_id)
--highlights  : BED4 file: chr start end color (hex color in col 4)
--regions-file: TSV file with columns: name, region [resolution, label, out]
                  name       - Region identifier (used in output filenames)
                  region     - Genomic coordinates (chr:start-end)
                  resolution - Optional Hi-C bin size override (default: 10000)
                  label      - Optional display label (default: name)
                  out        - Optional custom output filename

DEPENDENCIES
------------
    pip install cooler numpy matplotlib pysam scipy
    pip install hicstraw   # only needed for .hic files

COMPANION SCRIPTS
-----------------
This plotting module is designed to work alongside the following helpers:

  - generate_perm_nulls.py
      End-to-end permutation null orchestration for bearing_pvalue.py
      (per-sample and optional differential nulls).

  - rebin_qcat.py
      Rebins qcat.bgz (and optional p-value BigWig) to coarser bins for
      cleaner wide-region visualization.

  - bearing_hic_plot_pval_overlay.py
      Provides drop-in combined draw functions that overlay a compact
      jewel-purple -log10(p) track onto KL stacked bars using a twin axis
      and an opacity mask (significant bins full opacity, non-significant
      bins desaturated).

      Horizontal replacement:
          draw_epilogos_horizontal + draw_pval_horizontal
          -> draw_epilogos_with_pval_horizontal

      Vertical replacement:
          draw_epilogos_vertical + draw_pval_vertical
          -> draw_epilogos_with_pval_vertical

      Import example:
          from bearing_hic_plot_pval_overlay import (
              draw_epilogos_with_pval_horizontal,
              draw_epilogos_with_pval_vertical,
          )
"""

import argparse
import csv
import json
import math
import re
import sys
from pathlib import Path

import numpy as np
from bearing.utils import iter_non_comment_lines, parse_ucsc_region
from bearing.tsv import read_tsv_dict_rows
from bearing.plot_legend import draw_legend_panel
from bearing.plot_layout import create_main_figure
from bearing.plot_loaders import (
    cool_resolution_variant as _cool_resolution_variant_impl,
    detect_hic_format as _detect_hic_format_impl,
    resolve_contacts_for_region_resolution as _resolve_contacts_for_region_resolution_impl,
)
from bearing.hic_io import (
    load_contact_matrix as _load_contact_matrix_impl,
    load_cool_matrix as _load_cool_matrix_impl,
    load_hic_matrix as _load_hic_matrix_impl,
)
from bearing.track_primitives import (
    add_horizontal_highlights,
    add_vertical_highlights,
)
from bearing.plot_tracks import (
    add_vertical_tick_bars as _add_vertical_tick_bars_impl,
    draw_diff_horizontal as _draw_diff_horizontal_impl,
    draw_diff_vertical as _draw_diff_vertical_impl,
    draw_epilogos_horizontal as _draw_epilogos_horizontal_impl,
    draw_epilogos_vertical as _draw_epilogos_vertical_impl,
    draw_gene_track as _draw_gene_track_impl,
    draw_gene_track_vertical as _draw_gene_track_vertical_impl,
    draw_genomic_axis as _draw_genomic_axis_impl,
    draw_genomic_axis_vertical as _draw_genomic_axis_vertical_impl,
    draw_loops_horizontal as _draw_loops_horizontal_impl,
    draw_loops_vertical as _draw_loops_vertical_impl,
    draw_pval_diff_horizontal as _draw_pval_diff_horizontal_impl,
    draw_pval_diff_vertical as _draw_pval_diff_vertical_impl,
    draw_pval_horizontal as _draw_pval_horizontal_impl,
    draw_pval_vertical as _draw_pval_vertical_impl,
)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyArrowPatch, Arc, FancyBboxPatch

# ---------------------------------------------------------------------------
# Category definitions (must match bigwig_to_qcat.py)
# ---------------------------------------------------------------------------
ALL_CATEGORIES = [
    ("ATAC",                       "#be92e0"),
    ("RNAseq +",                   "#6495ed"),
    ("RNAseq -",                   "#1a3a8f"),
    ("CTCF",                       "#ff2200"),
    ("Cohesin",                    "#8b0000"),
    ("H3K27ac",                    "#00e676"),
    ("Bivalent/Poised TSS",        "#cd5c5c"),
    ("Flanking Bivalent TSS/Enh",  "#e9967a"),
    ("Bivalent Enhancer",          "#bdb76b"),
    ("Repressed PolyComb",         "#808080"),
    ("Weak Repressed PolyComb",    "#c0c0c0"),
    ("Quiescent/Low",              "#d0d0d0"),
]


def load_categories_yaml(path):
    """
    Load category definitions from a YAML or JSON file.

    Returns (categories_list, negative_strand_states) where:
      categories_list        : list of (name, color) tuples
      negative_strand_states : set of 1-based indices where negative_strand=true
    
    Supports two formats:
    
    1. Standard YAML/JSON (list of dicts):
       categories:
         - name: "ATAC"
           color: "#be92e0"
         - name: "CTCF"
           color: "#ff2200"
    
    2. Numeric-keyed JSON (from bearmon):
       {
         "categories": {
           "1": ["ATAC", "#be92e0"],
           "2": ["CTCF", "#ff2200"]
         }
       }
    """
    try:
        import yaml
    except ImportError:
        sys.exit(
            "ERROR: PyYAML is required for --categories.  "
            "Install with:  pip install pyyaml"
        )
    with open(path) as fh:
        doc = yaml.safe_load(fh)
    if not isinstance(doc, dict) or "categories" not in doc:
        sys.exit(f"ERROR: {path} must contain a top-level 'categories' key.")
    
    entries = doc["categories"]
    categories_list = []
    negative_strand_states = set()
    
    # Handle two different structures for the categories list
    
    if isinstance(entries, list):
        # Format 1: list of dicts
        for i, entry in enumerate(entries, start=1):
            if not isinstance(entry, dict):
                sys.exit(
                    f"ERROR: category {i} in {path} is not a dict (got {type(entry).__name__}: {repr(entry)}).\n"
                    f"       Each category entry must be a dict with 'name' and 'color' keys.\n"
                    f"       Expected format:\n"
                    f"         categories:\n"
                    f"           - name: 'ATAC'\n"
                    f"             color: '#be92e0'\n"
                    f"           - name: 'CTCF'\n"
                    f"             color: '#ff2200'"
                )
            name = entry.get("name")
            color = entry.get("color")
            if not name or not color:
                sys.exit(
                    f"ERROR: category {i} in {path} is missing 'name' or 'color'.\n"
                    f"       Got: {entry}"
                )
            categories_list.append((name, color))
            if entry.get("negative_strand", False):
                negative_strand_states.add(i)
    
    elif isinstance(entries, dict):
        # Format 2: numeric-key dict with [name, color] values (from bearmon JSON)
        # Sort by numeric key to ensure consistent ordering
        try:
            sorted_keys = sorted(entries.keys(), key=lambda x: int(x))
        except ValueError:
            # Keys are not numeric; just sort alphabetically
            sorted_keys = sorted(entries.keys())
        
        for idx, key in enumerate(sorted_keys, start=1):
            entry = entries[key]
            if isinstance(entry, (list, tuple)) and len(entry) >= 2:
                name = entry[0]
                color = entry[1]
                if not isinstance(name, str) or not isinstance(color, str):
                    sys.exit(
                        f"ERROR: category '{key}' in {path} has invalid types.\n"
                        f"       Expected [name_str, color_str], got {entry}"
                    )
                categories_list.append((name, color))
            elif isinstance(entry, dict):
                # Also support dict format within numeric keys
                name = entry.get("name")
                color = entry.get("color")
                if not name or not color:
                    sys.exit(
                        f"ERROR: category '{key}' in {path} is missing 'name' or 'color'.\n"
                        f"       Got: {entry}"
                    )
                categories_list.append((name, color))
                if entry.get("negative_strand", False):
                    negative_strand_states.add(idx)
            else:
                sys.exit(
                    f"ERROR: category '{key}' in {path} has unexpected format.\n"
                    f"       Expected [name, color] or {{'name': ..., 'color': ...}}\n"
                    f"       Got: {entry}"
                )
    else:
        sys.exit(
            f"ERROR: 'categories' in {path} must be a list or dict, "
            f"got {type(entries).__name__}."
        )
    
    if not categories_list:
        sys.exit(f"ERROR: no categories found in {path}")
    
    return categories_list, negative_strand_states


# ---------------------------------------------------------------------------
# Region parsing
# ---------------------------------------------------------------------------

def parse_region(region_str):
    """
    Parse a UCSC region string like chr6:50000000-52000000.
    Returns (chrom, start, end) as (str, int, int).
    """
    return parse_ucsc_region(region_str)


def _non_comment_lines(handle):
    """Yield TSV lines excluding blanks and lines starting with '#' (after trim)."""
    yield from iter_non_comment_lines(handle)


def load_regions_file(path):
    """
    Load genomic regions from a TSV file.
    
    Expected columns:
      - name (required): region identifier used in output filenames
      - region (required): genomic coordinates (chr:start-end)
      - resolution (optional): Hi-C bin resolution override
      - label (optional): custom display label
      - out (optional): custom output filename
    
    Returns list of dicts with keys: name, region, resolution, label, out
    """
    regions = []
    try:
        rows = read_tsv_dict_rows(Path(path), required_columns=["name", "region"])
    except ValueError as exc:
        msg = str(exc).replace("TSV", "regions file")
        sys.exit(f"ERROR: {msg}")

    for row in rows:
        if not row.get("name") or not row.get("region"):
            continue

        entry = {
            "name": row["name"],
            "region": row["region"],
            "resolution": int(row["resolution"]) if row.get("resolution") and row["resolution"].isdigit() else None,
            "label": row.get("label", ""),
            "out": row.get("out", ""),
        }
        regions.append(entry)

    if not regions:
        sys.exit("ERROR: no valid regions found in file")
    return regions


# ---------------------------------------------------------------------------
# Hi-C loading
# ---------------------------------------------------------------------------

def load_hic_matrix(hic_path, chrom, start, end, resolution):
    """
    Load a Hi-C contact matrix from a .hic file using hicstraw.
        Returns
        -------
        (mat, used_resolution)
            mat             : 2D numpy array (n_bins x n_bins)
            used_resolution : int resolution actually used
    """
    return _load_hic_matrix_impl(hic_path, chrom, start, end, resolution)


def load_cool_matrix(cool_path, chrom, start, end, resolution):
    """
    Load a Hi-C contact matrix from a .cool or .mcool file using cooler.

    For .mcool files the resolution must match one of the stored resolutions.
    For single-resolution .cool files the resolution argument is informational
    only (the file has exactly one resolution).

        Returns
        -------
        (mat, used_resolution)
            mat             : 2D numpy array (n_bins x n_bins)
            used_resolution : int resolution actually used

        Matrix values are log1p-transformed and NaN-zeroed. ICE balancing
        weights are applied when available.
    """
    return _load_cool_matrix_impl(cool_path, chrom, start, end, resolution)


def _detect_hic_format(path):
    """
    Return 'cool' for .cool/.mcool files, 'hic' for .hic files.
    Raises ValueError if the extension is unrecognised.
    """
    return _detect_hic_format_impl(path)


def load_contact_matrix(path, chrom, start, end, resolution, fmt=None):
    """
    Unified loader. fmt is 'cool' or 'hic'; if None it is auto-detected
    from the file extension.  .cool/.mcool is the default/preferred format.
    """
    return _load_contact_matrix_impl(path, chrom, start, end, resolution, fmt=fmt)


def _cool_resolution_variant(path, resolution):
    """
    If path points to a single-resolution .cool file whose filename embeds
    a bin size token, return a sibling filename adjusted to `resolution`.

    Supported patterns include:
      *_bs_10000.cool  -> *_bs_<resolution>.cool
      *_res_10000.cool -> *_res_<resolution>.cool
      *_bin_10000.cool -> *_bin_<resolution>.cool
    Returns None when no known resolution token is present.
    """
    return _cool_resolution_variant_impl(path, resolution)


def _resolve_contacts_for_region_resolution(contact_a, contact_b, resolution, hic_fmt=None):
    """
    For batch regions mode, try to switch both contact files to resolution-
    matched .cool siblings when available.

    Rule: only switch when BOTH A and B have matching existing .cool files.
    Otherwise keep original paths to avoid mixed-resolution inputs.

    Returns
    -------
    (resolved_a, resolved_b, note)
      note is a human-readable status string for logging.
    """
    return _resolve_contacts_for_region_resolution_impl(
        contact_a,
        contact_b,
        resolution,
        hic_fmt=hic_fmt,
    )


def make_split_hic(mat_a, mat_b):
    """
    Combine two Hi-C matrices into a split-diagonal view:
      - lower-left triangle = condition A
      - upper-right triangle = condition B
    The diagonal itself is averaged.
    Returns a single matrix of the same shape.
    """
    n = min(mat_a.shape[0], mat_b.shape[0])
    mat_a = mat_a[:n, :n]
    mat_b = mat_b[:n, :n]

    combined = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        for j in range(n):
            if i > j:           # lower-left: condition A
                combined[i, j] = mat_a[i, j]
            elif i < j:         # upper-right: condition B
                combined[i, j] = mat_b[i, j]
            else:               # diagonal: average
                combined[i, j] = 0.5 * (mat_a[i, j] + mat_b[i, j])
    return combined


def make_rgb_hic(mat_a, mat_b, vmax_percentile=98, palette="magenta-green"):
    """
        Build RGB Hi-C image for two-condition visualization.

        Palettes:
            - magenta-green (default): A=magenta, B=green, overlap=white, background=black
            - red-green: A=red, B=green, overlap=yellow, background=black
            - blue-red: A=blue, B=red, overlap=magenta, background=black
            - green-blue: A=green, B=blue, overlap=cyan, background=black
            - magenta-green-white: A=magenta, B=green, overlap=dark purple,
              background=white (cartographic convention). Uses bilinear interpolation
              across four corner colors to smoothly transition between conditions.

        Background is black for all palettes except magenta-green-white (white).
    """
    n = min(mat_a.shape[0], mat_b.shape[0])
    a = np.asarray(mat_a[:n, :n], dtype=np.float64)
    b = np.asarray(mat_b[:n, :n], dtype=np.float64)

    def _norm(x):
        pos = x[x > 0]
        if pos.size == 0:
            return np.zeros_like(x, dtype=np.float64)
        vmax = float(np.percentile(pos, vmax_percentile))
        if vmax <= 0:
            vmax = float(pos.max()) if pos.size else 1.0
        if vmax <= 0:
            vmax = 1.0
        y = np.clip(x / vmax, 0.0, 1.0)
        return y

    a_n = _norm(a)
    b_n = _norm(b)
    rgb = np.zeros((n, n, 3), dtype=np.uint8)

    a_u8 = (a_n * 255).astype(np.uint8)
    b_u8 = (b_n * 255).astype(np.uint8)

    if palette == "magenta-green":
        rgb[:, :, 0] = a_u8
        rgb[:, :, 1] = b_u8
        rgb[:, :, 2] = a_u8
    elif palette == "red-green":
        rgb[:, :, 0] = a_u8
        rgb[:, :, 1] = b_u8
        rgb[:, :, 2] = 0
    elif palette == "blue-red":
        rgb[:, :, 0] = b_u8
        rgb[:, :, 1] = 0
        rgb[:, :, 2] = a_u8
    elif palette == "green-blue":
        rgb[:, :, 0] = 0
        rgb[:, :, 1] = a_u8
        rgb[:, :, 2] = b_u8
    elif palette == "magenta-green-white":
        # Bivariate palette with white = no contact (cartographic convention).
        # Four corners:
        #   (a=0, b=0) → white (1.0, 1.0, 1.0) — no contact
        #   (a=1, b=0) → magenta (1.0, 0.0, 1.0) — high in A only
        #   (a=0, b=1) → green (0.0, 1.0, 0.0) — high in B only
        #   (a=1, b=1) → dark purple (0.2, 0.0, 0.3) — high in both
        # Bilinear interpolation: rgb = (1-a)(1-b)*c00 + a(1-b)*c10 + (1-a)b*c01 + ab*c11
        c00 = np.array([1.0, 1.0, 1.0], dtype=np.float64)  # white
        c10 = np.array([1.0, 0.0, 1.0], dtype=np.float64)  # magenta
        c01 = np.array([0.0, 1.0, 0.0], dtype=np.float64)  # green
        c11 = np.array([0.2, 0.0, 0.3], dtype=np.float64)  # dark purple
        
        # Compute bilinear interpolation for each pixel
        # Shape of a_n, b_n: (n, n); output shape: (n, n, 3)
        w00 = (1.0 - a_n) * (1.0 - b_n)  # shape (n, n)
        w10 = a_n * (1.0 - b_n)
        w01 = (1.0 - a_n) * b_n
        w11 = a_n * b_n
        
        rgb_float = (w00[:, :, None] * c00 +
                     w10[:, :, None] * c10 +
                     w01[:, :, None] * c01 +
                     w11[:, :, None] * c11)
        rgb = (rgb_float * 255).astype(np.uint8)
    else:
        raise ValueError(
            "Unknown RGB palette: "
            f"{palette}. Use one of: magenta-green, red-green, blue-red, green-blue, magenta-green-white"
        )
    return rgb


# ---------------------------------------------------------------------------
# qcat.bgz loading
# ---------------------------------------------------------------------------

def load_qcat_scores(qcat_path, chrom, start, end):
    """
    Load Bearing KL scores from a qcat.bgz file for a genomic region.

    Returns
    -------
    positions : list of (bin_start, bin_end) tuples within [start, end]
    score_mat : np.ndarray (n_bins, num_states) float32
    num_states: int
    """
    try:
        import pysam
    except ImportError:
        sys.exit("ERROR: pysam is required.  pip install pysam")

    tbx = pysam.TabixFile(str(qcat_path))
    positions = []
    rows = []
    num_states = 0

    # Resolve chromosome naming mismatches (chr2 vs 2) when possible.
    fetch_chrom = chrom
    contigs = set(tbx.contigs) if tbx.contigs is not None else set()
    if contigs and fetch_chrom not in contigs:
        alt = fetch_chrom[3:] if fetch_chrom.startswith("chr") else ("chr" + fetch_chrom)
        if alt in contigs:
            fetch_chrom = alt
        else:
            tbx.close()
            print(
                f"  WARNING: chromosome {chrom} not present in {qcat_path}. "
                f"Returning empty qcat track for this region."
            )
            return [], np.zeros((0, 1), dtype=np.float32), 1

    try:
        iterator = tbx.fetch(fetch_chrom, start, end)
    except Exception as e:
        tbx.close()
        print(
            f"  WARNING: could not fetch {fetch_chrom}:{start}-{end} from {qcat_path}: {e}. "
            "Returning empty qcat track for this region."
        )
        return [], np.zeros((0, 1), dtype=np.float32), 1

    for rec in iterator:
        parts = rec.split("\t")
        if len(parts) < 4:
            continue
        s, e = int(parts[1]), int(parts[2])
        meta = parts[3]
        if meta.startswith("{"):
            payload = json.loads(meta)
            pairs = payload.get("qcat", [])
        else:
            qcat_start = meta.find("qcat:")
            if qcat_start < 0:
                continue
            raw_start = meta.find(",raw:", qcat_start)
            if raw_start >= 0:
                qcat_payload = meta[qcat_start + 5:raw_start]
            else:
                qcat_payload = meta[qcat_start + 5:]
            pairs = json.loads(qcat_payload)
        max_state = max(int(p[1]) for p in pairs)
        num_states = max(num_states, max_state)
        positions.append((s, e))
        rows.append(pairs)

    tbx.close()

    if not positions:
        print(f"  WARNING: no qcat bins found in {fetch_chrom}:{start}-{end} "
              f"in {qcat_path}")
        return [], np.zeros((0, 1), dtype=np.float32), 1

    score_mat = np.zeros((len(positions), num_states), dtype=np.float32)
    for i, pairs in enumerate(rows):
        for score_val, state_idx in pairs:
            si = int(state_idx) - 1
            if 0 <= si < num_states:
                score_mat[i, si] = float(score_val)

    return positions, score_mat, num_states


# ---------------------------------------------------------------------------
# Annotation loaders
# ---------------------------------------------------------------------------

def load_loops(loops_path, chrom, start, end):
    """
    Load loop anchors from a BEDPE file within the region.
    Returns list of (s1, e1, s2, e2, score) tuples.
    """
    loops = []
    with open(loops_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue
            c1, s1, e1 = parts[0], int(parts[1]), int(parts[2])
            c2, s2, e2 = parts[3], int(parts[4]), int(parts[5])
            # Column 7+ varies by caller: some BEDPE put a numeric score there,
            # others (e.g. Mustache merged loops) put an anchor-name string like
            # 'AnchorA_3000_num365_'. Use the first parseable float at/after
            # col 6; fall back to 1.0 if none (loops are drawn unweighted).
            score = 1.0
            for tok in parts[6:]:
                try:
                    score = float(tok)
                    break
                except ValueError:
                    continue
            if c1 != chrom or c2 != chrom:
                continue
            # Keep if either anchor overlaps the region
            if e1 < start or s1 > end:
                continue
            loops.append((s1, e1, s2, e2, score))
    return loops


def load_genes(genes_path, chrom, start, end):
    """
    Load gene records from a BED6 file within the region.
    Returns list of (gene_start, gene_end, name, strand) tuples.
    """
    genes = []
    with open(genes_path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("track"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            c = parts[0]
            s, e = int(parts[1]), int(parts[2])
            name = parts[3] if len(parts) > 3 else ""
            strand = parts[5] if len(parts) > 5 else "+"
            if c != chrom:
                continue
            if e < start or s > end:
                continue
            genes.append((s, e, name, strand))
    return genes


def load_genes_gtf(gtf_path, chrom, start, end):
    """
    Load gene records from a GTF file within the region.

    Uses feature type 'gene' when available and labels entries by
    gene_name (fallback: gene_id). Coordinates are converted from GTF
    1-based inclusive to 0-based half-open for plotting compatibility.

    Returns list of (gene_start, gene_end, name, strand) tuples.
    """
    genes = []
    with open(gtf_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            c = parts[0]
            feature = parts[2]
            try:
                s = int(parts[3]) - 1
                e = int(parts[4])
            except ValueError:
                continue
            strand = parts[6] if parts[6] in {"+", "-"} else "+"
            attrs = parts[8]

            if c != chrom:
                continue
            if e < start or s > end:
                continue
            if feature != "gene":
                continue

            # Parse GTF key-value attributes: key "value";
            kv = {m.group(1): m.group(2)
                  for m in re.finditer(r'(\S+)\s+"([^"]*)"', attrs)}
            name = kv.get("gene_name") or kv.get("gene_id") or ""
            genes.append((max(0, s), e, name, strand))
    return genes


def load_highlights(highlights_path, chrom, start, end):
    """
    Load highlight regions from a BED4 file (col 4 = hex color).
    Returns list of (region_start, region_end, color) tuples.
    """
    highlights = []
    with open(highlights_path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("track"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            c = parts[0]
            s, e = int(parts[1]), int(parts[2])
            color = parts[3] if len(parts) > 3 else "#ffff00"
            if c != chrom:
                continue
            if e < start or s > end:
                continue
            highlights.append((s, e, color))
    return highlights


# ---------------------------------------------------------------------------
# Drawing helpers
# ---------------------------------------------------------------------------

HIC_CMAP = LinearSegmentedColormap.from_list(
    "hic_red", ["#ffffff", "#ffd0d0", "#ff6666", "#cc0000", "#660000"]
)


def genomic_to_ax(pos, region_start, region_end):
    """Convert a genomic coordinate to [0, 1] axis fraction."""
    return (pos - region_start) / (region_end - region_start)


def add_vertical_tick_bars(ax, region_start, region_end):
    """Draw unlabeled genomic tick bars on a vertical track's Y axis."""
    span = region_end - region_start
    for interval in [5e6, 2e6, 1e6, 5e5, 2e5, 1e5, 5e4, 1e4]:
        if 4 <= span / interval <= 12:
            break

    ticks_genomic = []
    pos = (region_start // int(interval) + 1) * int(interval)
    while pos < region_end:
        ticks_genomic.append(pos)
        pos += int(interval)

    ticks_ax = [genomic_to_ax(t, region_start, region_end)
                for t in ticks_genomic]
    ax.set_yticks(ticks_ax)
    ax.set_yticklabels([])
    ax.tick_params(axis="y", left=True, right=False,
                   length=2.5, width=0.6, colors="#777777")


def draw_epilogos_horizontal(ax, positions, score_mat, num_states,
                              region_start, region_end,
                              categories, highlights=None, label="",
                              y_max=None):
    """
    Draw a horizontal epilogos stacked bar track on ax.
    Each bin is a stacked bar of state scores, colored by category.
    Bins span the x axis from 0 to 1 (genomic fraction).
    """
    ax.set_xlim(0, 1)
    _y_max = y_max if y_max is not None else (score_mat.max() * 1.15 if score_mat.size > 0 else 1)
    ax.set_ylim(0, _y_max)

    # Highlights first (behind bars)
    add_horizontal_highlights(ax, highlights, region_start, region_end, genomic_to_ax)

    for i, (bs, be) in enumerate(positions):
        x0 = genomic_to_ax(bs, region_start, region_end)
        x1 = genomic_to_ax(be, region_start, region_end)
        width = max(x1 - x0, 0.0002)
        row = score_mat[i]
        # Stack states sorted ascending by index (consistent stacking order)
        bottom = 0.0
        for si in range(num_states):
            val = float(row[si])
            if val <= 0:
                continue
            color = (categories[si][1] if si < len(categories)
                     else "#cccccc")
            ax.bar(x0, val, width=width, bottom=bottom,
                   color=color, align="edge",
                   linewidth=0, zorder=2)
            bottom += val

    ax.set_ylabel(label, fontsize=7, labelpad=3)
    ax.tick_params(axis="both", labelsize=6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xticks([])


def draw_epilogos_vertical(ax, positions, score_mat, num_states,
                            region_start, region_end,
                            categories, highlights=None, label="",
                            y_max=None):
    """
    Draw a vertical epilogos track (rotated 90 degrees clockwise).
    Genomic position runs bottom-to-top on the y axis.
    Scores extend to the right on the x axis.
    """
    # Y axis inverted: top = region_start, bottom = region_end (matches Hi-C)
    # X axis inverted: bars extend left toward the Hi-C square
    ax.set_ylim(1, 0)
    _y_max = y_max if y_max is not None else (score_mat.max() * 1.15 if score_mat.size > 0 else 1)
    ax.set_xlim(_y_max, 0)

    add_vertical_highlights(ax, highlights, region_start, region_end, genomic_to_ax)

    for i, (bs, be) in enumerate(positions):
        y0 = genomic_to_ax(bs, region_start, region_end)
        y1 = genomic_to_ax(be, region_start, region_end)
        height = max(y1 - y0, 0.0002)
        row = score_mat[i]
        left = 0.0
        for si in range(num_states):
            val = float(row[si])
            if val <= 0:
                continue
            color = (categories[si][1] if si < len(categories)
                     else "#cccccc")
            ax.barh(y0, val, height=height, left=left,
                    color=color, align="edge",
                    linewidth=0, zorder=2)
            left += val

    ax.set_xlabel(label, fontsize=7, labelpad=3)
    ax.tick_params(axis="both", labelsize=6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    # Left spine acts as separator between Hi-C square and this track
    ax.spines["left"].set_visible(True)
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["left"].set_color("#aaaaaa")
    add_vertical_tick_bars(ax, region_start, region_end)


def draw_diff_horizontal(ax, positions, score_mat, num_states,
                          region_start, region_end,
                          categories, highlights=None, label="",
                          diff_max=None):
    """
    Draw a horizontal differential epilogos track.
    Y axis is inverted: negative scores (enriched in B) extend upward,
    positive scores (enriched in A) extend downward from zero.
    This mirrors the convention used in the vertical diff track where
    positive extends toward the Hi-C square (left/up = more in B,
    right/down = more in A).
    A zero reference line is drawn.
    """
    if score_mat.size == 0:
        ax.set_axis_off()
        return

    if diff_max is not None:
        abs_max = diff_max
    else:
        abs_max = np.abs(score_mat).max()
        if abs_max == 0:
            abs_max = 1.0
    ax.set_xlim(0, 1)
    # Inverted: negative values extend upward, positive downward
    ax.set_ylim(abs_max * 1.15, -abs_max * 1.15)

    # Zero reference line
    ax.axhline(0, color="#888888", linewidth=0.6, zorder=1)

    add_horizontal_highlights(ax, highlights, region_start, region_end, genomic_to_ax)

    for i, (bs, be) in enumerate(positions):
        x0 = genomic_to_ax(bs, region_start, region_end)
        x1 = genomic_to_ax(be, region_start, region_end)
        width = max(x1 - x0, 0.0002)
        row = score_mat[i]
        pos_bottom = 0.0
        neg_bottom = 0.0
        for si in range(num_states):
            val = float(row[si])
            if val == 0:
                continue
            color = (categories[si][1] if si < len(categories)
                     else "#cccccc")
            if val > 0:
                ax.bar(x0, val, width=width, bottom=pos_bottom,
                       color=color, align="edge", linewidth=0, zorder=2)
                pos_bottom += val
            else:
                ax.bar(x0, val, width=width, bottom=neg_bottom,
                       color=color, align="edge", linewidth=0, zorder=2)
                neg_bottom += val

    # Subtle gray baseline fill to make zero-crossing visible
    if positions:
        xs = [genomic_to_ax(bs, region_start, region_end) for bs, be in positions]
        totals = [sum(max(float(row[si]), 0) for si in range(num_states))
                  for row in score_mat]
        ax.fill_between(xs, 0, totals, alpha=0.06, color="#444444",
                        step="post", zorder=1)

    ax.set_ylabel(label, fontsize=7, labelpad=3)
    ax.tick_params(axis="both", labelsize=6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xticks([])


def draw_diff_vertical(ax, positions, score_mat, num_states,
                        region_start, region_end,
                        categories, highlights=None, label="",
                        diff_max=None):
    """
    Draw a vertical differential epilogos track.
    Genomic position runs top-to-bottom on the Y axis (matches Hi-C).
    X axis is inverted: negative scores (enriched in B) extend to the right
    toward the Hi-C square, positive scores (enriched in A) extend left.
    A zero reference line is drawn at x=0.
    """
    if score_mat.size == 0:
        ax.set_axis_off()
        return

    if diff_max is not None:
        abs_max = diff_max
    else:
        abs_max = np.abs(score_mat).max()
        if abs_max == 0:
            abs_max = 1.0

    # Y axis inverted: top = region_start, bottom = region_end (matches Hi-C)
    ax.set_ylim(1, 0)
    # X axis inverted: negative extends right (toward Hi-C), positive left
    ax.set_xlim(abs_max * 1.15, -abs_max * 1.15)

    # Zero reference line
    ax.axvline(0, color="#888888", linewidth=0.6, zorder=1)

    add_vertical_highlights(ax, highlights, region_start, region_end, genomic_to_ax)

    for i, (bs, be) in enumerate(positions):
        y0 = genomic_to_ax(bs, region_start, region_end)
        y1 = genomic_to_ax(be, region_start, region_end)
        height = max(y1 - y0, 0.0002)
        row = score_mat[i]
        pos_left = 0.0
        neg_left = 0.0
        for si in range(num_states):
            val = float(row[si])
            if val == 0:
                continue
            color = (categories[si][1] if si < len(categories)
                     else "#cccccc")
            if val > 0:
                ax.barh(y0, val, height=height, left=pos_left,
                        color=color, align="edge", linewidth=0, zorder=2)
                pos_left += val
            else:
                ax.barh(y0, val, height=height, left=neg_left,
                        color=color, align="edge", linewidth=0, zorder=2)
                neg_left += val

    ax.set_xlabel(label, fontsize=7, labelpad=3)
    ax.tick_params(axis="both", labelsize=6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["left"].set_color("#aaaaaa")
    ax.spines["bottom"].set_visible(False)
    add_vertical_tick_bars(ax, region_start, region_end)


def draw_hic_square(ax, mat_a, mat_b, combined_mat, resolution,
                    region_start, region_end,
                    loops=None, highlights=None,
                    label_a="Condition A", label_b="Condition B",
                    rgb_mode=False, rgb_palette="magenta-green"):
    """
    Draw the split Hi-C contact matrix on ax.
    Lower-left = condition A, upper-right = condition B.
    Diagonal separator line drawn in black.
    Optional loop arcs overlaid.
    """
    if rgb_mode:
        rgb = make_rgb_hic(mat_a, mat_b, palette=rgb_palette)
        n = rgb.shape[0]
        ax.imshow(
            rgb,
            origin="upper",
            extent=[0, 1, 1, 0],
            aspect="auto",
            interpolation="nearest",
            zorder=1,
        )
        vmax = None
    else:
        n = combined_mat.shape[0]
        vmax = np.percentile(combined_mat[combined_mat > 0], 95) if combined_mat.any() else 1
        ax.imshow(
            combined_mat,
            cmap=HIC_CMAP,
            vmin=0, vmax=vmax,
            origin="upper",
            extent=[0, 1, 1, 0],   # (left, right, bottom, top) in axis coords
            aspect="auto",
            interpolation="nearest",
            zorder=1,
        )

    # Diagonal separator
    ax.plot([0, 1], [0, 1], color="#222222", linewidth=0.8,
            linestyle="--", zorder=4)

    # Condition labels: A = lower-left, B = upper-right.
    # With origin="upper" and ylim(1,0): axes (0,1)=lower-left, (1,0)=upper-right.
    # bbox adds a semi-transparent white background for readability over the matrix.
    _label_bbox = dict(boxstyle="round,pad=0.2", facecolor="white",
                       alpha=0.7, edgecolor="none")
    ax.text(0.04, 0.04, label_a, transform=ax.transAxes,
            fontsize=7, color="#222222", ha="left", va="bottom",
            fontweight="bold", zorder=6, bbox=_label_bbox)
    ax.text(0.96, 0.96, label_b, transform=ax.transAxes,
            fontsize=7, color="#222222", ha="right", va="top",
            fontweight="bold", zorder=6, bbox=_label_bbox)

    # Highlights as vertical + horizontal bands
    if highlights:
        for hs, he, hcol in highlights:
            x0 = genomic_to_ax(max(hs, region_start), region_start, region_end)
            x1 = genomic_to_ax(min(he, region_end),   region_start, region_end)
            ax.axvspan(x0, x1, color=hcol, alpha=0.12, zorder=3)
            ax.axhspan(x0, x1, color=hcol, alpha=0.12, zorder=3)

    # Loop arcs
    if loops:
        max_score = max(lp[4] for lp in loops) if loops else 1.0
        for s1, e1, s2, e2, score in loops:
            mid1 = genomic_to_ax(0.5 * (s1 + e1), region_start, region_end)
            mid2 = genomic_to_ax(0.5 * (s2 + e2), region_start, region_end)
            # Draw as a point on the matrix (loop anchor intersection)
            alpha = 0.4 + 0.5 * (score / max_score)
            ax.plot(mid2, mid1, "o", color="#1a1a1a",
                    markersize=3.5, alpha=alpha, zorder=5)
            ax.plot(mid1, mid2, "o", color="#1a1a1a",
                    markersize=3.5, alpha=alpha, zorder=5)
            # Connect with a thin cross-hair line
            ax.plot([mid2, mid2], [mid1 - 0.01, mid1 + 0.01],
                    color="#1a1a1a", lw=0.6, alpha=alpha, zorder=5)
            ax.plot([mid2 - 0.01, mid2 + 0.01], [mid1, mid1],
                    color="#1a1a1a", lw=0.6, alpha=alpha, zorder=5)

    # --- Genomic ticks on top (x) and left (y) edges ---
    span = region_end - region_start
    for interval in [5e6, 2e6, 1e6, 5e5, 2e5, 1e5, 5e4, 1e4]:
        if 4 <= span / interval <= 12:
            break
    ticks_genomic = []
    pos = (region_start // int(interval) + 1) * int(interval)
    while pos < region_end:
        ticks_genomic.append(pos)
        pos += int(interval)
    ticks_ax = [genomic_to_ax(t, region_start, region_end)
                for t in ticks_genomic]
    tick_labels = [f"{t/1e6:.2f}" for t in ticks_genomic]

    # Top spine (x axis) -- labels above the matrix
    ax.set_xticks(ticks_ax)
    ax.set_xticklabels(tick_labels, fontsize=5, rotation=45, ha="left")
    ax.xaxis.set_label_position("top")
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", length=3, width=0.5, labelsize=5, top=True,
                   bottom=False, colors="#444444")

    # Left spine (y axis) -- labels to the left of the matrix
    ax.set_yticks(ticks_ax)
    ax.set_yticklabels(tick_labels, fontsize=5, ha="right")
    ax.yaxis.set_label_position("left")
    ax.yaxis.tick_left()
    ax.tick_params(axis="y", length=3, width=0.5, labelsize=5, left=True,
                   right=False, colors="#444444")

    ax.set_xlim(0, 1)
    ax.set_ylim(1, 0)
    ax.spines["top"].set_visible(True)
    ax.spines["top"].set_linewidth(0.5)
    ax.spines["top"].set_color("#444444")
    ax.spines["left"].set_visible(True)
    ax.spines["left"].set_linewidth(0.5)
    ax.spines["left"].set_color("#444444")
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)

    return vmax


def draw_gene_track(ax, genes, region_start, region_end,
                    highlights=None, label="Genes"):
    """
    Draw a simple gene body track below the epilogos.
    Genes shown as thick arrows indicating strand direction.
    """
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.5, 1.5)

    add_horizontal_highlights(ax, highlights, region_start, region_end, genomic_to_ax)

    y_levels = [0.2, 0.8]
    used = []

    for gs, ge, name, strand in genes:
        x0 = genomic_to_ax(max(gs, region_start), region_start, region_end)
        x1 = genomic_to_ax(min(ge, region_end),   region_start, region_end)
        if x1 <= x0:
            continue

        # Assign y level to avoid overlap
        yl = 0.2
        for prev_x0, prev_x1, prev_y in used:
            if not (x1 < prev_x0 or x0 > prev_x1) and prev_y == 0.2:
                yl = 0.8
                break
        used.append((x0, x1, yl))

        # Gene body
        ax.plot([x0, x1], [yl, yl], color="#555555", lw=2.5,
                solid_capstyle="butt", zorder=2)

        # Arrow for strand
        arrow_x = x1 if strand == "+" else x0
        dx = 0.01 if strand == "+" else -0.01
        ax.annotate("",
            xy=(arrow_x + dx, yl),
            xytext=(arrow_x, yl),
            arrowprops=dict(arrowstyle="->", color="#555555",
                           lw=1.0, mutation_scale=6),
            zorder=3,
        )

        # Label -- only draw if gene is wide enough to avoid clutter
        # Skip label if gene span < 0.4% of the region (avoids dense clusters)
        if (x1 - x0) >= 0.004:
            mid_x = 0.5 * (x0 + x1)
            ax.text(mid_x, yl + 0.15, name,
                    fontsize=6, ha="center", va="bottom",
                    color="#333333", style="italic",
                    rotation=0, zorder=3)

    ax.set_ylabel(label, fontsize=7, labelpad=3)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)


def draw_legend(ax, categories, num_states, hic_vmax=None,
                rgb_mode=False, label_a="Condition A", label_b="Condition B",
                rgb_palette="magenta-green"):
    """
    Draw a category color legend on a blank axes.
    If hic_vmax is provided, a Hi-C contact colorbar is drawn at the top
    of the legend panel above the chromatin state swatches.
    """
    draw_legend_panel(
        ax,
        categories,
        num_states,
        hic_cmap=HIC_CMAP,
        hic_vmax=hic_vmax,
        rgb_mode=rgb_mode,
        label_a=label_a,
        label_b=label_b,
        rgb_palette=rgb_palette,
    )


def draw_genomic_axis(ax, region_start, region_end, chrom):
    """
    Draw a simple genomic coordinate axis.
    """
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    span = region_end - region_start
    # Pick a sensible tick interval
    for interval in [5e6, 2e6, 1e6, 5e5, 2e5, 1e5, 5e4, 1e4]:
        n_ticks = span / interval
        if 4 <= n_ticks <= 12:
            break

    ticks_genomic = []
    pos = (region_start // int(interval) + 1) * int(interval)
    while pos < region_end:
        ticks_genomic.append(pos)
        pos += int(interval)

    ticks_ax = [genomic_to_ax(t, region_start, region_end)
                for t in ticks_genomic]

    ax.set_xticks(ticks_ax)
    ax.set_xticklabels(
        [f"{t/1e6:.2f}" for t in ticks_genomic],
        fontsize=6, rotation=45, ha="right",
    )
    ax.set_yticks([])
    ax.set_xlabel(f"{chrom} (Mb)", fontsize=7, labelpad=2)
    ax.tick_params(axis="x", length=3, width=0.6, labelsize=6)
    ax.spines["bottom"].set_linewidth(0.6)


def draw_gene_track_vertical(ax, genes, region_start, region_end,
                              highlights=None, label="Genes"):
    """
    Vertical gene track: genomic position on Y axis (bottom=start, top=end),
    gene bodies drawn as horizontal bands extending to the right.
    """
    # Y axis inverted: top = region_start, bottom = region_end (matches Hi-C)
    ax.set_ylim(1, 0)
    # Three x-levels for gene bodies, kept close to the right-side axis column.
    x_levels = [0.64, 0.79, 0.92]
    ax.set_xlim(0.56, 1.0)

    add_vertical_highlights(ax, highlights, region_start, region_end, genomic_to_ax)

    # used: list of (y0, y1, xl) for overlap detection across all levels
    used = []

    for gs, ge, name, strand in genes:
        y0 = genomic_to_ax(max(gs, region_start), region_start, region_end)
        y1 = genomic_to_ax(min(ge, region_end),   region_start, region_end)
        if y1 <= y0:
            continue

        # Assign the first x level that has no vertical overlap
        xl = x_levels[0]
        for lvl in x_levels:
            overlap = any(
                not (y1 < py0 or y0 > py1) and px == lvl
                for py0, py1, px in used
            )
            if not overlap:
                xl = lvl
                break

        used.append((y0, y1, xl))

        ax.plot([xl, xl], [y0, y1], color="#555555", lw=2.0,
                solid_capstyle="butt", zorder=2)

        # Arrow: + strand points downward (increasing genomic pos = down in axes)
        arrow_y = y1 if strand == "+" else y0
        dy = 0.01 if strand == "+" else -0.01
        ax.annotate("",
            xy=(xl, arrow_y + dy),
            xytext=(xl, arrow_y),
            arrowprops=dict(arrowstyle="->", color="#555555",
                           lw=1.0, mutation_scale=6),
            zorder=3,
        )

        mid_y = 0.5 * (y0 + y1)
        # Label only if gene is tall enough to avoid clutter
        if (y1 - y0) >= 0.004:
            ax.text(max(0.57, xl - 0.025), mid_y, name,
                fontsize=4, ha="right", va="center",
                    color="#333333", style="italic",
                    rotation=90, zorder=3)

    ax.set_xlabel(label, fontsize=7, labelpad=3)
    ax.set_xticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["left"].set_color("#aaaaaa")
    add_vertical_tick_bars(ax, region_start, region_end)


def draw_genomic_axis_vertical(ax, region_start, region_end, chrom):
    """
    Vertical genomic axis: coordinates run top (start) to bottom (end),
    matching the Hi-C matrix orientation (origin at top-left).
    Tick labels and ylabel are horizontal for readability.
    """
    # Inverted: top = region_start, bottom = region_end
    ax.set_ylim(1, 0)
    ax.set_xlim(0, 1)
    for sp in ["top", "right", "bottom"]:
        ax.spines[sp].set_visible(False)
    ax.spines["left"].set_linewidth(0.6)

    span = region_end - region_start
    for interval in [5e6, 2e6, 1e6, 5e5, 2e5, 1e5, 5e4, 1e4]:
        if 4 <= span / interval <= 12:
            break

    ticks_genomic = []
    pos = (region_start // int(interval) + 1) * int(interval)
    while pos < region_end:
        ticks_genomic.append(pos)
        pos += int(interval)

    ticks_ax = [genomic_to_ax(t, region_start, region_end)
                for t in ticks_genomic]

    ax.set_yticks(ticks_ax)
    ax.set_yticklabels(
        [f"{t/1e6:.2f}" for t in ticks_genomic],
        fontsize=5, rotation=45, ha="left",
    )
    ax.set_xticks([])
    ax.set_ylabel("Mb", fontsize=6, labelpad=2)
    ax.yaxis.set_label_position("right")
    ax.yaxis.tick_right()
    ax.tick_params(axis="y", length=3, width=0.6, labelsize=5, pad=1)


# ---------------------------------------------------------------------------
# P-value track helpers
# ---------------------------------------------------------------------------

def load_bigwig_values(bw_path, chrom, start, end, n_bins=None):
    """
    Load per-bin mean values from a BigWig file over a genomic region.
    Returns (positions, values) as 1-D numpy arrays.
    positions[i] is the genomic start of bin i; values[i] is the mean signal.
    Missing data (NaN) is returned as 0.0.
    """
    try:
        import pyBigWig
    except ImportError:
        sys.exit("ERROR: pyBigWig is required for --pval-a/--pval-b.\n"
                 "Install with:  pip install pyBigWig")
    bw = pyBigWig.open(str(bw_path))
    if n_bins is None:
        n_bins = max(1, (end - start) // 200)
    try:
        raw = bw.stats(chrom, start, end, type="mean", nBins=n_bins)
    except Exception:
        raw = [None] * n_bins
    bw.close()
    vals = np.array([v if v is not None else 0.0 for v in raw],
                    dtype=np.float64)
    bin_size = (end - start) / n_bins
    positions = np.array([start + i * bin_size for i in range(n_bins)])
    return positions, vals


def load_stats_tsv_values(tsv_path, chrom, start, end):
    """
    Load signed -log10(p-value) values from a stats TSV file.

    The TSV is expected to contain at least chrom, start, end, and pval
    columns. If a direction column is present, it is used to sign the values.
    Otherwise a signed bearing_score column is used when available.
    """
    candidates = {chrom}
    if chrom.startswith("chr"):
        candidates.add(chrom[3:])
    else:
        candidates.add("chr" + chrom)

    positions = []
    values = []
    with open(tsv_path, newline="") as fh:
        reader = csv.DictReader(iter_non_comment_lines(fh), delimiter="\t")
        if reader.fieldnames is None:
            return np.array([], dtype=np.float64), np.array([], dtype=np.float64)

        field_map = {
            (name.strip().lower() if name is not None else ""): name
            for name in reader.fieldnames
            if name is not None
        }
        required = ["chrom", "start", "end", "pval"]
        missing = [col for col in required if col not in field_map]
        if missing:
            raise SystemExit(
                f"ERROR: {tsv_path} is missing required columns for TSV p-value input: {', '.join(missing)}"
            )

        chrom_col = field_map["chrom"]
        start_col = field_map["start"]
        end_col = field_map["end"]
        pval_col = field_map["pval"]
        direction_col = field_map.get("direction")
        score_col = field_map.get("bearing_score")

        for row in reader:
            row_chrom = (row.get(chrom_col) or "").strip()
            if row_chrom not in candidates:
                continue

            try:
                row_start = int(float(row.get(start_col, "")))
                row_end = int(float(row.get(end_col, "")))
            except (TypeError, ValueError):
                continue

            if row_end <= start or row_start >= end:
                continue

            try:
                pval = float(row.get(pval_col, ""))
            except (TypeError, ValueError):
                continue

            if not np.isfinite(pval):
                continue

            safe_pval = max(pval, 1e-300)
            sign = 1.0
            if direction_col is not None:
                direction = (row.get(direction_col) or "").strip()
                if direction.startswith("-"):
                    sign = -1.0
                elif direction.startswith("+"):
                    sign = 1.0
            elif score_col is not None:
                try:
                    score = float(row.get(score_col, ""))
                except (TypeError, ValueError):
                    score = 0.0
                if score < 0:
                    sign = -1.0

            positions.append(row_start)
            values.append(sign * (-math.log10(safe_pval)))

    if not positions:
        return np.array([], dtype=np.float64), np.array([], dtype=np.float64)

    order = np.argsort(np.asarray(positions, dtype=np.float64))
    return (
        np.asarray(positions, dtype=np.float64)[order],
        np.asarray(values, dtype=np.float64)[order],
    )


def load_stats_tsv_with_categories(tsv_path, chrom, start, end):
    """
    Load p-values and per-category KL scores from a stats TSV file.
    
    Returns (positions, values, kl_scores, category_names) where:
      positions      : numpy array of bin start positions
      values         : numpy array of signed -log10(p-values)
      kl_scores      : numpy array of shape (n_bins, n_categories) with per-category KL values
      category_names : list of category names from kl_<name> columns
    """
    candidates = {chrom}
    if chrom.startswith("chr"):
        candidates.add(chrom[3:])
    else:
        candidates.add("chr" + chrom)

    positions = []
    values = []
    kl_rows = []
    category_names = []
    
    with open(tsv_path, newline="") as fh:
        reader = csv.DictReader(iter_non_comment_lines(fh), delimiter="\t")
        if reader.fieldnames is None:
            return (
                np.array([], dtype=np.float64),
                np.array([], dtype=np.float64),
                np.array([], dtype=np.float64).reshape(0, 0),
                [],
            )

        field_map = {
            (name.strip().lower() if name is not None else ""): name
            for name in reader.fieldnames
            if name is not None
        }
        
        # Extract category columns (kl_<name>)
        kl_cols = []
        for fname in reader.fieldnames:
            if fname is not None and fname.lower().startswith("kl_"):
                category_names.append(fname[3:])  # Strip "kl_" prefix
                kl_cols.append(fname)
        
        required = ["chrom", "start", "end", "pval"]
        missing = [col for col in required if col not in field_map]
        if missing:
            raise SystemExit(
                f"ERROR: {tsv_path} is missing required columns for TSV p-value input: {', '.join(missing)}"
            )

        chrom_col = field_map["chrom"]
        start_col = field_map["start"]
        end_col = field_map["end"]
        pval_col = field_map["pval"]
        direction_col = field_map.get("direction")
        score_col = field_map.get("bearing_score")

        for row in reader:
            row_chrom = (row.get(chrom_col) or "").strip()
            if row_chrom not in candidates:
                continue

            try:
                row_start = int(float(row.get(start_col, "")))
                row_end = int(float(row.get(end_col, "")))
            except (TypeError, ValueError):
                continue

            if row_end <= start or row_start >= end:
                continue

            try:
                pval = float(row.get(pval_col, ""))
            except (TypeError, ValueError):
                continue

            if not np.isfinite(pval):
                continue

            safe_pval = max(pval, 1e-300)
            sign = 1.0
            if direction_col is not None:
                direction = (row.get(direction_col) or "").strip()
                if direction.startswith("-"):
                    sign = -1.0
                elif direction.startswith("+"):
                    sign = 1.0
            elif score_col is not None:
                try:
                    score = float(row.get(score_col, ""))
                except (TypeError, ValueError):
                    score = 0.0
                if score < 0:
                    sign = -1.0

            positions.append(row_start)
            values.append(sign * (-math.log10(safe_pval)))
            
            # Load per-category KL scores
            kl_row = []
            for col in kl_cols:
                try:
                    kl_val = float(row.get(col, "0"))
                except (TypeError, ValueError):
                    kl_val = 0.0
                kl_row.append(kl_val)
            kl_rows.append(kl_row)

    if not positions:
        return (
            np.array([], dtype=np.float64),
            np.array([], dtype=np.float64),
            np.array([], dtype=np.float64).reshape(0, 0),
            category_names,
        )

    order = np.argsort(np.asarray(positions, dtype=np.float64))
    pos_array = np.asarray(positions, dtype=np.float64)[order]
    val_array = np.asarray(values, dtype=np.float64)[order]
    kl_array = np.asarray(kl_rows, dtype=np.float64)[order] if kl_rows else np.array([]).reshape(len(positions), 0)
    
    return pos_array, val_array, kl_array, category_names


def load_pval_track_values(track_path, chrom, start, end, n_bins=None):
    """Load p-value track values from either a BigWig or a stats TSV."""
    track_path = Path(track_path)
    if track_path.suffix.lower() == ".tsv":
        return load_stats_tsv_values(track_path, chrom, start, end)
    return load_bigwig_values(track_path, chrom, start, end, n_bins=n_bins)


def draw_pval_horizontal(ax, positions, values, region_start, region_end,
                         highlights=None, label="", color="#a43cca",
                         y_max=None, cutoff_value=None):
    """
    Horizontal -log10(p-value) track drawn as a filled area.
    Aligned with the epilogos track above it, with inverted y so values
    extend downward away from qcat bars.
    y_max (optional) forces the vertical scale; otherwise determined from
    values.
    """
    ax.set_xlim(0, 1)
    if y_max is None:
        vmax = float(values.max()) * 1.1 if values.size > 0 and values.max() > 0 else 1.0
    else:
        vmax = y_max
    ax.set_ylim(vmax, 0)

    add_horizontal_highlights(ax, highlights, region_start, region_end, genomic_to_ax)

    xs = np.array([genomic_to_ax(p, region_start, region_end)
                   for p in positions])
    ax.fill_between(xs, 0, values, color=color, alpha=0.65,
                    linewidth=0, zorder=2)
    ax.plot(xs, values, color=color, linewidth=0.6, zorder=3)

    if cutoff_value is not None and cutoff_value >= 0:
        ax.axhline(cutoff_value, color="#666666", linestyle="--",
                   linewidth=0.8, zorder=4)

    ax.set_ylabel(label, fontsize=7, labelpad=3)
    ax.tick_params(axis="both", labelsize=6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xticks([])


def draw_pval_vertical(ax, positions, values, region_start, region_end,
                       highlights=None, label="", color="#a43cca",
                       x_max=None, cutoff_value=None):
    """
    Vertical -log10(p-value) track (values extend right, away from epilogos).
    Genomic position runs top-to-bottom on the Y axis (matches Hi-C orientation).
    x_max (optional) forces horizontal scale.
    """
    ax.set_ylim(1, 0)
    if x_max is None:
        vmax = float(values.max()) * 1.1 if values.size > 0 and values.max() > 0 else 1.0
    else:
        vmax = x_max
    ax.set_xlim(0, vmax)

    add_vertical_highlights(ax, highlights, region_start, region_end, genomic_to_ax)

    ys = np.array([genomic_to_ax(p, region_start, region_end)
                   for p in positions])
    ax.fill_betweenx(ys, 0, values, color=color, alpha=0.65,
                     linewidth=0, zorder=2)
    ax.plot(values, ys, color=color, linewidth=0.6, zorder=3)

    if cutoff_value is not None and cutoff_value >= 0:
        ax.axvline(cutoff_value, color="#666666", linestyle="--",
                   linewidth=0.8, zorder=4)

    ax.set_xlabel(label, fontsize=7, labelpad=3)
    ax.tick_params(axis="both", labelsize=6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["left"].set_color("#aaaaaa")
    add_vertical_tick_bars(ax, region_start, region_end)


def draw_pval_diff_horizontal(ax, positions, values, region_start, region_end,
                              highlights=None, label="",
                              pos_color="#0c4a6e", neg_color="#d97706",
                              y_max=None, cutoff_value=None,
                              diff_score_positions=None, diff_score_matrix=None,
                              categories=None):
    return _draw_pval_diff_horizontal_impl(
        ax, positions, values, region_start, region_end,
        highlights=highlights, label=label,
        pos_color=pos_color, neg_color=neg_color,
        y_max=y_max, cutoff_value=cutoff_value,
        diff_score_positions=diff_score_positions,
        diff_score_matrix=diff_score_matrix,
        categories=categories,
    )


def draw_pval_diff_vertical(ax, positions, values, region_start, region_end,
                            highlights=None, label="",
                            pos_color="#0c4a6e", neg_color="#d97706",
                            x_max=None, cutoff_value=None,
                            diff_score_positions=None, diff_score_matrix=None,
                            categories=None):
    return _draw_pval_diff_vertical_impl(
        ax, positions, values, region_start, region_end,
        highlights=highlights, label=label,
        pos_color=pos_color, neg_color=neg_color,
        x_max=x_max, cutoff_value=cutoff_value,
        diff_score_positions=diff_score_positions,
        diff_score_matrix=diff_score_matrix,
        categories=categories,
    )


def draw_loops_horizontal(ax, loops, region_start, region_end,
                          highlights=None, label="", color="#4c78a8",
                          anchor_color="#4c78a8"):
    """Draw a horizontal loop-arc track with anchor-region highlighting."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    if highlights:
        for hs, he, hcol in highlights:
            x0 = genomic_to_ax(max(hs, region_start), region_start, region_end)
            x1 = genomic_to_ax(min(he, region_end),   region_start, region_end)
            ax.axvspan(x0, x1, color=hcol, alpha=0.18, zorder=0)

    if loops:
        max_dist = max(abs((0.5 * (s2 + e2)) - (0.5 * (s1 + e1)))
                       for s1, e1, s2, e2, _ in loops)
        max_dist = max(max_dist, 1.0)
        max_score = max(lp[4] for lp in loops)
        max_score = max(max_score, 1e-9)

        for s1, e1, s2, e2, score in loops:
            a1s = genomic_to_ax(max(s1, region_start), region_start, region_end)
            a1e = genomic_to_ax(min(e1, region_end),   region_start, region_end)
            a2s = genomic_to_ax(max(s2, region_start), region_start, region_end)
            a2e = genomic_to_ax(min(e2, region_end),   region_start, region_end)
            if a1e <= 0 or a2e <= 0 or a1s >= 1 or a2s >= 1:
                continue

            alpha = 0.35 + 0.55 * min(1.0, score / max_score)

            # Highlight anchor intervals
            ax.axvspan(a1s, a1e, color=anchor_color, alpha=alpha, zorder=1)
            ax.axvspan(a2s, a2e, color=anchor_color, alpha=alpha, zorder=1)

            m1 = genomic_to_ax(0.5 * (s1 + e1), region_start, region_end)
            m2 = genomic_to_ax(0.5 * (s2 + e2), region_start, region_end)
            left, right = sorted([m1, m2])
            width = right - left
            if width <= 0:
                continue
            dist_bp = abs((0.5 * (s2 + e2)) - (0.5 * (s1 + e1)))
            height = 0.10 + 0.65 * (dist_bp / max_dist)

            arc = Arc(((left + right) / 2.0, 0.0),
                      width=width, height=height,
                      angle=0, theta1=0, theta2=180,
                      lw=0.9, color=color, alpha=alpha, zorder=3)
            ax.add_patch(arc)

    ax.set_ylabel(label, fontsize=7, labelpad=3)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)


def draw_loops_vertical(ax, loops, region_start, region_end,
                        highlights=None, label="", color="#e05c3a",
                        anchor_color="#e05c3a"):
    """Draw a vertical loop-arc track with anchor-region highlighting."""
    ax.set_ylim(1, 0)
    ax.set_xlim(0, 1)

    if highlights:
        for hs, he, hcol in highlights:
            y0 = genomic_to_ax(max(hs, region_start), region_start, region_end)
            y1 = genomic_to_ax(min(he, region_end),   region_start, region_end)
            ax.axhspan(y0, y1, color=hcol, alpha=0.18, zorder=0)

    if loops:
        max_dist = max(abs((0.5 * (s2 + e2)) - (0.5 * (s1 + e1)))
                       for s1, e1, s2, e2, _ in loops)
        max_dist = max(max_dist, 1.0)
        max_score = max(lp[4] for lp in loops)
        max_score = max(max_score, 1e-9)

        for s1, e1, s2, e2, score in loops:
            a1s = genomic_to_ax(max(s1, region_start), region_start, region_end)
            a1e = genomic_to_ax(min(e1, region_end),   region_start, region_end)
            a2s = genomic_to_ax(max(s2, region_start), region_start, region_end)
            a2e = genomic_to_ax(min(e2, region_end),   region_start, region_end)
            if a1e <= 0 or a2e <= 0 or a1s >= 1 or a2s >= 1:
                continue

            alpha = 0.35 + 0.55 * min(1.0, score / max_score)

            # Highlight anchor intervals
            ax.axhspan(a1s, a1e, color=anchor_color, alpha=alpha, zorder=1)
            ax.axhspan(a2s, a2e, color=anchor_color, alpha=alpha, zorder=1)

            m1 = genomic_to_ax(0.5 * (s1 + e1), region_start, region_end)
            m2 = genomic_to_ax(0.5 * (s2 + e2), region_start, region_end)
            y1, y2 = sorted([m1, m2])
            if abs(y2 - y1) <= 0:
                continue
            dist_bp = abs((0.5 * (s2 + e2)) - (0.5 * (s1 + e1)))
            rad = 0.12 + 0.30 * (dist_bp / max_dist)

            arc = FancyArrowPatch(
                (0.02, y1),
                (0.02, y2),
                connectionstyle=f"arc3,rad={rad}",
                arrowstyle="-",
                linewidth=0.9,
                color=color,
                alpha=alpha,
                zorder=3,
            )
            ax.add_patch(arc)

    ax.set_xlabel(label, fontsize=7, labelpad=3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["left"].set_color("#aaaaaa")
    add_vertical_tick_bars(ax, region_start, region_end)
    ax.set_xticks([])


# ---------------------------------------------------------------------------
# Compatibility wrappers for extracted track renderers
# ---------------------------------------------------------------------------

def add_vertical_tick_bars(ax, region_start, region_end):
    return _add_vertical_tick_bars_impl(ax, region_start, region_end)


def draw_epilogos_horizontal(ax, positions, score_mat, num_states,
                             region_start, region_end,
                             categories, highlights=None, label="",
                             y_max=None):
    return _draw_epilogos_horizontal_impl(
        ax, positions, score_mat, num_states,
        region_start, region_end,
        categories, highlights=highlights, label=label,
        y_max=y_max,
    )


def draw_epilogos_vertical(ax, positions, score_mat, num_states,
                           region_start, region_end,
                           categories, highlights=None, label="",
                           y_max=None):
    return _draw_epilogos_vertical_impl(
        ax, positions, score_mat, num_states,
        region_start, region_end,
        categories, highlights=highlights, label=label,
        y_max=y_max,
    )


def draw_diff_horizontal(ax, positions, score_mat, num_states,
                         region_start, region_end,
                         categories, highlights=None, label="",
                         diff_max=None):
    return _draw_diff_horizontal_impl(
        ax, positions, score_mat, num_states,
        region_start, region_end,
        categories, highlights=highlights, label=label,
        diff_max=diff_max,
    )


def draw_diff_vertical(ax, positions, score_mat, num_states,
                       region_start, region_end,
                       categories, highlights=None, label="",
                       diff_max=None):
    return _draw_diff_vertical_impl(
        ax, positions, score_mat, num_states,
        region_start, region_end,
        categories, highlights=highlights, label=label,
        diff_max=diff_max,
    )


def draw_gene_track(ax, genes, region_start, region_end,
                    highlights=None, label="Genes"):
    return _draw_gene_track_impl(
        ax, genes, region_start, region_end,
        highlights=highlights, label=label,
    )


def draw_genomic_axis(ax, region_start, region_end, chrom):
    return _draw_genomic_axis_impl(ax, region_start, region_end, chrom)


def draw_gene_track_vertical(ax, genes, region_start, region_end,
                             highlights=None, label="Genes"):
    return _draw_gene_track_vertical_impl(
        ax, genes, region_start, region_end,
        highlights=highlights, label=label,
    )


def draw_genomic_axis_vertical(ax, region_start, region_end, chrom):
    return _draw_genomic_axis_vertical_impl(ax, region_start, region_end, chrom)


def draw_pval_horizontal(ax, positions, values, region_start, region_end,
                         highlights=None, label="", color="#a43cca",
                         y_max=None, cutoff_value=None):
    return _draw_pval_horizontal_impl(
        ax, positions, values, region_start, region_end,
        highlights=highlights, label=label, color=color,
        y_max=y_max, cutoff_value=cutoff_value,
    )


def draw_pval_vertical(ax, positions, values, region_start, region_end,
                       highlights=None, label="", color="#a43cca",
                       x_max=None, cutoff_value=None):
    return _draw_pval_vertical_impl(
        ax, positions, values, region_start, region_end,
        highlights=highlights, label=label, color=color,
        x_max=x_max, cutoff_value=cutoff_value,
    )


def draw_pval_diff_horizontal(ax, positions, values, region_start, region_end,
                              highlights=None, label="",
                              pos_color="#0c4a6e", neg_color="#d97706",
                              y_max=None, cutoff_value=None,
                              diff_score_positions=None, diff_score_matrix=None,
                              categories=None):
    return _draw_pval_diff_horizontal_impl(
        ax, positions, values, region_start, region_end,
        highlights=highlights, label=label,
        pos_color=pos_color, neg_color=neg_color,
        y_max=y_max, cutoff_value=cutoff_value,
        diff_score_positions=diff_score_positions,
        diff_score_matrix=diff_score_matrix,
        categories=categories,
    )


def draw_pval_diff_vertical(ax, positions, values, region_start, region_end,
                            highlights=None, label="",
                            pos_color="#0c4a6e", neg_color="#d97706",
                            x_max=None, cutoff_value=None,
                            diff_score_positions=None, diff_score_matrix=None,
                            categories=None):
    return _draw_pval_diff_vertical_impl(
        ax, positions, values, region_start, region_end,
        highlights=highlights, label=label,
        pos_color=pos_color, neg_color=neg_color,
        x_max=x_max, cutoff_value=cutoff_value,
        diff_score_positions=diff_score_positions,
        diff_score_matrix=diff_score_matrix,
        categories=categories,
    )


def draw_loops_horizontal(ax, loops, region_start, region_end,
                          highlights=None, label="", color="#4c78a8",
                          anchor_color="#4c78a8"):
    return _draw_loops_horizontal_impl(
        ax, loops, region_start, region_end,
        highlights=highlights, label=label, color=color,
        anchor_color=anchor_color,
    )


def draw_loops_vertical(ax, loops, region_start, region_end,
                        highlights=None, label="", color="#e05c3a",
                        anchor_color="#e05c3a"):
    return _draw_loops_vertical_impl(
        ax, loops, region_start, region_end,
        highlights=highlights, label=label, color=color,
        anchor_color=anchor_color,
    )


# ---------------------------------------------------------------------------
# Main figure assembly
# ---------------------------------------------------------------------------

def make_figure(
    hic_a_path, hic_b_path,
    qcat_a_path, qcat_b_path,
    region_str, resolution,
    out_path,
    loops_path=None,
    genes_path=None,
    gtf_path=None,
    highlights_path=None,
    diff_qcat_path=None,
    pval_a_path=None,
    pval_b_path=None,
    pval_diff_path=None,
    loops_a_path=None,
    loops_b_path=None,
    pval_cutoff=0.05,
    label_a="Condition A",
    label_b="Condition B",
    categories=None,
    hic_fmt=None,
    pval_overlay=False,
    rgb_hic=False,
    rgb_palette="magenta-green",
):
    chrom, region_start, region_end = parse_region(region_str)
    span = region_end - region_start
    print(f"\nRegion: {chrom}:{region_start:,}-{region_end:,}  ({span/1e6:.2f} Mb)")

    if categories is None:
        categories = ALL_CATEGORIES

    # --- Load data ---
    print("Loading Hi-C matrices...")
    mat_a, used_res_a = load_contact_matrix(
        hic_a_path, chrom, region_start, region_end,
        resolution, fmt=hic_fmt,
    )
    mat_b, used_res_b = load_contact_matrix(
        hic_b_path, chrom, region_start, region_end,
        resolution, fmt=hic_fmt,
    )
    combined = None if rgb_hic else make_split_hic(mat_a, mat_b)
    matrix_n = min(mat_a.shape[0], mat_b.shape[0])
    print(f"  Matrix size: {matrix_n} x {matrix_n} bins")
    if used_res_a == used_res_b:
        print(f"  Hi-C resolution used: {used_res_a:,} bp")
        resolution_note = f"Resolution: {used_res_a:,} bp"
    else:
        print(
            "  WARNING: Hi-C resolutions differ between A and B: "
            f"A={used_res_a:,} bp, B={used_res_b:,} bp"
        )
        resolution_note = (
            f"Resolution: A={used_res_a:,} bp, B={used_res_b:,} bp"
        )

    print("Loading Bearing epilogos scores...")
    pos_a, scores_a, ns_a = load_qcat_scores(
        qcat_a_path, chrom, region_start, region_end)
    pos_b, scores_b, ns_b = load_qcat_scores(
        qcat_b_path, chrom, region_start, region_end)
    num_states = max(ns_a, ns_b)
    # Pad score matrices to num_states if needed
    if scores_a.shape[1] < num_states:
        scores_a = np.pad(scores_a, ((0,0),(0,num_states-scores_a.shape[1])))
    if scores_b.shape[1] < num_states:
        scores_b = np.pad(scores_b, ((0,0),(0,num_states-scores_b.shape[1])))
    print(f"  Cond A: {len(pos_a)} bins, Cond B: {len(pos_b)} bins, "
          f"{num_states} states")

    loops = None
    if loops_path:
        print("Loading loops...")
        loops = load_loops(loops_path, chrom, region_start, region_end)
        print(f"  {len(loops)} loops in region")

    loops_a = None
    if loops_a_path:
        print("Loading condition A loops...")
        loops_a = load_loops(loops_a_path, chrom, region_start, region_end)
        print(f"  {len(loops_a)} loops in region")

    loops_b = None
    if loops_b_path:
        print("Loading condition B loops...")
        loops_b = load_loops(loops_b_path, chrom, region_start, region_end)
        print(f"  {len(loops_b)} loops in region")

    genes = None
    if genes_path:
        print("Loading genes...")
        genes = load_genes(genes_path, chrom, region_start, region_end)
        print(f"  {len(genes)} gene records in region")
    elif gtf_path:
        print("Loading genes from GTF...")
        genes = load_genes_gtf(gtf_path, chrom, region_start, region_end)
        print(f"  {len(genes)} gene records in region")

    highlights = None
    if highlights_path:
        print("Loading highlights...")
        highlights = load_highlights(highlights_path, chrom, region_start, region_end)
        print(f"  {len(highlights)} highlight regions")

    # --- Load differential qcat if provided ---
    pos_diff, scores_diff, ns_diff = [], np.zeros((0, 1), dtype=np.float32), 1
    has_diff = diff_qcat_path is not None
    if has_diff:
        print("Loading differential epilogos scores...")
        pos_diff, scores_diff, ns_diff = load_qcat_scores(
            diff_qcat_path, chrom, region_start, region_end)
        # Display convention in this figure is A-B across both diff panels.
        # Flip here so labels and rendered direction match user expectation.
        scores_diff = -scores_diff
        ns_diff_full = max(ns_diff, num_states)
        if scores_diff.shape[1] < ns_diff_full:
            scores_diff = np.pad(
                scores_diff,
                ((0, 0), (0, ns_diff_full - scores_diff.shape[1])),
            )
        print(f"  Diff: {len(pos_diff)} bins")

    # --- Load p-value BigWig tracks if provided ---
    pos_pval_a, vals_pval_a = None, None
    pos_pval_b, vals_pval_b = None, None
    pos_pval_diff, vals_pval_diff = None, None
    if pval_a_path is not None:
        print("Loading p-value track for condition A...")
        pos_pval_a, vals_pval_a = load_pval_track_values(
            pval_a_path, chrom, region_start, region_end)
        print(f"  {len(pos_pval_a)} bins")
    if pval_b_path is not None:
        print("Loading p-value track for condition B...")
        pos_pval_b, vals_pval_b = load_pval_track_values(
            pval_b_path, chrom, region_start, region_end)
        print(f"  {len(pos_pval_b)} bins")
    if pval_diff_path is not None:
        print("Loading signed differential p-value track...")
        pos_pval_diff, vals_pval_diff = load_pval_track_values(
            pval_diff_path, chrom, region_start, region_end)
        vals_pval_diff = -vals_pval_diff
        print(f"  {len(pos_pval_diff)} bins")

    pval_cutoff_value = None
    if pval_cutoff is not None and pval_cutoff > 0:
        pval_cutoff_value = -math.log10(pval_cutoff)

    # --- Figure layout ---
    #
    # 7-column x 7-row GridSpec:
    #
    #         col0(sq)   col1(epi_B)  col2(pval_B)  col3(loop_B)  col4(diff_B) col5(gene) col6(axis)
    #  row0   Hi-C       epi_B_v      pval_B_v      loop_B_v      diff_B_v     gene_v     axis_v
    #  row1   epi_A_h
    #  row2   pval_A_h
    #  row3   loop_A_h
    #  row4   diff_A_h
    #  row5   gene_h                                [  legend  spans col5:7, row5:7  ]
    #  row6   axis_h
    #
    # pval tracks only drawn when --pval-a/--pval-b are supplied.
    # diff tracks only drawn when has_diff=True.
    # All axes are always allocated so the layout is symmetric.

    has_genes  = genes is not None and len(genes) > 0
    has_pval_a = pos_pval_a is not None and len(pos_pval_a) > 0
    has_pval_b = pos_pval_b is not None and len(pos_pval_b) > 0
    has_pval_diff = pos_pval_diff is not None and len(pos_pval_diff) > 0
    has_loop_a = loops_a is not None and len(loops_a) > 0
    has_loop_b = loops_b is not None and len(loops_b) > 0

    use_pval_overlay = bool(pval_overlay and (has_pval_a or has_pval_b))
    draw_epi_with_pval_h = None
    draw_epi_with_pval_v = None
    if use_pval_overlay:
        try:
            from bearing_hic_plot_pval_overlay import (
                draw_epilogos_with_pval_horizontal as draw_epi_with_pval_h,
                draw_epilogos_with_pval_vertical as draw_epi_with_pval_v,
            )
            print("Using p-value overlay mode for epilogos tracks.")
        except ImportError as e:
            print(
                "WARNING: --pval-overlay was requested but "
                "bearing_hic_plot_pval_overlay.py is unavailable "
                f"({e}). Falling back to separate p-value panels."
            )
            use_pval_overlay = False

    fig, axes = create_main_figure(
        use_pval_overlay=use_pval_overlay,
        has_loops=has_loop_a or has_loop_b,
    )
    ax_hic = axes["hic"]
    ax_epi_a = axes["epi_a"]
    ax_pval_a = axes["pval_a"]
    ax_loop_a = axes["loop_a"]
    ax_diff_a = axes["diff_a"]
    ax_diff_pval_a = axes["diff_pval_a"]
    ax_gene_h = axes["gene_h"]
    ax_axis_h = axes["axis_h"]
    ax_epi_b = axes["epi_b"]
    ax_pval_b = axes["pval_b"]
    ax_loop_b = axes["loop_b"]
    ax_diff_b = axes["diff_b"]
    ax_diff_pval_b = axes["diff_pval_b"]
    ax_gene_v = axes["gene_v"]
    ax_axis_v = axes["axis_v"]
    ax_legend = axes["legend"]

    # --- Draw panels ---
    print("Rendering Hi-C panel...")
    hic_vmax = draw_hic_square(
        ax_hic, mat_a, mat_b, combined, resolution,
        region_start, region_end,
        loops=loops, highlights=highlights,
        label_a=label_a, label_b=label_b,
        rgb_mode=rgb_hic,
        rgb_palette=rgb_palette,
    )
    ax_hic.text(
        0.02, 0.98,
        resolution_note,
        transform=ax_hic.transAxes,
        fontsize=6,
        ha="left",
        va="top",
        color="#222222",
        zorder=7,
        bbox={"boxstyle": "round,pad=0.18", "facecolor": "white", "alpha": 0.7, "edgecolor": "none"},
    )

    # compute shared y-limit so each assay uses same vertical scale
    y_max_shared = 0.0
    if scores_a.size:
        y_max_shared = max(y_max_shared, float(scores_a.max()))
    if scores_b.size:
        y_max_shared = max(y_max_shared, float(scores_b.max()))
    if has_diff and scores_diff.size:
        # diff may have negative values; take absolute
        y_max_shared = max(y_max_shared, float(np.abs(scores_diff).max()))
    if y_max_shared == 0:
        y_max_shared = 1.0
    else:
        y_max_shared *= 1.05  # small headroom

    # p-value tracks should also share a common maximum among themselves
    y_max_shared_pval = 0.0
    if has_pval_a and vals_pval_a.size:
        y_max_shared_pval = max(y_max_shared_pval, float(vals_pval_a.max()))
    if has_pval_b and vals_pval_b.size:
        y_max_shared_pval = max(y_max_shared_pval, float(vals_pval_b.max()))
    if y_max_shared_pval == 0:
        y_max_shared_pval = 1.0
    else:
        y_max_shared_pval *= 1.05

    y_max_shared_pval_diff = 0.0
    if has_pval_diff and vals_pval_diff.size:
        y_max_shared_pval_diff = float(np.abs(vals_pval_diff).max())
    if y_max_shared_pval_diff == 0:
        y_max_shared_pval_diff = 1.0
    else:
        y_max_shared_pval_diff *= 1.05

    print("Rendering horizontal epilogos, pvalue and diff tracks...")
    if use_pval_overlay and has_pval_a:
        draw_epi_with_pval_h(
            ax_epi_a, pos_a, scores_a, num_states,
            region_start, region_end,
            categories[:num_states],
            pval_positions=pos_pval_a,
            pval_values=vals_pval_a,
            pval_alpha=pval_cutoff,
            highlights=highlights,
            label=label_a,
            y_max=y_max_shared,
        )
    else:
        draw_epilogos_horizontal(
            ax_epi_a, pos_a, scores_a, num_states,
            region_start, region_end,
            categories[:num_states],
            highlights=highlights,
            label=label_a,
            y_max=y_max_shared,
        )

    if use_pval_overlay:
        ax_pval_a.set_axis_off()
    elif has_pval_a:
        draw_pval_horizontal(
            ax_pval_a, pos_pval_a, vals_pval_a,
            region_start, region_end,
            highlights=highlights,
            label=label_a + " -log\u2081\u2080p",
            y_max=y_max_shared_pval,
            cutoff_value=pval_cutoff_value,
        )
    else:
        ax_pval_a.set_axis_off()

    if has_loop_a:
        draw_loops_horizontal(
            ax_loop_a, loops_a,
            region_start, region_end,
            highlights=highlights,
            label=label_a + " loops",
            color="#4c78a8",
            anchor_color="#4c78a8",
        )
    else:
        ax_loop_a.set_axis_off()

    if has_diff:
        draw_diff_horizontal(
            ax_diff_a, pos_diff, scores_diff, num_states,
            region_start, region_end,
            categories[:num_states],
            highlights=highlights,
            label=label_a + " - " + label_b,
        )
    else:
        ax_diff_a.set_axis_off()

    if has_pval_diff:
        draw_pval_diff_horizontal(
            ax_diff_pval_a, pos_pval_diff, vals_pval_diff,
            region_start, region_end,
            highlights=highlights,
            label=(label_a + " - " + label_b + " signed -log\u2081\u2080p"),
            y_max=y_max_shared_pval_diff,
            cutoff_value=pval_cutoff_value,
        )
    else:
        ax_diff_pval_a.set_axis_off()

    if has_genes:
        draw_gene_track(
            ax_gene_h, genes,
            region_start, region_end,
            highlights=highlights,
            label="Genes",
        )
    else:
        ax_gene_h.set_axis_off()

    draw_genomic_axis(ax_axis_h, region_start, region_end, chrom)

    print("Rendering vertical epilogos, pvalue and diff tracks...")
    if use_pval_overlay and has_pval_b:
        draw_epi_with_pval_v(
            ax_epi_b, pos_b, scores_b, num_states,
            region_start, region_end,
            categories[:num_states],
            pval_positions=pos_pval_b,
            pval_values=vals_pval_b,
            pval_alpha=pval_cutoff,
            highlights=highlights,
            label=label_b,
            y_max=y_max_shared,
        )
    else:
        draw_epilogos_vertical(
            ax_epi_b, pos_b, scores_b, num_states,
            region_start, region_end,
            categories[:num_states],
            highlights=highlights,
            label=label_b,
            y_max=y_max_shared,
        )

    if use_pval_overlay:
        ax_pval_b.set_axis_off()
    elif has_pval_b:
        draw_pval_vertical(
            ax_pval_b, pos_pval_b, vals_pval_b,
            region_start, region_end,
            highlights=highlights,
            label=label_b + " -log\u2081\u2080p",
            x_max=y_max_shared_pval,
            cutoff_value=pval_cutoff_value,
        )
    else:
        ax_pval_b.set_axis_off()

    if has_loop_b:
        draw_loops_vertical(
            ax_loop_b, loops_b,
            region_start, region_end,
            highlights=highlights,
            label=label_b + " loops",
            color="#e05c3a",
            anchor_color="#e05c3a",
        )
    else:
        ax_loop_b.set_axis_off()

    if has_diff:
        draw_diff_vertical(
            ax_diff_b, pos_diff, scores_diff, num_states,
            region_start, region_end,
            categories[:num_states],
            highlights=highlights,
            label=label_a + " - " + label_b,
        )
    else:
        ax_diff_b.set_axis_off()

    if has_pval_diff:
        draw_pval_diff_vertical(
            ax_diff_pval_b, pos_pval_diff, vals_pval_diff,
            region_start, region_end,
            highlights=highlights,
            label=(label_a + " - " + label_b + " signed -log\u2081\u2080p"),
            x_max=y_max_shared_pval_diff,
            cutoff_value=pval_cutoff_value,
        )
    else:
        ax_diff_pval_b.set_axis_off()

    if has_genes:
        draw_gene_track_vertical(
            ax_gene_v, genes,
            region_start, region_end,
            highlights=highlights,
            label="Genes",
        )
    else:
        ax_gene_v.set_axis_off()

    draw_genomic_axis_vertical(ax_axis_v, region_start, region_end, chrom)

    print("Rendering legend...")
    draw_legend(ax_legend, categories[:num_states], num_states,
                hic_vmax=hic_vmax,
                rgb_mode=rgb_hic,
                label_a=label_a,
                label_b=label_b,
                rgb_palette=rgb_palette)

    # --- Save ---
    print(f"Saving to {out_path}...")
    fig.savefig(str(out_path), dpi=150, bbox_inches="tight",
                facecolor="white")
    plt.close(fig)
    print(f"Done.  Figure: {out_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Composite Hi-C + Bearing epilogos comparison figure "
            "for two conditions at a genomic locus."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--contact-a", required=True, metavar="FILE",
                        help=("Hi-C contact file for condition A. "
                              "Accepts .cool or .mcool (default/preferred) "
                              "or .hic (Juicer format). "
                              "Format is auto-detected from extension."))
    parser.add_argument("--contact-b", required=True, metavar="FILE",
                        help="Hi-C contact file for condition B (same formats as --contact-a)")
    parser.add_argument("--format", dest="hic_fmt",
                        choices=["cool", "hic"], default=None, metavar="FMT",
                        help=("Force Hi-C format: 'cool' (.cool/.mcool) or 'hic' (.hic). "
                              "Default: auto-detect from file extension."))
    parser.add_argument("--qcat-a",  required=True, metavar="FILE",
                        help="Bearing qcat.bgz for condition A")
    parser.add_argument("--qcat-b",  required=True, metavar="FILE",
                        help="Bearing qcat.bgz for condition B")
    parser.add_argument("--diff-qcat", metavar="FILE", default=None,
                        help=(
                            "Differential qcat.bgz from compare_qcat.py "
                            "(A - B signed scores). When supplied, a diff "
                            "track is drawn below cond A epilogos and beside "
                            "cond B epilogos."
                        ))
    parser.add_argument("--pval-a", metavar="FILE", default=None,
                        help=(
                            "-log10(p-value) track for condition A. Accepts either "
                            "a BigWig (*.neglog10p.bw) or a stats TSV (*.stats.tsv) "
                            "with chrom/start/end/pval columns. Adds a horizontal "
                            "p-value track below the cond A epilogos."
                        ))
    parser.add_argument("--pval-b", metavar="FILE", default=None,
                        help=(
                            "-log10(p-value) track for condition B. Accepts either "
                            "a BigWig (*.neglog10p.bw) or a stats TSV (*.stats.tsv) "
                            "with chrom/start/end/pval columns. Adds a vertical "
                            "p-value track beside the cond B epilogos."
                        ))
    parser.add_argument("--pval-diff", metavar="FILE", default=None,
                        help=(
                            "Signed -log10(p-value) track from bearing_pvalue.py --diff. "
                            "Accepts either a BigWig (*.neglog10p.bw) or a stats TSV "
                            "(*.stats.tsv) with chrom/start/end/pval columns and an "
                            "optional direction column. When supplied, diff panels show "
                            "signed significance instead of differential epilogos."
                        ))
    parser.add_argument("--pval-cutoff", type=float, default=0.05,
                        metavar="P",
                        help=(
                            "P-value threshold drawn on p-value tracks as a dashed "
                            "dark gray line at -log10(P). Default: 0.05"
                        ))
    parser.add_argument("--pval-overlay", action="store_true",
                        help=(
                            "Overlay p-value significance directly on epilogos tracks "
                            "using bearing_hic_plot_pval_overlay.py. "
                            "Significant bins remain opaque; non-significant bins are "
                            "desaturated, with a compact jewel-purple -log10(p) twin axis."
                        ))
    parser.add_argument("--rgb-hic", action="store_true",
                        help=(
                            "Render Hi-C contacts as an RGB image: condition "
                            "A in red+blue (magenta), condition B in the green "
                            "channel. Co-enriched bins appear white; "
                            "condition-specific bins appear magenta or green. "
                            "Background is black. Replaces the split-diagonal "
                            "colourmap."
                        ))
    parser.add_argument(
        "--rgb-palette",
        choices=["magenta-green", "red-green", "blue-red", "green-blue", "magenta-green-white"],
        default="magenta-green",
        help=(
            "Color mapping for --rgb-hic. Default: magenta-green "
            "(A=magenta, B=green, overlap=white, background=black). "
            "Use magenta-green-white for white background with dark-purple overlap."
        ),
    )
    # Region specification (mutually exclusive)
    region_group = parser.add_mutually_exclusive_group(required=True)
    region_group.add_argument("--region", metavar="CHR:START-END",
                              help="Genomic region to plot, e.g. chr6:50000000-52000000")
    region_group.add_argument("--regions-file", metavar="TSV",
                              help="TSV file with columns: name, region [resolution, label, out]")
    
    parser.add_argument("--resolution", type=int, default=10000, metavar="BP",
                        help="Hi-C bin resolution in bp (default: 10000, overridden by regions file)")
    
    # Output specification
    parser.add_argument("--out", metavar="FILE",
                        help="Output PDF path (required for single region, ignored in batch mode)")
    parser.add_argument("--outdir", metavar="DIR", default=".",
                        help="Output directory for batch mode (default: current directory)")
    parser.add_argument("--label-a", default="Condition A", metavar="STR",
                        help="Label for condition A (default: 'Condition A')")
    parser.add_argument("--label-b", default="Condition B", metavar="STR",
                        help="Label for condition B (default: 'Condition B')")
    parser.add_argument("--loops",      metavar="BEDPE",
                        help="Loop anchors BEDPE file (optional)")
    parser.add_argument("--loops-a", metavar="BEDPE", default=None,
                        help=(
                            "Condition A loop BEDPE (e.g. Mustache output). "
                            "Drawn as a horizontal loop-arc track below "
                            "condition A p-value."
                        ))
    parser.add_argument("--loops-b", metavar="BEDPE", default=None,
                        help=(
                            "Condition B loop BEDPE (e.g. Mustache output). "
                            "Drawn as a vertical loop-arc track next to "
                            "condition B p-value."
                        ))
    parser.add_argument("--genes",      metavar="BED",
                        help="Gene annotation BED6 file (optional)")
    parser.add_argument("--gtf",      metavar="GTF",
                        help=(
                            "Gene annotation GTF file (optional). "
                            "Labels use gene_name when present, else gene_id."
                        ))
    parser.add_argument("--highlights", metavar="BED",
                        help="Region highlights BED4 file, col4=hex color (optional)")
    parser.add_argument(
        "--categories", metavar="YAML",
        help=(
            "YAML file defining category names and colors. "
            "Overrides the built-in 15-state mm10 definitions. "
            "Example files are in categories/."
        ),
    )

    args = parser.parse_args()

    if not (0 < args.pval_cutoff <= 1):
        parser.error("--pval-cutoff must be > 0 and <= 1")
    if args.genes and args.gtf:
        parser.error("Use only one of --genes or --gtf")

    # Load custom categories if provided
    cli_categories = None
    if args.categories:
        cli_categories, _ = load_categories_yaml(args.categories)
        print(f"Categories loaded from: {args.categories}  ({len(cli_categories)} states)")

    # Validate output arguments
    if args.region and not args.out:
        parser.error("--out is required when using --region")

    print(f"\n{'='*60}")
    print(f"  bearing_hic_plot.py")
    print(f"  Format:     {args.hic_fmt or 'auto-detect'}")
    print(f"  Cond A:     {args.label_a}  ({args.contact_a})")
    print(f"  Cond B:     {args.label_b}  ({args.contact_b})")
    if args.loops:      print(f"  Loops:      {args.loops}")
    if args.loops_a:    print(f"  Loops A:    {args.loops_a}")
    if args.loops_b:    print(f"  Loops B:    {args.loops_b}")
    if args.genes:      print(f"  Genes:      {args.genes}")
    if args.gtf:        print(f"  GTF:        {args.gtf}")
    if args.highlights: print(f"  Highlights: {args.highlights}")
    if args.diff_qcat:  print(f"  Diff qcat:  {args.diff_qcat}")
    if args.pval_a:     print(f"  P-val A:    {args.pval_a}")
    if args.pval_b:     print(f"  P-val B:    {args.pval_b}")
    if args.pval_diff:  print(f"  P-val diff: {args.pval_diff}")
    if args.pval_a or args.pval_b or args.pval_diff:
        print(f"  P cutoff:   {args.pval_cutoff} (-log10={-math.log10(args.pval_cutoff):.3f})")
    if args.pval_overlay:
        print("  P overlay:  enabled")
    if args.rgb_hic:
        print("  Hi-C mode:  RGB")
        print(f"  RGB palette: {args.rgb_palette}")
    print(f"{'='*60}\n")

    # Batch mode: process multiple regions from file
    if args.regions_file:
        regions = load_regions_file(args.regions_file)
        outdir = Path(args.outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        
        print(f"Batch mode: processing {len(regions)} region(s) from {args.regions_file}")
        print(f"Output directory: {outdir}\n")
        
        for i, reg in enumerate(regions, start=1):
            name = reg["name"]
            region_str = reg["region"]
            resolution = reg["resolution"] if reg["resolution"] else args.resolution
            label = reg["label"] if reg["label"] else name

            hic_a_path, hic_b_path, resolve_note = _resolve_contacts_for_region_resolution(
                args.contact_a,
                args.contact_b,
                resolution,
                hic_fmt=args.hic_fmt,
            )
            
            # Determine output path
            if reg["out"]:
                out_path = outdir / reg["out"]
            else:
                # Sanitize name for filename
                safe_name = "".join(
                    ch if ch.isalnum() or ch in {"-", "_"} else "_" 
                    for ch in name
                )
                out_path = outdir / f"{safe_name}_hic.pdf"
            
            print(f"[{i}/{len(regions)}] {name}: {region_str} @ {resolution:,} bp -> {out_path}")
            print(f"  Hi-C files: {hic_a_path} | {hic_b_path}")
            print(f"  Resolution match: {resolve_note}")
            
            try:
                make_figure(
                    hic_a_path      = hic_a_path,
                    hic_b_path      = hic_b_path,
                    qcat_a_path     = args.qcat_a,
                    qcat_b_path     = args.qcat_b,
                    region_str      = region_str,
                    resolution      = resolution,
                    out_path        = out_path,
                    loops_path      = args.loops,
                    loops_a_path    = args.loops_a,
                    loops_b_path    = args.loops_b,
                    genes_path      = args.genes,
                    gtf_path        = args.gtf,
                    highlights_path = args.highlights,
                    diff_qcat_path  = args.diff_qcat,
                    pval_a_path     = args.pval_a,
                    pval_b_path     = args.pval_b,
                    pval_diff_path  = args.pval_diff,
                    pval_cutoff     = args.pval_cutoff,
                    pval_overlay    = args.pval_overlay,
                    label_a         = args.label_a,
                    label_b         = args.label_b,
                    hic_fmt         = args.hic_fmt,
                    categories      = cli_categories,
                    rgb_hic         = args.rgb_hic,
                    rgb_palette     = args.rgb_palette,
                )
            except Exception as e:
                print(f"  ERROR: Failed to process region {name}: {e}")
                continue
        
        print(f"\nBatch complete. {len(regions)} figure(s) written to {outdir}/")
    
    # Single region mode
    else:
        print(f"Single region mode:")
        print(f"  Region:     {args.region}")
        print(f"  Resolution: {args.resolution:,} bp")
        print(f"  Output:     {args.out}\n")
        
        make_figure(
            hic_a_path      = args.contact_a,
            hic_b_path      = args.contact_b,
            qcat_a_path     = args.qcat_a,
            qcat_b_path     = args.qcat_b,
            region_str      = args.region,
            resolution      = args.resolution,
            out_path        = args.out,
            loops_path      = args.loops,
            loops_a_path    = args.loops_a,
            loops_b_path    = args.loops_b,
            genes_path      = args.genes,
            gtf_path        = args.gtf,
            highlights_path = args.highlights,
            diff_qcat_path  = args.diff_qcat,
            pval_a_path     = args.pval_a,
            pval_b_path     = args.pval_b,
            pval_diff_path  = args.pval_diff,
            pval_cutoff     = args.pval_cutoff,
            pval_overlay    = args.pval_overlay,
            label_a         = args.label_a,
            label_b         = args.label_b,
            hic_fmt         = args.hic_fmt,
            categories      = cli_categories,
            rgb_hic         = args.rgb_hic,
            rgb_palette     = args.rgb_palette,
        )


if __name__ == "__main__":
    main()
