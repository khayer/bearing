from __future__ import annotations

import re
from pathlib import Path


def detect_hic_format(path):
    """Return 'cool' for .cool/.mcool, 'hic' for .hic (case-insensitive); else ValueError."""
    s = str(path).lower()
    if s.endswith(".cool") or s.endswith(".mcool"):
        return "cool"
    if s.endswith(".hic"):
        return "hic"
    raise ValueError(
        f"Cannot detect Hi-C format from filename: {path}\n"
        f"Expected .cool, .mcool, or .hic"
    )


# Filename tokens that encode a .cool file's bin size, e.g. ``sample_bs_5000.cool``.
_COOL_RES_TOKEN = re.compile(
    r"(_(?:bs|res|resolution|bin|binsize)_)(\d+)(\.cool)$", flags=re.IGNORECASE
)


def cool_resolution_variant(path, resolution):
    """Return ``path`` with its filename resolution token swapped for ``resolution``.

    Returns None when the filename carries no recognised resolution token.
    """
    p = Path(path)
    m = _COOL_RES_TOKEN.search(p.name)
    if not m:
        return None
    new_name = p.name[: m.start()] + m.group(1) + str(int(resolution)) + m.group(3)
    return p.with_name(new_name)


def resolve_contacts_for_region_resolution(contact_a, contact_b, resolution, hic_fmt=None):
    """Swap a pair of .cool paths for their ``resolution``-specific siblings.

    Returns ``(path_a, path_b, message)``. The substitution only happens when
    both files are plain .cool files with a resolution token and both
    resolution-matched siblings exist; otherwise the originals are returned
    with a message explaining why.
    """
    path_a = Path(contact_a)
    path_b = Path(contact_b)
    skipped = "auto-resolution file matching skipped"

    fmt_a = hic_fmt or detect_hic_format(path_a)
    fmt_b = hic_fmt or detect_hic_format(path_b)
    if fmt_a != "cool" or fmt_b != "cool":
        return path_a, path_b, f"{skipped} (non-cool format)"

    names = (str(path_a).lower(), str(path_b).lower())
    if any(n.endswith(".mcool") for n in names):
        return path_a, path_b, f"{skipped} (.mcool uses internal resolutions)"
    if not all(n.endswith(".cool") for n in names):
        return path_a, path_b, f"{skipped} (unsupported extension)"

    cand_a = cool_resolution_variant(path_a, resolution)
    cand_b = cool_resolution_variant(path_b, resolution)
    if cand_a is None or cand_b is None:
        return path_a, path_b, f"{skipped} (no resolution token in filename)"

    has_a = cand_a.exists()
    has_b = cand_b.exists()
    if has_a and has_b:
        return cand_a, cand_b, "using resolution-matched .cool files"
    if has_a or has_b:
        return path_a, path_b, "found only one matching .cool file; keeping original pair"
    return path_a, path_b, "no matching .cool files found; keeping original pair"


def _parse_item_rgb(col9):
    """Parse BED itemRgb as 'R,G,B' or a packed 0xRRGGBB integer; None if invalid."""
    col9 = col9.strip()
    try:
        if "," in col9:
            rgb = tuple(int(x) for x in col9.split(","))
            return rgb if len(rgb) == 3 else None
        n = int(col9)
    except ValueError:
        return None
    return ((n >> 16) & 0xFF, (n >> 8) & 0xFF, n & 0xFF)


def _is_bed_header(line):
    return not line or line.startswith("#") or line.lower().startswith("track")


def load_bed_for_region(bed_path, chrom, region_start, region_end):
    """
    Load BED (BED4 to BED9) and return list of feature dicts overlapping region.

    Returns list of dicts with keys: chrom, start, end, name, score, strand,
    item_rgb (tuple or None), bed_columns (int)

    Missing files yield an empty list; header, malformed and off-chromosome
    lines are skipped. Overlap is inclusive of the region endpoints.
    """
    p = Path(bed_path)
    features = []
    if not p.exists():
        return features

    with open(p, "r") as fh:
        for line in fh:
            line = line.strip()
            if _is_bed_header(line):
                continue
            parts = line.split("\t")
            if len(parts) < 3 or parts[0] != chrom:
                continue
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue
            if end < region_start or start > region_end:
                continue

            n_cols = len(parts)
            features.append({
                "chrom": parts[0],
                "start": start,
                "end": end,
                "name": parts[3] if n_cols >= 4 else "",
                "score": parts[4] if n_cols >= 5 else "",
                "strand": parts[5] if n_cols >= 6 else None,
                "item_rgb": _parse_item_rgb(parts[8]) if n_cols >= 9 else None,
                "bed_columns": n_cols,
            })
    return features
