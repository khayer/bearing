from __future__ import annotations

import re
from pathlib import Path


def detect_hic_format(path):
    s = str(path).lower()
    if s.endswith(".cool") or s.endswith(".mcool"):
        return "cool"
    if s.endswith(".hic"):
        return "hic"
    raise ValueError(
        f"Cannot detect Hi-C format from filename: {path}\n"
        f"Expected .cool, .mcool, or .hic"
    )


def cool_resolution_variant(path, resolution):
    p = Path(path)
    s = p.name
    res_s = str(int(resolution))
    patterns = [
        r"(_bs_)(\d+)(\.cool)$",
        r"(_res_)(\d+)(\.cool)$",
        r"(_resolution_)(\d+)(\.cool)$",
        r"(_bin_)(\d+)(\.cool)$",
        r"(_binsize_)(\d+)(\.cool)$",
    ]
    for pat in patterns:
        m = re.search(pat, s, flags=re.IGNORECASE)
        if not m:
            continue
        repl = m.group(1) + res_s + m.group(3)
        new_name = re.sub(pat, repl, s, flags=re.IGNORECASE)
        return p.with_name(new_name)
    return None


def resolve_contacts_for_region_resolution(contact_a, contact_b, resolution, hic_fmt=None):
    path_a = Path(contact_a)
    path_b = Path(contact_b)

    fmt_a = hic_fmt or detect_hic_format(path_a)
    fmt_b = hic_fmt or detect_hic_format(path_b)
    if fmt_a != "cool" or fmt_b != "cool":
        return path_a, path_b, "auto-resolution file matching skipped (non-cool format)"

    if str(path_a).lower().endswith(".mcool") or str(path_b).lower().endswith(".mcool"):
        return path_a, path_b, "auto-resolution file matching skipped (.mcool uses internal resolutions)"
    if not str(path_a).lower().endswith(".cool") or not str(path_b).lower().endswith(".cool"):
        return path_a, path_b, "auto-resolution file matching skipped (unsupported extension)"

    cand_a = cool_resolution_variant(path_a, resolution)
    cand_b = cool_resolution_variant(path_b, resolution)
    if cand_a is None or cand_b is None:
        return path_a, path_b, "auto-resolution file matching skipped (no resolution token in filename)"

    has_a = cand_a.exists()
    has_b = cand_b.exists()
    if has_a and has_b:
        return cand_a, cand_b, "using resolution-matched .cool files"
    if has_a or has_b:
        return path_a, path_b, "found only one matching .cool file; keeping original pair"
    return path_a, path_b, "no matching .cool files found; keeping original pair"


def load_bed_for_region(bed_path, chrom, region_start, region_end):
    """
    Load BED (BED4 to BED9) and return list of feature dicts overlapping region.

    Returns list of dicts with keys: chrom, start, end, name, score, strand,
    item_rgb (tuple or None), bed_columns (int)
    """
    from pathlib import Path

    p = Path(bed_path)
    features = []
    if not p.exists():
        return features

    with open(p, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.lower().startswith("track"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            chrom_f = parts[0]
            try:
                start = int(parts[1])
                end = int(parts[2])
            except Exception:
                continue
            if chrom_f != chrom:
                continue
            if end < region_start or start > region_end:
                continue

            name = parts[3] if len(parts) >= 4 else ""
            score = parts[4] if len(parts) >= 5 else ""
            strand = parts[5] if len(parts) >= 6 else None
            item_rgb = None
            if len(parts) >= 9:
                col9 = parts[8].strip()
                # parse comma sep R,G,B
                try:
                    if "," in col9:
                        rgb = tuple(int(x) for x in col9.split(","))
                        if len(rgb) == 3:
                            item_rgb = rgb
                    else:
                        # integer rgb
                        n = int(col9)
                        r = (n >> 16) & 0xFF
                        g = (n >> 8) & 0xFF
                        b = n & 0xFF
                        item_rgb = (r, g, b)
                except Exception:
                    item_rgb = None

            features.append({
                "chrom": chrom_f,
                "start": start,
                "end": end,
                "name": name,
                "score": score,
                "strand": strand,
                "item_rgb": item_rgb,
                "bed_columns": len(parts),
            })
    return features
