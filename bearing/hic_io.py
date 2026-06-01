from __future__ import annotations

import sys

import numpy as np

from .plot_loaders import detect_hic_format


def load_hic_matrix(hic_path, chrom, start, end, resolution):
    """Load a Hi-C contact matrix from a .hic file using hicstraw."""
    try:
        import hicstraw
    except ImportError:
        sys.exit(
            "ERROR: hicstraw is required for .hic files.\n"
            "Install with:  pip install hicstraw"
        )

    hic = hicstraw.HiCFile(str(hic_path))
    available = hic.getResolutions()
    if resolution not in available:
        closest = min(available, key=lambda r: abs(r - resolution))
        print(f"  WARNING: resolution {resolution} not in .hic file. Using {closest} instead.")
        resolution = closest

    norm_order = ["KR", "SCALE", "VC_SQRT", "NONE"]
    result = None
    used_norm = None
    chrom_str = chrom.replace("chr", "")
    for norm in norm_order:
        try:
            candidate = hic.getMatrixZoomData(
                chrom_str, chrom_str, "observed", norm, "BP", resolution,
            )
            candidate.getRecords(start, start + resolution, start, start + resolution)
            result = candidate
            used_norm = norm
            break
        except Exception:
            continue

    if result is None:
        sys.exit(
            f"ERROR: Could not load Hi-C matrix for {chrom} at {resolution} bp with any normalization. "
            f"Check that this resolution exists in the .hic file."
        )

    if used_norm != "KR":
        print(f"  NOTE: KR normalization not available for {chrom} at {resolution} bp -- using {used_norm} instead.")

    records = result.getRecords(start, end, start, end)
    n_bins = (end - start) // resolution + 1
    mat = np.zeros((n_bins, n_bins), dtype=np.float64)

    for rec in records:
        i = (rec.binX - start) // resolution
        j = (rec.binY - start) // resolution
        if 0 <= i < n_bins and 0 <= j < n_bins:
            mat[i, j] = rec.counts
            mat[j, i] = rec.counts

    mat = np.where(np.isnan(mat), 0.0, mat)
    mat = np.log1p(mat)
    return mat, int(resolution)


def load_cool_matrix(cool_path, chrom, start, end, resolution):
    """Load a Hi-C contact matrix from a .cool/.mcool file using cooler."""
    try:
        import cooler
    except ImportError:
        sys.exit(
            "ERROR: cooler is required for .cool/.mcool files.\n"
            "Install with:  pip install cooler"
        )

    cool_path = str(cool_path)
    used_resolution = int(resolution)
    if cool_path.endswith(".mcool"):
        try:
            clr = cooler.Cooler(cool_path + "::/resolutions/" + str(resolution))
        except Exception:
            try:
                import h5py
                with h5py.File(cool_path, "r") as f:
                    available = sorted(int(r) for r in f.get("resolutions", {}).keys())
            except Exception:
                available = []
            if available:
                closest = min(available, key=lambda r: abs(r - resolution))
                print(f"  WARNING: resolution {resolution} not in .mcool. Using {closest} instead.")
                used_resolution = int(closest)
                clr = cooler.Cooler(cool_path + "::/resolutions/" + str(closest))
            else:
                sys.exit(
                    f"ERROR: Could not open {cool_path} at resolution {resolution}. "
                    f"Run 'cooler ls {cool_path}' to check available resolutions."
                )
    else:
        clr = cooler.Cooler(cool_path)
        if getattr(clr, "binsize", None) is not None:
            used_resolution = int(clr.binsize)

    region_str = f"{chrom}:{start}-{end}"
    try:
        mat = clr.matrix(balance=True).fetch(region_str)
    except Exception:
        print("  NOTE: ICE balancing weights not found, using raw counts.")
        mat = clr.matrix(balance=False).fetch(region_str)

    mat = np.where(np.isnan(mat), 0.0, mat).astype(np.float64)
    mat = np.log1p(mat)
    return mat, used_resolution


def load_contact_matrix(path, chrom, start, end, resolution, fmt=None):
    """Unified loader across cool/hic formats."""
    if fmt is None:
        fmt = detect_hic_format(path)
    if fmt == "cool":
        return load_cool_matrix(path, chrom, start, end, resolution)
    return load_hic_matrix(path, chrom, start, end, resolution)
