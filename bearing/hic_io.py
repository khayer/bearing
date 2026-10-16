from __future__ import annotations

import importlib
import sys

import numpy as np

from .plot_loaders import detect_hic_format

# Normalizations tried for .hic files, in order of preference.
HIC_NORM_ORDER = ("KR", "SCALE", "VC_SQRT", "NONE")


def _import_or_exit(module_name, file_kind):
    """Import an optional dependency or exit with an install hint."""
    try:
        return importlib.import_module(module_name)
    except ImportError:
        sys.exit(
            f"ERROR: {module_name} is required for {file_kind} files.\n"
            f"Install with:  pip install {module_name}"
        )


def _closest_resolution(available, resolution):
    return min(available, key=lambda r: abs(r - resolution))


def _log_contacts(mat):
    """Replace NaNs with 0 and return log1p(counts) as float64."""
    mat = np.where(np.isnan(mat), 0.0, mat).astype(np.float64)
    return np.log1p(mat)


def _hic_zoom_data(hic, chrom, start, resolution):
    """Return (zoom_data, norm) for the first normalization that yields records.

    Each candidate is probed with a single-bin query at ``start`` because
    missing normalization vectors only raise once records are requested.
    Returns (None, None) when no normalization works.
    """
    chrom_str = chrom.replace("chr", "")
    for norm in HIC_NORM_ORDER:
        try:
            candidate = hic.getMatrixZoomData(
                chrom_str, chrom_str, "observed", norm, "BP", resolution,
            )
            candidate.getRecords(start, start + resolution, start, start + resolution)
            return candidate, norm
        except Exception:
            continue
    return None, None


def load_hic_matrix(hic_path, chrom, start, end, resolution):
    """Load a Hi-C contact matrix from a .hic file using hicstraw."""
    hicstraw = _import_or_exit("hicstraw", ".hic")

    hic = hicstraw.HiCFile(str(hic_path))
    available = hic.getResolutions()
    if resolution not in available:
        closest = _closest_resolution(available, resolution)
        print(f"  WARNING: resolution {resolution} not in .hic file. Using {closest} instead.")
        resolution = closest

    zoom_data, used_norm = _hic_zoom_data(hic, chrom, start, resolution)
    if zoom_data is None:
        sys.exit(
            f"ERROR: Could not load Hi-C matrix for {chrom} at {resolution} bp with any normalization. "
            f"Check that this resolution exists in the .hic file."
        )
    if used_norm != "KR":
        print(f"  NOTE: KR normalization not available for {chrom} at {resolution} bp -- using {used_norm} instead.")

    n_bins = (end - start) // resolution + 1
    mat = np.zeros((n_bins, n_bins), dtype=np.float64)
    for rec in zoom_data.getRecords(start, end, start, end):
        i = (rec.binX - start) // resolution
        j = (rec.binY - start) // resolution
        if 0 <= i < n_bins and 0 <= j < n_bins:
            mat[i, j] = mat[j, i] = rec.counts

    return _log_contacts(mat), int(resolution)


def _mcool_resolutions(cool_path):
    """List resolutions stored in an .mcool, or [] if it cannot be read."""
    try:
        import h5py
        with h5py.File(cool_path, "r") as f:
            return sorted(int(r) for r in f.get("resolutions", {}).keys())
    except Exception:
        return []


def _open_mcool(cooler, cool_path, resolution):
    """Open an .mcool at ``resolution``, falling back to the closest stored one.

    Returns (Cooler, used_resolution); exits when nothing can be opened.
    """
    try:
        return cooler.Cooler(f"{cool_path}::/resolutions/{resolution}"), int(resolution)
    except Exception:
        pass
    available = _mcool_resolutions(cool_path)
    if not available:
        sys.exit(
            f"ERROR: Could not open {cool_path} at resolution {resolution}. "
            f"Run 'cooler ls {cool_path}' to check available resolutions."
        )
    closest = _closest_resolution(available, resolution)
    print(f"  WARNING: resolution {resolution} not in .mcool. Using {closest} instead.")
    return cooler.Cooler(f"{cool_path}::/resolutions/{closest}"), int(closest)


def load_cool_matrix(cool_path, chrom, start, end, resolution):
    """Load a Hi-C contact matrix from a .cool/.mcool file using cooler."""
    cooler = _import_or_exit("cooler", ".cool/.mcool")

    cool_path = str(cool_path)
    if cool_path.endswith(".mcool"):
        clr, used_resolution = _open_mcool(cooler, cool_path, resolution)
    else:
        clr = cooler.Cooler(cool_path)
        binsize = getattr(clr, "binsize", None)
        used_resolution = int(binsize) if binsize is not None else int(resolution)

    region_str = f"{chrom}:{start}-{end}"
    try:
        mat = clr.matrix(balance=True).fetch(region_str)
    except Exception:
        print("  NOTE: ICE balancing weights not found, using raw counts.")
        mat = clr.matrix(balance=False).fetch(region_str)

    return _log_contacts(mat), used_resolution


def load_contact_matrix(path, chrom, start, end, resolution, fmt=None):
    """Unified loader across cool/hic formats."""
    if fmt is None:
        fmt = detect_hic_format(path)
    if fmt == "cool":
        return load_cool_matrix(path, chrom, start, end, resolution)
    return load_hic_matrix(path, chrom, start, end, resolution)
