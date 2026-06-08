#!/usr/bin/env python3
"""
Batch runner for bearing_hic_plot.py:
- Compare each condition to one user-specified reference condition.
- Reuse a regions template TSV via --regions-file.
- Auto-wire qcat + pvalue tracks from a sample sheet.

Expected sample sheet columns:
  sample, condition, replicate, out

P-value track convention (from bearing_pvalue.py):
    <results-dir>/<sample>.stats.tsv preferred, falling back to
    <results-dir>/<sample>.neglog10p.bw
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
import gzip
from pathlib import Path
from typing import Dict, List, Tuple

from bearing.utils import (
    iter_non_comment_lines,
    parse_key_value_items,
    parse_ucsc_region,
    resolve_path,
    sanitize_token,
)
from bearing.tsv import read_tsv_dict_rows
from bearing.tsv import read_tsv_table
from bearing.diagnostics import require_files, require_mapping_paths
from bearing.runner import execute_command


def _parse_key_value(items: List[str], arg_name: str) -> Dict[str, str]:
    return parse_key_value_items(items, arg_name)


def _safe(s: str) -> str:
    return sanitize_token(s)


def _resolve(path_value: str, base_dir: Path) -> Path:
    return resolve_path(path_value, base_dir)


def _non_comment_lines(handle):
    yield from iter_non_comment_lines(handle)


def _parse_region(region_str: str) -> Tuple[str, int, int]:
    return parse_ucsc_region(region_str)


def _load_regions(path: Path) -> List[dict]:
    regions: List[dict] = []
    try:
        rows = read_tsv_dict_rows(path, required_columns=["name", "region"])
    except ValueError as exc:
        raise ValueError(str(exc).replace("TSV", "Regions file")) from exc

    for row in rows:
        name = row.get("name", "")
        region = row.get("region", "")
        if not name or not region:
            continue
        regions.append({"name": name, "region": region})
    return regions


def _strip_qcat_suffix(path: Path) -> str:
    name = path.name
    if name.endswith(".qcat.bgz"):
        return name[:-len(".qcat.bgz")]
    if name.endswith(".bgz"):
        return name[:-len(".bgz")]
    return path.stem


def _orient_diff_qcat(src: Path, dst: Path) -> Path:
    """Rewrite a signed diff qcat.bgz file so positive values mean ref - cond."""
    import json
    import os

    try:
        import pysam
    except ImportError as exc:
        raise RuntimeError("pysam is required to orient diff qcat files") from exc

    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp_tsv = str(dst) + ".tmp.tsv"
    with gzip.open(src, "rt") as fin, open(tmp_tsv, "w") as fout:
        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            meta = parts[3]
            if meta.startswith("{"):
                payload = json.loads(meta)
                pairs = payload.get("qcat", [])
            else:
                qcat_start = meta.find("qcat:")
                if qcat_start < 0:
                    fout.write(line + "\n")
                    continue
                raw_start = meta.find(",raw:", qcat_start)
                if raw_start >= 0:
                    prefix = meta[:qcat_start + 5]
                    qcat_payload = meta[qcat_start + 5:raw_start]
                    raw_suffix = meta[raw_start:]
                else:
                    prefix = meta[:qcat_start + 5]
                    qcat_payload = meta[qcat_start + 5:]
                    raw_suffix = ""

                pairs = json.loads(qcat_payload)
            flipped = [[float(f"{-float(score):.6g}"), int(state_idx)] for score, state_idx in pairs]
            if meta.startswith("{"):
                payload["qcat"] = flipped
                parts[3] = json.dumps(payload, separators=(",", ":"))
            else:
                parts[3] = prefix + json.dumps(flipped, separators=(",", ":")) + raw_suffix
            fout.write("\t".join(parts) + "\n")

    pysam.tabix_compress(tmp_tsv, str(dst), force=True)
    pysam.tabix_index(str(dst), preset="bed", force=True)
    os.remove(tmp_tsv)
    return dst


def _orient_diff_stats_tsv(src: Path, dst: Path) -> Path:
    """Rewrite a signed diff stats TSV so it reports the forward ref - cond direction."""
    dst.parent.mkdir(parents=True, exist_ok=True)
    with open(src, newline="") as fin, open(dst, "w", newline="") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{src} has no header")
        fieldnames = list(reader.fieldnames)
        lower = {name.lower(): name for name in fieldnames if name is not None}
        score_col = lower.get("bearing_score")
        tested_col = lower.get("bearing_score_tested")
        direction_col = lower.get("direction")

        writer = csv.DictWriter(fout, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in reader:
            if score_col and row.get(score_col):
                try:
                    row[score_col] = f"{-float(row[score_col]):.6f}"
                except ValueError:
                    pass
            if tested_col and row.get(tested_col):
                try:
                    row[tested_col] = f"{abs(float(row[tested_col])):.6f}"
                except ValueError:
                    pass
            if direction_col and row.get(direction_col):
                d = row[direction_col].strip()
                if d == "+":
                    row[direction_col] = "-"
                elif d == "-":
                    row[direction_col] = "+"
            writer.writerow(row)
    return dst


def _orient_diff_bigwig(src: Path, dst: Path) -> Path:
    """Rewrite a signed diff BigWig so positive values mean ref - cond."""
    try:
        import pyBigWig
    except ImportError as exc:
        raise RuntimeError("pyBigWig is required to orient diff BigWig files") from exc

    dst.parent.mkdir(parents=True, exist_ok=True)
    with pyBigWig.open(str(src)) as bw_in, pyBigWig.open(str(dst), "w") as bw_out:
        chroms = bw_in.chroms()
        bw_out.addHeader([(chrom, size) for chrom, size in chroms.items()])
        for chrom, size in chroms.items():
            intervals = bw_in.intervals(chrom)
            if not intervals:
                continue
            starts, ends, values = [], [], []
            for start, end, value in intervals:
                starts.append(int(start))
                ends.append(int(end))
                values.append(float(-value))
            bw_out.addEntries([chrom] * len(starts), starts, ends=ends, values=values)
    return dst


def _prepare_forward_oriented_diff_assets(cond: str, ref: str, comp_dir: Path, results_dir: Path, outdir: Path):
    """Return forward-oriented diff qcat and p-value assets for ref vs cond."""
    fwd = comp_dir / f"diff_{ref}_vs_{cond}.qcat.bgz"
    rev = comp_dir / f"diff_{cond}_vs_{ref}.qcat.bgz"
    oriented_dir = outdir / "oriented_diff_assets"
    oriented_dir.mkdir(parents=True, exist_ok=True)

    diff_qcat = None
    diff_pval = None
    source_base = None

    if fwd.exists():
        diff_qcat = fwd
        source_base = _strip_qcat_suffix(fwd)
    elif rev.exists():
        source_base = _strip_qcat_suffix(rev)
        diff_qcat = _orient_diff_qcat(rev, oriented_dir / f"diff_{ref}_vs_{cond}.qcat.bgz")
        print(
            f"[INFO] {cond}: oriented {rev.name} -> {diff_qcat.name} so diff tracks display {ref} - {cond}",
            file=sys.stderr,
        )

    if diff_qcat is not None:
        diff_base = _strip_qcat_suffix(diff_qcat)
        lookup_base = source_base or diff_base
        cand_diff_stats = results_dir / f"{lookup_base}.stats.tsv"
        cand_diff_pval = results_dir / f"{lookup_base}.neglog10p.bw"
        if cand_diff_stats.exists():
            if diff_qcat.parent == oriented_dir:
                diff_pval = _orient_diff_stats_tsv(
                    cand_diff_stats,
                    oriented_dir / f"{diff_base}.stats.tsv",
                )
                print(
                    f"[INFO] {cond}: oriented {cand_diff_stats.name} for {ref} - {cond}",
                    file=sys.stderr,
                )
            else:
                diff_pval = cand_diff_stats
        elif cand_diff_pval.exists():
            if diff_qcat.parent == oriented_dir:
                diff_pval = _orient_diff_bigwig(
                    cand_diff_pval,
                    oriented_dir / f"{diff_base}.neglog10p.bw",
                )
                print(
                    f"[INFO] {cond}: oriented {cand_diff_pval.name} for {ref} - {cond}",
                    file=sys.stderr,
                )
            else:
                diff_pval = cand_diff_pval
        else:
            print(
                f"[WARN] {cond}: diff qcat detected but no diff p-value track found at {cand_diff_stats} or {cand_diff_pval}",
                file=sys.stderr,
            )

    return diff_qcat, diff_pval


def _count_significant_bins_in_region(bw_path: Path, chrom: str, start: int, end: int,
                                      cutoff_neglog10: float) -> int:
    try:
        import pyBigWig
    except ImportError as exc:
        raise RuntimeError("pyBigWig is required for summary counting.") from exc

    count = 0
    with pyBigWig.open(str(bw_path)) as bw:
        chroms = bw.chroms()
        fetch_chrom = chrom
        if fetch_chrom not in chroms:
            alt = fetch_chrom[3:] if fetch_chrom.startswith("chr") else ("chr" + fetch_chrom)
            if alt in chroms:
                fetch_chrom = alt
            else:
                return 0

        intervals = bw.intervals(fetch_chrom, start, end)
        if not intervals:
            return 0
        for _, _, value in intervals:
            if abs(float(value)) >= cutoff_neglog10:
                count += 1
    return count


def _load_rows(sheet_path: Path) -> List[dict]:
    fields, raw_rows = read_tsv_table(sheet_path)
    required = ["condition", "out"]
    missing = [h for h in required if h not in fields]
    if missing:
        raise ValueError(
            f"Sample sheet missing required columns: {', '.join(missing)}"
        )
    if "sample" not in fields and "name" not in fields:
        raise ValueError("Sample sheet missing required columns: sample or name")

    rows: List[dict] = []
    for raw in raw_rows:
        row = {
            "sample": raw.get("sample", "") or raw.get("name", ""),
            "condition": raw.get("condition", ""),
            "replicate": raw.get("replicate", "0") or "0",
            "out": raw.get("out", ""),
        }
        if not row["sample"] or not row["condition"] or not row["out"]:
            continue
        rows.append(row)

    if not rows:
        raise ValueError("No usable rows found in sample sheet.")
    return rows


def _pick_representatives(rows: List[dict]) -> Dict[str, dict]:
    by_cond: Dict[str, List[dict]] = {}
    for r in rows:
        by_cond.setdefault(r["condition"], []).append(r)

    chosen: Dict[str, dict] = {}
    for cond, cond_rows in by_cond.items():
        def rep_key(r: dict) -> Tuple[int, str]:
            try:
                rep_num = int(r["replicate"])
            except Exception:
                rep_num = 10**9
            return rep_num, r["sample"]

        chosen[cond] = sorted(cond_rows, key=rep_key)[0]
    return chosen


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Batch-generate bearing_hic_plot figures: each condition vs a reference condition.",
    )
    ap.add_argument("--sheet", required=True, help="Sample sheet TSV with sample/condition/replicate/out columns.")
    ap.add_argument("--regions-file", required=True, help="Regions template TSV for bearing_hic_plot.py --regions-file.")
    ap.add_argument("--reference-condition", required=True, help="Condition used as the comparator baseline.")

    ap.add_argument(
        "--contact",
        action="append",
        default=[],
        metavar="COND=PATH",
        help="Condition-specific Hi-C file mapping. Repeat this flag.",
    )
    ap.add_argument(
        "--loops",
        action="append",
        default=[],
        metavar="COND=BEDPE",
        help="Optional condition-specific loop BEDPE mapping. Repeat this flag.",
    )
    ap.add_argument(
        "--label",
        action="append",
        default=[],
        metavar="COND=LABEL",
        help="Optional display label mapping for conditions. Defaults to condition names.",
    )

    ap.add_argument("--results-dir", default="results", help="Directory containing <sample>.stats.tsv or <sample>.neglog10p.bw files.")
    ap.add_argument(
        "--no-pvals",
        action="store_true",
        help="Disable p-value tracks entirely (do not require or pass --pval-* inputs).",
    )
    ap.add_argument("--comparison-dir", default=None, help="Directory containing compare_qcat diff_*.qcat.bgz files (optional).")
    ap.add_argument("--outdir", default="hic_batch", help="Top-level output directory for plots.")
    ap.add_argument("--triangle", action="store_true", help="Use bearing_hic_plot_triangle.py for generated figures.")
    ap.add_argument("--combined", action="store_true", help="Use bearing_hic_combined_plot.py (combined Hi-C + BEARING browser) for generated figures.")
    ap.add_argument("--diff-pvals-only", action="store_true",
                    help="Render the differential (B-A) p-value track but do NOT "
                         "require or pass per-sample p-value tracks. Use when the "
                         "pipeline only computes per-comparison diff p-values.")

    ap.add_argument("--highlights", default=None, help="Optional BED for highlighted regions.")
    ap.add_argument("--gtf", default=None, help="Optional GTF for genes.")
    ap.add_argument("--genes", default=None, help="Optional BED6 for genes.")
    ap.add_argument(
        "--bed",
        action="append",
        default=[],
        metavar="FILE",
        help=(
            "BED file to overlay as a track row when generating triangle figures. "
            "Repeat for multiple files."
        ),
    )
    ap.add_argument(
        "--bed-style",
        action="append",
        default=[],
        metavar="FILE=STYLE",
        help=(
            "Per-file rendering style override for --bed entries. STYLE is 'cbe' or 'itemRgb'. "
            "Repeat for each file you want to override."
        ),
    )
    ap.add_argument("--pval-cutoff", type=float, default=0.05, help="p-value cutoff for overlay/significance tracks.")
    ap.add_argument("--pval-overlay", action="store_true", help="Enable p-value overlay mode in bearing_hic_plot.py.")
    ap.add_argument("--rgb-hic", action="store_true", help="Enable RGB Hi-C mode in bearing_hic_plot.py.")
    ap.add_argument(
        "--rgb-palette",
        choices=["magenta-green", "red-green", "blue-red", "green-blue", "magenta-green-white"],
        default="magenta-green",
        help=(
            "Color mapping used when --rgb-hic is enabled. "
            "Default: magenta-green."
        ),
    )
    ap.add_argument(
        "--categories",
        default=None,
        help=("Optional categories file (YAML or JSON) to forward to bearing_hic_plot.py. "
              "If omitted, the script will auto-detect YAML files under the categories/ directory."),
    )
    ap.add_argument("--hic-format", choices=["auto", "cooler", "juicer"], default="auto")
    ap.add_argument("--python", default=sys.executable, help="Python executable to use for child calls.")
    ap.add_argument("--run", action="store_true", help="Execute commands. Without this, only print commands (dry-run).")
    ap.add_argument(
        "--summary-tsv",
        default=None,
        help=(
            "Optional output TSV path for per-region significant-bin summary. "
            "Default: <outdir>/significant_bins_summary.tsv"
        ),
    )

    args = ap.parse_args()

    script_dir = Path(__file__).resolve().parent
    sheet_path = _resolve(args.sheet, Path.cwd())
    regions_path = _resolve(args.regions_file, Path.cwd())
    results_dir = _resolve(args.results_dir, Path.cwd())
    outdir = _resolve(args.outdir, Path.cwd())
    bearing_plot = script_dir / (
        "bearing_hic_combined_plot.py" if getattr(args, "combined", False)
        else ("bearing_hic_plot_triangle.py" if args.triangle else "bearing_hic_plot.py"))

    if not sheet_path.exists():
        raise SystemExit(f"ERROR: missing sheet: {sheet_path}")
    if not regions_path.exists():
        raise SystemExit(f"ERROR: missing regions file: {regions_path}")
    if not bearing_plot.exists():
        raise SystemExit(f"ERROR: missing plot script next to this script: {bearing_plot}")

    try:
        contact_map = _parse_key_value(args.contact, "--contact")
        loops_map = _parse_key_value(args.loops, "--loops") if args.loops else {}
        label_map = _parse_key_value(args.label, "--label") if args.label else {}
    except ValueError as e:
        raise SystemExit(f"ERROR: {e}")

    # Early validation for explicit optional file arguments.
    optional_paths = []
    if args.highlights:
        optional_paths.append(_resolve(args.highlights, Path.cwd()))
    if args.gtf:
        optional_paths.append(_resolve(args.gtf, Path.cwd()))
    if args.genes:
        optional_paths.append(_resolve(args.genes, Path.cwd()))
    if args.bed:
        for b in args.bed:
            optional_paths.append(_resolve(b, Path.cwd()))
    try:
        require_files(optional_paths, context="optional input")
    except FileNotFoundError as exc:
        raise SystemExit(f"ERROR: {exc}")

    try:
        require_mapping_paths(
            contact_map,
            lambda p: _resolve(p, Path.cwd()),
            context="--contact input",
        )
    except FileNotFoundError as exc:
        raise SystemExit(f"ERROR: {exc}")

    try:
        require_mapping_paths(
            loops_map,
            lambda p: _resolve(p, Path.cwd()),
            context="--loops input",
        )
    except FileNotFoundError as exc:
        raise SystemExit(f"ERROR: {exc}")

    rows = _load_rows(sheet_path)
    reps = _pick_representatives(rows)

    ref = args.reference_condition
    if ref not in reps:
        raise SystemExit(
            f"ERROR: reference condition '{ref}' not present in sample sheet. "
            f"Available: {', '.join(sorted(reps))}"
        )

    missing_contacts = [c for c in reps if c not in contact_map]
    if missing_contacts:
        raise SystemExit(
            "ERROR: missing --contact mappings for conditions: "
            + ", ".join(sorted(missing_contacts))
        )

    outdir.mkdir(parents=True, exist_ok=True)

    def sample_qcat(sample_row: dict) -> Path:
        return _resolve(sample_row["out"], sheet_path.parent)

    def sample_pval(sample_row: dict) -> Path:
        stats_tsv = results_dir / f"{sample_row['sample']}.stats.tsv"
        if stats_tsv.exists():
            return stats_tsv
        return results_dir / f"{sample_row['sample']}.neglog10p.bw"

    ref_row = reps[ref]
    ref_qcat = sample_qcat(ref_row)
    ref_pval = sample_pval(ref_row)
    ref_contact = _resolve(contact_map[ref], Path.cwd())

    ref_required = [ref_qcat, ref_contact]
    if not args.no_pvals and not args.diff_pvals_only:
        ref_required.append(ref_pval)
    for p in ref_required:
        if not p.exists():
            raise SystemExit(f"ERROR: missing reference input: {p}")

    comparisons = [c for c in sorted(reps) if c != ref]
    if not comparisons:
        raise SystemExit("ERROR: no non-reference conditions found to compare.")

    print("=" * 72)
    print("batch_bearing_hic_plots.py")
    print(f"  Reference condition: {ref} ({ref_row['sample']})")
    print(f"  Regions file:        {regions_path}")
    print(f"  Comparisons:         {len(comparisons)}")
    print(f"  Plot mode:           {'TRIANGLE' if args.triangle else 'STANDARD'}")
    print(f"  Run mode:            {'EXECUTE' if args.run else 'DRY-RUN'}")
    print("=" * 72)

    n_ok = 0
    for cond in comparisons:
        cond_row = reps[cond]
        cond_qcat = sample_qcat(cond_row)
        cond_pval = sample_pval(cond_row)
        cond_contact = _resolve(contact_map[cond], Path.cwd())

        required_inputs = [cond_qcat, cond_contact]
        if not args.no_pvals and not args.diff_pvals_only:
            required_inputs.append(cond_pval)
        missing = [p for p in required_inputs if not p.exists()]
        if missing:
            print(f"\n[SKIP] {cond}: missing inputs", file=sys.stderr)
            for m in missing:
                print(f"  - {m}", file=sys.stderr)
            continue

        pair_outdir = outdir / f"{_safe(ref)}_vs_{_safe(cond)}"
        pair_outdir.mkdir(parents=True, exist_ok=True)

        diff_qcat = None
        diff_pval = None
        if args.comparison_dir:
            comp_dir = _resolve(args.comparison_dir, Path.cwd())
            diff_qcat, diff_pval = _prepare_forward_oriented_diff_assets(
                cond=cond,
                ref=ref,
                comp_dir=comp_dir,
                results_dir=results_dir,
                outdir=outdir,
            )

        # Select categories file to forward to the plot script.
        # Priority: explicit CLI --categories -> auto-detect YAML in categories/ -> none
        categories_yaml = None
        if args.categories:
            try:
                categories_yaml = _resolve(args.categories, Path.cwd())
                if not categories_yaml.exists():
                    raise FileNotFoundError(categories_yaml)
                print(f"[INFO] Using categories file from CLI: {categories_yaml}", file=sys.stderr)
            except Exception:
                raise SystemExit(f"ERROR: --categories path not found: {args.categories}")
        else:
            cats_dir = script_dir / "categories"
            if cats_dir.exists() and cats_dir.is_dir():
                yfiles = sorted([p for p in cats_dir.iterdir() if p.is_file() and p.suffix.lower() in (".yaml", ".yml")])
                if yfiles:
                    preferred = None
                    for name in ("mm10_15state.yaml", "mm10_15state.yml", "hg38_15state.yaml", "hg38_15state.yml"):
                        cand = cats_dir / name
                        if cand in yfiles:
                            preferred = cand
                            break
                    categories_yaml = preferred or yfiles[0]
                    print(f"[INFO] Using categories YAML: {categories_yaml}", file=sys.stderr)

        cmd = [
            args.python,
            str(bearing_plot),
            "--contact-a", str(ref_contact),
            "--contact-b", str(cond_contact),
            "--qcat-a", str(ref_qcat),
            "--qcat-b", str(cond_qcat),
            "--label-a", label_map.get(ref, ref),
            "--label-b", label_map.get(cond, cond),
            "--regions-file", str(regions_path),
            "--outdir", str(pair_outdir),
        ]

        if categories_yaml:
            cmd.extend(["--categories", str(categories_yaml)])

        if not args.no_pvals and not args.diff_pvals_only:
            cmd.extend([
                "--pval-a", str(ref_pval),
                "--pval-b", str(cond_pval),
                "--pval-cutoff", str(args.pval_cutoff),
            ])

        if args.hic_format != "auto":
            cmd.extend(["--format", args.hic_format])
        if args.pval_overlay and not args.no_pvals:
            cmd.append("--pval-overlay")
        elif args.pval_overlay and args.no_pvals:
            print("[WARN] --pval-overlay ignored because --no-pvals is set.", file=sys.stderr)
        if args.rgb_hic:
            cmd.append("--rgb-hic")
            cmd.extend(["--rgb-palette", args.rgb_palette])
        if diff_qcat is not None:
            cmd.extend(["--diff-qcat", str(diff_qcat)])
        if diff_pval is not None and not args.no_pvals:
            cmd.extend(["--pval-diff", str(diff_pval)])
        if args.highlights:
            cmd.extend(["--highlights", str(_resolve(args.highlights, Path.cwd()))])
        if args.gtf:
            cmd.extend(["--gtf", str(_resolve(args.gtf, Path.cwd()))])
        if args.genes:
            cmd.extend(["--genes", str(_resolve(args.genes, Path.cwd()))])
        if args.bed:
            for b in args.bed:
                cmd.extend(["--bed", str(_resolve(b, Path.cwd()))])
        if args.bed_style:
            for s in args.bed_style:
                cmd.extend(["--bed-style", s])
        if ref in loops_map:
            cmd.extend(["--loops-a", str(_resolve(loops_map[ref], Path.cwd()))])
        if cond in loops_map:
            cmd.extend(["--loops-b", str(_resolve(loops_map[cond], Path.cwd()))])

        print(f"\n# {ref} vs {cond}  (samples: {ref_row['sample']} vs {cond_row['sample']})")
        rc = execute_command(cmd, run=args.run, label=f"{ref} vs {cond}")
        if rc != 0:
            print(f"[ERROR] command failed for {ref} vs {cond} (exit {rc})", file=sys.stderr)
            continue

        n_ok += 1

    print(f"\nCompleted comparisons: {n_ok}/{len(comparisons)}")
    if not args.run:
        print("(dry-run only; rerun with --run to execute)")

    # --- Region x condition significance summary ---
    if args.no_pvals:
        print("\nSkipping significant-bin summary because --no-pvals was requested.")
        return 0

    try:
        regions = _load_regions(regions_path)
    except Exception as e:
        print(f"\n[WARN] Could not load regions for summary: {e}", file=sys.stderr)
        return 0

    cutoff_neglog10 = -math.log10(args.pval_cutoff)
    cond_order = [ref] + [c for c in sorted(reps) if c != ref]
    summary_rows = []

    for reg in regions:
        region_name = reg["name"]
        region_str = reg["region"]
        try:
            chrom, start, end = _parse_region(region_str)
        except Exception:
            print(f"[WARN] Skipping malformed region string: {region_str}", file=sys.stderr)
            continue

        row = {
            "region_name": region_name,
            "region": region_str,
        }
        for cond in cond_order:
            pval_bw = sample_pval(reps[cond])
            if not pval_bw.exists():
                row[cond] = "NA"
                continue
            try:
                row[cond] = str(_count_significant_bins_in_region(
                    pval_bw, chrom, start, end, cutoff_neglog10
                ))
            except Exception as e:
                print(
                    f"[WARN] Summary count failed for {cond} in {region_name}: {e}",
                    file=sys.stderr,
                )
                row[cond] = "NA"

        summary_rows.append(row)

    if summary_rows:
        header = ["region_name", "region"] + cond_order
        print("\n" + "=" * 72)
        print(
            "Significant-bin summary per region "
            f"(threshold: -log10(p) >= {cutoff_neglog10:.3f})"
        )
        print("\t".join(header))
        for r in summary_rows:
            print("\t".join(r.get(h, "") for h in header))

        summary_path = (
            _resolve(args.summary_tsv, Path.cwd())
            if args.summary_tsv else
            (outdir / "significant_bins_summary.tsv")
        )
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        with open(summary_path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=header, delimiter="\t")
            writer.writeheader()
            for r in summary_rows:
                writer.writerow({k: r.get(k, "") for k in header})
        print(f"Summary TSV: {summary_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
