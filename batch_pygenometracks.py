#!/usr/bin/env python3
"""
batch_pygenometracks.py
=======================
Run pyGenomeTracks for every sample in a sheet across one or more regions.

This avoids repeating commands like:
    pyGenomeTracks --tracks nina_tracks.ini --region chr6:... -o out.png

SAMPLE SHEET
------------
TSV with header. Required columns:
  - sample
And one of:
  - tracks_ini
  - qcat (auto-derives tracks ini as <base>_tracks.ini)

Optional columns:
  - enabled (1/0, true/false)

REGIONS
-------
Provide either:
  1) --region NAME=chr:start-end    (repeatable)
  2) --regions-file TSV with columns: name, region

OUTPUT NAMING
-------------
Default output file pattern:
  {sample}_{name}.png
Override with --out-template.

USAGE
-----
python batch_pygenometracks.py \
  --sheet samples.tsv \
  --region wide=chr6:30000000-50000000 \
  --region tcrb=chr6:40793981-41688054 \
  --region rc=chr6:41525453-41567200 \
  --outdir plots \
  --run

# Dry run (print commands only)
python batch_pygenometracks.py --sheet samples.tsv --region v1=chr6:40829239-41076570
"""

import argparse
import sys
from pathlib import Path

from bearing.utils import (
    iter_non_comment_lines,
    parse_named_region_item,
    sanitize_token,
)
from bearing.tsv import read_tsv_dict_rows
from bearing.tsv import read_tsv_table
from bearing.diagnostics import require_files
from bearing.runner import execute_commands


def _truthy(value):
    if value is None:
        return True
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def _non_comment_lines(handle):
    """Yield TSV lines excluding blanks and lines starting with '#' (after trim)."""
    yield from iter_non_comment_lines(handle)


def parse_sheet(path):
    rows = []
    try:
        fields, raw_rows = read_tsv_table(Path(path))
    except ValueError as exc:
        sys.exit(f"ERROR: {exc}")

    has_sample = "sample" in fields
    has_ini = "tracks_ini" in fields
    has_qcat = ("qcat" in fields) or ("out" in fields)
    if not has_sample:
        sys.exit("ERROR: sheet must include a 'sample' column")
    if not (has_ini or has_qcat):
        sys.exit("ERROR: sheet must include 'tracks_ini', 'qcat', or 'out' column")

    for row in raw_rows:
        if not row.get("sample"):
            continue
        if not _truthy(row.get("enabled", "1")):
            continue

        qcat_path = row.get("qcat") or row.get("out")

        if row.get("tracks_ini"):
            ini = Path(row["tracks_ini"])
        elif qcat_path:
            qcat = qcat_path
            base = qcat.replace(".qcat.bgz", "").replace(".bgz", "")
            ini = Path(base + "_tracks.ini")
        else:
            continue

        rows.append({
            "sample": row["sample"],
            "tracks_ini": ini,
            "pvalue_bw": row.get("pvalue_bw") or None,
        })

    if not rows:
        sys.exit("ERROR: no usable sample rows found in sheet")
    return rows


def parse_region_item(item):
    try:
        return parse_named_region_item(item, arg_name="region")
    except ValueError as exc:
        sys.exit(f"ERROR: {exc}")


def parse_regions_file(path):
    regions = []
    try:
        rows = read_tsv_dict_rows(Path(path), required_columns=["name", "region"])
    except ValueError as exc:
        msg = str(exc).replace("TSV", "regions file")
        sys.exit(f"ERROR: {msg}")

    for row in rows:
        if not row.get("name") or not row.get("region"):
            continue
        regions.append({"name": row["name"], "region": row["region"]})
    return regions


def sanitize(text):
    return sanitize_token(text, allowed="-_.")


def derive_ini_with_inverted_pvalue(source_ini, outdir, sample, pvalue_bw,
                                     pval_cutoff=3.0):
    """
    Read a tracks.ini, keep the leading Epilogos section intact, and replace
    any "epilogos inverted" / "Epilogos inverted" section with a bigwig
    track pointing at the per-sample p-value BigWig.

    If no inverted-Epilogos section is found, the new section is appended
    at the end before any genes/x-axis sections.

    Returns the path of the derived ini.
    """
    if not pvalue_bw:
        # No p-value BigWig provided for this sample - return source ini
        # unchanged (pyGenomeTracks will use it as-is).
        return Path(source_ini)

    source_ini = Path(source_ini)
    if not source_ini.exists():
        return source_ini  # caller will handle the missing-file error

    derived = Path(outdir) / f"{sample}_tracks_with_pval.ini"
    text = source_ini.read_text()

    new_section_lines = [
        "",
        f"[pvalue_inverted_{sample}]",
        f"file = {pvalue_bw}",
        f"title = -log10(p) (inverted)",
        "color = #707070",
        "alpha = 0.85",
        "type = fill",
        "orientation = inverted",
        "min_value = 0",
        "max_value = 8",
        "height = 1.5",
        "show_data_range = true",
        "number_of_bins = 700",
    ]
    if pval_cutoff and pval_cutoff > 0:
        new_section_lines.extend([
            "",
            f"[pvalue_cutoff_line_{sample}]",
            "file_type = hlines",
            f"y_values = {pval_cutoff}",
            "color = #d97706",
            "line_style = dashed",
            "line_width = 1.0",
            "overlay_previous = share-y",
        ])

    new_section = "\n".join(new_section_lines) + "\n"

    # Try to find an existing inverted-Epilogos section and replace only
    # that section, leaving the leading Epilogos row intact.
    import re
    pattern = re.compile(
        r"(?ms)^\[[^\]]*[Ii]nverted[^\]]*\]\n.*?(?=^\[|\Z)",
        re.DOTALL,
    )
    if pattern.search(text):
        modified = pattern.sub(new_section.rstrip() + "\n", text, count=1)
    else:
        # No inverted-Epilogos section found: insert after the first spacer
        # following the main Epilogos block, or before genes/x-axis/end.
        insert_pattern = re.compile(
            r"(?ms)^\[[^\]]*(?:genes|x-axis|x_axis)[^\]]*\]",
            re.IGNORECASE,
        )
        m = insert_pattern.search(text)
        if m:
            modified = text[:m.start()] + new_section + text[m.start():]
        else:
            modified = text.rstrip() + "\n" + new_section

    derived.write_text(modified)
    return derived


def main():
    parser = argparse.ArgumentParser(
        description="Batch-run pyGenomeTracks for each sample/region from TSV."
    )
    parser.add_argument("--sheet", required=True, help="Sample TSV")
    parser.add_argument(
        "--region", action="append", default=[], metavar="NAME=CHR:START-END",
        help="Region spec (repeatable), e.g. --region tcrb=chr6:40793981-41688054",
    )
    parser.add_argument(
        "--regions-file", default=None, metavar="TSV",
        help="Optional TSV with columns: name, region",
    )
    parser.add_argument(
        "--outdir", default=".", metavar="DIR",
        help="Directory for output images (default: current directory)",
    )
    parser.add_argument(
        "--out-template", default="{sample}_{name}.pdf", metavar="STR",
        help="Output filename template (default: {sample}_{name}.pdf). "
             "Use .png extension for raster output.",
    )
    parser.add_argument(
        "--dpi", type=int, default=300, metavar="N",
        help="Output resolution in DPI. pyGenomeTracks default is 72; "
             "300 is print-quality. Default: 300.",
    )
    parser.add_argument(
        "--pval-cutoff", type=float, default=3.0, metavar="-LOG10P",
        help="Significance line threshold (in -log10(p) units) for the "
             "inverted p-value track. Default: 3.0 (i.e. p < 0.001). "
             "Set to 0 to disable.",
    )
    parser.add_argument(
        "--pval-dir", default=None, metavar="DIR",
        help="Optional directory containing per-sample neglog10p.bw files. "
             "If provided, the script will look for <sample>.neglog10p.bw "
             "inside this directory and use it as the p-value BigWig.",
    )
    parser.add_argument(
        "--threads", type=int, default=1,
        help="Parallel jobs to run (default: 1)",
    )
    parser.add_argument(
        "--run", action="store_true",
        help="Execute commands. Without this flag, only print commands.",
    )
    parser.add_argument(
        "--skip-missing-ini", action="store_true",
        help="Skip rows whose tracks ini is missing instead of failing.",
    )
    args = parser.parse_args()

    samples = parse_sheet(args.sheet)

    # If a p-value directory is provided, try to auto-discover per-sample
    # neglog10p.bigwig files named <sample>.neglog10p.bw and set
    # the row's pvalue_bw if not already provided in the sheet.
    if args.pval_dir:
        pval_dir = Path(args.pval_dir)
        if not pval_dir.exists():
            print(f"WARNING: pval dir does not exist: {pval_dir}")
        else:
            for s in samples:
                if s.get("pvalue_bw"):
                    continue
                candidate = pval_dir / f"{s['sample']}.neglog10p.bw"
                if candidate.exists():
                    s["pvalue_bw"] = str(candidate)
                else:
                    s["pvalue_bw"] = None

    regions = [parse_region_item(r) for r in args.region]
    if args.regions_file:
        regions.extend(parse_regions_file(args.regions_file))
    if not regions:
        sys.exit("ERROR: provide at least one --region or --regions-file")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    commands = []
    for s in samples:
        sample = sanitize(s["sample"])
        ini = Path(s["tracks_ini"])

        if not ini.exists():
            msg = f"Missing tracks ini: {ini} (sample={sample})"
            if args.skip_missing_ini:
                print("WARNING: " + msg)
                continue
            sys.exit("ERROR: " + msg)

        try:
            require_files([ini], context="tracks ini")
        except FileNotFoundError as exc:
            if args.skip_missing_ini:
                print("WARNING: " + str(exc))
                continue
            sys.exit(f"ERROR: {exc}")

        # Derive a per-run ini that replaces any inverted Epilogos section
        # with an inverted p-value bigwig track for this sample (if provided).
        ini_for_run = derive_ini_with_inverted_pvalue(
            ini,
            outdir,
            sample,
            s.get("pvalue_bw"),
            pval_cutoff=args.pval_cutoff,
        )

        for reg in regions:
            name = sanitize(reg["name"])
            region = reg["region"]
            out_name = args.out_template.format(sample=sample, name=name, region=region)
            out_path = outdir / out_name
            cmd = [
                "pyGenomeTracks",
                "--tracks", str(ini_for_run),
                "--region", region,
                "--dpi", str(args.dpi),
                "-o", str(out_path),
            ]
            commands.append(cmd)

    if not commands:
        sys.exit("ERROR: no commands generated")

    print(f"Generated {len(commands)} pyGenomeTracks command(s).")
    total, failures = execute_commands(commands, run=args.run, threads=max(1, args.threads))
    if not args.run:
        print("\nDry run complete. Add --run to execute.")
        return
    if failures:
        sys.exit(f"ERROR: {failures}/{total} command(s) failed")
    print("All pyGenomeTracks commands completed successfully.")


if __name__ == "__main__":
    main()
