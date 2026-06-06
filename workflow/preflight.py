#!/usr/bin/env python3
"""
preflight.py

Validate that every input file referenced by a BEARING run exists BEFORE
launching the pipeline, so a missing file fails in seconds instead of hours in.

Checks:
  - every BigWig in the sheet's `bw` column (per sample)
  - every Hi-C file in cool/insul/pca1, expanded over config hic.resolutions
  - the reference files named in the config (chrom_sizes, genes, gtf, etc.)
  - the categories YAML and regions file

Usage:
  python3 workflow/preflight.py --configfile workflow/config/config.yaml
  python3 workflow/preflight.py --configfile workflow/config/config.yaml --core-only
  # resolve relative paths against a specific directory (e.g. your data dir):
  python3 workflow/preflight.py --configfile workflow/config/config.yaml --base-dir /path/to/run

Exit 0 = all required inputs present. Exit 1 = at least one missing (listed).
--core-only skips Hi-C inputs (use when running the no-Hi-C path).
ASCII-only.
"""

import argparse
import os
import sys

import yaml


def load_sheet(path):
    rows = []
    with open(path) as fh:
        header = None
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            parts = line.split("\t")
            if header is None:
                header = parts
                continue
            rows.append(dict(zip(header, [p.strip() for p in parts])))
    return rows


def resolve(path, base_dir):
    if not path:
        return path
    if os.path.isabs(path):
        return path
    return os.path.normpath(os.path.join(base_dir, path))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--configfile", required=True)
    ap.add_argument("--base-dir", default=None,
                    help="directory relative paths resolve against "
                         "(default: the directory you run this from)")
    ap.add_argument("--core-only", action="store_true",
                    help="skip Hi-C inputs (cool/insul/pca1)")
    args = ap.parse_args()

    with open(args.configfile) as fh:
        cfg = yaml.safe_load(fh)

    cfgdir = os.path.dirname(os.path.abspath(args.configfile))
    # config paths like "config/samples.tsv" are written relative to the
    # workflow/ directory, which is the PARENT of the configfile's dir
    # (workflow/config/). Resolve against workflow/.
    wfdir = os.path.dirname(cfgdir) if os.path.basename(cfgdir) == "config" else cfgdir
    base = args.base_dir or os.getcwd()

    def cfg_path(p):
        if not p:
            return ""
        return p if os.path.isabs(p) else os.path.normpath(os.path.join(wfdir, p))

    missing = []
    checked = 0

    def check(path, label):
        nonlocal checked
        if not path:
            return
        checked += 1
        if not os.path.exists(path):
            missing.append((label, path))

    # --- reference files from config ---
    for key in ["chrom_sizes", "genes_bed", "agr_genes_bed", "gtf",
                "blacklist_external", "categories", "regions_file"]:
        val = cfg.get(key, "")
        if val:
            check(cfg_path(val), "config:%s" % key)

    # --- sample sheet inputs ---
    sheet_path = cfg_path(cfg["samples_sheet"])
    if not os.path.exists(sheet_path):
        print("FATAL: sample sheet not found: %s" % sheet_path, file=sys.stderr)
        sys.exit(1)
    rows = load_sheet(sheet_path)

    hic_cfg = (cfg.get("hic", {}) or {})
    resolutions = hic_cfg.get("resolutions", []) or []
    # Per-column resolution requirements: a Hi-C column is only read at the
    # resolutions the rules that consume it actually use. Checking the full
    # global grid would flag files (e.g. pca1 at 10kb) that no rule reads.
    comp_res = hic_cfg.get("compartment_resolutions", resolutions) or resolutions
    tad_res = hic_cfg.get("tad_resolutions", resolutions) or resolutions
    contact_res = hic_cfg.get("contact_resolution", None)
    # cool matrices are read at: contact_resolution (figures / contact isolation
    # / crosslocus) and the O/E grid resolutions (25k, 50k, 100k). They are NOT
    # read at the coarse compartment resolutions, so don't demand those.
    cool_res = sorted({int(x) for x in (
        ([contact_res] if contact_res else []) + [25000, 50000, 100000]
    )}) or resolutions
    col_res = {"cool": cool_res, "insul": tad_res, "pca1": comp_res}

    n_bearing = 0
    n_hic = 0
    for r in rows:
        sample = r.get("sample", "?")
        bw = r.get("bw", "")
        if bw:
            n_bearing += 1
            for p in bw.split(","):
                check(resolve(p.strip(), base), "bw[%s]" % sample)
        if not args.core_only:
            # A row may override which resolutions it needs via an optional
            # `resolutions` column (comma/space-separated), e.g. Hi-C-only
            # backgrounds that exist at only a subset of resolutions. When set,
            # it applies to all that row's Hi-C columns; otherwise each column
            # is checked at the resolutions its consuming rule actually uses.
            row_res = r.get("resolutions", "").replace(",", " ").split()
            row_res = [x for x in row_res if x]
            for col in ["cool", "insul", "pca1"]:
                pat = r.get(col, "")
                if not pat:
                    continue
                if col == "cool":
                    n_hic += 1
                res_list = row_res or col_res.get(col, resolutions)
                if "{res}" in pat:
                    for res in res_list:
                        check(resolve(pat.replace("{res}", str(res)), base),
                              "%s[%s]@%s" % (col, sample, res))
                else:
                    check(resolve(pat, base), "%s[%s]" % (col, sample))

    # --- report ---
    print("Preflight: %d sample(s) with BigWigs, %d Hi-C condition(s)%s"
          % (n_bearing, n_hic, " (Hi-C skipped: --core-only)" if args.core_only else ""))
    print("Checked %d file path(s); base dir for relative paths: %s"
          % (checked, base))
    if missing:
        print("\nMISSING %d file(s):" % len(missing), file=sys.stderr)
        for label, path in missing:
            print("  [%s] %s" % (label, path), file=sys.stderr)
        print("\nFix these (or pass --base-dir / use absolute paths) before "
              "launching the pipeline.", file=sys.stderr)
        sys.exit(1)
    print("\nAll required inputs present. Safe to launch.")


if __name__ == "__main__":
    main()
