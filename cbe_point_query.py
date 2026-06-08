#!/usr/bin/env python3
"""
cbe_point_query.py
==================
Per-CBE point query. CBEs are ~18 bp CTCF-binding-element motifs, far shorter
than the 200 bp BEARING bins, so a regional-enrichment test over a CBE is
structurally empty (zero bins fall *within* an 18 bp region). The honest query
for a point feature is therefore: which 200 bp bin does each CBE fall in, and
what is that bin's BEARING differential score, p-value, FDR call, and per-track
decomposition?

For each CBE this reports the bin whose interval contains the CBE midpoint. If
no bin contains it (e.g. the bin was masked for low signal or blacklist), the
row is reported with bin_status='no_bin' and NA stats - itself informative.

Consumes one or more diff_<A>_vs_<B>.stats.tsv files (comparison parsed from
filename) and writes one long-format row per (comparison, CBE).

ASCII-only. Handles gzip-compressed stats files.
"""

import argparse
import csv
import gzip
import os
import re
import sys


def _open(path):
    with open(path, "rb") as fh:
        if fh.read(2) == b"\x1f\x8b":
            return gzip.open(path, "rt")
    return open(path, "r")


def parse_comparison_name(path):
    base = os.path.basename(path)
    base = re.sub(r"\.stats\.tsv(\.gz)?$", "", base)
    base = re.sub(r"^diff_", "", base)
    return base


def read_cbes(path):
    cbes = []
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            name = (p[3] if len(p) > 3 else "%s:%s-%s" % (p[0], p[1], p[2]))
            name = name.split("#")[0].strip()
            cbes.append((p[0], int(p[1]), int(p[2]), name))
    return cbes


def load_bins(path):
    """Return (rows, fields, fdr_col, kl_cols). rows are dicts with typed
    chrom/start/end retained as r['_chrom'], r['_start'], r['_end']."""
    fh = _open(path)
    reader = csv.DictReader(fh, delimiter="\t")
    fields = reader.fieldnames or []
    fdr_col = next((c for c in fields if c.startswith("significant_fdr")), None)
    kl_cols = [c for c in fields if c.startswith("kl_")]
    rows = []
    for r in reader:
        try:
            r["_chrom"] = r["chrom"]
            r["_start"] = int(r["start"])
            r["_end"] = int(r["end"])
        except (KeyError, ValueError):
            continue
        rows.append(r)
    fh.close()
    return rows, fields, fdr_col, kl_cols


def index_by_chrom(rows):
    """Group bin rows by chrom, each sorted by start, with a parallel starts
    list for binary search."""
    import bisect
    by = {}
    for r in rows:
        by.setdefault(r["_chrom"], []).append(r)
    index = {}
    for c, rs in by.items():
        rs.sort(key=lambda x: x["_start"])
        starts = [x["_start"] for x in rs]
        index[c] = (rs, starts)
    return index, bisect


def find_bin(index, bisect, chrom, pos):
    """Return the bin row whose [start,end) contains pos, or None."""
    if chrom not in index:
        return None
    rs, starts = index[chrom]
    i = bisect.bisect_right(starts, pos) - 1
    if 0 <= i < len(rs):
        r = rs[i]
        if r["_start"] <= pos < r["_end"]:
            return r
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--diffs", nargs="+", required=True,
                    help="diff_<A>_vs_<B>.stats.tsv files (one per comparison)")
    ap.add_argument("--cbe-bed", required=True, help="CBE BED4 (chrom,start,end,name)")
    ap.add_argument("--categories", default=None,
                    help="optional categories JSON to label kl_N columns by track name")
    ap.add_argument("--out", required=True, help="output long-format TSV")
    args = ap.parse_args()

    # kl_1, kl_2, ... (numeric suffix, file order) -> kl_<TrackName>, mapped
    # POSITIONALLY from the categories JSON so it is robust to 0- or 1-indexed
    # keys. Columns already named (kl_CTCF) are left as-is.
    kl_label = {}
    if args.categories and os.path.exists(args.categories):
        try:
            import json
            obj = json.load(open(args.categories))
            cats = obj.get("categories", obj) if isinstance(obj, dict) else {}
            ordered = []
            for k in sorted(cats, key=lambda x: int(x)):
                v = cats[k]
                ordered.append(v[0] if isinstance(v, (list, tuple)) else str(v))
            for i, name in enumerate(ordered):
                kl_label["kl_%d" % (i + 1)] = "kl_%s" % name
        except (ValueError, KeyError, TypeError):
            kl_label = {}

    cbes = read_cbes(args.cbe_bed)

    # Determine track set from the first readable diff (for consistent columns).
    kl_union = []
    out_rows = []
    for diff_path in args.diffs:
        comp = parse_comparison_name(diff_path)
        rows, fields, fdr_col, kl_cols = load_bins(diff_path)
        for k in kl_cols:
            if k not in kl_union:
                kl_union.append(k)
        index, bisect = index_by_chrom(rows)
        for (chrom, start, end, name) in cbes:
            mid = (start + end) // 2
            b = find_bin(index, bisect, chrom, mid)
            rec = {
                "comparison": comp, "cbe_name": name, "cbe_chrom": chrom,
                "cbe_start": start, "cbe_end": end,
            }
            if b is None:
                rec.update({"bin_status": "no_bin", "bin_start": "", "bin_end": "",
                            "bearing_score": "", "direction": "", "pval": "",
                            "pval_adj_bh": "", "fdr_significant": ""})
                for k in kl_cols:
                    rec[k] = ""
            else:
                rec.update({
                    "bin_status": "ok",
                    "bin_start": b["_start"], "bin_end": b["_end"],
                    "bearing_score": b.get("bearing_score", ""),
                    "direction": b.get("direction", ""),
                    "pval": b.get("pval", ""),
                    "pval_adj_bh": b.get("pval_adj_bh", ""),
                    "fdr_significant": b.get(fdr_col, "") if fdr_col else "",
                })
                for k in kl_cols:
                    rec[k] = b.get(k, "")
            out_rows.append(rec)

    base_cols = ["comparison", "cbe_name", "cbe_chrom", "cbe_start", "cbe_end",
                 "bin_status", "bin_start", "bin_end", "bearing_score",
                 "direction", "pval", "pval_adj_bh", "fdr_significant"]
    labeled_kl = [kl_label.get(k, k) for k in kl_union]
    header = base_cols + labeled_kl
    # rename per-row kl keys to their labeled names so DictWriter aligns
    if kl_label:
        for rec in out_rows:
            for raw in list(rec.keys()):
                if raw in kl_label:
                    rec[kl_label[raw]] = rec.pop(raw)
    with open(args.out, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=header, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for rec in out_rows:
            w.writerow(rec)

    # Compact stderr summary: how many CBEs are FDR-significant per comparison.
    sys.stderr.write("Wrote per-CBE point query: %s (%d rows)\n"
                     % (args.out, len(out_rows)))
    by_comp = {}
    for rec in out_rows:
        c = rec["comparison"]
        sig = str(rec.get("fdr_significant", "")) in ("1", "1.0")
        nobin = rec.get("bin_status") == "no_bin"
        d = by_comp.setdefault(c, [0, 0, 0])
        d[0] += 1
        d[1] += 1 if sig else 0
        d[2] += 1 if nobin else 0
    for c, (n, nsig, nnb) in sorted(by_comp.items()):
        sys.stderr.write("  %-16s %d CBEs | %d FDR-sig | %d with no bin\n"
                         % (c, n, nsig, nnb))


if __name__ == "__main__":
    main()
