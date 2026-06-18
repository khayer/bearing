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


def build_cbe_index(cbes):
    """{chrom: sorted [(mid, cbe_idx), ...]} for the CBE midpoints."""
    by = {}
    for i, (chrom, start, end, _name) in enumerate(cbes):
        mid = (start + end) // 2
        by.setdefault(chrom, []).append((mid, i))
    for c in by:
        by[c].sort()
    return by


def query_diff_for_cbes(path, mids_by_chrom):
    """Stream a diff stats TSV ONCE, keeping only the bins whose [start,end)
    contains a CBE midpoint. Returns (header, col, fdr_col, kl_cols, hits) where
    hits maps cbe_idx -> the matching row (list of cells). Memory is O(#CBEs),
    not O(#bins) -- the genome-wide table is never held in memory."""
    import bisect
    fh = _open(path)
    reader = csv.reader(fh, delimiter="\t", lineterminator="\n")
    header = next(reader)
    col = {name: i for i, name in enumerate(header)}
    fdr_col = next((c for c in header if c.startswith("significant_fdr")), None)
    kl_cols = [c for c in header if c.startswith("kl_")]
    try:
        ci_c, ci_s, ci_e = col["chrom"], col["start"], col["end"]
    except KeyError:
        fh.close()
        raise SystemExit("[ERROR] %s lacks chrom/start/end columns" % path)
    hits = {}
    for parts in reader:
        if len(parts) <= ci_e:
            continue
        mids = mids_by_chrom.get(parts[ci_c])
        if not mids:
            continue
        try:
            s = int(parts[ci_s]); e = int(parts[ci_e])
        except ValueError:
            continue
        # CBE midpoints with mid >= s, scanning upward while mid < e
        j = bisect.bisect_left(mids, (s, -1))
        while j < len(mids) and mids[j][0] < e:
            hits[mids[j][1]] = parts
            j += 1
    fh.close()
    return header, col, fdr_col, kl_cols, hits


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--diffs", nargs="+", required=True,
                    help="diff_<A>_vs_<B>.stats.tsv files (one per comparison)")
    ap.add_argument("--cbe-bed", required=True, help="CBE BED4 (chrom,start,end,name)")
    ap.add_argument("--categories", default=None,
                    help="optional categories JSON to label kl_N columns by track name")
    ap.add_argument("--out", required=True, help="output long-format TSV")
    ap.add_argument("--summary-out", default=None,
                    help="optional per-comparison significance summary TSV "
                         "(n CBEs, n FDR-sig, n no-bin, significant CBE names)")
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
    mids_by_chrom = build_cbe_index(cbes)

    # Stream each diff once, pulling only the CBE-containing bins.
    kl_union = []
    out_rows = []
    for diff_path in args.diffs:
        comp = parse_comparison_name(diff_path)
        header, col, fdr_col, kl_cols, hits = query_diff_for_cbes(
            diff_path, mids_by_chrom)
        for k in kl_cols:
            if k not in kl_union:
                kl_union.append(k)
        for ci, (chrom, start, end, name) in enumerate(cbes):
            rec = {
                "comparison": comp, "cbe_name": name, "cbe_chrom": chrom,
                "cbe_start": start, "cbe_end": end,
            }
            parts = hits.get(ci)
            if parts is None:
                rec.update({"bin_status": "no_bin", "bin_start": "", "bin_end": "",
                            "bearing_score": "", "direction": "", "pval": "",
                            "pval_adj_bh": "", "fdr_significant": ""})
                for k in kl_cols:
                    rec[k] = ""
            else:
                def cell(name, default=""):
                    i = col.get(name)
                    return parts[i] if i is not None and i < len(parts) else default
                rec.update({
                    "bin_status": "ok",
                    "bin_start": cell("start"), "bin_end": cell("end"),
                    "bearing_score": cell("bearing_score"),
                    "direction": cell("direction"),
                    "pval": cell("pval"),
                    "pval_adj_bh": cell("pval_adj_bh"),
                    "fdr_significant": cell(fdr_col) if fdr_col else "",
                })
                for k in kl_cols:
                    rec[k] = cell(k)
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
    sig_names = {}
    for rec in out_rows:
        c = rec["comparison"]
        sig = str(rec.get("fdr_significant", "")) in ("1", "1.0")
        nobin = rec.get("bin_status") == "no_bin"
        d = by_comp.setdefault(c, [0, 0, 0])
        d[0] += 1
        d[1] += 1 if sig else 0
        d[2] += 1 if nobin else 0
        if sig:
            sig_names.setdefault(c, []).append(rec["cbe_name"])
    for c, (n, nsig, nnb) in sorted(by_comp.items()):
        sys.stderr.write("  %-16s %d CBEs | %d FDR-sig | %d with no bin\n"
                         % (c, n, nsig, nnb))
    total_sig = sum(v[1] for v in by_comp.values())
    sys.stderr.write("TOTAL FDR-significant CBE x comparison cells: %d\n" % total_sig)

    if args.summary_out:
        with open(args.summary_out, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["comparison", "n_cbe", "n_fdr_sig", "n_no_bin",
                        "significant_cbes"])
            for c in sorted(by_comp):
                n, nsig, nnb = by_comp[c]
                w.writerow([c, n, nsig, nnb,
                            ";".join(sorted(sig_names.get(c, []))) or "."])
        sys.stderr.write("Wrote significance summary: %s\n" % args.summary_out)


if __name__ == "__main__":
    main()
