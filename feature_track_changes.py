#!/usr/bin/env python3
"""
feature_track_changes.py -- per-feature, per-track BEARING differential changes
(with p-values) for a BED of features, across one or more comparisons.

Reports, for each known CTCF element (CBE) or AgR gene, how each chromatin track
(CTCF, Cohesin/RAD21, ATAC, H3K27ac, RNAseq+/-) shifts between the two conditions
of each comparison, with significance. The per-track value is BEARING's SIGNED
per-track contribution kl_<track>: > 0 = enriched in the first condition (e.g.
DN), < 0 = enriched in the second. It is a compositional differential (how much
of the bin's local signal composition that track accounts for), NOT a raw
fold-change in coverage.

Two modes, auto-selected by feature width:
  * point  (feature <= --point-max bp, e.g. 18 bp CBEs): report the single bin
    that contains the feature midpoint -- its per-track kl + p-value/FDR.
  * region (larger features, e.g. AgR gene bodies): aggregate over every bin
    overlapping the feature -- per-track sum / mean / peak, n_bins, n significant
    bins, and the smallest pval_adj_bh in the feature.

Streams each diff stats file ONCE, keeping only feature-overlapping bins (never
holds the genome-wide table in memory), so it is safe on the full 13.6M-bin files.

  python feature_track_changes.py \
    --features annotations/cbe_mm10.bed \
    --diffs results/pvalue/diff_*.stats.tsv \
    --categories results/DN_rep1_cats.json \
    --out cbe_track_changes.tsv

ASCII only. Handles gzip diff files.
"""
import argparse
import bisect
import csv
import glob
import gzip
import json
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
    return re.sub(r"^diff_", "", base)


def clean_track(name):
    """Make a track name safe + readable as a column label."""
    n = name.strip().replace(" ", "")
    n = n.replace("+", "Pos").replace("-", "Neg").replace("/", "_")
    return n


def load_track_labels(categories):
    """kl_1.. -> clean track name, positionally from the categories JSON."""
    labels = {}
    if categories and os.path.exists(categories):
        try:
            obj = json.load(open(categories))
            cats = obj.get("categories", obj) if isinstance(obj, dict) else {}
            for i, k in enumerate(sorted(cats, key=lambda x: int(x))):
                v = cats[k]
                nm = v[0] if isinstance(v, (list, tuple)) else str(v)
                labels["kl_%d" % (i + 1)] = clean_track(nm)
        except (ValueError, KeyError, TypeError):
            labels = {}
    return labels


def read_features(path):
    feats = []
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            try:
                s, e = int(p[1]), int(p[2])
            except ValueError:
                continue
            name = (p[3].split("#")[0].strip() if len(p) > 3
                    else "%s:%d-%d" % (p[0], s, e))
            feats.append((p[0], s, e, name))
    return feats


def merged_by_chrom(feats):
    """Per chrom: (starts[], ends[]) of merged non-overlapping feature spans, for
    a fast 'does this bin touch any feature' pre-filter."""
    by = {}
    for (c, s, e, _n) in feats:
        by.setdefault(c, []).append((s, e))
    out = {}
    for c, iv in by.items():
        iv.sort()
        m = [list(iv[0])]
        for s, e in iv[1:]:
            if s <= m[-1][1]:
                m[-1][1] = max(m[-1][1], e)
            else:
                m.append([s, e])
        out[c] = ([x[0] for x in m], [x[1] for x in m])
    return out


def touches(merged, chrom, s, e):
    mc = merged.get(chrom)
    if not mc:
        return False
    starts, ends = mc
    i = bisect.bisect_right(starts, e - 1) - 1
    return i >= 0 and ends[i] > s


def feats_by_chrom_indexed(feats, mode, point_max):
    """Per chrom, features sorted by start with parallel starts[] for scanning.
    Each item: (fstart, fend, fmid, idx, is_point)."""
    by = {}
    for idx, (c, s, e, _n) in enumerate(feats):
        is_point = (mode == "point") or (mode == "auto" and (e - s) <= point_max)
        by.setdefault(c, []).append((s, e, (s + e) // 2, idx, is_point))
    out = {}
    for c, items in by.items():
        items.sort()
        out[c] = ([it[0] for it in items], items)
    return out


def assign(feats_idx, chrom, s, e):
    """Return list of (feature_idx, is_point) for features this bin serves:
    point features whose midpoint is in [s,e); region features overlapping [s,e)."""
    fc = feats_idx.get(chrom)
    if not fc:
        return []
    starts, items = fc
    hi = bisect.bisect_left(starts, e)  # features with fstart < e
    hits = []
    for j in range(hi):
        fs, fe, fmid, idx, is_point = items[j]
        if is_point:
            if s <= fmid < e:
                hits.append((idx, True))
        else:
            if fe > s:  # overlap
                hits.append((idx, False))
    return hits


def fnum(x):
    try:
        return float(x)
    except (ValueError, TypeError):
        return None


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--features", required=True, help="BED of features (CBE or AgR)")
    ap.add_argument("--diffs", nargs="+", required=True,
                    help="diff_<A>_vs_<B>.stats.tsv files (globs ok)")
    ap.add_argument("--categories", default=None, help="categories JSON for labels")
    ap.add_argument("--mode", choices=["auto", "point", "region"], default="auto")
    ap.add_argument("--point-max", type=int, default=200,
                    help="auto mode: features <= this width are point (default 200)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    diffs = []
    for d in args.diffs:
        diffs.extend(sorted(glob.glob(d)) if any(ch in d for ch in "*?[") else [d])
    if not diffs:
        sys.exit("[ERROR] no diff files matched")

    feats = read_features(args.features)
    if not feats:
        sys.exit("[ERROR] no features parsed from %s" % args.features)
    labels = load_track_labels(args.categories)
    merged = merged_by_chrom(feats)
    feats_idx = feats_by_chrom_indexed(feats, args.mode, args.point_max)
    is_point_feat = {}
    for c in feats_idx:
        for (_fs, _fe, _fm, idx, isp) in feats_idx[c][1]:
            is_point_feat[idx] = isp

    kl_order = []
    point_rows, region_rows = [], []

    for dp in diffs:
        comp = parse_comparison_name(dp)
        fh = _open(dp)
        rd = csv.reader(fh, delimiter="\t")
        header = next(rd)
        col = {n: i for i, n in enumerate(header)}
        fdr_col = next((c for c in header if c.startswith("significant_fdr")), None)
        kl_cols = [c for c in header if c.startswith("kl_")]
        for k in kl_cols:
            if k not in kl_order:
                kl_order.append(k)
        try:
            ci_c, ci_s, ci_e = col["chrom"], col["start"], col["end"]
        except KeyError:
            sys.exit("[ERROR] %s lacks chrom/start/end" % dp)

        # per-feature region accumulators for this comparison
        acc = {}     # idx -> dict
        point_hit = {}  # idx -> parts (for point features)

        for parts in rd:
            if len(parts) <= ci_e:
                continue
            chrom = parts[ci_c]
            if not touches(merged, chrom, int(parts[ci_s]) if parts[ci_s].lstrip("-").isdigit() else -1,
                           int(parts[ci_e]) if parts[ci_e].lstrip("-").isdigit() else -1):
                continue
            try:
                s = int(parts[ci_s]); e = int(parts[ci_e])
            except ValueError:
                continue
            for idx, is_point in assign(feats_idx, chrom, s, e):
                if is_point:
                    point_hit[idx] = parts
                    continue
                a = acc.setdefault(idx, {"n": 0, "nsig": 0, "minp": None,
                                         "bsum": 0.0, "tsum": {}, "tpeak": {}})
                a["n"] += 1
                if fdr_col and parts[col[fdr_col]] in ("1", "1.0"):
                    a["nsig"] += 1
                padj = fnum(parts[col["pval_adj_bh"]]) if "pval_adj_bh" in col else None
                if padj is not None:
                    a["minp"] = padj if a["minp"] is None else min(a["minp"], padj)
                bs = fnum(parts[col["bearing_score"]]) if "bearing_score" in col else None
                if bs is not None:
                    a["bsum"] += bs
                for k in kl_cols:
                    val = fnum(parts[col[k]])
                    if val is None:
                        continue
                    a["tsum"][k] = a["tsum"].get(k, 0.0) + val
                    if abs(val) > abs(a["tpeak"].get(k, 0.0)):
                        a["tpeak"][k] = val
        fh.close()

        # emit point rows
        for idx, (c, fs, fe, name) in enumerate(feats):
            if not is_point_feat.get(idx, False):
                continue
            rec = {"comparison": comp, "feature": name, "chrom": c,
                   "feat_start": fs, "feat_end": fe}
            parts = point_hit.get(idx)
            if parts is None:
                rec["bin_status"] = "no_bin"
                point_rows.append(rec)
                continue
            rec["bin_status"] = "ok"
            rec["bin_start"] = parts[col["start"]]
            rec["bin_end"] = parts[col["end"]]
            rec["bearing_score"] = parts[col["bearing_score"]] if "bearing_score" in col else ""
            rec["direction"] = parts[col["direction"]] if "direction" in col else ""
            rec["pval"] = parts[col["pval"]] if "pval" in col else ""
            rec["pval_adj_bh"] = parts[col["pval_adj_bh"]] if "pval_adj_bh" in col else ""
            rec["fdr_significant"] = parts[col[fdr_col]] if fdr_col else ""
            for k in kl_cols:
                rec[labels.get(k, k)] = parts[col[k]]
            point_rows.append(rec)

        # emit region rows
        for idx, (c, fs, fe, name) in enumerate(feats):
            if is_point_feat.get(idx, False):
                continue
            rec = {"comparison": comp, "feature": name, "chrom": c,
                   "feat_start": fs, "feat_end": fe}
            a = acc.get(idx)
            if a is None or a["n"] == 0:
                rec.update({"n_bins": 0, "n_sig_bins": 0, "min_pval_adj_bh": "",
                            "bearing_sum": ""})
                region_rows.append(rec)
                continue
            rec["n_bins"] = a["n"]
            rec["n_sig_bins"] = a["nsig"]
            rec["min_pval_adj_bh"] = "%.4g" % a["minp"] if a["minp"] is not None else ""
            rec["bearing_sum"] = "%.4f" % a["bsum"]
            for k in kl_cols:
                lab = labels.get(k, k)
                sm = a["tsum"].get(k, 0.0)
                rec["%s_sum" % lab] = "%.4f" % sm
                rec["%s_mean" % lab] = "%.4f" % (sm / a["n"])
                rec["%s_peak" % lab] = "%.4f" % a["tpeak"].get(k, 0.0)
            region_rows.append(rec)

    labeled = [labels.get(k, k) for k in kl_order]
    n_point = sum(1 for v in is_point_feat.values() if v)
    n_region = len(is_point_feat) - n_point

    wrote = []
    if point_rows:
        base = ["comparison", "feature", "chrom", "feat_start", "feat_end",
                "bin_status", "bin_start", "bin_end", "bearing_score", "direction",
                "pval", "pval_adj_bh", "fdr_significant"]
        out_p = args.out if not region_rows else _suffix(args.out, "_point")
        _write(out_p, base + labeled, point_rows)
        wrote.append((out_p, len(point_rows), "point"))
    if region_rows:
        base = ["comparison", "feature", "chrom", "feat_start", "feat_end",
                "n_bins", "n_sig_bins", "min_pval_adj_bh", "bearing_sum"]
        cols = base + ["%s_%s" % (lab, suf) for lab in labeled
                       for suf in ("sum", "mean", "peak")]
        out_r = args.out if not point_rows else _suffix(args.out, "_region")
        _write(out_r, cols, region_rows)
        wrote.append((out_r, len(region_rows), "region"))

    sys.stderr.write("Features: %d point, %d region | %d comparisons\n"
                     % (n_point, n_region, len(diffs)))
    for path, n, kind in wrote:
        sys.stderr.write("  wrote %s (%s, %d rows)\n" % (path, kind, n))


def _suffix(path, suf):
    root, ext = os.path.splitext(path)
    return root + suf + ext


def _write(path, cols, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


if __name__ == "__main__":
    main()
