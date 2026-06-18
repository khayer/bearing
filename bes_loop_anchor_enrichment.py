#!/usr/bin/env python3
"""
bes_loop_anchor_enrichment.py  (#2 -- anchor-stratified BES)
===========================================================
Test whether BEARING differential scores concentrate at chromatin-loop anchors
relative to non-anchor bins. This pulls BEARING into topology by asking: does
chromatin that participates in a called loop also show compositional change?

For one comparison (one diff_<A>_vs_<B>.stats.tsv) and a loop BEDPE, every
200 bp bin is labelled ANCHOR if it overlaps either anchor of any loop, else
NON-ANCHOR. We then compare the two strata on:
  - |BES| magnitude            -> Mann-Whitney U (two-sided)
  - FDR-significant fraction    -> Fisher exact (2x2)
and report per-track mean |kl| in each stratum (which track drives any
anchor enrichment -- expected CTCF/Cohesin if topology-linked).

Anchors are 5-25 kb wide (Mustache multi-resolution), so each spans many
200 bp bins; this is a regional stratum, not a point query. Loops are filtered
to the locus spanned by the diff table.

ASCII-only. scipy is optional; if unavailable, p-values are reported as NA and
only descriptive statistics are emitted (so the script never hard-fails in a
minimal environment -- the real run has scipy).
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
    return re.sub(r"^diff_", "", base)


def read_loops(path):
    """Return list of (chrom, start, end) anchor intervals (both A and B)."""
    anchors = []
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            try:
                anchors.append((p[0], int(p[1]), int(p[2])))
                anchors.append((p[3], int(p[4]), int(p[5])))
            except ValueError:
                continue
    return anchors


def build_anchor_lookup(anchors):
    """Per-chrom sorted (start,end) for interval overlap testing."""
    by = {}
    for (c, s, e) in anchors:
        by.setdefault(c, []).append((s, e))
    for c in by:
        by[c].sort()
    return by


def overlaps_anchor(lookup, chrom, bstart, bend):
    import bisect
    if chrom not in lookup:
        return False
    iv = lookup[chrom]
    starts = [s for (s, e) in iv]
    # find anchors whose start <= bend; check those for overlap
    i = bisect.bisect_right(starts, bend)
    for j in range(i - 1, -1, -1):
        s, e = iv[j]
        if e <= bstart:
            # since sorted by start, but ends not monotonic; scan a little
            if starts[j] < bstart - 10_000_000:
                break
            continue
        if s < bend and e > bstart:
            return True
    return False


def _track_namer(categories_path):
    """Return f(kl_col)->display name. Maps kl_1, kl_2, ... (numeric suffix, in
    file order) to track names POSITIONALLY from a categories JSON, so it is
    robust to whether the JSON keys are 0- or 1-indexed. Columns already named
    (kl_CTCF) or unmatched fall back to stripping the kl_ prefix."""
    mapping = {}
    if categories_path and os.path.exists(categories_path):
        try:
            import json
            obj = json.load(open(categories_path))
            cats = obj.get("categories", obj) if isinstance(obj, dict) else {}
            ordered = []
            for k in sorted(cats, key=lambda x: int(x)):
                v = cats[k]
                ordered.append(v[0] if isinstance(v, (list, tuple)) else str(v))
            for i, name in enumerate(ordered):
                mapping["kl_%d" % (i + 1)] = name
        except (ValueError, KeyError, TypeError):
            mapping = {}

    def namer(col):
        return mapping.get(col, col.replace("kl_", ""))
    return namer


def _parse_locus(s):
    """'chr6:40790000-41690000' -> (chrom, start, end), or None."""
    if not s:
        return None
    try:
        chrom, rng = s.split(":")
        a, b = rng.replace(",", "").split("-")
        return (chrom, int(a), int(b))
    except ValueError:
        raise SystemExit("ERROR: bad --locus '%s' (want chr:start-end)" % s)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--diff", required=True, help="diff_<A>_vs_<B>.stats.tsv")
    ap.add_argument("--loops", required=True, help="loop BEDPE for the condition")
    ap.add_argument("--out", required=True, help="output summary TSV (one row)")
    ap.add_argument("--categories", default=None,
                    help="optional categories JSON to label kl_N tracks by name")
    ap.add_argument("--locus", default=None,
                    help="restrict to a genomic window chr:start-end (e.g. the "
                         "Tcrb locus) instead of the genome-wide noise floor")
    ap.add_argument("--fdr-only", action="store_true",
                    help="restrict the |BES| comparison to FDR-significant bins")
    ap.add_argument("--bins-out", default=None,
                    help="optional per-bin labelled TSV (bin, stratum, BES, fdr)")
    args = ap.parse_args()

    namer = _track_namer(args.categories)
    locus = _parse_locus(args.locus)
    comp = parse_comparison_name(args.diff)
    lookup = build_anchor_lookup(read_loops(args.loops))

    fh = _open(args.diff)
    reader = csv.DictReader(fh, delimiter="\t", lineterminator="\n")
    fields = reader.fieldnames or []
    fdr_col = next((c for c in fields if c.startswith("significant_fdr")), None)
    kl_cols = [c for c in fields if c.startswith("kl_")]

    anchor = {"abs_bes": [], "sig": 0, "n": 0, "kl": {k: 0.0 for k in kl_cols}}
    nonanc = {"abs_bes": [], "sig": 0, "n": 0, "kl": {k: 0.0 for k in kl_cols}}

    bins_rows = []
    for r in reader:
        try:
            c = r["chrom"]; s = int(r["start"]); e = int(r["end"])
            bes = float(r.get("bearing_score", "nan"))
        except (KeyError, ValueError):
            continue
        sig = int(float(r.get(fdr_col, "0") or 0)) if fdr_col else 0
        if locus is not None and (c != locus[0] or e <= locus[1] or s >= locus[2]):
            continue
        if args.fdr_only and sig != 1:
            continue
        grp = anchor if overlaps_anchor(lookup, c, s, e) else nonanc
        grp["abs_bes"].append(abs(bes))
        grp["sig"] += sig
        grp["n"] += 1
        for k in kl_cols:
            try:
                grp["kl"][k] += abs(float(r.get(k, "0") or 0))
            except ValueError:
                pass
        if args.bins_out is not None:
            bins_rows.append((c, s, e,
                              "anchor" if grp is anchor else "non_anchor",
                              bes, sig))
    fh.close()

    # Statistics (scipy optional).
    mwu_p = fisher_p = "NA"
    try:
        from scipy.stats import mannwhitneyu, fisher_exact
        if anchor["abs_bes"] and nonanc["abs_bes"]:
            mwu_p = "%.4g" % mannwhitneyu(anchor["abs_bes"], nonanc["abs_bes"],
                                          alternative="two-sided").pvalue
        a_sig, a_ns = anchor["sig"], anchor["n"] - anchor["sig"]
        n_sig, n_ns = nonanc["sig"], nonanc["n"] - nonanc["sig"]
        fisher_p = "%.4g" % fisher_exact([[a_sig, a_ns], [n_sig, n_ns]])[1]
    except ImportError:
        sys.stderr.write("WARN: scipy unavailable; p-values reported as NA\n")

    def mean(xs):
        return (sum(xs) / len(xs)) if xs else 0.0

    a_dom = namer(max(anchor["kl"], key=anchor["kl"].get)) if (kl_cols and anchor["n"]) else ""

    with open(args.out, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(["comparison", "n_anchor_bins", "n_nonanchor_bins",
                    "mean_abs_bes_anchor", "mean_abs_bes_nonanchor",
                    "frac_fdr_sig_anchor", "frac_fdr_sig_nonanchor",
                    "mannwhitney_p_absbes", "fisher_p_fdrsig",
                    "dominant_track_anchor"])
        w.writerow([
            comp, anchor["n"], nonanc["n"],
            "%.4f" % mean(anchor["abs_bes"]),
            "%.4f" % mean(nonanc["abs_bes"]),
            "%.4f" % (anchor["sig"] / anchor["n"] if anchor["n"] else 0),
            "%.4f" % (nonanc["sig"] / nonanc["n"] if nonanc["n"] else 0),
            mwu_p, fisher_p, a_dom,
        ])

    if args.bins_out is not None:
        with open(args.bins_out, "w", newline="") as bf:
            bw = csv.writer(bf, delimiter="\t", lineterminator="\n")
            bw.writerow(["chrom", "start", "end", "stratum", "bearing_score",
                         "fdr_significant"])
            bw.writerows(bins_rows)

    sys.stderr.write(
        "%s: anchor n=%d (mean|BES|=%.3f, FDRsig=%d) vs non-anchor n=%d "
        "(mean|BES|=%.3f, FDRsig=%d) | MWU p=%s Fisher p=%s\n"
        % (comp, anchor["n"], mean(anchor["abs_bes"]), anchor["sig"],
           nonanc["n"], mean(nonanc["abs_bes"]), nonanc["sig"], mwu_p, fisher_p))


if __name__ == "__main__":
    main()
