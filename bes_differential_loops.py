#!/usr/bin/env python3
"""
bes_differential_loops.py  (#3 -- differential loops vs differential BES)
========================================================================
Ask whether loops that are GAINED or LOST between two conditions show a BEARING
compositional shift at their anchors -- the most direct "BEARING in topology"
question. If condition B gains a loop relative to A, do the anchor bins carry a
BEARING differential signal (and in which track -- CTCF/Cohesin expected)?

Inputs
  --loops-a / --loops-b : Mustache loop BEDPE for the two conditions (A=ref).
  --diff                : diff_<A>_vs_<B>.stats.tsv (BEARING differential).
Loop differential calling
  Two loops match if BOTH anchors overlap (reciprocal interval overlap) within
  --slop bp. A loop in B with no match in A is GAINED; a loop in A with no match
  in B is LOST; matched loops are SHARED. (Mustache merges resolutions, so
  anchors are 5-25 kb; default --slop 0 already tolerates that via interval
  overlap. --slop only pads further.)
For each loop class (gained/lost/shared) we collect the unique anchor bins and
report, from the BEARING diff table: n anchor bins, mean |BES|, FDR-significant
count, and the dominant per-track |kl|. A summary row per class is written.

ASCII-only. scipy optional (Fisher exact on gained-vs-shared FDR fraction);
NA if unavailable.
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
    """Return list of loops as ((cA,sA,eA),(cB,sB,eB))."""
    loops = []
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            try:
                loops.append(((p[0], int(p[1]), int(p[2])),
                              (p[3], int(p[4]), int(p[5]))))
            except ValueError:
                continue
    return loops


def iv_overlap(a, b, slop):
    (ca, sa, ea) = a
    (cb, sb, eb) = b
    if ca != cb:
        return False
    return (sa - slop) < eb and (ea + slop) > sb


def loops_match(l1, l2, slop):
    a1, b1 = l1
    a2, b2 = l2
    # anchors may be stored in either order; test both pairings
    direct = iv_overlap(a1, a2, slop) and iv_overlap(b1, b2, slop)
    swap = iv_overlap(a1, b2, slop) and iv_overlap(b1, a2, slop)
    return direct or swap


def classify(loops_a, loops_b, slop):
    matched_b = set()
    shared, lost = [], []
    for la in loops_a:
        hit = None
        for j, lb in enumerate(loops_b):
            if j in matched_b:
                continue
            if loops_match(la, lb, slop):
                hit = j
                break
        if hit is not None:
            matched_b.add(hit)
            shared.append(la)
        else:
            lost.append(la)
    gained = [lb for j, lb in enumerate(loops_b) if j not in matched_b]
    return {"gained": gained, "lost": lost, "shared": shared}


def anchor_intervals(loops):
    ivs = []
    for (a, b) in loops:
        ivs.append(a)
        ivs.append(b)
    return ivs


def _parse_locus(s):
    if not s:
        return None
    try:
        chrom, rng = s.split(":")
        a, b = rng.replace(",", "").split("-")
        return (chrom, int(a), int(b))
    except ValueError:
        raise SystemExit("ERROR: bad --locus '%s' (want chr:start-end)" % s)


def load_diff(path, locus=None, fdr_only=False):
    fh = _open(path)
    reader = csv.DictReader(fh, delimiter="\t")
    fields = reader.fieldnames or []
    fdr_col = next((c for c in fields if c.startswith("significant_fdr")), None)
    kl_cols = [c for c in fields if c.startswith("kl_")]
    rows = []
    for r in reader:
        try:
            r["_c"] = r["chrom"]; r["_s"] = int(r["start"]); r["_e"] = int(r["end"])
            r["_bes"] = float(r.get("bearing_score", "nan"))
            r["_sig"] = int(float(r.get(fdr_col, "0") or 0)) if fdr_col else 0
        except (KeyError, ValueError):
            continue
        if locus is not None and (r["_c"] != locus[0] or r["_e"] <= locus[1]
                                  or r["_s"] >= locus[2]):
            continue
        if fdr_only and r["_sig"] != 1:
            continue
        rows.append(r)
    fh.close()
    # index by chrom for overlap
    by = {}
    for r in rows:
        by.setdefault(r["_c"], []).append(r)
    for c in by:
        by[c].sort(key=lambda x: x["_s"])
    return by, kl_cols


def bins_in_intervals(by, kl_cols, ivs):
    import bisect
    seen = set()
    abs_bes = []
    sig = 0
    kl = {k: 0.0 for k in kl_cols}
    for (c, s, e) in ivs:
        if c not in by:
            continue
        rs = by[c]
        starts = [x["_s"] for x in rs]
        i0 = bisect.bisect_left(starts, s) - 1
        if i0 < 0:
            i0 = 0
        for r in rs[i0:]:
            if r["_s"] >= e:
                break
            if r["_e"] <= s:
                continue
            key = (r["_c"], r["_s"])
            if key in seen:
                continue
            seen.add(key)
            abs_bes.append(abs(r["_bes"]))
            sig += r["_sig"]
            for k in kl_cols:
                try:
                    kl[k] += abs(float(r.get(k, "0") or 0))
                except ValueError:
                    pass
    return abs_bes, sig, kl


def _track_namer(categories_path):
    """kl_1..6 (1-indexed) -> display name via categories JSON, else strip kl_."""
    mapping = {}
    if categories_path and os.path.exists(categories_path):
        try:
            import json
            cats = json.load(open(categories_path)).get("categories", {})
            for k, v in cats.items():
                name = v[0] if isinstance(v, (list, tuple)) else str(v)
                mapping["kl_%d" % (int(k) + 1)] = name
        except (ValueError, KeyError, TypeError):
            mapping = {}

    def namer(col):
        return mapping.get(col, col.replace("kl_", ""))
    return namer


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--loops-a", required=True, help="reference-condition loop BEDPE")
    ap.add_argument("--loops-b", required=True, help="comparison-condition loop BEDPE")
    ap.add_argument("--diff", required=True, help="diff_<A>_vs_<B>.stats.tsv")
    ap.add_argument("--slop", type=int, default=0,
                    help="extra bp padding for anchor overlap (default 0)")
    ap.add_argument("--categories", default=None,
                    help="optional categories JSON to label kl_N tracks by name")
    ap.add_argument("--locus", default=None,
                    help="restrict to a genomic window chr:start-end (e.g. Tcrb)")
    ap.add_argument("--fdr-only", action="store_true",
                    help="restrict BES tallies to FDR-significant bins")
    ap.add_argument("--out", required=True, help="output per-class summary TSV")
    args = ap.parse_args()

    namer = _track_namer(args.categories)
    comp = parse_comparison_name(args.diff)
    la = read_loops(args.loops_a)
    lb = read_loops(args.loops_b)
    classes = classify(la, lb, args.slop)
    by, kl_cols = load_diff(args.diff, locus=_parse_locus(args.locus),
                            fdr_only=args.fdr_only)

    def mean(xs):
        return (sum(xs) / len(xs)) if xs else 0.0

    stats = {}
    for cls in ("gained", "lost", "shared"):
        ivs = anchor_intervals(classes[cls])
        abs_bes, sig, kl = bins_in_intervals(by, kl_cols, ivs)
        dom = namer(max(kl, key=kl.get)) if (kl_cols and abs_bes) else ""
        stats[cls] = (len(classes[cls]), len(abs_bes), mean(abs_bes), sig, dom)

    fisher_p = "NA"
    try:
        from scipy.stats import fisher_exact
        g_n, g_bins, _, g_sig, _ = stats["gained"]
        s_n, s_bins, _, s_sig, _ = stats["shared"]
        if g_bins and s_bins:
            fisher_p = "%.4g" % fisher_exact(
                [[g_sig, g_bins - g_sig], [s_sig, s_bins - s_sig]])[1]
    except ImportError:
        sys.stderr.write("WARN: scipy unavailable; Fisher p reported as NA\n")

    with open(args.out, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["comparison", "loop_class", "n_loops", "n_anchor_bins",
                    "mean_abs_bes", "n_fdr_sig", "dominant_track",
                    "fisher_p_gained_vs_shared_fdr"])
        for cls in ("gained", "lost", "shared"):
            n_loops, n_bins, m_bes, sig, dom = stats[cls]
            w.writerow([comp, cls, n_loops, n_bins, "%.4f" % m_bes, sig, dom,
                        fisher_p if cls == "gained" else ""])

    sys.stderr.write(
        "%s: gained=%d lost=%d shared=%d loops | "
        "gained anchors mean|BES|=%.3f (FDRsig=%d, %s) | Fisher(gained vs shared)=%s\n"
        % (comp, len(classes["gained"]), len(classes["lost"]), len(classes["shared"]),
           stats["gained"][2], stats["gained"][3], stats["gained"][4], fisher_p))


if __name__ == "__main__":
    main()
