#!/usr/bin/env python3
"""
feature_changes_across_runs.py -- show whether the per-track differential change
at each CBE / AgR feature is STABLE across scoring regimes (default / qnorm /
seed / jsd / adaptive).

Reads each run's feature_track_changes.py output (cbe_track_changes.tsv for CBEs,
agr_track_changes.tsv for AgR genes), pulls the chosen track value, and builds a
(feature, comparison, track) x regime matrix. Flags cells whose SIGN flips across
regimes -- the instability an advisor will probe ("is this CTCF change real or an
artifact of how you scored it?").

  python feature_changes_across_runs.py --kind cbe \
    --query default=workflow/results \
    --query qnorm=workflow/results_qnorm \
    --query jsd=workflow/results_jsd \
    --query adaptive=workflow/results_adaptive \
    --tracks CTCF,Cohesin --out cbe_ctcf_cohesin_across_runs.tsv

PATH may be a run dir (then <dir>/regional/<file> is read or built on demand) or
a feature_track_changes.tsv directly. ASCII only.
"""
import argparse
import csv
import glob
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
EPS = 0.01  # |value| below this is treated as ~0 (neither + nor -) for sign tests

KIND = {
    "cbe": {"file": "cbe_track_changes.tsv", "mode": "point",
            "bed": os.path.join(HERE, "annotations", "cbe_mm10.bed")},
    "agr": {"file": "agr_track_changes.tsv", "mode": "region",
            "bed": os.path.join(HERE, "annotations", "AgRgenes_mm10_s.bed")},
}


def resolve(path, fname):
    return os.path.join(path, "regional", fname) if os.path.isdir(path) else path


def build(run_path, qpath, kind, categories, features):
    diffs = sorted(glob.glob(os.path.join(run_path, "pvalue", "diff_*.stats.tsv")))
    if not diffs:
        return False
    os.makedirs(os.path.dirname(qpath), exist_ok=True)
    cats = categories or os.path.join(run_path, "DN_rep1_cats.json")
    bed = features or KIND[kind]["bed"]
    cmd = [sys.executable, os.path.join(HERE, "feature_track_changes.py"),
           "--features", bed, "--diffs", *diffs,
           "--mode", KIND[kind]["mode"], "--out", qpath]
    if os.path.exists(cats):
        cmd += ["--categories", cats]
    sys.stderr.write("[build] %s -> %s\n"
                     % (os.path.basename(run_path.rstrip("/")), qpath))
    subprocess.run(cmd, check=True)
    return True


_POINT_BASE = {"comparison", "feature", "chrom", "feat_start", "feat_end",
               "bin_status", "bin_start", "bin_end", "bearing_score", "direction",
               "pval", "pval_adj_bh", "fdr_significant"}


def derive_tracks(path, kind, agg):
    """All track names present in a feature_track_changes file."""
    with open(path) as fh:
        cols = (csv.reader(fh, delimiter="\t").__next__())
    if kind == "cbe":
        return [c for c in cols if c not in _POINT_BASE]
    suf = "_" + agg
    return [c[:-len(suf)] for c in cols if c.endswith(suf)]


def load(path, tracks, kind, agg):
    """Return {(feature, comparison, track): value} for the requested tracks."""
    out = {}
    with open(path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        cols = rd.fieldnames or []
        for r in rd:
            feat, comp = r.get("feature", ""), r.get("comparison", "")
            for t in tracks:
                key = t if kind == "cbe" else "%s_%s" % (t, agg)
                if key not in cols:
                    continue
                v = r.get(key, "")
                try:
                    out[(feat, comp, t)] = float(v)
                except (ValueError, TypeError):
                    out[(feat, comp, t)] = None
    return out


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kind", choices=["cbe", "agr"], required=True)
    ap.add_argument("--query", action="append", required=True, metavar="LABEL=PATH")
    ap.add_argument("--tracks", default=None,
                    help="comma-separated tracks (default: ALL tracks in the file, "
                         "e.g. ATAC,RNAseqPos,RNAseqNeg,CTCF,Cohesin,H3K27ac)")
    ap.add_argument("--agg", choices=["sum", "mean", "peak"], default="sum",
                    help="AgR/region aggregate column to read (default sum)")
    ap.add_argument("--features", default=None,
                    help="override the default CBE/AgR bed for on-the-fly build")
    ap.add_argument("--categories", default=None)
    ap.add_argument("--no-build", action="store_true")
    ap.add_argument("--eps", type=float, default=EPS,
                    help="|value| below this counts as ~0 (neither + nor -) in the "
                         "sign-flip test (default %.2g)" % EPS)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    fname = KIND[args.kind]["file"]

    # Pass 1: resolve each run to its query file, building on demand.
    resolved = []
    for spec in args.query:
        if "=" not in spec:
            sys.exit("[ERROR] --query expects LABEL=PATH, got: %s" % spec)
        label, path = spec.split("=", 1)
        qpath = resolve(path, fname)
        if not os.path.exists(qpath):
            if args.no_build or not os.path.isdir(path):
                sys.exit("[ERROR] no %s for run '%s' at %s" % (fname, label, qpath))
            if not build(path, qpath, args.kind, args.categories, args.features):
                sys.exit("[ERROR] run '%s' has no pvalue/diff_*.stats.tsv at %s"
                         % (label, path))
        resolved.append((label, qpath))

    # Default to every track present in the file.
    if args.tracks:
        tracks = [t.strip() for t in args.tracks.split(",") if t.strip()]
    else:
        tracks = derive_tracks(resolved[0][1], args.kind, args.agg)
        sys.stderr.write("[tracks] using all: %s\n" % ", ".join(tracks))

    # Pass 2: load values.
    runs = [(label, load(qpath, tracks, args.kind, args.agg))
            for label, qpath in resolved]

    keys = []
    seen = set()
    for _lab, d in runs:
        for k in d:
            if k not in seen:
                seen.add(k)
                keys.append(k)
    keys.sort()

    labels = [lab for lab, _d in runs]
    agg_note = "" if args.kind == "cbe" else " (%s)" % args.agg
    print("=" * 84)
    print("Per-track change across run regimes -- %s features%s" % (args.kind.upper(), agg_note))
    print("  tracks: %s | regimes: %s" % (", ".join(tracks), ", ".join(labels)))
    print("-" * 84)

    n_flip = 0
    flips = []
    rows_out = []
    eps = args.eps
    for (feat, comp, track) in keys:
        vals = []
        for _lab, d in runs:
            vals.append(d.get((feat, comp, track)))
        present = [v for v in vals if v is not None]
        pos = [v for v in present if v > eps]
        neg = [v for v in present if v < -eps]
        flip = bool(pos) and bool(neg)
        # weak_side: magnitude of the smaller of the two opposing extremes, so
        # strong re-segmentation reversals rank above barely-over-floor noise.
        # Empty (NA) when the cell is not a flip.
        weak = min(max(pos), abs(min(neg))) if flip else None
        # flip_source: which regime dissents, among present above-eps regimes.
        #   <label>_only  : exactly one regime has the minority sign and the
        #                   majority (>=2) agrees -- e.g. adaptive_only is a
        #                   re-segmentation reversal; default_only often means a
        #                   stale base, not biology.
        #   mixed         : no single dissenter (even split / >1 in minority).
        #   incomplete    : a regime is missing or sub-eps, so not all regimes
        #                   contribute a sign -- attribution is undefined.
        fsource = ""
        if flip:
            signed = [(lab, 1 if v > 0 else -1)
                      for lab, v in zip(labels, vals)
                      if v is not None and abs(v) > eps]
            if len(signed) < len(runs):
                fsource = "incomplete"
            else:
                pos_l = [l for l, s in signed if s > 0]
                neg_l = [l for l, s in signed if s < 0]
                minority = pos_l if len(pos_l) <= len(neg_l) else neg_l
                majority = len(signed) - len(minority)
                fsource = ("%s_only" % minority[0]
                           if len(minority) == 1 and majority >= 2 else "mixed")
        if flip:
            n_flip += 1
            flips.append((feat, comp, track, vals, weak, fsource))
        rng = (max(present) - min(present)) if present else 0.0
        rows_out.append((feat, comp, track, vals, flip, rng, weak, fsource))

    print("  %d (feature x comparison x track) cells | %d flip sign across regimes"
          " (eps=%.3g)" % (len(keys), n_flip, eps))
    if flips:
        flips.sort(key=lambda x: x[4], reverse=True)   # strongest reversals first
        src_tally = {}
        for _f, _c, _t, _v, _w, fs in flips:
            src_tally[fs] = src_tally.get(fs, 0) + 1
        print("  sign-flipping cells, strongest first "
              "(regime values %s | weak_side | source):" % ",".join(labels))
        for feat, comp, track, vals, weak, fsource in flips[:25]:
            sv = ", ".join("%.3f" % v if v is not None else "NA" for v in vals)
            print("    %-22s %-14s %-8s [%s]  weak=%.3f  %s"
                  % (feat, comp, track, sv, weak, fsource))
        if len(flips) > 25:
            print("    ... and %d more" % (len(flips) - 25))
        print("  flip_source tally: "
              + ", ".join("%s=%d" % (k, src_tally[k]) for k in sorted(src_tally)))
    else:
        print("  No sign flips: every track change keeps a consistent direction "
              "across all\n  scoring regimes -- the changes are robust to scoring "
              "choice.")
    print("=" * 84)

    if args.out:
        with open(args.out, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["feature", "comparison", "track"] + labels
                       + ["sign_flip", "range", "weak_side", "flip_source"])
            for feat, comp, track, vals, flip, rng, weak, fsource in rows_out:
                w.writerow([feat, comp, track]
                           + ["%.4f" % v if v is not None else "" for v in vals]
                           + ["1" if flip else "0", "%.4f" % rng,
                              "%.4f" % weak if weak is not None else "", fsource])
        print("wrote %s" % args.out)


if __name__ == "__main__":
    main()
