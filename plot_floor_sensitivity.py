#!/usr/bin/env python3
"""
plot_floor_sensitivity.py -- Supplementary Figure S11 (differential-floor
sensitivity), with NO hard-coded data. Every plotted value is read from a
pipeline output or derived from one, and echoed into the --source-data TSV.

Three panels:
  A  What each floor can detect, in compositional terms (uses the REAL
     background Q from --q-json, not a constant).
  B  Calls are insensitive to the floor across the swept range.
  C  Why raising the floor raises the significance cutoff faster than it raises
     any single bin's own p-value -- so 0.5 is conservative, not hit-maximising.

ASCII only in source. Unicode only inside display strings. matplotlib (Agg)
only, no seaborn.
"""

import argparse
import hashlib
import json
import os
import sys

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

VERSION = "1.0"

FLOOR_COLORS = {0.1: "#9e9e9e", 0.25: "#66bb6a", 0.5: "#d32f2f",
                1.0: "#1976d2", 1.5: "#7b1fa2"}
CMP_COLORS = {"DN_vs_DP": "#1976d2", "DN_vs_EbKO": "#43a047",
              "DN_vs_ProB": "#8e24aa", "DN_vs_S3T3": "#ef6c00"}
PROD_FLOOR = 0.5


# ---------------------------------------------------------------------------
# IO
# ---------------------------------------------------------------------------

def sha256(path):
    h = hashlib.sha256()
    try:
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        return h.hexdigest()
    except OSError:
        return "NA"


def floor_tag(f):
    return str(f).replace(".", "p")


def read_sweep(path):
    """comparison -> floor -> {tested, null, sig, min_pval}."""
    out = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < len(header):
                continue
            cmp_name = p[idx["comparison"]]
            f = float(p[idx["floor"]])
            out.setdefault(cmp_name, {})[f] = {
                "tested": int(p[idx["bins_tested"]]),
                "null": int(p[idx["null_scores"]]),
                "sig": int(p[idx["bins_bh_significant"]]),
                "min_pval": float(p[idx["min_pval"]]),
            }
    return out


def read_stats(path):
    """(chrom,start) -> (bearing_score, pval) from a production stats.tsv."""
    out = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        need = ("chrom", "start", "pval", "bearing_score")
        if not all(k in idx for k in need):
            return out
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) <= idx["bearing_score"]:
                continue
            try:
                key = (p[idx["chrom"]], int(p[idx["start"]]))
                out[key] = (float(p[idx["bearing_score"]]),
                            float(p[idx["pval"]]))
            except ValueError:
                continue
    return out


def _open_maybe_gz(path):
    if path.endswith(".gz"):
        import gzip
        return gzip.open(path, "rt")
    return open(path)


def read_per_bin(per_bin_dir, cmp_name, floor):
    """(chrom,start) -> (bearing_score, pval, pval_adj_bh, significant)."""
    tag = floor_tag(floor)
    for cand in ("%s_floor%s.tsv.gz" % (cmp_name, tag),
                 "%s_floor%s.tsv" % (cmp_name, tag)):
        path = os.path.join(per_bin_dir, cand)
        if os.path.exists(path):
            break
    else:
        return {}
    out = {}
    with _open_maybe_gz(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < len(header):
                continue
            key = (p[idx["chrom"]], int(p[idx["start"]]))
            out[key] = (float(p[idx["bearing_score"]]), float(p[idx["pval"]]),
                        float(p[idx["pval_adj_bh"]]), int(p[idx["significant"]]))
    return out


def load_q(path):
    with open(path) as fh:
        d = json.load(fh)
    q_by_idx = {int(k): float(v) for k, v in d["Q"].items()}
    order = sorted(q_by_idx)
    q = np.array([q_by_idx[i] for i in order], dtype=np.float64)
    names = d.get("categories", {})
    return q, d.get("sample", "?"), names


def read_example_bins(path):
    """List of (chrom, start, label)."""
    out = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if p[0].lower() == "chrom":
                continue
            if len(p) < 3:
                continue
            out.append((p[0], int(p[1]), p[2]))
    return out


# ---------------------------------------------------------------------------
# Panel A math
# ---------------------------------------------------------------------------

def bes_scalar_q(P, q):
    """Clamped-KL BES with a single scalar background q for every state."""
    P = np.asarray(P, dtype=np.float64)
    mask = P > q
    if not np.any(mask):
        return 0.0
    return float(np.sum(P[mask] * np.log2(P[mask] / q)))


def detection_curve(q, xs):
    """y(x) = |ceiling - BES([1-x, x])| for a pure-track-A bin losing fraction
    x to track B, both tracks sharing background q. ceiling = -log2(q)."""
    ceiling = -np.log2(q)
    ys = np.array([abs(ceiling - bes_scalar_q([1.0 - x, x], q)) for x in xs])
    return ceiling, ys


def smallest_detectable_split(q, floor, xs):
    """Smallest x in xs with y(x) >= floor. Returns (x, ratio_str) or None."""
    _, ys = detection_curve(q, xs)
    hit = np.nonzero(ys >= floor)[0]
    if len(hit) == 0:
        return None
    x = xs[hit[0]]
    if x <= 0:
        return x, "100/0"
    return x, "%d/%d" % (round((1 - x) * 100), round(x * 100))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sweep-tsv", required=True)
    ap.add_argument("--per-bin-dir", required=True)
    ap.add_argument("--stats-dir", required=True)
    ap.add_argument("--q-json", required=True)
    ap.add_argument("--mechanism-comparison", default="DN_vs_DP")
    ap.add_argument("--example-bins", required=True)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--out-pdf", required=True)
    ap.add_argument("--out-png", required=True)
    ap.add_argument("--source-data", required=True)
    args = ap.parse_args()

    sweep = read_sweep(args.sweep_tsv)
    q, q_sample, q_names = load_q(args.q_json)
    mech = args.mechanism_comparison
    if mech not in sweep:
        sys.exit("ERROR: mechanism comparison %s not in sweep TSV" % mech)
    floors = sorted(sweep[mech].keys())

    src_rows = []   # (panel, series, x, y, note)

    # ----- figure -----
    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(15.0, 4.6))

    # ===================== Panel A =====================
    xs = np.arange(0.0, 0.5 + 1e-9, 0.001)
    q_bar = float(np.mean(q))
    q_min = float(np.min(q))
    q_max = float(np.max(q))
    ceiling, ys = detection_curve(q_bar, xs)
    # band across the per-track Q range (min q = most detectable)
    _, ys_lo = detection_curve(q_max, xs)   # highest q -> lowest curve
    _, ys_hi = detection_curve(q_min, xs)   # lowest q  -> highest curve
    axA.fill_between(xs, ys_lo, ys_hi, color="#cccccc", alpha=0.5, lw=0,
                     label="per-track Q range")
    axA.plot(xs, ys, color="#222222", lw=2.0,
             label="mean Q = %.4f (ceiling %.2f)" % (q_bar, ceiling))
    for f in floors:
        col = FLOOR_COLORS.get(f, "#000000")
        axA.axhline(f, ls="--", lw=1.0, color=col)
        res = smallest_detectable_split(q_bar, f, xs)
        if res is None:
            txt = "floor %g: even 50/50 undetectable" % f
            axA.annotate(txt, xy=(0.5, f), xytext=(0.30, f + 0.03),
                         fontsize=6.5, color=col, ha="left", va="bottom")
            src_rows.append(("A", "floor_marker", f, None, "not reached"))
        else:
            x, ratio = res
            axA.plot([x], [f], "o", color=col, ms=4)
            axA.annotate("floor %g -> %s" % (f, ratio), xy=(x, f),
                         xytext=(x + 0.01, f + 0.03), fontsize=6.5, color=col,
                         ha="left", va="bottom")
            src_rows.append(("A", "floor_marker", f, x,
                             "smallest detectable split %s" % ratio))
    for x, yv in zip(xs, ys):
        src_rows.append(("A", "detection_curve", float(x), float(yv), ""))
    axA.set_xlabel("fraction x moved to a second track")
    axA.set_ylabel("|delta BES|  (bits)")
    axA.set_title("A. What each floor can detect")
    axA.set_xlim(0, 0.5)
    axA.legend(fontsize=6.5, loc="upper left", frameon=False)

    # ===================== Panel B =====================
    max_dev = 0.0
    for cmp_name in sorted(sweep.keys()):
        d = sweep[cmp_name]
        if PROD_FLOOR not in d or d[PROD_FLOOR]["sig"] == 0:
            continue
        base = d[PROD_FLOOR]["sig"]
        fs = sorted(d.keys())
        ys_b = [100.0 * d[f]["sig"] / base for f in fs]
        col = CMP_COLORS.get(cmp_name, "#555555")
        axB.plot(fs, ys_b, "-o", ms=3, color=col, label=cmp_name)
        for f, yv in zip(fs, ys_b):
            src_rows.append(("B", cmp_name, f, yv, ""))
            if 0.1 <= f <= 0.5:
                max_dev = max(max_dev, abs(yv - 100.0))
    axB.axvline(PROD_FLOOR, color="#d32f2f", lw=1.0, ls=":")
    axB.axvspan(0.1, 0.5, color="#eeeeee", alpha=0.6, lw=0)
    axB.set_xscale("log")
    axB.set_xlabel("differential floor (log scale)")
    axB.set_ylabel("bins significant\n(% of floor 0.5)")
    axB.set_title("B. Calls are insensitive to the floor")
    axB.annotate("max deviation over 0.1-0.5:\n%.1f%%" % max_dev,
                 xy=(0.03, 0.03), xycoords="axes fraction", fontsize=7,
                 ha="left", va="bottom")
    axB.legend(fontsize=6.5, loc="upper left", frameon=False)
    src_rows.append(("B", "max_dev_0.1_0.5", None, max_dev, "percent"))

    # ===================== Panel C =====================
    d = sweep[mech]
    n_null = {f: d[f]["null"] for f in floors}
    pstar = {}
    for f in floors:
        t = d[f]["tested"]
        pstar[f] = (args.fdr * d[f]["sig"] / t) if t > 0 else np.nan
        src_rows.append(("C", "pstar", f, pstar[f], mech))
    fs_arr = np.array(floors)
    ps_arr = np.array([pstar[f] for f in floors])
    axC.plot(fs_arr, ps_arr, color="#000000", lw=2.5, label="BH cutoff p*(f)")
    axC.fill_between(fs_arr, 1e-12, ps_arr, color="#e0e0e0", alpha=0.6, lw=0)

    def trace_from_C(Cp1, s_abs, label, color):
        """Draw p(f) = Cp1 / (n_null_f + 1) over floors f <= |s|. Star the first
        floor where the trace dips below p*(f)."""
        tf = [f for f in floors if f <= s_abs + 1e-12]
        if not tf:
            return
        pf = [Cp1 / (n_null[f] + 1.0) for f in tf]
        axC.plot(tf, pf, "-o", ms=3, color=color, label=label)
        crossed = None
        for f, pv in zip(tf, pf):
            src_rows.append(("C", label, f, pv, "s=%.4g" % s_abs))
            if crossed is None and not np.isnan(pstar[f]) and pv <= pstar[f]:
                crossed = f
        if crossed is not None:
            axC.plot([crossed], [Cp1 / (n_null[crossed] + 1.0)], "*",
                     color=color, ms=13, mec="black", mew=0.4)

    # (1) strongest bin: C from min_pval at production floor
    n0 = n_null[PROD_FLOOR]
    Cp1_strong = d[PROD_FLOOR]["min_pval"] * (n0 + 1.0)
    trace_from_C(Cp1_strong, max(floors) + 1.0,
                 "strongest bin (C=%d)" % max(0, round(Cp1_strong - 1)),
                 "#111111")
    src_rows.append(("C", "strongest_bin_C", None, round(Cp1_strong - 1),
                     "from min_pval at floor %g" % PROD_FLOOR))

    # (2) real flipping bin: deterministic from per-bin sweep output
    pb_lo = read_per_bin(args.per_bin_dir, mech, PROD_FLOOR)
    pb_hi = read_per_bin(args.per_bin_dir, mech, max(floors))
    flip_key = None
    if pb_lo and pb_hi:
        cand = [(k, pb_lo[k][0]) for k in pb_lo
                if pb_lo[k][3] == 0 and k in pb_hi and pb_hi[k][3] == 1]
        if cand:
            cand.sort(key=lambda kv: abs(kv[1]))
            # median by |bearing_score|, deterministic
            flip_key = cand[len(cand) // 2][0]
    if flip_key is not None:
        score, pval_lo = pb_lo[flip_key][0], pb_lo[flip_key][1]
        Cp1_flip = pval_lo * (n0 + 1.0)
        trace_from_C(Cp1_flip, abs(score),
                     "flipping bin %s:%d" % flip_key, "#00838f")
        src_rows.append(("C", "flipping_bin", None, None,
                         "%s:%d score=%.4g" % (flip_key[0], flip_key[1], score)))
    else:
        src_rows.append(("C", "flipping_bin", None, None, "none found"))

    # (3) example bins from production stats.tsv
    stats_path = os.path.join(args.stats_dir, "diff_%s.stats.tsv" % mech)
    stats = read_stats(stats_path)
    ex_bins = read_example_bins(args.example_bins)
    ex_colors = ["#8e24aa", "#ef6c00", "#43a047"]
    for i, (chrom, start, label) in enumerate(ex_bins):
        rec = stats.get((chrom, start))
        if rec is None:
            src_rows.append(("C", "example", None, None,
                             "%s:%d NOT in stats" % (chrom, start)))
            continue
        score, pval = rec
        Cp1 = pval * (n0 + 1.0)
        trace_from_C(Cp1, abs(score), "ex: %s" % label,
                     ex_colors[i % len(ex_colors)])

    # rate annotation: n_null shrink vs p* rise, floor 0.5 -> max floor
    fmax = max(floors)
    null_ratio = (n_null[PROD_FLOOR] / n_null[fmax]) if n_null[fmax] else np.nan
    pstar_ratio = (pstar[fmax] / pstar[PROD_FLOOR]) if pstar[PROD_FLOOR] else np.nan
    axC.annotate(
        "%g -> %g:\nnull x%.1f smaller\np* x%.1f larger" % (
            PROD_FLOOR, fmax, null_ratio, pstar_ratio),
        xy=(0.97, 0.03), xycoords="axes fraction", fontsize=7,
        ha="right", va="bottom")
    src_rows.append(("C", "null_ratio_0.5_to_max", None, null_ratio, ""))
    src_rows.append(("C", "pstar_ratio_0.5_to_max", None, pstar_ratio, ""))
    axC.set_xscale("log")
    axC.set_yscale("log")
    axC.set_xlabel("differential floor (log scale)")
    axC.set_ylabel("p-value (log scale)")
    axC.set_title("C. Cutoff moves faster than the bin")
    axC.legend(fontsize=6.0, loc="lower left", frameon=False)

    fig.tight_layout()
    # CreationDate=None keeps the PDF byte-identical across runs (no timestamp).
    fig.savefig(args.out_pdf, metadata={"CreationDate": None})
    fig.savefig(args.out_png, dpi=150, metadata={"Software": None})
    plt.close(fig)

    # ----- source data -----
    inputs = {
        "sweep_tsv": args.sweep_tsv,
        "q_json": args.q_json,
        "example_bins": args.example_bins,
        "stats": stats_path,
    }
    with open(args.source_data, "w") as fh:
        fh.write("# SuppFig_S11 floor sensitivity source data\n")
        fh.write("# script_version=%s\n" % VERSION)
        fh.write("# numpy=%s matplotlib=%s\n" % (np.__version__, matplotlib.__version__))
        fh.write("# mechanism_comparison=%s fdr=%g\n" % (mech, args.fdr))
        for k, p in inputs.items():
            fh.write("# input %s=%s sha256=%s\n" % (k, p, sha256(p)))
        fh.write("# Q_source=%s Q=%s mean=%.6g min=%.6g max=%.6g\n" % (
            q_sample, ",".join("%.6g" % v for v in q), q_bar, q_min, q_max))
        if flip_key is not None:
            fh.write("# flipping_bin=%s:%d\n" % flip_key)
        fh.write("panel\tseries\tx\ty\tnote\n")
        for panel, series, x, y, note in src_rows:
            fh.write("%s\t%s\t%s\t%s\t%s\n" % (
                panel, series,
                "" if x is None else ("%.6g" % x),
                "" if y is None else ("%.6g" % y),
                note))

    print("Wrote %s, %s, %s" % (args.out_pdf, args.out_png, args.source_data))
    return 0


if __name__ == "__main__":
    sys.exit(main())
