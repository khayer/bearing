#!/usr/bin/env python3
# ----------------------------------------------------------------------
# File     : dev/rna_dominance_diagnostic.py
# Version  : 1.2.0
# Date     : 2026-11-08 (1.2.0: per-track tail statistics of the raw signal,
#            PREFIX.tail_stats.tsv, and --tail-only; scoring unchanged.
#            1.1.0 2026-11-07: data-dependent transform variants log1p_rna
#            and asinh_auto, known-region check, lead-share of top bins;
#            1.0.0 2026-11-05)
# Authors  : Katharina E. Hayer (katharinaehayer@gmail.com) and Claude
#            (Anthropic), co-created
# ----------------------------------------------------------------------
"""
rna_dominance_diagnostic.py
===========================
Scoring-only sensitivity check: how much of the RNA dominance of BEARING's
differential calls comes from HOW the panel is scored?

In the l1 production run 83-100% of the genome-wide BH regions are led by an
RNA strand and no region is led by CTCF. Candidate scoring causes:
  * the strand split: one transcribed strand of a gene is ~100% of a bin's
    signal, the most extreme composition BEARING can score;
  * dynamic range: RNA (CPM) and ChIP / CUT&RUN tracks are on different
    scales before the composition is formed (RNA is heavy-tailed: a few
    exonic bins carry values orders of magnitude above the ChIP tracks).
This script rescores two conditions on a few chromosomes under pre-specified
variants and reports which track leads the strongest differential bins and
regions in each:

  production        the production scorer (no normalisation)
  merge_strands     RNA+ and RNA- summed into one RNA track before scoring
  nzq               per-track nonzero-quantile normalisation (the existing
                    bigwig_to_qcat --normalize-method nonzero-quantile)
  merge_strands_nzq both
  log1p_rna         (1.1.0) RNA tracks only: x -> log1p(x / c_i)
  asinh_auto        (1.1.0) every track: x -> asinh(x / c_i)

c_i is DATA-DEPENDENT: the median of the nonzero raw values of track i, pooled
over all samples of both conditions and all chromosomes read (one cofactor per
track, shared by the samples, so real between-condition differences in a track
are kept). Values well below c_i stay nearly linear; values far above c_i are
compressed logarithmically. Both transforms map 0 to 0, so the low-signal mask
(applied, as in production, to the PRE-transform raw signal) and on/off bins are
unchanged; only the relative weight of large values changes. The cofactors are
written to PREFIX.cofactors.tsv.

Scoring uses the production functions of bigwig_to_qcat.py
(mean_signal_in_bins, negative-strand abs, blacklist zeroing,
normalize_signal_matrix, signals_to_prob, kl_scores_per_bin with the low-
signal mask on the PRE-normalisation signal and the composition prior on the
post-normalisation signal), the differential is mean(A reps) - mean(B reps)
per track (compare_qcat), the bin statistic is l1 = sum_i |d_i| and bins
with l1 >= --floor are "tested", as in bearing_pvalue --diff-stat l1.

What this is NOT: no permutation null, no p-values. The "strongest" bins are
the top --top-fracs of the tested bins by l1, a proxy for what a test would
call. Q (the background composition) is computed from the chromosomes read
(--chroms), not genome-wide; use more chromosomes for a closer match.

OUTPUT (--out PREFIX)
  PREFIX.lead_tracks.tsv  one row per variant x bin set x unit (bins or
                          regions): counts and shares of the lead track, and
                          (bins) the median lead share = |d_lead| / l1
  PREFIX.regions.tsv      (1.1.0) per variant x --check-regions region: tested
                          bins, bins in the top fractions, best l1 rank
                          percentile among tested bins, lead track
  PREFIX.cofactors.tsv    (1.1.0) c_i per track
  PREFIX.tail_stats.tsv   (1.2.0) per track x group (each condition and
                          both pooled; all replicates of a group together):
                          bins read, fraction of bins > 0, quantiles q50 /
                          q90 / q99 / q99.9 and max of the nonzero raw values
                          (after abs of negative-strand tracks and blacklist
                          zeroing, before any transform), q99/q50 and
                          q99.9/q50 (the dynamic range of the track). With
                          --tail-only only this file and the summary are
                          written (no scoring).
  PREFIX.background_q.tsv background composition Q per variant and sample
  PREFIX.summary.txt      human-readable

Bin sets per variant: all tested bins; top fractions of tested bins by l1;
and "prod_top<f>": the bins that are in the production top fraction, scored
under the variant (same bins, so the change in lead track is directly due to
the scoring). Regions: bins of a set merged when < --gap bp apart; lead
track = largest summed |d_i| over the region.

Usage (from workflow/, cluster env: module purge; conda activate bearing;
env -u PYTHONPATH):
  python ../dev/rna_dominance_diagnostic.py \\
    --samples config/samples.tsv --cond-a DN --cond-b S3T3 --chroms chr6 \\
    --base-dir /mnt/isilon/bassing_lab/integration_paper \\
    --chrom-sizes resources/mm10.chrom.sizes \\
    --categories ../categories/mm10_6track_panel.yaml \\
    --blacklist results_prior1_n100_l1/blacklist/merged_blacklist.bed \\
    --variants production,nzq,log1p_rna,asinh_auto \\
    --out results_prior1_n100_l1/table_sources/rna_dominance_v11_DN_vs_S3T3_chr6

ASCII only. Reads BigWigs; writes nothing else.
"""

import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)

import bigwig_to_qcat as b2q  # noqa: E402  production scorer

VERSION = "1.2.0"
VARIANTS = ("production", "merge_strands", "nzq", "merge_strands_nzq",
            "log1p_rna", "asinh_auto")
DEFAULT_REGIONS = ("DJ_RC=chr6:41500000-41551500,"
                   "Vbeta_cluster=chr6:41040000-41290000,"
                   "Trbv1=chr6:40880000-40905000,"
                   "Cd4=chr6:124825000-124950000")


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

def read_sample_sheet(path, base_dir):
    """{condition: [(sample, [bw paths in panel order]), ...]} from the
    BEARING sample sheet (rows with a non-empty bw column)."""
    out = {}
    header = None
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if header is None:
                header = f
                continue
            row = dict(zip(header, f + [""] * (len(header) - len(f))))
            bw = row.get("bw", "").strip()
            if not bw or not row.get("sample", "").strip():
                continue
            paths = [p if os.path.isabs(p) else os.path.normpath(os.path.join(base_dir, p))
                     for p in bw.split(",")]
            out.setdefault(row["condition"], []).append((row["sample"], paths))
    return out


def read_chrom_sizes(path):
    sizes = {}
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) >= 2:
                sizes[f[0]] = int(f[1])
    return sizes


def read_track_names(categories_path, n):
    names, neg = None, set()
    if categories_path:
        cats, neg = b2q.load_categories_yaml(categories_path)
        names = []
        for c in cats:
            if isinstance(c, dict):
                names.append(str(c.get("name", "")))
            elif isinstance(c, (list, tuple)) and c:
                names.append(str(c[0]))
            else:
                names.append(str(c))
    if not names or len(names) != n:
        names = ["track%d" % (i + 1) for i in range(n)]
    return names, set(neg)


def find_rna_strands(names, override=None):
    """0-based indices of RNA+ and RNA- from the track names, or --rna-cols."""
    if override:
        a, b = [int(x) - 1 for x in override.split(",")]
        return a, b
    plus = [i for i, n in enumerate(names) if "rna" in n.lower() and "+" in n]
    minus = [i for i, n in enumerate(names) if "rna" in n.lower() and "-" in n]
    if len(plus) != 1 or len(minus) != 1:
        sys.exit("ERROR: could not identify the two RNA strand tracks from the "
                 "category names %s; pass --rna-cols PLUS,MINUS (1-based)" % names)
    return plus[0], minus[0]


def read_sample(paths, chroms, sizes, neg_states, blacklist):
    """Per chromosome raw signal matrix (bins x tracks), as _compute_chrom_cache
    builds it before normalisation: mean per 200-bp bin, abs() of negative-
    strand tracks, blacklisted bins zeroed. Returns {chrom: (bins, matrix)}."""
    import pyBigWig
    out = {}
    for chrom in chroms:
        bins = b2q.bins_for_chrom(sizes[chrom])
        m = np.zeros((len(bins), len(paths)), dtype=np.float32)
        for si, p in enumerate(paths):
            with pyBigWig.open(p) as bw:
                m[:, si] = b2q.mean_signal_in_bins(bw, chrom, bins)
        for neg_idx in neg_states:
            if neg_idx - 1 < m.shape[1]:
                m[:, neg_idx - 1] = np.abs(m[:, neg_idx - 1])
        if blacklist is not None:
            m[b2q.bins_overlapping_blacklist(bins, blacklist, chrom), :] = 0.0
        out[chrom] = (bins, m)
    return out


def parse_regions(spec):
    """'name=chr:start-end,...' -> [(name, chrom, start, end)]"""
    out = []
    for item in [s for s in (spec or "").split(",") if s.strip()]:
        name, loc = item.split("=")
        chrom, rng = loc.split(":")
        a, b = rng.replace("_", "").split("-")
        out.append((name.strip(), chrom.strip(), int(a), int(b)))
    return out


# ---------------------------------------------------------------------------
# Variants and scoring
# ---------------------------------------------------------------------------

def track_cofactors(data):
    """Median of the nonzero raw values of each track, pooled over every
    sample and chromosome read (data = {sample: (cond, {chrom: (bins, m)})}).
    A track with no nonzero value gets cofactor 1."""
    pooled = None
    for _, (_, sd) in data.items():
        for chrom in sd:
            m = np.clip(sd[chrom][1].astype(np.float64), 0.0, None)
            if pooled is None:
                pooled = [[] for _ in range(m.shape[1])]
            for i in range(m.shape[1]):
                v = m[:, i]
                pooled[i].append(v[v > 0])
    cof = []
    for parts in pooled:
        v = np.concatenate(parts) if parts else np.array([])
        cof.append(float(np.median(v)) if v.size else 1.0)
    return np.array(cof, dtype=np.float64)


TAIL_QS = (0.5, 0.9, 0.99, 0.999)


def tail_stats(data, names, cond_a, cond_b):
    """Per track and group (cond_a, cond_b, pooled): distribution of the
    nonzero raw values (data = {sample: (cond, {chrom: (bins, m)})}).
    Returns a list of row dicts."""
    groups = ((cond_a, (cond_a,)), (cond_b, (cond_b,)), ("pooled", (cond_a, cond_b)))
    rows = []
    for i, name in enumerate(names):
        for gname, conds in groups:
            parts, n_all = [], 0
            for _, (cond, sd) in data.items():
                if cond not in conds:
                    continue
                for chrom in sd:
                    v = np.clip(sd[chrom][1][:, i].astype(np.float64), 0.0, None)
                    n_all += v.size
                    parts.append(v[v > 0])
            v = np.concatenate(parts) if parts else np.array([])
            r = dict(track=name, group=gname, bins=int(n_all),
                     frac_nonzero=(v.size / n_all) if n_all else float("nan"))
            if v.size:
                qs = np.quantile(v, TAIL_QS)
                for q, x in zip(TAIL_QS, qs):
                    r["q%g" % (100 * q)] = float(x)
                r["max"] = float(v.max())
                r["q99_over_q50"] = float(qs[2] / qs[0]) if qs[0] > 0 else float("nan")
                r["q99.9_over_q50"] = float(qs[3] / qs[0]) if qs[0] > 0 else float("nan")
            else:
                for q in TAIL_QS:
                    r["q%g" % (100 * q)] = float("nan")
                r["max"] = r["q99_over_q50"] = r["q99.9_over_q50"] = float("nan")
            rows.append(r)
    return rows


def tail_text(rows):
    out = ["raw-signal tails (nonzero values; before any transform):",
           "  %-12s %-8s %8s %10s %10s %10s %10s %8s %9s"
           % ("track", "group", "nonzero", "q50", "q90", "q99", "q99.9",
              "q99/q50", "q999/q50")]
    for r in rows:
        out.append("  %-12s %-8s %7.1f%% %10.4g %10.4g %10.4g %10.4g %8.3g %9.3g"
                   % (r["track"], r["group"], 100 * r["frac_nonzero"], r["q50"],
                      r["q90"], r["q99"], r["q99.9"], r["q99_over_q50"],
                      r["q99.9_over_q50"]))
    return out


def merge_strands(x, ip, im):
    """Sum the two RNA strand columns into one (placed at the RNA+ position);
    the RNA- column is dropped."""
    y = np.delete(x, im, axis=1)
    pos = ip if ip < im else ip - 1
    y[:, pos] = x[:, ip] + x[:, im]
    return y


def variant_matrix(raw, variant, ip, im, cof=None):
    """(signal used for P and the prior, pre-normalisation signal used for the
    low-signal mask) for one chromosome of one sample. For every variant the
    mask signal is the untransformed (clipped, strand-merged where relevant)
    raw signal, so the masked bins are the same as in production."""
    x = np.clip(raw.astype(np.float64), 0.0, None)
    if variant in ("merge_strands", "merge_strands_nzq"):
        x = merge_strands(x, ip, im)
    pre = x.copy()
    if variant in ("nzq", "merge_strands_nzq"):
        x = b2q.normalize_signal_matrix(x, method="nonzero-quantile")
    elif variant == "log1p_rna":
        for j in (ip, im):
            x[:, j] = np.log1p(x[:, j] / cof[j])
    elif variant == "asinh_auto":
        x = np.arcsinh(x / cof[np.newaxis, :])
    return x, pre


def score_sample(sample_data, variant, ip, im, min_signal, alpha, cof=None):
    """KL scores (bins x tracks, chromosomes concatenated) and Q, with Q = the
    mean of P over all bins read, as production computes it genome-wide."""
    mats, pres, Ps = [], [], []
    for chrom in sample_data:
        x, pre = variant_matrix(sample_data[chrom][1], variant, ip, im, cof)
        mats.append(x)
        pres.append(pre)
        Ps.append(b2q.signals_to_prob(x))
    P = np.concatenate(Ps)
    X = np.concatenate(mats)
    PRE = np.concatenate(pres)
    Q = P.mean(axis=0)
    s, _ = b2q.kl_scores_per_bin(P, Q, raw_signal_matrix=PRE, min_signal=min_signal,
                                 prior_strength=alpha, prior_source=X)
    return s.astype(np.float32), Q


def variant_names(names, variant, ip, im):
    if variant in ("merge_strands", "merge_strands_nzq"):
        out = list(names)
        out[ip] = "RNA (both strands)"
        del out[im]
        return out
    return list(names)


# ---------------------------------------------------------------------------
# Lead-track summaries
# ---------------------------------------------------------------------------

def lead_counts_bins(d, idx, n_tracks):
    if idx.size == 0:
        return np.zeros(n_tracks, dtype=int)
    lead = np.argmax(np.abs(d[idx]), axis=1)
    return np.bincount(lead, minlength=n_tracks)


def median_lead_share(d, idx):
    """Median over bins of |d_lead| / l1: 1.0 = a single track carries the
    whole differential (transform-invariant one-track bins)."""
    if idx.size == 0:
        return float("nan")
    a = np.abs(d[idx])
    l1 = a.sum(axis=1)
    ok = l1 > 0
    return float(np.median(a[ok].max(axis=1) / l1[ok])) if ok.any() else float("nan")


def regions_from_bins(idx, chrom_of, start_of, end_of, gap):
    """Merge sorted bin indices into regions (same chromosome, start - previous
    end <= gap). Returns a list of index arrays."""
    if idx.size == 0:
        return []
    idx = np.sort(idx)
    regs, cur = [], [idx[0]]
    for a, b in zip(idx[:-1], idx[1:]):
        if chrom_of[b] == chrom_of[a] and start_of[b] - end_of[a] <= gap:
            cur.append(b)
        else:
            regs.append(np.array(cur))
            cur = [b]
    regs.append(np.array(cur))
    return regs


def lead_counts_regions(d, regs, n_tracks):
    c = np.zeros(n_tracks, dtype=int)
    for r in regs:
        c[int(np.argmax(np.abs(d[r]).sum(axis=0)))] += 1
    return c


def region_check(d, l1, tested, order, top_sets, regions, chroms, chrom_of,
                 start_of, end_of, names):
    """Per known region: tested bins, bins in each top set, best rank
    percentile (0 = strongest tested bin) and the lead track."""
    rank = np.full(l1.size, -1, dtype=np.int64)
    rank[order] = np.arange(order.size)
    tested_mask = np.zeros(l1.size, dtype=bool)
    tested_mask[tested] = True
    out = []
    for name, chrom, a, b in regions:
        if chrom not in chroms:
            continue
        ci = chroms.index(chrom)
        sel = np.where((chrom_of == ci) & (end_of > a) & (start_of < b))[0]
        st = sel[tested_mask[sel]]
        r = dict(region=name, tested=int(st.size))
        for set_name, idx in top_sets:
            r["in_" + set_name] = int(np.isin(st, idx).sum())
        if st.size and order.size:
            r["best_rank_pct"] = 100.0 * rank[st].min() / order.size
            r["lead"] = names[int(np.argmax(np.abs(d[st]).sum(axis=0)))]
        else:
            r["best_rank_pct"] = float("nan")
            r["lead"] = ""
        out.append(r)
    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def build_parser():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samples", required=True, help="BEARING sample sheet (TSV)")
    ap.add_argument("--base-dir", default=None,
                    help="directory the sheet's relative BigWig paths resolve "
                         "against: the config's data_dir, as in the Snakefile "
                         "(default: current directory)")
    ap.add_argument("--cond-a", required=True)
    ap.add_argument("--cond-b", required=True)
    ap.add_argument("--chroms", default="chr6", help="comma-separated")
    ap.add_argument("--chrom-sizes", required=True)
    ap.add_argument("--categories", default=None, help="panel YAML (track names, "
                    "negative-strand tracks)")
    ap.add_argument("--rna-cols", default=None,
                    help="1-based PLUS,MINUS RNA columns if names do not show them")
    ap.add_argument("--blacklist", default=None, help="merged blacklist BED")
    ap.add_argument("--min-signal", type=float, default=0.1)
    ap.add_argument("--prior-strength", type=float, default=1.0)
    ap.add_argument("--floor", type=float, default=0.5,
                    help="testing floor on l1 (pvalue_min_signal)")
    ap.add_argument("--top-fracs", default="0.01,0.001")
    ap.add_argument("--gap", type=int, default=1000)
    ap.add_argument("--variants", default=",".join(VARIANTS))
    ap.add_argument("--check-regions", default=DEFAULT_REGIONS,
                    help="name=chr:start-end,... known regions to track across "
                         "variants (default: Tcrb DJ_RC, Vbeta cluster, Trbv1, Cd4)")
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--tail-only", action="store_true",
                    help="(1.2.0) write only the raw-signal tail statistics "
                         "(PREFIX.tail_stats.tsv + summary); no scoring")
    return ap


def _write_tsv(path, rows):
    cols = []
    for r in rows:
        for k in r:
            if k not in cols:
                cols.append(k)
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(("%.4f" % r[k]) if isinstance(r.get(k), float)
                               else str(r.get(k, "")) for k in cols) + "\n")


def main(argv=None):
    args = build_parser().parse_args(argv)
    base = args.base_dir or os.getcwd()
    sheet = read_sample_sheet(args.samples, base)
    for c in (args.cond_a, args.cond_b):
        if c not in sheet:
            sys.exit("ERROR: condition %s not in %s (have %s)"
                     % (c, args.samples, sorted(sheet)))
    chroms = [c for c in args.chroms.split(",") if c]
    sizes = read_chrom_sizes(args.chrom_sizes)
    missing = [c for c in chroms if c not in sizes]
    if missing:
        sys.exit("ERROR: chromosomes not in %s: %s" % (args.chrom_sizes, missing))
    n_tracks = len(sheet[args.cond_a][0][1])
    names, neg = read_track_names(args.categories, n_tracks)
    ip, im = find_rna_strands(names, args.rna_cols)
    blacklist = b2q.load_blacklist(args.blacklist) if args.blacklist else None
    variants = [v for v in args.variants.split(",") if v]
    bad = [v for v in variants if v not in VARIANTS]
    if bad or "production" not in variants:
        sys.exit("ERROR: --variants must include production and only %s" % (VARIANTS,))
    variants = ["production"] + [v for v in variants if v != "production"]
    fracs = [float(x) for x in args.top_fracs.split(",") if x]
    regions = parse_regions(args.check_regions)

    print("reading %s (%d reps) and %s (%d reps) on %s"
          % (args.cond_a, len(sheet[args.cond_a]), args.cond_b,
             len(sheet[args.cond_b]), ",".join(chroms)), file=sys.stderr)
    data = {}
    for cond in (args.cond_a, args.cond_b):
        for sample, paths in sheet[cond]:
            data[sample] = (cond, read_sample(paths, chroms, sizes, neg, blacklist))
    cof = track_cofactors(data)
    trows = tail_stats(data, names, args.cond_a, args.cond_b)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    _write_tsv(args.out + ".tail_stats.tsv", trows)
    if args.tail_only:
        txt = "\n".join(["rna_dominance_diagnostic.py %s  %s vs %s  chroms %s  "
                         "(--tail-only: raw-signal tails, no scoring)"
                         % (VERSION, args.cond_a, args.cond_b, ",".join(chroms)), ""]
                        + tail_text(trows)) + "\n"
        with open(args.out + ".summary.txt", "w") as fh:
            fh.write(txt)
        sys.stdout.write(txt)
        return [], [], []

    # bin coordinates (identical for every sample)
    any_sample = next(iter(data.values()))[1]
    chrom_of, start_of, end_of = [], [], []
    for ci, chrom in enumerate(chroms):
        bins = any_sample[chrom][0]
        chrom_of.append(np.full(len(bins), ci, dtype=np.int32))
        start_of.append(np.array([b[0] for b in bins], dtype=np.int64))
        end_of.append(np.array([b[1] for b in bins], dtype=np.int64))
    chrom_of = np.concatenate(chrom_of)
    start_of = np.concatenate(start_of)
    end_of = np.concatenate(end_of)

    rows, qrows, rrows, text = [], [], [], []
    text.append("cofactors c_i (median nonzero raw value, pooled): "
                + "  ".join("%s %.4g" % (names[i], cof[i]) for i in range(len(names))))
    text.append("")
    text.extend(tail_text(trows))
    text.append("")
    prod_top = {}
    for variant in variants:
        vn = variant_names(names, variant, ip, im)
        nt = len(vn)
        sc = {args.cond_a: [], args.cond_b: []}
        for sample, (cond, sd) in data.items():
            s, Q = score_sample(sd, variant, ip, im, args.min_signal,
                                args.prior_strength, cof)
            sc[cond].append(s)
            qrows.append(dict(variant=variant, sample=sample, condition=cond,
                              **{vn[i]: float(Q[i]) for i in range(nt)}))
        d = np.mean(sc[args.cond_a], axis=0) - np.mean(sc[args.cond_b], axis=0)
        l1 = np.abs(d).sum(axis=1)
        tested = np.where(l1 >= args.floor)[0]
        order = tested[np.argsort(-l1[tested], kind="stable")]
        sets = [("tested", tested)]
        top_sets = []
        for f in fracs:
            k = max(1, int(round(f * tested.size))) if tested.size else 0
            top = np.sort(order[:k])
            sets.append(("top%g" % f, top))
            top_sets.append(("top%g" % f, top))
            if variant == "production":
                prod_top[f] = top
            else:
                sets.append(("prod_top%g" % f, prod_top[f]))
        text.append("== %s  (%d tracks: %s)  tested bins (l1 >= %g): %d"
                    % (variant, nt, ", ".join(vn), args.floor, tested.size))
        for set_name, idx in sets:
            for unit in ("bins", "regions"):
                if unit == "bins":
                    c = lead_counts_bins(d, idx, nt)
                else:
                    c = lead_counts_regions(d, regions_from_bins(
                        idx, chrom_of, start_of, end_of, args.gap), nt)
                tot = int(c.sum())
                rna = [i for i, n in enumerate(vn) if "rna" in n.lower()]
                rna_share = (c[rna].sum() / tot) if tot else float("nan")
                r = dict(variant=variant, set=set_name, unit=unit, n=tot,
                         rna_share=rna_share)
                if unit == "bins":
                    r["median_lead_share"] = median_lead_share(d, idx)
                for i in range(nt):
                    r["n_" + vn[i]] = int(c[i])
                rows.append(r)
                extra = ("  lead share %.2f" % r["median_lead_share"]) if unit == "bins" else ""
                text.append("  %-14s %-7s n=%-7d RNA %5.1f%%%s  | %s"
                            % (set_name, unit, tot, 100 * rna_share if tot else float("nan"),
                               extra,
                               "  ".join("%s %d" % (vn[i], c[i]) for i in range(nt))))
        for r in region_check(d, l1, tested, order, top_sets, regions, chroms,
                              chrom_of, start_of, end_of, vn):
            r = dict(variant=variant, **r)
            rrows.append(r)
            text.append("  region %-14s tested %-5d %s best rank %.2f%%  lead %s"
                        % (r["region"], r["tested"],
                           "  ".join("%s %d" % (k, v) for k, v in r.items() if k.startswith("in_")),
                           r["best_rank_pct"], r["lead"]))
        text.append("")

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    _write_tsv(args.out + ".lead_tracks.tsv", rows)
    _write_tsv(args.out + ".background_q.tsv", qrows)
    _write_tsv(args.out + ".regions.tsv", rrows)
    _write_tsv(args.out + ".cofactors.tsv",
               [dict(track=names[i], cofactor=float(cof[i])) for i in range(len(names))])
    head = ["rna_dominance_diagnostic.py %s  %s vs %s  chroms %s  min_signal %g  "
            "prior %g  floor %g  gap %d  (scoring only; no p-values; Q from these "
            "chromosomes)" % (VERSION, args.cond_a, args.cond_b, ",".join(chroms),
                              args.min_signal, args.prior_strength, args.floor, args.gap),
            "lead track = largest |per-track differential| (bins) or largest summed "
            "|differential| over the region (regions); prod_top<f> = the production "
            "top-f bins scored under this variant; lead share = median |d_lead|/l1 "
            "(1 = one track carries the whole differential)", ""]
    txt = "\n".join(head + text) + "\n"
    with open(args.out + ".summary.txt", "w") as fh:
        fh.write(txt)
    sys.stdout.write(txt)
    return rows, qrows, rrows


if __name__ == "__main__":
    main()
