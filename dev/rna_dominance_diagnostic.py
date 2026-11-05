#!/usr/bin/env python3
# ----------------------------------------------------------------------
# File     : dev/rna_dominance_diagnostic.py
# Version  : 1.0.0
# Date     : 2026-11-05
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
    scales before the composition is formed.
This script rescores two conditions on a few chromosomes under four
pre-specified variants and reports which track leads the strongest
differential bins and regions in each:

  production        the production scorer (no normalisation)
  merge_strands     RNA+ and RNA- summed into one RNA track before scoring
  nzq               per-track nonzero-quantile normalisation (the existing
                    bigwig_to_qcat --normalize-method nonzero-quantile)
  merge_strands_nzq both

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
                          regions): counts and shares of the lead track
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
    --chrom-sizes resources/mm10.chrom.sizes \\
    --categories ../categories/mm10_6track_panel.yaml \\
    --blacklist results_prior1_n100_l1/blacklist/merged_blacklist.bed \\
    --out results_prior1_n100_l1/table_sources/rna_dominance_DN_vs_S3T3_chr6

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

VARIANTS = ("production", "merge_strands", "nzq", "merge_strands_nzq")


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


# ---------------------------------------------------------------------------
# Variants and scoring
# ---------------------------------------------------------------------------

def merge_strands(x, ip, im):
    """Sum the two RNA strand columns into one (placed at the RNA+ position);
    the RNA- column is dropped."""
    y = np.delete(x, im, axis=1)
    pos = ip if ip < im else ip - 1
    y[:, pos] = x[:, ip] + x[:, im]
    return y


def variant_matrix(raw, variant, ip, im):
    """(signal used for P and the prior, pre-normalisation signal used for the
    low-signal mask) for one chromosome of one sample."""
    x = np.clip(raw.astype(np.float64), 0.0, None)
    if variant in ("merge_strands", "merge_strands_nzq"):
        x = merge_strands(x, ip, im)
    pre = x.copy()
    if variant in ("nzq", "merge_strands_nzq"):
        x = b2q.normalize_signal_matrix(x, method="nonzero-quantile")
    return x, pre


def score_sample(sample_data, variant, ip, im, min_signal, alpha):
    """KL scores (bins x tracks, chromosomes concatenated) and Q, with Q = the
    mean of P over all bins read, as production computes it genome-wide."""
    mats, pres, Ps = [], [], []
    for chrom in sample_data:
        x, pre = variant_matrix(sample_data[chrom][1], variant, ip, im)
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


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def build_parser():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samples", required=True, help="BEARING sample sheet (TSV)")
    ap.add_argument("--base-dir", default=None,
                    help="directory the sheet's relative BigWig paths are relative "
                         "to (default: current directory, i.e. run from workflow/)")
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
    ap.add_argument("--out", required=True, help="output prefix")
    return ap


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
    fracs = [float(x) for x in args.top_fracs.split(",") if x]

    print("reading %s (%d reps) and %s (%d reps) on %s"
          % (args.cond_a, len(sheet[args.cond_a]), args.cond_b,
             len(sheet[args.cond_b]), ",".join(chroms)), file=sys.stderr)
    data = {}
    for cond in (args.cond_a, args.cond_b):
        for sample, paths in sheet[cond]:
            data[sample] = (cond, read_sample(paths, chroms, sizes, neg, blacklist))

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

    rows, qrows, text = [], [], []
    prod_top = {}
    for variant in variants:
        vn = variant_names(names, variant, ip, im)
        nt = len(vn)
        sc = {args.cond_a: [], args.cond_b: []}
        for sample, (cond, sd) in data.items():
            s, Q = score_sample(sd, variant, ip, im, args.min_signal,
                                args.prior_strength)
            sc[cond].append(s)
            qrows.append(dict(variant=variant, sample=sample, condition=cond,
                              **{vn[i]: float(Q[i]) for i in range(nt)}))
        d = np.mean(sc[args.cond_a], axis=0) - np.mean(sc[args.cond_b], axis=0)
        l1 = np.abs(d).sum(axis=1)
        tested = np.where(l1 >= args.floor)[0]
        order = tested[np.argsort(-l1[tested], kind="stable")]
        sets = [("tested", tested)]
        for f in fracs:
            k = max(1, int(round(f * tested.size))) if tested.size else 0
            sets.append(("top%g" % f, np.sort(order[:k])))
            if variant == "production":
                prod_top[f] = np.sort(order[:k])
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
                for i in range(nt):
                    r["n_" + vn[i]] = int(c[i])
                rows.append(r)
                text.append("  %-14s %-7s n=%-7d RNA %5.1f%%  | %s"
                            % (set_name, unit, tot, 100 * rna_share if tot else float("nan"),
                               "  ".join("%s %d" % (vn[i], c[i]) for i in range(nt))))
        text.append("")

    # write outputs (union of track columns; blank where a variant lacks one)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    all_cols = []
    for r in rows:
        for k in r:
            if k not in all_cols:
                all_cols.append(k)
    with open(args.out + ".lead_tracks.tsv", "w") as fh:
        fh.write("\t".join(all_cols) + "\n")
        for r in rows:
            fh.write("\t".join(("%.4f" % r[k]) if isinstance(r.get(k), float)
                               else str(r.get(k, "")) for k in all_cols) + "\n")
    qcols = []
    for r in qrows:
        for k in r:
            if k not in qcols:
                qcols.append(k)
    with open(args.out + ".background_q.tsv", "w") as fh:
        fh.write("\t".join(qcols) + "\n")
        for r in qrows:
            fh.write("\t".join(("%.4f" % r[k]) if isinstance(r.get(k), float)
                               else str(r.get(k, "")) for k in qcols) + "\n")
    head = ["rna_dominance_diagnostic.py 1.0.0  %s vs %s  chroms %s  min_signal %g  "
            "prior %g  floor %g  gap %d  (scoring only; no p-values; Q from these "
            "chromosomes)" % (args.cond_a, args.cond_b, ",".join(chroms),
                              args.min_signal, args.prior_strength, args.floor, args.gap),
            "lead track = largest |per-track differential| (bins) or largest summed "
            "|differential| over the region (regions); prod_top<f> = the production "
            "top-f bins scored under this variant", ""]
    txt = "\n".join(head + text) + "\n"
    with open(args.out + ".summary.txt", "w") as fh:
        fh.write(txt)
    sys.stdout.write(txt)
    return rows, qrows


if __name__ == "__main__":
    main()
