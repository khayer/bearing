#!/usr/bin/env python3
"""
top_diff_sites.py
=================
Find the highest BEARING differential sites genome-wide for each comparison.

IMPORTANT - ranking and the compositional ceiling:
  bearing_score is a KL divergence between NORMALIZED state distributions, so
  it discards absolute signal magnitude. A near-empty bin (e.g. a pseudogene
  with one stray antisense read) has a fully concentrated composition and pins
  the compositional ceiling, while a real multi-feature locus scores lower.
  Ranking by raw score therefore surfaces low-signal artifacts. Use the levers
  below to rank by something that respects real signal:

  --require-fdr-sig         keep only significant_fdr0.05 == 1 (cheap, no extra
                            inputs).
  --exclude-track NAME      drop sites whose dominant track is NAME (repeat);
                            e.g. --exclude-track "RNAseq -" to kill antisense
                            spikes.
  --rank-by pval            rank by adjusted p-value instead of score.
  --sheet + --reference + --min-signal FLOOR
                            absolute-signal floor read from the source bigwigs.
                            For each comparison, signal is summed across the
                            assay tracks, aggregated across ALL replicates of
                            the reference and of the comparison condition
                            (--replicate-aggregate, default max = keep a bin if
                            any replicate has signal), then combined across the
                            two conditions (--floor-combiner, default max =
                            keep on/off events present in only one condition).
                            Bins below FLOOR are dropped. This is the way to use
                            it now, since qcat files carry no raw: field.
  --rank-by signal_weighted rank by |score| * combined-signal instead of a hard
                            floor.
  --per-track               report a top-N for EACH dominant assay separately
                            instead of one global list, so CTCF / Cohesin / ATAC
                            / H3K27ac get their own rankings and are not crowded
                            out by expression (RNAseq). Candidates are chosen
                            within each track, and the floor still applies.

The legacy --qcat-a/--qcat-b path still works if your qcat files happen to
carry raw: fields (most do not).

Bins within --merge-distance are merged into one site. Output: coordinates,
peak signed score, dominant track (argmax |kl|), per-condition signal,
nearest gene (with --gtf), and any pval / significant_fdr columns found.

Usage (top-N per assay, no harsh global floor -- recommended):
  top_diff_sites.py \\
    --diff EbKO=results_v6_paper/diff_DN_vs_EbKO.stats.tsv \\
    --diff ProB=results_v6_paper/diff_DN_vs_ProB.stats.tsv \\
    --diff DP=results_v6_paper/diff_DN_vs_DP.stats.tsv \\
    --diff 3T3=results_v6_paper/diff_DN_vs_S3T3.stats.tsv \\
    --sheet samples_calib.tsv --reference DN \\
    --per-track --rank-by signal_weighted \\
    --signal-cache ./signal_cache \\
    --categories DN_rep1_cats.json \\
    --gtf ../bearing_score/gencode.vM23.annotation_modified_overlaps_removed_sorted.gtf \\
    --n 10 --out top_diff_sites_per_assay.tsv

Per-assay, no harsh floor: --per-track gives a top-N for each assay; ranking
within each assay by signal_weighted (|score| * combined signal) floats
signal-backed sites to the top and sinks low-signal ceiling junk, with no
absolute cutoff. If you prefer a hard cut, use a RELATIVE per-assay floor:
--min-signal-pctile 50 keeps the upper-signal half within each assay (scale-
appropriate), unlike one absolute --min-signal across all assays.

Speed: signal is read only for the top --signal-candidates highest-score bins
per assay, so runtime does not scale with the genome. With --signal-cache the
candidate set is independent of the floor, so sweeping is near-instant after
the first run.
"""

import argparse
import csv
import gzip
import hashlib
import json
import os
import re

import numpy as np
import pandas as pd

AUTOSOMES = {"chr{}".format(i) for i in range(1, 23)}


def load_categories_json(path):
    with open(path) as f:
        d = json.load(f)
    cats = d["categories"]
    return [cats[k][0] for k in sorted(cats.keys(), key=int)]


def detect_columns(path):
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
    kl = sorted([c for c in header if c.startswith("kl_")], key=lambda c: int(c.split("_")[1]))
    score_col = "bearing_score_tested" if "bearing_score_tested" in header else (
        "bearing_score" if "bearing_score" in header else None)
    signed_col = "bearing_score" if "bearing_score" in header else score_col
    extra = [c for c in header
             if c.lower().startswith("neglog10") or "pval" in c.lower()
             or c.startswith("significant_fdr") or c == "score_normalised"]
    padj = next((c for c in header if c.lower() in ("pval_adj_bh", "padj", "qval", "fdr")), None)
    fdr = next((c for c in header if c.startswith("significant_fdr")), None)
    return header, kl, score_col, signed_col, extra, padj, fdr


def _qcat_raw_list(meta):
    """Extract the raw per-state list from a qcat meta field, mirroring
    compare_qcat.parse_qcat_bgz._parse_qcat_meta. Returns a list or None."""
    meta = meta.strip()
    if meta.startswith("{"):
        try:
            return json.loads(meta).get("raw")
        except ValueError:
            return None
    q = meta.find("qcat:")
    if q < 0:
        return None
    r = meta.find(",raw:", q)
    if r < 0:
        return None
    frag = meta[r + 5:]            # everything after ",raw:"
    try:
        return json.loads(frag)    # raw is the last field in the canonical format
    except ValueError:
        pass
    # Fallback: raw list followed by other fields -> bracket-match.
    depth = 0
    for j, ch in enumerate(frag):
        if ch == "[":
            depth += 1
        elif ch == "]":
            depth -= 1
            if depth == 0:
                try:
                    return json.loads(frag[:j + 1])
                except ValueError:
                    return None
    return None


def _raw_sum(raw):
    """Sum a raw list whose entries are scalars or [value, state] pairs."""
    if raw is None:
        return None
    total = 0.0
    for v in raw:
        total += float(v[0]) if isinstance(v, (list, tuple)) else float(v)
    return total


def load_qcat_signal(path):
    """Return dict (chrom, start) -> total raw signal, scanning the whole qcat.
    Reads via pysam (like compare_qcat) with a gzip fallback, and reports
    clearly if no raw: signal field is present."""
    sig = {}
    n_lines = [0]
    n_raw = [0]
    first_meta = [None]

    def handle(parts):
        if len(parts) < 4:
            return
        n_lines[0] += 1
        if first_meta[0] is None:
            first_meta[0] = parts[3][:160]
        s = _raw_sum(_qcat_raw_list(parts[3]))
        if s is None:
            return
        n_raw[0] += 1
        try:
            sig[(parts[0], int(parts[1]))] = s
        except ValueError:
            pass

    used = None
    try:
        import pysam
        tbx = pysam.TabixFile(str(path))
        for chrom in tbx.contigs:
            for rec in tbx.fetch(chrom):
                handle(rec.split("\t"))
        tbx.close()
        used = "pysam"
    except Exception:
        opener = gzip.open if str(path).endswith((".gz", ".bgz")) else open
        with opener(path, "rt") as fh:
            for line in fh:
                if not line or line[0] == "#" or line.startswith("track"):
                    continue
                handle(line.rstrip("\n").split("\t"))
        used = "gzip"

    if n_raw[0] == 0:
        print("  WARNING: no raw: signal field found in {} "
              "({} lines scanned via {}).".format(path, n_lines[0], used))
        if first_meta[0]:
            print("  First meta seen: " + first_meta[0])
        print("  -> The absolute-signal floor needs raw: fields in the qcat. "
              "These files appear not to carry them; see the note from this run.")
    return sig


def parse_sample_sheet_bw(path):
    """Parse a sample sheet into {condition: [ [bigwig paths per track], ... ]},
    one inner list per replicate. Expects tab-separated columns including
    'condition' and 'bw' (comma-separated bigwig paths in track order)."""
    cond_reps = {}
    with open(path, newline="") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for row in rdr:
            cond = (row.get("condition") or "").strip()
            bw = (row.get("bw") or "").strip()
            if not cond or not bw:
                continue
            bws = [p.strip() for p in bw.split(",") if p.strip()]
            cond_reps.setdefault(cond, []).append(bws)
    return cond_reps


def open_condition_handles(bw_lists, base_dir):
    """Open bigwig handles for a condition: list (over replicates) of
    [pyBigWig handles per track]. Paths resolve relative to base_dir."""
    import pyBigWig
    rep_handles = []
    for bws in bw_lists:
        hs = []
        for p in bws:
            pp = p if os.path.isabs(p) else os.path.join(base_dir, p)
            hs.append(pyBigWig.open(pp))
        rep_handles.append(hs)
    return rep_handles


def close_handles(rep_handles):
    for hs in rep_handles:
        for h in hs:
            try:
                h.close()
            except Exception:
                pass


def signal_for_bins(rep_handles, bins, track_aggregate="sum", rep_aggregate="max"):
    """
    Total signal at each given bin, read directly from bigwigs.

    Only the supplied bins are queried (bounded work, independent of genome
    size). Within a replicate the per-track mean coverage is combined by
    track_aggregate ('sum' or 'max'); across replicates by rep_aggregate
    ('max' keeps a bin if any replicate has signal). Returns np.ndarray aligned
    to bins.
    """
    n = len(bins)
    if n == 0:
        return np.zeros(0, dtype=np.float64)
    rep_tot = []
    for hs in rep_handles:
        tot = np.zeros(n, dtype=np.float64)
        for h in hs:
            chd = h.chroms()
            col = np.zeros(n, dtype=np.float64)
            for i, (c, s, e) in enumerate(bins):
                if c not in chd:
                    continue
                try:
                    v = h.stats(c, int(s), int(e), type="mean")
                except RuntimeError:
                    v = None
                val = v[0] if (v and v[0] is not None) else 0.0
                if val > 0:
                    col[i] = val
            tot = np.maximum(tot, col) if track_aggregate == "max" else tot + col
        rep_tot.append(tot)
    stacked = np.vstack(rep_tot)
    if rep_aggregate == "mean":
        return stacked.mean(axis=0)
    if rep_aggregate == "min":
        return stacked.min(axis=0)
    return stacked.max(axis=0)


def load_gtf_genes(path):
    rows = {}
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "gene":
                continue
            try:
                s = int(p[3]) - 1
                e = int(p[4])
            except ValueError:
                continue
            m = re.search(r'gene_name\s+"([^"]*)"', p[8]) or re.search(r'gene_id\s+"([^"]*)"', p[8])
            rows.setdefault(p[0], []).append((s, e, m.group(1) if m else "?"))
    genes = {}
    for ch, lst in rows.items():
        lst.sort()
        genes[ch] = (np.array([r[0] for r in lst]), np.array([r[1] for r in lst]),
                     [r[2] for r in lst])
    return genes


def nearest_gene(genes, chrom, center):
    if chrom not in genes:
        return "", ""
    starts, ends, names = genes[chrom]
    inside = np.where((starts <= center) & (ends >= center))[0]
    if len(inside):
        return names[inside[0]], "0"
    dist = np.where(center > ends, center - ends, np.where(center < starts, starts - center, 0))
    i = int(np.argmin(np.abs(dist)))
    return names[i], str(int(abs(dist[i])))


def top_sites_for_diff(path, names, args, genes, signal_fn):
    header, kl_cols, score_col, signed_col, extra, padj, fdr = detect_columns(path)
    if score_col is None:
        raise SystemExit("{}: no bearing_score[_tested] column".format(path))
    want = list(dict.fromkeys(
        [c for c in ("chrom", "start", "end") if c in header]
        + [score_col, signed_col] + kl_cols + extra))
    df = pd.read_csv(path, sep="\t", usecols=want, low_memory=False).copy()
    if getattr(args, "chroms", None):
        df = df[df["chrom"].isin(set(args.chroms))]
    if not args.all_chroms:
        df = df[df["chrom"].isin(AUTOSOMES)]
    df = df.dropna(subset=[score_col]).reset_index(drop=True)

    # Dominant track per bin (argmax |kl|), used for output and --exclude-track.
    if kl_cols:
        kl_abs = df[kl_cols].abs().to_numpy()
        dom_idx = np.argmax(kl_abs, axis=0) if kl_abs.size == 0 else np.argmax(kl_abs, axis=1)
        df["_dominant"] = [names[i] if i < len(names) else kl_cols[i] for i in dom_idx]
    else:
        df["_dominant"] = ""

    # Candidate selection for signal: only the highest-score bins can win a
    # top-N ranking, so signal is read for those alone. In per-track mode,
    # candidates are picked WITHIN each dominant track, so weaker assays are
    # not starved by the expression-dominated bins.
    combiner = getattr(args, "floor_combiner", "max")
    sig_avail = signal_fn is not None
    per_track = getattr(args, "per_track", False)
    df["_sigA"] = np.nan
    df["_sigB"] = np.nan
    df["_sigcomb"] = np.nan
    if sig_avail and len(df):
        mag_all = df[score_col].abs().to_numpy()
        if per_track:
            dom = df["_dominant"].to_numpy()
            tracks_present = list(pd.unique(df["_dominant"]))
            per_m = max(args.n * 50,
                        int(getattr(args, "signal_candidates", 20000)) // max(1, len(tracks_present)))
            cand_set = set()
            for t in tracks_present:
                gi = np.where(dom == t)[0]
                if len(gi) > per_m:
                    gi = gi[np.argpartition(mag_all[gi], -per_m)[-per_m:]]
                cand_set.update(int(i) for i in gi)
            cand_idx = np.array(sorted(cand_set), dtype=int)
        else:
            K = max(int(getattr(args, "signal_candidates", 20000)), args.n)
            cand_idx = (np.argpartition(mag_all, -K)[-K:] if len(df) > K
                        else np.arange(len(df)))
        cand_bins = [(df.at[i, "chrom"], int(df.at[i, "start"]), int(df.at[i, "end"]))
                     for i in cand_idx]
        print("  resolving signal for {} candidate bins ...".format(len(cand_bins)),
              flush=True)
        sigA_arr, sigB_arr = signal_fn(cand_bins)
        df.loc[cand_idx, "_sigA"] = sigA_arr
        df.loc[cand_idx, "_sigB"] = sigB_arr
        a = df["_sigA"].to_numpy()
        b = df["_sigB"].to_numpy()
        if combiner == "sum":
            df["_sigcomb"] = a + b
        elif combiner == "min":
            df["_sigcomb"] = np.minimum(a, b)
        else:
            df["_sigcomb"] = np.maximum(a, b)

    # Filters.
    n0 = len(df)
    if args.require_fdr_sig and fdr is not None:
        df = df[df[fdr].astype(float) == 1]
    if args.exclude_track:
        df = df[~df["_dominant"].isin(set(args.exclude_track))]
    if args.min_signal > 0:
        if not sig_avail:
            raise SystemExit("--min-signal needs a signal source (--sheet or --qcat-a/--qcat-b)")
        df = df[df["_sigcomb"] >= args.min_signal]   # NaN (non-candidate) -> dropped
    df = df.reset_index(drop=True)
    print("  {} -> {} bins after filters".format(n0, len(df)))

    # Rank key (higher = better).
    mag = df[score_col].abs().to_numpy()
    if args.rank_by == "signal_weighted":
        if not sig_avail:
            raise SystemExit("--rank-by signal_weighted needs a signal source "
                             "(--sheet or --qcat-a/--qcat-b)")
        df["_rankkey"] = mag * np.nan_to_num(df["_sigcomb"].to_numpy(), nan=0.0)
    elif args.rank_by == "pval":
        if padj is None:
            raise SystemExit("--rank-by pval: no adjusted p-value column found")
        df["_rankkey"] = -df[padj].astype(float).to_numpy()
    else:
        df["_rankkey"] = mag
    df = df.dropna(subset=["_rankkey"]).reset_index(drop=True)

    def _merge_top(sub):
        ranked = sub.sort_values("_rankkey", ascending=False)
        chosen = []
        claimed = {}
        md = args.merge_distance
        for idx, row in ranked.iterrows():
            ch = row["chrom"]
            s, e = int(row["start"]), int(row["end"])
            spans = claimed.setdefault(ch, [])
            if any(not (e + md < cs or s - md > ce) for cs, ce in spans):
                continue
            spans.append((s, e))
            chosen.append(idx)
            if len(chosen) >= args.n:
                break
        return chosen

    def _record(rank, idx, track_label):
        row = df.loc[idx]
        rec = {
            "rank": rank, "track": track_label, "chrom": row["chrom"],
            "start": int(row["start"]), "end": int(row["end"]),
            "center": (int(row["start"]) + int(row["end"])) // 2,
            "score": "{:.6g}".format(float(row[score_col])),
            "signed_score": "{:.6g}".format(float(row[signed_col])),
            "dominant_track": row["_dominant"],
        }
        if sig_avail:
            rec["signal_A"] = "" if pd.isna(row["_sigA"]) else "{:.4g}".format(float(row["_sigA"]))
            rec["signal_B"] = "" if pd.isna(row["_sigB"]) else "{:.4g}".format(float(row["_sigB"]))
        for c in extra:
            v = row[c]
            if pd.isna(v):
                rec[c] = ""
            else:
                try:
                    rec[c] = "{:.6g}".format(float(v))
                except (ValueError, TypeError):
                    rec[c] = str(v)
        if genes is not None:
            g, dist = nearest_gene(genes, row["chrom"], rec["center"])
            rec["nearest_gene"] = g
            rec["gene_dist"] = dist
        return rec

    pctile = getattr(args, "min_signal_pctile", None)

    def _apply_pctile(sub, label):
        # Relative floor: keep bins at/above the P-th percentile of combined
        # signal WITHIN this group. Computed per group, so each assay gets a
        # scale-appropriate cut instead of one absolute number across all.
        if pctile is None or not sig_avail:
            return sub
        vals = sub["_sigcomb"].to_numpy()
        finite = vals[np.isfinite(vals)]
        if finite.size == 0:
            return sub
        thr = float(np.percentile(finite, pctile))
        kept = sub[sub["_sigcomb"] >= thr]
        print("  {}: signal pctile {} -> floor {:.4g}, {} -> {} bins".format(
            label, pctile, thr, len(sub), len(kept)))
        return kept

    out = []
    if per_track:
        present = set(df["_dominant"])
        order = [t for t in names if t in present]
        order += [t for t in pd.unique(df["_dominant"]) if t not in order]
        for t in order:
            sub = _apply_pctile(df[df["_dominant"] == t], t)
            chosen = _merge_top(sub)
            for rank, idx in enumerate(chosen, 1):
                out.append(_record(rank, idx, t))
    else:
        chosen = _merge_top(_apply_pctile(df, "all"))
        for rank, idx in enumerate(chosen, 1):
            out.append(_record(rank, idx, ""))
    return out, extra, sig_avail, (genes is not None)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--diff", action="append", required=True, metavar="LABEL=PATH")
    ap.add_argument("--sheet", default=None,
                    help="Sample sheet (tab-separated; columns include 'condition' "
                         "and 'bw' = comma-separated bigwig paths per track). When "
                         "given, the signal floor is read from the source bigwigs: "
                         "all replicates of the reference and of each comparison "
                         "condition. Paths are resolved relative to the sheet.")
    ap.add_argument("--reference", default="DN",
                    help="Reference condition name in the sheet (condition A). Default DN.")
    ap.add_argument("--qcat-a", dest="qcat_a", default=None,
                    help="(Legacy) reference qcat with raw: fields, if present.")
    ap.add_argument("--qcat-b", dest="qcat_b", default=None,
                    help="(Legacy) condition qcat with raw: fields (one --diff at a time).")
    ap.add_argument("--min-signal", type=float, default=0.0,
                    help="Absolute floor: keep bins whose combined signal (see "
                         "--floor-combiner) >= this. 0 = off. One absolute number "
                         "is harsh across assays of different scale; prefer "
                         "--min-signal-pctile with --per-track.")
    ap.add_argument("--min-signal-pctile", dest="min_signal_pctile", type=float,
                    default=None, metavar="P",
                    help="Relative floor: within each group keep bins at/above the "
                         "P-th percentile of combined signal (0-100). With "
                         "--per-track the percentile is computed PER ASSAY, so each "
                         "assay gets a scale-appropriate cut instead of one global "
                         "absolute number. e.g. 50 keeps the upper-signal half of "
                         "each assay's high-score bins.")
    ap.add_argument("--floor-combiner", dest="floor_combiner",
                    choices=["max", "sum", "min"], default="max",
                    help="Combine the two conditions' signal: 'max' keeps on/off "
                         "events (signal in only one condition); 'sum' adds; 'min' "
                         "requires signal in both. Default max.")
    ap.add_argument("--replicate-aggregate", dest="replicate_aggregate",
                    choices=["max", "mean", "min"], default="max",
                    help="Combine replicate signal within a condition (bigwig mode): "
                         "'max' keeps a bin if ANY replicate has signal (best chance "
                         "of keeping real loci); 'mean'/'min' are more conservative. "
                         "Default max.")
    ap.add_argument("--track-aggregate", dest="track_aggregate",
                    choices=["sum", "max"], default="sum",
                    help="Combine per-assay signal within a replicate: 'sum' = total "
                         "activity across assays; 'max' = strongest single assay. "
                         "Default sum.")
    ap.add_argument("--bin-size", dest="bin_size", type=int, default=200,
                    help="Bin size for bigwig signal; must match the qcat grid. Default 200.")
    ap.add_argument("--signal-candidates", dest="signal_candidates", type=int, default=20000,
                    help="Only read bigwig signal for the top-K highest-score bins "
                         "(per comparison). The top-N can only come from high-score "
                         "bins, so this bounds the work and avoids scanning every bin. "
                         "Raise if a strong floor leaves too few. Default 20000.")
    ap.add_argument("--signal-cache", dest="signal_cache", default=None, metavar="DIR",
                    help="Cache per-comparison candidate signal here. The candidate "
                         "set does not depend on --min-signal, so threshold sweeps "
                         "reuse the cache and are instant after the first run.")
    ap.add_argument("--chroms", nargs="+", default=None,
                    help="Restrict to these chromosomes (e.g. --chroms chr6 for a quick test).")
    ap.add_argument("--rank-by", choices=["score", "signal_weighted", "pval"],
                    default="score")
    ap.add_argument("--require-fdr-sig", action="store_true",
                    help="Keep only significant_fdr0.05 == 1 bins.")
    ap.add_argument("--exclude-track", action="append", default=[], metavar="NAME",
                    help="Drop sites whose dominant track is NAME (repeatable).")
    ap.add_argument("--categories", default=None)
    ap.add_argument("--gtf", default=None)
    ap.add_argument("--n", type=int, default=10,
                    help="Top sites to report (per dominant track when --per-track).")
    ap.add_argument("--per-track", dest="per_track", action="store_true",
                    help="Report a top-N for EACH dominant assay track separately, "
                         "instead of one global top-N. Gives the chromatin assays "
                         "(CTCF, Cohesin, ATAC, H3K27ac) their own ranked lists "
                         "instead of being crowded out by expression (RNAseq).")
    ap.add_argument("--merge-distance", type=int, default=1000)
    ap.add_argument("--all-chroms", action="store_true")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    names = load_categories_json(args.categories) if args.categories else []
    genes = load_gtf_genes(args.gtf) if args.gtf else None
    if genes is not None:
        print("Loaded genes for {} chromosomes".format(len(genes)))

    sheet_mode = bool(args.sheet)
    cond_reps = None
    base_dir = None
    ref_handles = None
    sigA_qcat = sigB_qcat = None

    if sheet_mode:
        cond_reps = parse_sample_sheet_bw(args.sheet)
        base_dir = os.path.dirname(os.path.abspath(args.sheet))
        if args.reference not in cond_reps:
            raise SystemExit("Reference condition '{}' not in sheet. Available: {}"
                             .format(args.reference, ", ".join(sorted(cond_reps))))
        ref_handles = open_condition_handles(cond_reps[args.reference], base_dir)
    elif args.qcat_a and args.qcat_b:
        print("Loading qcat signal A ...", flush=True)
        sigA_qcat = load_qcat_signal(args.qcat_a)
        print("  {} bins".format(len(sigA_qcat)))
        print("Loading qcat signal B ...", flush=True)
        sigB_qcat = load_qcat_signal(args.qcat_b)
        print("  {} bins".format(len(sigB_qcat)))
        if len(args.diff) > 1:
            raise SystemExit(
                "--qcat-b is condition-specific; run one --diff at a time when "
                "using qcat signal. (Use --sheet for multiple comparisons.)")
    elif args.qcat_a or args.qcat_b:
        raise SystemExit("Provide BOTH --qcat-a and --qcat-b, or use --sheet.")

    def _cache_lookup(tag, cand_bins, paths_a, paths_b, compute):
        if not args.signal_cache:
            return compute()
        params = (args.bin_size, args.track_aggregate, args.replicate_aggregate,
                  args.floor_combiner)
        h = hashlib.md5()
        h.update(repr([tag, sorted(paths_a), sorted(paths_b), params]).encode())
        h.update(np.asarray([s for _, s, _ in cand_bins], dtype=np.int64).tobytes())
        h.update("|".join(c for c, _, _ in cand_bins).encode())
        fpath = os.path.join(args.signal_cache, "sig_{}_{}.npz".format(tag, h.hexdigest()[:16]))
        if os.path.exists(fpath):
            d = np.load(fpath)
            print("  (cached signal: {})".format(fpath))
            return d["a"], d["b"]
        a, b = compute()
        try:
            os.makedirs(args.signal_cache, exist_ok=True)
            np.savez(fpath, a=a, b=b)
        except Exception as exc:
            print("  WARNING: could not write signal cache: {}".format(exc))
        return a, b

    def make_signal_fn(label):
        if sheet_mode:
            cond_handles = open_condition_handles(cond_reps[label], base_dir)
            ref_paths = [p for rep in cond_reps[args.reference] for p in rep]
            cond_paths = [p for rep in cond_reps[label] for p in rep]

            def fn(cand_bins):
                def compute():
                    a = signal_for_bins(ref_handles, cand_bins,
                                        args.track_aggregate, args.replicate_aggregate)
                    b = signal_for_bins(cond_handles, cand_bins,
                                        args.track_aggregate, args.replicate_aggregate)
                    return a, b
                a, b = _cache_lookup(label, cand_bins, ref_paths, cond_paths, compute)
                return a, b
            return fn, cond_handles
        if sigA_qcat is not None:
            def fn(cand_bins):
                a = np.array([sigA_qcat.get((c, s), 0.0) for c, s, _ in cand_bins])
                b = np.array([sigB_qcat.get((c, s), 0.0) for c, s, _ in cand_bins])
                return a, b
            return fn, None
        return None, None

    all_rows, extra_cols, had_sig, had_genes = [], [], False, False
    for item in args.diff:
        if "=" not in item:
            raise SystemExit("--diff expects LABEL=PATH: " + item)
        label, path = item.split("=", 1)
        if sheet_mode and label not in cond_reps:
            raise SystemExit("Comparison condition '{}' not in sheet. Available: {}"
                             .format(label, ", ".join(sorted(cond_reps))))
        signal_fn, cond_handles = make_signal_fn(label)
        print("\n# {}  (top {} by {}{})".format(
            label, args.n, args.rank_by,
            "; min-signal {} [{}]".format(args.min_signal, args.floor_combiner)
            if args.min_signal else ""))
        try:
            sites, extra, had_sig, had_genes = top_sites_for_diff(
                path, names, args, genes, signal_fn)
        finally:
            if cond_handles is not None:
                close_handles(cond_handles)
        extra_cols = extra
        cols = ["rank"]
        if args.per_track:
            cols += ["track"]
        cols += ["chrom", "start", "end", "score", "signed_score", "dominant_track"]
        if had_sig:
            cols += ["signal_A", "signal_B"]
        if had_genes:
            cols += ["nearest_gene", "gene_dist"]
        cols += extra_cols
        print("  " + "  ".join("{:>10}".format(c) for c in cols))
        for r in sites:
            print("  " + "  ".join("{:>10}".format(str(r.get(c, ""))) for c in cols))
            r["comparison"] = label
            all_rows.append(r)

        # Per-assay (or overall) significance summary.
        fdr_key = next((c for c in extra_cols if c.startswith("significant_fdr")), None)
        if fdr_key:
            def _is_sig(rec):
                v = rec.get(fdr_key, "")
                try:
                    return float(v) == 1.0
                except (ValueError, TypeError):
                    return str(v).lower() in ("true", "1")
            groups = {}
            for r in sites:
                groups.setdefault(r.get("track", "all"), []).append(r)
            order = [t for t in names if t in groups] + \
                    [t for t in groups if t not in names]
            parts = []
            for t in (order if args.per_track else ["all"]):
                grp = groups.get(t, []) if args.per_track else sites
                if not grp:
                    continue
                ns = sum(1 for r in grp if _is_sig(r))
                parts.append("{} {}/{} sig".format(t if args.per_track else "FDR0.05",
                                                    ns, len(grp)))
            if parts:
                print("  [significant at FDR 0.05]  " + " | ".join(parts))

    if ref_handles is not None:
        close_handles(ref_handles)

    base = ["comparison", "rank"]
    if args.per_track:
        base += ["track"]
    base += ["chrom", "start", "end", "center", "score", "signed_score", "dominant_track"]
    if had_sig:
        base += ["signal_A", "signal_B"]
    if had_genes:
        base += ["nearest_gene", "gene_dist"]
    header = base + extra_cols
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in all_rows:
            w.writerow(r)
    print("\nWrote", args.out, "({} rows)".format(len(all_rows)))


if __name__ == "__main__":
    main()
