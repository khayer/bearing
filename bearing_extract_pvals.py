#!/usr/bin/env python3
"""
bearing_extract_pvals.py

Return the V1P/4C FDR numbers for the manuscript Results paragraph (R3) WITHOUT
uploading the big file. Runs locally against the bearing_pvalue.py --diff output
for the DN-vs-dV1P comparison (the file with per-bin p-values and BH-adjusted
q-values).

It reports, per pre-specified Tcrb region and genome-wide:
  - number of scored/tested bins in the region
  - number significant at FDR < 0.05
  - minimum and median adjusted p-value among significant bins
  - the list of significant bins (chrom, start, end, padj)
and genome-wide:
  - total significant bins, and how many fall OUTSIDE the Tcrb locus
    (to support the claim that the V1P signal is focal, not scattered).

ASCII only. Standard library only (gzip handles .bgz because bgzip output is
gzip-readable). No repo import required.

USAGE
  python3 bearing_extract_pvals.py PVALFILE
  python3 bearing_extract_pvals.py PVALFILE --fdr 0.05
  # if column auto-detection fails, first inspect the header:
  python3 bearing_extract_pvals.py PVALFILE --describe

Default regions are the manuscript v5 Tcrb BED (Trbv1, Vbeta cluster, NEGCTL,
DJ_RC, Ebeta/Trbv31/CBE3). Override with --regions BEDFILE (chrom start end name).
"""
import argparse, gzip, io, sys, re, statistics

# ---- manuscript v5 Tcrb regions (mm10), chrom start end name -----------------
DEFAULT_REGIONS = [
    ("chr6", 40880000, 40905000, "Trbv1"),
    ("chr6", 41040000, 41290000, "Vbeta_cluster"),
    ("chr6", 41313087, 41486888, "NEGCTL"),
    ("chr6", 41500000, 41551500, "DJ_RC"),
    ("chr6", 41552000, 41562000, "Ebeta_Trbv31_CBE3"),
]
# whole Tcrb locus (focal zoom) used for the "focal, not scattered" statement
TCRB_LOCUS = ("chr6", 40790000, 41690000)


def opener(path):
    if path.endswith(".gz") or path.endswith(".bgz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def sniff_header(path):
    """Return (delimiter, header_list, first_data_line, header_present)."""
    with opener(path) as fh:
        first = None
        for line in fh:
            if line.strip() == "" or line.startswith("#"):
                continue
            first = line.rstrip("\n")
            break
    delim = "\t" if "\t" in first else ("," if "," in first else None)
    cols = first.split(delim) if delim else first.split()
    # header present if it has non-numeric tokens in the expected name slots
    looks_header = any(re.search(r"[A-Za-z_]", c) for c in cols[:4]) and not _is_float(cols[1] if len(cols) > 1 else "")
    return delim, cols, looks_header


def _is_float(s):
    try:
        float(s); return True
    except Exception:
        return False


def find_cols(header):
    """Map logical fields to column indices from a header row (case-insensitive)."""
    h = [c.strip().lower() for c in header]
    def pick(cands, required=True, contains=False):
        for i, name in enumerate(h):
            for c in cands:
                if (name == c) or (contains and c in name):
                    return i
        if required:
            return None
        return None
    idx = {}
    idx["chrom"] = pick(["chrom", "chr", "seqnames", "#chrom"], contains=True)
    idx["start"] = pick(["start", "chromstart", "pos", "bin_start"], contains=True)
    idx["end"]   = pick(["end", "chromend", "stop", "bin_end"], contains=True)
    # adjusted p / q value
    idx["padj"]  = pick(["pval_adj_bh", "padj", "p_adj", "qvalue", "q_value", "q", "fdr", "adj_pval", "pval_adj"], contains=True)
    # raw p (fallback / reporting)
    idx["pval"]  = pick(["pval", "p_value", "pvalue", "p"], required=False, contains=True)
    # explicit significance flag if present
    idx["sig"]   = pick(["significant_fdr0.05", "significant", "sig"], required=False, contains=True)
    return idx


def main():
    ap = argparse.ArgumentParser(description="Extract V1P/4C FDR numbers from a bearing_pvalue --diff table.")
    ap.add_argument("pvalfile")
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--regions", help="BED file: chrom start end name (overrides defaults)")
    ap.add_argument("--describe", action="store_true", help="print detected delimiter/header/columns and exit")
    ap.add_argument("--col-chrom", type=int, help="0-based column index override")
    ap.add_argument("--col-start", type=int)
    ap.add_argument("--col-end", type=int)
    ap.add_argument("--col-padj", type=int)
    args = ap.parse_args()

    delim, header, looks_header = sniff_header(args.pvalfile)
    if args.describe:
        print("delimiter:", repr(delim))
        print("header present:", looks_header)
        for i, c in enumerate(header):
            print(f"  col[{i}] = {c!r}")
        print("\nIf the padj/chrom/start/end columns are not auto-found, re-run with")
        print("  --col-chrom N --col-start N --col-end N --col-padj N")
        return

    idx = find_cols(header) if looks_header else {"chrom":0,"start":1,"end":2,"padj":None,"pval":None,"sig":None}
    for k in ("chrom","start","end","padj"):
        ov = getattr(args, "col_"+k, None)
        if ov is not None:
            idx[k] = ov
    if idx["chrom"] is None or idx["start"] is None or idx["padj"] is None:
        sys.exit("ERROR: could not locate chrom/start/padj columns. Run with --describe, then pass --col-* overrides.")

    regions = DEFAULT_REGIONS
    if args.regions:
        regions = []
        with open(args.regions) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                p = line.split()
                regions.append((p[0], int(p[1]), int(p[2]), p[3] if len(p) > 3 else f"{p[0]}:{p[1]}-{p[2]}"))

    # accumulators
    region_hits = {r[3]: [] for r in regions}
    region_n    = {r[3]: 0 for r in regions}
    gw_sig = 0
    tcrb_sig = 0
    total_rows = 0

    with opener(args.pvalfile) as fh:
        started = False
        for line in fh:
            if line.strip() == "" or line.startswith("#"):
                continue
            if looks_header and not started:
                started = True
                continue  # skip header row
            parts = line.rstrip("\n").split(delim) if delim else line.split()
            try:
                chrom = parts[idx["chrom"]]
                start = int(float(parts[idx["start"]]))
                end   = int(float(parts[idx["end"]])) if idx["end"] is not None else start + 200
                padj  = parts[idx["padj"]]
                if padj == "" or padj.lower() in ("na", "nan", "none"):
                    continue
                padj = float(padj)
            except (IndexError, ValueError):
                continue
            total_rows += 1
            is_sig = padj < args.fdr
            if is_sig:
                gw_sig += 1
                if chrom == TCRB_LOCUS[0] and start >= TCRB_LOCUS[1] and start < TCRB_LOCUS[2]:
                    tcrb_sig += 1
            for (rc, rs, re_, name) in regions:
                if chrom == rc and start >= rs and start < re_:
                    region_n[name] += 1
                    if is_sig:
                        region_hits[name].append((chrom, start, end, padj))

    # ---- report ----
    print(f"# bearing_extract_pvals.py  file={args.pvalfile}  FDR<{args.fdr}")
    print(f"# rows with a usable adjusted p-value: {total_rows}")
    print(f"# genome-wide significant bins: {gw_sig}")
    print(f"# significant bins inside the Tcrb focal locus "
          f"({TCRB_LOCUS[0]}:{TCRB_LOCUS[1]:,}-{TCRB_LOCUS[2]:,}): {tcrb_sig}")
    print(f"# significant bins OUTSIDE the Tcrb locus: {gw_sig - tcrb_sig}  "
          f"(focality: {'GOOD - concentrated at Tcrb' if tcrb_sig and (gw_sig - tcrb_sig) < tcrb_sig else 'inspect'})")
    print()
    print(f"{'region':22s} {'n_bins':>7s} {'n_sig':>6s} {'min_padj':>11s} {'median_padj_sig':>16s}")
    for (rc, rs, re_, name) in regions:
        hits = region_hits[name]
        if hits:
            padjs = [h[3] for h in hits]
            print(f"{name:22s} {region_n[name]:7d} {len(hits):6d} {min(padjs):11.3e} {statistics.median(padjs):16.3e}")
        else:
            print(f"{name:22s} {region_n[name]:7d} {0:6d} {'-':>11s} {'-':>16s}")
    print()
    # detailed Trbv1 list (the R3 headline)
    for (rc, rs, re_, name) in regions:
        if name == "Trbv1" and region_hits[name]:
            print("# Trbv1 significant bins (the R3 loop-loss recovery):")
            for (c, s, e, q) in sorted(region_hits[name], key=lambda x: x[3]):
                print(f"   {c}:{s:,}-{e:,}   padj={q:.3e}")
    print()
    print("# READY-TO-PASTE for R3 (fill the brackets from the numbers above):")
    tr = region_hits.get("Trbv1", [])
    rc = region_hits.get("DJ_RC", [])
    if tr:
        print(f'   "...recovered as {len(tr)} FDR-significant bins at Trbv1 '
              f'(minimum adjusted p = {min(h[3] for h in tr):.1e}), with the recombination center '
              f'itself unchanged ({len(rc)} significant bins), under the same circular-shift '
              f'permutation null and Benjamini-Hochberg FDR."')


if __name__ == "__main__":
    main()
