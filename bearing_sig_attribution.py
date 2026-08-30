#!/usr/bin/env python3
"""
bearing_sig_attribution.py

The R3-defensible number: per-track attribution restricted to the FDR-SIGNIFICANT
bins in each region (not all bins), by joining the differential qcat with the
p-value table. Answers: "of the bins that actually reach significance at Trbv1,
what fraction of the differential is the contact channel?" - and, as a control,
shows the same for the NEGCTL region so the platform floor is explicit.

Why: over ALL bins the O/E contact channels dominate |diff| everywhere in the
locus (the DN-vs-V1P maps differ by platform), so the all-bins attribution is not
specific. Restricting to FDR-significant bins is the honest readout.

Inputs:
  QCAT  = diff_DN_vs_V1P.qcat.bgz  (bgzip+tabix; contact channels in the panel)
  PVAL  = diff_DN_vs_V1P.stats.tsv (per-bin pval + BH-adjusted q)

ASCII only. Uses tabix (or pysam) for the qcat; streams the pval tsv once.

USAGE
  python3 bearing_sig_attribution.py diff_DN_vs_V1P.qcat.bgz diff_DN_vs_V1P.stats.tsv \
      --track-names "RNAseq +,RNAseq -,CTCF,Cohesin,NIPBL,H3K27ac,4C RC,4C Trbv1" \
      --contact-name "4C RC" --fdr 0.05
"""
import argparse, gzip, io, json, os, shutil, subprocess, sys

DEFAULT_REGIONS = [
    ("chr6", 40880000, 40905000, "Trbv1"),
    ("chr6", 41040000, 41290000, "Vbeta_cluster"),
    ("chr6", 41313087, 41486888, "NEGCTL"),
    ("chr6", 41500000, 41551500, "DJ_RC"),
    ("chr6", 41552000, 41562000, "Ebeta_Trbv31_CBE3"),
]


def opener(path):
    if path.endswith((".gz", ".bgz")):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def is_indexed(path):
    return path.endswith((".bgz", ".gz")) and (os.path.exists(path + ".tbi") or os.path.exists(path + ".csi"))


def region_lines(path, chrom, start, end):
    if is_indexed(path) and shutil.which("tabix"):
        p = subprocess.Popen(["tabix", path, f"{chrom}:{start+1}-{end}"], stdout=subprocess.PIPE, text=True)
        for line in p.stdout:
            yield line.rstrip("\n")
        p.stdout.close(); p.wait()
        return
    try:
        import pysam
        if is_indexed(path):
            for row in pysam.TabixFile(path).fetch(chrom, start, end):
                yield row
            return
    except Exception:
        pass
    with opener(path) as fh:  # slow fallback
        for line in fh:
            if not line or line[0] == "#":
                continue
            p = line.split("\t", 3)
            if len(p) >= 2 and p[0] == chrom:
                try:
                    s = int(p[1])
                except ValueError:
                    continue
                if start <= s < end:
                    yield line.rstrip("\n")


def parse_qcat_payload(payload):
    m = payload.find("qcat:")
    if m < 0:
        return None
    rest = payload[m + 5:]
    cut = rest.find(",raw:")
    arr = (rest[:cut] if cut >= 0 else rest).strip().rstrip(",")
    try:
        pairs = json.loads(arr)
    except Exception:
        return None
    out = []
    for p in pairs:
        try:
            out.append((int(p[1]), float(p[0])))
        except Exception:
            pass
    return out


def sniff_pval_cols(path):
    with opener(path) as fh:
        for line in fh:
            if line.strip() and not line.startswith("#"):
                cols = line.rstrip("\n").split("\t")
                low = [c.strip().lower() for c in cols]
                def find(cs):
                    for i, n in enumerate(low):
                        for c in cs:
                            if c in n:
                                return i
                    return None
                return {
                    "chrom": find(["chrom", "chr", "seqname"]),
                    "start": find(["start", "pos"]),
                    "padj": find(["pval_adj_bh", "padj", "qvalue", "q_value", "fdr", "adj"]),
                }, cols
    return None, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("qcat")
    ap.add_argument("pval")
    ap.add_argument("--track-names")
    ap.add_argument("--contact-name")
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--regions", help="BED override: chrom start end name (default = 5 Tcrb regions)")
    args = ap.parse_args()

    names = args.track_names.split(",") if args.track_names else None
    cols, header = sniff_pval_cols(args.pval)
    if cols is None or cols["chrom"] is None or cols["start"] is None or cols["padj"] is None:
        sys.exit(f"ERROR: could not find chrom/start/padj in pval header: {header}")

    regions = DEFAULT_REGIONS
    if args.regions:
        regions = []
        with open(args.regions) as fh:
            for line in fh:
                if line.strip() and not line.startswith("#"):
                    p = line.split()
                    regions.append((p[0], int(p[1]), int(p[2]), p[3] if len(p) > 3 else f"{p[0]}:{p[1]}"))

    # 1) collect ALL FDR-significant (chrom,start) genome-wide from the pval tsv
    #    (the significant set is small - a few thousand bins - so region-agnostic is fine)
    sig = set()
    with opener(args.pval) as fh:
        first = True
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            if first:
                first = False
                continue  # header
            p = line.rstrip("\n").split("\t")
            try:
                q = float(p[cols["padj"]])
            except (IndexError, ValueError):
                continue
            if q < args.fdr:
                try:
                    sig.add((p[cols["chrom"]], int(float(p[cols["start"]]))))
                except (IndexError, ValueError):
                    continue

    # base detection (0- vs 1-based track ids) from Tcrb qcat rows (file-global property)
    prescan = set()
    for line in region_lines(args.qcat, "chr6", 40790000, 41690000):
        pr = parse_qcat_payload(line.split("\t", 3)[3]) if line.count("\t") >= 3 else None
        if pr:
            for tid, _ in pr:
                prescan.add(tid)
        if len(prescan) >= 8:
            break
    base = 0 if (0 in prescan or not prescan) else 1

    def tname(tid):
        j = tid - base
        return names[j].strip() if (names and 0 <= j < len(names)) else f"track{tid}"

    key = (args.contact_name or "").lower()
    def is_contact(tid):
        nm = tname(tid).lower()
        return (key in nm) if key else any(k in nm for k in ("4c", "contact", "oe", "hic"))

    print(f"# sig-restricted attribution  FDR<{args.fdr}  (Tcrb-span significant bins: {len(sig)})")
    print(f"{'region':20s} {'n_sig':>6s} {'contact%':>9s}  dominant per-track over SIGNIFICANT bins")
    for (rc, rs, re_, name) in regions:
        acc = {}
        n = 0
        for line in region_lines(args.qcat, rc, rs, re_):
            parts = line.split("\t")
            try:
                s = int(parts[1])
            except (IndexError, ValueError):
                continue
            if (rc, s) not in sig:
                continue
            pairs = parse_qcat_payload(parts[3]) if len(parts) > 3 else None
            if not pairs:
                continue
            n += 1
            for tid, sc in pairs:
                acc[tid] = acc.get(tid, 0.0) + abs(sc)
        tot = sum(acc.values())
        if n == 0 or tot == 0:
            print(f"{name:20s} {n:6d} {'-':>9s}  (no significant bins)")
            continue
        cfrac = 100.0 * sum(v for t, v in acc.items() if is_contact(t)) / tot
        top = sorted(acc.items(), key=lambda kv: -kv[1])[:3]
        topstr = ", ".join(f"{tname(t)} {100.0*v/tot:.0f}%" for t, v in top)
        print(f"{name:20s} {n:6d} {cfrac:8.1f}%  {topstr}")
    print()
    print("# Read: compare Trbv1 (real loop anchor) vs NEGCTL (negative control). If the contact% and")
    print("# the per-track profile look the SAME in NEGCTL as in Trbv1, the signal is platform floor,")
    print("# not loop-specific -> lean on cross-condition insulation specificity + the ARIMA-matched run.")


if __name__ == "__main__":
    main()
