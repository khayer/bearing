#!/usr/bin/env python3
"""
cbe_across_runs.py -- consolidate per-CBE significance across run regimes.

Answers one question: does ANY known CTCF binding element (CBE) reach FDR
significance under ANY scoring regime (default / qnorm / seed / jsd / adaptive)
and ANY comparison? Reads each run's cbe_point_query.tsv (produced by the
cbe_point_query rule) and reports a combined matrix + verdict.

Each CBE is a ~18 bp point feature mapped to the bin that contains it (the
regime's bin grid -- 200 bp fixed, or variable-width adaptive), so this is a
point query, not a regional-enrichment test (an 18 bp region contains no whole
200 bp bin, so regional enrichment over CBEs is structurally empty).

  python cbe_across_runs.py \
    --query default=workflow/results \
    --query jsd=workflow/results_jsd \
    --query adaptive=workflow/results_adaptive \
    --out cbe_significance_across_runs.tsv

PATH may be a run directory (then <dir>/regional/cbe_point_query.tsv is read) or
a cbe_point_query.tsv file directly. ASCII only.
"""
import argparse
import csv
import os
import sys


def resolve(path):
    if os.path.isdir(path):
        return os.path.join(path, "regional", "cbe_point_query.tsv")
    return path


def load_query(path):
    """Return {(cbe_name, comparison): (is_sig, bin_status)} from a query TSV."""
    out = {}
    with open(path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            name = r.get("cbe_name", "")
            comp = r.get("comparison", "")
            sig = str(r.get("fdr_significant", "")).strip() in ("1", "1.0")
            out[(name, comp)] = (sig, r.get("bin_status", ""))
    return out


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--query", action="append", required=True, metavar="LABEL=PATH",
                    help="run label and its run dir or cbe_point_query.tsv "
                         "(repeatable)")
    ap.add_argument("--out", default=None, help="combined matrix TSV")
    args = ap.parse_args()

    runs = []
    for spec in args.query:
        if "=" not in spec:
            sys.exit("[ERROR] --query expects LABEL=PATH, got: %s" % spec)
        label, path = spec.split("=", 1)
        qpath = resolve(path)
        if not os.path.exists(qpath):
            sys.exit("[ERROR] no cbe_point_query.tsv for run '%s' at %s\n"
                     "        build it: snakemake .../%s/regional/cbe_point_query.tsv"
                     % (label, qpath, os.path.basename(path.rstrip("/"))))
        runs.append((label, load_query(qpath)))

    # universe of CBEs and comparisons
    cbes, comps = [], []
    for _lab, q in runs:
        for (name, comp) in q:
            if name not in cbes:
                cbes.append(name)
            if comp not in comps:
                comps.append(comp)
    cbes.sort()
    comps.sort()

    # per-run, per-comparison significant counts
    print("=" * 76)
    print("CBE significance across run regimes  (FDR-significant CBE x comparison)")
    print("-" * 76)
    print("  %d CBEs x %d comparisons x %d runs" % (len(cbes), len(comps), len(runs)))
    print("-" * 76)
    grand = 0
    per_run_total = {}
    for label, q in runs:
        n_sig = sum(1 for v in q.values() if v[0])
        per_run_total[label] = n_sig
        grand += n_sig
        by_comp = {}
        for (name, comp), (sig, _st) in q.items():
            if sig:
                by_comp[comp] = by_comp.get(comp, 0) + 1
        detail = ", ".join("%s:%d" % (c, by_comp[c]) for c in sorted(by_comp)) or "none"
        print("  %-12s  %3d significant   (%s)" % (label, n_sig, detail))
    print("-" * 76)

    # CBEs significant in >= 1 (run, comparison)
    ever_sig = sorted({name for _lab, q in runs
                       for (name, _c), (sig, _st) in q.items() if sig})
    print("  CBEs significant in >= 1 regime x comparison: %d" % len(ever_sig))
    if ever_sig:
        for nm in ever_sig:
            where = []
            for label, q in runs:
                hits = [c for (n, c), (sig, _st) in q.items() if n == nm and sig]
                if hits:
                    where.append("%s[%s]" % (label, ",".join(sorted(hits))))
            print("    %-24s %s" % (nm, "; ".join(where)))
        print("\n  VERDICT: at least one known CTCF site reaches significance.")
    else:
        print("\n  VERDICT: NO known CTCF site reaches FDR significance in any "
              "regime x comparison.\n  Consistent with BEARING's CTCF differential "
              "being compositional / low-coverage\n  (CTCF Q-drift ~0.002 across "
              "sample pairs); coverage-equalized adaptive\n  binning does not change "
              "the conclusion." if grand == 0 else "")
    print("=" * 76)

    if args.out:
        with open(args.out, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            cols = ["cbe_name"]
            for label, _q in runs:
                for comp in comps:
                    cols.append("%s|%s" % (label, comp))
            cols.append("n_sig_total")
            w.writerow(cols)
            for name in cbes:
                row = [name]
                tot = 0
                for label, q in runs:
                    for comp in comps:
                        v = q.get((name, comp))
                        if v is None:
                            row.append("")
                        elif v[1] == "no_bin":
                            row.append("nb")
                        else:
                            row.append("1" if v[0] else "0")
                            tot += int(v[0])
                row.append(str(tot))
                w.writerow(row)
        print("wrote %s" % args.out)


if __name__ == "__main__":
    main()
