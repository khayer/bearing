#!/usr/bin/env python3
"""
dchic_region_result.py -- pull a named region's result out of dcHiC output and
state it in a form that can be cited.

dcHiC writes, under
    DifferentialResult/<diffdir>/fdr_result/
        differential.intra_sample_combined.bedGraph          (all bins)
        differential.intra_sample_combined.Filtered.bedGraph (significant only)
with per-experiment mean PC values and an adjusted p-value per bin.

This script reports, for a region of interest:
  - which bins overlap it
  - the mean PC per experiment in those bins (sign = compartment call)
  - whether the compartment SIGN differs between two named experiments
  - the dcHiC adjusted p-value
  - how many bins the region actually spans, which is usually the crux

USAGE
-----
    python dchic_region_result.py \
        --bedgraph DifferentialResult/WT_vs_V1P/fdr_result/differential.intra_sample_combined.bedGraph \
        --region chr6:40872100-40908238 --label V1 \
        --exp-a WT --exp-b V1P \
        --out v1_compartment_result.tsv

ASCII only.
"""

import argparse
import sys


def parse_region(s):
    chrom, rest = s.split(":")
    lo, hi = rest.replace(",", "").split("-")
    return chrom, int(lo), int(hi)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bedgraph", required=True,
                    help="differential.intra_sample_combined.bedGraph from dcHiC")
    ap.add_argument("--region", required=True, help="chr6:40872100-40908238")
    ap.add_argument("--label", default="region")
    ap.add_argument("--exp-a", required=True, help="experiment prefix, e.g. WT")
    ap.add_argument("--exp-b", required=True, help="experiment prefix, e.g. V1P")
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    chrom, lo, hi = parse_region(a.region)

    with open(a.bedgraph) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [l.rstrip("\n").split("\t") for l in fh]
    print("columns: %s" % header)

    def col(name):
        for i, h in enumerate(header):
            if h.strip().lower() == name.strip().lower():
                return i
        return None

    ci, si, ei = col("chr") or 0, col("start") or 1, col("end") or 2
    ai, bi = col(a.exp_a), col(a.exp_b)
    pi = col("padj") or col("adj.pvalue") or col("qvalue") or col("padjust")
    if ai is None or bi is None:
        sys.exit("ERROR: could not find columns %r and %r in the bedGraph header.\n"
                 "  Available: %s" % (a.exp_a, a.exp_b, header))

    # Detect a multi-condition MFA run. dcHiC computes ONE padj per bin that tests
    # whether that bin varies across ALL experiments in the run, not the pairwise
    # exp_a-vs-exp_b difference. With 3+ experiments the padj is therefore a
    # GLOBAL (across-panel) significance, and a bin can be "significant" here while
    # exp_a and exp_b are indistinguishable -- the significance may be driven by a
    # third condition (e.g. a fibroblast). The SIGN CHANGE between exp_a and exp_b
    # is still read pairwise from the two PC columns and is reliable; the padj is
    # not pairwise. Label it so downstream readers of the TSV are not misled.
    _meta = {"chr", "start", "end", "replicate_wt", "sample_maha",
             "pval", "padj", "dist_clust", "padjust", "adj.pvalue", "qvalue"}
    exp_cols = [h for h in header
                if h.strip().lower() not in _meta
                and not re.search(r"_rep\d", h)]
    multi_condition = len(exp_cols) > 2
    padj_label = "padj_global" if multi_condition else "padj"
    if multi_condition:
        print("\n  NOTE: %d experiments in this run (%s). dcHiC's padj is a GLOBAL"
              % (len(exp_cols), ", ".join(exp_cols)))
        print("  test of variation across all of them, NOT the pairwise %s-vs-%s"
              % (a.exp_a, a.exp_b))
        print("  difference. Trust the SIGN CHANGE (read from the two PC columns)")
        print("  as the pairwise quantity; treat padj as panel-level significance.")
        print("  For a clean pairwise test, rerun dcHiC on just these two conditions.")
    if pi is None:
        print("WARNING: no adjusted-p column found; reporting PC values only",
              file=sys.stderr)

    hits = []
    for r in rows:
        try:
            if r[ci] != chrom:
                continue
            s, e = int(r[si]), int(r[ei])
        except (ValueError, IndexError):
            continue
        if s < hi and e > lo:
            hits.append(r)

    if not hits:
        sys.exit("No dcHiC bins overlap %s. Check the chromosome naming and the "
                 "resolution." % a.region)

    binsize = int(hits[0][ei]) - int(hits[0][si])
    print("\n=== %s  %s ===" % (a.label, a.region))
    print("region width : %s bp" % format(hi - lo, ","))
    print("bin size     : %s bp" % format(binsize, ","))
    print("bins covered : %d" % len(hits))
    if len(hits) <= 2:
        print("\n  NOTE: this region spans %d compartment bin(s). A compartment call")
        print("  here is a statement about %d bin(s), not a regional trend. dcHiC's own"
              % len(hits))
        print("  documentation warns that isolated significant bins can be misleading,")
        print("  and its --distclust / --numberclust filters exist for exactly this.")

    out_rows = []
    _sigtag = "(sig, global)" if multi_condition else "(significant)"
    print("\n  %-22s %10s %10s %12s %s" % ("bin", a.exp_a, a.exp_b, padj_label, "sign change?"))
    for r in hits:
        va, vb = float(r[ai]), float(r[bi])
        p = float(r[pi]) if pi is not None and r[pi] not in ("NA", "") else float("nan")
        flip = (va > 0) != (vb > 0)
        sig = (p <= a.alpha) if p == p else False
        print("  %s:%-12s %10.4f %10.4f %12.3g %s%s"
              % (r[ci], "%s-%s" % (r[si], r[ei]), va, vb, p,
                 "YES" if flip else "no", "  " + _sigtag if sig else ""))
        out_rows.append((r[ci], r[si], r[ei], "%.6f" % va, "%.6f" % vb,
                         "%.6g" % p, "1" if flip else "0",
                         "1" if sig else "0"))

    n_flip = sum(1 for r in out_rows if r[6] == "1")
    n_sig = sum(1 for r in out_rows if r[7] == "1")
    _sigword = "global-significant (panel-level)" if multi_condition else "significant"
    print("\n  VERDICT: %d/%d bin(s) change %s-vs-%s compartment sign; "
          "%d/%d %s at %s <= %g"
          % (n_flip, len(hits), a.exp_a, a.exp_b, n_sig, len(hits),
             _sigword, padj_label, a.alpha))
    if multi_condition:
        print("  (sign change is the pairwise quantity; %s is panel-level -- see NOTE above)"
              % padj_label)
    if n_sig == 0:
        print("  -> no significant compartment change at %s between %s and %s."
              % (a.label, a.exp_a, a.exp_b))
        print("     State this as a bounded negative: with this design, a compartment")
        print("     change of the size dcHiC can detect is not present. It is NOT")
        print("     evidence that nothing changed -- see the power note below.")

    if a.out:
        with open(a.out, "w") as fh:
            fh.write("# dchic_region_result.py\n")
            fh.write("# bedgraph=%s\n# region=%s label=%s\n" % (a.bedgraph, a.region, a.label))
            fh.write("# exp_a=%s exp_b=%s alpha=%s binsize=%d bins=%d\n"
                     % (a.exp_a, a.exp_b, a.alpha, binsize, len(hits)))
            fh.write("# run_experiments=%s multi_condition=%s\n"
                     % (",".join(exp_cols), multi_condition))
            if multi_condition:
                fh.write("# WARNING: padj is a GLOBAL across-panel test (all "
                         "experiments), NOT pairwise %s-vs-%s. sign_change is "
                         "pairwise and reliable; the 'significant' flag reflects "
                         "panel-level padj.\n" % (a.exp_a, a.exp_b))
            _pcol = "padj_global" if multi_condition else "padj"
            _scol = "global_significant" if multi_condition else "significant"
            fh.write("chr\tstart\tend\t%s_PC\t%s_PC\t%s\tsign_change\t%s\n"
                     % (a.exp_a, a.exp_b, _pcol, _scol))
            for r in out_rows:
                fh.write("\t".join(r) + "\n")
        print("\nwrote %s" % a.out)

    print("""
POWER NOTE -- read before writing this up
-----------------------------------------
A null result here is a statement about detectable effect size, not about
biology. Three things bound it:

  1. Bins. A region narrower than one compartment bin cannot produce a
     "regional" compartment result; the test is on one or two bins.
  2. Replicates. dcHiC uses replicate variance as an IHW covariate. With n=2
     per condition the variance estimate is minimal, and the covariate carries
     little information.
  3. Multiple testing. Correction is genome-wide over every bin, so a
     single-bin, single-locus effect must be extreme to survive.

The defensible sentence is of the form: "at <resolution> resolution, with n
replicates per genotype, dcHiC reports no significant compartment change at
<region> (padj = X)". Not: "the compartment does not change".
""")
    return 0


if __name__ == "__main__":
    sys.exit(main())
