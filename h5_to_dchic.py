#!/usr/bin/env python3
"""
h5_to_dchic.py -- convert HiCExplorer .h5 matrices to dcHiC sparse-matrix input.

WHY
---
dcHiC takes HiC-Pro style sparse matrices (<indexA> <indexB> <count>) plus a bed
file. It does NOT read .h5, .cool or .hic directly, and it does NOT accept
precomputed PC1 tracks -- it computes compartment scores itself, because its
statistics depend on doing so (MFA across all input maps, quantile normalization
of component scores, then a Mahalanobis distance per bin).

This reads .h5 directly through the `hicmatrix` package that ships with
HiCExplorer, so there is no intermediate .cool file.

TWO THINGS TO SETTLE BEFORE YOU TRUST THE OUTPUT
------------------------------------------------
1. BALANCED vs RAW. Files named `raw_corrected_KR_*` are KR-BALANCED, i.e. the
   values are corrected floats, not integer counts. dcHiC's input spec calls the
   third column "count". Compartment calling normally wants a balanced matrix
   (HiCExplorer's own hicPCA does), so balanced is probably right here -- but
   dcHiC also performs internal processing, and feeding it pre-balanced values
   may double-normalize. CHECK the dcHiC wiki/issues for their recommendation
   before publishing a number. This script reports which kind it detected and
   will not silently pretend floats are counts (see --value-mode).

2. RESOLUTION. Compartments are a large-scale property. 100 kb is the usual
   primary choice; 250 kb is a sensible robustness check. Do not run this at
   10 kb for a compartment question.

USAGE
-----
  # one call per replicate
  python h5_to_dchic.py \
      --h5 /path/raw_corrected_KR_dV1P_rep1_bs_100000.h5 \
      --sample dV1P_rep1 --experiment dV1P --outdir dchic_in_100kb

  # then build the dcHiC input file
  python h5_to_dchic.py --write-input dchic_in_100kb \
      --out dchic_in_100kb/input.txt --res-label 100kb

Requires: hicmatrix (installed with HiCExplorer), numpy, scipy. ASCII only.
"""

import argparse
import os
import re
import sys

# dcHiC's README: drop chrM, chrY and non-standard contigs; they break the PCA step.
DEFAULT_KEEP = ["chr%d" % i for i in range(1, 20)] + ["chrX"]


def load_h5(path):
    try:
        from hicmatrix import HiCMatrix as hm
    except ImportError:
        sys.exit("ERROR: the 'hicmatrix' package is required (ships with HiCExplorer).\n"
                 "  conda activate <your hicexplorer env>\n"
                 "  -- or convert first:\n"
                 "     hicConvertFormat --matrices in.h5 --outFileName out.cool \\\n"
                 "         --inputFormat h5 --outputFormat cool")
    print("  loading %s" % path)
    return hm.hiCMatrix(path)


def convert(h5_path, sample, experiment, outdir, keep_chroms, value_mode, blacklist_bed=None):
    import numpy as np
    from scipy.sparse import triu

    os.makedirs(outdir, exist_ok=True)
    hic = load_h5(h5_path)
    ci = hic.cut_intervals
    if not ci:
        sys.exit("ERROR: no cut_intervals in %s" % h5_path)
    binsize = int(ci[0][2]) - int(ci[0][1])
    print("    bins: %s | binsize: %s bp" % (format(len(ci), ","), format(binsize, ",")))
    if binsize < 50000:
        print("    WARNING: binsize %d bp is fine for a compartment analysis; "
              "100 kb is the usual primary resolution." % binsize, file=sys.stderr)

    keep = set(keep_chroms)
    orig2idx = {}
    kept = []
    for i, iv in enumerate(ci):
        chrom, start, end = str(iv[0]), int(iv[1]), int(iv[2])
        if chrom in keep:
            orig2idx[i] = len(kept) + 1          # dcHiC indices are 1-based
            kept.append((chrom, start, end))
    if not kept:
        sys.exit("ERROR: no bins on the requested chromosomes. Present: %s"
                 % sorted({str(x[0]) for x in ci})[:10])
    print("    kept %s bins on %d chromosomes" % (format(len(kept), ","),
                                                  len({k[0] for k in kept})))

    bl = set()
    if blacklist_bed:
        iv = {}
        with open(blacklist_bed) as fh:
            for line in fh:
                if line.startswith(("#", "track")):
                    continue
                f = line.split()
                iv.setdefault(f[0], []).append((int(f[1]), int(f[2])))
        for n, (chrom, s, e) in enumerate(kept, start=1):
            for bs, be in iv.get(chrom, []):
                if s < be and e > bs:
                    bl.add(n)
                    break
        print("    blacklisted %d bins" % len(bl))

    bed_path = os.path.join(outdir, "%s.bed" % sample)
    with open(bed_path, "w") as fh:
        for n, (chrom, s, e) in enumerate(kept, start=1):
            if blacklist_bed:
                fh.write("%s\t%d\t%d\t%d\t%d\n" % (chrom, s, e, n, 1 if n in bl else 0))
            else:
                fh.write("%s\t%d\t%d\t%d\n" % (chrom, s, e, n))

    chrom_of = {i: str(iv[0]) for i, iv in enumerate(ci)}
    m = triu(hic.matrix.tocoo(), k=0).tocoo()   # upper triangle incl. diagonal

    # decide how to render values
    vals = m.data
    is_float = not np.allclose(vals, np.round(vals))
    if value_mode == "auto":
        mode = "float" if is_float else "int"
    else:
        mode = value_mode
    print("    matrix values look like %s; writing as %s"
          % ("BALANCED floats" if is_float else "integer counts", mode))
    if is_float and mode == "int":
        print("    WARNING: rounding balanced values to integers loses information.",
              file=sys.stderr)

    mat_path = os.path.join(outdir, "%s.matrix" % sample)
    n_written = n_trans = 0
    with open(mat_path, "w") as fh:
        for a, b, v in zip(m.row, m.col, m.data):
            if v == 0:
                continue
            ia, ib = orig2idx.get(a), orig2idx.get(b)
            if ia is None or ib is None:
                continue
            if chrom_of[a] != chrom_of[b]:
                n_trans += 1
                continue                       # cis only
            if ia > ib:
                ia, ib = ib, ia
            if mode == "int":
                iv_ = int(round(float(v)))
                if iv_ <= 0:
                    continue
                fh.write("%d\t%d\t%d\n" % (ia, ib, iv_))
            else:
                fh.write("%d\t%d\t%.6g\n" % (ia, ib, float(v)))
            n_written += 1
    print("    wrote %s (%s cis pixels; %s trans pixels dropped)"
          % (os.path.basename(mat_path), format(n_written, ","), format(n_trans, ",")))
    print("    wrote %s" % os.path.basename(bed_path))
    return mat_path, bed_path


def derive_experiment(sample):
    """dV1P_rep1 -> dV1P ; WT_R2 -> WT ; ProB_rep10 -> ProB"""
    m = re.match(r'^(.*?)[_-](?:rep|Rep|REP|R)\d+$', sample)
    return m.group(1) if m else sample


def write_input(indir, out, res_label):
    rows = []
    for f in sorted(os.listdir(indir)):
        if not f.endswith(".matrix"):
            continue
        sample = f[:-len(".matrix")]
        bed = os.path.join(indir, sample + ".bed")
        if not os.path.exists(bed):
            print("  WARNING: no .bed for %s, skipping" % sample, file=sys.stderr)
            continue
        exp = derive_experiment(sample)
        rep_prefix = "%s_%s" % (sample, res_label)
        for bad in ("-", "."):
            if bad in rep_prefix or bad in exp:
                sys.exit("ERROR: dcHiC forbids '%s' in replicate/experiment prefixes.\n"
                         "  offending: %s / %s" % (bad, rep_prefix, exp))
        rows.append((os.path.join(indir, f), bed, rep_prefix, exp))
    if not rows:
        sys.exit("ERROR: no .matrix files in %s" % indir)
    with open(out, "w") as fh:
        for r in rows:
            fh.write("\t".join(r) + "\n")
    print("wrote %s (%d rows)\n" % (out, len(rows)))
    from collections import Counter
    per_exp = Counter(r[3] for r in rows)
    print("replicates per experiment:")
    single = False
    for k, v in sorted(per_exp.items()):
        flag = ""
        if v < 2:
            flag = "   <-- ONE REPLICATE: dcHiC drops to plain FDR, no IHW"
            single = True
        print("   %-16s %d%s" % (k, v, flag))
    if single:
        print("\nNOTE: with a single replicate in any condition, dcHiC cannot use "
              "replicate\nvariance as an IHW covariate for that condition. The run "
              "still works; the\nmultiple-testing correction is simply weaker.")
    if len(per_exp) < 2:
        print("\nERROR-ish: only one experiment found. dcHiC needs at least two "
              "conditions to\ncompare -- add the control (e.g. WT) before running.",
              file=sys.stderr)
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--h5", help="HiCExplorer .h5 matrix")
    ap.add_argument("--sample", help="replicate name, e.g. dV1P_rep1 (no dots/dashes)")
    ap.add_argument("--experiment", default=None,
                    help="condition name, e.g. dV1P (default: derived from --sample)")
    ap.add_argument("--outdir", default="dchic_in")
    ap.add_argument("--chroms", default=",".join(DEFAULT_KEEP))
    ap.add_argument("--blacklist", default=None)
    ap.add_argument("--value-mode", choices=["auto", "int", "float"], default="auto",
                    help="how to write the third matrix column. 'auto' keeps floats "
                         "for balanced matrices and integers for raw counts.")
    ap.add_argument("--write-input", metavar="DIR")
    ap.add_argument("--out", default="input.txt")
    ap.add_argument("--res-label", default="100kb")
    a = ap.parse_args()

    if a.write_input:
        write_input(a.write_input, a.out, a.res_label)
        return 0
    if not (a.h5 and a.sample):
        ap.error("--h5 and --sample are required (or use --write-input)")
    exp = a.experiment or derive_experiment(a.sample)
    for bad in ("-", "."):
        if bad in a.sample or bad in exp:
            sys.exit("ERROR: dcHiC forbids '%s' in replicate/experiment prefixes." % bad)
    print("sample=%s experiment=%s" % (a.sample, exp))
    convert(a.h5, a.sample, exp, a.outdir,
            [c for c in a.chroms.split(",") if c], a.value_mode, a.blacklist)
    return 0


if __name__ == "__main__":
    sys.exit(main())
