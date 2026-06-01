#!/usr/bin/env python3
"""
simulate_bearing_tracks.py
==========================
Ground-truth generator for the BEARING pipeline.

WHY THIS EXISTS
---------------
Biological recapitulation (BEARING reproduces known Tcrb/Igh biology) cannot,
on its own, rule out the objection that the method was tuned to find what was
already known to be there. This script supplies the one thing recapitulation
cannot: a KNOWN per-bin composition that BEARING must recover. We plant a
per-bin probability vector

    P = (1 - s) * Q + s * a

where Q is the per-sample background composition, a is a known allocation
vector over the six tracks (sums to 1) describing how the deviation from Q is
distributed, and s in (0, 1] is a known strength describing how far the bin is
pushed from Q toward a. NULL bins have P = Q and should score low.

The companion evaluator (evaluate_bearing_recovery.py) runs the REAL pipeline
on the BigWigs written here and checks that the recovered per-track allocation
matches the planted a, and that the recovered composition is stable across the
three free parameters (low-signal mask, zero-clamp epsilon/pseudocount,
normalization quantile). That parameter sweep is what answers whether a
quantitative claim such as the 55-71 percent cohesin/CTCF split is a real
compositional feature or an artifact of the mask, clamp, and normalization
choices.

TRACK ORDER (fixed, matches bigwig_to_qcat.py state order 1..6)
    0 ATAC, 1 RNA+, 2 RNA-, 3 CTCF, 4 RAD21 (a.k.a. Cohesin), 5 H3K27ac

REALISM MODEL
-------------
1. Each bin gets a total-coverage factor C drawn from a heavy-tailed
   (log-normal) distribution, so per-bin sequencing depth varies the way real
   coverage does.
2. The expected per-track signal is mu_track = C * snr * P_track (NULL bins use
   P = Q). The snr scalar sets how far above the per-track noise floor the
   planted signal sits.
3. Observed per-track signal is drawn from a negative binomial around mu with a
   tunable dispersion, then a small constant noise floor is added.
4. Each track is smoothed with a Gaussian kernel over a few bins to create
   fragment-length autocorrelation, so the circular-shift null used by BEARING
   has realistic local structure rather than i.i.d. noise. Smoothing bleeds
   composition across bin boundaries; block-edge bins are therefore flagged
   separately from block-interior bins in the truth table so the evaluator can
   score interiors cleanly.

OUTPUTS
-------
    <prefix>_atac.bw, _rnaplus.bw, _rnaminus.bw, _ctcf.bw, _rad21.bw,
        _h3k27ac.bw                 six BigWig tracks
    <prefix>.chrom.sizes            chrom-sizes header (one line)
    <prefix>_truth.tsv              per-bin truth table (see columns below)
    <prefix>_blocks.bed             BED of planted blocks (name = block_id)

TRUTH TABLE COLUMNS
    chrom, start, end, bin_index, block_id, bin_class, strength_s,
    a_ATAC, a_RNAplus, a_RNAminus, a_CTCF, a_RAD21, a_H3K27ac, is_null
  bin_class is one of: null, interior, edge
  For null bins, strength_s = 0 and the a_* columns repeat Q (the background).

DEPENDENCIES
    numpy, scipy, pyBigWig
"""

import argparse
import sys

import numpy as np

try:
    from scipy.ndimage import gaussian_filter1d
except Exception as exc:  # pragma: no cover
    sys.stderr.write("ERROR: scipy is required (scipy.ndimage.gaussian_filter1d)\n")
    raise

# -----------------------------------------------------------------------------
# Fixed six-track layout. Names match the first six categories used by
# bigwig_to_qcat.py (RAD21 is labelled "Cohesin" there). Order is load-bearing:
# the evaluator maps qcat state index 1..6 onto exactly this order.
# -----------------------------------------------------------------------------
TRACK_NAMES = ["ATAC", "RNA+", "RNA-", "CTCF", "RAD21", "H3K27ac"]
N_TRACKS = 6

# Short tokens used in output filenames, in track order.
TRACK_FILE_TOKENS = ["atac", "rnaplus", "rnaminus", "ctcf", "rad21", "h3k27ac"]

# Default mildly-skewed background composition Q over the six tracks
# (sums to 1). Deliberately non-uniform so KL contributions are non-trivial.
DEFAULT_Q = np.array([0.22, 0.14, 0.10, 0.18, 0.20, 0.16], dtype=np.float64)


def _normalize(vec):
    """Return a copy of vec normalized to sum to 1 (safe against zeros)."""
    vec = np.asarray(vec, dtype=np.float64)
    total = vec.sum()
    if total <= 0:
        raise ValueError("cannot normalize a non-positive vector")
    return vec / total


def named_allocations():
    """
    Return the two named allocation vectors mirroring the dual-scale claim,
    plus a dict for reference. Each vector sums to 1, six entries in TRACK
    order [ATAC, RNA+, RNA-, CTCF, RAD21, H3K27ac].

    architecture : RAD21/CTCF dominated (cohesin + insulator), small remainder
    transcription: RNA+ dominated, small remainder
    """
    small = 0.05  # remainder mass spread over the non-dominant tracks

    # architecture: 0.45 RAD21, 0.30 CTCF, remainder small and even
    arch = np.full(N_TRACKS, 0.0)
    arch[4] = 0.45  # RAD21
    arch[3] = 0.30  # CTCF
    rem_idx = [0, 1, 2, 5]
    rem_mass = 1.0 - arch.sum()
    for i in rem_idx:
        arch[i] = rem_mass / len(rem_idx)
    arch = _normalize(arch)

    # transcription: 0.70 RNA+, remainder small and even
    txn = np.full(N_TRACKS, 0.0)
    txn[1] = 0.70  # RNA+
    rem_idx = [0, 2, 3, 4, 5]
    rem_mass = 1.0 - txn.sum()
    for i in rem_idx:
        txn[i] = rem_mass / len(rem_idx)
    txn = _normalize(txn)

    return {"architecture": arch, "transcription": txn}


def grid_allocations():
    """
    A small panel of additional allocation vectors so recovered-vs-planted can
    be characterized over several distinct compositions (not just the two named
    ones). Each sums to 1, six entries in TRACK order.
    """
    allocs = {}

    # single-track-dominant vectors, one per track (tests rank recovery cleanly)
    for ti, name in enumerate(TRACK_NAMES):
        v = np.full(N_TRACKS, 0.04)
        v[ti] = 1.0 - 0.04 * (N_TRACKS - 1)
        allocs["dom_%s" % TRACK_FILE_TOKENS[ti]] = _normalize(v)

    # a couple of two-track mixes
    mix1 = np.array([0.05, 0.05, 0.05, 0.35, 0.45, 0.05])  # CTCF+RAD21 (insulator)
    allocs["mix_ctcf_rad21"] = _normalize(mix1)
    mix2 = np.array([0.30, 0.30, 0.05, 0.05, 0.05, 0.25])  # ATAC+RNA+ +H3K27ac (active)
    allocs["mix_active"] = _normalize(mix2)

    # Non-saturated architecture block: RAD21/CTCF dominant, but ATAC and H3K27ac
    # are ALSO pushed above background so they survive the zero-clamp. This keeps
    # the recovered RAD21+CTCF fraction in an informative ~0.73 range (not pinned
    # at 1.0), so the parameter-sensitivity panel can actually reveal whether the
    # mask/clamp/normalization perturb a NON-trivial composition. The value also
    # brackets the real 55-71% cohesin/CTCF split.
    arch_active = np.array([0.26, 0.05, 0.03, 0.24, 0.30, 0.20])
    allocs["arch_plus_active"] = _normalize(arch_active)

    return allocs


def build_block_plan(q, n_bins, bin_size, block_bins, gap_bins,
                     strengths, start_offset_bins):
    """
    Lay out planted blocks along the chromosome.

    Returns a list of block dicts:
        { "block_id", "alloc" (6,), "strength", "bin_start", "bin_end" }
    bin_end is exclusive. Blocks are placed left to right with `gap_bins` NULL
    bins between them. Raises if the plan does not fit in n_bins.

    Block panel:
      - the two named blocks at full-ish strength (s=0.8) so the named-block
        bar charts are unambiguous
      - the named allocations swept across the strength grid
      - the grid allocations each at a mid strength (s=0.6)
    """
    named = named_allocations()
    grid = grid_allocations()

    plan = []

    # 1. Two named blocks at s = 0.8.
    for name in ("architecture", "transcription"):
        plan.append({"block_id": "named_%s_s0.80" % name,
                     "alloc": named[name], "strength": 0.80})

    # 2. Named allocations across the strength grid.
    for name in ("architecture", "transcription"):
        for s in strengths:
            plan.append({"block_id": "%s_s%.2f" % (name, s),
                         "alloc": named[name], "strength": float(s)})

    # 3. Grid allocations at a single mid strength.
    for key, vec in grid.items():
        plan.append({"block_id": "%s_s0.60" % key,
                     "alloc": vec, "strength": 0.60})

    # 4. The non-saturated arch_plus_active block across the strength grid, so
    #    the parameter-sensitivity panel can sweep a composition whose RAD21+CTCF
    #    fraction is informative (~0.73) rather than pinned at 1.0.
    for s in strengths:
        plan.append({"block_id": "arch_plus_active_s%.2f" % s,
                     "alloc": grid["arch_plus_active"], "strength": float(s)})

    # Assign coordinates.
    cursor = start_offset_bins
    for blk in plan:
        bs = cursor
        be = bs + block_bins
        if be > n_bins:
            raise ValueError(
                "block plan does not fit: need bin %d but only %d bins. "
                "Increase --length or reduce --block-bins / number of blocks."
                % (be, n_bins))
        blk["bin_start"] = bs
        blk["bin_end"] = be
        cursor = be + gap_bins

    return plan


def assemble_P(q, n_bins, plan):
    """
    Build the per-bin planted probability matrix P (n_bins x 6), the per-bin
    block id array, and the bin_class array (null / interior / edge).

    Edge bins are the first and last bin of each block (they will be most
    contaminated by Gaussian smoothing); all other in-block bins are interior.
    """
    P = np.tile(q[np.newaxis, :], (n_bins, 1))  # NULL bins default to Q
    block_id = np.array(["null"] * n_bins, dtype=object)
    bin_class = np.array(["null"] * n_bins, dtype=object)
    strength = np.zeros(n_bins, dtype=np.float64)

    for blk in plan:
        a = blk["alloc"]
        s = blk["strength"]
        bs, be = blk["bin_start"], blk["bin_end"]
        planted = (1.0 - s) * q + s * a
        planted = planted / planted.sum()  # guard against fp drift
        P[bs:be, :] = planted[np.newaxis, :]
        block_id[bs:be] = blk["block_id"]
        strength[bs:be] = s
        bin_class[bs:be] = "interior"
        # mark the first and last bin of the block as edge
        bin_class[bs] = "edge"
        bin_class[be - 1] = "edge"

    return P, block_id, bin_class, strength


def simulate_signal(P, rng, coverage_log_mean, coverage_log_sigma,
                    snr, nb_dispersion, noise_floor, smooth_sigma_bins):
    """
    Turn the planted probability matrix P (n_bins x 6) into observed per-track
    signal (n_bins x 6), applying the realism model described in the header.

    Parameters
    ----------
    P                 : (n_bins, 6) planted probability matrix
    rng               : numpy Generator
    coverage_log_mean : mean of the log-normal per-bin coverage factor
    coverage_log_sigma: sigma of the log-normal per-bin coverage factor
    snr               : scalar multiplying the per-track expectation
    nb_dispersion     : negative-binomial dispersion (size r); smaller = noisier
    noise_floor       : constant added to every track after NB draw
    smooth_sigma_bins : Gaussian smoothing sigma (in bins) per track

    Returns
    -------
    signal : (n_bins, 6) float array of observed signal (>= 0)
    """
    n_bins = P.shape[0]

    # Per-bin total coverage factor (heavy-tailed, positive).
    coverage = rng.lognormal(mean=coverage_log_mean, sigma=coverage_log_sigma,
                             size=n_bins)  # (n_bins,)

    # Expected per-track signal: coverage * snr * P.
    mu = coverage[:, np.newaxis] * snr * P  # (n_bins, 6)
    mu = np.clip(mu, 1e-9, None)

    # Negative-binomial draw around mu. Parameterize NB by (r, p) with
    # mean = r (1 - p) / p. Fix r = nb_dispersion, solve p = r / (r + mu).
    r = float(nb_dispersion)
    p = r / (r + mu)
    # numpy's negative_binomial(n, p) returns number of failures before n
    # successes; its mean is n (1 - p) / p, matching the parameterization.
    signal = rng.negative_binomial(r, p).astype(np.float64)

    # Add a small constant noise floor so no track is exactly zero everywhere.
    signal = signal + noise_floor

    # Per-track Gaussian smoothing for fragment-length autocorrelation.
    if smooth_sigma_bins and smooth_sigma_bins > 0:
        for ti in range(signal.shape[1]):
            signal[:, ti] = gaussian_filter1d(
                signal[:, ti], sigma=float(smooth_sigma_bins), mode="wrap")

    signal = np.clip(signal, 0.0, None)
    return signal


def write_bigwigs(signal, chrom, chrom_len, bin_size, prefix):
    """
    Write six BigWig files, one per track, using fixed-width 200 bp intervals.
    Requires pyBigWig. Returns the list of written paths in track order.
    """
    try:
        import pyBigWig
    except ImportError:
        sys.stderr.write(
            "ERROR: pyBigWig is required to write BigWigs. "
            "Install with: pip install pyBigWig\n")
        raise

    n_bins = signal.shape[0]
    starts = (np.arange(n_bins, dtype=np.int64) * bin_size).astype(np.int64)
    ends = starts + bin_size
    ends[-1] = min(int(ends[-1]), chrom_len)  # clamp final interval to chrom end

    paths = []
    for ti, token in enumerate(TRACK_FILE_TOKENS):
        path = "%s_%s.bw" % (prefix, token)
        bw = pyBigWig.open(path, "w")
        bw.addHeader([(chrom, int(chrom_len))])
        vals = signal[:, ti].astype(np.float64)
        bw.addEntries(
            [chrom] * n_bins,
            [int(x) for x in starts],
            ends=[int(x) for x in ends],
            values=[float(v) for v in vals],
        )
        bw.close()
        paths.append(path)
    return paths


def write_chrom_sizes(chrom, chrom_len, prefix):
    path = "%s.chrom.sizes" % prefix
    with open(path, "w") as fh:
        fh.write("%s\t%d\n" % (chrom, int(chrom_len)))
    return path


def write_truth_tsv(prefix, chrom, bin_size, P, block_id, bin_class, strength,
                    is_null):
    path = "%s_truth.tsv" % prefix
    header = (["chrom", "start", "end", "bin_index", "block_id", "bin_class",
               "strength_s"]
              + ["a_%s" % t for t in
                 ["ATAC", "RNAplus", "RNAminus", "CTCF", "RAD21", "H3K27ac"]]
              + ["is_null"])
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        n_bins = P.shape[0]
        for i in range(n_bins):
            start = i * bin_size
            end = start + bin_size
            a = P[i]
            row = [chrom, str(start), str(end), str(i),
                   str(block_id[i]), str(bin_class[i]),
                   "%.6f" % strength[i]]
            row += ["%.6f" % a[j] for j in range(N_TRACKS)]
            row += ["1" if is_null[i] else "0"]
            fh.write("\t".join(row) + "\n")
    return path


def write_blocks_bed(prefix, chrom, bin_size, plan):
    path = "%s_blocks.bed" % prefix
    with open(path, "w") as fh:
        for blk in plan:
            start = blk["bin_start"] * bin_size
            end = blk["bin_end"] * bin_size
            # BED: chrom start end name score strand
            fh.write("%s\t%d\t%d\t%s\t%d\t.\n"
                     % (chrom, start, end, blk["block_id"],
                        int(round(1000 * blk["strength"]))))
    return path


def parse_q(arg, default):
    if arg is None:
        return _normalize(default)
    parts = [p for p in arg.replace(",", " ").split() if p]
    if len(parts) != N_TRACKS:
        raise ValueError("--background-q needs exactly %d values (got %d)"
                         % (N_TRACKS, len(parts)))
    return _normalize(np.array([float(x) for x in parts], dtype=np.float64))


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Generate synthetic six-track BigWigs with a known per-bin "
                    "composition for BEARING ground-truth recovery testing.")
    ap.add_argument("--prefix", default="bearing_sim",
                    help="output filename prefix (default: bearing_sim)")
    ap.add_argument("--chrom", default="chrSim",
                    help="synthetic chromosome name (default: chrSim)")
    ap.add_argument("--length", type=int, default=20_000_000,
                    help="chromosome length in bp (default: 20000000 = 20 Mb)")
    ap.add_argument("--bin-size", type=int, default=200,
                    help="bin size in bp (default: 200)")
    ap.add_argument("--background-q", default=None,
                    help="six comma- or space-separated values for Q "
                         "(default: mildly skewed %s)"
                         % ",".join("%.2f" % x for x in DEFAULT_Q))
    ap.add_argument("--block-bins", type=int, default=25,
                    help="number of bins per planted block (default: 25 = 5 kb)")
    ap.add_argument("--gap-bins", type=int, default=75,
                    help="NULL bins between planted blocks (default: 75)")
    ap.add_argument("--start-offset-bins", type=int, default=500,
                    help="leading NULL bins before the first block (default: 500)")
    ap.add_argument("--strength-grid", default="0.1,0.2,0.3,0.5,0.7,0.9",
                    help="comma-separated strengths s for the named-block sweep")
    # realism knobs
    ap.add_argument("--snr", type=float, default=40.0,
                    help="signal-to-noise scalar on per-track expectation "
                         "(default: 40)")
    ap.add_argument("--coverage-log-mean", type=float, default=0.0,
                    help="log-normal mean for per-bin coverage (default: 0.0)")
    ap.add_argument("--coverage-log-sigma", type=float, default=0.5,
                    help="log-normal sigma for per-bin coverage (default: 0.5)")
    ap.add_argument("--nb-dispersion", type=float, default=8.0,
                    help="negative-binomial dispersion r; smaller = noisier "
                         "(default: 8.0)")
    ap.add_argument("--noise-floor", type=float, default=0.5,
                    help="constant added to every track post-NB (default: 0.5)")
    ap.add_argument("--smooth-sigma-bins", type=float, default=1.5,
                    help="Gaussian smoothing sigma in bins for autocorrelation "
                         "(default: 1.5)")
    ap.add_argument("--seed", type=int, default=12345,
                    help="random seed (default: 12345)")
    ap.add_argument("--no-bigwig", action="store_true",
                    help="skip BigWig writing (truth table + BED only); useful "
                         "for testing without pyBigWig")
    args = ap.parse_args(argv)

    if args.length % args.bin_size != 0:
        sys.stderr.write(
            "WARNING: length %d not divisible by bin-size %d; trailing bp "
            "ignored.\n" % (args.length, args.bin_size))

    n_bins = args.length // args.bin_size
    q = parse_q(args.background_q, DEFAULT_Q)
    strengths = [float(x) for x in args.strength_grid.split(",") if x.strip()]
    rng = np.random.default_rng(args.seed)

    plan = build_block_plan(
        q, n_bins, args.bin_size,
        block_bins=args.block_bins, gap_bins=args.gap_bins,
        strengths=strengths, start_offset_bins=args.start_offset_bins)

    P, block_id, bin_class, strength = assemble_P(q, n_bins, plan)
    is_null = (bin_class == "null")

    signal = simulate_signal(
        P, rng,
        coverage_log_mean=args.coverage_log_mean,
        coverage_log_sigma=args.coverage_log_sigma,
        snr=args.snr,
        nb_dispersion=args.nb_dispersion,
        noise_floor=args.noise_floor,
        smooth_sigma_bins=args.smooth_sigma_bins)

    # outputs
    cs_path = write_chrom_sizes(args.chrom, args.length, args.prefix)
    truth_path = write_truth_tsv(args.prefix, args.chrom, args.bin_size,
                                 P, block_id, bin_class, strength, is_null)
    bed_path = write_blocks_bed(args.prefix, args.chrom, args.bin_size, plan)

    bw_paths = []
    if not args.no_bigwig:
        bw_paths = write_bigwigs(signal, args.chrom, args.length,
                                 args.bin_size, args.prefix)

    # summary
    n_planted = int((~is_null).sum())
    print("simulate_bearing_tracks.py summary")
    print("  chromosome      : %s  (%d bp, %d bins of %d bp)"
          % (args.chrom, args.length, n_bins, args.bin_size))
    print("  background Q    : %s"
          % ", ".join("%s=%.3f" % (TRACK_NAMES[i], q[i]) for i in range(N_TRACKS)))
    print("  planted blocks  : %d  (%d planted bins, %d null bins)"
          % (len(plan), n_planted, n_bins - n_planted))
    print("  block size      : %d bins (%d bp), gap %d bins"
          % (args.block_bins, args.block_bins * args.bin_size, args.gap_bins))
    print("  named blocks    : architecture (RAD21/CTCF), transcription (RNA+)")
    print("  strength grid   : %s" % strengths)
    print("  chrom.sizes     : %s" % cs_path)
    print("  truth table     : %s" % truth_path)
    print("  blocks BED      : %s" % bed_path)
    if bw_paths:
        print("  BigWigs         : %s" % ", ".join(bw_paths))
    else:
        print("  BigWigs         : (skipped, --no-bigwig)")


if __name__ == "__main__":
    main()
