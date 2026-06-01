#!/usr/bin/env bash
#SBATCH --job-name=bearing_sweep
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=6:00:00
#SBATCH --output=bearing_sweep_%j.log
# -----------------------------------------------------------------------------
# Submit with:   sbatch run_bearing_recovery_sweep.sh [OUTROOT]
# Run directly:  bash   run_bearing_recovery_sweep.sh [OUTROOT]
#
# Edit the two lines below for your environment (conda install path + env name,
# and the working directory) before submitting.
# -----------------------------------------------------------------------------
#if [ -n "${SLURM_JOB_ID:-}" ]; then
#    # running under SLURM: set up the environment
#    source ~/miniconda3/etc/profile.d/conda.sh
#    conda activate bearing
#    cd /mnt/isilon/bassing_lab/integration_paper/bearing_sim_test
#fi
# =============================================================================
# run_bearing_recovery_sweep.sh
# -----------------------------------------------------------------------------
# COARSE GRID robustness sweep for the BEARING ground-truth test.
#
# Sweeps every COMBINATION of a small grid over:
#   SIMULATION parameters (require re-simulating the synthetic data):
#     --snr            signal-to-noise of the planted signal
#     --nb-dispersion  negative-binomial dispersion (smaller = noisier)
#   SCORING parameter (changes only how BEARING scores fixed data):
#     --normalize-tracks   per-track (nonzero-quantile) normalization, off/on
#
# WHY THIS GRID (changed from the earlier pseudocount sweep):
#   The earlier sweep varied the zero-clamp pseudocount over {1e-7,1e-6,1e-5}
#   and every triple of rows was bit-identical: with six tracks summing to
#   O(1) per bin, a pseudocount that small never moves P where real signal
#   exists, and the low-signal mask already zeros sub-floor cells BEFORE
#   normalization. So that axis only RE-CONFIRMED insensitivity at 3x the
#   compute. We pin the pseudocount to a single representative value (1e-6) and
#   instead vary a knob that can actually move recovery: per-track
#   normalization. --normalize-tracks rescales the six assays to comparable
#   distributions before P is formed, which changes WHICH track dominates each
#   bin's composition -- exactly the thing the attribution recovers.
#
#   The expected clamped-KL allocation in the evaluator is analytic from the
#   planted (a, s, Q), so it does NOT change between the off/on arms; recovered-
#   vs-expected therefore stays a fair faithfulness test under normalization.
#
# Grid size = |SNR| x |DISP| x |NORMTRACKS|. With the defaults below that is
# 3 x 2 x 2 = 12 cells (vs 18 before, and all 12 are informative). For each
# (snr, disp) we simulate ONCE (the sim depends on neither scoring knob), then
# score+evaluate with normalization off and on.
#
# A master TSV collects, per cell: per-block-vs-expected slope and Pearson,
# rank recovery, and the informative-block RAD21+CTCF fraction. That table is
# the headline robustness result -- recovery should be stable across the grid,
# and the off-vs-on comparison shows whether track normalization preserves
# faithful attribution.
#
# OPTIONAL pseudocount bound (PSEUDO_PROBE=1): one extra off-grid run at a
# LARGE pseudocount on a single (snr, disp) cell, to bound where the
# pseudocount actually starts to bite. Turns "insensitive over 1e-7..1e-5"
# into "insensitive up to ~1e-3, degrades beyond" -- a stronger caveat
# sentence. Off by default.
#
# USAGE
#   bash run_bearing_recovery_sweep.sh [OUTROOT]
#
# Activate the env (numpy scipy matplotlib pyBigWig pysam) BEFORE running.
# =============================================================================
set -euo pipefail

# ---- paths to the scripts (edit if not in the current directory) ------------
SIM=./simulate_bearing_tracks.py
EVAL=./evaluate_bearing_recovery.py
SCORER=~/data/tools/scripts/bigwig_to_qcat/bigwig_to_qcat.py
VIZ=./visualize_bearing_tracks.py

PY=${PYTHON:-python3}
OUTROOT=${1:-bearing_sweep_new}
CHROM=chrSim
LENGTH=20000000          # 20 Mb full run; set 2000000 for a smoke test
SEED=12345

# ---- COARSE grids -----------------------------------------------------------
SNR_GRID="20 40 80"            # 3 values
DISP_GRID="4 16"               # 2 values
NORMTRACKS_GRID="off on"       # 2 values: per-track normalization off / on

# Pseudocount is PINNED (no longer swept); representative mid value.
PSEUDOCOUNT=1e-6

# Optional one-shot large-pseudocount probe (see header). Set PSEUDO_PROBE=1
# to enable; it runs on the cell named by PROBE_SNR/PROBE_DISP with normtracks
# OFF, at PROBE_PSEUDO.
PSEUDO_PROBE=${PSEUDO_PROBE:-0}
PROBE_SNR=${PROBE_SNR:-20}
PROBE_DISP=${PROBE_DISP:-4}
PROBE_PSEUDO=${PROBE_PSEUDO:-1e-2}

mkdir -p "$OUTROOT"
MASTER="$OUTROOT/sweep_master.tsv"
echo -e "snr\tdisp\tnormalize_tracks\tpseudocount\toutdir\tperblock_slope_vs_expected\tperblock_pearson_vs_expected\trank_recovery\tinformative_frac_min\tinformative_frac_max" > "$MASTER"

BW () { local p=$1; echo "${p}_atac.bw ${p}_rnaplus.bw ${p}_rnaminus.bw ${p}_ctcf.bw ${p}_rad21.bw ${p}_h3k27ac.bw"; }

# pull headline numbers from a recovery_summary.txt into the master TSV
collect () {
  local snr=$1 disp=$2 normtracks=$3 pseudo=$4 outdir=$5
  "$PY" - "$snr" "$disp" "$normtracks" "$pseudo" "$outdir" "$outdir/recovery_summary.txt" "$MASTER" <<'PYEOF'
import sys, re
snr, disp, normtracks, pseudo, outdir, sumf, master = sys.argv[1:8]
slope=pear=rank=fmin=fmax="NA"
try:
    txt=open(sumf).read()
    m=re.search(r"PER-BLOCK mean vs EXPECTED.*?slope=([0-9.]+)", txt)
    if m: slope=m.group(1)
    m=re.search(r"PER-BLOCK mean vs EXPECTED.*?Pearson=([0-9.]+)", txt)
    if m: pear=m.group(1)
    m=re.search(r"rank recovery\):\s*([0-9.]+)", txt)
    if m: rank=m.group(1)
    m=re.search(r"INFORMATIVE-block RAD21\+CTCF fraction range across all sweeps:\s*([0-9.]+)\s*-\s*([0-9.]+)", txt)
    if m: fmin, fmax = m.group(1), m.group(2)
except FileNotFoundError:
    pass
with open(master,"a") as fh:
    fh.write("\t".join([snr,disp,normtracks,pseudo,outdir,slope,pear,rank,fmin,fmax])+"\n")
PYEOF
}

# score+evaluate one cell; normtracks is "off" or "on"
evaluate_cell () {
  local snr=$1 disp=$2 normtracks=$3 pseudo=$4 prefix=$5 outdir=$6
  local nt_flag=()
  if [ "$normtracks" = "on" ]; then nt_flag=(--normalize-tracks); fi
  # --pseudocount is applied to EVERY scorer run inside the evaluator, so
  # slope/Pearson/rank are genuinely recomputed at this pseudocount and
  # normalization setting (not just panel e).
  "$PY" "$EVAL" \
        --truth "${prefix}_truth.tsv" \
        --bw $(BW "$prefix") \
        --chrom-sizes "${prefix}.chrom.sizes" \
        --scorer "$SCORER" \
        --pseudocount "$pseudo" \
        "${nt_flag[@]}" \
        --outdir "$outdir" >/dev/null
  collect "$snr" "$disp" "$normtracks" "$pseudo" "$outdir"
}

CELL=0
TOTAL=$(( $(wc -w <<<"$SNR_GRID") * $(wc -w <<<"$DISP_GRID") * $(wc -w <<<"$NORMTRACKS_GRID") ))

for snr in $SNR_GRID; do
  for disp in $DISP_GRID; do
    # simulate once per (snr, disp) -- neither scoring knob affects simulation
    prefix="$OUTROOT/sim_snr${snr}_disp${disp}"
    echo ">>> simulate snr=$snr disp=$disp -> $prefix"
    "$PY" "$SIM" --prefix "$prefix" --chrom "$CHROM" --length "$LENGTH" \
          --seed "$SEED" --snr "$snr" --nb-dispersion "$disp" >/dev/null

    for normtracks in $NORMTRACKS_GRID; do
      CELL=$((CELL+1))
      outdir="$OUTROOT/eval_snr${snr}_disp${disp}_nt-${normtracks}"
      echo ">>> [$CELL/$TOTAL] evaluate snr=$snr disp=$disp normalize_tracks=$normtracks"
      evaluate_cell "$snr" "$disp" "$normtracks" "$PSEUDOCOUNT" "$prefix" "$outdir"
    done

    # one overview+zoom track picture per (snr, disp) for the record
    if [ -f "$VIZ" ]; then
      "$PY" "$SCORER" --bw $(BW "$prefix") --out "${prefix}.qcat.bgz" \
            --chrom-sizes "${prefix}.chrom.sizes" --no-extras >/dev/null 2>&1 || true
      "$PY" "$VIZ" --prefix "$prefix" --qcat "${prefix}.qcat.bgz" \
            --out "$OUTROOT/tracks_snr${snr}_disp${disp}.png" >/dev/null 2>&1 || true
    fi
  done
done

# ---- optional large-pseudocount probe (bounds where pseudocount bites) ------
if [ "$PSEUDO_PROBE" = "1" ]; then
  probe_prefix="$OUTROOT/sim_snr${PROBE_SNR}_disp${PROBE_DISP}"
  if [ -f "${probe_prefix}_truth.tsv" ]; then
    probe_out="$OUTROOT/probe_snr${PROBE_SNR}_disp${PROBE_DISP}_nt-off_eps${PROBE_PSEUDO}"
    echo ">>> [probe] large pseudocount snr=$PROBE_SNR disp=$PROBE_DISP eps=$PROBE_PSEUDO (normtracks off)"
    evaluate_cell "$PROBE_SNR" "$PROBE_DISP" "off" "$PROBE_PSEUDO" "$probe_prefix" "$probe_out"
  else
    echo ">>> [probe] skipped: ${probe_prefix}_truth.tsv not found " \
         "(PROBE_SNR/PROBE_DISP must be in the swept grid)" >&2
  fi
fi

echo
echo "=============================================================="
echo "Coarse grid complete ($TOTAL cells). Master: $MASTER"
echo "Per-cell figures + summaries under: $OUTROOT/eval_*/"
echo "Pseudocount pinned at $PSEUDOCOUNT; normalize_tracks swept {off,on}."
[ "$PSEUDO_PROBE" = "1" ] && echo "Large-pseudocount probe row appended (eps=$PROBE_PSEUDO)."
echo "=============================================================="
column -t -s $'\t' "$MASTER" 2>/dev/null || cat "$MASTER"
