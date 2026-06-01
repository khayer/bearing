#!/usr/bin/env bash
#SBATCH --job-name=bearing_calib
#SBATCH --output=bearing_calib_%j.out
#SBATCH --error=bearing_calib_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=96G
#
# sbatch_replicate_calibration.sh
# ===============================
# SLURM wrapper for run_replicate_calibration.sh (genome-wide).
#
# WHY THE BIG --mem
#   compare_qcat.py parses each qcat genome-wide and (in the current build)
#   holds the raw six-track arrays per bin too (include_raw=True is hardcoded).
#   At 200 bp over mm10 that is ~1.3e7 bins per qcat as Python objects, several
#   GB each, two parsed at once in diff-only mode. bearing_calibration.py then
#   pools all permuted diff nulls (n_null = n_bins * n_perms). 96G is a safe
#   starting request for genome-wide + 3 perms; see the SIZING note below to
#   right-size from the first run's seff/sacct.
#
# WHAT IT RUNS
#   The genome-wide within-condition replicate-differential FDR calibration
#   across all conditions, reusing the pre-computed per-sample permutation
#   qcats. No re-scoring; this is diff + p-value only.
#
# SUBMIT
#   sbatch sbatch_replicate_calibration.sh
#
#   Override defaults at submit time via --export, e.g. fewer conditions or
#   more perms:
#     sbatch --export=ALL,N_PERMS=10,CONDITIONS="DN EbKO" \
#            sbatch_replicate_calibration.sh
#   Or change SLURM resources on the command line (these win over the
#   #SBATCH lines above):
#     sbatch --mem=128G --time=24:00:00 sbatch_replicate_calibration.sh
#
# SIZING (do this after the first run)
#   seff $SLURM_JOB_ID
#   sacct -j $SLURM_JOB_ID --format=JobID,State,ExitCode,MaxRSS,ReqMem,Elapsed,MaxVMSize
#   If MaxRSS << ReqMem, lower --mem next time; if the job was killed (OOM),
#   raise it. A bare "Killed" with a PID in the .err file is an OOM, not a
#   timeout (a timeout says CANCELLED ... DUE TO TIME LIMIT).

set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration -- edit these to your paths, or override via --export at submit
# ---------------------------------------------------------------------------
SCRIPTS_DIR="${SCRIPTS_DIR:-$HOME/data/tools/scripts/bigwig_to_qcat}"
RUNNER="${RUNNER:-$SCRIPTS_DIR/run_replicate_calibration.sh}"

SHEET="${SHEET:-samples_calib.tsv}"
OUT="${OUT:-calib_run}"
# Default template kept out of a ${:-...} expansion: the literal {sample}/{perm}
# tokens would otherwise be mangled by brace handling. Set PERM_TEMPLATE via
# --export to override.
DEFAULT_PERM_TEMPLATE='perm_v6_paper{perm}/{sample}/{sample}_perm{perm}.qcat.bgz'
PERM_TEMPLATE="${PERM_TEMPLATE:-$DEFAULT_PERM_TEMPLATE}"
N_PERMS="${N_PERMS:-3}"
FDR="${FDR:-0.05}"
MIN_SIGNAL="${MIN_SIGNAL:-0.0}"
OBS_QCAT_COL="${OBS_QCAT_COL:-out}"
CONDITIONS="${CONDITIONS:-}"          # empty = all conditions in the sheet
PYTHON_BIN="${PYTHON_BIN:-python}"
WORKDIR="${WORKDIR:-$SLURM_SUBMIT_DIR}"

# ---------------------------------------------------------------------------
# Optional environment module / conda activation. Uncomment and adapt.
# ---------------------------------------------------------------------------
# module load python/3.11
# source ~/miniconda3/etc/profile.d/conda.sh && conda activate bearing

# ---------------------------------------------------------------------------
# Preamble
# ---------------------------------------------------------------------------
echo "============================================================"
echo "  SLURM job        : ${SLURM_JOB_ID:-<interactive>}"
echo "  node             : $(hostname)"
echo "  start            : $(date)"
echo "  workdir          : $WORKDIR"
echo "  scripts-dir      : $SCRIPTS_DIR"
echo "  runner           : $RUNNER"
echo "  sheet            : $SHEET"
echo "  out              : $OUT"
echo "  perm template    : $PERM_TEMPLATE"
echo "  n-perms          : $N_PERMS"
echo "  conditions       : ${CONDITIONS:-<all>}"
echo "  mem requested    : ${SLURM_MEM_PER_NODE:-?} MB (node)"
echo "  cpus             : ${SLURM_CPUS_PER_TASK:-?}"
echo "============================================================"

cd "$WORKDIR"

# ---------------------------------------------------------------------------
# Preflight: fail before burning queue time on a missing file
# ---------------------------------------------------------------------------
fail=0
[[ -f "$RUNNER"        ]] || { echo "ERROR: runner not found: $RUNNER" >&2; fail=1; }
[[ -f "$SHEET"         ]] || { echo "ERROR: sheet not found: $SHEET (cwd=$PWD)" >&2; fail=1; }
for s in compare_qcat.py bearing_pvalue.py bearing_calibration.py; do
    [[ -f "$SCRIPTS_DIR/$s" ]] || { echo "ERROR: missing $SCRIPTS_DIR/$s" >&2; fail=1; }
done
command -v "$PYTHON_BIN" >/dev/null 2>&1 || { echo "ERROR: python not on PATH: $PYTHON_BIN" >&2; fail=1; }
[[ $fail -eq 0 ]] || { echo "Preflight failed; not submitting work." >&2; exit 1; }

# ---------------------------------------------------------------------------
# Assemble runner arguments
# ---------------------------------------------------------------------------
args=(
    --sheet         "$SHEET"
    --out           "$OUT"
    --perm-template "$PERM_TEMPLATE"
    --n-perms       "$N_PERMS"
    --fdr           "$FDR"
    --min-signal    "$MIN_SIGNAL"
    --obs-qcat-col  "$OBS_QCAT_COL"
    --scripts-dir   "$SCRIPTS_DIR"
)
[[ -n "$CONDITIONS" ]] && args+=( --conditions "$CONDITIONS" )

echo "+ PYTHON=$PYTHON_BIN bash $RUNNER ${args[*]}"
echo "------------------------------------------------------------"

# /usr/bin/time -v gives a peak-RSS line in the .err to cross-check seff.
TIME_BIN=""
command -v /usr/bin/time >/dev/null 2>&1 && TIME_BIN="/usr/bin/time -v"

PYTHON="$PYTHON_BIN" $TIME_BIN bash "$RUNNER" "${args[@]}"
rc=$?

echo "------------------------------------------------------------"
echo "  runner exit code : $rc"
echo "  end              : $(date)"
echo "  aggregate table  : $OUT/calibration_summary.tsv"
echo "============================================================"
exit $rc
