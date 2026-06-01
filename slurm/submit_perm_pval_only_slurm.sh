#!/usr/bin/env bash
set -euo pipefail

# Submit only the permutation p-value stage as a Slurm array.
# Assumes perm*/<sample>/*.qcat.bgz already exist.
#
# Usage:
#   ./submit_perm_pval_only_slurm.sh [manifest] [n_perms] [perm_prefix] [results_dir] \
#                                    [cpus] [mem] [time] [null_subsample]
#
# Example:
#   ./submit_perm_pval_only_slurm.sh .slurm_perm_manifest.tsv 4 perm results_march 2 32G 03:00:00 1500000

MANIFEST="${1:-.slurm_perm_manifest.tsv}"
N_PERMS="${2:-4}"
PERM_PREFIX="${3:-perm}"
RESULTS_DIR="${4:-results_march}"

PVAL_CPUS="${5:-2}"
PVAL_MEM="${6:-32G}"
PVAL_TIME="${7:-03:00:00}"
NULL_SUBSAMPLE="${8:-1500000}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BEARING_PVALUE="$SCRIPT_DIR/bearing_pvalue.py"

if [[ ! -f "$MANIFEST" ]]; then
  echo "ERROR: manifest not found: $MANIFEST" >&2
  exit 1
fi
if [[ ! -f "$BEARING_PVALUE" ]]; then
  echo "ERROR: missing helper script: $BEARING_PVALUE" >&2
  exit 1
fi

N_SAMPLES=$(wc -l < "$MANIFEST" | tr -d ' ')
if [[ "$N_SAMPLES" -lt 1 ]]; then
  echo "ERROR: manifest has no rows: $MANIFEST" >&2
  exit 1
fi

mkdir -p "$RESULTS_DIR" .slurm_logs

JOBID=$(sbatch --parsable \
  --job-name=perm_pval_only \
  --cpus-per-task="$PVAL_CPUS" \
  --mem="$PVAL_MEM" \
  --time="$PVAL_TIME" \
  --array="1-${N_SAMPLES}" \
  --output=".slurm_logs/perm_pval_only_%A_%a.out" \
  --export=ALL,MANIFEST="$MANIFEST",N_PERMS="$N_PERMS",PERM_PREFIX="$PERM_PREFIX",RESULTS_DIR="$RESULTS_DIR",BEARING_PVALUE="$BEARING_PVALUE",NULL_SUBSAMPLE="$NULL_SUBSAMPLE" \
  --wrap='
set -euo pipefail

line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MANIFEST")
sample=$(printf "%s\n" "$line" | cut -f1)
out=$(printf "%s\n" "$line" | cut -f3)

if [[ -z "$sample" || -z "$out" ]]; then
  echo "ERROR: malformed manifest row $SLURM_ARRAY_TASK_ID in $MANIFEST" >&2
  exit 2
fi

out_name=${out##*/}
base=${out_name%.qcat.bgz}

resolve_null() {
  local perm="$1"
  local dir="${PERM_PREFIX}${perm}/${sample}"
  local a="${dir}/${base}_perm${perm}.qcat.bgz"
  local b="${dir}/null_perm${perm}.qcat.bgz"
  local c="${dir}/${sample}_perm${perm}.qcat.bgz"

  if [[ -f "$a" ]]; then echo "$a"; return 0; fi
  if [[ -f "$b" ]]; then echo "$b"; return 0; fi
  if [[ -f "$c" ]]; then echo "$c"; return 0; fi

  mapfile -t found < <(find "$dir" -maxdepth 1 -name "*.qcat.bgz" | sort)
  if [[ ${#found[@]} -eq 1 ]]; then
    echo "${found[0]}"
    return 0
  fi

  return 1
}

nulls=()
for ((p=1; p<=N_PERMS; p++)); do
  q=$(resolve_null "$p") || {
    echo "ERROR: missing null qcat sample=$sample perm=$p" >&2
    exit 3
  }
  nulls+=("$q")
done

mkdir -p "$RESULTS_DIR"
cmd=(
  python "$BEARING_PVALUE"
  --qcat "$out"
  --null-qcat "${nulls[@]}"
  --out-prefix "$RESULTS_DIR/$sample"
  --fdr 0.05
  --score-plot
)

if [[ -n "${NULL_SUBSAMPLE:-}" && "${NULL_SUBSAMPLE}" != "0" ]]; then
  cmd+=(--null-subsample "$NULL_SUBSAMPLE")
fi

"${cmd[@]}"
')

echo "Submitted pvalue-only array."
echo "  Job ID: $JOBID"
echo "  Array: 1-$N_SAMPLES"
echo "  Null subsample per file: $NULL_SUBSAMPLE"
echo

echo "Monitor with:"
echo "  squeue -j $JOBID"
echo "  tail -f .slurm_logs/perm_pval_only_${JOBID}_*.out"
