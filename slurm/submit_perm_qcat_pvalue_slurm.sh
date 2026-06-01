#!/usr/bin/env bash
set -euo pipefail

# submit_perm_qcat_pvalue_slurm.sh  (v3)
# ============================================================
# Full permutation null pipeline as SLURM array jobs.
#
# Replaces serial execution in generate_perm_nulls.py with four
# Slurm-array stages chained via --dependency=afterok:
#
#   Stage 1: array over (perm, sample)    -> shift_bigwig.py    [NEW in v3]
#   Stage 2: array over (perm, sample)    -> bigwig_to_qcat.py  [existed in v2]
#   Stage 3: array over (perm)            -> compare_qcat.py    [NEW in v3, optional]
#   Stage 4: array over sample            -> bearing_pvalue.py  [existed in v2]
#   Stage 5: array over diff_comparison   -> bearing_pvalue.py  [NEW in v3, optional]
#
# Stage 5 replaces serial run_perm_diff_pvalue.sh by running each
# differential-comparison p-value calculation as an independent SLURM
# task. Only runs when --diff-sheet was provided AND a directory of
# real-data diff qcats is available (default: ./comparison_main; can be
# overridden with COMPARISON_DIR_OVERRIDE environment variable).
#
# Each (perm, sample) pair runs as an independent SLURM task, so the
# cluster scheduler runs as many in parallel as resources allow. This
# is dramatically faster than the serial outer-perm loop in
# generate_perm_nulls.py.
#
# Floors are loaded from --floors-tsv (precomputed once on real BigWigs),
# skipping the percentile-floor sampling pass per perm. Floors are
# invariant under circular shift so reuse is statistically sound.
#
# Usage:
#   ./submit_perm_qcat_pvalue_slurm.sh [sheet] [n_perms] [perm_prefix] [results_dir]
#                                      [shift_cpus] [shift_mem] [shift_time]
#                                      [qcat_cpus]  [qcat_mem]  [qcat_time]
#                                      [diff_cpus]  [diff_mem]  [diff_time]
#                                      [pval_cpus]  [pval_mem]  [pval_time]
#                                      [null_subsample] [blacklist_bed]
#                                      [floors_tsv] [min_signal]
#                                      [seed] [min_shift] [diff_sheet]
#
# Example (10 perms, full diff pipeline):
#   ./submit_perm_qcat_pvalue_slurm.sh samples_template.tsv 10 perm results \
#       4 16G 04:00:00 \
#       8 24G 04:00:00 \
#       4 24G 02:00:00 \
#       2 24G 02:00:00 \
#       1500000 detected_blacklist.bed \
#       floors.tsv 0.1 \
#       42 1000000 samples_template.tsv
#
# Quick test (no diff):
#   ./submit_perm_qcat_pvalue_slurm.sh samples_template.tsv 10 perm results \
#       4 16G 04:00:00 \
#       8 24G 04:00:00 \
#       4 24G 02:00:00 \
#       2 24G 02:00:00 \
#       1500000 detected_blacklist.bed \
#       floors.tsv 0.1

# ────────────────────────────────────────────────────────────
# Positional arguments
# ────────────────────────────────────────────────────────────
SHEET="${1:-samples_template.tsv}"
N_PERMS="${2:-10}"
PERM_PREFIX="${3:-perm}"
RESULTS_DIR="${4:-results}"

SHIFT_CPUS="${5:-4}"
SHIFT_MEM="${6:-16G}"
SHIFT_TIME="${7:-04:00:00}"

QCAT_CPUS="${8:-8}"
QCAT_MEM="${9:-24G}"
QCAT_TIME="${10:-04:00:00}"

DIFF_CPUS="${11:-4}"
DIFF_MEM="${12:-24G}"
DIFF_TIME="${13:-02:00:00}"

PVAL_CPUS="${14:-2}"
PVAL_MEM="${15:-24G}"
PVAL_TIME="${16:-02:00:00}"

NULL_SUBSAMPLE="${17:-1500000}"
BLACKLIST_BED="${18:-}"
FLOORS_TSV="${19:-}"
MIN_SIGNAL="${20:-}"
SEED="${21:-42}"
MIN_SHIFT="${22:-1000000}"
DIFF_SHEET="${23:-}"
CATEGORIES_YAML="${24:-}"     # optional: --categories YAML passed to bigwig_to_qcat AND compare_qcat
CHROMS_LIST="${25:-}"          # optional: space- or comma-separated list of chroms (passed through everywhere)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SHIFT_BW="$SCRIPT_DIR/shift_bigwig.py"
BIGWIG_TO_QCAT="$SCRIPT_DIR/bigwig_to_qcat.py"
COMPARE_QCAT="$SCRIPT_DIR/compare_qcat.py"
BEARING_PVALUE="$SCRIPT_DIR/bearing_pvalue.py"

for f in "$SHEET" "$SHIFT_BW" "$BIGWIG_TO_QCAT" "$BEARING_PVALUE"; do
  [[ -f "$f" ]] || { echo "ERROR: required file missing: $f" >&2; exit 1; }
done

if [[ -n "$DIFF_SHEET" ]]; then
  [[ -f "$COMPARE_QCAT" ]] || { echo "ERROR: --diff-sheet given but compare_qcat.py not found at: $COMPARE_QCAT" >&2; exit 1; }
  [[ -f "$DIFF_SHEET" ]] || { echo "ERROR: diff sheet not found: $DIFF_SHEET" >&2; exit 1; }
fi
[[ -z "$BLACKLIST_BED" || -f "$BLACKLIST_BED" ]] || { echo "ERROR: blacklist BED not found: $BLACKLIST_BED" >&2; exit 1; }
[[ -z "$FLOORS_TSV" || -f "$FLOORS_TSV" ]] || { echo "ERROR: floors TSV not found: $FLOORS_TSV" >&2; exit 1; }
[[ -z "$CATEGORIES_YAML" || -f "$CATEGORIES_YAML" ]] || { echo "ERROR: categories YAML not found: $CATEGORIES_YAML" >&2; exit 1; }

mkdir -p "$RESULTS_DIR" .slurm_logs

# Build manifest: sample, bw_csv, out, condition, replicate
MANIFEST=".slurm_perm_manifest.tsv"
awk -F '\t' '
  BEGIN { OFS="\t"; sample_i=0; bw_i=0; out_i=0; cond_i=0; rep_i=0 }
  /^[[:space:]]*#/ || NF==0 { next }
  sample_i==0 {
    for (i=1; i<=NF; i++) {
      h=$i; gsub(/^[[:space:]]+|[[:space:]]+$/, "", h); hl=tolower(h)
      if (hl=="sample" || hl=="name") sample_i=i
      else if (hl=="bw") bw_i=i
      else if (hl=="out") out_i=i
      else if (hl=="condition") cond_i=i
      else if (hl=="replicate") rep_i=i
    }
    if (sample_i==0 || bw_i==0 || out_i==0) {
      print "ERROR: sheet must contain columns: sample, bw, out" > "/dev/stderr"
      exit 2
    }
    next
  }
  {
    sample=$sample_i; bw=$bw_i; out=$out_i
    cond=(cond_i ? $cond_i : "")
    rep=(rep_i ? $rep_i : "")
    gsub(/^[[:space:]]+|[[:space:]]+$/, "", sample)
    gsub(/^[[:space:]]+|[[:space:]]+$/, "", bw)
    gsub(/^[[:space:]]+|[[:space:]]+$/, "", out)
    gsub(/^[[:space:]]+|[[:space:]]+$/, "", cond)
    gsub(/^[[:space:]]+|[[:space:]]+$/, "", rep)
    if (sample=="" || bw=="" || out=="") next
    print sample, bw, out, cond, rep
  }
' "$SHEET" > "$MANIFEST"

N_SAMPLES=$(wc -l < "$MANIFEST" | tr -d ' ')
[[ "$N_SAMPLES" -ge 1 ]] || { echo "ERROR: no valid sample rows in: $SHEET" >&2; exit 1; }

TOTAL_SHIFT_TASKS=$((N_SAMPLES * N_PERMS))
TOTAL_QCAT_TASKS=$((N_SAMPLES * N_PERMS))

echo "============================================================"
echo "Submission summary:"
echo "  Sheet:             $SHEET"
echo "  Samples:           $N_SAMPLES"
echo "  Permutations:      $N_PERMS"
echo "  Output prefix:     $PERM_PREFIX"
echo "  Results dir:       $RESULTS_DIR"
echo "  Shift tasks (S1):  $TOTAL_SHIFT_TASKS"
echo "  Qcat tasks (S2):   $TOTAL_QCAT_TASKS"
[[ -n "$DIFF_SHEET" ]] && echo "  Diff tasks (S3):   $N_PERMS  (one per perm)"
echo "  Pval tasks (S4):   $N_SAMPLES"
echo "  Seed:              $SEED"
echo "  Min shift:         $MIN_SHIFT bp"
[[ -n "$BLACKLIST_BED" ]] && echo "  Blacklist BED:     $BLACKLIST_BED"
[[ -n "$FLOORS_TSV" ]] && echo "  Floors TSV:        $FLOORS_TSV  (perms skip p5 pass)"
[[ -n "$MIN_SIGNAL" ]] && echo "  Min signal:        $MIN_SIGNAL"
[[ -n "$DIFF_SHEET" ]] && echo "  Diff sheet:        $DIFF_SHEET"
[[ -n "$CATEGORIES_YAML" ]] && echo "  Categories YAML:   $CATEGORIES_YAML"
[[ -n "$CHROMS_LIST" ]] && echo "  Chroms restricted: $CHROMS_LIST"
echo "============================================================"

# ────────────────────────────────────────────────────────────
# Stage 1: Shift
# ────────────────────────────────────────────────────────────
echo "Submitting Stage 1 (shift) as array 1-$TOTAL_SHIFT_TASKS"

SHIFT_JOBID=$(sbatch --parsable \
  --job-name=perm_shift \
  --cpus-per-task="$SHIFT_CPUS" \
  --mem="$SHIFT_MEM" \
  --time="$SHIFT_TIME" \
  --array="1-${TOTAL_SHIFT_TASKS}" \
  --output=".slurm_logs/perm_shift_%A_%a.out" \
  --export=ALL,MANIFEST="$MANIFEST",N_SAMPLES="$N_SAMPLES",N_PERMS="$N_PERMS",PERM_PREFIX="$PERM_PREFIX",SHIFT_BW="$SHIFT_BW",SEED="$SEED",MIN_SHIFT="$MIN_SHIFT" \
  --wrap='
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"

echo "PWD: $PWD"
echo "SLURM_SUBMIT_DIR: ${SLURM_SUBMIT_DIR:-(unset)}"

idx0=$((SLURM_ARRAY_TASK_ID - 1))
perm=$((idx0 / N_SAMPLES + 1))
sample_idx=$((idx0 % N_SAMPLES + 1))
perm_seed=$((SEED + perm))

line=$(sed -n "${sample_idx}p" "$MANIFEST")
sample=$(printf "%s\n" "$line" | cut -f1)
bw_csv=$(printf "%s\n" "$line" | cut -f2)

[[ -n "$sample" && -n "$bw_csv" ]] || { echo "ERROR: malformed manifest row $sample_idx" >&2; exit 2; }

out_dir="${PERM_PREFIX}${perm}/${sample}"
mkdir -p "$out_dir"

IFS="," read -r -a local_bws <<< "$bw_csv"
n_workers=${SLURM_CPUS_PER_TASK:-1}
n_bws=${#local_bws[@]}

echo "Stage 1: perm=${perm} sample=${sample} -- shifting ${n_bws} bw(s) with up to ${n_workers} workers"

job_list=$(mktemp)
trap "rm -f $job_list" EXIT
for bw in "${local_bws[@]}"; do
  fn=${bw##*/}
  stem=${fn%.bw}
  out_path="${out_dir}/${stem}_perm${perm}.bw"
  if [[ -f "$out_path" ]]; then
    echo "[skip exists] $out_path" >&2
    continue
  fi
  # Deterministic file-level seed. Uses cksum (POSIX, always available)
  # since sha256sum may not be in PATH on compute nodes and we want a
  # solution that does not depend on python being available.
  # Note: this seed differs from generate_perm_nulls.py'\''s
  # _deterministic_file_seed() (which uses SHA-256), but it'\''s still
  # deterministic and reproducible from the same inputs, which is what
  # matters for permutation runs.
  file_seed=$(printf "%s|%s|%s|%s" "$perm_seed" "$sample" "$bw" "$out_path" | cksum | cut -d" " -f1)
  if [[ -z "$file_seed" || ! "$file_seed" =~ ^[0-9]+$ ]]; then
    echo "ERROR: could not compute file_seed for $bw (got: \"$file_seed\")" >&2
    echo "       cksum availability: $(command -v cksum || echo NOT_FOUND)" >&2
    exit 3
  fi
  echo "$bw|$out_path|$file_seed"
done > "$job_list"

n_pending=$(wc -l < "$job_list" | tr -d " ")
if [[ "$n_pending" == "0" ]]; then
  echo "All shifted BigWigs already exist for perm=${perm} sample=${sample}"
  exit 0
fi

echo "Running ${n_pending} pending shift(s) with up to ${n_workers} parallel workers..."

# Use a simple for-loop with & for parallelism (avoids xargs quoting hazards).
# Throttle to n_workers concurrent processes using a semaphore via FIFO.
running=0
shift_pids=()
shift_failed=0
while IFS="|" read -r in_path out_path seed; do
  if [[ -z "$seed" ]]; then
    echo "ERROR: empty seed in job list for $in_path" >&2
    exit 4
  fi
  # Throttle
  while [[ $running -ge $n_workers ]]; do
    # Wait for ANY running job to finish
    if wait -n 2>/dev/null; then
      running=$((running - 1))
    else
      # wait -n returns non-zero if a job failed; track it but continue
      shift_failed=$((shift_failed + 1))
      running=$((running - 1))
    fi
  done
  # Launch shift
  echo "  shift: $in_path -> $out_path (seed=$seed)"
  python "$SHIFT_BW" --bw "$in_path" --out "$out_path" --seed "$seed" --min-shift "$MIN_SHIFT" &
  shift_pids+=($!)
  running=$((running + 1))
done < "$job_list"

# Wait for all remaining jobs
for pid in "${shift_pids[@]}"; do
  if ! wait "$pid"; then
    shift_failed=$((shift_failed + 1))
  fi
done

if [[ $shift_failed -gt 0 ]]; then
  echo "ERROR: $shift_failed shift task(s) failed for perm=${perm} sample=${sample}" >&2
  exit 5
fi

echo "Stage 1 complete: perm=${perm} sample=${sample}"
')

echo "  SHIFT job id:  $SHIFT_JOBID"

# ────────────────────────────────────────────────────────────
# Stage 2: Score
# ────────────────────────────────────────────────────────────
echo "Submitting Stage 2 (qcat) as array 1-$TOTAL_QCAT_TASKS  (afterok:$SHIFT_JOBID)"

QCAT_JOBID=$(sbatch --parsable \
  --job-name=perm_qcat \
  --cpus-per-task="$QCAT_CPUS" \
  --mem="$QCAT_MEM" \
  --time="$QCAT_TIME" \
  --array="1-${TOTAL_QCAT_TASKS}" \
  --dependency="afterok:${SHIFT_JOBID}" \
  --output=".slurm_logs/perm_qcat_%A_%a.out" \
  --export=ALL,MANIFEST="$MANIFEST",N_SAMPLES="$N_SAMPLES",N_PERMS="$N_PERMS",PERM_PREFIX="$PERM_PREFIX",BIGWIG_TO_QCAT="$BIGWIG_TO_QCAT",BLACKLIST_BED="$BLACKLIST_BED",FLOORS_TSV="$FLOORS_TSV",MIN_SIGNAL="$MIN_SIGNAL",CATEGORIES_YAML="$CATEGORIES_YAML",CHROMS_LIST="$CHROMS_LIST" \
  --wrap='
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"

idx0=$((SLURM_ARRAY_TASK_ID - 1))
perm=$((idx0 / N_SAMPLES + 1))
sample_idx=$((idx0 % N_SAMPLES + 1))

line=$(sed -n "${sample_idx}p" "$MANIFEST")
sample=$(printf "%s\n" "$line" | cut -f1)
bw_csv=$(printf "%s\n" "$line" | cut -f2)
out=$(printf "%s\n" "$line" | cut -f3)

[[ -n "$sample" && -n "$bw_csv" && -n "$out" ]] || { echo "ERROR: malformed manifest row $sample_idx" >&2; exit 2; }

out_name=${out##*/}
base=${out_name%.qcat.bgz}

IFS="," read -r -a local_bws <<< "$bw_csv"
shifted=()
for bw in "${local_bws[@]}"; do
  fn=${bw##*/}
  stem=${fn%.bw}
  shifted+=("${PERM_PREFIX}${perm}/${sample}/${stem}_perm${perm}.bw")
done

null_qcat="${PERM_PREFIX}${perm}/${sample}/${base}_perm${perm}.qcat.bgz"
mkdir -p "${PERM_PREFIX}${perm}/${sample}"

if [[ -f "$null_qcat" ]]; then
  echo "[skip exists] $null_qcat"
  exit 0
fi

qcat_cmd=(
  python "$BIGWIG_TO_QCAT"
  --bw "${shifted[@]}"
  --out "$null_qcat"
  --jobs "${SLURM_CPUS_PER_TASK:-1}"
  --no-extras
)

[[ -n "${BLACKLIST_BED:-}" ]] && qcat_cmd+=(--blacklist "$BLACKLIST_BED")
[[ -n "${FLOORS_TSV:-}" ]] && qcat_cmd+=(--floors-tsv "$FLOORS_TSV" --sample-name "$sample")
[[ -n "${MIN_SIGNAL:-}" ]] && qcat_cmd+=(--min-signal "$MIN_SIGNAL")
[[ -n "${CATEGORIES_YAML:-}" ]] && qcat_cmd+=(--categories "$CATEGORIES_YAML")
if [[ -n "${CHROMS_LIST:-}" ]]; then
  IFS=", " read -r -a chrom_arr <<< "$CHROMS_LIST"
  qcat_cmd+=(--chroms "${chrom_arr[@]}")
fi

echo "Running: ${qcat_cmd[*]}"
"${qcat_cmd[@]}"
echo "Stage 2 complete: perm=${perm} sample=${sample}"
')

echo "  QCAT job id:   $QCAT_JOBID"

# ────────────────────────────────────────────────────────────
# Stage 3 (optional): Diff
# ────────────────────────────────────────────────────────────
DIFF_JOBID=""
if [[ -n "$DIFF_SHEET" ]]; then
  echo "Submitting Stage 3 (diff) as array 1-$N_PERMS  (afterok:$QCAT_JOBID)"

  DIFF_JOBID=$(sbatch --parsable \
    --job-name=perm_diff \
    --cpus-per-task="$DIFF_CPUS" \
    --mem="$DIFF_MEM" \
    --time="$DIFF_TIME" \
    --array="1-${N_PERMS}" \
    --dependency="afterok:${QCAT_JOBID}" \
    --output=".slurm_logs/perm_diff_%A_%a.out" \
    --export=ALL,MANIFEST="$MANIFEST",N_SAMPLES="$N_SAMPLES",N_PERMS="$N_PERMS",PERM_PREFIX="$PERM_PREFIX",COMPARE_QCAT="$COMPARE_QCAT",CATEGORIES_YAML="$CATEGORIES_YAML",CHROMS_LIST="$CHROMS_LIST" \
    --wrap='
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"

perm=$SLURM_ARRAY_TASK_ID
diff_dir="${PERM_PREFIX}${perm}/diff_comparison"
mkdir -p "$diff_dir"

null_sheet="${diff_dir}/_null_sheet_perm${perm}.tsv"
printf "sample\tcondition\treplicate\tqcat\n" > "$null_sheet"

while IFS=$'\''\t'\'' read -r sample bw_csv out cond rep; do
  out_name=${out##*/}
  base=${out_name%.qcat.bgz}
  null_qcat="${PERM_PREFIX}${perm}/${sample}/${base}_perm${perm}.qcat.bgz"
  if [[ ! -f "$null_qcat" ]]; then
    echo "ERROR: missing null qcat for ${sample} perm${perm}: $null_qcat" >&2
    exit 2
  fi
  printf "%s\t%s\t%s\t%s\n" "$sample" "$cond" "$rep" "$null_qcat" >> "$null_sheet"
done < "$MANIFEST"

diff_cmd=(
  python "$COMPARE_QCAT"
  --sheet "$null_sheet"
  --out "$diff_dir"
  --diff-only
  --no-clip
  --skip-pca
  --skip-q-pair-jsd
  --workers "${SLURM_CPUS_PER_TASK:-1}"
)

[[ -n "${CATEGORIES_YAML:-}" ]] && diff_cmd+=(--categories "$CATEGORIES_YAML")
if [[ -n "${CHROMS_LIST:-}" ]]; then
  IFS=", " read -r -a chrom_arr <<< "$CHROMS_LIST"
  diff_cmd+=(--chroms "${chrom_arr[@]}")
fi

echo "Running: ${diff_cmd[*]}"
"${diff_cmd[@]}"
echo "Stage 3 complete: perm=${perm}"
')

  echo "  DIFF job id:   $DIFF_JOBID"
fi

# ────────────────────────────────────────────────────────────
# Stage 4: P-value
# ────────────────────────────────────────────────────────────
PVAL_DEPENDENCY="afterok:${QCAT_JOBID}"
[[ -n "$DIFF_JOBID" ]] && PVAL_DEPENDENCY="${PVAL_DEPENDENCY}:${DIFF_JOBID}"

echo "Submitting Stage 4 (pvalue) as array 1-$N_SAMPLES  (afterok:$PVAL_DEPENDENCY)"

PVAL_JOBID=$(sbatch --parsable \
  --job-name=perm_pval \
  --cpus-per-task="$PVAL_CPUS" \
  --mem="$PVAL_MEM" \
  --time="$PVAL_TIME" \
  --array="1-${N_SAMPLES}" \
  --dependency="${PVAL_DEPENDENCY}" \
  --output=".slurm_logs/perm_pval_%A_%a.out" \
  --export=ALL,MANIFEST="$MANIFEST",N_PERMS="$N_PERMS",PERM_PREFIX="$PERM_PREFIX",RESULTS_DIR="$RESULTS_DIR",BEARING_PVALUE="$BEARING_PVALUE",NULL_SUBSAMPLE="$NULL_SUBSAMPLE" \
  --wrap='
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"

line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MANIFEST")
sample=$(printf "%s\n" "$line" | cut -f1)
bw_csv=$(printf "%s\n" "$line" | cut -f2)
out=$(printf "%s\n" "$line" | cut -f3)

[[ -n "$sample" && -n "$bw_csv" && -n "$out" ]] || { echo "ERROR: malformed manifest row $SLURM_ARRAY_TASK_ID" >&2; exit 2; }

out_name=${out##*/}
base=${out_name%.qcat.bgz}

nulls=()
for ((p=1; p<=N_PERMS; p++)); do
  null_qcat="${PERM_PREFIX}${p}/${sample}/${base}_perm${p}.qcat.bgz"
  if [[ ! -f "$null_qcat" ]]; then
    dir="${PERM_PREFIX}${p}/${sample}"
    mapfile -t found < <(find "$dir" -maxdepth 1 -name "*.qcat.bgz" | sort)
    if [[ ${#found[@]} -eq 1 ]]; then
      null_qcat="${found[0]}"
    else
      echo "ERROR: could not resolve null qcat for sample=${sample} perm=${p}" >&2
      exit 2
    fi
  fi
  nulls+=("$null_qcat")
done

mkdir -p "$RESULTS_DIR"
pvalue_cmd=(
  python "$BEARING_PVALUE"
  --qcat "$out"
  --null-qcat "${nulls[@]}"
  --out-prefix "$RESULTS_DIR/$sample"
  --fdr 0.05
  --score-plot
)

if [[ -n "${NULL_SUBSAMPLE:-}" && "${NULL_SUBSAMPLE}" != "0" ]]; then
  pvalue_cmd+=(--null-subsample "$NULL_SUBSAMPLE")
fi

echo "Running: ${pvalue_cmd[*]}"
"${pvalue_cmd[@]}"
echo "Stage 4 complete: sample=${sample}"
')

echo "  PVAL job id:   $PVAL_JOBID"

# ────────────────────────────────────────────────────────────
# Stage 5 (optional): Diff p-value calculation per comparison.
# Only runs if --diff-sheet was provided (i.e., Stage 3 ran).
# Depends on Stage 3 + Stage 4. Replaces run_perm_diff_pvalue.sh.
# ────────────────────────────────────────────────────────────
DIFF_PVAL_JOBID=""
COMPARISON_DIR=""
if [[ -n "$DIFF_SHEET" ]]; then
  # Compute the number of pairwise comparisons from unique conditions in the
  # diff sheet. We use the same MANIFEST column 4 (condition) we built earlier.
  N_CONDITIONS=$(awk -F'\t' '$4 != ""' "$MANIFEST" | cut -f4 | sort -u | wc -l | tr -d ' ')
  N_DIFF_COMPS=$(( N_CONDITIONS * (N_CONDITIONS - 1) / 2 ))

  if [[ "$N_DIFF_COMPS" -lt 1 ]]; then
    echo "WARNING: cannot determine diff comparisons (need >=2 conditions in MANIFEST col 4)."
    echo "         Skipping Stage 5. Run run_perm_diff_pvalue.sh manually if needed."
  else
    # The real-data diff qcats need to live somewhere. By convention, they're
    # produced by compare_qcat.py on the REAL sample sheet into a directory.
    # We'll assume the user has produced these separately (e.g., into
    # comparison_vdj/ or similar) OR we let them pass via the 24th positional
    # argument. For now we look for $COMPARISON_DIR env var, else derive from
    # PERM_PREFIX naming convention.
    if [[ -n "${COMPARISON_DIR_OVERRIDE:-}" ]]; then
      COMPARISON_DIR="$COMPARISON_DIR_OVERRIDE"
    else
      # Default: assume real-data diff qcats are in ./comparison_main
      COMPARISON_DIR="comparison_main"
    fi

    if [[ ! -d "$COMPARISON_DIR" ]]; then
      echo "WARNING: real-data diff qcats directory not found: $COMPARISON_DIR"
      echo "         Stage 5 needs real-data diff qcats produced by"
      echo "         compare_qcat.py on the REAL sample sheet."
      echo "         Generate them with:"
      echo "           python compare_qcat.py --sheet $SHEET --out $COMPARISON_DIR --diff-only ..."
      echo "         Then re-submit Stage 5 separately via run_perm_diff_pvalue.sh"
      echo "         or set COMPARISON_DIR_OVERRIDE and re-run this script."
      echo "         Skipping Stage 5."
    else
      echo "Submitting Stage 5 (diff pvalue) as array 1-$N_DIFF_COMPS  (afterok:$QCAT_JOBID${DIFF_JOBID:+:$DIFF_JOBID})"

      DIFF_PVAL_DEPENDENCY="afterok:${QCAT_JOBID}"
      [[ -n "$DIFF_JOBID" ]] && DIFF_PVAL_DEPENDENCY="${DIFF_PVAL_DEPENDENCY}:${DIFF_JOBID}"

      DIFF_PVAL_JOBID=$(sbatch --parsable \
        --job-name=perm_diff_pval \
        --cpus-per-task="$PVAL_CPUS" \
        --mem="$PVAL_MEM" \
        --time="$PVAL_TIME" \
        --array="1-${N_DIFF_COMPS}" \
        --dependency="${DIFF_PVAL_DEPENDENCY}" \
        --output=".slurm_logs/perm_diff_pval_%A_%a.out" \
        --export=ALL,COMPARISON_DIR="$COMPARISON_DIR",N_PERMS="$N_PERMS",PERM_PREFIX="$PERM_PREFIX",RESULTS_DIR="$RESULTS_DIR",BEARING_PVALUE="$BEARING_PVALUE",NULL_SUBSAMPLE="$NULL_SUBSAMPLE" \
        --wrap='
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"

# Enumerate the diff qcat files at runtime (after Stages 2-3 finish, so the
# array index can map to a specific file).
mapfile -t diff_files < <(find "$COMPARISON_DIR" -maxdepth 1 -name "diff_*.qcat.bgz" | sort)

if [[ ${#diff_files[@]} -eq 0 ]]; then
  echo "ERROR: no diff_*.qcat.bgz files found in $COMPARISON_DIR" >&2
  exit 2
fi

if [[ "$SLURM_ARRAY_TASK_ID" -gt "${#diff_files[@]}" ]]; then
  echo "Skipping: array task ${SLURM_ARRAY_TASK_ID} exceeds ${#diff_files[@]} diff files"
  exit 0
fi

idx=$((SLURM_ARRAY_TASK_ID - 1))
main_qcat="${diff_files[$idx]}"
comp=$(basename "$main_qcat")
comp="${comp#diff_}"
comp="${comp%.qcat.bgz}"

echo "Stage 5: comparison=${comp}"

nulls=()
for ((p=1; p<=N_PERMS; p++)); do
  null_qcat="${PERM_PREFIX}${p}/diff_comparison/diff_${comp}.qcat.bgz"
  if [[ -f "$null_qcat" ]]; then
    nulls+=("$null_qcat")
  else
    echo "WARNING: null qcat not found for perm$p: $null_qcat" >&2
  fi
done

if [[ ${#nulls[@]} -eq 0 ]]; then
  echo "ERROR: no null qcats found for $comp" >&2
  exit 2
fi

echo "Found ${#nulls[@]} null qcats for $comp"

mkdir -p "$RESULTS_DIR"
out_prefix="$RESULTS_DIR/diff_${comp}"
pvalue_cmd=(
  python "$BEARING_PVALUE"
  --qcat "$main_qcat"
  --null-qcat "${nulls[@]}"
  --diff
  --out-prefix "$out_prefix"
  --fdr 0.05
  --score-plot
)

if [[ -n "${NULL_SUBSAMPLE:-}" && "${NULL_SUBSAMPLE}" != "0" ]]; then
  pvalue_cmd+=(--null-subsample "$NULL_SUBSAMPLE")
fi

echo "Running: ${pvalue_cmd[*]}"
"${pvalue_cmd[@]}"
echo "Stage 5 complete: comparison=${comp}"
')

      echo "  DIFF_PVAL job id: $DIFF_PVAL_JOBID"
    fi
  fi
fi

echo
echo "============================================================"
echo "Submitted multi-stage pipeline:"
echo "  Stage 1 (shift):     job $SHIFT_JOBID  [$TOTAL_SHIFT_TASKS array tasks]"
echo "  Stage 2 (qcat):      job $QCAT_JOBID  [$TOTAL_QCAT_TASKS array tasks, depends on $SHIFT_JOBID]"
[[ -n "$DIFF_JOBID" ]] && echo "  Stage 3 (diff):      job $DIFF_JOBID  [$N_PERMS array tasks, depends on $QCAT_JOBID]"
echo "  Stage 4 (pval):      job $PVAL_JOBID  [$N_SAMPLES array tasks, depends on Stage 2/3]"
[[ -n "$DIFF_PVAL_JOBID" ]] && echo "  Stage 5 (diff_pval): job $DIFF_PVAL_JOBID  [$N_DIFF_COMPS array tasks, depends on Stage 2/3]"
echo
JOBS_LIST="$SHIFT_JOBID,$QCAT_JOBID"
[[ -n "$DIFF_JOBID" ]] && JOBS_LIST="$JOBS_LIST,$DIFF_JOBID"
JOBS_LIST="$JOBS_LIST,$PVAL_JOBID"
[[ -n "$DIFF_PVAL_JOBID" ]] && JOBS_LIST="$JOBS_LIST,$DIFF_PVAL_JOBID"
echo "Monitor:"
echo "  squeue -j $JOBS_LIST"
echo
echo "Cancel everything:"
echo "  scancel $JOBS_LIST"
echo "============================================================"
