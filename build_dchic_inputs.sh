#!/usr/bin/env bash
# build_dchic_inputs.sh -- reproducible replacement for the manual mk() helper.
#
# Builds the dcHiC input subset files in each dchic_in_<RES>/ directory from the
# per-sample .matrix/.bed produced by run_dchic_convert.sh:
#   input_manuscript.txt  DN,DP,EbKO,ProB,3T3  (n=2 each; the 5-condition panel)
#   input_lymphoid.txt     DN,DP,EbKO,ProB      (no 3T3; the clean lymphoid test)
#
# VERIFY before first use: the sample ids below must match the .matrix/.bed
# basenames in dchic_in_<RES>/ EXACTLY (casing matters: s3T3 lowercase, EBKO
# upper, as in the KR h5 filenames). Adjust the lists if your files differ.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

# Parent dir for the dchic_in_<res> scratch (must match run_dchic_convert.sh).
# Defaults to the repo root for a standalone run; the Snakemake rule points it
# at the results outdir.
WORK="${DCHIC_WORKDIR:-$PWD}"

RESOLUTIONS="${RESOLUTIONS:-100kb 250kb}"
MANUSCRIPT_SAMPLES="DN_rep1 DN_rep2 DP_rep1 DP_rep2 EBKO_rep1 EBKO_rep2 ProB_rep1 ProB_rep2 s3T3_rep1 s3T3_rep2"
LYMPHOID_SAMPLES="DN_rep1 DN_rep2 DP_rep1 DP_rep2 EBKO_rep1 EBKO_rep2 ProB_rep1 ProB_rep2"

# columns: matrix_path  bed_path  sample_label  experiment
build_input() {
  local dir="$1"; shift
  local reslabel="${dir##*dchic_in_}"   # basename suffix, works for abs or rel dir
  for s in "$@"; do
    local exp="${s%_rep*}"
    printf "%s\t%s\t%s_%s\t%s\n" \
      "$dir/$s.matrix" "$dir/$s.bed" "$s" "$reslabel" "$exp"
  done
}

# Fail loud if any listed sample has no converted .matrix/.bed in $dir. Without
# this, a casing mismatch (e.g. EBKO vs EbKO) silently drops a condition and
# dcHiC runs on a smaller-but-valid-looking panel -- exactly the kind of error
# that produces a wrong table with no visible failure.
check_samples_present() {
  local dir="$1"; shift
  local missing=()
  for s in "$@"; do
    [ -s "$dir/$s.matrix" ] && [ -s "$dir/$s.bed" ] || missing+=("$s")
  done
  if [ ${#missing[@]} -gt 0 ]; then
    echo "FATAL: $dir is missing .matrix/.bed for: ${missing[*]}" >&2
    echo "       The sample ids in build_dchic_inputs.sh must match the converted" >&2
    echo "       basenames EXACTLY (casing: s3T3 lower, EBKO upper). Fix the lists" >&2
    echo "       or re-run run_dchic_convert.sh, then retry." >&2
    exit 1
  fi
}

for res in $RESOLUTIONS; do
  dir="$WORK/dchic_in_${res}"
  if [ ! -d "$dir" ]; then
    echo "SKIP: $dir not present (run run_dchic_convert.sh first)"; continue
  fi
  check_samples_present "$dir" $MANUSCRIPT_SAMPLES
  check_samples_present "$dir" $LYMPHOID_SAMPLES
  build_input "$dir" $MANUSCRIPT_SAMPLES > "$dir/input_manuscript.txt"
  build_input "$dir" $LYMPHOID_SAMPLES   > "$dir/input_lymphoid.txt"
  echo "wrote $dir/input_manuscript.txt ($(wc -l < "$dir/input_manuscript.txt") rows), input_lymphoid.txt ($(wc -l < "$dir/input_lymphoid.txt") rows)"
done
