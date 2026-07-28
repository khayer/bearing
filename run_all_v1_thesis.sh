#!/usr/bin/env bash
# V1 compartment question, THESIS version (DN gains rep3/rep4; n=4 for DN).
# This is the SENSITIVITY check on run_all_v1.sh: if the A->B switch at V1 holds
# at n=4, it is not an n=2 variance artifact. Unbalanced design (DN n=4, mutants
# n=2) -- state that in the writeup.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

export INPUT_PREFIX="input_v1_thesis"
export DIFF_PREFIX="v1_thesis"
export RESOLUTIONS="100kb 50kb"
source ./_dchic_run.sh

dchic_run_all
dchic_region "chr6:40872100-40908238" "V1" "DN" "dV1P"
dchic_region "chr6:40872100-40908238" "V1" "DN" "dV1CTCF"

echo ""
echo "V1 THESIS run complete. Primary = run_all_v1.sh (n=2, matches manuscript);"
echo "this n=4 run is the power/sensitivity check. If V1 still flips A->B here,"
echo "the result is robust to replicate count."
