#!/usr/bin/env bash
# V1 compartment question, manuscript-matched (n=2 per condition; DN is control).
# Runs cis/select/analyze/viz at 100kb and 50kb, then extracts the V1 region for
# both dV1P (promoter deletion) and dV1CTCF (CTCF-site deletion).
#
# Run from the directory that CONTAINS dchic_in_100kb/ and dchic_in_50kb/,
# with the dchic conda env active:
#     conda activate dchic
#     bash run_all_v1.sh 2>&1 | tee run_all_v1.out
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

export INPUT_PREFIX="input_v1"
export DIFF_PREFIX="v1"
export RESOLUTIONS="100kb 50kb"
source ./_dchic_run.sh

dchic_run_all

# V1 region: chr6:40,872,100-40,908,238. DN is the control (config: control: DN).
dchic_region "chr6:40872100-40908238" "V1" "DN" "dV1P"
dchic_region "chr6:40872100-40908238" "V1" "DN" "dV1CTCF"

echo ""
echo "V1 run complete. Compare padj and sign-change across 100kb and 50kb:"
echo "  paper/table_sources/V1_compartment_dchic_dV1P_100kb.tsv"
echo "  paper/table_sources/V1_compartment_dchic_dV1P_50kb.tsv"
