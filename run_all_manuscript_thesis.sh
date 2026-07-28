#!/usr/bin/env bash
# Manuscript compartment claim, THESIS version (DN and DP gain rep3/rep4).
# Sensitivity check on run_all_manuscript.sh. Unbalanced (DN/DP n=4, rest n=2).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

export INPUT_PREFIX="input_thesis"
export DIFF_PREFIX="manuscript_thesis"
export RESOLUTIONS="100kb 250kb"   # manuscript compartment ladder is 100/250/500 kb; 250 as the confirmation res
source ./_dchic_run.sh

dchic_run_all

echo ""
echo "Manuscript THESIS run complete (DN/DP n=4). Compare significant-bin counts"
echo "and the Tcrb-region calls against run_all_manuscript.sh to show the"
echo "stability claim does not depend on replicate count."
