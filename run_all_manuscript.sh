#!/usr/bin/env bash
# Manuscript compartment-STABILITY claim, tested (n=2 per condition, matches the
# merged matrices). Replaces the sign-flip count with dcHiC's tested calls.
# Positive control: S3T3-vs-any-lymphoid MUST be strongly significant. If it is
# not, the run is wrong -- stop and check before trusting lymphoid-vs-lymphoid.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

export INPUT_PREFIX="input_manuscript"
export DIFF_PREFIX="manuscript"
export RESOLUTIONS="100kb 250kb"   # manuscript compartment ladder is 100/250/500 kb; 250 as the confirmation res
source ./_dchic_run.sh

dchic_run_all

echo ""
echo "Manuscript run complete. The differential bedGraph per resolution is:"
echo "  dchic_in_<res>/DifferentialResult/manuscript_<res>/fdr_result/differential.intra_sample_combined.pcQnm.bedGraph"
echo ""
echo "For the Tcrb-locus compartment call across conditions, extract the Tcrb TAD:"
echo "  python3 dchic_region_result.py --bedgraph <that file> \\"
echo "      --region chr6:40400000-42400000 --label Tcrb --exp-a DN --exp-b S3T3 --out ..."
echo "Repeat --exp-b for DP, EBKO, ProB to tabulate the stability claim."
