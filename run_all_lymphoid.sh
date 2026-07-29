#!/usr/bin/env bash
# run_all_lymphoid.sh -- lymphoid-only dcHiC compartment run (no 3T3), the clean
# test of the manuscript's DP/EbKO stability + Pro-B RC shift. Mirrors
# run_all_manuscript.sh. Writes its region-extract TSVs to a separate directory
# (paper/table_sources/lymphoid) so they do not collide with the 5-condition run
# (dchic_region names files by LABEL+EXPB, which are identical across runs).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

export INPUT_PREFIX="input_lymphoid"
export DIFF_PREFIX="lymphoid"
export RESOLUTIONS="${RESOLUTIONS:-100kb}"          # 100 kb is the primary lymphoid test
export TABLE_OUT="${TABLE_OUT:-$(pwd)/paper/table_sources/lymphoid}"
mkdir -p "$TABLE_OUT"
source ./_dchic_run.sh

dchic_run_all

# region extracts: DN vs each lymphoid condition over the Tcrb window
dchic_region chr6:40400000-42400000 Tcrb DN DP
dchic_region chr6:40400000-42400000 Tcrb DN EBKO
dchic_region chr6:40400000-42400000 Tcrb DN ProB

echo ""
echo "Lymphoid run complete. TSVs in $TABLE_OUT :"
echo "  Tcrb_compartment_dchic_{DP,EBKO,ProB}_<res>.tsv"
