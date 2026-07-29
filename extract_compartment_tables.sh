#!/usr/bin/env bash
# extract_compartment_tables.sh -- (re)extract the Tcrb compartment TSVs from
# EXISTING dcHiC bedGraphs, without re-running dcHiC. Fast, idempotent. Produces
# the exact inputs build_compartment_supp_table.py reads:
#   paper/table_sources/Tcrb_compartment_dchic_{S3T3,DP,EBKO,ProB}_<res>.tsv  (5-cond)
#   paper/table_sources/lymphoid/Tcrb_compartment_dchic_{DP,EBKO,ProB}_<res>.tsv
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
RES="${RES:-100kb}"

# 5-condition (manuscript) extracts
export INPUT_PREFIX="input_manuscript" DIFF_PREFIX="manuscript" RESOLUTIONS="$RES"
export TABLE_OUT="$(pwd)/paper/table_sources"; mkdir -p "$TABLE_OUT"
source ./_dchic_run.sh
dchic_region chr6:40400000-42400000 Tcrb DN S3T3
dchic_region chr6:40400000-42400000 Tcrb DN DP
dchic_region chr6:40400000-42400000 Tcrb DN EBKO
dchic_region chr6:40400000-42400000 Tcrb DN ProB

# lymphoid-only extracts (separate directory)
export INPUT_PREFIX="input_lymphoid" DIFF_PREFIX="lymphoid" RESOLUTIONS="$RES"
export TABLE_OUT="$(pwd)/paper/table_sources/lymphoid"; mkdir -p "$TABLE_OUT"
dchic_region chr6:40400000-42400000 Tcrb DN DP
dchic_region chr6:40400000-42400000 Tcrb DN EBKO
dchic_region chr6:40400000-42400000 Tcrb DN ProB
echo "compartment TSVs written under paper/table_sources/ and paper/table_sources/lymphoid/"
