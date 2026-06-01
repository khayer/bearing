#!/usr/bin/env bash
#
# run_replicate_calibration.sh
# ============================
# Within-condition replicate-differential FDR calibration across every
# condition in a BEARING sample sheet, REUSING pre-computed per-sample
# permutation qcats (no re-scoring).
#
# Given:
#   * real per-sample qcats         (from the master sheet's qcat/out column)
#   * permuted per-sample qcats     (already scored; located by a path
#                                    template with {sample} and {perm})
# this script, per condition with exactly two replicates r1,r2:
#   1. diffs the real pair            (compare_qcat.py --diff-only) -> observed
#   2. diffs each permuted pair       (compare_qcat.py --diff-only) -> nulls
#   3. runs bearing_calibration.py    -> p-value uniformity diagnostics
# then aggregates lambda / floor% / BH-significant / verdict for all
# conditions into one table.
#
# Building a rep1-vs-rep2 diff only needs aligning two existing qcats and
# subtracting, which is exactly compare_qcat.py's job -- so no shift_bigwig.py,
# no bigwig_to_qcat.py, no generate_perm_nulls.py here. If you do NOT already
# have permuted per-sample qcats, generate them first (generate_perm_nulls.py
# or shift_bigwig.py + bigwig_to_qcat.py), then point --perm-template at them.
#
# Diff filenames from compare_qcat.py are token-sanitized, so this script
# GLOBS diff_*_vs_*.qcat.bgz rather than guessing the name. Each diff run uses
# a 2-row sheet, so the glob is unambiguous (it errors on zero or >1 matches).
#
# USAGE
#   ./run_replicate_calibration.sh \
#       --sheet samples.tsv \
#       --out   calib_run/ \
#       --perm-template 'perm_v6_paper{perm}/{sample}/{sample}_perm{perm}.qcat.bgz' \
#       --n-perms 3 \
#       --scripts-dir ~/data/tools/scripts/bigwig_to_qcat/ \
#       [--obs-qcat-col out] [--fdr 0.05] [--min-signal 0.0] \
#       [--conditions "DN EbKO"] [--force] [--dry-run]
#
# The real per-sample qcat for each replicate is read from the master sheet
# column named by --obs-qcat-col (default: out). The permuted per-sample qcat
# is --perm-template with {sample} and {perm} substituted.
#
# CLUSTER NOTE
#   This is fast (diff + p-value only). Conditions run serially; to fan out on
#   SLURM, submit one job per condition with --conditions "<COND>".
#
# ENV
#   PYTHON   python interpreter to use (default: python)

set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
SHEET=""
OUT=""
PERM_TEMPLATE=""
N_PERMS=3
OBS_QCAT_COL="out"
FDR=0.05
MIN_SIGNAL=0.0
CONDITIONS=""
SCRIPTS_DIR=""
FORCE=0
DRY_RUN=0

PY="${PYTHON:-python}"

usage() { sed -n '2,52p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }
expand_home() { case "$1" in "~"/*) printf '%s' "${HOME}/${1#\~/}";; *) printf '%s' "$1";; esac; }

# ---------------------------------------------------------------------------
# Parse args
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sheet)         SHEET="$2"; shift 2 ;;
        --out)           OUT="$2"; shift 2 ;;
        --perm-template) PERM_TEMPLATE="$2"; shift 2 ;;
        --n-perms)       N_PERMS="$2"; shift 2 ;;
        --obs-qcat-col)  OBS_QCAT_COL="$(printf '%s' "$2" | tr '[:upper:]' '[:lower:]')"; shift 2 ;;
        --fdr)           FDR="$2"; shift 2 ;;
        --min-signal)    MIN_SIGNAL="$2"; shift 2 ;;
        --conditions)    CONDITIONS="$2"; shift 2 ;;
        --scripts-dir)   SCRIPTS_DIR="$(expand_home "$2")"; shift 2 ;;
        --force)         FORCE=1; shift ;;
        --dry-run)       DRY_RUN=1; shift ;;
        -h|--help)       usage 0 ;;
        *)               echo "ERROR: unknown argument: $1" >&2; usage 2 ;;
    esac
done

[[ -z "$SHEET"         ]] && { echo "ERROR: --sheet is required" >&2; usage 2; }
[[ -z "$OUT"           ]] && { echo "ERROR: --out is required" >&2; usage 2; }
[[ -z "$PERM_TEMPLATE" ]] && { echo "ERROR: --perm-template is required" >&2; usage 2; }
[[ -f "$SHEET" ]] || { echo "ERROR: sheet not found: $SHEET" >&2; exit 1; }
case "$PERM_TEMPLATE" in
    *'{sample}'*'{perm}'*|*'{perm}'*'{sample}'*) : ;;
    *) echo "ERROR: --perm-template must contain both {sample} and {perm}" >&2; exit 2 ;;
esac

if [[ -z "$SCRIPTS_DIR" ]]; then
    SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
COMPARE_QCAT="$SCRIPTS_DIR/compare_qcat.py"
BEARING_PVALUE="$SCRIPTS_DIR/bearing_pvalue.py"
BEARING_CALIB="$SCRIPTS_DIR/bearing_calibration.py"
for s in "$COMPARE_QCAT" "$BEARING_PVALUE" "$BEARING_CALIB"; do
    [[ -f "$s" ]] || { echo "ERROR: missing helper script: $s" >&2; exit 1; }
done

mkdir -p "$OUT"
SHEETS_DIR="$OUT/sheets"; mkdir -p "$SHEETS_DIR"
AGG="$OUT/calibration_summary.tsv"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
run() { echo "+ $*"; if [[ $DRY_RUN -eq 0 ]]; then "$@"; fi; }

expand_template() {  # template sample perm
    local t="$1"; t="${t//\{sample\}/$2}"; t="${t//\{perm\}/$3}"; printf '%s' "$t"
}

extract() {  # label file sed_expr -> value or "NA"; never fails
    local v; v=$(grep -m1 -- "$1" "$2" 2>/dev/null | sed -E "$3" 2>/dev/null) || true
    [[ -z "$v" ]] && v="NA"; printf '%s' "$v"
}

write_qcat_sheet() {  # outfile s1 q1 s2 q2  (condition := sample so the diff is rep-vs-rep)
    { printf 'sample\tcondition\treplicate\tqcat\n'
      printf '%s\t%s\t1\t%s\n' "$2" "$2" "$3"
      printf '%s\t%s\t1\t%s\n' "$4" "$4" "$5"
    } > "$1"
}

locate_one_diff() {  # dir -> path of the single diff_*_vs_*.qcat.bgz, or empty
    local d="$1" g
    shopt -s nullglob; g=( "$d"/diff_*_vs_*.qcat.bgz ); shopt -u nullglob
    if [[ ${#g[@]} -eq 1 ]]; then printf '%s' "${g[0]}"; fi
    return 0
}

# Distinct conditions in sheet order (real condition column).
mapfile -t ALL_CONDS < <(
    awk -F'\t' '
        NR==1{ for(i=1;i<=NF;i++) if(tolower($i)=="condition") ci=i
               if(!ci){print "ERROR: no condition column">"/dev/stderr"; exit 2} next }
        $0!~/[^[:space:]]/ || $0 ~ /^[[:space:]]*#/ { next }
        { c=$ci; if(!(c in seen)){seen[c]=1; print c} }
    ' "$SHEET"
)

declare -a CONDS=()
if [[ -n "$CONDITIONS" ]]; then
    for c in $CONDITIONS; do
        if printf '%s\n' "${ALL_CONDS[@]}" | grep -qx -- "$c"; then CONDS+=("$c")
        else echo "WARNING: condition '$c' not in sheet; skipping" >&2; fi
    done
else
    CONDS=("${ALL_CONDS[@]}")
fi
[[ ${#CONDS[@]} -gt 0 ]] || { echo "ERROR: no conditions to process" >&2; exit 1; }

echo "============================================================"
echo "  run_replicate_calibration.sh  (reusing pre-computed perms)"
echo "  sheet         : $SHEET"
echo "  out           : $OUT"
echo "  conditions    : ${CONDS[*]}"
echo "  n-perms       : $N_PERMS"
echo "  perm template : $PERM_TEMPLATE"
echo "  obs qcat col  : $OBS_QCAT_COL"
echo "  fdr           : $FDR    min-signal: $MIN_SIGNAL"
echo "  scripts-dir   : $SCRIPTS_DIR"
[[ $DRY_RUN -eq 1 ]] && echo "  MODE          : DRY RUN"
echo "============================================================"

printf 'condition\tn_bins\tlambda\tpct_at_floor\tbh_significant\tverdict\n' > "$AGG"

# ---------------------------------------------------------------------------
# Per-condition pipeline
# ---------------------------------------------------------------------------
for cond in "${CONDS[@]}"; do
    echo ""
    echo ">>> condition: $cond"

    # Pull the two replicate rows (sample, real qcat) for this condition.
    if ! pair=$(awk -F'\t' -v cond="$cond" -v qcol="$OBS_QCAT_COL" '
            NR==1{ for(i=1;i<=NF;i++){h=tolower($i)
                     if(h=="sample"||h=="name")si=i; else if(h=="condition")ci=i; else if(h==qcol)qi=i}
                   if(!si||!ci||!qi) exit 3; next }
            $0!~/[^[:space:]]/ || $0 ~ /^[[:space:]]*#/ {next}
            $ci==cond { n++; if(n==1){s1=$si;q1=$qi} else if(n==2){s2=$si;q2=$qi} }
            END{ if(n!=2) exit 4; printf "%s\t%s\t%s\t%s\n", s1,q1,s2,q2 }' "$SHEET"); then
        echo "  WARNING: $cond lacks exactly 2 replicates with a '$OBS_QCAT_COL' column; skipping" >&2
        continue
    fi
    IFS=$'\t' read -r s1 q1 s2 q2 <<< "$pair"
    echo "  replicates: $s1 ($q1) , $s2 ($q2)"

    cond_dir="$OUT/$cond"; obs_dir="$cond_dir/observed"; perm_root="$cond_dir/perm"
    mkdir -p "$obs_dir" "$perm_root"

    # --- observed diff (real rep1 vs real rep2) ------------------------------
    if [[ $DRY_RUN -eq 0 && ( ! -f "$q1" || ! -f "$q2" ) ]]; then
        echo "  ERROR: real qcat missing ($q1 or $q2); skipping $cond" >&2; continue
    fi
    existing_obs="$(locate_one_diff "$obs_dir")"
    if [[ -n "$existing_obs" && $FORCE -eq 0 ]]; then
        echo "  [observed] exists, skipping (use --force)"
    else
        write_qcat_sheet "$SHEETS_DIR/obs_${cond}.tsv" "$s1" "$q1" "$s2" "$q2"
        run "$PY" "$COMPARE_QCAT" --sheet "$SHEETS_DIR/obs_${cond}.tsv" --out "$obs_dir" --diff-only
    fi

    # --- permuted diffs (perm rep1 vs perm rep2, per round) ------------------
    declare -a nulls=()
    for ((k=1; k<=N_PERMS; k++)); do
        pq1="$(expand_template "$PERM_TEMPLATE" "$s1" "$k")"
        pq2="$(expand_template "$PERM_TEMPLATE" "$s2" "$k")"
        if [[ $DRY_RUN -eq 0 && ( ! -f "$pq1" || ! -f "$pq2" ) ]]; then
            echo "  WARNING: perm $k qcat missing ($pq1 or $pq2); skipping this round" >&2
            continue
        fi
        pk_dir="$perm_root/perm${k}"; mkdir -p "$pk_dir"
        existing_pk="$(locate_one_diff "$pk_dir")"
        if [[ -n "$existing_pk" && $FORCE -eq 0 ]]; then
            echo "  [perm $k] exists, skipping (use --force)"
        else
            write_qcat_sheet "$SHEETS_DIR/perm${k}_${cond}.tsv" "$s1" "$pq1" "$s2" "$pq2"
            run "$PY" "$COMPARE_QCAT" --sheet "$SHEETS_DIR/perm${k}_${cond}.tsv" --out "$pk_dir" --diff-only
        fi
        if [[ $DRY_RUN -eq 0 ]]; then
            d="$(locate_one_diff "$pk_dir")"
            if [[ -n "$d" ]]; then nulls+=("$d"); else echo "  WARNING: no diff produced for perm $k of $cond" >&2; fi
        fi
    done

    if [[ $DRY_RUN -eq 1 ]]; then
        echo "  [calibration] (dry-run) would diff observed + up to $N_PERMS perm rounds, then run bearing_calibration.py"
        continue
    fi

    obs_diff="$(locate_one_diff "$obs_dir")"
    [[ -n "$obs_diff" ]] || { echo "  ERROR: no observed diff for $cond; skipping" >&2; continue; }
    [[ ${#nulls[@]} -gt 0 ]] || { echo "  ERROR: no permuted diff nulls for $cond; skipping" >&2; continue; }
    [[ ${#nulls[@]} -eq $N_PERMS ]] || echo "  WARNING: using ${#nulls[@]} perm nulls (expected $N_PERMS) for $cond" >&2

    # --- calibration ---------------------------------------------------------
    out_prefix="$cond_dir/$cond"
    run "$PY" "$BEARING_CALIB" \
        --diff-qcat "$obs_diff" \
        --null-qcat "${nulls[@]}" \
        --out-prefix "$out_prefix" \
        --min-signal "$MIN_SIGNAL" \
        --fdr "$FDR" \
        --label "$cond rep1 vs rep2" \
        --bearing-pvalue "$BEARING_PVALUE"

    summ="${out_prefix}_summary.txt"
    if [[ -f "$summ" ]]; then
        nbins=$(extract 'observed bins'        "$summ" 's/^[^:]*:[[:space:]]*//')
        lam=$(extract 'genomic-inflation lam'  "$summ" 's/^[^:]*:[[:space:]]*//')
        floorpct=$(extract 'bins at floor'     "$summ" 's/.*\(([0-9.]+)%\).*/\1/')
        bhsig=$(extract 'BH-significant'       "$summ" 's/^[^:]*:[[:space:]]*([0-9,]+).*/\1/')
        verdict=$(extract '^VERDICT:'          "$summ" 's/^VERDICT:[[:space:]]*([A-Z]+).*/\1/')
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$cond" "$nbins" "$lam" "$floorpct" "$bhsig" "$verdict" >> "$AGG"
    else
        printf '%s\tNA\tNA\tNA\tNA\tNO_SUMMARY\n' "$cond" >> "$AGG"
    fi
    unset nulls
done

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
echo ""
echo "============================================================"
echo "  Calibration summary (all conditions)"
echo "============================================================"
if command -v column >/dev/null 2>&1; then column -t -s $'\t' "$AGG"; else cat "$AGG"; fi
echo ""
echo "Per-condition outputs: $OUT/<condition>/<condition>_calibration.pdf (+ _perbin.tsv.gz, _summary.txt)"
echo "Aggregate table:       $AGG"
if [[ $DRY_RUN -eq 1 ]]; then echo "(dry run: no commands were executed)"; fi
exit 0
