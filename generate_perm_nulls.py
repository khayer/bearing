#!/usr/bin/env python3
"""
generate_perm_nulls.py
======================
Orchestrate the full permutation null pipeline for bearing_pvalue.py.

Runs shift_bigwig.py, bigwig_to_qcat.py, and optionally compare_qcat.py
across N permutation rounds to produce null qcat.bgz files for both
standard (per-sample) and differential (condition A vs B) significance
testing.

WHAT THIS SCRIPT DOES
---------------------
For each permutation round i = 1..N:

  Standard null (per sample):
    1. Circularly shift all BigWig tracks (independently per track).
    2. Score the shifted tracks with bigwig_to_qcat.py.
    3. The resulting null_permi.qcat.bgz files are used directly by
       bearing_pvalue.py --null-qcat for per-sample significance.

  Differential null (--diff-sheet):
    After steps 1-2, also run compare_qcat.py on the permuted scores
    to generate permuted diff qcats. These are used by bearing_pvalue.py
    --diff --null-qcat for differential significance.

HOW MANY PERMUTATIONS?
----------------------
The minimum achievable p-value is 1/(n_null_bins + 1).
With mm10 at 200 bp (~11M scored bins per sample per permutation):

  N=3 : min p ~ 3e-8   (sufficient for FDR 0.05 in most cases)
  N=5 : min p ~ 2e-8   (comfortable for FDR 0.01)
  N=10: min p ~ 9e-9   (robust; recommended for publication)

Computational cost scales linearly with N and with the number of samples.
For 4 samples, N=3 requires 12 bigwig_to_qcat.py runs (~2-6 hours total
depending on hardware).

USAGE
-----
  # Standard per-sample nulls only:
  python generate_perm_nulls.py \\
    --sheet samples.tsv \\
    --n-perms 5 \\
    --out-dir perm_nulls/

  # Standard + differential nulls:
  python generate_perm_nulls.py \\
    --sheet samples.tsv \\
    --diff-sheet samples.tsv \\
    --n-perms 5 \\
    --out-dir perm_nulls/

  # Dry run (print commands without executing):
  python generate_perm_nulls.py \\
    --sheet samples.tsv \\
    --n-perms 3 \\
    --out-dir perm_nulls/ \\
    --dry-run

AFTER RUNNING
-------------
The script prints the exact bearing_pvalue.py commands for every sample
and every condition pair once all permutations are complete.

Example output structure for 3 permutations, 2 conditions:

  perm_nulls/
    perm1/
      DN_rep1/
        atac_perm1.bw  ctcf_perm1.bw  ...    <- shifted BigWigs
        null_perm1.qcat.bgz                  <- permuted scoring
      DN_rep2/  3T3_rep1/  3T3_rep2/  ...
      diff_comparison/
        diff_DN_vs_3T3.qcat.bgz              <- permuted diff (if --diff-sheet)
    perm2/  perm3/

  bearing_pvalue.py commands emitted:
    --qcat DN_rep1.qcat.bgz
    --null-qcat perm_nulls/perm1/DN_rep1/null_perm1.qcat.bgz
                perm_nulls/perm2/DN_rep1/null_perm2.qcat.bgz
                perm_nulls/perm3/DN_rep1/null_perm3.qcat.bgz

SAMPLE SHEET FORMAT
-------------------
Same TSV format as bigwig_to_qcat.py and shift_bigwig.py:

  sample    condition   replicate   bw                          out
  DN_rep1   DN          1           atac.bw,ctcf.bw,...         DN_rep1.qcat.bgz
  DN_rep2   DN          2           atac2.bw,ctcf2.bw,...       DN_rep2.qcat.bgz
  3T3_rep1  3T3         1           atac_3t3.bw,...             3T3_rep1.qcat.bgz

For --diff-sheet, the same sheet is typically used (it reads the condition
column to determine which diffs to compute).

DEPENDENCIES
------------
  pip install pyBigWig numpy
  (bigwig_to_qcat.py, compare_qcat.py, shift_bigwig.py must be in PATH or
   same directory)
"""

import argparse
import concurrent.futures
import csv
import datetime
import hashlib
import json
import subprocess
import sys
from itertools import combinations
from pathlib import Path

from bearing.sheet import extract_bw_values, sample_name_from_row
from bearing.tsv import read_tsv_table
from bearing.runner import format_command


# ---------------------------------------------------------------------------
# Sheet parsing
# ---------------------------------------------------------------------------

def load_sheet(path):
    """
    Parse TSV. Returns list of dicts with keys:
      sample, condition, replicate, bw_paths (list of Path), out (str or None)
    """
    samples = []
    fields, rows = read_tsv_table(Path(path))
    if "bw" not in fields and not any(
        h.startswith("bw") and h[2:].isdigit() for h in fields
    ):
        sys.exit(
            "ERROR: sheet must have a 'bw' column (comma-separated paths) or bw1, bw2, ... columns."
        )

    for row in rows:
        sample = sample_name_from_row(row)
        condition = row.get("condition", "")
        replicate = row.get("replicate", "0") or "0"
        out = row.get("out") or None
        bw_list = extract_bw_values(row)

        if not bw_list:
            print(
                f"  WARNING: no BigWig paths for sample '{sample}', skipping.",
                file=sys.stderr,
            )
            continue

        samples.append({
            "sample": sample,
            "condition": condition,
            "replicate": replicate,
            "bw_paths": [Path(p) for p in bw_list],
            "out": out,
        })

    if not samples:
        sys.exit("ERROR: no usable rows in sheet.")
    return samples


def qcat_name_from_row(row):
    """Derive the original qcat filename for a sample row."""
    if row["out"]:
        return row["out"]
    if row["sample"]:
        return row["sample"] + ".qcat.bgz"
    return "sample.qcat.bgz"


def null_qcat_name(row, suffix):
    """Derive the null qcat filename for a permuted sample."""
    base = qcat_name_from_row(row)
    for ext in (".qcat.bgz", ".bgz"):
        if base.endswith(ext):
            base = base[:-len(ext)]
            break
    return Path(base).name + suffix + ".qcat.bgz"


def _null_qcat_candidates(row, suffix):
    """Return likely null qcat basenames for one sample/perm suffix."""
    sample = row.get("sample", "")
    suffix_clean = suffix[1:] if suffix.startswith("_") else suffix
    preferred = null_qcat_name(row, suffix)

    candidates = [
        preferred,
        f"null_{suffix_clean}.qcat.bgz",
        f"{sample}_{suffix_clean}.qcat.bgz" if sample else "",
        f"{sample}{suffix}.qcat.bgz" if sample else "",
    ]
    # De-duplicate while preserving order.
    out = []
    seen = set()
    for c in candidates:
        if c and c not in seen:
            out.append(c)
            seen.add(c)
    return out


def resolve_null_qcat_path(row, perm_dir, suffix):
    """Find an existing null qcat path for one sample/perm if available."""
    sample = row["sample"]
    sample_perm_dir = perm_dir / sample

    for name in _null_qcat_candidates(row, suffix):
        p = sample_perm_dir / name
        if p.exists():
            return p

    # Last-resort fallback: if exactly one qcat exists in the sample perm dir,
    # use it to keep skip-scoring workflows robust across naming conventions.
    matches = sorted(sample_perm_dir.glob("*.qcat.bgz"))
    if len(matches) == 1:
        return matches[0]

    return sample_perm_dir / null_qcat_name(row, suffix)


def _expected_shifted_bws(row, perm_dir, suffix):
    sample_perm_dir = perm_dir / row["sample"]
    return [sample_perm_dir / (p.stem + suffix + ".bw") for p in row["bw_paths"]]


def _all_shifted_exist(samples, perm_dir, suffix):
    return all(
        all(p.exists() for p in _expected_shifted_bws(row, perm_dir, suffix))
        for row in samples
    )


def _all_sample_nulls_exist(samples, perm_dir, suffix):
    for row in samples:
        p = resolve_null_qcat_path(row, perm_dir, suffix)
        if not p.exists():
            return False
    return True


def _all_diff_nulls_exist(diff_samples, diff_perm_dir):
    conditions = sorted(set(s["condition"] for s in diff_samples))
    for cond_A, cond_B in combinations(conditions, 2):
        safe_A = cond_A.replace(" ", "_")
        safe_B = cond_B.replace(" ", "_")
        diff_path = diff_perm_dir / f"diff_{safe_A}_vs_{safe_B}.qcat.bgz"
        if not diff_path.exists():
            return False
    return True


def _deterministic_file_seed(perm_seed, sample, in_path, out_path):
    raw = f"{perm_seed}|{sample}|{in_path}|{out_path}".encode("utf-8")
    return int(hashlib.sha256(raw).hexdigest()[:8], 16)


def _utc_now_iso():
    return datetime.datetime.now(datetime.timezone.utc).replace(microsecond=0).isoformat()


def _count_shifted_existing(samples, perm_dir, suffix):
    return sum(
        1
        for row in samples
        for p in _expected_shifted_bws(row, perm_dir, suffix)
        if p.exists()
    )


def _count_null_qcat_existing(samples, perm_dir, suffix):
    return sum(
        1
        for row in samples
        if resolve_null_qcat_path(row, perm_dir, suffix).exists()
    )


def _count_diff_qcat_existing(diff_samples, diff_perm_dir):
    if not diff_samples:
        return 0
    n = 0
    conditions = sorted(set(s["condition"] for s in diff_samples))
    for cond_A, cond_B in combinations(conditions, 2):
        safe_A = cond_A.replace(" ", "_")
        safe_B = cond_B.replace(" ", "_")
        diff_path = diff_perm_dir / f"diff_{safe_A}_vs_{safe_B}.qcat.bgz"
        if diff_path.exists():
            n += 1
    return n


def _write_progress_state(progress_json, state):
    progress_json.parent.mkdir(parents=True, exist_ok=True)
    progress_json.write_text(json.dumps(state, indent=2, sort_keys=True) + "\n")


def _append_progress_log(progress_log, msg):
    progress_log.parent.mkdir(parents=True, exist_ok=True)
    with open(progress_log, "a", encoding="utf-8") as fh:
        fh.write(f"[{_utc_now_iso()}] {msg}\n")


# ---------------------------------------------------------------------------
# Command runner
# ---------------------------------------------------------------------------

def run_cmd(cmd, dry_run=False, label=""):
    """Print and optionally execute a shell command."""
    cmd_str = format_command(cmd)
    if label:
        print(f"\n  [{label}]")
    print(f"  $ {cmd_str}")
    if not dry_run:
        result = subprocess.run(cmd)
        if result.returncode != 0:
            sys.exit(f"ERROR: command failed (exit {result.returncode}): "
                     f"{cmd_str}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Orchestrate permutation null generation for bearing_pvalue.py.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument("--sheet", required=True, metavar="TSV",
                    help="Sample sheet TSV (same format as bigwig_to_qcat.py).")
    ap.add_argument("--diff-sheet", default=None, metavar="TSV",
                    help="If given, also run compare_qcat.py on permuted "
                         "scores to generate null diff qcats. Typically the "
                         "same TSV as --sheet (reads the condition column).")
    ap.add_argument("--n-perms", type=int, default=5, metavar="N",
                    help="Number of permutation rounds (default: 5). "
                         "See docstring for guidance on choosing N.")
    ap.add_argument("--out-dir", required=True, metavar="DIR",
                    help="Root output directory. Permutation rounds are "
                         "written to out-dir/perm1/, out-dir/perm2/, etc.")
    ap.add_argument("--seed", type=int, default=42, metavar="INT",
                    help="Master random seed (default: 42). Each permutation "
                         "round uses a derived seed for reproducibility.")
    ap.add_argument("--perm-index", type=int, default=None, metavar="I",
                    help="Run ONLY permutation round I (1-based) instead of "
                         "1..n_perms. Uses the same derived seed (seed+I) and "
                         "writes to out-dir/permI/, so independent jobs produce "
                         "identical output to a single batch run. Enables "
                         "per-permutation parallelism across cluster jobs.")
    ap.add_argument("--only-sample", default=None, metavar="NAME",
                    help="Process ONLY this sample (by sample name). Combined "
                         "with --perm-index, lets one cluster job shift+score a "
                         "single sample for a single permutation, for maximum "
                         "parallelism. The differential step still needs all "
                         "samples present, so run it (without --only-sample) "
                         "after every sample's observed null exists.")
    ap.add_argument("--min-shift", type=int, default=1_000_000, metavar="BP",
                    help="Minimum circular shift in bp (default: 1,000,000).")
    ap.add_argument("--chroms", nargs="+", metavar="CHR", default=None,
                    help="Restrict scoring to specific chromosomes "
                         "(passed to bigwig_to_qcat.py --chroms). "
                         "Useful for quick tests.")
    ap.add_argument("--chrom-sizes", default=None, metavar="FILE",
                    help="Custom chrom.sizes file (passed to bigwig_to_qcat.py).")
    ap.add_argument("--categories", default=None, metavar="YAML",
                    help="YAML category file (passed to bigwig_to_qcat.py "
                         "and compare_qcat.py).")
    ap.add_argument("--blacklist", default=None, metavar="BED",
                    help="Blacklist BED passed to bigwig_to_qcat.py --blacklist during perm scoring.")
    ap.add_argument("--floors-tsv", default=None, metavar="PATH",
                    help="Pre-computed per-track noise floors TSV (from a "
                         "prior real-data --write-floors-tsv run). Passed "
                         "through to bigwig_to_qcat.py --floors-tsv for each "
                         "permuted sample, which skips the percentile-floor "
                         "sampling pass (~4 min saved per perm scoring run). "
                         "Floors are invariant under circular shift so reuse "
                         "is statistically sound.")
    ap.add_argument("--min-signal", type=float, default=None, metavar="FLOAT",
                    help="Passed through to bigwig_to_qcat.py --min-signal "
                         "during perm scoring (default: bigwig_to_qcat.py's "
                         "own default).")
    # Normalization must match the OBSERVED scoring so the null is on the same
    # scale. Circular shifting preserves each track's value multiset, so
    # quantile mapping of a shifted track equals the shifted normalized track;
    # applying the same normalization here is exact, not approximate.
    ap.add_argument("--normalize-tracks", action="store_true",
                    help="Pass --normalize-tracks to bigwig_to_qcat.py during "
                         "perm scoring (must match the observed run).")
    ap.add_argument("--normalize-method", default=None, metavar="METHOD",
                    help="Pass through to bigwig_to_qcat.py --normalize-method "
                         "during perm scoring (e.g. cohort-quantile).")
    ap.add_argument("--cohort-reference", default=None, metavar="NPZ",
                    help="Pass through to bigwig_to_qcat.py --cohort-reference "
                         "during perm scoring (required for cohort-quantile).")
    ap.add_argument("--score-method", default=None, choices=["kl", "jsd"],
                    help="Pass through to bigwig_to_qcat.py --score-method during "
                         "perm scoring, so the null matches the observed scoring.")
    ap.add_argument("--no-extras", action="store_true",
                    help="Pass --no-extras to bigwig_to_qcat.py during perm "
                         "scoring (skips cats.json, tracks.ini, plots). "
                         "Recommended for permutation runs; saves ~30 s per "
                         "scoring run plus disk I/O.")
    ap.add_argument("--jobs", type=int, default=1, metavar="N",
                    help="Parallel jobs for bigwig_to_qcat.py (default: 1).")
    ap.add_argument("--score-workers", type=int, default=1, metavar="N",
                    help="Parallel sample processes for bigwig_to_qcat.py per permutation (default: 1).")
    ap.add_argument("--shift-workers", type=int, default=1, metavar="N",
                    help="Parallel workers for shift_bigwig.py (default: 1).")
    ap.add_argument("--compare-workers", type=int, default=1, metavar="N",
                    help="Parallel workers for compare_qcat.py during "
                         "differential null generation (default: 1).")
    ap.add_argument("--compare-consensus-q", action="store_true",
                    help="Pass --consensus-q to compare_qcat.py when "
                         "generating permuted differential qcats.")
    ap.add_argument("--dry-run", action="store_true",
                    help="Print commands without executing them.")
    ap.add_argument("--skip-scoring", action="store_true",
                    help="Skip BigWig shifting and scoring; only run "
                         "compare_qcat.py on already-generated null qcats. "
                         "Useful if shifting is done but diffing was missed.")
    ap.add_argument("--no-resume", action="store_true",
                    help="Disable resume behavior. By default, completed "
                         "permutations and existing outputs are skipped.")
    args = ap.parse_args()
    resume = not args.no_resume

    # Resolve script locations (same directory as this script, or PATH)
    script_dir = Path(__file__).parent
    def _script(name):
        local = script_dir / name
        return str(local) if local.exists() else name

    SHIFT_BW  = _script("shift_bigwig.py")
    SCORE_BW  = _script("bigwig_to_qcat.py")
    COMPARE   = _script("compare_qcat.py")
    PVALUE    = _script("bearing_pvalue.py")

    out_root = Path(args.out_dir)
    out_root.mkdir(parents=True, exist_ok=True)

    # When running a single permutation (--perm-index, used by the per-perm
    # Snakemake jobs), write per-index progress files so concurrent jobs never
    # write the same file. A shared progress.json/.log would be clobbered by
    # 10 parallel jobs and can hit partial-write races on shared filesystems.
    if args.perm_index is not None:
        _ptag = "perm%d" % args.perm_index
        if args.only_sample:
            _ptag += "." + args.only_sample
        progress_json = out_root / ("progress.%s.json" % _ptag)
        progress_log = out_root / ("progress.%s.log" % _ptag)
    else:
        progress_json = out_root / "progress.json"
        progress_log = out_root / "progress.log"

    samples = load_sheet(args.sheet)
    if args.only_sample:
        samples = [s for s in samples if s["sample"] == args.only_sample]
        if not samples:
            sys.exit("ERROR: --only-sample '%s' not found in sheet"
                     % args.only_sample)
    diff_samples = load_sheet(args.diff_sheet) if args.diff_sheet else None

    # Collect null qcat paths per sample for final summary
    # null_qcats[sample_name] = [path1, path2, ...]
    null_qcats = {s["sample"]: [] for s in samples}
    # null_diff_qcats[(cond_A, cond_B)] = [path1, path2, ...]
    null_diff_qcats = {}
    if diff_samples:
        conditions = sorted(set(s["condition"] for s in diff_samples))
        for cond_A, cond_B in combinations(conditions, 2):
            null_diff_qcats[(cond_A, cond_B)] = []

    total_shifted_expected = sum(len(s["bw_paths"]) for s in samples)
    total_qcat_expected = len(samples)
    total_diff_expected = len(null_diff_qcats)
    progress_state = {
        "started_at": _utc_now_iso(),
        "sheet": args.sheet,
        "diff_sheet": args.diff_sheet,
        "out_dir": str(out_root),
        "n_perms": args.n_perms,
        "resume_enabled": resume,
        "permutation_template": {
            "shift_expected": total_shifted_expected,
            "qcat_expected": total_qcat_expected,
            "diff_expected": total_diff_expected,
        },
        "permutations": {
            str(i): {
                "status": "pending",
                "shift_done": 0,
                "qcat_done": 0,
                "diff_done": 0,
                "updated_at": _utc_now_iso(),
            }
            for i in range(1, args.n_perms + 1)
        },
        "last_message": "initialized",
        "updated_at": _utc_now_iso(),
    }

    def _update_perm_progress(perm_i, *, status=None, message=None, phase=None, target=None):
        suffix_local = f"_perm{perm_i}"
        perm_dir_local = out_root / f"perm{perm_i}"
        diff_dir_local = perm_dir_local / "diff_comparison"
        node = progress_state["permutations"][str(perm_i)]
        node["shift_done"] = _count_shifted_existing(samples, perm_dir_local, suffix_local)
        node["qcat_done"] = _count_null_qcat_existing(samples, perm_dir_local, suffix_local)
        node["diff_done"] = _count_diff_qcat_existing(diff_samples, diff_dir_local)
        if status is not None:
            node["status"] = status
        if phase is not None:
            node["phase"] = phase
        if target is not None:
            node["target"] = target
        node["updated_at"] = _utc_now_iso()
        if message:
            progress_state["last_message"] = message
            _append_progress_log(progress_log, message)
        progress_state["updated_at"] = _utc_now_iso()
        _write_progress_state(progress_json, progress_state)

    _append_progress_log(progress_log, "run initialized")
    _write_progress_state(progress_json, progress_state)

    # ── Banner ────────────────────────────────────────────────────────────
    print("=" * 64)
    print("  generate_perm_nulls.py")
    print(f"  Sheet:       {args.sheet}")
    print(f"  Samples:     {len(samples)}")
    print(f"  Permutations:{args.n_perms}")
    print(f"  Output:      {out_root}/")
    print(f"  Score workers: {max(1, args.score_workers)}")
    print(f"  Shift workers: {max(1, args.shift_workers)}")
    if args.floors_tsv:
        print(f"  Floors TSV (reused for all perms): {args.floors_tsv}")
    if args.no_extras:
        print(f"  --no-extras: enabled for perm scoring")
    print(f"  Progress JSON: {progress_json}")
    print(f"  Progress log:  {progress_log}")
    if args.diff_sheet:
        print(f"  Diff sheet:  {args.diff_sheet}")
        print(f"  Compare workers: {max(1, args.compare_workers)}")
        if args.compare_consensus_q:
            print("  Compare consensus Q: enabled")
        if null_diff_qcats:
            pairs_str = ", ".join(
                f"{a} vs {b}" for a, b in null_diff_qcats)
            print(f"  Diff pairs:  {pairs_str}")
    if args.dry_run:
        print("  MODE:        DRY RUN (no commands executed)")
    print("=" * 64)

    # ── Permutation rounds ────────────────────────────────────────────────
    if args.perm_index is not None:
        if args.perm_index < 1 or args.perm_index > args.n_perms:
            sys.exit("ERROR: --perm-index %d out of range 1..%d"
                     % (args.perm_index, args.n_perms))
        _perm_rounds = [args.perm_index]
    else:
        _perm_rounds = list(range(1, args.n_perms + 1))
    for perm_i in _perm_rounds:
        suffix = f"_perm{perm_i}"
        perm_seed = args.seed + perm_i  # deterministic per-round seed
        perm_dir = out_root / f"perm{perm_i}"
        perm_dir.mkdir(parents=True, exist_ok=True)
        diff_perm_dir = perm_dir / "diff_comparison"

        print(f"\n{'='*64}")
        print(f"  Permutation {perm_i} / {args.n_perms}  "
              f"(seed={perm_seed}, suffix={suffix})")
        print(f"{'='*64}")
        _update_perm_progress(
            perm_i,
            status="running",
            phase="initializing",
            target=None,
            message=f"perm{perm_i} started",
        )

        if resume:
            score_done = _all_sample_nulls_exist(samples, perm_dir, suffix)
            diff_done = True
            if args.diff_sheet:
                diff_done = diff_perm_dir.exists() and _all_diff_nulls_exist(
                    diff_samples, diff_perm_dir
                )
            if score_done and (not args.diff_sheet or diff_done):
                print(f"\n  [resume] perm{perm_i} already complete; skipping.")
                for row in samples:
                    sample = row["sample"]
                    null_qcats[sample].append(
                        resolve_null_qcat_path(row, perm_dir, suffix)
                    )
                if args.diff_sheet and diff_samples:
                    conditions = sorted(set(s["condition"] for s in diff_samples))
                    for cond_A, cond_B in combinations(conditions, 2):
                        safe_A = cond_A.replace(" ", "_")
                        safe_B = cond_B.replace(" ", "_")
                        diff_path = (diff_perm_dir /
                                     f"diff_{safe_A}_vs_{safe_B}.qcat.bgz")
                        null_diff_qcats[(cond_A, cond_B)].append(diff_path)
                _update_perm_progress(
                    perm_i,
                    status="complete",
                    phase="complete",
                    message=f"perm{perm_i} resume-skip complete",
                )
                continue

        if not args.skip_scoring:
            # ── Step 1: shift all BigWigs via sheet ───────────────────────
            print(f"\nStep 1/{2 + bool(args.diff_sheet)}: "
                  f"Shifting BigWigs...")
            shift_missing = []
            for row in samples:
                shifted_paths = _expected_shifted_bws(row, perm_dir, suffix)
                for in_path, out_path in zip(row["bw_paths"], shifted_paths):
                    if not out_path.exists():
                        shift_missing.append((row["sample"], in_path, out_path))

            if resume and not shift_missing:
                print("  [resume] all shifted BigWigs already exist; skipping shift step.")
                _update_perm_progress(
                    perm_i,
                    status="running",
                    phase="shifting",
                    message=f"perm{perm_i} shift step skipped (already complete)",
                )
            elif resume and shift_missing:
                total_expected = sum(len(r["bw_paths"]) for r in samples)
                if len(shift_missing) < total_expected:
                    print(
                        f"  [resume] shifting only missing files: "
                        f"{len(shift_missing)}/{total_expected}."
                    )

                    def _run_shift_missing(task):
                        sample, in_path, out_path = task
                        out_path.parent.mkdir(parents=True, exist_ok=True)
                        file_seed = _deterministic_file_seed(
                            perm_seed, sample, in_path, out_path
                        )
                        cmd = [
                            "python", SHIFT_BW,
                            "--bw", str(in_path),
                            "--out", str(out_path),
                            "--seed", str(file_seed),
                            "--min-shift", str(args.min_shift),
                        ]
                        if args.dry_run:
                            run_cmd(cmd, dry_run=True, label=f"shift {sample}")
                            return 0, format_command(cmd)
                        result = subprocess.run(cmd)
                        return result.returncode, format_command(cmd)

                    if args.dry_run:
                        for task in shift_missing:
                            _run_shift_missing(task)
                    else:
                        failures = []
                        max_workers = min(max(1, int(args.shift_workers)), len(shift_missing))
                        if max_workers == 1:
                            for task in shift_missing:
                                code, cmd_str = _run_shift_missing(task)
                                if code != 0:
                                    failures.append((code, cmd_str))
                                    break
                        else:
                            with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as ex:
                                futures = [ex.submit(_run_shift_missing, t) for t in shift_missing]
                                for fut in concurrent.futures.as_completed(futures):
                                    code, cmd_str = fut.result()
                                    if code != 0:
                                        failures.append((code, cmd_str))
                            
                        if failures:
                            code, cmd_str = failures[0]
                            sys.exit(
                                f"ERROR: command failed (exit {code}) during partial resume shifting: {cmd_str}"
                            )
                    _update_perm_progress(
                        perm_i,
                        status="running",
                        phase="shifting",
                        message=f"perm{perm_i} partial resume shifting complete",
                    )
                else:
                    shift_cmd = [
                        "python", SHIFT_BW,
                        "--sheet", args.sheet,
                        "--out-dir", str(perm_dir),
                        "--suffix", suffix,
                        "--seed", str(perm_seed),
                        "--min-shift", str(args.min_shift),
                        "--workers", str(max(1, args.shift_workers)),
                    ]
                    run_cmd(shift_cmd, dry_run=args.dry_run,
                            label=f"shift perm{perm_i}")
                    _update_perm_progress(
                        perm_i,
                        status="running",
                        phase="shifting",
                        message=f"perm{perm_i} shift step complete",
                    )
            else:
                shift_cmd = [
                    "python", SHIFT_BW,
                    "--sheet", args.sheet,
                    "--out-dir", str(perm_dir),
                    "--suffix", suffix,
                    "--seed", str(perm_seed),
                    "--min-shift", str(args.min_shift),
                    "--workers", str(max(1, args.shift_workers)),
                ]
                run_cmd(shift_cmd, dry_run=args.dry_run,
                        label=f"shift perm{perm_i}")
                _update_perm_progress(
                    perm_i,
                    status="running",
                    phase="shifting",
                    message=f"perm{perm_i} shift step complete",
                )

            # ── Step 2: score each permuted sample ────────────────────────
            print(f"\nStep 2/{2 + bool(args.diff_sheet)}: "
                  f"Scoring permuted BigWigs...")
            score_tasks = []
            for row in samples:
                sample = row["sample"]
                sample_perm_dir = perm_dir / sample
                # Shifted BigWigs land in perm_dir/<sample>/
                shifted_bws = [
                    sample_perm_dir / (p.stem + suffix + ".bw")
                    for p in row["bw_paths"]
                ]
                null_qcat = sample_perm_dir / null_qcat_name(row, suffix)

                score_cmd = [
                    "python", SCORE_BW,
                    "--bw"] + [str(p) for p in shifted_bws] + [
                    "--out", str(null_qcat),
                ]
                if args.chroms:
                    score_cmd += ["--chroms"] + args.chroms
                if args.chrom_sizes:
                    score_cmd += ["--chrom-sizes", args.chrom_sizes]
                if args.categories:
                    score_cmd += ["--categories", args.categories]
                if args.blacklist:
                    score_cmd += ["--blacklist", args.blacklist]
                if args.floors_tsv:
                    score_cmd += ["--floors-tsv", args.floors_tsv,
                                  "--sample-name", sample]
                if args.min_signal is not None:
                    score_cmd += ["--min-signal", str(args.min_signal)]
                if args.normalize_tracks:
                    score_cmd += ["--normalize-tracks"]
                if args.normalize_method:
                    score_cmd += ["--normalize-method", args.normalize_method]
                if args.cohort_reference:
                    score_cmd += ["--cohort-reference", args.cohort_reference]
                if args.score_method:
                    score_cmd += ["--score-method", args.score_method]
                if args.no_extras:
                    score_cmd += ["--no-extras"]
                if args.jobs > 1:
                    score_cmd += ["--jobs", str(args.jobs)]

                label = f"score {sample} perm{perm_i}"

                if resume:
                    existing = resolve_null_qcat_path(row, perm_dir, suffix)
                    if existing.exists():
                        print(f"  [resume] existing null qcat found for {sample}; skipping score.")
                        null_qcats[sample].append(existing)
                        continue

                score_tasks.append((sample, null_qcat, score_cmd, label))

            _update_perm_progress(
                perm_i,
                status="running",
                phase="scoring",
                target=f"{len(score_tasks)}/{len(samples)} samples queued" if score_tasks else "no samples queued",
                message=f"perm{perm_i} scoring queued: {len(score_tasks)} sample(s)",
            )

            if args.dry_run:
                for sample, null_qcat, score_cmd, label in score_tasks:
                    run_cmd(score_cmd, dry_run=True, label=label)
                    null_qcats[sample].append(null_qcat)
                    _update_perm_progress(
                        perm_i,
                        status="running",
                        phase="scoring",
                        target=f"{len(null_qcats[sample])} completed for {sample}",
                        message=f"perm{perm_i} score complete: {sample}",
                    )
            else:
                score_workers = max(1, int(args.score_workers))
                if score_workers == 1 or len(score_tasks) == 1:
                    for sample, null_qcat, score_cmd, label in score_tasks:
                        run_cmd(score_cmd, dry_run=False, label=label)
                        null_qcats[sample].append(null_qcat)
                        _update_perm_progress(
                            perm_i,
                            status="running",
                            phase="scoring",
                            target=f"completed {sample}",
                            message=f"perm{perm_i} score complete: {sample}",
                        )
                else:
                    for _, _, score_cmd, label in score_tasks:
                        run_cmd(score_cmd, dry_run=True, label=label)

                    def _run_score(task):
                        sample, null_qcat, score_cmd, label = task
                        result = subprocess.run(score_cmd)
                        return sample, null_qcat, label, result.returncode, format_command(score_cmd)

                    failures = []
                    max_workers = min(score_workers, len(score_tasks))
                    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as ex:
                        futures = [ex.submit(_run_score, t) for t in score_tasks]
                        for fut in concurrent.futures.as_completed(futures):
                            sample, null_qcat, label, code, cmd_str = fut.result()
                            if code != 0:
                                failures.append((label, code, cmd_str))
                            else:
                                null_qcats[sample].append(null_qcat)
                                _update_perm_progress(
                                    perm_i,
                                    status="running",
                                    phase="scoring",
                                    target=f"completed {sample}",
                                    message=f"perm{perm_i} score complete: {sample}",
                                )

                    if failures:
                        label, code, cmd_str = failures[0]
                        sys.exit(
                            f"ERROR: command failed (exit {code}) in {label}: {cmd_str}"
                        )
            _update_perm_progress(
                perm_i,
                status="running",
                phase="scoring",
                message=f"perm{perm_i} score step complete",
            )
        else:
            # In skip-scoring mode, discover existing null qcats so final
            # bearing_pvalue command suggestions are still emitted.
            for row in samples:
                sample = row["sample"]
                resolved = resolve_null_qcat_path(row, perm_dir, suffix)
                null_qcats[sample].append(resolved)
            _update_perm_progress(
                perm_i,
                status="running",
                phase="scoring",
                message=f"perm{perm_i} skip-scoring mode: discovered existing qcats",
            )

        # ── Step 3 (optional): generate permuted diff qcats ───────────────
        if args.diff_sheet:
            diff_perm_dir.mkdir(parents=True, exist_ok=True)

            if resume and _all_diff_nulls_exist(diff_samples, diff_perm_dir):
                print("\n  [resume] all permuted diff qcats already exist; skipping diff step.")
                conditions = sorted(set(s["condition"] for s in diff_samples))
                for cond_A, cond_B in combinations(conditions, 2):
                    safe_A = cond_A.replace(" ", "_")
                    safe_B = cond_B.replace(" ", "_")
                    diff_path = (diff_perm_dir /
                                 f"diff_{safe_A}_vs_{safe_B}.qcat.bgz")
                    null_diff_qcats[(cond_A, cond_B)].append(diff_path)
                _update_perm_progress(
                    perm_i,
                    status="complete",
                    phase="complete",
                    message=f"perm{perm_i} diff step skipped (already complete)",
                )
                continue

            # Write a temporary sheet pointing to the null qcats for this perm
            tmp_sheet = perm_dir / f"null_sheet_perm{perm_i}.tsv"
            _write_null_sheet(
                diff_samples, perm_dir, suffix, tmp_sheet,
                dry_run=args.dry_run)

            step_n = 3
            print(f"\nStep {step_n}/{2 + bool(args.diff_sheet)}: "
                  f"Generating permuted diff qcats...")
            diff_cmd = [
                "python", COMPARE,
                "--sheet", str(tmp_sheet),
                "--out", str(diff_perm_dir),
                "--diff-only",
                "--no-clip",
                "--skip-pca",
                "--skip-q-pair-jsd",
                "--workers", str(max(1, args.compare_workers)),
            ]
            if args.chroms:
                diff_cmd += ["--chroms"] + args.chroms
            if args.categories:
                diff_cmd += ["--categories", args.categories]
            if args.compare_consensus_q:
                diff_cmd += ["--consensus-q"]

            run_cmd(diff_cmd, dry_run=args.dry_run,
                    label=f"diff perm{perm_i}")

            # Collect produced diff qcat paths
            if diff_samples:
                conditions = sorted(set(s["condition"] for s in diff_samples))
                for cond_A, cond_B in combinations(conditions, 2):
                    safe_A = cond_A.replace(" ", "_")
                    safe_B = cond_B.replace(" ", "_")
                    diff_path = (diff_perm_dir /
                                 f"diff_{safe_A}_vs_{safe_B}.qcat.bgz")
                    null_diff_qcats[(cond_A, cond_B)].append(diff_path)

            _update_perm_progress(
                perm_i,
                status="complete",
                phase="complete",
                message=f"perm{perm_i} diff step complete",
            )
        else:
            _update_perm_progress(
                perm_i,
                status="complete",
                phase="complete",
                message=f"perm{perm_i} complete",
            )

    # ── Final summary: print bearing_pvalue.py commands ───────────────────
    print(f"\n{'='*64}")
    print("  All permutations complete.")
    print(f"{'='*64}")
    progress_state["completed_at"] = _utc_now_iso()
    progress_state["last_message"] = "all permutations complete"
    progress_state["updated_at"] = _utc_now_iso()
    _append_progress_log(progress_log, "all permutations complete")
    _write_progress_state(progress_json, progress_state)
    print("\nRun bearing_pvalue.py for each sample:\n")

    for row in samples:
        sample = row["sample"]
        orig_qcat = qcat_name_from_row(row)
        null_paths = null_qcats[sample]
        if not null_paths and not args.dry_run:
            continue
        if args.dry_run:
            # In dry run, show expected paths
            null_paths = [
                out_root / f"perm{i}" / sample /
                null_qcat_name(row, f"_perm{i}")
                for i in range(1, args.n_perms + 1)
            ]

        null_args = " \\\n    ".join(str(p) for p in null_paths)
        print(f"# {sample}")
        print(f"python {PVALUE} \\")
        print(f"  --qcat {orig_qcat} \\")
        print(f"  --null-qcat {null_args} \\")
        print(f"  --out-prefix results/{sample} \\")
        print(f"  --fdr 0.05 --score-plot")
        print()

    if args.diff_sheet and null_diff_qcats:
        print("Run bearing_pvalue.py --diff for each condition pair:\n")
        conditions = sorted(set(s["condition"] for s in diff_samples))
        for cond_A, cond_B in combinations(conditions, 2):
            safe_A = cond_A.replace(" ", "_")
            safe_B = cond_B.replace(" ", "_")
            orig_diff = f"diff_{safe_A}_vs_{safe_B}.qcat.bgz"
            null_paths = null_diff_qcats.get((cond_A, cond_B), [])
            if not null_paths and not args.dry_run:
                continue
            if args.dry_run:
                null_paths = [
                    out_root / f"perm{i}" / "diff_comparison" /
                    f"diff_{safe_A}_vs_{safe_B}.qcat.bgz"
                    for i in range(1, args.n_perms + 1)
                ]
            null_args = " \\\n    ".join(str(p) for p in null_paths)
            print(f"# {cond_A} vs {cond_B}")
            print(f"python {PVALUE} \\")
            print(f"  --qcat {orig_diff} \\")
            print(f"  --null-qcat {null_args} \\")
            print(f"  --diff \\")
            print(f"  --out-prefix results/diff_{safe_A}_vs_{safe_B} \\")
            print(f"  --fdr 0.05 --score-plot")
            print()


def _write_null_sheet(samples, perm_dir, suffix, out_path, dry_run=False):
    """
    Write a temporary sample sheet pointing to the permuted null qcats
    for use by compare_qcat.py --sheet. Only includes the qcat column
    (no bw paths needed since qcats already exist).
    """
    if dry_run:
        print(f"  [dry run] would write null sheet -> {out_path}")
        return

    with open(out_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["sample", "condition", "replicate", "qcat"])
        for row in samples:
            sample = row["sample"]
            null_qcat = resolve_null_qcat_path(row, perm_dir, suffix)
            writer.writerow([
                sample,
                row["condition"],
                row["replicate"],
                str(null_qcat),
            ])


if __name__ == "__main__":
    main()
