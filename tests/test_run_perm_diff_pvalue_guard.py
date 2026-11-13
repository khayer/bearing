# ----------------------------------------------------------------------
# File     : tests/test_run_perm_diff_pvalue_guard.py
# Version  : 1.0.0 (2026-11-11)
# Authors  : Katharina E. Hayer (katharinaehayer@gmail.com) and Claude
#            (Anthropic), co-created
# Purpose  : slurm/run_perm_diff_pvalue.sh must exit non-zero when a
#            comparison ends without a stats.tsv (hollow pvalue.done guard).
# ----------------------------------------------------------------------
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "slurm" / "run_perm_diff_pvalue.sh"

STUB_OK = """import sys
a = sys.argv
out = a[a.index("--out-prefix") + 1]
open(out + ".stats.tsv", "w").write("chrom\\tstart\\n")
"""
STUB_FAIL = "import sys\nsys.exit(1)\n"
STUB_FAIL_DP = """import sys
a = sys.argv
q = a[a.index("--qcat") + 1]
if "DN_vs_DP" in q:
    sys.exit(1)
out = a[a.index("--out-prefix") + 1]
open(out + ".stats.tsv", "w").write("chrom\\tstart\\n")
"""


def _setup(tmp_path, stub, comps=("DN_vs_DP", "DN_vs_EbKO"), n_perms=2):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    shutil.copy(SCRIPT, bin_dir / "run_perm_diff_pvalue.sh")
    (bin_dir / "bearing_pvalue.py").write_text(stub)
    cmp_dir = tmp_path / "compare"
    cmp_dir.mkdir()
    for c in comps:
        (cmp_dir / f"diff_{c}.qcat.bgz").write_text("x")
        for p in range(1, n_perms + 1):
            d = tmp_path / "perm" / f"perm{p}" / "diff_comparison"
            d.mkdir(parents=True, exist_ok=True)
            (d / f"diff_{c}.qcat.bgz").write_text("x")
    return bin_dir


def _run(tmp_path, bin_dir, env_extra=None, n_perms=2):
    env = dict(os.environ)
    env.pop("ALLOW_SKIP", None)
    env.pop("FORCE_PVALUE", None)
    # make "python" in the script resolve to this interpreter
    shim = tmp_path / "shim"
    shim.mkdir(exist_ok=True)
    if not (shim / "python").exists():
        (shim / "python").symlink_to(sys.executable)
    env["PATH"] = f"{shim}{os.pathsep}{env['PATH']}"
    if env_extra:
        env.update(env_extra)
    return subprocess.run(
        ["bash", str(bin_dir / "run_perm_diff_pvalue.sh"), str(tmp_path / "compare"),
         str(n_perms), str(tmp_path / "perm" / "perm"), str(tmp_path / "pvalue"),
         "kl", "", "l1"],
        env=env, capture_output=True, text=True)


@pytest.mark.skipif(shutil.which("bash") is None, reason="needs bash")
def test_all_written_exits_zero(tmp_path):
    r = _run(tmp_path, _setup(tmp_path, STUB_OK))
    assert r.returncode == 0, r.stderr
    assert "2 of 2" in r.stdout


@pytest.mark.skipif(shutil.which("bash") is None, reason="needs bash")
def test_all_failed_exits_nonzero(tmp_path):
    r = _run(tmp_path, _setup(tmp_path, STUB_FAIL))
    assert r.returncode == 1
    assert "no stats.tsv for: DN_vs_DP DN_vs_EbKO" in r.stderr


@pytest.mark.skipif(shutil.which("bash") is None, reason="needs bash")
def test_one_failed_exits_nonzero_and_allow_skip(tmp_path):
    b = _setup(tmp_path, STUB_FAIL_DP)
    r = _run(tmp_path, b)
    assert r.returncode == 1
    assert "1 of 2" in r.stdout and "DN_vs_DP" in r.stderr
    r2 = _run(tmp_path, b, {"ALLOW_SKIP": "1"})
    assert r2.returncode == 0, r2.stderr


@pytest.mark.skipif(shutil.which("bash") is None, reason="needs bash")
def test_missing_nulls_exits_nonzero(tmp_path):
    b = _setup(tmp_path, STUB_OK)
    shutil.rmtree(tmp_path / "perm")
    r = _run(tmp_path, b)
    assert r.returncode == 1
    assert "no null qcats found" in r.stderr


@pytest.mark.skipif(shutil.which("bash") is None, reason="needs bash")
def test_resume_counts_as_present(tmp_path):
    b = _setup(tmp_path, STUB_FAIL)
    (tmp_path / "pvalue").mkdir()
    for c in ("DN_vs_DP", "DN_vs_EbKO"):
        (tmp_path / "pvalue" / f"diff_{c}.stats.tsv").write_text("chrom\n")
    r = _run(tmp_path, b)
    assert r.returncode == 0, r.stderr
    assert "[resume]" in r.stdout
