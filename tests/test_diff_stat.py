# ----------------------------------------------------------------------
# File     : tests/test_diff_stat.py
# Version  : 1.0.0
# Date     : 2026-10-12
# Authors  : Katharina E. Hayer (katharinaehayer@gmail.com) and Claude
#            (Anthropic), co-created
# ----------------------------------------------------------------------
"""bearing_pvalue.py --diff-stat: the default (sum) is unchanged, l1 tests the
sum of |per-track differentials| on the observed bins, the floor and the null,
and a planted composition switch is found by l1 and not by sum. ASCII only."""

import csv
import gzip
import json
import os
import subprocess
import sys

import numpy as np
import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
bp = pytest.importorskip("bearing_pvalue")


def test_diff_statistic_values():
    per = {1: 1.2, 2: -1.2, 3: 0.1}
    assert bp.diff_statistic(sum(per.values()), per, "sum") == pytest.approx(0.1)
    assert bp.diff_statistic(sum(per.values()), per, "l1") == pytest.approx(2.5)
    with pytest.raises(ValueError):
        bp.diff_statistic(0.0, per, "max")


def _write(path, n, rng, switch_every=0):
    with gzip.open(path, "wt") as fh:
        for i in range(n):
            d = rng.normal(0, 0.5, 6)
            if switch_every and i % switch_every == 0:
                d[1] += 2.0
                d[3] -= 2.0          # one track up, another down: sum cancels
            pairs = sorted([[float("%.6g" % v), k + 1] for k, v in enumerate(d)],
                           key=lambda x: -abs(x[0]))
            fh.write("chr1\t%d\t%d\tid:%d,qcat:%s\n"
                     % (i * 200, i * 200 + 200, i + 1,
                        json.dumps(pairs, separators=(",", ":"))))


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    pytest.importorskip("pyBigWig")
    d = tmp_path_factory.mktemp("diffstat")
    rng = np.random.default_rng(3)
    _write(str(d / "obs.qcat.gz"), 4000, rng, switch_every=40)
    nulls = []
    for k in range(3):
        p = str(d / ("null%d.qcat.gz" % k))
        _write(p, 4000, rng)
        nulls.append(p)
    (d / "cs.txt").write_text("chr1\t1000000\n")
    out = {}
    for st in ("default", "sum", "l1"):
        cmd = [sys.executable, os.path.join(REPO, "bearing_pvalue.py"),
               "--qcat", str(d / "obs.qcat.gz"), "--null-qcat"] + nulls + [
               "--out-prefix", str(d / st), "--diff", "--chrom-sizes",
               str(d / "cs.txt"), "--null-min-bytes", "0"]
        if st != "default":
            cmd += ["--diff-stat", st]
        r = subprocess.run(cmd, capture_output=True, text=True)
        assert r.returncode == 0, r.stderr[-2000:]
        with open(str(d / (st + ".stats.tsv"))) as fh:
            out[st] = list(csv.DictReader(fh, delimiter="\t"))
    return d, out


def test_default_is_sum_and_unchanged(runs):
    d, out = runs
    a = open(str(d / "default.stats.tsv")).read()
    b = open(str(d / "sum.stats.tsv")).read()
    assert a == b
    for r in out["sum"]:
        assert float(r["bearing_score_tested"]) == pytest.approx(
            abs(float(r["bearing_score"])), abs=1e-6)


def test_l1_tested_column_and_floor(runs):
    _, out = runs
    for r in out["l1"]:
        kl = [float(v) for k, v in r.items() if k.startswith("kl_")]
        assert float(r["bearing_score_tested"]) == pytest.approx(
            sum(abs(x) for x in kl), abs=1e-5)
        assert float(r["bearing_score_tested"]) >= 0.5
        # direction stays the sign of the net differential
        net = float(r["bearing_score"])
        assert r["direction"] == ("+" if net > 0 else "-" if net < 0 else "0")


def test_switch_found_by_l1_not_sum(runs):
    _, out = runs

    def called_switch(rows):
        return sum(1 for r in rows if r["significant_fdr0.05"] == "1"
                   and (int(r["start"]) // 200) % 40 == 0)
    assert called_switch(out["sum"]) == 0
    assert called_switch(out["l1"]) >= 50          # of 100 planted switches
