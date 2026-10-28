# ----------------------------------------------------------------------
# File     : tests/test_fig7e_regional_empirical.py
# Version  : 1.0.0
# Date     : 2026-10-28
# Authors  : Katharina E. Hayer (katharinaehayer@gmail.com) and Claude
#            (Anthropic), co-created
# ----------------------------------------------------------------------
"""Fig 7E renderer: reads regional_empirical tables, writes pdf/png/source."""
import csv
import os
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(REPO, "paper", "figures", "fig7e_regional_empirical.py")
COLS = ["locus", "comparison", "region", "chrom", "start", "end", "tested_bins",
        "k", "k_pos", "k_neg", "n_locus", "p_spatial", "p_directional",
        "p_combined", "g_dir", "testable", "emp_n", "emp_count", "emp_p",
        "q_testable", "q_all"]


def _row(cmp, region, nb, emp, q):
    testable = nb >= 5
    return dict(locus="tcrb", comparison=cmp, region=region, chrom="chr6",
                start=0, end=1, tested_bins=nb, k=2, k_pos=2, k_neg=0, n_locus=4,
                p_spatial=0.1, p_directional=0.1, p_combined=0.1, g_dir=0.5,
                testable="yes" if testable else "no",
                emp_n=2000 if testable else "", emp_count=0 if testable else "",
                emp_p=emp if testable else "", q_testable=q if testable else "",
                q_all=q if testable else 1)


def test_fig7e_writes_outputs(tmp_path):
    tsv = tmp_path / "regional_empirical_tcrb.tsv"
    rows = [_row("DN_vs_A", "DJ_RC", 50, 0.0005, 0.002),
            _row("DN_vs_A", "Trbv1", 3, None, None),
            _row("DN_vs_B", "DJ_RC", 40, 0.3, 0.6),
            _row("DN_vs_B", "Trbv1", 8, 1.0, 1.0)]
    with open(tsv, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    prefix = str(tmp_path / "out" / "fig7e")
    subprocess.run([sys.executable, SCRIPT, "--tables", str(tsv),
                    "--comparisons", "DN_vs_A", "DN_vs_B",
                    "--out-prefix", prefix], check=True)
    for ext in (".pdf", ".png", ".source.tsv"):
        assert os.path.getsize(prefix + ext) > 0
    with open(prefix + ".source.tsv") as fh:
        src = list(csv.DictReader(fh, delimiter="\t"))
    assert len(src) == 4
    untest = [r for r in src if r["testable"] == "no"]
    assert len(untest) == 1 and untest[0]["region"] == "Trbv1"
    call = [r for r in src if r["q_testable"] and float(r["q_testable"]) < 0.05]
    assert [(r["comparison"], r["region"]) for r in call] == [("DN_vs_A", "DJ_RC")]
