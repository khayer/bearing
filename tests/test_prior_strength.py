# ----------------------------------------------------------------------
# File     : tests/test_prior_strength.py
# Version  : 1.1.0
# Date     : 2026-10-07 (1.0.0: 2026-09-28)
# Authors  : Katharina E. Hayer (katharinaehayer@gmail.com) and Claude
#            (Anthropic), co-created
# Purpose  : pin the composition prior (--prior-strength): off = identical
#            scoring and provenance digest; on = no hard-floor cliff and no
#            single-track ceiling at barely-detectable bins.
# ----------------------------------------------------------------------
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import bigwig_to_qcat as b2q  # noqa: E402
from score_provenance import score_provenance_signature  # noqa: E402

Q = np.array([0.30, 0.20, 0.068, 0.15, 0.12, 0.162])


def _score(x, **kw):
    x = np.asarray(x, dtype=float)
    P = b2q.signals_to_prob(x)
    s, _ = b2q.kl_scores_per_bin(P, Q, raw_signal_matrix=x, **kw)
    return s.sum(axis=1)


def test_prior_off_is_identical():
    x = np.array([[0, 0, 0.15, 0, 0, 0], [5, 3, 1, 8, 2, 4.0]])
    a = _score(x, min_signal=0.1)
    b = _score(x, min_signal=0.1, prior_strength=0.0, prior_source=x)
    assert np.array_equal(a, b)


def test_hard_floor_cliff_and_ceiling_exist_without_prior():
    above = _score([[0, 0, 0.15, 0, 0, 0]], min_signal=0.1)[0]
    below = _score([[0, 0, 0.09, 0, 0, 0]], min_signal=0.1)[0]
    assert below == 0.0
    assert abs(above - np.log2(1.0 / Q[2])) < 1e-3   # single-track ceiling


def test_prior_removes_cliff_and_keeps_strong_signal():
    x = np.array([[0, 0, 0.15, 0, 0, 0],
                  [0, 0, 0.09, 0, 0, 0],
                  [0, 0, 20.0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0.0]])
    s = _score(x, min_signal=0.0, prior_strength=1.0, prior_source=x)
    assert abs(s[0] - s[1]) < 0.2            # no cliff between 0.15 and 0.09
    assert s[2] > 0.9 * np.log2(1.0 / Q[2])  # strong single track keeps score
    assert s[3] == 0.0                       # empty bin: P = Q -> 0
    assert np.all(np.diff(_score(np.array([[0, 0, v, 0, 0, 0] for v in
                  (0.05, 0.1, 0.5, 2, 10)]), min_signal=0.0,
                  prior_strength=1.0, prior_source=np.array(
                  [[0, 0, v, 0, 0, 0] for v in (0.05, 0.1, 0.5, 2, 10)])))
                  > 0)                       # score grows smoothly with signal


def test_prior_requires_source():
    try:
        _score([[1, 1, 1, 1, 1, 1]], prior_strength=1.0)
    except ValueError:
        return
    raise AssertionError("expected ValueError without prior_source")


def test_provenance_digest_unchanged_when_prior_off():
    args = (False, None, "kl", 0.1, "cats.yaml", None, None)
    old = score_provenance_signature(*args)
    zero = score_provenance_signature(*args, prior_strength=0.0)
    on = score_provenance_signature(*args, prior_strength=1.0)
    assert old == zero
    assert on[0] != old[0] and "prior_strength=1" in on[1]


def test_consensus_rescore_matches_scorer_with_prior():
    # compare_qcat --consensus-q must apply the same prior as bigwig_to_qcat
    # (1.1.0: it previously ignored prior_strength and used min_signal 0.01).
    import compare_qcat as cq
    x = np.array([[0, 0, 0.15, 0, 0, 0], [5, 3, 1, 8, 2, 4.0],
                  [0.02, 0, 0, 0, 0, 0.03]])
    raw = {i: x[i] for i in range(len(x))}
    for alpha in (0.0, 1.0):
        got, _, _ = cq.rescore_bins_with_consensus_q(
            raw, Q, min_signal=0.1, prior_strength=alpha)
        want = _score(x, min_signal=0.1, prior_strength=alpha, prior_source=x)
        for i in range(len(x)):
            assert np.isclose(float(got[i].sum()), want[i], atol=1e-5)
