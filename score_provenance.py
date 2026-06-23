"""Scoring-provenance signature for BEARING qcats (standard-library only).

Single source of truth shared by bigwig_to_qcat.py (which stamps the signature
of the settings it ACTUALLY used, next to each qcat) and assert_score_provenance.py
(which recomputes the EXPECTED signature from the active config and halts the
pipeline on any mismatch).

It exists because the scoring-determining settings (track normalization,
divergence method, min-signal floor, category panel, adaptive bins, cohort
reference) are passed to the scorer as parameters, not as tracked file inputs.
Under Snakemake rerun-triggers=mtime a change to one of those parameters does
NOT re-trigger scoring, so a populated regime tree can silently reuse qcats that
were scored under DIFFERENT settings and produce false results. Stamping the
signature at score time and verifying it before the p-value layer catches that.

The signature is intentionally path-form robust: file fields contribute only
their basename (so absolute-vs-relative differences do not cause false
mismatches), scalar flags contribute their value. It captures the axes that
change the per-bin score. It does NOT detect a same-path different-content edit
of one of those files -- it guards the identity of the FLAGS and inputs, which
is the regime-confusion failure mode it is built for.
"""

import hashlib
import os


def _base(path):
    """Basename of a path field, or the literal 'none' when unset/empty."""
    if not path:
        return "none"
    return os.path.basename(str(path))


def score_provenance_signature(normalize_tracks, normalize_method, score_method,
                               min_signal, categories_path, bins_bed_path,
                               cohort_reference_path=None):
    """Return (digest, payload) for the scoring-determining settings.

    digest  -- sha256 hex of the canonical payload (written as the .sig first line)
    payload -- the human-readable field listing (written as comment lines)

    normalize_method is forced to 'none' when normalize_tracks is False, so a
    run that does not normalize has a canonical signature regardless of the
    method default that argparse may have filled in.
    """
    norm_on = bool(normalize_tracks)
    fields = [
        "norm_tracks=%d" % (1 if norm_on else 0),
        "norm_method=%s" % ((normalize_method or "none") if norm_on else "none"),
        "score_method=%s" % (score_method or "kl"),
        "min_signal=%s" % ("%.6g" % float(min_signal)),
        "categories=%s" % _base(categories_path),
        "bins_bed=%s" % _base(bins_bed_path),
        "cohort_ref=%s" % _base(cohort_reference_path),
    ]
    payload = "\n".join(fields)
    digest = hashlib.sha256(payload.encode("ascii")).hexdigest()
    return digest, payload
