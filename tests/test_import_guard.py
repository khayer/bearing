"""Import / entry-point guard.

Purpose
-------
This test does not check any scientific behaviour. It is a structural safety
net: it asserts that every canonical module still imports and that every
``bearing-*`` console entry point still resolves to a callable ``main``.

Why it exists
-------------
The repository is being cleaned up (scratch scripts moved to ``dev/``,
duplicate files removed). The danger in that reorg is moving or deleting a
file that another module imports at run time -- for example
``bearing_hic_plot_pval_overlay`` (imported by three Hi-C plotters) or
``score_provenance`` (imported by ``bigwig_to_qcat``). Those modules are never
named in the docs, so a name-based cleanup would flag them as orphans. If one
is moved, this test goes red immediately instead of the breakage surfacing
only when a reviewer runs that specific plot.

Maintenance contract
---------------------
``CANONICAL_MODULES`` is the explicit list of top-level scripts that must
always import. When a script is intentionally retired to ``dev/``, remove it
from this list in the same commit. Adding a new pipeline script? Add it here.
The list is the machine-checkable definition of "these are load-bearing."

The Hi-C scripts import heavy optional deps (hicstraw, cooler) lazily inside
functions, so importing the module itself succeeds with only the core install.
If that ever changes, wrap the relevant entries with importorskip rather than
deleting them.
"""

import importlib
from importlib.metadata import entry_points

import pytest

# --- Canonical top-level pipeline modules (Section A + B of the cleanup) ----
CANONICAL_MODULES = [
    # core scoring / comparison / statistics
    "bigwig_to_qcat",
    "compare_qcat",
    "bearing_pvalue",
    "generate_perm_nulls",
    "regional_enrichment",
    "consolidate_regional_enrichment",
    "shift_bigwig",
    "rebin_qcat",
    "bearing_calibration",
    "diffuse_architecture_top1pct",
    "build_cohort_quantile_reference",
    "compare_adaptive_vs_default",
    "compare_default_vs_qnorm",
    "crosslocus_fdr_counts",
    "significant_bins_summary",
    "cbe_point_query",
    "feature_track_changes",
    "suppfig_s1_jsd_decomposition",
    "plot_calibration_summary",
    # BES / Hi-C analyses
    "bes_hic_correlation",
    "bes_differential_loops",
    "bes_loop_anchor_enrichment",
    "bearing_tad_compare",
    "sanity_check_tad_boundaries",
    "bearing_region_decomposition",
    "bearing_decomposition_report",
    # plotting entry points + plotters
    "bearing_hic_combined_plot",
    "bearing_hic_plot",
    "bearing_hic_plot_triangle",
    "bearing_hic_kl_track_plot",
    "bearing_kl_track_plot",
    "batch_bearing_hic_plots",
    "batch_pygenometracks",
    # CI assertion helpers
    "assert_qcat_nonempty",
    "assert_score_provenance",
]

# --- Modules that LOOK orphaned by name but are imported at run time. -------
#     Each maps to the canonical modules that import it. Moving one of these
#     without updating its importers is exactly the failure this guards.
LOAD_BEARING_HELPERS = {
    "bearing_hic_plot_pval_overlay": [
        "bearing_hic_combined_plot",
        "bearing_hic_plot_triangle",
        "bearing_hic_plot",
    ],
    "score_provenance": [
        "bigwig_to_qcat",
        "assert_score_provenance",
    ],
}

# --- bearing/ package submodules -------------------------------------------
PACKAGE_MODULES = [
    "bearing",
    "bearing.diagnostics",
    "bearing.hic_io",
    "bearing.plot_layout",
    "bearing.plot_legend",
    "bearing.plot_loaders",
    "bearing.plot_tracks",
    "bearing.runner",
    "bearing.sheet",
    "bearing.track_primitives",
    "bearing.tsv",
    "bearing.utils",
    "bearing.validate",
]

# --- Console entry points declared in pyproject.toml [project.scripts] ------
EXPECTED_ENTRY_POINTS = {
    "bearing-score": "bigwig_to_qcat:main",
    "bearing-compare": "compare_qcat:main",
    "bearing-pvalue": "bearing_pvalue:main",
    "bearing-regional": "regional_enrichment:main",
    "bearing-perm": "generate_perm_nulls:main",
    "bearing-shift": "shift_bigwig:main",
    "bearing-rebin": "rebin_qcat:main",
}


@pytest.mark.parametrize("modname", CANONICAL_MODULES)
def test_canonical_module_imports(modname):
    """Every canonical top-level script must import cleanly."""
    importlib.import_module(modname)


@pytest.mark.parametrize("modname", PACKAGE_MODULES)
def test_package_module_imports(modname):
    """Every bearing/ package submodule must import cleanly."""
    importlib.import_module(modname)


@pytest.mark.parametrize("helper,importers", sorted(LOAD_BEARING_HELPERS.items()))
def test_load_bearing_helper_present(helper, importers):
    """Helper modules that are imported (not documented) must stay importable,
    and their importers must still reference them.

    If a cleanup moves `helper` out of the top-level import path, importing it
    here fails. If a cleanup removes the import from a dependent, the source
    check below fails. Either way the reorg is caught before it ships.
    """
    importlib.import_module(helper)
    for importer in importers:
        mod = importlib.import_module(importer)
        src = ""
        if getattr(mod, "__file__", None):
            with open(mod.__file__, "r", encoding="utf-8") as fh:
                src = fh.read()
        assert helper in src, (
            "%s no longer imports %s -- update LOAD_BEARING_HELPERS if this is "
            "intentional" % (importer, helper)
        )


def test_entry_points_registered():
    """Every bearing-* console script from pyproject is installed."""
    installed = {
        ep.name: ep.value
        for ep in entry_points(group="console_scripts")
        if ep.name.startswith("bearing-")
    }
    missing = set(EXPECTED_ENTRY_POINTS) - set(installed)
    assert not missing, "missing console entry points: %s" % sorted(missing)
    for name, target in EXPECTED_ENTRY_POINTS.items():
        assert installed[name] == target, (
            "%s points at %s, expected %s" % (name, installed[name], target)
        )


@pytest.mark.parametrize("name", sorted(EXPECTED_ENTRY_POINTS))
def test_entry_point_loads_callable(name):
    """Each entry point must resolve to a callable main()."""
    (ep,) = [
        e
        for e in entry_points(group="console_scripts")
        if e.name == name
    ]
    fn = ep.load()
    assert callable(fn), "%s did not resolve to a callable" % name
