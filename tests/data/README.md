# Simulated test data

Tiny, hand-built/simulated fixtures used by `tests/test_core_characterization.py`.
Regenerate-free: they are small and committed as-is.

| File | Purpose |
|---|---|
| `sim_features.bed` | BED3–BED9 features on `chrSim` incl. comment/`track` lines, malformed rows, comma and packed-integer itemRgb, invalid RGB, and an off-chromosome feature. |
| `sim_samples.tsv` | Sample sheet with mixed-case headers, comment and blank lines, padded cells and empty `bwN` columns. |
| `sim_regions.tsv` | Named regions with comma-formatted coordinates and an embedded comment. |
| `sim_contacts.tsv` | Sparse upper-triangle Hi-C records (seed 7, 100 bp bins over `chrSim:0-1000`) plus one out-of-region record; fed to fake `hicstraw` objects. |
