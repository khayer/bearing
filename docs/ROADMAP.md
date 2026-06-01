# Bearing — Publication Roadmap
Last updated: 2026-03-07

> **⚠️ NOTE**: This roadmap has been split into focused TODO files:  
> - **[code_todo.md](code_todo.md)** — Software architecture, packaging, testing  
> - **[manuscript_todo.md](manuscript_todo.md)** — Literature review, experiments, writing  
> 
> This file is kept for historical reference.

## Recent visual update

- 2026-05-02: Differential p-value tracks in triangle plots now support a
  low-alpha filled background plus per-bin lollipop markers colored by the
  dominant differential track, with nearest-neighbor alignment for slightly
  offset p-value and diff-score bins.

---

## 1. Software architecture (pre-publication requirements)

### 1.1 YAML-based category configuration
- Replace the hardcoded ALL_CATEGORIES list in bigwig_to_qcat.py with a
  user-supplied YAML config file
- Schema should define per-track: name, color, strand (pos/neg/abs), optional
  normalization group, optional display label
- Provide a default categories.yaml for the current 15-state mm10 setup so
  existing users are not broken
- CLI flag: --categories categories.yaml
- Validation: check that YAML keys match the number of BigWig files supplied,
  emit clear errors for missing fields
- This is the single biggest blocker for community adoption -- people need to
  plug in their own assay combinations without editing source code

### 1.2 Modular pipeline structure
- Split bigwig_to_qcat.py into importable modules:
    bearing/
      __init__.py
      io.py          -- BigWig reading, chrom sizes, BED parsing
      score.py       -- signals_to_prob(), kl_scores_per_bin(), quantile_normalize
      compare.py     -- JSD, Spearman, PCA, differential qcat (from compare_qcat.py)
      viz.py         -- write_ini(), write_debug_bigwigs()
      config.py      -- YAML loading, category validation
      cli.py         -- argparse entry points
- Expose bearing-score and bearing-compare as console_scripts in setup.cfg
- This allows people to import bearing.score in their own notebooks

### 1.3 Packaging
- pyproject.toml / setup.cfg with pinned dependency ranges
- PyPI package: pip install bearing-epigenomics (or similar, check name)
- Conda recipe for bioconda
- Docker / Singularity container for HPC use

### 1.4 Testing
- Synthetic BigWig fixtures (small, deterministic) for unit tests
- pytest suite covering: score correctness, YAML validation, qcat round-trip,
  JSD symmetry, diff qcat sign convention
- CI via GitHub Actions on push

---

## 2. Literature review (TODO)

### 2.1 Multi-omics integration methods -- comprehensive search needed
- Key question: who has done genomic bin-level integration of >3 continuous
  signal tracks without discretization?
- Search terms to cover:
    - "multi-omics integration" + "epigenomics" + "continuous signal"
    - "chromatin accessibility" + "histone modification" + "joint scoring"
    - "information theory" + "epigenomics" + "KL divergence" + "genomic bins"
    - MOFA, MEFISTO, Cobolt, totalVI (single-cell side -- note differences)
    - cisTopic, ArchR, EpiScanpy (single-cell ATAC -- different modality)
    - IDEAS, RFECS, Segway (continuous HMM approaches -- closest relatives)
    - ChromDiff, CSREP, epilogos diff (differential state -- already reviewed)

### 2.2 Specific lab / preprint to check
- kay nieselt lab, University of Tuebingen, Germany
  - Look for recent preprints on chromatin / epigenome integration
  - Check bioRxiv author search: nieselt + Tuebingen
  - May have overlap with the multi-assay scoring concept or the biological
    application -- confirm before finalizing the positioning statement
  - NOTE: biological application details to be provided separately
    (see Section 4 below)

### 2.3 4C-seq integration literature
- 4C-seq scoring in the context of multi-assay frameworks is unusual --
  document why we include it and who else has done similar
- Relevant: virtual 4C from Hi-C, NG Capture-C, Tiled-C
- Key question: has anyone scored 4C contact strength jointly with chromatin
  accessibility and histone marks at the bin level?

### 2.4 Hi-C / 3D genome integration
- Arc visualization in genome browsers: already supported by pyGenomeTracks
- Scoring angle: weighting KL scores by loop anchor contact frequency
  (no one has done this -- potential methods novelty)
- Search: "Hi-C" + "chromatin state" + "loop anchor" + "scoring" 2023-2026
- Relevant tools: Mustache, HiCCUPS, FitHiC2 for loop calling
  (Bearing would consume their output BED files)

---

## 3. GitHub repository structure

    bearing/
      README.md
      LICENSE
      pyproject.toml
      categories/
        mm10_15state.yaml       -- default config (current hardcoded categories)
        mm10_5state_minimal.yaml -- minimal example for new users
        hg38_15state.yaml       -- human equivalent (to be added)
      bearing/
        __init__.py
        io.py
        score.py
        compare.py
        viz.py
        config.py
        cli.py
      tests/
        fixtures/               -- small synthetic BigWig files
        test_score.py
        test_io.py
        test_compare.py
      docs/
        method.md               -- scoring_method.md (existing, move here)
        abstract.md             -- manuscript abstract draft (existing)
      use_cases/                -- see Section 4
        README.md
        mm10_enhancer_study/
          ...

---

## 4. Use case / biological application (vignette)

### 4.1 Overview
- A self-contained worked example in use_cases/mm10_enhancer_study/
- Someone should be able to clone the repo, download the example data, run
  three commands, and reproduce the key figure
- Format: step-by-step README.md + shell script + expected output screenshots

### 4.2 Biological context
- To be filled in with details from the lab
- Involves mouse mm10 data, DJ1/4C viewpoints, enhancer activity
- Likely centers on a specific locus where multi-assay Bearing score reveals
  coordinated regulatory activity not visible from individual tracks


### 4.3 What the vignette should demonstrate
- Running bearing-score on the full assay set
- Using --normalize-tracks and --debug-tracks to inspect intermediate signal
- Running bearing-compare on replicates + conditions
- Adding Hi-C arcs to the INI and rendering a final figure
- Interpreting the differential qcat output between two conditions

---

## 5. Manuscript plan

### 5.1 Target
- Primary: Nature Methods or Genome Research
- Fallback: Bioinformatics, NAR (if scoped as software)

### 5.2 Key figures (draft plan)
1. Method schematic: BigWigs -> probability landscape -> KL scoring -> qcat
2. Genome browser screenshot: Bearing track + individual BigWigs + Hi-C arcs
   at the application locus
3. Multi-sample comparison panel: JSD heatmap + PCA showing replicate
   clustering and condition separation
4. Differential epilogos track: condition A vs B at key loci, with state
   breakdown
5. Benchmarking: comparison to ChromHMM + epilogos on same data --
   show what Bearing recovers that discretization loses

### 5.3 Remaining experiments needed
- Benchmarking against discrete ChromHMM pipeline on same mm10 data
- Runtime and memory profiling (genome-wide, 15 tracks, typical HPC node)
- At least one additional organism or dataset to show generalizability
- Replicate QC validation: show JSD / Spearman metrics predict concordance

### 5.4 Author list and contributions
- To be decided
- Bassing, Allyn, glynn, Wu, Lee, Vincent, Clement, Oltz, Sacan , (me Hayer)

---

## 6. Immediate next steps (prioritized)

1. YAML category config (Section 1.1) -- unblock community use
2. Biological application details from lab (Section 4.2) -- needed for
   abstract anchor and use case vignette
3. nieselt lab literature search (Section 2.2) -- positioning
4. Modularize codebase into bearing/ package (Section 1.2)
5. Write use case vignette (Section 4.3)
6. Benchmarking figure against ChromHMM (Section 5.3)
