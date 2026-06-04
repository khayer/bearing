# BEARING Statistics: Significance, Null Model, and Differential Testing

This document explains the statistical layer of BEARING: what the per-bin
p-value means, how the permutation null is built, what the differential
("diff") score is, and what the differential p-value tests. It is the companion
to [`scoring_method.md`](scoring_method.md), which covers how the per-bin score
itself is computed. Read that one first if you want the score mechanics; read
this one for what the numbers let you *conclude*.

---

## TL;DR -- The Intuition in Plain Language

At every 200 bp window of the genome, BEARING looks at the *mix* of the six
signals (ATAC, RNAseq+, RNAseq-, CTCF, cohesin/Rad21, H3K27ac) -- how much of
the activity here is openness, transcription, CTCF, and so on. Most windows
have an unremarkable, near-average mix. The **BEARING score** flags windows
whose mix stands out from the genome-wide average: places where the assays
pile up together in an unusual way.

But a window could look unusual just by luck. So we play a shuffling game: we
slide each track around its chromosome by a large random offset, which keeps
each signal looking like itself (same peakiness, same coverage) but scrambles
where the different tracks line up relative to each other. We re-score many of
these shuffled genomes. If a real window scores higher than almost any shuffled
window ever does, that is our evidence that the marks are genuinely
co-localized there, not lined up by accident. That evidence is the **p-value**.

The **diff** is the same idea applied to two conditions: at each window, which
marks went up or down going from, say, DN to DP. The **diff p-value** asks
whether that change is bigger than the shuffling noise -- a real, directional
rearrangement of chromatin rather than a coincidence.

One-line analogy: the score finds the genome's unusually busy intersections;
the shuffle test confirms the traffic is really coordinated rather than cars
happening to pass at once; the diff shows which intersections got busier or
quieter between two conditions.

---

## What the BEARING score is (and is not)

For each bin, the six per-track signals are normalized into a **composition**
P (the fraction of this bin's total activity contributed by each assay). The
genome-wide average composition is the **background** Q. The score sums the
per-track epilogos S1 / KL contributions `P_i * log2(P_i / Q_i)` (negatives
clamped to zero). A bin scores high when one or more assays are enriched
*relative to what is typical for that assay in this sample*.

- **It is:** a per-bin, per-assay measure of compositional distinctiveness --
  how unusual this bin's mix of chromatin marks is compared to the rest of the
  sample. It is reference-model-free: Q is the sample's own genome-wide
  average, not an external annotated reference or a trained model.
- **It is not:** an absolute signal level (a bin can have high raw signal but a
  typical composition and score low), a peak call, or a between-condition
  comparison (that is the diff). On its own it carries no significance -- a
  high score could in principle arise by chance, which is what the p-value is
  for.

---

## The statistical question

Two questions, one per analysis mode:

- **Single-sample:** *Where in the genome is a sample's combination of
  chromatin marks more spatially co-localized than chance would produce?*
- **Differential:** *Where does that structure change between two conditions,
  and in which direction?*

---

## The null hypothesis and the permutation null

The null hypothesis is **spatial independence of the assays**: the tracks have
no real positional relationship to one another, so any apparent co-localization
is coincidental.

We realize this null by **circular-shift permutation**. For each input BigWig,
apply a large independent random circular shift along the chromosome, then
re-score. Shifting preserves each track's own statistical character -- its
peakiness, coverage, and autocorrelation are unchanged -- but destroys the
*alignment between tracks*: CTCF peaks no longer sit where they did relative to
ATAC, and so on. The shifted genome is therefore a world in which each assay
still looks normal but their co-occurrence is random. This is exactly the null
we want: it isolates cross-track co-localization as the only thing being
tested.

Each permutation round (`n_perms` in the config) produces a full shuffled,
re-scored genome. Pooling all shuffled bins gives the **null score
distribution**.

Why permutation rather than a parametric model: a 2-component Gamma mixture was
evaluated and **fails on the flat score distributions** seen in broadly active
cell types (e.g. DN thymocytes, S3T3 fibroblasts), where there is no clean
low-score background peak to fit. The empirical permutation null makes no
assumption about the distribution's shape and is valid for any assay
combination, any cell type, and any score distribution.

---

## What the p-value means

The empirical one-sided p-value for a bin is the fraction of null bins scoring
at least as high as the observed bin:

    p = P(null score >= observed score)

Small p means a score this high essentially does not occur when the assays are
spatially decoupled -- so the real bin reflects genuine co-localized chromatin
structure, not coincidence. P-values are corrected across all bins with
Benjamini-Hochberg FDR (default q = 0.05).

Two practical notes:

- The minimum achievable p-value is `1 / (n_null_bins + 1)`, so more
  permutations buy finer resolution at the significant tail.
- The output is a BigWig of `-log10(p)` per bin plus a per-bin statistics TSV,
  so significance can be browsed as a genome track.

---

## What the diff is

The single-sample score says where one sample's composition is distinctive.
The **diff** compares two conditions, computed per bin per assay on the
replicate-averaged scores:

    diff_i = (mean KL in condition A) - (mean KL in condition B)

It is **signed**: positive means assay *i* is more enriched in A at this bin,
negative means more enriched in B. So `diff_DN_vs_DP` shows which marks shifted,
in which direction, between the two stages.

- **It is:** a signed, per-assay change in compositional enrichment between two
  conditions.
- **It is not:** a fold-change in raw signal. And on its own it is only a
  magnitude -- you do not yet know whether a given difference exceeds the noise
  between replicates and shuffles.

---

## What the diff p-value means

Same logic as the single-sample p-value, applied to the differences and
testing **both directions** via the absolute value:

    p = P(|null diff| >= |observed diff|)

The null here is the differential of the *shifted* data: diff qcats built from
the circularly-shifted tracks, pooled across permutations. A significant diff
p-value means the change between the two conditions at this bin is larger than
what spatially-decoupled assays produce -- a real, directional rearrangement of
chromatin rather than coincidental assay overlap.

The differential output BigWig is **signed**: positive `-log10(p)` where A > B,
negative where B > A, so significance and direction are visible in one track.
FDR-corrected at q = 0.05.

---

## Summary table

| Quantity            | Question                                                      | Null hypothesis                                            | Output                                  |
|---------------------|---------------------------------------------------------------|------------------------------------------------------------|-----------------------------------------|
| BEARING score       | How distinctive is this bin's chromatin composition?          | (descriptive; no test on its own)                          | per-bin score, qcat track               |
| Single-sample p     | Is the score higher than chance co-localization would give?   | The assays are spatially independent (alignment is random) | `-log10(p)` BigWig, stats TSV, FDR 0.05 |
| Diff                | How does the composition change between two conditions?        | (descriptive; signed magnitude)                            | signed per-assay diff qcat              |
| Diff p              | Is the change larger than shuffling noise, and which way?      | The between-condition difference is no larger than under spatial independence | signed `-log10(p)` BigWig, FDR 0.05     |

---

## A note on wording for write-ups

What the permutation null literally tests is **spatial co-localization /
alignment of the tracks** -- shifting destroys cross-track positional
relationships. It does not test mechanistic or causal coordination. The
defensible phrasing is "non-random spatial co-occurrence of the assays." The
biological motivation (recombination centers, enhancer-promoter modules) can
be stated separately as *why* such co-occurrence is interesting, without
claiming the test itself proves mechanism.

---

## See also

- [`scoring_method.md`](scoring_method.md) -- how the per-bin score is computed
  (KL mechanics, normalization, alternatives).
- `bearing_pvalue.py` -- p-value computation (empirical permutation, with
  Gamma fallbacks documented but not recommended).
- `generate_perm_nulls.py` -- builds the circular-shift permutation nulls
  (observed and differential).
