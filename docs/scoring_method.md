# Scoring Method: KL Divergence for Multi-Signal Genomic Data

## TL;DR -- The Intuition in Plain Language

Imagine you are standing in a library. Most shelves are a mix of fiction,
science, and history -- that is the background, the "average shelf." Now you
walk past a shelf that is almost entirely neuroscience textbooks. That shelf
stands out. The more different it is from the average, the more it tells you
something specific is happening there.

Bearing does the same thing for the genome.

**Your tracks are the books.** Each BigWig file is one type of signal --
ATAC-seq, H3K27ac, CTCF, 4C, and so on. At every 200 bp bin across the
genome, Bearing reads the signal level from each track.

**Each bin becomes a shelf.** The signal values are converted into a
probability vector: "of all the signal present in this bin, what fraction
comes from each track?" A bin where 80% of the signal is ATAC and 10% is
CTCF looks very different from a bin where 60% is H3K27ac and 20% is RNAseq.

**The background is the average shelf.** Bearing computes Q -- the
genome-wide average signal composition. Most bins look roughly like Q.
A typical piece of chromatin is neither particularly open nor particularly
modified; it is just... average.

**The score measures how surprising a bin is.** The KL divergence score for
each bin asks: how much information would you need to update your expectation
if you saw this bin after only knowing the genome-wide average? A bin that
looks just like the average scores near zero. A bin with a dramatic spike in
CTCF and cohesin -- with very little of anything else -- scores high for those
two tracks, because that combination is rare and specific.

**Each track gets its own bar in the stacked visualization.** The score is
decomposed per track, so you can see not just that a bin is unusual, but
*which signal is making it unusual*. A high CTCF bar means CTCF is the unusual
thing here. A tall H3K27ac bar on top of an ATAC bar means active enhancer.
This is exactly like a sequence logo, where the height of each letter tells
you how informative that nucleotide position is, and which base dominates.

### A worked example

Suppose you have 5 tracks: ATAC, CTCF, H3K27ac, RNAseq+, RNAseq-.
The genome-wide average (Q) is roughly:

```bash
ATAC      15%
CTCF       8%
H3K27ac   30%
RNAseq+   25%
RNAseq-   22%
```

Now consider three bins:

**Bin 1 -- a CTCF insulator site:**

```bash
ATAC      12%   (similar to background -- no score)
CTCF      65%   (8x higher than background -- high score)
H3K27ac    8%   (lower than background -- clamped to 0)
RNAseq+    8%   (lower than background -- clamped to 0)
RNAseq-    7%   (lower than background -- clamped to 0)
```

Result: a tall red CTCF bar. The bin is dominated by CTCF, which is unusual.

**Bin 2 -- an active enhancer:**

```bash
ATAC      45%   (3x higher than background -- medium score)
CTCF       5%   (similar to background -- no score)
H3K27ac   45%   (1.5x higher than background -- medium score)
RNAseq+    3%   (lower -- clamped to 0)
RNAseq-    2%   (lower -- clamped to 0)
```

Result: a green H3K27ac bar on top of a teal ATAC bar. Co-occurrence of
accessibility and H3K27ac -- the classic enhancer signature.

**Bin 3 -- an unremarkable intergenic region:**

```bash
ATAC      14%   (similar to background)
CTCF       9%   (similar to background)
H3K27ac   28%   (similar to background)
RNAseq+   26%   (similar to background)
RNAseq-   23%   (similar to background)
```

Result: near-zero score for all tracks. Nothing interesting here.

### Why not just look at each track separately?

You could. But with 15 tracks, that means 15 separate browser tracks to scan
simultaneously, mentally asking "is CTCF high here? is H3K27ac also high? is
ATAC open?" for every locus. Bearing collapses that cognitive load into a
single stacked bar where the height tells you how unusual the locus is, and
the color tells you what is driving it. Think of it as a single dashboard
instead of 15 separate gauges.

### What does the differential track add?

When comparing two conditions (e.g. control vs. knockout), the differential
track shows `score_A - score_B` per bin per track. A positive bar means that
track is more active in condition A at this locus. A negative bar means it is
more active in condition B. The stacked bars still use the same color scheme,
so a red positive bar means CTCF gained signal in A, and a green negative bar
means H3K27ac was lost.

---

## Overview

This document describes the information-theoretic scoring approach used in
`bigwig_to_qcat.py` to convert continuous BigWig signal tracks into per-bin
saliency scores for visualization as epilogos stacked bar tracks. It explains
the mathematical foundations, how the method was adapted from the original
epilogos framework for discrete chromatin states to continuous signals, and
why this approach is preferable to common alternatives such as raw signal
visualization, fold enrichment, or peak calling.

---

## Background: What Are We Trying to Score?

Given multiple genomic signal tracks measured across the same genome (for
example, ATAC-seq, RNA-seq, CTCF ChIP-seq, and histone modifications), the
goal is to assign a score to every 200 bp genomic bin that answers the
question:

> "How unusual is the combination of signals in this bin, and which signal
> is driving that unusualness?"

This is the core idea behind epilogos (Quon, Reynolds et al., bioRxiv 2025):
to reduce a complex multi-track dataset into a single intuitive track that
captures the dominant biology at each genomic position. The approach is
explicitly inspired by sequence logos, which summarize the information content
of aligned DNA or protein sequences, and applies the same information-theoretic
framework to genomic annotations.

---

## The Scoring Method: S1 KL Divergence

### 1. Binning and Signal Extraction

The genome is tiled into non-overlapping 200 bp bins. For each bin b and each
of N signal tracks (states), the mean BigWig signal value is extracted:

```bash
x_{b,i}  for bin b, state i in {1, ..., N}
```

This produces a raw signal matrix of shape (num_bins x N).

### 2. Normalization to a Probability Distribution

Raw signal values are not directly comparable across tracks because each
assay has its own scale and dynamic range. To put all signals on a common
footing, the raw signal vector for each bin is normalized to a probability
distribution P:

```bash
Step 1: clip to zero   -- x_{b,i} = max(x_{b,i}, 0)
Step 2: add pseudocount -- x_{b,i} = x_{b,i} + epsilon   (epsilon = 1e-6)
Step 3: normalize      -- P_{b,i} = x_{b,i} / sum_j(x_{b,j})
```

After normalization, the values P_{b,1} ... P_{b,N} sum to 1 for every bin b.
The pseudocount prevents division by zero and log(0) errors in bins where all
signals are zero. Note that this introduces a subtle artifact in genuinely
zero-signal bins -- see section 5 (Low-Signal Masking) for how this is handled. The vector P_{b} can be interpreted as: "given that a signal
is present in this bin, what fraction of it comes from each track?"

### 3. Computing the Background Distribution Q

The background distribution Q represents the genome-wide average signal
composition -- in other words, what a "typical" bin looks like:

```latex
Q_i = (1 / num_bins) * sum_b( P_{b,i} )
```

Q is computed as the mean of all per-bin probability vectors across the entire
genome (or across the supplied regions in testing mode). It encodes how much
each track contributes on average. Tracks that are broadly active genome-wide
(such as RNAseq or H3K27ac) will have a higher Q than tracks with sparse,
localized signal (such as a 4C viewpoint).

### 4. Per-State KL Divergence Score

For each bin b and state i, the S1 epilogos score is the per-state
contribution to the Kullback-Leibler (KL) divergence between the observed
distribution P_{b} and the background Q:

```latex
score_{b,i} = P_{b,i} * log2( P_{b,i} / Q_i )
```

Negative values -- which arise when P_{b,i} < Q_i, meaning the state is
suppressed relative to background -- are clamped to zero:

```latex
score_{b,i} = max( 0,  P_{b,i} * log2( P_{b,i} / Q_i ) )
```

This is intentional and consistent with the epilogos convention: only the
positive (salient) contributions are visualized. The full KL divergence across
all states for a bin sums the positive contributions:

```latex
KL(P_b || Q) = sum_i [ P_{b,i} * log2( P_{b,i} / Q_i ) ]
```

High total KL divergence means the bin has an unusual signal composition
compared to the genome-wide average. Each state's individual score tells you
how much of that unusualness is attributable to that particular signal track.

### 5. Low-Signal Masking

Before the scores are written to the output file, bins where the total
pre-pseudocount signal across all tracks is below a minimum threshold are
explicitly zeroed out:

```latex
if sum_i( x_{b,i} ) < min_signal:
    score_{b,i} = 0  for all i
```

The default threshold is 0.01 and can be adjusted with the --min-signal
argument (set to 0 to disable entirely).

This step is necessary because the pseudocount added in step 2 introduces
a subtle artifact in zero-signal bins. Without masking, a bin where every
BigWig track reads zero gets a uniform P = [1/N, ..., 1/N] after
normalization. For focal tracks with Q_i < 1/N (e.g. CTCF, 4C viewpoints),
the ratio P_i / Q_i is greater than 1, producing a small but nonzero positive
KL score -- even though there is no real signal. In repeat-masked or
unmappable regions this generates spurious scores that are purely an artifact
of the pseudocount, not biology.

The masking is applied on the same non-negative matrix used to build P,
before pseudocount addition. In default mode (no cross-track normalization),
this is the clipped raw signal (after absolute-value correction for configured
negative-strand tracks). If `--normalize-tracks` is enabled, masking is
applied after the chosen cross-track normalization step.

### 6. Sorting and Output

Within each bin, the (state, score) pairs are sorted in descending order of
score. The dominant state -- the one most responsible for the bin's unusualness
-- is placed first. This is what pyGenomeTracks uses to draw the stacked bar:
the tallest bar segment appears at the bottom, creating an intuitive visual
where the "lead signal" is always immediately visible.

---

## Why KL Divergence? Comparison With Alternative Approaches

### Raw Signal Visualization

The most naive approach is to display each BigWig track as a separate line or
bar track. This works well for a small number of tracks but becomes unworkable
when integrating many signals simultaneously. Comparing signal heights across
tracks requires mental arithmetic by the viewer, different assays have
incomparable scales, and there is no single summary that captures which
genomic positions are "interesting" across all tracks simultaneously.

KL scoring collapses N tracks into a single scored track while preserving
which state is dominant -- the best of both worlds.

### Fold Enrichment Over Input

ChIP-seq analysis commonly reports fold enrichment: the signal in the IP
sample divided by the signal in a matched input control. This is a ratio that
accounts for local GC bias and mappability. However, fold enrichment:

- Requires a matched input control, which is not available for all assay types
  (RNA-seq, 4C, ATAC-seq typically do not use it)
- Is computed per-track, not across tracks simultaneously
- Does not provide a way to compare or combine signals from different assays
- Has a scale that varies between experiments and is not directly interpretable
  as probability

The KL approach does not require a control track. Instead, the genome-wide
background Q plays the role of the "expected" signal, computed directly from
the data itself.

### Peak Calling (MACS2, SEACR, etc.)

Peak calling is the standard approach for identifying enriched regions in
ChIP-seq and ATAC-seq data. Tools such as MACS2 model the read distribution
using a Poisson or negative binomial null model and report genomic intervals
that exceed a significance threshold (Comparative analysis of ChIP-seq peak
callers, PubMed 33412758). SEACR uses the global signal distribution as a
threshold for CUT&RUN data (Meers et al., Epigenetics & Chromatin 2019).

Peak calling has several well-known limitations in the multi-track integration
context:

- It produces a binary result (peak or no peak) rather than a continuous score,
  discarding information about the magnitude and relative contribution of
  different signals
- Thresholds are set per-track, making it difficult to combine signals from
  assays with very different signal-to-noise ratios
- It is not designed to summarize which of multiple concurrent signals is
  dominant at a given position
- Peak callers are optimized for specific assay types and are not directly
  transferable across all data types

KL scoring is continuous (not binary), does not require threshold setting,
and naturally handles signals from multiple assay types simultaneously.

### Simple Normalization (CPM, RPKM, z-score)

Counts per million (CPM) or RPKM normalization corrects for sequencing depth
but does not address the question of which genomic bins are unusual. A z-score
over each track flags positions that deviate from the per-track mean, but:

- Z-scores are computed independently per track, so they do not capture
  co-enrichment patterns across tracks
- They are sensitive to the assumption of normality, which is often violated
  by highly skewed ChIP-seq or ATAC-seq signal distributions
- They do not provide a natural way to weight the contribution of each track
  proportionally

The KL approach is non-parametric: it makes no assumptions about the
distribution of signal values. It only requires that signals can be
normalized to a probability distribution, which is always possible given
non-negative values.

### Entropy-Based Methods

Shannon entropy H(P) = -sum_i P_i * log2(P_i) measures the uniformity of the
signal distribution at a bin. A bin where all signals are equal has maximum
entropy; a bin dominated by a single signal has low entropy. Entropy could in
principle be used as a saliency score, but it has a key disadvantage: it is
symmetric. A bin where only ATAC is active and a bin where only RNAseq is
active would get the same entropy score, even though they represent completely
different biology.

KL divergence is asymmetric: it measures divergence from the specific
background Q, not just internal uniformity. This means it captures whether a
particular state is unusually high relative to its own genome-wide expectation.
A bin with high ATAC and low RNAseq gets a large positive score for ATAC (if
ATAC is typically low genome-wide) and no contribution from RNAseq.

---

## Adaptation from Discrete to Continuous Signals

The original epilogos S1 metric (Quon, Reynolds et al., bioRxiv 2025) is
defined for discrete chromatin state calls: at each genomic bin, each
biosample is assigned exactly one chromatin state, and P_{b,i} is the
fraction of biosamples assigned to state i at bin b. The metric is therefore
designed for a matrix of discrete labels, not continuous values.

This pipeline adapts the same framework for continuous BigWig signals by
treating the normalized signal vector as P. The key properties that make this
valid are:

1. Normalization: after clipping and pseudocount addition, the signal vector
   is normalized to sum to 1, giving it the same mathematical properties as
   a probability distribution.

2. Background: the genome-wide mean Q is computed in the same way as for
   discrete states, and plays the same role as the expected frequency.

3. Interpretation: the per-state KL score still has the same meaning --
   "how much does state i contribute to the unusualness of this bin?" -- but
   now "state" refers to a signal track rather than a discrete chromatin label.

The main difference in interpretation is that in the original epilogos, P
directly represents the fraction of samples in each state. Here, P represents
the relative contribution of each assay to the total signal in the bin. Both
are legitimate probability distributions and the KL computation is identical.

The adaptation is not without limitations. Because Q is computed from the
specific set of samples provided, it reflects the composition of this
particular experiment rather than a reference genome-wide background from a
large cohort. In the original epilogos, Q is estimated from hundreds of
biosamples, which makes it more robust.

A more significant issue arises from the dynamic range imbalance between assay
types. When forming the joint probability vector P, the normalization step
divides each bin's raw signal vector by its total signal. In bins where a
broad mark (e.g. H3K27ac or RNAseq) is active, that track will claim a large
fraction of the normalized vector -- not because it is unusually high at that
specific locus, but simply because it is broadly active genome-wide. As a
consequence, Q_i for broad tracks becomes large, which means P_{b,i} / Q_i
stays close to 1 even at genuinely enriched loci, suppressing their KL score.
Focal tracks (e.g. CTCF, 4C viewpoints) have low Q_i and therefore score
strongly wherever they are active.

The result is a systematic bias: broad tracks are chronically underscored and
focal tracks are overscored relative to what a biologically neutral comparison
would suggest.

### Track Normalization: Correcting the Dynamic Range Imbalance

To address this, the pipeline offers a `--normalize-tracks` option that applies
cross-track quantile normalization before forming P.

Current implementation details in `bigwig_to_qcat.py`:

- Default method: `--normalize-method nonzero-quantile`
- Legacy method: `--normalize-method quantile`

Both methods are rank-based and preserve exact zeros in each track. The
difference is how the shared reference distribution is built.

Default `nonzero-quantile` procedure:

1. For each track, collect and sort only non-zero values.
2. Build a shared non-zero reference distribution by averaging matched
  quantiles across tracks.
3. Map each track's non-zero values to that reference by within-track rank;
  zeros stay zero.

Legacy `quantile` procedure:

1. Sort full columns (including zeros) for all tracks.
2. Build a shared reference as the mean at each rank.
3. Map non-zero values by within-track rank to that full reference; zeros stay
  zero.

Because zeros participate in the legacy reference, very weak non-zero values
can be mapped down to zero. The default `nonzero-quantile` method was added to
reduce that behavior and better preserve low non-zero signal.

After normalization, tracks are brought onto a comparable marginal scale. A
value at the 90th percentile of CTCF is made comparable to a value at the 90th
percentile of H3K27ac, reducing between-assay dynamic range imbalance before
the tracks compete in the joint probability vector.

The KL computation then proceeds identically. The key effect is on Q: because
all tracks now have the same marginal distribution, Q_i will be approximately
equal across states (each track claims roughly 1/N of the total signal on
average). This means the ratio P_{b,i} / Q_i genuinely reflects per-bin
enrichment relative to that track's own typical signal level, rather than
reflecting the genome-wide breadth of the assay.

**When to use `--normalize-tracks`:**
Use it when your set of BigWig files mixes assay types with substantially
different signal breadth -- for example, a broad histone mark alongside a
focal transcription factor or a 4C viewpoint. If all tracks are of similar
breadth (e.g. all ChIP-seq for histone marks), normalization has less impact.

**Limitation:** Cross-track quantile normalization assumes all tracks should be equally
"important" at a genome-wide level. If one track genuinely is more broadly
active (e.g. a permissive histone mark in a highly transcribed cell type),
normalization will suppress that global difference and only preserve relative
enrichment patterns. This is usually desirable for visualization purposes but
should be considered when interpreting absolute score magnitudes.

---

## The Three Epilogos Metrics: S1, S2, S3

The epilogos paper defines three saliency metrics of increasing complexity
(Quon, Reynolds et al., bioRxiv 2025):

**S1** is the standard per-state KL divergence described above. It is the
metric implemented in this pipeline. The authors note it provides a good
balance of computational efficiency and biological interpretability.

**S2** extends S1 by incorporating pairwise co-occurrence patterns between
chromatin states. It accounts for the fact that some states are
cell-type-specific (like enhancer states) while others are broadly active
(like promoter states). Observing multiple enhancer states co-occurring at
the same position is more surprising than observing multiple promoter states.
S2 captures this by computing a pairwise co-occurrence matrix and comparing
it to a genome-wide expected co-occurrence background.

**S3** further extends S2 by modeling between-biosample similarities,
accounting for the fact that some biosamples are more similar to each other
than to the population at large.

S2 and S3 are not implemented in this pipeline. They require a large number
of biosamples (ideally dozens or hundreds) and substantially more computation.
For single-experiment multi-track data of the kind this pipeline is designed
for, S1 is the appropriate choice.

---

## Biological Validation of the KL Approach

The original epilogos paper demonstrates that S1 saliency scores correlate
with multiple independent lines of biological evidence, including evolutionary
sequence conservation (a strong indicator of functional relevance), the
density of trait-associated genetic variants from GWAS studies, and regulatory
sequence content. In their analysis, the genome was divided into 100 equally
sized groups by S1 saliency score, and higher-scoring bins showed
progressively stronger enrichment for all three biological signals.

This gives confidence that the KL divergence framework captures genuine
biological signal and not just technical variation. The same reasoning applies
when adapting the framework to continuous BigWig signals: bins that receive
high scores should represent positions where the combination of assay signals
is genuinely unusual relative to the genome-wide background.

---

## Summary

| Property | Peak calling | Fold enrichment | Raw signal | KL (no norm.) | KL + normalization (`nonzero-quantile` default; `quantile` legacy) |
| ---------- | -------------- | ----------------- | ------------ | --------------- | -------------------------------------------------------------------- |
| Continuous score | No (binary) | Yes | Yes | Yes | Yes |
| Multi-track integration | No | No | No | Yes | Yes |
| Requires control track | Sometimes | Yes | No | No | No |
| Threshold-free | No | No | Yes | Yes | Yes |
| Per-state attribution | No | No | No | Yes | Yes |
| Accounts for background | Via null model | Via input | No | Yes (genome-wide Q) | Yes (genome-wide Q) |
| Assay-agnostic | No | No | Yes | Partial | Yes |
| Corrects dynamic range imbalance | No | No | No | No | Yes |
| Handles low-signal/repeat regions | Via peak threshold | N/A | No | Partial (pseudocount artifact) | Yes (min-signal mask) |

The KL divergence approach is uniquely suited to the task of integrating
multiple heterogeneous genomic signal tracks into a single interpretable
visualization, which is why it forms the basis of the epilogos framework and
this pipeline.

---

## References

Quon J, Reynolds A, Tripician N, Rynes E, Teodosiadis A, Kellis M, Meuleman W.
Epilogos: information-theoretic navigation of multi-tissue functional genomic
annotations. bioRxiv (2025).
<https://doi.org/10.1101/2025.06.18.660301>

Lopez-Delisle L, Rabbani L, Wolff J, Bhardwaj V, Backofen R, Gruning B,
Ramirez F, Manke T. pyGenomeTracks: reproducible plots for multivariate
genomic datasets. Bioinformatics (2021) 37(3):422-423.
<https://doi.org/10.1093/bioinformatics/btaa692>

Zhang Y, Liu T, Meyer CA, et al. Model-based analysis of ChIP-Seq (MACS2).
Genome Biology (2008) 9:R137.
<https://doi.org/10.1186/gb-2008-9-9-r137>

Meers MP, Tenenbaum D, Henikoff S. Peak calling by Sparse Enrichment Analysis
for CUT&RUN chromatin profiling. Epigenetics & Chromatin (2019) 12:42.
<https://doi.org/10.1186/s13072-019-0287-4>

Corces MR et al. Chromatin accessibility profiling methods. Nature Reviews
Methods Primers (2022).
<https://doi.org/10.1038/s43586-022-00182-6>

Ernst J, Kellis M. ChromHMM: automating chromatin-state discovery and
characterization. Nature Methods (2012) 9:215-216.
<https://doi.org/10.1038/nmeth.1906>
