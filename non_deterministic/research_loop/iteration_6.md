# Iteration 6: Layer 23 Residual Decomposition — The Uncertainty Axis

**Direction**: Decompose layer 23's residual contribution (layer_23 - layer_22) to find which dimensions encode uncertainty.

---

## The Headline Finding

Layer 23's uncertainty encoding is a **single principal axis that explains 83.5% of variance** across all 50 prompts. A single dimension (dim 609, r=+0.9956) can perfectly classify all 25 experimental and 25 control prompts. The signal is so strong that **274 out of 1024 dimensions** (26.8%) individually correlate with fragility at |r| > 0.9.

The BOS attention sink (dim 487) and the uncertainty encoding are **completely independent circuits** — zero overlap in their top-10 active dimensions.

---

## Part 1: The Uncertainty Subspace Is Actually 1D

PCA of the layer-23 residual (averaged over generated positions) across all 50 prompts:

| PCs | Variance Explained |
|-----|-------------------|
| 1 | **83.5%** |
| 2 | 92.1% |
| 3 | 95.7% |
| 5 | 98.3% |
| 10 | 99.6% |

PC1 separates experimental from control with Cohen's d = **-4.41**. The uncertainty signal is essentially one-dimensional: a single direction in 1024-space that varies smoothly from "confident" to "uncertain."

## Part 2: Per-Dimension Correlations — Pervasive But Low-Rank

While the signal lives on one axis, it's distributed across many dimensions:

- **274 dims** (26.8%) have |r| > 0.9 with fragility
- **596 dims** (58.2%) have |r| > 0.7
- **761 dims** (74.3%) have |r| > 0.5

Top uncertainty dimensions (more active when fragile):

| Dim | r with fragility | Exp mean | Ctrl mean | Ratio |
|-----|-----------------|----------|-----------|-------|
| 609 | **+0.9956** | 1.021 | 0.333 | 3.1x |
| 791 | +0.9904 | 2.013 | 1.056 | 1.9x |
| 455 | +0.9904 | 2.833 | 1.341 | 2.1x |
| 190 | +0.9901 | 1.750 | 0.649 | 2.7x |
| 218 | +0.9901 | 1.763 | 0.595 | 3.0x |
| 737 | +0.9897 | 9.651 | 3.248 | **3.0x** |

Top confidence dimensions (more active when confident):

| Dim | r with fragility | Exp mean | Ctrl mean | Ratio |
|-----|-----------------|----------|-----------|-------|
| 751 | -0.9736 | 1.123 | 2.157 | 0.52x |
| 576 | -0.9712 | 0.301 | 0.439 | 0.68x |
| 898 | -0.9687 | 0.471 | 0.933 | 0.50x |
| 788 | -0.9677 | 1.387 | 2.713 | 0.51x |

The interpretation: most dimensions participate in the uncertainty encoding, but they all move together along a single axis. This is why PCA shows 83.5% on PC1 despite the signal being distributed across 274+ dims — they're all correlated.

## Part 3: Perfect Classification From One Dimension

**100% classification accuracy (25/25 exp, 25/25 ctrl) using any single one of the top correlated dimensions.** Even dim 609 alone, with a simple median-split, perfectly separates the groups. This holds for top-1, top-3, top-5, top-10, top-20, and top-50 feature sets.

**Practical implication**: You could build a real-time "fragility probe" by:
1. Hooking layer 23's residual (layer_23 - layer_22)
2. Projecting onto the PC1 direction (or just reading dim 609)
3. The scalar value predicts the prompt's fragility with near-perfect accuracy

## Part 4: Residual Dynamics Over Generation

The temporal pattern of layer-23 residual norms reveals a **cooling effect**:

**Most fragile experimental (frag=0.625):**
- Start: 87.3 → Q1: 64.9 → Mid: 68.8 → End: 65.7
- Starts high, stays high — layer 23 works hard throughout

**Most confident control (frag=0.004):**
- Start: 96.9 → Q1: 43.1 → Mid: 37.2 → End: 33.6
- Starts even HIGHER than fragile, then collapses to half

Confident prompts actually trigger LARGER layer-23 residuals at the very start of generation (positions 0-7), but this rapidly drops as the model enters its confident repetitive loop. Fragile prompts sustain high residuals because the model never achieves confidence.

## Part 5: PC1 Loadings — The Uncertainty Direction

The uncertainty axis (PC1) is dominated by:
- **Dim 16**: loading = -0.380 (the strongest contributor)
- **Dim 737**: loading = -0.338
- **Dim 876**: loading = -0.227
- **Dim 453**: loading = -0.197
- **Dim 464**: loading = -0.151

(Negative loading + negative PC1 score for experimental = these dims are LARGER for uncertain prompts.)

These are NOT the same dimensions as the BOS attention sink (dim 487). The BOS sink and the uncertainty encoding are completely independent:
- Dim 487 fragility correlation: r = +0.18 (rank 936/1024 — nearly uncorrelated)
- Dim 487 PC1 loading: -0.004 (negligible)
- Top-10 overlap between BOS dims and uncertainty dims: **ZERO**

## Part 6: Concentration vs Distribution

The difference between experimental and control mean residuals is distributed:
- Top 1 dim: 1.9% of total |diff|
- Top 10 dims: 9.2%
- Top 50 dims: 24.1%
- Top 100 dims: 36.9%

This seems to contradict the PCA finding (83.5% on PC1). The resolution: the MEAN difference is spread across dimensions, but the VARIANCE structure is low-rank. All 274 high-correlation dimensions co-vary — when one goes up, they all go up — creating a single effective degree of freedom.

## Theoretical Interpretation

Layer 23's uncertainty encoding works like this:
1. **A single direction** in 1024-space encodes the fragility spectrum
2. **Hundreds of dimensions** participate in this direction (it's not one "uncertainty neuron")
3. The direction is **orthogonal** to the BOS attention sink — separate circuits
4. The signal is **temporally sustained** for fragile prompts but **decays** for confident ones
5. The decay time constant is the bifurcation: confident prompts "cool" layer 23's residual within ~64 tokens, while fragile prompts never cool

This suggests that Pythia-410M allocates a substantial portion of its final layer's capacity to a single purpose: encoding confidence/uncertainty. The fact that 274/1024 dims correlate at |r| > 0.9 means roughly 25% of layer 23's representational capacity is dedicated to this signal.

## Scripts
- `research_loop/iter6_layer23_decomp.py`
