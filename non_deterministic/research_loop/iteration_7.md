# Iteration 7: The Uncertainty Axis Through All 24 Layers

**Direction**: Trace the uncertainty signal from embedding through layer 23. Does it rotate, flip, or emerge suddenly?

---

## The Headline

The uncertainty signal exists at **every single layer** and achieves 100% classification accuracy from layer 0 onward. But it undergoes continuous **rotation through different subspaces** — each layer's uncertainty direction is nearly orthogonal to distant layers. The previously reported "sign flip" at layer 22-23 is purely a norm phenomenon, not a direction change.

---

## Two Metrics Tell Different Stories

### Norm-Fragility Correlation (L2 magnitude of layer output)
```
Layers 0-8:   oscillating, weak (r = -0.6 to +0.5)
Layers 9-21:  strongly NEGATIVE (r = -0.70 to -0.98, peak at L17: -0.976)
Layer 22:     FLIPS to +0.69
Layer 23:     POSITIVE +0.97
Final_norm:   flips BACK to -0.96
```

### PCA PC1-Fragility Correlation (direction of maximum variation)
```
ALL layers:   POSITIVE, ranging from +0.72 (embed) to +0.99 (L23)
No sign flip ever occurs.
```

### The Reconciliation

The uncertainty signal is a **direction** in 1024-space. What changes across layers is its **magnitude relative to the background**:

- **Layers 9-21**: Background representations are larger for confident prompts (distributed, high-norm). The uncertainty signal is present but small. Total norm is dominated by background → negative norm-frag correlation.
- **Layer 23**: The uncertainty component grows huge (1.6x) for fragile prompts. It now dominates the total norm → positive norm-frag correlation.
- **Final_norm**: Layer normalization divides by norm, compressing the large L23 norms → norm flips back negative. Direction preserved.

**The "correction layer" story is reframed**: Layer 23 doesn't create or flip the uncertainty signal. It **amplifies** it until it overwhelms the background.

## The Uncertainty Axis Rotates Through Layers

Cosine similarity of each layer's PC1 direction with layer 23's PC1:

| Layer | cos(PC1, L23) | Interpretation |
|-------|--------------|----------------|
| embed | -0.05 | Orthogonal |
| layer_5 | -0.02 | Orthogonal |
| layer_10 | +0.04 | Orthogonal |
| layer_15 | +0.18 | Slight alignment |
| layer_17 | +0.15 | Slight alignment |
| layer_20 | +0.20 | Slight alignment |
| layer_22 | **+0.38** | Moderate alignment |
| layer_23 | **+1.00** | Reference |
| final_norm | **+0.72** | Strong preservation |

**The uncertainty direction at layer 23 is largely ORTHOGONAL to the direction at layers 0-17.** The model represents uncertainty using different features at different depths. Even though PC1 always correlates with fragility, it points in a DIFFERENT direction at each layer.

Pairwise cosine similarities between consecutive groups:
- layer_5 ↔ layer_10: 0.48 (moderate rotation)
- layer_10 ↔ layer_15: 0.52 (moderate)
- layer_15 ↔ layer_17: **0.77** (stable)
- layer_17 ↔ layer_20: **0.85** (very stable — these layers share a subspace)
- layer_20 ↔ layer_22: 0.37 (sharp rotation begins)
- layer_22 ↔ layer_23: 0.38 (rotation continues)

The uncertainty representation goes through three phases:
1. **Layers 0-8**: Uncertainty encoded in rapidly rotating subspaces (low inter-layer cos)
2. **Layers 9-20**: A relatively stable subspace (cos 0.77-0.85 between neighboring layers)
3. **Layers 21-23**: Rapid rotation into the final subspace used for output

## Emergence: The Signal Is There From the Start

| Layer | |r(PC1, frag)| | Classification Accuracy | PC1 Variance Explained |
|-------|-----------------|------------------------|----------------------|
| embed | 0.77 | 96% | 43% |
| layer_0 | 0.86 | **100%** | 46% |
| layer_7 | 0.91 | 100% | 64% |
| layer_12 | 0.94 | 100% | 71% |
| layer_17 | **0.97** | 100% | **82%** |
| layer_21 | 0.98 | 100% | **86%** |
| layer_23 | **0.99** | 100% | 81% |

The uncertainty signal is already strong enough for perfect classification at **layer 0** (the very first transformer block). What subsequent layers do is:
1. **Increase |r|**: from 0.77 at embed to 0.99 at L23 (noise reduction)
2. **Increase variance explained**: PC1 captures more of the total variation
3. **Rotate** the encoding direction into the final output subspace

## Layer Residual Analysis

Each layer's CONTRIBUTION to the uncertainty signal (residual = layer_k - layer_{k-1}):

| Layer residual | PC1 r(frag) | PC1 VarExp | Interpretation |
|---------------|-------------|-----------|----------------|
| layer_7 | +0.94 | 59% | Strong uncertainty contributor |
| layer_12 | +0.94 | 67% | Strong contributor |
| **layer_16** | **+0.41** | 77% | **Weak!** Not uncertainty-related |
| **layer_17** | **+0.29** | 67% | **Very weak!** Doing something else |
| layer_18 | +0.94 | 55% | Strong contributor resumes |
| **layer_23** | **+0.99** | **85%** | **The dominant contributor** |

**Layers 16-17 are doing something OTHER than uncertainty encoding.** Their residuals contribute almost nothing to the fragility signal (r = +0.29-0.41). This is despite the fact that their OUTPUTS strongly correlate with fragility (r > 0.93). These layers pass through the accumulated uncertainty signal without modifying it — they're working on other aspects of language processing.

## Theoretical Picture

The transformer processes uncertainty like a relay race:
1. **Layer 0**: The uncertainty signal is created — already enough for 100% classification
2. **Layers 1-20**: Each layer rotates the signal into a new subspace while passing it forward. Most layers add to the signal (r > 0.9), but layers 16-17 pause to do non-uncertainty work
3. **Layers 21-23**: Rapid rotation and amplification into the final output subspace, where it dominates the norm

The signal is never lost, never flipped, and never created from scratch. It's a gradually refined encoding that rotates through the model's representational space while maintaining its information content.

## Scripts
- `research_loop/iter7_axis_across_layers.py`
