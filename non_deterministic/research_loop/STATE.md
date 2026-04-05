# Hardware-Induced Divergence: Research Loop State

## Iteration: 13
## Timestamp: 2026-04-05 ~05:00 PDT

---

## Last Iteration (13)

**Direction investigated**: Can fragility be predicted from prompt text alone (no inference)?

**Key findings**:
- **65% of fragility variance is predictable from 34 text features** (linear regression R-sq=0.645)
- **A single binary feature ("is this arithmetic?") explains 50%** (R-sq=0.503)
- Simple rules (category + type) get R-sq=0.578 — nearly as good as full regression
- `max_number_digits` is the second-best continuous predictor (r=+0.49): bigger numbers → more fragile
- Exp vs ctrl classification from text alone: 67.2% accuracy (vs 100% from activations)
- Quartile prediction: 91% within-1-quartile accuracy
- 35% of variance is MODEL-SPECIFIC: same prompt type but different generation paths

**Detailed report**: `research_loop/iteration_13.md`

---

## Established Findings (verified)

### Fragility Landscape
- 1999 prompts scored, fragility range 0.4%-62.5%, median 9.4%, mean 11.7%
- Category ranking: reasoning (18.3%) >> continuation (13.6%) > factual (9.3%) > conversational (8.8%) > open_ended (8.3%)
- Within reasoning: multiplication (33.7%) > addition (30.7%) > subtraction (26.8%) > percentage (22.0%) >> sequences (5.4%) > word problems (4.0%)
- ~~Prompt length vs fragility: r=-0.82 within reasoning~~ **80% CONFOUNDED (Iter 12)**: the correlation is driven by problem TYPE (arithmetic vs word problems), not length. Within-type, r drops to -0.17. Addition and subtraction show POSITIVE within-type correlations.

### Position Fragility (Iter 2)
- Position 0: 17.0% (newline), Position 1: 23.5% (newline), **Position 2: 82.9%** (peak — "The" vs "I")
- Smooth monotonic decline after position 2: plateau ~50% through pos 18, then ~28% at pos 31, ~7% at pos 64-127, ~3% at pos 128+
- Per-category at position 2: continuation 94.5% > open_ended 87.8% > reasoning 85.3% > factual 79.3% > conversational 68.0%

### The "The" vs "I" Contest (Iter 2)
- 42.4% of position-2 fragility: choosing between sentence-starting words (The/I/A/")
- Stylistic choice (3rd vs 1st person), not content — but changes entire response structure

### Simulated Cross-Platform Divergence (Iter 3, CALIBRATED Iter 11)

**Divergence rate vs noise level (with calibrated platform mapping):**

| ε (logit noise) | % diverging | Exp:Ctrl ratio | Realistic scenario |
|-----------------|-------------|----------------|--------------------|
| 1e-4 | 0.3% | — | Same-vendor GPUs, FP32 (lower bound) |
| **1e-3** | **4.4%** | **9.0x** | **NVIDIA vs NVIDIA, FP32** |
| **5e-3** | **19.7%** | **4.2x** | **NVIDIA vs Apple Silicon, FP32** |
| 1e-2 | 31.7% | 3.9x | BF16 cross-platform (lower bound) |
| 5e-2 | 74.5% | 2.3x | BF16 cross-platform (upper bound) |
| 1e-1 | 90.4% | 1.5x | Extreme / quantized models |

**Semantic severity of first divergence (all noise levels):**
- 82-90% content-changing (word or digit swap)
- 0% cosmetic — every flip changes meaning
- Structural (sentence-start) flips increase at higher noise: 3.4% → 11.2%

**First divergence position distribution (ε=1e-2):**
- Positions 8-31 account for 47.3% of first flips (first 1-2 sentences)
- Only 5.5% at positions 0-2, 2.4% at positions 128+

### Cascade Statistics at Threshold 1.0 (NEW — Iter 3)
- Total runs: 35,136. Mean length: 1.70. Max: 16.
- 66.4% isolated (single fragile position), vs 93.6% at threshold 0.1
- 54.4% of prompts have max run >= 5, 5.2% >= 10
- Geometric decay: ~50% fewer runs at each length increment

### "Unflippable" Prompts (Iter 3 + REVISED Iter 4)
- 70 control prompts remain stable even at ε=1e-1 (all 256 margins > 0.1)
- **All 70 generate repetitive text** — echo prompt or repeat phrases. This is the mechanism.
- Category: reasoning 60% unflippable, factual 47.5%, open_ended 35%, conversational 25%, continuation 7.5%
- Reasoning is bimodal: word problems/sequences = very confident, arithmetic = very fragile
- Unflippable start with "The" 64% of time; margin curve ramps monotonically from 1.9 → 8.4
- Position 2 does NOT distinguish groups (both ~83-87% fragile) — it's the post-pos-2 trajectory that matters

### Two Attractor States (Iter 4, REVISED Iter 10)
- **Confident repetition**: model echoes/rephrases prompt or a simple pattern → all margins high → unflippable
- **Uncertain computation**: model attempts arithmetic or creative branching → margins collapse → fragile
- ~~Bifurcation point: by position ~32, the trajectory is determined~~ **REFUTED (Iter 10)**
- **REVISED**: The separation is GRADUAL, not a bifurcation. Predictive power: r=0.08 at pos 32, r=0.59 at pos 64, r=0.77 at pos 96
- 75.7% of prompts have simple monotonic-decrease trajectories. Only 0.6% U-shaped.
- The early32 distribution is unimodal (no bimodal gap) — attractor states are a continuous spectrum of decline rates
- Most flippable prompts show a margin DIP at positions 16-31 on AVERAGE, but this doesn't reliably predict individual prompt trajectories

### Semantic Taxonomy of Fragile Positions (Iter 1, threshold 0.1)
- 83.7% meaning-changing, 3.0% cosmetic, 13.2% ambiguous
- word-vs-word dominates early (63%), digit-vs-digit dominates late (71%)
- Crossover at ~position 55 (driven by reasoning prompts)

### Digit Confusion Structure (Iter 1)
- 98.6% in reasoning prompts. Adjacent digit confusion: 49.6% (2.5x baseline)
- Top pair '1' vs '3' is a repetitive loop artifact (93 prompts)

### Antonym/Opposition at Fragile Positions (NEW — Iter 8)
- **5.47% of fragile positions are oppositions** (3,267/59,744)
- Dominated by 'a' vs 'the' (1,404 = 83% of curated antonyms) — grammatical, not meaning-inverting
- True meaning inversions: man/woman (100), black/white (60), large/small (41), he/she (22), +/- (1,354 in reasoning)
- Opposition rate is flat across margin sizes (4-7%) — NOT concentrated at smallest margins
- At eps=1e-2: 4.6% of first divergences are oppositions; ~1% are truly meaning-inverting
- Per-category: reasoning 9.4% (sign flips), factual 5.7%, open_ended 3.8%, continuation 3.5%, conversational 1.7%
- Most spectacular concrete cases: black<->white coat color, out<->in direction, man<->woman character

### Activation Architecture (from 50-prompt subset, REVISED Iter 7)
- **NORM correlation** oscillates: negative at layers 9-21 (peak L17: r=-0.98), positive at L22-23 (L23: r=+0.97)
- **PC1 correlation** is ALWAYS positive: from r=+0.77 (embed) to r=+0.99 (L23)
- The "sign flip" is a norm phenomenon, not a direction flip. Background magnitude changes, not the signal.
- Layer 23 AMPLIFIES an existing uncertainty signal until it dominates total norm

### Dimension 487 / BOS Sink
- Single dim creates 50x norm anomaly at position 0, layers 5-21. Resolved by layer 23.
- **Completely independent** of the uncertainty encoding — zero overlap in top-10 dims

### Uncertainty Encoding
- Uncertain: compress to PR=22 at layer 20, explode at layer 23
- Confident: maintain PR=114 throughout. Cohen's d up to -10.6.

### Layer-23 Uncertainty Axis (Iter 6)
- **A single PC explains 83.5% of variance** in layer-23 residuals across 50 prompts
- The uncertainty signal is essentially 1D: a direction in 1024-space
- **274/1024 dims** correlate with fragility at |r| > 0.9 (they all co-vary along the axis)
- Top loading dims: 16 (-0.38), 737 (-0.34), 876 (-0.23), 453 (-0.20)
- Single-dim median-split gives **100% classification accuracy** (exp vs ctrl)
- Could build a real-time fragility probe: hook layer 23 residual, project onto PC1
- Temporal dynamics: confident prompts trigger LARGE initial residuals then cool rapidly; fragile prompts sustain
- The uncertainty encoding uses ~25% of layer 23's capacity (274/1024 dims)
- Completely orthogonal to the BOS attention sink (dim 487: r=0.18, rank 936/1024)

### Uncertainty Axis Across Layers (NEW — Iter 7)
- **100% classification from layer 0** — the model "knows" its confidence from the first block
- PC1-frag correlation climbs: embed 0.77 → L7 0.91 → L17 0.97 → L23 0.99
- **The axis ROTATES through different subspaces**: cos(L5, L23) ~ 0, cos(L15, L23) ~ 0.18
- Three phases: layers 0-8 (rapid rotation), 9-20 (stable, cos ~0.8), 21-23 (rapid rotation to output)
- **Layers 16-17 are "non-uncertainty" layers**: residual PC1 r=0.29-0.41 (they pass the signal without modifying it)
- Layer 23 residual has the strongest fragility signal of any layer (r=+0.99)
- The signal is a relay: never lost, never created from scratch, just rotated and amplified

### Prompt-Based Fragility Prediction (NEW — Iter 13)
- **65% of fragility variance predictable from prompt text alone** (34 features, linear regression)
- `is_arithmetic` alone: R-sq=0.503 (50% — the single most important predictor)
- Category alone: R-sq=0.209. Simple rules: R-sq=0.578. Full regression: R-sq=0.645.
- max_number_digits: r=+0.49 (bigger numbers → more fragile)
- Quartile prediction: 91% within-1-quartile from text alone
- Exp/ctrl text-only classification: 67.2% accuracy (vs 100% from layer-23 activations)
- The remaining 35% is model-specific: depends on generation path, not just prompt properties

### Quality-Stability Tradeoff (NEW -- Iter 5)
- **No sweet spot. The tradeoff is monotonic for non-reasoning prompts (r=+0.74).**
- Fragile text = diverse text. Stable text = repetitive text. They are mathematically opposed.
- Per-category: fragile third has 2.5-8.4x better composite quality than stable third
- Reasoning is an exception (U-shaped) due to bimodal structure
- The "most stable" text is degenerate loops; the "best quality" fragile text is math quizzes
- This is inherent to greedy decoding: diversity needs small margins, stability needs large ones

## Unverified Claims (need double-checking)
- ~~The r=-0.82 length-fragility correlation~~ **RESOLVED Iter 12**: 80% confounded by problem type
- Layer 23 "amplification layer" interpretation (revised from "correction") — Pythia-specific or general?
- Dimension 487 — truly a single dimension or are there others?
- 100% classification accuracy on 50 prompts — could be overfit (only have activations for 50)
- ~~The "bifurcation by position 32" claim~~ **REFUTED in Iter 10**
- ~~100% repetition rate in unflippable~~ **RESOLVED Iter 12**: real but universal (98-99% of all model output is repetitive)

## Open Questions
1. **Downstream cascade simulation**: After a token flips at position P, how does the divergence propagate? We can't fully simulate without the model, but can we estimate from cascade structure?
2. **Antonym/contradiction pairs**: Systematic search for semantic opposites at fragile positions — the most consequential divergences
3. **Attention head analysis in layer 23**: Which heads drive the correction behavior?
4. **Cross-platform noise calibration**: What is the actual logit-space noise magnitude for CUDA vs MPS vs A100?
5. **Does position of first divergence correlate with activation patterns?** Can layer norms predict it?
6. ~~Can fragility be predicted from the prompt alone?~~ **ANSWERED Iter 13**: 65% predictable, 50% from `is_arithmetic` alone
7. **Does the quality-stability tradeoff hold for sampling-based decoding?** (We only tested greedy)
8. **Is the Pythia-410M degenerate behavior (repetitive loops) specific to this model size, or does it scale?**

## Suggested Directions for Next Iteration
1. **What are layers 16-17 doing?** They contribute almost nothing to the uncertainty signal. Are they working on syntax? Grammar?
2. **Fragility probe validation**: Run the model on new prompts, hook layer 23, project onto PC1, verify fragility prediction on held-out data
3. **Generate paper figures**: Actually produce the 8 proposed figures from the synthesis
4. **Updated synthesis**: Revise the iteration 9 synthesis to incorporate iter 10 (bifurcation refutation) and iter 11 (noise calibration)
5. **Verify the "a" vs "the" pattern**: 1,404 cases dominate opposition counts. Is the definite/indefinite choice random or driven by specific prompt features?

## Data Locations
- Logit traces: `non_deterministic/experiment_0/logit_traces/prompt_NNNN.pt`
  - Keys: input_ids, generated_ids, top5_values (n_steps, 5), top5_indices (n_steps, 5)
- Fragility scores: `non_deterministic/prompts/fragility_scores.csv`
- Prompts: `non_deterministic/prompts/candidates_2000.jsonl`
- Experimental/control sets: `non_deterministic/prompts/experimental_200.jsonl`, `control_200.jsonl`
- Exp 1A outputs: `non_deterministic/experiment_1/outputs/h_nvidia_run1.jsonl`
- Activations (50 prompts): `non_deterministic/experiment_1/activations/h_nvidia/prompt_NNNN.pt`
- Iteration 1 scripts: `research_loop/analyze_competing_tokens.py`, `verify_position16.py`, `deep_semantic_analysis.py`, `digit_analysis.py`, `position_shift_analysis.py`, `verify_1vs3.py`
- Iteration 2 scripts: `research_loop/iter2_position_analysis.py`, `iter2_pos2_deep.py`
- Iteration 3 scripts: `research_loop/iter3_simulated_divergence.py`
- Iteration 4 scripts: `research_loop/iter4_unflippable.py`
- Iteration 5 scripts: `research_loop/iter5_quality_stability.py`
- Iteration 6 scripts: `research_loop/iter6_layer23_decomp.py`
- Iteration 7 scripts: `research_loop/iter7_axis_across_layers.py`
- Iteration 8 scripts: `research_loop/iter8_antonym_mining.py`
- Iteration 9 scripts: `research_loop/iter9_synthesis.py`
- Iteration 10 scripts: `research_loop/iter10_bifurcation.py`
- Iteration 11: `research_loop/iteration_11.md` (literature synthesis, no code)
- Iteration 12 scripts: `research_loop/iter12_verify_claims.py`
- Iteration 13 scripts: `research_loop/iter13_prompt_prediction.py`
