# Iteration 10: The Bifurcation Claim Is Wrong

**Direction**: Verify the claim "by position ~32, the trajectory is determined" (Iteration 4).

**Verdict: REFUTED.** Early fragility (positions 0-31) does NOT predict late fragility. The correlation at cutpoint 32 is r=0.08 — essentially zero. The "bifurcation" is not a sharp transition but a gradual divergence that only becomes reliable at position ~96.

---

## The Data

### Early vs Late Fragility Correlation (how well does early predict late?)

| Cutpoint | r(early, late) |
|----------|---------------|
| 8 | 0.017 |
| 16 | -0.055 |
| **32** | **0.080** |
| 48 | 0.167 |
| 64 | 0.191 |
| **96** | **0.331** |

At position 32, the correlation is r=0.08 — essentially random.

### Cumulative Prefix Prediction of Overall Fragility

| Prefix length | r(prefix, overall) | Quartile accuracy |
|--------------|-------------------|-------------------|
| 1 | 0.029 | 50.5% |
| **2** | **0.359** | 50.5% |
| 8 | 0.091 | 59.3% |
| 16 | 0.059 | 58.6% |
| **32** | **0.303** | **79.3%** |
| 64 | 0.592 | 91.2% |
| **96** | **0.773** | **98.8%** |
| 128 | 0.875 | 99.9% |

**Key pattern**: Position 2 alone has r=0.36 (a genuine signal from the "The" vs "I" choice). But positions 3-16 ADD NOISE — the universal high-fragility plateau dilutes the signal. Starting at ~position 20, each new position adds genuine information. Reliable prediction (r > 0.75) requires ~96 tokens.

## Why the Claim Was Wrong

The "bifurcation by position 32" came from Iteration 4's observation that:
- Unflippable prompts have margins that ramp up from positions 0-32
- Fragile prompts have a margin DIP at positions 16-31

This is true as an **average pattern across groups**. But when applied to **individual prompts**, positions 0-31 are dominated by the universal position-2 spike and opening-phrase uncertainty. This affects ALL prompts similarly. The per-prompt differentiation only emerges gradually.

## State-Switching Prompts

- **Calm→Fragile** (early < 0.3, late > 0.15): **43 prompts** (2.2%). All reasoning — these are word problems that start by echoing the question (calm) then enter computation (fragile).
- **Fragile→Calm** (early > 0.5, late < 0.05): **383 prompts** (19.2%). Dominated by continuation (the model's opening uncertainty resolves into confident narrative).

The asymmetry is dramatic: 383 vs 43. Most prompts START fragile (position 2) then calm down. Very few get MORE fragile.

## Trajectory Shape Distribution

| Shape | Count | % |
|-------|-------|---|
| Monotonic decrease | 1,513 | **75.7%** |
| Irregular | 463 | 23.2% |
| Flat | 12 | 0.6% |
| U-shape | 11 | **0.6%** |

76% of prompts simply get less fragile over time. No trajectory reversal, no bifurcation — just a steady decline from the position-2 peak.

## The Early32 Distribution Is Not Bimodal

If there were a true bifurcation at position 32, we'd expect a bimodal distribution of early fragility (two clusters). Instead:

```
early32 histogram:
  0.00-0.05:    1
  0.10-0.15:    8
  0.15-0.20:   64
  0.25-0.30:  139
  0.30-0.35:  232     <<<< mode region
  0.40-0.45:  277     <<<< peak
  0.50-0.55:  283     <<<< peak
  0.55-0.60:  287     <<<< peak
  0.65-0.70:  176
  0.75-0.80:   47
```

The distribution is unimodal and continuous. There is no gap, no bifurcation, no two states at position 32.

## Revised Picture

**Old claim**: "By position ~32, the trajectory is determined. There are two attractor states."

**New picture**: There is a **gradual divergence of slopes**, not a sharp bifurcation:
1. Positions 0-31: HIGH fragility for almost ALL prompts (universal sentence-start uncertainty)
2. Positions 32-64: Trajectories begin to separate — reasoning stays high, others decline
3. Positions 64-128: Clear separation emerges (r=0.59-0.77)
4. Positions 128+: Fully separated (r=0.88)

The "two attractor states" from Iteration 4 are better understood as a **continuous spectrum of decline rates**: some prompts decline fast (confident repetition), others decline slowly or not at all (arithmetic loops). But there's no discrete boundary at position 32.

## What Survives From the Attractor State Theory

The observation that:
- Confident prompts produce repetitive text (all margins high) → confirmed
- Fragile prompts enter computation loops (margins collapse) → confirmed
- Reasoning is bimodal → confirmed

But the TIMING claim ("bifurcation by position 32") is wrong. It's a gradient, not a step function. And the critical threshold for prediction is ~96 tokens, not 32.

## Scripts
- `research_loop/iter10_bifurcation.py`
