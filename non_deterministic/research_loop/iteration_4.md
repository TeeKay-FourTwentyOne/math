# Iteration 4: Anatomy of the "Unflippable" Prompts

**Direction**: Characterize the 70 control prompts that remain stable even at ε=1e-1 (all 256 margins > 0.1).

---

## The Punchline

Unflippable prompts are ones where **the model enters a highly repetitive generation loop** — either echoing the prompt text or repeating a simple phrase structure. Repetition makes each token predictable → large margins → stability. This is the mirror image of the fragile arithmetic loops: both are repetitive, but unflippable loops repeat *known text* while fragile loops repeat *uncertain computations*.

---

## Category Distribution: Reasoning Dominates Both Extremes

| Category | Control (40 each) | Unflippable | % Unflippable |
|----------|-------------------|-------------|---------------|
| **reasoning** | 40 | **24** | **60.0%** |
| factual | 40 | 19 | 47.5% |
| open_ended | 40 | 14 | 35.0% |
| conversational | 40 | 10 | 25.0% |
| **continuation** | 40 | **3** | **7.5%** |

Reasoning is simultaneously the MOST fragile category (highest mean fragility) and the MOST likely to be unflippable. This confirms its bimodal nature:
- Short arithmetic (multiplication, addition): extremely fragile
- Word problems and sequences: extremely confident (echo the question)

Continuation prompts are almost never unflippable (only 7.5%) because narrative generation always encounters uncertain moments — which branch to take, which word to choose.

## Margin Curve Shape: The Distinguishing Feature

Average margin by generation window:

| Window | Unflippable | Flippable ctrl | Most flippable exp |
|--------|-------------|----------------|-------------------|
| pos 0-7 | 1.90 | 1.56 | **2.49** |
| pos 8-15 | 2.63 | 1.97 | 1.77 |
| pos 16-31 | 4.09 | 3.30 | **1.23** |
| pos 32-63 | 5.75 | 5.05 | 1.57 |
| pos 64-127 | 7.23 | 6.61 | 3.72 |
| pos 128-191 | 7.92 | 7.43 | 4.30 |
| pos 192-255 | 8.35 | 7.92 | 5.90 |

**Key pattern**: All three groups show monotonically increasing margins in late generation (pos 64+). The difference is in the opening:
- Unflippable and flippable control both ramp UP from position 0
- **Most flippable experimental actually STARTS confident (2.49) then DROPS to 1.23** at positions 16-31 before slowly recovering

This dip at positions 16-31 is the signature of entering an arithmetic computation zone where the model loses confidence. It's the inverse of the unflippable pattern.

## What Unflippable Generated Text Looks Like

**100% of unflippable prompts generate repetitive text.** Concrete examples:

1. **"Invent and describe a pillow that records your dreams"** (min_margin=0.498):
   > The pillow is a pillow that records your dreams.\n\nThe pillow is a pillow that records your dreams.

2. **"What comes next: 45, 56, 67, 78, 89, ?"** (min_margin=0.490):
   > \#\#\# **45**\n\n**45**\n\n**45**\n\n**45**\n\n**45**...

3. **"A car travels at 87 mph for 11 hours. How far does it go?"** (min_margin=0.403):
   > The car travels at 87 miles per hour for 11 hours. How far does it go?\n\nThe car travels at 87 miles...

The model is **confidently stuck** — it found a comfortable pattern and repeats it with high certainty. The text is low-quality but high-confidence.

## Generated Text Properties

| Property | Unflippable | Flippable ctrl | Most flippable |
|----------|-------------|----------------|----------------|
| Gen text length | 926 chars | 894 chars | 766 chars |
| Has paragraph breaks | 87% | 68% | 41% |
| Has repetition | **100%** | 99% | 59% |
| Starts with "The" | **64%** (45/70) | 27% (35/130) | 17% (12/70) |

Unflippable prompts overwhelmingly start with "The" (64%) — the model's most confident sentence-start. When the model starts with "The" followed by echoing the prompt, it locks into a repetitive pattern with consistently high margins.

## Position 2 Does NOT Distinguish Them

Both groups are equally fragile at position 2:
- Unflippable: 82.9% fragile (frac < 1.0)
- Flippable ctrl: 86.9% fragile

The sentence-start ambiguity is universal. What makes prompts unflippable is everything AFTER position 2 — whether the model enters a confident repetitive loop or faces uncertain territory.

## Simple Predictors

| Feature | Unflippable | Flippable ctrl | Cohen's d |
|---------|-------------|----------------|-----------|
| Mean margin | 6.99 | 6.44 | **0.77** |
| Median margin | 7.12 | 6.59 | 0.66 |
| Prompt text length | 59 chars | 76 chars | **-0.51** |
| Prompt word count | 11.6 | 13.9 | -0.38 |

Mean margin is the best predictor (d=0.77), but prompt length also matters (d=-0.51): shorter prompts are more likely to be unflippable because their responses tend to be more formulaic.

## The Barely-Unflippable Edge Cases

The 10 prompts closest to the unflippable threshold (min_margin 0.104-0.117):
- "In what year did the fall of the Roman Empire occur?" (factual, min_margin=0.104)
- "Good morning! I need your help with something." (conversational, 0.108)
- "Write a paragraph about chaos." (open_ended, 0.113)

These are prompts that are mostly confident but have ONE moment of near-uncertainty. A slightly different prompt phrasing could push them below the threshold.

## Theoretical Implication

The model has two "attractor states" for text generation:
1. **Confident repetition**: echo/rephrase the prompt or a simple pattern → all margins high → unflippable
2. **Uncertain computation**: attempt arithmetic or creative branching → some margins tiny → fragile

The bifurcation between these states happens early in generation (positions 2-15). By position 32, the model's trajectory is largely determined — either in a confident loop or not.

## Scripts
- `research_loop/iter4_unflippable.py`
