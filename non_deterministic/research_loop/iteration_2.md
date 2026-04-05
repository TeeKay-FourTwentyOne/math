# Iteration 2: The Position-2 Sentence-Start Ambiguity

**Direction**: Verify the "position 16 universal decision point" claim and map the fine-grained position-fragility curve.

**Verdict**: The position-16 claim is technically correct (49.6% at threshold 1.0) but **deeply misleading**. Position 16 is not a decision point — it's an unremarkable point on a smooth decline. The real story is **position 2**, which hits 82.9% fragility.

---

## The \n\n Prefix Pattern

Pythia-410M almost always starts its response with two newlines:
- **Position 0**: `\n` (newline) in 83% of prompts. Fragility: 17%.
- **Position 1**: `\n` again in 72.5% of prompts. Fragility: 23.5%.
- **Position 2**: The first actual content token. Fragility: **82.9%** — the true peak.

The double-newline prefix is a training artifact: the model learned that answers follow prompts with a blank line separator.

## The "The" vs "I" Contest

At position 2, the model faces a near-coin-flip between sentence-starting words. The top competing pairs (margin < 1.0):

| Rank | Token 1 | Token 2 | Count | % of fragile |
|------|---------|---------|-------|-------------|
| 1 | 'The' | 'I' | 208 | 12.5% |
| 2 | 'I' | 'The' | 190 | 11.5% |
| 3 | 'I' | 'A' | 120 | 7.2% |
| 4 | 'The' | 'A' | 101 | 6.1% |
| 5 | '"' | 'The' | 84 | 5.1% |
| 6 | 'A' | 'The' | 60 | 3.6% |

**42.4% of all position-2 fragility** comes from the model choosing between "The...", "I...", "A...", and '"...' as sentence starters. This is a stylistic/structural choice (third-person vs first-person vs indefinite vs dialogue), not a content decision.

## Per-Category Position-2 Fragility

| Category | Fragility at pos 2 | Interpretation |
|----------|-------------------|----------------|
| continuation | **94.5%** | Many valid ways to start a narrative |
| open_ended | 87.8% | Creative responses have open structure |
| reasoning | 85.3% | Echo vs rephrase decision |
| factual | 79.3% | Somewhat more constrained answers |
| conversational | 68.0% | Most constrained response patterns |

## The True Position Curve (threshold=1.0)

There is no cliff at position 16. The curve is:
- **Position 0**: 17.0% (newline is dominant)
- **Position 2**: **82.9%** (sentence-start ambiguity — the true peak)
- **Largest single drop**: pos 2→3 (-23.1 pp, from 82.9% to 59.8%)
- **Positions 3-18**: Noisy plateau at ~46-63% (building the opening phrase)
- **Positions 19-31**: Gradual decline to 28%
- **Positions 32-63**: ~12-34%
- **Positions 64-127**: ~7%
- **Positions 128+**: ~3-4%

The decline is smooth and monotonic (after the pos-2 spike). Position 16 (49.6%) sits at roughly the geometric midpoint of the curve, which is why it seemed like a threshold — but it has no structural significance.

## Position 0 Conditional Analysis

Position 2 fragility depends on what position 0 was:

| Position 0 token | N prompts | Pos 2 fragility |
|-------------------|-----------|-----------------|
| `\n` | 1659 | 81.4% |
| ` she` | 11 | **100.0%** |
| `It` | 29 | 96.6% |
| ` he` | 174 | 95.4% |
| `I` | 20 | 95.0% |
| ` a` | 49 | 87.8% |
| ` the` | 12 | 83.3% |
| ` I` | 30 | 60.0% |

When position 0 is NOT a newline (the model starts with actual content), position 2 fragility is even HIGHER (95%+). The only exception is ` I` (60%), where the pronoun already constrains the sentence structure.

## Concrete Examples

**Most fragile at position 2** (margin = 0.000393):
- Prompt: "Who wrote the novel 'The Picture of Dorian Gray'?"
- Generated: `\n` `\n` **`The`** `novel` `was` `written`
- Contest: `'The'` vs `'I'` (margin 0.0004 — effectively a coin flip)
- On a different platform, this prompt would start with "I believe..." instead of "The novel was..."

**Least fragile at position 2** (margin = 7.04):
- Prompt: "What is 315 divided by 15?"
- Generated: `\n` `\n` ` 315` ` divided` ` by` ` 15`
- The model confidently echoes the question before computing.

## Increasing Fragility (rare counter-examples)

Only **3 out of 1999 prompts** show fragility that INCREASES during generation:
1. "If 3 workers can finish a job in 7 days, how many days would 7 workers take?"
2. "How many hours does it take to travel 189 miles at 29 miles per hour?"
3. "How many hours does it take to travel 316 miles at 22 miles per hour?"

All three are reasoning prompts that initially echo the question (low fragility) then enter an arithmetic computation phase (high fragility). The fragility curve is U-shaped: [0.22, 0.09, 0.16, 0.38, 0.53, 0.38, 0.31, 0.19] over 8 windows of 32 tokens.

## Implications for Cross-Platform Experiments

1. **Position 2 is the highest-risk single position**: 83% of prompts will produce the first divergent token here if cross-platform FP noise exceeds the margin. Since the contest is typically "The" vs "I", the divergence will be SEMANTICALLY SIGNIFICANT — entire response styles will differ.

2. **The divergence cascade**: A flip at position 2 from "The" to "I" changes the entire sentence structure. This means early-position divergences produce maximally different outputs, while late-position divergences (which are more likely to be digit-for-digit swaps) produce minimal semantic change.

3. **Category targeting**: Continuation prompts (94.5% fragile at pos 2) will show the most position-2 divergence. Conversational prompts (68%) will show the least.

## Previous Iteration Error Corrected

The verify_position16.py script from iteration 1 used threshold 0.1 (finding 7.3% fragility at pos 16) to "debunk" the 49.6% claim. This was a threshold mismatch — the correct threshold is 1.0, under which 49.6% is indeed correct. The claim was not wrong in its number, but wrong in its interpretation: position 16 is not special.
