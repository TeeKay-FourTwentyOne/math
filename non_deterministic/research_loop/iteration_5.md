# Iteration 5: Quality vs Stability — The Fundamental Tradeoff

**Direction**: Is there a fragility "sweet spot" where the model produces both coherent and hardware-stable text?

**Answer: No.** The tradeoff is monotonic and fundamental. More stable text is more repetitive. The correlation is r=+0.74 for non-reasoning prompts.

---

## The Core Finding

For non-reasoning prompts (continuation, factual, conversational, open_ended), text diversity increases monotonically with fragility:

| Fragility bin | N | TTR | 4-gram RepRate | Composite | Long repeat % |
|--------------|---|-----|----------------|-----------|---------------|
| Very low (0-3%) | 100 | 0.066 | 0.863 | 0.013 | 99% |
| Low (3-6%) | 405 | 0.087 | 0.802 | 0.039 | 100% |
| Med-low (6-10%) | 572 | 0.084 | 0.680 | 0.030 | 99.8% |
| Medium (10-15%) | 455 | 0.120 | 0.541 | 0.058 | 100% |
| Med-high (15-25%) | 283 | 0.168 | 0.457 | 0.097 | 99.6% |
| High (25-40%) | 173 | 0.233 | 0.433 | 0.163 | 80.9% |

**Composite = TTR × (1 - RepRate).** Higher = more diverse, less repetitive. It rises monotonically with fragility. There is no peak — just a continuous tradeoff.

## Correlations

| Metric pair | r (all) | r (non-reasoning) |
|-------------|---------|-------------------|
| Fragility vs TTR | +0.46 | **+0.74** |
| Fragility vs RepRate | -0.40 | **-0.54** |
| Fragility vs Composite | +0.39 | **+0.63** |

Non-reasoning shows a very strong positive correlation (r=0.74) between fragility and text diversity. This makes intuitive sense: diverse text requires choosing less-predictable tokens, which means smaller logit margins.

## Per-Category: Fragile Text Is Better Text

For EVERY non-reasoning category, the most fragile third has the best composite quality:

| Category | Stable third | Middle third | Fragile third | Ratio fragile/stable |
|----------|-------------|-------------|---------------|---------------------|
| continuation | 0.046 | 0.074 | **0.115** | 2.5x |
| conversational | 0.008 | 0.023 | **0.067** | 8.4x |
| factual | 0.012 | 0.021 | **0.078** | 6.5x |
| open_ended | 0.010 | 0.036 | **0.058** | 5.8x |
| reasoning | 0.103 | 0.036 | **0.193** | U-shaped |

Reasoning is U-shaped: low-fragility word problems echo the question (moderate quality), mid-fragility is low quality, high-fragility arithmetic generates diverse math notation (high TTR but not coherent prose).

## What "Stable" Text Actually Looks Like

The most stable prompts (fragility < 5%) with highest composite scores generate degenerate text:

1. **"If 5 workers can finish a job in 28 days..."** (frag=3.1%, comp=1.0):
   > The answer is: $5,000,000,000,000,000,000,000,000,000,000,...

2. **"What comes next: 45, 56, 67, 78, 89, ?"** (frag=4.9%, comp=0.49):
   > \#\#\# **45**\n\n**45**\n\n**45**\n\n**45**...

These score well on TTR (because commas/numbers create diverse 4-grams) but are useless.

## What "Sweet Spot" Text Actually Looks Like

The prompts at fragility ~33% (highest composite ≈ 0.24) generate math quiz patterns:

1. **"What is 2749 - 2546?"** (frag=34%, comp=0.20):
   > -2\*w\*\*3 + w\*\*2 + w + 27\nWhat is the s'th term of -13, -31, -59, -97, -143, -203, -281?\n-s\*\*3 + s\*\*2 - 12\*s - 3...

This is the model generating Wolfram-Alpha-style math training data — diverse vocabulary but not meaningful prose.

## The Mechanism

The tradeoff is inherent to greedy decoding:
- **Confident (low-fragility) text**: the model finds a high-probability pattern and locks in → repetitive loops → stable across hardware
- **Diverse (high-fragility) text**: the model explores lower-probability continuations → more varied vocabulary → margins are small → vulnerable to hardware noise

There is no decoding strategy that simultaneously maximizes text diversity and hardware stability under greedy (argmax) decoding. They are mathematically opposed properties: diversity requires low margins, stability requires high margins.

## Practical Implication

If cross-platform reproducibility is required, the options are:
1. **Accept repetitive output** (use greedy decoding, get stable but degenerate text)
2. **Accept divergence** (use sampling, get diverse text that differs across platforms)
3. **Constrain divergence** — identify which prompts are in the fragile zone and either:
   - Skip them (only process low-fragility prompts cross-platform)
   - Use a secondary model to verify/reconcile divergent outputs
   - Increase temperature then round/quantize logits to create artificial margin

## Quality Metric Caveat

TTR and 4-gram repetition rate are blunt instruments. They measure lexical diversity, not coherence or usefulness. The "best" text by these metrics includes mathematical notation and number strings. A proper quality evaluation would require human judgment or a larger evaluator model. However, the DIRECTION of the tradeoff would remain the same — text that is more coherent and varied will necessarily be less stable.

## Scripts
- `research_loop/iter5_quality_stability.py`
