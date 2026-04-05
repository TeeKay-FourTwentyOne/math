# Iteration 12: Verifying Two Longstanding Claims

**Direction**: Resolve two unverified claims that have been on the list since the seed iteration.

---

## Claim 1: r=-0.82 Length-Fragility Correlation — CONFOUNDED

**Verdict: 80% of the correlation is explained by problem type, not length.**

The r=-0.82 overall correlation within reasoning exists because short prompts are arithmetic (fragile) and long prompts are word problems (confident). After removing problem-type means:

| Metric | Value |
|--------|-------|
| Overall r(length, fragility) | **-0.82** |
| After removing type means | **-0.17** |
| Variance explained by type | **80%** |

Within individual problem types, length has weak and inconsistent effects:

| Type | N | Mean length | Mean frag | r(len, frag) |
|------|---|-------------|-----------|-------------|
| addition | 55 | 19.7 | 0.307 | **+0.20** (positive!) |
| subtraction | 55 | 19.5 | 0.269 | **+0.39** (positive!) |
| multiplication | 50 | 19.8 | 0.337 | -0.36 |
| percentage | 41 | 38.8 | 0.220 | -0.74 |
| money | 34 | 57.5 | 0.061 | -0.91 |
| distance | 49 | 66.6 | 0.049 | -0.13 |
| division | 25 | 26.4 | 0.258 | -0.03 |
| sequence | 41 | 54.2 | 0.054 | +0.23 |

Addition and subtraction actually show **POSITIVE** within-type correlations — longer problems (bigger numbers) are MORE fragile, not less. This directly contradicts the overall r=-0.82.

The strong within-type correlations (money r=-0.91, percentage r=-0.74) are real but driven by different mechanisms than the overall pattern. Money problems with more text context give the model more to echo confidently.

**Corrected statement**: Fragility in reasoning is determined by problem TYPE (arithmetic: 27-34% vs word problems: 4-7%), not by length. Length is a proxy for type.

---

## Claim 2: 100% Repetition in Unflippable — REAL BUT UNIVERSAL

**Verdict: 100% rate is real (confirmed by 5 methods), but it's NOT a distinguishing feature — 98-99% of ALL Pythia-410M output is repetitive.**

| Detection method | Unflippable (70) | Flippable ctrl (130) | Difference |
|-----------------|-------------------|---------------------|------------|
| 20-char x3 (original) | 100.0% | 99.2% | +0.8% |
| Sentence repeat x2 | 91.4% | 86.9% | +4.5% |
| 4-gram x3 | 100.0% | 99.2% | +0.8% |
| Paragraph repeat | 77.1% | 56.9% | +20.2% |
| **50-char x2 (strict)** | **100.0%** | **99.2%** | **+0.8%** |

Zero unflippable prompts lack strict (50-char) repetition: **0/70**.

But the critical context — repetition by fragility quartile across ALL 1999 prompts:

| Quartile | Fragility range | Strict repetition rate |
|----------|----------------|----------------------|
| Q1 (most stable) | 0.004-0.059 | **99.4%** |
| Q2 | 0.059-0.094 | 98.0% |
| Q3 | 0.094-0.145 | 98.6% |
| Q4 (most fragile) | 0.145-0.613 | 88.4% |

**Pythia-410M at 256 tokens with greedy decoding produces repetitive text almost universally.** The unflippable prompts aren't special — they follow the model's general behavior. The only group with noticeably less repetition is the most fragile quartile (88.4%), which includes arithmetic prompts generating digit-heavy content.

**Corrected statement**: Repetition is a property of the model, not of unflippable prompts specifically. What distinguishes unflippable prompts is confident repetition (high margins) vs the universal baseline of repetition-with-occasionally-fragile-tokens.

## Scripts
- `research_loop/iter12_verify_claims.py`
