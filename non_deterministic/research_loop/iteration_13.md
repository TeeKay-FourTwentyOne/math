# Iteration 13: Predicting Fragility From Prompt Text Alone

**Direction**: Can we predict which prompts will diverge across platforms without running inference?

**Answer: Partially — 65% of variance is predictable, dominated by one question: "is this arithmetic?"**

---

## Key Results

### Hierarchy of Predictors (R-squared)

| Predictor | R-squared | What it means |
|-----------|-----------|---------------|
| Category alone | 0.209 (21%) | Knowing reasoning vs continuation vs factual gives modest signal |
| **`is_arithmetic` alone** | **0.503 (50%)** | **A single binary feature explains HALF of all fragility variance** |
| Simple rules (category + type) | 0.578 (58%) | Hand-crafted rules near the ceiling |
| Linear regression (34 features) | **0.645 (65%)** | All features together — only +7% over simple rules |

### Top 5 Individual Feature Correlations

| Feature | r with fragility | Interpretation |
|---------|-----------------|----------------|
| **is_arithmetic** | **+0.71** | Simple arithmetic in reasoning = fragile |
| max_number_digits | +0.49 | Bigger numbers = harder = more fragile |
| has_times | +0.42 | Multiplication marker |
| cat_reasoning | +0.40 | Reasoning category |
| has_digits | +0.40 | Presence of numbers |

**Weakest**: `has_question_mark` (r=0.001 — literally zero signal).

### Practical Triage (from prompt text alone)

| Prompt type | Expected fragility | Divergence at eps=1e-3 |
|------------|-------------------|----------------------|
| Arithmetic (What is X + Y?) | ~28% | High risk |
| Continuation (narrative) | ~14% | Medium risk |
| Reasoning (word problems) | ~5% | Low risk |
| Factual / Conversational / Open-ended | ~8-9% | Low risk |

### Exp vs Ctrl Classification From Text Alone

- Accuracy: **67.2%** (131/200 exp, 138/200 ctrl correctly classified)
- Much worse than the 100% achievable from activation data
- The text-only predictor compresses the predicted range (mean predicted 0.15 vs actual 0.24 for exp)

## What The 35% Unexplained Variance Contains

The remaining 35% not predictable from text includes:
1. **Model-specific generation path**: which token sequences the model produces
2. **Interaction between prompt and model**: same category but different random seeds
3. **Within-type variance**: two addition prompts of the same length can have very different fragility (e.g., "What is 816 + 721?" at 0.26 vs "What is 3099 + 2140?" at 0.59)

The `max_number_digits` feature captures some of this (bigger numbers → more fragile), but much per-prompt variation remains.

## The `max_number_digits` Effect

This is a genuinely interesting micro-finding: prompts with larger numbers produce more fragile output, with r=+0.49. The mechanism: 4-digit addition requires more uncertain digit computations than 2-digit addition. Pythia-410M is more confident about "24 + 37" than "3947 + 2116" because the larger numbers require carrying digits that the model hasn't memorized.

## Linear Regression Coefficients (Top 10)

| Feature | Coefficient | Direction |
|---------|------------|-----------|
| is_arithmetic | +0.121 | Arithmetic → 12pp more fragile |
| has_times | +0.112 | Multiplication → 11pp more fragile |
| has_digits | -0.093 | Having digits (non-arithmetic) → actually less fragile |
| cat_continuation | +0.089 | Continuation → 9pp more fragile |
| kw_buy | -0.064 | Money word problems → less fragile |
| cat_reasoning | +0.044 | Reasoning (non-arithmetic) → 4pp more fragile |
| max_number_digits | +0.031 | Each additional digit → 3pp more fragile |
| has_plus | +0.031 | Addition → 3pp more fragile |
| kw_capital | -0.025 | "Capital of X" → less fragile |

Note the `has_digits: -0.093`: digits in NON-arithmetic prompts (e.g., "What year did...") are associated with LOWER fragility. The model is confident about repeating numbers from the prompt.

## Quartile Prediction

| Method | Exact quartile | Within 1 quartile |
|--------|---------------|-------------------|
| Simple rules | 46.3% | 91.5% |
| Linear regression | 45.7% | 90.4% |

91% within-1-quartile accuracy means you can roughly rank prompts by fragility from text alone, but fine-grained ordering requires inference.

## Scripts
- `research_loop/iter13_prompt_prediction.py`
