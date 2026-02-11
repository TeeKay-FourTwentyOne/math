# External Input from Gemini Deep Research & Claude Opus 4.6

**Date**: 2026-02-10
**Source**: Gemini Deep Research (literature) + Claude Opus 4.6 (web UI analysis)

---

## 1. For Computational Progress (n=32, m=63)

**Reference**: Đoković-Kotsireas, "Compression of periodic complementary sequences
and applications," Designs, Codes and Cryptography, 2015.

**Key idea**: Exploit multiplier group orbits for composite m to compress the search
space. Their sequence compression methods are specifically designed for periodic
complementary sequences over composite moduli.

**Action**: Implement their compression technique, apply to m=63 = 9×7.
This is an engineering task — the algorithm is described in the paper.

**Why this matters**: n=32 (m=63) is our first unsolved case. SA gets stuck at cost=4
(2 violated positions). The compression might reveal structure that makes the search
tractable, or prove infeasibility.

---

## 2. For General Proof (p ≡ 3 mod 4)

**Reference**: Szekeres, "Tournaments and Hadamard matrices," Enseignement
Mathématique 15, pp. 269–278, 1969.

**Key idea**: Szekeres' cyclotomic splitting technique constructs supplementary
difference sets by splitting cyclotomic classes. The technique works over Z_{(p-1)/2}.

**Research question**: Can the cyclotomic splitting be adapted from Z_{(p-1)/2} to Z_p
with asymmetric threshold constraints (our problem)?

**Connection to our findings**: The k=n-2 universality (Theorem 5) tells us exactly
what |D11| to target. This is a crucial input for adapting the Szekeres construction.

**Nature of task**: Mathematical investigation, not coding. Needs careful reading of
the construction and analysis of whether character sum estimates transfer to our
setting.

---

## 3. CRITICAL PROOF INSIGHT: Energy/Parseval Argument

**From Opus 4.6, responding to our gap/std findings:**

> The proof needs to go through an energy/Parseval argument, not a pointwise
> probability argument. The right framework is: bound the total excess
> Σ_{d ∈ D₁₁} max(0, A(d)+B(d)−(n−2)) using the total energy constraint,
> not bound each Pr[A(d)+B(d) > n−2] independently.

**Why this changes the proof strategy:**

Our pointwise approach requires Lemma L6 (joint probability bound over all d ∈ D11),
which requires handling correlations (the hardest step). The energy/Parseval approach
BYPASSES this entirely by working with the total violation.

**The mechanism (Opus 4.6):**

> Choosing |D₁₁| = n−2 (fewer tight positions) and |D₂₂| = n (more loose positions)
> creates an asymmetry where the Parseval constraint works in your favor: there's more
> room to absorb excess on the loose side than deficit needed on the tight side.

**Parseval identity in our context:**

For any (D11, D12):
  Σ_{d=1}^{p-1} (A(d) + B(d)) = |D11|(|D11|-1) + |D12|(|D12|-1) = FIXED

At k = n-2:
- D11 positions (|D11| = n-2): threshold n-2, average A+B ≈ n-3
- D22 positions (|D22| = n): threshold n-3, average A+B ≈ n-3

The total sum is fixed. If we can show the total EXCESS at D11 positions is bounded,
then the existence of a valid pair follows.

**Key advantage**: The energy argument works GLOBALLY, not pointwise. It naturally
handles correlations because it deals with the SUM, not individual events.

**How to formalize:**

Define excess E = Σ_{d ∈ D11} max(0, A(d)+B(d)-(n-2)).
If E = 0, all D11 constraints are satisfied.
By Parseval/total energy:
  Σ_{d ∈ D11} A(d)+B(d) = fixed - Σ_{d ∈ D22} A(d)+B(d)
If the D22 values don't exceed their (loose) threshold much, the D11 sum is controlled.
Combined with variance bounds, this might show E = 0 for some (D11, D12).

**Status**: This is the most promising proof approach identified so far.
Needs formalization using the exact moments (E[B] = (p-3)/4, Var[B] = p/16).

---

## Integration with Our Findings

| Our finding | How it connects |
|-------------|----------------|
| k=n-2 universality (Theorem 5) | Tells Szekeres adaptation what |D11| to target |
| Var[B] = p/16 (exact_moments.py) | Quantifies fluctuation for energy argument |
| E[valid pairs] ≈ 2^{p-2} | Confirms first moment works; energy argument is cleaner |
| Positive association (corr analysis) | Energy argument doesn't need this |
| Gap = -1 exactly | The small gap is what makes the energy argument possible |

---

## Recommended Next Steps (Priority Order)

1. **Formalize the energy/Parseval argument** using k=n-2 formulation and exact moments
2. **Read Szekeres 1969** for cyclotomic splitting technique
3. **Look up Đoković-Kotsireas 2015** for m=63 compression
4. **Combine** energy argument with Szekeres construction for complete proof
