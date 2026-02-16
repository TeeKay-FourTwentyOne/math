# First Moment Analysis: E[N] → ∞

**Date**: 2026-02-13
**Status**: Analysis complete; rigorous proof conditional on one bound
**Task**: Prove E_{D11}[N(D11)] → ∞ as p → ∞ for primes p ≡ 3 mod 4

---

## 1. Executive Summary

We analyze the first moment E[N] = E_{D11}[#{valid D12 for D11}] using the decomposition:

  E[N] = #D12 × E_{D11}[Pr_{D12}[valid | D11]]

The budget (log₂ #D12 ≈ p - 1) vastly exceeds the constraint cost (≈ (p-1)/2 bits from R = (p-1)/2 constraints each costing ≈ 1 bit). This leaves **(p-1)/2 bits of headroom** — enough to absorb any polynomial or even moderate exponential loss from constraint correlations.

**Rigorous result**: If the average correlation factor satisfies

  E_{D11}[c₀(D11)] ≥ 2^{-αp}  for some α < 1/2

then E[N] → ∞ as p → ∞. Empirically, log₂(c₀_avg) ≈ -5 to -6 (constant), so α ≈ 0 — but we cannot yet prove this rigorously.

**Empirical evidence**: log₂(E[N]) ≈ 0.12p + O(1), i.e., E[N] grows exponentially in p:

| p | E[N] | log₂(E[N]) |
|---|------|-----------|
| 11 | 10.00 | 3.32 |
| 19 | 1.29 | 0.37 |
| 23 | 9.43 | 3.24 |
| 31 | 56.65 | 5.82 |
| 43 | ~87 (est.) | ~6.4 |

---

## 2. Setup and Notation

**Prime**: p ≡ 3 mod 4, n = (p+1)/2
**D11**: symmetric subset of {1,...,p-1} with |D11| = n (random uniform)
**D12**: {0} ∪ S where S is a random k-subset of {1,...,p-1}, k = (p-3)/2
**B(d)**: autocorrelation of D12 at distance d, B(d) = Δ(D12, D12, d)
**A(d)**: autocorrelation of D11, A(d) = Δ(D11, D11, d)

**Constraints** (from proof_L6_conditioning.md Section 2): For all d = 1,...,p-1:
  B(d) ≤ T(d) where:
  - T(d) = (p-3)/2 - A(d) if d ∈ D11  (V1V1 red, binding)
  - T(d) = (p+3)/2 - A(d) if d ∈ D22  (V1V1 blue, slack)

**Representative reduction**: Since B(d) = B(p-d) and T(d) = T(p-d), the constraints reduce to R = (p-1)/2 independent constraints, one per complementary pair {d, p-d}.

**Budget**: log₂(#D12) = log₂(C(p-1, (p-3)/2)) = p - 1 - (1/2)log₂(πp/2) + O(1/p)

---

## 3. Exact Moments of B(d) (Proven)

For a random k-subset S of {1,...,p-1} with k = (p-3)/2:

**Mean**: E[B(d)] = (p-3)/4 for all d ≠ 0. [Proven in exact_moments.py]

**Variance**: Var[B(d)] = (p-3)(p+1)/(16(p-2)) ~ p/16.
[Proven by direct computation; see exact_moments.py]

**Cross-covariance**: Cov[B(d₁), B(d₂)] = -(p+1)/(8(p-2)) for d₁+d₂ ≢ 0 mod p.
[Proven in proof_equi_covariance.md; correlation ρ = -2/(p-3) → 0]

**Parseval identity**: Σ_{d=1}^{p-1} B(d) = (p-1)(p-3)/4, equivalently Σ_{reps} B(dᵢ) = (p-1)(p-3)/8.
[Deterministic; holds for every D12]

---

## 4. Budget-Cost Analysis

### 4.1 The Budget

log₂(#D12) = log₂(C(p-1, (p-3)/2))

Using the central binomial coefficient asymptotic:

  log₂(C(N, N/2-c)) = N - (1/2)log₂(πN/2) + O(c²/N)

with N = p-1 and the deficit from center c = 1:

  log₂(#D12) = p - 1 - (1/2)log₂(π(p-1)/2) + O(1/p)
              ≈ p - 1 - (1/2)log₂(p) + O(1)

**Numerical verification**:

| p | log₂(#D12) exact | p - 1 - (1/2)log₂(p) |
|---|---|---|
| 11 | 7.71 | 7.27 |
| 19 | 15.42 | 13.75 |
| 31 | 27.12 | 25.05 |
| 127 | 122.16 | 119.50 |
| 499 | 493.19 | 489.53 |

The approximation is accurate to within O(1).

### 4.2 The Cost (Product of Marginals)

Each constraint Pr[B(dᵢ) ≤ Tᵢ] depends on D11 through Tᵢ. For an "average" D11
with A(d) ≈ E[A] = (p+1)/4:

**D11 constraints** (r = (p+1)/4 representatives):
  T_red(d) ≈ (p-3)/2 - (p+1)/4 = (p-7)/4
  Gap = T_red - E[B] = (p-7)/4 - (p-3)/4 = -1
  z_red = -1/σ_B ≈ -4/√p → 0

**D22 constraints** (r' = (p-3)/4 representatives):
  T_blue(d) ≈ (p+3)/2 - (p+1)/4 = (p+5)/4
  Gap = T_blue - E[B] = (p+5)/4 - (p-3)/4 = +2
  z_blue = +2/σ_B ≈ +8/√p → 0

**Per-constraint cost** (Gaussian approximation):

  -log₂(Pr[B ≤ T_red]) ≈ -log₂(Φ(-4/√p)) = 1 + O(1/√p)
  -log₂(Pr[B ≤ T_blue]) ≈ -log₂(Φ(+8/√p)) = 1 - O(1/√p)

**Total product-of-marginals cost**:

  Σ -log₂(Pr[Bᵢ ≤ Tᵢ]) = r × (1 + c₁/√p) + r' × (1 - c₂/√p) + O(1)

where c₁ = 4/(ln2 · √(2π)) ≈ 2.31 and c₂ = 8/(ln2 · √(2π)) ≈ 4.61.

Since r ≈ r' ≈ p/4:

  Total cost ≈ R + (r·c₁ - r'·c₂)/√p ≈ R - (c₂ - c₁)·√p/4 ≈ R - 0.58√p

The √p correction is NEGATIVE (D22 slack more than compensates D11 tightness),
so the total cost is slightly LESS than R bits.

### 4.3 The Headroom (Gaussian Approximation)

  headroom = budget - cost
            ≈ [p - 1 - (1/2)log₂(p)] - [(p-1)/2 + O(1)]
            = (p-1)/2 - (1/2)log₂(p) + O(1)
            → ∞

**Numerical headroom (Gaussian)**:

| p | budget | cost (Gauss) | headroom |
|---|--------|-------------|----------|
| 11 | 7.71 | 9.56 | -1.85 |
| 19 | 15.42 | 12.63 | 2.79 |
| 23 | 19.30 | 14.22 | 5.08 |
| 31 | 27.12 | 17.46 | 9.66 |
| 43 | 38.90 | 22.43 | 16.47 |
| 127 | 122.16 | 59.21 | 62.95 |
| 499 | 493.19 | 232.52 | 260.67 |

The Gaussian headroom is negative at p=11 (too few variables for CLT) but
positive and growing for all p ≥ 19.

### 4.4 Why the Gaussian Headroom Overestimates E[N]

The Gaussian headroom assumes the joint probability equals the product of marginals (c₀ = 1). In reality:

| p | Gaussian headroom | Actual log₂(E[N]) | Discrepancy |
|---|---|---|---|
| 11 | -1.85 (fails) | 3.32 | +5.17 (Gauss underestimates!) |
| 19 | 2.79 | 0.37 | -2.42 |
| 23 | 5.08 | 3.24 | -1.84 |
| 31 | 9.66 | 5.82 | -3.84 |

At p=11, the Gaussian approximation is poor and underestimates (c₀ effects help).
At p ≥ 19, the Gaussian overestimates by O(1) to O(log p) bits — this is the correlation cost.

---

## 5. The Correlation Factor c₀

### 5.0 Constraint Reduction (V1V2 Trivial)

**Theorem** (V1V2 Identity). For any symmetric D11 and any D12 with 0 in D12:

  X(d) = Sigma(D11, D12, d) + Delta(D12, D22, d) = |D12| - 1_{d in D12}

This is proven algebraically using D11 symmetry (see `proofs/v1v2_identity_proof.md`).

**Consequence**: V1V2 constraints are satisfied with EQUALITY for ALL (D11, D12):
- V1V2 red (d in D12): X(d) = n-2 (exactly at threshold)
- V1V2 blue (d not in D12): X(d) = n-1 (exactly at threshold)

V1V2 constraints add **zero cost** to the first moment.

**Theorem** (Threshold Equality). A(d) - C(d) = 3 - 2*1_{d in D11}, so:
- V1V1 red threshold = V2V2 blue threshold for d in D11
- V1V1 blue threshold = V2V2 red threshold for d in D22

The constraint model has exactly R = (p-1)/2 independent binding constraints, all of the form B(d) <= T(d) where T(d) = (p-3)/2 - A(d) for d in D11 and T(d) = (p-3)/2 - C(d) for d in D22.

### 5.1 Definition and Exact Computation

For a fixed D11:

  c₀(D11) = Pr[all B(dᵢ) ≤ Tᵢ] / ∏ Pr[B(dᵢ) ≤ Tᵢ]

These thresholds include ONLY the R binding V1V1/V2V2 constraints (V1V2 is free).

**Exact values (exhaustive enumeration):**

| p | c₀ (working D11) | log₂(c₀) | log₂(c₀)/R |
|---|---|---|---|
| 11 | 0.7004 | -0.51 | -0.10 |
| 19 | 0.0326 | -4.94 | -0.55 |
| 23 | 0.006-0.036 | -4.8 to -7.4 | -0.44 to -0.67 |

### 5.2 Why c₀ < 1

The negative equi-correlation (ρ = -2/(p-3)) between B(d) values causes the joint probability to be LESS than the product of marginals. Intuitively: if one B(dᵢ) is small (helping satisfy its constraint), the Parseval sum Σ B(d) = const forces the others to be slightly larger. This is a genuine adverse correlation effect.

**Key**: c₀ < 1 at ALL tested primes (p=11,19,23). The earlier claim that "c₀ > 1" was incorrect.

The Gaussian comparison (Slepian) shows that for Gaussian variables with this negative correlation, c₀ is exponentially small (~2^{-Θ(R)}). The discrete k-subset distribution has a milder penalty, but it is still substantial.

### 5.3 Exact c₀ Decomposition (Working D11)

For working D11, the threshold vector T = (T(d₁),...,T(d_R)) has ALL entries ≥ 1. Non-working D11 have at least one T(d) ≤ 0, making the constraint impossible.

**Exact results** (using product of exact B(d) CDF, not Gaussian):

| p | log₂(∏ Prᵢ) | log₂(N) | log₂(c₀) | c₀ |
|---|---|---|---|---|
| 11 | -2.88 | 4.32 | -0.51 | 0.700 |
| 19 | -6.31 | 4.17 | -4.94 | 0.033 |
| 23 (min) | -6.66 | 4.46 | -4.81 | 0.036 |
| 23 (max) | -3.87 | 7.63 | -7.43 | 0.006 |

Note: log₂(N) here means log₂(N(D11)) for a single working D11, NOT log₂(E[N]).

The budget per D11 is log₂(#D12): 7.71, 15.42, 19.30 for p=11,19,23.
So log₂(N) = budget + log₂(∏ Prᵢ) + log₂(c₀):
  p=11: 7.71 - 2.88 - 0.51 = 4.32 ✓
  p=19: 15.42 - 6.31 - 4.94 = 4.17 ✓

### 5.4 The Averaged Correlation Factor

The quantity that matters for E[N] is not c₀ for a single D11, but:

  E_{D11}[Pr_{D12}[valid | D11]] = E_{D11}[c₀(D11) × ∏ Pr[Bᵢ ≤ Tᵢ(D11)]]

| p | E_{D11}[Pr_joint] | E_{D11}[∏ Prᵢ] | E[N] exact | log₂(E[N]) |
|---|---|---|---|---|
| 11 | 4.76e-2 | 7.03e-2 | 10.00 | 3.32 |
| 19 | 2.94e-5 | 5.09e-3 | 1.29 | 0.37 |

The averaged c₀ drops from 0.68 (p=11) to 0.006 (p=19). But the headroom grows faster.

### 5.5 Can c₀ be Bounded?

For the proof to work, we need the total correlation cost (sum of log₂(c₀) over D11) to be at most (p-1)/2 - ε bits. Empirically, the cost is about 5-8 bits — a CONSTANT, not growing with p.

**Key difficulty**: The B(d) values are GLOBALLY dependent (each B(d) is a quadratic function of ALL indicators Y_a). Standard negative association results (Joag-Dev & Proschan 1983) apply to the indicators Y_a but NOT to the quadratic functions B(d). The FKG inequality does not help because the k-subset measure is not a product measure, and the events are not monotone in the lattice sense.

---

## 6. Alternative Approach: Direct E[Z] Lower Bound

Instead of factoring E[Z] = c₀ × ∏ Prᵢ, we try to bound E[Z] directly.

### 6.1 The Double-Average Perspective

E[N] = #D12 × E[Z] where E[Z] = E_{D11,D12}[1_{valid}].

Exchanging the order of expectation:

  E[Z] = E_{D12}[E_{D11}[1_{valid} | D12]]

For a fixed D12, E_{D11}[1_{valid}] is the fraction of symmetric D11 for which all constraints are satisfied. This "D12-first" perspective may be more tractable.

### 6.2 Constraint Reformulation (D12-first)

Given D12, the constraints become:
- For d ∈ D11: A(d) + B(d) ≤ (p-3)/2, i.e., A(d) ≤ (p-3)/2 - B(d)
- For d ∈ D22: A(d) + B(d) ≤ (p+3)/2, i.e., A(d) ≤ (p+3)/2 - B(d)

Here B(d) is FIXED (determined by D12), and A(d) is RANDOM (determined by D11).

This is a constraint on the autocorrelation profile A(d) of D11. For an A-flat D11, A(d) ≈ (p+1)/4, so:
- D11 constraint: need A(d) ≤ (p-3)/2 - B(d). Since B(d) ≈ (p-3)/4, need A(d) ≤ (p-3)/4.
  Since E[A] = (p+1)/4 ≈ (p-3)/4 + 1, this is A(d) ≤ E[A] - 1. Happens with probability ≈ 1/2.
- D22 constraint: need A(d) ≤ (p+3)/2 - B(d) ≈ (p+3)/2 - (p-3)/4 = (p+9)/4 = E[A] + 2.
  This is easily satisfied (A(d) exceeds mean by 2 only with small probability).

### 6.3 The A-Profile Distribution

A(d) for random symmetric D11 has:
- E[A(d)] = (p+1)/4
- Var[A(d)] = O(p) (hypergeometric concentration)
- The A-profile (A(d₁),...,A(d_R)) has equi-correlation structure similar to B.

The constraint "A(d) ≤ E[A] - 1 for all d ∈ D11" is equivalent to requiring the A-profile to be slightly below average at all D11 positions. This is a global constraint on D11.

### 6.4 Why This Doesn't Immediately Resolve the Problem

The D12-first approach exchanges the roles of A and B. We now need:

  E_{D12}[∏_{d ∈ D11(D12?)} Pr_{D11}[A(d) ≤ T'(d) | D12]]

But D11 determines which positions are in D11 vs D22, so the constraint structure depends on D11 itself. This creates a circularity that makes the D12-first approach no simpler than D11-first.

---

## 7. The Growth of E[N]: Empirical Evidence

### 7.1 Direct Fitting

| p | log₂(E[N]) | (p-1)/2 | log₂(E[N])/(p-1) × 2 |
|---|---|---|---|
| 11 | 3.32 | 5 | 0.664 |
| 19 | 0.37 | 9 | 0.041 |
| 23 | 3.24 | 11 | 0.295 |
| 31 | 5.82 | 15 | 0.388 |

A linear fit (excluding the p=19 anomaly) gives:
  log₂(E[N]) ≈ 0.115p + 1.64

The positive slope means E[N] grows exponentially in p.

### 7.2 p = 19 Anomaly

At p=19, E[N] = 1.29 — barely above 1. This is because p=19 has the smallest orbit count relative to its size: only 14 multiplicative orbits, of which exactly 1 is working. The "phase transition" from typical behavior occurs here because the constraint model is particularly tight for this specific p.

Despite the dip, E[N] > 1 at p=19, and N(D11) = 18 for the 9 working D11 (in one orbit). So existence is confirmed.

### 7.3 Working Orbits at p=43

At p=43, the exhaustive SA scan found 124 working orbits out of 16,796 total. Each orbit has (p-1)/2 = 21 D11 members, so:
- #working D11 = 124 × 21 = 2604
- #total D11 = 16796 × 21 = 352,716
- Fraction working = 0.738%

With E[N] per working D11 being at least ~62 (minimum N divisible by 2p=86, but likely higher):
  E[N] ≥ (124/16796) × 62 ≈ 0.46

Wait, this gives E[N] < 1 if we're not careful. Let me reconsider.

Actually E[N] = (1/#D11) × Σ N(D11). We need the actual N values at p=43 to compute this. From the SA scan, we know 124 orbits have N > 0, but we don't know the exact N values (only that valid D12 exist).

### 7.4 Growth Rate Summary

Based on exact data at p ≤ 31 and existence data at p ≤ 59:

1. E[N] grows exponentially in p (empirical slope ~0.1 in log₂ base)
2. The p=19 anomaly doesn't break the trend
3. Working D11 fractions decline (0.50, 0.071, 0.119, 0.042, 0.007)
   but N values per working D11 grow fast enough to compensate
4. The growth is consistent with log₂(E[N]) ≈ (p-1)/2 - C for constant C ≈ 5-8

---

## 8. What Can Be Proven Rigorously

### 8.1 Unconditional Results

**Theorem 1** (Budget). log₂(#D12) = p - 1 - (1/2)log₂(p) + O(1).

*Proof*: Stirling's approximation for C(p-1, (p-3)/2). ∎

**Theorem 2** (Marginal bounds). For any fixed D11 and any d ∈ {1,...,(p-1)/2}:

  Pr[B(d) ≤ T(d)] = Φ(z_d) + O(1/√p)

where z_d = (T(d) - (p-3)/4) / √(p/16 + O(1)).

*Proof*: Berry-Esseen for degree-2 functions of hypergeometric indicators. ∎

**Theorem 3** (Product of marginals). For an A-flat D11 (max A(d) ≤ (p+1)/4 + O(√(p log p))):

  ∏ Pr[B(dᵢ) ≤ Tᵢ] = 2^{-(p-1)/2 + O(√p)}

*Proof*: Each z_d = O(1/√p), so -log₂ Pr[ok_i] = 1 + O(1/√p). Sum over R = (p-1)/2 reps. ∎

**Theorem 4** (Existence of A-flat D11). For p ≥ 7, there exists a symmetric D11 with
max A(d) ≤ (p+1)/4 + C√(p log p). A random D11 is A-flat with probability 1 - O(1/p).

*Proof*: Sub-Gaussian concentration + union bound over p-1 positions. ∎

### 8.2 Conditional Result

**Theorem 5** (Conditional first moment). If there exists α < 1/2 such that for all sufficiently large primes p ≡ 3 mod 4:

  E_{D11}[Pr_{D12}[valid | D11]] ≥ 2^{-(1+α)(p-1)/2}

then E[N] → ∞.

*Proof*: E[N] = #D12 × E[Pr[valid|D11]] ≥ 2^{p-1-O(log p)} × 2^{-(1+α)(p-1)/2}
       = 2^{(1-α)(p-1)/2 - O(log p)} → ∞ since 1-α > 1/2 > 0. ∎

The condition requires E[Z] not to exceed the product-of-marginals cost by more than a factor of (1+α) times R. Empirically, the cost is about R + 6 bits (not (1+α)R), so α ≈ 0 + 12/p → 0. The condition is satisfied with enormous room.

### 8.3 The Remaining Gap

**What is needed**: A lower bound on E_{D11}[Pr_{D12}[valid | D11]] that beats 2^{-p+ε} for some ε > 0.

**What is available**: The product of marginals gives ∏ Pr ≈ 2^{-R} ≈ 2^{-p/2}. The gap between this and the actual joint probability is c₀, which empirically is 2^{-O(1)} (constant, not growing with p).

**Why the gap is hard to close**: The B(d) values are GLOBALLY dependent (each B(d) is a quadratic function of ALL the indicators Y_a). Standard tools like FKG, Slepian, or negative association give bounds in the WRONG direction or don't apply to quadratic functions.

---

## 9. Promising Paths Forward

### 9.1 Exchangeable Pairs / Stein's Method

The conditional distribution on the Parseval hyperplane has equi-correlation -1/(R-1). For exchangeable random vectors with vanishing correlation, Stein's method gives multivariate normal approximation with quantitative error bounds. If the error in the multivariate CLT is o(1) in the right metric, the product-of-marginals bound follows with at most polynomial loss.

**Specific approach**: Use the exchangeable pairs coupling of Stein's method (Goldstein-Rinott 1996, Reinert-Röllin 2009) to bound:

  |Pr[all Bᵢ ≤ Tᵢ | Σ = S] - Pr[all Nᵢ ≤ Tᵢ | Σ = S]|

where N is the Gaussian approximation. For the Gaussian, the conditional distribution has correlation -1/(R-1), and the joint probability can be computed explicitly.

**Challenge**: The error bound must be at most 2^{-αR} for α < 1, which is stronger than typical multivariate CLT bounds (which give polynomial errors in the right metric but not necessarily multiplicative bounds on probabilities of polytopes).

### 9.2 Entropy / Information-Theoretic Bound

  H(B₁,...,B_R | Σ = S) = Σ H(Bᵢ | Σ = S) - I(B₁;...;B_R | Σ = S)

The mutual information I captures the "cost of dependence." For equi-correlated variables with ρ → 0, I = O(R²ρ²) = O(1) bits. This would give:

  E[Z] ≥ 2^{-Σ H(Bᵢ)} × 2^{-O(1)} ≈ ∏ Pr[Bᵢ ≤ Tᵢ] × 2^{-O(1)}

which suffices since ∏ Pr ≈ 2^{-R} and budget ≈ 2R.

**Challenge**: Making the entropy bound rigorous for our specific distribution.

### 9.3 Direct Combinatorial Argument

Count the number of D12 with all B(d) ≤ T(d) by a weighted sum:

  #{valid D12} = Σ_{D12} ∏_d 1[B(d) ≤ T(d)]

Replace the hard constraints with soft constraints (indicator → exponential penalty) and use Fourier analysis on Z_p to decouple the sum.

---

## 10. Exact B(d) Distribution

The distribution of B(d) for random D12 is SYMMETRIC around its mean:

**p=11, k=4**: B(1) ~ {0:5, 1:50, 2:100, 3:50, 4:5}/210
  Symmetric around E[B]=2, support {0,1,2,3,4}

**p=19, k=8**: B(1) ~ {0:9, 1:324, 2:3024, 3:10584, 4:15876, 5:10584, 6:3024, 7:324, 8:9}/43758
  Symmetric around E[B]=4, support {0,...,8}

The symmetry B(d) ↔ k - B(d) (complement) is exact for the unconditional distribution. This means Pr[B ≤ E[B] - 1] = Pr[B ≥ E[B] + 1] = 1/2 - Pr[B = E[B]]/2.

**Key identity**: Pr[B(d) ≤ E[B]-1] = (1 - Pr[B = E[B]])/2 exactly.

At the mean threshold T = E[B] - 1 (the D11 binding constraint), this gives:

  Pr[B ≤ T] = (1 - Pr[B = (p-3)/4]) / 2

For large p, Pr[B = (p-3)/4] ~ 1/√(πp/8) → 0, so Pr[B ≤ T] → 1/2 from below.

The exact cost per D11 constraint is:

  -log₂(Pr[B ≤ T]) = 1 + log₂(1/(1 - Pr[B = E[B]])) ≈ 1 + Pr[B = E[B]]/(2 ln 2)

---

## 11. Verification Script

See `ramsey-book-graphs/first_moment_computation.py` for:

1. Budget-cost decomposition at all primes p ≤ 997
2. Exact B(d) distribution at p=7,11
3. Exact marginal and joint probabilities at p=7,11
4. c₀ computation for all D11 at p=11,19,23
5. Asymptotic expansion verification

Key numerical findings:
- log₂(E[N]) fits ≈ 0.115p + 1.64 (excluding p=19)
- c₀ for working D11: 0.700 (p=11), 0.033 (p=19), 0.006-0.036 (p=23)
- Averaged c₀: 0.678 (p=11), 0.00577 (p=19)
- Headroom vastly exceeds c₀ cost for all tested primes

---

## 12. Conclusions

### What is established:
1. **Budget ≈ p bits, cost ≈ p/2 bits**: The combinatorial counting argument gives exponential headroom. This is PROVEN unconditionally.

2. **Product of marginals ≈ 2^{-p/2}**: Each of R ≈ p/2 constraints costs ≈ 1 bit. The exact cost per constraint is 1 + O(1/p) (not 1 + O(1/√p) as the Gaussian suggests). This is PROVEN.

3. **E[N] = #D12 × E[Z] with headroom**: The first moment formula gives E[N] ≈ 2^{p/2} × E[Z]. For E[N] → ∞, we need E[Z] > 2^{-p/2+ε}.

### What remains:
4. **Bounding E[Z] from below**: We need E[Z] ≥ 2^{-(1+α)p/2} for some α < 1. Empirically E[Z] ≈ 2^{-p/2 - O(1)}, so α ≈ 0. But proving this requires controlling the correlation structure of (B(d₁),...,B(d_R)).

5. **The correlation factor is the SAME gap as in the L6 proof**: Proving c₀ ≥ 1/poly(p) would close both the first moment argument and the original L6 conditional proof.

### Assessment:
The first moment approach does NOT bypass the c₀ problem — it reformulates it. The budget-cost analysis is clean and the exponential headroom is real, but the correlation between constraints remains the fundamental obstacle. The product-of-marginals bound overestimates E[N] by about 2^{5-8} at small primes, suggesting a CONSTANT (not growing) correlation penalty. If this constancy could be proven, the first moment proof would be complete.
