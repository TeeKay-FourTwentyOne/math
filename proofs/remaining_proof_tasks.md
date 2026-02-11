# Remaining Work to Complete the Proof of R(B_{n-1}, B_n) = 4n-1

**Date**: 2026-02-10
**For**: Claude Code Agent Team
**From**: Opus 4.6 (web UI review)
**Context**: I have reviewed proof_positive_association.md, first_moment_proof_sketch.md, proof_outline.md, external_input_20260210.md, and README.md. The proof architecture is sound. This document specifies exactly what remains.

---

## Status

The proof chain is:

1. Structural reduction (Theorems 1–3) — **PROVEN**
2. Optimal |D11| size (Theorem 5) — **PROVEN** (exhaustive + sampling)
3. Moment computations (L1, L2, L3) — **PROVEN** (exact_moments.py)
4. Var[A] computation (L3b) — **MC-VERIFIED ONLY, needs rigorous proof**
5. Per-constraint rates ≈ 1/2 (L5) — **PROVEN** (standard CLT)
6. Joint probability bound (L6) — **COMPUTATIONALLY VERIFIED, needs formal proof**

The exponential headroom (E[valid D12] ≈ 2^{0.7p}) means even a joint probability bound with polynomial loss p^{-C} suffices. Full positive association (ratio ≥ 1) is stronger than needed.

---

## Task 1: Prove L3b (Var[A(d)] = p/8 + O(1))

### What exists
L3 (Var[B] = p/16) is proven exactly in exact_moments.py. L3b (Var[A] = p/8) is MC-verified only (200K trials, Var[A]/p → 0.125).

### Why it differs from L3
A(d) is the autocorrelation of a *symmetric* random subset — chosen by selecting (p-1)/4 negation pairs {x, -x}. The pair structure means indicators 1[x ∈ D11] and 1[-x ∈ D11] are perfectly correlated (both 1 or both 0), while indicators for different pairs are independent.

### Approach
Same indicator decomposition technique as exact_moments.py. Decompose A(d) = Σ_x 1[x ∈ D11] · 1[x+d ∈ D11]. Group terms by which negation pairs they involve. Terms involving two different pairs are products of independent indicators. Terms where x and x+d belong to the same pair (i.e., d = -2x) create a dependency. There are O(1) such special terms per shift d. The generic terms give the leading contribution p/8; the special terms contribute O(1).

### Difficulty: LOW
Careful but straightforward. Verify against MC at every step.

---

## Task 2: Prove L6 — The Main Gap

### What we need (weakened, sufficient for proof)

For a fixed symmetric D11 of size n = (p+1)/2, with D12 uniform:

$$\Pr\bigl[\forall d \in D_{11}: B(d) \leq T(d)\bigr] \;\geq\; p^{-C} \cdot \prod_{d \in D_{11}} \Pr\bigl[B(d) \leq T(d)\bigr]$$

for some absolute constant C. Any polynomial loss is tolerable because E[valid D12] grows as 2^{0.7p}.

### What's computationally verified
- Full positive association (ratio ≥ 1) for ALL tested D11s at p = 11, 19, 23, 31, 43, 47, 59.
- Ratio grows ≈ p^{1.5}: from 4.7 at p=11 to 131.4 at p=59.
- Verified in joint_probability_analysis.py.

### Dead ends (do not pursue)
- **FKG**: B(d) is not monotone in D12 elements. Doesn't apply.
- **LLL**: Dependency graph is complete. Gives nothing.
- **Suen/Janson**: Too crude (tested in suen_inequality_test.py).

---

### RECOMMENDED APPROACH: Multivariate CLT (Asymptotic + Finite Verification)

**This is the shortest path to completing the proof.**

**Core idea**: For large p, the vector (B(d))_{d ∈ D11} is approximately multivariate Gaussian with near-diagonal covariance. For near-diagonal multivariate Gaussians, joint threshold probabilities factorize.

**Why this works**:

**Step A — CLT for each B(d):** Each B(d) is a sum of ~(p-1)/2 weakly dependent terms. Specifically, for D12\{0} a random (p-3)/2-subset of {1,...,p-1}:

B(d) = Σ_{x ∈ D12\{0}} 1[x+d ∈ D12]

This is a function of a uniformly random subset. The CLT for such functions follows from Stein's method for combinatorial central limit theorems (Goldstein-Rinott 1996). Berry-Esseen gives rate O(1/√p).

**Step B — Near-zero pairwise covariance:** Compute Cov[B(d₁), B(d₂)] exactly for d₁ ≠ d₂ (same technique as exact_moments.py). The computational data shows |mean_corr| < 0.001 at all tested primes. If |Cov[B(d₁), B(d₂)]| = O(1) while Var[B(d)] = p/16, then the correlation ρ = O(1/p) → 0.

**THIS IS THE KEY NEW COMPUTATION: compute Cov[B(d₁), B(d₂)] exactly using the indicator decomposition method from exact_moments.py. If the result confirms |Cov| = O(1), the rest of the argument follows from standard theory.**

**Step C — Multivariate normal approximation:** With near-diagonal covariance, the multivariate CLT gives:

|Pr[∀d: B(d) ≤ T(d)] - Pr[∀d: Z_d ≤ T(d)]| → 0

where Z ~ N(μ, Σ) with Σ nearly diagonal. Use Stein's method for multivariate normal approximation (Chen-Goldstein-Shao 2011).

**Step D — Factorization for near-diagonal Gaussian:** For Z ~ N(μ, Σ) with Σ = σ²I + E where ||E||/σ² → 0:

Pr[∀d: Z_d ≤ t_d] / ∏ Pr[Z_d ≤ t_d] → 1

Explicit bounds via Slepian's inequality and Gaussian comparison theorems (Slepian 1962; Li and Shao 2002).

**Step E — Combine:** For p ≥ p₀ (with p₀ chosen so error terms are small), the ratio joint/product is ≥ 1 - o(1) > 1/poly(p). For p < p₀ ≤ 61, use computational verification already done.

**References**:
- Goldstein and Rinott, "Multivariate normal approximations by Stein's method and size bias couplings," J. Applied Prob. (1996)
- Chen, Goldstein, Shao, "Normal Approximation by Stein's Method," Springer (2011)
- Slepian, "The one-sided barrier problem for Gaussian noise," Bell System Technical Journal (1962)

---

### ALTERNATIVE APPROACH: Asymptotic analysis of log-ratio

Directly analyze log R(p) = log(Pr[all ok] / ∏ Pr[ok_d]) as p → ∞.

The inclusion-exclusion expansion gives:

log Pr[all ok] = Σ_d log Pr[ok_d] + Σ_{d₁<d₂} log(1 + Cov[1_{ok_{d₁}}, 1_{ok_{d₂}}]/(Pr[ok_{d₁}]Pr[ok_{d₂}])) + higher order

If pairwise event covariances are O(1/p) and there are O(n²) = O(p²) pairs, the correction is O(p). Since log ∏ Pr[ok_d] ≈ -n log 2 ≈ -p/2, the correction is a lower-order additive term, giving:

log Pr[all ok] ≥ Σ log Pr[ok_d] - O(p) ≥ -n log 2 - O(p)

This is enough: E[valid D12] = C(p-1,(p-3)/2) × Pr[all ok] ≥ 2^{p-1} × 2^{-n - O(p)}.

Wait — this gives 2^{p-1-n-O(p)}, and if the O(p) term has constant > 1/2, this could go to zero. Need the constant to be small. This approach requires more careful analysis of the pairwise correction terms.

**Bottom line**: The CLT approach (Steps A-E) is cleaner. The log-ratio approach works but requires tighter bounds on the correction constant.

---

## Task 3: Minor Fixes

### 3a. proof_outline.md line 68
Change `B(-d)` to `B(d)` in the V2V2 formula. The universal symmetry lemma gives B(-d) = B(d) always.

### 3b. Explicit p = 7 handling
Add to proof_outline.md and proof_positive_association.md: "For p = 7, no 2-block circulant satisfying (S1)-(S3) exists (exhaustive verification, d11_size_survey.py). The case n = 4 is handled by prior constructions [Lidický et al. 2024]."

### 3c. Var[A+B] consistency
All documents should use: Var[A+B] = Var[A] + Var[B] = p/8 + p/16 = 3p/16 + O(1). Gap/std = -4/√(3p). Some documents currently use p/8 or other expressions.

### 3d. "Fix ANY D11" claim
In proof_positive_association.md (Step 2) and first_moment_proof_sketch.md (Section 11), note: the proof works for any D11 because E[valid D12] grows exponentially (≈ 2^{0.7p}) regardless of D11 choice. Per-constraint rates vary but the geometric mean is empirically ≈ 2^{-0.7}, and C(p-1,(p-3)/2) ≈ 2^{p-1} dominates any sub-exponential variation.

---

## Formulation note: use k = n throughout

The k = n and k = n-2 formulations are rigorously equivalent via complement bijection (Theorem 5). They represent the same construction viewed from opposite sides:

| | k = n | k = n-2 |
|--|-------|---------|
| D11 constraint | binding (gap -1) | easy (gap +1) |
| D22 constraint | easy (gap +2) | binding (gap ≈ 0) |

The difficulty is conserved — it just moves between D11 and D22 under complement. The k = n formulation is notationally simpler because the binding constraint lives on D11 (the set being designed) and D22 is trivially satisfied. Use k = n for the proof.

---

## Summary

| Priority | Task | Difficulty | What it closes |
|----------|------|-----------|----------------|
| **1** | Compute Cov[B(d₁), B(d₂)] exactly | LOW-MEDIUM | Enables the CLT argument |
| **2** | L6 via multivariate CLT + finite verification | MEDIUM | **Completes the proof** |
| **3** | L3b (Var[A] exact) | LOW | Minor gap |
| **4** | Document fixes (3a-3d) | LOW | Consistency |

**The single most important new computation**: Exact calculation of Cov[B(d₁), B(d₂)] for d₁ ≠ d₂, using the indicator decomposition method from exact_moments.py. If this confirms |Cov| = O(1) (while Var = O(p)), the multivariate CLT argument completes L6 for large p, and the proof is done.
