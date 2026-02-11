# L6 Recovery: Task Specification After Slepian Sign Error

**Date**: 2026-02-10
**For**: Claude Code Agent Team
**Priority**: CRITICAL — this is the only remaining gap in the proof

---

## What Happened

The proof of L6 via multivariate CLT + Slepian's inequality is **invalid**. The sign is backwards.

### The error
Slepian's inequality states: if Cov(Zᵢ, Zⱼ) ≤ 0 for all i ≠ j, then

  Pr[all Zᵢ ≤ tᵢ] **≤** ∏ Pr[Zᵢ ≤ tᵢ]

NOT ≥. Negative correlation makes "all below threshold" events LESS likely than independence, not more. The proof claimed the opposite.

### Concrete verification
Z₁, Z₂ standard normal, ρ = -0.5, threshold t = 0:
- Pr[Z₁ ≤ 0, Z₂ ≤ 0] = 1/4 + arcsin(-0.5)/(2π) = 1/6
- ∏ Pr[Zᵢ ≤ 0] = 1/4
- Ratio = 2/3 < 1. Negative association, not positive.

### What this means
The Gaussian approximation to (B(d₁),...,B(dᵣ)) has ρ = -2/(p-3) < 0. In the Gaussian limit, the joint probability is LESS than the product of marginals. But the computational evidence shows ratios of 4.7× to 131× — massive positive association. Therefore:

**The observed positive association is a non-Gaussian phenomenon.** It cannot be captured by any argument that goes through pairwise covariances and a Gaussian approximation. The higher-order structure of the discrete k-subset measure creates the positive association, and the two-moment (Gaussian) approximation destroys it.

---

## What Still Holds

Everything else in the proof is correct:

- ✅ Structural reduction (Theorems 1-3)
- ✅ Optimal |D11| size (Theorem 5)
- ✅ Moment computations: E[B] = (p-3)/4, Var[B] = p/16 + O(1)
- ✅ Equi-covariance: Cov[B(d₁), B(d₂)] = -(p+1)/(8(p-2)) for non-complementary pairs
- ✅ Equi-covariance proof works for ALL odd m, not just primes (no primality used)
- ✅ Computational verification of positive association for p = 11, 19, 23, 31, 43, 47, 59
- ✅ E[valid D12] ≈ 2^{0.7p} (the first moment count is correct)
- ✅ If L6 holds (even with polynomial loss), the proof is complete

---

## The Mechanism Behind Positive Association

Understanding WHY positive association holds despite negative pairwise correlation is essential for proving it.

### The hard sum constraint
Σ_{d=1}^{m-1} B(d) = |D12|(|D12|-1) = EXACT CONSTANT

This is not a soft/approximate constraint — it's an algebraic identity. Every D12 satisfies it exactly.

### How this creates positive association
Consider the D11 events {B(d) ≤ T(d)} for d ∈ D11.

When several B(d) values for d ∈ D11 are simultaneously low (below threshold), the hard sum constraint forces the B(d) values for d ∈ D22 to be high. This "energy redistribution" pushes excess to D22 positions, where the threshold is LOOSER (n+1 vs n-2 at k=n). So the configurations where D11 events are jointly satisfied are exactly the configurations where energy has been successfully redistributed to D22 — and the asymmetric threshold structure means there's room for this.

In contrast, the Gaussian model treats the sum constraint as a soft effect (captured by the covariance matrix), losing the hard-constraint structure. The Gaussian predicts that when one B(d) is low, others are slightly more likely to be high (negative correlation) — but it misses that the TOTAL redistribution to D22 is what matters, and the D22 threshold has room.

### Key insight
The positive association is NOT between the random variables B(d). It's between the EVENTS {B(d) ≤ T(d)}. These are different things. Variables can be negatively correlated while the corresponding threshold events are positively associated, when there is a latent structure (the sum constraint + asymmetric thresholds) that couples the events favorably.

---

## Viable Approaches to Proving L6

### Approach 1: Conditional on the sum constraint (MOST PROMISING)

**Idea**: Condition on the total energy at D11 positions: S₁ = Σ_{d ∈ D11} B(d).

When S₁ is low (say S₁ ≤ (n-2) × n, the maximum if all D11 positions are at threshold), then there's a good chance ALL individual B(d) ≤ n-2 (since the average is at most n-2, and individual values concentrate around their conditional mean).

When S₁ is high, the event fails.

So: Pr[all D11 ok] ≈ Pr[S₁ ≤ target] × Pr[all individual ≤ n-2 | S₁ ≤ target].

The first factor involves a SINGLE random variable S₁ (not a joint probability over many events). The second factor is a max-vs-average bound conditioned on the sum.

**Concrete steps**:
1. Compute E[S₁] and Var[S₁] exactly. S₁ = Σ_{d ∈ D11} B(d). Using equi-covariance:
   - E[S₁] = n × E[B] = n(p-3)/4
   - Var[S₁] = n × Var[B] + n(n-1) × Cov_non-comp (all D11 pairs are non-complementary since D11 is symmetric: if d ∈ D11 then p-d ∈ D11, so d and p-d are both in D11, and their B values are identical — they don't contribute independent constraints)

   WAIT — this needs care. The n constraints for D11 include complementary pairs. Since D11 is symmetric, d ∈ D11 implies p-d ∈ D11. So the constraints come in complementary pairs, and B(d) = B(p-d) identically. There are only r = n/2 independent constraints (approximately). The sum S₁ double-counts: S₁ = 2 × Σ_{independent d in D11} B(d) (approximately).

   This needs to be worked out carefully.

2. Bound Pr[S₁ ≤ n(n-2)] using the mean and variance of S₁.

3. Conditioned on S₁ = s ≤ n(n-2), bound Pr[max_{d ∈ D11} B(d) > T(d)].
   The conditional distribution of (B(d))_{d ∈ D11} given S₁ = s has mean s/n ≤ n-2 per coordinate.
   Need to bound how much any single B(d) can deviate above its conditional mean.
   The equi-covariance structure means the conditional variance of each B(d) given S₁ is
   Var[B] × (1 - n × Cov² / (Var × Var_S₁)) or similar — work this out.

**Task for the team**: Formalize this conditioning argument. Compute Var[S₁], then bound the conditional max.

### Approach 2: Direct second moment on valid D12 count

**Idea**: For fixed D11, let N = #{valid D12}. We know E[N] ≈ 2^{0.7p}. Compute E[N²] and apply Paley-Zygmund: Pr[N > 0] ≥ E[N]²/E[N²].

E[N²] = Σ_{D12, D12'} Pr[D12 valid AND D12' valid].

For two independent D12, D12' drawn uniformly: Pr[both valid] = Pr[D12 valid]² (independence). So E[N²] = E[N]² when averaging over independent pairs. But N² counts ORDERED pairs from the SAME probability space. The correlation comes from overlap.

**Concretely**: E[N²] = Σ_s C(|D12|, s) × C(m - |D12|, |D12| - s) × Pr[both valid | intersection size s] where s = |D12 ∩ D12'|.

When s is close to |D12| (high overlap), the two events are nearly identical: Pr[both valid] ≈ Pr[one valid].
When s is close to the expected intersection of two random sets, the events are nearly independent.

**Task for the team**: For p = 11 and p = 19, compute E[N²] exactly:
1. Enumerate all valid D12 sets for a fixed good D11.
2. For each pair of valid D12 sets, record their intersection size.
3. Compute E[N²] = (# valid pairs) and E[N]² = (# valid)².
4. Report E[N²]/E[N]² — if this is polynomial in p, the proof works.

Do this for multiple D11 choices to check robustness.

### Approach 3: Entropy / information-theoretic bound

**Idea**: Use the conditional entropy chain rule. The events {B(d) ≤ T(d)} are determined by D12. Reveal D12 one element at a time. Each element changes each B(d) by at most 2. Use a martingale or entropy argument to bound the joint probability from below.

**Specifically**: Let x₁, ..., x_k be the non-zero elements of D12 revealed sequentially. After revealing x₁, ..., xᵢ, the conditional probability of "all ok" changes by a bounded amount (each new element affects each B(d) by at most 2, and there are n constraints).

If the conditional probability never drops too fast, the final probability is bounded below.

**Task for the team**: Implement and test the following for p = 11, 19:
1. For a fixed good D11, enumerate all valid D12.
2. Fix an ordering of elements.
3. Track Pr[all ok | first i elements revealed] as i goes from 0 to k.
4. Measure the maximum single-step drop in log-probability.
5. If the max drop is O(1/p) per step, then after k = O(p) steps, the total drop is O(1), giving Pr[all ok] ≥ Ω(1) × E[Pr[all ok]] — which is enough.

### Approach 4: Numerical/asymptotic study of the positive association ratio

**Idea**: Even without proving L6 analytically, understand the mechanism well enough to identify a provable proxy.

**Task for the team**:
1. For each tested prime, compute the joint probability Pr[all ok] and the product ∏ Pr[ok_d] to high precision.
2. Compute the LOG of the ratio: log(Pr[all ok]) - Σ log(Pr[ok_d]).
3. Decompose this into pairwise and higher-order contributions using inclusion-exclusion or cumulant expansion.
4. Determine which order of interaction dominates the positive association.
5. If the positive association is driven primarily by a single low-order effect (e.g., the sum constraint), this suggests which proof technique will work.

---

## Priority Order

1. **Approach 1** (conditioning on sum): Most likely to yield a clean proof. Start here.
2. **Approach 2** (second moment on N): Computational first, may give the proof directly. Run the E[N²] computation for p = 11, 19 immediately.
3. **Approach 4** (numerical decomposition): Will inform which approach is correct. Run in parallel.
4. **Approach 3** (entropy/martingale): Backup if 1 and 2 fail.

---

## What Success Looks Like

Any of the following would complete the proof:

**(a)** Prove Pr[all D11 ok | D11] ≥ 2^{-(p-1)} for any fixed D11.
   (This is the threshold for E[valid D12] ≥ 1.)

**(b)** Prove Pr[all D11 ok | D11] ≥ p^{-C} × ∏ Pr[ok_d] for some constant C.
   (Polynomial loss version — sufficient due to exponential headroom.)

**(c)** Prove E_{D11}[#{valid D12 for D11}] ≥ 1 by any method.
   (Averaging over D11 is also fine — we just need existence of ONE valid pair.)

**(d)** Prove E[N²]/E[N]² = O(poly(p)) for fixed D11.
   (Paley-Zygmund then gives Pr[N > 0] → 1.)

The bar is LOW because E[valid D12] ≈ 2^{0.7p}. We have massive exponential headroom. We just need to avoid losing more than 0.7p bits in the joint probability bound.
