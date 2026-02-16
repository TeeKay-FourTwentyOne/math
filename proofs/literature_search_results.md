# Literature Search Results: Tools for Proving L6

**Date**: 2026-02-13
**Goal**: Find mathematical tools that could prove: for every sufficiently large prime p ≡ 3 mod 4, there exist (D11, D12) in Z_p satisfying all ~p Ramsey constraints simultaneously.

---

## Executive Summary

**Most promising tools** (ranked by applicability):

1. **Negative Association + Reverse FKG** — The uniform k-subset measure is negatively associated (Joag-Dev & Proschan 1983). Our constraints {B(d) ≤ T} are *decreasing* events in the D12 indicators. For NA measures, Pr[∩ decreasing events] ≥ ∏ Pr[decreasing event_i]. This gives **exactly the product lower bound** we need — and we already know the product exceeds the true probability (our data shows positive association). **This is the wrong direction** — it gives a *lower* bound that we already have numerically, but doesn't help prove the lower bound is positive.

2. **Kim-Vu Polynomial Concentration on the Slice** — B(d) is a degree-2 polynomial in D12 indicators. Kim-Vu (2000) + Filmus lifting (2016) give sub-Gaussian tail bounds for low-degree polynomials on slices. Could bound Pr[B(d) > T] for individual constraints, feeding into Janson or union bound.

3. **Janson Inequality (Extended)** — Classic Janson bounds Pr[no bad event] ≥ exp(-μ - Δ) where μ = Σ Pr[bad_i] and Δ = Σ_{i~j} Pr[bad_i ∩ bad_j]. In our setting, the "bad events" are {B(d) > T}. If μ → 0 and Δ = o(μ²), then Pr[all ok] → 1. The challenge: our events are NOT of the standard "subgraph containment" form that Janson was designed for.

4. **Stein's Method / Multivariate Normal Approximation** — The B(d) vector is approximately multivariate normal (by CLT on the slice). Reinert & Röllin (2009) give quantitative bounds. Could prove Pr[all B(d) ≤ T] ≈ Pr[all Gaussian ≤ T] and then compute the Gaussian probability.

5. **Entropy/Shearer's Lemma** — Potentially useful for counting arguments but seems indirect for our existence problem.

---

## 1. Hypercontractivity on the Johnson Graph / Boolean Slice

### Key Papers
- **Filmus, O'Donnell, Wu (2022)**: "Log-Sobolev inequality for the multislice" — Electronic Journal of Probability 27, 1-30. [arXiv:1809.03546](https://arxiv.org/abs/1809.03546)
- **Filmus (2016)**: "An orthogonal basis for functions over a slice of the Boolean hypercube" — Electronic Journal of Combinatorics 23(1). [arXiv:1406.0142](https://arxiv.org/abs/1406.0142)
- **Sambale & Sinulis (2022)**: "Concentration inequalities on the multislice and for sampling without replacement" — Journal of Theoretical Probability 35(4), 2712-2737. [arXiv:2010.16289](https://arxiv.org/abs/2010.16289)

### What These Give Us
- **Log-Sobolev constant**: For the slice {x ∈ {0,1}^p : Σx_i = k} with k = (p-1)/2, the log-Sobolev constant ρ satisfies ρ^{-1} ≤ Cp for a constant C. This implies hypercontractivity: ||T_t f||_q ≤ ||f||_2 for appropriate noise operator T_t and q > 2.
- **Filmus basis**: Explicit orthogonal basis for functions on the slice, analogous to Fourier characters on the cube. B(d) is a degree-2 function in this basis.
- **Concentration**: For Lipschitz functions on the slice, sub-Gaussian concentration with variance proxy ~1/p (from the finite-sampling correction factor 1 - n/N).
- **Talagrand convex distance**: Holds on the multislice (Sambale-Sinulis).

### Applicability Assessment
**(a) Does it apply?** YES — our D12 is a uniform k-subset of Z_p, exactly the slice setting. B(d) is a degree-2 polynomial in the D12 indicator variables.

**(b) What would it give us?**
- Individual tail bounds: Pr[B(d) > T] via hypercontractivity for degree-2 functions.
- Potentially: a "level-d inequality" bounding Pr[∩ {B(d) ≤ T}] using spectral information.

**(c) Obstacles:**
- Hypercontractivity gives *concentration* (tail bounds for a single function), not directly *joint* tail bounds for ~p correlated functions.
- The standard approach would bound Pr[B(d) > T] for each d, then use union bound — but union bound over ~p events may lose too much.
- Need to combine with a method that handles the joint structure (Janson, second moment, etc.).

**VERDICT: Useful building block (individual tail bounds), but not sufficient alone.**

---

## 2. Positive Association and FKG-type Inequalities

### Key Papers
- **Joag-Dev & Proschan (1983)**: "Negative association of random variables, with applications" — Annals of Statistics 11(1), 286-295.
- **Borcea, Brändén, Liggett (2009)**: "Negative dependence and the geometry of polynomials" — JAMS 22(2), 521-567. [arXiv:0707.2340](https://arxiv.org/abs/0707.2340)
- **Gladkov (2024)**: "A strong FKG inequality for multiple events" — Bulletin of the London Mathematical Society 56, 2794-2801. [arXiv:2305.02653](https://arxiv.org/abs/2305.02653)
- **Pemantle (2000)**: "Towards a theory of negative dependence" — J. Math. Phys. 41, 1371-1390. [arXiv:math/0404095](https://arxiv.org/abs/math/0404095)
- **Dubhashi & Ranjan (1998)**: "Negative dependence through the FKG inequality" — BRICS Report RS-96-27.

### The Key Structural Fact

**The uniform k-subset measure is negatively associated (NA).** This is a classical result (Joag-Dev & Proschan 1983): sampling without replacement gives NA random variables.

For NA random variables: if f_1, ..., f_m are *nondecreasing* functions of disjoint subsets of variables, then Cov(f_i, f_j) ≤ 0.

**Critical observation for our problem:**

Our constraint events are {B(d) ≤ T} where B(d) = |{(i,j) : i,j ∈ D12, i-j ≡ d mod p}|.

- B(d) is an *increasing* function of D12 (adding an element can only increase or maintain B(d)).
- So {B(d) ≤ T} is a *decreasing* event.
- For NA measures, decreasing events are *positively* correlated: Pr[∩_d {B(d) ≤ T}] ≥ ∏_d Pr[B(d) ≤ T].

**This confirms our numerical observation that Pr[all ok] > ∏ Pr[ok_d]** (the actual probability exceeds the independence prediction).

### Borcea-Brändén-Liggett: Strong Rayleigh

The uniform k-subset measure is **strongly Rayleigh** (the strongest form of negative dependence). This implies:
- Negative association
- Ultra-log-concavity of marginals
- Concentration inequalities (Pemantle-Peres 2014)
- All previously studied notions of negative dependence

### Gladkov 2024: Strong FKG for Multiple Events

Gladkov's result gives a strong FKG inequality for multiple events under product measures. However, our measure (uniform k-subset) is NOT a product measure, so this doesn't directly apply. The strong Rayleigh property is the correct framework.

### Applicability Assessment
**(a) Does it apply?** YES — uniform k-subset is NA/strongly Rayleigh, and our events are decreasing.

**(b) What would it give us?**
- Pr[all B(d) ≤ T] ≥ ∏ Pr[B(d) ≤ T] — a **product lower bound**.
- If we can show ∏ Pr[B(d) ≤ T] > 0 (equivalently, each Pr[B(d) ≤ T] > 0, which is trivially true since T ≥ 0), this gives existence.
- But the product bound is trivially positive! The real question is whether it grows fast enough.

**(c) Obstacles:**
- The product bound ∏ Pr[B(d) ≤ T] can be exponentially small (each factor < 1, ~p factors).
- We need to show this product is positive, which is trivially true, OR we need E[N] > 0 which requires ∏ Pr > 0 over all D11.
- **The NA lower bound is nontrivial only if we need to go from "probably ≥ 0" to "definitely > 0".** Since E[N] = Σ_{D11} ∏ Pr[B(d) ≤ T | D11] × (binomial coefficient normalization), and each term is positive, E[N] > 0 is automatic. The real question is whether E[N] → ∞.

**VERDICT: Confirms positive association but doesn't directly close the gap. The NA/product lower bound is the FLOOR, not the ceiling. We need E[N] → ∞, which is a first-moment computation, not a correlation inequality.**

---

## 3. Anti-concentration of Quadratic Forms on Subsets

### Key Papers
- **Costello, Tao, Vu (2006)**: "Random symmetric matrices are almost surely non-singular" — Duke Math J. 135(2), 395-413. [arXiv:math/0505156](https://arxiv.org/abs/math/0505156)
- **Kwan & Sudakov (2020)**: "Anticoncentration for subgraph statistics" — Annals of Probability 48(3), 1371-1405. [arXiv:1807.05202](https://arxiv.org/abs/1807.05202)
- **Fox, Kwan, Sudakov (2020)**: "Anticoncentration for subgraph counts in random graphs" — [arXiv:1905.12749](https://arxiv.org/abs/1905.12749)
- **Fox, Kwan, Sauermann (2020)**: "Combinatorial anti-concentration inequalities, with applications" — Mathematical Proceedings of the Cambridge Philosophical Society. [arXiv:1905.12142](https://arxiv.org/abs/1905.12142)
- **Kim & Vu (2000)**: "Concentration of multivariate polynomials and its applications" — Combinatorica 20(3), 417-434.

### Relevance

B(d) = Σ_{i,j ∈ D12, i-j≡d} 1 is a degree-2 polynomial (U-statistic) in the D12 indicators. The anti-concentration question asks: how spread out is B(d)?

**Costello-Tao-Vu** developed a quadratic Littlewood-Offord theory for random symmetric matrices. Their result: for random ±1 variables ξ_i, the quadratic form Q = Σ a_{ij} ξ_i ξ_j satisfies Pr[Q = v] ≤ C n^{-1/8} for any fixed v, under non-degeneracy conditions. Recent work (2023) achieves the optimal n^{-1/2} bound.

**Kwan-Sudakov** prove anti-concentration for subgraph statistics. For a graph G on n vertices, the number of edges spanned by a random k-subset is anti-concentrated: Pr[X = ℓ] ≤ O(k^{-1/2} log^O(1) k).

### Applicability Assessment
**(a) Does it apply?** PARTIALLY. B(d) IS a degree-2 polynomial in the D12 indicator variables on the slice. However:
- Costello-Tao-Vu works for independent ±1 variables, not the slice (k-subset) measure.
- Kwan-Sudakov works for k-subsets but gives *anti-concentration* (upper bound on point probabilities), not *concentration* in our desired range.

**(b) What would it give us?**
- Anti-concentration: Pr[B(d) = v] = O(1/√k) for any fixed v. This is useful for showing B(d) is spread out.
- But we need Pr[B(d) ≤ T] for our specific threshold T, which is a *tail bound*, not an anti-concentration bound.
- Could be useful combined with normal approximation: if B(d) is approximately normal with known mean and variance, anti-concentration + CLT gives Pr[B(d) ≤ T].

**(c) Obstacles:**
- Doesn't directly give joint bounds.
- Tail bounds on the slice for degree-2 polynomials need the Filmus-Sambale framework (Section 1).

**VERDICT: Useful for understanding individual B(d) distribution, but not directly applicable to the joint existence problem.**

---

## 4. Janson Inequality / Extended Janson

### Key Papers
- **Janson (1990)**: "Poisson approximation for large deviations" — Random Structures & Algorithms 1, 221-229.
- **Suen (1990)**: "A correlation inequality and a Poisson limit theorem for noninteracting random variables" — Annals of Probability 18, 1305-1315.
- **Janson, Luczak, Ruciński (2000)**: Book "Random Graphs" — Chapter on Janson inequalities.
- **Janson (2016)**: "The lower tail: Poisson approximation revisited" — Random Structures & Algorithms 48(2), 219-246. [arXiv:1406.1248](https://arxiv.org/abs/1406.1248)
- **Scott (2012)**: "The Janson inequalities for general up-sets" — [arXiv:1203.1024](https://arxiv.org/abs/1203.1024)

### Standard Janson Setup
Consider independent random variables, events A_1, ..., A_m where A_i depends on a subset S_i of the variables. Define:
- μ = Σ Pr[A_i]
- Δ = Σ_{i~j} Pr[A_i ∩ A_j] (sum over dependent pairs i,j where S_i ∩ S_j ≠ ∅)

Then: **Pr[∩ Ā_i] ≥ exp(-μ - Δ)**

If Δ = o(μ), this gives Pr[∩ Ā_i] ≥ exp(-μ(1+o(1))).

### Application to Our Problem

Define bad events A_d = {B(d) > T} for each d in our constraint set. Then:
- We want Pr[∩ Ā_d] = Pr[all B(d) ≤ T] > 0.
- μ = Σ Pr[B(d) > T] — the expected number of violated constraints.
- Δ = Σ_{d~d'} Pr[B(d) > T, B(d') > T].

**If μ → 0** (each constraint has vanishing violation probability), and Δ = o(μ), then Pr[all ok] → 1 for that D11.

### Critical Obstacle: Non-independence of Base Variables

**Janson's inequality is designed for INDEPENDENT base random variables.** Our D12 is drawn from the uniform k-subset measure, which has dependent indicators.

**Possible workaround**:
1. Use the "Poissonized" version: sample each element independently with probability k/p = (p-1)/(2p) ≈ 1/2, condition on the total being exactly k. The conditioning introduces dependencies but the "pre-conditioned" measure is a product measure. Janson + de-Poissonization might work.
2. Use Scott (2012): Janson inequalities for general up-sets, which extends to non-product settings using positive correlation.
3. Use the strongly Rayleigh property to transfer Janson-type bounds from the product to the slice.

### Applicability Assessment
**(a) Does it apply?** POTENTIALLY, with modifications. The standard Janson requires independent base variables, but:
- Scott's extension works for up-sets (our bad events are UP-sets since B(d) is increasing in D12).
- De-Poissonization is standard in probabilistic combinatorics.

**(b) What would it give us?**
- If μ → 0 for appropriate D11: Pr[all B(d) ≤ T | D11] → 1, hence N(D11) → C(p, (p-1)/2).
- This is MUCH stronger than we need.
- More realistically: if μ = O(log p), then Pr[all ok] ≥ exp(-O(log p)) = poly(1/p), giving E[N] ≥ poly(p) × C(p,k) × poly(1/p) = poly growth.

**(c) Obstacles:**
- We need to compute μ = Σ Pr[B(d) > T | D11] — requires tail bounds for B(d) on the slice.
- We need to compute Δ — requires joint tail bounds.
- The dependency structure is FULL (every pair (d,d') has correlated B-values), so Δ involves all pairs.
- The equi-correlation Cov[B(d),B(d')] = -(p+1)/(8(p-2)) is NEGATIVE, which is favorable (Scott's extension uses positive correlation of Ā_d events).

**VERDICT: PROMISING. The Janson approach via Scott's extension or de-Poissonization could work if we can compute μ and show Δ/μ → 0. The negative correlation between B(d) values is FAVORABLE. This is a leading candidate.**

---

## 5. Talagrand's Inequality on the Slice

### Key Papers
- **Sambale & Sinulis (2022)**: Already cited above. Proves Talagrand's convex distance inequality for the multislice.
- **Talagrand (1995)**: "Concentration of measure and isoperimetric inequalities in product spaces" — Publications Mathématiques de l'IHÉS 81, 73-205.

### What It Gives
Talagrand's convex distance inequality on the slice: for a set A in the slice with measure at least 1/2,
Pr[d_T(x, A) ≥ t] ≤ C exp(-t²/4)
where d_T is the Talagrand (convex) distance.

### Applicability Assessment
**(a) Does it apply?** The slice setting matches perfectly.

**(b) What would it give?** Concentration of Lipschitz (or convex) functions. But our question is about the probability that ALL ~p constraints are satisfied simultaneously — this is about a complex event, not a single function being concentrated.

**(c) Obstacles:** Talagrand concentration is for a *single* function. We need a *joint* statement about ~p degree-2 polynomials.

**VERDICT: Useful supporting tool (e.g., to prove N is concentrated around E[N]), but not directly applicable to the existence question.**

---

## 6. Recent Results on R(B_n, B_m) and Diagonal Ramsey

### Key Papers
- **Campos, Griffiths, Morris, Sahasrabudhe (2023)**: "An exponential improvement for diagonal Ramsey" — Annals of Mathematics (to appear). [arXiv:2303.09521](https://arxiv.org/abs/2303.09521)
  - Proved R(k) ≤ (4-ε)^k, first exponential improvement since Erdős-Szekeres (1935).

- **Balister, Bollobás, et al. (2024)**: "Upper bounds for multicolour Ramsey numbers" — [arXiv:2410.17197](https://arxiv.org/abs/2410.17197)
  - Extended to r-colour case: R_r(k) ≤ e^{-δk} r^{rk}.

- **Book-cycle Ramsey** (2024-2025): Exact values of R(B_n^{(2)}, C_m) for large n. Not directly relevant.

### Relevance to Our Problem
The Campos et al. breakthrough uses a new "book" algorithm with a probabilistic argument about graph colorings. Their key lemma involves controlling the codegree structure in random graphs. However:
- Their method is for the UPPER bound on diagonal Ramsey numbers.
- Our problem is a LOWER bound (constructing a 2-coloring) for off-diagonal book Ramsey numbers.
- The techniques don't directly transfer.

**VERDICT: Inspirational but not directly applicable. Different problem direction (upper vs lower bounds).**

---

## 7. Supplementary Difference Sets in Z_p

### Key Papers
- **Wallis (1972)**: Supplementary difference sets constructions.
- **Lehmer (1974)**: Cyclotomic difference sets.
- **Sun (2025)**: "Cyclotomic matrices and power difference sets" — [arXiv:2511.13613](https://arxiv.org/abs/2511.13613)
  - Studies cyclotomic matrices and their spectral/determinant properties.

### Connection to Our Problem
Our (D11, D12) pair resembles a 2-{p; (p+1)/2, (p-1)/2; λ} supplementary difference set. The spectral complementarity condition (verified numerically) says |hat(1_{D11})|² + |hat(1_{D12})|² ≈ constant — exactly the SDS condition.

### Known Results on SDS Existence
- For p ≡ 1 mod 4: Paley construction gives SDS (and our Ramsey construction).
- For p ≡ 3 mod 4: NO known cyclotomic SDS construction for the exact parameters we need.
- Golay pairs (the strongest form of SDS) DON'T exist for p ≡ 3 mod 4 (EKS theorem).
- The "near-SDS" condition (spectral complementarity with small error) may suffice, but existence for all p ≡ 3 mod 4 is open.

### Applicability Assessment
**(a) Does it apply?** Directly — our problem IS (essentially) constructing near-SDS pairs.

**(b) What would it give?** An algebraic construction guaranteeing existence.

**(c) Obstacles:** The exact SDS parameters we need are blocked by EKS for Golay pairs. Cyclotomic constructions work for p ≡ 1 mod 4 but not p ≡ 3 mod 4.

**VERDICT: Our problem IS an SDS problem in disguise, but the algebraic approach for p ≡ 3 mod 4 remains stuck. Probabilistic methods seem necessary.**

---

## 8. Entropy / Shearer's Lemma

### Key Papers
- **Galvin (2014)**: "Three tutorial lectures on entropy and counting"
- **Zhao (2022)**: MIT OCW lectures on entropy methods in combinatorics.
- **Radhakrishnan (2003)**: "Entropy and counting"

### Shearer's Lemma
For random variables X = (X_1, ..., X_n) and a covering family S of subsets of [n]:
k·H[X] ≤ Σ_{S ∈ S} H[X_S]
where k is the minimum coverage.

### Applicability Assessment
**(a) Does it apply?** Potentially, for counting valid (D11, D12) pairs.

**(b) What would it give?** An upper bound on the entropy of the valid pair distribution, which could give a lower bound on the number of valid pairs.

**(c) Obstacles:** Shearer's lemma typically gives UPPER bounds on counts (via entropy), while we need LOWER bounds. The standard entropy method would need to be applied in reverse, which is non-standard.

**VERDICT: Not directly useful. Entropy methods are better for upper bounds on counts.**

---

## 9. Stein's Method for Multivariate Normal Approximation

### Key Papers
- **Reinert & Röllin (2009)**: "Multivariate normal approximation with Stein's method of exchangeable pairs under a general linearity condition" — Annals of Probability 37(6), 2150-2173. [arXiv:0711.1082](https://arxiv.org/abs/0711.1082)
- **Chatterjee & Meckes (2008)**: "Multivariate normal approximation using exchangeable pairs"

### The Idea
The vector (B(d₁), ..., B(d_m)) should be approximately multivariate normal by CLT on the slice. Stein's method can quantify the approximation error:
||Law(B) - N(μ, Σ)||_TV ≤ error(p)

If error(p) → 0, then:
Pr[all B(d) ≤ T] ≈ Pr[all Gaussian ≤ T]

And the Gaussian probability can be computed from the known covariance structure.

### Applicability Assessment
**(a) Does it apply?** YES — the B(d) vector is a sum of weakly dependent terms (each pair (i,j) in D12 contributes to B(d) for d = i-j). The exchangeable pair construction is natural: swap one element of D12 with one outside.

**(b) What would it give us?**
- Quantitative CLT: distance between B-vector distribution and multivariate normal.
- Then: Pr[all B(d) ≤ T] ≈ Gaussian probability, which we can compute.
- If the Gaussian probability is bounded below by a positive quantity (for suitable D11), existence follows.

**(c) Obstacles:**
- Need the approximation error to be smaller than the target probability.
- The covariance matrix is nearly proportional to identity on the Parseval hyperplane — this is FAVORABLE (the Gaussian probability for equi-correlated normals with small correlation is well-studied).
- The number of constraints grows with p (~p constraints), so the target probability may be exponentially small, requiring exponentially good normal approximation.

**VERDICT: PROMISING as a quantitative tool. If we can show the multivariate CLT error is o(target probability), this closes the gap. The equi-correlated structure with near-zero correlation is very favorable for Gaussian approximation.**

---

## 10. Negative Dependence (Strong Rayleigh) + Decoupling

### Key Insight (Novel Combination)

The uniform k-subset measure is strongly Rayleigh. This gives us:
1. **Negative association**: decreasing events are positively correlated (Section 2).
2. **Concentration**: Lipschitz functions concentrate (Pemantle-Peres 2014).
3. **Stochastic domination**: The k-subset measure is stochastically dominated by independent Bernoulli(k/p) indicators (in a certain precise sense).

**Decoupling idea**:
- Compute Pr[all B(d) ≤ T] under the PRODUCT measure (independent Bernoulli with Pr[i ∈ D12] = k/p).
- Use strong Rayleigh properties to transfer bounds to the slice measure.
- The product measure makes the B(d)'s into sums of independent random variables, where standard Janson/concentration inequalities apply directly.

### Key Reference for Decoupling
- **de la Peña & Giné (1999)**: "Decoupling: From Dependence to Independence" — Springer.
- **Filmus lifting**: The Filmus basis provides explicit lifting from slice to cube functions, preserving L² norms. Could be used to "decouple" the slice constraint to a product-space constraint.

**VERDICT: This is potentially the most powerful combination. Use Filmus lifting to transfer the problem to the product space, apply standard Janson inequality there, then transfer back. The loss in the transfer should be controlled by the log-Sobolev constant.**

---

## Recommended Strategy

Based on this literature search, the most promising proof path for L6 combines:

### Path A: Janson via De-Poissonization
1. Work with product (Bernoulli 1/2) measure.
2. Compute μ = Σ_d Pr_product[B(d) > T] using standard tail bounds for sums of independent Bernoulli products.
3. Compute Δ = Σ_{d,d'} Pr_product[B(d) > T ∧ B(d') > T].
4. Apply Janson: Pr_product[all B(d) ≤ T] ≥ exp(-μ - Δ).
5. De-Poissonize: transfer back to uniform k-subset measure.
6. The de-Poissonization costs at most a factor of O(√p) (standard).

### Path B: Multivariate Normal Approximation
1. Use Stein's method (exchangeable pairs on the slice) to prove quantitative multivariate CLT for (B(d))_d.
2. Compute Pr[all Gaussian(μ,Σ) ≤ T] using the known covariance structure.
3. Show the approximation error is smaller than the Gaussian probability.

### Path C: First Moment + Second Moment (Current Approach, Enhanced)
1. Prove E[N] → ∞ (first moment: expected number of valid D12 over all D11).
2. Prove E[N²]/E[N]² = O(poly(p)) (second moment: bounded variance ratio).
3. Conclude N > 0 with positive probability by Paley-Zygmund.
4. **Enhancement**: Use negative association to tighten the second moment bound.
   - Since {B(d) ≤ T} events are positively correlated (NA of base measure), the second moment E[N²] is bounded above by... well, positive correlation means E[N²] ≥ E[N]² (we get a LOWER bound on E[N²], which is the wrong direction for Paley-Zygmund).
   - **Wait**: positive association gives E[N²] ≥ E[N]², which means the second moment ratio E[N²]/E[N]² ≥ 1 — this is automatic and unhelpful.
   - For Paley-Zygmund, we need E[N²]/E[N]² to be bounded (not too large). The NA property doesn't help with this.

### Conclusion

**Path A (Janson via de-Poissonization) is the most promising new approach.** It avoids the second moment entirely and goes directly to an exponential lower bound on Pr[all constraints satisfied | D11]. The key computation needed is μ (expected number of violated constraints under the product measure for appropriate D11).

**Path B (Stein's method)** is a strong alternative if the quantitative CLT can be made tight enough.

**Path C (enhanced first+second moment)** remains viable if E[N] → ∞ can be established; the literature search hasn't found a magic tool for the second moment ratio.
