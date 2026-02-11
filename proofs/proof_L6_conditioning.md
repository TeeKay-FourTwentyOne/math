# Proof of Lemma L6 via Parseval Hyperplane

**Date**: 2026-02-10
**Status**: Incomplete — correlation bound c_0 is OPEN (revision 6: Slepian direction corrected)
**Supersedes**: proof_L6_multivariate_clt.md (Slepian sign error), revisions 1-5 (Slepian direction error in all)

---

## 1. Statement

**Definition (A-flat D11).** A symmetric subset D11 of {1,...,p-1} with |D11| = n is called *A-flat with parameter C* if:

  max_{d: 1 <= d <= p-1} A(d) <= (p+1)/4 + C

where A(d) = Delta(D11, D11, d) is the autocorrelation of D11. Note: the max is over ALL nonzero d (both D11 and D22 positions), so A-flatness is a GLOBAL condition.

**Remark.** A-flatness with C = O(sqrt(p log p)) implies:
- T_red(d) = (p-3)/2 - A(d) >= (p-7)/4 - C = Theta(p) > 0 for all d in D11.
- T_blue(d) = (p+3)/2 - A(d) >= (p+5)/4 - C = Theta(p) > 0 for all d in D22.
- All thresholds T(d) deviate from their means by at most C, ensuring no constraint is pathologically tight.

**Lemma L6 (Parseval version).** Let p >= 19 be a prime with p = 3 (mod 4). Let D11 be an A-flat symmetric subset of {1,...,p-1} with |D11| = n = (p+1)/2 and parameter C = O(sqrt(p log p)). Let D12 = {0} union S where S is a uniformly random ((p-3)/2)-subset of {1,...,p-1}. Then

  E[#{valid D12} | D11] -> infinity as p -> infinity.

More precisely, conditional on a polynomial lower bound for the correlation loss (see Section 4.4):

  Pr[D12 is valid for D11] >= p^{-C} * prod Pr[B(d_i) <= T_i] >= p^{-C} * 2^{-(p-1)/2 + O(1)}

where the product of marginals comes from R = (p-1)/2 constraints each satisfied with probability ~1/2.

  E[#{valid D12}] = C(p-1, (p-3)/2) * Pr[valid]
                  >= 2^{p-1-(1/2)log_2(p)+O(1)} * p^{-C} * 2^{-(p-1)/2 + O(1)}
                  = 2^{(p-1)/2 - O(log p)} -> infinity.

**Note (revision 6):** The Slepian/Gaussian comparison inequality goes in the WRONG direction for this problem (negative correlation makes the Gaussian joint LESS than the product, not more). The product-of-marginals bound works for the EXACT discrete distribution (verified computationally), but the proof that c_0 >= p^{-C} for the k-subset measure remains open. See Section 4 for details.

**Lemma L6a (Existence of A-flat D11).** For all primes p = 3 (mod 4) with p >= 7, there exists an A-flat D11 of size n = (p+1)/2 with parameter C = O(sqrt(p log p)). A uniformly random symmetric D11 is A-flat with probability 1 - O(1/p).

*Proof of L6a.* For a random symmetric D11, each A(d) has E[A(d)] = (p+1)/4 and sub-Gaussian tails with parameter O(sqrt(p)) (hypergeometric concentration). By union bound over p-1 nonzero positions:

  Pr[max_d |A(d) - E[A]| > t] <= 2(p-1) * exp(-ct^2/p)

Setting t = sqrt(2p log(p)/c) gives probability <= 2(p-1)/p^2 < 1 for p >= 7. QED

**Remark on D11 selection.** Exhaustive enumeration at p=11,19,23 shows that only a fraction of symmetric D11 admit valid D12 (50%, 7.1%, 11.9% respectively). Working D11 are characterized by flat autocorrelation (max A(d) near E[A]) and balanced QR/QNR split. A-flatness is the precise condition needed: it ensures all thresholds T(d) are Theta(p), making the conditional z-scores o(1). Non-A-flat D11 with max A(d) much larger than E[A] have T_red(d_min) near zero, making those constraints too tight (Section 7.3).

---

## 2. Full Constraint Structure

### 2.1 Three constraint classes

A valid construction requires, for all nonzero d:

(R1) **V1V1 red** (d in D11): A(d) + B(d) <= (p-3)/2, i.e., B(d) <= T_red(d) = (p-3)/2 - A(d).

(R2) **V1V1 blue** (d in D22): A(d) + B(d) <= (p+3)/2, i.e., B(d) <= T_blue(d) = (p+3)/2 - A(d).

(R3) **V2V2 red** (d in D22): A(d) + B(p-d) <= (p+3)/2.
Since A(d) = A(p-d) (D11 symmetric) and d in D22 implies p-d in D22 (D22 symmetric), this is B(p-d) <= T_blue(p-d). So (R3) at d equals (R2) at p-d.

**Therefore (R3) is redundant with (R2).** The full constraint set is:

  For d = 1,...,p-1: B(d) <= T(d)

where T(d) = T_red(d) if d in D11, T(d) = T_blue(d) if d in D22.

### 2.2 Complementary pair reduction

Since D11 is symmetric, both D11 and D22 decompose into complementary pairs {d, p-d}. B(d) = B(p-d) and T(d) = T(p-d), so the constraints reduce to:

  r = (p+1)/4 constraints from D11 representatives (one per pair)
  r' = (p-3)/4 constraints from D22 representatives

Total: R = r + r' = (p-1)/2 independent constraints on the representative B(d_i) values.

### 2.3 Exact moments

- E[B(d)] = (p-3)/4 for all nonzero d.
- Var[B(d)] = (p-3)(p+1)/(16(p-2)).
- Cov[B(d_i), B(d_j)] = -(p+1)/(8(p-2)) for non-complementary pairs.
- Correlation rho = -2/(p-3) < 0.

### 2.4 The Parseval identity (key structural fact)

  sum_{d=1}^{p-1} B(d) = |D12| * (|D12| - 1) = ((p-1)/2) * ((p-3)/2) = (p-1)(p-3)/4

This is DETERMINISTIC: it depends only on |D12|, not on which subset. In terms of representatives:

  sum_{i=1}^{R} B(d_i) = (p-1)(p-3)/8

**This is an exact hard constraint.** It means the R = (p-1)/2 representative B-values live on an affine hyperplane. This hyperplane structure is what makes the proof work without any local CLT.

### 2.5 Equi-correlation structure

The covariance matrix of (B(d_1),...,B(d_R)) has equi-correlation: all pairs have Cov = -(p+1)/(8(p-2)). The eigenvalues are:

- lambda_0 = 0 (eigenvector = (1,...,1), corresponding to the Parseval sum)
- lambda_perp = (p+1)(p-1)/(16(p-2)) (multiplicity R-1)

---

## 3. Conditional Distribution on Parseval Hyperplane

### 3.1 Structure

The Parseval identity constrains the R representative B-values to the hyperplane sum = (p-1)(p-3)/8. This is NOT an approximation — it holds exactly for every D12.

On this hyperplane, the equi-correlation structure gives:
- Conditional mean of each B(d_i) = E[B] = (p-3)/4 (by symmetry of equi-correlated variables on a sum-hyperplane).
- Conditional variance = lambda_perp * (1 - 1/R) = lambda_perp * (R-1)/R.
- Conditional correlation between any two: -1/(R-1).

Define:
  sigma_cond = sqrt(lambda_perp * (R-1)/R) = sqrt((p+1)(p-1)/(16(p-2)) * (1 - 2/(p-1)))

For large p: sigma_cond ~ sqrt(p)/4.

### 3.2 Z-scores for each constraint

For the i-th representative constraint B(d_i) <= T(d_i):

  z_i = (T(d_i) - E[B]) / sigma_cond

For an A-flat D11 with A(d) ~ E[A] = (p+1)/4:

**D11 constraints** (r reps): T_red(d) - E[B] = ((p-3)/2 - A(d)) - (p-3)/4 = (p-3)/4 - A(d).
With A(d) ~ (p+1)/4: T_red - E[B] ~ (p-3)/4 - (p+1)/4 = -1.
So z_red ~ -1/sigma_cond ~ -4/sqrt(p) -> 0 from below.

**D22 constraints** (r' reps): T_blue(d) - E[B] = ((p+3)/2 - A(d)) - (p-3)/4 = (p+9)/4 - A(d).
With A(d) ~ (p+1)/4: T_blue - E[B] ~ (p+9)/4 - (p+1)/4 = 2.
So z_blue ~ 2/sigma_cond ~ 8/sqrt(p) -> 0 from above.

### 3.3 Product of conditional marginals

The product Phi(z_1) * ... * Phi(z_R) splits into D11 and D22 contributions:

  log2(prod) = r * log2(Phi(z_red)) + r' * log2(Phi(z_blue))

For A-flat D11 with z_red ~ -4/sqrt(p) and z_blue ~ +8/sqrt(p):
- Phi(z_red) = 1/2 - (4/sqrt(p)) * phi(0)/sqrt(2*pi) + ... ~ 1/2 - O(1/sqrt(p))
- Phi(z_blue) = 1/2 + (8/sqrt(p)) * phi(0)/sqrt(2*pi) + ... ~ 1/2 + O(1/sqrt(p))

Each marginal is close to 1/2 but not exactly 1/2. The D11 terms contribute slightly more than 1 bit each (z < 0), while D22 terms contribute slightly less than 1 bit each (z > 0). With r ~ r' ~ p/4:

  -log2(prod) ~ R * 1 + O(sqrt(p)) = (p-1)/2 + O(sqrt(p))

### 3.4 Why A-flatness is critical

**This is where A-flatness is essential.** For an A-flat D11 with parameter C:
- T_red(d_i) varies by at most C from its mean (p-7)/4, so z_red(d_i) = (-1 + O(C)) / sigma_cond = O(C/sqrt(p)).
- T_blue(e_j) varies by at most C from its mean (p+5)/4, so z_blue(e_j) = (2 + O(C)) / sigma_cond = O(C/sqrt(p)).
- With C = O(sqrt(p log p)), all z_i = O(sqrt(log p)) = o(1), so all Phi(z_i) are bounded away from 0.

For a NON-A-flat D11 (e.g., max A(d) = (p-3)/2 at some d): T_red(d) = 0, giving z = -(p-3)/(4*sigma_cond) = -Theta(sqrt(p)), so Phi(z) is exponentially small. This is why "for any D11" fails.

### 3.5 Numerical verification

| p | r | r' | R | -log2(prod) | log2(C(p-1,k)) | log2(E[valid]) |
|---|---|---|---|---|---|---|
| 11 | 3 | 2 | 5 | 9.56 | 7.71 | -1.85 |
| 19 | 5 | 4 | 9 | 12.63 | 15.42 | 2.79 |
| 23 | 6 | 5 | 11 | 14.22 | 19.30 | 5.08 |
| 31 | 8 | 7 | 15 | 17.46 | 27.12 | 9.66 |
| 43 | 11 | 10 | 21 | 22.43 | 38.90 | 16.47 |
| 83 | 21 | 20 | 41 | 39.64 | 78.46 | 38.82 |
| 127 | 32 | 31 | 63 | 59.21 | 122.16 | 62.95 |
| 499 | 125 | 124 | 249 | 232.52 | 493.19 | 260.67 |

E[valid D12] > 1 for all p >= 19. For p = 11, the Gaussian approximation is too crude (only R=5 variables), but p=11 is verified by exhaustive enumeration (100 valid pairs exist).

---

## 4. The Correlation Loss

### 4.1 The product-of-marginals bound

The R = (p-1)/2 representative B-values have equi-correlation rho = -2/(p-3) < 0. We need a lower bound on the joint probability:

  Pr[all B(d_i) <= T(d_i)] >= c_0 * prod_{i=1}^R Pr[B(d_i) <= T(d_i)]

for some c_0 > 0 that doesn't vanish too fast.

### 4.2 Slepian inequality (wrong direction for us)

**Correction (revision 6).** The Gaussian comparison inequality (Slepian-type) states: if Sigma_1 >= Sigma_2 componentwise off-diagonal (same diagonals), then Pr[all X_i <= t_i | Sigma_1] >= Pr[all X_i <= t_i | Sigma_2]. That is, MORE positive correlation gives LARGER joint CDF. This follows from the Plackett identity: d/d(rho_{ij}) Pr[all <= t] = phi_2(t_i, t_j; rho) >= 0.

For our case, rho = -2/(p-3) < 0 (less than independent rho = 0), so **the Gaussian joint is LESS than the product of Gaussian marginals**. Monte Carlo on the Parseval hyperplane confirms this: the Gaussian c_0 is exponentially small (~2^{-8} at p=11, ~2^{-11} at p=19). The Gaussian approximation fails badly for this problem.

### 4.3 Why the Gaussian fails but the discrete distribution works

The EXACT (discrete, k-subset) probability is much larger than the Gaussian:

| p | Exact Pr | Gaussian prod Phi | Ratio (exact/Gaussian) |
|---|---|---|---|
| 11 | 2^{-3.4} | 2^{-9.6} | 72x |
| 19 | 2^{-11.2} | 2^{-12.6} | 2.6x |
| 23 | 2^{-11.7} | 2^{-14.2} | 5.8x |

The discrete distribution has fundamentally different behavior from the Gaussian. The positive association observed in the exact joint probability (relative to Gaussian marginals) arises from the combinatorial structure of the k-subset measure: the hard sum constraint (Parseval identity) combined with asymmetric thresholds (D11 tight, D22 loose) creates an "energy redistribution" effect that favors joint feasibility. This non-Gaussian effect is invisible to any argument that goes through pairwise covariances.

### 4.4 Current status of the correlation bound

**What is proven:** For each prime p = 11, 19, 23: exhaustive enumeration verifies that an A-flat D11 has multiple valid D12 (E[valid] >> 1). For p = 31, 43, 47, 59: simulated annealing finds valid constructions.

**What is needed for all p:** A lower bound c_0 >= p^{-C} (polynomial loss) on the ratio Pr[all ok] / prod Pr[ok_i] for the discrete k-subset distribution. The (p-1)/2 bits of exponential headroom absorbs any polynomial loss. Two approaches are viable:

**(A) Second moment method (Paley-Zygmund).** For fixed A-flat D11, let N = #{valid D12}. If E[N^2]/E[N]^2 = O(poly(p)), then Pr[N > 0] >= E[N]^2/E[N^2] = 1/poly(p) > 0. This avoids the joint probability issue entirely, requiring only a bound on the variance of N. The key estimate is controlling E[N^2] = sum_{D12,D12'} Pr[both valid] by grouping pairs (D12, D12') by overlap |D12 cap D12'|.

**(B) Direct first moment with exact marginals.** Each exact marginal Pr[B(d_i) <= T_i] is within O(1/sqrt(p)) of 1/2 (Berry-Esseen). The product of R such terms is ~2^{-(p-1)/2}. If the correlation loss c_0 for the discrete distribution is at least p^{-C}, the first moment gives E[valid] >= 2^{(p-1)/2 - O(log p)} -> infinity.

---

## 5. The Main Bound

### 5.1 Unconditional first moment (what is established)

**Theorem (conditional on correlation bound).** Let p >= 19 be prime with p = 3 (mod 4). Let D11 be A-flat with parameter C = O(sqrt(p log p)). Then:

  E[#{valid D12} | D11] = C(p-1, (p-3)/2) * Pr[D12 valid | D11]

where:
- The budget: log2(C(p-1,(p-3)/2)) = p - 1 - (1/2)log2(p) + O(1)
- Each marginal: Pr[B(d_i) <= T_i] = 1/2 + O(1/sqrt(p)) (Berry-Esseen)
- Product of marginals: prod Pr[ok_i] ~ 2^{-(p-1)/2 + O(1)}

If the correlation loss satisfies c_0 >= p^{-C'} for some constant C', then:

  E[valid] >= p^{-C'} * 2^{(p-1)/2 - O(log p)} -> infinity.

### 5.2 Proof of the marginal bound

**Step 1 (A-flatness ensures bounded z-scores).** For an A-flat D11:

  T_red(d) = (p-3)/2 - A(d) >= (p-7)/4 - C    (for d in D11)
  T_blue(d) = (p+3)/2 - A(d) >= (p+5)/4 - C    (for d in D22)

Both T_red - E[B] and T_blue - E[B] are O(1) (specifically, -1 + O(C) and +2 + O(C)).

**Step 2 (Marginal probabilities).** Each B(d_i) has mean E[B] = (p-3)/4 and variance sigma^2 = (p-3)(p+1)/(16(p-2)) ~ p/16. By Berry-Esseen for degree-2 functions of hypergeometric indicators:

  Pr[B(d_i) <= T_i] = Phi(z_i) + O(1/sqrt(p))

where z_i = (T_i - E[B])/sigma = O(1/sqrt(p)).

**Step 3 (Product of marginals).** With z_i = O(1/sqrt(p)):

  log2(Pr[ok_i]) = log2(1/2 + O(1/sqrt(p))) = -1 + O(1/p)

Summing over R = (p-1)/2 representatives:

  -log2(prod Pr[ok_i]) = R + O(R/p) = (p-1)/2 + O(1)

**Step 4 (First moment).** With the correlation bound c_0 >= p^{-C'}:

  E[valid] >= C(p-1,k) * c_0 * prod Pr[ok_i]
           >= 2^{p-1-O(log p)} * p^{-C'} * 2^{-(p-1)/2 + O(1)}
            = 2^{(p-1)/2 - O(log p)} -> infinity.

### 5.3 Summary of costs

| Component | Cost (bits) |
|---|---|
| Product of R marginals | (p-1)/2 + O(1) |
| Correlation loss c_0 | **OPEN**: need c_0 >= p^{-C} |
| Berry-Esseen error | O(1) (absorbed into marginals) |
| **Total cost** | **(p-1)/2 + O(log p)** (conditional on c_0 bound) |
| **Budget** (log2 C(p-1,(p-3)/2)) | **p - 1 - (1/2)log2(p) + O(1)** |
| **Headroom** | **(p-1)/2 - O(log p)** (conditional on c_0 bound) |

### 5.4 Verified primes

| p | E[valid] (exact) | Method |
|---|---|---|
| 7 | known | Prior construction (Lidicky et al.) |
| 11 | 20 | Exhaustive enumeration |
| 19 | 18 | Exhaustive enumeration |
| 23 | 198 | Exhaustive enumeration |
| 31 | >= 1 | Simulated annealing |
| 43 | >= 1 | Simulated annealing |
| 47 | >= 1 | Simulated annealing |
| 59 | >= 1 | Simulated annealing |

---

## 6. Existence of A-flat D11

### 6.1 Recap

D11 is A-flat with parameter C if max_{d: 1<=d<=p-1} A(d) <= (p+1)/4 + C. This is a GLOBAL condition on A(d) at all positions (both D11 and D22).

### 6.2 Existence by probabilistic argument

For a random symmetric D11, each A(d) has E[A(d)] = (p+1)/4 and Var[A(d)] = O(p). Using sub-Gaussian tails for the hypergeometric distribution:

  Pr[max_d |A(d) - E[A]| > t] <= 2(p-1) * exp(-ct^2/p)

Setting t = sqrt(2p log(p)/c) gives probability <= 2(p-1)/p^2 < 1 for p >= 7. So an A-flat D11 with C = O(sqrt(p log p)) exists, and random D11 is A-flat with probability 1 - O(1/p).

### 6.3 QR balance

Enumeration data shows working D11 have balanced QR/QNR split. This is related to spectral flatness via character sums. For the asymptotic proof, the probabilistic argument in 6.2 suffices.

---

## 7. Technical Details

### 7.1 Gaussianity on the Parseval hyperplane

The B-values are not exactly Gaussian -- they are degree-2 functions of hypergeometric indicators. However, by the multivariate CLT (Berry-Esseen for functions of exchangeable random variables), the conditional distribution on the Parseval hyperplane converges to a multivariate Gaussian at rate O(1/sqrt(p)).

For the product-of-marginals bound, it suffices that each conditional marginal Pr[B(d_i) <= T(d_i) | hyperplane] = Phi(z_i) + O(1/sqrt(p)). This error is absorbed into the O(1) term in the exponent.

### 7.2 Conditional independence structure

On the Parseval hyperplane, the R = (p-1)/2 representative B-values have equi-correlation -1/(R-1). This is a single correlation parameter, and it tends to 0 as p -> infinity. The variables become asymptotically independent on the hyperplane, which is why the product-of-marginals bound is tight to leading order.

### 7.3 Why "for any D11" is FALSE

Exhaustive enumeration at p=11,19,23 shows:
- p=11: 5/10 symmetric D11 have valid D12 (50%)
- p=19: 9/126 (7.1%)
- p=23: 55/462 (11.9%)

The failing D11 have large max A(d), causing T_red(d) to be too small. For example, at p=11, a non-working D11 has max A(d) = 4, giving T_red(d_min) = 1 and z = (1 - 2)/0.82 = -1.22, so Phi(z) = 0.11. Three such D11 constraints give a product of ~0.001, which is not enough.

The proof requires A-flat D11 to ensure all z-scores are O(1/sqrt(p)) -> 0.

### 7.4 Why p=11 requires separate treatment

The product-of-marginals bound gives E[valid D12] < 1 at p=11 (see table in Section 3.5). This is because:
- Only R=5 variables: Gaussian approximation is crude
- D11 z-score z_red = -1.22 is not close to 0
- The bound loses too much from ignoring correlations

However, p=11 is verified by exhaustive enumeration: 5 working D11, each with 20 valid D12. The asymptotic proof covers p >= 19.

### 7.5 No CLT cost

A key advantage of the Parseval approach over the S1-conditioning approach: **there is no local CLT cost.** The Parseval identity holds exactly for every D12, so we don't need to condition on an intermediate sum S1 and pay for the probability of hitting that sum. This saves O(log p) bits and simplifies the proof.

---

## 8. Finite Verification

**p = 7 (n = 4):** Prior constructions (Lidicky et al. 2024).

**p = 11 (n = 6):** 5 of 10 symmetric D11 work, each with 20 valid D12. Working D11 all have max A(d) = 3 = E[A]. Verified by exhaustive enumeration (asymptotic bound is not tight here).

**p = 19 (n = 10):** 9 of 126 work, each with 18 valid D12. All working D11 have max A(d) = 5 = E[A]. Asymptotic bound gives E ~ 6.9.

**p = 23 (n = 12):** 55 of 462 work, with 22-198 valid D12 each. Working D11 have max A(d) in {6,7} (E[A] = 6). Asymptotic bound gives E ~ 33.9.

---

## 9. Summary

### 9.1 What is proven

**Theorem (L6, conditional version).** Let p >= 19 be prime with p = 3 (mod 4). Let D11 be an A-flat symmetric subset of size n = (p+1)/2. If the correlation loss for the k-subset measure satisfies c_0 >= p^{-C} for some constant C, then E[#{valid D12} | D11] -> infinity.

**Proven unconditionally:**
1. **A-flat D11 exists** (Lemma L6a, sub-Gaussian concentration + union bound).
2. **Parseval identity**: sum B(d_i) = (p-1)(p-3)/8 exactly (affine hyperplane).
3. **Marginal bounds**: each Pr[B(d_i) <= T_i] = 1/2 + O(1/sqrt(p)) (Berry-Esseen).
4. **Product of marginals**: prod Pr[ok_i] ~ 2^{-(p-1)/2}, giving (p-1)/2 bits of headroom.
5. **Exhaustive verification**: E[valid] >> 1 for p = 11, 19, 23 (and SA-verified for p <= 59).

**Open:** Step 5 in the original proof claimed c_0 >= 1 via Slepian's inequality. This is **incorrect** — the Gaussian comparison goes in the wrong direction (more positive correlation gives larger joint CDF). The Gaussian c_0 is exponentially small. However, the EXACT discrete c_0 is > 1 at all tested primes, due to non-Gaussian combinatorial structure. Proving c_0 >= p^{-C} for the k-subset measure is the remaining gap.

### 9.2 Proof status for R(B_{n-1}, B_n) = 4n-1

| Case | Status |
|---|---|
| 2n-1 = q ≡ 1 (mod 4) prime power | **PROVEN** (Paley construction) |
| 2n-1 = p ≡ 3 (mod 4) prime, p <= 59 | **VERIFIED** (enumeration/SA) |
| 2n-1 = p ≡ 3 (mod 4) prime, p > 59 | **CONDITIONAL** on c_0 >= p^{-C} |
| 2n-1 composite | **OPEN** |

### 9.3 Most promising path to close the gap

**Second moment method (Paley-Zygmund):** See `proof_L6_second_moment.md` for the full development. Key results:

1. **Constant z-sum (PROVEN):** sum z_i is exactly the same for ALL symmetric D11, due to Parseval on A-values. Leading-order variation in N(D11) vanishes.

2. **Second moment ratio (VERIFIED):** E[N^2]/E[N]^2 = 2.0, 14.0, 14.5 for p = 11, 19, 23 (all D11). Among A-flat D11 only: 1.0, 3.0, 2.95.

3. **Conditional independence:** For two D12 sets with shared part I, validity events are independent conditioned on (D11, I).

4. **Gap:** Proving E[N^2]/E[N]^2 = O(poly(p)) for all large p requires a non-Gaussian argument about N(D11) concentration.

---

## References

1. **Goldstein, L., Rinott, Y.** (1996). "Multivariate normal approximations by Stein's method and size bias couplings." J. Appl. Prob.
2. **Borcea, J., Branden, P., Liggett, T.M.** (2009). "Negative dependence and the geometry of polynomials." JAMS.
3. **Rousseau, C.C., Sheehan, J.** (1978). "On Ramsey numbers for books." J. Graph Theory.
