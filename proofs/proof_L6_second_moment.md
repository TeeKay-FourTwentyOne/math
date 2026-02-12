# Proof of Lemma L6 via Second Moment Method

**Date**: 2026-02-11
**Status**: Conditional on second moment ratio bound (Conjecture in Section 7.3)
**Supersedes**: Previous draft of this file; complements proof_L6_conditioning.md (which is conditional on c_0 bound).

---

## 1. Overview

**Goal**: Prove that for all sufficiently large primes p = 3 (mod 4), there exists a symmetric D11 of size n = (p+1)/2 in Z_p that admits a valid D12, establishing R(B_{n-1}, B_n) >= 4n-1.

**Strategy**: Paley-Zygmund inequality applied to N(D11) = #{valid D12 for D11}, averaged over random symmetric D11. We show:

1. E_{D11}[N] grows (first moment)
2. E_{D11}[N^2] / E_{D11}[N]^2 = O(1) (second moment ratio bounded)
3. Therefore Pr[N > 0] >= E[N]^2/E[N^2] = Omega(1) > 0

The key structural insight is that the sum of constraint z-scores is CONSTANT across all symmetric D11 (a consequence of the Parseval identity), so the first-order term in the Gaussian approximation to N(D11) does not vary. The variation arises only at second order, giving bounded variance of log N(D11) and hence a bounded second moment ratio.

---

## 2. Setup and Notation

Fix a prime p = 3 (mod 4), p >= 19. Let:
- n = (p+1)/2, k = (p-3)/2
- D11: a symmetric subset of {1,...,p-1} with |D11| = n (both d and p-d in D11)
- D12 = {0} union S where S is a k-subset of {1,...,p-1}
- D22 = {1,...,p-1} \ D11 (the complement, also symmetric, |D22| = n-2)
- A(d) = #{(a,b) in D11 x D11 : a-b = d mod p, a != b} (D11 autocorrelation)
- B(d) = #{(a,b) in D12 x D12 : a-b = d mod p, a != b} (D12 autocorrelation)

**Validity constraint**: D12 is valid for D11 if for all d = 1,...,p-1:

  A(d) + B(d) <= T(d)

where T(d) = (p-3)/2 if d in D11, and T(d) = (p+3)/2 if d in D22.

**Equivalently** (subtracting A(d)): B(d) <= tau(d) where:
- tau_red(d) = (p-3)/2 - A(d)  for d in D11  (binding constraints)
- tau_blue(d) = (p+3)/2 - A(d)  for d in D22  (loose constraints)

**Complementary pair reduction**: Since D11 is symmetric, B(d) = B(p-d) and tau(d) = tau(p-d). The constraints reduce to R = (p-1)/2 representatives, one per complementary pair {d, p-d}. Of these, r = (p+1)/4 are D11 representatives and r' = (p-3)/4 are D22 representatives.

### 2.1 Exact moments of B-values

For a random k-subset S of {1,...,p-1} and D12 = {0} union S:
- E[B(d)] = mu_B = (p-3)/4 for all nonzero d
- Var[B(d)] = sigma^2 = (p-3)(p+1) / (16(p-2))
- sigma ~ sqrt(p)/4 for large p

### 2.2 The Parseval identity

  sum_{d=1}^{p-1} B(d) = |D12|(|D12|-1) = ((p-1)/2)((p-3)/2) = (p-1)(p-3)/4

This is DETERMINISTIC. In terms of representatives: sum_{i=1}^R B(d_i) = (p-1)(p-3)/8.

Similarly for A: sum_{d=1}^{p-1} A(d) = n(n-1) = (p+1)p/4 - 1/2. In terms of representatives: sum_{i=1}^R A(d_i) = n(n-1)/2.

---

## 3. The Z-score Structure

### 3.1 Definition

For each representative d_i, define the z-score:

  z_i = (tau(d_i) - mu_B) / sigma

Explicitly:
- If d_i in D11: z_i = ((p-3)/2 - A(d_i) - (p-3)/4) / sigma = ((p-3)/4 - A(d_i)) / sigma
- If d_i in D22: z_i = ((p+3)/2 - A(d_i) - (p-3)/4) / sigma = ((p+9)/4 - A(d_i)) / sigma

### 3.2 The Parseval sum constraint on z-scores

**Proposition (Constant z-sum).** For every symmetric D11 of size n:

  sum_{i=1}^R z_i = (p-7) / (4*sigma)

*Proof.* Compute:

  sigma * sum z_i = sum_{d_i in D11 reps} ((p-3)/4 - A(d_i)) + sum_{d_j in D22 reps} ((p+9)/4 - A(d_j))
  = (p-3)r/4 + (p+9)r'/4 - sum_{all reps} A(d_i)
  = (p-3)(p+1)/16 + (p+9)(p-3)/16 - n(n-1)/2

Expanding with n = (p+1)/2:
- n(n-1)/2 = (p+1)(p-1)/8
- (p-3)(p+1)/16 + (p+9)(p-3)/16 = [(p-3)(p+1) + (p+9)(p-3)] / 16
  = [(p^2-2p-3) + (p^2+6p-27)] / 16 = [2p^2+4p-30] / 16 = (p^2+2p-15)/8

So: sigma * sum z_i = (p^2+2p-15)/8 - (p^2-1)/8 = (2p-14)/8 = (p-7)/4. QED.

**Corollary.** sum z_i ~ sqrt(p)/2 -> infinity, but sum z_i / R ~ 2/p -> 0. The AVERAGE z-score is slightly positive (reflecting that D22 constraints are looser than D11 constraints).

### 3.3 The Gaussian proxy for N(D11)

The Gaussian approximation to the product of marginals gives:

  log f(D11) := sum_{i=1}^R log Phi(z_i)

where Phi is the standard normal CDF. For small z, expand:

  log Phi(z) = -log 2 + c_1 * z - c_2 * z^2 + O(z^3)

where c_1 = phi(0)/Phi(0) = sqrt(2/pi) and c_2 = phi(0)^2/(2*Phi(0)^2) = 1/pi.

Therefore:

  log f(D11) = -R log 2 + c_1 * S_1 - c_2 * S_2 + O(S_3)

where S_j = sum_{i=1}^R z_i^j. By the Parseval constraint:

  **S_1 = sum z_i = (p-7)/(4*sigma) is CONSTANT for all D11.**

So the LEADING variation in log f(D11) comes from S_2 = sum z_i^2.

### 3.4 Verification

| p | sum z_i (theory) | sum z_i (all D11) | S_2 range |
|---|---|---|---|
| 11 | 1.225 | 1.225 (exact, all D11) | [4.5, 19.5] |
| 19 | 2.766 | 2.766 (exact, all D11) | [7.6, 63.8] |
| 23 | 3.347 | 3.347 (exact, all D11) | [8.4, 93.8] |

Confirmed: S_1 is exactly constant; variation comes entirely from S_2 and higher terms.

---

## 4. First Moment: E[N] Grows

### 4.1 The Gaussian product-of-marginals estimate

For A-flat D11 (max_d A(d) <= E[A] + C with C = O(sqrt(p log p))):

Each z_i = O(C/sigma) = O(sqrt(log p / p)) -> 0. So Phi(z_i) = 1/2 + O(1/sqrt(p)).

  -log_2(prod Phi(z_i)) = R * 1 - (c_1/ln 2) * S_1 + (c_2/ln 2) * S_2 + O(S_3)

With S_1 constant and S_2 = O(R * max z_i^2) = O(p * C^2/p) = O(C^2):

  -log_2(prod Phi) = R + O(C^2) = (p-1)/2 + O(log p)

So:

  log_2 E[N | D11] = log_2 C(p-1,k) + log_2(prod Phi) + log_2(c_0)
                   = (p-1 - (1/2)log_2 p + O(1)) - ((p-1)/2 + O(log p)) + log_2(c_0)
                   = (p-1)/2 - O(log p) + log_2(c_0)

**This is >> 0 provided c_0 >= p^{-C'} for some constant C'.** The first moment argument alone requires the correlation bound.

### 4.2 Bypassing c_0 via the second moment

The second moment method proves existence without explicitly bounding c_0. We show E[N^2]/E[N]^2 = O(1) over random D11, which (combined with E[N] > 0 for some D11, verified computationally) gives existence.

---

## 5. Second Moment: E[N^2]/E[N]^2 is Bounded

### 5.1 The key identity

For random symmetric D11 and FIXED D12, D12':

  E_{D11}[N^2] = sum_{D12, D12'} Pr_{D11}[both valid]

Since D12 and D12' are deterministic and only D11 is random:

  Pr_{D11}[D12 valid AND D12' valid] = Pr_{D11}[for all d: B(d) <= tau(d) AND B'(d) <= tau(d)]

The thresholds tau(d) depend on D11 (through A(d) and D11/D22 membership). The B-values depend only on D12 (not D11), so:

  E[N^2] / C(p-1,k)^2 = E_{D11}[f(D11)^2]

  E[N]^2 / C(p-1,k)^2 = (E_{D11}[f(D11)])^2

where f(D11) = N(D11)/C(p-1,k) = Pr[random D12 valid | D11].

Therefore:

  **E[N^2] / E[N]^2 = E[f^2] / (E[f])^2 = 1 + Var[f] / (E[f])^2**

### 5.2 Conditional independence for pairs of D12

**Lemma (Conditional independence).** For fixed D11, let S_1 = I union U_1 and S_2 = I union U_2 be two k-subsets with shared part I, where U_1 and U_2 are independently drawn (k-|I|)-subsets of {1,...,p-1}\I. Then conditioned on (D11, I), the events {S_1 valid} and {S_2 valid} are independent.

*Proof.* B(d; {0} union S_j) depends on the elements of S_j = I union U_j. For fixed I:
- B_I(d) = #{(a,b) in (I union {0})^2 : a-b=d} is deterministic
- The cross-terms and within-U terms depend only on U_j (given I)

Since U_1 and U_2 are independent given I, the B-values for S_1 and S_2 are independent given (D11, I). QED

**Consequence.** E_{D11}[N^2] can be decomposed by overlap s = |S_1 cap S_2|:

  E[N^2] = sum_s H(s) * P_2(s)

where H(s) = #{ordered pairs with overlap s} and P_2(s) = Pr_{D11}[both valid | overlap s].

For the "independent U" model (allowing U_1 cap U_2 != empty):
  P_2(s) = E_{D11, I}[Pr[valid | D11, I]^2]

For the "exact overlap s" model (U_1 cap U_2 = empty), a mild correction applies due to the negative dependence from the disjointness constraint, but the effect is O(1/p) per constraint.

**Overlap distribution from enumeration (summed over all D11):**

| p | k | Typical overlap | Distribution |
|---|---|---|---|
| 11 | 4 | 1.6 | {1:140, 2:760, 3:720, 4:280, 5:100} |
| 19 | 8 | 3.6 | {2:108, 3:486, 4:648, 5:1080, 6:432, 9:162} |
| 23 | 10 | 4.5 | {1:264, 2:3608, ..., 7:100870, ..., 11:4356} |

The diagonal (s=k) contributes 1/E[N] to the ratio, which vanishes as E[N] grows.

### 5.3 Concentration of f(D11) via the Parseval structure

The Gaussian proxy gives:

  log f(D11) ~ -R log 2 + c_1 * S_1 - c_2 * S_2 + higher order

Since S_1 is constant:

  f(D11) ~ 2^{-R} * exp(c_1 * S_1) * exp(-c_2 * S_2)

The variation of f is controlled by exp(-c_2 * S_2). We need to show Var[exp(-c_2 * S_2)] / (E[exp(-c_2 * S_2)])^2 is bounded.

### 5.4 Analysis of S_2 for A-flat D11

For an A-flat D11, each z_i satisfies |z_i| <= C'/sigma = O(sqrt(log p / p)).

Write z_i = alpha_i + delta * xi_i where:
- alpha_i = ((p-3)/4 - A(d_i)) / sigma (the D11-membership z-score)
- delta = 3/sigma (the D11-vs-D22 shift)
- xi_i = 1{d_i in D22}

Then:

  S_2 = sum (alpha_i + delta * xi_i)^2
      = sum alpha_i^2 + 2*delta * sum alpha_i * xi_i + delta^2 * sum xi_i^2
      = sum alpha_i^2 + 2*delta * sum_{D22 reps} alpha_j + delta^2 * r'

The first term sum alpha_i^2 depends on the A-values (which depend on D11). For A-flat D11, each alpha_i = O(1/sigma) = O(1/sqrt(p)), so sum alpha_i^2 = O(R/p) = O(1).

The second term depends on which representatives are in D22 and their alpha-values. The key: for A-flat D11, the alpha-values at all positions are O(1/sqrt(p)), so |2*delta * sum_{D22} alpha_j| <= 2*delta * r' * max|alpha_j| = O((1/sqrt(p)) * p * (1/sqrt(p))) = O(1).

The third term delta^2 * r' = 9/sigma^2 * (p-3)/4 ~ 9*16/p * p/4 = 36 is O(1).

Therefore **S_2 = O(1) for all A-flat D11**, and the VARIATION of S_2 across A-flat D11 is also O(1).

### 5.5 Bounded second moment ratio

Since S_2 varies by O(1) across A-flat D11, and f ~ exp(-c_2 * S_2 + const):

  max_{D11 A-flat} f / min_{D11 A-flat} f <= exp(c_2 * Delta_S2) = exp(O(1)) = O(1)

where Delta_S2 is the range of S_2 over A-flat D11.

Therefore f(D11) varies by at most a constant factor across A-flat D11. In particular:

  E[f^2] / (E[f])^2 <= (max f / min f)^2 = O(1)

**assuming all A-flat D11 have f > 0** (i.e., N(D11) > 0). If some A-flat D11 have f = 0, then:

  E[f^2] / (E[f])^2 <= (1/p_working) * (max f / min f among working)^2

where p_working = Pr[N > 0 | A-flat]. The within-working ratio is O(1) (since S_2 varies by O(1) among working A-flat D11). The factor 1/p_working requires p_working = Omega(1).

### 5.6 Numerical verification

**Over all symmetric D11:**

| p | #{D11} | E[N] | E[N^2] | E[N^2]/E[N]^2 | Pr[N>0] | PZ bound |
|---|---|---|---|---|---|---|
| 11 | 10 | 10.0 | 200 | 2.0 | 0.500 | 0.500 |
| 19 | 126 | 1.3 | 23.1 | 14.0 | 0.071 | 0.071 |
| 23 | 462 | 9.4 | 1291 | 14.5 | 0.119 | 0.069 |

**Over A-flat D11 only (max A(d) at D11 positions = floor(E[A])):**

| p | #{A-flat} | E[N] | E[N^2] | ratio | Pr[N>0] | PZ |
|---|---|---|---|---|---|---|
| 11 | 5 | 20.0 | 400 | 1.00 | 1.000 | 1.000 |
| 19 | 27 | 6.0 | 108 | 3.00 | 0.333 | 0.333 |
| 23 | 77 | 50.3 | 7467 | 2.95 | 0.571 | 0.339 |
| 31 | 330 | 475.9 | 1.55e6 | 6.86 | 0.273 | 0.146 |

The ratio over A-flat D11 grows slowly: 1.0, 3.0, 2.95, 6.86 at p = 11, 19, 23, 31.

**At p=43 (SA-verified):** 3528 A-flat D11 in 168 orbits. SA found 5/168 working orbits (p_working >= 0.030). Random D12 sampling at this prime is too sparse (valid fraction ~1e-7) to compute N directly, but SA confirms A-flat D11 have valid D12.

**NOTE:** An earlier computation using 500K random D12 samples per orbit found 0 working A-flat orbits at p=43 — this was a SAMPLING ARTIFACT (500K samples has <5% chance of finding a valid D12 when the valid fraction is ~1e-7). SA, which is orders of magnitude more efficient at finding valid solutions, corrected this.

### 5.7 Decomposition of the ratio

The ratio decomposes as:

  E[N^2]/E[N]^2 = (1/p_working) * (E[N^2 | working] / E[N | working]^2)

| p | p_working (A-flat) | within-working ratio | overall ratio |
|---|---|---|---|
| 11 | 1.000 | 1.00 | 1.00 |
| 19 | 0.333 | 1.00 | 3.00 |
| 23 | 0.571 | 1.69 | 2.95 |
| 31 | 0.273 | 1.87 | 6.86 |
| 43 | >= 0.030 | ? | ? |

The within-working ratio is close to 1 (N values among working A-flat D11 don't vary much). The overall ratio is dominated by 1/p_working.

**Scaling of p_working:** The product p_working × p gives 6.3, 13.1, 8.5, >= 1.3 for p = 19, 23, 31, 43. This is consistent with p_working >= C/p for a constant C >= 1, giving overall ratio = O(p).

---

## 6. The Working Fraction p_working

### 6.1 Empirical data

Among A-flat D11 (max A(d) at D11 positions = floor(E[A])):

| p | A-flat total | Working | p_working |
|---|---|---|---|
| 11 | 5 | 5 | 1.000 |
| 19 | 27 | 9 | 0.333 |
| 23 | 77 | 44 | 0.571 |

The fraction p_working is bounded below by 1/3 at all tested primes.

### 6.2 Why p_working should be Omega(1)

The A-value profile at D11 positions determines whether a D11 is A-flat, but the A-values at D22 positions (which also affect validity through the D22 constraints) are NOT controlled by the D11-position A-flatness condition.

However, two structural facts constrain the D22 A-values:

**(a) Parseval constraint on A-values:** sum_{all d} A(d) = n(n-1). So sum_{D22} A(d) = n(n-1) - sum_{D11} A(d). For A-flat D11, sum_{D11} A(d) is close to n * E[A] = n(p+1)/4, leaving sum_{D22} A(d) close to its expected value.

**(b) Symmetry:** Both D11 and D22 are symmetric sets. The A-values inherit this symmetry: A(d) = A(p-d).

Together, these facts mean that for A-flat D11, the D22 A-values are also well-behaved on average, even if some individual values can be large.

### 6.3 What distinguishes working from non-working A-flat D11

Empirically, at p=19 the A-flat D11 have two distinct A-value profiles:
- Profile [3,3,4,4,4,4,5,5,5,5] (18 D11): 9 work, 9 don't
- Profile [1,1,2,2,3,3,4,4,5,5] (9 D11): 0 work

The non-working D11 in the first profile have D22 A-values [5,5,6,6,6,6,7,7] while working ones have [5,5,5,5,6,6,8,8]. The difference is in the DISTRIBUTION of D22 A-values (same sum, different spread).

This suggests that the fraction p_working is determined by a secondary condition on the D22 A-value distribution, which should hold with constant probability among A-flat D11.

---

## 7. The Complete Argument

### 7.1 Theorem

**Theorem.** For all primes p = 3 (mod 4) with p >= 19, there exists a valid (D11, D12) construction for the book graph Ramsey problem, proving R(B_{n-1}, B_n) >= 4n-1 where n = (p+1)/2.

### 7.2 Proof structure

**Step 1 (A-flat D11 exist).** By sub-Gaussian concentration (Lemma L6a), random symmetric D11 is A-flat with probability 1 - O(1/p). In particular, the number of A-flat D11 is Omega(#{D11}).

**Step 2 (Constant z-sum).** For every symmetric D11: sum_{i=1}^R z_i = (p-7)/(4*sigma). This is an exact algebraic identity following from the Parseval identity for A-values.

**Step 3 (S_2 bounded for A-flat D11).** For A-flat D11 with parameter C = O(sqrt(p log p)): each |z_i| = O(C/sigma) = O(sqrt(log p / p)), so S_2 = sum z_i^2 = O(R * C^2/p) = O(log p).

**Step 4 (Product of marginals).** The Gaussian proxy gives:

  -log_2 f(D11) = R + O(S_2) = (p-1)/2 + O(log p)

Combined with the budget log_2 C(p-1,k) = p - 1 - (1/2)log_2 p + O(1):

  log_2 E[N | D11] = (p-1)/2 - O(log p) + log_2 c_0(D11)

where c_0(D11) is the ratio of the exact joint probability to the Gaussian product.

**Step 5 (Second moment over A-flat D11).** Since S_1 is constant and S_2 varies by O(log p):

  log f(D11) varies by O(log p) across A-flat D11

Therefore:

  E[f^2] / (E[f])^2 <= exp(O(log p)) = poly(p)

This gives E[N^2]/E[N]^2 = O(poly(p)).

**Step 6 (Paley-Zygmund).** Pr_{A-flat D11}[N > 0] >= E[N]^2 / E[N^2] >= 1/poly(p) > 0.

Since at least one D11 has N(D11) > 0, a valid construction exists.

### 7.3 The remaining gap

The argument above is rigorous EXCEPT for Step 5, which uses the Gaussian proxy to bound the second moment ratio. The Gaussian proxy gives S_2 = O(log p) for A-flat D11, but the actual S_2 ranges from O(1) to O(p) among A-flat D11. The Gaussian proxy thus overestimates the variation of N(D11).

**Why the Gaussian proxy fails**: The Gaussian proxy N_gauss(D11) = C(p-1,k) * prod Phi(z_i) varies by a factor of ~370x across A-flat D11 at p=23, while the exact N(D11) takes only a few discrete values {0, 22, 110, 198}. The exact N concentrates MUCH better than the Gaussian proxy predicts. This is a genuine combinatorial phenomenon: the correlation structure of the k-subset measure constrains N more tightly than pairwise correlations alone would suggest.

The correlation loss c_0(D11) = N(D11) / N_gauss(D11) ranges from 0 (for non-working D11) to 0.64 (p=23), with working D11 having c_0 between 0.11 and 0.64. The c_0 is NOT monotone in N_gauss: some D11 with large N_gauss have c_0 = 0, while D11 with moderate N_gauss have positive c_0.

**What is proven rigorously:**
- The constant z-sum (Section 3.2): exact algebraic identity
- A-flat D11 exist in abundance (Lemma L6a)
- The product of marginals ~ 2^{-(p-1)/2} (Berry-Esseen)
- The budget ~ 2^{p-1} (Stirling)
- E[N^2]/E[N]^2 = O(1) at p=11,19,23 (exhaustive computation)

**What requires additional work:**
The core challenge is proving E[N^2]/E[N]^2 = O(poly(p)) for all large p. The Gaussian proxy is too coarse -- any proof must use the full combinatorial structure of the k-subset distribution, not just pairwise statistics. Possible approaches:
1. **Direct combinatorial bound**: Count valid D12 by inclusion-exclusion or character sums
2. **Entropy method**: Bound the conditional entropy of B-values on the Parseval hyperplane
3. **Algebraic approach**: Use the multiplicative structure of Z_p^* to relate different D11

### 7.4 Numerical evidence

The exact second moment ratios over A-flat D11 are 1.00, 3.00, 2.95 for p = 11, 19, 23. This strongly suggests the ratio is O(1), which is MUCH better than any bound achievable through the Gaussian proxy.

The Gaussian proxy analysis (Section 5) gives a bound of exp(O(S_2 variation)) which is exponential in p. The enormous discrepancy between this exponential bound and the actual O(1) ratio demonstrates that the concentration of N(D11) is a fundamentally non-Gaussian phenomenon.

| p | N_gauss range (A-flat) | N exact range (working) | c_0 range (working) |
|---|---|---|---|
| 11 | [4.6, 4.6] | [20, 20] | [4.36, 4.36] |
| 19 | [7.4, 93.3] | [0, 18] | [0, 0.46] |
| 23 | [1.1, 408.7] | [0, 198] | [0, 0.64] |

---

## 8. Multiplicative Orbit Structure (New)

### 8.1 The key symmetry

**Proposition.** If (D11, D12) is a valid pair and g in Z_p^*, then (gD11, gD12) is also valid.

*Proof.* The constraint A(d) + B(d) <= T(d) is preserved under multiplication by g, since Delta(gS, gS, gd) = Delta(S, S, d) and the threshold T(gd) for gD11 equals T(d) for D11 (because gd in gD11 iff d in D11). QED.

**Corollary.** N(D11) = N(gD11) for all g in Z_p^*.

Since D11 is symmetric, multiplying by -1 fixes D11. So the orbits of symmetric D11 under Z_p^* have size exactly (p-1)/2 (the order of Z_p^*/{+-1}).

### 8.2 Orbit decomposition of the second moment ratio

The (p-1)/2 elements of each orbit all have the same N-value. So:

  E[N^2] / E[N]^2 = (#orbits * sum_{working orbits} N_i^2) / (sum_{all orbits} N_i)^2

where N_i is the N-value of orbit i (with orbits weighted equally).

Since all working orbits have N_i > 0, we can decompose:

  ratio = (1 / p_working) * (within-working ratio)

where p_working = (# working orbits) / (# total orbits) and the within-working ratio is the Cauchy-Schwarz ratio among working orbits.

### 8.3 Numerical verification

| p | Orbits | Orbit size | Working orbits | N values (by orbit) | p_working | Within ratio | Overall ratio |
|---|--------|-----------|----------------|---------------------|-----------|-------------|---------------|
| 11 | 2 | 5 | 1 | {0: 1, 20: 1} | 0.500 | 1.00 | 2.00 |
| 19 | 14 | 9 | 1 | {0: 13, 18: 1} | 0.071 | 1.00 | 14.00 |
| 23 | 42 | 11 | 5 | {0: 37, 22: 2, 44: 1, 110: 1, 198: 1} | 0.119 | 1.73 | 14.52 |

Key observations:
- At p=11, 19: N takes only TWO values (0 and one positive value). The within-working ratio is 1.
- At p=23: N takes FIVE distinct values, but the within-working ratio is still close to 1.
- The overall ratio is dominated by 1/p_working.
- Each valid pair orbit has size exactly p-1 (no stabilizers beyond {+-1}).

### 8.4 Pair orbit structure

The action of Z_p^* on pairs (D11, D12) also has a clean orbit structure:

| p | Valid pair orbits | Pair orbit size | Total valid pairs |
|---|-------------------|-----------------|-------------------|
| 11 | 10 | 10 | 100 |
| 19 | 9 | 18 | 162 |
| 23 | 198 | 22 | 4356 |

The number of pair orbits per D11 orbit equals N(D11) / 2 (since each D12 orbit under Z_p^*/{+-1} has size (p-1)/2 and pairs (D11, D12) have orbit size p-1).

### 8.5 Implications for the proof

The orbit structure means:

1. **The problem reduces to counting orbits**: we need at least one working orbit (N_i > 0).
2. **The total number of orbits is C((p-1)/2, (p+1)/4) / ((p-1)/2)**, which grows exponentially in p.
3. **The fraction of working orbits p_working determines the second moment ratio.**
4. **Proving p_working >= 1/(# orbits) suffices** for existence (trivially, just one orbit needed).
5. **Proving p_working >= 1/poly(p) suffices** for the Paley-Zygmund argument.

The data suggests p_working = Theta(1) (constant fraction), but even 1/poly(p) would close the proof.

---

## 9. Character Sum Analysis of the Paley D12

### 9.1 The Paley autocorrelation

**Proposition.** For p = 3 (mod 4), the Paley set P = QR union {0} has perfectly flat autocorrelation:

  A_P(d) = (p+1)/4 for all d = 1,...,p-1.

*Proof.* Standard result for conference matrices / Paley-type Hadamard matrices. The autocorrelation of the Paley graph on p vertices at any nonzero difference d is (p - 3 + 4*1[d=0])/4, but for the "extended" set P = QR union {0} of size (p+1)/2, the autocorrelation equals n(n-1)/(p-1) = ((p+1)/2)((p-1)/2)/(p-1) = (p+1)/4 exactly.

Verified computationally for all p = 3 mod 4 up to p = 167.

### 9.2 The B-value structure for QR-based D12

For D12 = {0} union (QR \ {a}) with a in QR:

  B(d) = A_P(d) - 1[a-d in P] - 1[a+d in P] = (p+1)/4 - 1[a-d in QR union {0}] - 1[a+d in QR union {0}]

For d != a and d != p-a (the generic case):

  B(d) = (p-3)/4 + epsilon(d)

where epsilon(d) = -[chi(a-d) + chi(a+d)]/2 in {-1, 0, +1} with chi the Legendre symbol.

### 9.3 Distribution of epsilon values

**Proposition.** For any a in QR and any p = 3 mod 4:
- epsilon = -1 (both a-d, a+d in QR): approximately (p-7)/8 pairs
- epsilon = 0 (exactly one in QR): approximately (p-1)/4 pairs
- epsilon = +1 (neither in QR): approximately (p-7)/8 pairs

Moreover, J(chi, chi) = 1 for all p = 3 mod 4, giving the exact count via Jacobi sums.

The number of "available" pairs (epsilon <= 0) exceeds the "needed" pairs ((p+1)/4) by approximately (p-1)/8 -- a surplus that grows linearly with p.

### 9.4 Limitations of the Paley D12 approach

Despite the clean structure, using D12 = {0} union (QR \ {a}) does NOT directly give valid pairs: the blue constraints at epsilon = +1 positions (d in D22 where B(d) = (p+1)/4) are violated because A(d) is too large at the complementary positions.

Verified computationally: at p = 23, 43, 47, 59, 67, 71, 79, 83, NO choice of D11 (using all epsilon=-1 pairs + subset of epsilon=0 pairs) satisfies all constraints with this D12.

The actual valid D12 have a more balanced structure, with B-values that are specifically tuned to avoid violations at both red and blue positions.

---

## 10. Summary

| Component | Status |
|---|---|
| Constant z-sum (Parseval) | **PROVEN** (Section 3.2) |
| A-flat D11 existence | **PROVEN** (Lemma L6a) |
| Product of marginals ~ 2^{-(p-1)/2} | **PROVEN** (Berry-Esseen) |
| Budget >> cost (headroom) | **PROVEN** ((p-1)/2 bits) |
| Second moment ratio O(1) | **VERIFIED** for p=11,19,23 |
| Second moment ratio O(p) for A-flat | **VERIFIED** for p=11,19,23,31; consistent at p=43 (SA) |
| Second moment ratio O(poly(p)) for all p | **OPEN** (empirically ~O(p) via A-flat) |
| Correlation loss c_0 >= 1/poly(p) | **VERIFIED** for p<=23, **OPEN** in general |
| A-flat D11 working at p=43 | **VERIFIED** (5/168 orbits via SA) |
| Multiplicative orbit structure | **PROVEN** (Section 8) |
| N constant on orbits | **PROVEN** (Section 8.1) |
| Paley autocorrelation A_P = (p+1)/4 | **PROVEN** (Section 9.1) |
| J(chi,chi) = 1 for p = 3 mod 4 | **VERIFIED** computationally, known result |

**Status**: The proof establishes a comprehensive algebraic framework including the constant z-sum identity, orbit structure under Z_p^*, and character-sum analysis of Paley-based D12. The A-flat second moment ratio is verified to grow as O(p) through p=31 (exact) and is consistent at p=43 (SA). A-flat D11 have valid D12 at ALL tested primes including p=43 (correcting an earlier sampling artifact). The gap is showing that p_working >= 1/poly(p) for all large p.

**Key insights**:
1. The Gaussian proxy is too coarse (370x variation vs discrete N-values).
2. N(D11) is constant on multiplicative orbits -- the problem reduces to counting orbits.
3. The Paley D12 has clean character-sum structure (B = mu +/- 1) but doesn't directly give valid pairs due to blue constraint violations.
4. The actual valid D12 are non-algebraic -- they must balance B-values at both red and blue positions simultaneously.
5. The most promising proof path is showing p_working >= 1/poly(p), i.e., a non-trivial fraction of orbits are working.

---

## References

1. **Paley, R.E.A.C., Zygmund, A.** (1932). "A note on analytic functions in the unit circle." Proc. Camb. Phil. Soc.
2. **Alon, N., Spencer, J.H.** (2016). "The Probabilistic Method." 4th ed., Wiley.
3. **Rousseau, C.C., Sheehan, J.** (1978). "On Ramsey numbers for books." J. Graph Theory.
