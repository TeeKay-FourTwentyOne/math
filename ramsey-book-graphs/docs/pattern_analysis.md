# Deep Pattern Analysis of Known Ramsey Book Graph Constructions

**Author**: researcher agent
**Date**: 2026-02-10
**Data**: Known constructions for n=6,8,10,12,14,16,18,20,22 (all m = 3 mod 4 cases)

---

## 0. Universal Structural Properties

Every known construction satisfies:

| Property | Value | Verified |
|----------|-------|----------|
| D22 = {1,...,m-1} \ D11 | Always | All 9 cases |
| D11 symmetric (x in D11 iff -x in D11) | Always | All 9 cases |
| D22 symmetric | Always (follows from D11 symmetric + complement) | All 9 cases |
| D12 NOT symmetric | Always | All 9 cases |
| 0 in D12 | Always | All 9 cases |
| \|D12\| = n-1 | Always | All 9 cases |

### Size Patterns

For n >= 12 (the "normal" regime):

| n | m | \|D11\| | \|D12\| | \|D22\| | d1 | d2 |
|---|---|---------|---------|---------|----|----|
| 12 | 23 | 12 (=n) | 11 (=n-1) | 10 (=n-2) | 23 (=m) | 21 (=m-2) |
| 14 | 27 | 14 | 13 | 12 | 27 | 25 |
| 16 | 31 | 16 | 15 | 14 | 31 | 29 |
| 18 | 35 | 18 | 17 | 16 | 35 | 33 |
| 20 | 39 | 20 | 19 | 18 | 39 | 37 |
| 22 | 43 | 22 | 21 | 20 | 43 | 41 |

Pattern: **|D11| = n, |D12| = n-1, |D22| = n-2, d1 = m, d2 = m-2**.

For small n (6,8,10), the roles of "large block" and "small block" are flipped:
|D11| = n-2, |D12| = n-1, |D22| = n, giving d1 = m-2 and d2 = m.

The transition happens at n=12 where |D11| first equals n.

---

## 1. Theorem: V1V2 Constraints Are Automatic

**This is the most important structural discovery.**

**Theorem.** If D22 = {1,...,m-1} \ D11, D11 is symmetric, and |D12| = n-1, then for ALL d in {0,...,m-1}:

```
Sigma(D11, D12, d) + Delta(D12, D22, d) = |D12| - [d in D12]
```

Specifically:
- For d in D12 (red edge): common neighbors = n-2 (exactly at red threshold)
- For d not in D12 (blue edge): common neighbors = n-1 (exactly at blue threshold)

**Proof.**

Write f(d) = Sigma(D11, D12, d) + Delta(D12, D22, d).

Since D22 = {1,...,m-1} \ D11, we expand Delta(D12, D22, d):

```
Delta(D12, D22, d) = |{b in D12 : (b-d) mod m != 0 and (b-d) mod m not in D11}|
                   = (|D12| - [d in D12]) - Delta(D12, D11, d)
```

(The (b-d)=0 case subtracts [d in D12], and then we subtract those b where (b-d) in D11.)

Note that 0 is not in D11, so Delta(D12, D11, d) = |{b in D12 : (b-d) in D11}|, which includes b=d only if 0 in D11 (it's not), so no special handling needed.

Therefore:
```
f(d) = Sigma(D11, D12, d) - Delta(D12, D11, d) + |D12| - [d in D12]
```

Now the key:
- **Delta(D12, D11, d) = Sigma(D12, D11, d)** because D11 is symmetric: 1_{D11}(b-d) = 1_{D11}(d-b).
- **Sigma(D11, D12, d) = Sigma(D12, D11, d)** by commutativity of convolution.

Therefore Sigma(D11, D12, d) - Delta(D12, D11, d) = 0 for all d, and:

```
f(d) = |D12| - [d in D12]
```

Setting |D12| = n-1 gives f(d) = n-2 for d in D12 and n-1 for d not in D12. QED.

**Consequence**: The V1V2 block is completely free -- any D11 (symmetric), D12 with |D12|=n-1 and 0 in D12, and D22 = complement(D11) will automatically satisfy all V1V2 constraints. The problem reduces to satisfying V1V1 and V2V2 only.

---

## 2. Theorem: V2V2 is Algebraically Determined by V1V1

**Theorem.** With D22 = complement(D11), define:
- A(d) = Delta(D11, D11, d, m)
- B(d) = Delta(D12, D12, d, m)
- V1V1(d) = A(d) + B(d)
- V2V2(d) = Delta(D22, D22, d, m) + Delta(D12^T, D12^T, d, m)

Then:
```
V2V2(d) = A(d) + B(m-d) + (m - 2 - 2|D11|) + 2*[d in D11]
```

**Proof.** Two identities:

1. Delta(D22, D22, d) = A(d) + (m-2-2|D11|) + 2[d in D11], derived from 1_{D22} = 1_{S} - 1_{D11} where S = {1,...,m-1}, expanding the product.

2. Delta(D12^T, D12^T, d) = Delta(D12, D12, m-d), because substituting a = -a' transforms the counting.

Both identities verified computationally for all 9 known constructions.

**Key structural consequences:**
- d in D11 (red in V1V1) iff d not in D22 (blue in V2V2) -- **colorings are flipped**
- The V2V2 red constraint on d (where d in D22, i.e., d not in D11) involves A(d) which is already constrained by V1V1 blue
- When D12 is "nearly symmetric" (B(d) ~ B(m-d)), V2V2 constraints are nearly redundant given V1V1

---

## 3. Spectral / Fourier Analysis of D11

For each known D11, I computed the character sums chi_k(D11) = sum_{d in D11} omega^{kd}.

### Key findings:

| n | m | \|D11\| | avg \|chi_k\|^2 | max \|chi_k\|^2 | min \|chi_k\|^2 | flatness ratio |
|---|---|---------|------------------|------------------|------------------|----------------|
| 6 | 11 | 4 | 2.80 | 6.32 | 0.06 | 113.0 |
| 8 | 15 | 6 | 3.86 | 9.00 | 0.11 | 78.7 |
| 10 | 19 | 8 | 4.89 | 11.84 | 0.12 | 102.9 |
| 12 | 23 | 12 | 6.00 | 25.86 | 0.11 | 228.4 |
| 14 | 27 | 14 | 7.00 | 35.14 | 0.003 | 13907.2 |
| 16 | 31 | 16 | 8.00 | 43.88 | 0.37 | 120.2 |
| 18 | 35 | 18 | 9.00 | 47.59 | 0.13 | 364.9 |
| 20 | 39 | 20 | 10.00 | 49.36 | 0.05 | 1066.9 |
| 22 | 43 | 22 | 11.00 | 41.53 | 0.00 | inf |

**Observations:**
- D11 is **far from a difference set** (flatness ratio >> 1). A perfect difference set would have |chi_k|^2 constant.
- The average |chi_k|^2 = |D11|(|D11|-1)/(m-1) + |D11|/m ... actually avg = (m*|D11| - |D11|^2)/(m-1) = |D11|(m-|D11|)/(m-1) when |D11| divides appropriately. For n >= 12, avg = |D11|. Intriguing.
- The spectrum is highly non-flat: a few large peaks and many near-zero values.
- For n=22 (m=43): two frequencies have |chi_k|^2 = 0 exactly (k=20 and k=23). This means D11 is an "annihilator" for those characters.

**Interpretation**: D11 is NOT chosen to minimize spectral fluctuation (as a difference set would). Instead, the large spectral peaks correspond to directions where Delta(D11,D11,d) deviates most from its mean. The construction works not because D11 is spectrally flat, but because the fluctuations are correctly aligned with the D11/D22 partition.

---

## 4. Cyclotomic Class Analysis

### Prime m cases (m = 3 mod 4)

For prime m with m = 3 mod 4, the quadratic residues (QR) and non-residues (QNR) each have size (m-1)/2.

**D11 always splits evenly between QR and QNR:**

| n | m | \|D11 cap QR\| | \|D11 cap QNR\| | QR fraction |
|---|---|----------------|-----------------|-------------|
| 6 | 11 | 2 | 2 | 0.500 |
| 10 | 19 | 4 | 4 | 0.500 |
| 12 | 23 | 6 | 6 | 0.500 |
| 16 | 31 | 8 | 8 | 0.500 |
| 22 | 43 | 11 | 11 | 0.500 |

This is **exact** for all prime m cases: |D11 cap QR| = |D11 cap QNR| = |D11|/2.

This makes sense because D11 is symmetric (x in D11 iff -x in D11), and for prime p = 3 mod 4, -1 is a QNR. So negation swaps QR and QNR. Since D11 is closed under negation, it must contain equal numbers of QRs and QNRs.

**D12 does NOT split evenly** -- the QR fraction varies (0.357 to 0.625). D12 is not symmetric, so no such constraint applies.

### Composite m cases

For composite m (n=8,14,18,20), the QR structure is more complex:
- m=15 (3*5): |QR|=5, |QNR|=9 (unequal because Z_15 is not a field)
- m=27 (3^3): |QR|=10, |QNR|=16
- m=35 (5*7): |QR|=11, |QNR|=23
- m=39 (3*13): |QR|=13, |QNR|=25

D11 still has the property that |D11 cap QR| / |D11| ~ |QR|/(m-1), but the exact split depends on the specific composite structure.

### Higher cyclotomic classes

For prime m, I checked 3-class and 6-class cyclotomic decompositions. D11 does NOT concentrate in specific classes -- coverage is roughly proportional to class size. D11 is determined by a combination of conditions beyond simple cyclotomic class membership.

---

## 5. D12 Structure Analysis

### What D12 is NOT:
- D12 is NOT a translate of D11 (D12 != D11 + c for any c)
- D12 is NOT a translate of D22
- D12 is NOT a scalar multiple of D11 (D12 != t*D11 for any unit t)
- D12 is NOT determined by a Delta threshold (D12 != {d : Delta(D11,D11,d) >= k})
- D12 is NOT determined by a Delta+Sigma threshold
- D12 is NOT symmetric

### What D12 IS:
- |D12| = n-1 (always)
- 0 in D12 (always)
- |D12 cap D12^T| varies but D12 shares substantial overlap with its transpose
- |D11 cap D12| ~ |D11|/2 (roughly half of D12 elements are in D11)

### Delta(D11,D11,d) distribution for D12 vs complement

The values Delta(D11,D11,d) for d in D12 and d not in D12 are **interleaved** -- there is no clean threshold that separates them. D12 membership involves BOTH the D11 autocorrelation AND additional structure.

### Key observation for D12

D12 must satisfy the V1V1 constraint (through B(d) = Delta(D12,D12,d)) and the V2V2 constraint (through B(m-d)). This is a non-trivial balancing act. The fact that D12 is not symmetric means B(d) != B(m-d), which provides the extra degrees of freedom needed to satisfy both V1V1 and V2V2 constraints simultaneously.

---

## 6. Multiplier / Automorphism Analysis

For every known construction:
- Only multipliers fixing D11: {1, m-1} (i.e., identity and negation)
- Only multiplier fixing D12: {1} (identity only)
- No non-trivial multiplier fixes both

This means the constructions have **minimal automorphism group** within the circulant structure. They are not built from any obvious "multiplier symmetry" -- the sets are genuinely irregular, chosen to balance the Delta constraints rather than to maximize algebraic structure.

---

## 7. Delta Fluctuation Bounds

### V1V1 block

| n | m | avg | std | range | norm_std (std/sqrt(m)) | gap/std |
|---|---|-----|-----|-------|------------------------|---------|
| 6 | 11 | 3.20 | 0.79 | 2 | 0.238 | -1.01 |
| 8 | 15 | 5.14 | 0.66 | 2 | 0.171 | -1.29 |
| 10 | 19 | 7.11 | 0.58 | 2 | 0.134 | -1.53 |
| 12 | 23 | 11.00 | 1.57 | 4 | 0.328 | 0.64 |
| 14 | 27 | 13.00 | 1.50 | 4 | 0.288 | 0.67 |
| 16 | 31 | 15.00 | 1.89 | 5 | 0.340 | 0.53 |
| 18 | 35 | 17.00 | 1.81 | 5 | 0.306 | 0.55 |
| 20 | 39 | 19.00 | 1.86 | 6 | 0.298 | 0.54 |
| 22 | 43 | 21.00 | 1.74 | 5 | 0.265 | 0.58 |

**Key observations:**

1. **Normalized std ~ 0.26-0.34 * sqrt(m)** for n >= 12. This is consistent with pseudo-random behavior of the difference sets.

2. **gap/std ~ 0.5-0.7** for n >= 12: the red edges need their lambda values to be only about 0.5-0.7 standard deviations below the mean. This is a very mild requirement!

3. For n=6,8,10 (small cases where |D11| < |D22|), the gap is negative, meaning the average is already below the threshold -- it's easy to satisfy the red constraint.

4. **The challenge grows slowly**: gap/std stays around 0.5-0.6 as n increases, which means the problem doesn't get fundamentally harder for larger n (at least for the V1V1 constraint).

### V1V2 block

As proven in Section 1, the V1V2 Delta values are EXACTLY flat: the standard deviation is approximately 0.5 (coming solely from the binary distinction red/blue, not from actual fluctuation).

### Implications for proof strategy

The normalized standard deviation of ~0.3*sqrt(m) suggests that for random symmetric subsets D11 of Z_m with |D11| = (m+1)/2, the autocorrelation Delta(D11,D11,d) fluctuates by O(sqrt(m)) around its mean. The construction only needs these fluctuations to exceed a threshold of about 1 unit, which is < sqrt(m) for large m. This suggests:

**A probabilistic/counting argument might work**: show that among all symmetric D11 with correct size, a positive fraction have Delta fluctuations sufficient to support a valid D12.

---

## 8. Categorization of n <= 50

| Category | Count | Values |
|----------|-------|--------|
| Paley (m = prime power, 1 mod 4) | 15 | 3,5,7,9,13,15,19,21,25,27,31,37,41,45,49 |
| Prime m = 3 mod 4 | 12 | 4,6,10,12,16,22,24,30,34,36,40,42 |
| Prime power m = 3 mod 4 (not prime) | 1 | 14 (m=27=3^3) |
| Composite m | 20 | 8,11,17,18,20,23,26,28,29,32,33,35,38,39,43,44,46,47,48,50 |

### Immediate targets:
- **n=24, m=47** (prime, 3 mod 4): Should be solvable by SA with structural constraints
- **n=23, m=45** (3^2 * 5): Hardest near-term case -- composite, requires new approach
- **n=25, m=49** (7^2, 1 mod 4): Paley construction over GF(49) should work

### Hard composite cases for the general proof:
The composite cases with m = 1 mod 4 are: n=11(m=21), n=17(m=33), n=23(m=45), n=29(m=57), n=33(m=65), n=35(m=69), n=39(m=77), n=43(m=85), n=47(m=93).

The composite cases with m = 3 mod 4 are: n=8(m=15), n=18(m=35), n=20(m=39), n=26(m=51), n=28(m=55), n=32(m=63), n=38(m=75), n=44(m=87), n=46(m=91), n=48(m=95), n=50(m=99).

Note that n=8,18,20 are already solved (by the existing known constructions from Lidicky et al.), so composite m = 3 mod 4 is not an insurmountable obstacle.

---

## 9. Composite m=45 (n=23) Analysis

### Ring structure
Z_45 = Z_9 x Z_5 via CRT.

- **Units**: phi(45) = 24 elements
- **QR mod 45**: 11 elements = {1,4,9,10,16,19,25,31,34,36,40}
- **QR mod 9**: {1,4,7} (3 elements)
- **QR mod 5**: {1,4} (2 elements)

### Multiplier orbits in Z_45
Under multiplication by units:
- Orbit 0: {0} (trivial)
- Orbit 1: 24 elements (all units)
- Orbit 2: {3,6,12,21,24,33,39,42} (8 elements, divisible by 3 but not 5 or 9)
- Orbit 3: {5,10,20,25,35,40} (6 elements, divisible by 5 but not 3)
- Orbit 4: {9,18,27,36} (4 elements, divisible by 9 but not 5)
- Orbit 5: {15,30} (2 elements, divisible by both 3 and 5)

### Symmetry constraints
Since m=45 is odd and has no element equal to its own negation (except 0), D11 must have **even** size.

For the standard construction: |D11| = n = 23 would require an ODD number of elements -- impossible for a symmetric set! So we need either:
- |D11| = 22, |D22| = 22, |D12| = 22 (with d1 = 44, d2 = 44)
- |D11| = 24, |D22| = 20, |D12| = 22 (with d1 = 46 = m+1... check this)

Wait -- the standard pattern for even n has |D11|=n, but n=23 is odd! For the first time we encounter odd n in the "large D11" regime. This means |D11| cannot equal n exactly; it must be n-1=22 or n+1=24.

**This is a fundamentally new challenge for n=23**. The solver should try both |D11|=22 and |D11|=24 configurations.

### Product structure
D11 cannot be a pure product set A x B in Z_9 x Z_5 (no factorization of 22 or 24 fits). It must be a non-product subset, likely a union of cosets or a more complex combinatorial object.

### Recommended approach for solver
1. Try both |D11| = 22 (11 symmetric pairs) and |D11| = 24 (12 pairs)
2. Use the CRT structure: elements at (a,b) in Z_9 x Z_5, try D11 that respects the orbit structure
3. Start with D11 containing entire orbits where possible
4. |D12| = 22 (= n-1), with 0 in D12

---

## 10. Key Constraints for Solver (Actionable Summary)

For any n with m = 2n-1:

1. **D22 = {1,...,m-1} \ D11** (always)
2. **D11 must be symmetric** (always)
3. **|D12| = n-1, 0 in D12** (always)
4. **V1V2 constraints are automatically satisfied** (free, no search needed)
5. **Only V1V1 and V2V2 constraints need to be checked**
6. V1V1(d) = Delta(D11,D11,d) + Delta(D12,D12,d) <= n-2 for d in D11
7. V2V2(d) = A(d) + B(m-d) + (m-2-2|D11|) + 2[d in D11] <= n-2 for d in D22

For n=23 (m=45):
- |D11| must be even: try 22 or 24
- |D12| = 22
- Try multiplier-orbit-aware initialization
- CRT decomposition Z_45 = Z_9 x Z_5 may help structure the search

For n=24 (m=47, prime):
- Standard approach: |D11| = 24 (symmetric), |D12| = 23
- D11 should have 12 QR + 12 QNR (forced by symmetry)
- Use SA with informed initialization

---

## 11. Fourier-Analytic Reformulation (Phase 2)

### Fourier representation of constraints

For a set S in Z_m, define the Fourier coefficient: hat{S}(k) = sum_{s in S} omega^{ks} where omega = e^{2*pi*i/m}.

The autocorrelation Delta(S,S,d) has Fourier representation:
```
Delta(S, S, d) = (1/m) * sum_k |hat{S}(k)|^2 * omega^{-kd}
```

Since D11 is symmetric, hat{D11}(k) is REAL for all k (verified computationally).

### Combined power spectrum

Define P(k) = |hat{D11}(k)|^2 + |hat{D12}(k)|^2. Then:
```
V1V1(d) = A(d) + B(d) = (1/m) * sum_k P(k) * omega^{-kd}
```

The k=0 term gives the average: P(0)/m = (|D11|^2 + |D12|^2)/m.

### Exact values (standard regime, n >= 12)

With |D11| = n, |D12| = n-1, m = 2n-1:
```
P(0)/m = (n^2 + (n-1)^2)/(2n-1) = n - (n-1)/(2n-1) ~ n - 1/2
```

The fluctuation F(d) = V1V1(d) - P(0)/m must satisfy:
```
F(d) <= (n-2) - P(0)/m = -2 + (n-1)/(2n-1) ~ -3/2   for d in D11 (red)
F(d) <= (n+1) - P(0)/m = +1 + (n-1)/(2n-1) ~ +3/2   for d not in D11 (blue)
```

### Remarkable tightness of constructions

For ALL known constructions with n >= 12, the V1V1 red lambda values A(d)+B(d) for d in D11 are:
- Exactly at threshold (n-2) for most red edges
- At most 2 below threshold for a few edges
- NEVER above threshold

Examples:
```
n=12: red values = [9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10], threshold=10
n=14: red values = [11, 11, 12, 12, 12, ..., 12, 12], threshold=12
n=16: red values = [12, 12, 12, 12, 13, 13, 14, 14, ..., 14], threshold=14
n=22: red values = [18, 18, 19, 19, ..., 20, 20, ..., 20], threshold=20
```

This extreme tightness suggests the constructions are essentially unique (up to multiplier equivalence) -- there's very little room to perturb them.

### Negative correlation property

For n >= 12: sum_{d in D11} V1V1(d) < |D11| * avg(V1V1). This means D11 elements have LOWER-than-average A+B values. The red/blue partition is "anti-correlated" with the fluctuation -- exactly what's needed for the constraint.

### Red slack and blue slack

| n | red_slack | blue_slack |
|---|-----------|------------|
| 12 | 2 | 6 |
| 14 | 2 | 8 |
| 16 | 10 | 2 |
| 18 | 10 | 4 |
| 20 | 4 | 12 |
| 22 | 12 | 6 |

Slack = budget - actual sum. Both red and blue constraints have modest slack that grows roughly proportional to n.

---

## 12. Proof Strategy Implications

The V1V2 theorem, V2V2 algebraic identity, and Fourier analysis together suggest:

### Approach A: Existential proof via counting / probabilistic method
For prime m = 3 mod 4:
1. The constraint requires F(d) <= -3/2 for half the elements and F(d) <= +3/2 for the other half
2. The fluctuation std is ~0.3*sqrt(m) ~ 0.3*sqrt(2n) ~ 0.4*sqrt(n)
3. The required deviation is only ~1.5 (constant!), which is ~3.5/sqrt(n) standard deviations
4. As n grows, this becomes EASIER in normalized terms
5. A second-moment or Lovász Local Lemma argument might prove existence for all large n

### Approach B: Explicit construction for composite m
For composite m = p1^a1 * ... * pk^ak:
1. Use CRT: Z_m = Z_{p1^a1} x ... x Z_{pk^ak}
2. Construct D11 component-wise where possible
3. Use "lifting" from solutions for each prime power component
4. Handle the interaction terms (D11 cannot be a pure product)

### Approach C: Fourier-analytic construction
The Fourier formulation reduces the problem to: given |hat{D11}(k)|^2, find D12 with |D12|=n-1 such that |hat{D12}(k)|^2 "corrects" the spectrum:
```
|hat{D11}(k)|^2 + |hat{D12}(k)|^2 has the right Fourier inverse
```
This is a kind of "spectral complementation" problem.

### Approach D: Direct algebraic construction
The fact that the V1V2 theorem holds for ALL D11 and D12 (given complement and size constraints) means we might be able to construct D12 as a function of D11 using some algebraic recipe. The relationship between the V1V1 and V2V2 constraints through D12's "flip" (B(d) vs B(m-d)) suggests that D12 should be "half-symmetric" in some precise sense.

### Most promising path
The Fourier approach (C) combined with the probabilistic bound (A) seems most promising:
1. Prove that for prime p = 3 mod 4, every symmetric D11 with |D11| = (p+1)/2 admits a valid D12
2. For prime powers, use finite field constructions
3. For composite m, show CRT-based constructions work using the orbit structure

---

## 13. Computational Results: New Constructions (Session 2)

### Paley family (q = 1 mod 4)

Proved and verified: For ALL prime powers q = 1 mod 4, the construction
D11 = D12 = QR(GF(q)*), D22 = QNR(GF(q)*) is valid.

Computationally verified for q = 5,9,13,17,25,29,37,41,49,53,61,73,81,89,97.

This covers n = 3,5,7,9,13,15,19,21,25,27,31,37,41,45,49 (and infinitely many more).

Key file: `paley_general.py` (includes theoretical proof and GF(p^k) implementation).

### SA-discovered constructions

| n | m | method | |D11| | |D12| | time | seed |
|---|---|--------|-------|-------|------|------|
| 23 | 45 (=9x5) | fast-joint-SA | 22 | 22 | 74s | 0 |
| 24 | 47 (prime, 3 mod 4) | fast-joint-SA | 22 | 23 | 128s | 0 |
| 26 | 51 (=3x17) | fast-joint-SA | 24 | 25 | 1421s | 2 |

All hit thresholds exactly: max_red = n-2, max_blue = n-1.

### n=23 analysis

The n=23 solution has:
- A(d)+B(d) histogram: {18:2, 19:2, 20:8, 21:14, 22:18}
- 28 of 44 constraints are tight (10 red, 18 blue)
- No clear CRT algebraic pattern; D11 draws roughly evenly from each character class
- Average A+B = 21.00 exactly

### Key insight: Fast SA methodology

The breakthrough was replacing `verify_construction()` (O(m * m^2) per call) with:
1. Precompute Delta arrays: O(m^2) per set
2. Compute cost from Delta arrays: O(m) per evaluation
3. On swap, recompute only the changed Delta arrays: O(m^2) per swap

This makes each SA iteration O(m^2) instead of O(m^3), enabling millions of iterations
in minutes for m ~ 50.

### Coverage summary (n <= 50)

All n <= 31 are now covered:
- Paley proven: n = 3,5,7,9,13,15,19,21,25,27,31 (and infinitely many more)
- SA verified: n = 4,6,8,10,11,12,14,16,17,18,20,22,23,24,26,28,29,30
- Prior work: n <= 21

n=32 (m=63 = 9x7) is the first OPEN case. See Section 14 below.

---

## 14. The n=32 Obstacle: Deep Analysis (Session 3)

### V2V2 = V1V1 Identity (General m)

**Theorem**: For ANY m (prime, prime power, or composite), the V2V2 constraints are algebraically identical to V1V1 when expressed in terms of A(d)+B(d).

**Proof**: B(d) := Delta(D12, D12, d) = |{(a,b) in D12xD12 : a-b = d mod m}|.
Then B(m-d) = |{(a,b) : a-b = -d}| = |{(b,a) : b-a = d}| = B(d).
Also, Delta(D12T, D12T, d) = Delta(D12, D12, m-d) = B(m-d) = B(d).
Therefore V2V2(d) = A(d) + B(d) + constant terms, matching V1V1(d).

**Consequence**: The "8 violations" pattern observed in near-solutions is actually 4 independent violations, each counted twice (once in V1V1, once in V2V2).

### Tight counting constraint (General)

For ALL n >= 3, the average A(d)+B(d) for d=1..m-1 equals approximately the binding threshold:

```
avg = (|D11|^2 - |D11| + |D12|^2 - |D12|) / (m-1)
```

With the "natural" |D11| ~ (m-1)/2 and |D12| = n-1 = (m-1)/2:
```
avg ~ ((m-1)/2)^2 * 2 / (m-1) = (m-1)/2 ~ n-1
```

The binding threshold is n-2, so avg exceeds threshold by ~1. This is universal.

### n=32 Specific Findings

For m = 63 = 9 x 7:

1. **|D11|=32**: avg A+B = 31.0, threshold = 30, slack = -1.0
2. **|D11|=30**: avg A+B = 29.03, threshold = 29 (blue), slack = -0.03
3. **|D11|=28,34**: Much worse (slack -2.2 to -3.1)

Near-solutions (cost=4, i.e., 2 real violations):
- All violations have excess exactly 1
- Violated differences form complementary pairs {d, m-d}
- Exhaustive D12 fix always fails (D11 is the root cause)
- 100+ seeds across 5 solver variants all converge to same barrier

CRT decomposition (Z_63 = Z_9 x Z_7) classes:
- Negation maps QR <-> QNR in both Z_9* and Z_7*
- Symmetric pairs: RR<->NN (18 elements), RN<->NR (18), R0<->N0 (6), 0R<->0N (6), ZR<->ZN (12), Z0<->Z0 (2)
- All CRT-aware algebraic initializations converge to same cost=4 barrier

**Open question**: Does a valid 2-block circulant construction exist for m=63? The counting argument says the budget is sufficient (capacity 1950 vs need 1922), but the autocorrelation structure may prevent achieving cost=0.
