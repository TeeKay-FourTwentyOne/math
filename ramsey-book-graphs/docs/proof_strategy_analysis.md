# Proof Strategy Analysis: R(B_{n-1}, B_n) = 4n-1

**Date**: 2026-02-10
**Status**: Framework established, key obstacles identified

---

## 1. What We Have

### Complete structural reduction (Theorems 1-3)

The problem of showing R(B_{n-1}, B_n) >= 4n-1 via 2-block circulants reduces to:

**Find symmetric D11 ⊂ {1,...,m-1} and D12 ⊂ Z_m with |D12| = n-1 such that:**

```
A(d) + B(d) ≤ n-2         for all d ∈ D11     (binding)
A(d) + B(d) ≤ 2|D11|-n+1  for all d ∈ D22     (loose)
```

where A(d) = C(D11,D11;d), B(d) = C(D12,D12;d), and D22 = {1,...,m-1}\D11.

V1V2 edges are automatically satisfied. V2V2 constraints are algebraically identical to V1V1.

### Paley case: DONE

For m = prime power ≡ 1 mod 4: D11 = D12 = QR(GF(m)) works. Proven and verified for all q ≤ 100.

### Computational verification

Valid constructions found via SA for all n ≤ 31, including:
- 9 primes p ≡ 3 mod 4: p = 7, 11, 19, 23, 31, 43, 47, 59
- 8 composites: m = 15, 21, 33, 35, 39, 45, 51, 55, 57

---

## 2. The Central Obstacle

### Spectral compensation is essential

**Key finding**: Setting D12 = QR(p) gives perfectly constant B(d) = (p-3)/4 for p ≡ 3 mod 4. But this FAILS because:

1. With constant B, the A(d) constraint becomes A(d) ≤ avg_A - 1 at D11 positions
2. The autocorrelation of symmetric sets is too concentrated — A(d) at D11 positions frequently exceeds avg_A
3. Even when A at D11 passes, A at D22 overflows (pushing excess from D11 to D22)

**What actually works**: The known solutions have B(d) varying (std ~1.2-1.5), with B(d) negatively correlated with A(d) at D11 positions:
- Correlation: ρ = -0.81 (p=43), -0.78 (p=47), -0.94 (p=59)
- Where A peaks, B dips, keeping A+B at threshold

This is **spectral compensation**: D12's Fourier spectrum |D̂12(k)|² is negatively correlated with D11's |D̂11(k)|².

### Why simple algebraic constructions fail

Tested for all p ≡ 3 mod 4 up to 83:
- Interval in dlog space: FAILS at every starting position
- Biquadratic residues: FAILS
- Alternating even dlogs: FAILS

The max A(d) at D11 positions exceeds threshold by ~40%. The autocorrelation of sets with algebraic structure in Z_p tends to have large peaks.

### LP-feasible D11 exist but are irregular

Exhaustive search (p=23, 31) and sampling (p=43) shows:
- 7.1% of random symmetric D11 are LP-feasible at p=23
- 3.7% at p=31
- ~4.5% at p=43

These D11 have no obvious algebraic structure in the cyclotomic framework.

---

## 3. Viable Proof Strategies

### Strategy A: Fourier LP + spectral realization

**Approach**: Show that for any p ≡ 3 mod 4, there exist D11 and D12 such that:
1. The Fourier LP is feasible (a target spectrum exists)
2. A binary set D12 of size n-1 can realize (or approximate) that spectrum

**Current status**:
- Step 1 verified computationally for all tested primes
- Step 2 is the key open problem: "spectral realization of binary sequences"

**Difficulty**: The spectral realization problem asks whether a nonneg function on Z_m can be the power spectrum of a (0,1)-sequence of given weight. This is related to the "turnstile problem" and is NP-hard in general, but may have structure-specific solutions for our parameters.

### Strategy B: Probabilistic existence (joint D11, D12)

**Approach**: Choose both D11 and D12 from some distribution and show P[all constraints satisfied] > 0.

**Key statistics**:
- E[A(d) + B(d)] = (|D11|² + |D12|² - |D11| - |D12|) / (m-1) for d ≠ 0
- The binding threshold is exactly E[A+B] - 1 (one unit below average)
- Need ALL ~n/2 D11 positions to have A+B below average

**Difficulty**: The events {A(d)+B(d) > threshold} for different d are positively correlated (higher-order, Phase 3 analysis: 50-133× excess over independence). Standard union bound / LLL both fail by large margins.

**Possible approach**: Use a conditioned second moment. Condition on D11 being LP-feasible (probability ~5%), then show that D12 satisfying the constraints exists with probability bounded away from 0 given this D11.

### Strategy C: Character sum approach

**Approach**: For p ≡ 3 mod 4, construct D11 using specific character-theoretic properties and bound A(d) + B(d) using Weil bounds.

**Structure**:
1. D11 = union of symmetric pairs selected by a character condition
2. D12 = some modification of QR (e.g., QR ∪ {0} with one element removed)
3. Bound A(d)+B(d) using: A(d) = (1/m)Σ|D̂11(k)|²ω^{kd}, B(d) ≈ (p-3)/4 + perturbation

**Difficulty**: The perturbation from modifying QR to get 0 ∈ D12 creates O(1) changes to B(d), which may be enough. But bounding max_d A(d) requires bounding max_k |D̂11(k)|², which for random sets is O(p log p) (too large).

### Strategy D: Composite case via product construction

**Approach**: For composite m = ab, use CRT Z_m ≅ Z_a × Z_b.
Build D11, D12 as products of sets in the components.
The autocorrelation decomposes: C(A×B, A×B; (d_a, d_b)) = C(A,A;d_a) · C(B,B;d_b).

**Difficulty**: Direct products give multiplicative autocorrelation structure, which doesn't match the additive threshold. But tensor constructions over product groups might work for specific factorizations.

---

## 4. Recommended Next Steps (prioritized for proof progress)

### Priority 1: Conditioned existence argument

For p ≡ 3 mod 4:
1. Fix D11 to be any LP-feasible D11 (existence verified computationally)
2. Show that for this D11, a random D12 of size n-1 satisfies constraints with positive probability
3. Key tool: the Fourier LP guarantees a feasible spectrum exists; need to show a random set approximates it

This avoids the hardest part (constructing D11) by using existence + concentration.

### Priority 2: D12 near QR(p)

D12 = QR gives constant B = (p-3)/4. A small perturbation (swap one element) changes B(d) by O(1) at specific shifts. Can we:
1. Start with D12 = QR
2. Apply a bounded number of swaps to create the needed compensation
3. Bound the resulting B(d) using character sums

This is more tractable than the full construction problem.

### Priority 3: Understanding the LP-feasible D11 structure

Why do ~5% of random D11 pass the LP? Is there a clean characterization?
- Test: is LP-feasibility related to the flatness of |D̂11(k)|²?
- Test: do LP-feasible D11 form a connected component under pair-swaps?

---

## 5. Key Quantities

| Quantity | Formula | Value for |D11|=n |
|----------|---------|-------------------|
| Average A+B | (|D11|²+|D12|²-|D11|-|D12|)/(m-1) | (n²-n + (n-1)²-(n-1))/(2n-2) ≈ n-1 |
| D11 threshold | n-2 | avg - 1 |
| D22 threshold | 2n-n+1 = n+1 | avg + 2 |
| Total slack | n-4 | grows linearly |
| B(d) if D12=QR | (p-3)/4 | constant |
| A(d) average | (p+1)/4 | = n/2 |
| A(d) D11 threshold (QR case) | (p-3)/4 | = avg_A - 1 |

---

## 6. Connection to Known Open Problems

The spectral realization question connects to:
1. **Turán power sum problem**: Minimizing max |Σ z^n_k| over unit-modulus sequences
2. **Circulant Hadamard conjecture**: Does a circulant Hadamard matrix of order > 4 exist?
3. **Flat polynomial problem**: For which n does there exist a polynomial Σ a_k z^k with a_k ∈ {0,1} and nearly flat |P(z)|² on the unit circle?

The connection to (3) is closest: we need D12 such that |D̂12(k)|² "fills in the gaps" of |D̂11(k)|² to create a nearly flat total power spectrum P(k).
