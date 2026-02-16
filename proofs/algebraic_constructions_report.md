# Algebraic/Cyclotomic Constructions for p ≡ 3 mod 4

**Date**: 2026-02-14
**Purpose**: Systematically explore explicit algebraic constructions of (D11, D12) for primes p ≡ 3 mod 4.

---

## 1. Framework: Coset Decomposition

**Key structural observation**: Every symmetric D11 of size (p+1)/2 in Z_p can be uniquely described as a union of (p+1)/4 cosets from the (p-1)/2-th power residue decomposition.

Let g be a primitive root mod p. Define:
- e = (p-1)/2 (odd, since p ≡ 3 mod 4)
- C_j = {g^j, g^{j+e}} = {g^j, -g^j} for j = 0, 1, ..., e-1

Each C_j has size 2 and is self-symmetric (since g^{j+e} = -g^j). There are e = (p-1)/2 such cosets, partitioning {1, ..., p-1}.

**Any symmetric D11 with |D11| = (p+1)/2 is a union of (p+1)/4 such cosets.**

The number of D11 candidates is C(e, (p+1)/4), and multiplication by g cyclically shifts coset indices, giving orbits of size e (generically).

| p | e = (p-1)/2 | target = (p+1)/4 | C(e, target) | #orbits | #working | working fraction |
|---|---|---|---|---|---|---|
| 11 | 5 | 3 | 10 | 2 | 1 | 50.0% |
| 19 | 9 | 5 | 126 | 14 | 1 | 7.1% |
| 23 | 11 | 6 | 462 | 42 | 5 | 11.9% |
| 31 | 15 | 8 | 6,435 | 429 | 18 | 4.2% |
| 43 | 21 | 11 | 352,716 | 16,796 | 124 | 0.74% |

---

## 2. Constructions Tested

### 2.1 Consecutive Cosets: D11 = C_0 ∪ C_1 ∪ ... ∪ C_{k-1}

**Result**: Works ONLY at p=11. Fails at p=19,23,31,43 with rapidly growing violation cost.

| p | SA cost |
|---|---------|
| 11 | 0 (VALID) |
| 19 | 8 |
| 23 | 16 |
| 31 | 48 |
| 43 | 116 |

**Reason**: Consecutive cosets in the g-power ordering correspond to consecutive powers of g, which cluster in Z_p. This gives poor spectral properties (high autocorrelation peaks).

### 2.2 Center-Based: D11 = elements closest to p/2

**Result**: Same as consecutive cosets (equivalent construction). Only works at p=11.

### 2.3 QR-Guided Pair Selection

Seven strategies tested (QR-low, QR-high, smallest QR element, etc.).

**Result at p=11**: The "closest to center" strategy works (matches the working orbit). All others fail.
**Result at p≥19**: No QR-guided strategy consistently produces working D11.

### 2.4 Symmetrized QR Union/Intersection

- QR ∩ neg(QR) = empty set for all p ≡ 3 mod 4 (since χ(-1) = -1)
- QR ∪ neg(QR) = Z_p* (too large, needs trimming)
- Various trimming strategies tested. No consistent working construction found.

### 2.5 Power Residue Coset Unions

For each divisor e | (p-1), tried all possible unions of negation-paired cosets of e-th power residues that give the correct size.

**Result**: At p=11, found 20 working constructions from 5th power residue cosets (all in the same orbit). At p=19, found all 9 working constructions from 9th power residue cosets (all in same orbit). These match the general coset framework from Section 1.

### 2.6 Half-System Criteria

Tested criteria based on:
- Products of consecutive Legendre symbols: χ(d)·χ(d+1) = ±1
- Linear conditions: ad + b mod p in a half-interval
- Various other algebraic criteria

**Result**: Some produce valid D11 at individual primes, but NO criterion works across multiple primes.

### 2.7 Index-Based Algebraic Conditions

Tried to identify algebraic properties of the working coset INDEX set (viewed as a subset of Z_e):
- Difference representation (autocorrelation of index set)
- Gap structure (circular gaps between consecutive indices)
- Modular residue patterns

**Result**: No simple algebraic condition on the index set separates working from non-working orbits at p ≥ 23. Specifically:
- Most-flat difference representation (range=0) does NOT imply working (p=23: both range=0 orbits are non-working)
- Gap variance does not separate (working orbits span the full range of gap variances)
- No modular pattern found

---

## 3. Structural Analysis of Known Solutions

### 3.1 QR Balance

All working D11 have a **perfect 50/50 split** between QR and QNR elements:
- |D11 ∩ QR| = |D11 ∩ QNR| = (p+1)/4

This is automatic: since D11 is symmetric and χ(-1) = -1 for p ≡ 3 mod 4, every d ∈ QR ∩ D11 is paired with p-d ∈ QNR ∩ D11.

### 3.2 DFT Flatness

At p=11: working D11 have flatness 4.115, non-working have 1.917.
At larger p: flatness alone does not separate (threshold shifts with p).

### 3.3 Multiplicative Orbit Structure (Already Known)

N(D11) is constant on orbits under multiplication by the primitive root g. The orbit structure is:
- p=11: 2 orbits of size 5 each
- p=19: 14 orbits of size 9 each
- p=23: 42 orbits of size 11 each
- p=31: 429 orbits of size 15 each
- p=43: 16,796 orbits of size 21 each

### 3.4 Known Solution Index Patterns

For the p=43 known working solution:
- D11 indices (log base g): show periodic structure with period e=(p-1)/2
- Index differences: [3,1,2,2,1,1,1,5,1,2,2,3,1,2,2,1,1,1,5,1,2] (period 21)
- Evenly balanced across all small moduli

For the p=47 solution: similar balanced structure.

---

## 4. Key Negative Results

1. **No explicit algebraic construction works across all primes p ≡ 3 mod 4**. Every construction we tested either:
   - Works only at p=11 (too small to be representative)
   - Produces non-working D11 at some prime

2. **No simple algebraic invariant separates working from non-working orbits** at p ≥ 23. This was already known from the exhaustive invariant search, and our coset-level analysis confirms it extends to the coset index structure.

3. **The working orbit fraction p_working is declining faster than 1/p**: p × p_working goes 5.5, 1.36, 2.74, 1.30, 0.32 for p = 11, 19, 23, 31, 43. This means a random coset selection has a vanishingly small probability of being working.

---

## 5. Positive Observations

1. **The coset framework provides a clean parametrization**: D11 construction reduces to choosing (p+1)/4 indices from {0,...,(p-3)/2} in Z_{(p-1)/2}. This is a compact, well-structured search space.

2. **Working solutions exist at every tested prime**: Despite the declining fraction, the absolute count of working orbits grows (1, 1, 5, 18, 124 for p = 11, 19, 23, 31, 43). The total N = Σ N(orbit) also grows.

3. **The QR balance and spectral complementarity constraints are necessary but not sufficient**: They prune the search space but don't identify the working constructions.

4. **The total number of valid (D11, D12) pairs grows rapidly**: 100, 162, 4356, and much more at larger p. This growth is what the first moment proof needs to capture.

---

## 6. Implications for L6 Proof

**An explicit algebraic construction approach to L6 appears unviable.** No algebraic formula for D11 works across primes, and the working orbits lack a consistent algebraic characterization.

The proof must instead use:
1. **Probabilistic/counting arguments**: Show that E[N] > 0 over random D11 (or random coset selections)
2. **Spectral/moment bounds**: Use the DFT structure (Fourier analysis on Z_p) to bound constraint correlations
3. **Asymptotic arguments**: Show that the total count of valid pairs grows, even though the fraction of working orbits shrinks

The coset framework may still be useful as a structured domain for probabilistic arguments (e.g., random union of cosets has well-controlled DFT properties).

---

## 7. Files

- Construction code: `ramsey-book-graphs/algebraic_constructions.py`
- This report: `proofs/algebraic_constructions_report.md`
