# Exhaustive D12 Enumeration Results

**Date**: 2026-02-10
**Purpose**: Verify existence of valid D12 for symmetric D11 at small primes.

---

## Summary

| p | n | #D11 | #D12/D11 | Working D11 | Total valid pairs | Min ratio | Time |
|---|---|------|----------|-------------|-------------------|-----------|------|
| 11 | 6 | 10 | 210 | 5 (50%) | 100 | 2.674 | 0.0s |
| 19 | 10 | 126 | 43,758 | 9 (7.1%) | 162 | 0.217 | 0.4s |
| 23 | 12 | 462 | 646,646 | 55 (11.9%) | 4,356 | 0.088 | 33.4s |

**Key finding**: NOT all symmetric D11 have a valid D12. However, a substantial fraction do, and the total number of valid (D11, D12) pairs grows rapidly with p. The proof only needs existence of at least one valid pair per prime, which is confirmed for all three primes.

---

## Critical Structural Observation: A-flatness

The single most predictive feature of whether a D11 admits a valid D12 is the **maximum autocorrelation** max_{d in D11} A(d):

| p | Working D11: max A(d) | Non-working D11: max A(d) |
|---|----------------------|--------------------------|
| 11 | {3} | {4} |
| 19 | {5} | {5, 6, 7, 8} |
| 23 | {6, 7} | {6, 7, 8, 9, 10} |

For p=11, the separation is perfect: max_A = 3 works, max_A = 4 does not.
For p=19, working D11 all have max_A = 5 (the minimum possible).
For p=23, working D11 have max_A in {6, 7}; all with max_A >= 8 fail.

**Interpretation**: The binding threshold is (p-3)/2. Each constraint needs A(d) + B(d) <= (p-3)/2. Since E[B(d)] = (p-3)/4, the slack for B(d) is (p-3)/2 - A(d) - (p-3)/4 = (p-3)/4 - A(d). Large A(d) leaves no room for B(d) fluctuations.

---

## QR Structure

For p = 3 mod 4, the quadratic residues QR are NOT symmetric (since -1 is a QNR). Therefore QR cannot serve directly as D11.

All working D11 have a **perfect 50/50 split** between QR and QNR elements:
- p=11: all working D11 have |D11 ∩ QR| = |D11 ∩ QNR| = 3
- p=19: all working D11 have |D11 ∩ QR| = |D11 ∩ QNR| = 5
- p=23: all working D11 have |D11 ∩ QR| = |D11 ∩ QNR| = 6

This is consistent with maximizing spectral flatness: mixing QR and QNR elements evenly distributes the autocorrelation mass.

---

## Joint vs Product of Marginals

The "ratio" column gives P[all constraints satisfied jointly] / Product of P[each constraint satisfied].

- p=11: ratio = 2.674 for all working D11 (strong positive association)
- p=19: ratio = 0.217 for all working D11 (negative association, but joint event still occurs)
- p=23: ratio ranges from 0.088 to 0.278 across working D11

**For the proof**: The Slepian-based argument requires ratio >= p^{-C} for some constant C. Even the worst ratio 0.088 at p=23 far exceeds 23^{-C} for any reasonable C. The polynomial loss claimed in L6 is amply satisfied.

---

## S1 Statistics for Working D11

For each working D11, S1 = sum_{d in D11} B(d) among valid D12:

### p = 11 (threshold per d: 4, so S1 <= 24 needed across 6 positions)
- All working D11: S1 mean ~ 12.6, range [8, 16]

### p = 19 (threshold per d: 8, so S1 <= 80 needed across 10 positions)
- All working D11: S1 mean ~ 39.0, range [33, 45]

### p = 23 (threshold per d: 10, so S1 <= 120 needed across 12 positions)
- Working D11 with 198 valid D12: S1 mean ~ 54-58, range varies

---

## Constraints

For each (D11, D12) pair with m = p, |D11| = (p+1)/2, |D12| = (p-1)/2, 0 in D12:
- **Binding** (d in D11): A(d) + B(d) <= (p-3)/2
- **Loose** (d in D22):   A(d) + B(d) <= (p+3)/2

where A(d) = Delta(D11,D11,d), B(d) = Delta(D12,D12,d).

---

## Per-Prime Details

### p = 11 (n = 6)

10 symmetric D11 (choosing 3 pairs from 5 available). 5 work, 5 don't.

| D11 | #Valid D12 | max A(d) | A values at D11 positions |
|-----|------------|----------|--------------------------|
| [1,2,3,8,9,10] | 0 | 4 | {1:4, 2:2, 3:2, 8:2, 9:2, 10:4} |
| **[1,2,4,7,9,10]** | **20** | **3** | {1:2, 2:3, 4:1, 7:1, 9:3, 10:2} |
| **[1,2,5,6,9,10]** | **20** | **3** | {1:3, 2:1, 5:2, 6:2, 9:1, 10:3} |
| [1,3,4,7,8,10] | 0 | 4 | {1:1, 3:4, 4:2, 7:2, 8:4, 10:1} |
| **[1,3,5,6,8,10]** | **20** | **3** | {1:1, 3:2, 5:3, 6:3, 8:2, 10:1} |
| [1,4,5,6,7,10] | 0 | 4 | {1:2, 4:4, 5:1, 6:1, 7:4, 10:2} |
| **[2,3,4,7,8,9]** | **20** | **3** | {2:2, 3:1, 4:3, 7:3, 8:1, 9:2} |
| [2,3,5,6,8,9] | 0 | 4 | {2:1, 3:2, 5:4, 6:4, 8:2, 9:1} |
| [2,4,5,6,7,9] | 0 | 4 | {2:4, 4:2, 5:2, 6:2, 7:2, 9:4} |
| **[3,4,5,6,7,8]** | **20** | **3** | {3:3, 4:2, 5:1, 6:1, 7:2, 8:3} |

### p = 19 (n = 10)

126 symmetric D11 (choosing 5 pairs from 9). 9 work, 117 don't.
All working D11 have max_A = 5 (the minimum possible for |D11|=10 in Z_19).
Each working D11 has exactly 18 valid D12.

### p = 23 (n = 12)

462 symmetric D11 (choosing 6 pairs from 11). 55 work, 407 don't.
Working D11 have max_A in {6, 7}. Valid D12 counts: 22, 44, 110, or 198.

---

## Implications for the Proof

1. **L6 does NOT hold in the "for all D11" form**. Many symmetric D11 have zero valid D12.

2. **L6 DOES hold in the "there exists" form**: For each tested prime, multiple D11 have valid D12, confirming R(B_{n-1}, B_n) >= 4n-1 by construction.

3. **The first moment argument must be refined**: Rather than showing "for random D12, E[valid] > 0 for ALL D11", the proof should show "for random (D11, D12), E[#valid pairs] > 0". The growing count of valid pairs (100, 162, 4356 for p=11,19,23) supports this.

4. **A-flatness is the correct D11 selection criterion**: D11 with flat autocorrelation (small max A(d)) are the ones that work. This connects to the spectral flatness discussion in the proof outline.

---

## Complement Symmetry Verification (Theorem 5)

Exhaustive enumeration confirms the complement bijection between |D11|=n and |D11|=n-2:

| p | |D11|=n: working / total / valid | |D11|=n-2: working / total / valid |
|---|--------------------------------|-----------------------------------|
| 11 | 5/10 / 100 | 5/10 / 100 |
| 19 | 9/126 / 162 | 9/126 / 162 |
| 23 | 55/462 / 4356 | 55/462 / 4356 |

Counts are **identical** in both formulations. The bijection maps (D11, D12) to (complement(D11), negate(D12)).

Thresholds differ:
- |D11|=n: D11 thresh=(p-3)/2 (binding), D22 thresh=(p+3)/2 (loose, 3 units above)
- |D11|=n-2: D11 thresh=(p-3)/2 (binding), D22 thresh=(p-5)/2 (also binding, 1 unit below D11)

The |D11|=n-2 formulation makes both constraints tight, revealing the true bottleneck structure.

---

## Enhanced Sampling for Larger Primes (p=31,43,47,59)

### D11 Selection Strategy Comparison

Three strategies tested with 1M D12 samples per D11:

| Strategy | Description | p=31 full hits |
|----------|------------|---------------|
| Greedy D11-only | Minimize max A(d) over d in D11 | 1/5 D11 |
| Greedy global | Minimize max A(d) over ALL d=1..p-1 | 0/10 D11 |
| QR-balanced | Random among |D11 ∩ QR| = |D11 ∩ QNR| | 4/10 D11 |

**Key finding**: For p=31, we found ~35 valid (D11, D12) pairs by random sampling. The winning D11 have max_A(D11) = 8 AND max_A(all) <= 11 AND QR-balanced.

Neither extremal strategy is optimal:
- Greedy D11-only achieves max_A(D11) = 8 but max_A(all) = 15-29 (D22 bottleneck)
- Greedy global achieves max_A(all) = 9-14 but max_A(D11) = 9-18 (binding fails)

### Binding Ratio Growth (Positive Association)

For known-good D11 from SA solutions:

| p | Binding rate | Binding ratio | log2(E[bind-valid]) |
|---|-------------|--------------|-------------------|
| 31 | 1.94% | 11.5 | 21.4 |
| 43 | 0.026% | 47.0 | 27.0 |
| 47 | 0.010% | 69.4 | 29.6 |
| 59 | 0.018% | 118.4 | 42.2 |

The binding ratio grows as ~p^1.5 and E[binding-valid] grows as ~2^(0.75p).

### D22 Bottleneck Analysis

| p | max_A(D22) | D22 thresh | E[B] | slack (thresh - maxA - E[B]) |
|---|-----------|-----------|------|-----|
| 31 | 12 | 17 | 7.2 | -2.2 |
| 43 | 14 | 23 | 10.2 | -1.2 |
| 47 | 17 | 25 | 11.2 | -3.2 |
| 59 | 21 | 31 | 14.3 | -4.3 |

D22 slack is negative: the tightest D22 position has A(d) + E[B] > threshold. The D22 constraints pass only when B(d) falls in the left tail, which is helped by positive association from the binding constraints.

### Gaussian First-Moment Estimates

| p | Best Gaussian log(E[valid]) | E[valid] |
|---|---------------------------|----------|
| 31 | -4.30 | 0.014 |
| 43 | -6.02 | 0.002 |
| 47 | -6.60 | 0.001 |
| 59 | -8.49 | 0.0002 |

The Gaussian (independence) estimate is pessimistic - it predicts E[valid] < 1 for all p >= 31. The actual positive association ratio (growing as p^1.5) compensates, but the proof must capture this.

---

## Raw Data

Full per-D11 data saved to:
- `enumeration_data/enumeration_p{11,19,23}.json` (exhaustive)
- `enumeration_data/sampling_aflat_p{31,43,47,59}.json` (sampling with A-flat/QR-balanced D11)
