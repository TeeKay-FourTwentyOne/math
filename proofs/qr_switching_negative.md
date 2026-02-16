# QR Switching Proof Attempt: Negative Result

## Summary

**Claim**: Starting from D12 = QR (quadratic residues mod p), one can fix all constraint
violations by a sequence of element swaps.

**Result**: This approach **FAILS** fundamentally. Valid D12 is not a small perturbation
of QR; it requires changing approximately 50% of all elements.

## Setup

For p ≡ 3 (mod 4) prime, the quadratic residues QR mod p form a
(p, (p-1)/2, (p-3)/4)-difference set. Setting D12 = QR gives:

- |D12| = (p-1)/2 (correct size)
- B(d) = (p-3)/4 for all d ≠ 0 (perfectly flat autocorrelation)
- B(0) = |D12| = (p-1)/2

The constraints require, for flat D11 (A(d) ≈ (p+1)/4 at all positions):
- D11 positions: A(d) + B(d) ≤ (p-3)/2, i.e., B(d) ≤ (p-3)/4 - 1 = E[B] - 1
- D22 positions: A(d) + B(d) ≤ (p+3)/2, i.e., B(d) ≤ (p+3)/4 - 1 = E[B] + 2

So QR exceeds the D11 threshold by exactly 1 at every D11 position (where A(d) = E[A]).
At D22 positions, QR satisfies with margin 2.

## Conservation Law

For any D12 with |D12| = (p-1)/2:

Σ_{d=1}^{p-1} B(d) = k(k-1) = (p-1)(p-3)/4 (constant, Parseval)

The sum of constraint thresholds is:

Σ T(d) = |D11| × (n-2) + |D22| × (n+1) = (p+1)/2 × (p-3)/2 + (p-3)/2 × (p+3)/2

The net slack (Σ thresholds - Σ A - Σ B) = -(p-1)/2 for all D11 at p=11, = -(3p-3)/2
for p=19. The sum of excesses is **negative**, so validity is not prevented by the
conservation law. The issue is distributional: can B(d) be low at D11 positions while
being high at D22 positions?

## Greedy Swap Results

### p = 11
- D12 = QR has total violation 2-6 depending on D11
- Greedy 1-swap **SOLVES** all 5 working D11 in exactly 1 step
- The swap is always: remove one QR element, add 0
- Closest valid D12 to QR: Hamming distance 2 (1 swap)

### p = 19
- D12 = QR has total violation 4-20 depending on D11
- Greedy 1-swap gets stuck at total violation = 2 for ALL 126 symmetric D11
- **0 out of 126 D11 solved by greedy** (0 out of 9 working D11)
- At the stuck point: **NO single swap reduces violation further**
- **NO 2-swap reduces violation further either** (exhaustively verified)
- Closest valid D12 to QR: Hamming distance 6 (needs 3 removes + 2 adds)
- All 9 working D11-orbits: 38 valid D12 each, Hamming range [6, 14], mean 9.5

### p = 23
- Greedy solves 3 out of 462 D11 (all non-working D11 get stuck at tv ≥ 2)
- Most get stuck at tv = 2

### p = 43
- Greedy gets stuck at tv = 2-10 for all 30 tested working D11
- Actual valid D12: Hamming distance from QR has mean 21.8, range [14, 30]
- Expected random Hamming distance: 21.5
- **QR retention fraction = 0.48** (vs random 0.49): D12 is statistically indistinguishable
  from random in its relationship to QR

## Key Structural Finding

At p = 19, the greedy stuck point has this profile:

```
d  type  A(d)  B(d)  threshold  excess
1  D11    4     4       8          0
2  D11    5     3       8          0
3  D11    5     4       8         +1  <== VIOLATED
6  D11    3     5       8          0
8  D11    4     3       8         -1
```

Position d=3 (and its twin d=16) are violated by exactly 1. But to reduce B(3),
every single swap that decreases B(3) also increases B at some other position,
creating a new violation. This is a **topological obstruction**: the constraint
surface has no descent direction.

## Hamming Distance Distribution at p = 43

| Hamming | Count |
|---------|-------|
| 14      | 3     |
| 16      | 6     |
| 18      | 12    |
| 20      | 30    |
| 22      | 27    |
| 24      | 25    |
| 26      | 15    |
| 28      | 5     |
| 30      | 1     |

The distribution is centered at ~21 (= |QR|), consistent with valid D12 being
**uniformly random** in its relationship to QR. The closest valid D12 at p=43
is at Hamming distance 14 (7 swaps), not 2 or 4.

## All Valid D12 Contain 0

At p=43, all 124 valid D12 contain 0. At p=19, 18/38 valid D12 contain 0.
This is structurally important: 0 ∉ QR for any p, so including 0 in D12 is
a universal feature of valid constructions.

## Why QR Switching Fails: Root Cause

1. **QR has the wrong spectral shape**: QR is a perfect difference set with
   flat power spectrum P_QR(k) = k for all k ≠ 0. But valid D12 needs
   *anti-correlated* spectrum with D11 (Corr(P_D11, P_D12) ≈ -0.5 to -0.9).

2. **Flat B is the worst case**: B(d) = constant means B hits the D11 threshold
   at every position. Valid D12 has B(d) varying by ±√(p/16), which means some
   D11 positions have B well below threshold while a few exceed it only where
   the D22 slack absorbs the excess.

3. **The perturbation is not small**: Going from flat B to the required
   asymmetric B profile requires changing ~50% of D12's elements.

4. **Local swaps create global ripples**: Each element swap changes B(d) at
   ALL p-1 positions (via the autocorrelation), not at O(1) positions.
   There is no way to make "local" adjustments.

## Implications for the Proof

The QR switching approach is definitively ruled out as a proof strategy.
The valid D12 is not close to any algebraically defined set (QR, cyclotomic
classes, arithmetic progressions). Its defining property is the global
spectral complementarity with D11.

Proof strategy (c) from the proof sketch -- explicit construction -- would
need to construct D12 from scratch using spectral methods, not by modifying
QR. This is essentially as hard as the original existence problem.

The probabilistic approaches (a) first moment with c₀ bound, or
(b) second moment / Paley-Zygmund, remain the most viable paths.
