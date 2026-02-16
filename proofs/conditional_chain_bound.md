# Conditional Chain Bound for c₀_reps

## Goal
Prove: for flat D11, log₂(c₀_reps) ≥ -Cp for some C < 1/2.

Combined with headroom_reps ≈ p/2 + O(√p), this gives:
log₂(E[N]) = headroom_reps + log₂(c₀_reps) ≥ (1/2 - C)p + O(√p) → ∞

## Empirical Data (flattest working D11 at each prime)

| p  | N    | budget | cost_r | HR_r  | HR_f  | lg₂c₀r | lg₂c₀f | lg₂E[N] |
|----|------|--------|--------|-------|-------|---------|---------|---------|
| 11 | 44   | 8.85   | 2.88   | 5.97  | 3.09  | -0.51   | 2.37    | 5.46    |
| 19 | 38   | 16.50  | 6.31   | 10.19 | 3.88  | -4.94   | 1.37    | 5.25    |
| 23 | 46   | 20.37  | 7.41   | 12.95 | 5.54  | -7.43   | -0.02   | 5.52    |
| 31 | 1426 | 28.16  | 10.96  | 17.20 | 6.24  | -6.73   | 4.23    | 10.48   |

**Linear fit: log₂(c₀_reps) ≈ -0.323×p + 1.87**
**Linear fit: log₂(E[N]) ≈ 0.244×p + 1.55**

Key ratio: |log₂(c₀_reps)|/p ≈ 0.05, 0.26, 0.32, 0.22 — non-monotone but < 0.5.

## Equicorrelation Structure

The B(d) variables for d ∈ reps = {1,...,(p-1)/2} satisfy:
- E[B(d)] = μ = k(k-1)/(p-1) = (p-3)/4
- Var(B(d)) = σ² = (p-3)(p+1)/(16(p-2)) ≈ p/16
- Cov(B(d), B(d')) = -σ²ρ where ρ = -2/(p-3) (for d' ∉ {d, p-d})
- Σ_{reps} B(d) = k(k-1)/2 = constant (exact sum constraint)

The (m-1)-dimensional conditional distribution (given the sum) has ρ = -2/(p-3) ≈ -2/p.

## Conditional Chain Decomposition

c₀_reps = ∏_{j=1}^{m} r_j where r_j = Pr[E_j | E₁∩...∩E_{j-1}] / Pr[E_j]

### Gaussian Heuristic

Under equicorrelated Gaussian with ρ < 0, conditioning on E₁,...,E_{j-1}:
- Shifts E[B(d_j)] upward by ≈ j × |ρ| × σ × √(2/π)
  (because each conditioning event removes the lower tail)
- The shift in z-score: Δz ≈ j × |ρ| × √(2/π) ≈ j × (2/p) × 0.80

Change in marginal probability: ΔF ≈ φ(z) × Δz where z ≈ 0 (threshold near median).
r_j ≈ 1 - ΔF/F(z) ≈ 1 - 2Δz/√(2π) ≈ 1 - 4j|ρ|/π

Total:
log₂(c₀_reps) ≈ Σ log₂(1 - 4j|ρ|/π) ≈ -(4|ρ|/π) × m(m+1)/(2 ln 2)

With |ρ| = 2/(p-3) and m = (p-1)/2:
≈ -(8/(π(p-3))) × ((p-1)/2)²/(2 ln 2) ≈ -p/(π ln 2) ≈ **-0.459p**

Predicted log₂(c₀_reps) by Gaussian heuristic:

| p  | Predicted | Actual | Ratio |
|----|-----------|--------|-------|
| 11 | -6.89     | -0.51  | 13.5  |
| 19 | -10.33    | -4.94  | 2.1   |
| 23 | -12.12    | -7.43  | 1.6   |
| 31 | -15.74    | -6.73  | 2.3   |

Gaussian overestimates by 2-14×. The ratio stabilizes around 2× for p ≥ 19.

### Critical Observation

Even the pessimistic Gaussian gives C ≈ 0.459 < 0.5:
log₂(E[N]) ≈ (0.5 - 0.459)p + O(√p) = 0.041p + O(√p) → ∞

The actual C ≈ 0.32 gives much more comfortable margin:
log₂(E[N]) ≈ (0.5 - 0.32)p = **0.18p → ∞**

## What Needs to be Proven

**Theorem (to prove)**: For p ≡ 3 mod 4 prime, sufficiently large, and D11 symmetric
with max_d |A(d) - (p+1)/4| ≤ O(√(p log p)):

log₂(c₀_reps) ≥ -p/(π ln 2) + O(√p)

equivalently, c₀_reps ≥ 2^{-p/(π ln 2) + O(√p)}

Since headroom_reps ≈ p/2 + 1.15√p > p/(π ln 2) ≈ 0.459p for large p,
this gives E[N] → ∞.

### Proof Strategy

1. **Normal approximation**: Show (B(d_1),...,B(d_{m-1})) is approximately
   multivariate normal with the equicorrelated covariance (on the (m-1)-simplex).
   - Berry-Esseen type: error O(1/√k) = O(1/√p) per marginal
   - Multivariate: use Stein's method on the Johnson scheme

2. **Conditional chain under normal**: For the (m-1)-variate Gaussian on the simplex
   {Σ x_i = const}, compute the exact conditional chain ratios.
   - Result: c₀_reps,Gauss = exp(-m²ρ²σ² × f(z)/2) where f depends on thresholds

3. **Transfer**: The error from discrete→Gaussian adds at most O(m/√p) = O(√p)
   to log₂(c₀_reps), which is absorbed in the O(√p) term.

### Gap in the Proof

The main difficulty is step 2: computing c₀ for the multivariate Gaussian on the
simplex. The covariance is SINGULAR (rank m-1 instead of m), and standard Slepian
doesn't apply directly.

However, projecting to the (m-1)-dimensional subspace orthogonal to 1 gives a
non-degenerate Gaussian with covariance:
Σ_proj = (1-ρ)I_{m-1} + ... (details needed)

In this subspace, the events E_d become half-spaces, and the joint probability
can be computed by the orthant probability formula.

For equicorrelated Gaussian in the non-degenerate (m-1)-dim projection:
Pr[all E_d] = Pr[X ∈ rectangular region] where X ~ N(0, Σ_proj)

This is a classical problem studied by Plackett (1954), Slepian (1962), etc.

## Alternative Approach: Direct Construction via Difference Sets

QR gives B(d) = (p-3)/4 exactly for all d. But the D11 threshold is (p-7)/4 = E[B] - 1.
So QR violates every D11 constraint by exactly 1.

If we could find D12 with B(d) = (p-3)/4 - 1 for d ∈ D11 and B(d) = (p-3)/4 + δ for d ∈ D22, the constraint would be satisfied.

This requires "spectrally complementary" D12 — exactly the pattern observed in SA solutions.

## Next Steps

1. Compute the (m-1)-dim projected Gaussian orthant probability exactly
2. Compare with actual c₀_reps values
3. Bound the discrete-to-Gaussian error using Stein's method
4. Alternatively: use the overlap decomposition to bound E[N²]/E[N]² directly
