# V1V2 Identity and Constraint Simplification

## Theorem (V1V2 Identity)

For any prime p, any symmetric D11 ⊂ {1,...,p-1} with |D11| = (p+1)/2, and any D12 ⊂ {0,...,p-1} with 0 ∈ D12, |D12| = (p-1)/2:

$$X(d) := \Sigma(D_{11}, D_{12}, d, p) + \Delta(D_{12}, D_{22}, d, p) = |D_{12}| - \mathbf{1}_{d \in D_{12}}$$

for all d ∈ {0,...,p-1}, where D22 = {1,...,p-1} \ D11.

## Proof

Recall:
- Σ(D11, D12, d, p) = #{a ∈ D11 : (d-a) mod p ∈ D12}
- Δ(D12, D22, d, p) = #{b ∈ D12 : (b-d) mod p ∈ D22}

**Claim**: Each element i ∈ D12 with i ≠ d contributes exactly 1 to X(d).

For i ∈ D12, i ≠ d: the element (d-i) mod p ∈ {1,...,p-1} (since i ≠ d, and working mod p).

Since D11 ∪ D22 = {1,...,p-1} (disjoint), exactly one of two cases holds:

**Case A**: (d-i) mod p ∈ D11. Setting a = (d-i) mod p, we have a ∈ D11 and (d-a) mod p = i ∈ D12. This contributes +1 to Σ(D11, D12, d, p).

**Case B**: (d-i) mod p ∈ D22. Then (i-d) mod p = (-(d-i)) mod p ∈ D22 since D22 is symmetric (because D11 is symmetric and D22 = {1,...,p-1} \ D11). So b = i satisfies b ∈ D12 and (b-d) mod p ∈ D22. This contributes +1 to Δ(D12, D22, d, p).

**Element d itself**: If d ∈ D12, then (d-d) mod p = 0, and 0 ∉ D11 ∪ D22 (since D11 ⊂ {1,...,p-1} and D22 ⊂ {1,...,p-1}). So d is NOT counted in either sum.

**No double counting**: Each i ∈ D12 \ {d} contributes to exactly one of Σ or Δ (by the disjointness of cases A and B). Conversely, every term counted in Σ or Δ corresponds to some i ∈ D12 (the Σ term at a corresponds to i = (d-a) mod p ∈ D12, and the Δ term at b corresponds to i = b ∈ D12).

Therefore: X(d) = |D12 \ {d}| = |D12| - 1_{d ∈ D12}. □

## Consequence: V1V2 Constraints are Automatically Satisfied

With n = (p+1)/2, N = 2p, d₁ = |D11| + |D12| = p, d₂ = |D22| + |D12| = p-2:

**Red V1V2 edges** (d ∈ D12): Common neighbors = X(d) = |D12| - 1 = (p-1)/2 - 1 = (p-3)/2 = n-2.
Threshold: ≤ n-2 = (p-3)/2. **Satisfied with equality.** ✓

**Blue V1V2 edges** (d ∉ D12): Blue common neighbors = (N-2) - d₁ - d₂ + X(d) = (2p-2) - p - (p-2) + |D12| = |D12| = (p-1)/2 = n-1.
Threshold: ≤ n-1 = (p-1)/2. **Satisfied with equality.** ✓

This holds for ANY choice of (D11, D12) with the correct sizes. No optimization or constraint satisfaction is needed for V1V2.

## Consequence: Reduced Constraint System

After removing V1V2, the Ramsey constraints R(B_{n-1}, B_n) ≥ 4n-1 reduce to:

For each d ∈ {1,...,p-1}, define B(d) = Δ(D12, D12, d, p) = #{(a,b) ∈ D12² : a-b ≡ d mod p}.

**V1V1 constraints:**
- d ∈ D11 (red): A(d) + B(d) ≤ n-2 = (p-3)/2
- d ∈ D22 (blue): A(d) + B(d) ≤ n+1 = (p+3)/2

**V2V2 constraints** (using Δ(D12^T, D12^T, d, p) = B(p-d)):
- d ∈ D22 (red): C(d) + B(p-d) ≤ n-2 = (p-3)/2
- d ∈ D11 (blue): C(d) + B(p-d) ≤ n-3 = (p-5)/2

Combining for each d ∈ {1,...,p-1}, the threshold on B(d) is:
- **d ∈ D11**: T(d) = min(n-2-A(d), n-3-C(d))
- **d ∈ D22**: T(d) = min(n+1-A(d), n-2-C(d))

## Empirical Finding: V1V1 Dominance

**At p=11, 19, 23**: For ALL working D11, the V1V1 constraint is the binding one:
- d ∈ D11: T(d) = n-2-A(d) (V1V1 red binds, since A(d) ≥ C(d)+1)
- d ∈ D22: T(d) = n+1-A(d) (V1V1 blue binds, since A(d) ≥ C(d)+3)

This is explained by the size asymmetry: |D11|² - |D22|² = 2(p-1), so on average A(d) - C(d) = 2(p-2)/(p-1) ≈ 2.

**Implication**: The effective constraint system simplifies to:
$$\forall d \in \{1,...,p-1\}: A(d) + B(d) \leq \begin{cases} (p-3)/2 & d \in D_{11} \\ (p+3)/2 & d \in D_{22} \end{cases}$$

This is just the V1V1 codegree condition!

## Computational Verification

V1V2 identity verified at p = 7, 11, 19, 23, 31 (multiple random D11, D12 each).
Cross-validation against `ramsey_core.verify_construction`: 200/200 agreement at p = 11, 19, 23.
Binding type analysis: 100% V1V1 dominance for all 69 working D11 tested (p=11,19,23).

## Files
- Proof + computation: `ramsey-book-graphs/v1v2_identity_and_c0.py`
- Results: `ramsey-book-graphs/v1v2_c0_results.json`
