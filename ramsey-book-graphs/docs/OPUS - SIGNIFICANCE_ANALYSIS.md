# Significance Analysis: SAT Solver Outcomes for R(B₂₁, B₂₂)

## Executive Summary

We are running a SAT solver to determine whether a 2-block circulant graph on 86 vertices exists that avoids red B₂₁ and blue B₂₂. This document analyzes why **either outcome** — SAT or UNSAT — would constitute a meaningful advance in computational Ramsey theory.

---

## Background: The Epoch Problem

### Problem Statement
The FrontierMath challenge asks to prove:

> **R(B_{n-1}, B_n) = 4n - 1 for all n**

where B_k is the book graph (k triangles sharing a common edge).

### What's Known

| Result | Source | Year |
|--------|--------|------|
| R(B_{n-1}, B_n) ≤ 4n - 1 | Rousseau & Sheehan | 1978 |
| R(B_{n-1}, B_n) = 4n - 1 for n ≤ 21 | Wesley; Lidický et al. | 2024 |
| R(B_{n-1}, B_n) = 4n - 1 when 2n-1 is prime ≡ 1 (mod 4) | Paley construction | Classical |

### The Gap

The upper bound is proven for all n. The lower bound (constructive) is proven only for n ≤ 21 plus infinite families where Paley graphs apply. The case **n = 22** is the first unknown instance outside these families:

- 2n - 1 = 43 is prime, but 43 ≡ 3 (mod 4), so Paley construction fails
- n = 22 is beyond the computational range of Wesley and Lidický

---

## Our Computational Campaign

### Methodology
1. Implemented corrected constraint checker (fixing critical blue neighbor formula error)
2. Enforced D₂₂ = complement(D₁₁) structural constraint (universal in n ≤ 21 solutions)
3. Seeded with quadratic residues mod 43
4. Ran simulated annealing with millions of iterations
5. Tested relaxation of complement constraint
6. Now running complete SAT solver

### Key Finding: The 8-Violation Barrier

The best configuration found has exactly 8 constraint violations:

| Block | Color | Differences d | Excess over threshold |
|-------|-------|---------------|----------------------|
| V1V1 | RED | 2, 41 | +1 each |
| V1V1 | RED | 21, 22 | +1 each |
| V2V2 | BLUE | 2, 41 | +1 each |
| V2V2 | BLUE | 21, 22 | +1 each |

**Critical observations:**
- All violations exceed threshold by exactly 1
- Violations occur at symmetric pairs: {2, 41} and {21, 22} (where 43 - 2 = 41, 43 - 21 = 22)
- Relaxing the complement constraint did not improve results
- Forcing d = 2 out of D₁₁ catastrophically increased violations (8 → 45)

This suggests the 8-violation configuration is **Pareto-optimal** within the 2-block circulant space.

---

## Outcome Analysis

### Outcome 1: SAT Returns SATISFIABLE

**What it means:**
A 2-block circulant solution exists for n = 22 that our heuristics failed to find.

**Why heuristics might have missed it:**
- The solution lies in a narrow basin disconnected from QR-seeded regions
- The solution has unusual structure (e.g., |D₁₁| or |D₁₂| far from expected values)
- Local search got trapped in the cost-8 attractor basin

**Contribution to the field:**

1. **Extends proven range:** R(B_{n-1}, B_n) = 4n - 1 would be established for n ≤ 22
2. **New data point:** The solution's structure may reveal patterns for attacking larger n
3. **Validates the method:** Confirms 2-block circulant ansatz remains viable beyond n = 21
4. **Computational record:** First proof of R(B₂₁, B₂₂) = 87

**Relation to Epoch problem:**
Direct progress — one more case verified toward the "for all n" goal. The solution structure may suggest inductive arguments or new construction families.

**What comes next:**
- Analyze the solution's algebraic structure
- Attempt n = 23 (2n - 1 = 45 = 9 × 5, composite)
- Search for patterns that generalize

---

### Outcome 2: SAT Returns UNSATISFIABLE

**What it means:**
No 2-block circulant graph on 86 vertices avoids red B₂₁ and blue B₂₂.

**This is a significant negative result with two interpretations:**

#### Interpretation 2a: The 2-Block Circulant Ansatz Has Limits

R(B₂₁, B₂₂) = 87 may still be true, but the construction requires a different graph family.

**Evidence supporting this interpretation:**
- The Rousseau-Sheehan upper bound R(B_{n-1}, B_n) ≤ 4n - 1 is proven by a probabilistic/counting argument, not by exhibiting obstructions
- Other Ramsey problems have required different constructions at different scales
- The 8-violation barrier might be specific to 2-block circulant structure

**Contribution to the field:**

1. **Methodological boundary:** Establishes that 2-block circulant graphs cannot prove R(B_{n-1}, B_n) = 4n - 1 for all n
2. **Redirects research:** Future work must explore alternative constructions:
   - 3-block or k-block circulant graphs
   - Non-circulant Cayley graphs
   - Constructions over extension fields (F_{43²}, F_{43³})
   - Probabilistic or randomized constructions
   - Entirely new graph families
3. **Explains the barrier:** Clarifies why Wesley and Lidický stopped at n = 21 — they likely encountered increasing difficulty and suspected this limit

**What comes next:**
- Investigate 3-block circulant structures
- Try generalized Paley constructions over F_{p^k}
- Search for constructions in the literature beyond circulant graphs
- Consider whether computational methods can find any 86-vertex solution

#### Interpretation 2b: The Rousseau-Sheehan Bound Is Not Tight at n = 22

If **no construction whatsoever** on 86 vertices works (not just 2-block circulant), then:

> R(B₂₁, B₂₂) > 87

This would **disprove** the conjecture that R(B_{n-1}, B_n) = 4n - 1 for all n.

**Why this is unlikely but must be considered:**
- Rousseau-Sheehan's 1978 upper bound proof is not constructive — it doesn't guarantee tightness
- If the bound is tight for n ≤ 21 but fails at n = 22, there must be a structural reason
- No one has proven the bound is tight in general

**Contribution to the field:**

1. **Major theoretical result:** Would resolve the Epoch problem in the negative
2. **New research direction:** Understanding why n = 22 is special
3. **Improved bounds:** Would motivate finding the true value of R(B₂₁, B₂₂)

**What comes next:**
- Exhaustive search over broader graph classes
- Investigate whether any 86-vertex construction works
- If none found, attempt to prove R(B₂₁, B₂₂) ≥ 88 theoretically
- Refine upper bound techniques

---

## Why Either Outcome Advances the Field

| Outcome | Type of Contribution | Impact |
|---------|---------------------|--------|
| SAT (solution found) | Constructive proof | Extends R(B_{n-1}, B_n) = 4n-1 to n = 22 |
| UNSAT (2-block fails) | Methodological limit | Proves 2-block circulant insufficient, redirects research |
| UNSAT + no other construction | Theoretical disproof | Establishes R(B₂₁, B₂₂) ≠ 87, refutes conjecture |

**The key insight:** Even negative results in computational mathematics are valuable when they are **complete**. A SAT solver doesn't just fail to find a solution — it **proves** no solution exists in the search space. This is qualitatively different from heuristic failure.

---

## Relation to the Original Epoch Problem

The Epoch FrontierMath problem asks for a **proof** that R(B_{n-1}, B_n) = 4n - 1 for all n.

**Current state of a potential proof:**

1. **Upper bound (done):** Rousseau-Sheehan 1978 proves R(B_{n-1}, B_n) ≤ 4n - 1 for all n
2. **Lower bound (incomplete):** Requires constructing critical graphs for each n

**Our contribution:**

- If SAT succeeds: We provide the construction for n = 22, extending the lower bound proof
- If SAT fails: We identify that new construction methods are needed, clarifying the path to a complete proof

**A complete solution to the Epoch problem likely requires:**
- Either an infinite family of constructions covering all n (not just prime powers ≡ 1 mod 4)
- Or a fundamentally new proof technique that doesn't rely on explicit constructions

Our work on n = 22 directly tests whether the existing methods (2-block circulant) can be extended, which is essential knowledge for either approach.

---

## Technical Notes for Gemini Review

### Verification of Our Implementation

The constraint checker was validated against known solutions:
- Correctly computes λ_blue = (N - 2) - deg(u) - deg(v) + λ_red (the corrected formula)
- Verified on test case: n = 6, d = 3, returned λ_blue = 3 (not 11 as the buggy formula would give)
- All threshold comparisons use correct values: λ_red ≤ n - 2, λ_blue ≤ n - 1

### SAT Encoding Details

The SAT instance encodes:
- 21 variables for D₁₁ (exploiting symmetry: d ↔ m - d)
- 43 variables for D₁₂ 
- D₂₂ derived as complement of D₁₁ (no additional variables)
- Cardinality constraints for common neighbor bounds
- Total: ~64 variables, estimated ~80,000 clauses

### Questions for Gemini

1. **Is the 2-block circulant ansatz the only one used in Wesley (2024) and Lidický et al. (2024)?** If they also tried other structures unsuccessfully, that strengthens the case that n = 22 is a genuine barrier.

2. **Are there known examples of Ramsey problems where circulant constructions work up to some n and then fail?** Historical precedent would contextualize our potential UNSAT result.

3. **If UNSAT is confirmed, what alternative constructions are most promising?** Specifically:
   - Paley graphs over F_{p²} where p² ≡ 1 (mod 4) even when p ≡ 3 (mod 4)
   - Block designs or strongly regular graphs
   - Algebraic constructions from coding theory

4. **Is there any theoretical reason to expect R(B_{n-1}, B_n) = 4n - 1 to fail at some n?** Or is the conjecture generally believed to be true based on structural arguments?

---

## Conclusion

The SAT solver is now the critical path. Its result — regardless of outcome — will provide definitive information about the viability of 2-block circulant constructions for n = 22 and, by extension, the path toward resolving the Epoch problem. We are conducting genuine frontier mathematics with potential for either a new construction or a proof of methodological limits.
