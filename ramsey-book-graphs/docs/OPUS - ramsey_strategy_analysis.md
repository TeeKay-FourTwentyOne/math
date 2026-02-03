# Critical Analysis: Ramsey Book Graph Implementation Strategy

**Date:** January 31, 2026  
**Subject:** Review of Gemini's Technical Specification for R(B_{n-1}, B_n) Search  
**Status:** Contains critical errors requiring immediate correction

---

## Executive Summary

The Gemini specification provides a reasonable foundation but contains **one critical mathematical error** and several strategic gaps that will waste computational resources if not addressed. This document details the issues and provides corrections based on the latest research (Wesley 2024, Lidický et al. 2024).

---

## 1. CRITICAL BUG: Blue Common Neighbor Formula

### The Error (Section 2.2)

The specification states:
> If edge is absent (BLUE), require λ_blue < n. Note: λ_blue = (N - 2) - λ_red.

**This formula is mathematically incorrect.**

### Correct Formula

For a non-edge (u,v) in the red graph, a vertex w is a common BLUE neighbor if and only if BOTH (u,w) AND (v,w) are non-edges in the red graph. This is NOT the same as "total vertices minus red common neighbors."

Let N_red(u) denote the red neighborhood of u. Then:

```
λ_blue(u,v) = |{w : w ∉ N_red(u) ∧ w ∉ N_red(v)}|
            = N - 2 - |N_red(u) ∪ N_red(v)|
            = N - 2 - deg_red(u) - deg_red(v) + λ_red(u,v)
```

For 2-block circulant graphs where degrees are constant within blocks:
- Let d₁ = |D₁₁| + |D₁₂| (degree of vertices in V₁)
- Let d₂ = |D₂₂| + |D₁₂| (degree of vertices in V₂)

Then for a blue edge (u,v):
- If u,v ∈ V₁: λ_blue = N - 2 - 2d₁ + λ_red
- If u,v ∈ V₂: λ_blue = N - 2 - 2d₂ + λ_red  
- If u ∈ V₁, v ∈ V₂: λ_blue = N - 2 - d₁ - d₂ + λ_red

**Impact:** Using the wrong formula will accept invalid constructions or reject valid ones.

---

## 2. KEY FINDING: D₂₂ = complement(D₁₁) in All Known Solutions

### Research Evidence

From Lidický et al. (2024):
> "The polycirculant graph constructions we found all have the extra property that the first circulant block (the subgraph induced by the first vertex-orbit of the 2-polycirculant automorphism) is isomorphic to the complement of the second circulant block."

From Wesley (2024):
> "In all of our cases we are showing bounds of the form R(B_r, B_s) ≥ 2m + 1, and we have D₂₂ = D̄₁₁ \ {0}, that is the complement of D₁₁ in Z_m \ {0}."

### Implication

**The search space can be dramatically reduced.** Instead of searching over independent (D₁₁, D₁₂, D₂₂), we should:

1. Fix D₂₂ = (Z_m \ {0}) \ D₁₁
2. Search only over (D₁₁, D₁₂)

This reduces the search space from ~3m/2 variables to ~m variables—a factor of ~2^(m/2) reduction in the SAT encoding.

### Recommended Implementation Change

```python
# Instead of independent D22
D22 = set(range(1, m)) - D11  # Complement of D11 in {1,...,m-1}
```

---

## 3. WRONG TARGET PRIORITIZATION

### The Error

The specification targets "composite n" specifically n = 18, 22, 27.

### Analysis

The actual mathematical criterion is whether **2n - 1 is a prime power ≡ 1 (mod 4)**, NOT whether n is composite.

| n | 2n-1 | Prime Power? | ≡ 1 (mod 4)? | Paley Works? | Status |
|---|------|--------------|--------------|--------------|--------|
| 18 | 35 = 5×7 | No | N/A | No | **Actually needs search** |
| 22 | 43 | Yes (prime) | 43 ≡ 3 (mod 4) | **No** | **Needs search** |
| 27 | 53 | Yes (prime) | 53 ≡ 1 (mod 4) | **Yes** | Should work with Paley! |

**n = 27 should NOT be a target** — the Paley construction should already work.

### Correct Priority Targets

The FIRST unknown case is n = 22, where verification was stopped ("generation simply started to take a lot of time" — Lidický et al.).

Correct priority order:
1. **n = 22** (first unknown, 2n-1 = 43 is prime ≡ 3 mod 4)
2. **n = 23** (2n-1 = 45 = 9×5, not prime power)
3. **n = 24** (2n-1 = 47 is prime ≡ 3 mod 4)
4. **n = 25** (2n-1 = 49 = 7², prime power but ≡ 1 mod 4 — should work!)
5. **n = 26** (2n-1 = 51 = 3×17, not prime power)

---

## 4. MISSING SYMMETRY BREAKING

### The Problem

The SAT encoding has massive symmetry that will cause the solver to explore equivalent solutions:

1. **Cyclic rotation:** If (D₁₁, D₁₂, D₂₂) is valid, so is (D₁₁ + k, D₁₂ + k, D₂₂ + k) for any k
2. **Reflection:** If D is valid, so is -D (mod m)
3. **Block swap:** If D₂₂ = D̄₁₁, swapping blocks gives the same graph

### Required Additions

Add lexicographic ordering constraints:
```
# Fix smallest element of D11 to be present (break rotation)
x_1 = True

# If d ∈ D11, require d ≤ m-d (break reflection)  
For all d > m/2: x_d → x_{m-d}
```

---

## 5. MISSING ALGEBRAIC SEEDING

### The Problem

The heuristic solver uses random initialization. This ignores all known structure.

### Better Approach

When 2n-1 = p is prime (but p ≡ 3 mod 4, so standard Paley fails):

1. **Try Paley anyway:** D₁₁ = quadratic residues mod p. Even if it doesn't satisfy all constraints, it may be close.

2. **Try twisted Paley:** D₁₁ = QR ∪ {some non-residues}, tuned to balance the common neighbor counts.

3. **Lift from nearby construction:** If a solution exists for n-1 or n+1, try adapting it.

### Specific Seeds for n = 22 (m = 43)

```python
# Quadratic residues mod 43
QR_43 = {1, 4, 6, 9, 10, 11, 13, 14, 15, 16, 17, 21, 23, 24, 25, 31, 35, 36, 38, 40, 41}
# |QR_43| = 21

# Try D11 = QR_43, D12 = QR_43, D22 = complement(D11)
# Then perturb from this starting point
```

---

## 6. PATTERN EXTRACTION SHOULD PRECEDE SEARCH

### Current Gap

The specification jumps directly to search without analyzing known solutions.

### Required Pre-Processing

Before any search, download and analyze ALL known constructions from:
- https://github.com/gwen-mckinley/ramsey-books-wheels
- https://github.com/Steven-VO/circulant-Ramsey/tree/master/RamseyGraphs/Books

Extract and tabulate:
1. |D₁₁|/m ratio for each n
2. |D₁₂|/m ratio for each n  
3. Verify D₂₂ = complement(D₁₁) holds universally
4. Check if D₁₂ ⊆ D₁₁ or D₁₂ ⊇ D₁₁ or neither
5. Look for algebraic structure in D₁₁ (cosets? quadratic residues? difference sets?)

---

## 7. FORMULA VERIFICATION

### The Wesley Lemma (Correct as Stated)

The common neighbor counting formulas in Section 2.1 appear correct, based on Wesley (2024) Lemma 7:

For u, v ∈ V₁ with difference d:
```
λ(u,v) = Δ(D₁₁, D₁₁, d) + Δ(D₁₂, D₁₂, d)
```

For u ∈ V₁, v ∈ V₂ with difference d:
```
λ(u,v) = Σ(D₁₁, D₁₂, d) + Δ(D₁₂, D₂₂, d)
```

Where:
- Δ(A, B, d) = |{(a,b) ∈ A×B : a - b ≡ d (mod m)}|
- Σ(A, B, d) = |{(a,b) ∈ A×B : a + b ≡ d (mod m)}|

### Implementation Note

These can be computed efficiently in O(|A|) time using indicator arrays:
```python
def delta(A, B, d, m):
    count = 0
    B_indicator = [False] * m
    for b in B:
        B_indicator[b] = True
    for a in A:
        if B_indicator[(a - d) % m]:
            count += 1
    return count
```

---

## 8. RECOMMENDED REVISED EXECUTION PLAN

### Phase 0: Data Acquisition (NEW)
1. Clone both GitHub repositories
2. Parse all known constructions for n ≤ 21
3. Tabulate structural properties
4. Verify D₂₂ = complement(D₁₁) hypothesis

### Phase 1: Core Implementation (REVISED)
1. Implement ramsey_core.py with **correct** blue neighbor formula
2. **Enforce D₂₂ = complement(D₁₁)** to reduce search space
3. Add symmetry breaking constraints
4. Verify against known solutions for n = 6, 10, 14, 21

### Phase 2: Targeted Search (REVISED)
1. **Start with n = 22**, not n = 18 (which is already solved)
2. Seed heuristic search with Paley-like structures
3. Use SAT with reduced variables (only D₁₁, D₁₂)
4. Implement incremental solving: if n works, warm-start n+1

### Phase 3: Analysis & Generalization
1. If solutions found, analyze for patterns
2. Attempt to formulate general algebraic construction
3. Test hypothesized construction on larger n

---

## 9. RESOURCE ESTIMATES

### SAT Encoding Size (Corrected)

For n = 22 (m = 43):
- Variables for D₁₁: ⌊43/2⌋ = 21 (exploiting symmetry)
- Variables for D₁₂: 43
- Variables for D₂₂: **0** (if we enforce D₂₂ = complement)
- **Total: 64 variables** (down from ~107)

Constraints per difference value d:
- Red constraint: atmost(n-2) on ~m² terms → O(m²) clauses
- Blue constraint: atmost(n-1) on ~m² terms → O(m²) clauses
- Total: O(m³) clauses ≈ 80,000 clauses for n = 22

This should be tractable for modern SAT solvers.

---

## 10. OPEN QUESTIONS FOR INVESTIGATION

1. **Why does D₂₂ = complement(D₁₁) work?** Is there a theoretical reason, or is it just empirically observed?

2. **What determines D₁₂?** In some constructions D₁₂ = D₁₁. Is this always optimal?

3. **Can we use number-theoretic sieves?** If 2n-1 has special factorization properties, can we exploit them?

4. **Is there a spectral characterization?** Paley graphs have nice eigenvalue properties. Can we formulate the book-free condition spectrally?

---

## Summary of Required Changes

| Issue | Severity | Fix |
|-------|----------|-----|
| Blue neighbor formula | **CRITICAL** | Use λ_blue = N - 2 - deg(u) - deg(v) + λ_red |
| Wrong targets | HIGH | Focus on n=22,23,24,26 not n=18,27 |
| Missing D₂₂ constraint | HIGH | Enforce D₂₂ = complement(D₁₁) |
| No symmetry breaking | MEDIUM | Add lexicographic constraints |
| Random initialization | MEDIUM | Seed with Paley-like structures |
| No pattern extraction | MEDIUM | Analyze known solutions first |

---

## References

1. Wesley, W.J. (2024). "Lower Bounds for Book Ramsey Numbers." arXiv:2410.03625v2.

2. Lidický, B., McKinley, G., Pfender, F., Van Overberghe, S. (2024). "Small Ramsey numbers for books, wheels, and generalizations." arXiv:2407.07285v2.

3. Rousseau, C.C., Sheehan, J. (1978). "On Ramsey numbers for books." J. Graph Theory 2:77-87.

4. Radziszowski, S.P. (2024). "Small Ramsey Numbers." Electronic J. Combinatorics, Dynamic Survey DS1.
