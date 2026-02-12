# p_working Scaling Analysis: Fraction of Working Orbits

## Statement
The proof of R(B_{n-1}, B_n) = 4n-1 for all n reduces to showing that the fraction p_working of symmetric D11 (in multiplicative orbits) with N(D11) > 0 satisfies p_working ≥ C/p for some constant C > 0 and all sufficiently large primes p ≡ 3 mod 4.

## Definitive Table

| p | #orbits | #working | p_working | p × p_working | Source |
|---|---------|----------|-----------|---------------|--------|
| 11 | 2 | 1 | 0.500 | **5.50** | Exact |
| 19 | 14 | 1 | 0.071 | **1.36** | Exact |
| 23 | 42 | 5 | 0.119 | **2.74** | Exact |
| 31 | 429 | 18 | 0.042 | **1.30** | Exact |
| 43 | 16,796 | ≥8 | ~0.010 | **~0.43** | SA (3/300 random) |
| 47 | 58,786 | ≥1 | ? | ? | Known sol exists |
| 59 | 2,674,440 | ≥1 | ? | ? | Known sol exists |

**Orbit count formula**: #orbits = C((p-1)/2, (p+1)/4) / ((p-1)/2), grows as ~2^{p/2}/sqrt(p).

## Key Corrections
- **p=31: 18 working orbits** (not 7 as in prior findings). The 7 was just the count with minimum N=62. Full exact exhaustive C enumeration found 18/429 = 0.042 working.
- **p=43 prior "6/268"**: The 268 was a targeted sample (e.g., A-flat D11), not random. Random SA scan of 300 orbits from the full 16,796 found 3, of which 1 matched the known 6 and 2 were NEW working orbits.
- **p=47, p=59**: SA solutions in the registry used |D11| = (p-3)/2 (not the standard (p+1)/2). These convert to standard form by swapping V1↔V2: new D11 = old D22, new D12 = -D12 mod p. Verified valid.

## SA Reliability
- C SA solver tested on all 6+1 known working D11 at p=43 and p=47: **100% success rate** (0-2 trials, 167K-794K iterations).
- C SA orbit scan: 300 random orbits per prime, 12 trials × 2M iterations each.
- At p=23, SA found 19/200 working (exact: 5/42 = 0.119, SA: 0.095). Mild underestimate (~20% miss rate across orbits due to random sampling, not SA failure).

## 95% Wilson CI for p=43
- p_working ∈ [0.0034, 0.029]
- p × p_working ∈ [0.15, 1.25]
- Estimated #working orbits: 57 – 486 (out of 16,796)

## Assessment

### Exact data (p ≤ 31): Supports p_working = Θ(1/p)
- p × p_working = {5.50, 1.36, 2.74, 1.30} — oscillates but stays in [1, 6].
- The oscillation (1.36 → 2.74 → 1.30) suggests non-monotonic behavior.

### SA data (p = 43): Inconclusive
- Point estimate p × p_working = 0.43, below the [1, 6] range of exact primes.
- But 95% CI includes 1.25, so consistent with stabilization near 1.
- Would need ~3000 orbits (not 300) to determine if p × p_working < 1 with confidence.

### p = 47, 59: Existence confirmed, fraction unknown
- Working orbits exist (at least 1 each), confirming R(B_{n-1}, B_n) = 4n-1 at these primes.
- Random orbit sampling is infeasible (58K – 2.7M total orbits; too few working to hit randomly).
- Targeted search (using spectral/algebraic criteria to filter candidate D11) would be needed.

### Bottom line
The conjecture p_working = Θ(1/p) is **supported by exact data for p ≤ 31** but **unresolved for p ≥ 43**. The SA estimate at p=43 is marginally consistent with the conjecture but could also indicate a steeper decline. Proving the conjecture likely requires structural arguments about the distribution of working D11 among multiplicative orbits, not computational verification.

## Files
- SA orbit scanner: `ramsey-book-graphs/sa_orbit_scan.c`
- SA reliability test: `ramsey-book-graphs/sa_orbit_test.c`
- Exact p=31 data: `ramsey-book-graphs/exact_N_p31_results.jsonl`
- Algebraic invariants: `ramsey-book-graphs/algebraic_invariants{,_large}.json`
