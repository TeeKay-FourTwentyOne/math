# Plan: Prove R(B_{n-1}, B_n) = 4n-1 for Primes p = 3 mod 4

## Context

We have proven R(B_{n-1}, B_n) = 4n-1 for all n where m = 2n-1 is a prime power ≡ 1 mod 4 (Paley construction). The verified range extends continuously through n ≤ 31. The next step toward proving the conjecture for ALL n is to handle primes p ≡ 3 mod 4, which would give a second infinite family. Combined with Paley, this covers all n where 2n-1 is prime — a positive-density set.

**Why this path over extending verified range:** Each SA verification is independent and doesn't generalize. The key lemma, if proven, handles infinitely many cases at once. The empirical data strongly suggests it's true: the normalized difficulty *decreases* with n.

**Literature:** Wesley (2025) independently proved the same Paley family. The p ≡ 3 mod 4 case is **open in the literature**. Our V1V2 theorem and V2V2 reduction go beyond published results.

## The Key Lemma

For prime p ≡ 3 mod 4, there exist:
- Symmetric D11 ⊂ {1,...,p-1} with |D11| = (p+1)/2
- D12 ⊂ Z_p with |D12| = (p-1)/2, 0 ∈ D12

Such that all four constraints hold (where A(d) = Delta(D11,D11,d), B(d) = Delta(D12,D12,d)):

| Constraint | Condition | Threshold | Binding? |
|-----------|-----------|-----------|----------|
| V1V1 red | d ∈ D11 | A(d)+B(d) ≤ (p-3)/2 | YES |
| V1V1 blue | d ∉ D11 | A(d)+B(d) ≤ (p+1)/2 | loose |
| V2V2 blue | d ∈ D11 | A(d)+B(p-d) ≤ (p-3)/2 | YES |
| V2V2 red | d ∉ D11 | A(d)+B(p-d) ≤ (p+3)/2 | loose |

The binding constraints require A(d)+B(d) and A(d)+B(p-d) both ≤ (p-3)/2 for d ∈ D11. The average A(d)+B(d) ≈ (p-1)/2, so this asks for values ~1 below average — but with std ~0.3√p growing, this becomes proportionally easier for larger p.

## The Core Difficulty

A naive probabilistic argument (random D12) fails because:
- E[B(d)] ≈ (p-3)/4, and A(d) ≈ (p+1)/4 on average
- So E[A(d)+B(d)] ≈ (p-1)/2 = threshold + 1
- Each constraint violates with probability ≈ 1/2 - O(1/√p)
- Expected violations ≈ p/4 → union bound diverges

The proof must exploit **structure**: either the negative correlation between D11 membership and A(d), or an algebraic/spectral property of a carefully chosen (D11, D12) pair.

---

## Implementation Plan

### Phase 1: Exhaustive Enumeration for Small Primes

**File:** `ramsey-book-graphs/enumerate_solutions.py`

For p = 7, 11, 19 (n = 4, 6, 10), the search space is small enough for complete enumeration of ALL valid (D11, D12) pairs.

**Algorithm:**
- p=7: 3 symmetric D11 choices × 10 D12 choices = 30 total → trivial
- p=11: 10 D11 × 210 D12 = 2,100 → trivial
- p=19: 126 D11 × 43,758 D12 ≈ 5.5M → feasible with fast Delta evaluation

**Data to collect per valid pair:**
- Power spectrum P(k) = |D̂11(k)|² + |D̂12(k)|²
- How close P(k) is to the "Legendre pair" value p (flat spectrum)
- Phase relationship between D̂11(k) and D̂12(k)
- Cyclotomic class distribution of D11 and D12
- |D12 ∩ D12^T| (asymmetry measure)

**Key question answered:** For each D11, how many valid D12 exist? Do ALL symmetric D11 admit a valid D12, or only special ones? This determines whether we prove the universal or existential version of the lemma.

### Phase 2: Statistical Sampling for Medium Primes

**File:** `ramsey-book-graphs/sample_d12_space.py`

For p = 23, 31, 43, 47, 59 (where we already have SA solutions):

1. **Fix the known D11** → sample 10⁶ random D12 of correct size containing 0
   - Measure: fraction with zero violations, violation cost distribution
   - Answer: is valid D12 "rare" (1 in 10⁶) or "common" (1 in 100)?

2. **Random D11 survey** → sample 1000 random symmetric D11, for each sample 10⁴ D12
   - Answer: does D11 choice matter? Are some D11 much better than others?

3. **Spectral correction analysis** → for the known valid pairs, compute the "residual spectrum" R(k) = P(k) - |D̂11(k)|² that D12 provides
   - Answer: does D12's spectrum "correct" D11's spectrum toward flatness?

### Phase 3: Correlation Structure Analysis

**File:** `ramsey-book-graphs/correlation_analysis.py`

The make-or-break question for the probabilistic proof: are the bad events (B_d for d ∈ D11) positively or negatively correlated?

For each prime p up to ~100:
1. Fix a known D11
2. Compute the exact covariance matrix Cov(B(d), B(d')) for d, d' ∈ D11 under uniformly random D12
3. Compute the joint probability Pr[all constraints satisfied] using:
   - Multivariate normal approximation with the exact covariance
   - Direct Monte Carlo (for small p)
4. Compare to what union bound / LLL would predict

**If negative correlation:** The bad events tend NOT to co-occur → proof via FKG inequality or second moment method should work.

**If positive correlation:** Need to understand the dependency graph for LLL, or abandon the probabilistic approach.

### Phase 4: Proof Attempt

Based on findings from Phases 1-3, pursue the most promising approach:

**Approach A (probabilistic, if correlations are favorable):**

File: `ramsey-book-graphs/prob_existence.py`

1. Choose D12 via a **structured random construction** (not uniform):
   - Define D12 as a "spectral corrector" — choose D12 so |D̂12(k)|² ≈ T(k) where T(k) is computed from D11
   - This biases D12 toward satisfying the constraints while remaining random enough for concentration arguments
2. Prove tail bound: Pr[B(d) > threshold - A(d)] ≤ exp(-c·f(p)) for each d ∈ D11
3. Union bound over |D11| = (p+1)/2 constraints: works if c·f(p) > log(p)
4. This gives "for all p ≥ p₀" existence; SA covers p < p₀

**Approach B (algebraic, if Phase 1 reveals structure):**

File: `ramsey-book-graphs/algebraic_p3mod4.py`

1. For p ≡ 3 mod 4, construct D11 using cyclotomic classes:
   - D11 consists of (p+1)/4 negation pairs {d, p-d} from {1,...,p-1}
   - Choose pairs based on index (discrete log) mod small numbers
2. Derive D12 algebraically from D11 (recipe from Phase 1 enumeration patterns)
3. Verify computationally for all p ≡ 3 mod 4 up to p ~ 500
4. Prove the construction works using character sum estimates (Weil bound, etc.)

**Approach C (Legendre pair connection):**

File: `ramsey-book-graphs/legendre_relaxation.py`

1. Our constraint is a relaxation of the Legendre pair condition: we need |D̂11(k)|² + |D̂12(k)|² to produce bounded fluctuations, not exact constancy
2. Start with "approximate Legendre pairs" from the literature (Djokovic-Kotsireas)
3. Show the approximation error is within our tolerance
4. This connects to supplementary difference sets — existing constructions may apply

### Phase 5: Verification & Assembly

**SA extension:** Push verified range to cover all primes p ≡ 3 mod 4 below whatever p₀ the proof requires (run `fast_sa_general.py` for p = 67, 71, 79, 83).

**Independent verification:** Every construction (SA or algebraic) verified by `validate_construction.py`.

**Proof document:** Update `docs/proof_outline.md` with the new theorem covering p ≡ 3 mod 4.

---

## Execution Order

| Step | Phase | Description | Dependencies |
|------|-------|-------------|-------------|
| 1 | 1 | Exhaustive enumeration (p=7,11,19) | None |
| 2 | 2 | D12 sampling for known solutions (p=23,31,43,47,59) | None |
| 3 | 2 | Random D11 survey | None |
| 4 | 3 | Correlation structure analysis | Steps 1-2 inform interpretation |
| 5 | 4 | Proof attempt (approach chosen based on 1-4) | Steps 1-4 |
| 6 | 5 | SA extension + verification | Step 5 determines p₀ |

Steps 1-3 can run **in parallel**. Step 4 can start as soon as 1-2 produce data. The proof approach (step 5) is chosen based on empirical results.

## Critical Files

| File | Role |
|------|------|
| `ramsey_core.py` | Delta, Sigma, verify_construction — foundation |
| `fourier_constraints.py` | Fourier analysis, character sums — template for spectral code |
| `fast_sa_general.py` | SA solver — extend for new primes |
| `validate_construction.py` | Independent brute-force validator |
| `docs/proof_outline.md` | Proof document to update |

## Success Criteria

- **Minimum:** Identify the correct proof approach with strong computational evidence
- **Target:** Prove the key lemma for all p ≡ 3 mod 4 with p ≥ p₀, verify SA constructions for p < p₀
- **Stretch:** Find an algebraic construction that works for ALL p ≡ 3 mod 4 (no SA needed)
