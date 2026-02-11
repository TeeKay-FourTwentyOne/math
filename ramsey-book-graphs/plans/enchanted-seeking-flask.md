# Plan: QR Perturbation Approach for R(B_{n-1}, B_n) = 4n-1

**Reference**: See `docs/proof_strategy_analysis.md` for broader context and alternative strategies.

## Context

We want a general proof that R(B_{n-1}, B_n) = 4n-1 for all n. The upper bound is known (Rousseau & Sheehan 1978); we need the lower bound. For p = 2n-1 a prime power with p ≡ 1 mod 4, the Paley construction works. The remaining gap is **p ≡ 3 mod 4** (and composites).

The complete structural reduction (Theorem 3, `docs/book_ramsey_algebraic_reduction.md`) shows all constraints collapse to:

| For d in | Constraint |
|----------|-----------|
| D11      | A(d) + B(d) ≤ n-2 |
| D22      | A(d) + B(d) ≤ 2k-n+1 (k=\|D11\|) |

where A(d) = C(D11,D11;d), B(d) = C(D12,D12;d). Note: B(d) = B(-d) always (Universal Symmetry Lemma), so the four original constraints merge into two.

**Key fact**: For p ≡ 3 mod 4, D12 = QR(p) gives B(d) = (p-3)/4 for all nonzero d (constant). But this fails because A(d) must be below its average at every D11 position simultaneously — too demanding.

**This plan**: Perturb D12 = (QR ∪ {0}) \ {r} for r ∈ QR to introduce controlled variation in B(d) that compensates for A(d) peaks.

---

## Step 1: Derive and verify the perturbation formula

**Formula** (for d ≠ 0, d ≠ r, d ≠ p-r):

```
B_new(d) = (p-3)/4 - (chi(r+d) + chi(r-d)) / 2
```

where chi is the Legendre symbol (x/p). This gives B_new(d) ∈ {(p-7)/4, (p-3)/4, (p+1)/4}.

**Derivation sketch**: D12_new = (QR ∪ {0}) \ {r}. Relative to C(QR,QR;d):
- Gain from adding 0: pairs (0,d) and (-d,0) contribute +1_QR(d) + 1_QR(-d) = 1 (since exactly one of d, -d is QR for p ≡ 3 mod 4)
- Loss from removing r: pairs (r, r+d) and (r-d, r) contribute -1_QR(r+d) - 1_QR(r-d)
- Net: B_new(d) = (p-3)/4 + 1 - 1_QR(r+d) - 1_QR(r-d)
- Using 1_QR(x) = (1+chi(x))/2: = (p-3)/4 - (chi(r+d)+chi(r-d))/2

**Edge cases at d = r and d = p-r**: slight correction terms (compute directly).

**Verification**: Implement direct computation of B_new(d) for all d, compare with formula. This confirms correctness before using it.

**Files**: Create `ramsey-book-graphs/qr_perturbation_test.py`

---

## Step 2: Enumerate (D11, r) pairs for known SA solutions

For each prime p ≡ 3 mod 4 in {7, 11, 19, 23, 31, 43, 47, 59}:

1. Load the known SA-found D11 from `solutions_registry.json`
2. For each r ∈ QR(p) ((p-1)/2 choices):
   - Compute D12 = (QR ∪ {0}) \ {r}
   - Compute B_new(d) for all nonzero d (use direct computation, verified against formula)
   - Check all constraints: A(d) + B(d) ≤ n-2 for d ∈ D11; A(d) + B(d) ≤ 2k-n+1 for d ∈ D22
3. Record: which r values work, success fraction, which constraints are tightest

**Key metric**: success_rate(p) = |{r ∈ QR : (D11_known, r) is valid}| / |QR|

If success_rate is bounded away from 0 as p grows, an averaging argument over r could close the proof for p ≡ 3 mod 4.

---

## Step 3: Exhaustive D11 search (small primes)

For p = 23 (462 symmetric D11s) and p = 31 (6435 symmetric D11s):

1. Enumerate all symmetric D11 of size n (or n-1 if n is odd)
2. For each D11, try all r ∈ QR
3. Record: total (D11, r) pairs that work, joint success rate

This tells us whether the perturbation approach works for *some* D11 even if the SA-found D11 fails.

---

## Step 4: Analyze and extend

Based on results:
- **If single-swap works for most primes**: Analyze the success rate trend. If bounded away from 0, formalize the averaging argument.
- **If single-swap fails for large primes**: Extend to multi-swap: D12 = (QR ∪ {0} ∪ QNR_add) \ QR_remove with |QR_remove| = |QNR_add| + 1. Search over small swap counts k = 2, 3.
- **Character sum analysis**: For working (D11, r) pairs, analyze which chi(r+d) + chi(r-d) patterns produce the needed compensation. Connect to Weil bounds on character sums.

---

## Implementation details

**New file**: `ramsey-book-graphs/qr_perturbation_test.py`

Reuse from existing codebase:
- `quadratic_residues(p)` from `qr_d12_analysis.py`
- `autocorrelation(S, m)` from `qr_d12_analysis.py`
- Known D11 solutions from `solutions_registry.json`

Functions to implement:
1. `legendre_symbol(a, p)` — Euler criterion
2. `compute_b_perturbation(r, p)` — exact B_new(d) for all d via direct computation
3. `verify_formula(r, p)` — compare formula vs direct computation
4. `check_validity(D11, B, p, n)` — check all constraints
5. `sweep_all_r(D11, p, n)` — try all r ∈ QR, report success rate
6. `main()` — loop over primes, report results table

**Output format**:
```
p=7:  n=4,  |QR|=3,   success: 2/3 (66.7%)
p=11: n=6,  |QR|=5,   success: 3/5 (60.0%)
...
p=59: n=30, |QR|=29,  success: ?/29 (?%)
```

---

## Verification

1. Run `python qr_perturbation_test.py` — should produce the success rate table
2. For any working (D11, r) pair, verify by constructing the full 2-block circulant and checking all edge constraints directly (cross-check against the SA validator)
3. Compare B_new(d) from formula vs direct C(D12,D12;d) computation — should match exactly

---

## Success criteria

- **Strong signal**: success rate ≥ 10% for all primes p ≤ 59, not decreasing with p
- **Proof path**: if averaging argument works (expected success rate > 0 for random r), formalize as: "for any LP-feasible D11, a random r ∈ QR satisfies all constraints with probability ≥ c > 0"
- **Weak signal**: works only for specific D11 — still useful for case-by-case computation but not a general proof
