# Proof Outline: R(B_{n-1}, B_n) = 4n - 1

**Date**: 2026-02-10
**Status**: Verified for all n <= 31. Proven for infinitely many n (Paley family).

---

## 1. Problem Statement

**Book graph** B_n = K_2 + K̄_n: n triangles sharing a common edge (2n+2 vertices).

**Goal**: Prove R(B_{n-1}, B_n) = 4n - 1 for all n >= 3.

**Known**: Upper bound R(B_{n-1}, B_n) <= 4n - 1 (Rousseau & Sheehan, 1978).

**Our contribution**: Lower bound constructions proving R(B_{n-1}, B_n) >= 4n - 1.

---

## 2. Construction Framework

All constructions use **2-block circulant graphs** on N = 4n - 2 = 2m vertices (m = 2n - 1):

- Vertex set V = V_1 ∪ V_2 where |V_1| = |V_2| = m
- Adjacency defined by difference sets D11, D12, D22 ⊂ Z_m (or GF(q) for Paley)
- Red edge (u,v) iff their difference is in the corresponding set; blue otherwise

**Universal structural constraints** (verified for all known constructions):

| Property | Constraint |
|----------|-----------|
| D22 = {1,...,m-1} \ D11 | D22 is the complement of D11 |
| D11 symmetric | x ∈ D11 ⟺ -x ∈ D11 |
| 0 ∈ D12 | Always (for m ≡ 3 mod 4 cases) |
| \|D12\| = n - 1 | Cross-block set has exactly n-1 elements |
| \|D11\| even | Required since D11 is symmetric in Z_m (m odd) |

**Validity conditions**: The graph avoids red B_{n-1} and blue B_n iff:
- Every red edge has ≤ n-2 red common neighbors
- Every blue edge has ≤ n-1 blue common neighbors

---

## 3. Key Theorems

### Theorem 1: V1V2 Auto-Satisfaction (NEW)

**If** D22 = complement(D11), D11 is symmetric, and |D12| = n-1, **then** for all d ∈ {0,...,m-1}:

```
V1V2 common neighbors at d = |D12| - [d ∈ D12]
```

- Red V1V2 edges: exactly n-2 common neighbors (at threshold)
- Blue V1V2 edges: exactly n-1 common neighbors (at threshold)

**Proof**: Sigma(D11, D12, d) = Delta(D12, D11, d) by symmetry of D11 and commutativity of convolution. Terms cancel exactly. See `docs/pattern_analysis.md` for full proof.

**Consequence**: V1V2 constraints are automatically satisfied. The problem reduces to V1V1 and V2V2 only.

**Computational verification**: Confirmed for all 9 known constructions (n=6,8,10,12,14,16,18,20,22) via independent brute-force (`verify_v1v2_theorem.py`).

### Theorem 2: V2V2 Algebraic Determination (NEW)

With A(d) = Delta(D11, D11, d) and B(d) = Delta(D12, D12, d):

```
V2V2(d) = A(d) + B(-d) + (m - 2 - 2|D11|) + 2·[d ∈ D11]
```

The V1V1 red edges correspond exactly to the V2V2 blue edges and vice versa. V2V2 is fully determined by the V1V1 autocorrelation values and the D12 autocorrelation at negated arguments.

### Theorem 3: Paley Construction (PROVEN)

**For any prime power q ≡ 1 (mod 4)**, the construction over the Cayley graph on (GF(q)+, GF(q)+):

```
D11 = D12 = QR(GF(q)*),  D22 = QNR(GF(q)*)
```

is a valid construction proving R(B_{n-1}, B_n) = 4n - 1 where n = (q+1)/2.

**Proof**: QR in GF(q) is an almost-difference set with Delta(QR, QR, d) = (q-5)/4 for d ∈ QR and (q-1)/4 for d ∈ QNR. This gives V1V1(d) = 2·Delta(QR,QR,d) which satisfies:
- Red (d ∈ QR): (q-5)/2 ≤ (q-3)/2 = n-2 ✓
- Blue (d ∈ QNR): blue_common = (q-1)/2 = n-1 ✓

V1V2 is free by Theorem 1. V2V2 follows by symmetry since -1 ∈ QR when q ≡ 1 (mod 4).

**Computational verification**: All 15 prime powers q ≡ 1 (mod 4) up to q = 100 independently verified via brute-force (`validate_paley_bruteforce.py`), including the non-trivial case GF(3^4) = GF(81).

---

## 4. Coverage Map for n ≤ 50

| n | m = 2n-1 | Factorization | Type | Status | Method |
|---|----------|---------------|------|--------|--------|
| 3 | 5 | prime | pp ≡ 1(4) | **PROVEN** | Paley |
| 4 | 7 | prime | p ≡ 3(4) | verified | prior work |
| 5 | 9 | 3² | pp ≡ 1(4) | **PROVEN** | Paley |
| 6 | 11 | prime | p ≡ 3(4) | verified | prior work |
| 7 | 13 | prime | pp ≡ 1(4) | **PROVEN** | Paley |
| 8 | 15 | 3×5 | composite | verified | prior work |
| 9 | 17 | prime | pp ≡ 1(4) | **PROVEN** | Paley |
| 10 | 19 | prime | p ≡ 3(4) | verified | prior work |
| 11 | 21 | 3×7 | composite | verified | prior work |
| 12 | 23 | prime | p ≡ 3(4) | verified | prior work |
| 13 | 25 | 5² | pp ≡ 1(4) | **PROVEN** | Paley |
| 14 | 27 | 3³ | p ≡ 3(4) | verified | prior work |
| 15 | 29 | prime | pp ≡ 1(4) | **PROVEN** | Paley |
| 16 | 31 | prime | p ≡ 3(4) | verified | prior work |
| 17 | 33 | 3×11 | composite | verified | prior work |
| 18 | 35 | 5×7 | composite | verified | prior work |
| 19 | 37 | prime | pp ≡ 1(4) | **PROVEN** | Paley |
| 20 | 39 | 3×13 | composite | verified | prior work |
| 21 | 41 | prime | pp ≡ 1(4) | **PROVEN** | Paley |
| 22 | 43 | prime | p ≡ 3(4) | verified | SA-informed |
| 23 | 45 | 3²×5 | composite | **verified** | fast SA (NEW) |
| 24 | 47 | prime | p ≡ 3(4) | **verified** | fast SA (NEW) |
| 25 | 49 | 7² | pp ≡ 1(4) | **PROVEN** | Paley (NEW) |
| 26 | 51 | 3×17 | composite | **verified** | fast SA (NEW) |
| 27 | 53 | prime | pp ≡ 1(4) | **PROVEN** | Paley |
| 28 | 55 | 5×11 | composite | **verified** | fast SA (NEW) |
| 29 | 57 | 3×19 | composite | **verified** | fast SA (NEW) |
| 30 | 59 | prime | p ≡ 3(4) | **verified** | fast SA (NEW) |
| 31 | 61 | prime | pp ≡ 1(4) | **PROVEN** | Paley |
| 32 | 63 | 3²×7 | composite | OPEN | n=32 hard |
| 33 | 65 | 5×13 | composite | OPEN | |
| 34 | 67 | prime | p ≡ 3(4) | OPEN | SA feasible |
| 35 | 69 | 3×23 | composite | OPEN | |
| 36 | 71 | prime | p ≡ 3(4) | OPEN | SA feasible |
| 37 | 73 | prime | pp ≡ 1(4) | **PROVEN** | Paley |
| 38 | 75 | 3×5² | composite | OPEN | |
| 39 | 77 | 7×11 | composite | OPEN | |
| 40 | 79 | prime | p ≡ 3(4) | OPEN | SA feasible |
| 41 | 81 | 3⁴ | pp ≡ 1(4) | **PROVEN** | Paley |
| 42 | 83 | prime | p ≡ 3(4) | OPEN | SA feasible |
| 43 | 85 | 5×17 | composite | OPEN | |
| 44 | 87 | 3×29 | composite | OPEN | |
| 45 | 89 | prime | pp ≡ 1(4) | **PROVEN** | Paley |
| 46 | 91 | 7×13 | composite | OPEN | |
| 47 | 93 | 3×31 | composite | OPEN | |
| 48 | 95 | 5×19 | composite | OPEN | |
| 49 | 97 | prime | pp ≡ 1(4) | **PROVEN** | Paley |
| 50 | 99 | 3²×11 | composite | OPEN | |

**Summary for n ≤ 50**:
- **PROVEN** (Paley infinite family): 15 values (n = 3,5,7,9,13,15,19,21,25,27,31,37,41,45,49)
- **Computationally verified** (SA): 16 values (n = 4,6,8,10,11,12,14,16,17,18,20,22,23,24,26,28,29,30)
- **OPEN**: 19 values (n = 32,33,34,35,36,38,39,40,42,43,44,46,47,48,50)

**Continuous verified range**: n = 3 through 31 (all proven or computationally verified).

---

## 5. Proof Methods

### Case 1: m = prime power ≡ 1 (mod 4) — PROVEN

The Paley construction over GF(q) gives an explicit, algebraic construction. This is a complete proof for infinitely many n, including:
- All n = (p+1)/2 where p is prime ≡ 1 (mod 4)
- All n = (p^k + 1)/2 where p^k ≡ 1 (mod 4) (e.g., GF(9), GF(25), GF(49), GF(81))

The proof uses the SRG property of the Paley graph: QR is an almost-difference set giving perfectly controlled autocorrelation values.

### Case 2: m = prime ≡ 3 (mod 4) — COMPUTATIONAL

For primes m ≡ 3 (mod 4), constructions are found via simulated annealing with structural constraints. Verified cases: m = 7,11,19,23,27,31,43,47,59.

**Towards a proof**: The Fourier analysis shows:
- The required fluctuation from average is ~1.5 (constant, independent of n)
- The standard deviation of Delta values grows as ~√n
- So the normalized difficulty DECREASES as n grows
- A probabilistic existence argument should work for sufficiently large n

**Key lemma needed**: For symmetric D11 ⊂ Z_p with |D11| = (p+1)/2, there exists D12 with |D12| = (p-1)/2 and 0 ∈ D12 such that A(d) + B(d) ≤ (p-1)/2 for all d ∈ D11.

### Case 3: m composite — COMPUTATIONAL

Composite cases are hardest due to lack of prime field structure. However, the fast SA approach works effectively:
- n=23 (m=45 = 9×5): solved in 74s
- n=26 (m=51 = 3×17): solved in 1421s
- n=28 (m=55 = 5×11): solved in 3259s
- n=29 (m=57 = 3×19): solved in 685s

The CRT decomposition Z_m ≅ Z_a × Z_b provides additional algebraic structure that could guide construction, but so far the "blind" SA search with universal structural constraints has been sufficient.

**Difficulty wall at n=32 (m=63 = 9×7)**: This case has resisted extensive SA search. It may require CRT-aware search strategies or a fundamentally different construction method.

---

## 6. Observations

### All constructions are extremally tight
Every verified construction hits both thresholds exactly:
- max red common neighbors = n - 2 (red threshold)
- max blue common neighbors = n - 1 (blue threshold)

There is zero slack. This suggests constructions are essentially unique (up to automorphism).

### SA search time scaling
For cases solved by fast SA, time roughly scales as O(m³) with high variance:

| m | Time (s) | Type |
|---|----------|------|
| 45 | 74 | composite |
| 47 | 128 | prime |
| 51 | 1,421 | composite |
| 55 | 3,259 | composite |
| 57 | 685 | composite |
| 59 | 4,918 | prime |

### Parity constraint
For odd m, symmetric D11 must have even cardinality. When the "standard" |D11| = n would be odd, the solver must use |D11| = n ± 1 instead. This affects cases where n is odd and m ≡ 1 (mod 4).

---

## 7. Tools Produced

| Tool | Description |
|------|-------------|
| `sa_solver.py` | Generalized SA solver for any n (parameterized from n=22 solver) |
| `fast_sa_general.py` | Fast joint SA with O(m²) incremental Delta updates |
| `validate_construction.py` | Standalone brute-force validator (no ramsey_core.py dependency) |
| `paley_general.py` | Paley construction for any GF(q) with q ≡ 1 (mod 4) |
| `validate_paley_bruteforce.py` | Standalone Paley validator with GF(q) arithmetic |
| `verify_v1v2_theorem.py` | Computational verification of V1V2 theorem |
| `n23_candidates.py` | D11 candidate analysis for composite m cases |
| `solutions_registry.json` | Complete registry of all verified constructions |

---

## 8. Open Problems

1. **General proof for primes ≡ 3 (mod 4)**: The Fourier analysis suggests a probabilistic existence argument should work. The key lemma (Section 5, Case 2) remains unproven.

2. **General construction for composite m**: No algebraic construction family is known. SA works computationally but doesn't constitute a proof. CRT-based constructions over product rings are a promising direction.

3. **The n=32 barrier**: m = 63 = 9×7 is the first case to resist our SA approach. Understanding why may reveal structural insights about composite cases with multiple small prime factors.

4. **Complete proof for all n**: Requires covering three cases — Paley (done), primes ≡ 3 mod 4 (needs existence argument), and composite m (needs construction family or case analysis).

---

## 9. References

- Rousseau & Sheehan (1978): Upper bound R(B_{n-1}, B_n) ≤ 4n - 1
- Lidicky et al. (2024): Verified range through n ≤ 21, known constructions for even n
- This work: Extended verified range to n ≤ 31, proved Paley infinite family, discovered V1V2 theorem
