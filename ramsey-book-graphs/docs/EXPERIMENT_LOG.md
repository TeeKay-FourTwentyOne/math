# Experiment Log: R(B_{n-1}, B_n) = 4n - 1

**Last updated:** February 8, 2026
**Target audience:** AI systems continuing this research

---

## Problem Statement

**Book graph B_n**: K_2 + K_n_bar (n triangles sharing a common edge, equivalently K_2 joined to an independent set of size n).

**Conjecture**: R(B_{n-1}, B_n) = 4n - 1 for all n >= 2.

**Known results**:
- Upper bound R(B_{n-1}, B_n) <= 4n - 1 proven by Rousseau & Sheehan (1978)
- Lower bound verified for n <= 21 (Wesley; Lidicky et al., 2024)
- Lower bound verified for all n where 2n-1 is a prime power with 2n-1 = 1 (mod 4) via Paley graph construction

**Current frontier**: n = 22 is the first open case. We need to show R(B_21, B_22) >= 87, i.e., find a 2-coloring of K_86 with no red B_21 and no blue B_22.

---

## Construction Framework: 2-Block Circulant Graphs

All known constructions for small n use **2-block circulant graphs** on N = 4n-2 vertices.

### Structure
- Vertex set V = V1 union V2, each of size m = 2n-1
- Vertices within and between blocks are connected by circulant structure
- Three difference sets parameterize the graph:
  - **D11** subset {1,...,m-1}: edges within V1 (must be symmetric: d in D11 iff m-d in D11)
  - **D12** subset {0,...,m-1}: edges from V1 to V2 (no symmetry requirement)
  - **D22** subset {1,...,m-1}: edges within V2 (must be symmetric)

### Degree structure
- All vertices in V1 have degree d1 = |D11| + |D12|
- All vertices in V2 have degree d2 = |D22| + |D12|
- Symmetry of D11 forces |D11| to be even; same for |D22|

### Common neighbor formulas
For a red edge (u,v), the number of common red neighbors (lambda_red) is:

- **V1V1 edge at difference d**: lambda_red = Delta(D11, D11, d) + Delta(D12, D12, d)
- **V2V2 edge at difference d**: lambda_red = Delta(D22, D22, d) + Delta(D12^T, D12^T, d)
- **V1V2 edge at difference d**: lambda_red = Sigma(D11, D12, d) + Delta(D12, D22, d)

Where:
- Delta(A, B, d) = |{a in A : (a-d) mod m in B}| (difference counting)
- Sigma(A, B, d) = |{a in A : (d-a) mod m in B}| (sum counting)
- D12^T = {-x mod m : x in D12}

### Constraint thresholds (for book-avoidance)
- **Red edges must avoid B_{n-1}**: lambda_red <= n-2 (= 20 for n=22)
- **Blue edges must avoid B_n**: lambda_blue <= n-1 (= 21 for n=22)

### Critical blue formula (CORRECTED)
Blue common neighbors are computed by inclusion-exclusion:

    lambda_blue = (N - 2) - deg_u - deg_v + lambda_red

This means for blue edges, lambda_red must satisfy:
- V1V1 blue: lambda_red <= 2*d1 - (N-2) + (n-1)
- V2V2 blue: lambda_red <= 2*d2 - (N-2) + (n-1)
- V1V2 blue: lambda_red <= d1 + d2 - (N-2) + (n-1)

**WARNING**: An early version of the code had lambda_blue = (N-2) - lambda_red (wrong). The correct formula is the one above. This was fixed before any of the experiments below.

### Complement constraint (historical)
All known solutions for n <= 21 use D22 = {1,...,m-1} \ D11 (complement). This forces d1 + d2 = (m-1) + 2*|D12|.

---

## For n=22: Key Parameters

- m = 43 (prime, 43 = 3 mod 4)
- N = 86 vertices
- red_threshold = 20 (lambda_red <= 20 for red edges)
- blue_threshold = 21 (lambda_blue <= 21 for blue edges)
- D11 has 21 pair slots (elements come in symmetric pairs {d, 43-d})
- D12 has 43 individual slots (0 through 42)
- D22 has 21 pair slots

### Quadratic residues mod 43
QR_43 = {1, 4, 6, 9, 10, 11, 13, 14, 15, 16, 17, 21, 23, 24, 25, 31, 35, 36, 38, 40, 41}

Properties:
- |QR_43| = 21
- Delta(QR_43, QR_43, d) = 10 for ALL d != 0 (perfectly uniform)
- QR_43 is NOT symmetric since 43 = 3 (mod 4), meaning -1 is a quadratic non-residue
- 12 of 21 pairs {d, 43-d} have their small element in QR

Since 43 = 3 (mod 4), the Paley graph (adjacency = QR) is a *tournament* (directed), not an undirected graph. The standard Paley construction for Ramsey lower bounds requires p = 1 (mod 4).

---

## Experiment 1: Complement SAT Solver (sat_n22_v3.py)

### Setup
- D22 = complement(D11), so only 21 + 43 = 64 primary variables
- Tests 7 configurations: |D11| in {18,20,22}, |D12| in {20,21,22}
- PySAT with CaDiCaL solver, cardinality constraints via sequential counter encoding
- Conditional encoding: for each difference d, encode "if d in D11 then lambda_red <= 20" and "if d not in D11 then lambda_red <= blue_threshold"

### Results
- Ran for 24+ hours on (|D11|=20, |D12|=21) configuration
- Terminated without returning a result (presumed UNSAT based on time)
- Problem size: ~80,000 clauses, ~64 primary variables but many auxiliary variables

### Assessment
The complement constraint forces degrees into a narrow band (d1 ~ 41-43, d2 ~ 41-43) which is the "density dead zone" identified below. Long runtime is consistent with exploring a large search tree before proving infeasibility.

---

## Experiment 2: Gemini Density Trap Critique

An external review (Gemini model) argued that the complement constraint D22 = complement(D11) forces the graph into a density regime where lambda_red ~ 21, barely exceeding the red threshold of 20. The critique recommended:

1. Decoupling D22 from D11 (add 21 more variables)
2. Targeting asymmetric degrees: d1 ~ 45, d2 ~ 40
3. Exploring skew-Hadamard constructions for p = 3 (mod 4)

**Verdict**: The critique was *directionally correct* (the density is indeed the problem) but the *proposed fix is insufficient* (see Experiments 3-5 below).

---

## Experiment 3: Skew-Hadamard Algebraic Exploration (skew_search.py)

### Setup
Systematic exploration of algebraically-motivated configurations:
- **Series A**: D12 = QR_43, D22 = complement(D11), varying |D11|
- **Series B**: D12 = QR_43, D11 = D22 (equal blocks)
- **Series C**: D12 = QR_43, independent D11 and D22 at various densities
- **Series D**: D12 = QNR_43, complement D22
- **Series E**: D12 = QR_43 union {0}, complement D22
- **Series F**: D12 = QNR_43 union {0}, complement D22

### Results
Best algebraic configuration: 28 violations, achieved at |D11|=22, |D12|=21, |D22|=20.

No algebraic seed produced fewer than 28 violations. The QR-based D12 gives perfectly uniform Delta(D12,D12,d) = 10, but this doesn't help D11/D22 enough because those must be symmetric while QR_43 is not.

### Key insight
The algebraic regularity of QR_43 (constant Delta) is helpful for cross-block terms but cannot compensate for the within-block density obstruction. The problem is not the choice of D12 — it's the fundamental density constraints on D11 and D22.

---

## Experiment 4: Simulated Annealing Comparison

### 4a: Complement D22 (D22 = complement(D11))
- 8 trials x 500,000 iterations, SA with T_init=10, cooling=0.99997
- **Best: cost=16, 16 violations**
- Configuration: |D11|=20, |D12|=21, |D22|=22, d1=41, d2=43
- All violations have excess exactly 1
- Violation pattern: 8 at V1V1 shifts, 8 at V2V2 shifts (mirrored)
- Violating shifts: {1, 9, 18, 20, 23, 25, 34, 42}

Note: A prior run with different parameters achieved cost=8 with only 8 violations at shifts {2, 21, 22, 41} — this appears to be the true global optimum for complement D22.

### 4b: Independent D22 (free variables)
- Same SA parameters
- **Best: cost=44, 33 violations** — MUCH WORSE
- Degrees still converge to d1~42, d2~42

### 4c: D12 = QR_43 fixed, independent D22
- **Best: cost=36, 20 violations** — intermediate

### Summary table

| D22 mode       | Best cost | Violations | Degrees      |
|----------------|-----------|------------|--------------|
| Complement     | 8*        | 8          | d1=43, d2=41 |
| Complement     | 16        | 16         | d1=41, d2=43 |
| Independent    | 44        | 33         | d1=42, d2=42 |
| QR + independent| 36       | 20         | d1=43, d2=41 |

*From earlier session with longer search.

### Key finding
Independent D22 makes things WORSE, not better. The complement constraint actually helps by correlating D11 and D22 in a way that partially cancels violations between the two blocks. With independent D22, the search space is larger but the obstruction is the same, and the correlation benefit is lost.

### Universal convergence
All SA runs, regardless of initialization or D22 mode, converge to d1 in {41, 42, 43} and d2 in {41, 42, 43}. This is not a heuristic artifact — it reflects the only feasible degree regime (see density analysis below).

---

## Experiment 5: Decoupled SAT Solver (sat_decoupled.py)

### Setup
- Independent D11, D12, D22 variables: 21 + 43 + 21 = 85 primary variables
- Fixed cardinalities per run to make degrees constant
- Blue thresholds become constants (not conditional on which elements are chosen)
- Sweeps configurations sorted by distance from (20, 21, 20)

### Results (partial)
- Configuration (20, 21, 20): 459,222 clauses, 222,495 variables
- Solver started but terminated before completing (process killed after several hours)
- No SAT result obtained for any configuration

### Assessment
The decoupled SAT has ~3x more clauses than the complement SAT due to independent D22 constraints. Combined with the average-case impossibility (Experiment 6), a complete run would likely return UNSAT.

---

## Experiment 6: Average-Case Impossibility Analysis

### The density obstruction (mathematical proof)

For V1V1 edges, the expected lambda_red over random shifts is:

    E[lambda_red] ~ (|D11|^2 + |D12|^2) / m

This must satisfy BOTH:
1. **Red bound**: E[lambda_red] <= 20 (for d in D11)
2. **Blue bound**: E[lambda_red] <= 2*d1 - 63 (for d not in D11)

Let s = |D11| + |D12| = d1 and p = |D11| * |D12|.
Then |D11|^2 + |D12|^2 = s^2 - 2p.

- Red requires: p >= (s^2 - 860) / 2
- Blue requires: p >= (s^2 - 86*s + 2709) / 2
- Maximum achievable: p_max = floor(s/2) * ceil(s/2)

| s (=d1) | Red: p >= | Blue: p >= | p_max | Red OK? | Blue OK? |
|---------|-----------|------------|-------|---------|----------|
| 39      | 310       | 385        | 380   | Yes     | No       |
| 40      | 370       | 406        | 400   | Yes     | No       |
| 41      | 411       | 432        | 420   | Yes     | No       |
| 42      | 452       | 462        | 441   | No      | No       |
| 43      | 495       | 495        | 462   | No      | No       |
| 44      | 538       | 532        | 484   | No      | No       |
| 45      | 583       | 573        | 506   | No      | No       |

**For EVERY integer value of d1, at least one constraint is infeasible on average.**

The red upper bound requires s <= 41 (approximately), while the blue lower bound requires s >= 42 (approximately). These ranges don't overlap.

### Comprehensive density sweep
A sweep over all feasible (|D11|, |D12|, |D22|) triples in range 14-28 x 16-27 x 14-28 found **zero configurations** where expected lambda_red is within all bounds simultaneously.

### Why difference sets can't save us
The only escape from the average-case argument is if Delta(D, D, d) is constant (not fluctuating around the mean). This requires a **difference set**.

For m = 43:
- (43, 21, 10) difference set exists (QR_43), but |QR| = 21 is odd, and D11 must have even cardinality (symmetric pairs). Also QR_43 is NOT symmetric since 43 = 3 (mod 4).
- (43, 20, ?): would need 20*19 = 42*lambda, lambda = 380/42 is not an integer. No such difference set exists.
- No symmetric difference set exists at any useful parameter for m = 43.

### Root cause
For primes p = 3 (mod 4), -1 is a quadratic non-residue, so the quadratic residues are never symmetric. The QR set is the only family of difference sets for primes, and it can't be used for D11 or D22 (which require symmetry). This is why the Paley construction works for p = 1 (mod 4) but fails for p = 3 (mod 4).

---

## CRITICAL UPDATE: Average-Case Analysis Was Misleading

### The average-case argument does NOT rule out solutions

Analysis of ALL known m = 3 (mod 4) constructions (n=6,8,10,12,14,16,18,20) reveals that **every single one has the same average-case "impossibility"**:
- The dense block always has avg_lambda exceeding red_threshold by ~1.5
- The sparse block always has avg_lambda exceeding blue_bound by ~0.6
- Yet valid constructions exist for ALL of these cases

The constructions succeed because Delta(D,D,d) fluctuates around the mean. The difference set D11 is chosen so red edges (d in D11) land at differences where lambda is below the threshold, while blue edges (d not in D11) land at differences where lambda is above average (which is fine for blue).

### Universal structural pattern for m = 3 (mod 4)

ALL known constructions for m = 3 (mod 4) share:
- D22 = complement(D11) (ALWAYS)
- 0 in D12 (ALWAYS)
- Degrees: {d1, d2} = {m, m-2}
- Sizes: |D11| = (m+1)/2 or (m-3)/2; |D12| = (m-1)/2
- Max red = n-2 exactly (tight), max blue = n-1 exactly (tight)

### Extrapolation to n=22

The excess and fluctuation range both scale predictably:

| n  | excess over red | V1V1 range | V2V2 range | V1V2 range |
|----|-----------------|------------|------------|------------|
| 12 | 1.52            | 4          | 3          | 1          |
| 14 | 1.52            | 4          | 3          | 1          |
| 16 | 1.52            | 5          | 3          | 1          |
| 18 | 1.51            | 5          | 3          | 1          |
| 20 | 1.51            | 6          | 6          | 1          |
| 22 | 1.51            | ~6-7?      | ~6-7?      | 1?         |

**A solution for n=22 is likely feasible** with the parameters |D11|=22, |D12|=21, |D22|=20, d1=43, d2=41 (or the flipped version).

### Revised conclusions

1. **The average-case "proof" of impossibility was wrong** — it proved the average is infeasible, but ALL prior cases are also average-infeasible yet have solutions
2. **The complement constraint IS correct** — every known construction uses it
3. **The structural pattern is universal** — and it predicts exactly the configuration the SA converges to
4. **The 8-violation barrier is a search failure, not a mathematical impossibility** — the SA gets close but can't close the gap
5. **A more powerful search (SAT, or improved heuristic) may find the solution**

### Previous conclusions (partially retracted)

~~The average-case density analysis proves no cardinality configuration is feasible on average.~~ This is technically true but misleading — all prior cases are also average-infeasible.

The complement constraint is indeed better than independent D22 (Experiment 4 confirms this). The Gemini critique's recommendation to decouple D22 was wrong.

---

## Experiment 7: Analysis of Known Constructions (NEW)

### Data source
Constructions from Steven-VO/circulant-Ramsey GitHub and Lidicky et al. (2024).
All m = 3 (mod 4) cases n=6,8,10,12,14,16,18,20 verified.

### Universal pattern discovered
ALL known m = 3 (mod 4) constructions share:
- D22 = complement(D11)
- 0 in D12
- Degrees {d1, d2} = {m, m-2}
- |D11| = (m+1)/2 or (m-3)/2, |D12| = (m-1)/2
- Max red = n-2 (tight), max blue = n-1 (tight)

### Average-case analysis: all cases exceed threshold
| n  | dense_avg | red_thresh | excess | V1V1_range | valid |
|----|-----------|------------|--------|------------|-------|
| 6  | 3.73*     | 4          | -0.27* | 2          | YES   |
| 8  | 5.67*     | 6          | -0.33* | 2          | YES   |
| 10 | 7.63*     | 8          | -0.37* | 2          | YES   |
| 12 | 11.52     | 10         | +1.52  | 4          | YES   |
| 14 | 13.52     | 12         | +1.52  | 4          | YES   |
| 16 | 15.52     | 14         | +1.52  | 5          | YES   |
| 18 | 17.51     | 16         | +1.51  | 5          | YES   |
| 20 | 19.51     | 18         | +1.51  | 6          | YES   |
| 22 | 21.51     | 20         | +1.51  | ???        | ???   |

*For n<=10, the "dense" block is V2V2 (d2=m), and V1V1 average is below threshold.
For n>=12, the pattern flips: V1V1 is the dense block.

### Key finding
The +1.51 excess is CONSTANT across all n >= 12. The fluctuation range grows
slowly (4, 4, 5, 5, 6, ...). All prior cases succeed despite the average
exceeding the threshold. **There is no reason to believe n=22 is fundamentally
different.**

### SOLUTION FOUND (Experiment 7b)

The informed SA (sa_n22_informed.py) found a valid construction on its FIRST trial
after only 369,287 iterations (~5 minutes)!

**The solution:**
```
D11 = [1, 2, 5, 10, 11, 13, 16, 17, 18, 19, 20, 23, 24, 25, 26, 27, 30, 32, 33, 38, 41, 42]
D12 = [0, 2, 5, 6, 8, 11, 15, 16, 20, 24, 25, 27, 28, 31, 32, 34, 35, 36, 37, 39, 41]
D22 = [3, 4, 6, 7, 8, 9, 12, 14, 15, 21, 22, 28, 29, 31, 34, 35, 36, 37, 39, 40]
```

Properties:
- |D11| = 22, |D12| = 21, |D22| = 20 (matches universal pattern exactly)
- d1 = 43, d2 = 41 (matches universal pattern exactly)
- D22 = complement(D11) (as predicted)
- 0 in D12 (as predicted)
- Max red common neighbors: 20 (threshold 20) - TIGHT
- Max blue common neighbors: 21 (threshold 21) - TIGHT
- Zero violations

**Verified by complete brute-force:** all 3,655 vertex pairs checked.

### Why the solution was found now but not before

Key improvements from analyzing known constructions:
1. **Fixed cardinalities**: swap moves maintain |D11|=22, |D12|=21 exactly (no time wasted exploring wrong densities)
2. **Fixed 0 in D12**: eliminates one degree of freedom that was never useful
3. **Correct structural constraints**: complement D22, symmetric D11 pairs

---

## Code Inventory

| File | Purpose | Status |
|------|---------|--------|
| `ramsey_core.py` | Core library: BlockCirculantGraph, Delta/Sigma, verify_construction | Verified correct |
| `solver_heuristic.py` | Tabu search and SA with parallel execution | Working |
| `solver_sat.py` | General SAT solver (red constraints only, no blue) | Incomplete |
| `sat_n22_v3.py` | SAT with complement D22, fixed cardinalities | Working, likely UNSAT |
| `sat_decoupled.py` | SAT with independent D11/D12/D22 | Working, untested to completion |
| `skew_search.py` | Algebraic seeds + SA for skew-Hadamard exploration | Complete |
| `analyze_known.py` | Verify known constructions, analyze Delta patterns | Complete |
| `sat_n22_targeted.py` | Targeted SAT with pattern-informed constraints | Superseded (SA found solution first) |
| `sa_n22_informed.py` | Fixed-cardinality SA with swap moves | Complete (found solution) |
| `validate_n22_full.py` | Standalone brute-force validator (no dependencies) | Complete (PASS) |

---

## Next Directions (prioritized)

### Priority 1: Study existing constructions for small m = 3 (mod 4)
For n <= 21, some values have m = 2n-1 = 3 (mod 4): m=3,7,19,23,31. How were those cases solved? If by exhaustive search on small instances, that doesn't help n=22. If by a specific algebraic construction, we may be able to generalize.

### Priority 2: Cayley graphs on non-cyclic groups
The dihedral group D_43 (order 86) is the main candidate. It has non-commutative elements and different algebraic structure from Z_86. A Cayley graph Cay(D_43, S) has the same number of vertices but a fundamentally different adjacency structure. The SA/SAT framework can be adapted.

### Priority 3: Multi-block or non-circulant structures
86 = 2 x 43 has no 3-way factorization, so k-block circulant for k >= 3 is not directly applicable. However, other structured graph families (e.g., based on finite geometry) might work.

### Priority 4: Extension field constructions
F_{43^2} has 1849 elements. It's unclear how to extract a graph on exactly 86 vertices, but substructures (e.g., certain cosets or subfields) might yield useful constructions.

### De-prioritized: Different n values
While solving easier n values (where 2n-1 = 1 mod 4) is straightforward via Paley, it doesn't advance the frontier since those cases are already covered by the known infinite family. The challenge is specifically the p = 3 (mod 4) cases.

---

## Experiment 8: Final Standalone Validation

### Purpose
Independent brute-force verification with zero dependencies on `ramsey_core.py`.
The script `validate_n22_full.py` builds the full 86x86 adjacency matrix from the
difference sets and directly counts common neighbors for every pair.

### Method
1. Build adjacency matrix: V1={0..42}, V2={43..85}
   - V1-V1: edge iff (v-u) mod 43 in D11
   - V2-V2: edge iff (v-u) mod 43 in D22
   - V1-V2: edge iff (v-43-u) mod 43 in D12
2. Verify structural properties (symmetry, complement, cardinalities, degrees)
3. For each of C(86,2)=3655 pairs, count common neighbors in the edge's color
   by iterating all 84 other vertices (no formulas, no shortcuts)

### Output
```
============================================================
STANDALONE VALIDATOR: R(B_21, B_22) = 87
============================================================

--- Structural Properties ---
|D11| = 22 (expected 22): OK
|D12| = 21 (expected 21): OK
|D22| = 20 (expected 20): OK
D11 symmetric: OK
D22 symmetric: OK
D22 = complement(D11): OK
0 in D12: OK
d1 = |D11| + |D12| = 43
d2 = |D22| + |D12| = 41

--- Building adjacency matrix ---
V1 degrees: {43} (expected {43})
V2 degrees: {41} (expected {41})

--- Checking all pairs ---
Total pairs checked: 3655
Red edges: 1806
Blue edges: 1849
Max red common neighbors:  20 (threshold 20)
Max blue common neighbors: 21 (threshold 21)
Violations: 0

============================================================
RESULT: PASS

The 2-block circulant graph on 86 vertices contains
no red B_21 and no blue B_22.
Therefore R(B_21, B_22) >= 87.
Combined with the upper bound R(B_21, B_22) <= 87,
this proves R(B_21, B_22) = 87.
============================================================
```

### Summary
- 3655 pairs checked exhaustively
- 1806 red edges, 1849 blue edges
- Max red common neighbors: 20 (= threshold, tight)
- Max blue common neighbors: 21 (= threshold, tight)
- Zero violations
- **PASS**: R(B_21, B_22) = 87 is confirmed

---

## Appendix: File Locations

- Problem statement: `ramsey-book-graphs.pdf`
- Code: `ramsey-book-graphs/*.py`
- Documentation: `ramsey-book-graphs/docs/*.md`
- Standalone validator: `ramsey-book-graphs/validate_n22_full.py`
- This file: `ramsey-book-graphs/docs/EXPERIMENT_LOG.md`
