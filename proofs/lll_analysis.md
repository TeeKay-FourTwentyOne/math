# Lovász Local Lemma Analysis for Book Graph Ramsey Numbers

**Date**: 2026-02-15
**Purpose**: Assess whether the Lovász Local Lemma (LLL) or related probabilistic methods can prove existence of valid (D11, D12) pairs for all sufficiently large primes p ≡ 3 mod 4.

## 1. Problem Setup

### 1.1 Goal
Prove R(B_{n-1}, B_n) = 4n-1 for all n where m = 2n-1 = p is prime, p ≡ 3 mod 4.

### 1.2 Construction
- 2-block circulant graph on 2m = 4n-2 vertices
- Partition: V1 (size p), V2 (size p-2)
- Difference sets:
  - D11 ⊂ {1,...,p-1}, symmetric, |D11| = (p+1)/2
  - D12 ⊂ {0,...,p-1}, |D12| = (p-1)/2
  - D22 determined by D11 via diagonal complement

### 1.3 Constraints
For each d ∈ {1,...,p-1}:
- If d ∈ D11: A(d) + B(d) ≤ n-2 AND C(d) + B(p-d) ≤ n-3
- If d ∈ D22: C(d) + B(p-d) ≤ n-2 AND A(d) + B(d) ≤ n+1

Where:
- A(d) = |{x ∈ D11 : x+d mod p ∈ D11}|
- B(d) = |{x ∈ D12 : x+d mod p ∈ D12}|
- C(d) = |{x ∈ D22 : x+d mod p ∈ D22}|

For symmetric D11 with d ∈ D11: C(d) = A(d) - 1, so constraint reduces to:
**A(d) + B(d) ≤ n-2** for all d ∈ D11

Plus V1V2 cross-constraints (p total differences including 0).

### 1.4 Key Properties
- Expected A(d) = (p+1)/4 ≈ n/2 for random symmetric D11
- Expected B(d) ≈ (p-3)/4 ≈ n/2 - 1 for random D12
- Threshold: n-2 ≈ p/2 - 1
- Expected A(d) + B(d) ≈ p/2 - 1/2, **just barely below threshold**
- Variance: σ²(A+B) ≈ p/4, so typical fluctuation ~ √p/2 ~ n^(1/2)

**CRITICAL BARRIER**: Constraint violation probability ≈ 1/2 for random (D11, D12).

## 2. Standard LLL Conditions

### 2.1 Symmetric LLL (Erdős-Lovász 1975)
**Condition**: If each bad event E_i has:
- Pr[E_i] ≤ p for all i
- E_i is mutually independent of all but d other events
- e·p·(d+1) ≤ 1

Then Pr[⋂ Ē_i] > 0 (no bad events occur).

### 2.2 Application to This Problem

**Bad events**: E_d = "constraint at d is violated"
- Number of events: ≈ p (one per difference)
- Pr[E_d] ≈ 1/2 (expected A+B ≈ threshold)
- Dependency: each E_d depends on ALL other events (through global D11, D12)

**LLL check**:
- e·p·(d+1) ≈ 2.718 · (1/2) · (p+1) ≈ 1.36p >> 1 for large p

**VERDICT**: Standard symmetric LLL **FAILS** — condition violated by factor ~p.

### 2.3 Shearer's Optimal Condition (1985)
For a dependency graph G and probability vector p, the optimal condition is:

Pr[⋂ Ē_i] > 0 if q_S(p₁,...,p_n) > 0 for all S ⊆ [n]

where q_S is the alternating-sign independence polynomial.

For dense dependency graphs (like ours, where everything depends on everything), this doesn't improve the bound significantly.

**Source**: [Directed Lovász Local Lemma and Shearer's Lemma](https://arxiv.org/abs/1611.00502)

## 3. Alternative LLL Variants

### 3.1 Lopsided LLL (Erdős-Spencer)
**Idea**: Replace mutual independence with negative dependency.

**Negative dependency graph**: Edge (E_i, E_j) exists if occurrence of E_i can INCREASE Pr[E_j].

**Advantage**: Often sparser than mutual dependency graph.

**Application to this problem**:
- Fix D11, vary D12
- Question: Are constraints negatively dependent?
- For d, d' ∈ D11: both constraints involve B(d), B(d')
- If B(d) is large, does that make B(d') more or less likely to be large?

**Analysis**: The B(d) values are correlated through the global constraint Σ B(d) = |D12|(|D12|-1). However:
- If B(d) is large at one d, the "budget" for other d' is reduced
- This suggests **negative correlation** — potentially helpful!

**Challenge**: Even with negative dependency, need:
- Each event has small negative dependency neighborhood
- Our problem: B(d) values all interact through the single global D12 choice

**VERDICT**: Lopsided LLL might help reduce dependency graph density, but unlikely to overcome p ≈ 1/2 barrier with dense dependencies.

**Sources**:
- [Probabilistic Methods in Combinatorics — Yufei Zhao](https://yufeizhao.com/pm/6.pdf)
- [Lecture on Lopsided LLL](https://www.cs.ubc.ca/~nickhar/W15/Lecture22Notes.pdf)

### 3.2 Cluster Expansion LLL (Bissacot et al. 2011)
**Improvement**: For dependency graph G, restrict Shearer's condition to independent sets in G.

**When it helps**: Dependency graph has clique structure (neighborhoods form cliques).

**Application to this problem**:
- Our dependency graph is nearly complete (almost every constraint depends on every other)
- No significant clique structure to exploit

**VERDICT**: Cluster expansion unlikely to provide substantial improvement.

**Sources**:
- [An Improvement of the Lovász Local Lemma via Cluster Expansion](https://arxiv.org/abs/0910.1824)
- [Stanford Lecture on Cluster Expansion](https://theory.stanford.edu/~jvondrak/MATH233A-2018/Math233-lec06.pdf)

## 4. Entropy Compression (Moser-Tardos 2010)

### 4.1 Method
Algorithmic approach to LLL: resample violated constraints until all satisfied.

**Key insight**: If the "information compression" rate is positive, the algorithm terminates.

**Condition**: Same as standard LLL (e·p·(d+1) ≤ 1), but provides constructive algorithm.

### 4.2 Application to This Problem
**Resampling strategy**:
1. Choose random (D11, D12)
2. If constraint at d violated, resample... what?
   - D11? Changes A(d) at d, but also A(d') for all other d'
   - D12? Changes B(d) at d, but also B(d') for all other d'

**Challenge**: Each resample affects ALL constraints, not just the violated one. This is the "global coupling" problem.

**VERDICT**: Entropy compression doesn't overcome the fundamental barrier (p ≈ 1/2, dense dependencies). It's an algorithmic framework, not a relaxation of the LLL conditions.

**Sources**:
- [Moser's entropy compression argument (Terence Tao)](https://terrytao.wordpress.com/2009/08/05/mosers-entropy-compression-argument/)
- [Moser-Tardos resampling algorithm](https://arxiv.org/abs/2001.00880)

## 5. Alternative Probability Spaces

### 5.1 Fix D11, Apply LLL to D12
**Strategy**: Instead of randomizing both D11 and D12:
1. Fix a "good" D11 (e.g., from a working orbit)
2. Randomize only D12
3. Apply LLL to find valid D12

**Bad events**: E_d = "B(d) is too large" (for d where A(d) is already determined)

**Analysis**:
- Number of constraints: p
- For d ∈ D11 with A(d) = a: need B(d) ≤ n-2-a
- Expected B(d) ≈ (p-3)/4
- If a is near its expectation (p+1)/4, then threshold is ≈ p/4 - 5/4
- Constraint violation probability ≈ 1/2 (still!)

**Dependency**:
- B(d) values are independent if we condition on size |D12|
- But |D12| = (p-1)/2 is fixed!
- So B(d) are negatively correlated (budget constraint)

**Problem**: Even with negative correlation, p ≈ 1/2 is too large for LLL.

**Possible improvement**: Choose D11 such that A(d) values are ALL small (say, ≤ n/4). Then thresholds for B(d) are larger, reducing violation probabilities.

**Obstruction**: Constant z-sum constraint: Σ_{d∈D11} A(d) = n(n-1)/2 exactly for symmetric D11. If all A(d) are small, we need many of them — but |D11| = n, so average A(d) = (n-1)/2. Can't make all A(d) small!

### 5.2 Fix D12, Apply LLL to D11
**Strategy**:
1. Fix a "good" D12
2. Randomize symmetric D11
3. Apply LLL to find valid D11

**Challenge**: Symmetric D11 is not a simple k-subset — it's a k/2-subset of multiplicative orbits {d, p-d}.

**Analysis**:
- For each d ∈ {1,...,(p-1)/2}, decide whether to include orbit {d, p-d} in D11
- This gives 2^((p-1)/2) ≈ 2^(n/2) possible symmetric D11
- Each d requires: A(d) + B(d) ≤ n-2 if d ∈ D11
- A(d) depends on ALL orbit inclusion decisions (global coupling)

**Dependency graph**: Still dense (nearly complete).

**VERDICT**: Doesn't improve the situation.

### 5.3 Random D11 from a Restricted Family
**Strategy**: Randomize D11 from a "nice" family (e.g., QR-based, low DFT flatness).

**Obstruction**: Prior work shows:
- QR-guided D11 construction produces non-working orbits (algebraic invariants don't separate)
- No single algebraic invariant separates working from non-working at p ≥ 19
- Fixed thresholds on (dft_flatness, a_d11_var) fail across primes

**VERDICT**: No known restricted family with guaranteed good properties.

## 6. Alternative Probabilistic Methods

### 6.1 Second Moment Method
**Already explored**: Double-averaging analysis shows:
- E[N] grows (3.32, 0.36, 3.24, 5.82 bits for p = 11, 19, 23, 31)
- E[N²]/E[N]² grows (2.0, 14.0, 14.5, 71.25)
- Paley-Zygmund bound: Pr[N > 0] ≥ 1/(E[N²]/E[N]²)

**Challenge**: Second moment ratio grows roughly linearly with p, so PZ bound decays as ~1/p. This gives existence for each specific p, but:
- Need ratio = O(poly(p)) for all large p (not proven)
- Empirical data stops at p=31 (ratio = 71.25)
- Product-measure approach fails at p ≈ 230 (margin/p → -1)

**Verdict**: Second moment method works empirically up to p=31, but proof of ratio = O(poly(p)) for all p is the OPEN PROBLEM (exactly the L6 gap!).

### 6.2 Rödl Nibble / Semi-Random Method
**Idea**: Build the object incrementally, "nibble" by "nibble," carefully controlling properties.

**Classic application**: Rödl (1985) proved Erdős-Hanani conjecture on approximate Steiner systems.

**Application to this problem**:
- Build D12 incrementally: start with {0}, add elements one at a time
- At each step, ensure constraints remain satisfiable

**Challenge**:
- Constraint marginals are coupled through D11
- Adding one element to D12 affects B(d) for all d simultaneously
- Hard to maintain "locally good" properties while growing D12

**Verdict**: Rödl nibble is powerful for hypergraph problems with LOCAL constraints. Our constraints are GLOBAL (coupling through D11, D12). Unclear if nibble strategy applies.

**Sources**:
- [Ramsey Theory and the Semi-Random Method](https://arxiv.org/html/2510.19978)
- [Random Ramsey Theorem](https://www.cambridge.org/core/journals/combinatorics-probability-and-computing/article/abs/short-proof-of-the-random-ramsey-theorem/1F00D51A100FEEC7CD7814EC2F76AA83)

### 6.3 Variance Method / Azuma-Hoeffding
**Idea**: Use concentration inequalities to show that constraint violations are rare.

**Condition**: Need exponential concentration (e.g., Pr[A+B > n-2] ≤ e^(-c·p)).

**Analysis**:
- A + B is a sum of dependent 0-1 variables
- Variance: Var(A+B) ≈ p/4
- Standard deviation: σ(A+B) ≈ √p/2
- Threshold distance: (n-2) - E[A+B] ≈ -1/2 < 0 (in wrong direction!)

**Problem**: Threshold is BELOW expectation, so we need LOWER tail bound, not upper tail. But lower tail is fat (not exponentially concentrated).

**Verdict**: Concentration inequalities don't help when threshold ≈ expectation and we need lower tail.

### 6.4 Algebraic / Finite Field Methods
**Classic applications**:
- Paley construction for p ≡ 1 mod 4 (already proven!)
- Graph power constructions (Alon-Roichman, Boppana)
- Polynomial method for Ramsey bounds

**Application to this problem**:
- For p ≡ 3 mod 4, QR is not a difference set
- Character sum analysis shows no algebraic invariant separates working from non-working
- DFT analysis shows spectral complementarity (anti-correlated D11, D12 spectra), but doesn't yield constructive proof

**Verdict**: No known algebraic construction for p ≡ 3 mod 4. This is why we resort to probabilistic/computational search.

**Sources**:
- [Algebraic Methods for Ramsey Lower Bounds](https://arxiv.org/pdf/math/0608013)
- [Finite Geometry and Ramsey Theory Mini-course](https://anuragbishnoi.wordpress.com/minicourse/)

## 7. Why LLL Fails: Root Cause Analysis

### 7.1 The Fundamental Barrier
The LLL requires **p·(d+1) ≤ 1/e**, where:
- p = Pr[single constraint violated]
- d = max dependency degree

For our problem:
- p ≈ 1/2 (threshold ≈ expectation)
- d ≈ p (nearly complete dependency graph)
- Product: p·(d+1) ≈ p/2 >> 1

**Gap**: We're off by a factor of ~p/2 ≈ n/4, which grows linearly with problem size!

### 7.2 Why Violation Probability is 1/2
The constraint A(d) + B(d) ≤ n-2 has:
- E[A(d)] ≈ n/2 (for random symmetric D11)
- E[B(d)] ≈ n/2 - 1 (for random D12)
- Threshold: n-2

**Slack**: E[A+B] - (n-2) ≈ -1/2

This tiny negative slack (only -1/2, not -Θ(n)) means the constraint is "barely satisfied on average." With Var(A+B) ≈ n/2, the standard deviation is ~√n/2, which is >> 1/2. So Pr[violation] ≈ Φ(1/√(2n)) ≈ 1/2 - O(1/√n) ≈ 1/2.

**Why this is tight**:
- Can't increase slack by adjusting expectations (constant z-sum constraint fixes E[A])
- Can't decrease variance significantly (coupling through global D11, D12)

### 7.3 Why Dependency Graph is Dense
Each constraint A(d) + B(d) ≤ n-2 involves:
- The global D11 (affects A(d) for all d)
- The global D12 (affects B(d) for all d)

Changing D11 to satisfy constraint at d affects A(d') at all other d'. Similarly for D12 and B(d').

**Implication**: The dependency graph is nearly COMPLETE — every constraint depends on every other constraint.

### 7.4 Product of Barriers
LLL fails because:
1. High violation probability (p ≈ 1/2) AND
2. Dense dependencies (d ≈ p)
3. Product p·(d+1) ≈ p/2 grows linearly

**To succeed, need**:
- Either p = o(1/d) ≈ o(1/p) → p = o(1/√p) → impossible!
- Or reduce dependency graph density significantly

## 8. Empirical vs. Theoretical Gap

### 8.1 What We Know Empirically
From p=11,19,23,31,43,47,59 data:
- Valid (D11, D12) pairs exist at all tested primes
- p_working = fraction of working orbits declines (5.50, 1.36, 2.74, 1.30, 0.32 for p×p_working)
- But at p=43: still 124/16796 working orbits (0.74% of all orbits)
- Each working orbit has (p-1)/2 = 21 distinct D11 choices
- So at p=43: 124 · 21 = 2604 working D11 sets (out of C(21, 11) = 352,716 symmetric D11)

**Implication**: Existence is not rare — thousands of solutions exist even at p=43!

### 8.2 What LLL Predicts
If we treat constraints as independent with p = 1/2:
- Pr[all p constraints satisfied] ≈ (1/2)^p = 2^(-p) → 0 exponentially fast

**Reality**:
- Pr[N > 0] = 0.5, 0.0714, 0.119, 0.042 for p = 11, 19, 23, 31
- Decays polynomially (roughly ~1/p), not exponentially

**Key difference**: Constraints are NOT independent! The negative correlations (budget constraints, spectral complementarity) provide substantial "savings" over the independence assumption.

### 8.3 The Missing Theory
We need a theoretical framework that:
1. Accounts for negative correlations between constraints
2. Proves Pr[all constraints satisfied] ≥ c/p^α for constants c > 0, α < ∞
3. Scales to all sufficiently large p

**LLL doesn't do this** because it assumes WORST-CASE dependencies (or at best, negative independence graph). Our problem has:
- Negative correlations (budget constraints)
- But also positive correlations (coupling through global D11, D12)
- Net effect: polynomial decay, not exponential

**No existing probabilistic method captures this mixed dependency structure.**

## 9. Recommendations

### 9.1 LLL and Variants: Not Promising
**Conclusion**: Standard LLL, Lopsided LLL, Cluster Expansion LLL, and Entropy Compression all fail because:
- Violation probability p ≈ 1/2 is too high
- Dependency graph is too dense
- Product p·(d+1) ≈ p/2 >> 1

**Do NOT pursue**: Further investigation of LLL-based approaches is unlikely to succeed.

### 9.2 Second Moment Method: Promising but Incomplete
**Current status**:
- Works empirically up to p=31 (E[N²]/E[N]² = 71.25 = O(p))
- Product-measure approach fails at p ≈ 230

**Open problem**: Prove E[N²]/E[N]² = O(poly(p)) for all large p.

**Key challenge**: Controlling the second moment requires understanding:
- Overlap structure of valid D12 sets for each D11
- Correlation between N(D11) values across D11 orbits
- How spectral complementarity affects the variance

**This is exactly the L6 gap identified in prior work.**

### 9.3 Alternative Directions

#### 9.3.1 Refined Second Moment Analysis
**Goal**: Decompose E[N²] to identify dominant terms.

**Approach**:
- E[N²] = E[N] + Σ_{D12 ≠ D12'} Pr[both valid | D11]
- Overlap decomposition: classify pairs by |D12 ∩ D12'|
- Use DFT/spectral properties to bound correlations

**Precedent**: Prior work on "overlap_decomposition.py" (see MEMORY.md).

**Next steps**:
- Compute overlap correlation structure at p=31,43
- Identify whether high-overlap pairs contribute disproportionately to E[N²]
- Prove structural bounds on overlap-weighted sums

#### 9.3.2 Spectral Complementarity as a Constructive Principle
**Observation**: Valid pairs have anti-correlated spectra (corr(|D̂11|², |D̂12|²) = -0.82 to -0.87 for p=47,59).

**Goal**: Turn this into a constructive algorithm.

**Approach**:
1. Given D11, compute its spectrum |D̂11|²
2. Search for D12 with spectrum |D̂12|² that:
   - Has negative correlation with |D̂11|²
   - Satisfies |D̂11|² + |D̂12|² ≈ constant (combined flatness)
3. Use iterative refinement (gradient descent, simulated annealing)

**Challenge**: Spectrum-to-set mapping is non-injective (phase problem). Need algorithmic way to navigate this.

#### 9.3.3 Conditional Variance Method
**Goal**: Condition on "easy" structure, prove variance is small in residual.

**Approach**:
- Partition D11 space into "nice" and "generic" parts
- For nice D11, prove directly that N(D11) > 0 (e.g., algebraic construction)
- For generic D11, prove E[N² | generic] / E[N | generic]² = O(poly(p))

**Challenge**: Defining "nice" partition that:
- Covers enough probability mass
- Admits clean variance bounds

#### 9.3.4 Algorithmic Proof via SA Certification
**Radical idea**: Instead of probabilistic existence proof, use simulated annealing as the algorithm, then CERTIFY its success probability.

**Approach**:
1. Prove that SA with parameters (T₀, cooling schedule, iterations) succeeds with probability ≥ δ > 0
2. Use landscape analysis: show energy barriers are O(poly(p)), acceptance ratio stays Ω(1/poly(p))
3. Combine with empirical data (100% success rate at p=43 for known working orbits)

**Precedent**: "Algorithmic LLL" (Moser-Tardos) provides constructive existence via resampling algorithm. Could we do similar for SA?

**Challenge**: Proving convergence of SA for this problem requires analyzing the constraint landscape, which is complex.

### 9.4 Literature Gaps
Based on web search, no existing work directly addresses:
- LLL with p ≈ 1/2 and dense dependency graphs
- Circulant graph constructions for Ramsey numbers at p ≡ 3 mod 4
- Negative dependency exploitation when dependencies are still dense
- Second moment bounds for highly correlated counting variables with global budget constraints

**Implication**: This problem may require NEW TECHNIQUES, not just application of existing probabilistic methods.

## 10. Summary

| Method | Condition | Status for This Problem | Verdict |
|--------|-----------|------------------------|---------|
| Standard LLL | e·p·(d+1) ≤ 1 | ~1.36·p >> 1 | **FAILS** |
| Lopsided LLL | Negative dependency | Still dense dependencies | **Unlikely** |
| Cluster Expansion LLL | Clique structure | Dependency graph nearly complete | **Unlikely** |
| Entropy Compression | Same as LLL | e·p·(d+1) >> 1 | **FAILS** |
| Fix D11, LLL on D12 | Reduce dependencies | p ≈ 1/2 still too high | **FAILS** |
| Second Moment Method | E[N²]/E[N]² = O(poly) | Empirical: true up to p=31; general: OPEN | **L6 GAP** |
| Rödl Nibble | Incremental construction | Global constraints, unclear if applicable | **Unknown** |
| Algebraic Construction | Finite field structure | No construction known for p ≡ 3 mod 4 | **Not Available** |

**Bottom Line**:
- LLL and standard probabilistic methods **definitively fail** due to p ≈ 1/2 + dense dependencies.
- Second moment method is the **most promising** approach but requires proving E[N²]/E[N]² = O(poly(p)).
- **This is exactly the L6 gap** — no existing technique resolves it.
- May require **new probabilistic/combinatorial methods** that exploit:
  - Spectral complementarity
  - Negative budget correlations
  - Overlap structure of valid pairs
  - Circulant symmetry

## 11. Sources

### Lovász Local Lemma - General
- [Constructive proof of the general Lovász local lemma | Journal of the ACM](https://dl.acm.org/doi/10.1145/1667053.1667060)
- [Probabilistic Methods in Combinatorics — Yufei Zhao](https://yufeizhao.com/pm/6.pdf)
- [Lovász local lemma - Wikipedia](https://en.wikipedia.org/wiki/Lov%C3%A1sz_local_lemma)
- [Lecture 2. The Lovász Local Lemma - Stanford CS Theory](https://theory.stanford.edu/~jvondrak/MATH233A-2018/Math233-lec02.pdf)

### Entropy Compression / Moser-Tardos
- [Moser's entropy compression argument | What's new (Terence Tao)](https://terrytao.wordpress.com/2009/08/05/mosers-entropy-compression-argument/)
- [Moser-Tardos resampling algorithm, entropy compression method and the subset gas](https://arxiv.org/abs/2001.00880)
- [Entropy compression - Wikipedia](https://en.wikipedia.org/wiki/Entropy_compression)

### Lopsided LLL
- [An Algorithmic Proof of the Lopsided Lovász Local Lemma](https://www.cs.ubc.ca/~nickhar/W15/Lecture22Notes.pdf)
- [MIT Lectures 10–14: Lovász Local Lemma](https://ocw.mit.edu/courses/18-226-probabilistic-methods-in-combinatorics-fall-2022/mit18_226_f22_lec10-14.pdf)

### Cluster Expansion
- [An Improvement of the Lovász Local Lemma via Cluster Expansion](https://arxiv.org/abs/0910.1824)
- [Lecture 6. The Cluster Expansion Lemma - Stanford CS Theory](https://theory.stanford.edu/~jvondrak/MATH233A-2018/Math233-lec06.pdf)
- [The Lovász Local Lemma: constructive aspects, stronger variants](https://theory.stanford.edu/~jvondrak/data/LLL-combsem.pdf)

### Shearer's Condition
- [Directed Lovász Local Lemma and Shearer's Lemma](https://arxiv.org/abs/1611.00502)
- [Variable Version Lovász Local Lemma: Beyond Shearer's Bound](https://arxiv.org/abs/1709.05143)

### Semi-Random Method / Rödl Nibble
- [Refined Absorption: A New Proof of the Existence Conjecture](https://arxiv.org/html/2510.19978)
- [A Short Proof of the Random Ramsey Theorem](https://www.cambridge.org/core/journals/combinatorics-probability-and-computing/article/abs/short-proof-of-the-random-ramsey-theorem/1F00D51A100FEEC7CD7814EC2F76AA83)

### Algebraic Methods
- [New lower bounds for Ramsey numbers of graphs and hypergraphs](http://homepages.math.uic.edu/~mubayi/papers/hypramsub.pdf)
- [Algebraic Methods for Ramsey Lower Bounds](https://arxiv.org/pdf/math/0608013)
- [Mini-course in Finite Geometry and Ramsey Theory](https://anuragbishnoi.wordpress.com/minicourse/)

### Second Moment Method
- [Second moment method - Wikipedia](https://en.wikipedia.org/wiki/Second_moment_method)
- [Lecture notes (MIT 18.226) Probabilistic Methods in Combinatorics](https://yufeizhao.com/pm/probmethod_notes.pdf)
- [First and Second Moment Methods](https://cse.buffalo.edu/~hungngo/classes/2011/Spring-694/lectures/sm.pdf)

### Paley Construction
- [Paley graph - Wikipedia](https://en.wikipedia.org/wiki/Paley_graph)
- [Paley construction - Wikipedia](https://en.wikipedia.org/wiki/Paley_construction)
- [Transitive Subtournaments of k-th Power Paley Digraphs](https://link.springer.com/article/10.1007/s00373-024-02792-7)

### Recent Ramsey Bounds
- [An exponential improvement for Ramsey lower bounds (2025)](http://staff.ustc.edu.cn/~jiema/ramsey-lower-bound_arXiv.pdf)
- [Some recent results in Ramsey theory (2026)](https://arxiv.org/html/2601.05221)
