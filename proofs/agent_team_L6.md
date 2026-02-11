# Agent Team: Solve L6 for Book Ramsey Proof

## Context (READ FIRST)

We are proving R(B_{n-1}, B_n) = 4n-1. The proof is complete EXCEPT for one lemma (L6). The previous proof of L6 via Slepian's inequality was WRONG — Slepian goes the wrong direction for negative correlations and "all below threshold" events. See L6_recovery_tasks.md for the full error analysis.

### What we need to prove
For fixed symmetric D11 ⊂ {1,...,m-1} with |D11| = n = (m+1)/2, and D12 a uniformly random ((m-1)/2)-subset of Z_m containing 0:

  Pr[∀d ∈ D11: B(d) ≤ T(d)] ≥ 2^{-(m-1)}

where B(d) = Δ(D12, D12; d) and T(d) = (m-3)/2 - A(d).

This is a LOW bar. E[valid D12] ≈ 2^{0.7m}, so we have ~0.7m bits of exponential headroom. We just need the joint probability to not collapse worse than 2^{-m}.

### Key facts already proven
- E[B(d)] = (m-3)/4
- Var[B(d)] ≈ m/16
- All non-complementary Cov[B(d1), B(d2)] = -(m+1)/(8(m-2)), i.e., ρ = -2/(m-3)
- Σ B(d) = |D12|(|D12|-1) EXACTLY (hard constraint, not approximate)
- Each per-constraint probability Pr[B(d) ≤ T(d)] ≈ 1/2
- Computational evidence: joint/product ratio ranges from 4.7× (p=11) to 131× (p=59), GROWING with p

### Why the old approach failed
Pairwise correlation is NEGATIVE (ρ < 0). For Gaussians with negative correlation, Pr[all ≤ t] < ∏Pr[≤ t]. But the ACTUAL discrete distribution shows massive positive association. The positive association is a higher-order effect driven by the hard sum constraint + asymmetric threshold structure, invisible to any Gaussian/two-moment approximation.

---

## Agent 1: Second Moment Method (COMPUTATIONAL — START IMMEDIATELY)

### Goal
Compute E[N²]/E[N]² for the count N = #{valid D12} for fixed D11. If this ratio is polynomial in p, Paley-Zygmund gives Pr[N>0] → 1 and L6 follows.

### Tasks

**Task 1a: Exact computation for p=11**
```python
# For p=11, m=11, k=4, N=10, n=6
# D12 is a random 4-subset of {1,...,10} plus {0}
# For a known-good D11 (e.g., from solutions_registry.json or SA):
# 1. Enumerate ALL C(10,4) = 210 possible D12 sets
# 2. For each, compute A(d)+B(d) for all d in D11
# 3. Mark which D12 are valid (all constraints satisfied)
# 4. N = count of valid D12
# 5. E[N] = N (uniform distribution, each has prob 1/210, so E[N] = N)
#    Actually: E[N] = (# valid) since we're counting over uniform draw
# 6. For E[N²], count ordered pairs of valid D12:
#    E[N²] = #{(D12, D12') : both valid} where both drawn uniformly
#    = (# valid)² when draws are independent
#    BUT N² for a single draw is just N² = (# valid in the single draw)²
#
# CORRECTION: N is a random variable = 1[D12 is valid].
# No — N is the NUMBER of valid D12 sets. It's a fixed integer, not random.
# 
# Let me reframe. For fixed D11, define:
#   X(D12) = 1[D12 satisfies all constraints]
#   E_uniform[X] = #{valid D12} / C(10,4)
#
# What we actually want is just: #{valid D12} ≥ 1.
# So just COUNT valid D12 for each D11.

# ACTUAL USEFUL COMPUTATION:
# For EACH symmetric D11 of size 6 in Z_11:
#   Count how many of the 210 D12 sets are valid.
#   Report the minimum count over all D11.
# If minimum count ≥ 1, the theorem holds for p=11.
```

**Task 1b: Same for p=19, p=23**
For p=19: C(18,8) = 43758 D12 sets. Enumerate all.
For p=23: C(22,10) = 646646. Enumerate all. (May take minutes, feasible.)

**Task 1c: If direct enumeration proves existence for p ≤ 23, do sampling for p=31, 43, 47, 59**
For each prime, fix 50+ random symmetric D11 of size n. For each:
- Sample 10M random D12 sets
- Count how many are valid
- Report estimated E[valid D12]
- If even ONE valid D12 is found, record it explicitly

**Deliverable**: A table showing #{valid D12} for each tested (p, D11) pair. If we find valid D12 for ALL D11 at each p ≤ 23 by exhaustive enumeration, that's a partial result worth recording.

---

## Agent 2: Conditioning on Partial Sums

### Goal
Prove L6 by conditioning on S₁ = Σ_{d ∈ D11} B(d), then bounding the conditional max.

### Tasks

**Task 2a: Compute Var[S₁] exactly**

S₁ = Σ_{d ∈ D11} B(d). But B(d) = B(m-d) identically, and D11 is symmetric (d ∈ D11 ⟹ m-d ∈ D11). So if D11 = {d₁, m-d₁, d₂, m-d₂, ...}, then S₁ = 2 Σ_{i=1}^{r} B(dᵢ) where r = n/2 (approximately; handle the case d = m-d carefully — this happens when 2d ≡ 0 mod m, i.e. d = m/2, which only exists when m is even, so never for odd m).

So S₁ = 2 S₁', where S₁' = Σ_{i=1}^{r} B(dᵢ) over r = (p-1)/4 independent representatives.

Var[S₁'] = r × Var[B] + r(r-1) × Cov_non-comp
         = r × Var[B] + r(r-1) × (-2Var/(m-3))

Compute this exactly. Verify numerically for p=11, 19.

**Task 2b: Distribution of S₁ — numerical study**

For p = 11, 19, 23: enumerate all D12 and compute the distribution of S₁. Plot it. Is it approximately normal? Compute its mean, variance, skewness, kurtosis.

**Task 2c: Conditional max deviation**

Conditioned on S₁ = s, what is the distribution of max_{d ∈ D11} B(d)?

For p = 11, 19: for each value of s, compute:
- The conditional distribution of (B(d₁), ..., B(dᵣ)) given S₁' = s/2
- The conditional max
- Pr[max > T | S₁' = s/2]

Plot Pr[all ok | S₁ = s] vs s. Identify the critical value of s where this transitions from ~1 to ~0.

**Task 2d: Prove the bound**

The goal is to show:
  Pr[all ok] = Σ_s Pr[S₁=s] × Pr[all ok | S₁=s] ≥ 2^{-(m-1)}

Strategy: split into s ≤ s* (favorable, Pr[all ok|s] is not too small) and s > s* (unfavorable, but Pr[S₁>s*] might be small enough). The hard sum constraint means S₁ + S₂ = const, so controlling S₁ controls S₂ too.

### Deliverable
Exact formula for Var[S₁], numerical distributions, and a proof sketch showing the conditioning argument works (or a clear explanation of why it doesn't).

---

## Agent 3: Strong Rayleigh / Negative Dependence Literature

### Goal
Search the literature for theorems about quadratic functions of strongly Rayleigh (SR) measures. The k-subset measure is SR. We need positive association of THRESHOLD events for QUADRATIC functions of SR indicators.

### Tasks

**Task 3a: Literature search**

Search for and read these papers:
1. Borcea-Brändén-Liggett 2009 "Negative dependence and geometry of polynomials"
   - What does SR imply for quadratic functions? Linear functions of SR indicators are NA (negatively associated). But B(d) is quadratic.

2. Pemantle 2000 "Towards a theory of negative dependence"
   - Survey of NA, CNA, ULC properties. Any results for quadratic forms?

3. Anari-Gharan-Vinzant 2021 "Log-concave polynomials IV"
   - Extensions of SR theory. Anything about polynomial functions?

4. Wagner 2008 "Negatively correlated random variables and Mason's conjecture"

5. Any paper on "positive association of threshold events for negatively dependent random variables"

**Task 3b: Check if B(d) ≤ T(d) events are monotone**

An event is INCREASING if adding an element to S can only make it more likely. An event is DECREASING if adding an element can only make it less likely.

Is the event {B(d) ≤ T(d)} monotone in either direction?

Adding element x to D12:
- B(d) changes by: +1[x-d ∈ D12] + 1[x+d ∈ D12] (roughly — work out the exact formula)
- So B(d) can increase or decrease depending on which other elements are present.

The event {B(d) ≤ T(d)} is therefore NOT monotone. This means FKG-type inequalities (which require monotone events) don't apply directly.

BUT: is there a monotone coupling? Can we write {B(d) ≤ T(d)} as an intersection or union of monotone events?

**Task 3c: Check the "stochastic covering" property**

For SR measures, Feder-Mihail 1992 showed the "stochastic covering" property: conditioning on an element being in S makes other elements LESS likely to be in S (negative correlation). This is for indicators, not for quadratic functions. But it might help build a coupling argument.

### Deliverable
A summary of what the SR/NA literature says about our specific setting. Identify any theorem that could apply, or prove that none of the standard theorems apply and explain why.

---

## Agent 4: Energy Redistribution Mechanism

### Goal
Quantify the energy redistribution mechanism precisely enough to prove L6 directly.

### Tasks

**Task 4a: Partition function approach**

Define Z(D11) = #{D12 : all D11 constraints satisfied}. We want Z ≥ 1.

The total energy is fixed: Σ B(d) = C. The D11 constraints require Σ_{d ∈ D11} B(d) ≤ n(n-2). The D22 threshold is n+1 per position. Think of this as a constrained partition function.

For p = 11, 19: compute Z exactly. Then compute Z with various relaxations:
- Z_relaxed = #{D12 : S₁ ≤ n(n-2)} (only the sum constraint, not individual positions)
- Compare Z / Z_relaxed = Pr[all individual ok | sum ok]

**Task 4b: Entropy of valid D12 sets**

For p = 11, 19: enumerate all valid D12 for a fixed D11. Compute:
- How many valid D12 exist
- The entropy of the uniform distribution on valid D12 sets
- How this compares to log₂ C(m-1, k) (the entropy of uniformly random D12)
- The "information cost" = log₂ C(m-1, k) - log₂ Z

This gives a direct measurement of how many bits the constraints cost. If the information cost is < m-1 bits (i.e., Z > 1), we're done for that p.

**Task 4c: Pairwise overlap of valid D12 sets**

For p = 11, 19: for all pairs of valid D12 sets, compute |D12 ∩ D12'|. What's the distribution? Are valid D12 sets clustered (high overlap) or spread out? This informs whether the second moment method (Agent 1) will work.

### Deliverable
Tables of Z(D11) for all D11 at p=11,19. Overlap distribution. Information cost analysis.

---

## Agent 5: Composite m Extension

### Goal  
Verify that the entire proof framework (not just equi-covariance) extends to composite m.

### Tasks

**Task 5a: Verify equi-covariance for composite m**

The proof in proof_equi_covariance.md never uses primality. Verify this claim by:
1. Reading the proof and confirming no step requires m prime
2. Running the computational verification for composite m = 15, 21, 25, 27, 33, 35
   - For each: pick random D11 of size (m+1)/2
   - Compute all pairwise Cov[B(d1), B(d2)] by MC sampling
   - Check they match -(m+1)/(8(m-2)) for non-complementary pairs

**Task 5b: Check structural reduction for composite m**

Theorems 1-3 (structural reduction) were stated for prime p. Re-derive for general odd m:
- Does the 2-block circulant construction still work on Z_m when m is composite?
- Are the adjacency constraints the same?
- Does the complement structure D22 = {1,...,m-1} \ D11 still give the right threshold?

Key concern: the construction builds a graph on 2m vertices with vertex sets V1, V2 each of size m, with edges determined by circulant structure on Z_m. This works for any m, not just prime.

**Task 5c: Optimal |D11| for composite m**

For composite m = 15, 21, 25, 27, 33, 35:
- Run the D11 size survey (as in d11_size_survey.py) to find optimal |D11|
- Is |D11| = (m+1)/2 still optimal?
- If not, what is, and does the first moment argument still close?

### Deliverable
Verification tables for composite m. Clear statement of which parts of the proof carry over unchanged and which need modification.

---

## Coordination Notes

- Agents 1 and 4 need the same enumerations (all D12 for fixed D11 at small p). Share this computation.
- Agent 2 depends on Agent 1's enumeration data for the conditional distribution study.
- Agent 3 is pure literature review and can run fully independently.
- Agent 5 is independent of the L6 question and can run in parallel.
- If Agent 1 finds that #{valid D12} ≥ 1 for ALL D11 at p ≤ 23 by exhaustive search, this is already a publishable partial result (extends computational verification from n ≤ 21 to n ≤ 12 by proof).
- The MOST IMPORTANT single computation is Agent 1, Task 1a. If we can exhaustively verify existence for small p, and Agent 2 or 3 provides an asymptotic argument for large p, the proof is complete.

## Success Criteria

The proof is complete if ANY of these are achieved:

1. Direct enumeration shows ≥1 valid D12 exists for every symmetric D11 of size n, for all p ≤ p₀, PLUS an asymptotic argument for p > p₀.

2. Paley-Zygmund via E[N²]/E[N]² = O(poly(p)).

3. Conditioning argument: Pr[all ok | S₁ ≤ s*] × Pr[S₁ ≤ s*] ≥ 2^{-(m-1)}.

4. Any other valid lower bound on Pr[all D11 ok] that exceeds 2^{-(m-1)}.
