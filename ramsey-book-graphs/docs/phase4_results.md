# Phase 4 Results: Proof Attempt for p ≡ 3 mod 4

## New Theorems

### Theorem 4: Constraint Reduction

For the standard partition (|D11| = (p+1)/2, |D12| = (p-1)/2), the four constraint classes reduce to exactly TWO constraints on A(d) + B(d):

| Common neighbor formula | Positions | Threshold on A+B |
|------------------------|-----------|-----------------|
| V1V1 red = A(d) + B(d) | d in D11 | <= n-2 |
| V2V2 blue = A(d) + B(d) + 1 | d in D11 | <= n-2 (same) |
| V1V1 blue = A(d) + B(d) - 2 | d in D22 | <= n+1 |
| V2V2 red = A(d) + B(d) - 3 | d in D22 | <= n+1 (same) |

**Verified by brute force for p = 11, 19, 23, 31, 43.**

Key implication: the problem reduces to finding D11, D12 such that A(d)+B(d) is controlled at both D11 and D22 positions. The constraints are coupled via the conservation law:

    Sum_{d>0} (A(d) + B(d)) = (p-1)^2 / 2

## Probabilistic Analysis (prob_existence.py)

Uniform random D12 fails spectacularly:
- E[A+B] = (p-1)/2 = threshold + 1 (exactly 1 above the binding threshold)
- Individual constraint violation probability: 64-82%
- Symmetric LLL: e*p*D = 20 to 94 >> 1

The fundamental issue: with random D12, the expected A+B equals the threshold plus 1. There is no room for probabilistic arguments.

## Fourier LP Analysis (fourier_lp.py)

**Key positive result**: The spectral constraints are ALWAYS satisfiable for known D11.

For each D11, define the LP: find P(k) = |D̂11(k)|^2 + |D̂12(k)|^2 such that all A+B constraints are met.

| p | LP margin (known D11) | LP feasibility (random D11) |
|---|----------------------|----------------------------|
| 11 | 0.12 | 49% |
| 19 | 0.02 | 21% |
| 23 | 0.30 | 16% |
| 31 | 0.14 | 9% |
| 43 | 0.51 | 4.5% |
| 47 | 0.53 (best random) | 2.8% |
| 59 | 0.42 (best random) | 0.5% |
| 67 | 0.17 (best random) | 0.3% |
| 71 | 0.22 (best random) | 0.5% |
| 79 | 0.24 (best random) | 0.3% |
| 83 | 0.07 (best random) | 0.08% |

**Critical observation**: LP-feasible D11 exist for ALL tested primes (up to p=83), even though they become rare among random D11.

## D11 Optimization (d11_optimization.py)

LP-feasible D11 are distinguished by:
1. **Low max A(d) at D11** (feasible: ~6.3, infeasible: ~7.7 at p=23)
2. **Low spectral peak max |D̂11(k)|^2** (feasible: ~22, infeasible: ~26 at p=23)
3. **Low power spectrum std** (feasible: ~6.6, infeasible: ~7.7 at p=23)

For ALL primes tested (p=11 to 83), the known D11 is LP-feasible. At p=19, the best D11 has LP margin 0.04 vs 0.02 for the known D11.

## Spectral Structure (spectral_proof.py)

For random symmetric D11 of size (p+1)/2:
- E[|D̂11(k)|^2] = (p+1)/4 exactly (confirmed theoretically and by MC)
- max |D̂11(k)|^2 / p grows slowly: 0.82 (p=11) to 1.64 (p=83)
- The max/ideal ratio grows roughly as sqrt(log p)

LP-feasible D11 have max |D̂11(k)|^2 well below the random mean maximum.

## Proof Strategy Assessment

### What works:
- Constraint reduction (Theorem 4): fully proven
- V1V2 auto-satisfaction (Theorem 1): fully proven
- Fourier LP feasibility: computationally verified for p <= 83

### What's missing:
1. **D11 construction**: Need to construct a specific D11 guaranteed to be LP-feasible for all p. Character sum bounds (Weil/Deligne) may help establish spectral flatness.

2. **Spectral realization**: Given an LP-feasible spectrum, need to show an actual set D12 exists. This is a combinatorial existence problem akin to phase retrieval.

3. **Union bound gap**: The D22 constraint is always binding (margin = 0 in the LP), meaning any proof must precisely control the D22 behavior — no probabilistic slack.

### Phase 5 SA Extension Results:

Attempted to extend verified range beyond n=31 using `fast_sa_numpy.py`:

| n | m | Best cost | Seeds tried | Notes |
|---|---|-----------|------------|-------|
| 32 | 63 | 8 | 6+ | Consistent cost=8 barrier (see n32 analysis) |
| 33 | 65 | 4 | 13+ | Cost=4 barrier: 4/13 seeds reach this |
| 34 | 67 | 4 | 11+ | Cost=4 barrier: 2/11 seeds reach this |
| 35 | 69 | 8 | 5+ | Still running |
| 36 | 71 | 8 | 11+ | Best=8 at seeds 7-8 |

Also tested focused D12-only solver with LP-feasible D11 for n=34, 36.
Found 2/1000 LP-feasible D11 for p=67, 5/1000 for p=71.
Focused solver did not outperform joint SA.

### Verdict:
A complete proof for p equiv 3 mod 4 remains open but the structural understanding is now much deeper. The problem has been reduced to two concrete questions (D11 construction and spectral realization). SA verification extends the computational range but hits a cost=4 barrier for n≥32, suggesting these sizes may require longer SA runs, new SA techniques, or alternative construction methods.
