#!/usr/bin/env python3
"""
Extend the hyperplane conditioning approach to composite m = 2n-1.

For composite m, the construction uses Z_m (cyclic group of order m).
The key question: does the bound give positive margin?

Differences from prime p:
- Marginals Pr[B(d) = j] may depend on d (not uniform)
- Multiplicative orbits have varying sizes
- But the hyperplane Σ B(d) = S still holds
"""

import numpy as np
from math import comb, log2, log, exp, gcd
from collections import Counter
import time


def is_prime(n):
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0 or n % 3 == 0: return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0: return False
        i += 6
    return True


def is_prime_power(n):
    """Check if n is a prime power q^k with q ≡ 1 mod 4."""
    for p in range(2, int(n**0.5) + 1):
        if n % p == 0:
            k = 0
            m = n
            while m % p == 0:
                m //= p
                k += 1
            if m == 1:
                return p % 4 == 1  # prime power with base ≡ 1 mod 4
            return False
    # n itself is prime
    return n % 4 == 1


def compute_margin_composite(m, delta=0, n_samples=500_000, seed=42):
    """
    Estimate the hyperplane conditioning margin for Z_m.

    For composite m, we can't compute the exact cycle PMF analytically.
    Instead, we estimate marginals and the convolution by Monte Carlo.
    """
    k = (m - 1) // 2
    n = (m + 1) // 2
    S = k * (k - 1) // 2
    rng = np.random.default_rng(seed)

    # D11: symmetric set, first (m+1)/4 pairs (if m odd)
    # For composite m, D11 can still be symmetric
    D11_set = set()
    n_pairs = (m + 1) // 4  # Need |D11| = (m+1)/2

    # Simple D11: consecutive pairs
    for d in range(1, n_pairs + 1):
        D11_set.add(d)
        D11_set.add(m - d)

    actual_d11_size = len(D11_set)
    target_d11_size = (m + 1) // 2

    # May need to adjust for even/odd
    if actual_d11_size < target_d11_size:
        # Add middle element if needed
        if m // 2 not in D11_set:
            D11_set.add(m // 2)

    # Compute A-profile for this D11
    D11_vec = np.zeros(m, dtype=int)
    for d in D11_set:
        D11_vec[d] = 1

    # Representatives: d = 1,...,k
    reps = list(range(1, k + 1))

    A_vals = {}
    for d in reps:
        A_vals[d] = sum(D11_vec[i] & D11_vec[(i + d) % m] for i in range(m))

    # Thresholds
    T = {}
    for d in reps:
        if d in D11_set:
            T[d] = (m - 3) // 2 - A_vals[d] + delta
        else:
            T[d] = (m + 3) // 2 - A_vals[d] - delta

    # Check feasibility
    infeasible = sum(1 for d in reps if T[d] < 0)
    if infeasible > 0:
        return None

    # Monte Carlo: sample D12's and estimate P[E]
    # Also compute marginal distributions for each d
    marginals = {d: Counter() for d in reps}
    n_valid = 0
    n_checked = 0

    for _ in range(n_samples):
        # Random k-subset of Z_m
        D12 = set(rng.choice(m, k, replace=False))

        all_ok = True
        for d in reps:
            b = sum(1 for x in D12 if (x + d) % m in D12)
            marginals[d][b] += 1
            if b > T[d]:
                all_ok = False
        if all_ok:
            n_valid += 1
        n_checked += 1

    # Estimate Pr[E]
    pr_E = n_valid / n_checked if n_checked > 0 else 0

    # Budget
    budget = sum(log2((m - i) / (i + 1)) for i in range(k))

    # Check marginal uniformity
    # For each d, compute Pr[B(d) ≤ T(d)]
    trunc_probs = {}
    for d in reps:
        total_marg = sum(marginals[d].values())
        prob = sum(marginals[d][j] for j in range(T[d] + 1)) / total_marg
        trunc_probs[d] = prob

    trunc_cost = -sum(log2(max(p, 1e-15)) for p in trunc_probs.values())

    # Check if marginals are approximately identical across d
    marg_variance = np.var([trunc_probs[d] for d in reps])

    # Estimate margin
    if pr_E > 0:
        log2_EN = budget + log2(pr_E)
    else:
        log2_EN = -1e10  # Lower bound from product measure

    # Product-measure estimate
    naive_margin = budget - trunc_cost

    return {
        'm': m, 'k': k, 'n': n,
        'budget': budget,
        'trunc_cost': trunc_cost,
        'naive_margin': naive_margin,
        'pr_E_sampled': pr_E,
        'log2_EN': log2_EN if pr_E > 0 else 'N/A',
        'n_valid': n_valid,
        'n_samples': n_checked,
        'marg_uniformity': marg_variance < 0.001,
        'marg_var': marg_variance,
        'min_trunc_prob': min(trunc_probs.values()),
        'max_trunc_prob': max(trunc_probs.values()),
        'infeasible': infeasible,
    }


def main():
    print("COVERAGE ANALYSIS: which n values are proven?")
    print("=" * 80)

    proven = {}
    gaps = []

    for n in range(2, 130):
        m = 2 * n - 1
        if is_prime(m) and m % 4 == 1:
            proven[n] = f"Paley (p={m} ≡ 1 mod 4)"
        elif is_prime(m) and m % 4 == 3 and m <= 227:
            proven[n] = f"HC proof (p={m} ≡ 3 mod 4)"
        elif is_prime_power(m):
            proven[n] = f"Paley (q={m} prime power ≡ 1 mod 4)"
        elif n <= 21:
            proven[n] = f"Verified (n ≤ 21)"
        else:
            gaps.append((n, m))

    print(f"\nProven: {len(proven)} values of n")
    print(f"Gaps: {len(gaps)} values of n")
    print()

    # List the gaps
    print("GAPS (n, m=2n-1):")
    for n, m in gaps[:40]:
        factors = []
        temp = m
        for p in range(2, int(m**0.5) + 1):
            while temp % p == 0:
                factors.append(p)
                temp //= p
        if temp > 1:
            factors.append(temp)

        factorization = " × ".join(str(f) for f in factors)
        print(f"  n={n:3d}, m={m:3d} = {factorization}")

    print(f"\n... total {len(gaps)} gaps in n=2..129")

    # Try hyperplane conditioning for small composite m
    print("\n\nHYPERPLANE CONDITIONING FOR COMPOSITE m")
    print("=" * 80)

    composite_m = [m for n, m in gaps if m <= 100]

    for m in composite_m[:8]:
        n = (m + 1) // 2
        print(f"\n--- m={m}, n={n} ---")
        t0 = time.time()
        result = compute_margin_composite(m, delta=0, n_samples=200_000)
        elapsed = time.time() - t0

        if result is None:
            print("  INFEASIBLE")
            continue

        print(f"  Budget: {result['budget']:.1f} bits")
        print(f"  Trunc cost: {result['trunc_cost']:.1f} bits")
        print(f"  Naive margin: {result['naive_margin']:.1f} bits")
        print(f"  Marginals uniform: {result['marg_uniformity']} "
              f"(var={result['marg_var']:.6f})")
        print(f"  Pr[E] sampled: {result['pr_E_sampled']:.6f} "
              f"({result['n_valid']}/{result['n_samples']})")
        if isinstance(result['log2_EN'], float):
            print(f"  log2(E[N]) estimate: {result['log2_EN']:.2f}")
        print(f"  [{elapsed:.1f}s]")


if __name__ == '__main__':
    main()
