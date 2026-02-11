"""
Fourier LP feasibility analysis for the constraint reduction problem.

Given D11, determine whether a SPECTRUM P(k) exists that satisfies all constraints.
This is a LINEAR PROGRAM in the Fourier domain.

FORMULATION:
  Variables: P(k) for k = 1, ..., (p-1)/2  (power spectrum, symmetric: P(k) = P(p-k))

  Constraints:
    (i)   P(k) >= |D̂11(k)|²        for all k   [non-negativity of |D̂12(k)|²]
    (ii)  2·Σ P(k) = S             where S = p(|D11|+|D12|) - P(0)  [Parseval]
    (iii) P(0)/p + (2/p) Σ P(k) cos(2πkd/p) ≤ n-2  for d ∈ D11   [binding]
    (iv)  P(0)/p + (2/p) Σ P(k) cos(2πkd/p) ≤ n+1  for d ∈ D22   [loose]

If feasible: a target spectrum exists → existence proof may be possible.
If infeasible: no D12 can work for this D11 → need a different D11 or approach.

Also computes the MARGIN: how much slack exists in the LP optimal solution.
"""

import sys
import os
import json
import time
import math
import numpy as np
from scipy.optimize import linprog

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ramsey_core import BlockCirculantGraph, verify_construction, Delta


def autocorrelation_fft(indicator, m):
    fft_val = np.fft.fft(indicator)
    autocorr = np.fft.ifft(np.abs(fft_val) ** 2).real
    return np.round(autocorr).astype(int)


KNOWN = {
    11: {"D11": {1, 2, 4, 7, 9, 10}, "D12": {0, 4, 5, 7, 10}},
    19: {"D11": {1, 2, 3, 6, 8, 11, 13, 16, 17, 18},
         "D12": {0, 2, 6, 8, 9, 12, 13, 17, 18}},
    23: {"D11": {5, 6, 7, 8, 9, 11, 12, 14, 15, 16, 17, 18},
         "D12": {0, 1, 2, 6, 10, 13, 14, 16, 18, 20, 21}},
    31: {"D11": {6, 7, 8, 10, 11, 12, 14, 15, 16, 17, 19, 20, 21, 23, 24, 25},
         "D12": {0, 1, 2, 3, 8, 11, 12, 13, 15, 18, 20, 21, 24, 27, 29}},
    43: {"D11": {1, 2, 5, 10, 11, 13, 16, 17, 18, 19, 20, 23, 24, 25, 26, 27,
                 30, 32, 33, 38, 41, 42},
         "D12": {0, 2, 5, 6, 8, 11, 15, 16, 20, 24, 25, 27, 28, 31, 32, 34,
                 35, 36, 37, 39, 41}},
}


def fourier_lp(p, D11, minimize_max=True):
    """
    Solve the Fourier LP for given D11.

    Variables: P(k) for k = 1, ..., F  where F = (p-1)/2.

    Returns feasibility status and optimal margin.
    """
    n = (p + 1) // 2
    F = (p - 1) // 2  # number of independent Fourier modes
    D11_set = set(D11)
    D22 = set(range(1, p)) - D11_set

    d11_size = len(D11)
    d12_size = (p - 1) // 2

    # Compute |D̂11(k)|² for k = 1, ..., F
    ind11 = np.zeros(p, dtype=np.float64)
    for j in D11:
        ind11[j] = 1.0
    dft11 = np.fft.fft(ind11)  # D̂11(k) for k=0,...,p-1
    power11 = np.abs(dft11) ** 2  # |D̂11(k)|²

    # P(0) = |D11|² + |D12|²
    P0 = d11_size ** 2 + d12_size ** 2

    # Parseval: 2·Σ_{k=1}^F P(k) = p(|D11|+|D12|) - P(0)
    S = p * (d11_size + d12_size) - P0

    # Cosine matrix: cos(2πkd/p) for k=1..F, d=1..p-1
    cos_matrix = np.zeros((p - 1, F))
    for k in range(1, F + 1):
        for d_idx, d in enumerate(range(1, p)):
            cos_matrix[d_idx, k - 1] = math.cos(2 * math.pi * k * d / p)

    # Objective: minimize max(A(d)+B(d) at D11) or just check feasibility
    # A(d)+B(d) = P(0)/p + (2/p) Σ_{k=1}^F P(k) cos(2πkd/p)

    # For the LP, we minimize t subject to:
    # (2/p) Σ P(k) cos(2πkd/p) ≤ t - P(0)/p  for d ∈ D11 (with t ≤ n-2)
    # (2/p) Σ P(k) cos(2πkd/p) ≤ (n+1) - P(0)/p  for d ∈ D22
    # 2·Σ P(k) = S
    # P(k) ≥ power11(k)  for k = 1..F

    # Variables: [P(1), P(2), ..., P(F), t]  (F+1 variables)
    # We want to minimize t.

    c = np.zeros(F + 1)
    c[F] = 1.0  # minimize t

    # Inequality constraints: A_ub @ x <= b_ub
    A_ub_list = []
    b_ub_list = []

    # For d ∈ D11: (2/p) Σ P(k) cos(2πkd/p) - t ≤ -P(0)/p + (n-2)
    # Wait, we want: P(0)/p + (2/p) Σ P(k) cos(2πkd/p) ≤ t
    # So: (2/p) Σ P(k) cos(2πkd/p) - t ≤ -P(0)/p
    for d in sorted(D11):
        row = np.zeros(F + 1)
        for k in range(1, F + 1):
            row[k - 1] = (2.0 / p) * math.cos(2 * math.pi * k * d / p)
        row[F] = -1.0  # coefficient of t
        A_ub_list.append(row)
        b_ub_list.append(-P0 / p)

    # Also need t ≤ n-2 (binding threshold)
    row_t = np.zeros(F + 1)
    row_t[F] = 1.0
    A_ub_list.append(row_t)
    b_ub_list.append(n - 2)

    # For d ∈ D22: P(0)/p + (2/p) Σ P(k) cos(2πkd/p) ≤ n+1
    for d in sorted(D22):
        row = np.zeros(F + 1)
        for k in range(1, F + 1):
            row[k - 1] = (2.0 / p) * math.cos(2 * math.pi * k * d / p)
        A_ub_list.append(row)
        b_ub_list.append(n + 1 - P0 / p)

    A_ub = np.array(A_ub_list)
    b_ub = np.array(b_ub_list)

    # Equality constraint: 2·Σ P(k) = S
    A_eq = np.zeros((1, F + 1))
    for k in range(F):
        A_eq[0, k] = 2.0
    b_eq = np.array([S])

    # Bounds: P(k) ≥ power11(k) for k=1..F, t ≥ -∞
    bounds = []
    for k in range(1, F + 1):
        bounds.append((power11[k], None))
    bounds.append((None, None))  # t is unbounded

    # Solve
    result = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                     bounds=bounds, method='highs')

    if result.success:
        optimal_t = result.x[F]
        P_optimal = result.x[:F]

        # Compute A(d)+B(d) at all positions with optimal P
        ab_values = {}
        for d in range(1, p):
            ab = P0 / p
            for k in range(1, F + 1):
                ab += (2.0 / p) * P_optimal[k - 1] * math.cos(2 * math.pi * k * d / p)
            ab_values[d] = ab

        # Check margins
        d11_ab = [ab_values[d] for d in sorted(D11)]
        d22_ab = [ab_values[d] for d in sorted(D22)]

        d11_margin = min(n - 2 - ab for ab in d11_ab)
        d22_margin = min(n + 1 - ab for ab in d22_ab)

        # Compute optimal |D̂12(k)|² = P(k) - |D̂11(k)|²
        d12_power = [P_optimal[k] - power11[k + 1] for k in range(F)]

        return {
            "feasible": True,
            "optimal_t": optimal_t,
            "d11_margin": d11_margin,
            "d22_margin": d22_margin,
            "max_ab_d11": max(d11_ab),
            "max_ab_d22": max(d22_ab),
            "P_optimal": P_optimal.tolist(),
            "d12_power_min": min(d12_power),
            "d12_power_max": max(d12_power),
        }
    else:
        return {
            "feasible": False,
            "status": result.message,
        }


def check_known_solution(p, D11, D12):
    """Verify the known solution satisfies the LP constraints (sanity check)."""
    n = (p + 1) // 2
    F = (p - 1) // 2
    d11_size = len(D11)
    d12_size = len(D12)

    ind11 = np.zeros(p, dtype=np.float64)
    for j in D11:
        ind11[j] = 1.0
    ind12 = np.zeros(p, dtype=np.float64)
    for j in D12:
        ind12[j] = 1.0

    dft11 = np.fft.fft(ind11)
    dft12 = np.fft.fft(ind12)
    P = np.abs(dft11) ** 2 + np.abs(dft12) ** 2

    A = autocorrelation_fft(ind11, p)
    B = autocorrelation_fft(ind12, p)

    # Check A+B at all positions
    max_d11 = max(int(A[d]) + int(B[d]) for d in D11)
    max_d22 = max(int(A[d]) + int(B[d]) for d in set(range(1, p)) - set(D11))

    return {
        "max_ab_d11": max_d11,
        "max_ab_d22": max_d22,
        "P_values": P[:F + 1].tolist(),
        "P_range": (float(min(P[1:])), float(max(P[1:]))),
    }


def random_d11_lp(p, num_trials=100):
    """Test LP feasibility for random D11."""
    F = (p - 1) // 2
    pairs = []
    seen = set()
    for d in range(1, p):
        if d not in seen:
            pairs.append((d, (-d) % p))
            seen.add(d)
            seen.add((-d) % p)

    num_pairs_needed = (p + 1) // 4
    rng = np.random.default_rng(42)

    feasible_count = 0
    margins = []

    for trial in range(num_trials):
        chosen = rng.choice(len(pairs), size=num_pairs_needed, replace=False)
        D11 = set()
        for i in chosen:
            D11.add(pairs[i][0])
            D11.add(pairs[i][1])

        result = fourier_lp(p, D11)
        if result["feasible"]:
            feasible_count += 1
            margins.append(result["d11_margin"])

    return {
        "feasible_count": feasible_count,
        "total": num_trials,
        "feasible_rate": feasible_count / num_trials,
        "margins": margins,
    }


def main():
    primes = [11, 19, 23, 31, 43]

    print("=" * 90)
    print("FOURIER LP FEASIBILITY ANALYSIS")
    print("=" * 90)
    print("""
For each D11, solve the LP: minimize max(A(d)+B(d)) over d ∈ D11
subject to:
  - P(k) ≥ |D̂11(k)|²  (non-negativity of |D̂12|²)
  - 2·Σ P(k) = S       (Parseval)
  - A(d)+B(d) ≤ n+1     for d ∈ D22

If the optimal value ≤ n-2, then a valid SPECTRUM exists.
""")

    for p in primes:
        n = (p + 1) // 2
        D11 = KNOWN[p]["D11"]
        D12 = KNOWN[p]["D12"]
        D22 = set(range(1, p)) - D11

        print(f"\n{'='*80}")
        print(f"p = {p}, n = {n}, F = {(p-1)//2} Fourier modes")
        print(f"|D11| = {len(D11)}, |D12| = {len(D12)}")
        print(f"Tight threshold: {n-2}, Loose threshold: {n+1}")
        print(f"{'='*80}")

        # Sanity check: known solution
        known = check_known_solution(p, D11, D12)
        print(f"\n  Known solution check:")
        print(f"    max A+B at D11: {known['max_ab_d11']} (threshold {n-2})")
        print(f"    max A+B at D22: {known['max_ab_d22']} (threshold {n+1})")
        print(f"    P(k) range: {known['P_range'][0]:.2f} to {known['P_range'][1]:.2f}")

        # LP with known D11
        t0 = time.time()
        lp = fourier_lp(p, D11)
        elapsed = time.time() - t0

        print(f"\n  LP RESULT (known D11):")
        print(f"    Feasible: {lp['feasible']}")
        if lp["feasible"]:
            print(f"    Optimal max(A+B at D11): {lp['optimal_t']:.4f}")
            print(f"    D11 margin: {lp['d11_margin']:.4f}")
            print(f"    D22 margin: {lp['d22_margin']:.4f}")
            print(f"    |D̂12|² range: [{lp['d12_power_min']:.4f}, {lp['d12_power_max']:.4f}]")
        else:
            print(f"    Status: {lp.get('status', 'unknown')}")
        print(f"    Time: {elapsed:.3f}s")

        # Random D11 LP analysis
        num_random = min(200, max(50, 10000 // p))
        print(f"\n  Random D11 LP analysis ({num_random} trials):")
        t0 = time.time()
        rand_lp = random_d11_lp(p, num_trials=num_random)
        elapsed = time.time() - t0

        print(f"    LP-feasible D11: {rand_lp['feasible_count']}/{rand_lp['total']} "
              f"({rand_lp['feasible_rate']:.2%})")
        if rand_lp["margins"]:
            print(f"    D11 margins: min={min(rand_lp['margins']):.4f}, "
                  f"max={max(rand_lp['margins']):.4f}, "
                  f"mean={sum(rand_lp['margins'])/len(rand_lp['margins']):.4f}")
        print(f"    Time: {elapsed:.1f}s")

    # Test larger primes (without known D11)
    print(f"\n{'='*90}")
    print("LARGER PRIMES (random D11 only)")
    print(f"{'='*90}")

    larger_primes = [47, 59, 67, 71, 79, 83]
    for p in larger_primes:
        n = (p + 1) // 2
        num_random = min(100, max(20, 5000 // p))
        print(f"\n  p={p}, n={n}: Testing {num_random} random D11...")
        t0 = time.time()
        rand_lp = random_d11_lp(p, num_trials=num_random)
        elapsed = time.time() - t0

        print(f"    LP-feasible: {rand_lp['feasible_count']}/{rand_lp['total']} "
              f"({rand_lp['feasible_rate']:.2%})")
        if rand_lp["margins"]:
            print(f"    Margins: min={min(rand_lp['margins']):.4f}, "
                  f"max={max(rand_lp['margins']):.4f}, "
                  f"mean={sum(rand_lp['margins'])/len(rand_lp['margins']):.4f}")
        print(f"    Time: {elapsed:.1f}s")

    # Summary
    print(f"\n{'='*90}")
    print("KEY INSIGHT")
    print(f"{'='*90}")
    print("""
If the LP is feasible for most random D11:
  → The spectral constraints are NOT the bottleneck
  → The difficulty is realizing the target spectrum with an actual set D12
  → A probabilistic argument may work in the SPECTRAL domain

If the LP is feasible for ALL tested D11:
  → ANY symmetric D11 admits a valid spectrum
  → The proof reduces to showing spectral realizability
  → This is a strong form of the conjecture

If the LP is infeasible for most D11:
  → D11 must be specially chosen
  → The proof must specify which D11 works
""")


if __name__ == "__main__":
    main()
