#!/usr/bin/env python3
"""Complete survey of |D11| sizes for p ≡ 3 mod 4 primes.

For each prime p and each feasible even |D11| = k, exhaustively count
all valid (D11, D12) pairs satisfying the algebraic constraints:
  A(d) + B(d) ≤ n-2      for d ∈ D11
  A(d) + B(d) ≤ 2k-n+1   for d ∈ D22

This reveals the optimal |D11| size and whether it follows a pattern.
"""

import json
import numpy as np
from itertools import combinations
from collections import Counter, defaultdict
import time
import os
from math import comb


def quadratic_residues(p):
    return sorted(set((x * x) % p for x in range(1, p)))


def autocorrelation(S, m):
    indicator = np.zeros(m)
    for x in S:
        indicator[x % m] = 1
    F = np.fft.fft(indicator)
    power = np.abs(F) ** 2
    result = np.real(np.fft.ifft(power))
    return np.rint(result).astype(int)


def find_all_symmetric_d11(p, k):
    """Enumerate all symmetric D11 ⊂ {1,...,p-1} of size k."""
    pairs = [(x, p - x) for x in range(1, (p + 1) // 2)]
    num_pairs = k // 2
    result = []
    for combo in combinations(range(len(pairs)), num_pairs):
        D11 = []
        for i in combo:
            D11.append(pairs[i][0])
            D11.append(pairs[i][1])
        result.append(sorted(D11))
    return result


def precompute_all_b(p, d12_size):
    """Precompute B(d) for ALL D12 ⊂ Z_p of given size.

    Returns (B_matrix[N, p], d12_list) where N = C(p, d12_size).
    """
    B_list = []
    d12_list = []
    for D12 in combinations(range(p), d12_size):
        D12 = list(D12)
        B = autocorrelation(D12, p)
        B_list.append(B)
        d12_list.append(D12)
    return np.array(B_list, dtype=np.int32), d12_list


def survey_single_k(all_d11, B_matrix, p, n, k):
    """Count valid (D11, D12) pairs for a given |D11| = k.

    Returns dict with statistics.
    """
    d11_thresh = n - 2
    d22_thresh = 2 * k - n + 1

    if d22_thresh < 0:
        return {
            'k': k, 'num_d11': len(all_d11), 'num_d12': len(B_matrix),
            'total_valid': 0, 'd11_with_valid': 0,
            'feasible': False, 'd22_thresh': d22_thresh,
            'd11_thresh': d11_thresh,
            'best_d11': None, 'best_count': 0,
        }

    total_valid = 0
    d11_with_valid = 0
    best_d11 = None
    best_count = 0
    best_d12_indices = []

    for D11 in all_d11:
        D11_set = set(D11)
        A = autocorrelation(D11, p)

        # Build threshold array
        threshold = np.zeros(p, dtype=np.int32)
        for d in range(1, p):
            threshold[d] = d11_thresh if d in D11_set else d22_thresh

        # Vectorized check
        sum_ab = B_matrix[:, 1:] + A[np.newaxis, 1:]
        thresh_row = threshold[np.newaxis, 1:]
        valid = np.all(sum_ab <= thresh_row, axis=1)
        cnt = int(np.sum(valid))

        total_valid += cnt
        if cnt > 0:
            d11_with_valid += 1
        if cnt > best_count:
            best_count = cnt
            best_d11 = D11
            best_d12_indices = list(np.where(valid)[0])

    return {
        'k': k, 'num_d11': len(all_d11), 'num_d12': len(B_matrix),
        'total_valid': total_valid, 'd11_with_valid': d11_with_valid,
        'feasible': True, 'd22_thresh': d22_thresh,
        'd11_thresh': d11_thresh,
        'best_d11': best_d11, 'best_count': best_count,
        'best_d12_indices': best_d12_indices[:5],
    }


def run_survey(p):
    """Complete survey for a single prime."""
    n = (p + 1) // 2
    d12_size = n - 1
    num_pairs = (p - 1) // 2  # symmetric pairs available

    total_d12 = comb(p, d12_size)
    print(f"\n{'='*80}")
    print(f"SURVEY: p={p}, n={n}, |D12|={d12_size}, "
          f"C({p},{d12_size})={total_d12} D12 candidates")
    print(f"{'='*80}")

    # Precompute all B arrays
    t0 = time.time()
    B_matrix, d12_list = precompute_all_b(p, d12_size)
    precomp_time = time.time() - t0
    print(f"  Precomputed {len(d12_list)} B arrays in {precomp_time:.1f}s")

    # Survey each even k from 2 to p-1
    results = []
    qr_set = set(quadratic_residues(p))

    for k in range(2, p, 2):
        all_d11 = find_all_symmetric_d11(p, k)
        if not all_d11:
            continue

        t1 = time.time()
        stats = survey_single_k(all_d11, B_matrix, p, n, k)
        elapsed = time.time() - t1

        feasible_str = "" if stats['feasible'] else " [INFEASIBLE: D22 thresh < 0]"
        valid_str = f"{stats['total_valid']:>8d}" if stats['feasible'] else "       -"

        print(f"  k={k:3d}: {len(all_d11):6d} D11s, "
              f"D11_thr={stats['d11_thresh']:3d}, D22_thr={stats['d22_thresh']:3d}, "
              f"valid={valid_str}, "
              f"D11s_ok={stats['d11_with_valid']:5d}"
              f"{feasible_str} [{elapsed:.1f}s]")

        # Show best D11 and example D12
        if stats['total_valid'] > 0 and stats['best_d11']:
            d11 = stats['best_d11']
            # Analyze D12 structure
            example_d12_idx = stats['best_d12_indices'][0]
            example_d12 = d12_list[example_d12_idx]
            d12_set = set(example_d12)
            qr_overlap = len(d12_set & qr_set)
            has_0 = 0 in d12_set
            print(f"         Best D11: {d11} "
                  f"({stats['best_count']} D12s work)")
            print(f"         Example D12: {example_d12} "
                  f"(0∈D12:{has_0}, |∩QR|={qr_overlap}/{len(qr_set)})")

        results.append(stats)

    return results, d12_list


def sa_solution_analysis():
    """Analyze |D11| choices in SA-found solutions."""
    print("\n" + "=" * 80)
    print("SA SOLUTION |D11| ANALYSIS")
    print("=" * 80)

    script_dir = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(script_dir, 'solutions_registry.json')
    with open(path, 'r') as f:
        data = json.load(f)

    print(f"\n  {'n':>3s} {'m':>4s} {'m%4':>4s} {'prime':>6s} "
          f"{'|D11|':>5s} {'|D12|':>5s} {'|D22|':>5s} "
          f"{'k-n':>4s} {'k/(m-1)':>8s} {'solver':>18s}")
    print("  " + "-" * 75)

    for sol in data['solutions']:
        m = sol['m']
        n_val = sol['n']
        d11 = sol['D11']
        d12 = sol['D12']
        k = len(d11)
        d12_len = len(d12)
        d22_len = m - 1 - k

        from sympy import isprime as _ip
        try:
            prime = _ip(m)
        except ImportError:
            prime = all(m % i != 0 for i in range(2, int(m**0.5)+1))

        k_minus_n = k - n_val
        k_ratio = k / (m - 1)

        print(f"  {n_val:3d} {m:4d} {m%4:4d} {'yes' if prime else 'no':>6s} "
              f"{k:5d} {d12_len:5d} {d22_len:5d} "
              f"{k_minus_n:+4d} {k_ratio:8.3f} {sol['solver']:>18s}")

    print("\n  Key: k-n shows |D11|-n. Paley uses k=n-1, SA often uses k=n or k=n-2.")


def summarize_results(all_results):
    """Print summary across all primes."""
    print("\n" + "=" * 80)
    print("CROSS-PRIME SUMMARY")
    print("=" * 80)

    print(f"\n  Optimal |D11| by prime (maximizing total valid pairs):")
    print(f"  {'p':>4s} {'n':>4s} {'best k':>7s} {'k-n':>5s} "
          f"{'k/(m-1)':>8s} {'valid':>10s} {'D11s_ok':>8s}")
    print("  " + "-" * 55)

    for p, results in all_results.items():
        n = (p + 1) // 2
        best = max(results, key=lambda r: r['total_valid'])
        if best['total_valid'] > 0:
            k = best['k']
            print(f"  {p:4d} {n:4d} {k:7d} {k-n:+5d} "
                  f"{k/(p-1):8.3f} {best['total_valid']:10d} "
                  f"{best['d11_with_valid']:8d}")
        else:
            print(f"  {p:4d} {n:4d}    NONE     -        - "
                  f"         0        0")


def main():
    print("=" * 80)
    print("COMPLETE |D11| SIZE SURVEY")
    print("For primes p ≡ 3 mod 4: which |D11| maximizes valid constructions?")
    print("=" * 80)

    # SA analysis first (quick)
    sa_solution_analysis()

    # Full surveys
    all_results = {}
    for p in [7, 11, 19]:
        results, _ = run_survey(p)
        all_results[p] = results

    # p=23: check if feasible
    p = 23
    n = 12
    d12_size = 11
    total_d12 = comb(p, d12_size)
    print(f"\n  p=23: C(23,11)={total_d12} D12 candidates")

    if total_d12 <= 2000000:
        results, _ = run_survey(p)
        all_results[p] = results
    else:
        print(f"  Skipping full enumeration (too many D12s)")

    summarize_results(all_results)

    # Deeper analysis of the optimal k
    print("\n" + "=" * 80)
    print("ANALYSIS: What makes the optimal k work?")
    print("=" * 80)

    for p in sorted(all_results):
        n = (p + 1) // 2
        results = all_results[p]
        best = max(results, key=lambda r: r['total_valid'])
        if best['total_valid'] == 0:
            print(f"\n  p={p}: No valid pairs at any k")
            continue

        k = best['k']
        qr = quadratic_residues(p)

        print(f"\n  p={p} (n={n}): optimal k={k} (k-n={k-n:+d})")
        print(f"    D11 threshold: A+B ≤ {n-2}")
        print(f"    D22 threshold: A+B ≤ {2*k-n+1}")
        print(f"    Average A(d) for |D11|={k}: "
              f"{k*(k-1)/(p-1):.2f}")
        print(f"    Average B(d) for |D12|={n-1}: "
              f"{(n-1)*(n-2)/(p-1):.2f}")
        print(f"    Average A+B: "
              f"{(k*(k-1)+(n-1)*(n-2))/(p-1):.2f}")
        print(f"    D11 slack (thresh - avg): "
              f"{n-2 - (k*(k-1)+(n-1)*(n-2))/(p-1):.2f}")

        # Distribution over k values
        print(f"\n    Valid pairs by k:")
        for r in results:
            if r['total_valid'] > 0:
                bar = "#" * min(50, r['total_valid'])
                print(f"      k={r['k']:3d}: {r['total_valid']:8d} {bar}")


if __name__ == '__main__':
    main()
