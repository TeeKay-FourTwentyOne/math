#!/usr/bin/env python3
"""Probabilistic estimation of valid (D11, D12) pair density.

Key finding from exhaustive survey: optimal |D11| = n-2 universally.
This script estimates whether valid pairs continue to exist at larger p
using Monte Carlo sampling of random D12 sets.

For each prime, we estimate:
  1. Per-D12 success rate for known/best D11
  2. Per-D11 success rate (fraction of random D11s admitting valid D12)
  3. Growth trend across primes
"""

import json
import numpy as np
from itertools import combinations
from math import comb
import time
import os


def quadratic_residues(p):
    return sorted(set((x * x) % p for x in range(1, p)))


def autocorrelation_batch(indicators, p):
    """Batch FFT autocorrelation for M indicator vectors.

    indicators: (M, p) float array
    Returns: (M, p) int array of autocorrelation values
    """
    F = np.fft.fft(indicators, axis=1)
    power = np.abs(F) ** 2
    result = np.real(np.fft.ifft(power, axis=1))
    return np.rint(result).astype(np.int32)


def autocorrelation(S, m):
    indicator = np.zeros(m)
    for x in S:
        indicator[x % m] = 1
    F = np.fft.fft(indicator)
    power = np.abs(F) ** 2
    result = np.real(np.fft.ifft(power))
    return np.rint(result).astype(int)


def generate_random_d12_batch(p, d12_size, M):
    """Generate M uniform random D12 subsets of Z_p with given size.

    Returns indicator matrix of shape (M, p).
    """
    # Use argpartition for fast uniform subset selection
    rand = np.random.rand(M, p)
    sorted_indices = np.argpartition(rand, d12_size, axis=1)[:, :d12_size]
    indicators = np.zeros((M, p), dtype=np.float64)
    rows = np.repeat(np.arange(M), d12_size)
    cols = sorted_indices.ravel()
    indicators[rows, cols] = 1.0
    return indicators


def generate_random_symmetric_d11(p, k, M):
    """Generate M uniform random symmetric D11 ⊂ {1,...,p-1} of size k.

    Returns list of D11 sets.
    """
    pairs = [(x, p - x) for x in range(1, (p + 1) // 2)]
    num_pairs = k // 2
    total_pairs = len(pairs)

    d11_list = []
    for _ in range(M):
        selected = np.random.choice(total_pairs, num_pairs, replace=False)
        D11 = set()
        for i in selected:
            D11.add(pairs[i][0])
            D11.add(pairs[i][1])
        d11_list.append(D11)
    return d11_list


def check_d11_with_random_d12s(D11, p, n, M_d12):
    """Check M_d12 random D12s against a fixed D11.

    Returns number of valid D12s found.
    """
    D11_set = set(D11)
    k = len(D11_set)
    d12_size = n - 1

    A = autocorrelation(list(D11_set), p)

    d11_thresh = n - 2
    d22_thresh = 2 * k - n + 1

    if d22_thresh < 0:
        return 0

    threshold = np.zeros(p, dtype=np.int32)
    for d in range(1, p):
        threshold[d] = d11_thresh if d in D11_set else d22_thresh

    # Process in batches to control memory
    batch_size = min(M_d12, 100000)
    total_valid = 0

    for start in range(0, M_d12, batch_size):
        end = min(start + batch_size, M_d12)
        M_batch = end - start

        indicators = generate_random_d12_batch(p, d12_size, M_batch)
        B_all = autocorrelation_batch(indicators, p)

        sum_ab = B_all[:, 1:] + A[np.newaxis, 1:]
        thresh_row = threshold[np.newaxis, 1:]
        valid = np.all(sum_ab <= thresh_row, axis=1)
        total_valid += int(np.sum(valid))

    return total_valid


def load_registry():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(script_dir, 'solutions_registry.json')
    try:
        with open(path, 'r') as f:
            return json.load(f)
    except FileNotFoundError:
        return {'solutions': []}


def is_prime(n):
    if n < 2:
        return False
    for i in range(2, int(n**0.5) + 1):
        if n % i == 0:
            return False
    return True


def main():
    print("=" * 80)
    print("PROBABILISTIC ESTIMATION OF VALID (D11, D12) PAIR DENSITY")
    print("=" * 80)

    data = load_registry()

    # Build D11 lookup for primes with known solutions
    known_d11 = {}
    for sol in data['solutions']:
        m = sol['m']
        if is_prime(m) and m % 4 == 3:
            known_d11[m] = {
                'D11': sol['D11'],
                'D22': sol.get('D22', []),
                'n': sol['n'],
            }

    primes = [7, 11, 19, 23, 31, 43, 47, 59]
    M_d12 = 500000  # D12 samples per D11

    # ==================== Phase 1: Known D11 estimation ====================
    print("\n" + "=" * 80)
    print("PHASE 1: Per-D12 success rate with known/best D11 at k=n-2")
    print("=" * 80)

    print(f"\n  Sampling {M_d12:,} random D12s per D11")
    print(f"\n  {'p':>4s} {'n':>4s} {'k':>4s} {'valid':>8s} "
          f"{'rate':>10s} {'95% CI':>16s} {'source':>20s}")
    print("  " + "-" * 70)

    # For exact results (small primes), use exhaustive counts
    exact_counts = {
        11: {'k': 4, 'best_d11': [1, 2, 9, 10], 'valid': 44,
             'total': 462, 'source': 'exhaustive'},
        19: {'k': 8, 'best_d11': [1, 2, 3, 9, 10, 16, 17, 18],
             'valid': 38, 'total': 92378, 'source': 'exhaustive'},
        23: {'k': 10,
             'best_d11': [1, 2, 3, 4, 10, 13, 19, 20, 21, 22],
             'valid': 414, 'total': 1352078, 'source': 'exhaustive'},
    }

    results = {}

    for p in primes:
        n = (p + 1) // 2
        k_opt = n - 2  # optimal |D11|
        d12_size = n - 1

        if p in exact_counts:
            ec = exact_counts[p]
            D11 = ec['best_d11']
            valid = ec['valid']
            total = ec['total']
            rate = valid / total
            # Wilson score interval
            z = 1.96
            denom = 1 + z**2 / total
            center = (rate + z**2 / (2*total)) / denom
            margin = z * np.sqrt((rate*(1-rate) + z**2/(4*total))/total) / denom
            ci_lo, ci_hi = max(0, center - margin), center + margin
            source = ec['source']
        elif p == 7:
            D11 = [1, 2, 5, 6]
            valid = 0
            total = 35
            rate = 0
            ci_lo, ci_hi = 0, 0.082  # upper bound
            source = 'exhaustive (none)'
        else:
            # Use known D11 or its complement at k=n-2
            if p in known_d11:
                sol = known_d11[p]
                if len(sol['D11']) == k_opt:
                    D11 = sol['D11']
                    source = 'registry D11'
                elif sol['D22'] and len(sol['D22']) == k_opt:
                    D11 = sol['D22']
                    source = 'registry D22→D11'
                else:
                    # Generate a random D11 at k=n-2
                    d11_list = generate_random_symmetric_d11(p, k_opt, 50)
                    best_d11 = None
                    best_valid = -1
                    for d11 in d11_list:
                        v = check_d11_with_random_d12s(
                            d11, p, n, min(M_d12, 50000))
                        if v > best_valid:
                            best_valid = v
                            best_d11 = d11
                    D11 = sorted(best_d11)
                    source = f'best of 50 random'
            else:
                # No registry; try random D11s
                d11_list = generate_random_symmetric_d11(p, k_opt, 100)
                best_d11 = None
                best_valid = -1
                for d11 in d11_list:
                    v = check_d11_with_random_d12s(d11, p, n, 50000)
                    if v > best_valid:
                        best_valid = v
                        best_d11 = d11
                D11 = sorted(best_d11)
                source = f'best of 100 random'

            # Full sampling with selected D11
            t0 = time.time()
            valid = check_d11_with_random_d12s(D11, p, n, M_d12)
            elapsed = time.time() - t0
            total = M_d12
            rate = valid / total
            z = 1.96
            margin = z * np.sqrt(rate * (1 - rate) / total) if rate > 0 else 0
            ci_lo = max(0, rate - margin)
            ci_hi = rate + margin
            if valid == 0:
                ci_hi = 3.0 / total  # rule of 3
            source += f' [{elapsed:.1f}s]'

        k = len(D11) if isinstance(D11, list) else k_opt

        ci_str = f"[{ci_lo:.2e}, {ci_hi:.2e}]"
        print(f"  {p:4d} {n:4d} {k:4d} {valid:8d} "
              f"{rate:10.4e} {ci_str:>16s} {source:>20s}")

        results[p] = {
            'n': n, 'k': k, 'valid': valid, 'total': total,
            'rate': rate, 'ci': (ci_lo, ci_hi), 'D11': D11,
        }

    # ==================== Phase 2: Per-D11 success rate ====================
    print("\n" + "=" * 80)
    print("PHASE 2: Per-D11 success rate (random D11 at k=n-2)")
    print("=" * 80)

    M_d11 = 200   # D11s to sample
    M_d12_quick = 100000  # D12s per D11

    print(f"\n  Sampling {M_d11} random D11s, "
          f"{M_d12_quick:,} D12s per D11")
    print(f"\n  {'p':>4s} {'n':>4s} {'D11s sampled':>12s} "
          f"{'D11s w/valid':>12s} {'rate':>8s} "
          f"{'avg valid D12':>14s}")
    print("  " + "-" * 60)

    for p in primes:
        if p == 7:
            print(f"  {p:4d} {4:4d}          3            0     0.0%              0")
            continue

        n = (p + 1) // 2
        k = n - 2

        if k < 2:
            continue

        t0 = time.time()
        d11_list = generate_random_symmetric_d11(p, k, M_d11)

        d11_with_valid = 0
        total_valid = 0

        for i, D11 in enumerate(d11_list):
            valid = check_d11_with_random_d12s(D11, p, n, M_d12_quick)
            total_valid += valid
            if valid > 0:
                d11_with_valid += 1

            if (i + 1) % 50 == 0:
                elapsed = time.time() - t0
                print(f"    [{elapsed:.0f}s] p={p}: {i+1}/{M_d11} D11s, "
                      f"{d11_with_valid} with valid D12",
                      flush=True)

        elapsed = time.time() - t0
        d11_rate = d11_with_valid / M_d11 * 100
        avg_valid = total_valid / M_d11

        print(f"  {p:4d} {n:4d} {M_d11:12d} {d11_with_valid:12d} "
              f"{d11_rate:7.1f}% {avg_valid:14.1f}  [{elapsed:.0f}s]")

    # ==================== Phase 3: Growth trend ====================
    print("\n" + "=" * 80)
    print("PHASE 3: Growth trend analysis")
    print("=" * 80)

    # Exact data from survey
    exact_data = [
        # (p, n, num_d11_at_k, num_d12, total_valid, d11s_with_valid)
        (11, 6, 10, 462, 220, 5),
        (19, 10, 126, 92378, 342, 9),
        (23, 12, 462, 1352078, 9108, 55),
    ]

    print(f"\n  Exact counts from exhaustive survey (k=n-2):")
    print(f"  {'p':>4s} {'n':>4s} "
          f"{'C(p,n-1)':>12s} {'valid pairs':>12s} "
          f"{'per-pair rate':>14s} "
          f"{'D11s ok':>8s} {'D11 rate':>10s} "
          f"{'avg D12/D11':>12s}")
    print("  " + "-" * 80)

    for p, n, num_d11, num_d12, total_valid, d11s_ok in exact_data:
        pair_rate = total_valid / (num_d11 * num_d12)
        d11_rate = d11s_ok / num_d11 * 100
        avg_d12 = total_valid / num_d11

        print(f"  {p:4d} {n:4d} "
              f"{num_d12:12,d} {total_valid:12,d} "
              f"{pair_rate:14.2e} "
              f"{d11s_ok:8d} {d11_rate:9.1f}% "
              f"{avg_d12:12.1f}")

    print(f"\n  Key metric: avg valid D12s per random D11")
    print(f"    p=11: 22.0  (growing)")
    print(f"    p=19:  2.7  (dropped)")
    print(f"    p=23: 19.7  (recovered!)")
    print(f"\n  If this metric stays bounded away from 0,")
    print(f"  existence of valid pairs is guaranteed for large p.")

    # Theoretical analysis
    print(f"\n  Theoretical average A+B analysis (k=n-2):")
    for p in [11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83]:
        n = (p + 1) // 2
        k = n - 2
        avg_ab = (k*(k-1) + (n-1)*(n-2)) / (p - 1)
        d11_thresh = n - 2
        d22_thresh = 2*k - n + 1
        d11_slack = d11_thresh - avg_ab
        d22_slack = d22_thresh - avg_ab
        print(f"    p={p:3d}: avg A+B = {avg_ab:.2f}, "
              f"D11 thresh={d11_thresh} (slack {d11_slack:+.2f}), "
              f"D22 thresh={d22_thresh} (slack {d22_slack:+.2f})")


if __name__ == '__main__':
    main()
