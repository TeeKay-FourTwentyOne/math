#!/usr/bin/env python3
"""
For each prime p ≡ 3 mod 4, find the symmetric D11 of size (p+1)/2
that minimizes max_d A(d) over ALL d in {1,...,p-1}.

Report the gap = min_{D11} max_d A(d) - E[A] where E[A] = (p+1)/4.

If this gap stays bounded (≤ 2 or 3) as p grows, the proof closes.
"""

import numpy as np
from itertools import combinations
from math import comb, floor
from collections import defaultdict
import time


def get_symmetric_pairs(p):
    pairs = []
    seen = set()
    for d in range(1, p):
        if d not in seen:
            pairs.append((d, (-d) % p))
            seen.add(d)
            seen.add((-d) % p)
    return pairs


def enumerate_symmetric_D11(p, size):
    pairs = get_symmetric_pairs(p)
    num_pairs = size // 2
    for chosen_pairs in combinations(range(len(pairs)), num_pairs):
        D11 = set()
        for i in chosen_pairs:
            d, neg_d = pairs[i]
            D11.add(d)
            D11.add(neg_d)
        yield frozenset(D11)


def batch_autocorrelation_max(p, size, sample_limit=None):
    """Compute max_d A(d) for all (or sampled) symmetric D11.

    Returns list of (max_A, D11_set) sorted by max_A.
    """
    pairs = get_symmetric_pairs(p)
    num_pairs = len(pairs)
    num_choose = size // 2
    total = comb(num_pairs, num_choose)

    print(f"  Total symmetric D11: {total}")

    if sample_limit and total > sample_limit:
        print(f"  Sampling {sample_limit} random D11...")
        return sample_autocorrelation_max(p, size, pairs, sample_limit)

    # Enumerate all
    results = []
    batch_size = min(5000, total)
    batch_indicators = np.zeros((batch_size, p), dtype=np.float64)
    batch_d11s = []
    batch_idx = 0

    for d11_idx, D11 in enumerate(enumerate_symmetric_D11(p, size)):
        for j in D11:
            batch_indicators[batch_idx, j] = 1.0
        batch_d11s.append(D11)
        batch_idx += 1

        if batch_idx == batch_size or d11_idx == total - 1:
            # Process batch
            actual_batch = batch_idx
            mat = batch_indicators[:actual_batch]
            fft_vals = np.fft.fft(mat, axis=1)
            autocorr = np.fft.ifft(np.abs(fft_vals) ** 2, axis=1).real
            autocorr = np.round(autocorr).astype(int)

            # max over d=1,...,p-1
            max_A = autocorr[:, 1:].max(axis=1)

            for i in range(actual_batch):
                results.append((int(max_A[i]), batch_d11s[i]))

            # Reset batch
            batch_indicators[:] = 0
            batch_d11s = []
            batch_idx = 0

        if (d11_idx + 1) % 100000 == 0:
            print(f"    {d11_idx+1}/{total}...")

    results.sort(key=lambda x: x[0])
    return results


def sample_autocorrelation_max(p, size, pairs, num_samples):
    """Sample random symmetric D11 and compute max_d A(d)."""
    num_pairs = len(pairs)
    num_choose = size // 2

    results = []
    batch_size = 5000

    for start in range(0, num_samples, batch_size):
        actual = min(batch_size, num_samples - start)
        indicators = np.zeros((actual, p), dtype=np.float64)
        d11s = []

        for i in range(actual):
            chosen = np.random.choice(num_pairs, num_choose, replace=False)
            D11 = set()
            for idx in chosen:
                d, neg_d = pairs[idx]
                D11.add(d)
                D11.add(neg_d)
            d11s.append(frozenset(D11))
            for j in D11:
                indicators[i, j] = 1.0

        fft_vals = np.fft.fft(indicators, axis=1)
        autocorr = np.fft.ifft(np.abs(fft_vals) ** 2, axis=1).real
        autocorr = np.round(autocorr).astype(int)
        max_A = autocorr[:, 1:].max(axis=1)

        for i in range(actual):
            results.append((int(max_A[i]), d11s[i]))

    results.sort(key=lambda x: x[0])
    return results


def analyze_best_d11(p, D11_set):
    """Detailed analysis of a specific D11."""
    D11 = set(D11_set)
    D22 = set(range(1, p)) - D11

    indicator = np.zeros(p, dtype=np.float64)
    for j in D11:
        indicator[j] = 1.0
    fft_val = np.fft.fft(indicator)
    A = np.round(np.fft.ifft(np.abs(fft_val) ** 2).real).astype(int)

    A_at_D11 = sorted([int(A[d]) for d in D11])
    A_at_D22 = sorted([int(A[d]) for d in D22])
    max_A_all = max(int(A[d]) for d in range(1, p))
    max_A_D11 = max(int(A[d]) for d in D11)
    max_A_D22 = max(int(A[d]) for d in D22)

    QR = {(x * x) % p for x in range(1, p)}
    qr_count = len(D11 & QR)

    return {
        "max_A_all": max_A_all,
        "max_A_D11": max_A_D11,
        "max_A_D22": max_A_D22,
        "A_profile_D11": A_at_D11,
        "A_profile_D22": A_at_D22,
        "qr_count": qr_count,
        "qnr_count": len(D11) - qr_count,
    }


def main():
    primes = [11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83]

    print("=" * 80)
    print("GLOBAL FLATNESS: min_{D11} max_d A(d) for symmetric D11")
    print("=" * 80)

    summary = []

    for p in primes:
        n = (p + 1) // 2
        E_A = (p + 1) / 4

        pairs = get_symmetric_pairs(p)
        total = comb(len(pairs), n // 2)

        # Decide enumeration vs sampling
        sample_limit = None
        if total > 2_000_000:
            sample_limit = 500_000

        print(f"\n{'=' * 60}")
        print(f"p = {p}, n = {n}, E[A] = {E_A}")

        t0 = time.time()
        results = batch_autocorrelation_max(p, n, sample_limit)
        elapsed = time.time() - t0

        # Statistics
        all_max = [r[0] for r in results]
        min_max_A = min(all_max)
        gap = min_max_A - E_A

        # Distribution of max_A values
        hist = defaultdict(int)
        for m in all_max:
            hist[m] += 1

        # Count how many achieve the minimum
        num_best = hist[min_max_A]

        print(f"  Time: {elapsed:.1f}s")
        print(f"  min(max_d A(d)) = {min_max_A}")
        print(f"  E[A] = {E_A}")
        print(f"  Gap = min - E[A] = {gap}")
        print(f"  #D11 achieving min = {num_best} / {len(results)}")
        print(f"  Fraction achieving min: {num_best/len(results)*100:.2f}%")
        if sample_limit:
            print(f"  (Based on {sample_limit} random samples)")

        # Distribution
        print(f"  max_A distribution:")
        for val in sorted(hist.keys())[:10]:
            pct = hist[val] / len(results) * 100
            print(f"    max_A = {val}: {hist[val]} ({pct:.1f}%)")
        if len(hist) > 10:
            print(f"    ... ({len(hist)} distinct values total)")

        # Analyze the best D11
        best_D11 = results[0][1]
        info = analyze_best_d11(p, best_D11)
        print(f"  Best D11 analysis:")
        print(f"    D11 = {sorted(best_D11)}")
        print(f"    max_A_D11 = {info['max_A_D11']}, max_A_D22 = {info['max_A_D22']}")
        print(f"    A@D11 = {info['A_profile_D11']}")
        print(f"    A@D22 = {info['A_profile_D22']}")
        print(f"    QR/QNR = {info['qr_count']}/{info['qnr_count']}")

        # Also show a few more best D11 to see if they're all in the same orbit
        if num_best > 1:
            # Check if all best D11 have the same A-profile
            profiles = set()
            for max_val, D11_set in results[:num_best+5]:
                if max_val == min_max_A:
                    info2 = analyze_best_d11(p, D11_set)
                    profiles.add(tuple(info2['A_profile_D11'] + info2['A_profile_D22']))
            print(f"    Distinct (D11,D22) profiles among best: {len(profiles)}")

        summary.append({
            "p": p,
            "n": n,
            "E_A": E_A,
            "min_max_A": min_max_A,
            "gap": gap,
            "num_best": num_best,
            "total": len(results),
            "sampled": sample_limit is not None,
        })

    # Final summary table
    print(f"\n{'=' * 80}")
    print("SUMMARY TABLE")
    print(f"{'=' * 80}")
    print(f"{'p':>5} {'n':>4} {'E[A]':>7} {'min max_A':>10} {'gap':>6} {'#best':>8} {'total':>10} {'method':>10}")
    print("-" * 70)
    for s in summary:
        method = "sample" if s["sampled"] else "exact"
        print(f"{s['p']:>5} {s['n']:>4} {s['E_A']:>7.1f} {s['min_max_A']:>10} {s['gap']:>6.1f} {s['num_best']:>8} {s['total']:>10} {method:>10}")

    # Key question
    print(f"\n{'=' * 80}")
    print("KEY QUESTION: Does the gap stay bounded?")
    print(f"{'=' * 80}")
    gaps = [s["gap"] for s in summary]
    print(f"Gaps: {[s['gap'] for s in summary]}")
    if max(gaps) <= 3:
        print("YES — gap ≤ 3 for all tested primes. Proof may close!")
    elif max(gaps) <= 5:
        print("MAYBE — gap ≤ 5, grows slowly. Need to check larger primes.")
    else:
        print("UNCLEAR — gap is growing. Need theoretical analysis.")


if __name__ == "__main__":
    main()
