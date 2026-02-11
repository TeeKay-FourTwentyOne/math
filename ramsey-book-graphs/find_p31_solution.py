#!/usr/bin/env python3
"""Find valid (D11, D12) pairs for p=31 at k=n-2=14.

Sampling found that ~2% of random D11s admit valid D12s.
This script finds them explicitly and analyzes their structure.
"""

import numpy as np
from itertools import combinations
import time


def autocorrelation(S, m):
    indicator = np.zeros(m)
    for x in S:
        indicator[x % m] = 1
    F = np.fft.fft(indicator)
    power = np.abs(F) ** 2
    result = np.real(np.fft.ifft(power))
    return np.rint(result).astype(int)


def autocorrelation_batch(indicators, p):
    F = np.fft.fft(indicators, axis=1)
    power = np.abs(F) ** 2
    result = np.real(np.fft.ifft(power, axis=1))
    return np.rint(result).astype(np.int32)


def generate_random_d12_batch(p, d12_size, M):
    rand = np.random.rand(M, p)
    sorted_indices = np.argpartition(rand, d12_size, axis=1)[:, :d12_size]
    indicators = np.zeros((M, p), dtype=np.float64)
    rows = np.repeat(np.arange(M), d12_size)
    cols = sorted_indices.ravel()
    indicators[rows, cols] = 1.0
    return indicators, sorted_indices


def quadratic_residues(p):
    return sorted(set((x * x) % p for x in range(1, p)))


def search_p31():
    p = 31
    n = 16
    k = 14  # n-2
    d12_size = 15  # n-1
    d11_thresh = n - 2  # 14
    d22_thresh = 2 * k - n + 1  # 13

    print(f"p={p}, n={n}, k={k} (=n-2)")
    print(f"D11 threshold: A+B <= {d11_thresh}")
    print(f"D22 threshold: A+B <= {d22_thresh}")
    print(f"|D12| = {d12_size}")

    # Generate symmetric D11 candidates
    pairs = [(x, p - x) for x in range(1, (p + 1) // 2)]
    # k=14 means 7 pairs out of 15
    num_pairs_needed = k // 2  # 7
    total_d11 = int(np.math.factorial(len(pairs))
                    / np.math.factorial(num_pairs_needed)
                    / np.math.factorial(len(pairs) - num_pairs_needed))
    print(f"Symmetric D11 candidates: C(15,7) = {total_d11}")

    M_d12 = 500000  # D12 samples per D11
    batch_size = 100000

    found_pairs = []
    d11_tested = 0

    print(f"\nSearching with {M_d12:,} random D12s per D11...")
    print()

    t0 = time.time()

    for combo in combinations(range(len(pairs)), num_pairs_needed):
        D11 = []
        for i in combo:
            D11.append(pairs[i][0])
            D11.append(pairs[i][1])
        D11 = sorted(D11)
        D11_set = set(D11)

        A = autocorrelation(D11, p)

        threshold = np.zeros(p, dtype=np.int32)
        for d in range(1, p):
            threshold[d] = d11_thresh if d in D11_set else d22_thresh

        total_valid = 0
        first_valid_d12 = None

        for start in range(0, M_d12, batch_size):
            end = min(start + batch_size, M_d12)
            M_batch = end - start

            indicators, indices = generate_random_d12_batch(
                p, d12_size, M_batch)
            B_all = autocorrelation_batch(indicators, p)

            sum_ab = B_all[:, 1:] + A[np.newaxis, 1:]
            thresh_row = threshold[np.newaxis, 1:]
            valid = np.all(sum_ab <= thresh_row, axis=1)
            cnt = int(np.sum(valid))
            total_valid += cnt

            if cnt > 0 and first_valid_d12 is None:
                idx = np.argmax(valid)
                first_valid_d12 = sorted(indices[idx].tolist())

        d11_tested += 1

        if total_valid > 0:
            found_pairs.append({
                'D11': D11,
                'D12': first_valid_d12,
                'valid_count': total_valid,
                'rate': total_valid / M_d12,
            })
            elapsed = time.time() - t0
            print(f"  FOUND #{len(found_pairs)}: D11={D11}")
            print(f"    D12={first_valid_d12}")
            print(f"    {total_valid}/{M_d12} valid D12s "
                  f"(rate={total_valid/M_d12:.2e})")
            print(f"    [{elapsed:.0f}s, {d11_tested}/{total_d11} D11s tested]")

        if d11_tested % 500 == 0:
            elapsed = time.time() - t0
            print(f"  [{elapsed:.0f}s] {d11_tested}/{total_d11} D11s, "
                  f"{len(found_pairs)} found")

        # Stop after finding enough or testing all
        if len(found_pairs) >= 20:
            print(f"\n  Found {len(found_pairs)} valid D11s, stopping search.")
            break

    elapsed = time.time() - t0
    print(f"\n{'='*70}")
    print(f"RESULTS: {len(found_pairs)} valid D11s found in "
          f"{d11_tested}/{total_d11} tested ({elapsed:.0f}s)")
    print(f"{'='*70}")

    if found_pairs:
        # Analyze structure
        qr = set(quadratic_residues(p))
        qnr = set(range(1, p)) - qr

        print(f"\nStructural analysis of valid D11s:")
        for i, pair in enumerate(found_pairs):
            D11 = set(pair['D11'])
            d11_in_qr = len(D11 & qr)
            d11_in_qnr = len(D11 & qnr)

            D12 = set(pair['D12'])
            d12_in_qr = len(D12 & qr)
            has_0 = 0 in D12

            print(f"\n  #{i+1}: D11={pair['D11']}")
            print(f"       |D11 n QR|={d11_in_qr}, "
                  f"|D11 n QNR|={d11_in_qnr}")
            print(f"       D12={pair['D12']}")
            print(f"       |D12 n QR|={d12_in_qr}, 0 in D12={has_0}")
            print(f"       valid D12 rate: {pair['rate']:.2e}")

        # Verify one solution fully
        print(f"\nFull verification of solution #1:")
        pair = found_pairs[0]
        D11 = pair['D11']
        D12 = pair['D12']
        D11_set = set(D11)
        A = autocorrelation(D11, p)
        B = autocorrelation(D12, p)

        print(f"  A(d) + B(d) at D11 positions (threshold {d11_thresh}):")
        d11_vals = [(d, int(A[d]+B[d])) for d in sorted(D11_set)]
        print(f"    {d11_vals}")
        max_d11 = max(v for _, v in d11_vals)
        print(f"    max = {max_d11} {'OK' if max_d11 <= d11_thresh else 'FAIL'}")

        D22_set = set(range(1, p)) - D11_set
        print(f"  A(d) + B(d) at D22 positions (threshold {d22_thresh}):")
        d22_vals = [(d, int(A[d]+B[d])) for d in sorted(D22_set)]
        print(f"    {d22_vals}")
        max_d22 = max(v for _, v in d22_vals)
        print(f"    max = {max_d22} {'OK' if max_d22 <= d22_thresh else 'FAIL'}")

        if max_d11 <= d11_thresh and max_d22 <= d22_thresh:
            print(f"\n  VERIFIED: Valid construction for R(B_15, B_16) >= 63")


if __name__ == '__main__':
    search_p31()
