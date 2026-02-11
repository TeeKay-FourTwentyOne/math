#!/usr/bin/env python3
"""Fast search for valid (D11, D12) pairs at p=31.

Uses k=n (|D11|=(p+1)/2=16) formulation where only D11 constraints bind.
This is complement-equivalent to k=n-2=14 by Theorem 5.

Strategy: For each random symmetric D11, test many random D12s.
When a valid pair is found, analyze its structure in detail.
"""

import numpy as np
from math import comb
import time
import json


def autocorrelation_fft(indicator, p):
    F = np.fft.fft(indicator)
    power = np.abs(F) ** 2
    return np.round(np.fft.ifft(power).real).astype(int)


def quadratic_residues(p):
    return sorted(set((x * x) % p for x in range(1, p)))


def generate_random_symmetric_d11(p, rng):
    pairs = [(x, p - x) for x in range(1, (p + 1) // 2)]
    num_pairs = (p + 1) // 4
    chosen = rng.choice(len(pairs), size=num_pairs, replace=False)
    D11 = set()
    for i in chosen:
        D11.add(pairs[i][0])
        D11.add(pairs[i][1])
    return D11


def search_p31():
    p = 31
    n = 16
    d11_size = n  # (p+1)/2 = 16
    d12_size = n - 1  # (p-1)/2 = 15
    d11_thresh = n - 2  # 14
    d22_thresh = n + 1  # 17 (loose)

    print(f"p={p}, n={n}")
    print(f"|D11|={d11_size} (k=n formulation)")
    print(f"|D12|={d12_size}")
    print(f"D11 threshold: A+B <= {d11_thresh}")
    print(f"D22 threshold: A+B <= {d22_thresh} (loose)")
    print(f"# symmetric D11 = C({(p-1)//2},{(p+1)//4}) = "
          f"{comb((p-1)//2, (p+1)//4)}")

    rng = np.random.default_rng(42)
    qr = set(quadratic_residues(p))
    qnr = set(range(1, p)) - qr

    M_d11 = 2000
    M_d12 = 200000
    batch_size = 50000

    found = []
    d11_tested = 0

    t0 = time.time()
    print(f"\nTesting {M_d11} random D11s with {M_d12:,} D12s each...")

    for trial in range(M_d11):
        D11 = generate_random_symmetric_d11(p, rng)
        D11_sorted = sorted(D11)

        ind11 = np.zeros(p, dtype=np.float64)
        for x in D11:
            ind11[x] = 1.0
        A = autocorrelation_fft(ind11, p)

        threshold = np.zeros(p, dtype=np.int32)
        for d in range(1, p):
            threshold[d] = d11_thresh if d in D11 else d22_thresh

        total_valid = 0
        first_valid_d12 = None

        for start in range(0, M_d12, batch_size):
            end = min(start + batch_size, M_d12)
            M_batch = end - start

            rand_vals = rng.random((M_batch, p - 1))
            sorted_idx = np.argpartition(rand_vals, d12_size - 1,
                                         axis=1)[:, :d12_size - 1]

            indicators = np.zeros((M_batch, p), dtype=np.float64)
            indicators[:, 0] = 1.0
            rows = np.repeat(np.arange(M_batch), d12_size - 1)
            cols = sorted_idx.ravel() + 1
            indicators[rows, cols] = 1.0

            F = np.fft.fft(indicators, axis=1)
            power = np.abs(F) ** 2
            B_all = np.round(np.fft.ifft(power, axis=1).real).astype(np.int32)

            sum_ab = B_all[:, 1:] + A[np.newaxis, 1:]
            thresh_row = threshold[np.newaxis, 1:]
            valid = np.all(sum_ab <= thresh_row, axis=1)
            cnt = int(np.sum(valid))
            total_valid += cnt

            if cnt > 0 and first_valid_d12 is None:
                idx = np.argmax(valid)
                d12_elements = sorted(np.where(indicators[idx] > 0.5)[0].tolist())
                first_valid_d12 = d12_elements

        d11_tested += 1

        if total_valid > 0:
            found.append({
                'D11': D11_sorted,
                'D12': first_valid_d12,
                'valid_count': total_valid,
                'rate': total_valid / M_d12,
            })
            elapsed = time.time() - t0
            print(f"  FOUND #{len(found)}: "
                  f"{total_valid}/{M_d12} valid D12s "
                  f"(rate={total_valid/M_d12:.2e}) "
                  f"[{elapsed:.0f}s, {d11_tested} D11s]")

        if d11_tested % 200 == 0:
            elapsed = time.time() - t0
            print(f"  [{elapsed:.0f}s] {d11_tested}/{M_d11} D11s, "
                  f"{len(found)} found")

    elapsed = time.time() - t0
    print(f"\n{'='*70}")
    print(f"RESULTS: {len(found)}/{d11_tested} D11s have valid D12 "
          f"({len(found)/d11_tested*100:.1f}%) [{elapsed:.0f}s]")
    print(f"{'='*70}")

    if not found:
        print("\nNo valid pairs found. Try more samples or different strategy.")
        return

    # Analyze structure
    print(f"\nStructural analysis of {len(found)} valid D11s:")

    for i, pair in enumerate(found[:10]):
        D11 = set(pair['D11'])
        D12 = set(pair['D12'])

        d11_in_qr = len(D11 & qr)
        d11_in_qnr = len(D11 & qnr)
        d12_in_qr = len(D12 & qr)
        d12_has_0 = 0 in D12

        print(f"\n  #{i+1}: D11={pair['D11']}")
        print(f"       |D11∩QR|={d11_in_qr}, |D11∩QNR|={d11_in_qnr}")
        print(f"       D12={pair['D12']}")
        print(f"       |D12∩QR|={d12_in_qr}, 0∈D12={d12_has_0}")
        print(f"       valid D12 rate: {pair['rate']:.2e}")

    # Full verification of first solution
    print(f"\n{'='*70}")
    print(f"FULL VERIFICATION of solution #1")
    print(f"{'='*70}")

    pair = found[0]
    D11 = pair['D11']
    D12 = pair['D12']
    D11_set = set(D11)

    ind11 = np.zeros(p, dtype=np.float64)
    for x in D11:
        ind11[x] = 1.0
    A = autocorrelation_fft(ind11, p)

    ind12 = np.zeros(p, dtype=np.float64)
    for x in D12:
        ind12[x] = 1.0
    B = autocorrelation_fft(ind12, p)

    print(f"\n  A(d) + B(d) at D11 positions (threshold {d11_thresh}):")
    d11_vals = [(d, int(A[d] + B[d])) for d in sorted(D11_set)]
    for d, v in d11_vals:
        status = "OK" if v <= d11_thresh else "FAIL"
        print(f"    d={d:2d}: A+B={v:3d}  {status}")
    max_d11 = max(v for _, v in d11_vals)

    D22_set = set(range(1, p)) - D11_set
    print(f"\n  A(d) + B(d) at D22 positions (threshold {d22_thresh}):")
    d22_vals = [(d, int(A[d] + B[d])) for d in sorted(D22_set)]
    for d, v in d22_vals:
        status = "OK" if v <= d22_thresh else "FAIL"
        print(f"    d={d:2d}: A+B={v:3d}  {status}")
    max_d22 = max(v for _, v in d22_vals)

    if max_d11 <= d11_thresh and max_d22 <= d22_thresh:
        print(f"\n  VERIFIED: Valid construction for R(B_15, B_16) >= 63")

        # Also verify complement (k=n-2) form
        D11_comp = sorted(set(range(1, p)) - D11_set)
        D12_comp = sorted({(-x) % p for x in D12})
        print(f"\n  Complement form (k=n-2=14):")
        print(f"    D11' = {D11_comp}")
        print(f"    D12' = {D12_comp}")

    # Save results
    output = {
        'p': p,
        'n': n,
        'd11_tested': d11_tested,
        'found': len(found),
        'd11_success_rate': len(found) / d11_tested,
        'solutions': [{
            'D11': s['D11'],
            'D12': s['D12'],
            'valid_d12_count': s['valid_count'],
            'valid_d12_rate': s['rate'],
        } for s in found[:20]],
    }

    output_path = 'p31_solutions.json'
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\nResults saved to {output_path}")


if __name__ == '__main__':
    search_p31()
