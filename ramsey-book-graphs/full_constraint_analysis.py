"""
Full constraint analysis: diagnose why "loose" V1V1 blue / V2V2 red constraints
are actually binding in practice.

Key finding from Phases 1-3: The correlation analysis (which checks only V1V1 red)
finds ~2% hit rate, but the sampling (which checks all constraints) finds ~0%.
The gap comes from D22 constraints being violated.

This script:
1. For each prime, shows A(d) distribution at D11 vs D22 positions
2. Computes effective B(d) threshold at each position
3. Runs MC checking ALL constraints
4. Compares "good" D11 (QR-balanced) vs "bad" D11
"""

import sys
import os
import json
import time
import math
import numpy as np
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ramsey_core import BlockCirculantGraph, verify_construction, Delta


def autocorrelation_fft(indicator, m):
    """Compute Delta(S,S,d) for all d via FFT."""
    fft_val = np.fft.fft(indicator)
    autocorr = np.fft.ifft(np.abs(fft_val) ** 2).real
    return np.round(autocorr).astype(int)


def get_symmetric_pairs(p):
    pairs = []
    seen = set()
    for d in range(1, p):
        if d not in seen:
            pairs.append((d, (-d) % p))
            seen.add(d)
            seen.add((-d) % p)
    return pairs


def random_symmetric_d11(p, size, rng):
    pairs = get_symmetric_pairs(p)
    chosen = rng.choice(len(pairs), size=size // 2, replace=False)
    D11 = set()
    for i in chosen:
        D11.add(pairs[i][0])
        D11.add(pairs[i][1])
    return frozenset(D11)


def analyze_d11(p, D11):
    """Analyze A(d) distribution for a given D11."""
    indicator = np.zeros(p, dtype=np.float64)
    for j in D11:
        indicator[j] = 1.0
    A = autocorrelation_fft(indicator, p)

    D22 = set(range(1, p)) - D11
    threshold_tight = (p - 3) // 2
    threshold_loose = (p + 3) // 2

    # A(d) at D11 vs D22 positions
    A_d11 = [int(A[d]) for d in sorted(D11)]
    A_d22 = [int(A[d]) for d in sorted(D22)]

    # Effective B threshold at each position
    B_thresh_d11 = [threshold_tight - int(A[d]) for d in sorted(D11)]
    B_thresh_d22 = [threshold_loose - int(A[d]) for d in sorted(D22)]

    # Min effective threshold (most constrained position)
    min_B_d11 = min(B_thresh_d11) if B_thresh_d11 else float("inf")
    min_B_d22 = min(B_thresh_d22) if B_thresh_d22 else float("inf")
    overall_min_B = min(min_B_d11, min_B_d22)
    bottleneck = "D11" if min_B_d11 <= min_B_d22 else "D22"

    # QR distribution
    QR = set()
    for x in range(1, p):
        QR.add((x * x) % p)
    d11_qr = len(D11 & QR)
    d11_qnr = len(D11) - d11_qr

    return {
        "A_d11": A_d11,
        "A_d22": A_d22,
        "mean_A_d11": sum(A_d11) / len(A_d11) if A_d11 else 0,
        "mean_A_d22": sum(A_d22) / len(A_d22) if A_d22 else 0,
        "max_A_d22": max(A_d22) if A_d22 else 0,
        "B_thresh_d11": B_thresh_d11,
        "B_thresh_d22": B_thresh_d22,
        "min_B_d11": min_B_d11,
        "min_B_d22": min_B_d22,
        "overall_min_B": overall_min_B,
        "bottleneck": bottleneck,
        "d11_qr": d11_qr,
        "d11_qnr": d11_qnr,
        "A_array": A,
    }


def mc_full_check(p, D11, A, num_mc=200000):
    """Monte Carlo with ALL constraints checked."""
    d12_size = (p - 1) // 2
    threshold_tight = (p - 3) // 2
    threshold_loose = (p + 3) // 2

    D11_set = set(D11)
    rng = np.random.default_rng(42)

    valid_count = 0
    d11_only_valid = 0  # valid by binding constraints only
    d22_violators = 0  # valid on D11 but fail on D22

    for _ in range(num_mc):
        # Random D12
        others = rng.choice(range(1, p), size=d12_size - 1, replace=False)
        d12_ind = np.zeros(p, dtype=np.float64)
        d12_ind[0] = 1.0
        for x in others:
            d12_ind[x] = 1.0
        B = autocorrelation_fft(d12_ind, p)

        # Check binding constraints (d ∈ D11)
        d11_ok = True
        for d in D11:
            if int(A[d]) + int(B[d]) > threshold_tight:
                d11_ok = False
                break

        if d11_ok:
            d11_only_valid += 1
            # Check D22 constraints
            d22_ok = True
            for d in range(1, p):
                if d not in D11_set:
                    if int(A[d]) + int(B[d]) > threshold_loose:
                        d22_ok = False
                        break
            if d22_ok:
                valid_count += 1
            else:
                d22_violators += 1

    return {
        "valid_all": valid_count,
        "valid_binding_only": d11_only_valid,
        "d22_violators": d22_violators,
        "total": num_mc,
        "rate_all": valid_count / num_mc,
        "rate_binding": d11_only_valid / num_mc,
        "rate_d22_given_d11": (
            d22_violators / d11_only_valid if d11_only_valid > 0 else float("nan")
        ),
    }


# Known solutions
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


def main():
    primes = [11, 19, 23, 31, 43]

    for p in primes:
        n = (p + 1) // 2
        print(f"\n{'=' * 80}")
        print(f"FULL CONSTRAINT ANALYSIS: p = {p}, n = {n}")
        print(f"Binding threshold (D11): {(p - 3) // 2}")
        print(f"Loose threshold (D22):   {(p + 3) // 2}")
        print(f"|D11| = {(p + 1) // 2}, |D22| = {(p - 3) // 2}")
        print(f"{'=' * 80}")

        if p in KNOWN:
            D11 = KNOWN[p]["D11"]
            D12 = KNOWN[p]["D12"]
        else:
            continue

        info = analyze_d11(p, D11)

        print(f"\n  A(d) distribution for KNOWN D11:")
        print(f"    D11 positions: mean={info['mean_A_d11']:.2f}, "
              f"values={info['A_d11']}")
        print(f"    D22 positions: mean={info['mean_A_d22']:.2f}, "
              f"max={info['max_A_d22']}, values={info['A_d22']}")
        print(f"    QR balance: {info['d11_qr']} QR, {info['d11_qnr']} QNR")

        print(f"\n  Effective B(d) thresholds:")
        print(f"    D11: {info['B_thresh_d11']} (min={info['min_B_d11']})")
        print(f"    D22: {info['B_thresh_d22']} (min={info['min_B_d22']})")
        print(f"    BOTTLENECK: {info['bottleneck']}")
        print(f"    Overall min B threshold: {info['overall_min_B']}")

        # Verify known solution
        d12_ind = np.zeros(p, dtype=np.float64)
        for j in D12:
            d12_ind[j] = 1.0
        B_known = autocorrelation_fft(d12_ind, p)

        print(f"\n  Known D12 B(d) values:")
        B_at_d11 = [int(B_known[d]) for d in sorted(D11)]
        B_at_d22 = [int(B_known[d]) for d in sorted(set(range(1, p)) - D11)]
        print(f"    D11: {B_at_d11}")
        print(f"    D22: {B_at_d22}")

        # Check margins for known solution
        margins_d11 = [info['B_thresh_d11'][i] - B_at_d11[i]
                       for i in range(len(B_at_d11))]
        margins_d22 = [info['B_thresh_d22'][i] - B_at_d22[i]
                       for i in range(len(B_at_d22))]
        print(f"    Margins at D11: {margins_d11} (min={min(margins_d11)})")
        print(f"    Margins at D22: {margins_d22} (min={min(margins_d22)})")

        # Monte Carlo with FULL constraint check
        num_mc = min(500000, max(50000, 10000000 // p))
        print(f"\n  Monte Carlo ({num_mc:,} trials) with FULL constraint check:")
        t0 = time.time()
        mc = mc_full_check(p, D11, info["A_array"], num_mc=num_mc)
        elapsed = time.time() - t0

        print(f"    Valid (binding only): {mc['valid_binding_only']}/{mc['total']} "
              f"= {mc['rate_binding']:.6f}")
        print(f"    D22 violators (pass D11, fail D22): {mc['d22_violators']}")
        print(f"    Valid (ALL): {mc['valid_all']}/{mc['total']} "
              f"= {mc['rate_all']:.6f}")
        if mc["valid_binding_only"] > 0:
            cond_rate = mc["valid_all"] / mc["valid_binding_only"]
            print(f"    P[D22 ok | D11 ok]: {cond_rate:.4f}")
            print(f"    P[D22 fail | D11 ok]: {mc['rate_d22_given_d11']:.4f}")
        print(f"    Time: {elapsed:.1f}s")

        # Compare QR-balanced vs random D11
        print(f"\n  D11 QUALITY COMPARISON:")
        rng = np.random.default_rng(99)
        d11_size = (p + 1) // 2
        QR = set()
        for x in range(1, p):
            QR.add((x * x) % p)

        # Sample random D11 and analyze
        good_count = 0  # QR-balanced
        bad_count = 0
        good_mc_rates = []
        bad_mc_rates = []

        for trial in range(50):
            rand_D11 = random_symmetric_d11(p, d11_size, rng)
            rand_info = analyze_d11(p, rand_D11)
            is_balanced = abs(rand_info["d11_qr"] - rand_info["d11_qnr"]) <= 2

            # Quick MC (5000 trials)
            mc_result = mc_full_check(p, rand_D11, rand_info["A_array"], num_mc=5000)

            if is_balanced:
                good_count += 1
                good_mc_rates.append(mc_result["rate_all"])
            else:
                bad_count += 1
                bad_mc_rates.append(mc_result["rate_all"])

        print(f"    QR-balanced D11 ({good_count}/50):")
        if good_mc_rates:
            print(f"      Hit rates: min={min(good_mc_rates):.6f}, "
                  f"max={max(good_mc_rates):.6f}, "
                  f"mean={sum(good_mc_rates) / len(good_mc_rates):.6f}")
        print(f"    Non-balanced D11 ({bad_count}/50):")
        if bad_mc_rates:
            print(f"      Hit rates: min={min(bad_mc_rates):.6f}, "
                  f"max={max(bad_mc_rates):.6f}, "
                  f"mean={sum(bad_mc_rates) / len(bad_mc_rates):.6f}")

    # Summary table
    print(f"\n{'=' * 80}")
    print("KEY INSIGHT SUMMARY")
    print(f"{'=' * 80}")
    print("""
The "loose" D22 constraints are NOT loose in practice because:
1. A(d) at D22 positions is systematically HIGHER than at D11 positions
2. The effective B(d) threshold at D22 can be TIGHTER than at D11
3. Constraints are coupled: pushing A(d)+B(d) down at D11 pushes it up at D22

For a proof, we need:
- D11 chosen to balance A(d) between D11 and D22 positions
- The QR-balanced D11 achieve this (all working D11 at p=11,19 are QR-balanced)
- Then a random D12 satisfies all constraints with probability ~1/poly(p)
""")


if __name__ == "__main__":
    main()
