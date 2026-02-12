#!/usr/bin/env python3
"""
For p = 11, 19, 23, 31: find ALL D11 with min possible max_d A(d),
then compute N(D11) for each to check if globally flat D11 have valid D12.
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


def run_prime(p):
    """Find globally flattest D11 and compute their N values."""
    n = (p + 1) // 2
    k = (p - 3) // 2
    thresh_binding = (p - 3) // 2
    thresh_loose = (p + 3) // 2
    E_A = (p + 1) / 4

    pairs = get_symmetric_pairs(p)
    total_d11 = comb(len(pairs), n // 2)
    total_d12 = comb(p - 1, k)

    print(f"\n{'=' * 70}")
    print(f"p = {p}, n = {n}, E[A] = {E_A}, |D12| = {total_d12}")
    print(f"Total symmetric D11 = {total_d11}")

    t0 = time.time()

    # Step 1: Find all D11 and their max_A_all
    print(f"  Computing max_A for all D11...")
    d11_data = []
    for D11 in enumerate_symmetric_D11(p, n):
        indicator = np.zeros(p, dtype=np.float64)
        for j in D11:
            indicator[j] = 1.0
        A = np.round(np.fft.ifft(np.abs(np.fft.fft(indicator)) ** 2).real).astype(int)
        max_A_all = int(max(A[d] for d in range(1, p)))
        max_A_D11 = int(max(A[d] for d in D11))
        D22 = set(range(1, p)) - D11
        max_A_D22 = int(max(A[d] for d in D22))
        d11_data.append({
            "D11": D11,
            "max_A_all": max_A_all,
            "max_A_D11": max_A_D11,
            "max_A_D22": max_A_D22,
            "A": A,
        })

    # Find minimum max_A
    min_max = min(d["max_A_all"] for d in d11_data)
    flattest = [d for d in d11_data if d["max_A_all"] == min_max]
    # Also get those with max_A_all = min_max + 1
    near_flat = [d for d in d11_data if d["max_A_all"] == min_max + 1]

    print(f"  Min max_A_all = {min_max} (gap = {min_max - E_A})")
    print(f"  #{len(flattest)} D11 achieve this minimum")
    print(f"  #{len(near_flat)} D11 have max_A_all = {min_max + 1}")

    # Step 2: Precompute all D12 autocorrelations
    print(f"  Precomputing {total_d12} D12 autocorrelations...")
    d12_matrix = np.zeros((total_d12, p), dtype=np.float64)
    d12_matrix[:, 0] = 1.0
    for i, chosen in enumerate(combinations(range(1, p), k)):
        for j in chosen:
            d12_matrix[i, j] = 1.0

    fft_vals = np.fft.fft(d12_matrix, axis=1)
    B_matrix = np.round(np.fft.ifft(np.abs(fft_vals) ** 2, axis=1).real).astype(np.int32)

    # Step 3: Compute N for flattest D11
    print(f"  Computing N for {len(flattest)} globally flattest D11...")
    for d in flattest:
        D11 = d["D11"]
        D22 = set(range(1, p)) - D11
        A = d["A"]

        d11_indices = np.array(sorted(D11), dtype=np.int32)
        d22_indices = np.array(sorted(D22), dtype=np.int32)
        A_d11 = np.array([int(A[i]) for i in d11_indices], dtype=np.int32)
        A_d22 = np.array([int(A[i]) for i in d22_indices], dtype=np.int32)

        F_binding = A_d11[np.newaxis, :] + B_matrix[:, d11_indices]
        valid_binding = np.all(F_binding <= thresh_binding, axis=1)

        F_loose = A_d22[np.newaxis, :] + B_matrix[:, d22_indices]
        valid_loose = np.all(F_loose <= thresh_loose, axis=1)

        valid_mask = valid_binding & valid_loose
        N = int(valid_mask.sum())
        d["N"] = N

    # Step 4: Also compute N for near-flat D11
    print(f"  Computing N for {len(near_flat)} near-flat D11...")
    for d in near_flat:
        D11 = d["D11"]
        D22 = set(range(1, p)) - D11
        A = d["A"]

        d11_indices = np.array(sorted(D11), dtype=np.int32)
        d22_indices = np.array(sorted(D22), dtype=np.int32)
        A_d11 = np.array([int(A[i]) for i in d11_indices], dtype=np.int32)
        A_d22 = np.array([int(A[i]) for i in d22_indices], dtype=np.int32)

        F_binding = A_d11[np.newaxis, :] + B_matrix[:, d11_indices]
        valid_binding = np.all(F_binding <= thresh_binding, axis=1)

        F_loose = A_d22[np.newaxis, :] + B_matrix[:, d22_indices]
        valid_loose = np.all(F_loose <= thresh_loose, axis=1)

        valid_mask = valid_binding & valid_loose
        N = int(valid_mask.sum())
        d["N"] = N

    elapsed = time.time() - t0

    # Report
    print(f"\n  Results (elapsed: {elapsed:.1f}s):")
    print(f"\n  --- Globally flattest (max_A_all = {min_max}, gap = {min_max - E_A}) ---")

    N_values_flat = defaultdict(int)
    for d in flattest:
        N_values_flat[d["N"]] += 1

    print(f"  N-value distribution:")
    for N_val in sorted(N_values_flat.keys()):
        print(f"    N = {N_val}: {N_values_flat[N_val]} D11")

    working_flat = sum(1 for d in flattest if d["N"] > 0)
    print(f"  Working: {working_flat}/{len(flattest)}")

    # Show a few examples
    for d in flattest[:5]:
        D11_sorted = sorted(d["D11"])
        A_profile = sorted(int(d["A"][i]) for i in range(1, p))
        print(f"    D11={D11_sorted[:6]}..., max_A_D11={d['max_A_D11']}, max_A_D22={d['max_A_D22']}, N={d['N']}")

    # Near-flat results
    print(f"\n  --- Near-flat (max_A_all = {min_max + 1}, gap = {min_max + 1 - E_A}) ---")

    N_values_near = defaultdict(int)
    for d in near_flat:
        N_values_near[d["N"]] += 1

    print(f"  N-value distribution:")
    for N_val in sorted(N_values_near.keys()):
        print(f"    N = {N_val}: {N_values_near[N_val]} D11")

    working_near = sum(1 for d in near_flat if d["N"] > 0)
    print(f"  Working: {working_near}/{len(near_flat)}")

    # Comparison with A-flat (D11-only flatness)
    a_flat_count = sum(1 for d in d11_data if d["max_A_D11"] <= floor(E_A))
    a_flat_working = sum(1 for d in d11_data if d["max_A_D11"] <= floor(E_A) and d.get("N") is not None and d["N"] > 0)

    # We need to compute N for all D11 to get this comparison correctly
    # For now, just note the counts
    print(f"\n  Comparison:")
    print(f"    Globally flat (gap={min_max - E_A}): {len(flattest)} D11, {working_flat} working")
    print(f"    Near-flat (gap={min_max + 1 - E_A}): {len(near_flat)} D11, {working_near} working")

    return {
        "p": p,
        "min_max": min_max,
        "gap": min_max - E_A,
        "num_flattest": len(flattest),
        "working_flattest": working_flat,
        "N_dist_flat": dict(N_values_flat),
        "num_near": len(near_flat),
        "working_near": working_near,
        "N_dist_near": dict(N_values_near),
    }


def main():
    print("=" * 70)
    print("DO GLOBALLY FLAT D11 HAVE VALID D12?")
    print("=" * 70)

    results = []
    for p in [11, 19, 23]:
        r = run_prime(p)
        results.append(r)

    # Summary
    print(f"\n{'=' * 70}")
    print("SUMMARY: Global flatness → validity?")
    print(f"{'=' * 70}")
    print(f"{'p':>5} {'gap':>5} {'#flat':>7} {'#working':>9} {'fraction':>10} {'N values (flat)':>30}")
    print("-" * 70)
    for r in results:
        frac = r["working_flattest"] / r["num_flattest"] if r["num_flattest"] > 0 else 0
        n_vals = sorted(r["N_dist_flat"].keys())
        print(f"{r['p']:>5} {r['gap']:>5.0f} {r['num_flattest']:>7} {r['working_flattest']:>9} {frac:>10.2%} {str(n_vals):>30}")


if __name__ == "__main__":
    main()
