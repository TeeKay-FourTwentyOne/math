#!/usr/bin/env python3
"""Global A-flatness and full-constraint marginal analysis for p=11,19,23.

For each D11 from exhaustive enumeration:
1. Compute A(d) for ALL d=1..p-1 (not just D11)
2. Compute D22 per-constraint marginals
3. Compute product of ALL marginals (D11 + D22)
4. Compare joint prob vs product
5. Compute log2(E[valid]) headroom

Output: detailed tables for the proof fix.
"""

import numpy as np
import json
import os
import sys
from math import comb, log2


def autocorrelation_full(D11, p):
    """Compute A(d) = Delta(D11, D11, d) for ALL d=0..p-1."""
    indicator = np.zeros(p, dtype=np.float64)
    for j in D11:
        indicator[j] = 1.0
    fft_val = np.fft.fft(indicator)
    autocorr = np.fft.ifft(np.abs(fft_val) ** 2).real
    return np.round(autocorr).astype(int)


def analyze_prime(p, data):
    """Full analysis for one prime."""
    summary = data["summary"]
    per_d11 = data["per_d11"]
    n = summary["n"]
    k = summary["k"]
    threshold_binding = summary["threshold_binding"]
    threshold_loose = summary["threshold_loose"]
    num_d12 = summary["num_d12"]

    print(f"\n{'='*80}")
    print(f"p = {p}, n = {n}, k = {k}")
    print(f"Binding threshold: A(d)+B(d) <= {threshold_binding}  (d in D11)")
    print(f"Loose threshold:   A(d)+B(d) <= {threshold_loose}  (d in D22)")
    print(f"Total D12: {num_d12}")
    print(f"log2(C({p-1},{k})) = {log2(comb(p-1,k)):.2f}")
    print(f"{'='*80}")

    # Precompute all D12 autocorrelations for exact marginal computation
    all_d12_B = precompute_all_B(p, k)  # (num_d12, p) array

    results = []

    for entry in per_d11:
        D11 = set(entry["D11"])
        D22 = set(range(1, p)) - D11
        num_valid = entry["num_valid_d12"]
        joint_prob = num_valid / num_d12

        # Full autocorrelation at ALL positions
        A = autocorrelation_full(D11, p)

        max_A_d11 = max(int(A[d]) for d in D11)
        max_A_d22 = max(int(A[d]) for d in D22) if D22 else 0
        max_A_all = max(int(A[d]) for d in range(1, p))

        # Compute D22 per-constraint marginals from precomputed B values
        d22_marginal_counts = {}
        for d in sorted(D22):
            Ad = int(A[d])
            T_d22 = threshold_loose - Ad
            # Count D12 with B(d) <= T_d22
            count = int(np.sum(all_d12_B[:, d] <= T_d22))
            d22_marginal_counts[d] = count

        # D11 marginals (already in data, but let's recompute for consistency)
        d11_marginal_counts = {}
        for d in sorted(D11):
            Ad = int(A[d])
            T_d11 = threshold_binding - Ad
            count = int(np.sum(all_d12_B[:, d] <= T_d11))
            d11_marginal_counts[d] = count

        # Products
        prod_d11 = 1.0
        log2_prod_d11 = 0.0
        for d in sorted(D11):
            p_marg = d11_marginal_counts[d] / num_d12
            if p_marg > 0:
                prod_d11 *= p_marg
                log2_prod_d11 += log2(p_marg)
            else:
                prod_d11 = 0.0
                log2_prod_d11 = float('-inf')
                break

        prod_d22 = 1.0
        log2_prod_d22 = 0.0
        for d in sorted(D22):
            p_marg = d22_marginal_counts[d] / num_d12
            if p_marg > 0:
                prod_d22 *= p_marg
                log2_prod_d22 += log2(p_marg)
            else:
                prod_d22 = 0.0
                log2_prod_d22 = float('-inf')
                break

        prod_all = prod_d11 * prod_d22
        log2_prod_all = log2_prod_d11 + log2_prod_d22

        # Ratios
        ratio_d11 = joint_prob / prod_d11 if prod_d11 > 0 else (float('inf') if num_valid > 0 else 0.0)
        ratio_all = joint_prob / prod_all if prod_all > 0 else (float('inf') if num_valid > 0 else 0.0)

        # log2(E[valid])
        log2_E = log2(num_valid) if num_valid > 0 else float('-inf')

        # Number of constraints
        r1 = len(D11)  # binding constraints
        r2 = len(D22)  # loose constraints

        r = {
            "D11": sorted(D11),
            "num_valid": num_valid,
            "joint_prob": joint_prob,
            "max_A_d11": max_A_d11,
            "max_A_d22": max_A_d22,
            "max_A_all": max_A_all,
            "A_at_d11": {d: int(A[d]) for d in sorted(D11)},
            "A_at_d22": {d: int(A[d]) for d in sorted(D22)},
            "r1": r1,
            "r2": r2,
            "r_total": r1 + r2,
            "log2_prod_d11": log2_prod_d11,
            "log2_prod_d22": log2_prod_d22,
            "log2_prod_all": log2_prod_all,
            "ratio_d11": ratio_d11,
            "ratio_all": ratio_all,
            "log2_E": log2_E,
            "d22_marginal_counts": d22_marginal_counts,
            "d22_marginal_min": min(d22_marginal_counts[d] / num_d12 for d in D22) if D22 else 1.0,
            "d22_marginal_max": max(d22_marginal_counts[d] / num_d12 for d in D22) if D22 else 1.0,
        }
        results.append(r)

    return results


def precompute_all_B(p, k):
    """Precompute B(d) for all D12 = {0} ∪ (k elements from {1..p-1}).

    Returns (num_d12, p) integer array.
    """
    from itertools import combinations
    candidates = list(range(1, p))
    d12_list = []
    for combo in combinations(candidates, k):
        D12 = {0} | set(combo)
        d12_list.append(D12)

    num_d12 = len(d12_list)
    indicator_matrix = np.zeros((num_d12, p), dtype=np.float64)
    for i, D12 in enumerate(d12_list):
        for j in D12:
            indicator_matrix[i, j] = 1.0

    # Batch FFT
    fft_vals = np.fft.fft(indicator_matrix, axis=1)
    autocorr = np.fft.ifft(np.abs(fft_vals) ** 2, axis=1).real
    return np.round(autocorr).astype(np.int32)


def print_table(p, results, num_d12):
    """Print the key table."""
    working = [r for r in results if r["num_valid"] > 0]
    non_working = [r for r in results if r["num_valid"] == 0]

    print(f"\n--- WORKING D11 ({len(working)}) ---")
    print(f"{'D11':>30s} | {'#val':>5s} | {'maxA_D11':>8s} | {'maxA_D22':>8s} | {'maxA_all':>8s} | "
          f"{'r1+r2':>5s} | {'lg2_prod_D11':>12s} | {'lg2_prod_all':>12s} | "
          f"{'ratio_D11':>10s} | {'ratio_all':>10s} | {'lg2_E':>8s}")
    print("-" * 160)
    for r in sorted(working, key=lambda x: -x["num_valid"]):
        d11_str = str(r["D11"])
        if len(d11_str) > 30:
            d11_str = d11_str[:27] + "..."
        print(f"{d11_str:>30s} | {r['num_valid']:5d} | {r['max_A_d11']:8d} | {r['max_A_d22']:8d} | {r['max_A_all']:8d} | "
              f"{r['r_total']:5d} | {r['log2_prod_d11']:12.2f} | {r['log2_prod_all']:12.2f} | "
              f"{r['ratio_d11']:10.4f} | {r['ratio_all']:10.4f} | {r['log2_E']:8.2f}")

    print(f"\n--- NON-WORKING D11 (top 10 by max_A_all) ---")
    non_sorted = sorted(non_working, key=lambda x: x["max_A_all"])
    for r in non_sorted[:10]:
        d11_str = str(r["D11"])
        if len(d11_str) > 30:
            d11_str = d11_str[:27] + "..."
        print(f"{d11_str:>30s} | {r['num_valid']:5d} | {r['max_A_d11']:8d} | {r['max_A_d22']:8d} | {r['max_A_all']:8d} | "
              f"{r['r_total']:5d} | {r['log2_prod_d11']:12.2f} | {r['log2_prod_all']:12.2f} | "
              f"{'N/A':>10s} | {'N/A':>10s} | {'N/A':>8s}")

    # D22 marginal analysis
    print(f"\n--- D22 MARGINAL ANALYSIS ---")
    print(f"For WORKING D11:")
    for r in sorted(working, key=lambda x: -x["num_valid"]):
        d11_str = str(r["D11"])
        if len(d11_str) > 30:
            d11_str = d11_str[:27] + "..."
        d22_margs = {d: r["d22_marginal_counts"][d] / num_d12 for d in sorted(r["d22_marginal_counts"])}
        min_m = min(d22_margs.values()) if d22_margs else 1.0
        max_m = max(d22_margs.values()) if d22_margs else 1.0
        A_at_worst = None
        worst_d = None
        for d, m in d22_margs.items():
            if m == min_m:
                worst_d = d
                A_at_worst = r["A_at_d22"][d]
                break
        print(f"  {d11_str:>30s}: D22 marg range [{min_m:.4f}, {max_m:.4f}], "
              f"worst d={worst_d} A(d)={A_at_worst}, "
              f"lg2(prod_D22)={r['log2_prod_d22']:.2f}")

    # Correlation: max_A_d22 vs working
    print(f"\n--- max_A_D22 vs WORKING ---")
    from collections import Counter
    w_counts = Counter(r["max_A_d22"] for r in working)
    nw_counts = Counter(r["max_A_d22"] for r in non_working)
    all_vals = sorted(set(list(w_counts.keys()) + list(nw_counts.keys())))
    for v in all_vals:
        wc = w_counts.get(v, 0)
        nwc = nw_counts.get(v, 0)
        tot = wc + nwc
        pct = 100 * wc / tot if tot > 0 else 0
        print(f"  max_A_D22={v}: {wc} working, {nwc} non-working ({pct:.0f}% working)")

    # Headroom: log2(E[valid]) for working D11
    print(f"\n--- HEADROOM: log2(E[valid D12]) ---")
    log2_C = log2(comb(p - 1, (p - 3) // 2))
    for r in sorted(working, key=lambda x: -x["num_valid"]):
        # E[valid] = num_valid (exhaustive)
        # headroom = log2(E) = log2(num_valid)
        # Compare with (p+1)/2
        target = (p + 1) / 2
        print(f"  D11={r['D11']}: log2(E)={r['log2_E']:.2f}, target (p+1)/2={target:.1f}, "
              f"ratio log2(E)/target = {r['log2_E']/target:.3f}")


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(script_dir, "..", "enumeration_data")

    print("=" * 80)
    print("GLOBAL A-FLATNESS AND FULL-CONSTRAINT MARGINAL ANALYSIS")
    print("=" * 80)
    sys.stdout.flush()

    all_prime_results = {}

    for p in [11, 19, 23]:
        path = os.path.join(data_dir, f"enumeration_p{p}.json")
        with open(path) as f:
            data = json.load(f)

        results = analyze_prime(p, data)
        num_d12 = data["summary"]["num_d12"]
        print_table(p, results, num_d12)
        all_prime_results[p] = results
        sys.stdout.flush()

    # Cross-prime headroom summary
    print(f"\n{'='*80}")
    print("CROSS-PRIME HEADROOM SUMMARY")
    print(f"{'='*80}")
    print(f"{'p':>4s} {'n':>4s} {'#work':>6s} {'best_lg2E':>10s} {'(p+1)/2':>8s} {'lg2C':>8s} "
          f"{'best_ratio_all':>15s} {'best_ratio_D11':>15s}")
    for p, results in all_prime_results.items():
        n = (p + 1) // 2
        working = [r for r in results if r["num_valid"] > 0]
        if not working:
            continue
        best_E = max(r["log2_E"] for r in working)
        best_ratio_all = max(r["ratio_all"] for r in working
                             if isinstance(r["ratio_all"], (int, float)) and r["ratio_all"] < float('inf'))
        best_ratio_d11 = max(r["ratio_d11"] for r in working
                              if isinstance(r["ratio_d11"], (int, float)) and r["ratio_d11"] < float('inf'))
        log2_C = log2(comb(p - 1, (p - 3) // 2))
        print(f"  {p:4d} {n:4d} {len(working):6d} {best_E:10.2f} {(p+1)/2:8.1f} {log2_C:8.1f} "
              f"{best_ratio_all:15.4f} {best_ratio_d11:15.4f}")

    # Save detailed results
    outpath = os.path.join(data_dir, "global_aflat_analysis.json")
    save_data = {}
    for p, results in all_prime_results.items():
        save_data[str(p)] = []
        for r in results:
            r_copy = dict(r)
            # Convert sets to lists for JSON
            r_copy["d22_marginal_counts"] = {str(k): v for k, v in r_copy["d22_marginal_counts"].items()}
            r_copy["A_at_d11"] = {str(k): v for k, v in r_copy["A_at_d11"].items()}
            r_copy["A_at_d22"] = {str(k): v for k, v in r_copy["A_at_d22"].items()}
            save_data[str(p)].append(r_copy)

    with open(outpath, "w") as f:
        json.dump(save_data, f, indent=2, default=str)
    print(f"\nDetailed results saved to {outpath}")


if __name__ == "__main__":
    main()
