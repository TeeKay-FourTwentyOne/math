#!/usr/bin/env python3
"""Cross-prime analysis of sampling results.

Computes key metrics across all tested primes to identify growth trends.
"""

import json
import os
import sys
from math import log, log2, comb, sqrt, exp
import numpy as np


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(script_dir, "..", "enumeration_data")

    print("=" * 100)
    print("CROSS-PRIME ANALYSIS")
    print("=" * 100)

    # Exhaustive data for small primes
    exhaustive = {}
    for p in [11, 19, 23]:
        path = os.path.join(data_dir, f"enumeration_p{p}.json")
        if os.path.exists(path):
            with open(path) as f:
                data = json.load(f)
            n = (p + 1) // 2
            k = (p - 3) // 2
            total_d12 = comb(p - 1, k)
            total_d11 = data["summary"]["num_d11"]
            total_valid = sum(d["num_valid_d12"] for d in data["per_d11"])
            working_d11 = sum(1 for d in data["per_d11"] if d["num_valid_d12"] > 0)

            # Get best D11 stats
            best = max(data["per_d11"], key=lambda d: d["num_valid_d12"])

            exhaustive[p] = {
                "n": n,
                "total_d11": total_d11,
                "working_d11": working_d11,
                "total_valid_pairs": total_valid,
                "total_d12": total_d12,
                "best_d11_valid": best["num_valid_d12"],
                "log2_total_d12": log2(total_d12),
            }

    # Sampling data for larger primes
    sampling = {}
    for p in [31, 43, 47, 59]:
        path = os.path.join(data_dir, f"sampling_aflat_p{p}.json")
        if not os.path.exists(path):
            continue
        with open(path) as f:
            data = json.load(f)

        n = data["n"]
        k = (p - 3) // 2
        total_d12 = comb(p - 1, k)

        # Find best Gaussian estimate
        best_gauss = None
        best_binding_rate = 0
        known_result = None

        for r in data["results"]:
            if r.get("gaussian") and isinstance(r["gaussian"].get("log_E_valid"), (int, float)):
                le = r["gaussian"]["log_E_valid"]
                if best_gauss is None or le > best_gauss:
                    best_gauss = le
            if r["binding_rate"] > best_binding_rate:
                best_binding_rate = r["binding_rate"]
            if r["type"] == "known":
                known_result = r

        sampling[p] = {
            "n": n,
            "total_d12": total_d12,
            "log2_total_d12": log2(total_d12),
            "best_gauss_log_E": best_gauss,
            "best_binding_rate": best_binding_rate,
            "known_result": known_result,
        }

    # Print exhaustive summary
    print("\n--- EXHAUSTIVE RESULTS (p=11,19,23) ---")
    print(f"{'p':>4s} {'n':>4s} {'#D11':>8s} {'working':>8s} {'%work':>8s} "
          f"{'#valid':>10s} {'log2(C)':>8s} {'best/D11':>10s}")
    for p, e in exhaustive.items():
        pct = 100 * e["working_d11"] / e["total_d11"]
        print(f"  {p:4d} {e['n']:4d} {e['total_d11']:8d} {e['working_d11']:8d} "
              f"{pct:7.1f}% {e['total_valid_pairs']:10d} "
              f"{e['log2_total_d12']:8.1f} {e['best_d11_valid']:10d}")

    # Print sampling summary
    print("\n--- SAMPLING RESULTS (p=31,43,47,59) ---")
    print(f"{'p':>4s} {'n':>4s} {'log2(C)':>8s} {'known_bind':>12s} "
          f"{'known_ratio':>12s} {'gauss_logE':>12s} {'gauss_E':>15s}")
    for p, s in sampling.items():
        kr = s["known_result"]
        if kr:
            kb = f"{kr['binding_rate']:.6f}"
            ratio = f"{kr['ratio_binding']:.2f}" if isinstance(kr.get("ratio_binding"), (int, float)) else "N/A"
        else:
            kb = "N/A"
            ratio = "N/A"
        ge = s["best_gauss_log_E"]
        ge_str = f"{ge:.2f}" if ge else "N/A"
        ev = f"{exp(ge):.4e}" if ge and ge > -500 else "N/A"
        print(f"  {p:4d} {s['n']:4d} {s['log2_total_d12']:8.1f} "
              f"{kb:>12s} {ratio:>12s} {ge_str:>12s} {ev:>15s}")

    # Growth analysis
    print("\n--- GROWTH ANALYSIS ---")
    print("\nKey metrics vs p:")
    print(f"{'p':>4s} {'n':>4s} {'log2(C)':>8s} {'budget':>8s} {'cost_est':>10s} {'margin':>10s}")

    for p in sorted(list(exhaustive.keys()) + list(sampling.keys())):
        n = (p + 1) // 2
        k = (p - 3) // 2
        log2C = log2(comb(p - 1, k))

        # Exact moments
        mu_B = (p - 3) * (p - 1) / (4 * (p - 2))
        var_B = (p - 3) * (p + 1) / (16 * (p - 2))
        sigma_B = sqrt(var_B)
        threshold = (p - 3) / 2

        # For "ideal" D11 with A(d) = (n-1)*(n-2)/(p-1) for all d (perfect flatness)
        # A(d) ~ n^2 / (4p) ~ p/4
        ideal_A = (n * (n - 1)) / (p - 1)  # = |D11|*(|D11|-1)/(p-1)
        slack = threshold - ideal_A - mu_B
        z = slack / sigma_B

        # Cost per constraint ~ -log2(Phi(z))
        from scipy.stats import norm
        cost_per = -norm.logcdf(z) / log(2) if z < 5 else 0
        total_cost = cost_per * n  # n binding constraints
        margin = log2C - total_cost

        # For known solution D11 (if available)
        if p in sampling and sampling[p]["best_gauss_log_E"]:
            actual_margin = sampling[p]["best_gauss_log_E"] / log(2)
        elif p in exhaustive:
            tv = exhaustive[p]["total_valid_pairs"]
            actual_margin = log2(tv) if tv > 0 else float('-inf')
        else:
            actual_margin = None

        am_str = f"{actual_margin:.1f}" if actual_margin else "N/A"

        print(f"  {p:4d} {n:4d} {log2C:8.1f} {log2C:8.1f} {total_cost:10.1f} {margin:10.1f}  (actual: {am_str})")

    # Binding ratio growth
    print("\n--- BINDING RATIO (positive association) ---")
    print("For known-good D11:")
    for p in [31, 43, 47, 59]:
        if p in sampling and sampling[p]["known_result"]:
            kr = sampling[p]["known_result"]
            print(f"  p={p}: binding ratio = {kr['ratio_binding']:.2f}")

    # What fraction of D12 space is binding-valid for best D11?
    print("\n--- BINDING-ONLY PASS RATE ---")
    for p in [31, 43, 47, 59]:
        if p in sampling and sampling[p]["known_result"]:
            kr = sampling[p]["known_result"]
            rate = kr["binding_rate"]
            total = comb(p - 1, (p - 3) // 2)
            E_binding = rate * total
            print(f"  p={p}: rate={rate:.6f}, E[binding-valid] = {E_binding:.0f} "
                  f"(log2 = {log2(E_binding):.1f})" if E_binding > 0 else
                  f"  p={p}: rate=0")

    # Max A(d) at D22 positions for known D11
    print("\n--- D22 BOTTLENECK ---")
    print("For known-good D11:")
    for p in [31, 43, 47, 59]:
        if p in sampling and sampling[p]["known_result"]:
            kr = sampling[p]["known_result"]
            threshold_loose = (p + 3) // 2
            max_A_d22 = kr.get("max_A_d22", "?")
            mu_B = (p - 3) * (p - 1) / (4 * (p - 2))
            slack = threshold_loose - max_A_d22 - mu_B if isinstance(max_A_d22, (int, float)) else None
            slack_str = f"{slack:.2f}" if slack else "?"
            print(f"  p={p}: max_A(D22)={max_A_d22}, threshold={threshold_loose}, "
                  f"E[B]={mu_B:.2f}, slack={slack_str}")


if __name__ == "__main__":
    main()
