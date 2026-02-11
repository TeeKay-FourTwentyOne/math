#!/usr/bin/env python3
"""Sampling verification for p=31,43,47,59: estimate fraction of valid D12 per D11.

For each prime p:
1. Generate 100+ random symmetric D11 of size n=(p+1)/2
2. For each D11, sample random D12 sets in batches (batch FFT for speed)
3. Also test known-good D11 from SA solutions
4. Report estimated E[valid D12] per D11 and overall statistics

Optimized: batch FFT computation of B(d) for multiple D12 simultaneously.
"""

import numpy as np
import json
import time
import os
import sys
from math import comb
from collections import defaultdict


# Known solutions (|D11| = n = (p+1)/2 convention)
KNOWN_SOLUTIONS = {
    31: {
        "D11": {6, 7, 8, 10, 11, 12, 14, 15, 16, 17, 19, 20, 21, 23, 24, 25},
        "D12": {0, 1, 2, 3, 8, 11, 12, 13, 15, 18, 20, 21, 24, 27, 29},
    },
    43: {
        "D11": {1, 2, 5, 10, 11, 13, 16, 17, 18, 19, 20, 23, 24, 25, 26, 27,
                30, 32, 33, 38, 41, 42},
        "D12": {0, 2, 5, 6, 8, 11, 15, 16, 20, 24, 25, 27, 28, 31, 32, 34,
                35, 36, 37, 39, 41},
    },
    47: None,
    59: None,
}


def load_registry():
    """Load solutions from solutions_registry.json."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(script_dir, "solutions_registry.json")
    if not os.path.exists(path):
        return
    with open(path) as f:
        registry = json.load(f)
    for sol in registry["solutions"]:
        p = sol["m"]
        if p not in {47, 59}:
            continue
        D11_raw = set(sol["D11"])
        D12_raw = set(sol["D12"])
        target = (p + 1) // 2
        if len(D11_raw) != target:
            D22 = set(range(1, p)) - D11_raw
            D12T = {(-x) % p for x in D12_raw}
            D11_raw, D12_raw = D22, D12T
        KNOWN_SOLUTIONS[p] = {"D11": D11_raw, "D12": D12_raw}


def batch_autocorrelation(indicator_matrix):
    """Compute Delta(S,S,d) for all sets via batch FFT.
    indicator_matrix: (batch_size, p) float
    Returns: (batch_size, p) int
    """
    fft_vals = np.fft.fft(indicator_matrix, axis=1)
    autocorr = np.fft.ifft(np.abs(fft_vals) ** 2, axis=1).real
    return np.round(autocorr).astype(np.int32)


def autocorrelation_single(S, p):
    """Compute A(d) for a single set."""
    indicator = np.zeros(p, dtype=np.float64)
    for j in S:
        indicator[j] = 1.0
    fft_val = np.fft.fft(indicator)
    autocorr = np.fft.ifft(np.abs(fft_val) ** 2).real
    return np.round(autocorr).astype(np.int32)


def random_symmetric_D11(p, rng):
    """Generate a random symmetric D11 of size (p+1)/2."""
    n = (p + 1) // 2
    num_pairs = n // 2
    pairs = [(d, (p - d) % p) for d in range(1, (p + 1) // 2)]
    chosen = rng.choice(len(pairs), size=num_pairs, replace=False)
    D11 = set()
    for i in chosen:
        D11.add(pairs[i][0])
        D11.add(pairs[i][1])
    return D11


def generate_random_d12_batch(p, k, batch_size, rng):
    """Generate a batch of random D12 = {0} ∪ S indicator matrices.

    Returns (batch_size, p) float indicator matrix.
    """
    matrix = np.zeros((batch_size, p), dtype=np.float64)
    matrix[:, 0] = 1.0  # 0 always in D12
    candidates = np.arange(1, p)
    for i in range(batch_size):
        S = rng.choice(candidates, size=k, replace=False)
        matrix[i, S] = 1.0
    return matrix


def sample_d12_for_d11(D11, p, A, num_samples, batch_size, rng):
    """Sample D12 in batches, count valid ones.

    Returns dict with stats.
    """
    n = (p + 1) // 2
    k = (p - 3) // 2
    threshold_binding = (p - 3) // 2
    threshold_loose = (p + 3) // 2

    d11_arr = np.array(sorted(D11), dtype=np.int32)
    d22_arr = np.array(sorted(set(range(1, p)) - D11), dtype=np.int32)

    A_at_d11 = A[d11_arr]
    A_at_d22 = A[d22_arr]

    num_valid = 0
    marginal_hits = np.zeros(len(d11_arr), dtype=np.int64)
    total_sampled = 0

    while total_sampled < num_samples:
        bs = min(batch_size, num_samples - total_sampled)
        d12_batch = generate_random_d12_batch(p, k, bs, rng)
        B_batch = batch_autocorrelation(d12_batch)  # (bs, p)

        # Binding constraints: A(d) + B(d) <= threshold_binding for d in D11
        F_binding = A_at_d11[np.newaxis, :] + B_batch[:, d11_arr]
        binding_ok = F_binding <= threshold_binding  # (bs, |D11|)

        # Per-constraint marginal counts
        marginal_hits += binding_ok.sum(axis=0)

        # Joint binding check
        all_binding_ok = np.all(binding_ok, axis=1)

        # Loose constraints: A(d) + B(d) <= threshold_loose for d in D22
        if d22_arr.size > 0:
            F_loose = A_at_d22[np.newaxis, :] + B_batch[:, d22_arr]
            all_loose_ok = np.all(F_loose <= threshold_loose, axis=1)
            valid_mask = all_binding_ok & all_loose_ok
        else:
            valid_mask = all_binding_ok

        num_valid += int(valid_mask.sum())
        total_sampled += bs

    rate = num_valid / total_sampled
    marginal_rates = marginal_hits / total_sampled
    prod_marginals = float(np.prod(marginal_rates)) if np.all(marginal_rates > 0) else 0.0
    ratio = rate / prod_marginals if prod_marginals > 0 else (float('inf') if num_valid > 0 else 0.0)

    num_d12_total = comb(p - 1, k)
    E_valid = rate * num_d12_total

    return {
        "num_valid": num_valid,
        "total_sampled": total_sampled,
        "rate": rate,
        "E_valid": E_valid,
        "ratio": ratio,
        "marginal_rate_min": float(marginal_rates.min()),
        "marginal_rate_max": float(marginal_rates.max()),
        "marginal_rate_mean": float(marginal_rates.mean()),
    }


def run_prime(p, num_random_d11=100, num_d12_samples=500000, batch_size=10000):
    """Run sampling verification for a single prime."""
    n = (p + 1) // 2
    k = (p - 3) // 2
    num_d12_total = comb(p - 1, k)

    print(f"\n{'=' * 80}")
    print(f"PRIME p = {p}, n = {n}")
    print(f"|D11| = {n}, |D12| = {n-1}, k = {k}")
    print(f"C({p-1},{k}) = {num_d12_total:,.0f}")
    print(f"Sampling {num_d12_samples:,} D12 per D11, {num_random_d11} random D11")
    print(f"Batch size: {batch_size}")
    print(f"{'=' * 80}")

    rng = np.random.default_rng(42)
    t0 = time.time()

    all_results = []

    # Test known D11 first
    known = KNOWN_SOLUTIONS.get(p)
    if known:
        D11 = known["D11"]
        A = autocorrelation_single(D11, p)
        max_A = max(int(A[d]) for d in D11)
        print(f"\nKnown D11: max_A = {max_A}")

        stats = sample_d12_for_d11(D11, p, A, num_d12_samples, batch_size, rng)
        print(f"  Valid: {stats['num_valid']}/{stats['total_sampled']} = {stats['rate']:.6f}")
        print(f"  E[valid D12]: {stats['E_valid']:.1f}")
        print(f"  Marginals: min={stats['marginal_rate_min']:.4f}, "
              f"max={stats['marginal_rate_max']:.4f}")
        print(f"  Ratio: {stats['ratio']:.4f}")

        all_results.append({
            "type": "known",
            "D11": sorted(D11),
            "max_A": max_A,
            **stats,
        })

    # Random D11
    print(f"\nSampling random D11...")
    max_A_dist = defaultdict(lambda: {"total": 0, "working": 0})
    working_count = 0

    for i in range(num_random_d11):
        D11 = random_symmetric_D11(p, rng)
        A = autocorrelation_single(D11, p)
        max_A = max(int(A[d]) for d in D11)
        max_A_dist[max_A]["total"] += 1

        stats = sample_d12_for_d11(D11, p, A, num_d12_samples, batch_size, rng)

        if stats["num_valid"] > 0:
            working_count += 1
            max_A_dist[max_A]["working"] += 1

        all_results.append({
            "type": "random",
            "D11": sorted(D11),
            "max_A": max_A,
            **stats,
        })

        if (i + 1) % 10 == 0:
            elapsed = time.time() - t0
            print(f"  {i+1}/{num_random_d11}: "
                  f"working={working_count}, "
                  f"this: max_A={max_A}, valid={stats['num_valid']}, "
                  f"{elapsed:.1f}s")

    elapsed = time.time() - t0

    # Summary
    print(f"\nCompleted in {elapsed:.1f}s")
    print(f"Working D11: {working_count}/{num_random_d11} "
          f"({100*working_count/num_random_d11:.1f}%)")

    print(f"\nmax_A distribution:")
    for mA in sorted(max_A_dist):
        info = max_A_dist[mA]
        pct = 100 * info["working"] / info["total"] if info["total"] > 0 else 0
        print(f"  max_A={mA}: {info['total']} D11, {info['working']} working ({pct:.1f}%)")

    working_results = [r for r in all_results if r["num_valid"] > 0]
    if working_results:
        rates = [r["rate"] for r in working_results]
        E_valids = [r["E_valid"] for r in working_results]
        ratios = [r["ratio"] for r in working_results
                  if isinstance(r["ratio"], (int, float)) and r["ratio"] < float('inf')]
        print(f"\nWorking D11 statistics:")
        print(f"  Rate: min={min(rates):.6f}, max={max(rates):.6f}")
        print(f"  E[valid D12]: min={min(E_valids):.1f}, max={max(E_valids):.1f}")
        if ratios:
            print(f"  Ratio: min={min(ratios):.6f}, max={max(ratios):.6f}")
        print(f"  max_A values: {sorted(set(r['max_A'] for r in working_results))}")

    return {
        "p": p,
        "n": n,
        "num_random_d11": num_random_d11,
        "num_d12_samples": num_d12_samples,
        "working_count": working_count,
        "elapsed": elapsed,
        "results": all_results,
        "max_A_dist": {str(k): v for k, v in max_A_dist.items()},
    }


def main():
    print("=" * 80)
    print("SAMPLING VERIFICATION FOR LARGER PRIMES")
    print("=" * 80)

    load_registry()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.join(script_dir, "..", "enumeration_data")
    os.makedirs(output_dir, exist_ok=True)

    all_summaries = {}

    for p in [31, 43, 47, 59]:
        # Adjust parameters: fewer samples for larger p (batch FFT memory)
        if p <= 31:
            num_d11, num_d12, bs = 100, 500000, 10000
        elif p <= 47:
            num_d11, num_d12, bs = 60, 200000, 5000
        else:
            num_d11, num_d12, bs = 40, 100000, 2000

        data = run_prime(p, num_random_d11=num_d11, num_d12_samples=num_d12, batch_size=bs)
        all_summaries[str(p)] = {
            "p": p, "n": data["n"],
            "working": data["working_count"],
            "total_d11": data["num_random_d11"],
            "fraction_working": data["working_count"] / data["num_random_d11"],
            "elapsed": data["elapsed"],
        }

        # Save results
        outpath = os.path.join(output_dir, f"sampling_p{p}.json")
        working_only = [r for r in data["results"] if r["num_valid"] > 0]
        save_data = {
            "summary": all_summaries[str(p)],
            "known_result": next((r for r in data["results"] if r["type"] == "known"), None),
            "working_results": working_only[:30],
            "max_A_dist": data["max_A_dist"],
        }
        with open(outpath, "w") as f:
            json.dump(save_data, f, indent=2, default=str)
        print(f"\nData saved to {outpath}")

    # Cross-prime summary
    print(f"\n{'=' * 80}")
    print("CROSS-PRIME SUMMARY")
    print(f"{'=' * 80}")
    print(f"{'p':>4s} {'n':>4s} {'working':>8s} {'total':>8s} {'frac':>8s} {'time':>8s}")
    print(f"  {'-'*42}")
    for p_str, s in all_summaries.items():
        print(f"  {s['p']:4d} {s['n']:4d} {s['working']:8d} {s['total_d11']:8d} "
              f"{s['fraction_working']:8.4f} {s['elapsed']:8.1f}s")


if __name__ == "__main__":
    main()
