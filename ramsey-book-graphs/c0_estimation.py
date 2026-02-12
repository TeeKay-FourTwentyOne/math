#!/usr/bin/env python3
"""
Estimate c₀(D11) = Pr[D12 valid | D11] / ∏_d Pr[B(d) ≤ T(d)]

c₀ captures how much the Parseval hyperplane constraint (negative correlation
among B(d) values) affects joint validity relative to independent marginals.

For the proof:
  c₀ > 1 for working D11 ↔ negative correlation HELPS (dependence beneficial)
  c₀ ~ O(1) as p → ∞ ↔ second moment method works

Method:
  1. Estimate F(T) = Pr[B(1) ≤ T] by sampling random D12 (universal by Z_p symmetry)
  2. For each D11 orbit: compute T(d) = threshold - A(d) for all positions d
  3. c₀ = (N_estimate / total_D12) / ∏_d F(T(d))
"""

import numpy as np
import json
import time
import os
from itertools import combinations
from math import comb, sqrt, floor, log, exp
from collections import defaultdict


def get_symmetric_pairs(p):
    """Get all negation pairs {d, p-d} from {1,...,p-1}."""
    pairs = []
    seen = set()
    for d in range(1, p):
        if d not in seen:
            pairs.append((d, (-d) % p))
            seen.add(d)
            seen.add((-d) % p)
    return pairs


def enumerate_symmetric_D11(p, size):
    """Enumerate all symmetric subsets of {1,...,p-1} with given even size."""
    pairs = get_symmetric_pairs(p)
    num_pairs = size // 2
    for chosen_pairs in combinations(range(len(pairs)), num_pairs):
        D11 = set()
        for i in chosen_pairs:
            d, neg_d = pairs[i]
            D11.add(d)
            D11.add(neg_d)
        yield frozenset(D11)


def primitive_root(p):
    """Find a primitive root mod p."""
    for g in range(2, p):
        seen = set()
        val = 1
        for _ in range(p - 1):
            val = (val * g) % p
            seen.add(val)
        if len(seen) == p - 1:
            return g
    return None


def canonical_orbit_rep(D11_set, p):
    """Find lexicographically smallest orbit member under multiplication by Z_p*."""
    g = primitive_root(p)
    best = tuple(sorted(D11_set))
    multiplier = 1
    for _ in range((p - 1) // 2):
        multiplier = (multiplier * g) % p
        transformed = frozenset((d * multiplier) % p for d in D11_set)
        candidate = tuple(sorted(transformed))
        if candidate < best:
            best = candidate
    return best


def autocorrelation_single(D_set, p):
    """Compute autocorrelation A(d) for a single set via FFT."""
    indicator = np.zeros(p, dtype=np.float64)
    for j in D_set:
        indicator[j] = 1.0
    fft_val = np.fft.fft(indicator)
    autocorr = np.fft.ifft(np.abs(fft_val) ** 2).real
    return np.round(autocorr).astype(int)


def estimate_marginal_cdf(p, k, num_samples=5_000_000, seed=42):
    """
    Estimate F(T) = Pr[B(1) ≤ T] for a random k-subset D12 of {1,...,p-1} union {0}.

    D12 = {0} ∪ S where S is a (k-1)-subset of {1,...,p-1}, so |D12| = k.
    But wait -- from the code, D12 doesn't necessarily contain 0.

    Actually, from invariant_hunt.py: D12 = {0} ∪ S where S is a (k)-subset of {1,...,p-1}.
    So |D12| = k+1 = (p-1)/2.

    B(d) = autocorrelation of D12 at shift d = #{(a,b) ∈ D12×D12 : a-b ≡ d (mod p)}.
    By Z_p symmetry, B(1) has the same distribution as B(d) for any nonzero d.

    Returns: dict mapping T -> Pr[B(1) <= T] (CDF values)
    """
    rng = np.random.default_rng(seed)
    d12_size = (p - 1) // 2  # |D12|
    k_choose = d12_size - 1  # size of S in {1,...,p-1}

    # Sample random D12 sets and compute B(1) for each
    population = np.arange(1, p)  # {1, ..., p-1}
    b1_values = np.zeros(num_samples, dtype=np.int32)

    for i in range(num_samples):
        # D12 = {0} ∪ random (d12_size-1)-subset of {1,...,p-1}
        S = rng.choice(population, size=k_choose, replace=False)
        D12 = np.empty(d12_size, dtype=np.int32)
        D12[0] = 0
        D12[1:] = S

        # B(1) = #{(a,b) ∈ D12×D12 : a-b ≡ 1 (mod p)}
        # For each a ∈ D12, check if (a-1) mod p ∈ D12
        D12_set = set(D12)
        b1 = sum(1 for a in D12 if (a - 1) % p in D12_set)
        b1_values[i] = b1

    # Build CDF
    max_b = int(b1_values.max())
    cdf = {}
    for T in range(max_b + 2):  # include one beyond max
        cdf[T] = float(np.mean(b1_values <= T))

    return cdf, b1_values


def estimate_marginal_cdf_fast(p, num_samples=5_000_000, seed=42):
    """
    Faster vectorized estimation of F(T) = Pr[B(1) ≤ T].

    Uses batch FFT approach for autocorrelation at position 1.
    """
    rng = np.random.default_rng(seed)
    d12_size = (p - 1) // 2
    k_choose = d12_size - 1

    population = np.arange(1, p)

    # Process in batches for memory efficiency
    batch_size = min(100_000, num_samples)
    n_batches = (num_samples + batch_size - 1) // batch_size

    all_b1 = []

    for batch_idx in range(n_batches):
        current_batch = min(batch_size, num_samples - batch_idx * batch_size)
        if current_batch <= 0:
            break

        # Build indicator matrix for this batch
        indicators = np.zeros((current_batch, p), dtype=np.float64)
        indicators[:, 0] = 1.0  # 0 is always in D12

        for i in range(current_batch):
            S = rng.choice(population, size=k_choose, replace=False)
            indicators[i, S] = 1.0

        # Compute B(1) via FFT: B(d) = IFFT(|FFT(indicator)|^2)[d]
        fft_vals = np.fft.fft(indicators, axis=1)
        autocorr = np.fft.ifft(np.abs(fft_vals) ** 2, axis=1).real
        b1_batch = np.round(autocorr[:, 1]).astype(np.int32)
        all_b1.append(b1_batch)

    b1_values = np.concatenate(all_b1)

    # Build CDF
    max_b = int(b1_values.max())
    cdf = {}
    for T in range(max_b + 2):
        cdf[T] = float(np.mean(b1_values <= T))

    return cdf, b1_values


def compute_thresholds(D11_set, A_values, p):
    """
    Compute threshold T(d) for each nonzero position d.

    For d ∈ D11: A(d) + B(d) ≤ thresh_binding  → T(d) = thresh_binding - A(d)
    For d ∈ D22: A(d) + B(d) ≤ thresh_loose    → T(d) = thresh_loose - A(d)

    thresh_binding = (p-3)/2
    thresh_loose = (p+3)/2
    """
    thresh_binding = (p - 3) // 2
    thresh_loose = (p + 3) // 2

    thresholds = {}
    for d in range(1, p):
        A_d = int(A_values[d])
        if d in D11_set:
            thresholds[d] = thresh_binding - A_d
        else:
            thresholds[d] = thresh_loose - A_d

    return thresholds


def compute_c0_for_orbit(D11_set, p, cdf, N_estimate, total_d12):
    """
    Compute c₀ for a given D11.

    c₀ = Pr[valid | D11] / ∏_d F(T(d))
       = (N_estimate / total_d12) / ∏_d F(T(d))

    Work in log space to avoid underflow.
    """
    A = autocorrelation_single(D11_set, p)
    thresholds = compute_thresholds(D11_set, A, p)

    # Compute log(∏ F(T(d))) = Σ log(F(T(d)))
    log_marginal_product = 0.0
    threshold_values = []

    max_cdf_T = max(cdf.keys())

    for d in range(1, p):
        T_d = thresholds[d]
        threshold_values.append(T_d)

        if T_d >= max_cdf_T:
            # F(T) = 1 for large T, log(1) = 0
            pass
        elif T_d < 0:
            # Impossible constraint: F(T) = 0
            log_marginal_product = float('-inf')
            break
        else:
            f_t = cdf.get(T_d, 0.0)
            if f_t <= 0:
                log_marginal_product = float('-inf')
                break
            log_marginal_product += log(f_t)

    # Compute log(Pr[valid])
    if N_estimate <= 0:
        return {
            "c0": 0.0,
            "log_c0": float('-inf'),
            "log_pr_valid": float('-inf'),
            "log_marginal_product": log_marginal_product,
            "threshold_values": sorted(threshold_values),
            "threshold_spread": max(threshold_values) - min(threshold_values),
            "threshold_var": float(np.var(threshold_values)),
        }

    log_pr_valid = log(N_estimate) - log(total_d12)
    log_c0 = log_pr_valid - log_marginal_product

    try:
        c0 = exp(log_c0)
    except OverflowError:
        c0 = float('inf')

    return {
        "c0": c0,
        "log_c0": log_c0,
        "log_pr_valid": log_pr_valid,
        "log_marginal_product": log_marginal_product,
        "threshold_values": sorted(threshold_values),
        "threshold_spread": max(threshold_values) - min(threshold_values),
        "threshold_var": float(np.var(threshold_values)),
    }


def run_p31(num_cdf_samples=5_000_000):
    """Run c₀ estimation at p=31."""
    p = 31
    n = 16
    d11_size = n  # (p+1)/2 = 16
    d12_size = n - 1  # (p-1)/2 = 15
    k_choose = d12_size - 1  # 14
    total_d12 = comb(p - 1, k_choose)  # C(30, 14) = 145422675

    print(f"\n{'='*80}")
    print(f"c₀ ESTIMATION: p = {p}")
    print(f"{'='*80}")
    print(f"  n = {n}, |D11| = {d11_size}, |D12| = {d12_size}")
    print(f"  Total D12 = {total_d12}")
    print(f"  Binding threshold = {(p-3)//2}, Loose threshold = {(p+3)//2}")

    # Step 1: Estimate F(T) via sampling
    print(f"\n  Step 1: Estimating marginal CDF F(T) with {num_cdf_samples:,} samples...")
    t0 = time.time()
    cdf, b1_values = estimate_marginal_cdf_fast(p, num_samples=num_cdf_samples)
    t1 = time.time()
    print(f"  Done in {t1-t0:.1f}s")
    print(f"  B(1) stats: mean={np.mean(b1_values):.3f}, std={np.std(b1_values):.3f}")
    print(f"  B(1) range: [{b1_values.min()}, {b1_values.max()}]")
    print(f"  CDF: {dict(sorted(cdf.items()))}")

    # Step 2: Load orbit data from existing results
    print(f"\n  Step 2: Loading orbit data...")
    results_path = os.path.join(os.path.dirname(__file__), "near_flat_p31_results.json")
    with open(results_path) as f:
        orbit_data = json.load(f)

    orbits = orbit_data["orbit_details"]
    print(f"  Loaded {len(orbits)} orbits (A-flat + near-flat)")

    # Step 3: Compute c₀ for each orbit
    print(f"\n  Step 3: Computing c₀ for each orbit...")
    t2 = time.time()

    results = []
    for orb in orbits:
        D11 = frozenset(orb["canonical"])
        N_est = orb["N_estimate"]
        num_valid = orb["num_valid"]
        total_sampled = orb["total_sampled"]

        c0_data = compute_c0_for_orbit(D11, p, cdf, N_est, total_d12)

        result = {
            "canonical": orb["canonical"],
            "max_A": orb["max_A"],
            "is_aflat": orb["is_aflat"],
            "num_valid": num_valid,
            "total_sampled": total_sampled,
            "N_estimate": N_est,
            "rate": orb["rate"],
            "c0": c0_data["c0"],
            "log_c0": c0_data["log_c0"],
            "log_pr_valid": c0_data["log_pr_valid"],
            "log_marginal_product": c0_data["log_marginal_product"],
            "threshold_spread": c0_data["threshold_spread"],
            "threshold_var": c0_data["threshold_var"],
            "threshold_values": c0_data["threshold_values"],
        }
        results.append(result)

    t3 = time.time()
    print(f"  Done in {t3-t2:.1f}s")

    # Step 4: Analysis
    print(f"\n  Step 4: Analysis")
    print(f"  {'='*70}")

    working = [r for r in results if r["num_valid"] > 0]
    non_working = [r for r in results if r["num_valid"] == 0]
    aflat_working = [r for r in working if r["is_aflat"]]
    aflat_non_working = [r for r in non_working if r["is_aflat"]]
    nf_working = [r for r in working if not r["is_aflat"]]

    print(f"\n  Working orbits ({len(working)}):")
    for r in sorted(working, key=lambda x: -x["c0"]):
        print(f"    D11={r['canonical'][:6]}... | max_A={r['max_A']} | "
              f"N_est={r['N_estimate']:.0f} | c₀={r['c0']:.4f} | "
              f"spread={r['threshold_spread']} | var={r['threshold_var']:.2f} | "
              f"aflat={r['is_aflat']}")

    if working:
        c0_vals = [r["c0"] for r in working]
        print(f"\n  Working orbits c₀: min={min(c0_vals):.4f}, max={max(c0_vals):.4f}, "
              f"mean={np.mean(c0_vals):.4f}, median={np.median(c0_vals):.4f}")

    if aflat_working:
        c0_vals_af = [r["c0"] for r in aflat_working]
        print(f"  A-flat working c₀: min={min(c0_vals_af):.4f}, max={max(c0_vals_af):.4f}, "
              f"mean={np.mean(c0_vals_af):.4f}")

    # Check non-working orbits - what's their marginal product?
    print(f"\n  Non-working orbits ({len(non_working)}): checking marginal products...")
    finite_lmp = [r["log_marginal_product"] for r in non_working
                  if r["log_marginal_product"] > float('-inf')]
    if finite_lmp:
        print(f"    log(∏F): min={min(finite_lmp):.2f}, max={max(finite_lmp):.2f}, "
              f"mean={np.mean(finite_lmp):.2f}")
    impossible = sum(1 for r in non_working if r["log_marginal_product"] == float('-inf'))
    print(f"    {impossible} orbits have impossible constraints (some T(d) < 0)")

    # Schur-convexity test: does more threshold spread → higher c₀?
    if working:
        print(f"\n  Schur-convexity test (threshold spread vs c₀):")
        for r in sorted(working, key=lambda x: x["threshold_spread"]):
            print(f"    spread={r['threshold_spread']:2d} var={r['threshold_var']:7.2f} "
                  f"→ c₀={r['c0']:.4f}")

    # A-flat vs non-A-flat comparison
    print(f"\n  A-flat working: {len(aflat_working)}, non-A-flat working: {len(nf_working)}")
    print(f"  A-flat non-working: {len(aflat_non_working)}, non-A-flat non-working: {len([r for r in non_working if not r['is_aflat']])}")

    # Summary stats
    summary = {
        "p": p,
        "num_cdf_samples": num_cdf_samples,
        "total_d12": total_d12,
        "cdf": {str(k): v for k, v in sorted(cdf.items())},
        "b1_mean": float(np.mean(b1_values)),
        "b1_std": float(np.std(b1_values)),
        "num_orbits": len(results),
        "num_working": len(working),
        "num_aflat_working": len(aflat_working),
        "num_nf_working": len(nf_working),
        "c0_working_stats": {
            "min": min(c0_vals) if working else None,
            "max": max(c0_vals) if working else None,
            "mean": float(np.mean(c0_vals)) if working else None,
            "median": float(np.median(c0_vals)) if working else None,
        } if working else {},
        "c0_aflat_working_stats": {
            "min": min(c0_vals_af) if aflat_working else None,
            "max": max(c0_vals_af) if aflat_working else None,
            "mean": float(np.mean(c0_vals_af)) if aflat_working else None,
        } if aflat_working else {},
        "orbit_results": results,
    }

    return summary


def run_p43(num_cdf_samples=5_000_000):
    """Run c₀ estimation at p=43 (stretch goal)."""
    p = 43
    n = 22
    d11_size = n  # 22
    d12_size = n - 1  # 21
    k_choose = d12_size - 1  # 20
    total_d12 = comb(p - 1, k_choose)  # C(42, 20) = very large

    print(f"\n{'='*80}")
    print(f"c₀ ESTIMATION: p = {p}")
    print(f"{'='*80}")
    print(f"  n = {n}, |D11| = {d11_size}, |D12| = {d12_size}")
    print(f"  Total D12 = {total_d12}")
    print(f"  Binding threshold = {(p-3)//2}, Loose threshold = {(p+3)//2}")

    # Step 1: Estimate F(T) via sampling
    print(f"\n  Step 1: Estimating marginal CDF F(T) with {num_cdf_samples:,} samples...")
    t0 = time.time()
    cdf, b1_values = estimate_marginal_cdf_fast(p, num_samples=num_cdf_samples, seed=123)
    t1 = time.time()
    print(f"  Done in {t1-t0:.1f}s")
    print(f"  B(1) stats: mean={np.mean(b1_values):.3f}, std={np.std(b1_values):.3f}")
    print(f"  B(1) range: [{b1_values.min()}, {b1_values.max()}]")
    print(f"  CDF: {dict(sorted(cdf.items()))}")

    # Step 2: Load orbit data
    print(f"\n  Step 2: Loading orbit data...")
    results_path = os.path.join(os.path.dirname(__file__), "near_flat_p43_results.json")
    with open(results_path) as f:
        orbit_data = json.load(f)

    orbits = orbit_data["orbit_results"]
    working_orbits = [o for o in orbits if o.get("working", False)]
    print(f"  Loaded {len(orbits)} orbits, {len(working_orbits)} working")

    # For p=43, we don't have exact N. For working orbits found by SA,
    # we can only say N >= 1. We'll estimate Pr[valid] from SA success.
    # Actually, we can compute the marginal product for all orbits and see
    # if the working ones have larger marginal product.
    print(f"\n  Step 3: Computing marginal products for all orbits...")
    t2 = time.time()

    results = []
    for orb in orbits:
        D11 = frozenset(orb["canonical"])

        # We don't have N for p=43, so compute marginal product only
        A = autocorrelation_single(D11, p)
        thresholds = compute_thresholds(D11, A, p)

        max_cdf_T = max(cdf.keys())
        log_marginal_product = 0.0
        threshold_values = []

        for d in range(1, p):
            T_d = thresholds[d]
            threshold_values.append(T_d)

            if T_d >= max_cdf_T:
                pass  # F(T) ≈ 1
            elif T_d < 0:
                log_marginal_product = float('-inf')
                break
            else:
                f_t = cdf.get(T_d, 0.0)
                if f_t <= 0:
                    log_marginal_product = float('-inf')
                    break
                log_marginal_product += log(f_t)

        result = {
            "canonical": orb["canonical"],
            "max_A": orb["max_A"],
            "is_aflat": orb["is_aflat"],
            "working": orb.get("working", False),
            "sa_cost": orb.get("sa_cost", None),
            "log_marginal_product": log_marginal_product,
            "threshold_spread": max(threshold_values) - min(threshold_values) if threshold_values else 0,
            "threshold_var": float(np.var(threshold_values)) if threshold_values else 0,
        }
        results.append(result)

    t3 = time.time()
    print(f"  Done in {t3-t2:.1f}s")

    # Analysis
    print(f"\n  Step 4: Analysis")
    working_results = [r for r in results if r["working"]]
    non_working_results = [r for r in results if not r["working"]]

    print(f"\n  Working orbits ({len(working_results)}):")
    for r in sorted(working_results, key=lambda x: -x["log_marginal_product"]):
        print(f"    D11={r['canonical'][:6]}... | max_A={r['max_A']} | "
              f"log(∏F)={r['log_marginal_product']:.2f} | "
              f"spread={r['threshold_spread']} | var={r['threshold_var']:.2f}")

    print(f"\n  Non-working orbits log(∏F) stats:")
    finite_lmp = [r["log_marginal_product"] for r in non_working_results
                  if r["log_marginal_product"] > float('-inf')]
    if finite_lmp:
        print(f"    min={min(finite_lmp):.2f}, max={max(finite_lmp):.2f}, mean={np.mean(finite_lmp):.2f}")

    if working_results:
        w_lmp = [r["log_marginal_product"] for r in working_results
                 if r["log_marginal_product"] > float('-inf')]
        if w_lmp:
            print(f"\n  Working orbits log(∏F): mean={np.mean(w_lmp):.2f}")
        if finite_lmp:
            print(f"  Non-working orbits log(∏F): mean={np.mean(finite_lmp):.2f}")
            print(f"  Gap: {np.mean(w_lmp) - np.mean(finite_lmp):.2f}" if w_lmp else "")

    summary = {
        "p": p,
        "num_cdf_samples": num_cdf_samples,
        "total_d12": total_d12,
        "cdf": {str(k): v for k, v in sorted(cdf.items())},
        "b1_mean": float(np.mean(b1_values)),
        "b1_std": float(np.std(b1_values)),
        "num_orbits": len(results),
        "num_working": len(working_results),
        "orbit_results": results,
    }

    return summary


if __name__ == "__main__":
    all_results = {}

    # Run p=31 (main target)
    p31_results = run_p31(num_cdf_samples=5_000_000)
    all_results["p31"] = p31_results

    # Run p=43 (stretch goal)
    try:
        p43_results = run_p43(num_cdf_samples=5_000_000)
        all_results["p43"] = p43_results
    except Exception as e:
        print(f"\n  p=43 failed: {e}")
        import traceback
        traceback.print_exc()

    # Save results
    out_path = os.path.join(os.path.dirname(__file__), "c0_estimation_results.json")

    # Convert any non-serializable values
    def sanitize(obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        elif isinstance(obj, (np.floating,)):
            return float(obj)
        elif isinstance(obj, float) and (obj == float('inf') or obj == float('-inf')):
            return str(obj)
        elif isinstance(obj, dict):
            return {k: sanitize(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [sanitize(x) for x in obj]
        return obj

    with open(out_path, 'w') as f:
        json.dump(sanitize(all_results), f, indent=2)

    print(f"\n  Results saved to {out_path}")
