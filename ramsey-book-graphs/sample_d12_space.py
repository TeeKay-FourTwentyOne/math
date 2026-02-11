"""
Phase 2: Statistical sampling of D12 space for medium primes p ≡ 3 mod 4.

For p = 23, 31, 43, 47, 59 (where we already have SA solutions), samples random D12
to measure how common valid solutions are, whether D11 choice matters, and what
spectral properties distinguish valid D12.

Three experiments per prime:
1. Fix known D11, sample random D12 → hit rate
2. Random D11 survey → does D11 choice matter?
3. Spectral correction analysis → does D12 correct D11's spectrum?

Convention: |D11| = (p+1)/2, |D12| = (p-1)/2, 0 ∈ D12.
"""

import sys
import os
import json
import time
import math
import cmath
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ramsey_core import BlockCirculantGraph, verify_construction, Delta

# Known solutions for primes ≡ 3 mod 4.
# Stored in the "standard" convention: |D11| = (p+1)/2.
# Solutions that were stored with |D11| = (p-3)/2 are relabeled (V1↔V2 swap).
KNOWN_P3MOD4 = {
    # From fourier_constraints.py KNOWN dict (some need relabeling)
    # n=6, p=11: original |D11|=4, relabeled to |D11|=6
    11: {
        "n": 6,
        "D11": {1, 2, 4, 7, 9, 10},  # was D22
        "D12": {0, 4, 5, 7, 10},  # D12^T of original
    },
    # n=10, p=19: original |D11|=8, relabeled to |D11|=10
    19: {
        "n": 10,
        "D11": {1, 2, 3, 6, 8, 11, 13, 16, 17, 18},  # was D22
        "D12": {0, 2, 6, 8, 9, 12, 13, 17, 18},  # D12^T of original
    },
    # n=12, p=23: original |D11|=12 = (p+1)/2, no relabeling needed
    23: {
        "n": 12,
        "D11": {5, 6, 7, 8, 9, 11, 12, 14, 15, 16, 17, 18},
        "D12": {0, 1, 2, 6, 10, 13, 14, 16, 18, 20, 21},
    },
    # n=16, p=31: |D11|=16 = (p+1)/2, no relabeling
    31: {
        "n": 16,
        "D11": {6, 7, 8, 10, 11, 12, 14, 15, 16, 17, 19, 20, 21, 23, 24, 25},
        "D12": {0, 1, 2, 3, 8, 11, 12, 13, 15, 18, 20, 21, 24, 27, 29},
    },
    # n=22, p=43: |D11|=22 = (p+1)/2, no relabeling
    43: {
        "n": 22,
        "D11": {1, 2, 5, 10, 11, 13, 16, 17, 18, 19, 20, 23, 24, 25, 26, 27,
                30, 32, 33, 38, 41, 42},
        "D12": {0, 2, 5, 6, 8, 11, 15, 16, 20, 24, 25, 27, 28, 31, 32, 34,
                35, 36, 37, 39, 41},
    },
}


def normalize_and_load_solution(p, D11_raw, D12_raw):
    """Normalize solution to |D11| = (p+1)/2 convention."""
    target = (p + 1) // 2
    D11 = set(D11_raw)
    D12 = set(D12_raw)
    if len(D11) == target:
        return D11, D12
    # Need to swap V1 ↔ V2
    D22 = set(range(1, p)) - D11
    D12T = {(-x) % p for x in D12}
    assert len(D22) == target, f"Swap failed: |D22|={len(D22)} != {target}"
    return D22, D12T


def load_registry_solutions():
    """Load solutions from solutions_registry.json that need relabeling."""
    registry_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "solutions_registry.json"
    )
    if not os.path.exists(registry_path):
        return
    with open(registry_path) as f:
        registry = json.load(f)

    primes_3mod4 = {47, 59}  # Primes in registry that need relabeling
    for sol in registry["solutions"]:
        p = sol["m"]
        if p not in primes_3mod4:
            continue
        D11_raw = sol["D11"]
        D12_raw = sol["D12"]
        D11, D12 = normalize_and_load_solution(p, D11_raw, D12_raw)
        n = (p + 1) // 2
        KNOWN_P3MOD4[p] = {"n": n, "D11": D11, "D12": D12}


def autocorrelation_fft(indicator, m):
    """Compute Delta(S,S,d) for all d via FFT."""
    fft_val = np.fft.fft(indicator)
    autocorr = np.fft.ifft(np.abs(fft_val) ** 2).real
    return np.round(autocorr).astype(int)


def random_d12(p, size, rng):
    """Generate a random D12 of given size containing 0."""
    candidates = list(range(1, p))
    chosen = rng.choice(candidates, size=size - 1, replace=False)
    return frozenset({0} | set(chosen.tolist()))


def random_symmetric_d11(p, size, rng):
    """Generate a random symmetric D11 of given even size from {1,...,p-1}."""
    # Build list of negation pairs
    pairs = []
    seen = set()
    for d in range(1, p):
        if d not in seen:
            pairs.append((d, (-d) % p))
            seen.add(d)
            seen.add((-d) % p)
    num_pairs = size // 2
    chosen_idx = rng.choice(len(pairs), size=num_pairs, replace=False)
    D11 = set()
    for i in chosen_idx:
        D11.add(pairs[i][0])
        D11.add(pairs[i][1])
    return frozenset(D11)


def check_constraints_fast(A, B, D11_mask, p):
    """Check all constraints using precomputed A(d) and B(d) arrays.

    Returns (valid, max_violation, violation_count).
    """
    threshold_tight = (p - 3) // 2
    threshold_loose = (p + 3) // 2

    F = A + B
    # Binding: A(d) + B(d) <= threshold_tight for d in D11
    d11_violations = 0
    max_viol = 0
    for d in range(1, p):
        if D11_mask[d]:
            excess = F[d] - threshold_tight
            if excess > 0:
                d11_violations += 1
                max_viol = max(max_viol, excess)
        else:
            # Loose V1V1 blue
            excess = F[d] - threshold_loose
            if excess > 0:
                d11_violations += 1
                max_viol = max(max_viol, excess)

    # Also check V2V2 loose: A(d) + B(p-d) for d not in D11
    for d in range(1, p):
        if not D11_mask[d]:
            bd_neg = int(B[(p - d) % p])
            excess_v2 = int(A[d]) + bd_neg - threshold_loose
            if excess_v2 > 0:
                d11_violations += 1
                max_viol = max(max_viol, excess_v2)

    return d11_violations == 0, max_viol, d11_violations


def experiment1_fixed_d11(p, known_D11, known_D12, num_samples=1_000_000):
    """Experiment 1: Fix known D11, sample random D12."""
    print(f"\n  EXPERIMENT 1: Fixed known D11, sample {num_samples:,} random D12")
    n = (p + 1) // 2
    d12_size = (p - 1) // 2

    # Precompute A(d) for fixed D11
    d11_indicator = np.zeros(p, dtype=np.float64)
    for j in known_D11:
        d11_indicator[j] = 1.0
    A = autocorrelation_fft(d11_indicator, p)

    D11_mask = np.zeros(p, dtype=bool)
    for j in known_D11:
        D11_mask[j] = True

    # Verify known solution first
    d12_indicator = np.zeros(p, dtype=np.float64)
    for j in known_D12:
        d12_indicator[j] = 1.0
    B_known = autocorrelation_fft(d12_indicator, p)
    valid, max_viol, _ = check_constraints_fast(A, B_known, D11_mask, p)
    print(f"    Known D12 valid: {valid} (max_violation={max_viol})")

    # Sample random D12
    rng = np.random.default_rng(42)
    valid_count = 0
    violation_costs = []
    violation_counts_hist = defaultdict(int)
    t0 = time.time()

    batch_size = min(10000, num_samples)
    num_batches = (num_samples + batch_size - 1) // batch_size

    for batch in range(num_batches):
        actual_batch = min(batch_size, num_samples - batch * batch_size)
        # Generate batch of random D12 indicators
        d12_batch = np.zeros((actual_batch, p), dtype=np.float64)
        for i in range(actual_batch):
            D12 = random_d12(p, d12_size, rng)
            for j in D12:
                d12_batch[i, j] = 1.0

        # Batch autocorrelation
        fft_vals = np.fft.fft(d12_batch, axis=1)
        B_batch = np.fft.ifft(np.abs(fft_vals) ** 2, axis=1).real
        B_batch = np.round(B_batch).astype(int)

        # Check each
        for i in range(actual_batch):
            B = B_batch[i]
            is_valid, max_v, num_v = check_constraints_fast(A, B, D11_mask, p)
            if is_valid:
                valid_count += 1
            violation_counts_hist[num_v] += 1
            if not is_valid:
                violation_costs.append(max_v)

        if (batch + 1) % max(1, num_batches // 5) == 0:
            elapsed = time.time() - t0
            checked = (batch + 1) * batch_size
            rate = checked / elapsed
            print(f"    {checked:,}/{num_samples:,} checked, "
                  f"{valid_count} valid, {rate:.0f}/sec")

    elapsed = time.time() - t0
    hit_rate = valid_count / num_samples

    print(f"    Results: {valid_count}/{num_samples:,} valid "
          f"({hit_rate:.6f} = 1 in {1 / hit_rate:.0f})" if hit_rate > 0
          else f"    Results: 0/{num_samples:,} valid")
    print(f"    Time: {elapsed:.1f}s")

    # Violation distribution
    print(f"    Violation count distribution (top 5):")
    for vc, cnt in sorted(violation_counts_hist.items(), key=lambda x: -x[1])[:5]:
        print(f"      {vc} violations: {cnt} ({100 * cnt / num_samples:.1f}%)")

    if violation_costs:
        print(f"    Max violation: min={min(violation_costs)}, "
              f"max={max(violation_costs)}, "
              f"mean={sum(violation_costs) / len(violation_costs):.2f}")

    return {"valid_count": valid_count, "total": num_samples, "hit_rate": hit_rate}


def experiment2_random_d11(p, num_d11=1000, num_d12_per=10000):
    """Experiment 2: Random D11 survey."""
    print(f"\n  EXPERIMENT 2: Sample {num_d11} random D11, "
          f"{num_d12_per} D12 each")
    n = (p + 1) // 2
    d11_size = (p + 1) // 2
    d12_size = (p - 1) // 2

    rng = np.random.default_rng(123)
    d11_results = []
    t0 = time.time()

    for d11_idx in range(num_d11):
        D11 = random_symmetric_d11(p, d11_size, rng)

        d11_indicator = np.zeros(p, dtype=np.float64)
        for j in D11:
            d11_indicator[j] = 1.0
        A = autocorrelation_fft(d11_indicator, p)

        D11_mask = np.zeros(p, dtype=bool)
        for j in D11:
            D11_mask[j] = True

        valid_count = 0
        total_max_viol = 0
        for _ in range(num_d12_per):
            D12 = random_d12(p, d12_size, rng)
            d12_indicator = np.zeros(p, dtype=np.float64)
            for j in D12:
                d12_indicator[j] = 1.0
            B = autocorrelation_fft(d12_indicator, p)
            is_valid, max_v, _ = check_constraints_fast(A, B, D11_mask, p)
            if is_valid:
                valid_count += 1
            total_max_viol += max_v

        d11_results.append({
            "valid": valid_count,
            "avg_max_viol": total_max_viol / num_d12_per,
        })

        if (d11_idx + 1) % max(1, num_d11 // 5) == 0:
            elapsed = time.time() - t0
            print(f"    D11 {d11_idx + 1}/{num_d11}, {elapsed:.1f}s")

    valid_counts = [r["valid"] for r in d11_results]
    avg_viols = [r["avg_max_viol"] for r in d11_results]

    d11_with_valid = sum(1 for v in valid_counts if v > 0)
    total_valid = sum(valid_counts)
    overall_rate = total_valid / (num_d11 * num_d12_per)

    print(f"\n    D11 with at least 1 valid D12: {d11_with_valid}/{num_d11}")
    print(f"    Overall hit rate: {overall_rate:.6f}")
    print(f"    Valid per D11: min={min(valid_counts)}, max={max(valid_counts)}, "
          f"mean={sum(valid_counts) / len(valid_counts):.1f}")
    print(f"    Avg max violation: min={min(avg_viols):.2f}, "
          f"max={max(avg_viols):.2f}, mean={sum(avg_viols) / len(avg_viols):.2f}")

    # Check if D11 choice matters
    std_valid = (sum((v - sum(valid_counts) / num_d11) ** 2
                     for v in valid_counts) / num_d11) ** 0.5
    print(f"    Std of valid count across D11: {std_valid:.2f}")
    print(f"    Does D11 choice matter? "
          f"{'YES (high variance)' if std_valid > sum(valid_counts) / num_d11 * 0.5 else 'NO (low variance)'}")

    return {
        "d11_with_valid": d11_with_valid,
        "total_d11": num_d11,
        "overall_rate": overall_rate,
        "valid_counts": valid_counts[:50],  # sample
    }


def experiment3_spectral(p, known_D11, known_D12):
    """Experiment 3: Spectral correction analysis."""
    print(f"\n  EXPERIMENT 3: Spectral correction analysis")
    m = p

    # Compute power spectra
    d11_indicator = np.zeros(p, dtype=np.float64)
    for j in known_D11:
        d11_indicator[j] = 1.0
    d12_indicator = np.zeros(p, dtype=np.float64)
    for j in known_D12:
        d12_indicator[j] = 1.0

    fft_d11 = np.fft.fft(d11_indicator)
    fft_d12 = np.fft.fft(d12_indicator)

    ps_d11 = np.abs(fft_d11) ** 2
    ps_d12 = np.abs(fft_d12) ** 2
    P = ps_d11 + ps_d12

    # Ideal: P(k) = p for k > 0 (Legendre pair condition)
    P_nonzero = P[1:]
    residual = P_nonzero - p
    rms_residual = np.sqrt(np.mean(residual ** 2))

    print(f"    |hat{{D11}}(k)|^2 range (k>0): [{ps_d11[1:].min():.2f}, "
          f"{ps_d11[1:].max():.2f}]")
    print(f"    |hat{{D12}}(k)|^2 range (k>0): [{ps_d12[1:].min():.2f}, "
          f"{ps_d12[1:].max():.2f}]")
    print(f"    P(k) range (k>0): [{P_nonzero.min():.2f}, {P_nonzero.max():.2f}]")
    print(f"    RMS deviation from Legendre (P(k)=p): {rms_residual:.3f}")
    print(f"    P(0) = {P[0]:.1f} = |D11|^2 + |D12|^2 = "
          f"{len(known_D11)}^2 + {len(known_D12)}^2 = "
          f"{len(known_D11) ** 2 + len(known_D12) ** 2}")

    # Correlation between D11 and D12 spectra
    corr = np.corrcoef(ps_d11[1:], ps_d12[1:])[0, 1]
    print(f"    Correlation(|D11|^2, |D12|^2): {corr:.4f}")
    print(f"    {'NEGATIVE correlation: D12 compensates D11 spectrum' if corr < -0.3 else ''}")
    print(f"    {'POSITIVE correlation: spectra reinforce' if corr > 0.3 else ''}")

    # Phase analysis: for D11 symmetric, hat{D11}(k) is real.
    # What about hat{D12}(k)?
    d12_imag_energy = np.sum(np.imag(fft_d12[1:]) ** 2)
    d12_real_energy = np.sum(np.real(fft_d12[1:]) ** 2)
    print(f"    D12 Fourier: real energy={d12_real_energy:.1f}, "
          f"imag energy={d12_imag_energy:.1f}")
    print(f"    D12 Fourier imaginary fraction: "
          f"{d12_imag_energy / (d12_real_energy + d12_imag_energy):.4f}")

    # "Ideal D12 spectrum": what |hat{D12}(k)|^2 would need to be for P(k) = p?
    ideal_d12_ps = p - ps_d11[1:]
    negative_count = np.sum(ideal_d12_ps < 0)
    print(f"\n    Ideal |D12|^2(k) = p - |D11|^2(k):")
    print(f"      Negative values (infeasible): {negative_count}/{p - 1}")
    if negative_count == 0:
        print(f"      Range: [{ideal_d12_ps.min():.2f}, {ideal_d12_ps.max():.2f}]")
        # How close is actual D12 to ideal?
        actual_d12_ps = ps_d12[1:]
        deviation = np.sqrt(np.mean((actual_d12_ps - ideal_d12_ps) ** 2))
        print(f"      RMS deviation of actual from ideal: {deviation:.3f}")

    return {
        "rms_legendre_residual": float(rms_residual),
        "spectral_correlation": float(corr),
        "d12_imag_fraction": float(
            d12_imag_energy / (d12_real_energy + d12_imag_energy)
        ),
        "infeasible_ideal_count": int(negative_count),
    }


def main():
    load_registry_solutions()

    # Determine sample sizes based on prime
    sample_sizes = {
        23: (500_000, 500, 5000),
        31: (200_000, 200, 2000),
        43: (100_000, 100, 1000),
        47: (100_000, 100, 1000),
        59: (50_000, 50, 500),
    }

    all_results = {}
    for p in sorted(KNOWN_P3MOD4.keys()):
        if p < 23:
            continue  # Skip small primes (covered by enumeration)

        info = KNOWN_P3MOD4[p]
        n = info["n"]
        D11 = info["D11"]
        D12 = info["D12"]

        # Verify sizes
        assert len(D11) == (p + 1) // 2, f"p={p}: |D11|={len(D11)} != {(p + 1) // 2}"
        assert len(D12) == (p - 1) // 2, f"p={p}: |D12|={len(D12)} != {(p - 1) // 2}"
        assert 0 in D12, f"p={p}: 0 not in D12"

        print(f"\n{'=' * 80}")
        print(f"PRIME p = {p}, n = {n}")
        print(f"|D11| = {len(D11)}, |D12| = {len(D12)}")
        print(f"{'=' * 80}")

        num_d12, num_d11, num_d12_per = sample_sizes.get(p, (50000, 50, 500))

        r1 = experiment1_fixed_d11(p, D11, D12, num_samples=num_d12)
        r2 = experiment2_random_d11(p, num_d11=num_d11, num_d12_per=num_d12_per)
        r3 = experiment3_spectral(p, D11, D12)

        all_results[str(p)] = {"exp1": r1, "exp2": r2, "exp3": r3}

    # Summary
    print(f"\n{'=' * 80}")
    print("SUMMARY ACROSS ALL PRIMES")
    print(f"{'=' * 80}")
    for p_str, r in sorted(all_results.items(), key=lambda x: int(x[0])):
        p = int(p_str)
        exp1 = r["exp1"]
        exp2 = r["exp2"]
        exp3 = r["exp3"]
        hr = exp1["hit_rate"]
        print(f"  p={p}: hit_rate={hr:.6f}"
              + (f" (1 in {1 / hr:.0f})" if hr > 0 else " (none found)")
              + f", D11_work={exp2['d11_with_valid']}/{exp2['total_d11']}"
              + f", spec_corr={exp3['spectral_correlation']:.3f}")

    # Save
    output_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "sampling_results.json"
    )
    with open(output_path, "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\nResults saved to {output_path}")


if __name__ == "__main__":
    main()
