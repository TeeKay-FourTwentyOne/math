"""
Phase 1: Exhaustive enumeration of valid (D11, D12) pairs for small primes p ≡ 3 mod 4.

For p = 7, 11, 19 (n = 4, 6, 10), exhaustively enumerate ALL symmetric D11 of size
(p+1)/2 and ALL D12 containing 0 of size (p-1)/2, checking which pairs satisfy all
Ramsey constraints.

Key question: Do ALL symmetric D11 admit a valid D12, or only special ones?

Convention: Uses Delta from ramsey_core (a-d convention). The binding constraint is
    Delta(D11,D11,d) + Delta(D12,D12,d) <= (p-3)/2  for all d in D11
"""

import sys
import os
import cmath
import math
import time
import json
from itertools import combinations
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ramsey_core import BlockCirculantGraph, verify_construction, Delta


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


def build_d12_indicator_matrix(p, d12_size):
    """Build indicator matrix for all D12 candidates (containing 0, size d12_size).

    Returns (matrix, list_of_D12_sets).
    """
    remaining = d12_size - 1  # 0 is always included
    candidates = list(range(1, p))
    d12_list = []
    for chosen in combinations(candidates, remaining):
        d12_list.append(frozenset({0} | set(chosen)))

    # Build indicator matrix
    matrix = np.zeros((len(d12_list), p), dtype=np.float64)
    for i, D12 in enumerate(d12_list):
        for j in D12:
            matrix[i, j] = 1.0
    return matrix, d12_list


def batch_autocorrelation(indicator_matrix):
    """Compute Delta(S,S,d) for all sets and all d using FFT.

    indicator_matrix: (num_sets, m) where each row is an indicator vector.
    Returns: (num_sets, m) integer matrix where [i,d] = Delta(S_i, S_i, d).
    """
    fft_vals = np.fft.fft(indicator_matrix, axis=1)
    autocorr = np.fft.ifft(np.abs(fft_vals) ** 2, axis=1).real
    return np.round(autocorr).astype(int)


def autocorrelation_single(indicator, m):
    """Compute Delta(S,S,d) for a single set using FFT."""
    fft_val = np.fft.fft(indicator)
    autocorr = np.fft.ifft(np.abs(fft_val) ** 2).real
    return np.round(autocorr).astype(int)


def fourier_coeff(S, k, m):
    """hat{S}(k) = sum_{s in S} omega^{ks}."""
    omega = cmath.exp(2j * cmath.pi / m)
    return sum(omega ** (k * s) for s in S)


def power_spectrum(S, m):
    """Compute |hat{S}(k)|^2 for k = 0,...,m-1."""
    return [abs(fourier_coeff(S, k, m)) ** 2 for k in range(m)]


def analyze_valid_pair(p, D11, D12, A_vals, B_vals):
    """Collect spectral and structural data on a valid (D11, D12) pair."""
    m = p

    ps_d11 = power_spectrum(D11, m)
    ps_d12 = power_spectrum(D12, m)

    # Combined power spectrum
    P = [ps_d11[k] + ps_d12[k] for k in range(m)]

    # Legendre pair comparison: ideal P(k) = p for all k > 0
    P_nonzero = P[1:]
    legendre_deviation = math.sqrt(sum((pk - p) ** 2 for pk in P_nonzero) / len(P_nonzero))

    # D12 asymmetry: |D12 ∩ D12^T|
    D12T = frozenset((-x) % m for x in D12)
    asymmetry = len(D12 & D12T)

    # Cyclotomic class distribution (quadratic residues mod p)
    QR = set()
    for x in range(1, p):
        QR.add((x * x) % p)
    QNR = set(range(1, p)) - QR

    d11_qr = len(D11 & QR)
    d11_qnr = len(D11 & QNR)
    d12_qr = len(D12 & QR)
    d12_qnr = len(D12 & QNR)

    # Constraint margins at each d in D11
    threshold = (p - 3) // 2
    margins = []
    for d in sorted(D11):
        margins.append(threshold - (A_vals[d] + B_vals[d]))

    return {
        "legendre_deviation": legendre_deviation,
        "mean_P_nonzero": sum(P_nonzero) / len(P_nonzero),
        "asymmetry": asymmetry,
        "D11_QR": d11_qr,
        "D11_QNR": d11_qnr,
        "D12_QR": d12_qr,
        "D12_QNR": d12_qnr,
        "margins": margins,
        "min_margin": min(margins),
        "max_margin": max(margins),
        "D12_is_symmetric": D12 == D12T,
    }


def run_prime(p):
    """Run exhaustive enumeration for prime p ≡ 3 mod 4."""
    n = (p + 1) // 2
    d11_size = (p + 1) // 2  # = n
    d12_size = (p - 1) // 2  # = n - 1
    threshold_red = (p - 3) // 2  # = n - 2
    threshold_blue = (p + 3) // 2  # loose

    pairs = get_symmetric_pairs(p)
    num_d11 = math.comb(len(pairs), d11_size // 2)
    num_d12 = math.comb(p - 1, d12_size - 1)
    total = num_d11 * num_d12

    print(f"\n{'=' * 80}")
    print(f"PRIME p = {p}, n = {n}, m = {p}")
    print(f"|D11| = {d11_size}, |D12| = {d12_size}")
    print(f"Binding threshold: A(d)+B(d) <= {threshold_red} for d in D11")
    print(f"Search space: {num_d11} D11 x {num_d12} D12 = {total:,}")
    print(f"{'=' * 80}")

    t0 = time.time()

    # Precompute all D12 indicator vectors and autocorrelations
    print("Precomputing D12 autocorrelations via batch FFT...")
    d12_matrix, d12_list = build_d12_indicator_matrix(p, d12_size)
    B_matrix = batch_autocorrelation(d12_matrix)  # (num_d12, p)
    print(f"  {len(d12_list)} D12 candidates processed in {time.time() - t0:.1f}s")

    # Also precompute B(p-d) for V2V2 constraints
    neg_indices = np.array([(p - d) % p for d in range(p)])
    B_neg_matrix = B_matrix[:, neg_indices]  # B_neg_matrix[:, d] = B(p-d)

    valid_pairs = []
    d11_valid_counts = {}
    d11_list = list(enumerate_symmetric_D11(p, d11_size))

    for d11_idx, D11 in enumerate(d11_list):
        # Compute A(d) for this D11
        d11_indicator = np.zeros(p, dtype=np.float64)
        for j in D11:
            d11_indicator[j] = 1.0
        A = autocorrelation_single(d11_indicator, p)

        # D11 mask for binding constraint
        d11_mask = np.zeros(p, dtype=bool)
        for j in D11:
            d11_mask[j] = True

        # Complement mask (d not in D11, d > 0)
        complement_mask = ~d11_mask & (np.arange(p) > 0)

        # Vectorized constraint check: A[d] + B[d] for all D12 candidates
        # F[i, d] = A(d) + B_i(d)
        F_matrix = A[np.newaxis, :] + B_matrix  # (num_d12, p)

        # Binding: max over d in D11 must be <= threshold_red
        if d11_mask.any():
            d11_max = F_matrix[:, d11_mask].max(axis=1)
            valid_mask = d11_max <= threshold_red
        else:
            valid_mask = np.ones(len(d12_list), dtype=bool)

        # Loose V1V1 blue: max over d not in D11 must be <= threshold_blue
        if complement_mask.any():
            comp_max = F_matrix[:, complement_mask].max(axis=1)
            valid_mask &= comp_max <= threshold_blue

        # Loose V2V2 red: A(d) + B(p-d) for d not in D11
        G_matrix = A[np.newaxis, :] + B_neg_matrix  # A(d) + B(p-d)
        if complement_mask.any():
            comp_max_v2 = G_matrix[:, complement_mask].max(axis=1)
            valid_mask &= comp_max_v2 <= threshold_blue

        count = int(valid_mask.sum())
        d11_valid_counts[D11] = count

        if count > 0:
            valid_indices = np.where(valid_mask)[0]
            for idx in valid_indices:
                D12 = d12_list[idx]
                B_vals = {d: int(B_matrix[idx, d]) for d in range(p)}
                A_vals = {d: int(A[d]) for d in range(p)}
                valid_pairs.append((D11, D12, A_vals, B_vals))

        if (d11_idx + 1) % max(1, len(d11_list) // 10) == 0:
            elapsed = time.time() - t0
            print(f"  D11 {d11_idx + 1}/{len(d11_list)}, "
                  f"valid so far: {len(valid_pairs)}, {elapsed:.1f}s")

    elapsed = time.time() - t0
    print(f"\nCompleted in {elapsed:.1f}s: {total:,} pairs checked, "
          f"{len(valid_pairs)} valid")

    # Summary: which D11 admit valid D12?
    print(f"\nD11 ANALYSIS ({len(d11_list)} symmetric D11 of size {d11_size}):")
    for D11 in d11_list:
        count = d11_valid_counts[D11]
        pct = 100 * count / num_d12 if num_d12 > 0 else 0
        print(f"  D11 = {sorted(D11)}: {count}/{num_d12} valid D12 ({pct:.2f}%)")

    admits_all = all(v > 0 for v in d11_valid_counts.values())
    admits_none = all(v == 0 for v in d11_valid_counts.values())
    print(f"\nDo ALL symmetric D11 admit a valid D12? {admits_all}")
    if not admits_all and not admits_none:
        working = sum(1 for v in d11_valid_counts.values() if v > 0)
        print(f"  {working}/{len(d11_list)} D11 have at least one valid D12")

    # Verify a sample of valid pairs using full verify_construction
    if valid_pairs:
        print(f"\nVERIFICATION SAMPLE (checking {min(5, len(valid_pairs))} pairs "
              f"with full verify_construction):")
        for i, (D11, D12, A_vals, B_vals) in enumerate(valid_pairs[:5]):
            D22 = frozenset(range(1, p)) - D11
            G = BlockCirculantGraph(n=n, D11=set(D11), D12=set(D12), D22=set(D22))
            result = verify_construction(G)
            status = "PASS" if result.valid else "FAIL"
            print(f"  Pair {i}: {status} (max_red={result.max_red_common}, "
                  f"max_blue={result.max_blue_common})")
            if not result.valid:
                print(f"    VIOLATIONS: {result.violations[:5]}")

    # Spectral analysis
    if valid_pairs:
        print(f"\nSPECTRAL ANALYSIS of {len(valid_pairs)} valid pairs:")
        analyses = []
        # Analyze up to 200 pairs (spectral analysis is slow for large sets)
        for D11, D12, A_vals, B_vals in valid_pairs[:200]:
            analysis = analyze_valid_pair(p, D11, D12, A_vals, B_vals)
            analyses.append(analysis)

        leg_devs = [a["legendre_deviation"] for a in analyses]
        mean_Ps = [a["mean_P_nonzero"] for a in analyses]
        asyms = [a["asymmetry"] for a in analyses]
        margins = [a["min_margin"] for a in analyses]

        print(f"  Legendre deviation: min={min(leg_devs):.3f}, "
              f"max={max(leg_devs):.3f}, mean={sum(leg_devs) / len(leg_devs):.3f}")
        print(f"  Mean P(k>0): min={min(mean_Ps):.3f}, "
              f"max={max(mean_Ps):.3f} (ideal: {p})")
        print(f"  D12∩D12^T size: {sorted(set(asyms))}")
        print(f"  Min margin (slack): {sorted(set(margins))}")

        # QR distribution
        qr_patterns = defaultdict(int)
        for a in analyses:
            pattern = (a["D11_QR"], a["D11_QNR"], a["D12_QR"], a["D12_QNR"])
            qr_patterns[pattern] += 1
        print(f"\n  QR distribution (D11_QR, D11_QNR, D12_QR, D12_QNR):")
        for pattern, count in sorted(qr_patterns.items()):
            print(f"    {pattern}: {count} pairs")

        # D12 symmetry
        sym_count = sum(1 for a in analyses if a["D12_is_symmetric"])
        print(f"\n  D12 is symmetric: {sym_count}/{len(analyses)}")

    return {
        "p": p,
        "n": n,
        "total_checked": total,
        "num_valid": len(valid_pairs),
        "all_D11_work": admits_all,
        "d11_counts": {str(sorted(k)): v for k, v in d11_valid_counts.items()},
        "valid_pairs": [
            (sorted(D11), sorted(D12))
            for D11, D12, _, _ in valid_pairs[:100]
        ],
    }


def main():
    primes = [7, 11, 19]
    all_results = {}

    for p in primes:
        result = run_prime(p)
        all_results[str(p)] = result

    # Cross-prime summary
    print(f"\n{'=' * 80}")
    print("CROSS-PRIME SUMMARY")
    print(f"{'=' * 80}")
    for p_str, r in all_results.items():
        p = int(p_str)
        print(f"  p={p}: {r['num_valid']} valid pairs out of {r['total_checked']:,}, "
              f"all D11 work: {r['all_D11_work']}")

    # Save results
    output_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "enumeration_results.json"
    )
    with open(output_path, "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\nResults saved to {output_path}")


if __name__ == "__main__":
    main()
