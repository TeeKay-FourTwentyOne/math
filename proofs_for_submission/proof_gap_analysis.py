#!/usr/bin/env python3
"""
Analyze two potential gaps in the HC proof:

Gap 1: The P ≥ Q bound uses Pr_Q[E∩H] which includes non-achievable profiles.
        Correct bound uses Σ_{achievable b ∈ E∩H} Q(b) ≤ Pr_Q[E∩H].

Gap 2: The "flat A-profile" assumed in the margin computation may not be
        achievable by any actual D11.

For small primes (p ≤ 23), enumerate exhaustively.
For larger primes, estimate the non-achievable mass correction.
"""

from math import comb, log2, sqrt, pi
from itertools import combinations
import numpy as np
import sys


def cycle_pmf(p):
    k = (p - 1) // 2
    total = comb(p, k)
    pmf = []
    for j in range(k):
        num = p * comb(k - 1, j) * comb(p - k - 1, k - 1 - j)
        denom = (k - j) * total
        pmf.append(num / denom)
    return pmf


def pmf_variance(pmf):
    """Compute variance of a PMF."""
    mean = sum(j * p for j, p in enumerate(pmf))
    var = sum(j**2 * p for j, p in enumerate(pmf)) - mean**2
    return mean, var


def estimate_non_achievable_fraction(p, pmf):
    """Estimate fraction of Q-mass on non-achievable profiles.

    Non-achievable = DFT of autocorrelation has negative values.
    Under Q (independent B(d)), hat(B)(j) = k + 2 Σ B(d) cos(2πjd/p).

    E[hat(B)(j)] = (p+1)/4
    Var[hat(B)(j)] = 4 Var(B) × Σ cos²(2πjd/p)
    """
    k = (p - 1) // 2
    mean_B, var_B = pmf_variance(pmf)

    # E[hat(B)(j)] for j ≠ 0
    # = k + 2 * mean_B * Σ cos(2πjd/p) for d=1..k
    # Σ_{d=1}^k cos(2πjd/p) = -1/2 for all j ≠ 0 (mod p)
    E_hat = k + 2 * mean_B * (-0.5)  # = k - mean_B

    # Var[hat(B)(j)] = 4 var_B * Σ cos²(2πjd/p)
    # Σ cos² ≈ k/2 for generic j
    cos2_sum = sum(np.cos(2 * pi * 1 * d / p)**2 for d in range(1, k+1))
    Var_hat = 4 * var_B * cos2_sum

    if Var_hat <= 0:
        return 0.0, E_hat, 0.0

    sigma_hat = sqrt(Var_hat)
    z = E_hat / sigma_hat

    # Pr[hat(B)(j) < 0] ≈ Φ(-z) for each j
    from scipy.stats import norm
    prob_neg_one = norm.cdf(-z)

    # Union bound: Pr[∃j: hat(B)(j) < 0] ≤ k × prob_neg_one
    # (there are k independent-ish spectral values j=1,...,k)
    delta_ub = min(1.0, k * prob_neg_one)

    return delta_ub, z, prob_neg_one


def compute_exact_gaps(p):
    """For small p, exact computation of gaps."""
    k = (p - 1) // 2
    n = (p + 1) // 2
    S = k * (k - 1) // 2
    R = k
    n_D12 = comb(p, k)

    pmf = cycle_pmf(p)
    f_mode = max(pmf)
    log2_fmode = log2(f_mode)
    log2_Cpk = log2(n_D12)

    print(f"\n{'='*60}")
    print(f"p = {p}, k = {k}, n = {n}, R = {R}, C(p,k) = {n_D12}")
    print(f"f_mode = {f_mode:.6f}, log2(f_mode) = {log2_fmode:.4f}")
    print(f"Trivial Bound condition = {log2_Cpk + R * log2_fmode:.3f}")
    print(f"{'='*60}")

    # Enumerate all D12 and their B-profiles
    D12_profiles = {}
    for D12 in combinations(range(p), k):
        D12_set = set(D12)
        bp = []
        for d in range(1, k+1):
            bp.append(sum(1 for x in D12_set if (x+d) % p in D12_set))
        D12_profiles[D12] = tuple(bp)

    # Profile -> count
    from collections import Counter
    profile_counts = Counter(D12_profiles.values())
    achievable = set(profile_counts.keys())

    # All symmetric D11
    pairs = [(d, p-d) for d in range(1, k+1)]
    n_pairs = n // 2

    best_N = 0
    best_info = None
    all_results = []

    for pair_sel in combinations(range(k), n_pairs):
        D11 = set()
        for idx in pair_sel:
            D11.add(pairs[idx][0])
            D11.add(pairs[idx][1])

        # A-profile
        A = []
        for d in range(1, k+1):
            A.append(sum(1 for x in D11 if (x+d) % p in D11))

        # Thresholds
        T = []
        for i in range(k):
            d = i + 1
            if d in D11 or (p-d) in D11:
                T.append(n - 2 - A[i])
            else:
                T.append(n + 1 - A[i])

        # Count valid D12
        valid_profiles = set()
        N = 0
        for D12, bp in D12_profiles.items():
            if all(bp[i] <= T[i] for i in range(R)):
                N += 1
                valid_profiles.add(bp)

        # Q values
        def Q(b):
            v = 1.0
            for bi in b:
                if bi < 0 or bi >= len(pmf):
                    return 0.0
                v *= pmf[bi]
            return v

        # Achievable Q-mass in E∩H
        ach_Q = sum(Q(bp) for bp in valid_profiles)

        # Full Pr_Q[E∩H] - need to enumerate all b in E∩H
        # For efficiency, only do this for best D11

        all_results.append({
            'D11': sorted(D11), 'A': A, 'T': T, 'N': N,
            'ach_Q': ach_Q, 'n_ach_profiles': len(valid_profiles)
        })

        if N > best_N:
            best_N = N
            best_info = all_results[-1]

    # Print summary
    n_working = sum(1 for r in all_results if r['N'] > 0)
    print(f"\nSymmetric D11 sets: {len(all_results)}")
    print(f"With N > 0: {n_working}")
    print(f"Average N: {sum(r['N'] for r in all_results) / len(all_results):.2f}")

    if best_info:
        r = best_info
        print(f"\nBest D11: {r['D11']}")
        print(f"  A-profile: {r['A']}")
        print(f"  Thresholds: {r['T']}")
        print(f"  N(D11) = {r['N']}")
        print(f"  log2(N) = {log2(r['N']):.3f}")
        print(f"  # achievable profiles in E = {r['n_ach_profiles']}")
        print(f"  Σ_ach Q in E∩H = {r['ach_Q']:.6e}")
        if r['ach_Q'] > 0:
            corrected_std = log2_Cpk + log2(r['ach_Q'])
            corrected_imp = log2(r['ach_Q']) - R * log2_fmode
            print(f"  Corrected standard margin = {corrected_std:.3f}")
            print(f"  Corrected improved margin = {corrected_imp:.3f}")

        # Now compute Pr_Q[E∩H] for best D11's thresholds
        T = r['T']
        total_Q_EH = 0.0
        ach_Q_EH = 0.0

        # Recursive enumeration
        def enum_Q(idx, rem, profile):
            nonlocal total_Q_EH, ach_Q_EH
            if idx == R:
                if rem == 0:
                    bp = tuple(profile)
                    qv = Q(bp)
                    if qv > 0:
                        total_Q_EH += qv
                        if bp in achievable:
                            ach_Q_EH += qv
                return
            max_v = min(T[idx], k-1, rem)
            rest_max = sum(min(T[j], k-1) for j in range(idx+1, R))
            min_v = max(0, rem - rest_max)
            for v in range(min_v, max_v + 1):
                profile.append(v)
                enum_Q(idx + 1, rem - v, profile)
                profile.pop()

        enum_Q(0, S, [])

        print(f"\n  Pr_Q[E∩H] = {total_Q_EH:.6e}")
        print(f"  Σ_ach Q in E∩H = {ach_Q_EH:.6e}")
        non_ach = total_Q_EH - ach_Q_EH
        print(f"  Non-achievable Q = {non_ach:.6e}")
        if total_Q_EH > 0:
            frac = ach_Q_EH / total_Q_EH
            print(f"  Achievable fraction = {frac:.4f}")
            print(f"  Gap correction = {log2(frac):.3f} bits" if frac > 0 else "  GAP: zero achievable mass!")

            std_with_full_Q = log2_Cpk + log2(total_Q_EH)
            std_with_ach_Q = log2_Cpk + log2(ach_Q_EH) if ach_Q_EH > 0 else float('-inf')
            imp_with_full_Q = log2(total_Q_EH) - R * log2_fmode
            imp_with_ach_Q = log2(ach_Q_EH) - R * log2_fmode if ach_Q_EH > 0 else float('-inf')

            print(f"\n  Standard margin (full Q):     {std_with_full_Q:.3f}")
            print(f"  Standard margin (ach only):   {std_with_ach_Q:.3f}")
            print(f"  Improved margin (full Q):     {imp_with_full_Q:.3f}")
            print(f"  Improved margin (ach only):   {imp_with_ach_Q:.3f}")

    # Also check the "flat profile" computation
    print(f"\n--- Flat Profile Analysis ---")
    A_flat_val = n * (n-1) // (2 * R)  # average A value
    remainder = n * (n-1) // 2 - R * A_flat_val
    print(f"  Average A = {n*(n-1)/2/R:.2f}")
    print(f"  Closest integer A = {A_flat_val}, remainder = {remainder}")

    # Check if any D11 has all A ≤ A_flat_val + 1 and ≥ A_flat_val - 1
    near_flat = []
    for r in all_results:
        max_dev = max(abs(a - n*(n-1)/2/R) for a in r['A'])
        if max_dev <= 1.5:
            near_flat.append((max_dev, r))

    if near_flat:
        near_flat.sort(key=lambda x: x[0])
        print(f"  D11 with max |A - mean| ≤ 1.5: {len(near_flat)}")
        md, r = near_flat[0]
        print(f"  Flattest: A={r['A']}, max_dev={md:.2f}, N={r['N']}")
    else:
        print(f"  WARNING: No D11 with near-flat A-profile!")
        min_dev = min(max(abs(a - n*(n-1)/2/R) for a in r['A']) for r in all_results)
        print(f"  Minimum max deviation = {min_dev:.2f}")


def estimate_gaps_large(p):
    """For larger primes, estimate the gaps without enumeration."""
    k = (p - 1) // 2
    n = (p + 1) // 2
    R = k
    pmf = cycle_pmf(p)
    f_mode = max(pmf)
    mean_B, var_B = pmf_variance(pmf)

    log2_fmode = log2(f_mode)
    log2_Cpk = sum(log2((p - i) / (i + 1)) for i in range(k))

    # Estimate non-achievable fraction
    try:
        delta_ub, z, prob_neg_one = estimate_non_achievable_fraction(p, pmf)
    except ImportError:
        from math import erfc
        E_hat = k - mean_B
        cos2_sum = sum(np.cos(2 * pi * 1 * d / p)**2 for d in range(1, k+1))
        Var_hat = 4 * var_B * cos2_sum
        sigma_hat = sqrt(Var_hat) if Var_hat > 0 else 1
        z = E_hat / sigma_hat
        prob_neg_one = 0.5 * erfc(z / sqrt(2))
        delta_ub = min(1.0, k * prob_neg_one)

    # Improved margin from the computation
    # (Use saddle-point from improved_margin_scaled.py)
    from improved_margin_scaled import compute_improved_scaled
    result = compute_improved_scaled(p)

    if result:
        imp = result['improved']
        correction = -log2(1 - delta_ub) if delta_ub < 1 else float('inf')
        corrected = imp - correction if correction != float('inf') else float('-inf')

        print(f"p={p:5d} k={k:4d} | improved={imp:8.2f} | z={z:.3f} "
              f"δ_ub={delta_ub:.4f} correction={correction:.2f} | "
              f"corrected={corrected:.2f} | "
              f"{'OK' if corrected > 0 else 'FAIL'}")
        return corrected
    else:
        print(f"p={p:5d} FAILED to compute")
        return None


def main():
    # Part 1: Exact analysis for small primes
    print("PART 1: EXACT ANALYSIS")
    print("=" * 60)

    for p in [11, 19]:
        compute_exact_gaps(p)

    # Part 2: Estimated gaps for larger primes
    print("\n\nPART 2: ESTIMATED NON-ACHIEVABLE MASS CORRECTION")
    print("=" * 70)
    print(f"{'p':>5} {'k':>4} | {'improved':>8} | {'z':>5} "
          f"{'δ_ub':>6} {'correction':>10} | {'corrected':>9} | status")
    print("-" * 70)

    sys.path.insert(0, 'ramsey-book-graphs')

    primes_3mod4 = [p for p in range(11, 1000)
                    if all(p % d != 0 for d in range(2, int(p**0.5)+1))
                    and p % 4 == 3]

    all_positive = True
    for p in primes_3mod4:
        c = estimate_gaps_large(p)
        if c is not None and c <= 0:
            all_positive = False

    if all_positive:
        print("\nAll corrected margins positive!")
    else:
        print("\nSOME corrected margins are NEGATIVE — gap is significant!")


if __name__ == '__main__':
    main()
