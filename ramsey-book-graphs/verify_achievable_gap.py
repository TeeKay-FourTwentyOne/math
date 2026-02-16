#!/usr/bin/env python3
"""
Verify the gap between Pr_Q[E∩H] and Σ_{achievable b ∈ E∩H} Q(b).

For each prime p ≡ 3 mod 4, enumerate ALL D12 subsets, find all achievable
B-profiles, and compare the Q-mass on achievable vs all profiles in E∩H.
"""

from math import comb, log2, prod
from itertools import combinations
import numpy as np
from collections import Counter


def cycle_pmf(p):
    """Compute cycle PMF f(j) for j = 0, ..., k-1."""
    k = (p - 1) // 2
    total = comb(p, k)
    pmf = []
    for j in range(k):
        num = p * comb(k - 1, j) * comb(p - k - 1, k - 1 - j)
        denom = (k - j) * total
        pmf.append(num / denom)
    return pmf


def compute_B_profile(D12_set, p):
    """Compute B(d) for d = 1, ..., k where k = (p-1)/2."""
    k = (p - 1) // 2
    D12 = set(D12_set)
    profile = []
    for d in range(1, k + 1):
        count = sum(1 for x in D12 if (x + d) % p in D12)
        profile.append(count)
    return tuple(profile)


def flat_D11(p):
    """Build the 'flat' symmetric D11 with delta=1 asymmetry."""
    k = (p - 1) // 2
    n = (p + 1) // 2
    # D11 = symmetric subset of {1,...,p-1} with |D11| = n
    # Take pairs {d, p-d} for d = 1, ..., k
    # Need n/2 = (p+1)/4 pairs (for p ≡ 3 mod 4, n is even)
    n_pairs = n // 2
    # All pairs: {1,p-1}, {2,p-2}, ..., {k, k+1}
    # Take first n_pairs pairs
    D11 = set()
    for d in range(1, n_pairs + 1):
        D11.add(d)
        D11.add(p - d)
    return D11


def compute_A_profile(D11_set, p):
    """Compute A(d) for d = 1, ..., k."""
    k = (p - 1) // 2
    D11 = set(D11_set)
    profile = []
    for d in range(1, k + 1):
        count = sum(1 for x in D11 if (x + d) % p in D11)
        profile.append(count)
    return profile


def main():
    for p in [11, 19, 23]:
        k = (p - 1) // 2
        n = (p + 1) // 2
        S = k * (k - 1) // 2
        R = k  # number of representative distances

        pmf = cycle_pmf(p)
        f_mode = max(pmf)

        # Build D11 and compute A-profile, thresholds
        D11 = flat_D11(p)
        D22 = set(range(1, p)) - D11
        A_prof = compute_A_profile(D11, p)

        # Thresholds: for d ∈ D11: T(d) = (p-3)/2 - A(d)
        #             for d ∈ D22: T(d) = (p+3)/2 - C(d) where C(d) = A(d) - 1 (D11) or A(d) - 3 (D22)
        # Actually from the constraint model:
        # d ∈ D11: tightest is V2V2 blue: C(d) + B(p-d) ≤ n-3, C(d) = A(d)-1
        #   so B(d) ≤ n-3 - C(d) = n-3 - (A(d)-1) = n-2-A(d) = (p-3)/2 - A(d)
        # d ∈ D22: tightest is V1V1 red: A(d) + B(d) ≤ n-2
        #   But A(d) for d ∈ D22... hmm, need to be careful about which constraint is tightest

        # For flat D11: all A(d_i) for d ∈ D11 reps are approximately (p+1)/4
        # Let's compute thresholds correctly
        thresholds = []
        for d_idx in range(k):
            d = d_idx + 1
            A_d = A_prof[d_idx]
            if d in D11 or (p - d) in D11:
                # d ∈ D11: V2V2 blue constraint: B(d) ≤ n-2-A(d)
                T_d = n - 2 - A_d
            else:
                # d ∈ D22: V2V2 red constraint: C(d)+B(p-d)≤n-2
                # C(d) = A(d)-3 for d ∈ D22
                C_d = A_d - 3
                T_d = n - 2 - C_d
            thresholds.append(T_d)

        print(f"\n=== p = {p}, k = {k}, n = {n}, S = {S}, R = {R} ===")
        print(f"  D11 = {sorted(D11)}")
        print(f"  A-profile = {A_prof}")
        print(f"  Thresholds = {thresholds}")
        print(f"  f_mode = {f_mode:.6f}")

        # Enumerate all D12 subsets
        all_elements = list(range(p))
        n_D12 = comb(p, k)
        print(f"  Enumerating {n_D12} D12 subsets...")

        # Collect all achievable B-profiles
        profile_counts = Counter()  # B-profile -> count of D12 with this profile
        for D12 in combinations(all_elements, k):
            bp = compute_B_profile(D12, p)
            profile_counts[bp] += 1

        achievable_profiles = set(profile_counts.keys())
        print(f"  # achievable profiles on H: {len(achievable_profiles)}")

        # Enumerate ALL integer profiles on H ∩ E (with entries in PMF support and sum = S)
        # Profile b ∈ {0,...,k-1}^R with Σ b_i = S and b_i ≤ T_i for all i
        def Q_value(profile):
            """Compute Q(b) = Π f(b_i)."""
            val = 1.0
            for b_i in profile:
                if b_i < 0 or b_i >= len(pmf):
                    return 0.0
                val *= pmf[b_i]
            return val

        # Generate all profiles in E ∩ H by recursion
        max_vals = [min(t, k - 1) for t in thresholds]

        total_Q_EH = 0.0  # Pr_Q[E∩H]
        achievable_Q_EH = 0.0  # Σ_{achievable b ∈ E∩H} Q(b)
        non_achievable_Q_EH = 0.0
        n_profiles_EH = 0
        n_achievable_EH = 0
        n_non_achievable_EH = 0

        def enumerate_profiles(idx, remaining_sum, current_profile):
            nonlocal total_Q_EH, achievable_Q_EH, non_achievable_Q_EH
            nonlocal n_profiles_EH, n_achievable_EH, n_non_achievable_EH

            if idx == R:
                if remaining_sum == 0:
                    bp = tuple(current_profile)
                    q_val = Q_value(bp)
                    if q_val > 0:
                        n_profiles_EH += 1
                        total_Q_EH += q_val
                        if bp in achievable_profiles:
                            achievable_Q_EH += q_val
                            n_achievable_EH += 1
                        else:
                            non_achievable_Q_EH += q_val
                            n_non_achievable_EH += 1
                return

            max_val = min(max_vals[idx], remaining_sum)
            # Min: need enough remaining slots to fill remaining_sum
            remaining_slots = R - idx - 1
            if remaining_slots > 0:
                max_possible_from_rest = sum(max_vals[idx+1:])
                min_val = max(0, remaining_sum - max_possible_from_rest)
            else:
                min_val = remaining_sum

            for v in range(min_val, max_val + 1):
                current_profile.append(v)
                enumerate_profiles(idx + 1, remaining_sum - v, current_profile)
                current_profile.pop()

        enumerate_profiles(0, S, [])

        print(f"  # profiles in E∩H (Q>0): {n_profiles_EH}")
        print(f"    achievable: {n_achievable_EH}")
        print(f"    non-achievable: {n_non_achievable_EH}")
        print(f"  Pr_Q[E∩H] = {total_Q_EH:.6e}")
        print(f"  Σ_{{achievable}} Q = {achievable_Q_EH:.6e}")
        print(f"  Σ_{{non-achiev}} Q = {non_achievable_Q_EH:.6e}")
        if total_Q_EH > 0:
            ratio = achievable_Q_EH / total_Q_EH
            print(f"  Achievable fraction of Q-mass: {ratio:.6f}")
            print(f"  Gap in log2: {log2(ratio) if ratio > 0 else '-inf':.3f} bits")

        # Compute both margins
        log2_Cpk = log2(comb(p, k))
        log2_PrQ = log2(total_Q_EH) if total_Q_EH > 0 else float('-inf')
        log2_ach = log2(achievable_Q_EH) if achievable_Q_EH > 0 else float('-inf')
        R_val = R
        log2_fmode = log2(f_mode)

        standard_margin = log2_Cpk + log2_PrQ
        improved_margin = log2_PrQ - R_val * log2_fmode
        corrected_standard = log2_Cpk + log2_ach
        corrected_improved = log2_ach - R_val * log2_fmode

        print(f"\n  Standard margin (using Pr_Q[E∩H]):  {standard_margin:.3f}")
        print(f"  Corrected standard (achievable only): {corrected_standard:.3f}")
        print(f"  Improved margin (using Pr_Q[E∩H]):   {improved_margin:.3f}")
        print(f"  Corrected improved (achievable only): {corrected_improved:.3f}")

        # Verify P ≥ Q for achievable profiles in E∩H
        min_PQ = float('inf')
        for bp, count in profile_counts.items():
            # Check if in E
            in_E = all(bp[i] <= thresholds[i] for i in range(R))
            if not in_E:
                continue
            P_val = count / n_D12
            Q_val = Q_value(bp)
            if Q_val > 0:
                ratio = P_val / Q_val
                if ratio < min_PQ:
                    min_PQ = ratio
        print(f"  min P/Q on achievable ∩ E∩H: {min_PQ:.4f}")


if __name__ == '__main__':
    main()
