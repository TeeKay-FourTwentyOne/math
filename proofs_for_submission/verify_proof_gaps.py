#!/usr/bin/env python3
"""
Rigorous verification of the HC proof for small primes.

For each p, enumerate ALL symmetric D11, ALL D12, and verify:
1. Which D11 give valid D12 (N(D11) > 0)?
2. Does the improved margin bound hold?
3. Is the "flat A-profile" assumption valid?
4. How much Q-mass falls on non-achievable profiles?
"""

from math import comb, log2
from itertools import combinations
from collections import Counter
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


def compute_B_profile(D12_set, p):
    k = (p - 1) // 2
    D12 = set(D12_set)
    profile = []
    for d in range(1, k + 1):
        count = sum(1 for x in D12 if (x + d) % p in D12)
        profile.append(count)
    return tuple(profile)


def compute_A_profile(D11_set, p):
    k = (p - 1) // 2
    D11 = set(D11_set)
    profile = []
    for d in range(1, k + 1):
        count = sum(1 for x in D11 if (x + d) % p in D11)
        profile.append(count)
    return profile


def get_thresholds(A_prof, D11_set, p):
    """Compute binding thresholds T(d) for each representative d."""
    k = (p - 1) // 2
    n = (p + 1) // 2
    thresholds = []
    for d_idx in range(k):
        d = d_idx + 1
        A_d = A_prof[d_idx]
        if d in D11_set or (p - d) in D11_set:
            # d ∈ D11: tightest constraint
            # V1V1 red: A(d)+B(d) ≤ n-2 → B(d) ≤ n-2-A(d)
            # V2V2 blue: C(d)+B(p-d) ≤ n-3, C(d)=A(d)-1 → B(d) ≤ n-2-A(d)
            T_d = n - 2 - A_d
        else:
            # d ∈ D22
            # V2V2 red: C(d)+B(p-d) ≤ n-2, C(d)=A(d)-3 → B(d) ≤ n+1-A(d)
            # V1V1 blue: A(d)+B(d) ≤ n+1 (slack)
            # The binding constraint for D22 is V2V2 red
            T_d = n + 1 - A_d
        thresholds.append(T_d)
    return thresholds


def is_valid_D12(B_prof, thresholds):
    return all(B_prof[i] <= thresholds[i] for i in range(len(B_prof)))


def all_symmetric_D11(p):
    """Generate all symmetric D11 subsets of {1,...,p-1} with |D11| = (p+1)/2."""
    k = (p - 1) // 2
    n = (p + 1) // 2
    n_pairs = n // 2  # number of pairs to select
    # Pairs are {d, p-d} for d = 1, ..., k
    pairs = [(d, p - d) for d in range(1, k + 1)]
    for selected in combinations(range(k), n_pairs):
        D11 = set()
        for idx in selected:
            d, pd = pairs[idx]
            D11.add(d)
            D11.add(pd)
        yield frozenset(D11)


def main():
    primes = [11, 19, 23]
    if len(sys.argv) > 1:
        primes = [int(x) for x in sys.argv[1:]]

    for p in primes:
        k = (p - 1) // 2
        n = (p + 1) // 2
        S = k * (k - 1) // 2
        R = k
        n_D12 = comb(p, k)

        pmf = cycle_pmf(p)
        f_mode = max(pmf)
        log2_fmode = log2(f_mode)
        log2_Cpk = log2(comb(p, k))

        print(f"\n{'='*70}")
        print(f"p = {p}, k = {k}, n = {n}, S = {S}, R = {R}")
        print(f"C(p,k) = {comb(p,k)}, f_mode = {f_mode:.6f}")
        print(f"log2(C(p,k)) = {log2_Cpk:.3f}, R*log2(f_mode) = {R*log2_fmode:.3f}")
        print(f"Condition = {log2_Cpk + R*log2_fmode:.3f}")
        print(f"{'='*70}")

        # Precompute all B-profiles
        print(f"Enumerating {n_D12} D12 subsets...")
        all_B_profiles = {}
        for D12 in combinations(range(p), k):
            bp = compute_B_profile(D12, p)
            all_B_profiles[D12] = bp

        # For each symmetric D11
        n_D11_total = comb(k, n // 2)
        print(f"Checking {n_D11_total} symmetric D11 sets...\n")

        best_N = 0
        best_D11 = None
        best_improved = float('-inf')
        best_improved_D11 = None

        for D11 in all_symmetric_D11(p):
            D11_set = set(D11)
            A_prof = compute_A_profile(D11_set, p)
            thresholds = get_thresholds(A_prof, D11_set, p)

            # Count valid D12
            N = 0
            for D12, bp in all_B_profiles.items():
                if is_valid_D12(bp, thresholds):
                    N += 1

            # Compute Q-mass in E∩H
            # Q(b) = Π f(b_i)
            def Q_value(profile):
                val = 1.0
                for b_i in profile:
                    if b_i < 0 or b_i >= len(pmf):
                        return 0.0
                    val *= pmf[b_i]
                return val

            # Sum Q over all achievable profiles in E∩H
            achievable_in_E = set()
            for D12, bp in all_B_profiles.items():
                if is_valid_D12(bp, thresholds):
                    achievable_in_E.add(bp)

            ach_Q = sum(Q_value(bp) for bp in achievable_in_E)

            # Compute improved margin using achievable Q-mass
            if ach_Q > 0:
                corrected_improved = log2(ach_Q) - R * log2_fmode
            else:
                corrected_improved = float('-inf')

            # Track best
            if N > best_N:
                best_N = N
                best_D11 = (D11, A_prof, thresholds, N, ach_Q, corrected_improved)
            if corrected_improved > best_improved:
                best_improved = corrected_improved
                best_improved_D11 = (D11, A_prof, thresholds, N, ach_Q, corrected_improved)

            # Print if N > 0
            if N > 0:
                print(f"  D11={sorted(D11_set)}")
                print(f"    A={A_prof}, T={thresholds}")
                print(f"    N(D11)={N}, #ach_profiles_in_E={len(achievable_in_E)}")
                print(f"    Σ_ach Q = {ach_Q:.6e}")
                print(f"    Corrected improved margin = {corrected_improved:.3f}")
                print()

        print(f"\n--- Summary for p={p} ---")
        if best_D11:
            D11, A, T, N, aQ, ci = best_D11
            print(f"  Best N = {N} (D11={sorted(D11)})")
            print(f"    A={A}, T={T}")
        if best_improved_D11:
            D11, A, T, N, aQ, ci = best_improved_D11
            print(f"  Best corrected improved = {ci:.3f} (D11={sorted(D11)})")
            print(f"    N={N}")

        # Also check: average N over all D11
        total_N = 0
        n_D11_count = 0
        for D11 in all_symmetric_D11(p):
            D11_set = set(D11)
            A_prof = compute_A_profile(D11_set, p)
            thresholds = get_thresholds(A_prof, D11_set, p)
            N = sum(1 for D12, bp in all_B_profiles.items()
                    if is_valid_D12(bp, thresholds))
            total_N += N
            n_D11_count += 1

        avg_N = total_N / n_D11_count
        print(f"  Average N over {n_D11_count} D11 sets: {avg_N:.2f}")
        print(f"  D11 sets with N > 0: {sum(1 for D11 in all_symmetric_D11(p) for _ in [1] if sum(1 for D12, bp in all_B_profiles.items() if is_valid_D12(bp, get_thresholds(compute_A_profile(set(D11), p), set(D11), p))) > 0)}")


if __name__ == '__main__':
    main()
