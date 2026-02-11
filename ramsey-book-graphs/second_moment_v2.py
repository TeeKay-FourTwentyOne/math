#!/usr/bin/env python3
"""
Second moment computation for Paley-Zygmund.

For fixed A-flat D11, let X_S = 1[{0} union S is valid], N = sum_S X_S.

E[N] = C(p-1,k) * p_0

E[N^2] = sum_{S,S'} Pr[both valid]
       = C(p-1,k) * sum_s C(k,s)*C(p-1-k,k-s) * p(s)

where p(s) = Pr[S valid AND S' valid | |S cap S'| = s].

The ratio E[N^2]/E[N]^2 = E_s[p(s)] / p_0^2
where E_s is over hypergeometric distribution.

Paley-Zygmund: Pr[N>0] >= E[N]^2/E[N^2].
"""

import numpy as np
from math import comb, log2, sqrt
from itertools import combinations
from collections import defaultdict


def compute_B(S_set, p):
    """Compute B(d) for d=1,...,p-1."""
    B = [0] * p
    for a in S_set:
        for b in S_set:
            if a != b:
                d = (a - b) % p
                B[d] += 1
    return B


def compute_A(D11_set, p):
    """Compute A(d) for d=1,...,p-1."""
    A = [0] * p
    for a in D11_set:
        for b in D11_set:
            if a != b:
                d = (a - b) % p
                A[d] += 1
    return A


def find_best_d11(p):
    """Find an A-flat D11 by trying QR-based constructions."""
    n = (p + 1) // 2

    # Generate all symmetric subsets of size n
    pairs = []
    seen = set()
    for d in range(1, p):
        if d not in seen:
            pairs.append((d, p - d))
            seen.add(d)
            seen.add(p - d)

    num_pairs = len(pairs)
    reps_needed = n // 2

    best_D11 = None
    best_max_A = p
    all_good = []

    for combo in combinations(range(num_pairs), reps_needed):
        D11 = []
        for idx in combo:
            d1, d2 = pairs[idx]
            D11.append(d1)
            if d1 != d2:
                D11.append(d2)
        if len(D11) != n:
            continue

        D11_set = set(D11)
        A = compute_A(D11_set, p)
        max_A = max(A[d] for d in range(1, p))
        if max_A < best_max_A:
            best_max_A = max_A
            best_D11 = sorted(D11)
        if max_A <= (p + 1) // 4:  # perfectly flat
            all_good.append(sorted(D11))

    return best_D11, best_max_A, all_good


def second_moment_exhaustive(p, D11):
    """Compute E[N], E[N^2], and E[N^2]/E[N]^2 by exhaustive enumeration."""
    n = (p + 1) // 2
    k = (p - 3) // 2

    D11_set = set(D11)
    A = compute_A(D11_set, p)

    # Thresholds
    T = [0] * p
    for d in range(1, p):
        if d in D11_set:
            T[d] = (p - 3) // 2 - A[d]
        else:
            T[d] = (p + 3) // 2 - A[d]

    # Enumerate all valid S
    elements = list(range(1, p))
    valid_list = []

    for S in combinations(elements, k):
        S_set = set(S)
        B = compute_B(S_set, p)
        ok = True
        for d in range(1, p):
            if B[d] > T[d]:
                ok = False
                break
        if ok:
            valid_list.append(frozenset(S))

    N_val = len(valid_list)
    total = comb(p - 1, k)
    p_0 = N_val / total

    print(f"  D11 = {D11}")
    print(f"  max A(d) = {max(A[d] for d in range(1,p))}, E[A] = {(p+1)/4:.1f}")
    print(f"  N = #{'{'}valid D12{'}'} = {N_val}")
    print(f"  C(p-1,k) = {total}")
    print(f"  p_0 = Pr[valid] = {p_0:.6f} = 2^{{{log2(p_0):.2f}}}")
    print(f"  E[N] = {N_val}, log2(E[N]) = {log2(N_val) if N_val > 0 else '-inf':.2f}")

    if N_val == 0:
        print("  No valid D12!")
        return 0, float('inf')

    # Compute E[N^2] by overlap
    # E[N^2] = sum_{S,S' valid} 1
    #        = sum_s #{valid pairs with overlap s}
    overlap_hist = defaultdict(int)
    for i in range(len(valid_list)):
        for j in range(len(valid_list)):
            s = len(valid_list[i] & valid_list[j])
            overlap_hist[s] += 1

    E_N2 = sum(overlap_hist.values())  # = N_val^2
    assert E_N2 == N_val ** 2

    # But for the probabilistic formulation:
    # E[N^2] = sum_{S,S'} Pr[both valid]
    # where (S, S') ranges over all ORDERED pairs of k-subsets
    # = sum_{S,S'} X_S * X_{S'}
    # This equals N_val^2 since X_S is deterministic.
    #
    # The Paley-Zygmund ratio is:
    # E[N^2] / E[N]^2 where E is over the RANDOM choice
    #
    # Wait -- N = sum_S X_S is deterministic for fixed D11.
    # So Paley-Zygmund doesn't apply directly.
    #
    # The correct approach: think of N as a random variable over choice of D11.
    # OR: use the "conditional second moment" formulation.
    #
    # Actually, the standard approach is:
    # N is the number of valid D12 for a fixed D11.
    # We already know N > 0 iff there exists a valid D12.
    # The first moment E_{D12}[1[valid]] = N / C(p-1,k) = p_0.
    # The second moment E_{D12,D12'}[1[both valid]] = ???
    #
    # For random D12: E[N] = C(p-1,k) * p_0.
    # By second moment: Pr[N >= 1] >= E[N]^2 / E[N^2]
    # where E[N^2] = E[(sum X_S)^2] = sum_{S,S'} E[X_S X_{S'}]
    # But X_S is deterministic! X_S = 1[S valid for fixed D11].
    #
    # OK so the second moment method doesn't make sense for FIXED D11.
    # It makes sense for RANDOM D11, where X_S depends on both D11 and S.

    # Let me instead compute p(s) = Pr[S valid AND S' valid | overlap s]
    # where S, S' are random k-subsets with |S cap S'| = s.
    # This is: #{valid pairs with overlap s} / #{all pairs with overlap s}

    print(f"\n  Overlap analysis:")
    print(f"  {'s':>4} {'#valid_pairs':>14} {'#all_pairs':>14} {'p(s)':>12} {'p(s)/p0^2':>12}")

    E_s = k * k / (p - 1)
    total_ratio = 0
    total_weight = 0

    for s in sorted(overlap_hist.keys()):
        n_valid = overlap_hist[s]
        # Number of all ordered pairs with overlap s:
        # Choose S: C(p-1,k). Choose which s elements of S overlap: C(k,s).
        # Choose remaining k-s elements of S' from outside S: C(p-1-k, k-s).
        n_all = total * comb(k, s) * comb(p - 1 - k, k - s)
        p_s = n_valid / n_all if n_all > 0 else 0
        ratio = p_s / p_0**2 if p_0 > 0 else 0

        # Weight in the hypergeometric sum
        weight = comb(k, s) * comb(p - 1 - k, k - s) / total
        total_ratio += weight * ratio
        total_weight += weight

        print(f"  {s:4d} {n_valid:14d} {n_all:14d} {p_s:12.6f} {ratio:12.4f}")

    print(f"\n  E_s[p(s)/p0^2] = {total_ratio:.4f}")
    print(f"  (This is E[N^2]/E[N]^2 in the random formulation)")
    print(f"  Paley-Zygmund: Pr[N >= 1] >= 1/{total_ratio:.4f} = {1/total_ratio:.6f}")
    print(f"  E[s] (hypergeometric mean) = {E_s:.1f}")

    return N_val, total_ratio


def main():
    for p in [11, 19]:
        if p % 4 != 3:
            continue

        print("=" * 70)
        print(f"p = {p}")
        print("=" * 70)

        best_D11, best_max_A, good_D11s = find_best_d11(p)
        E_A = (p + 1) / 4

        print(f"\nBest D11: max_A = {best_max_A}, E[A] = {E_A}")
        if good_D11s:
            print(f"Found {len(good_D11s)} perfectly flat D11 (max A = {int(E_A)})")

        # Try ALL D11 that have valid D12, starting with flattest
        if good_D11s:
            print(f"\nUsing first perfectly flat D11:")
            for D11 in good_D11s[:2]:
                result = second_moment_exhaustive(p, D11)
                if result[0] > 0:
                    break
                print()
        else:
            # Try best and also try all D11s to find one with valid D12
            print(f"\nBest D11 has max_A = {best_max_A} > E[A] = {E_A}")
            # Search all D11s for ones with valid D12
            pairs = []
            seen = set()
            for d in range(1, p):
                if d not in seen:
                    pairs.append((d, p - d))
                    seen.add(d)
                    seen.add(p - d)
            num_pairs = len(pairs)
            n = (p + 1) // 2
            reps_needed = n // 2

            found = False
            for combo in combinations(range(num_pairs), reps_needed):
                D11 = []
                for idx in combo:
                    d1, d2 = pairs[idx]
                    D11.append(d1)
                    if d1 != d2:
                        D11.append(d2)
                if len(D11) != n:
                    continue
                D11 = sorted(D11)
                D11_set = set(D11)
                A = compute_A(D11_set, p)
                max_A = max(A[d] for d in range(1, p))
                if max_A <= int(E_A):
                    print(f"\nTrying flat D11 with max_A = {max_A}:")
                    result = second_moment_exhaustive(p, D11)
                    if result[0] > 0:
                        found = True
                        break
            if not found:
                print(f"\nNo flat D11 with valid D12 found. Trying best D11:")
                second_moment_exhaustive(p, best_D11)
        print()


if __name__ == '__main__':
    main()
