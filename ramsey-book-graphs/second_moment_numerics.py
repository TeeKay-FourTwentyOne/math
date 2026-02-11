#!/usr/bin/env python3
"""
Compute E[N^2] for the second moment method.

For fixed A-flat D11, N = #{valid D12}.
E[N^2] = sum_{D12, D12'} Pr[both valid]
       = sum_s W(s) * Pr[both valid | overlap s]

where W(s) = C(p-1,k) * C(k,s) * C(p-1-k,k-s) is the number of
ordered pairs (D12, D12') with |D12 cap D12'| = s + 1 (the +1 for {0}).

Actually: D12 = {0} union S, D12' = {0} union S', so the overlap in the
random parts is s = |S cap S'|, ranging from max(0, 2k-p+1) to k.
"""

import numpy as np
from math import comb, log2, sqrt
from itertools import combinations
import sys


def compute_B(S, p):
    """Compute B(d) for all d given S subset of {1,...,p-1}.
    B(d) = |{(a,b) in S x S : a - b = d mod p}| for d != 0.
    """
    S_set = set(S)
    B = {}
    for d in range(1, p):
        count = 0
        for a in S_set:
            b = (a - d) % p
            if b in S_set and b != 0:  # S is subset of {1,...,p-1}
                count += 1
        B[d] = count
    return B


def check_valid(B, T, p):
    """Check if B(d) <= T(d) for all d = 1,...,p-1."""
    for d in range(1, p):
        if B[d] > T[d]:
            return False
    return True


def compute_A(D11, p):
    """Compute A(d) = Delta(D11, D11, d) for all d."""
    D11_set = set(D11)
    A = {}
    for d in range(1, p):
        count = 0
        for a in D11_set:
            b = (a - d) % p
            if b in D11_set:
                count += 1
        A[d] = count
    return A


def exhaustive_second_moment(p):
    """Exhaustive computation of E[N], E[N^2] for small p."""
    n = (p + 1) // 2
    k = (p - 3) // 2  # |S| = |D12| - 1

    # Find a good (A-flat) D11
    # Use QR for p = 3 mod 4
    QR = set()
    for x in range(1, p):
        QR.add((x * x) % p)

    # D11 = QR union {0}... no, D11 is subset of {1,...,p-1}
    # Try D11 = QR (quadratic residues)
    D11 = sorted(QR)
    if len(D11) != n:
        # For p = 3 mod 4, |QR| = (p-1)/2. We need |D11| = (p+1)/2.
        # D11 should be symmetric (d in D11 iff p-d in D11)
        # Let's enumerate and find best D11
        pass

    # For small p, just find the best A-flat D11 by enumeration
    # Generate all symmetric subsets of size n
    # A symmetric subset: for each pair {d, p-d}, either both in or both out
    pairs = []
    seen = set()
    for d in range(1, p):
        if d not in seen:
            pairs.append((d, p - d))
            seen.add(d)
            seen.add(p - d)

    num_pairs = len(pairs)  # = (p-1)/2
    reps_needed = n // 2  # number of pairs to include (n is even since p = 3 mod 4)

    best_D11 = None
    best_max_A = p

    for combo in combinations(range(num_pairs), reps_needed):
        D11 = []
        for idx in combo:
            d1, d2 = pairs[idx]
            D11.append(d1)
            if d1 != d2:
                D11.append(d2)
        D11 = sorted(D11)
        if len(D11) != n:
            continue

        A = compute_A(D11, p)
        max_A = max(A.values())
        if max_A < best_max_A:
            best_max_A = max_A
            best_D11 = D11

    D11 = best_D11
    A = compute_A(D11, p)
    print(f"p = {p}, n = {n}, k = {k}")
    print(f"Best D11: {D11}")
    print(f"max A(d) = {max(A.values())}, E[A] = {(p+1)/4:.1f}")

    # Compute thresholds
    D11_set = set(D11)
    T = {}
    for d in range(1, p):
        if d in D11_set:
            T[d] = (p - 3) // 2 - A[d]
        else:
            T[d] = (p + 3) // 2 - A[d]

    # Enumerate all k-subsets S of {1,...,p-1}
    elements = list(range(1, p))
    valid_count = 0
    total = 0
    valid_sets = []

    for S in combinations(elements, k):
        total += 1
        B = compute_B(S, p)
        if check_valid(B, T, p):
            valid_count += 1
            valid_sets.append(set(S))

    E_N = valid_count
    print(f"E[N] = {E_N} out of C({p-1},{k}) = {comb(p-1,k)}")
    print(f"Pr[valid] = {E_N/comb(p-1,k):.6f}")

    if E_N == 0:
        print("No valid D12 found!")
        return

    # Compute E[N^2] by checking all pairs of valid D12
    E_N2 = 0
    overlap_counts = {}
    for i in range(len(valid_sets)):
        for j in range(len(valid_sets)):
            s = len(valid_sets[i] & valid_sets[j])
            E_N2 += 1
            overlap_counts[s] = overlap_counts.get(s, 0) + 1

    print(f"E[N^2] = {E_N2}")
    print(f"E[N]^2 = {E_N**2}")
    print(f"E[N^2]/E[N]^2 = {E_N2/E_N**2:.4f}")
    print(f"Paley-Zygmund: Pr[N>0] >= E[N]^2/E[N^2] = {E_N**2/E_N2:.6f}")
    print()

    # Overlap distribution
    print("Overlap distribution among valid pairs:")
    print(f"  s : count : fraction")
    for s in sorted(overlap_counts.keys()):
        print(f"  {s:3d}: {overlap_counts[s]:6d} : {overlap_counts[s]/E_N2:.4f}")

    # Expected overlap for random pairs
    E_s = k * k / (p - 1)
    print(f"\nExpected overlap (random): {E_s:.1f}")

    return E_N, E_N2, E_N2 / E_N**2


def main():
    for p in [11, 19]:
        if p % 4 != 3:
            continue
        print("=" * 60)
        result = exhaustive_second_moment(p)
        print()

    # p=23 might be too slow for exhaustive -- check
    p = 23
    k = (p - 3) // 2  # = 10
    total_subsets = comb(22, 10)
    print(f"p=23: C(22,10) = {total_subsets} -- {'feasible' if total_subsets < 1e7 else 'too large'}")


if __name__ == '__main__':
    main()
