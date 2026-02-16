#!/usr/bin/env python3
"""
Verify: is P(b) >= Q(b) for ALL achievable B-profiles on the hyperplane?
If so, c₀ >= Pr_Q[H|E] for any constraint event E.
This is the KEY step for the proof.

Also check: does supp(P) cover all of E ∩ H?
"""

import numpy as np
from math import comb, log2
from itertools import combinations
from collections import Counter
import json
import time


def verify(p):
    k = (p - 1) // 2
    n = (p + 1) // 2
    R = k
    reps = list(range(1, R + 1))
    S = k * (k - 1) // 2
    total = comb(p, k)

    t0 = time.time()
    print(f"\np={p}: k={k}, R={R}, S={S}, C(p,k)={total}")

    # Enumerate all D12, compute B-vectors
    joint_counts = Counter()  # B-vector -> count
    marginal_counts = [Counter() for _ in range(R)]

    for D12_tuple in combinations(range(p), k):
        D12_set = set(D12_tuple)
        bvec = tuple(
            sum(1 for i in D12_set if (i + d) % p in D12_set)
            for d in reps
        )
        joint_counts[bvec] += 1
        for j, b in enumerate(bvec):
            marginal_counts[j][b] += 1

    num_profiles = len(joint_counts)
    print(f"  Distinct B-profiles on H (=supp(P)): {num_profiles}")

    # Verify all on hyperplane
    for bvec, cnt in joint_counts.items():
        assert sum(bvec) == S, f"Profile {bvec} has sum {sum(bvec)} != {S}"

    # Compute P(b) and Q(b) for each profile
    min_ratio = float('inf')
    max_ratio = 0
    ratios = []
    below_1_count = 0

    for bvec, cnt in joint_counts.items():
        P_b = cnt / total
        Q_b = 1.0
        for j in range(R):
            Q_b *= marginal_counts[j][bvec[j]] / total
        ratio = P_b / Q_b
        ratios.append(ratio)
        min_ratio = min(min_ratio, ratio)
        max_ratio = max(max_ratio, ratio)
        if ratio < 1.0:
            below_1_count += 1

    print(f"\n  P/Q ratios on supp(P) (= on H ∩ achievable):")
    print(f"    min P/Q = {min_ratio:.8f}")
    print(f"    max P/Q = {max_ratio:.8f}")
    print(f"    mean P/Q = {np.mean(ratios):.6f}")
    print(f"    Profiles with P/Q < 1: {below_1_count}/{num_profiles}")
    print(f"    min P/Q > 1: {'YES ✓' if min_ratio > 1.0 else 'NO ✗'}")

    if min_ratio > 1:
        # min P/Q = Pr_Q[H] × min P/Q_H
        # Since Q_H(b) = Q(b)/Pr_Q[H], min P/Q_H = min(P/Q) / Pr_Q[H]
        # Actually min P(b)/Q(b) directly gives us what we need
        print(f"\n  *** P(b) >= Q(b) for ALL achievable b on H ***")
        print(f"  This means: c₀_reps >= Pr_Q[H|E] for ANY constraint event E")

    # Now check: how many integer points on E∩H exist but are NOT achievable?
    # For each working D11, count points in E∩H with vs without P-support
    half = (p - 1) // 2
    pairs = [(i + 1, p - (i + 1)) for i in range(half)]
    n_choose = (p + 1) // 4
    support = sorted(marginal_counts[0].keys())
    max_b = max(support)

    # Check a few D11
    d11_list = []
    for chosen in combinations(range(half), n_choose):
        d11 = set()
        for i in chosen:
            d11.add(pairs[i][0])
            d11.add(pairs[i][1])
        d11_list.append(frozenset(d11))

    # Check coverage for working D11 (take first few)
    print(f"\n  Checking supp(P) coverage of E∩H for working D11:")
    for d11 in d11_list[:5]:
        A = {}
        for d in range(1, p):
            A[d] = sum(1 for x in d11 if (x + d) % p in d11)

        T_reps = []
        for d in reps:
            if d in d11:
                T_reps.append((n - 2) - A[d])
            else:
                T_reps.append((n + 1) - A[d])

        if any(t < 0 for t in T_reps):
            continue

        # Count profiles in E∩H∩supp(P) and also in E∩H
        in_E_and_H_and_supp = 0
        N_valid = 0  # subsets
        for bvec, cnt in joint_counts.items():
            if all(bvec[j] <= T_reps[j] for j in range(R)):
                in_E_and_H_and_supp += 1
                N_valid += cnt

        if N_valid == 0:
            continue

        # Count integer lattice points in E∩H (may not be in supp(P))
        # This is computationally expensive for large R...
        # Use the convolution approach to count integer points
        from itertools import product as iproduct
        if p <= 19:
            # Brute force for small R
            in_E_and_H_total = 0
            for bvec in iproduct(*(range(T_reps[j]+1) for j in range(R))):
                if sum(bvec) == S:
                    in_E_and_H_total += 1

            print(f"    D11={sorted(d11)[:5]}...: E∩H lattice points={in_E_and_H_total}, "
                  f"of which in supp(P)={in_E_and_H_and_supp}, N_valid={N_valid}")
            if in_E_and_H_total == in_E_and_H_and_supp:
                print(f"    *** supp(P) COVERS all of E∩H! ***")
        else:
            print(f"    D11={sorted(d11)[:5]}...: supp(P)∩E∩H has {in_E_and_H_and_supp} profiles, "
                  f"N_valid={N_valid}")
        break  # Just check first working D11

    t1 = time.time()
    print(f"\n  ({t1-t0:.1f}s)")

    return {
        'p': p,
        'num_profiles': num_profiles,
        'min_PQ': min_ratio,
        'max_PQ': max_ratio,
        'below_1': below_1_count,
        'min_PQ_gt_1': min_ratio > 1.0
    }


def main():
    results = []
    for p in [11, 19, 23]:
        r = verify(p)
        results.append(r)

    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    for r in results:
        print(f"p={r['p']:2d}: min P/Q = {r['min_PQ']:.6f}, max = {r['max_PQ']:.4f}, "
              f"P/Q≥1: {'✓' if r['min_PQ_gt_1'] else '✗'}, "
              f"below_1: {r['below_1']}/{r['num_profiles']}")

    with open('ramsey-book-graphs/min_pq_results.json', 'w') as f:
        json.dump(results, f, indent=2)


if __name__ == '__main__':
    main()
