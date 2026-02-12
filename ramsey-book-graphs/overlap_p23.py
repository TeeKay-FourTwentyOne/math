#!/usr/bin/env python3
"""
Overlap decomposition at p=23 using batch FFT for speed.
Decomposes E[N²] by overlap class, separating same-orbit vs cross-orbit pairs.
"""

import numpy as np
from itertools import combinations
from collections import defaultdict
import time
import sys
import os
import json
from math import comb

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ramsey_core import Delta, Sigma, verify_construction, BlockCirculantGraph

P = 23
N_BOOK = (P + 1) // 2  # 12
K = (P - 3) // 2  # 10


def batch_autocorr(indicator_matrix):
    """Compute Delta(S,S,d) for all sets via FFT. Shape: (num_sets, p)."""
    fft = np.fft.fft(indicator_matrix, axis=1)
    ac = np.fft.ifft(np.abs(fft) ** 2, axis=1).real
    return np.round(ac).astype(int)


def autocorr_single(indicator, p):
    """Compute Delta(S,S,d) for a single set."""
    fft = np.fft.fft(indicator)
    ac = np.fft.ifft(np.abs(fft) ** 2).real
    return np.round(ac).astype(int)


def full_validity_check(d11_set, d12_set, p):
    """Full constraint check including V1V2 (slower but exact)."""
    n = (p + 1) // 2
    m = p
    d22_set = set(range(1, p)) - set(d11_set)
    G = BlockCirculantGraph(n=n, D11=set(d11_set), D12=set(d12_set),
                            D22=d22_set)
    result = verify_construction(G)
    return result.valid


def main():
    print(f"Overlap decomposition at p={P}, n={N_BOOK}, k={K}")
    t0 = time.time()

    # Load working D11 from enumeration data
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    with open(os.path.join(base_dir, 'enumeration_data', 'enumeration_p23.json')) as f:
        data = json.load(f)

    all_d11_data = data['per_d11']
    num_d11 = len(all_d11_data)
    working = [(i, set(e['D11']), e['num_valid_d12'])
               for i, e in enumerate(all_d11_data) if e['num_valid_d12'] > 0]
    print(f"Working D11: {len(working)} out of {num_d11}")

    # Build D12 indicator matrix and batch autocorrelation
    print(f"Building D12 matrix (C({P-1},{K}) = {comb(P-1,K)})...")
    t1 = time.time()
    d12_list = []
    d12_matrix = []
    for s in combinations(range(1, P), K):
        ind = np.zeros(P, dtype=np.float64)
        ind[0] = 1.0
        for x in s:
            ind[x] = 1.0
        d12_matrix.append(ind)
        d12_list.append(frozenset([0] + list(s)))
    d12_matrix = np.array(d12_matrix)  # (num_d12, P)
    num_d12 = len(d12_list)
    print(f"  {num_d12} D12 built in {time.time()-t1:.1f}s")

    print("Computing batch autocorrelation (FFT)...")
    t1 = time.time()
    B_matrix = batch_autocorr(d12_matrix)  # (num_d12, P)
    # B_neg[i,d] = B_i(P-d) for V2V2 constraints
    neg_idx = np.array([(P - d) % P for d in range(P)])
    B_neg = B_matrix[:, neg_idx]
    print(f"  Done in {time.time()-t1:.1f}s")

    # For each working D11, find valid D12 using vectorized V1V1 + full check
    print(f"\nFinding valid D12 for {len(working)} working D11...")
    valid_d12_by_d11 = {}

    for idx, (d11_idx, d11_set, expected_N) in enumerate(working):
        t1 = time.time()

        # Precompute A(d) and C(d) for this D11
        d11_ind = np.zeros(P, dtype=np.float64)
        for x in d11_set:
            d11_ind[x] = 1.0
        A = autocorr_single(d11_ind, P)

        d22_set = set(range(1, P)) - d11_set
        d22_ind = np.zeros(P, dtype=np.float64)
        for x in d22_set:
            d22_ind[x] = 1.0
        C = autocorr_single(d22_ind, P)

        d11_mask = np.zeros(P, dtype=bool)
        for x in d11_set:
            d11_mask[x] = True
        d22_mask = np.zeros(P, dtype=bool)
        for x in d22_set:
            d22_mask[x] = True

        n = N_BOOK

        # V1V1 red: A(d)+B(d) <= n-2 for d in D11
        F = A[np.newaxis, :] + B_matrix  # (num_d12, P)
        if d11_mask.any():
            v1v1_ok = F[:, d11_mask].max(axis=1) <= n - 2
        else:
            v1v1_ok = np.ones(num_d12, dtype=bool)

        # V2V2 blue (tightest): C(d)+B(p-d) <= n-3 for d in D11
        G = C[np.newaxis, :] + B_neg  # (num_d12, P)
        if d11_mask.any():
            v2v2_blue_ok = G[:, d11_mask].max(axis=1) <= n - 3
        else:
            v2v2_blue_ok = np.ones(num_d12, dtype=bool)

        # V2V2 red: C(d)+B(p-d) <= n-2 for d in D22
        if d22_mask.any():
            v2v2_red_ok = G[:, d22_mask].max(axis=1) <= n - 2
        else:
            v2v2_red_ok = np.ones(num_d12, dtype=bool)

        # V1V1 blue: A(d)+B(d) <= n+1 for d in D22 (slack, usually auto-satisfied)
        if d22_mask.any():
            v1v1_blue_ok = F[:, d22_mask].max(axis=1) <= n + 1
        else:
            v1v1_blue_ok = np.ones(num_d12, dtype=bool)

        # Combine V1V1 + V2V2 fast filter
        fast_ok = v1v1_ok & v2v2_blue_ok & v2v2_red_ok & v1v1_blue_ok
        candidates = np.where(fast_ok)[0]

        # Full validation (including V1V2) for candidates
        valid_indices = []
        for ci in candidates:
            if full_validity_check(d11_set, d12_list[ci], P):
                valid_indices.append(ci)

        valid_d12_by_d11[d11_idx] = valid_indices
        t2 = time.time()
        passed_fast = len(candidates)
        passed_full = len(valid_indices)
        print(f"  D11[{d11_idx}]: fast={passed_fast}, full={passed_full} "
              f"(expected {expected_N}) [{t2-t1:.1f}s]")
        assert passed_full == expected_N, f"Mismatch: {passed_full} != {expected_N}"

    # =============================================
    # B-PROFILE ANALYSIS
    # =============================================
    print(f"\n{'='*60}")
    print("B-PROFILE ORBIT ANALYSIS")
    print(f"{'='*60}")

    d12_B_profiles = {}
    for d11_idx, d11_set, expected_N in working:
        for vi in valid_d12_by_d11[d11_idx]:
            if vi not in d12_B_profiles:
                d12_B_profiles[vi] = tuple(B_matrix[vi, 1:])  # B(1),...,B(p-1)

    for d11_idx, d11_set, expected_N in working:
        vis = valid_d12_by_d11[d11_idx]
        bp_groups = defaultdict(list)
        for vi in vis:
            bp_groups[d12_B_profiles[vi]].append(vi)

        n_orbits = len(bp_groups)
        print(f"\nD11[{d11_idx}] (N={expected_N}): {n_orbits} B-profile orbits")
        for bp, members in sorted(bp_groups.items(), key=lambda x: -len(x[1])):
            print(f"  size={len(members)}, B[:6]={list(bp[:6])}")

    # =============================================
    # OVERLAP DECOMPOSITION
    # =============================================
    print(f"\n{'='*60}")
    print("OVERLAP DECOMPOSITION OF E[N²]")
    print(f"{'='*60}")

    EN = sum(e['num_valid_d12'] for e in all_d11_data) / num_d11
    EN2 = sum(e['num_valid_d12']**2 for e in all_d11_data) / num_d11

    overlap_total = defaultdict(float)
    overlap_same = defaultdict(float)
    overlap_cross = defaultdict(float)
    overlap_count = defaultdict(int)

    for d11_idx, d11_set, expected_N in working:
        vis = valid_d12_by_d11[d11_idx]
        for vi_a in vis:
            bp_a = d12_B_profiles[vi_a]
            for vi_b in vis:
                bp_b = d12_B_profiles[vi_b]
                s = len(d12_list[vi_a] & d12_list[vi_b])
                overlap_total[s] += 1.0 / num_d11
                overlap_count[s] += 1
                if bp_a == bp_b:
                    overlap_same[s] += 1.0 / num_d11
                else:
                    overlap_cross[s] += 1.0 / num_d11

    print(f"\nE[N] = {EN:.4f}, E[N²] = {EN2:.4f}, ratio = {EN2/EN**2:.4f}")
    print(f"\n{'s':>5} {'Total':>10} {'Same-orb':>10} {'Cross-orb':>10} "
          f"{'#Trips':>8} {'Cross%':>7}")
    print(f"{'-'*5} {'-'*10} {'-'*10} {'-'*10} {'-'*8} {'-'*7}")

    for s in sorted(overlap_total.keys()):
        t = overlap_total[s]
        sa = overlap_same.get(s, 0)
        cr = overlap_cross.get(s, 0)
        cnt = overlap_count[s]
        cp = cr / t * 100 if t > 0 else 0
        print(f"{s:>5} {t:>10.2f} {sa:>10.2f} {cr:>10.2f} {cnt:>8} {cp:>6.1f}%")

    total_c = sum(overlap_total.values())
    total_s = sum(overlap_same.values())
    total_x = sum(overlap_cross.values())
    print(f"{'TOT':>5} {total_c:>10.2f} {total_s:>10.2f} {total_x:>10.2f}")
    print(f"\nVerify: total={total_c:.2f}, E[N²]={EN2:.2f}")
    print(f"Same-orbit: {total_s/total_c:.1%}, Cross-orbit: {total_x/total_c:.1%}")

    # =============================================
    # P₂/(P₁·P₁') RATIOS
    # =============================================
    print(f"\n{'='*60}")
    print("P₂/(P₁·P₁') BY ORBIT TYPE AND OVERLAP")
    print(f"{'='*60}")

    # Compute P₁
    P1 = defaultdict(float)
    for d11_idx, _, _ in working:
        for vi in valid_d12_by_d11[d11_idx]:
            P1[vi] += 1.0 / num_d11

    # Compute P₂
    P2 = defaultdict(float)
    for d11_idx, _, _ in working:
        vis = valid_d12_by_d11[d11_idx]
        for a in vis:
            for b in vis:
                P2[(a, b)] += 1.0 / num_d11

    # Group ratios
    ratios_same = defaultdict(list)
    ratios_cross = defaultdict(list)

    for (a, b), p2v in P2.items():
        s = len(d12_list[a] & d12_list[b])
        p1p = P1[a] * P1[b]
        r = p2v / p1p if p1p > 0 else float('inf')
        if d12_B_profiles[a] == d12_B_profiles[b]:
            ratios_same[s].append(r)
        else:
            ratios_cross[s].append(r)

    print(f"\nSame-orbit P₂/(P₁·P₁'):")
    print(f"  {'s':>5} {'N':>6} {'Mean':>8} {'Min':>8} {'Max':>8}")
    for s in sorted(ratios_same.keys()):
        r = ratios_same[s]
        print(f"  {s:>5} {len(r):>6} {np.mean(r):>8.2f} {np.min(r):>8.2f} "
              f"{np.max(r):>8.2f}")

    if ratios_cross:
        print(f"\nCross-orbit P₂/(P₁·P₁'):")
        print(f"  {'s':>5} {'N':>6} {'Mean':>8} {'Min':>8} {'Max':>8}")
        for s in sorted(ratios_cross.keys()):
            r = ratios_cross[s]
            print(f"  {s:>5} {len(r):>6} {np.mean(r):>8.2f} {np.min(r):>8.2f} "
                  f"{np.max(r):>8.2f}")

    # =============================================
    # MARGINAL ANALYSIS: CROSS-ORBIT PAIR
    # =============================================
    print(f"\n{'='*60}")
    print("PER-POSITION MARGINALS: CROSS-ORBIT EXAMPLE")
    print(f"{'='*60}")

    # Find D11 with most B-profile orbits
    best = None
    for d11_idx, d11_set, expected_N in working:
        vis = valid_d12_by_d11[d11_idx]
        bps = set(d12_B_profiles[vi] for vi in vis)
        if best is None or len(bps) > best[0]:
            best = (len(bps), d11_idx, d11_set, expected_N, vis)

    n_orbits, d11_idx, d11_set, expected_N, vis = best
    bp_groups = defaultdict(list)
    for vi in vis:
        bp_groups[d12_B_profiles[vi]].append(vi)
    bp_list = sorted(bp_groups.keys(), key=lambda x: -len(bp_groups[x]))

    print(f"\nD11[{d11_idx}] (N={expected_N}, {n_orbits} orbits)")

    # Pick one D12 from each of the first 2 orbits
    vi_a = bp_groups[bp_list[0]][0]
    vi_b = bp_groups[bp_list[1]][0]
    B_a = np.array(d12_B_profiles[vi_a])
    B_b = np.array(d12_B_profiles[vi_b])
    B_max = np.maximum(B_a, B_b)
    tight = B_max - B_a

    print(f"D12_a (orbit 0, size {len(bp_groups[bp_list[0]])}): B={list(B_a[:11])}")
    print(f"D12_b (orbit 1, size {len(bp_groups[bp_list[1]])}): B={list(B_b[:11])}")
    print(f"max(B,B'): {list(B_max[:11])}")
    print(f"Tightening: {list(tight[:11])}")
    print(f"  {np.sum(tight>0)}/{P-1} positions tightened, max_tight={np.max(tight)}, "
          f"sum_tight={np.sum(tight)}")

    # Compute per-position marginals (V1V1+V2V2 only)
    reps = list(range(1, (P + 1) // 2))  # d=1,...,11

    marginals_j = np.zeros(P)  # joint
    marginals_a = np.zeros(P)  # single a
    marginals_b = np.zeros(P)  # single b

    for entry in all_d11_data:
        d11_e = set(entry['D11'])
        d22_e = set(range(1, P)) - d11_e
        d11_i = np.zeros(P, dtype=np.float64)
        d22_i = np.zeros(P, dtype=np.float64)
        for x in d11_e: d11_i[x] = 1.0
        for x in d22_e: d22_i[x] = 1.0
        A_e = autocorr_single(d11_i, P)
        C_e = autocorr_single(d22_i, P)

        for d in range(1, P):
            pd = P - d
            Ba, Bb = int(B_a[d-1]), int(B_b[d-1])
            Bm = max(Ba, Bb)
            Ba_pd, Bb_pd = int(B_a[pd-1]), int(B_b[pd-1])
            Bm_pd = max(Ba_pd, Bb_pd)
            Ad, Cd = int(A_e[d]), int(C_e[d])
            n = N_BOOK

            if d in d11_e:
                ok_a = (Ad + Ba <= n-2) and (Cd + Ba_pd <= n-3)
                ok_b = (Ad + Bb <= n-2) and (Cd + Bb_pd <= n-3)
                ok_j = (Ad + Bm <= n-2) and (Cd + Bm_pd <= n-3)
            else:
                ok_a = (Cd + Ba_pd <= n-2) and (Ad + Ba <= n+1)
                ok_b = (Cd + Bb_pd <= n-2) and (Ad + Bb <= n+1)
                ok_j = (Cd + Bm_pd <= n-2) and (Ad + Bm <= n+1)

            marginals_a[d] += ok_a
            marginals_b[d] += ok_b
            marginals_j[d] += ok_j

    marginals_a /= num_d11
    marginals_b /= num_d11
    marginals_j /= num_d11

    # Products over representatives
    prod_a = np.prod(marginals_a[reps])
    prod_b = np.prod(marginals_b[reps])
    prod_j = np.prod(marginals_j[reps])

    p2_ab = P2.get((vi_a, vi_b), 0)
    p1_a = P1[vi_a]
    p1_b = P1[vi_b]

    print(f"\nPer-rep marginals (d=1..11):")
    print(f"  Joint:    {[f'{marginals_j[d]:.3f}' for d in reps]}")
    print(f"  Single A: {[f'{marginals_a[d]:.3f}' for d in reps]}")
    print(f"  Single B: {[f'{marginals_b[d]:.3f}' for d in reps]}")
    print(f"  Tighten:  {[f'{marginals_a[d]-marginals_j[d]:.3f}' for d in reps]}")

    print(f"\nProducts over representatives:")
    print(f"  ∏ joint   = {prod_j:.6e}")
    print(f"  ∏ single A= {prod_a:.6e}")
    print(f"  ∏ single B= {prod_b:.6e}")
    print(f"  Actual P₂ = {p2_ab:.6e}")
    print(f"  P₁(A)     = {p1_a:.6e}")
    print(f"  P₁(B)     = {p1_b:.6e}")
    print(f"  P₁·P₁'   = {p1_a*p1_b:.6e}")

    if prod_j > 0:
        print(f"  c₀² = P₂/∏joint = {p2_ab/prod_j:.4f}")
    if prod_a > 0:
        print(f"  c₀(A) = P₁/∏A   = {p1_a/prod_a:.4f}")
    if prod_b > 0:
        print(f"  c₀(B) = P₁/∏B   = {p1_b/prod_b:.4f}")
    if p1_a * p1_b > 0:
        print(f"  P₂/(P₁·P₁')     = {p2_ab/(p1_a*p1_b):.4f}")

    # Compare same-orbit and cross-orbit c₀²
    print(f"\n--- Same-orbit pair for comparison ---")
    vi_a2 = bp_groups[bp_list[0]][1] if len(bp_groups[bp_list[0]]) > 1 else vi_a
    p2_same = P2.get((vi_a, vi_a2), 0)
    print(f"  P₂(same-orbit) = {p2_same:.6e}")
    print(f"  P₂(cross-orbit) = {p2_ab:.6e}")
    print(f"  Ratio cross/same: {p2_ab/p2_same:.4f}" if p2_same > 0 else "  Same=0")

    t_end = time.time()
    print(f"\nTotal time: {t_end-t0:.1f}s")


if __name__ == '__main__':
    main()
