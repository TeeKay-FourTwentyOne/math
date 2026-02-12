#!/usr/bin/env python3
"""
Second moment overlap decomposition for R(B_{n-1}, B_n).

Decomposes E[N²] = Σ_{s} Σ_{|D12∩D12'|=s} P₂(D12,D12') by overlap class s,
where P₂(D12,D12') = Pr_{D11}[both D12 and D12' valid for random D11].

Also computes per-position marginals and their products, comparing to actual P₂
to measure the correlation loss structure.

Steps:
1. For fixed (D12, D12'), tightened threshold at d: T(d) - max(B(d), B'(d))
2. Per-position marginals → product baseline → ratio P₂/product = c₀²
3. Sum over overlap classes → show Σ = E[N]² × poly(p)
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
from ramsey_core import Delta, Sigma, BlockCirculantGraph, verify_construction


def gen_symmetric_d11(p):
    """Generate all symmetric D11 ⊂ {1,...,p-1} with |D11| = (p+1)/2."""
    n = (p + 1) // 2
    pairs = [(d, p - d) for d in range(1, (p + 1) // 2)]
    r = n // 2  # = (p+1)/4 pairs to choose
    result = []
    for chosen in combinations(range(len(pairs)), r):
        d11 = set()
        for i in chosen:
            d11.add(pairs[i][0])
            d11.add(pairs[i][1])
        result.append(frozenset(d11))
    return result


def gen_all_d12(p):
    """Generate all D12 = {0} ∪ S, S a k-subset of {1,...,p-1}, k=(p-3)/2."""
    k = (p - 3) // 2
    result = []
    for s in combinations(range(1, p), k):
        result.append(frozenset([0] + list(s)))
    return result


def check_valid(d11_set, d12_set, p):
    """Check validity of (D11, D12) using the full constraint model."""
    n = (p + 1) // 2
    m = p
    d22_set = set(range(1, p)) - set(d11_set)
    d12t = frozenset((-x) % m for x in d12_set)
    N_total = 2 * m
    d1 = len(d11_set) + len(d12_set)
    d2 = len(d22_set) + len(d12_set)

    for d in range(1, m):
        A_d = Delta(d11_set, d11_set, d, m)
        B_d = Delta(d12_set, d12_set, d, m)
        C_d = Delta(d22_set, d22_set, d, m)
        Bt_d = Delta(d12t, d12t, d, m)  # = B(p-d)

        if d in d11_set:
            # V1V1 red: A(d)+B(d) ≤ n-2
            if A_d + B_d > n - 2:
                return False
            # V2V2 blue (tightest): C(d)+B(p-d) ≤ n-3
            if C_d + Bt_d > n - 3:
                return False
        else:
            # V2V2 red: C(d)+Bt(d) ≤ n-2
            if C_d + Bt_d > n - 2:
                return False
            # V1V1 blue: A(d)+B(d) ≤ n+1
            if A_d + B_d > n + 1:
                return False

    # V1V2 constraints
    for d in range(m):
        sig = Sigma(d11_set, d12_set, d, m)
        dlt = Delta(d12_set, d22_set, d, m)
        X = sig + dlt

        if d in d12_set:
            if X > n - 2:
                return False
        else:
            blue_cross = (N_total - 2) - d1 - d2 + X
            if blue_cross > n - 1:
                return False

    return True


def compute_B_profile(d12_set, p):
    """Compute B(d) = Delta(D12, D12, d, p) for all d in {1,...,p-1}."""
    return {d: Delta(d12_set, d12_set, d, p) for d in range(1, p)}


def compute_per_position_marginal(all_d11, d12_set, d12p_set, p):
    """
    For each position d, compute marginal Pr_{D11}[constraint at d satisfied
    for both D12 and D12'].

    Returns: dict d -> marginal probability (over V1V1/V2V2 constraints only),
             and the single-D12 marginals for comparison.
    """
    n = (p + 1) // 2
    m = p
    num_d11 = len(all_d11)

    B = compute_B_profile(d12_set, p)
    Bp = compute_B_profile(d12p_set, p)

    # For each d, count how many D11 satisfy the V1V1/V2V2 constraint at d for both
    joint_marginal = {}
    single_marginal = {}  # for D12 alone
    single_marginal_p = {}  # for D12' alone

    for d in range(1, p):
        count_joint = 0
        count_single = 0
        count_single_p = 0

        for d11 in all_d11:
            d22 = set(range(1, p)) - set(d11)
            A_d = Delta(d11, d11, d, m)
            C_d = Delta(d22, d22, d, m)

            B_max = max(B[d], Bp[d])
            B_pd_max = max(B.get((p - d) % p or p, 0), Bp.get((p - d) % p or p, 0))
            # Handle d and p-d for B values
            pd = p - d  # since d in {1,...,p-1}, p-d is also in {1,...,p-1}
            B_pd_max = max(B[pd], Bp[pd])

            if d in d11:
                # V1V1 red + V2V2 blue (tightest)
                ok_joint = (A_d + B_max <= n - 2) and (C_d + B_pd_max <= n - 3)
                ok_single = (A_d + B[d] <= n - 2) and (C_d + B[pd] <= n - 3)
                ok_single_p = (A_d + Bp[d] <= n - 2) and (C_d + Bp[pd] <= n - 3)
            else:
                # d in D22: V2V2 red + V1V1 blue
                ok_joint = (C_d + B_pd_max <= n - 2) and (A_d + B_max <= n + 1)
                ok_single = (C_d + B[pd] <= n - 2) and (A_d + B[d] <= n + 1)
                ok_single_p = (C_d + Bp[pd] <= n - 2) and (A_d + Bp[d] <= n + 1)

            count_joint += ok_joint
            count_single += ok_single
            count_single_p += ok_single_p

        joint_marginal[d] = count_joint / num_d11
        single_marginal[d] = count_single / num_d11
        single_marginal_p[d] = count_single_p / num_d11

    return joint_marginal, single_marginal, single_marginal_p


def run_analysis(p):
    """Full overlap decomposition analysis for prime p."""
    print(f"\n{'='*70}")
    print(f"  SECOND MOMENT OVERLAP DECOMPOSITION: p = {p}")
    print(f"{'='*70}")

    n = (p + 1) // 2
    k = (p - 3) // 2
    m = p

    t0 = time.time()

    # Generate all D11 and D12
    print(f"\nGenerating D11 (symmetric, size {n}) and D12 (size {(p-1)//2})...")
    all_d11 = gen_symmetric_d11(p)
    all_d12 = gen_all_d12(p)
    num_d11 = len(all_d11)
    num_d12 = len(all_d12)
    print(f"  #D11 = {num_d11}, #D12 = {num_d12}")
    print(f"  #D12 pairs = {num_d12**2}")

    # Find valid D12 for each D11
    print(f"\nComputing validity matrix...")
    t1 = time.time()
    valid_d12_for = {}  # d11_idx -> list of d12_idx
    N_values = {}

    for i, d11 in enumerate(all_d11):
        valid = []
        for j, d12 in enumerate(all_d12):
            if check_valid(d11, d12, p):
                valid.append(j)
        valid_d12_for[i] = valid
        N_values[i] = len(valid)
        if (i + 1) % max(1, num_d11 // 10) == 0:
            print(f"  D11 {i+1}/{num_d11}: N={len(valid)}")

    t2 = time.time()
    print(f"  Validity check: {t2-t1:.2f}s")

    # Basic statistics
    EN = sum(N_values.values()) / num_d11
    EN2 = sum(v**2 for v in N_values.values()) / num_d11
    ratio = EN2 / EN**2 if EN > 0 else float('inf')
    pr_positive = sum(1 for v in N_values.values() if v > 0) / num_d11

    print(f"\n--- Basic Statistics ---")
    print(f"  E[N] = {EN:.4f}")
    print(f"  E[N²] = {EN2:.4f}")
    print(f"  E[N²]/E[N]² = {ratio:.4f}")
    print(f"  Pr[N>0] = {pr_positive:.4f}")

    # ============================================================
    # STEP 1: Overlap decomposition of E[N²]
    # ============================================================
    print(f"\n{'='*70}")
    print(f"  STEP 1: Overlap decomposition of E[N²]")
    print(f"{'='*70}")

    # For each working D11, enumerate pairs of valid D12 and bin by overlap
    overlap_contrib = defaultdict(float)   # s -> contribution to E[N²]
    overlap_pair_count = defaultdict(int)   # s -> number of (D11, D12, D12') triples

    for i in range(num_d11):
        V = valid_d12_for[i]
        for j1 in V:
            for j2 in V:
                s = len(all_d12[j1] & all_d12[j2])
                overlap_contrib[s] += 1.0 / num_d11
                overlap_pair_count[s] += 1

    # Expected number of pairs at each overlap (for reference)
    # H(t) = C(p-1,k) * C(k,t) * C(p-1-k, k-t) where t = s-1 (since 0 is always shared)
    print(f"\n  {'Overlap s':>10} {'Contrib to E[N²]':>18} {'#Triples':>10} "
          f"{'H(s)':>12} {'Avg P₂':>12} {'Frac of E[N²]':>15}")
    print(f"  {'-'*10} {'-'*18} {'-'*10} {'-'*12} {'-'*12} {'-'*15}")

    total_contrib = sum(overlap_contrib.values())
    for s in sorted(overlap_contrib.keys()):
        c = overlap_contrib[s]
        cnt = overlap_pair_count[s]
        # H(s): number of ordered pairs with overlap exactly s
        # D12 = {0} ∪ S, overlap = 1 + |S∩S'|, so t = s-1
        t = s - 1
        H_s = comb(k, t) * comb(p - 1 - k, k - t) * num_d12  # ordered pairs
        avg_p2 = c / H_s if H_s > 0 else 0
        frac = c / total_contrib
        print(f"  {s:>10} {c:>18.4f} {cnt:>10} {H_s:>12} {avg_p2:>12.2e} {frac:>15.4f}")

    print(f"\n  Total = {total_contrib:.4f} (should be {EN2:.4f})")

    # ============================================================
    # STEP 2: Compute P₁(D12) for each D12
    # ============================================================
    print(f"\n{'='*70}")
    print(f"  STEP 2: P₁ analysis and P₂/(P₁·P₁') ratios")
    print(f"{'='*70}")

    # P₁(D12) = Pr_{D11}[D12 valid] = #{D11: valid} / #D11
    P1 = np.zeros(num_d12)
    for i in range(num_d11):
        for j in valid_d12_for[i]:
            P1[j] += 1.0 / num_d11

    P1_mean = P1.mean()
    P1_nonzero = P1[P1 > 0]
    print(f"\n  P₁ statistics:")
    print(f"    Mean P₁ = E[N]/#D12 = {P1_mean:.6e}")
    print(f"    #D12 with P₁>0: {len(P1_nonzero)} / {num_d12}")
    if len(P1_nonzero) > 0:
        print(f"    P₁ range (nonzero): [{P1_nonzero.min():.6e}, {P1_nonzero.max():.6e}]")

    # ============================================================
    # STEP 2b: P₂/(P₁·P₁') ratio by overlap class
    # ============================================================

    # For each working D11, compute P₂ and P₁·P₁' for each pair
    ratio_by_overlap = defaultdict(list)  # s -> list of P₂/(P₁·P₁')
    p2_by_overlap = defaultdict(list)

    for i in range(num_d11):
        V = valid_d12_for[i]
        for j1 in V:
            for j2 in V:
                s = len(all_d12[j1] & all_d12[j2])
                # This triple contributes 1/num_d11 to P₂(j1, j2)
                # But P₂(j1, j2) gets contributions from ALL D11, not just this one
                pass  # We need to accumulate P₂ first

    # Compute full P₂ matrix for pairs that appear (sparse)
    P2_dict = defaultdict(float)  # (j1, j2) -> P₂
    for i in range(num_d11):
        V = valid_d12_for[i]
        for j1 in V:
            for j2 in V:
                P2_dict[(j1, j2)] += 1.0 / num_d11

    # Now group by overlap and compute ratios
    for (j1, j2), p2_val in P2_dict.items():
        s = len(all_d12[j1] & all_d12[j2])
        p1p1 = P1[j1] * P1[j2]
        if p1p1 > 0:
            ratio_by_overlap[s].append(p2_val / p1p1)
        p2_by_overlap[s].append(p2_val)

    print(f"\n  P₂/(P₁·P₁') ratio by overlap (pairs with nonzero P₂):")
    print(f"  {'Overlap s':>10} {'#Pairs':>8} {'Mean ratio':>12} "
          f"{'Min ratio':>12} {'Max ratio':>12} {'Mean P₂':>12}")
    print(f"  {'-'*10} {'-'*8} {'-'*12} {'-'*12} {'-'*12} {'-'*12}")

    for s in sorted(ratio_by_overlap.keys()):
        rats = ratio_by_overlap[s]
        p2s = p2_by_overlap[s]
        print(f"  {s:>10} {len(rats):>8} {np.mean(rats):>12.4f} "
              f"{np.min(rats):>12.4f} {np.max(rats):>12.4f} {np.mean(p2s):>12.6f}")

    # ============================================================
    # STEP 3: Per-position marginal analysis for representative pairs
    # ============================================================
    print(f"\n{'='*70}")
    print(f"  STEP 3: Per-position marginal analysis")
    print(f"{'='*70}")

    # Pick representative pairs at different overlap levels
    # We'll use pairs from the highest-N working D11
    best_d11_idx = max(N_values, key=N_values.get)
    V_best = valid_d12_for[best_d11_idx]
    print(f"\n  Using D11[{best_d11_idx}] with N={N_values[best_d11_idx]}")
    print(f"  D11 = {sorted(all_d11[best_d11_idx])}")

    # Find pairs at different overlap levels
    overlap_examples = {}
    for j1 in V_best:
        for j2 in V_best:
            if j1 >= j2:
                continue
            s = len(all_d12[j1] & all_d12[j2])
            if s not in overlap_examples:
                overlap_examples[s] = (j1, j2)

    print(f"\n  Representative pairs at each overlap level:")
    for s in sorted(overlap_examples.keys()):
        j1, j2 = overlap_examples[s]
        d12_1 = all_d12[j1]
        d12_2 = all_d12[j2]

        # Compute B profiles
        B1 = compute_B_profile(d12_1, p)
        B2 = compute_B_profile(d12_2, p)

        # Compute max B
        B_max = {d: max(B1[d], B2[d]) for d in range(1, p)}

        # Per-position marginals
        joint_marg, single_marg1, single_marg2 = compute_per_position_marginal(
            all_d11, d12_1, d12_2, p)

        # Products of marginals (over representative pairs {d, p-d})
        # Use positions d = 1,...,(p-1)/2 as representatives
        reps = list(range(1, (p + 1) // 2))

        prod_joint = 1.0
        prod_single1 = 1.0
        prod_single2 = 1.0
        for d in reps:
            prod_joint *= joint_marg[d]
            prod_single1 *= single_marg1[d]
            prod_single2 *= single_marg2[d]

        # Actual P₂ and P₁ values
        actual_p2 = P2_dict.get((j1, j2), 0) + P2_dict.get((j2, j1), 0)
        # P₂ should be symmetric; use the dict value
        actual_p2 = P2_dict.get((j1, j2), 0)
        actual_p1_1 = P1[j1]
        actual_p1_2 = P1[j2]

        print(f"\n  Overlap s={s}: D12[{j1}] ∩ D12[{j2}]")
        print(f"    |D12 ∩ D12'| = {s}")
        print(f"    B(d):  {[B1[d] for d in range(1, p)]}")
        print(f"    B'(d): {[B2[d] for d in range(1, p)]}")
        print(f"    max(B,B'): {[B_max[d] for d in range(1, p)]}")
        print(f"    Per-rep marginals (joint):  {[f'{joint_marg[d]:.4f}' for d in reps]}")
        print(f"    Per-rep marginals (single): {[f'{single_marg1[d]:.4f}' for d in reps]}")
        print(f"    Product of joint marginals:  {prod_joint:.6e}")
        print(f"    Product of single marginals: {prod_single1:.6e}")
        print(f"    Actual P₂ = {actual_p2:.6e}")
        print(f"    Actual P₁(D12) = {actual_p1_1:.6e}")
        print(f"    Actual P₁(D12') = {actual_p1_2:.6e}")
        print(f"    P₁·P₁' = {actual_p1_1 * actual_p1_2:.6e}")
        if prod_joint > 0:
            print(f"    c₀² = P₂/∏marginals = {actual_p2 / prod_joint:.4f}")
        if actual_p1_1 * actual_p1_2 > 0:
            print(f"    P₂/(P₁·P₁') = {actual_p2 / (actual_p1_1 * actual_p1_2):.4f}")
        if prod_single1 > 0:
            print(f"    c₀(D12) = P₁/∏marginals = {actual_p1_1 / prod_single1:.4f}")

    # ============================================================
    # STEP 4: Polynomial bound verification
    # ============================================================
    print(f"\n{'='*70}")
    print(f"  STEP 4: Polynomial bound verification")
    print(f"{'='*70}")

    # For each overlap s, compute: contribution_s / (H(s) * E[P₁]²)
    # If this ratio is bounded by poly(p), we're done
    E_P1_sq = P1_mean ** 2

    print(f"\n  E[P₁] = {P1_mean:.6e}")
    print(f"  E[P₁]² = {E_P1_sq:.6e}")
    print(f"  E[N]² = {EN**2:.4f}")

    print(f"\n  {'Overlap s':>10} {'Contrib':>12} {'H(s)':>12} "
          f"{'Contrib/H(s)':>14} {'/ E[P₁]²':>12} {'Cumulative':>12}")
    print(f"  {'-'*10} {'-'*12} {'-'*12} {'-'*14} {'-'*12} {'-'*12}")

    cum = 0
    for s in sorted(overlap_contrib.keys()):
        c = overlap_contrib[s]
        t = s - 1
        H_s = comb(k, t) * comb(p - 1 - k, k - t) * num_d12
        avg = c / H_s if H_s > 0 else 0
        norm = avg / E_P1_sq if E_P1_sq > 0 else 0
        cum += c
        print(f"  {s:>10} {c:>12.4f} {H_s:>12} {avg:>14.2e} {norm:>12.4f} {cum:>12.4f}")

    print(f"\n  Σ contributions = {cum:.4f}")
    print(f"  E[N]² = {EN**2:.4f}")
    print(f"  Ratio Σ/E[N]² = {cum/EN**2:.4f} (= E[N²]/E[N]²)")

    t_end = time.time()
    print(f"\n  Total time: {t_end-t0:.2f}s")

    return {
        'p': p, 'EN': EN, 'EN2': EN2, 'ratio': ratio,
        'overlap_contrib': dict(overlap_contrib),
        'overlap_pair_count': dict(overlap_pair_count),
    }


if __name__ == '__main__':
    results = {}
    for p in [11, 19]:
        results[p] = run_analysis(p)

    # Summary comparison
    print(f"\n{'='*70}")
    print(f"  CROSS-PRIME SUMMARY")
    print(f"{'='*70}")
    for p in [11, 19]:
        r = results[p]
        print(f"\n  p={p}: E[N]={r['EN']:.2f}, E[N²]={r['EN2']:.2f}, "
              f"ratio={r['ratio']:.2f}")
        print(f"    Overlap distribution of E[N²]:")
        for s in sorted(r['overlap_contrib'].keys()):
            frac = r['overlap_contrib'][s] / r['EN2']
            print(f"      s={s}: {r['overlap_contrib'][s]:.2f} ({frac:.1%})")
