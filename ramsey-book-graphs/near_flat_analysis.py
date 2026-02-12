#!/usr/bin/env python3
"""Near-flat D11 second moment analysis for p=31 and p=43.

"Near-flat" = max_A(D11) <= floor(E[A]) + 1.
This is a relaxation of "A-flat" (max_A <= floor(E[A])).

At p=43, ALL A-flat D11 have N=0 (no valid D12), but the known solution
has max_A=12 = floor(E[A])+1, so near-flat is the right class.

TASK 1 (p=31): Enumerate all symmetric D11, sample D12 for near-flat class,
compute second moment ratio for near-flat vs A-flat.

TASK 2 (p=43): Enumerate, classify, then use SA to find working near-flat orbits.
For the known D11, heavy-sample to estimate N.

Usage: python -u near_flat_analysis.py
"""

import numpy as np
import json
import time
import sys
import os
import math
import random
from math import comb, log2
from itertools import combinations
from collections import defaultdict, Counter

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


# =============================================================================
# Core utilities (from p31_orbit_analysis.py)
# =============================================================================

def symmetric_pairs(p):
    """Return list of (d, p-d) pairs for d=1..p-1, each pair listed once."""
    return [(d, p - d) for d in range(1, (p + 1) // 2)]


def enumerate_all_symmetric_d11(p):
    """Enumerate all symmetric D11 of size (p+1)/2."""
    n = (p + 1) // 2
    num_pairs = n // 2
    pairs = symmetric_pairs(p)
    all_d11 = []
    for chosen in combinations(range(len(pairs)), num_pairs):
        d11 = set()
        for idx in chosen:
            d, comp = pairs[idx]
            d11.add(d)
            d11.add(comp)
        all_d11.append(frozenset(d11))
    return all_d11, pairs


def batch_autocorrelation_fft(all_d11, p, batch_size=5000):
    """Compute autocorrelation A(d) for all D11 via batch FFT."""
    results = {}
    d11_list = list(all_d11)
    for start in range(0, len(d11_list), batch_size):
        end = min(start + batch_size, len(d11_list))
        batch = d11_list[start:end]
        bs = len(batch)
        indicator_matrix = np.zeros((bs, p), dtype=np.float64)
        for i, d11 in enumerate(batch):
            for j in d11:
                indicator_matrix[i, j] = 1.0
        fft_vals = np.fft.fft(indicator_matrix, axis=1)
        autocorr = np.fft.ifft(np.abs(fft_vals) ** 2, axis=1).real
        A_batch = np.round(autocorr).astype(np.int32)
        for i, d11 in enumerate(batch):
            results[d11] = A_batch[i]
    return results


def multiplicative_orbit_canonical(d11_set, p):
    """Compute canonical form of D11 under multiplicative action of Z_p*."""
    canonical = tuple(sorted(d11_set))
    for g in range(2, p):
        transformed = frozenset((g * d) % p for d in d11_set)
        t = tuple(sorted(transformed))
        if t < canonical:
            canonical = t
    return canonical


def group_into_orbits(all_d11, p):
    """Group D11 sets by multiplicative orbit. Returns dict: canon -> list."""
    orbits = defaultdict(list)
    for d11 in all_d11:
        canon = multiplicative_orbit_canonical(d11, p)
        orbits[canon].append(d11)
    return orbits


def batch_sample_d12(p, A, d11_arr, d22_arr, num_samples, batch_size, rng):
    """Sample random D12 and count valid ones.

    Returns (num_valid, total_sampled).
    """
    k = (p - 3) // 2
    threshold_binding = (p - 3) // 2
    threshold_loose = (p + 3) // 2
    A_at_d11 = A[d11_arr]
    A_at_d22 = A[d22_arr]
    candidates = np.arange(1, p)
    num_valid = 0
    total_sampled = 0

    while total_sampled < num_samples:
        bs = min(batch_size, num_samples - total_sampled)
        indicator_matrix = np.zeros((bs, p), dtype=np.float64)
        indicator_matrix[:, 0] = 1.0
        for i in range(bs):
            chosen = rng.choice(candidates, size=k, replace=False)
            indicator_matrix[i, chosen] = 1.0

        fft_vals = np.fft.fft(indicator_matrix, axis=1)
        B_batch = np.round(np.fft.ifft(np.abs(fft_vals) ** 2, axis=1).real).astype(np.int32)

        F_binding = A_at_d11[np.newaxis, :] + B_batch[:, d11_arr]
        all_binding_ok = np.all(F_binding <= threshold_binding, axis=1)

        if d22_arr.size > 0:
            F_loose = A_at_d22[np.newaxis, :] + B_batch[:, d22_arr]
            all_loose_ok = np.all(F_loose <= threshold_loose, axis=1)
            valid_mask = all_binding_ok & all_loose_ok
        else:
            valid_mask = all_binding_ok

        num_valid += int(valid_mask.sum())
        total_sampled += bs

    return num_valid, total_sampled


# =============================================================================
# SA solver for finding valid D12 given fixed D11
# =============================================================================

def sa_find_d12(p, D11_set, max_iter=500000, n_trials=5, T_init=3.0,
                cooling=0.999995, verbose=False):
    """Use SA to find a valid D12 for a given D11.

    Returns (D12, cost) if found (cost=0), else (best_D12, best_cost).
    """
    m = p
    n = (p + 1) // 2
    N = 2 * m
    k = (p - 3) // 2  # |D12| - 1 (excluding 0)
    d12_size = k + 1
    D22_set = set(range(1, p)) - D11_set
    red_thresh = n - 2
    blue_thresh = n - 1

    D11_sorted = sorted(D11_set)
    D22_sorted = sorted(D22_set)
    d11_arr = np.array(D11_sorted, dtype=np.int32)
    d22_arr = np.array(D22_sorted, dtype=np.int32)

    # Precompute A(d) for D11
    indicator_d11 = np.zeros(p, dtype=np.float64)
    for j in D11_set:
        indicator_d11[j] = 1.0
    fft_d11 = np.fft.fft(indicator_d11)
    A = np.round(np.fft.ifft(np.abs(fft_d11) ** 2).real).astype(np.int32)

    A_at_d11 = A[d11_arr]
    A_at_d22 = A[d22_arr]

    def compute_cost(D12):
        indicator_d12 = np.zeros(p, dtype=np.float64)
        for j in D12:
            indicator_d12[j] = 1.0
        fft_d12 = np.fft.fft(indicator_d12)
        B = np.round(np.fft.ifft(np.abs(fft_d12) ** 2).real).astype(np.int32)

        # Also need D12^T autocorrelation for V2V2 constraints
        D12T = {(-x) % p for x in D12}
        indicator_d12t = np.zeros(p, dtype=np.float64)
        for j in D12T:
            indicator_d12t[j] = 1.0
        fft_d12t = np.fft.fft(indicator_d12t)
        B_T = np.round(np.fft.ifft(np.abs(fft_d12t) ** 2).real).astype(np.int32)

        # D22 autocorrelation
        indicator_d22 = np.zeros(p, dtype=np.float64)
        for j in D22_set:
            indicator_d22[j] = 1.0
        fft_d22 = np.fft.fft(indicator_d22)
        C = np.round(np.fft.ifft(np.abs(fft_d22) ** 2).real).astype(np.int32)

        d1 = len(D11_set) + len(D12)
        d2 = len(D22_set) + len(D12)

        cost = 0
        # V1V1: A(d) + B(d) constraints
        for d in range(1, p):
            common = int(A[d]) + int(B[d])
            if d in D11_set:
                if common > red_thresh:
                    cost += common - red_thresh
            else:
                blue_common = (N - 2) - 2 * d1 + common
                if blue_common > blue_thresh:
                    cost += blue_common - blue_thresh

        # V2V2: C(d) + B_T(d) constraints
        for d in range(1, p):
            common = int(C[d]) + int(B_T[d])
            if d in D22_set:
                if common > red_thresh:
                    cost += common - red_thresh
            else:
                blue_common = (N - 2) - 2 * d2 + common
                if blue_common > blue_thresh:
                    cost += blue_common - blue_thresh

        # V1V2: Sigma(D11, D12, d) + Delta(D12, D22, d) constraints
        # Use convolution for Sigma
        indicator_d12_rev = np.zeros(p, dtype=np.float64)
        for j in D12:
            indicator_d12_rev[(-j) % p] = 1.0
        sigma_vals = np.round(np.fft.ifft(fft_d11 * np.fft.fft(indicator_d12_rev)).real).astype(np.int32)
        # No wait - Sigma counts a+b=d, which is convolution, not correlation
        # Sigma(D11, D12, d) = sum_{a in D11} 1_{(d-a) in D12}
        # = (indicator_D11 * indicator_D12)(d) where * is convolution
        sigma_vals_correct = np.round(np.fft.ifft(fft_d11 * np.fft.fft(indicator_d12)).real).astype(np.int32)

        # Delta(D12, D22, d) = correlation
        delta_d12_d22 = np.round(np.fft.ifft(np.conj(np.fft.fft(indicator_d12)) * fft_d22).real).astype(np.int32)
        # Wait, Delta(A,B,d) = #{(a,b): a-b=d} = #{a in A: a-d in B}
        # = sum_a indicator_A(a) * indicator_B(a-d)
        # This is cross-correlation of A and B at lag d
        # = IFFT(FFT(A) * conj(FFT(B)))[d]

        indicator_d22_np = np.zeros(p, dtype=np.float64)
        for j in D22_set:
            indicator_d22_np[j] = 1.0
        delta_d12_d22 = np.round(np.fft.ifft(np.fft.fft(indicator_d12) * np.conj(np.fft.fft(indicator_d22_np))).real).astype(np.int32)

        for d in range(p):
            common = int(sigma_vals_correct[d]) + int(delta_d12_d22[d])
            if d in D12:
                if common > red_thresh:
                    cost += common - red_thresh
            else:
                blue_common = (N - 2) - d1 - d2 + common
                if blue_common > blue_thresh:
                    cost += blue_common - blue_thresh

        return cost

    best_global_cost = float('inf')
    best_global_d12 = None

    for trial in range(n_trials):
        # Random init: D12 = {0} + k random from {1,...,p-1}
        rest = list(range(1, p))
        random.shuffle(rest)
        D12 = {0} | set(rest[:k])

        cost = compute_cost(D12)
        best_cost = cost
        best_d12 = D12.copy()

        temp = T_init
        for it in range(max_iter):
            if cost == 0:
                return D12, 0

            # Swap move: remove one element (not 0), add one
            out_list = [x for x in D12 if x != 0]
            in_list = list(set(range(1, p)) - D12)
            if not out_list or not in_list:
                break
            d_out = random.choice(out_list)
            d_in = random.choice(in_list)

            D12.discard(d_out)
            D12.add(d_in)
            new_cost = compute_cost(D12)
            delta = new_cost - cost

            if delta < 0 or random.random() < math.exp(-delta / max(temp, 0.001)):
                cost = new_cost
                if cost < best_cost:
                    best_cost = cost
                    best_d12 = D12.copy()
            else:
                D12.discard(d_in)
                D12.add(d_out)

            temp *= cooling

        if best_cost < best_global_cost:
            best_global_cost = best_cost
            best_global_d12 = best_d12.copy()

        if verbose:
            print(f"    SA trial {trial}: best_cost={best_cost}")
            sys.stdout.flush()

    return best_global_d12, best_global_cost


def sa_find_d12_fast(p, D11_set, max_iter=500000, n_trials=5, T_init=3.0,
                     cooling=0.999995, verbose=False):
    """Faster SA using the ramsey_core verifier for cost."""
    sys.path.insert(0, SCRIPT_DIR)
    from ramsey_core import BlockCirculantGraph, verify_construction

    m = p
    n = (p + 1) // 2
    k = (p - 3) // 2
    D22_set = set(range(1, p)) - D11_set

    best_global_cost = float('inf')
    best_global_d12 = None

    for trial in range(n_trials):
        rest = list(range(1, p))
        random.shuffle(rest)
        D12 = {0} | set(rest[:k])

        G = BlockCirculantGraph(n=n, D11=D11_set, D12=D12, D22=D22_set)
        result = verify_construction(G)
        cost = sum(ex for _, _, ex in result.violations)

        best_cost = cost
        best_d12 = D12.copy()
        temp = T_init

        for it in range(max_iter):
            if cost == 0:
                return D12, 0

            out_list = [x for x in D12 if x != 0]
            in_list = list(set(range(1, p)) - D12)
            if not out_list or not in_list:
                break
            d_out = random.choice(out_list)
            d_in = random.choice(in_list)

            D12.discard(d_out)
            D12.add(d_in)

            G = BlockCirculantGraph(n=n, D11=D11_set, D12=D12, D22=D22_set)
            result = verify_construction(G)
            new_cost = sum(ex for _, _, ex in result.violations)

            delta = new_cost - cost
            if delta < 0 or random.random() < math.exp(-delta / max(temp, 0.001)):
                cost = new_cost
                if cost < best_cost:
                    best_cost = cost
                    best_d12 = D12.copy()
            else:
                D12.discard(d_in)
                D12.add(d_out)

            temp *= cooling

        if best_cost < best_global_cost:
            best_global_cost = best_cost
            best_global_d12 = best_d12.copy()

        if verbose:
            print(f"    SA trial {trial}: best_cost={best_cost}")
            sys.stdout.flush()

    return best_global_d12, best_global_cost


# =============================================================================
# Analysis for a single prime
# =============================================================================

def analyze_near_flat(p, num_samples_per_orbit=500000, batch_size=10000,
                      sa_iters=300000, sa_trials=3, heavy_samples=0):
    """Full near-flat analysis for prime p.

    Returns dict with all results.
    """
    t_start = time.time()
    n = (p + 1) // 2
    k = (p - 3) // 2
    total_d12 = comb(p - 1, k)
    a_flat_threshold = (p + 1) // 4
    near_flat_threshold = a_flat_threshold + 1
    threshold_binding = (p - 3) // 2

    num_pairs = (p - 1) // 2
    num_choose = n // 2
    total_d11_count = comb(num_pairs, num_choose)
    orbit_size = (p - 1) // 2  # for p prime, generic orbit

    print(f"\n{'='*80}")
    print(f"NEAR-FLAT ANALYSIS FOR p = {p}")
    print(f"{'='*80}")
    print(f"  n={n}, |D11|={n}, |D12|={k+1}")
    print(f"  E[A] = (p+1)/4 = {(p+1)/4:.1f}")
    print(f"  A-flat threshold: max_A <= {a_flat_threshold}")
    print(f"  Near-flat threshold: max_A <= {near_flat_threshold}")
    print(f"  Binding threshold: A(d)+B(d) <= {threshold_binding}")
    print(f"  C({num_pairs},{num_choose}) = {total_d11_count:,} symmetric D11")
    print(f"  C({p-1},{k}) = {total_d12:,} total D12")
    print(f"  Generic orbit size: {orbit_size}")
    sys.stdout.flush()

    # Step 1: Enumerate all D11 and compute autocorrelations
    print(f"\n--- Step 1: Enumerate and compute A-profiles ---")
    t0 = time.time()
    all_d11, pairs = enumerate_all_symmetric_d11(p)
    print(f"  Enumerated {len(all_d11):,} D11 in {time.time()-t0:.1f}s")
    sys.stdout.flush()

    t0 = time.time()
    A_dict = batch_autocorrelation_fft(all_d11, p)
    print(f"  Computed A-profiles in {time.time()-t0:.1f}s")
    sys.stdout.flush()

    # Step 2: Classify by max_A
    print(f"\n--- Step 2: Classification ---")
    maxA_to_d11 = defaultdict(list)
    for d11 in all_d11:
        A = A_dict[d11]
        max_A = max(int(A[d]) for d in d11)
        maxA_to_d11[max_A].append(d11)

    max_A_dist = {m: len(lst) for m, lst in sorted(maxA_to_d11.items())}
    print(f"  max_A distribution:")
    num_aflat = 0
    num_near_flat = 0
    for val in sorted(max_A_dist.keys()):
        count = max_A_dist[val]
        tags = []
        if val <= a_flat_threshold:
            tags.append("A-flat")
            num_aflat += count
            num_near_flat += count
        elif val <= near_flat_threshold:
            tags.append("near-flat-only")
            num_near_flat += count
        elif val > threshold_binding:
            tags.append("INFEASIBLE")
        tag_str = f" [{', '.join(tags)}]" if tags else ""
        print(f"    max_A={val:3d}: {count:6d} D11{tag_str}")

    num_near_flat_only = num_near_flat - num_aflat
    print(f"\n  A-flat D11: {num_aflat}")
    print(f"  Near-flat-only D11 (max_A={near_flat_threshold}): {num_near_flat_only}")
    print(f"  Total near-flat D11: {num_near_flat}")
    sys.stdout.flush()

    # Step 3: Group near-flat D11 by orbit
    print(f"\n--- Step 3: Orbit decomposition ---")
    t0 = time.time()

    # We only need orbits of near-flat D11
    near_flat_d11 = []
    for max_A_val in sorted(maxA_to_d11.keys()):
        if max_A_val <= near_flat_threshold:
            near_flat_d11.extend(maxA_to_d11[max_A_val])

    orbits_nf = group_into_orbits(near_flat_d11, p)
    print(f"  Near-flat orbits: {len(orbits_nf)} (in {time.time()-t0:.1f}s)")

    # Classify orbits by max_A
    aflat_orbits = {}
    nf_only_orbits = {}
    for canon, members in orbits_nf.items():
        rep = members[0]
        A = A_dict[rep]
        max_A = max(int(A[d]) for d in rep)
        if max_A <= a_flat_threshold:
            aflat_orbits[canon] = {"members": members, "max_A": max_A, "rep": rep}
        else:
            nf_only_orbits[canon] = {"members": members, "max_A": max_A, "rep": rep}

    print(f"  A-flat orbits: {len(aflat_orbits)}")
    print(f"  Near-flat-only orbits (max_A={near_flat_threshold}): {len(nf_only_orbits)}")
    sys.stdout.flush()

    # Step 4: Check for known solutions
    print(f"\n--- Step 4: Known solutions ---")
    known_d11_sets = {}
    known = {
        31: {6, 7, 8, 10, 11, 12, 14, 15, 16, 17, 19, 20, 21, 23, 24, 25},
        43: {1, 2, 5, 10, 11, 13, 16, 17, 18, 19, 20, 23, 24, 25, 26, 27,
             30, 32, 33, 38, 41, 42},
    }
    if p in known:
        D11_known = known[p]
        D11_fs = frozenset(D11_known)
        if D11_fs in A_dict:
            A = A_dict[D11_fs]
            mA = max(int(A[d]) for d in D11_known)
            known_d11_sets["hardcoded"] = (D11_fs, mA)
            print(f"  Known D11: max_A={mA}, A-flat={mA <= a_flat_threshold}, "
                  f"near-flat={mA <= near_flat_threshold}")
        else:
            print(f"  Known D11 NOT in enumeration (wrong size?)")

    # Also check solution files
    for fname in [f"solution_n{n}_sa.json", f"p{p}_solutions.json"]:
        fpath = os.path.join(SCRIPT_DIR, fname)
        if os.path.exists(fpath):
            try:
                with open(fpath) as f:
                    data = json.load(f)
                if isinstance(data, dict) and "D11" in data:
                    D11_raw = set(data["D11"])
                elif isinstance(data, dict) and "solutions" in data:
                    for sol in data["solutions"][:1]:
                        D11_raw = set(sol["D11"])
                        break
                elif isinstance(data, list) and len(data) > 0:
                    D11_raw = set(data[0].get("D11", []))
                else:
                    continue
                target = (p + 1) // 2
                if len(D11_raw) != target:
                    D11_raw = set(range(1, p)) - D11_raw
                D11_fs = frozenset(D11_raw)
                if D11_fs in A_dict:
                    A = A_dict[D11_fs]
                    mA = max(int(A[d]) for d in D11_raw)
                    known_d11_sets[fname] = (D11_fs, mA)
                    print(f"  {fname}: max_A={mA}")
            except Exception as e:
                print(f"  {fname}: error loading: {e}")
    sys.stdout.flush()

    # Step 5: D12 sampling or SA for orbit representatives
    rng = np.random.default_rng(seed=2024)
    random.seed(2024)

    orbit_results = []

    use_sampling = (p <= 31)  # Sampling feasible for p<=31

    if use_sampling:
        print(f"\n--- Step 5: D12 sampling for near-flat orbit reps ---")
        print(f"  {num_samples_per_orbit:,} samples per orbit rep")
        sys.stdout.flush()

        all_nf_orbits = list(orbits_nf.items())
        all_nf_orbits.sort(key=lambda x: (max(int(A_dict[x[1][0]][d]) for d in x[1][0]), x[0]))
        # Sort by max_A of representative

        for oi, (canon, members) in enumerate(all_nf_orbits):
            rep = members[0]
            A = A_dict[rep]
            max_A = max(int(A[d]) for d in rep)
            d11_arr = np.array(sorted(rep), dtype=np.int32)
            D22 = frozenset(range(1, p)) - rep
            d22_arr = np.array(sorted(D22), dtype=np.int32)

            t0 = time.time()
            num_valid, total_sampled = batch_sample_d12(
                p, A, d11_arr, d22_arr, num_samples_per_orbit, batch_size, rng
            )
            elapsed = time.time() - t0

            rate = num_valid / total_sampled if total_sampled > 0 else 0
            N_est = rate * total_d12

            orbit_results.append({
                "canonical": list(canon),
                "orbit_size": len(members),
                "max_A": max_A,
                "is_aflat": max_A <= a_flat_threshold,
                "is_near_flat_only": a_flat_threshold < max_A <= near_flat_threshold,
                "num_valid": num_valid,
                "total_sampled": total_sampled,
                "rate": rate,
                "N_estimate": N_est,
                "D11": sorted(rep),
            })

            class_tag = "A-flat" if max_A <= a_flat_threshold else "NF-only"
            print(f"  [{oi+1:3d}/{len(all_nf_orbits)}] {class_tag:7s} max_A={max_A} "
                  f"valid={num_valid:6d}/{total_sampled:,} "
                  f"N_est={N_est:12.0f} ({elapsed:.1f}s)")
            sys.stdout.flush()

    else:
        # p >= 43: Use SA to find valid D12
        print(f"\n--- Step 5: SA search for valid D12 (p={p}) ---")
        sys.stdout.flush()

        # Heavy sampling for known D11
        for name, (D11_fs, mA) in known_d11_sets.items():
            if heavy_samples > 0:
                print(f"\n  Heavy sampling known D11 ({name}, max_A={mA}): "
                      f"{heavy_samples:,} samples...")
                sys.stdout.flush()
                A_val = A_dict[D11_fs]
                d11_arr = np.array(sorted(D11_fs), dtype=np.int32)
                D22 = frozenset(range(1, p)) - D11_fs
                d22_arr = np.array(sorted(D22), dtype=np.int32)
                t0 = time.time()
                num_valid, total_sampled = batch_sample_d12(
                    p, A_val, d11_arr, d22_arr, heavy_samples, batch_size, rng
                )
                elapsed = time.time() - t0
                rate = num_valid / total_sampled if total_sampled > 0 else 0
                N_est = rate * total_d12
                print(f"    hits={num_valid}/{total_sampled:,} rate={rate:.3e} "
                      f"N_est={N_est:.3e} ({elapsed:.1f}s)")
                sys.stdout.flush()

        # SA for orbit reps
        all_nf_orbits = list(orbits_nf.items())
        total_nf_orbits = len(all_nf_orbits)
        print(f"\n  Running SA on {total_nf_orbits} near-flat orbit reps...")
        print(f"  SA params: {sa_iters:,} iters, {sa_trials} trials per orbit")
        sys.stdout.flush()

        num_working = 0
        for oi, (canon, members) in enumerate(all_nf_orbits):
            rep = members[0]
            A = A_dict[rep]
            max_A = max(int(A[d]) for d in rep)
            D11_set = set(rep)

            t0 = time.time()
            best_d12, best_cost = sa_find_d12_fast(
                p, D11_set, max_iter=sa_iters, n_trials=sa_trials,
                T_init=3.0, cooling=0.999995, verbose=False
            )
            elapsed = time.time() - t0

            working = (best_cost == 0)
            if working:
                num_working += 1

            orbit_results.append({
                "canonical": list(canon),
                "orbit_size": len(members),
                "max_A": max_A,
                "is_aflat": max_A <= a_flat_threshold,
                "is_near_flat_only": a_flat_threshold < max_A <= near_flat_threshold,
                "sa_cost": best_cost,
                "working": working,
                "D11": sorted(rep),
                "D12": sorted(best_d12) if working else None,
            })

            class_tag = "A-flat" if max_A <= a_flat_threshold else "NF-only"
            status = "FOUND" if working else f"cost={best_cost}"
            if (oi + 1) % 20 == 0 or working or oi == 0:
                print(f"  [{oi+1:3d}/{total_nf_orbits}] {class_tag:7s} max_A={max_A} "
                      f"{status} ({elapsed:.1f}s) [working so far: {num_working}]")
                sys.stdout.flush()

    # Step 6: Compute second moment ratios
    print(f"\n--- Step 6: Second moment ratio computation ---")
    sys.stdout.flush()

    if use_sampling:
        # Sampling-based ratio: E[N^2] / E[N]^2
        # Over all near-flat orbits (each orbit rep is equally likely within class)

        # A-flat class
        af_N = [r["N_estimate"] for r in orbit_results if r["is_aflat"]]
        af_count = sum(r["orbit_size"] for r in orbit_results if r["is_aflat"])

        # Near-flat-only class
        nf_N = [r["N_estimate"] for r in orbit_results if r["is_near_flat_only"]]
        nf_count = sum(r["orbit_size"] for r in orbit_results if r["is_near_flat_only"])

        # Combined near-flat
        all_nf_N = af_N + nf_N

        def compute_ratio(N_values, weights=None):
            """Compute E[N^2]/E[N]^2 with optional weights."""
            if not N_values:
                return float('inf'), 0, 0
            arr = np.array(N_values)
            if weights is not None:
                w = np.array(weights, dtype=np.float64)
                w = w / w.sum()
                E_N = float(np.sum(w * arr))
                E_N2 = float(np.sum(w * arr**2))
            else:
                E_N = float(np.mean(arr))
                E_N2 = float(np.mean(arr**2))
            ratio = E_N2 / (E_N**2) if E_N > 0 else float('inf')
            return ratio, E_N, E_N2

        # Weight by orbit size (each D11 equally likely)
        af_weights = [r["orbit_size"] for r in orbit_results if r["is_aflat"]]
        nf_weights = [r["orbit_size"] for r in orbit_results if r["is_near_flat_only"]]
        all_nf_weights = af_weights + nf_weights

        ratio_af, EN_af, EN2_af = compute_ratio(af_N, af_weights)
        ratio_nf_only, EN_nf, EN2_nf = compute_ratio(nf_N, nf_weights)
        ratio_all_nf, EN_all, EN2_all = compute_ratio(all_nf_N, all_nf_weights)

        # Working counts
        af_working = sum(1 for r in orbit_results if r["is_aflat"] and r["N_estimate"] > 0)
        nf_working = sum(1 for r in orbit_results if r["is_near_flat_only"] and r["N_estimate"] > 0)

        print(f"\n  A-flat class:")
        print(f"    D11 count: {af_count}")
        print(f"    Orbits: {len(af_N)} ({af_working} working)")
        print(f"    E[N] = {EN_af:.4e}")
        print(f"    E[N^2] = {EN2_af:.4e}")
        print(f"    Ratio E[N^2]/E[N]^2 = {ratio_af:.4f}")

        print(f"\n  Near-flat-only class (max_A={near_flat_threshold}):")
        print(f"    D11 count: {nf_count}")
        print(f"    Orbits: {len(nf_N)} ({nf_working} working)")
        print(f"    E[N] = {EN_nf:.4e}")
        print(f"    E[N^2] = {EN2_nf:.4e}")
        print(f"    Ratio E[N^2]/E[N]^2 = {ratio_nf_only:.4f}")

        print(f"\n  Combined near-flat (max_A <= {near_flat_threshold}):")
        print(f"    D11 count: {af_count + nf_count}")
        print(f"    Orbits: {len(all_nf_N)} ({af_working + nf_working} working)")
        print(f"    E[N] = {EN_all:.4e}")
        print(f"    E[N^2] = {EN2_all:.4e}")
        print(f"    Ratio E[N^2]/E[N]^2 = {ratio_all_nf:.4f}")

    else:
        # SA-based: compute p_working and estimate ratio
        af_orbits_total = sum(1 for r in orbit_results if r["is_aflat"])
        af_working = sum(1 for r in orbit_results if r["is_aflat"] and r.get("working", False))
        nf_orbits_total = sum(1 for r in orbit_results if r["is_near_flat_only"])
        nf_working = sum(1 for r in orbit_results if r["is_near_flat_only"] and r.get("working", False))

        af_count = sum(r["orbit_size"] for r in orbit_results if r["is_aflat"])
        nf_count = sum(r["orbit_size"] for r in orbit_results if r["is_near_flat_only"])

        p_working_af = af_working / af_orbits_total if af_orbits_total > 0 else 0
        p_working_nf = nf_working / nf_orbits_total if nf_orbits_total > 0 else 0
        p_working_all = (af_working + nf_working) / (af_orbits_total + nf_orbits_total) if (af_orbits_total + nf_orbits_total) > 0 else 0

        ratio_af = float('inf')
        ratio_nf_only = float('inf')
        ratio_all_nf = float('inf')

        print(f"\n  A-flat class:")
        print(f"    D11 count: {af_count}")
        print(f"    Orbits: {af_orbits_total} ({af_working} working)")
        print(f"    p_working = {p_working_af:.4f}")

        print(f"\n  Near-flat-only class (max_A={near_flat_threshold}):")
        print(f"    D11 count: {nf_count}")
        print(f"    Orbits: {nf_orbits_total} ({nf_working} working)")
        print(f"    p_working = {p_working_nf:.4f}")

        print(f"\n  Combined near-flat (max_A <= {near_flat_threshold}):")
        print(f"    D11 count: {af_count + nf_count}")
        print(f"    Orbits: {af_orbits_total + nf_orbits_total} ({af_working + nf_working} working)")
        print(f"    p_working = {p_working_all:.4f}")

    # Step 7: Summary table
    print(f"\n{'='*80}")
    print(f"SUMMARY for p={p}")
    print(f"{'='*80}")

    summary = {
        "p": p,
        "n": n,
        "total_d11": total_d11_count,
        "total_d12": total_d12,
        "a_flat_threshold": a_flat_threshold,
        "near_flat_threshold": near_flat_threshold,
        "max_A_distribution": max_A_dist,
        "num_aflat_d11": num_aflat,
        "num_near_flat_d11": num_near_flat,
        "num_near_flat_only_d11": num_near_flat_only,
    }

    if use_sampling:
        summary.update({
            "aflat_orbits": len(af_N),
            "aflat_working": af_working,
            "aflat_ratio": ratio_af,
            "aflat_EN": EN_af,
            "aflat_EN2": EN2_af,
            "nf_only_orbits": len(nf_N),
            "nf_only_working": nf_working,
            "nf_only_ratio": ratio_nf_only,
            "nf_only_EN": EN_nf,
            "nf_only_EN2": EN2_nf,
            "combined_nf_orbits": len(all_nf_N),
            "combined_nf_working": af_working + nf_working,
            "combined_nf_ratio": ratio_all_nf,
            "combined_nf_EN": EN_all,
            "combined_nf_EN2": EN2_all,
        })
    else:
        summary.update({
            "aflat_orbits": af_orbits_total,
            "aflat_working": af_working,
            "aflat_p_working": p_working_af,
            "nf_only_orbits": nf_orbits_total,
            "nf_only_working": nf_working,
            "nf_only_p_working": p_working_nf,
            "combined_nf_orbits": af_orbits_total + nf_orbits_total,
            "combined_nf_working": af_working + nf_working,
            "combined_nf_p_working": p_working_all,
        })

    elapsed_total = time.time() - t_start
    summary["elapsed_seconds"] = elapsed_total
    summary["orbit_results"] = orbit_results

    print(f"\n  Total elapsed: {elapsed_total:.1f}s")
    sys.stdout.flush()

    return summary


def main():
    print("=" * 80)
    print("NEAR-FLAT D11 SECOND MOMENT ANALYSIS")
    print("R(B_{n-1}, B_n) = 4n-1 for p = 3 mod 4 primes")
    print("=" * 80)
    print(f"Date: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    sys.stdout.flush()

    results = {}

    # =========================================================================
    # p = 31: Full sampling analysis
    # =========================================================================
    r31 = analyze_near_flat(
        p=31,
        num_samples_per_orbit=500000,
        batch_size=10000,
    )
    results["31"] = r31

    # =========================================================================
    # p = 43: SA-based analysis
    # =========================================================================
    r43 = analyze_near_flat(
        p=43,
        num_samples_per_orbit=0,  # not used for SA
        batch_size=10000,
        sa_iters=300000,
        sa_trials=3,
        heavy_samples=50_000_000,  # 50M for known D11
    )
    results["43"] = r43

    # =========================================================================
    # Combined summary table
    # =========================================================================
    print(f"\n{'='*80}")
    print("COMBINED RESULTS TABLE")
    print(f"{'='*80}")
    print(f"\n{'p':>4s} | {'Class':>10s} | {'#D11':>8s} | {'#Orbits':>8s} | "
          f"{'#Working':>8s} | {'p_work':>8s} | {'Ratio':>10s}")
    print("-" * 75)

    for p_str in ["31", "43"]:
        r = results[p_str]
        p_val = r["p"]

        # A-flat row
        if "aflat_ratio" in r:
            ratio_str = f"{r['aflat_ratio']:.2f}"
            pwork = f"{r['aflat_working']}/{r['aflat_orbits']}"
        else:
            ratio_str = "N/A (SA)"
            pwork = f"{r['aflat_working']}/{r['aflat_orbits']}"

        print(f"{p_val:4d} | {'A-flat':>10s} | {r['num_aflat_d11']:8d} | "
              f"{r.get('aflat_orbits',0):8d} | {r.get('aflat_working',0):8d} | "
              f"{pwork:>8s} | {ratio_str:>10s}")

        # Near-flat-only row
        if "nf_only_ratio" in r:
            ratio_str = f"{r['nf_only_ratio']:.2f}"
            pwork = f"{r['nf_only_working']}/{r['nf_only_orbits']}"
        else:
            ratio_str = "N/A (SA)"
            pwork = f"{r['nf_only_working']}/{r['nf_only_orbits']}"

        print(f"{p_val:4d} | {'NF-only':>10s} | {r['num_near_flat_only_d11']:8d} | "
              f"{r.get('nf_only_orbits',0):8d} | {r.get('nf_only_working',0):8d} | "
              f"{pwork:>8s} | {ratio_str:>10s}")

        # Combined row
        if "combined_nf_ratio" in r:
            ratio_str = f"{r['combined_nf_ratio']:.2f}"
            pwork = f"{r['combined_nf_working']}/{r['combined_nf_orbits']}"
        else:
            ratio_str = "N/A (SA)"
            pwork = f"{r['combined_nf_working']}/{r['combined_nf_orbits']}"

        print(f"{p_val:4d} | {'Combined':>10s} | {r['num_near_flat_d11']:8d} | "
              f"{r.get('combined_nf_orbits',0):8d} | {r.get('combined_nf_working',0):8d} | "
              f"{pwork:>8s} | {ratio_str:>10s}")
        print("-" * 75)

    # Save results
    output_path = os.path.join(SCRIPT_DIR, "near_flat_results.json")

    # Convert for JSON serialization
    def make_serializable(obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, frozenset):
            return sorted(obj)
        if isinstance(obj, set):
            return sorted(obj)
        return obj

    class NumpyEncoder(json.JSONEncoder):
        def default(self, obj):
            result = make_serializable(obj)
            if result is not obj:
                return result
            return super().default(obj)

    with open(output_path, "w") as f:
        json.dump(results, f, indent=2, cls=NumpyEncoder, default=str)

    print(f"\nResults saved to {output_path}")
    sys.stdout.flush()


if __name__ == "__main__":
    main()
