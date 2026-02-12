#!/usr/bin/env python3
"""Fast SA search for valid D12 for near-flat D11 orbits at p=43.

Uses NumPy FFT for fast cost evaluation.
Samples a manageable subset of the 1847 near-flat-only orbits.

Usage: python -u near_flat_p43_sa.py
"""

import numpy as np
import json
import time
import sys
import os
import random
from math import comb
from itertools import combinations
from collections import defaultdict

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def symmetric_pairs(p):
    return [(d, p - d) for d in range(1, (p + 1) // 2)]


def enumerate_all_symmetric_d11(p):
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
    canonical = tuple(sorted(d11_set))
    for g in range(2, p):
        transformed = frozenset((g * d) % p for d in d11_set)
        t = tuple(sorted(transformed))
        if t < canonical:
            canonical = t
    return canonical


def group_into_orbits(all_d11, p):
    orbits = defaultdict(list)
    for d11 in all_d11:
        canon = multiplicative_orbit_canonical(d11, p)
        orbits[canon].append(d11)
    return orbits


def fft_cost(D11_set, D12_set, p):
    """Compute violation cost using FFT. Returns total excess."""
    n = (p + 1) // 2
    m = p
    N = 2 * m
    D22_set = set(range(1, p)) - D11_set
    red_thresh = n - 2
    blue_thresh = n - 1

    # Indicators
    ind_d11 = np.zeros(p, dtype=np.float64)
    for j in D11_set: ind_d11[j] = 1.0
    ind_d12 = np.zeros(p, dtype=np.float64)
    for j in D12_set: ind_d12[j] = 1.0
    ind_d22 = np.zeros(p, dtype=np.float64)
    for j in D22_set: ind_d22[j] = 1.0

    D12T = {(-x) % p for x in D12_set}
    ind_d12t = np.zeros(p, dtype=np.float64)
    for j in D12T: ind_d12t[j] = 1.0

    # FFTs
    fft_d11 = np.fft.fft(ind_d11)
    fft_d12 = np.fft.fft(ind_d12)
    fft_d22 = np.fft.fft(ind_d22)
    fft_d12t = np.fft.fft(ind_d12t)

    # Autocorrelations: Delta(A,A,d) = IFFT(|FFT(A)|^2)
    A_auto = np.round(np.fft.ifft(np.abs(fft_d11)**2).real).astype(np.int64)
    B_auto = np.round(np.fft.ifft(np.abs(fft_d12)**2).real).astype(np.int64)
    C_auto = np.round(np.fft.ifft(np.abs(fft_d22)**2).real).astype(np.int64)
    BT_auto = np.round(np.fft.ifft(np.abs(fft_d12t)**2).real).astype(np.int64)

    # Cross: Sigma(D11, D12, d) = convolution = IFFT(FFT(D11) * FFT(D12))
    sigma = np.round(np.fft.ifft(fft_d11 * fft_d12).real).astype(np.int64)

    # Cross: Delta(D12, D22, d) = IFFT(FFT(D12) * conj(FFT(D22)))
    delta_d12_d22 = np.round(np.fft.ifft(fft_d12 * np.conj(fft_d22)).real).astype(np.int64)

    d1 = len(D11_set) + len(D12_set)
    d2 = len(D22_set) + len(D12_set)

    cost = 0

    # V1V1 edges (d=1..p-1)
    for d in range(1, p):
        common = int(A_auto[d]) + int(B_auto[d])
        if d in D11_set:
            if common > red_thresh:
                cost += common - red_thresh
        else:
            blue_common = (N - 2) - 2 * d1 + common
            if blue_common > blue_thresh:
                cost += blue_common - blue_thresh

    # V2V2 edges (d=1..p-1)
    for d in range(1, p):
        common = int(C_auto[d]) + int(BT_auto[d])
        if d in D22_set:
            if common > red_thresh:
                cost += common - red_thresh
        else:
            blue_common = (N - 2) - 2 * d2 + common
            if blue_common > blue_thresh:
                cost += blue_common - blue_thresh

    # V1V2 edges (d=0..p-1)
    for d in range(p):
        common = int(sigma[d]) + int(delta_d12_d22[d])
        if d in D12_set:
            if common > red_thresh:
                cost += common - red_thresh
        else:
            blue_common = (N - 2) - d1 - d2 + common
            if blue_common > blue_thresh:
                cost += blue_common - blue_thresh

    return cost


def sa_find_d12_fft(p, D11_set, max_iter=80000, n_trials=2,
                    T_init=3.0, cooling=0.999970):
    """SA to find valid D12 for given D11, using FFT cost."""
    k = (p - 3) // 2  # |D12| - 1

    best_global_cost = float('inf')
    best_global_d12 = None

    for trial in range(n_trials):
        rest = list(range(1, p))
        random.shuffle(rest)
        D12 = {0} | set(rest[:k])

        cost = fft_cost(D11_set, D12, p)
        best_cost = cost
        best_d12 = D12.copy()
        temp = T_init

        for it in range(max_iter):
            if cost == 0:
                return D12, 0

            out_list = [x for x in D12 if x != 0]
            in_list = list(set(range(1, p)) - D12)
            d_out = random.choice(out_list)
            d_in = random.choice(in_list)

            D12.discard(d_out)
            D12.add(d_in)
            new_cost = fft_cost(D11_set, D12, p)
            delta = new_cost - cost

            if delta < 0 or random.random() < np.exp(-delta / max(temp, 0.001)):
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

    return best_global_d12, best_global_cost


def main():
    P = 43
    n = (P + 1) // 2  # 22
    k = (P - 3) // 2  # 20
    a_flat_threshold = (P + 1) // 4  # 11
    near_flat_threshold = a_flat_threshold + 1  # 12

    print("=" * 80)
    print(f"NEAR-FLAT SA ANALYSIS FOR p = {P}")
    print("=" * 80)
    print(f"  n={n}, |D11|={n}, |D12|={k+1}")
    print(f"  A-flat threshold: {a_flat_threshold}, near-flat: {near_flat_threshold}")
    sys.stdout.flush()

    # Known D11
    known_d11 = frozenset({1, 2, 5, 10, 11, 13, 16, 17, 18, 19, 20, 23, 24, 25, 26, 27,
                           30, 32, 33, 38, 41, 42})

    # Step 1: Enumerate
    print(f"\n--- Enumerating D11 ---")
    t0 = time.time()
    all_d11, pairs = enumerate_all_symmetric_d11(P)
    print(f"  {len(all_d11):,} D11 in {time.time()-t0:.1f}s")
    sys.stdout.flush()

    print(f"  Computing A-profiles...")
    t0 = time.time()
    A_dict = batch_autocorrelation_fft(all_d11, P)
    print(f"  Done in {time.time()-t0:.1f}s")
    sys.stdout.flush()

    # Classify
    maxA_to_d11 = defaultdict(list)
    for d11 in all_d11:
        A = A_dict[d11]
        max_A = max(int(A[d]) for d in d11)
        maxA_to_d11[max_A].append(d11)

    max_A_dist = {m_: len(lst) for m_, lst in sorted(maxA_to_d11.items())}
    print(f"\n  max_A distribution:")
    for val in sorted(max_A_dist.keys()):
        tags = []
        if val <= a_flat_threshold: tags.append("A-flat")
        elif val <= near_flat_threshold: tags.append("NF-only")
        tag_str = f" [{', '.join(tags)}]" if tags else ""
        print(f"    max_A={val:3d}: {max_A_dist[val]:6d}{tag_str}")

    # Collect near-flat D11
    near_flat_d11 = []
    for max_A_val in sorted(maxA_to_d11.keys()):
        if max_A_val <= near_flat_threshold:
            near_flat_d11.extend(maxA_to_d11[max_A_val])

    # Check known D11
    if known_d11 in A_dict:
        A = A_dict[known_d11]
        mA = max(int(A[d]) for d in known_d11)
        print(f"\n  Known D11 max_A = {mA}")
    sys.stdout.flush()

    # Step 2: Group near-flat into orbits
    print(f"\n--- Orbit decomposition ---")
    t0 = time.time()
    orbits_nf = group_into_orbits(near_flat_d11, P)
    print(f"  {len(orbits_nf)} near-flat orbits in {time.time()-t0:.1f}s")

    # Classify
    aflat_orbits = []
    nf_only_orbits = []
    for canon, members in orbits_nf.items():
        rep = members[0]
        A = A_dict[rep]
        max_A = max(int(A[d]) for d in rep)
        entry = {"canon": canon, "members": members, "max_A": max_A, "rep": rep}
        if max_A <= a_flat_threshold:
            aflat_orbits.append(entry)
        else:
            nf_only_orbits.append(entry)

    print(f"  A-flat orbits: {len(aflat_orbits)}")
    print(f"  NF-only orbits: {len(nf_only_orbits)}")
    sys.stdout.flush()

    # Step 3: SA search
    # For 2015 orbits, full SA is expensive. Strategy:
    # - Run all 168 A-flat orbits (quick check)
    # - Sample 200 of the 1847 NF-only orbits
    # - Always include the known D11's orbit

    random.seed(42)

    SA_ITERS = 80000
    SA_TRIALS = 2

    # Find known D11's orbit
    known_canon = multiplicative_orbit_canonical(known_d11, P)
    known_orbit_idx = None

    results_list = []

    print(f"\n--- SA on A-flat orbits ({len(aflat_orbits)} total) ---")
    print(f"  {SA_ITERS:,} iters x {SA_TRIALS} trials each")
    sys.stdout.flush()

    af_working = 0
    for oi, entry in enumerate(aflat_orbits):
        D11_set = set(entry["rep"])
        t0 = time.time()
        best_d12, best_cost = sa_find_d12_fft(
            P, D11_set, max_iter=SA_ITERS, n_trials=SA_TRIALS
        )
        elapsed = time.time() - t0
        working = (best_cost == 0)
        if working:
            af_working += 1

        results_list.append({
            "canonical": list(entry["canon"]),
            "max_A": entry["max_A"],
            "orbit_size": len(entry["members"]),
            "is_aflat": True,
            "sa_cost": best_cost,
            "working": working,
        })

        if (oi + 1) % 10 == 0 or working:
            print(f"  [{oi+1:3d}/{len(aflat_orbits)}] max_A={entry['max_A']} "
                  f"{'FOUND' if working else f'cost={best_cost}'} "
                  f"({elapsed:.1f}s) [af_working={af_working}]")
            sys.stdout.flush()

    print(f"\n  A-flat: {af_working}/{len(aflat_orbits)} working")
    sys.stdout.flush()

    # NF-only orbits: sample 200 + always include known
    print(f"\n--- SA on NF-only orbit sample ---")

    # First check: is the known D11 in a NF-only orbit?
    known_in_nf = False
    known_entry_idx = None
    for i, entry in enumerate(nf_only_orbits):
        if entry["canon"] == known_canon:
            known_in_nf = True
            known_entry_idx = i
            break

    SAMPLE_SIZE = min(100, len(nf_only_orbits))
    sample_indices = set()
    if known_entry_idx is not None:
        sample_indices.add(known_entry_idx)
    while len(sample_indices) < SAMPLE_SIZE:
        sample_indices.add(random.randint(0, len(nf_only_orbits) - 1))
    sample_indices = sorted(sample_indices)

    print(f"  Sampling {len(sample_indices)} of {len(nf_only_orbits)} NF-only orbits")
    if known_in_nf:
        print(f"  Known D11 orbit included (idx={known_entry_idx})")
    sys.stdout.flush()

    nf_working = 0
    for si, oi in enumerate(sample_indices):
        entry = nf_only_orbits[oi]
        D11_set = set(entry["rep"])
        t0 = time.time()
        best_d12, best_cost = sa_find_d12_fft(
            P, D11_set, max_iter=SA_ITERS, n_trials=SA_TRIALS
        )
        elapsed = time.time() - t0
        working = (best_cost == 0)
        if working:
            nf_working += 1

        is_known = (entry["canon"] == known_canon)
        results_list.append({
            "canonical": list(entry["canon"]),
            "max_A": entry["max_A"],
            "orbit_size": len(entry["members"]),
            "is_aflat": False,
            "sa_cost": best_cost,
            "working": working,
            "is_known_orbit": is_known,
        })

        if (si + 1) % 20 == 0 or working or is_known:
            tag = " ***KNOWN***" if is_known else ""
            print(f"  [{si+1:3d}/{len(sample_indices)}] max_A={entry['max_A']} "
                  f"{'FOUND' if working else f'cost={best_cost}'} "
                  f"({elapsed:.1f}s) [nf_working={nf_working}]{tag}")
            sys.stdout.flush()

    print(f"\n  NF-only sample: {nf_working}/{len(sample_indices)} working")

    # Extrapolate
    p_working_nf_sample = nf_working / len(sample_indices) if sample_indices else 0
    estimated_nf_working_total = int(p_working_nf_sample * len(nf_only_orbits))

    print(f"  p_working (sample) = {p_working_nf_sample:.4f}")
    print(f"  Estimated total NF-only working: ~{estimated_nf_working_total}")
    sys.stdout.flush()

    # Summary
    print(f"\n{'='*80}")
    print(f"SUMMARY FOR p={P}")
    print(f"{'='*80}")
    p_af = af_working / len(aflat_orbits) if aflat_orbits else 0
    p_nf = p_working_nf_sample
    total_working = af_working + estimated_nf_working_total
    total_orbits = len(aflat_orbits) + len(nf_only_orbits)
    p_combined = total_working / total_orbits if total_orbits > 0 else 0

    print(f"  A-flat orbits:  {len(aflat_orbits):5d}, working: {af_working:4d}, "
          f"p_working = {p_af:.4f}")
    print(f"  NF-only orbits: {len(nf_only_orbits):5d}, working: ~{estimated_nf_working_total:4d} "
          f"(sample {nf_working}/{len(sample_indices)}), p_working = {p_nf:.4f}")
    print(f"  Combined:       {total_orbits:5d}, working: ~{total_working:4d}, "
          f"p_working = {p_combined:.4f}")

    # Save
    output = {
        "p": P,
        "n": n,
        "a_flat_threshold": a_flat_threshold,
        "near_flat_threshold": near_flat_threshold,
        "total_d11": len(all_d11),
        "num_aflat_d11": sum(len(maxA_to_d11[m_]) for m_ in maxA_to_d11 if m_ <= a_flat_threshold),
        "num_nf_only_d11": len(maxA_to_d11.get(near_flat_threshold, [])),
        "num_aflat_orbits": len(aflat_orbits),
        "num_nf_only_orbits": len(nf_only_orbits),
        "aflat_working": af_working,
        "nf_sample_size": len(sample_indices),
        "nf_sample_working": nf_working,
        "p_working_aflat": p_af,
        "p_working_nf_sample": p_nf,
        "estimated_nf_working_total": estimated_nf_working_total,
        "sa_params": {"max_iter": SA_ITERS, "n_trials": SA_TRIALS},
        "orbit_results": results_list,
        "max_A_distribution": max_A_dist,
    }

    output_path = os.path.join(SCRIPT_DIR, "near_flat_p43_results.json")
    with open(output_path, "w") as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\n  Saved to {output_path}")
    sys.stdout.flush()


if __name__ == "__main__":
    main()
