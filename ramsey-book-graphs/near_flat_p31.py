#!/usr/bin/env python3
"""Near-flat D11 second moment analysis for p=31.

Computes E[N^2]/E[N]^2 for A-flat and near-flat D11 classes.

Usage: python -u near_flat_p31.py
"""

import numpy as np
import json
import time
import sys
import os
from math import comb
from itertools import combinations
from collections import defaultdict

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
P = 31
N_VAL = (P + 1) // 2       # 16
K_D12 = (P - 3) // 2       # 14
A_FLAT = (P + 1) // 4      # 8
NEAR_FLAT = A_FLAT + 1     # 9
THRESH_BIND = (P - 3) // 2 # 14
THRESH_LOOSE = (P + 3) // 2 # 17
TOTAL_D12 = comb(P - 1, K_D12)
NUM_SAMPLES = 500_000
BATCH = 10000


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


def sample_d12(p, A, d11_arr, d22_arr, num_samples, batch_size, rng):
    k = (p - 3) // 2
    tb = (p - 3) // 2
    tl = (p + 3) // 2
    A_d11 = A[d11_arr]
    A_d22 = A[d22_arr]
    candidates = np.arange(1, p)
    nv = 0
    ns = 0
    while ns < num_samples:
        bs = min(batch_size, num_samples - ns)
        ind = np.zeros((bs, p), dtype=np.float64)
        ind[:, 0] = 1.0
        for i in range(bs):
            chosen = rng.choice(candidates, size=k, replace=False)
            ind[i, chosen] = 1.0
        fv = np.fft.fft(ind, axis=1)
        B = np.round(np.fft.ifft(np.abs(fv)**2, axis=1).real).astype(np.int32)
        ok_bind = np.all(A_d11[np.newaxis, :] + B[:, d11_arr] <= tb, axis=1)
        if d22_arr.size > 0:
            ok_loose = np.all(A_d22[np.newaxis, :] + B[:, d22_arr] <= tl, axis=1)
            ok = ok_bind & ok_loose
        else:
            ok = ok_bind
        nv += int(ok.sum())
        ns += bs
    return nv, ns


def main():
    t_start = time.time()
    print("=" * 80)
    print(f"NEAR-FLAT D11 ANALYSIS FOR p = {P}")
    print("=" * 80)
    print(f"  n={N_VAL}, |D11|={N_VAL}, |D12|={K_D12+1}")
    print(f"  A-flat: max_A <= {A_FLAT}, near-flat: max_A <= {NEAR_FLAT}")
    print(f"  C({(P-1)//2},{N_VAL//2}) = {comb((P-1)//2, N_VAL//2):,} D11")
    print(f"  C({P-1},{K_D12}) = {TOTAL_D12:,} D12")
    print(f"  Samples per orbit: {NUM_SAMPLES:,}")
    sys.stdout.flush()

    # Enumerate
    all_d11, pairs = enumerate_all_symmetric_d11(P)
    A_dict = batch_autocorrelation_fft(all_d11, P)

    # Classify
    maxA_to_d11 = defaultdict(list)
    for d11 in all_d11:
        A = A_dict[d11]
        max_A = max(int(A[d]) for d in d11)
        maxA_to_d11[max_A].append(d11)

    print(f"\n  max_A distribution:")
    for val in sorted(maxA_to_d11.keys()):
        count = len(maxA_to_d11[val])
        tags = []
        if val <= A_FLAT: tags.append("A-flat")
        elif val <= NEAR_FLAT: tags.append("NF-only")
        tag_str = f" [{', '.join(tags)}]" if tags else ""
        print(f"    max_A={val:3d}: {count:6d}{tag_str}")

    # Near-flat D11
    nf_d11 = []
    for mA in sorted(maxA_to_d11.keys()):
        if mA <= NEAR_FLAT:
            nf_d11.extend(maxA_to_d11[mA])

    # Orbits
    orbits = group_into_orbits(nf_d11, P)
    print(f"\n  Near-flat orbits: {len(orbits)}")

    af_orbits = []
    nf_only_orbits = []
    for canon, members in orbits.items():
        rep = members[0]
        A = A_dict[rep]
        max_A = max(int(A[d]) for d in rep)
        entry = {"canon": canon, "members": members, "max_A": max_A, "rep": rep}
        if max_A <= A_FLAT:
            af_orbits.append(entry)
        else:
            nf_only_orbits.append(entry)

    print(f"  A-flat orbits: {len(af_orbits)}, NF-only: {len(nf_only_orbits)}")
    sys.stdout.flush()

    # Sample D12 for all orbits
    rng = np.random.default_rng(seed=2024)
    all_entries = sorted(af_orbits + nf_only_orbits, key=lambda e: (e["max_A"], e["canon"]))

    print(f"\n--- Sampling ---")
    for i, entry in enumerate(all_entries):
        rep = entry["rep"]
        A = A_dict[rep]
        d11_arr = np.array(sorted(rep), dtype=np.int32)
        D22 = frozenset(range(1, P)) - rep
        d22_arr = np.array(sorted(D22), dtype=np.int32)
        t0 = time.time()
        nv, ns = sample_d12(P, A, d11_arr, d22_arr, NUM_SAMPLES, BATCH, rng)
        elapsed = time.time() - t0
        rate = nv / ns
        N_est = rate * TOTAL_D12
        entry["num_valid"] = nv
        entry["total_sampled"] = ns
        entry["rate"] = rate
        entry["N_estimate"] = N_est

        tag = "AF" if entry["max_A"] <= A_FLAT else "NF"
        if (i + 1) % 20 == 0 or nv > 0:
            print(f"  [{i+1:3d}/{len(all_entries)}] {tag} max_A={entry['max_A']} "
                  f"valid={nv:6d} N_est={N_est:12.0f} ({elapsed:.1f}s)")
            sys.stdout.flush()

    # Compute ratios
    print(f"\n--- Second moment ratios ---")

    def compute_ratio(entries):
        if not entries:
            return float('inf'), 0, 0, 0, 0
        weights = np.array([e["orbit_size"] for e in entries], dtype=np.float64)
        weights /= weights.sum()
        N_arr = np.array([e["N_estimate"] for e in entries])
        E_N = float(np.sum(weights * N_arr))
        E_N2 = float(np.sum(weights * N_arr**2))
        ratio = E_N2 / (E_N**2) if E_N > 0 else float('inf')
        n_working = sum(1 for e in entries if e["N_estimate"] > 0)
        return ratio, E_N, E_N2, n_working, len(entries)

    # Add orbit_size
    for e in all_entries:
        e["orbit_size"] = len(e["members"])

    af_entries = [e for e in all_entries if e["max_A"] <= A_FLAT]
    nf_entries = [e for e in all_entries if e["max_A"] > A_FLAT]

    r_af, en_af, en2_af, w_af, n_af = compute_ratio(af_entries)
    r_nf, en_nf, en2_nf, w_nf, n_nf = compute_ratio(nf_entries)
    r_all, en_all, en2_all, w_all, n_all = compute_ratio(all_entries)

    af_d11_count = sum(e["orbit_size"] for e in af_entries)
    nf_d11_count = sum(e["orbit_size"] for e in nf_entries)
    all_d11_count = af_d11_count + nf_d11_count

    print(f"\n  A-flat (max_A <= {A_FLAT}):")
    print(f"    #D11 = {af_d11_count}, #orbits = {n_af}, working = {w_af}")
    print(f"    E[N] = {en_af:.4e}, E[N^2] = {en2_af:.4e}")
    print(f"    Ratio = {r_af:.4f}")

    print(f"\n  NF-only (max_A = {NEAR_FLAT}):")
    print(f"    #D11 = {nf_d11_count}, #orbits = {n_nf}, working = {w_nf}")
    print(f"    E[N] = {en_nf:.4e}, E[N^2] = {en2_nf:.4e}")
    print(f"    Ratio = {r_nf:.4f}")

    print(f"\n  Combined near-flat (max_A <= {NEAR_FLAT}):")
    print(f"    #D11 = {all_d11_count}, #orbits = {n_all}, working = {w_all}")
    print(f"    E[N] = {en_all:.4e}, E[N^2] = {en2_all:.4e}")
    print(f"    Ratio = {r_all:.4f}")

    elapsed_total = time.time() - t_start
    print(f"\n  Total time: {elapsed_total:.1f}s")

    # Save
    output = {
        "p": P,
        "n": N_VAL,
        "a_flat_threshold": A_FLAT,
        "near_flat_threshold": NEAR_FLAT,
        "total_d11": len(all_d11),
        "total_d12": TOTAL_D12,
        "num_aflat_d11": af_d11_count,
        "num_nf_only_d11": nf_d11_count,
        "num_aflat_orbits": n_af,
        "num_nf_only_orbits": n_nf,
        "aflat_working": w_af,
        "nf_only_working": w_nf,
        "aflat_ratio": r_af,
        "nf_only_ratio": r_nf,
        "combined_ratio": r_all,
        "aflat_EN": en_af,
        "aflat_EN2": en2_af,
        "nf_only_EN": en_nf,
        "nf_only_EN2": en2_nf,
        "combined_EN": en_all,
        "combined_EN2": en2_all,
        "max_A_distribution": {str(k): len(v) for k, v in sorted(maxA_to_d11.items())},
        "orbit_details": [
            {
                "canonical": list(e["canon"]),
                "max_A": e["max_A"],
                "orbit_size": e["orbit_size"],
                "is_aflat": e["max_A"] <= A_FLAT,
                "num_valid": e["num_valid"],
                "total_sampled": e["total_sampled"],
                "rate": e["rate"],
                "N_estimate": e["N_estimate"],
            }
            for e in all_entries
        ],
        "elapsed_seconds": elapsed_total,
    }

    output_path = os.path.join(SCRIPT_DIR, "near_flat_p31_results.json")
    with open(output_path, "w") as f:
        json.dump(output, f, indent=2, default=str)
    print(f"  Saved to {output_path}")
    sys.stdout.flush()


if __name__ == "__main__":
    main()
