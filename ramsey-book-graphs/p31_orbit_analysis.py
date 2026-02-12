#!/usr/bin/env python3
"""Orbit analysis of ALL symmetric D11 at p=31 for Ramsey book graph construction.

p=31: n=16, |D11|=16 = 8 symmetric pairs from 15 total pairs.
C(15,8) = 6435 total symmetric D11.
|D12| = 14, E[A] = 8.0, binding threshold = 14, loose threshold = 17.

Exhaustive D12 enumeration impossible (C(30,14) ~ 145M).
Instead, ESTIMATE N(D11) via random sampling of D12.

Steps:
1. Enumerate all 6435 symmetric D11
2. Compute A(d) for each via FFT
3. Identify A-flat D11 (max A on D11 <= 8)
4. Group by multiplicative orbit (D11 ~ g*D11 mod p)
5. For each A-flat orbit rep, sample 500K random D12, check constraints
6. Report structural quantities and N estimates

Usage: python -u p31_orbit_analysis.py
"""

import numpy as np
import json
import time
import sys
import os
from math import comb, log2
from itertools import combinations
from collections import defaultdict


# =============================================================================
# Parameters for p = 31
# =============================================================================
P = 31
N = (P + 1) // 2          # 16
K_D12 = (P - 3) // 2      # 14  (elements of D12 from {1,...,30})
THRESHOLD_BINDING = (P - 3) // 2  # 14
THRESHOLD_LOOSE = (P + 3) // 2    # 17
A_FLAT_THRESHOLD = (P + 1) // 4   # 8
NUM_SAMPLES = 500_000
BATCH_SIZE = 10_000


def symmetric_pairs(p):
    """Return list of (d, p-d) pairs for d=1..p-1, each pair listed once."""
    pairs = []
    seen = set()
    for d in range(1, p):
        if d not in seen:
            comp = p - d
            pairs.append((min(d, comp), max(d, comp)))
            seen.add(d)
            seen.add(comp)
    return pairs


def enumerate_all_symmetric_d11(p):
    """Enumerate all symmetric D11 of size (p+1)/2.

    Each D11 is a union of (p+1)/4 symmetric pairs from the (p-1)/2 pairs.
    Returns list of frozensets.
    """
    n = (p + 1) // 2
    num_pairs = n // 2  # 8 pairs for p=31
    pairs = symmetric_pairs(p)
    assert len(pairs) == (p - 1) // 2, f"Expected {(p-1)//2} pairs, got {len(pairs)}"

    all_d11 = []
    for chosen in combinations(range(len(pairs)), num_pairs):
        d11 = set()
        for idx in chosen:
            d, comp = pairs[idx]
            d11.add(d)
            d11.add(comp)
        all_d11.append(frozenset(d11))

    return all_d11, pairs


def autocorrelation_fft(S, p):
    """Compute A(d) = #{(a,b) in S x S : a-b = d mod p, a != b} for all d=0..p-1."""
    indicator = np.zeros(p, dtype=np.float64)
    for j in S:
        indicator[j] = 1.0
    fft_val = np.fft.fft(indicator)
    autocorr = np.fft.ifft(np.abs(fft_val) ** 2).real
    return np.round(autocorr).astype(np.int32)


def primitive_root(p):
    """Find a primitive root mod p."""
    for g in range(2, p):
        seen = set()
        val = 1
        for _ in range(p - 1):
            seen.add(val)
            val = (val * g) % p
        if len(seen) == p - 1:
            return g
    raise ValueError(f"No primitive root found for p={p}")


def multiplicative_orbit(d11_set, p):
    """Compute the multiplicative orbit of a D11 set.

    Two D11 are in the same orbit if one = {g*d mod p : d in D11} for some g
    with gcd(g, p) = 1.

    Returns canonical form (smallest lexicographic sorted tuple).
    """
    # For p prime, the multiplicative group Z_p^* has order p-1.
    # We check all g in {1,...,p-1}.
    canonical = tuple(sorted(d11_set))

    for g in range(2, p):
        transformed = frozenset((g * d) % p for d in d11_set)
        t = tuple(sorted(transformed))
        if t < canonical:
            canonical = t

    return canonical


def group_by_orbit(all_d11, p):
    """Group D11 sets by multiplicative orbit.

    Returns dict: canonical_form -> list of D11 sets in that orbit.
    """
    orbits = defaultdict(list)
    d11_to_canonical = {}

    for d11 in all_d11:
        canon = multiplicative_orbit(d11, p)
        orbits[canon].append(d11)
        d11_to_canonical[d11] = canon

    return orbits, d11_to_canonical


def batch_sample_d12(p, k, A, d11_arr, d22_arr, A_at_d11, A_at_d22,
                     threshold_binding, threshold_loose,
                     num_samples, batch_size, rng):
    """Sample random D12 in batches and check both constraints.

    D12 = {0} union (k elements from {1,...,p-1}).
    Returns (num_valid_full, num_valid_binding_only, total_sampled).
    """
    candidates = np.arange(1, p)
    num_valid_full = 0
    num_valid_binding = 0
    total_sampled = 0

    while total_sampled < num_samples:
        bs = min(batch_size, num_samples - total_sampled)

        # Generate batch of random D12 indicator vectors
        indicator_matrix = np.zeros((bs, p), dtype=np.float64)
        indicator_matrix[:, 0] = 1.0  # 0 always in D12

        for i in range(bs):
            chosen = rng.choice(candidates, size=k, replace=False)
            indicator_matrix[i, chosen] = 1.0

        # Batch FFT to compute B(d) for all D12 simultaneously
        fft_vals = np.fft.fft(indicator_matrix, axis=1)
        B_batch = np.round(np.fft.ifft(np.abs(fft_vals) ** 2, axis=1).real).astype(np.int32)

        # Check binding constraints: A(d) + B(d) <= threshold_binding for d in D11
        F_binding = A_at_d11[np.newaxis, :] + B_batch[:, d11_arr]  # (bs, |D11|)
        all_binding_ok = np.all(F_binding <= threshold_binding, axis=1)
        num_valid_binding += int(all_binding_ok.sum())

        # Check loose constraints: A(d) + B(d) <= threshold_loose for d in D22
        if d22_arr.size > 0:
            F_loose = A_at_d22[np.newaxis, :] + B_batch[:, d22_arr]  # (bs, |D22|)
            all_loose_ok = np.all(F_loose <= threshold_loose, axis=1)
            valid_full = all_binding_ok & all_loose_ok
        else:
            valid_full = all_binding_ok

        num_valid_full += int(valid_full.sum())
        total_sampled += bs

    return num_valid_full, num_valid_binding, total_sampled


def main():
    t_start = time.time()

    print("=" * 80)
    print(f"ORBIT ANALYSIS OF SYMMETRIC D11 AT p = {P}")
    print("=" * 80)
    print(f"Parameters:")
    print(f"  p = {P}, n = {N}")
    print(f"  |D11| = {N} (= {N//2} symmetric pairs from {(P-1)//2} total)")
    print(f"  |D12| = {K_D12 + 1} (= {{0}} + {K_D12} from {{1,...,{P-1}}})")
    print(f"  C({(P-1)//2},{N//2}) = {comb((P-1)//2, N//2)} symmetric D11")
    print(f"  C({P-1},{K_D12}) = {comb(P-1, K_D12):,} total D12 (too many to enumerate)")
    print(f"  Binding threshold (d in D11): A(d)+B(d) <= {THRESHOLD_BINDING}")
    print(f"  Loose threshold   (d in D22): A(d)+B(d) <= {THRESHOLD_LOOSE}")
    print(f"  A-flat threshold: max A(d) on D11 <= {A_FLAT_THRESHOLD}")
    print(f"  D12 samples per orbit rep: {NUM_SAMPLES:,}")
    print(f"  Batch size: {BATCH_SIZE:,}")
    sys.stdout.flush()

    # =========================================================================
    # Step 1: Enumerate all symmetric D11
    # =========================================================================
    print(f"\n--- Step 1: Enumerating all symmetric D11 ---")
    sys.stdout.flush()

    all_d11, pairs = enumerate_all_symmetric_d11(P)
    print(f"  Total symmetric D11: {len(all_d11)}")
    assert len(all_d11) == comb((P - 1) // 2, N // 2)
    sys.stdout.flush()

    # =========================================================================
    # Step 2: Compute A(d) for each D11 and classify
    # =========================================================================
    print(f"\n--- Step 2: Computing autocorrelation profiles ---")
    sys.stdout.flush()

    d11_data = {}  # d11 (frozenset) -> dict of properties

    for idx, d11 in enumerate(all_d11):
        A = autocorrelation_fft(d11, P)
        D22 = frozenset(range(1, P)) - d11

        d11_sorted = sorted(d11)
        d22_sorted = sorted(D22)

        A_on_D11 = sorted([int(A[d]) for d in d11_sorted])
        A_on_D22 = sorted([int(A[d]) for d in d22_sorted])

        max_A_D11 = max(A_on_D11)
        max_A_D22 = max(A_on_D22)

        S_D11 = sum(int(A[d]) for d in d11_sorted)
        S_D22 = sum(int(A[d]) for d in d22_sorted)

        d11_data[d11] = {
            "D11": d11_sorted,
            "D22": d22_sorted,
            "A_profile_D11": A_on_D11,
            "A_profile_D22": A_on_D22,
            "max_A_D11": max_A_D11,
            "max_A_D22": max_A_D22,
            "S_D11": S_D11,
            "S_D22": S_D22,
            "min_slack_D11": THRESHOLD_BINDING - max_A_D11,
            "min_slack_D22": THRESHOLD_LOOSE - max_A_D22,
            "is_aflat": max_A_D11 <= A_FLAT_THRESHOLD,
            "A_full": A,  # keep for sampling
        }

    num_aflat = sum(1 for v in d11_data.values() if v["is_aflat"])
    print(f"  A-flat D11 (max A on D11 <= {A_FLAT_THRESHOLD}): {num_aflat} / {len(all_d11)}")

    # Distribution of max_A_D11
    from collections import Counter
    max_A_dist = Counter(v["max_A_D11"] for v in d11_data.values())
    print(f"  Distribution of max A(d) on D11:")
    for val in sorted(max_A_dist.keys()):
        marker = " <-- A-flat" if val <= A_FLAT_THRESHOLD else ""
        print(f"    max_A = {val:2d}: {max_A_dist[val]:5d} D11{marker}")
    sys.stdout.flush()

    # =========================================================================
    # Step 3: Group by multiplicative orbit
    # =========================================================================
    print(f"\n--- Step 3: Grouping by multiplicative orbit ---")
    sys.stdout.flush()

    orbits, d11_to_canonical = group_by_orbit(all_d11, P)

    print(f"  Total orbits: {len(orbits)}")

    # Classify orbits
    orbit_info = {}
    for canon, members in orbits.items():
        rep = members[0]  # pick first as representative
        data = d11_data[rep]
        is_aflat = data["is_aflat"]
        orbit_info[canon] = {
            "representative": rep,
            "size": len(members),
            "is_aflat": is_aflat,
            "max_A_D11": data["max_A_D11"],
            "max_A_D22": data["max_A_D22"],
        }

    aflat_orbits = {k: v for k, v in orbit_info.items() if v["is_aflat"]}
    non_aflat_orbits = {k: v for k, v in orbit_info.items() if not v["is_aflat"]}

    print(f"  A-flat orbits: {len(aflat_orbits)}")
    print(f"  Non-A-flat orbits: {len(non_aflat_orbits)}")

    # Orbit sizes
    orbit_sizes = Counter(v["size"] for v in orbit_info.values())
    print(f"  Orbit size distribution:")
    for sz in sorted(orbit_sizes.keys()):
        print(f"    size {sz:3d}: {orbit_sizes[sz]:4d} orbits")

    # Verify: sum of orbit sizes = total D11
    total_from_orbits = sum(v["size"] for v in orbit_info.values())
    assert total_from_orbits == len(all_d11), \
        f"Orbit sizes sum to {total_from_orbits}, expected {len(all_d11)}"
    print(f"  Verification: sum of orbit sizes = {total_from_orbits} (correct)")
    sys.stdout.flush()

    # =========================================================================
    # Step 4: Primitive root and orbit structure
    # =========================================================================
    g = primitive_root(P)
    print(f"\n  Primitive root mod {P}: g = {g}")
    print(f"  Multiplicative group order: {P - 1}")

    # The orbit of D11 under multiplication by g generates all multiplications
    # since g generates Z_p^*. So orbit size divides p-1.
    print(f"  All orbit sizes divide {P - 1} = {P - 1}")
    for sz in sorted(orbit_sizes.keys()):
        divides = (P - 1) % sz == 0
        print(f"    {sz} divides {P - 1}: {divides}")
    sys.stdout.flush()

    # =========================================================================
    # Step 5: Sample D12 for each A-flat orbit representative
    # =========================================================================
    print(f"\n--- Step 5: Sampling D12 for {len(aflat_orbits)} A-flat orbit reps ---")
    print(f"  ({NUM_SAMPLES:,} samples each, batch size {BATCH_SIZE:,})")
    sys.stdout.flush()

    rng = np.random.default_rng(seed=2024)

    # Sort A-flat orbits by max_A_D11 then by canonical form
    aflat_sorted = sorted(aflat_orbits.items(), key=lambda x: (x[1]["max_A_D11"], x[0]))

    sampling_results = {}

    for orbit_idx, (canon, oinfo) in enumerate(aflat_sorted):
        rep = oinfo["representative"]
        data = d11_data[rep]

        d11_arr = np.array(data["D11"], dtype=np.int32)
        d22_arr = np.array(data["D22"], dtype=np.int32)
        A = data["A_full"]
        A_at_d11 = A[d11_arr]
        A_at_d22 = A[d22_arr]

        t0 = time.time()
        num_valid_full, num_valid_binding, total_sampled = batch_sample_d12(
            P, K_D12, A, d11_arr, d22_arr, A_at_d11, A_at_d22,
            THRESHOLD_BINDING, THRESHOLD_LOOSE,
            NUM_SAMPLES, BATCH_SIZE, rng
        )
        elapsed = time.time() - t0

        # Estimate N(D11)
        total_d12 = comb(P - 1, K_D12)
        rate_full = num_valid_full / total_sampled
        rate_binding = num_valid_binding / total_sampled
        N_estimate_full = rate_full * total_d12
        N_estimate_binding = rate_binding * total_d12

        sampling_results[canon] = {
            "num_valid_full": num_valid_full,
            "num_valid_binding": num_valid_binding,
            "total_sampled": total_sampled,
            "rate_full": rate_full,
            "rate_binding": rate_binding,
            "N_estimate_full": N_estimate_full,
            "N_estimate_binding": N_estimate_binding,
            "elapsed": elapsed,
        }

        # Print progress
        d11_str = str(data["D11"])
        if len(d11_str) > 50:
            d11_str = d11_str[:47] + "..."
        print(f"  [{orbit_idx+1:3d}/{len(aflat_sorted)}] orbit_size={oinfo['size']:2d} "
              f"max_A_D11={data['max_A_D11']} max_A_D22={data['max_A_D22']} "
              f"valid_full={num_valid_full:6d} valid_bind={num_valid_binding:6d} "
              f"N_est={N_estimate_full:12.0f} ({elapsed:.1f}s)")
        sys.stdout.flush()

    # =========================================================================
    # Step 6: Also compute stats for non-A-flat orbits (no sampling)
    # =========================================================================
    print(f"\n--- Step 6: Non-A-flat orbit statistics ---")
    non_aflat_sorted = sorted(non_aflat_orbits.items(),
                               key=lambda x: (x[1]["max_A_D11"], x[0]))

    # Just show distribution
    non_aflat_maxA = Counter(v["max_A_D11"] for v in non_aflat_orbits.values())
    print(f"  Non-A-flat orbits by max_A_D11:")
    for val in sorted(non_aflat_maxA.keys()):
        print(f"    max_A_D11 = {val}: {non_aflat_maxA[val]} orbits")
    sys.stdout.flush()

    # =========================================================================
    # Step 7: Build results table
    # =========================================================================
    print(f"\n{'='*120}")
    print(f"RESULTS TABLE: A-flat orbits sorted by estimated N(D11)")
    print(f"{'='*120}")

    # Merge data
    table_rows = []
    for canon, oinfo in aflat_sorted:
        rep = oinfo["representative"]
        data = d11_data[rep]
        samp = sampling_results[canon]

        row = {
            "canonical": list(canon),
            "orbit_size": oinfo["size"],
            "D11": data["D11"],
            "D22": data["D22"],
            "A_profile_D11": data["A_profile_D11"],
            "A_profile_D22": data["A_profile_D22"],
            "max_A_D11": data["max_A_D11"],
            "max_A_D22": data["max_A_D22"],
            "S_D11": data["S_D11"],
            "S_D22": data["S_D22"],
            "min_slack_D11": data["min_slack_D11"],
            "min_slack_D22": data["min_slack_D22"],
            "num_valid_full": samp["num_valid_full"],
            "num_valid_binding": samp["num_valid_binding"],
            "total_sampled": samp["total_sampled"],
            "rate_full": samp["rate_full"],
            "rate_binding": samp["rate_binding"],
            "N_estimate_full": samp["N_estimate_full"],
            "N_estimate_binding": samp["N_estimate_binding"],
        }
        table_rows.append(row)

    # Sort by N_estimate_full descending
    table_rows.sort(key=lambda r: -r["N_estimate_full"])

    # Print table
    print(f"{'#':>3s} | {'orb_sz':>6s} | {'maxA_D11':>8s} | {'maxA_D22':>8s} | "
          f"{'S_D11':>5s} | {'S_D22':>5s} | {'slack11':>7s} | {'slack22':>7s} | "
          f"{'hits_full':>9s} | {'hits_bind':>9s} | {'N_est_full':>12s} | {'N_est_bind':>12s} | "
          f"{'A_profile_D11'}")
    print("-" * 160)

    for i, row in enumerate(table_rows):
        aprof = str(row["A_profile_D11"])
        if len(aprof) > 50:
            aprof = aprof[:47] + "..."

        print(f"{i+1:3d} | {row['orbit_size']:6d} | {row['max_A_D11']:8d} | {row['max_A_D22']:8d} | "
              f"{row['S_D11']:5d} | {row['S_D22']:5d} | {row['min_slack_D11']:7d} | {row['min_slack_D22']:7d} | "
              f"{row['num_valid_full']:9d} | {row['num_valid_binding']:9d} | "
              f"{row['N_estimate_full']:12.0f} | {row['N_estimate_binding']:12.0f} | "
              f"{aprof}")

    # =========================================================================
    # Step 8: Pattern analysis
    # =========================================================================
    print(f"\n{'='*80}")
    print("PATTERN ANALYSIS")
    print(f"{'='*80}")

    # Q1: Does D22 A-profile predict N?
    print(f"\n--- Q1: D22 A-profile vs N(D11) ---")
    working = [r for r in table_rows if r["num_valid_full"] > 0]
    non_working = [r for r in table_rows if r["num_valid_full"] == 0]

    print(f"  Working A-flat orbits (N > 0): {len(working)} / {len(table_rows)}")
    print(f"  Non-working A-flat orbits:     {len(non_working)} / {len(table_rows)}")

    if working:
        print(f"\n  Working orbits max_A_D22 values: {sorted(set(r['max_A_D22'] for r in working))}")
    if non_working:
        print(f"  Non-working orbits max_A_D22 values: {sorted(set(r['max_A_D22'] for r in non_working))}")

    # Correlation: max_A_D22 vs working
    print(f"\n  max_A_D22 vs working:")
    all_maxA_D22 = sorted(set(r["max_A_D22"] for r in table_rows))
    for val in all_maxA_D22:
        w = sum(1 for r in working if r["max_A_D22"] == val)
        nw = sum(1 for r in non_working if r["max_A_D22"] == val)
        tot = w + nw
        pct = 100 * w / tot if tot > 0 else 0
        print(f"    max_A_D22 = {val}: {w:4d} working, {nw:4d} non-working ({pct:.0f}% working)")

    # Q2: S_D11 constant?
    print(f"\n--- Q2: Is S_D11 = sum A(d) for d in D11 constant across all D11? ---")
    s_d11_values = set(v["S_D11"] for v in d11_data.values())
    print(f"  Distinct S_D11 values: {sorted(s_d11_values)}")
    if len(s_d11_values) == 1:
        print(f"  YES! S_D11 is constant = {s_d11_values.pop()}")
    else:
        s_d11_counter = Counter(v["S_D11"] for v in d11_data.values())
        print(f"  NO, S_D11 takes {len(s_d11_values)} values:")
        for val in sorted(s_d11_counter.keys()):
            print(f"    S_D11 = {val}: {s_d11_counter[val]} D11")

    # Q3: Relationship between min_slack_D22 and N
    print(f"\n--- Q3: min_slack_D22 vs N ---")
    for r in working[:10]:
        print(f"  N_est={r['N_estimate_full']:10.0f} min_slack_D22={r['min_slack_D22']} "
              f"max_A_D22={r['max_A_D22']} A_prof_D22={r['A_profile_D22']}")

    # Q4: Top vs bottom by S_D22
    if working:
        print(f"\n--- Q4: S_D22 for top N vs bottom N ---")
        top5 = working[:5] if len(working) >= 5 else working
        bottom5 = working[-5:] if len(working) >= 5 else working
        print(f"  Top {len(top5)} by N: S_D22 = {[r['S_D22'] for r in top5]}")
        print(f"  Bottom {len(bottom5)} by N: S_D22 = {[r['S_D22'] for r in bottom5]}")

    # Q5: Best orbit details
    if working:
        print(f"\n--- Q5: Best orbit (highest estimated N) ---")
        best = working[0]
        print(f"  D11 = {best['D11']}")
        print(f"  D22 = {best['D22']}")
        print(f"  A profile on D11: {best['A_profile_D11']}")
        print(f"  A profile on D22: {best['A_profile_D22']}")
        print(f"  max_A_D11 = {best['max_A_D11']}, max_A_D22 = {best['max_A_D22']}")
        print(f"  S_D11 = {best['S_D11']}, S_D22 = {best['S_D22']}")
        print(f"  min_slack_D11 = {best['min_slack_D11']}, min_slack_D22 = {best['min_slack_D22']}")
        print(f"  Estimated N(D11) = {best['N_estimate_full']:.0f}")
        print(f"  Sample hits: {best['num_valid_full']} / {best['total_sampled']}")
        print(f"  Orbit size: {best['orbit_size']}")

    # =========================================================================
    # Summary statistics
    # =========================================================================
    print(f"\n{'='*80}")
    print("SUMMARY")
    print(f"{'='*80}")
    total_d12 = comb(P - 1, K_D12)
    print(f"  p = {P}, n = {N}")
    print(f"  Total symmetric D11: {len(all_d11)}")
    print(f"  Total orbits: {len(orbits)}")
    print(f"  A-flat D11: {num_aflat}")
    print(f"  A-flat orbits: {len(aflat_orbits)}")
    print(f"  Working A-flat orbits (N_est > 0): {len(working)}")
    print(f"  Total D12 space: C({P-1},{K_D12}) = {total_d12:,}")
    print(f"  log2(C({P-1},{K_D12})) = {log2(total_d12):.2f}")

    if working:
        best_N = max(r["N_estimate_full"] for r in working)
        total_N = sum(r["N_estimate_full"] * r["orbit_size"] for r in working)
        print(f"  Best estimated N(D11): {best_N:,.0f}")
        print(f"  Total estimated solutions (sum over working orbits, weighted): {total_N:,.0f}")

    elapsed_total = time.time() - t_start
    print(f"\n  Total elapsed time: {elapsed_total:.1f}s")
    sys.stdout.flush()

    # =========================================================================
    # Save JSON
    # =========================================================================
    output_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               "p31_orbit_results.json")

    # Build non-A-flat orbit list
    non_aflat_list = []
    for canon, oinfo in non_aflat_sorted:
        rep = oinfo["representative"]
        data = d11_data[rep]
        non_aflat_list.append({
            "canonical": list(canon),
            "orbit_size": oinfo["size"],
            "D11": data["D11"],
            "A_profile_D11": data["A_profile_D11"],
            "A_profile_D22": data["A_profile_D22"],
            "max_A_D11": data["max_A_D11"],
            "max_A_D22": data["max_A_D22"],
            "S_D11": data["S_D11"],
            "S_D22": data["S_D22"],
            "min_slack_D11": data["min_slack_D11"],
            "min_slack_D22": data["min_slack_D22"],
        })

    save_data = {
        "p": P,
        "n": N,
        "parameters": {
            "threshold_binding": THRESHOLD_BINDING,
            "threshold_loose": THRESHOLD_LOOSE,
            "a_flat_threshold": A_FLAT_THRESHOLD,
            "d12_size": K_D12 + 1,
            "total_d12": total_d12,
            "num_samples": NUM_SAMPLES,
        },
        "summary": {
            "total_symmetric_d11": len(all_d11),
            "total_orbits": len(orbits),
            "num_aflat_d11": num_aflat,
            "num_aflat_orbits": len(aflat_orbits),
            "num_working_aflat_orbits": len(working),
            "max_A_D11_distribution": {str(k): v for k, v in sorted(max_A_dist.items())},
            "s_d11_values": sorted(s_d11_values) if len(s_d11_values) > 1 else list(s_d11_values),
            "orbit_sizes": {str(k): v for k, v in sorted(orbit_sizes.items())},
        },
        "aflat_orbits_by_N": table_rows,
        "non_aflat_orbits": non_aflat_list,
        "elapsed_seconds": elapsed_total,
    }

    with open(output_path, "w") as f:
        json.dump(save_data, f, indent=2, default=str)

    print(f"\n  Results saved to {output_path}")
    sys.stdout.flush()


if __name__ == "__main__":
    main()
