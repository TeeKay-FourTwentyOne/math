#!/usr/bin/env python3
"""
Exact computation of c_0(D11) = Pr[all ok] / prod_d Pr[ok_d]
for the Ramsey book graph construction, testing Schur-convexity conjecture.

For each prime p in {11, 19, 23} (p ≡ 3 mod 4):
  - Enumerate ALL symmetric D11
  - For each D11, enumerate ALL D12 candidates
  - Compute joint validity and per-position marginals
  - Compute c_0 = joint / product_of_marginals
  - Report relationship between c_0 and Var(T)
"""

import numpy as np
import json
import time
from itertools import combinations
from math import comb, sqrt, floor
from collections import defaultdict
from fractions import Fraction


def get_symmetric_pairs(p):
    """Get all negation pairs {d, p-d} from {1,...,p-1}."""
    pairs = []
    seen = set()
    for d in range(1, p):
        if d not in seen:
            pairs.append((d, (-d) % p))
            seen.add(d)
            seen.add((-d) % p)
    return pairs


def enumerate_symmetric_D11(p, size):
    """Enumerate all symmetric subsets of {1,...,p-1} with given even size."""
    pairs = get_symmetric_pairs(p)
    num_pairs = size // 2
    for chosen_pairs in combinations(range(len(pairs)), num_pairs):
        D11 = set()
        for i in chosen_pairs:
            d, neg_d = pairs[i]
            D11.add(d)
            D11.add(neg_d)
        yield frozenset(D11)


def autocorrelation_single(D_set, p):
    """Compute autocorrelation A(d) for a single set via FFT."""
    indicator = np.zeros(p, dtype=np.float64)
    for j in D_set:
        indicator[j] = 1.0
    fft_val = np.fft.fft(indicator)
    autocorr = np.fft.ifft(np.abs(fft_val) ** 2).real
    return np.round(autocorr).astype(int)


def batch_autocorrelation(matrix, p):
    """Compute autocorrelation for many sets at once via FFT."""
    fft_vals = np.fft.fft(matrix, axis=1)
    autocorr = np.fft.ifft(np.abs(fft_vals) ** 2, axis=1).real
    return np.round(autocorr).astype(np.int32)


def build_d12_candidates(p, k):
    """Build indicator matrix for all D12 = {0} ∪ S, S a k-subset of {1,...,p-1}."""
    num_d12 = comb(p - 1, k)
    matrix = np.zeros((num_d12, p), dtype=np.float64)
    matrix[:, 0] = 1.0
    for i, chosen in enumerate(combinations(range(1, p), k)):
        for j in chosen:
            matrix[i, j] = 1.0
    return matrix


def primitive_root(p):
    """Find a primitive root mod p."""
    for g in range(2, p):
        seen = set()
        val = 1
        for _ in range(p - 1):
            val = (val * g) % p
            seen.add(val)
        if len(seen) == p - 1:
            return g
    return None


def canonical_orbit_rep(D11_set, p):
    """Find lexicographically smallest orbit member under multiplication by Z_p*."""
    g = primitive_root(p)
    best = tuple(sorted(D11_set))
    multiplier = 1
    for _ in range((p - 1) // 2):
        multiplier = (multiplier * g) % p
        transformed = frozenset((d * multiplier) % p for d in D11_set)
        candidate = tuple(sorted(transformed))
        if candidate < best:
            best = candidate
    return best


def run_prime(p):
    """Run the full c_0 computation for prime p."""
    print(f"\n{'=' * 80}")
    print(f"SCHUR-CONVEXITY EXACT: p = {p}")
    print(f"{'=' * 80}")

    n = (p + 1) // 2
    d11_size = n                   # |D11| = (p+1)/2
    d12_size = n - 1               # |D12| = (p-1)/2
    k = d12_size - 1               # |S| = (p-3)/2 (elements added to {0})
    thresh_binding = (p - 3) // 2  # A(d)+B(d) <= this for d in D11
    thresh_loose = (p + 3) // 2    # A(d)+B(d) <= this for d in D22
    E_A = (p + 1) / 4
    floor_E_A = floor(E_A)
    mu_B = (p - 3) / 4
    sigma = sqrt((p - 3) * (p + 1) / (16 * (p - 2)))

    pairs = get_symmetric_pairs(p)
    num_pairs = len(pairs)
    num_d11 = comb(num_pairs, d11_size // 2)
    num_d12 = comb(p - 1, k)

    print(f"  n = {n}, |D11| = {d11_size}, |D12| = {d12_size}, k = {k}")
    print(f"  E[A] = {E_A}, floor(E[A]) = {floor_E_A}")
    print(f"  mu_B = {mu_B}, sigma = {sigma:.6f}")
    print(f"  Binding threshold: A(d)+B(d) <= {thresh_binding}")
    print(f"  Loose threshold:   A(d)+B(d) <= {thresh_loose}")
    print(f"  #symmetric D11 = {num_d11}")
    print(f"  #D12 candidates = {num_d12}")

    t0 = time.time()

    # Precompute all D12 autocorrelations
    print(f"\n  Precomputing {num_d12} D12 autocorrelations...")
    d12_matrix = build_d12_candidates(p, k)
    B_matrix = batch_autocorrelation(d12_matrix, p)
    print(f"  Done in {time.time() - t0:.1f}s")

    # Verify that marginal CDF F(T) is the same for all positions d
    # (i.e., Pr[B(d) <= T] depends only on T, not on d)
    print(f"\n  Verifying marginal uniformity...")
    B_at_pos = {}
    for d in range(1, p):
        B_at_pos[d] = B_matrix[:, d].copy()

    # Check: for positions with same threshold, the CDFs should match
    # Pick two different positions and compare their B-value distributions
    test_positions = list(range(1, min(p, 6)))
    for d in test_positions:
        vals, counts = np.unique(B_at_pos[d], return_counts=True)
        dist = dict(zip(vals.tolist(), counts.tolist()))
    # Compare distributions across all positions
    ref_d = 1
    ref_sorted = np.sort(B_at_pos[ref_d])
    uniform = True
    for d in range(2, p):
        d_sorted = np.sort(B_at_pos[d])
        if not np.array_equal(ref_sorted, d_sorted):
            uniform = False
            print(f"    WARNING: B distribution at d={d} differs from d={ref_d}")
            break
    if uniform:
        print(f"    VERIFIED: B(d) has identical distribution for all d in {{1,...,{p-1}}}")
    else:
        print(f"    NOTE: B(d) distributions differ across positions")

    # Precompute the universal marginal CDF: for each threshold T, Pr[B(d) <= T]
    # Since all positions have the same B distribution, we only need one
    B_vals_ref = B_at_pos[1]
    max_B = int(B_vals_ref.max())
    marginal_cdf = {}  # T -> count of D12 with B(d) <= T
    for T in range(0, max_B + 5):
        marginal_cdf[T] = int(np.sum(B_vals_ref <= T))

    # Enumerate all symmetric D11
    print(f"\n  Enumerating symmetric D11...")
    d11_list = list(enumerate_symmetric_D11(p, d11_size))
    assert len(d11_list) == num_d11, f"Expected {num_d11}, got {len(d11_list)}"

    results = []

    print(f"\n  Computing c_0 for all {num_d11} symmetric D11...")
    t1 = time.time()

    for idx, D11 in enumerate(d11_list):
        D11_sorted = sorted(D11)
        D22 = sorted(set(range(1, p)) - D11)

        # Autocorrelation of D11
        A = autocorrelation_single(D11, p)

        # Threshold vector T(d) for each constrained position
        # For d in D11: T(d) = thresh_binding - A(d)
        # For d in D22: T(d) = thresh_loose - A(d)
        T_D11 = [thresh_binding - int(A[d]) for d in D11_sorted]
        T_D22 = [thresh_loose - int(A[d]) for d in D22]

        A_at_D11 = [int(A[d]) for d in D11_sorted]
        A_at_D22 = [int(A[d]) for d in D22]
        max_A_D11 = max(A_at_D11)
        is_a_flat = (max_A_D11 <= floor_E_A)

        # Count joint validity: D12 satisfying ALL constraints
        d11_indices = np.array(D11_sorted, dtype=np.int32)
        d22_indices = np.array(D22, dtype=np.int32)
        A_d11_arr = np.array(A_at_D11, dtype=np.int32)
        A_d22_arr = np.array(A_at_D22, dtype=np.int32)

        # Check binding constraints: A(d) + B(d) <= thresh_binding for d in D11
        F_binding = A_d11_arr[np.newaxis, :] + B_matrix[:, d11_indices]
        valid_binding = np.all(F_binding <= thresh_binding, axis=1)

        # Check loose constraints: A(d) + B(d) <= thresh_loose for d in D22
        F_loose = A_d22_arr[np.newaxis, :] + B_matrix[:, d22_indices]
        valid_loose = np.all(F_loose <= thresh_loose, axis=1)

        valid_mask = valid_binding & valid_loose
        N = int(valid_mask.sum())

        # Per-position marginals
        # For each d in D11: count how many D12 satisfy B(d) <= T_D11[i]
        marginals_D11 = []
        for i, d in enumerate(D11_sorted):
            T = T_D11[i]
            count = int(np.sum(B_matrix[:, d] <= T))
            marginals_D11.append(count)
            # Verify against universal CDF
            if T in marginal_cdf:
                assert count == marginal_cdf[T], \
                    f"Marginal mismatch at d={d}, T={T}: got {count}, expected {marginal_cdf[T]}"

        # For each d in D22: count how many D12 satisfy B(d) <= T_D22[i]
        marginals_D22 = []
        for i, d in enumerate(D22):
            T = T_D22[i]
            count = int(np.sum(B_matrix[:, d] <= T))
            marginals_D22.append(count)
            if T in marginal_cdf:
                assert count == marginal_cdf[T], \
                    f"Marginal mismatch at d={d}, T={T}: got {count}, expected {marginal_cdf[T]}"

        # Compute product of marginal probabilities (as fractions for exactness)
        # Pr[ok_d] = marginal_d / num_d12
        # product = prod(marginal_d) / num_d12^(num_positions)
        # c_0 = (N / num_d12) / (prod(marginals) / num_d12^num_positions)
        #      = N * num_d12^(num_positions - 1) / prod(marginals)

        all_marginals = marginals_D11 + marginals_D22
        all_T = T_D11 + T_D22
        num_positions = len(all_marginals)

        # Use logarithms for the product (exact fractions would overflow)
        log_prod_marginals = sum(np.log(m) for m in all_marginals if m > 0)
        log_joint = np.log(N) if N > 0 else float('-inf')
        log_num_d12 = np.log(num_d12)

        # c_0 = (N/num_d12) / prod(marginal_d/num_d12)
        #      = N / num_d12 * num_d12^num_positions / prod(marginals)
        #      = N * num_d12^(num_positions-1) / prod(marginals)
        if N > 0 and all(m > 0 for m in all_marginals):
            log_c0 = log_joint + (num_positions - 1) * log_num_d12 - log_prod_marginals
            c0 = np.exp(log_c0)
        elif N == 0:
            c0 = 0.0
        else:
            c0 = float('inf')  # joint > 0 but some marginal = 0 (shouldn't happen)

        # Threshold statistics
        T_all = all_T
        mean_T = np.mean(T_all)
        var_T = np.var(T_all)
        range_T = max(T_all) - min(T_all)

        # Orbit representative
        orbit_rep = canonical_orbit_rep(D11, p)

        # z-scores at D11 positions
        z_D11 = [(T - mu_B) / sigma for T in T_D11]
        z_D22 = [(T - mu_B) / sigma for T in T_D22]

        entry = {
            "D11": D11_sorted,
            "D22": D22,
            "A_at_D11": A_at_D11,
            "A_at_D22": A_at_D22,
            "T_D11": T_D11,
            "T_D22": T_D22,
            "max_A_D11": max_A_D11,
            "is_a_flat": is_a_flat,
            "N": N,
            "num_d12": num_d12,
            "marginals_D11": marginals_D11,
            "marginals_D22": marginals_D22,
            "c0": float(c0),
            "log_c0": float(log_c0) if N > 0 and all(m > 0 for m in all_marginals) else None,
            "Var_T": float(var_T),
            "range_T": int(range_T),
            "mean_T": float(mean_T),
            "T_sorted": sorted(T_all),
            "orbit_rep": list(orbit_rep),
            "z_D11": [round(z, 6) for z in z_D11],
            "z_D22": [round(z, 6) for z in z_D22],
        }
        results.append(entry)

        if (idx + 1) % 50 == 0 or idx == len(d11_list) - 1:
            elapsed = time.time() - t1
            print(f"    [{idx+1}/{num_d11}] elapsed={elapsed:.1f}s, last N={N}, c0={c0:.6f}")

    elapsed = time.time() - t0
    print(f"\n  Total time for p={p}: {elapsed:.1f}s")

    # Group by orbit
    orbit_groups = defaultdict(list)
    for e in results:
        orbit_groups[tuple(e["orbit_rep"])].append(e)

    # Verify: all members of an orbit have same N and same c0
    print(f"\n  Orbit consistency check:")
    orbit_consistent = True
    for orbit, members in orbit_groups.items():
        N_vals = set(m["N"] for m in members)
        c0_vals = set(round(m["c0"], 8) for m in members)
        if len(N_vals) > 1 or len(c0_vals) > 1:
            print(f"    INCONSISTENT orbit {orbit[:4]}...: N={N_vals}, c0={c0_vals}")
            orbit_consistent = False
    if orbit_consistent:
        print(f"    All {len(orbit_groups)} orbits consistent (same N, same c0 within orbit)")

    # Summary table
    print(f"\n{'=' * 120}")
    print(f"SUMMARY TABLE: p = {p}")
    print(f"{'=' * 120}")
    header = f"{'orbit':>30s} {'|orb|':>5s} {'maxA':>5s} {'A-flat':>6s} {'N':>8s} " \
             f"{'c0':>12s} {'Var(T)':>10s} {'range(T)':>8s} {'Pr[all]':>14s} {'prod Pr':>14s}"
    print(header)
    print("-" * 120)

    orbit_summaries = []
    for orbit in sorted(orbit_groups.keys()):
        members = orbit_groups[orbit]
        e = members[0]  # representative
        orbit_label = str(list(orbit)[:5]) + ("..." if len(orbit) > 5 else "")
        pr_all = e["N"] / num_d12 if num_d12 > 0 else 0
        pr_prod = pr_all / e["c0"] if e["c0"] > 0 else 0

        print(f"{orbit_label:>30s} {len(members):>5d} {e['max_A_D11']:>5d} "
              f"{'Y' if e['is_a_flat'] else 'N':>6s} {e['N']:>8d} "
              f"{e['c0']:>12.6f} {e['Var_T']:>10.4f} {e['range_T']:>8d} "
              f"{pr_all:>14.8f} {pr_prod:>14.8f}")

        orbit_summaries.append({
            "orbit_rep": list(orbit),
            "orbit_size": len(members),
            "max_A_D11": e["max_A_D11"],
            "is_a_flat": e["is_a_flat"],
            "N": e["N"],
            "c0": e["c0"],
            "Var_T": e["Var_T"],
            "range_T": e["range_T"],
            "T_sorted": e["T_sorted"],
            "Pr_all": pr_all,
            "Pr_prod": pr_prod,
            "A_profile_D11": sorted(e["A_at_D11"]),
            "A_profile_D22": sorted(e["A_at_D22"]),
        })

    # Schur-convexity tests
    print(f"\n{'=' * 80}")
    print(f"SCHUR-CONVEXITY TESTS: p = {p}")
    print(f"{'=' * 80}")

    # Test 1: c_0 vs Var(T) monotonicity
    c0_varT_pairs = [(e["Var_T"], e["c0"], e["N"]) for e in orbit_summaries if e["N"] > 0]
    c0_varT_pairs.sort()
    print(f"\n  Test 1: c_0 vs Var(T) (among N > 0 orbits)")
    monotone = True
    for i in range(len(c0_varT_pairs) - 1):
        v1, c1, _ = c0_varT_pairs[i]
        v2, c2, _ = c0_varT_pairs[i + 1]
        if c2 < c1 - 1e-10:
            monotone = False
            print(f"    VIOLATION: Var(T)={v1:.4f} -> c0={c1:.6f}, Var(T)={v2:.4f} -> c0={c2:.6f}")
    if monotone:
        print(f"    MONOTONE: c_0 is non-decreasing in Var(T)")
    for v, c, n in c0_varT_pairs:
        print(f"      Var(T)={v:>10.4f}  c0={c:>12.6f}  N={n}")

    # Test 2: c_0 < 1 for A-flat (uniform T)?
    print(f"\n  Test 2: c_0 for A-flat D11 (most uniform T)")
    a_flat_orbits = [e for e in orbit_summaries if e["is_a_flat"]]
    for e in a_flat_orbits:
        status = "< 1" if e["c0"] < 1 else "> 1" if e["c0"] > 1 else "= 1"
        print(f"    orbit={str(e['orbit_rep'][:5])+'...':<30s} c0={e['c0']:.6f} ({status}), N={e['N']}")

    # Test 3: c_0 > 1 for spread T?
    print(f"\n  Test 3: c_0 for non-A-flat D11 (spread T)")
    non_flat_orbits = [e for e in orbit_summaries if not e["is_a_flat"] and e["N"] > 0]
    for e in non_flat_orbits:
        status = "> 1" if e["c0"] > 1 else "< 1" if e["c0"] < 1 else "= 1"
        print(f"    orbit={str(e['orbit_rep'][:5])+'...':<30s} c0={e['c0']:.6f} ({status}), "
              f"Var(T)={e['Var_T']:.4f}, N={e['N']}")

    # Test 4: Majorization check
    print(f"\n  Test 4: Majorization pairs")
    # For pairs with same sum of T but one majorizes the other
    for i, e1 in enumerate(orbit_summaries):
        for j, e2 in enumerate(orbit_summaries):
            if i >= j:
                continue
            if e1["N"] == 0 or e2["N"] == 0:
                continue
            T1 = e1["T_sorted"]
            T2 = e2["T_sorted"]
            if len(T1) != len(T2):
                continue
            sum1 = sum(T1)
            sum2 = sum(T2)
            if sum1 != sum2:
                continue
            # Check if T2 majorizes T1 (T2 more spread)
            # T2 majorizes T1 if partial sums of sorted-desc T2 >= T1
            T1_desc = sorted(T1, reverse=True)
            T2_desc = sorted(T2, reverse=True)
            partial_sums_1 = [sum(T1_desc[:k+1]) for k in range(len(T1_desc))]
            partial_sums_2 = [sum(T2_desc[:k+1]) for k in range(len(T2_desc))]
            if all(s2 >= s1 for s1, s2 in zip(partial_sums_1, partial_sums_2)):
                # T2 majorizes T1
                if T1_desc != T2_desc:
                    check = "OK" if e2["c0"] >= e1["c0"] - 1e-10 else "FAIL"
                    print(f"    T1={T1} -> c0={e1['c0']:.6f}")
                    print(f"    T2={T2} -> c0={e2['c0']:.6f}  (T2 majorizes T1) [{check}]")
            elif all(s1 >= s2 for s1, s2 in zip(partial_sums_1, partial_sums_2)):
                # T1 majorizes T2
                if T1_desc != T2_desc:
                    check = "OK" if e1["c0"] >= e2["c0"] - 1e-10 else "FAIL"
                    print(f"    T2={T2} -> c0={e2['c0']:.6f}")
                    print(f"    T1={T1} -> c0={e1['c0']:.6f}  (T1 majorizes T2) [{check}]")

    return {
        "p": p,
        "n": n,
        "num_d11": num_d11,
        "num_d12": num_d12,
        "num_orbits": len(orbit_groups),
        "orbit_summaries": orbit_summaries,
        "marginal_cdf": {str(k): v for k, v in sorted(marginal_cdf.items())},
    }


def main():
    print("=" * 80)
    print("SCHUR-CONVEXITY EXACT TEST")
    print("c_0(D11) = Pr[all ok] / prod_d Pr[ok_d]")
    print("=" * 80)

    all_results = {}

    for p in [11, 19, 23]:
        results = run_prime(p)
        all_results[str(p)] = results

    # Cross-prime summary
    print(f"\n{'=' * 80}")
    print(f"CROSS-PRIME SUMMARY")
    print(f"{'=' * 80}")
    for p_str, res in all_results.items():
        p = res["p"]
        orbits = res["orbit_summaries"]
        n_pos = len([o for o in orbits if o["N"] > 0])
        c0_vals = [o["c0"] for o in orbits if o["N"] > 0]
        var_vals = [o["Var_T"] for o in orbits if o["N"] > 0]
        print(f"\n  p = {p}: {res['num_orbits']} orbits, {n_pos} with N>0")
        if c0_vals:
            print(f"    c0 range: [{min(c0_vals):.6f}, {max(c0_vals):.6f}]")
            print(f"    Var(T) range: [{min(var_vals):.4f}, {max(var_vals):.4f}]")
            # Correlation
            if len(c0_vals) > 1:
                corr = np.corrcoef(var_vals, c0_vals)[0, 1]
                print(f"    Correlation(Var(T), c0) = {corr:.6f}")

    # Save results
    outpath = "/Users/stephenpadgett/Projects/math/ramsey-book-graphs/schur_convexity_exact.json"
    save_data = {}
    for p_str, res in all_results.items():
        save_entry = {
            "p": res["p"],
            "n": res["n"],
            "num_d11": res["num_d11"],
            "num_d12": res["num_d12"],
            "num_orbits": res["num_orbits"],
            "marginal_cdf": res["marginal_cdf"],
            "orbit_summaries": res["orbit_summaries"],
        }
        save_data[p_str] = save_entry

    with open(outpath, "w") as f:
        json.dump(save_data, f, indent=2, default=str)
    print(f"\nResults saved to {outpath}")


if __name__ == "__main__":
    main()
