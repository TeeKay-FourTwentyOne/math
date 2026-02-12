#!/usr/bin/env python3
"""
Invariant hunt: identify what structural invariant of D11 determines N(D11).

For primes p in {19, 23} (p ≡ 3 mod 4), enumerate all symmetric D11,
filter to A-flat, exhaustively count valid D12 for each, and compute
a battery of structural invariants to find what predicts N(D11).
"""

import numpy as np
import json
import time
import os
from itertools import combinations
from math import comb, sqrt, floor
from collections import defaultdict


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


def quadratic_residues(p):
    """Return set of quadratic residues mod p."""
    return {(x * x) % p for x in range(1, p)}


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
    """Find lexicographically smallest orbit member under multiplication by Z_p*.

    Since D11 is symmetric (closed under negation), the orbit under Z_p*
    has size at most (p-1)/2 (since negation is already in the stabilizer).
    """
    g = primitive_root(p)
    best = tuple(sorted(D11_set))
    # Generate orbit under multiplication by g^2 (since g^((p-1)/2) = -1 already in stabilizer)
    # Actually, orbit under full Z_p* but D11 is symmetric so g and -g give same set
    multiplier = 1
    for _ in range((p - 1) // 2):
        multiplier = (multiplier * g) % p
        transformed = frozenset((d * multiplier) % p for d in D11_set)
        candidate = tuple(sorted(transformed))
        if candidate < best:
            best = candidate
    return best


def run_prime(p):
    """Run the full invariant hunt for prime p."""
    print(f"\n{'=' * 80}")
    print(f"INVARIANT HUNT: p = {p}")
    print(f"{'=' * 80}")

    n = (p + 1) // 2
    d11_size = n              # |D11| = (p+1)/2
    d12_size = n - 1          # |D12| = (p-1)/2
    k = d12_size - 1          # |S| = (p-3)/2
    thresh_binding = (p - 3) // 2   # A(d)+B(d) <= this for d in D11
    thresh_loose = (p + 3) // 2     # A(d)+B(d) <= this for d in D22
    E_A = (p + 1) / 4
    floor_E_A = floor(E_A)
    mu_B = (p - 3) / 4
    sigma = sqrt((p - 3) * (p + 1) / (16 * (p - 2)))

    QR = quadratic_residues(p)
    QNR = set(range(1, p)) - QR

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
    print(f"  #QR = {len(QR)}, QR = {sorted(QR)}")
    print(f"  #QNR = {len(QNR)}, QNR = {sorted(QNR)}")

    t0 = time.time()

    # Precompute all D12 autocorrelations
    print(f"\n  Precomputing {num_d12} D12 autocorrelations...")
    d12_matrix = build_d12_candidates(p, k)
    B_matrix = batch_autocorrelation(d12_matrix, p)
    print(f"  Done in {time.time() - t0:.1f}s")

    # Enumerate all symmetric D11
    print(f"  Enumerating symmetric D11...")
    d11_list = list(enumerate_symmetric_D11(p, d11_size))
    assert len(d11_list) == num_d11

    # Compute A-values and filter to A-flat
    all_d11_data = []
    a_flat_data = []

    for d11_idx, D11 in enumerate(d11_list):
        A = autocorrelation_single(D11, p)
        D11_sorted = sorted(D11)
        D22 = set(range(1, p)) - D11
        D22_sorted = sorted(D22)

        A_at_D11 = [int(A[d]) for d in D11_sorted]
        A_at_D22 = [int(A[d]) for d in D22_sorted]
        max_A_D11 = max(A_at_D11)
        max_A_D22 = max(A_at_D22)
        max_A_all = max(max_A_D11, max_A_D22)

        is_a_flat = (max_A_D11 <= floor_E_A)

        entry = {
            "D11": D11_sorted,
            "D22": D22_sorted,
            "A_at_D11": A_at_D11,
            "A_at_D22": A_at_D22,
            "A_profile_D11": sorted(A_at_D11),
            "A_profile_D22": sorted(A_at_D22),
            "max_A_D11": max_A_D11,
            "max_A_D22": max_A_D22,
            "max_A_all": max_A_all,
            "sum_A_D11": sum(A_at_D11),
            "sum_A_D22": sum(A_at_D22),
            "is_a_flat": is_a_flat,
        }

        all_d11_data.append(entry)
        if is_a_flat:
            a_flat_data.append(entry)

    print(f"  Total D11: {len(all_d11_data)}, A-flat: {len(a_flat_data)}")

    # For A-flat D11: compute N(D11) and all invariants
    print(f"\n  Computing N(D11) for {len(a_flat_data)} A-flat D11...")
    t1 = time.time()

    for idx, entry in enumerate(a_flat_data):
        D11 = set(entry["D11"])
        D22 = set(entry["D22"])
        d11_indices = np.array(entry["D11"], dtype=np.int32)
        d22_indices = np.array(entry["D22"], dtype=np.int32)

        # Count valid D12
        F_binding = entry["A_at_D11"] + B_matrix[:, d11_indices]  # broadcasting
        # A_at_D11 is a list, need to convert
        A_d11_arr = np.array(entry["A_at_D11"], dtype=np.int32)
        A_d22_arr = np.array(entry["A_at_D22"], dtype=np.int32)

        F_binding = A_d11_arr[np.newaxis, :] + B_matrix[:, d11_indices]
        valid_binding = np.all(F_binding <= thresh_binding, axis=1)

        F_loose = A_d22_arr[np.newaxis, :] + B_matrix[:, d22_indices]
        valid_loose = np.all(F_loose <= thresh_loose, axis=1)

        valid_mask = valid_binding & valid_loose
        N = int(valid_mask.sum())

        entry["N"] = N

        # QR partition
        qr_in_D11 = len(D11 & QR)
        qnr_in_D11 = len(D11 & QNR)
        entry["qr_in_D11"] = qr_in_D11
        entry["qnr_in_D11"] = qnr_in_D11

        # Number of "high" D22 positions
        for thresh_name, thresh_val in [("E_A", floor_E_A), ("E_A+1", floor_E_A + 1),
                                         ("E_A+2", floor_E_A + 2), ("E_A+3", floor_E_A + 3)]:
            entry[f"num_D22_above_{thresh_name}"] = sum(1 for a in entry["A_at_D22"] if a > thresh_val)

        # Orbit representative
        orbit_rep = canonical_orbit_rep(D11, p)
        entry["orbit_rep"] = list(orbit_rep)

        # z-scores
        z_scores = []
        for d in entry["D11"]:
            T_d = thresh_binding - entry["A_at_D11"][entry["D11"].index(d)]
            z = (T_d - mu_B) / sigma
            z_scores.append(z)
        for d in entry["D22"]:
            T_d = thresh_loose - entry["A_at_D22"][entry["D22"].index(d)]
            z = (T_d - mu_B) / sigma
            z_scores.append(z)

        # z-scores at D11 positions only (these are the binding constraints)
        z_D11 = []
        for i, d in enumerate(entry["D11"]):
            T_d = thresh_binding - entry["A_at_D11"][i]
            z = (T_d - mu_B) / sigma
            z_D11.append(z)

        # z-scores at D22 positions only
        z_D22 = []
        for i, d in enumerate(entry["D22"]):
            T_d = thresh_loose - entry["A_at_D22"][i]
            z = (T_d - mu_B) / sigma
            z_D22.append(z)

        # S1 = sum of ALL z_i (both D11 and D22 representatives)
        # But since D11 is symmetric: each pair {d, p-d} gives same z,
        # so representative sum = half of sum over all
        # Actually let's compute S1 over all (p-1) positions
        S1_all = sum(z_scores)
        S1_D11 = sum(z_D11)
        S1_D22 = sum(z_D22)

        # S2 = sum of z_i^2
        S2_all = sum(z ** 2 for z in z_scores)
        S2_D11 = sum(z ** 2 for z in z_D11)
        S2_D22 = sum(z ** 2 for z in z_D22)

        # Expected S1 = (p-7)/(4*sigma) from constant z-sum theorem
        expected_S1 = (p - 7) / (4 * sigma)

        entry["z_D11"] = [round(z, 6) for z in z_D11]
        entry["z_D22"] = [round(z, 6) for z in z_D22]
        entry["S1_all"] = round(S1_all, 6)
        entry["S1_D11"] = round(S1_D11, 6)
        entry["S1_D22"] = round(S1_D22, 6)
        entry["S2_all"] = round(S2_all, 6)
        entry["S2_D11"] = round(S2_D11, 6)
        entry["S2_D22"] = round(S2_D22, 6)
        entry["expected_S1"] = round(expected_S1, 6)

        # min z at D11 positions (tightest binding constraint)
        entry["min_z_D11"] = round(min(z_D11), 6)
        entry["min_z_D22"] = round(min(z_D22), 6)

        # Product of Phi(z_i) at D11 positions (Gaussian proxy for N)
        from scipy.stats import norm as norm_dist
        log_prod_phi_D11 = sum(norm_dist.logcdf(z) for z in z_D11)
        entry["log_prod_phi_D11"] = round(log_prod_phi_D11, 6)

        if (idx + 1) % 10 == 0 or idx == len(a_flat_data) - 1:
            print(f"    [{idx+1}/{len(a_flat_data)}] D11={entry['D11'][:4]}... N={N}")

    elapsed = time.time() - t0
    print(f"  Total time: {elapsed:.1f}s")

    # =========================================================================
    # ANALYSIS
    # =========================================================================
    print(f"\n{'=' * 80}")
    print(f"ANALYSIS FOR p = {p}")
    print(f"{'=' * 80}")

    # Group by N-value
    n_groups = defaultdict(list)
    for entry in a_flat_data:
        n_groups[entry["N"]].append(entry)

    print(f"\n--- A-flat D11 grouped by N(D11) ---")
    print(f"{'N':>8} {'count':>6} {'QR range':>12} {'maxA_D22':>10} {'S2_D11':>12} {'min_z_D11':>12} {'log_prod_phi':>14}")
    print("-" * 80)
    for N_val in sorted(n_groups.keys()):
        group = n_groups[N_val]
        qr_vals = [e["qr_in_D11"] for e in group]
        maxA_D22_vals = [e["max_A_D22"] for e in group]
        S2_D11_vals = [e["S2_D11"] for e in group]
        min_z_vals = [e["min_z_D11"] for e in group]
        log_phi_vals = [e["log_prod_phi_D11"] for e in group]

        qr_range = f"[{min(qr_vals)},{max(qr_vals)}]"
        mA_range = f"[{min(maxA_D22_vals)},{max(maxA_D22_vals)}]"
        S2_range = f"[{min(S2_D11_vals):.2f},{max(S2_D11_vals):.2f}]"
        mz_range = f"[{min(min_z_vals):.3f},{max(min_z_vals):.3f}]"
        lp_range = f"[{min(log_phi_vals):.2f},{max(log_phi_vals):.2f}]"

        print(f"{N_val:>8} {len(group):>6} {qr_range:>12} {mA_range:>10} {S2_range:>12} {mz_range:>12} {lp_range:>14}")

    # Detailed profiles per group
    print(f"\n--- D11 A-profiles and D22 A-profiles per N-group ---")
    for N_val in sorted(n_groups.keys()):
        group = n_groups[N_val]
        print(f"\n  N = {N_val} ({len(group)} D11):")
        for e in group:
            print(f"    D11={e['D11']}")
            print(f"      A@D11={e['A_profile_D11']}, sum={e['sum_A_D11']}")
            print(f"      A@D22={e['A_profile_D22']}, sum={e['sum_A_D22']}, max={e['max_A_D22']}")
            print(f"      QR={e['qr_in_D11']}, QNR={e['qnr_in_D11']}")
            print(f"      S1_D11={e['S1_D11']}, S2_D11={e['S2_D11']}")
            print(f"      min_z_D11={e['min_z_D11']}, log_prod_phi={e['log_prod_phi_D11']:.4f}")
            print(f"      z@D11={e['z_D11']}")
            # High D22 positions
            high_str = ", ".join(f">{t}:{e[f'num_D22_above_{t}']}"
                                 for t in ["E_A", "E_A+1", "E_A+2", "E_A+3"])
            print(f"      D22 high: {high_str}")

    # =========================================================================
    # CROSS-TABULATIONS
    # =========================================================================
    print(f"\n{'=' * 80}")
    print(f"CROSS-TABULATIONS FOR p = {p}")
    print(f"{'=' * 80}")

    # 1. QR count vs N
    print(f"\n--- QR count in D11 vs N ---")
    qr_n_table = defaultdict(lambda: defaultdict(int))
    for e in a_flat_data:
        qr_n_table[e["qr_in_D11"]][e["N"]] += 1
    all_N = sorted(set(e["N"] for e in a_flat_data))
    all_QR = sorted(set(e["qr_in_D11"] for e in a_flat_data))
    header = f"{'QR':>6}" + "".join(f"{'N='+str(n):>8}" for n in all_N) + f"{'total':>8}"
    print(header)
    for qr in all_QR:
        row = f"{qr:>6}"
        total = 0
        for n in all_N:
            c = qr_n_table[qr][n]
            total += c
            row += f"{c:>8}"
        row += f"{total:>8}"
        print(row)

    # 2. max A at D22 vs N
    print(f"\n--- Max A(d) at D22 positions vs N ---")
    maxA_n_table = defaultdict(lambda: defaultdict(int))
    for e in a_flat_data:
        maxA_n_table[e["max_A_D22"]][e["N"]] += 1
    all_maxA = sorted(set(e["max_A_D22"] for e in a_flat_data))
    header = f"{'maxA':>6}" + "".join(f"{'N='+str(n):>8}" for n in all_N) + f"{'total':>8}"
    print(header)
    for ma in all_maxA:
        row = f"{ma:>6}"
        total = 0
        for n in all_N:
            c = maxA_n_table[ma][n]
            total += c
            row += f"{c:>8}"
        row += f"{total:>8}"
        print(row)

    # 3. D22 A-profile vs N
    print(f"\n--- D22 A-profile (sorted) vs N ---")
    profile_n_table = defaultdict(lambda: defaultdict(int))
    for e in a_flat_data:
        key = tuple(e["A_profile_D22"])
        profile_n_table[key][e["N"]] += 1
    print(f"{'D22 profile':>40}" + "".join(f"{'N='+str(n):>8}" for n in all_N))
    for profile in sorted(profile_n_table.keys()):
        row = f"{str(list(profile)):>40}"
        for n in all_N:
            row += f"{profile_n_table[profile][n]:>8}"
        print(row)

    # 4. D11 A-profile vs N
    print(f"\n--- D11 A-profile (sorted) vs N ---")
    d11_profile_n_table = defaultdict(lambda: defaultdict(int))
    for e in a_flat_data:
        key = tuple(e["A_profile_D11"])
        d11_profile_n_table[key][e["N"]] += 1
    print(f"{'D11 profile':>50}" + "".join(f"{'N='+str(n):>8}" for n in all_N))
    for profile in sorted(d11_profile_n_table.keys()):
        row = f"{str(list(profile)):>50}"
        for n in all_N:
            row += f"{d11_profile_n_table[profile][n]:>8}"
        print(row)

    # 5. S2_D11 (rounded) vs N
    print(f"\n--- S2_D11 (rounded to 0.1) vs N ---")
    s2_n_table = defaultdict(lambda: defaultdict(int))
    for e in a_flat_data:
        key = round(e["S2_D11"], 1)
        s2_n_table[key][e["N"]] += 1
    all_s2 = sorted(s2_n_table.keys())
    header = f"{'S2_D11':>8}" + "".join(f"{'N='+str(n):>8}" for n in all_N)
    print(header)
    for s2 in all_s2:
        row = f"{s2:>8.1f}"
        for n in all_N:
            row += f"{s2_n_table[s2][n]:>8}"
        print(row)

    # 6. num_D22_above_E_A vs N
    print(f"\n--- #D22 positions with A(d) > floor(E_A) vs N ---")
    high_n_table = defaultdict(lambda: defaultdict(int))
    for e in a_flat_data:
        key = e["num_D22_above_E_A"]
        high_n_table[key][e["N"]] += 1
    all_high = sorted(high_n_table.keys())
    header = f"{'#high':>6}" + "".join(f"{'N='+str(n):>8}" for n in all_N)
    print(header)
    for h in all_high:
        row = f"{h:>6}"
        for n in all_N:
            row += f"{high_n_table[h][n]:>8}"
        print(row)

    # 7. num_D22_above_E_A+1 vs N
    print(f"\n--- #D22 positions with A(d) > floor(E_A)+1 vs N ---")
    high_n_table2 = defaultdict(lambda: defaultdict(int))
    for e in a_flat_data:
        key = e["num_D22_above_E_A+1"]
        high_n_table2[key][e["N"]] += 1
    all_high2 = sorted(high_n_table2.keys())
    header = f"{'#high':>6}" + "".join(f"{'N='+str(n):>8}" for n in all_N)
    print(header)
    for h in all_high2:
        row = f"{h:>6}"
        for n in all_N:
            row += f"{high_n_table2[h][n]:>8}"
        print(row)

    # 8. COMBINED: (D11 profile, D22 profile) vs N -- perfect predictor test
    print(f"\n--- (D11 profile, D22 profile) vs N ---")
    combo_n_table = defaultdict(lambda: defaultdict(int))
    for e in a_flat_data:
        key = (tuple(e["A_profile_D11"]), tuple(e["A_profile_D22"]))
        combo_n_table[key][e["N"]] += 1
    print(f"  Distinct (D11,D22) profile pairs: {len(combo_n_table)}")
    print(f"  Perfect predictor? {all(len(v) == 1 for v in combo_n_table.values())}")
    for combo in sorted(combo_n_table.keys()):
        n_counts = combo_n_table[combo]
        if len(n_counts) > 1:
            print(f"  AMBIGUOUS: D11={list(combo[0])}, D22={list(combo[1])} -> N values: {dict(n_counts)}")

    # 9. Full A-value vector (A(1),...,A(p-1)) profile (sorted) vs N
    print(f"\n--- Full A-profile sorted vs N ---")
    full_profile_n_table = defaultdict(lambda: defaultdict(int))
    for e in a_flat_data:
        full = sorted(e["A_profile_D11"] + e["A_profile_D22"])
        key = tuple(full)
        full_profile_n_table[key][e["N"]] += 1
    print(f"  Distinct full profiles: {len(full_profile_n_table)}")
    print(f"  Perfect predictor? {all(len(v) == 1 for v in full_profile_n_table.values())}")

    # 10. Orbit analysis
    print(f"\n--- Orbit analysis ---")
    orbit_n_table = defaultdict(lambda: defaultdict(int))
    for e in a_flat_data:
        key = tuple(e["orbit_rep"])
        orbit_n_table[key][e["N"]] += 1
    print(f"  Distinct orbits among A-flat D11: {len(orbit_n_table)}")
    for orbit in sorted(orbit_n_table.keys()):
        n_counts = orbit_n_table[orbit]
        if len(n_counts) > 1:
            print(f"  MULTI-N orbit: rep={list(orbit)[:6]}... -> {dict(n_counts)}")
    print(f"  All orbits have unique N? {all(len(v) == 1 for v in orbit_n_table.values())}")

    # =========================================================================
    # SUMMARY: Best candidate invariant
    # =========================================================================
    print(f"\n{'=' * 80}")
    print(f"SUMMARY: BEST CANDIDATE INVARIANTS FOR p = {p}")
    print(f"{'=' * 80}")

    # Check each candidate for whether it perfectly predicts N
    candidates = {
        "QR count": lambda e: e["qr_in_D11"],
        "max A at D22": lambda e: e["max_A_D22"],
        "D22 A-profile": lambda e: tuple(e["A_profile_D22"]),
        "D11 A-profile": lambda e: tuple(e["A_profile_D11"]),
        "(D11,D22) profiles": lambda e: (tuple(e["A_profile_D11"]), tuple(e["A_profile_D22"])),
        "Full A-profile": lambda e: tuple(sorted(e["A_profile_D11"] + e["A_profile_D22"])),
        "S2_D11": lambda e: round(e["S2_D11"], 4),
        "#D22 above E_A": lambda e: e["num_D22_above_E_A"],
        "#D22 above E_A+1": lambda e: e["num_D22_above_E_A+1"],
        "min z at D11": lambda e: round(e["min_z_D11"], 4),
        "log prod phi D11": lambda e: round(e["log_prod_phi_D11"], 2),
        "orbit rep": lambda e: tuple(e["orbit_rep"]),
    }

    for name, func in candidates.items():
        inv_to_N = defaultdict(set)
        for e in a_flat_data:
            inv_to_N[func(e)].add(e["N"])
        perfect = all(len(v) == 1 for v in inv_to_N.values())
        n_distinct = len(inv_to_N)
        n_N_vals = len(set(e["N"] for e in a_flat_data))
        max_ambiguity = max(len(v) for v in inv_to_N.values())
        print(f"  {name:>25s}: distinct={n_distinct:>4}, perfect={perfect}, max_ambiguity={max_ambiguity}")
        if not perfect:
            # Show ambiguous cases
            for inv_val, n_vals in inv_to_N.items():
                if len(n_vals) > 1:
                    inv_str = str(inv_val) if len(str(inv_val)) < 50 else str(inv_val)[:47] + "..."
                    print(f"      AMBIG: {inv_str} -> N in {sorted(n_vals)}")

    # Check: does N>0 have a simple separator?
    print(f"\n--- Separating N>0 from N=0 ---")
    for name, func in candidates.items():
        vals_pos = set(func(e) for e in a_flat_data if e["N"] > 0)
        vals_zero = set(func(e) for e in a_flat_data if e["N"] == 0)
        overlap = vals_pos & vals_zero
        if not overlap:
            print(f"  {name}: PERFECTLY separates N>0 from N=0")
        else:
            print(f"  {name}: overlap = {len(overlap)} values")

    return {
        "p": p,
        "n": n,
        "num_d11": num_d11,
        "num_d12": num_d12,
        "num_a_flat": len(a_flat_data),
        "a_flat_data": a_flat_data,
        "n_groups": {str(k): len(v) for k, v in n_groups.items()},
    }


def main():
    print("=" * 80)
    print("INVARIANT HUNT: What structural invariant of D11 determines N(D11)?")
    print("=" * 80)

    all_results = {}

    for p in [19, 23]:
        results = run_prime(p)
        all_results[str(p)] = results

    # Save results
    outpath = "/Users/stephenpadgett/Projects/math/ramsey-book-graphs/invariant_hunt_results.json"

    # Prepare serializable version
    save_data = {}
    for p_str, res in all_results.items():
        save_entry = {
            "p": res["p"],
            "n": res["n"],
            "num_d11": res["num_d11"],
            "num_d12": res["num_d12"],
            "num_a_flat": res["num_a_flat"],
            "n_groups": res["n_groups"],
            "a_flat_entries": [],
        }
        for e in res["a_flat_data"]:
            save_e = {k: v for k, v in e.items()
                      if k not in ("D22",)}  # skip large redundant fields
            save_entry["a_flat_entries"].append(save_e)
        save_data[p_str] = save_entry

    with open(outpath, "w") as f:
        json.dump(save_data, f, indent=2, default=str)
    print(f"\nResults saved to {outpath}")


if __name__ == "__main__":
    main()
