#!/usr/bin/env python3
"""Enhanced sampling for p=31,43,47,59 with A-flat and QR-balanced D11.

Key improvements over sampling_d12_larger_primes.py:
1. Greedy construction of A-flat D11 (minimize max A(d))
2. QR-balanced D11 (|D11 ∩ QR| = |D11 ∩ QNR|)
3. Separate binding-only vs full-constraint tracking
4. D22 marginal tracking
5. Known solutions from registry for all primes
6. Gaussian first-moment estimates even when 0 hits

Usage: python -u sampling_aflat_d11.py [--primes 31 43 47 59] [--samples 1000000]
"""

import numpy as np
import json
import time
import os
import sys
from math import comb, log, log2, pi, sqrt, exp
from collections import defaultdict
import argparse


def load_all_known_solutions():
    """Load all known solutions from the registry + hardcoded ones."""
    solutions = {
        31: {
            "D11": {6, 7, 8, 10, 11, 12, 14, 15, 16, 17, 19, 20, 21, 23, 24, 25},
            "D12": {0, 1, 2, 3, 8, 11, 12, 13, 15, 18, 20, 21, 24, 27, 29},
        },
    }
    script_dir = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(script_dir, "solutions_registry.json")
    if not os.path.exists(path):
        return solutions
    with open(path) as f:
        registry = json.load(f)
    for sol in registry["solutions"]:
        p = sol["m"]
        D11 = set(sol["D11"])
        D12 = set(sol["D12"])
        target = (p + 1) // 2
        # Normalize to |D11| = (p+1)/2 convention
        if len(D11) != target:
            D22 = set(range(1, p)) - D11
            D12T = {(-x) % p for x in D12}
            D11, D12 = D22, D12T
        solutions[p] = {"D11": D11, "D12": D12}
    return solutions


def quadratic_residues(p):
    """Return QR mod p (excluding 0)."""
    return set((x * x) % p for x in range(1, p)) - {0}


def autocorrelation_single(S, p):
    """Compute A(d) for a single set S in Z_p."""
    indicator = np.zeros(p, dtype=np.float64)
    for j in S:
        indicator[j] = 1.0
    fft_val = np.fft.fft(indicator)
    autocorr = np.fft.ifft(np.abs(fft_val) ** 2).real
    return np.round(autocorr).astype(np.int32)


def batch_autocorrelation(indicator_matrix):
    """Compute Delta(S,S,d) for batch of sets via FFT."""
    fft_vals = np.fft.fft(indicator_matrix, axis=1)
    autocorr = np.fft.ifft(np.abs(fft_vals) ** 2, axis=1).real
    return np.round(autocorr).astype(np.int32)


def symmetric_pairs(p):
    """Return list of {d, p-d} pairs for d=1..p-1."""
    pairs = []
    seen = set()
    for d in range(1, p):
        if d not in seen:
            comp = (p - d) % p
            pairs.append((d, comp))
            seen.add(d)
            seen.add(comp)
    return pairs


def greedy_aflat_d11(p, rng, num_attempts=50):
    """Construct A-flat symmetric D11 by greedy pair selection.

    Try multiple random orderings and pick the one with smallest max A(d).
    """
    n = (p + 1) // 2
    num_pairs = n // 2
    pairs = symmetric_pairs(p)

    best_D11 = None
    best_max_A = float('inf')

    for _ in range(num_attempts):
        # Random permutation of pairs
        order = rng.permutation(len(pairs))
        D11 = set()
        chosen = []

        for idx in order:
            if len(chosen) >= num_pairs:
                break
            d, comp = pairs[idx]
            # Tentatively add this pair
            trial = D11 | {d, comp}
            A = autocorrelation_single(trial, p)
            max_A = max(int(A[x]) for x in trial)
            # Accept if max_A is reasonable (greedy criterion)
            chosen.append(idx)
            D11 = trial

        A = autocorrelation_single(D11, p)
        max_A = max(int(A[d]) for d in D11)
        if max_A < best_max_A:
            best_max_A = max_A
            best_D11 = D11

    return best_D11, best_max_A


def greedy_aflat_d11_incremental(p, rng, num_attempts=100):
    """Construct A-flat D11 by incremental greedy: at each step, add the pair
    that minimizes the resulting max A(d) over D11 positions only."""
    n = (p + 1) // 2
    num_pairs = n // 2
    pairs = symmetric_pairs(p)

    best_D11 = None
    best_max_A = float('inf')

    for attempt in range(num_attempts):
        D11 = set()
        available = list(range(len(pairs)))
        rng.shuffle(available)

        for step in range(num_pairs):
            if len(available) > 20:
                candidates = rng.choice(available, size=min(20, len(available)), replace=False).tolist()
            else:
                candidates = list(available)

            best_pair_idx = None
            best_pair_maxA = float('inf')

            for idx in candidates:
                d, comp = pairs[idx]
                trial = D11 | {d, comp}
                A = autocorrelation_single(trial, p)
                mA = max(int(A[x]) for x in trial)
                if mA < best_pair_maxA:
                    best_pair_maxA = mA
                    best_pair_idx = idx

            d, comp = pairs[best_pair_idx]
            D11 = D11 | {d, comp}
            available.remove(best_pair_idx)

        A = autocorrelation_single(D11, p)
        max_A = max(int(A[d]) for d in D11)
        if max_A < best_max_A:
            best_max_A = max_A
            best_D11 = D11

    return best_D11, best_max_A


def greedy_global_aflat_d11(p, rng, num_attempts=100):
    """Construct D11 minimizing max A(d) over ALL d=1..p-1, not just D11.

    This is the critical metric: D22 positions with high A(d) are the bottleneck.
    """
    n = (p + 1) // 2
    num_pairs = n // 2
    pairs = symmetric_pairs(p)

    best_D11 = None
    best_max_A_all = float('inf')

    for attempt in range(num_attempts):
        D11 = set()
        available = list(range(len(pairs)))
        rng.shuffle(available)

        for step in range(num_pairs):
            if len(available) > 20:
                candidates = rng.choice(available, size=min(20, len(available)), replace=False).tolist()
            else:
                candidates = list(available)

            best_pair_idx = None
            best_pair_score = float('inf')

            for idx in candidates:
                d, comp = pairs[idx]
                trial = D11 | {d, comp}
                A = autocorrelation_single(trial, p)
                # Score: max A(d) over ALL nonzero positions
                score = max(int(A[d]) for d in range(1, p))
                if score < best_pair_score:
                    best_pair_score = score
                    best_pair_idx = idx

            d, comp = pairs[best_pair_idx]
            D11 = D11 | {d, comp}
            available.remove(best_pair_idx)

        A = autocorrelation_single(D11, p)
        max_A_all = max(int(A[d]) for d in range(1, p))
        if max_A_all < best_max_A_all:
            best_max_A_all = max_A_all
            best_D11 = D11

    return best_D11, best_max_A_all


def qr_balanced_d11(p, rng, num_attempts=50):
    """Generate QR-balanced symmetric D11 (|D11 ∩ QR| = |D11 ∩ QNR|).

    For p ≡ 3 mod 4, QR is not symmetric, so we need to construct D11
    that is symmetric AND has equal QR/QNR representation.
    """
    n = (p + 1) // 2
    num_pairs = n // 2
    pairs = symmetric_pairs(p)
    QR = quadratic_residues(p)

    # Classify pairs by their QR composition
    # Each pair {d, p-d}: both QR, both QNR, or mixed
    # For p ≡ 3 mod 4: d is QR iff p-d is QNR (since -1 is QNR)
    # So every pair is mixed! This means any symmetric D11 is QR-balanced.

    # For p ≡ 1 mod 4: d is QR iff p-d is QR (since -1 is QR)
    # So pairs are either both-QR or both-QNR.

    best_D11 = None
    best_max_A = float('inf')

    for _ in range(num_attempts):
        order = rng.permutation(len(pairs))
        D11 = set()
        count = 0
        for idx in order:
            if count >= num_pairs:
                break
            d, comp = pairs[idx]
            D11.add(d)
            D11.add(comp)
            count += 1

        # Check QR balance
        qr_count = len(D11 & QR)
        qnr_count = len(D11 - QR)
        if qr_count != qnr_count:
            continue  # Skip non-balanced (only relevant for p ≡ 1 mod 4)

        A = autocorrelation_single(D11, p)
        max_A = max(int(A[d]) for d in D11)
        if max_A < best_max_A:
            best_max_A = max_A
            best_D11 = D11

    if best_D11 is None:
        # Fallback: for p ≡ 3 mod 4, all symmetric D11 are QR-balanced
        D11 = set()
        order = rng.permutation(len(pairs))
        for idx in order[:num_pairs]:
            d, comp = pairs[idx]
            D11.add(d)
            D11.add(comp)
        A = autocorrelation_single(D11, p)
        best_max_A = max(int(A[d]) for d in D11)
        best_D11 = D11

    return best_D11, best_max_A


def random_symmetric_D11(p, rng):
    """Generate a random symmetric D11 of size (p+1)/2."""
    n = (p + 1) // 2
    num_pairs = n // 2
    pairs = symmetric_pairs(p)
    chosen = rng.choice(len(pairs), size=num_pairs, replace=False)
    D11 = set()
    for i in chosen:
        D11.add(pairs[i][0])
        D11.add(pairs[i][1])
    return D11


def generate_random_d12_batch(p, k, batch_size, rng):
    """Generate batch of random D12 = {0} ∪ k elements from {1,..,p-1}."""
    matrix = np.zeros((batch_size, p), dtype=np.float64)
    matrix[:, 0] = 1.0
    candidates = np.arange(1, p)
    for i in range(batch_size):
        S = rng.choice(candidates, size=k, replace=False)
        matrix[i, S] = 1.0
    return matrix


def gaussian_first_moment(A, D11, D22, p):
    """Estimate E[#valid D12] using Gaussian approximation.

    For each constraint position d, B(d) ~ N(mu, sigma^2) approximately.
    E[B(d)] = (p-3)(p-1)/(4(p-2)) for |D12|=(p-1)/2
    Var[B(d)] = (p-3)(p+1)/(16(p-2))
    """
    k = (p - 3) // 2  # |D12 \ {0}|
    n12 = k + 1  # |D12| = (p-1)/2

    # Exact moments
    mu_B = (p - 3) * (p - 1) / (4 * (p - 2))
    var_B = (p - 3) * (p + 1) / (16 * (p - 2))
    sigma_B = sqrt(var_B)

    threshold_binding = (p - 3) // 2
    threshold_loose = (p + 3) // 2

    # For each constraint, P[A(d) + B(d) <= threshold] ~ Phi((threshold - A(d) - mu_B) / sigma_B)
    from scipy.stats import norm

    log_prob_sum = 0.0
    details = {"binding": [], "loose": []}

    for d in sorted(D11):
        Ad = int(A[d])
        slack = threshold_binding - Ad - mu_B
        z = slack / sigma_B
        lp = norm.logcdf(z)
        log_prob_sum += lp
        details["binding"].append({"d": d, "A": Ad, "slack": float(slack),
                                    "z": float(z), "log_p": float(lp)})

    for d in sorted(D22):
        Ad = int(A[d])
        slack = threshold_loose - Ad - mu_B
        z = slack / sigma_B
        lp = norm.logcdf(z)
        log_prob_sum += lp
        details["loose"].append({"d": d, "A": Ad, "slack": float(slack),
                                  "z": float(z), "log_p": float(lp)})

    num_d12_total = comb(p - 1, k)
    log_E_valid = log_prob_sum + log(num_d12_total)

    return {
        "mu_B": float(mu_B),
        "sigma_B": float(sigma_B),
        "log_prod_marginals": float(log_prob_sum),
        "log_C_pk": float(log(num_d12_total)),
        "log_E_valid": float(log_E_valid),
        "E_valid_approx": float(exp(log_E_valid)) if log_E_valid < 500 else "exp({:.1f})".format(log_E_valid),
        "binding_details": details["binding"],
        "loose_details": details["loose"],
    }


def sample_d12_for_d11(D11, p, A, num_samples, batch_size, rng):
    """Sample D12, track binding-only and full constraints separately."""
    k = (p - 3) // 2
    threshold_binding = (p - 3) // 2
    threshold_loose = (p + 3) // 2

    d11_arr = np.array(sorted(D11), dtype=np.int32)
    D22 = set(range(1, p)) - D11
    d22_arr = np.array(sorted(D22), dtype=np.int32)

    A_at_d11 = A[d11_arr]
    A_at_d22 = A[d22_arr]

    num_valid_full = 0
    num_valid_binding = 0
    binding_marginal_hits = np.zeros(len(d11_arr), dtype=np.int64)
    loose_marginal_hits = np.zeros(len(d22_arr), dtype=np.int64)
    total_sampled = 0

    while total_sampled < num_samples:
        bs = min(batch_size, num_samples - total_sampled)
        d12_batch = generate_random_d12_batch(p, k, bs, rng)
        B_batch = batch_autocorrelation(d12_batch)

        # Binding constraints
        F_binding = A_at_d11[np.newaxis, :] + B_batch[:, d11_arr]
        binding_ok = F_binding <= threshold_binding
        binding_marginal_hits += binding_ok.sum(axis=0)
        all_binding_ok = np.all(binding_ok, axis=1)
        num_valid_binding += int(all_binding_ok.sum())

        # Loose constraints
        if d22_arr.size > 0:
            F_loose = A_at_d22[np.newaxis, :] + B_batch[:, d22_arr]
            loose_ok = F_loose <= threshold_loose
            loose_marginal_hits += loose_ok.sum(axis=0)
            all_loose_ok = np.all(loose_ok, axis=1)
            valid_mask = all_binding_ok & all_loose_ok
        else:
            valid_mask = all_binding_ok

        num_valid_full += int(valid_mask.sum())
        total_sampled += bs

    num_d12_total = comb(p - 1, k)

    binding_rate = num_valid_binding / total_sampled
    full_rate = num_valid_full / total_sampled

    binding_marginals = binding_marginal_hits / total_sampled
    loose_marginals = loose_marginal_hits / total_sampled

    prod_binding = float(np.prod(binding_marginals)) if np.all(binding_marginals > 0) else 0.0
    prod_loose = float(np.prod(loose_marginals)) if np.all(loose_marginals > 0) else 0.0
    prod_all = prod_binding * prod_loose if prod_binding > 0 and prod_loose > 0 else 0.0

    ratio_binding = binding_rate / prod_binding if prod_binding > 0 else float('inf')
    ratio_full = full_rate / prod_all if prod_all > 0 else (float('inf') if num_valid_full > 0 else 0.0)

    return {
        "num_valid_full": num_valid_full,
        "num_valid_binding": num_valid_binding,
        "total_sampled": total_sampled,
        "binding_rate": binding_rate,
        "full_rate": full_rate,
        "E_valid_binding": binding_rate * num_d12_total,
        "E_valid_full": full_rate * num_d12_total,
        "binding_marginal_min": float(binding_marginals.min()),
        "binding_marginal_max": float(binding_marginals.max()),
        "binding_marginal_mean": float(binding_marginals.mean()),
        "loose_marginal_min": float(loose_marginals.min()) if len(loose_marginals) > 0 else 1.0,
        "loose_marginal_max": float(loose_marginals.max()) if len(loose_marginals) > 0 else 1.0,
        "loose_marginal_mean": float(loose_marginals.mean()) if len(loose_marginals) > 0 else 1.0,
        "ratio_binding": ratio_binding,
        "ratio_full": ratio_full,
        "log2_prod_binding": float(np.sum(np.log2(binding_marginals))) if np.all(binding_marginals > 0) else float('-inf'),
        "log2_prod_loose": float(np.sum(np.log2(loose_marginals))) if np.all(loose_marginals > 0) else float('-inf'),
    }


def analyze_d11(D11, p):
    """Compute structural properties of a D11."""
    QR = quadratic_residues(p)
    A = autocorrelation_single(D11, p)
    D22 = set(range(1, p)) - D11

    max_A_d11 = max(int(A[d]) for d in D11)
    max_A_d22 = max(int(A[d]) for d in D22) if D22 else 0
    max_A_all = max(max_A_d11, max_A_d22)
    mean_A_d11 = float(np.mean([int(A[d]) for d in D11]))
    mean_A_d22 = float(np.mean([int(A[d]) for d in D22])) if D22 else 0.0

    qr_count = len(D11 & QR)
    qnr_count = len(D11 - QR)

    return {
        "max_A_d11": max_A_d11,
        "max_A_d22": max_A_d22,
        "max_A_all": max_A_all,
        "mean_A_d11": mean_A_d11,
        "mean_A_d22": mean_A_d22,
        "qr_count": qr_count,
        "qnr_count": qnr_count,
        "qr_balanced": qr_count == qnr_count,
        "A_values_d11": {d: int(A[d]) for d in sorted(D11)},
    }


def run_prime(p, num_d12_samples=1000000, batch_size=10000, num_greedy=20,
              num_qr_balanced=10, num_random=20):
    """Run enhanced sampling for a single prime."""
    n = (p + 1) // 2
    k = (p - 3) // 2
    num_d12_total = comb(p - 1, k)

    print(f"\n{'=' * 80}")
    print(f"PRIME p = {p}, n = {n}")
    print(f"|D11| = {n}, |D12| = {n-1}, k = {k}")
    print(f"C({p-1},{k}) = {num_d12_total:,.0f} ({log2(num_d12_total):.1f} bits)")
    print(f"Sampling {num_d12_samples:,} D12 per D11")
    print(f"D11 sources: {num_greedy} greedy A-flat, {num_qr_balanced} QR-balanced, {num_random} random")
    print(f"{'=' * 80}")
    sys.stdout.flush()

    rng = np.random.default_rng(42)
    t0 = time.time()
    all_results = []

    # Load known solutions
    known_solutions = load_all_known_solutions()
    known = known_solutions.get(p)

    QR = quadratic_residues(p)
    print(f"p mod 4 = {p % 4}")
    print(f"|QR| = {len(QR)}")
    sys.stdout.flush()

    # 1. Test known solution D11
    if known:
        D11 = known["D11"]
        props = analyze_d11(D11, p)
        A = autocorrelation_single(D11, p)
        print(f"\n--- Known solution D11 ---")
        print(f"  max_A(D11)={props['max_A_d11']}, max_A(D22)={props['max_A_d22']}, "
              f"max_A(all)={props['max_A_all']}")
        print(f"  QR: {props['qr_count']}/{props['qnr_count']} "
              f"{'(balanced)' if props['qr_balanced'] else '(UNBALANCED)'}")
        sys.stdout.flush()

        stats = sample_d12_for_d11(D11, p, A, num_d12_samples, batch_size, rng)
        print(f"  Binding: {stats['num_valid_binding']}/{stats['total_sampled']} "
              f"= {stats['binding_rate']:.6f} (ratio={stats['ratio_binding']:.4f})")
        print(f"  Full:    {stats['num_valid_full']}/{stats['total_sampled']} "
              f"= {stats['full_rate']:.6f}")
        print(f"  Binding marginals: min={stats['binding_marginal_min']:.4f}, "
              f"max={stats['binding_marginal_max']:.4f}")
        print(f"  Loose marginals:   min={stats['loose_marginal_min']:.4f}, "
              f"max={stats['loose_marginal_max']:.4f}")
        print(f"  log2(prod_binding)={stats['log2_prod_binding']:.2f}, "
              f"log2(prod_loose)={stats['log2_prod_loose']:.2f}")

        # Gaussian estimate
        try:
            gauss = gaussian_first_moment(A, D11, set(range(1, p)) - D11, p)
            print(f"  Gaussian E[valid]: log={gauss['log_E_valid']:.2f}, "
                  f"E={gauss['E_valid_approx']}")
        except ImportError:
            gauss = None
            print(f"  (scipy not available for Gaussian estimate)")

        sys.stdout.flush()

        all_results.append({
            "type": "known",
            "D11": sorted(D11),
            **props,
            **stats,
            "gaussian": gauss,
        })

    # 2. Greedy A-flat D11 (D11-only optimization)
    print(f"\n--- Greedy A-flat D11 (D11-only, {num_greedy} attempts) ---")
    sys.stdout.flush()

    aflat_d11s = []
    for i in range(num_greedy):
        D11, max_A = greedy_aflat_d11_incremental(p, rng, num_attempts=20)
        aflat_d11s.append((D11, max_A))

    aflat_d11s.sort(key=lambda x: x[1])
    seen = set()
    unique_aflat = []
    for D11, max_A in aflat_d11s:
        key = tuple(sorted(D11))
        if key not in seen:
            seen.add(key)
            unique_aflat.append((D11, max_A))
    unique_aflat = unique_aflat[:5]  # Top 5

    print(f"  Best max_A(D11) values: {[mA for _, mA in unique_aflat]}")
    sys.stdout.flush()

    for i, (D11, max_A) in enumerate(unique_aflat):
        props = analyze_d11(D11, p)
        A = autocorrelation_single(D11, p)
        stats = sample_d12_for_d11(D11, p, A, num_d12_samples, batch_size, rng)

        print(f"  [{i}] max_A(D11)={max_A}, max_A(all)={props['max_A_all']}, "
              f"QR={'bal' if props['qr_balanced'] else 'no'}, "
              f"binding={stats['num_valid_binding']} (ratio={stats['ratio_binding']:.4f}), "
              f"full={stats['num_valid_full']}")

        try:
            gauss = gaussian_first_moment(A, D11, set(range(1, p)) - D11, p)
            print(f"      Gaussian log(E[valid])={gauss['log_E_valid']:.2f}")
        except ImportError:
            gauss = None
        sys.stdout.flush()

        all_results.append({
            "type": "greedy_d11_only",
            "D11": sorted(D11),
            **props,
            **stats,
            "gaussian": gauss,
        })

    # 2b. Greedy GLOBAL A-flat D11 (minimize max A(d) over ALL positions)
    print(f"\n--- Greedy GLOBAL A-flat D11 ({num_greedy} attempts) ---")
    sys.stdout.flush()

    gflat_d11s = []
    for i in range(num_greedy):
        D11, max_A_all = greedy_global_aflat_d11(p, rng, num_attempts=20)
        gflat_d11s.append((D11, max_A_all))

    gflat_d11s.sort(key=lambda x: x[1])
    seen = set()
    unique_gflat = []
    for D11, max_A_all in gflat_d11s:
        key = tuple(sorted(D11))
        if key not in seen:
            seen.add(key)
            unique_gflat.append((D11, max_A_all))
    unique_gflat = unique_gflat[:10]  # Top 10

    print(f"  Best max_A(all) values: {[mA for _, mA in unique_gflat]}")
    sys.stdout.flush()

    for i, (D11, _) in enumerate(unique_gflat):
        props = analyze_d11(D11, p)
        A = autocorrelation_single(D11, p)
        stats = sample_d12_for_d11(D11, p, A, num_d12_samples, batch_size, rng)

        print(f"  [{i}] max_A(D11)={props['max_A_d11']}, max_A(all)={props['max_A_all']}, "
              f"QR={'bal' if props['qr_balanced'] else 'no'}, "
              f"binding={stats['num_valid_binding']} (ratio={stats['ratio_binding']:.4f}), "
              f"full={stats['num_valid_full']}")

        try:
            gauss = gaussian_first_moment(A, D11, set(range(1, p)) - D11, p)
            print(f"      Gaussian log(E[valid])={gauss['log_E_valid']:.2f}")
        except ImportError:
            gauss = None
        sys.stdout.flush()

        all_results.append({
            "type": "greedy_global_aflat",
            "D11": sorted(D11),
            **props,
            **stats,
            "gaussian": gauss,
        })

    # 3. QR-balanced D11 (with A-flatness optimization)
    print(f"\n--- QR-balanced D11 ({num_qr_balanced} attempts) ---")
    sys.stdout.flush()

    for i in range(num_qr_balanced):
        D11, max_A = qr_balanced_d11(p, rng, num_attempts=50)
        props = analyze_d11(D11, p)
        A = autocorrelation_single(D11, p)
        stats = sample_d12_for_d11(D11, p, A, num_d12_samples, batch_size, rng)

        print(f"  [{i}] max_A(D11)={max_A}, max_A(all)={props['max_A_all']}, "
              f"QR={props['qr_count']}/{props['qnr_count']}, "
              f"binding={stats['num_valid_binding']}, full={stats['num_valid_full']}")

        try:
            gauss = gaussian_first_moment(A, D11, set(range(1, p)) - D11, p)
            print(f"      Gaussian log(E[valid])={gauss['log_E_valid']:.2f}")
        except ImportError:
            gauss = None
        sys.stdout.flush()

        all_results.append({
            "type": "qr_balanced",
            "D11": sorted(D11),
            **props,
            **stats,
            "gaussian": gauss,
        })

    # 4. Random D11 (baseline)
    print(f"\n--- Random D11 ({num_random} samples) ---")
    sys.stdout.flush()

    for i in range(num_random):
        D11 = random_symmetric_D11(p, rng)
        props = analyze_d11(D11, p)
        A = autocorrelation_single(D11, p)
        max_A = props["max_A_d11"]
        stats = sample_d12_for_d11(D11, p, A, num_d12_samples, batch_size, rng)

        if (i + 1) % 5 == 0:
            print(f"  [{i}] max_A(D11)={max_A}, max_A(all)={props['max_A_all']}, "
                  f"binding={stats['num_valid_binding']}, full={stats['num_valid_full']}")
            sys.stdout.flush()

        try:
            gauss = gaussian_first_moment(A, D11, set(range(1, p)) - D11, p)
        except ImportError:
            gauss = None

        all_results.append({
            "type": "random",
            "D11": sorted(D11),
            **props,
            **stats,
            "gaussian": gauss,
        })

    elapsed = time.time() - t0

    # Summary
    print(f"\n{'=' * 60}")
    print(f"SUMMARY for p = {p} ({elapsed:.1f}s)")
    print(f"{'=' * 60}")

    for dtype in ["known", "greedy_d11_only", "greedy_global_aflat", "qr_balanced", "random"]:
        subset = [r for r in all_results if r["type"] == dtype]
        if not subset:
            continue
        working_binding = sum(1 for r in subset if r["num_valid_binding"] > 0)
        working_full = sum(1 for r in subset if r["num_valid_full"] > 0)
        max_As = [r["max_A_d11"] for r in subset]
        print(f"\n  {dtype}: {len(subset)} D11")
        print(f"    max_A(D11): {min(max_As)}-{max(max_As)}")
        print(f"    Working (binding): {working_binding}/{len(subset)}")
        print(f"    Working (full): {working_full}/{len(subset)}")

        if any(r.get("gaussian") for r in subset):
            log_Es = [r["gaussian"]["log_E_valid"] for r in subset if r.get("gaussian")]
            print(f"    Gaussian log(E[valid]): {min(log_Es):.2f} to {max(log_Es):.2f}")

        # Binding ratio stats
        ratios = [r["ratio_binding"] for r in subset
                  if r["num_valid_binding"] > 0 and isinstance(r["ratio_binding"], (int, float))
                  and r["ratio_binding"] < float('inf')]
        if ratios:
            print(f"    Binding ratio: {min(ratios):.4f} to {max(ratios):.4f}")

    sys.stdout.flush()

    return {
        "p": p, "n": n,
        "elapsed": elapsed,
        "results": all_results,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--primes", nargs="+", type=int, default=[31, 43, 47, 59])
    parser.add_argument("--samples", type=int, default=1000000)
    parser.add_argument("--batch", type=int, default=10000)
    args = parser.parse_args()

    print("=" * 80)
    print("ENHANCED SAMPLING: A-FLAT AND QR-BALANCED D11")
    print("=" * 80)
    sys.stdout.flush()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.join(script_dir, "..", "enumeration_data")
    os.makedirs(output_dir, exist_ok=True)

    cross_prime = {}

    for p in args.primes:
        if p <= 31:
            ng, nqr, nr = 20, 10, 20
            bs = args.batch
        elif p <= 47:
            ng, nqr, nr = 15, 8, 15
            bs = min(args.batch, 5000)
        else:
            ng, nqr, nr = 10, 5, 10
            bs = min(args.batch, 3000)

        data = run_prime(p, num_d12_samples=args.samples, batch_size=bs,
                         num_greedy=ng, num_qr_balanced=nqr, num_random=nr)

        # Save
        outpath = os.path.join(output_dir, f"sampling_aflat_p{p}.json")
        save_results = []
        for r in data["results"]:
            r_copy = dict(r)
            # Trim verbose fields for JSON
            if r_copy.get("gaussian"):
                g = r_copy["gaussian"]
                r_copy["gaussian"] = {
                    "log_E_valid": g["log_E_valid"],
                    "E_valid_approx": g["E_valid_approx"],
                    "mu_B": g["mu_B"],
                    "sigma_B": g["sigma_B"],
                }
            if "A_values_d11" in r_copy:
                del r_copy["A_values_d11"]  # Too verbose for JSON
            save_results.append(r_copy)

        save_data = {
            "p": p, "n": data["n"],
            "elapsed": data["elapsed"],
            "num_d12_samples": args.samples,
            "results": save_results,
        }
        with open(outpath, "w") as f:
            json.dump(save_data, f, indent=2, default=str)
        print(f"\nSaved to {outpath}")

        # Cross-prime summary entry
        best_gaussian = None
        for r in data["results"]:
            if r.get("gaussian"):
                le = r["gaussian"]["log_E_valid"]
                if best_gaussian is None or le > best_gaussian:
                    best_gaussian = le
        cross_prime[p] = {
            "p": p, "n": data["n"],
            "best_gaussian_log_E": best_gaussian,
            "any_binding_hit": any(r["num_valid_binding"] > 0 for r in data["results"]),
            "any_full_hit": any(r["num_valid_full"] > 0 for r in data["results"]),
            "best_max_A_d11": min(r["max_A_d11"] for r in data["results"]),
        }

    # Cross-prime table
    print(f"\n{'=' * 80}")
    print("CROSS-PRIME SUMMARY")
    print(f"{'=' * 80}")
    print(f"{'p':>4s} {'n':>4s} {'best_maxA':>10s} {'binding_hit':>12s} {'full_hit':>10s} {'gauss_logE':>12s}")
    for p, s in cross_prime.items():
        gle = f"{s['best_gaussian_log_E']:.1f}" if s['best_gaussian_log_E'] else "N/A"
        print(f"  {p:4d} {s['n']:4d} {s['best_max_A_d11']:10d} "
              f"{'YES' if s['any_binding_hit'] else 'no':>12s} "
              f"{'YES' if s['any_full_hit'] else 'no':>10s} "
              f"{gle:>12s}")
    sys.stdout.flush()


if __name__ == "__main__":
    main()
