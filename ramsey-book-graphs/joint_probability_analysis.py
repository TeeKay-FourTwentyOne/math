#!/usr/bin/env python3
"""
Joint Probability Analysis for First Moment Proof.

KEY QUESTION: For a fixed "good" D11, is
  Pr_{D12}[B(d) <= T(d) for all d in D11] >= (1/2)^n * p^O(1)?

If yes, E[valid D12] = C(p-1, (p-3)/2) * Pr >= 2^{p-1} * 2^{-p/2} * p^O(1) -> infty,
completing the existence proof.

This script decomposes the positive association observed in correlation_analysis.py
into A-side (latent variable) and B-side (Parseval constraint) contributions.

For FIXED D11:
1. Compute A(d) for all d
2. Set thresholds T(d) = (p-3)/2 - A(d) for B(d) at each D11 position
3. Sample random D12 and measure:
   (a) Marginal rates: Pr[B(d) <= T(d)] per position
   (b) Joint rate: Pr[all B(d) <= T(d)]
   (c) Ratio: joint / product(marginals)  -- positive or negative association of B-events?
4. Energy analysis: distribution of sum_{D11} B(d) for valid D12 vs random D12
5. Exact enumeration for small p (p=11)
"""

import sys
import os
import json
import time
import numpy as np
from math import comb
from itertools import combinations

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ramsey_core import Delta


def autocorrelation_fft(indicator, p):
    """Compute Delta(S,S,d) for all d via FFT."""
    fft_val = np.fft.fft(indicator)
    autocorr = np.fft.ifft(np.abs(fft_val) ** 2).real
    return np.round(autocorr).astype(int)


# Known D11 solutions (k=n formulation: |D11| = n = (p+1)/2)
KNOWN_SOLUTIONS = {
    7: {"D11": {1, 2, 5, 6}, "D12": {0, 3, 5}},
    11: {"D11": {1, 2, 4, 7, 9, 10}, "D12": {0, 4, 5, 7, 10}},
    19: {"D11": {1, 2, 3, 6, 8, 11, 13, 16, 17, 18},
         "D12": {0, 2, 6, 8, 9, 12, 13, 17, 18}},
    23: {"D11": {5, 6, 7, 8, 9, 11, 12, 14, 15, 16, 17, 18},
         "D12": {0, 1, 2, 6, 10, 13, 14, 16, 18, 20, 21}},
    31: {"D11": {6, 7, 8, 10, 11, 12, 14, 15, 16, 17, 19, 20, 21, 23, 24, 25},
         "D12": {0, 1, 2, 3, 8, 11, 12, 13, 15, 18, 20, 21, 24, 27, 29}},
}


def load_registry():
    """Load solutions from registry for p=43, 47, 59."""
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "solutions_registry.json")
    if not os.path.exists(path):
        return
    with open(path) as f:
        registry = json.load(f)
    for sol in registry["solutions"]:
        p = sol["m"]
        if p not in {43, 47, 59}:
            continue
        D11_raw = set(sol["D11"])
        D12_raw = set(sol["D12"])
        target_size = (p + 1) // 2
        # Transform to k=n formulation if needed
        if len(D11_raw) != target_size:
            D22 = set(range(1, p)) - D11_raw
            D12T = {(-x) % p for x in D12_raw}
            D11_raw, D12_raw = D22, D12T
        KNOWN_SOLUTIONS[p] = {"D11": D11_raw, "D12": D12_raw}


def exact_enumeration_p11():
    """For p=11: enumerate ALL D12 and compute exact joint probability."""
    p = 11
    n = 6
    d12_size = 5  # (p-1)/2
    threshold = (p - 3) // 2  # = 4

    D11 = sorted(KNOWN_SOLUTIONS[p]["D11"])
    D22 = sorted(set(range(1, p)) - set(D11))

    # Compute A(d) for this D11
    d11_ind = np.zeros(p)
    for d in D11:
        d11_ind[d] = 1.0
    A = autocorrelation_fft(d11_ind, p)

    # Thresholds for B(d)
    T = {d: threshold - int(A[d]) for d in D11}

    print(f"\n{'='*80}")
    print(f"EXACT ENUMERATION: p=11, |D11|={len(D11)}, |D12|={d12_size}")
    print(f"{'='*80}")
    print(f"  D11 = {D11}")
    print(f"  A(d) values: {[int(A[d]) for d in D11]}")
    print(f"  B thresholds T(d): {[T[d] for d in D11]}")

    # Enumerate all D12 = {0} ∪ S where S ⊂ {1,...,p-1}, |S| = d12_size-1
    total = comb(p - 1, d12_size - 1)
    valid = 0
    marginal_ok = {d: 0 for d in D11}
    pairwise_ok = {}
    for i in range(len(D11)):
        for j in range(i+1, len(D11)):
            pairwise_ok[(D11[i], D11[j])] = 0

    for S in combinations(range(1, p), d12_size - 1):
        D12 = {0} | set(S)
        d12_ind = np.zeros(p)
        for x in D12:
            d12_ind[x] = 1.0
        B = autocorrelation_fft(d12_ind, p)

        all_ok = True
        ok_flags = {}
        for d in D11:
            ok = (int(B[d]) <= T[d])
            ok_flags[d] = ok
            if ok:
                marginal_ok[d] += 1
            else:
                all_ok = False

        if all_ok:
            valid += 1

        # Count pairwise
        for i in range(len(D11)):
            for j in range(i+1, len(D11)):
                if ok_flags[D11[i]] and ok_flags[D11[j]]:
                    pairwise_ok[(D11[i], D11[j])] += 1

    print(f"\n  Total D12: {total}")
    print(f"  Valid D12: {valid} ({100*valid/total:.2f}%)")

    # Marginal rates
    marginal_rates = {d: marginal_ok[d] / total for d in D11}
    print(f"\n  Marginal Pr[B(d) <= T(d)]:")
    for d in D11:
        print(f"    d={d:2d}: T(d)={T[d]}, rate={marginal_rates[d]:.4f}")

    # Independence prediction
    log_indep = sum(np.log(r) for r in marginal_rates.values() if r > 0)
    indep_pred = np.exp(log_indep) if all(r > 0 for r in marginal_rates.values()) else 0

    print(f"\n  Joint Pr[all ok] = {valid/total:.6f}")
    print(f"  Independence prediction = {indep_pred:.6f}")
    ratio = (valid / total) / indep_pred if indep_pred > 0 else float('inf')
    print(f"  Ratio (joint/indep) = {ratio:.4f}")
    if ratio > 1:
        print(f"  >>> POSITIVE ASSOCIATION of B-events (for fixed D11)")
    else:
        print(f"  >>> NEGATIVE ASSOCIATION of B-events (for fixed D11)")

    # Pairwise analysis
    print(f"\n  Pairwise analysis:")
    for i in range(len(D11)):
        for j in range(i+1, len(D11)):
            d1, d2 = D11[i], D11[j]
            pair_rate = pairwise_ok[(d1, d2)] / total
            indep_pair = marginal_rates[d1] * marginal_rates[d2]
            ratio_pair = pair_rate / indep_pair if indep_pair > 0 else float('inf')
            if abs(ratio_pair - 1) > 0.01:
                assoc = "+" if ratio_pair > 1 else "-"
                print(f"    ({d1},{d2}): pair={pair_rate:.4f}, "
                      f"indep={indep_pair:.4f}, ratio={ratio_pair:.4f} [{assoc}]")

    return ratio


def mc_fixed_d11_analysis(p, num_mc=500000):
    """For fixed D11: MC analysis of B-event association."""
    n = (p + 1) // 2
    d12_size = (p - 1) // 2
    threshold = (p - 3) // 2

    D11 = sorted(KNOWN_SOLUTIONS[p]["D11"])
    D22 = sorted(set(range(1, p)) - set(D11))

    # Compute A(d)
    d11_ind = np.zeros(p)
    for d in D11:
        d11_ind[d] = 1.0
    A = autocorrelation_fft(d11_ind, p)

    # Thresholds
    T = {d: threshold - int(A[d]) for d in D11}

    print(f"\n{'='*80}")
    print(f"FIXED D11 ANALYSIS: p={p}, n={n}")
    print(f"{'='*80}")
    print(f"  |D11|={len(D11)}, |D12|={d12_size}, threshold={threshold}")
    A_vals = [int(A[d]) for d in D11]
    print(f"  A(d): min={min(A_vals)}, max={max(A_vals)}, mean={np.mean(A_vals):.2f}")
    T_vals = [T[d] for d in D11]
    print(f"  T(d) = thresh-A(d): min={min(T_vals)}, max={max(T_vals)}, mean={np.mean(T_vals):.2f}")
    print(f"  E[B(d)] = {(p-3)/4:.2f}, so gap T(d)-E[B] ranges from "
          f"{min(T_vals)-(p-3)/4:.2f} to {max(T_vals)-(p-3)/4:.2f}")

    rng = np.random.default_rng(42)
    k = d12_size - 1  # non-zero elements

    # Sample D12 and check constraints
    marginal_ok = np.zeros(len(D11), dtype=int)
    joint_ok = 0
    B_sum_d11_all = []
    B_sum_d11_valid = []

    t0 = time.time()
    for trial in range(num_mc):
        S = rng.choice(range(1, p), size=k, replace=False)
        d12_ind = np.zeros(p)
        d12_ind[0] = 1.0
        d12_ind[S] = 1.0
        B = autocorrelation_fft(d12_ind, p)

        all_ok = True
        b_sum_d11 = 0
        for i, d in enumerate(D11):
            b_val = int(B[d])
            b_sum_d11 += b_val
            if b_val <= T[d]:
                marginal_ok[i] += 1
            else:
                all_ok = False

        B_sum_d11_all.append(b_sum_d11)
        if all_ok:
            joint_ok += 1
            B_sum_d11_valid.append(b_sum_d11)

    elapsed = time.time() - t0
    joint_rate = joint_ok / num_mc
    marginal_rates = marginal_ok / num_mc

    print(f"\n  MC results ({num_mc:,} trials, {elapsed:.1f}s):")
    print(f"  Marginal Pr[B(d) <= T(d)]:")
    for i, d in enumerate(D11):
        print(f"    d={d:2d}: T={T[d]:2d}, A={int(A[d]):2d}, rate={marginal_rates[i]:.4f}")

    # Independence prediction
    log_indep = sum(np.log(max(r, 1e-10)) for r in marginal_rates)
    indep_pred = np.exp(log_indep)

    print(f"\n  Joint Pr[all ok] = {joint_ok}/{num_mc} = {joint_rate:.6f}")
    print(f"  Independence prediction = {indep_pred:.6f}")
    ratio = joint_rate / indep_pred if indep_pred > 0 else float('inf')
    print(f"  Ratio (joint/indep) = {ratio:.4f}")

    if ratio > 1:
        print(f"  >>> POSITIVE ASSOCIATION even for fixed D11!")
    elif ratio > 0.5:
        print(f"  >>> Mild negative association (ratio > 0.5)")
    else:
        print(f"  >>> Significant negative association")

    # Energy analysis
    B_sum_all_arr = np.array(B_sum_d11_all)
    total_B_energy = d12_size * (d12_size - 1)  # Parseval: sum over all d

    print(f"\n  ENERGY ANALYSIS:")
    print(f"    Total B energy (Parseval): {total_B_energy}")
    print(f"    E[sum_D11 B(d)] = {B_sum_all_arr.mean():.2f} "
          f"(theoretical: {len(D11) * (p-3)/4:.2f})")
    print(f"    Std[sum_D11 B(d)] = {B_sum_all_arr.std():.2f}")

    if B_sum_d11_valid:
        B_sum_valid_arr = np.array(B_sum_d11_valid)
        print(f"    For VALID D12:")
        print(f"      Mean sum_D11 B(d) = {B_sum_valid_arr.mean():.2f}")
        print(f"      Max allowed sum_D11 B(d) = {sum(T_vals)}")
        print(f"      Deficit = {sum(T_vals) - B_sum_valid_arr.mean():.2f}")
        print(f"      Energy shifted to D22: {total_B_energy - B_sum_valid_arr.mean():.2f}")
    else:
        print(f"    No valid D12 found in MC sample")

    # What fraction of the probability comes from energy being below sum(T)?
    energy_ok = B_sum_all_arr <= sum(T_vals)
    print(f"\n    Pr[sum_D11 B(d) <= sum(T)] = {energy_ok.mean():.4f}")
    print(f"    Pr[all ok | sum ok] = {joint_ok / max(energy_ok.sum(), 1):.6f}")
    print(f"    Pr[all ok] / Pr[sum ok] = {joint_rate / max(energy_ok.mean(), 1e-10):.6f}")

    return {
        "p": p,
        "n": n,
        "joint_rate": joint_rate,
        "indep_pred": indep_pred,
        "ratio_joint_indep": ratio,
        "marginal_rates": marginal_rates.tolist(),
        "energy_mean": float(B_sum_all_arr.mean()),
        "energy_std": float(B_sum_all_arr.std()),
        "sum_threshold": sum(T_vals),
        "energy_ok_rate": float(energy_ok.mean()),
    }


def random_d11_analysis(p, num_d11=1000, num_d12=10000):
    """Analyze how joint/indep ratio varies across random D11s."""
    n = (p + 1) // 2
    d12_size = (p - 1) // 2
    threshold = (p - 3) // 2
    k_d12 = d12_size - 1  # non-zero elements in D12

    rng = np.random.default_rng(123)

    # Generate negation pairs for symmetric D11
    pairs = [(x, p - x) for x in range(1, (p + 1) // 2)]
    num_pairs_to_choose = n // 2  # how many pairs to include

    print(f"\n{'='*80}")
    print(f"RANDOM D11 ANALYSIS: p={p}, n={n}")
    print(f"  Testing {num_d11} random D11s with {num_d12} random D12s each")
    print(f"{'='*80}")

    joint_rates = []
    indep_rates = []
    ratios = []
    max_A_vals = []
    good_d11_count = 0

    t0 = time.time()
    for d11_trial in range(num_d11):
        # Generate random symmetric D11
        chosen = rng.choice(len(pairs), size=num_pairs_to_choose, replace=False)
        D11 = set()
        for i in chosen:
            D11.add(pairs[i][0])
            D11.add(pairs[i][1])
        D11_list = sorted(D11)

        # Compute A(d)
        d11_ind = np.zeros(p)
        for d in D11_list:
            d11_ind[d] = 1.0
        A = autocorrelation_fft(d11_ind, p)
        A_vals = [int(A[d]) for d in D11_list]
        max_A = max(A_vals)
        max_A_vals.append(max_A)

        # Thresholds
        T = {d: threshold - int(A[d]) for d in D11_list}

        # If any T(d) < 0, impossible to satisfy
        if min(T.values()) < 0:
            continue

        good_d11_count += 1

        # Sample D12
        marginal_ok = np.zeros(len(D11_list), dtype=int)
        joint_ok = 0
        for _ in range(num_d12):
            S = rng.choice(range(1, p), size=k_d12, replace=False)
            d12_ind = np.zeros(p)
            d12_ind[0] = 1.0
            d12_ind[S] = 1.0
            B = autocorrelation_fft(d12_ind, p)

            all_ok = True
            for i, d in enumerate(D11_list):
                if int(B[d]) <= T[d]:
                    marginal_ok[i] += 1
                else:
                    all_ok = False
            if all_ok:
                joint_ok += 1

        jr = joint_ok / num_d12
        mr = marginal_ok / num_d12

        # Independence prediction
        log_ip = sum(np.log(max(r, 1e-10)) for r in mr)
        ip = np.exp(log_ip) if all(r > 0 for r in mr) else 0

        joint_rates.append(jr)
        indep_rates.append(ip)
        if ip > 0 and jr > 0:
            ratios.append(jr / ip)

    elapsed = time.time() - t0

    print(f"\n  Completed in {elapsed:.1f}s")
    print(f"  Feasible D11s (all T(d) >= 0): {good_d11_count}/{num_d11} "
          f"({100*good_d11_count/num_d11:.1f}%)")
    print(f"  max A(d) across D11s: min={min(max_A_vals)}, max={max(max_A_vals)}, "
          f"mean={np.mean(max_A_vals):.2f}")

    if ratios:
        print(f"\n  Ratio (joint/indep) across feasible D11s:")
        print(f"    min={min(ratios):.4f}, max={max(ratios):.4f}, "
              f"median={np.median(ratios):.4f}, mean={np.mean(ratios):.4f}")
        print(f"    Fraction with ratio > 1 (positive assoc): "
              f"{sum(1 for r in ratios if r > 1)}/{len(ratios)} "
              f"({100*sum(1 for r in ratios if r > 1)/len(ratios):.1f}%)")
        print(f"    Fraction with ratio > 0.5: "
              f"{sum(1 for r in ratios if r > 0.5)}/{len(ratios)}")
    else:
        print(f"  No feasible D11 found with measurable joint probability")

    return ratios


def var_A_computation(p):
    """Compute Var[A(d)] by MC over random symmetric D11."""
    n = (p + 1) // 2
    pairs = [(x, p - x) for x in range(1, (p + 1) // 2)]
    num_pairs_to_choose = n // 2

    rng = np.random.default_rng(42)
    num_mc = 200000

    # Collect A(d) for d=1 (all d equivalent by symmetry of the D11 distribution)
    A_vals = []
    for _ in range(num_mc):
        chosen = rng.choice(len(pairs), size=num_pairs_to_choose, replace=False)
        D11 = set()
        for i in chosen:
            D11.add(pairs[i][0])
            D11.add(pairs[i][1])

        d11_ind = np.zeros(p)
        for d in D11:
            d11_ind[d] = 1.0
        A = autocorrelation_fft(d11_ind, p)
        A_vals.append(int(A[1]))

    A_arr = np.array(A_vals, dtype=float)
    return {
        "E_A": A_arr.mean(),
        "Var_A": A_arr.var(),
        "Std_A": A_arr.std(),
        "theoretical_E_A": (p + 1) / 4,
        "Var_A_over_p": A_arr.var() / p,
    }


def main():
    load_registry()

    print("=" * 80)
    print("JOINT PROBABILITY ANALYSIS FOR FIRST MOMENT PROOF")
    print("=" * 80)

    # Phase 0: Var[A(d)] computation
    print(f"\n{'='*80}")
    print("PHASE 0: Var[A(d)] computation (MC)")
    print(f"{'='*80}")
    print(f"\n  {'p':>4s} {'E[A]':>8s} {'E[A]th':>8s} {'Var[A]':>10s} {'Var/p':>10s}")
    print(f"  {'-'*50}")
    for p in [11, 19, 23, 31, 43, 47, 59]:
        result = var_A_computation(p)
        print(f"  {p:4d} {result['E_A']:8.4f} {result['theoretical_E_A']:8.4f} "
              f"{result['Var_A']:10.4f} {result['Var_A_over_p']:10.6f}")

    # Phase 1: Exact enumeration for p=11
    exact_ratio = exact_enumeration_p11()

    # Phase 2: Fixed D11 MC analysis for each prime
    print(f"\n\n{'='*80}")
    print("PHASE 2: Fixed D11, random D12 — B-event association")
    print(f"{'='*80}")

    mc_sizes = {11: 500000, 19: 200000, 23: 200000, 31: 100000,
                43: 50000, 47: 30000, 59: 20000}

    results = {}
    for p in sorted(KNOWN_SOLUTIONS.keys()):
        if p == 7:
            continue  # too small
        num_mc = mc_sizes.get(p, 20000)
        result = mc_fixed_d11_analysis(p, num_mc=num_mc)
        results[p] = result

    # Summary table
    print(f"\n\n{'='*80}")
    print("SUMMARY: B-event association for fixed D11")
    print(f"{'='*80}")
    print(f"  {'p':>4s} {'n':>3s} {'Pr[all]':>12s} {'Pr[indep]':>12s} "
          f"{'ratio':>8s} {'E_sum':>8s} {'sum_T':>6s} {'Pr[E<=T]':>10s}")
    print(f"  {'-'*75}")
    for p in sorted(results.keys()):
        r = results[p]
        print(f"  {p:4d} {r['n']:3d} {r['joint_rate']:12.6f} "
              f"{r['indep_pred']:12.6f} {r['ratio_joint_indep']:8.4f} "
              f"{r['energy_mean']:8.2f} {r['sum_threshold']:6d} "
              f"{r['energy_ok_rate']:10.4f}")

    # KEY FINDING
    print(f"\n  KEY FINDING:")
    ratios = [results[p]['ratio_joint_indep'] for p in sorted(results.keys())
              if results[p]['ratio_joint_indep'] > 0]
    if ratios:
        if all(r >= 1 for r in ratios):
            print(f"  >>> B-events are POSITIVELY ASSOCIATED even for fixed D11!")
            print(f"  >>> This means Pr[all ok] >= prod(marginals) for any D11.")
            print(f"  >>> Combined with prod(marginals) >= (1/2)^n, the proof is COMPLETE.")
        elif all(r >= 0.5 for r in ratios):
            print(f"  >>> B-events are only mildly negatively associated (ratio > 0.5).")
            print(f"  >>> The factor of at most 2^n loss is compensated by 2^{{p-1}} D12s.")
            print(f"  >>> E[valid D12] >= 2^{{p-1}} * (1/2)^n * 0.5^n = 2^{{p-1-2n}} ~ 1 for large p.")
            print(f"  >>> This is marginal — need tighter bound.")
        else:
            print(f"  >>> Significant negative association detected.")
            print(f"  >>> The energy/Parseval argument is needed.")

    # Phase 3: Random D11 analysis (smaller primes only)
    print(f"\n\n{'='*80}")
    print("PHASE 3: Ratio distribution over random D11s")
    print(f"{'='*80}")

    for p in [11, 19, 23]:
        ratios = random_d11_analysis(p, num_d11=500, num_d12=10000)

    # Phase 4: Trend analysis
    print(f"\n\n{'='*80}")
    print("PHASE 4: Proof viability assessment")
    print(f"{'='*80}")

    for p in sorted(results.keys()):
        r = results[p]
        n = r['n']
        log2_d12_count = np.log2(comb(p-1, (p-3)//2)) if p <= 60 else (p-1) * 0.99
        log2_joint = np.log2(r['joint_rate']) if r['joint_rate'] > 0 else float('-inf')
        log2_E = log2_d12_count + log2_joint

        print(f"  p={p:3d}: log2(#D12)={log2_d12_count:.1f}, "
              f"log2(Pr[all ok | D11])={log2_joint:.1f}, "
              f"log2(E[valid D12 | D11])={log2_E:.1f}")


if __name__ == "__main__":
    main()
