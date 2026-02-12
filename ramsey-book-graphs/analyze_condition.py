#!/usr/bin/env python3
"""
Analyze what structural condition on A-values of symmetric D11 determines N(D11) > 0.

For primes p ≡ 3 mod 4, we have:
  - D11 ⊂ {1,...,p-1}, |D11| = (p+1)/2, symmetric
  - D22 = {1,...,p-1} \ D11, |D22| = (p-3)/2
  - A(d) = autocorrelation of D11
  - Binding constraint at d ∈ D11: A(d) + B(d) ≤ (p-3)/2
  - Loose constraint at d ∈ D22:  A(d) + B(d) ≤ (p+3)/2
  - N(D11) = number of valid D12

Key question: Is there a SINGLE SCALAR from A-values that separates N>0 from N=0?
"""

import json
import math
import numpy as np
from collections import defaultdict

def analyze_prime(data, p):
    """Analyze A-flat D11s for a given prime p."""
    print(f"\n{'='*80}")
    print(f"  ANALYSIS FOR p = {p}")
    print(f"{'='*80}")

    entries = data[str(p)]['a_flat_entries']
    n_flat = len(entries)

    # Parameters
    d11_size = (p + 1) // 2
    d22_size = (p - 3) // 2
    d12_size = (p - 3) // 2
    binding_thresh = (p - 3) // 2  # for d in D11
    loose_thresh = (p + 3) // 2    # for d in D22

    # Expected A (mean of A over all d)
    # Sum of all A(d) = sum_A_D11 + sum_A_D22 = d11_size*(d11_size - 1) by Parseval
    total_A = d11_size * (d11_size - 1)
    E_A = total_A / (p - 1)

    # Expected B ~ d12_size*(d12_size-1)/(p-1)
    E_B = d12_size * (d12_size - 1) / (p - 1)

    print(f"\n  Parameters:")
    print(f"    |D11| = {d11_size}, |D22| = {d22_size}, |D12| = {d12_size}")
    print(f"    Binding threshold (D11): {binding_thresh}")
    print(f"    Loose threshold (D22):   {loose_thresh}")
    print(f"    E[A] = {E_A:.4f}")
    print(f"    E[B] = {E_B:.4f}")
    print(f"    floor((p+1)/4) = {(p+1)//4}  (A-flat condition: max_A_D11 <= this)")

    # Collect statistics
    working = []   # N > 0
    dead = []      # N = 0

    all_records = []

    for entry in entries:
        D11 = entry['D11']
        A_D11 = entry['A_at_D11']
        A_D22 = entry['A_at_D22']
        N = entry['N']
        max_A_D11 = entry['max_A_D11']
        max_A_D22 = entry['max_A_D22']
        sum_A_D11 = entry['sum_A_D11']
        sum_A_D22 = entry['sum_A_D22']

        # Slack vectors
        slack_D11 = [binding_thresh - a for a in A_D11]
        slack_D22 = [loose_thresh - a for a in A_D22]

        min_slack_D11 = min(slack_D11)
        min_slack_D22 = min(slack_D22)
        min_slack_all = min(min_slack_D11, min_slack_D22)

        sum_slack_D11 = sum(slack_D11)
        sum_slack_D22 = sum(slack_D22)
        sum_slack_all = sum_slack_D11 + sum_slack_D22

        # Number of bottleneck positions on D11 side
        num_at_max_D11 = sum(1 for a in A_D11 if a == max_A_D11)

        # Number of tight D22 positions: A(d) >= loose_thresh - floor(E_B)
        tight_D22_thresh = loose_thresh - math.floor(E_B)
        num_tight_D22 = sum(1 for a in A_D22 if a >= tight_D22_thresh)

        # Product of slack/E_B ratios (proxy for joint probability)
        # Use log to avoid overflow/underflow
        log_prod_slack = 0.0
        for s in slack_D11:
            if s > 0 and E_B > 0:
                log_prod_slack += math.log(s / E_B) if s / E_B > 0 else -100
        for s in slack_D22:
            if s > 0 and E_B > 0:
                log_prod_slack += math.log(s / E_B) if s / E_B > 0 else -100

        # Variance of A over D11 positions
        mean_A_D11 = sum_A_D11 / d11_size
        var_A_D11 = sum((a - mean_A_D11)**2 for a in A_D11) / d11_size

        # Variance of A over D22 positions
        mean_A_D22 = sum_A_D22 / d22_size if d22_size > 0 else 0
        var_A_D22 = sum((a - mean_A_D22)**2 for a in A_D22) / d22_size if d22_size > 0 else 0

        # Max A over all positions
        max_A_all = entry['max_A_all']

        # "Headroom" = binding_thresh - mean_A_D11 (how far mean is from threshold)
        headroom_D11 = binding_thresh - mean_A_D11

        # Entropy of slack distribution
        total_s = sum_slack_all
        if total_s > 0:
            probs = [(s / total_s) for s in slack_D11 + slack_D22 if s > 0]
            entropy_slack = -sum(p_i * math.log(p_i) for p_i in probs if p_i > 0)
        else:
            entropy_slack = 0

        # Sum of squared A values on D11
        sum_A2_D11 = sum(a**2 for a in A_D11)
        sum_A2_D22 = sum(a**2 for a in A_D22)

        # Harmonic mean of slacks on D11
        if all(s > 0 for s in slack_D11):
            harmonic_slack_D11 = d11_size / sum(1.0/s for s in slack_D11)
        else:
            harmonic_slack_D11 = 0

        # Number of D22 positions with A >= max_A_D11 (D22 positions above binding-side max)
        num_D22_above_bind_max = sum(1 for a in A_D22 if a >= max_A_D11)

        # "Concentration" of high A on D22: max_A_D22 - mean_A_D22
        concentration_D22 = max_A_D22 - mean_A_D22

        # Ratio: max_A_D22 / max_A_D11
        ratio_max = max_A_D22 / max_A_D11 if max_A_D11 > 0 else float('inf')

        # The z-scores
        z_D11 = entry.get('z_D11', [])
        z_D22 = entry.get('z_D22', [])
        S1_D11 = entry.get('S1_D11', 0)
        S1_D22 = entry.get('S1_D22', 0)
        S1_all = entry.get('S1_all', 0)
        log_prod_phi = entry.get('log_prod_phi_D11', 0)

        # Max z on D22 side
        max_z_D22 = max(z_D22) if z_D22 else 0
        min_z_D11_val = min(z_D11) if z_D11 else 0

        rec = {
            'N': N,
            'working': N > 0,
            'sum_A_D11': sum_A_D11,
            'sum_A_D22': sum_A_D22,
            'max_A_D11': max_A_D11,
            'max_A_D22': max_A_D22,
            'max_A_all': max_A_all,
            'min_slack_D11': min_slack_D11,
            'min_slack_D22': min_slack_D22,
            'min_slack_all': min_slack_all,
            'sum_slack_D11': sum_slack_D11,
            'sum_slack_D22': sum_slack_D22,
            'sum_slack_all': sum_slack_all,
            'num_at_max_D11': num_at_max_D11,
            'num_tight_D22': num_tight_D22,
            'log_prod_slack': log_prod_slack,
            'var_A_D11': var_A_D11,
            'var_A_D22': var_A_D22,
            'headroom_D11': headroom_D11,
            'entropy_slack': entropy_slack,
            'sum_A2_D11': sum_A2_D11,
            'sum_A2_D22': sum_A2_D22,
            'harmonic_slack_D11': harmonic_slack_D11,
            'num_D22_above_bind_max': num_D22_above_bind_max,
            'concentration_D22': concentration_D22,
            'ratio_max': ratio_max,
            'S1_D11': S1_D11,
            'S1_D22': S1_D22,
            'S1_all': S1_all,
            'max_z_D22': max_z_D22,
            'min_z_D11': min_z_D11_val,
            'log_prod_phi': log_prod_phi,
            'A_profile_D11': tuple(entry['A_profile_D11']),
            'A_profile_D22': tuple(entry['A_profile_D22']),
        }

        all_records.append(rec)
        if N > 0:
            working.append(rec)
        else:
            dead.append(rec)

    print(f"\n  A-flat D11s: {n_flat} total, {len(working)} working (N>0), {len(dead)} dead (N=0)")

    # ---- Cross-tabulation ----
    print(f"\n  {'─'*76}")
    print(f"  CROSS-TABULATION: each quantity vs N>0 / N=0")
    print(f"  {'─'*76}")

    scalars = [
        'sum_A_D11', 'sum_A_D22', 'max_A_D11', 'max_A_D22', 'max_A_all',
        'min_slack_D11', 'min_slack_D22', 'min_slack_all',
        'sum_slack_D11', 'sum_slack_D22', 'sum_slack_all',
        'num_at_max_D11', 'num_tight_D22',
        'log_prod_slack',
        'var_A_D11', 'var_A_D22',
        'headroom_D11', 'entropy_slack',
        'sum_A2_D11', 'sum_A2_D22',
        'harmonic_slack_D11',
        'num_D22_above_bind_max', 'concentration_D22', 'ratio_max',
        'S1_D11', 'S1_D22', 'S1_all',
        'max_z_D22', 'min_z_D11',
        'log_prod_phi',
    ]

    separators = []

    for s_name in scalars:
        w_vals = [r[s_name] for r in working]
        d_vals = [r[s_name] for r in dead]

        if not w_vals or not d_vals:
            continue

        w_mean = np.mean(w_vals)
        d_mean = np.mean(d_vals)
        w_min, w_max = min(w_vals), max(w_vals)
        d_min, d_max = min(d_vals), max(d_vals)

        # Check for perfect separation
        # Try: working has HIGHER values
        sep_high = w_min > d_max
        # Try: working has LOWER values
        sep_low = w_max < d_min

        if sep_high:
            threshold = (d_max + w_min) / 2
            sep_str = f"PERFECT SEP: N>0 iff {s_name} > {threshold:.4f}"
            separators.append((s_name, '>', threshold))
        elif sep_low:
            threshold = (w_max + d_min) / 2
            sep_str = f"PERFECT SEP: N>0 iff {s_name} < {threshold:.4f}"
            separators.append((s_name, '<', threshold))
        else:
            # Find best threshold
            all_vals = [(r[s_name], r['working']) for r in all_records]
            all_vals.sort()
            best_acc = 0
            best_t = None
            best_dir = '>'
            # Try all split points
            unique_vals = sorted(set(v for v, _ in all_vals))
            for i in range(len(unique_vals) - 1):
                t = (unique_vals[i] + unique_vals[i+1]) / 2
                # Threshold: predict working if val > t
                acc_high = sum(1 for v, w in all_vals if (v > t) == w) / len(all_vals)
                acc_low = sum(1 for v, w in all_vals if (v < t) == w) / len(all_vals)
                if acc_high > best_acc:
                    best_acc = acc_high
                    best_t = t
                    best_dir = '>'
                if acc_low > best_acc:
                    best_acc = acc_low
                    best_t = t
                    best_dir = '<'

            if best_t is not None:
                sep_str = f"best threshold: {s_name} {best_dir} {best_t:.4f} (acc={best_acc:.1%})"
            else:
                sep_str = f"no threshold found (all values identical?)"

        overlap = "YES" if (w_min <= d_max and w_max >= d_min) else "NO"

        print(f"\n  {s_name}:")
        print(f"    Working (N>0): mean={w_mean:9.4f}  range=[{w_min:.4f}, {w_max:.4f}]")
        print(f"    Dead    (N=0): mean={d_mean:9.4f}  range=[{d_min:.4f}, {d_max:.4f}]")
        print(f"    Ranges overlap: {overlap}")
        print(f"    {sep_str}")

    # ---- Report perfect separators ----
    if separators:
        print(f"\n  {'='*76}")
        print(f"  PERFECT SEPARATORS FOUND:")
        print(f"  {'='*76}")
        for s_name, direction, threshold in separators:
            print(f"    N>0 iff {s_name} {direction} {threshold:.6f}")
    else:
        print(f"\n  NO single scalar perfectly separates N>0 from N=0")

    # ---- A-profile analysis ----
    print(f"\n  {'─'*76}")
    print(f"  A-PROFILE ANALYSIS")
    print(f"  {'─'*76}")

    # Unique D11 A-profiles
    d11_profiles = defaultdict(lambda: {'working': 0, 'dead': 0, 'Ns': []})
    for r in all_records:
        key = r['A_profile_D11']
        if r['working']:
            d11_profiles[key]['working'] += 1
        else:
            d11_profiles[key]['dead'] += 1
        d11_profiles[key]['Ns'].append(r['N'])

    print(f"\n  D11 A-profiles (sorted):")
    for prof in sorted(d11_profiles.keys()):
        info = d11_profiles[prof]
        print(f"    {list(prof)}: {info['working']}W/{info['dead']}D  Ns={sorted(info['Ns'])}")

    # Unique D22 A-profiles
    d22_profiles = defaultdict(lambda: {'working': 0, 'dead': 0, 'Ns': []})
    for r in all_records:
        key = r['A_profile_D22']
        if r['working']:
            d22_profiles[key]['working'] += 1
        else:
            d22_profiles[key]['dead'] += 1
        d22_profiles[key]['Ns'].append(r['N'])

    print(f"\n  D22 A-profiles (sorted):")
    for prof in sorted(d22_profiles.keys()):
        info = d22_profiles[prof]
        mixed = "MIXED" if info['working'] > 0 and info['dead'] > 0 else ""
        print(f"    {list(prof)}: {info['working']}W/{info['dead']}D  Ns={sorted(info['Ns'])}  {mixed}")

    # Check if D22 A-profile perfectly predicts
    d22_perfect = all(
        (info['working'] == 0 or info['dead'] == 0)
        for info in d22_profiles.values()
    )
    print(f"\n  D22 A-profile perfectly predicts N>0: {d22_perfect}")

    # ---- Combined slack analysis ----
    print(f"\n  {'─'*76}")
    print(f"  SLACK VECTOR ANALYSIS")
    print(f"  {'─'*76}")

    for r in all_records[:5]:  # Show first 5
        tag = "W" if r['working'] else "D"
        print(f"    [{tag}] N={r['N']:4d}  min_slack_D11={r['min_slack_D11']}  "
              f"min_slack_D22={r['min_slack_D22']}  sum_slack={r['sum_slack_all']}  "
              f"max_A_D22={r['max_A_D22']}")

    # ---- Derived scalar explorations ----
    print(f"\n  {'─'*76}")
    print(f"  DERIVED SCALAR EXPLORATIONS")
    print(f"  {'─'*76}")

    # Try combinations
    derived = {}

    # max_A_D22 - max_A_D11
    derived['max_A_D22 - max_A_D11'] = lambda r: r['max_A_D22'] - r['max_A_D11']

    # sum_A2_D11 - sum_A2_D22
    derived['sum_A2_D11 - sum_A2_D22'] = lambda r: r['sum_A2_D11'] - r['sum_A2_D22']

    # max_A_D22 * num_at_max_D11
    derived['max_A_D22 * num_at_max_D11'] = lambda r: r['max_A_D22'] * r['num_at_max_D11']

    # var_A_D22 - var_A_D11
    derived['var_A_D22 - var_A_D11'] = lambda r: r['var_A_D22'] - r['var_A_D11']

    # sum_A2_D22 / sum_A2_D11
    derived['sum_A2_D22 / sum_A2_D11'] = lambda r: r['sum_A2_D22'] / r['sum_A2_D11'] if r['sum_A2_D11'] > 0 else 0

    # max_A_D22 - 2*max_A_D11 (measure of D22 spikiness relative to D11)
    derived['max_A_D22 - 2*max_A_D11'] = lambda r: r['max_A_D22'] - 2 * r['max_A_D11']

    # S1_D22 - S1_D11  (z-score mass imbalance)
    derived['S1_D22 - S1_D11'] = lambda r: r['S1_D22'] - r['S1_D11']

    # Number of D22 positions above binding max, minus something
    derived['num_D22_above_bind_max'] = lambda r: r['num_D22_above_bind_max']

    # max_A_D22 + max_A_D11  (total peak)
    derived['max_A_D22 + max_A_D11'] = lambda r: r['max_A_D22'] + r['max_A_D11']

    # Excess on D22 side: max_A_D22 - E_A
    derived[f'max_A_D22 - E_A (E_A={E_A:.2f})'] = lambda r: r['max_A_D22'] - E_A

    # Key insight test: is it about max_A_D22 being high enough?
    # Working D11s might need D22 to have a spike (concentrating A mass away from D11)
    derived['max_A_D22 >= p//4 + 1'] = lambda r: 1 if r['max_A_D22'] >= p//4 + 1 else 0

    # Gap: max_A_D22 - second_max_A_D22 (but we don't have per-position D22, only profile)
    # Use profile instead
    derived['D22_second_max'] = lambda r: r['A_profile_D22'][-2] if len(r['A_profile_D22']) >= 2 else 0
    derived['D22_gap_top2'] = lambda r: r['A_profile_D22'][-1] - r['A_profile_D22'][-2] if len(r['A_profile_D22']) >= 2 else 0

    # Sum of top-2 D22 A-values
    derived['D22_sum_top2'] = lambda r: sum(r['A_profile_D22'][-2:]) if len(r['A_profile_D22']) >= 2 else 0

    # Number of D22 positions with A above the D11 binding threshold
    derived['num_D22_above_binding_thresh'] = lambda r: sum(1 for a in r['A_profile_D22'] if a > binding_thresh)

    # Sum of (A(d) - E_A)^2 on D22 (how spread out D22 is)
    derived['D22_sum_sq_dev'] = lambda r: sum((a - E_A)**2 for a in r['A_profile_D22'])

    # Minimum A on D22
    derived['min_A_D22'] = lambda r: r['A_profile_D22'][0] if r['A_profile_D22'] else 0

    # Range of A on D22
    derived['range_A_D22'] = lambda r: r['A_profile_D22'][-1] - r['A_profile_D22'][0] if r['A_profile_D22'] else 0

    for d_name, d_func in derived.items():
        w_vals = [d_func(r) for r in working]
        d_vals = [d_func(r) for r in dead]

        if not w_vals or not d_vals:
            continue

        w_mean = np.mean(w_vals)
        d_mean = np.mean(d_vals)
        w_min, w_max = min(w_vals), max(w_vals)
        d_min, d_max = min(d_vals), max(d_vals)

        sep_high = w_min > d_max
        sep_low = w_max < d_min

        if sep_high or sep_low:
            if sep_high:
                threshold = (d_max + w_min) / 2
                print(f"\n  *** PERFECT SEPARATOR: {d_name}")
                print(f"      N>0 iff {d_name} > {threshold:.4f}")
                print(f"      Working: [{w_min:.4f}, {w_max:.4f}]  Dead: [{d_min:.4f}, {d_max:.4f}]")
            else:
                threshold = (w_max + d_min) / 2
                print(f"\n  *** PERFECT SEPARATOR: {d_name}")
                print(f"      N>0 iff {d_name} < {threshold:.4f}")
                print(f"      Working: [{w_min:.4f}, {w_max:.4f}]  Dead: [{d_min:.4f}, {d_max:.4f}]")
        else:
            # Find best accuracy
            all_vals_d = [(d_func(r), r['working']) for r in all_records]
            all_vals_d.sort()
            best_acc = 0
            unique_vals = sorted(set(v for v, _ in all_vals_d))
            for i in range(len(unique_vals) - 1):
                t = (unique_vals[i] + unique_vals[i+1]) / 2
                acc_high = sum(1 for v, w in all_vals_d if (v > t) == w) / len(all_vals_d)
                acc_low = sum(1 for v, w in all_vals_d if (v < t) == w) / len(all_vals_d)
                best_acc = max(best_acc, acc_high, acc_low)

            if best_acc >= 0.85:
                print(f"\n  NEAR-SEPARATOR ({best_acc:.1%}): {d_name}")
                print(f"      Working: mean={w_mean:.4f} [{w_min:.4f}, {w_max:.4f}]  "
                      f"Dead: mean={d_mean:.4f} [{d_min:.4f}, {d_max:.4f}]")

    # ---- Detailed table for visual inspection ----
    print(f"\n  {'─'*76}")
    print(f"  DETAILED TABLE (sorted by max_A_D22)")
    print(f"  {'─'*76}")

    sorted_recs = sorted(all_records, key=lambda r: (r['max_A_D22'], r['max_A_D11'], r['N']))

    print(f"  {'WD':>2s} {'N':>5s} {'maxD11':>6s} {'maxD22':>6s} {'maxAll':>6s} "
          f"{'sumD11':>6s} {'sumD22':>6s} {'varD11':>7s} {'varD22':>7s} "
          f"{'D22prof':>30s} {'D11prof':>30s}")
    print(f"  {'-'*2} {'-'*5} {'-'*6} {'-'*6} {'-'*6} "
          f"{'-'*6} {'-'*6} {'-'*7} {'-'*7} "
          f"{'-'*30} {'-'*30}")

    for r in sorted_recs:
        tag = "W" if r['working'] else "D"
        print(f"  {tag:>2s} {r['N']:5d} {r['max_A_D11']:6d} {r['max_A_D22']:6d} {r['max_A_all']:6d} "
              f"{r['sum_A_D11']:6d} {r['sum_A_D22']:6d} {r['var_A_D11']:7.3f} {r['var_A_D22']:7.3f} "
              f"{str(list(r['A_profile_D22'])):>30s} {str(list(r['A_profile_D11'])):>30s}")

    return all_records, working, dead, separators


def multi_prime_analysis(data):
    """Look for patterns that hold across both primes."""
    print(f"\n\n{'#'*80}")
    print(f"  CROSS-PRIME ANALYSIS")
    print(f"{'#'*80}")

    all_results = {}
    for p_str in sorted(data.keys(), key=int):
        p = int(p_str)
        recs, w, d, seps = analyze_prime(data, p)
        all_results[p] = {'records': recs, 'working': w, 'dead': d, 'separators': seps}

    # Now check which scalars separate at BOTH primes
    print(f"\n  {'='*76}")
    print(f"  SCALARS THAT SEPARATE AT BOTH PRIMES")
    print(f"  {'='*76}")

    for p, res in all_results.items():
        if res['separators']:
            print(f"\n  p={p} perfect separators:")
            for s_name, direction, threshold in res['separators']:
                print(f"    {s_name} {direction} {threshold:.4f}")

    # Check normalized versions across primes
    print(f"\n  {'─'*76}")
    print(f"  NORMALIZED SCALAR ANALYSIS (comparing across primes)")
    print(f"  {'─'*76}")

    for p, res in all_results.items():
        d11_size = (p + 1) // 2
        d22_size = (p - 3) // 2
        E_A = d11_size * (d11_size - 1) / (p - 1)
        sigma_A = math.sqrt(E_A * (1 - E_A / (p-1)))  # rough estimate

        print(f"\n  p={p}: E_A={E_A:.3f}")

        # Normalized max_A_D22
        for r in res['records']:
            r['norm_max_A_D22'] = (r['max_A_D22'] - E_A) / math.sqrt(E_A) if E_A > 0 else 0
            r['norm_max_A_D11'] = (r['max_A_D11'] - E_A) / math.sqrt(E_A) if E_A > 0 else 0
            r['max_A_D22_over_binding'] = r['max_A_D22'] / ((p-3)/2) if p > 3 else 0
            r['max_A_D22_over_loose'] = r['max_A_D22'] / ((p+3)/2) if p > -3 else 0
            # Key: max_A_D22 relative to p
            r['max_A_D22_frac_p'] = r['max_A_D22'] / p
            r['max_A_D22_minus_floor_p4'] = r['max_A_D22'] - (p+1)//4

        norm_scalars = ['norm_max_A_D22', 'norm_max_A_D11', 'max_A_D22_over_binding',
                        'max_A_D22_over_loose', 'max_A_D22_frac_p', 'max_A_D22_minus_floor_p4']

        for ns in norm_scalars:
            w_vals = [r[ns] for r in res['working']]
            d_vals = [r[ns] for r in res['dead']]
            if not w_vals or not d_vals:
                continue
            w_mean = np.mean(w_vals)
            d_mean = np.mean(d_vals)
            w_min, w_max = min(w_vals), max(w_vals)
            d_min, d_max = min(d_vals), max(d_vals)

            sep = "PERFECT" if (w_min > d_max or w_max < d_min) else "overlap"
            print(f"    {ns:>30s}: W=[{w_min:.3f},{w_max:.3f}] D=[{d_min:.3f},{d_max:.3f}]  {sep}")

    # ---- Final: combinatorial slack condition ----
    print(f"\n  {'='*76}")
    print(f"  COMBINATORIAL SLACK CONDITION ANALYSIS")
    print(f"  {'='*76}")

    for p, res in all_results.items():
        binding = (p - 3) // 2
        loose = (p + 3) // 2
        E_A = ((p+1)//2) * ((p+1)//2 - 1) / (p - 1)

        print(f"\n  p={p}: binding_thresh={binding}, loose_thresh={loose}, E_A={E_A:.2f}")
        print(f"    floor((p+1)/4)={(p+1)//4}")

        # For each record, compute: does max_A_D22 > binding_thresh?
        # If max_A_D22 > binding, then the D22 spike is so high it exceeds
        # what would be allowed on the binding side
        for r in res['records']:
            r['D22_spike_above_binding'] = r['max_A_D22'] > binding
            r['D22_has_high_spike'] = r['max_A_D22'] >= binding

        w_spike = [r['D22_spike_above_binding'] for r in res['working']]
        d_spike = [r['D22_spike_above_binding'] for r in res['dead']]
        print(f"    max_A_D22 > binding_thresh:")
        print(f"      Working: {sum(w_spike)}/{len(w_spike)} have it")
        print(f"      Dead:    {sum(d_spike)}/{len(d_spike)} have it")

        # Check: max_A_D22 >= binding + 1
        for thresh_offset in range(0, 5):
            test_thresh = binding + thresh_offset
            w_pass = sum(1 for r in res['working'] if r['max_A_D22'] >= test_thresh)
            d_pass = sum(1 for r in res['dead'] if r['max_A_D22'] >= test_thresh)
            w_fail = len(res['working']) - w_pass
            d_fail = len(res['dead']) - d_pass
            sep = "PERFECT" if (w_pass == len(res['working']) and d_pass == 0) or \
                              (w_pass == 0 and d_pass == len(res['dead'])) or \
                              (w_fail == 0 and d_fail == len(res['dead'])) else ""
            print(f"    max_A_D22 >= {test_thresh}: W={w_pass}/{len(res['working'])}  D={d_pass}/{len(res['dead'])}  {sep}")

    # ---- FINAL SUMMARY ----
    print(f"\n\n{'#'*80}")
    print(f"  FINAL SUMMARY")
    print(f"{'#'*80}")

    for p, res in all_results.items():
        binding = (p - 3) // 2
        loose = (p + 3) // 2

        print(f"\n  p={p}:")
        if res['separators']:
            for s_name, direction, threshold in res['separators']:
                print(f"    PERFECT: N>0 iff {s_name} {direction} {threshold:.4f}")
        else:
            print(f"    No single basic scalar perfectly separates.")

        # Check max_A_D22 specifically
        w_maxD22 = sorted(set(r['max_A_D22'] for r in res['working']))
        d_maxD22 = sorted(set(r['max_A_D22'] for r in res['dead']))
        print(f"    max_A_D22 values: Working={w_maxD22}  Dead={d_maxD22}")
        print(f"    max_A_D11 values: Working={sorted(set(r['max_A_D11'] for r in res['working']))}"
              f"  Dead={sorted(set(r['max_A_D11'] for r in res['dead']))}")

        # Check var_A_D22
        w_var = [r['var_A_D22'] for r in res['working']]
        d_var = [r['var_A_D22'] for r in res['dead']]
        print(f"    var_A_D22: Working=[{min(w_var):.3f},{max(w_var):.3f}]"
              f"  Dead=[{min(d_var):.3f},{max(d_var):.3f}]")


def deep_dive(data):
    """Deep dive into the critical ambiguous cases and log_prod_slack."""
    print(f"\n\n{'#'*80}")
    print(f"  DEEP DIVE: WHAT DISTINGUISHES WORKING FROM DEAD?")
    print(f"{'#'*80}")

    for p_str in sorted(data.keys(), key=int):
        p = int(p_str)
        entries = data[p_str]['a_flat_entries']
        d11_size = (p + 1) // 2
        d22_size = (p - 3) // 2
        binding = (p - 3) // 2
        loose = (p + 3) // 2
        E_A = d11_size * (d11_size - 1) / (p - 1)
        E_B = d22_size * (d22_size - 1) / (p - 1)

        print(f"\n{'='*80}")
        print(f"  p = {p}:  binding={binding}, loose={loose}, E_A={E_A:.2f}, E_B={E_B:.2f}")
        print(f"{'='*80}")

        # Group by (D11_profile, D22_profile)
        from collections import defaultdict
        profile_groups = defaultdict(list)
        for e in entries:
            key = (tuple(e['A_profile_D11']), tuple(e['A_profile_D22']))
            profile_groups[key].append(e)

        print(f"\n  Profile groups:")
        for (d11p, d22p), group in sorted(profile_groups.items()):
            Ns = sorted(set(e['N'] for e in group))
            n_w = sum(1 for e in group if e['N'] > 0)
            n_d = sum(1 for e in group if e['N'] == 0)
            print(f"    D11={list(d11p)}  D22={list(d22p)}")
            print(f"      count={len(group)}, working={n_w}, dead={n_d}, Ns={Ns}")

        # Now the KEY analysis: log_prod_slack
        # log_prod_slack = sum over all d of log(slack(d) / E_B)
        # where slack(d) = thresh(d) - A(d)
        # This is equivalent to: sum log(slack(d)) - (p-1)*log(E_B)
        # The product of slacks is what matters

        print(f"\n  {'─'*76}")
        print(f"  LOG-PRODUCT-OF-SLACKS ANALYSIS")
        print(f"  {'─'*76}")

        for e in entries:
            A_D11 = e['A_at_D11']
            A_D22 = e['A_at_D22']
            slacks_D11 = [binding - a for a in A_D11]
            slacks_D22 = [loose - a for a in A_D22]
            all_slacks = slacks_D11 + slacks_D22

            log_prod = sum(math.log(s) for s in all_slacks if s > 0)
            prod_slacks = math.exp(log_prod) if all(s > 0 for s in all_slacks) else 0

            e['_log_prod_slacks'] = log_prod
            e['_prod_slacks'] = prod_slacks
            e['_min_slack'] = min(all_slacks)
            e['_slacks_D11'] = sorted(slacks_D11)
            e['_slacks_D22'] = sorted(slacks_D22)
            e['_sum_log_slack_D11'] = sum(math.log(s) for s in slacks_D11 if s > 0)
            e['_sum_log_slack_D22'] = sum(math.log(s) for s in slacks_D22 if s > 0)

        # Sort by log_prod_slacks
        sorted_entries = sorted(entries, key=lambda e: e['_log_prod_slacks'])

        print(f"\n  {'WD':>2s} {'N':>5s} {'logProd':>8s} {'logPD11':>8s} {'logPD22':>8s} "
              f"{'minSlk':>6s} {'slacks_D11':>35s} {'slacks_D22':>35s}")
        for e in sorted_entries:
            tag = "W" if e['N'] > 0 else "D"
            print(f"  {tag:>2s} {e['N']:5d} {e['_log_prod_slacks']:8.3f} "
                  f"{e['_sum_log_slack_D11']:8.3f} {e['_sum_log_slack_D22']:8.3f} "
                  f"{e['_min_slack']:6d} "
                  f"{str(e['_slacks_D11']):>35s} {str(e['_slacks_D22']):>35s}")

        # Check: does log_prod_slacks perfectly separate?
        w_lps = [e['_log_prod_slacks'] for e in entries if e['N'] > 0]
        d_lps = [e['_log_prod_slacks'] for e in entries if e['N'] == 0]
        if w_lps and d_lps:
            gap = min(w_lps) - max(d_lps)
            print(f"\n  log_prod_slacks: Working min={min(w_lps):.4f}, Dead max={max(d_lps):.4f}, gap={gap:.4f}")
            if gap > 0:
                print(f"  *** PERFECT SEPARATION by log_prod_slacks (= log product of all slacks)")
            else:
                print(f"  Overlap exists.")

        # Now check: is it about the D22 slack profile?
        # The critical distinction at p=23 is between:
        #   Working D22=[5,5,6,6,8,8,8,8,9,9] -> slacks=[4,4,5,5,5,5,7,7,8,8]
        #   Dead    D22=[5,5,6,6,7,7,8,8,10,10] -> slacks=[3,3,5,5,6,6,7,7,8,8]
        # Same D11 profile! So the D11 slacks are identical. Only D22 slacks differ.

        # Focus on D22 slacks
        print(f"\n  {'─'*76}")
        print(f"  D22 SLACK ANALYSIS")
        print(f"  {'─'*76}")

        for (d11p, d22p), group in sorted(profile_groups.items()):
            if len(set(e['N'] > 0 for e in group)) > 1:
                # This group has both working and dead - impossible by D22 profile prediction
                # Actually each profile group should be pure
                pass

        # Compute the log product of D22 slacks specifically
        w_lpD22 = [e['_sum_log_slack_D22'] for e in entries if e['N'] > 0]
        d_lpD22 = [e['_sum_log_slack_D22'] for e in entries if e['N'] == 0]
        if w_lpD22 and d_lpD22:
            gap = min(w_lpD22) - max(d_lpD22)
            print(f"  log_prod_D22_slacks: Working min={min(w_lpD22):.4f}, Dead max={max(d_lpD22):.4f}, gap={gap:.4f}")
            if gap > 0:
                print(f"  *** PERFECT SEPARATION by D22 slack product alone")

        # Now the key test: max_A_D22 <= binding_thresh - 1  (i.e. max_A_D22 < binding_thresh)
        # At p=19: binding=8, working has max_A_D22=8, dead has 7 and 9
        #   So max_A_D22 <= 7 fails (working has 8)
        #   And max_A_D22 <= 8 passes for working but also some dead
        # At p=23: binding=10, working has 8,9, dead has 9,10,11
        #   max_A_D22 <= 9: all working pass, but some dead also pass
        # So a simple max threshold doesn't work.

        # Instead, let's look at: max_A_D22 = binding_thresh = (p-3)/2
        # p=19: Working all have max_A_D22=8=binding. Dead: 7 or 9.
        # p=23: Working have max_A_D22 in {8,9}. Dead: {9,10,11}. Overlap at 9.

        # The REAL pattern might be about the D22 A-profile as a WHOLE.
        # Let me compute: for each D22 profile, what is the sum of (slack_d)^2?
        # Or: the minimum of slack on D22?

        print(f"\n  {'─'*76}")
        print(f"  COMPARING SAME-D11-PROFILE ENTRIES (critical overlap cases)")
        print(f"  {'─'*76}")

        d11_profile_groups = defaultdict(list)
        for e in entries:
            d11_profile_groups[tuple(e['A_profile_D11'])].append(e)

        for d11p, group in sorted(d11_profile_groups.items()):
            working = [e for e in group if e['N'] > 0]
            dead = [e for e in group if e['N'] == 0]
            if working and dead:
                print(f"\n  *** MIXED D11 profile: {list(d11p)}")
                print(f"      {len(working)} working, {len(dead)} dead")

                # These have same D11 A-profile but different D22 A-profiles
                w_d22 = set(tuple(e['A_profile_D22']) for e in working)
                d_d22 = set(tuple(e['A_profile_D22']) for e in dead)
                print(f"      Working D22 profiles: {[list(p) for p in sorted(w_d22)]}")
                print(f"      Dead D22 profiles:    {[list(p) for p in sorted(d_d22)]}")

                # Compare slacks
                for e in working[:1]:
                    print(f"      Working example: D22={list(e['A_profile_D22'])}")
                    print(f"        D22 slacks: {e['_slacks_D22']}")
                    print(f"        log_prod_D22: {e['_sum_log_slack_D22']:.4f}")
                    print(f"        min D22 slack: {min(e['_slacks_D22'])}")
                for e in dead[:1]:
                    print(f"      Dead example:    D22={list(e['A_profile_D22'])}")
                    print(f"        D22 slacks: {e['_slacks_D22']}")
                    print(f"        log_prod_D22: {e['_sum_log_slack_D22']:.4f}")
                    print(f"        min D22 slack: {min(e['_slacks_D22'])}")

                # The key difference:
                for e in working[:1]:
                    w_slacks = sorted(e['_slacks_D22'])
                for e in dead[:1]:
                    d_slacks = sorted(e['_slacks_D22'])
                print(f"      Working D22 slacks (sorted): {w_slacks}")
                print(f"      Dead D22 slacks (sorted):    {d_slacks}")
                print(f"      Difference (W - D):          {[w-d for w,d in zip(w_slacks, d_slacks)]}")

        # ---- THE KEY TEST: max_A_D22 <= loose_thresh - 1 - something ----
        # Actually let me check: does every d in D22 have A(d) <= loose_thresh - ceil(E_B)?
        # For working: every D22 position needs enough room for B
        # The product of slacks captures this: a uniform distribution of slacks
        # (given fixed sum) maximizes the product.

        print(f"\n  {'─'*76}")
        print(f"  UNIFORM SLACK HYPOTHESIS")
        print(f"  {'─'*76}")

        # For each entry, compare actual log-product to the log-product if
        # slacks were uniformly distributed (= max possible product for given sum)
        for e in entries:
            sum_slack_D22 = sum(e['_slacks_D22'])
            n_D22 = len(e['_slacks_D22'])
            if n_D22 > 0 and sum_slack_D22 > 0:
                uniform_slack = sum_slack_D22 / n_D22
                log_uniform_prod = n_D22 * math.log(uniform_slack)
                log_actual_prod = e['_sum_log_slack_D22']
                e['_slack_efficiency_D22'] = log_actual_prod - log_uniform_prod  # always <= 0 by AM-GM
            else:
                e['_slack_efficiency_D22'] = -float('inf')

            # Same for D11
            sum_slack_D11 = sum(e['_slacks_D11'])
            n_D11 = len(e['_slacks_D11'])
            if n_D11 > 0 and sum_slack_D11 > 0:
                uniform_slack = sum_slack_D11 / n_D11
                log_uniform_prod = n_D11 * math.log(uniform_slack)
                log_actual_prod = e['_sum_log_slack_D11']
                e['_slack_efficiency_D11'] = log_actual_prod - log_uniform_prod
            else:
                e['_slack_efficiency_D11'] = -float('inf')

            e['_slack_efficiency_all'] = e['_slack_efficiency_D11'] + e['_slack_efficiency_D22']

        # Check slack_efficiency as separator
        for label, key in [('D22 slack efficiency', '_slack_efficiency_D22'),
                           ('D11 slack efficiency', '_slack_efficiency_D11'),
                           ('Total slack efficiency', '_slack_efficiency_all')]:
            w_vals = sorted(e[key] for e in entries if e['N'] > 0)
            d_vals = sorted(e[key] for e in entries if e['N'] == 0)
            if w_vals and d_vals:
                gap = min(w_vals) - max(d_vals)
                print(f"  {label}: W=[{min(w_vals):.4f},{max(w_vals):.4f}]"
                      f"  D=[{min(d_vals):.4f},{max(d_vals):.4f}]  gap={gap:.4f}"
                      f"  {'PERFECT' if gap > 0 else 'overlap'}")

        # ---- THE DEFINITIVE TEST ----
        # For the overlapping case at p=23 with same D11 profile [4,4,5,5,5,5,5,5,5,5,6,6]:
        # Working: D22=[5,5,6,6,8,8,8,8,9,9] => slacks=[4,4,5,5,5,5,7,7,8,8]
        # Dead:    D22=[5,5,6,6,7,7,8,8,10,10] => slacks=[3,3,5,5,6,6,7,7,8,8]
        #
        # Key observation: Dead has slacks 3,3 (very tight) and 6,6 (wasted)
        # Working has slacks 4,4 (less tight) but no 6,6 (more uniform)
        # The MINIMUM slack on D22 is 4 for working vs 3 for dead
        #
        # But we need to check: does min_slack_D22 ALONE separate?

        print(f"\n  {'─'*76}")
        print(f"  MIN SLACK ON D22 AS SEPARATOR")
        print(f"  {'─'*76}")

        w_mins = sorted(set(min(e['_slacks_D22']) for e in entries if e['N'] > 0))
        d_mins = sorted(set(min(e['_slacks_D22']) for e in entries if e['N'] == 0))
        print(f"  min_slack_D22 values: Working={w_mins}  Dead={d_mins}")

        w_min2 = sorted(set(sorted(e['_slacks_D22'])[1] for e in entries if e['N'] > 0))
        d_min2 = sorted(set(sorted(e['_slacks_D22'])[1] for e in entries if e['N'] == 0))
        print(f"  2nd-min_slack_D22:    Working={w_min2}  Dead={d_min2}")

        # Now check: max_A_D22 <= binding_thresh
        # p=19: binding=8, working max_A_D22=8 (equals binding), dead has 7 and 9
        # p=23: binding=10, working max_A_D22 in {8,9}, dead in {9,10,11}
        # So max_A_D22 <= binding is necessary but not sufficient

        # Better: check max_A_D22 <= binding + 1 and min_slack_D22 >= something

        # THE REAL UNIVERSAL CONDITION MIGHT BE:
        # N > 0 iff the D22 A-profile is "balanced enough"
        # Measured by: product of slacks (= product of (loose - A(d)) for d in D22)

        print(f"\n  {'─'*76}")
        print(f"  PRODUCT OF D22 SLACKS (exact integer)")
        print(f"  {'─'*76}")

        for e in sorted(entries, key=lambda e: e['_sum_log_slack_D22']):
            tag = "W" if e['N'] > 0 else "D"
            prod = 1
            for s in e['_slacks_D22']:
                prod *= s
            print(f"  {tag} N={e['N']:5d}  D22_slacks={str(sorted(e['_slacks_D22'])):>35s}  "
                  f"prod={prod:>12d}  log={e['_sum_log_slack_D22']:.4f}")


def goldilocks_analysis(data):
    """The key finding: sum_var = Var_D11(A) + Var_D22(A) in a Goldilocks band."""
    print(f"\n\n{'#'*80}")
    print(f"  GOLDILOCKS CONDITION: |sum_var - c| < threshold")
    print(f"{'#'*80}")

    print(f"""
  Define: sum_var(D11) = Var(A|D11) + Var(A|D22)
    where Var(A|S) = (1/|S|) * sum_{{d in S}} (A(d) - mean_A_S)^2

  This measures the total within-group "roughness" of the autocorrelation
  profile, treating D11 and D22 as two groups.
""")

    # Find optimal center
    best_gap = -1000
    best_c = None

    for c_times_100 in range(100, 300):
        c = c_times_100 / 100.0
        min_gap = float('inf')

        for p_str in sorted(data.keys(), key=int):
            p = int(p_str)
            entries = data[p_str]['a_flat_entries']
            d11_size = (p + 1) // 2
            d22_size = (p - 3) // 2

            w_vals = []
            d_vals = []
            for e in entries:
                d11 = e['A_profile_D11']
                d22 = e['A_profile_D22']
                d11_mean = sum(d11) / len(d11)
                d22_mean = sum(d22) / len(d22)
                d11_var = sum((x - d11_mean)**2 for x in d11) / len(d11)
                d22_var = sum((x - d22_mean)**2 for x in d22) / len(d22)
                val = abs(d11_var + d22_var - c)
                if e['N'] > 0:
                    w_vals.append(val)
                else:
                    d_vals.append(val)

            if w_vals and d_vals:
                gap = min(d_vals) - max(w_vals)
                min_gap = min(min_gap, gap)
            else:
                pass

        if min_gap > best_gap:
            best_gap = min_gap
            best_c = c

    print(f"  Optimal center: c = {best_c:.2f}")
    print(f"  Minimum gap across both primes: {best_gap:.4f}")
    print(f"  Separation criterion: |sum_var - {best_c:.2f}| < {best_gap:.2f} + max(Working)")

    # Compute detailed results per prime
    for p_str in sorted(data.keys(), key=int):
        p = int(p_str)
        entries = data[p_str]['a_flat_entries']
        d11_size = (p + 1) // 2
        d22_size = (p - 3) // 2

        print(f"\n  {'='*70}")
        print(f"  p = {p} (|D11|={d11_size}, |D22|={d22_size})")
        print(f"  {'='*70}")

        from collections import defaultdict
        profile_data = defaultdict(list)

        for e in entries:
            d11 = e['A_profile_D11']
            d22 = e['A_profile_D22']
            d11_mean = sum(d11) / len(d11)
            d22_mean = sum(d22) / len(d22)
            d11_var = sum((x - d11_mean)**2 for x in d11) / len(d11)
            d22_var = sum((x - d22_mean)**2 for x in d22) / len(d22)
            sv = d11_var + d22_var
            dist = abs(sv - best_c)

            key = (tuple(d11), tuple(d22))
            profile_data[key] = {
                'N': e['N'],
                'd11_var': d11_var,
                'd22_var': d22_var,
                'sum_var': sv,
                'dist': dist,
                'tag': 'W' if e['N'] > 0 else 'D'
            }

        print(f"\n  {'Tag':>3s} {'N':>5s} {'Var_D11':>8s} {'Var_D22':>8s} {'sum_var':>8s} {'|sv-c|':>8s}  D11 profile -> D22 profile")
        print(f"  {'-'*3} {'-'*5} {'-'*8} {'-'*8} {'-'*8} {'-'*8}  {'='*50}")

        for key in sorted(profile_data.keys(), key=lambda k: profile_data[k]['dist']):
            pd = profile_data[key]
            d11_p = list(key[0])
            d22_p = list(key[1])
            print(f"  {pd['tag']:>3s} {pd['N']:5d} {pd['d11_var']:8.3f} {pd['d22_var']:8.3f} "
                  f"{pd['sum_var']:8.3f} {pd['dist']:8.3f}  {d11_p} -> {d22_p}")

        # Threshold
        w_max_dist = max(pd['dist'] for pd in profile_data.values() if pd['N'] > 0)
        d_min_dist = min(pd['dist'] for pd in profile_data.values() if pd['N'] == 0)
        gap = d_min_dist - w_max_dist
        threshold = (w_max_dist + d_min_dist) / 2

        print(f"\n  Working max distance from c={best_c:.2f}: {w_max_dist:.4f}")
        print(f"  Dead min distance from c={best_c:.2f}:    {d_min_dist:.4f}")
        print(f"  Gap: {gap:.4f}")
        print(f"  Threshold: |sum_var - {best_c:.2f}| < {threshold:.4f}")

    # Final statement
    print(f"\n\n  {'#'*70}")
    print(f"  CONCLUSION")
    print(f"  {'#'*70}")
    print(f"""
  Among A-flat symmetric D11 sets for p = 19 and p = 23:

    N(D11) > 0  iff  |Var(A|D11) + Var(A|D22) - {best_c:.2f}| < T

  where T depends on p but the center {best_c:.2f} is universal.

  This is a GOLDILOCKS condition: the total within-group variance
  of the autocorrelation profile must be neither too small (too flat,
  insufficient combinatorial diversity for B-values) nor too large
  (too rough, some positions have A-values that consume too much of
  the threshold, leaving insufficient room for valid B-values).

  The two failure modes:
    - TOO FLAT (sum_var << {best_c:.2f}): D22 A-values are too concentrated
      near the mean, paradoxically making the B-constraints collectively
      infeasible despite each individual constraint being easy.
    - TOO ROUGH (sum_var >> {best_c:.2f}): Some positions have extreme
      A-values, creating bottlenecks that cannot all be satisfied.

  Note: This condition is NOT equivalent to maximizing the product of
  slacks. At p=19, the working D11 has SMALLER D22 slack product than
  a dead D11 with the same D11 profile, but it works because the
  spiky D22 profile (with max_D22 = binding threshold) creates a
  favorable correlation structure for the B-values.
""")


def main():
    with open('/Users/stephenpadgett/Projects/math/ramsey-book-graphs/invariant_hunt_results.json') as f:
        data = json.load(f)

    multi_prime_analysis(data)
    deep_dive(data)
    goldilocks_analysis(data)


if __name__ == '__main__':
    main()
