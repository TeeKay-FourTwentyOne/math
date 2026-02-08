"""
Analyze known constructions for m = 3 (mod 4) cases.

Verifies known solutions, analyzes Delta fluctuation patterns,
and extracts insights for attacking n=22.
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ramsey_core import BlockCirculantGraph, verify_construction, Delta, Sigma

# Known constructions for m = 3 (mod 4) cases
# Source: Steven-VO/circulant-Ramsey GitHub + Lidicky et al. 2024
KNOWN = {
    6: {
        "D11": {3, 5, 6, 8},
        "D12": {0, 1, 4, 6, 7},
        "D22": {1, 2, 4, 7, 9, 10},
    },
    8: {
        "D11": {3, 6, 7, 8, 9, 12},
        "D12": {0, 1, 4, 6, 8, 9, 13},
        "D22": {1, 2, 4, 5, 10, 11, 13, 14},
    },
    10: {
        "D11": {4, 5, 7, 9, 10, 12, 14, 15},
        "D12": {0, 1, 2, 6, 7, 10, 11, 13, 17},
        "D22": {1, 2, 3, 6, 8, 11, 13, 16, 17, 18},
    },
    12: {
        "D11": {5, 6, 7, 8, 9, 11, 12, 14, 15, 16, 17, 18},
        "D12": {0, 1, 2, 6, 10, 13, 14, 16, 18, 20, 21},
        "D22": {1, 2, 3, 4, 10, 13, 19, 20, 21, 22},
    },
    14: {
        "D11": {5, 7, 8, 9, 10, 11, 13, 14, 16, 17, 18, 19, 20, 22},
        "D12": {0, 1, 2, 7, 8, 10, 13, 14, 17, 18, 21, 23, 25},
        "D22": {1, 2, 3, 4, 6, 12, 15, 21, 23, 24, 25, 26},
    },
    16: {
        "D11": {6, 7, 8, 10, 11, 12, 14, 15, 16, 17, 19, 20, 21, 23, 24, 25},
        "D12": {0, 1, 2, 3, 8, 11, 12, 13, 15, 18, 20, 21, 24, 27, 29},
        "D22": {1, 2, 3, 4, 5, 9, 13, 18, 22, 26, 27, 28, 29, 30},
    },
    18: {
        "D11": {6, 8, 9, 10, 11, 13, 14, 15, 17, 18, 20, 21, 22, 24, 25, 26, 27, 29},
        "D12": {0, 1, 2, 3, 8, 9, 10, 12, 14, 15, 19, 20, 22, 24, 25, 28, 32},
        "D22": {1, 2, 3, 4, 5, 7, 12, 16, 19, 23, 28, 30, 31, 32, 33, 34},
    },
    20: {
        "D11": {7, 8, 9, 10, 12, 13, 14, 16, 18, 19, 20, 21, 23, 25, 26, 27, 29, 30, 31, 32},
        "D12": {0, 1, 2, 5, 6, 8, 12, 14, 15, 17, 18, 22, 23, 25, 27, 28, 33, 36, 37},
        "D22": {1, 2, 3, 4, 5, 6, 11, 15, 17, 22, 24, 28, 33, 34, 35, 36, 37, 38},
    },
}


def analyze_construction(n_param):
    """Analyze a known construction in detail."""
    if n_param not in KNOWN:
        print(f"No known construction for n={n_param}")
        return

    data = KNOWN[n_param]
    m = 2 * n_param - 1
    N = 4 * n_param - 2

    G = BlockCirculantGraph(n=n_param, D11=data["D11"], D12=data["D12"], D22=data["D22"])
    result = verify_construction(G)

    d1 = len(G.D11) + len(G.D12)
    d2 = len(G.D22) + len(G.D12)

    print(f"\n{'='*70}")
    print(f"n={n_param}, m={m}, N={N}")
    print(f"|D11|={len(G.D11)}, |D12|={len(G.D12)}, |D22|={len(G.D22)}")
    print(f"d1={d1}, d2={d2}")
    print(f"Valid: {result.valid}")
    print(f"Max red common: {result.max_red_common} (threshold {result.red_threshold})")
    print(f"Max blue common: {result.max_blue_common} (threshold {result.blue_threshold})")

    if not result.valid:
        print(f"VIOLATIONS ({len(result.violations)}):")
        for vtype, d, excess in result.violations:
            print(f"  {vtype} d={d} excess={excess}")

    # Check complement relationship
    all_nonzero = set(range(1, m))
    is_complement = (G.D22 == all_nonzero - G.D11)
    print(f"D22 = complement(D11): {is_complement}")

    # Check 0 in D12
    print(f"0 in D12: {0 in G.D12}")

    # Average-case analysis (same as what we did for n=22)
    red_thresh = n_param - 2
    blue_thresh = n_param - 1

    # V1V1 average lambda_red
    avg_lambda_V1V1 = (len(G.D11)**2 + len(G.D12)**2) / m
    blue_bound_V1V1 = 2 * d1 - (N - 2) + blue_thresh

    print(f"\nV1V1 average lambda_red: {avg_lambda_V1V1:.2f}")
    print(f"  Red threshold: {red_thresh}")
    print(f"  Blue bound (max lambda_red for blue): {blue_bound_V1V1}")
    print(f"  Avg exceeds red by: {avg_lambda_V1V1 - red_thresh:.2f}")
    print(f"  Avg exceeds blue by: {avg_lambda_V1V1 - blue_bound_V1V1:.2f}")

    # V2V2 average lambda_red
    D12T = {(-x) % m for x in G.D12}
    avg_lambda_V2V2 = (len(G.D22)**2 + len(G.D12)**2) / m
    blue_bound_V2V2 = 2 * d2 - (N - 2) + blue_thresh

    print(f"\nV2V2 average lambda_red: {avg_lambda_V2V2:.2f}")
    print(f"  Red threshold: {red_thresh}")
    print(f"  Blue bound: {blue_bound_V2V2}")
    print(f"  Avg exceeds red by: {avg_lambda_V2V2 - red_thresh:.2f}")
    print(f"  Avg exceeds blue by: {avg_lambda_V2V2 - blue_bound_V2V2:.2f}")

    # Actual Delta values for V1V1
    print(f"\nV1V1 Delta values (lambda_red for each difference d):")
    deltas_V1V1 = {}
    for d in range(1, m):
        lam = Delta(G.D11, G.D11, d, m) + Delta(G.D12, G.D12, d, m)
        deltas_V1V1[d] = lam

    red_vals = [deltas_V1V1[d] for d in range(1, m) if d in G.D11]
    blue_vals = [deltas_V1V1[d] for d in range(1, m) if d not in G.D11]

    print(f"  Red edges (d in D11): min={min(red_vals)}, max={max(red_vals)}, "
          f"threshold={red_thresh}")
    print(f"  Blue edges (d not in D11): min={min(blue_vals)}, max={max(blue_vals)}")
    blue_lambda_blue = [(N - 2) - 2*d1 + v for v in blue_vals]
    print(f"  Blue lambda_blue: min={min(blue_lambda_blue)}, max={max(blue_lambda_blue)}, "
          f"threshold={blue_thresh}")

    # Delta fluctuation range
    all_vals = list(deltas_V1V1.values())
    print(f"  All lambda_red: min={min(all_vals)}, max={max(all_vals)}, "
          f"range={max(all_vals)-min(all_vals)}")
    print(f"  Average: {sum(all_vals)/len(all_vals):.2f}")

    # Key question: how much deviation from average is needed?
    needed_for_red = avg_lambda_V1V1 - red_thresh
    print(f"\n  Deviation needed to satisfy red: {needed_for_red:.2f} below avg")
    print(f"  Achieved: red edges have lambda in [{min(red_vals)}, {max(red_vals)}], "
          f"vs avg {avg_lambda_V1V1:.2f}")

    # V2V2 Delta values
    deltas_V2V2 = {}
    for d in range(1, m):
        lam = Delta(G.D22, G.D22, d, m) + Delta(D12T, D12T, d, m)
        deltas_V2V2[d] = lam

    red_vals_22 = [deltas_V2V2[d] for d in range(1, m) if d in G.D22]
    blue_vals_22 = [deltas_V2V2[d] for d in range(1, m) if d not in G.D22]

    print(f"\nV2V2 Delta values:")
    all_22 = list(deltas_V2V2.values())
    print(f"  Red edges: min={min(red_vals_22)}, max={max(red_vals_22)}, threshold={red_thresh}")
    print(f"  All: min={min(all_22)}, max={max(all_22)}, range={max(all_22)-min(all_22)}")
    print(f"  Average: {sum(all_22)/len(all_22):.2f}")

    # V1V2 Delta values
    deltas_V1V2 = {}
    for d in range(m):
        lam = Sigma(G.D11, G.D12, d, m) + Delta(G.D12, G.D22, d, m)
        deltas_V1V2[d] = lam

    red_12 = [deltas_V1V2[d] for d in range(m) if d in G.D12]
    blue_12 = [deltas_V1V2[d] for d in range(m) if d not in G.D12]

    print(f"\nV1V2 Delta values:")
    all_12 = list(deltas_V1V2.values())
    print(f"  Red edges: min={min(red_12)}, max={max(red_12)}, threshold={red_thresh}")
    print(f"  All: min={min(all_12)}, max={max(all_12)}, range={max(all_12)-min(all_12)}")
    print(f"  Average: {sum(all_12)/len(all_12):.2f}")

    return {
        "n": n_param, "m": m, "valid": result.valid,
        "avg_V1V1": avg_lambda_V1V1,
        "avg_V2V2": avg_lambda_V2V2,
        "range_V1V1": max(all_vals) - min(all_vals),
        "range_V2V2": max(all_22) - min(all_22),
        "range_V1V2": max(all_12) - min(all_12),
        "d1": d1, "d2": d2,
    }


def extrapolate_to_n22():
    """What would n=22 need?"""
    m = 43
    N = 86
    n_param = 22
    red_thresh = 20
    blue_thresh = 21

    print(f"\n{'='*70}")
    print(f"EXTRAPOLATION TO n=22 (m=43)")
    print(f"{'='*70}")

    # Based on pattern: |D11|=22, |D12|=21, |D22|=20, d1=43, d2=41
    d11_size = 22
    d12_size = 21
    d22_size = 20
    d1 = d11_size + d12_size  # 43
    d2 = d22_size + d12_size  # 41

    avg_V1V1 = (d11_size**2 + d12_size**2) / m
    avg_V2V2 = (d22_size**2 + d12_size**2) / m
    blue_V1V1 = 2 * d1 - (N - 2) + blue_thresh  # 2*43 - 84 + 21 = 23
    blue_V2V2 = 2 * d2 - (N - 2) + blue_thresh  # 2*41 - 84 + 21 = 19
    blue_V1V2 = d1 + d2 - (N - 2) + blue_thresh  # 43+41 - 84 + 21 = 21

    print(f"|D11|={d11_size}, |D12|={d12_size}, |D22|={d22_size}")
    print(f"d1={d1}, d2={d2}")
    print()
    print(f"V1V1: avg_lambda = {avg_V1V1:.2f}")
    print(f"  Red edges need lambda <= {red_thresh}, excess = {avg_V1V1 - red_thresh:.2f}")
    print(f"  Blue edges need lambda <= {blue_V1V1}, excess = {avg_V1V1 - blue_V1V1:.2f}")
    print()
    print(f"V2V2: avg_lambda = {avg_V2V2:.2f}")
    print(f"  Red edges need lambda <= {red_thresh}, excess = {avg_V2V2 - red_thresh:.2f}")
    print(f"  Blue edges need lambda <= {blue_V2V2}, excess = {avg_V2V2 - blue_V2V2:.2f}")
    print()
    print(f"V1V2: avg_lambda = (Sigma+Delta terms, harder to compute without sets)")
    print(f"  Red edges need lambda <= {red_thresh}")
    print(f"  Blue edges need lambda <= {blue_V1V2}")

    # What fluctuation range is needed?
    # For V1V1: need red values <= 20 when avg is ~21.49
    # For V2V2: need red values <= 20 when avg is ~19.58, BUT blue needs <= 19
    print()
    print(f"V1V1: need Delta fluctuation of at least {avg_V1V1 - red_thresh:.2f} below avg for red edges")
    print(f"V2V2: avg {avg_V2V2:.2f} is below red threshold {red_thresh} (OK for red!)")
    print(f"  But blue bound is {blue_V2V2}, need lambda <= {blue_V2V2}")
    print(f"  Avg exceeds blue bound by {avg_V2V2 - blue_V2V2:.2f}")

    # Also try the flipped config
    print(f"\n--- Flipped: |D11|=20, |D12|=21, |D22|=22, d1=41, d2=43 ---")
    d11_size2 = 20
    d22_size2 = 22
    d1_2 = 41
    d2_2 = 43
    avg2_V1V1 = (d11_size2**2 + d12_size**2) / m
    avg2_V2V2 = (d22_size2**2 + d12_size**2) / m
    blue2_V1V1 = 2*d1_2 - 84 + 21  # 19
    blue2_V2V2 = 2*d2_2 - 84 + 21  # 23

    print(f"V1V1: avg = {avg2_V1V1:.2f}, red <= {red_thresh}, blue <= {blue2_V1V1}")
    print(f"  Avg exceeds red by {avg2_V1V1 - red_thresh:.2f}")
    print(f"  Avg exceeds blue by {avg2_V1V1 - blue2_V1V1:.2f}")
    print(f"V2V2: avg = {avg2_V2V2:.2f}, red <= {red_thresh}, blue <= {blue2_V2V2}")
    print(f"  Avg exceeds red by {avg2_V2V2 - red_thresh:.2f}")
    print(f"  Avg exceeds blue by {avg2_V2V2 - blue2_V2V2:.2f}")


def summary_table():
    """Print summary comparison across all cases."""
    print(f"\n{'='*70}")
    print(f"SUMMARY: Average-case excess and fluctuation ranges")
    print(f"{'='*70}")
    print(f"{'n':>3} {'m':>3} {'d1':>3} {'d2':>3} {'avg_11':>7} {'red_t':>5} "
          f"{'excess':>7} {'range_11':>8} {'range_22':>8} {'range_12':>8} {'valid':>5}")
    print("-" * 80)

    for n_param in sorted(KNOWN.keys()):
        info = analyze_construction(n_param)
        if info:
            print(f"{info['n']:3d} {info['m']:3d} {info['d1']:3d} {info['d2']:3d} "
                  f"{info['avg_V1V1']:7.2f} {n_param-2:5d} "
                  f"{info['avg_V1V1']-(n_param-2):7.2f} "
                  f"{info['range_V1V1']:8d} {info['range_V2V2']:8d} "
                  f"{info['range_V1V2']:8d} {'Y' if info['valid'] else 'N':>5}")


if __name__ == "__main__":
    # Analyze each known construction
    results = {}
    for n_param in sorted(KNOWN.keys()):
        results[n_param] = analyze_construction(n_param)

    # Extrapolate to n=22
    extrapolate_to_n22()

    # Summary
    summary_table()
