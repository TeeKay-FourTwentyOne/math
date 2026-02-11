"""
Find the OPTIMAL D11 for each prime: the one maximizing LP margin.
Then analyze what structural property distinguishes good D11 from bad ones.

For small primes, enumerate all D11. For larger primes, sample.
"""

import sys
import os
import json
import time
import math
import numpy as np
from scipy.optimize import linprog
from itertools import combinations

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ramsey_core import Delta


def autocorrelation_fft(indicator, m):
    fft_val = np.fft.fft(indicator)
    autocorr = np.fft.ifft(np.abs(fft_val) ** 2).real
    return np.round(autocorr).astype(int)


def get_pairs(p):
    pairs = []
    seen = set()
    for d in range(1, p):
        if d not in seen:
            pairs.append((d, (-d) % p))
            seen.add(d)
            seen.add((-d) % p)
    return pairs


def lp_margin(p, D11):
    """Compute the LP margin for a given D11. Returns margin (positive = feasible)."""
    n = (p + 1) // 2
    F = (p - 1) // 2
    D11_set = set(D11)

    d11_size = len(D11)
    d12_size = (p - 1) // 2

    ind11 = np.zeros(p, dtype=np.float64)
    for j in D11:
        ind11[j] = 1.0
    dft11 = np.fft.fft(ind11)
    power11 = np.abs(dft11) ** 2

    P0 = d11_size ** 2 + d12_size ** 2
    S = p * (d11_size + d12_size) - P0

    # Variables: P(1)...P(F), t
    c = np.zeros(F + 1)
    c[F] = 1.0

    A_ub_list = []
    b_ub_list = []

    # D11 constraints: (2/p) Σ P(k) cos(2πkd/p) - t ≤ -P0/p
    for d in sorted(D11):
        row = np.zeros(F + 1)
        for k in range(1, F + 1):
            row[k - 1] = (2.0 / p) * math.cos(2 * math.pi * k * d / p)
        row[F] = -1.0
        A_ub_list.append(row)
        b_ub_list.append(-P0 / p)

    # t ≤ n-2
    row_t = np.zeros(F + 1)
    row_t[F] = 1.0
    A_ub_list.append(row_t)
    b_ub_list.append(n - 2)

    # D22 constraints
    D22 = set(range(1, p)) - D11_set
    for d in sorted(D22):
        row = np.zeros(F + 1)
        for k in range(1, F + 1):
            row[k - 1] = (2.0 / p) * math.cos(2 * math.pi * k * d / p)
        A_ub_list.append(row)
        b_ub_list.append(n + 1 - P0 / p)

    A_ub = np.array(A_ub_list)
    b_ub = np.array(b_ub_list)

    A_eq = np.zeros((1, F + 1))
    for k in range(F):
        A_eq[0, k] = 2.0
    b_eq = np.array([S])

    bounds = [(power11[k], None) for k in range(1, F + 1)] + [(None, None)]

    result = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                     bounds=bounds, method='highs', options={'presolve': True})

    if result.success:
        return (n - 2) - result.x[F]  # margin = threshold - optimal_t
    else:
        return -999  # infeasible


def d11_properties(p, D11):
    """Compute structural properties of D11."""
    D11_set = set(D11)
    ind11 = np.zeros(p, dtype=np.float64)
    for j in D11:
        ind11[j] = 1.0
    A = autocorrelation_fft(ind11, p)
    dft = np.fft.fft(ind11)
    power = np.abs(dft) ** 2

    # A(d) stats at D11 positions
    A_d11 = [int(A[d]) for d in D11]
    A_d22 = [int(A[d]) for d in range(1, p) if d not in D11_set]

    # |D̂11(k)|² stats
    pk = power[1:(p + 1) // 2].real

    # QR distribution (for reference)
    QR = set()
    for x in range(1, p):
        QR.add((x * x) % p)

    return {
        "max_A_d11": max(A_d11),
        "min_A_d11": min(A_d11),
        "mean_A_d11": sum(A_d11) / len(A_d11),
        "range_A_d11": max(A_d11) - min(A_d11),
        "max_A_d22": max(A_d22),
        "max_power": float(max(pk)),
        "min_power": float(min(pk)),
        "power_range": float(max(pk) - min(pk)),
        "power_std": float(np.std(pk)),
        "d11_qr": len(D11_set & QR),
        "d11_qnr": len(D11_set) - len(D11_set & QR),
    }


def main():
    print("=" * 90)
    print("D11 OPTIMIZATION: Finding best D11 for each prime")
    print("=" * 90)

    # Known solutions for comparison
    KNOWN_D11 = {
        11: {1, 2, 4, 7, 9, 10},
        19: {1, 2, 3, 6, 8, 11, 13, 16, 17, 18},
        23: {5, 6, 7, 8, 9, 11, 12, 14, 15, 16, 17, 18},
        31: {6, 7, 8, 10, 11, 12, 14, 15, 16, 17, 19, 20, 21, 23, 24, 25},
        43: {1, 2, 5, 10, 11, 13, 16, 17, 18, 19, 20, 23, 24, 25, 26, 27,
             30, 32, 33, 38, 41, 42},
    }

    # p=11: enumerate ALL D11
    for p in [11, 19]:
        n = (p + 1) // 2
        pairs = get_pairs(p)
        num_pairs_needed = (p + 1) // 4
        total = math.comb(len(pairs), num_pairs_needed)

        print(f"\n{'='*80}")
        print(f"p = {p}, n = {n}")
        print(f"Total symmetric D11: {total}")
        print(f"{'='*80}")

        if total > 50000:
            print(f"  Too many to enumerate, sampling instead.")
            # Sample
            rng = np.random.default_rng(42)
            best_margin = -999
            best_d11 = None
            num_feasible = 0
            num_tested = 5000

            for trial in range(num_tested):
                chosen = rng.choice(len(pairs), size=num_pairs_needed, replace=False)
                D11 = set()
                for i in chosen:
                    D11.add(pairs[i][0])
                    D11.add(pairs[i][1])

                margin = lp_margin(p, D11)
                if margin > -900:
                    num_feasible += 1
                if margin > best_margin:
                    best_margin = margin
                    best_d11 = D11.copy()

            print(f"  Sampled {num_tested}: {num_feasible} feasible ({100*num_feasible/num_tested:.1f}%)")
        else:
            best_margin = -999
            best_d11 = None
            margins = []
            num_feasible = 0

            for idx, combo in enumerate(combinations(range(len(pairs)), num_pairs_needed)):
                D11 = set()
                for i in combo:
                    D11.add(pairs[i][0])
                    D11.add(pairs[i][1])

                margin = lp_margin(p, D11)
                margins.append(margin)
                if margin > -900:
                    num_feasible += 1
                if margin > best_margin:
                    best_margin = margin
                    best_d11 = D11.copy()

            print(f"  Enumerated all {total}: {num_feasible} feasible ({100*num_feasible/total:.1f}%)")
            if margins:
                feas_margins = [m for m in margins if m > -900]
                if feas_margins:
                    print(f"  Margin distribution: min={min(feas_margins):.4f}, "
                          f"max={max(feas_margins):.4f}, mean={sum(feas_margins)/len(feas_margins):.4f}")

        if best_d11:
            props = d11_properties(p, best_d11)
            print(f"\n  BEST D11: {sorted(best_d11)}")
            print(f"  LP margin: {best_margin:.4f}")
            print(f"  A(d) at D11: max={props['max_A_d11']}, min={props['min_A_d11']}, "
                  f"range={props['range_A_d11']}")
            print(f"  |D̂11(k)|²: max={props['max_power']:.2f}, min={props['min_power']:.2f}, "
                  f"std={props['power_std']:.2f}")
            print(f"  QR balance: {props['d11_qr']} QR, {props['d11_qnr']} QNR")

        # Known D11 properties
        known = KNOWN_D11[p]
        known_margin = lp_margin(p, known)
        known_props = d11_properties(p, known)
        print(f"\n  KNOWN D11: {sorted(known)}")
        print(f"  LP margin: {known_margin:.4f}")
        print(f"  A(d) at D11: max={known_props['max_A_d11']}, min={known_props['min_A_d11']}, "
              f"range={known_props['range_A_d11']}")
        print(f"  |D̂11(k)|²: max={known_props['max_power']:.2f}, min={known_props['min_power']:.2f}, "
              f"std={known_props['power_std']:.2f}")
        print(f"  QR balance: {known_props['d11_qr']} QR, {known_props['d11_qnr']} QNR")

        # Compare best vs known
        if best_d11 and best_d11 != known:
            print(f"\n  Is known D11 optimal? {'YES' if abs(best_margin - known_margin) < 1e-6 else 'NO'}")
            if abs(best_margin - known_margin) > 1e-6:
                print(f"  Best margin {best_margin:.4f} vs known margin {known_margin:.4f}")

    # Correlation analysis: what predicts LP margin?
    print(f"\n{'='*90}")
    print("PROPERTY ANALYSIS: What predicts LP feasibility?")
    print(f"{'='*90}")

    for p in [23, 31, 43]:
        n = (p + 1) // 2
        pairs = get_pairs(p)
        num_pairs_needed = (p + 1) // 4
        rng = np.random.default_rng(42)

        print(f"\n  p={p}: Sampling 500 D11, analyzing properties...")
        num_sample = 500
        data = []

        for trial in range(num_sample):
            chosen = rng.choice(len(pairs), size=num_pairs_needed, replace=False)
            D11 = set()
            for i in chosen:
                D11.add(pairs[i][0])
                D11.add(pairs[i][1])

            margin = lp_margin(p, D11)
            props = d11_properties(p, D11)
            data.append({
                "margin": margin,
                "feasible": margin > -900,
                **props,
            })

        feasible = [d for d in data if d["feasible"]]
        infeasible = [d for d in data if not d["feasible"]]

        print(f"    Feasible: {len(feasible)}/{num_sample}")

        if feasible and infeasible:
            for key in ["max_A_d11", "range_A_d11", "max_power", "power_std", "max_A_d22"]:
                f_vals = [d[key] for d in feasible]
                i_vals = [d[key] for d in infeasible]
                print(f"    {key:15s}: feasible={sum(f_vals)/len(f_vals):.2f}, "
                      f"infeasible={sum(i_vals)/len(i_vals):.2f}")

        # Check known D11
        known = KNOWN_D11[p]
        known_margin = lp_margin(p, known)
        known_props = d11_properties(p, known)
        print(f"    Known D11: margin={known_margin:.4f}, "
              f"max_A={known_props['max_A_d11']}, "
              f"power_std={known_props['power_std']:.2f}")

    # Large prime test: find best random D11
    print(f"\n{'='*90}")
    print("LARGE PRIME SEARCH: Best D11 from random sampling")
    print(f"{'='*90}")

    for p in [47, 59, 67, 71, 79, 83]:
        n = (p + 1) // 2
        pairs = get_pairs(p)
        num_pairs_needed = (p + 1) // 4
        rng = np.random.default_rng(42)

        num_sample = min(2000, max(500, 100000 // p))
        best_margin = -999
        best_d11 = None
        num_feasible = 0

        t0 = time.time()
        for trial in range(num_sample):
            chosen = rng.choice(len(pairs), size=num_pairs_needed, replace=False)
            D11 = set()
            for i in chosen:
                D11.add(pairs[i][0])
                D11.add(pairs[i][1])

            margin = lp_margin(p, D11)
            if margin > -900:
                num_feasible += 1
            if margin > best_margin:
                best_margin = margin
                best_d11 = D11.copy()

        elapsed = time.time() - t0
        print(f"\n  p={p}: {num_feasible}/{num_sample} feasible ({elapsed:.1f}s)")
        if best_d11 and best_margin > -900:
            props = d11_properties(p, best_d11)
            print(f"    Best margin: {best_margin:.4f}")
            print(f"    Best D11: max_A_d11={props['max_A_d11']}, "
                  f"power_std={props['power_std']:.2f}")
        else:
            print(f"    No feasible D11 found in sample")
            # Report the best (least infeasible)
            if best_d11:
                props = d11_properties(p, best_d11)
                print(f"    Least infeasible: margin={best_margin:.4f}, "
                      f"max_A_d11={props['max_A_d11']}")


if __name__ == "__main__":
    main()
