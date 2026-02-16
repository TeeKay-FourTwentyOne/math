#!/usr/bin/env python3
"""
Optimize D11 A-profile to maximize hyperplane bound margin.
CORRECT VERSION: uses A = (p+1)/4 for flat D11 (not E_B).

T_d11 = (p-3)/2 - A_d11
T_d22 = (p+3)/2 - A_d22

Constant sum: Σ A(d_i) = (p²-1)/8 over all reps
"""

import numpy as np
from math import comb, log2, log, exp, pi, sqrt
import json


def cycle_pmf(p):
    k = (p - 1) // 2
    pmf = np.zeros(k, dtype=np.float64)
    total = comb(p, k)
    for j in range(k):
        num = p * comb(k - 1, j) * comb(p - k - 1, k - 1 - j)
        denom = (k - j) * total
        pmf[j] = num / denom
    return pmf


def compute_margin(p, pmf, T_vals, S):
    """Compute margin for a list of thresholds."""
    budget = log2(comb(p, (p-1)//2))

    trunc_dists = []
    cost = 0
    for T in T_vals:
        if T < 0:
            return -1e10
        td = pmf[:min(T + 1, len(pmf))]
        Z = td.sum()
        if Z <= 0 or Z < 1e-300:
            return -1e10
        cost -= log2(Z)
        trunc_dists.append(td)

    HR = budget - cost

    # Convolution
    conv = trunc_dists[0].copy()
    for i in range(1, len(trunc_dists)):
        conv = np.convolve(conv, trunc_dists[i])

    pr_Q_HE = conv[S] if S < len(conv) else 0.0
    pr_Q_E = 1.0
    for td in trunc_dists:
        pr_Q_E *= td.sum()

    if pr_Q_E <= 0 or pr_Q_HE <= 0:
        return -1e10

    pr_Q_H_given_E = pr_Q_HE / pr_Q_E
    log2_pr = log2(pr_Q_H_given_E)

    return HR + log2_pr


def main():
    primes = []
    for n in range(11, 500, 2):
        if n % 4 == 3 and all(n % d != 0 for d in range(2, int(n**0.5) + 1)):
            primes.append(n)

    print(f"{'p':>5} {'flat':>8} {'δ=1':>8} {'δ=2':>8} {'δ=3':>8} {'best':>8} {'best_δ':>6} {'best/p':>8}")

    results = []
    for p in primes:
        k = (p - 1) // 2
        n = (p + 1) // 2
        S = k * (k - 1) // 2
        n_d11 = (p + 1) // 4
        n_d22 = (p - 3) // 4
        pmf = cycle_pmf(p)

        E_A = (p + 1) / 4  # Correct average A-value

        # For integer A-values, we need: n_d11 * A_d11 + n_d22 * A_d22 = (p²-1)/8
        # Flat: A_d11 = A_d22 = round((p+1)/4)
        # Check: for p ≡ 3 mod 4, (p+1)/4 is an integer
        A_flat = (p + 1) // 4
        assert n_d11 * A_flat + n_d22 * A_flat == (p * p - 1) // 8, \
            f"Sum check failed: {n_d11 * A_flat + n_d22 * A_flat} != {(p*p-1)//8}"

        # Flat thresholds
        T_d11_flat = (p - 3) // 2 - A_flat  # = (p-7)/4
        T_d22_flat = (p + 3) // 2 - A_flat  # = (p+5)/4 - (p+1)/4 = 1... wait

        # Recompute: T_d22 = (p+3)/2 - A = (p+3)/2 - (p+1)/4 = (2(p+3)-(p+1))/4 = (p+5)/4
        T_d22_flat = (p + 5) // 4

        T_flat = [T_d11_flat] * n_d11 + [T_d22_flat] * n_d22
        margin_flat = compute_margin(p, pmf, T_flat, S)

        margins = [margin_flat]

        best_margin = margin_flat
        best_delta = 0

        for delta in range(1, min(A_flat, 8)):
            A_d11 = A_flat - delta
            # A_d22 = (sum_total - n_d11 * A_d11) / n_d22
            sum_total = (p * p - 1) // 8
            remaining = sum_total - n_d11 * A_d11
            if n_d22 == 0:
                if remaining != 0:
                    margins.append(-1e10)
                    continue
                A_d22_val = 0
            else:
                if remaining % n_d22 != 0:
                    # Non-integer A_d22 — use mixed thresholds
                    A_d22_base = remaining // n_d22
                    A_d22_extra = remaining - n_d22 * A_d22_base
                    # A_d22_extra reps have A = A_d22_base + 1, rest have A_d22_base
                    T_d22_base = (p + 3) // 2 - A_d22_base
                    T_d22_high = (p + 3) // 2 - (A_d22_base + 1)
                    T_d11_val = (p - 3) // 2 - A_d11

                    T_vals = [T_d11_val] * n_d11 + \
                             [T_d22_high] * A_d22_extra + \
                             [T_d22_base] * (n_d22 - A_d22_extra)
                    m = compute_margin(p, pmf, T_vals, S)
                    margins.append(m)
                    if m > best_margin:
                        best_margin = m
                        best_delta = delta
                    continue

                A_d22_val = remaining // n_d22

            T_d11_val = (p - 3) // 2 - A_d11
            T_d22_val = (p + 3) // 2 - A_d22_val

            if T_d22_val < 0:
                margins.append(-1e10)
                continue

            T_vals = [T_d11_val] * n_d11 + [T_d22_val] * n_d22
            m = compute_margin(p, pmf, T_vals, S)
            margins.append(m)
            if m > best_margin:
                best_margin = m
                best_delta = delta

        # Format output
        m_strs = []
        for i, m in enumerate(margins[:4]):
            m_strs.append(f"{m:8.3f}" if m > -1e9 else "     N/A")

        while len(m_strs) < 4:
            m_strs.append("     N/A")

        print(f"{p:5d} {m_strs[0]} {m_strs[1]} {m_strs[2]} {m_strs[3]} "
              f"{best_margin:8.3f} {best_delta:6d} {best_margin/p:8.5f}")

        results.append({
            'p': p,
            'margin_flat': margin_flat,
            'best_margin': best_margin,
            'best_delta': best_delta,
            'best_margin_over_p': best_margin / p,
        })

    # Summary
    print("\n\nCRITICAL PRIMES")
    print("=" * 60)
    for i in range(1, len(results)):
        r0, r1 = results[i-1], results[i]
        if r0['margin_flat'] > 0 and r1['margin_flat'] <= 0:
            print(f"  FLAT margin crosses zero between p={r0['p']} and p={r1['p']}")
        if r0['best_margin'] > 0 and r1['best_margin'] <= 0:
            print(f"  BEST margin crosses zero between p={r0['p']} and p={r1['p']}")

    # Trend analysis
    large = [r for r in results if r['p'] >= 100]
    if len(large) > 5:
        ps = np.array([r['p'] for r in large], dtype=float)
        flat_m = np.array([r['margin_flat'] for r in large])
        best_m = np.array([r['best_margin'] for r in large])
        inv_sqrt_p = 1.0 / np.sqrt(ps)
        A = np.column_stack([np.ones_like(inv_sqrt_p), inv_sqrt_p, inv_sqrt_p**2])

        c_flat = np.linalg.lstsq(A, flat_m / ps, rcond=None)[0]
        c_best = np.linalg.lstsq(A, best_m / ps, rcond=None)[0]

        print(f"\n  flat margin/p → {c_flat[0]:.6f}")
        print(f"  best margin/p → {c_best[0]:.6f}")

    with open('ramsey-book-graphs/optimize_correct_results.json', 'w') as f:
        json.dump(results, f, indent=2, default=float)


if __name__ == '__main__':
    main()
