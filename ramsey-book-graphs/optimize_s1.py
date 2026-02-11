#!/usr/bin/env python3
"""
Optimize the conditioning point s1 to minimize total cost.

Total cost = -log2(Pr[S1 = s1]) + sum(-log2(Phi(z_i)))

The first term is minimized at s1 = E[S1] but the second term needs s1 small
(to make D11 constraints loose). The optimum balances these.

We also need s1 in the feasible range [s1_min, s1_max] from pigeonhole.
"""

import numpy as np
from math import sqrt, log2, pi
from scipy.stats import norm
from fractions import Fraction


def optimize_s1_prime(p):
    """Find optimal s1' (sum of D11 reps) for conditioning."""
    n = (p + 1) // 2
    r = n // 2
    r2 = (p - 1 - n) // 2
    R = r + r2

    # Exact moments
    v = Fraction((p-3)*(p+1), 16*(p-2))
    c_nc = Fraction(-(p+1), 8*(p-2))
    E_B = Fraction(p-3, 4)
    E_A = Fraction(p+1, 4)

    # S1' = sum of D11 reps
    E_S1p = r * E_B
    Var_S1p = float(r * v + r * (r-1) * c_nc)

    # Regression coefficients
    cov_d11_s1 = float(v + (r-1) * c_nc)
    cov_d22_s1 = float(r * c_nc)
    beta_d11 = cov_d11_s1 / Var_S1p
    beta_d22 = cov_d22_s1 / Var_S1p

    # Conditional variances
    var_cond_d11 = float(v) - cov_d11_s1**2 / Var_S1p
    var_cond_d22 = float(v) - cov_d22_s1**2 / Var_S1p
    sig_d11 = sqrt(var_cond_d11)
    sig_d22 = sqrt(var_cond_d22)

    # Thresholds for flat D11 (A(d) = E[A] for all d)
    T_red = float(Fraction(p-5, 4))  # (p-3)/2 - (p+1)/4
    T_blue = float(Fraction(p+5, 4))  # (p+3)/2 - (p+1)/4

    # Pigeonhole
    s1p_max = r * T_red   # sum of D11 thresholds / 2
    # D22 pigeonhole: sum of D22 thresholds
    S_full = float(Fraction((p-1)*(p-3), 4))  # Parseval sum
    s2p_max = r2 * T_blue  # sum of D22 thresholds / 2
    s1p_min = S_full / 2 - s2p_max

    E_S1p_f = float(E_S1p)

    print(f"\np = {p}, r = {r}, r' = {r2}")
    print(f"  E[S1'] = {E_S1p_f:.1f}, Std[S1'] = {sqrt(Var_S1p):.3f}")
    print(f"  s1'_max = {s1p_max:.1f} (D11 pigeonhole)")
    print(f"  s1'_min = {s1p_min:.1f} (D22 pigeonhole)")

    # Scan over s1' values
    best_margin = -1e10
    best_s1p = None
    results = []

    # Scan a range centered around E[S1'] with focus on the feasible region
    for s1p in np.linspace(max(0, s1p_min), min(S_full/2, s1p_max), 1000):
        delta = s1p - E_S1p_f

        # Local CLT cost (for S1')
        z_clt = delta / sqrt(Var_S1p)
        # Pr[S1' = s1'] ~ 1/sqrt(2*pi*Var) * exp(-z^2/2)
        log2_pr = -0.5 * log2(2 * pi * Var_S1p) - z_clt**2 / (2 * np.log(2))

        # Conditional means
        mu_d11 = float(E_B) + beta_d11 * 2 * delta  # S1 = 2*S1'
        mu_d22 = float(E_B) + beta_d22 * 2 * delta

        # z-scores for marginals
        z_red = (T_red - mu_d11) / sig_d11
        z_blue = (T_blue - mu_d22) / sig_d22

        if norm.cdf(z_red) <= 0 or norm.cdf(z_blue) <= 0:
            continue

        # Product of marginals
        log2_prod = r * log2(norm.cdf(z_red)) + r2 * log2(norm.cdf(z_blue))

        # Total margin = budget - cost = (p-1) - (-log2_pr - log2_prod)
        # = (p-1) + log2_pr + log2_prod
        margin = (p - 1) + log2_pr + log2_prod

        results.append((s1p, margin, log2_pr, log2_prod, z_clt, z_red, z_blue))

        if margin > best_margin:
            best_margin = margin
            best_s1p = s1p

    if best_s1p is not None:
        # Find best result
        for s1p, margin, log2_pr, log2_prod, z_clt, z_red, z_blue in results:
            if s1p == best_s1p:
                print(f"\n  OPTIMAL s1' = {best_s1p:.2f}")
                print(f"    z_CLT = {z_clt:.4f}")
                print(f"    -log2(Pr[S1'=s1']) = {-log2_pr:.2f} bits")
                print(f"    z_red = {z_red:.4f}, Phi = {norm.cdf(z_red):.6f}")
                print(f"    z_blue = {z_blue:.4f}, Phi = {norm.cdf(z_blue):.6f}")
                print(f"    -log2(prod Phi) = {-log2_prod:.2f} bits")
                print(f"    Total cost = {-log2_pr - log2_prod:.2f} bits")
                print(f"    Budget = {p-1} bits")
                print(f"    MARGIN = {margin:.2f} bits")
                break

    # Also show what happens at a few interesting points
    print(f"\n  Comparison at key s1' values:")
    print(f"  {'s1p':>10} {'z_CLT':>8} {'CLT cost':>10} {'z_red':>8} {'z_blue':>8} {'prod cost':>10} {'total':>10} {'margin':>8}")
    for s1p, margin, log2_pr, log2_prod, z_clt, z_red, z_blue in results:
        # Show a few representative points
        if abs(s1p - E_S1p_f) < 0.5 or abs(s1p - s1p_max) < 0.5 or abs(s1p - best_s1p) < 0.5 or abs(z_clt - (-1.0)) < 0.05 or abs(z_clt - (-2.0)) < 0.05:
            print(f"  {s1p:10.1f} {z_clt:8.3f} {-log2_pr:10.2f} {z_red:8.3f} {z_blue:8.3f} {-log2_prod:10.2f} {-log2_pr-log2_prod:10.2f} {margin:8.2f}")

    return best_margin


def main():
    print("=" * 70)
    print("OPTIMIZING CONDITIONING POINT s1")
    print("=" * 70)

    for p in [11, 19, 23, 31, 43, 59, 83, 127, 199, 499]:
        if p % 4 != 3:
            continue
        optimize_s1_prime(p)


if __name__ == '__main__':
    main()
