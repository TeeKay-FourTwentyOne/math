#!/usr/bin/env python3
"""
Verify the simplified Parseval-hyperplane approach (no S1 conditioning).

All (p-1)/2 representative B-values live on the hyperplane:
  sum_{reps} B(d_i) = (p-1)(p-3)/8

Conditional on this hyperplane, each B(d_i) has:
  - conditional mean = E[B] = (p-3)/4
  - conditional variance = lambda_perp * (1 - 2/(p-1))

where lambda_perp = (p+1)(p-1)/(16(p-2)).

The z-scores are:
  z_red(d) = (T_red(d) - E[B]) / sigma_cond   for d in D11
  z_blue(d) = (T_blue(d) - E[B]) / sigma_cond  for d in D22

For A-flat D11 with A(d) = E[A] = (p+1)/4:
  T_red = (p-7)/4, T_blue = (p+5)/4, E[B] = (p-3)/4
  z_red = ((p-7)/4 - (p-3)/4) / sigma_cond = -1/sigma_cond
  z_blue = ((p+5)/4 - (p-3)/4) / sigma_cond = 2/sigma_cond
"""

import numpy as np
from math import sqrt, log2
from scipy.stats import norm


def analyze_parseval_approach(p):
    """Analyze the simplified approach using Parseval hyperplane only."""
    n = (p + 1) // 2
    r = n // 2          # D11 reps
    r2 = (p - 1 - n) // 2  # D22 reps
    R = r + r2          # total reps = (p-1)/2

    # Exact moments
    E_B = (p - 3) / 4
    lam_perp = (p + 1) * (p - 1) / (16 * (p - 2))

    # Conditional variance on Parseval hyperplane (one linear constraint on R variables)
    var_cond = lam_perp * (1 - 1/R)  # conditioning removes 1 dof from R equi-correlated vars
    sigma_cond = sqrt(var_cond)

    # For perfectly A-flat D11
    T_red = (p - 7) / 4   # D11 threshold - E[B] = -1
    T_blue = (p + 5) / 4  # D22 threshold - E[B] = +2

    z_red = (T_red - E_B) / sigma_cond   # = -1/sigma_cond
    z_blue = (T_blue - E_B) / sigma_cond  # = 2/sigma_cond

    phi_red = norm.cdf(z_red)
    phi_blue = norm.cdf(z_blue)

    # Product of marginals (r D11 + r2 D22 constraints)
    # But one constraint is redundant (the sum), so effective constraints = R-1 = (p-3)/2
    # Actually no -- the sum is already enforced by Parseval, so ALL R constraints are non-trivial
    # The product of marginals overcounts by the hyperplane constraint
    log2_prod = r * log2(phi_red) + r2 * log2(phi_blue)

    # Budget: log2(C(p-1, (p-3)/2)) ~ p-1 - (1/2)log2(p)
    budget = p - 1  # approximate

    # NO CLT cost -- Parseval is exact!
    total_cost = -log2_prod
    margin = budget - total_cost

    print(f"p={p:4d}: r={r:3d}, r'={r2:3d}, R={R:3d}")
    print(f"  sigma_cond = {sigma_cond:.4f}")
    print(f"  z_red = {z_red:.4f}, Phi = {phi_red:.6f}, cost = {-log2(phi_red):.4f} bits")
    print(f"  z_blue = {z_blue:.4f}, Phi = {phi_blue:.6f}, cost = {-log2(phi_blue):.4f} bits")
    print(f"  Product cost: {r}*{-log2(phi_red):.4f} + {r2}*{-log2(phi_blue):.4f} = {total_cost:.2f} bits")
    print(f"  Budget: {budget} bits")
    print(f"  Margin: {margin:.2f} bits")
    print(f"  E[valid D12] ~ 2^{{{margin:.1f}}} / sqrt(p)")
    print()

    return margin


def compare_approaches(p):
    """Compare S1-conditioning vs Parseval-only approach."""
    n = (p + 1) // 2
    r = n // 2
    r2 = (p - 1 - n) // 2
    R = r + r2

    E_B = (p - 3) / 4
    lam_perp = (p + 1) * (p - 1) / (16 * (p - 2))
    var_cond = lam_perp * (1 - 1/R)
    sigma_cond = sqrt(var_cond)

    T_red = (p - 7) / 4
    T_blue = (p + 5) / 4

    z_red = (T_red - E_B) / sigma_cond
    z_blue = (T_blue - E_B) / sigma_cond

    log2_prod_parseval = r * log2(norm.cdf(z_red)) + r2 * log2(norm.cdf(z_blue))
    margin_parseval = (p - 1) + log2_prod_parseval

    return p, r, r2, -log2_prod_parseval, margin_parseval


def main():
    print("=" * 70)
    print("PARSEVAL-HYPERPLANE APPROACH (no S1 conditioning)")
    print("=" * 70)
    print()

    for p in [11, 19, 23, 31, 43, 59, 83, 127, 199, 499]:
        if p % 4 != 3:
            continue
        analyze_parseval_approach(p)

    print()
    print("=" * 70)
    print("COMPARISON TABLE: Parseval-only vs S1-conditioning")
    print("=" * 70)
    print(f"{'p':>5} {'r':>4} {'r2':>4} {'prod cost':>10} {'margin':>8}")
    print("-" * 35)
    for p in [11, 19, 23, 31, 43, 59, 83, 127, 199, 499]:
        if p % 4 != 3:
            continue
        _, r, r2, cost, margin = compare_approaches(p)
        print(f"{p:5d} {r:4d} {r2:4d} {cost:10.2f} {margin:8.2f}")


if __name__ == '__main__':
    main()
