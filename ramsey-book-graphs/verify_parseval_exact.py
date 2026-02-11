#!/usr/bin/env python3
"""
Verify the Parseval-hyperplane approach with exact budgets.

Budget = log2(C(p-1, (p-3)/2)), not just p-1.
Also check: does E[valid D12] > 1?
"""

import numpy as np
from math import sqrt, log2, comb
from scipy.stats import norm


def verify_exact(p):
    """Verify with exact budget."""
    n = (p + 1) // 2
    k = (p - 3) // 2  # |D12| - 1 = size of random subset S
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

    phi_red = norm.cdf(z_red)
    phi_blue = norm.cdf(z_blue)

    log2_prod = r * log2(phi_red) + r2 * log2(phi_blue)

    # Exact budget
    budget_exact = log2(comb(p - 1, k))
    budget_approx = p - 1

    # E[valid D12] = C(p-1,k) * product of marginals
    # (ignoring correlation loss, which is >= some constant)
    log2_E = budget_exact + log2_prod
    E_valid = 2 ** log2_E

    margin_exact = budget_exact + log2_prod
    margin_approx = budget_approx + log2_prod

    print(f"p={p:4d}: r={r:3d}, r'={r2:3d}")
    print(f"  z_red = {z_red:.4f}, z_blue = {z_blue:.4f}")
    print(f"  Phi_red = {phi_red:.6f}, Phi_blue = {phi_blue:.6f}")
    print(f"  -log2(prod) = {-log2_prod:.2f} bits")
    print(f"  Budget (exact) = log2(C({p-1},{k})) = {budget_exact:.2f}")
    print(f"  Budget (approx) = p-1 = {budget_approx}")
    print(f"  Margin (exact) = {margin_exact:.2f} bits")
    print(f"  E[valid D12] = 2^{{{log2_E:.2f}}} = {E_valid:.2f}")
    print(f"  E > 1? {'YES' if E_valid > 1 else 'NO'}")
    print()


def main():
    print("=" * 70)
    print("PARSEVAL APPROACH - EXACT BUDGET VERIFICATION")
    print("=" * 70)
    print()

    for p in [11, 19, 23, 31, 43, 59, 83, 127, 199, 499]:
        if p % 4 != 3:
            continue
        verify_exact(p)

    # Also check: what's the minimum p where E[valid] > 1?
    print("=" * 70)
    print("SCANNING FOR MINIMUM p WHERE E[valid D12] > 1")
    print("=" * 70)
    for p in range(11, 200, 4):
        if p % 4 != 3:
            continue
        # Check if prime
        if not all(p % i != 0 for i in range(2, int(sqrt(p)) + 1)):
            continue

        n = (p + 1) // 2
        k = (p - 3) // 2
        r = n // 2
        r2 = (p - 1 - n) // 2
        R = r + r2

        E_B = (p - 3) / 4
        lam_perp = (p + 1) * (p - 1) / (16 * (p - 2))
        var_cond = lam_perp * (1 - 1/R)
        sigma_cond = sqrt(var_cond)

        z_red = ((p - 7) / 4 - E_B) / sigma_cond
        z_blue = ((p + 5) / 4 - E_B) / sigma_cond

        log2_prod = r * log2(norm.cdf(z_red)) + r2 * log2(norm.cdf(z_blue))
        budget = log2(comb(p - 1, k))
        log2_E = budget + log2_prod
        E_valid = 2 ** log2_E

        status = "OK" if E_valid > 1 else "FAIL"
        print(f"  p={p:4d}: E[valid] = {E_valid:12.2f}  [{status}]")


if __name__ == '__main__':
    main()
