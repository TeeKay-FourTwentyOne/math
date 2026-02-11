#!/usr/bin/env python3
"""
Compute the Gaussian correlation loss c_0 on the Parseval hyperplane.

c_0 = Pr[all B(d_i) <= T(d_i) | hyperplane] / prod Phi(z_i)

For equi-correlated Gaussians with rho = -1/(R-1) on the hyperplane.
"""

import numpy as np
from scipy.stats import norm
from math import sqrt, log2

def compute_c0(p, num_samples=10_000_000):
    """Compute c_0 via Monte Carlo for the Gaussian case."""
    n = (p + 1) // 2
    r = (p + 1) // 4      # D11 reps
    r_prime = (p - 3) // 4  # D22 reps
    R = r + r_prime          # total reps = (p-1)/2

    E_A = (p + 1) / 4
    E_B = (p - 3) / 4
    lam_perp = (p + 1) * (p - 1) / (16 * (p - 2))
    sigma_perp = sqrt(lam_perp)
    sigma_cond = sqrt(lam_perp * (R - 1) / R)

    # Thresholds for perfectly A-flat D11 (A(d) = E[A] for all d)
    T_red = (p - 3) / 2 - E_A   # = (p-7)/4
    T_blue = (p + 3) / 2 - E_A  # = (p+5)/4

    # z-scores
    z_red = (T_red - E_B) / sigma_cond   # ~ -4/sqrt(p)
    z_blue = (T_blue - E_B) / sigma_cond  # ~ +8/sqrt(p)

    # Product of marginals
    log2_prod = r * log2(norm.cdf(z_red)) + r_prime * log2(norm.cdf(z_blue))

    # w_i = (T_i - E[B]) / sigma_perp  (threshold in Y-space)
    w_red = (T_red - E_B) / sigma_perp
    w_blue = (T_blue - E_B) / sigma_perp

    # Monte Carlo: generate Y_i = U_i - U_bar on hyperplane
    rng = np.random.default_rng(42)

    # Process in batches to avoid memory issues
    batch_size = min(num_samples, 1_000_000)
    num_batches = num_samples // batch_size
    total_ok = 0

    for _ in range(num_batches):
        U = rng.standard_normal((batch_size, R))
        U_bar = U.mean(axis=1, keepdims=True)
        Y = U - U_bar  # On hyperplane sum=0, corr = -1/(R-1)

        # Check all constraints
        # First r columns are D11 (threshold w_red), next r' are D22 (threshold w_blue)
        ok_red = np.all(Y[:, :r] <= w_red, axis=1)
        ok_blue = np.all(Y[:, r:] <= w_blue, axis=1)
        total_ok += np.sum(ok_red & ok_blue)

    pr_mc = total_ok / num_samples

    # c_0 = Pr_MC / prod Phi(z_i)
    prod_phi = 2 ** log2_prod
    c0 = pr_mc / prod_phi if prod_phi > 0 else float('inf')

    print(f"p={p:4d}: R={R:3d}, r={r:2d}, r'={r_prime:2d}")
    print(f"  z_red={z_red:.4f}, z_blue={z_blue:.4f}")
    print(f"  w_red={w_red:.4f}, w_blue={w_blue:.4f}")
    print(f"  -log2(prod Phi) = {-log2_prod:.2f}")
    print(f"  Pr_MC = {pr_mc:.6e} ({total_ok}/{num_samples})")
    if pr_mc > 0:
        print(f"  -log2(Pr_MC) = {-log2(pr_mc):.2f}")
        print(f"  c_0 = {c0:.4f}")
        print(f"  log2(c_0) = {log2(c0):.4f}")
    else:
        print(f"  c_0: no events observed (need more samples)")
    print()
    return c0, pr_mc, log2_prod


if __name__ == '__main__':
    print("=" * 60)
    print("GAUSSIAN c_0 ON PARSEVAL HYPERPLANE")
    print("=" * 60)

    for p in [11, 19, 23, 31, 43, 59, 83, 127]:
        if p % 4 != 3:
            continue
        n_samp = 10_000_000 if p <= 43 else 50_000_000
        compute_c0(p, num_samples=n_samp)
