#!/usr/bin/env python3
"""
Verify the conditional distribution structure when conditioning on S1.

Key question: when the FULL vector (B at D11 reps, B at D22 reps) has
equi-correlation rho = -2/(p-3), and we condition on S1 = sum of D11 reps,
what is the conditional distribution?

This is NOT the same as conditioning on the full sum. The full sum is
deterministic (Parseval), but S1 alone is random.
"""

import numpy as np
from math import sqrt, log2
from scipy.stats import norm
from fractions import Fraction


def analyze_conditioning(p):
    """Verify conditional structure for prime p."""
    n = (p + 1) // 2
    r = n // 2    # D11 reps
    r2 = (p - 1 - n) // 2  # D22 reps
    R = r + r2    # total reps = (p-1)/2

    # Exact moments
    v = Fraction((p-3)*(p+1), 16*(p-2))  # Var[B(d)]
    c_nc = Fraction(-(p+1), 8*(p-2))     # Cov[B(di), B(dj)] non-complementary
    lam_perp = v - c_nc  # = (p+1)(p-1)/(16(p-2))

    # Full covariance matrix: Sigma = v*I + c_nc*(J - I) = (v-c_nc)*I + c_nc*J
    # where J is the all-ones matrix, size R x R

    # Partition: X1 = (B(d1),...,B(dr)), X2 = (B(e1),...,B(e_{r2}))
    # Sigma11 = v*I_r + c_nc*(J_r - I_r) (among D11 reps)
    # Sigma22 = v*I_{r2} + c_nc*(J_{r2} - I_{r2}) (among D22 reps)
    # Sigma12 = c_nc * J_{r x r2} (cross covariance)

    # S1 = 1^T X1 (sum of D11 reps)
    # Var[S1] = 1^T Sigma11 1 = r*v + r*(r-1)*c_nc
    Var_S1 = r * v + r * (r-1) * c_nc
    print(f"\np = {p}, r = {r}, r' = {r2}, R = {R}")
    print(f"  Var[S1] = {float(Var_S1):.6f}")
    print(f"  Std[S1] = {sqrt(float(Var_S1)):.4f}")

    # Cov[X1_i, S1] = Sigma11[i,:] @ 1 = v + (r-1)*c_nc
    cov_x1_s1 = v + (r-1) * c_nc
    print(f"  Cov[B(d_i in D11), S1] = {float(cov_x1_s1):.6f}")

    # Cov[X2_j, S1] = Sigma12[j,:] @ 1 = r * c_nc
    cov_x2_s1 = r * c_nc
    print(f"  Cov[B(e_j in D22), S1] = {float(cov_x2_s1):.6f}")

    # Conditional mean of B(d_i) given S1 = s1:
    # E[B(d_i) | S1 = s1] = E[B] + (cov_x1_s1 / Var_S1) * (s1 - E[S1])
    beta_d11 = cov_x1_s1 / Var_S1
    print(f"  beta_D11 = Cov[B(d_i), S1] / Var[S1] = {float(beta_d11):.6f}")
    print(f"    => E[B(d_i) | S1=s1] = E[B] + {float(beta_d11):.6f} * (s1 - E[S1])")

    # Conditional mean of B(e_j) given S1 = s1:
    beta_d22 = cov_x2_s1 / Var_S1
    print(f"  beta_D22 = Cov[B(e_j), S1] / Var[S1] = {float(beta_d22):.6f}")
    print(f"    => E[B(e_j) | S1=s1] = E[B] + {float(beta_d22):.6f} * (s1 - E[S1])")

    # Conditional variance of B(d_i) given S1:
    # Var[B(d_i) | S1] = Var[B(d_i)] - Cov[B(d_i),S1]^2 / Var[S1]
    var_cond_d11 = v - cov_x1_s1**2 / Var_S1
    var_cond_d22 = v - cov_x2_s1**2 / Var_S1
    print(f"  Var[B(d_i) | S1] (D11) = {float(var_cond_d11):.6f}, sigma = {sqrt(float(var_cond_d11)):.4f}")
    print(f"  Var[B(e_j) | S1] (D22) = {float(var_cond_d22):.6f}, sigma = {sqrt(float(var_cond_d22)):.4f}")

    # Check: with s1 = s*_red (pigeonhole), what are the conditional means?
    E_A = Fraction(p + 1, 4)
    E_B = Fraction(p - 3, 4)
    E_S1 = r * E_B  # since S1' = sum of r reps

    # For flat D11: sum A(d_i in D11 reps) = r * E_A (approximately)
    s_star_red = r * (Fraction(p-3, 2) - E_A)  # = r * (p-5)/4

    print(f"\n  E[S1'] = {float(E_S1):.2f} (sum over D11 reps)")
    print(f"  s*_red/2 = {float(s_star_red):.2f} (D11 pigeonhole for S1')")

    delta_s = s_star_red - E_S1  # = r * ((p-5)/4 - (p-3)/4) = -r/2
    print(f"  delta = s*_red/2 - E[S1'] = {float(delta_s):.2f}")

    # Conditional means at s1' = s*_red/2:
    mu_d11 = E_B + beta_d11 * 2 * delta_s  # factor 2 because S1 = 2*S1'
    mu_d22 = E_B + beta_d22 * 2 * delta_s
    print(f"\n  At s1 = s*_red (optimal):")
    print(f"    E[B(d_i) | S1=s*_red] = {float(mu_d11):.4f}")
    print(f"    E[B(e_j) | S1=s*_red] = {float(mu_d22):.4f}")

    # z-scores for D11 and D22 constraints
    T_red = Fraction(p-5, 4)  # for flat D11
    T_blue = Fraction(3*p + 5, 4)  # for flat D11 (but wrong -- D22 thresh is (p+3)/2 - A(d))

    # Actually T_blue for d in D22: (p+3)/2 - A(d). For flat D11, A(d) at D22 positions
    # is ALSO ~ E[A] = (p+1)/4. So T_blue = (p+3)/2 - (p+1)/4 = (p+5)/4.
    T_blue = Fraction(p + 5, 4)

    sigma_d11 = sqrt(float(var_cond_d11))
    sigma_d22 = sqrt(float(var_cond_d22))

    z_red = float(T_red - mu_d11) / sigma_d11
    z_blue = float(T_blue - mu_d22) / sigma_d22

    print(f"    T_red = {float(T_red):.2f}, z_red = {z_red:.4f}, Phi(z_red) = {norm.cdf(z_red):.6f}")
    print(f"    T_blue = {float(T_blue):.2f}, z_blue = {z_blue:.4f}, Phi(z_blue) = {norm.cdf(z_blue):.6f}")

    # Cost per constraint
    cost_red = -log2(norm.cdf(z_red)) if norm.cdf(z_red) > 0 else 100
    cost_blue = -log2(norm.cdf(z_blue)) if norm.cdf(z_blue) > 0 else 100

    total = r * cost_red + r2 * cost_blue
    budget = p - 1

    print(f"\n    D11 cost: {r} * {cost_red:.4f} = {r * cost_red:.2f} bits")
    print(f"    D22 cost: {r2} * {cost_blue:.4f} = {r2 * cost_blue:.2f} bits")
    print(f"    Total: {total:.2f} bits")
    print(f"    Budget: {budget} bits")
    print(f"    Margin: {budget - total:.2f} bits")

    # Also check the z-score of s_opt in the local CLT
    z_s = float(2 * delta_s) / sqrt(float(4 * Var_S1))  # z-score of s*_red - E[S1]
    # Wait, S1 = 2*S1', Var[S1] = 4*Var[S1']
    # z = (s*_red - E[S1]) / Std[S1] = (2*s_star_red - 2*E_S1) / sqrt(4*Var_S1)
    #   = (2*delta_s) / (2*sqrt(Var_S1)) = delta_s / sqrt(Var_S1)
    z_s = float(delta_s) / sqrt(float(Var_S1))
    print(f"\n    z-score of s_opt in local CLT: {z_s:.4f}")
    print(f"    Phi(z_s) = {norm.cdf(z_s):.6f}")
    print(f"    Local CLT density factor: ~ exp(-{z_s**2/2:.2f}) = {np.exp(-z_s**2/2):.4f}")


def main():
    print("=" * 70)
    print("CONDITIONAL DISTRIBUTION STRUCTURE")
    print("=" * 70)

    for p in [11, 19, 23, 31, 43, 59, 83, 127, 199]:
        if p % 4 != 3:
            continue
        analyze_conditioning(p)


if __name__ == '__main__':
    main()
