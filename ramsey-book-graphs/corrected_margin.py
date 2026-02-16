#!/usr/bin/env python3
"""
Compute corrected improved margin accounting for non-achievable Q-mass.

The correction uses the DFT non-negativity bound:
  achievable_fraction ≈ Φ(z)^k where z = E[hat(B)(j)] / std(hat(B)(j))

Corrected improved margin = improved_margin + k × log₂(Φ(z))
"""

from math import comb, log2, sqrt, pi, erf, exp
import numpy as np
import sys

sys.path.insert(0, 'ramsey-book-graphs')


def cycle_pmf(p):
    k = (p - 1) // 2
    total = comb(p, k)
    pmf = np.zeros(k, dtype=np.float64)
    for j in range(k):
        num = p * comb(k - 1, j) * comb(p - k - 1, k - 1 - j)
        denom = (k - j) * total
        pmf[j] = num / denom
    return pmf


def phi(x):
    """Standard normal CDF."""
    return 0.5 * (1 + erf(x / sqrt(2)))


def compute_z_score(p, pmf):
    """Compute z-score for DFT non-negativity."""
    k = (p - 1) // 2
    mean_B = sum(j * pmf[j] for j in range(len(pmf)))
    var_B = sum(j**2 * pmf[j] for j in range(len(pmf))) - mean_B**2

    # E[hat(B)(j)] = k + 2 * mean_B * (-1/2) = k - mean_B
    E_hat = k - mean_B

    # Var[hat(B)(j)] = 4 * var_B * Σ_{d=1}^k cos²(2πd/p) for j=1
    cos2_sum = sum(np.cos(2 * pi * d / p)**2 for d in range(1, k + 1))
    Var_hat = 4 * var_B * cos2_sum

    if Var_hat <= 0:
        return float('inf'), var_B

    z = E_hat / sqrt(Var_hat)
    return z, var_B


def compute_improved_margin(p, delta=1):
    """Compute improved margin using saddle-point (from improved_margin_scaled.py)."""
    from improved_margin_scaled import compute_improved_scaled
    r = compute_improved_scaled(p, delta)
    if r is None:
        return None
    return r


def main():
    max_p = int(sys.argv[1]) if len(sys.argv) > 1 else 500

    primes_3mod4 = [p for p in range(11, max_p + 1)
                    if all(p % d != 0 for d in range(2, int(p**0.5) + 1))
                    and p % 4 == 3]

    # Exact alpha values for small primes
    exact_alpha = {11: 1.94, 19: 3.78, 23: 6.71, 31: 51.98}

    print("CORRECTED IMPROVED MARGIN (accounting for non-achievable Q-mass)")
    print("=" * 100)
    print(f"{'p':>5} {'k':>4} {'z':>6} {'Φ(z)':>6} {'Φ(z)^k':>10} {'log2':>7} "
          f"{'improved':>9} {'corrected':>10} {'α-boost':>8} {'final':>8} {'status':>6}")
    print("-" * 100)

    all_ok = True
    first_fail = None

    for p in primes_3mod4:
        k = (p - 1) // 2
        pmf = cycle_pmf(p)
        f_mode = float(np.max(pmf))
        log2_Cpk = sum(log2((p - i) / (i + 1)) for i in range(k))
        R = k

        # z-score
        z, var_B = compute_z_score(p, pmf)
        phi_z = phi(z)
        ach_frac = phi_z ** k  # estimated achievable fraction
        log2_ach_frac = k * log2(phi_z) if phi_z > 0 else float('-inf')

        # Improved margin
        result = compute_improved_margin(p)
        if result is None:
            print(f"{p:5d} FAILED")
            continue

        improved = result['improved']
        condition = result['condition']

        # Corrected improved margin = improved + log₂(achievable fraction)
        corrected = improved + log2_ach_frac

        # For small primes with exact α: additional boost from P/Q > 1
        alpha_boost = 0.0
        if p in exact_alpha:
            alpha_boost = log2(exact_alpha[p])

        # Final margin:
        # For p ≤ 31 with exact α: use standard margin + α boost + ach correction
        # For p ≥ 43: use improved margin + ach correction
        if p in exact_alpha:
            # N ≥ C(p,k) × α × Σ_ach Q = 2^(standard + log2(α) + log2_ach_frac)
            # = 2^(improved + condition + log2(α) + log2_ach_frac)
            # Hmm, for p ≤ 31, condition > 0 (TB doesn't hold)
            # Standard = improved + condition = log2(C(p,k)) + log2(Pr_Q[E∩H])
            # With α: standard + log2(α)
            # With ach correction: standard + log2(α) + log2_ach_frac
            standard = result['standard']
            final = standard + alpha_boost + log2_ach_frac
        else:
            final = corrected

        status = "OK" if final > 0 else "FAIL"
        if final <= 0 and first_fail is None:
            first_fail = p
            all_ok = False

        if (p <= 200 or p % 100 < 10 or p == first_fail or
                p in [227, 499, 991] or final <= 2):
            print(f"{p:5d} {k:4d} {z:6.3f} {phi_z:6.4f} {ach_frac:10.2e} "
                  f"{log2_ach_frac:7.2f} {improved:9.2f} {corrected:10.2f} "
                  f"{alpha_boost:8.2f} {final:8.2f} {status:>6}")

    print()
    print("=" * 100)
    if all_ok:
        print("ALL corrected margins POSITIVE!")
    else:
        print(f"FIRST FAILURE at p = {first_fail}")

    # Asymptotic analysis
    print("\nAsymptotic: corrected ≈ R × ((1/2)log₂(πk) - c** ) where")
    print("  c** = c* + |log₂(Φ(1))| = 3.264 + 0.250 = 3.514")
    print("  Positive when πk > 2^{2×3.514} = 2^{7.028} = 130.5")
    print("  i.e., k > 41.5, p > 85")
    print(f"  (z → 1 as p → ∞, so Φ(z) → Φ(1) = {phi(1):.4f})")
    print(f"  log₂(Φ(1)) = {log2(phi(1)):.4f}")


if __name__ == '__main__':
    main()
