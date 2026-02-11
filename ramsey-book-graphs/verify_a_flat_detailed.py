#!/usr/bin/env python3
"""
Detailed analysis of A-flatness and the conditioning proof.

The key insight: the product of conditional marginals is
  log2(prod) = sum_i log2(Phi(z_i))
where z_i = (T(d_i) - mu_cond) / sigma_cond.

If ALL z_i are near -2/sqrt(p), each Phi(z_i) ~ 1/2 and log2(prod) ~ -r.
But if one z_i is O(1) (i.e., T(d_i) = mu_cond - O(sigma_cond)),
then Phi(z_i) can be much smaller, and the contribution to log2(prod)
can be O(1) bits worse.

We need: log2(prod) >= -(p-1) - O(log p) for the proof to work.
Since log2(prod) >= -r - (sum of extra losses) and r = (p+1)/4,
we need: sum of extra losses <= 3(p-1)/4 - (p+1)/4 = (p-2)/2.

The extra loss per constraint is:
  loss_i = -log2(Phi(z_i)) - 1  (excess over the baseline -1 per constraint)

For z_i = -2/sqrt(p): loss_i ~ 2*sqrt(2/(pi*p))/ln(2) ~ 2.3/sqrt(p)
Total for r constraints: r * 2.3/sqrt(p) ~ 0.58*sqrt(p)

For z_i = -C (constant): loss_i = -log2(Phi(-C)) - 1

Key: even one constraint with T(d_i) = 0 (meaning A(d_i) = thresh = (p-3)/2)
gives z_i = (0 - (p-3)/4) / (sqrt(p)/4) = -(p-3)/sqrt(p) ~ -sqrt(p).
Then Phi(-sqrt(p)) ~ exp(-p/2) and loss_i ~ p/2.
This ALONE exceeds the budget of (p-2)/2!

So the proof DOES need the A-flat condition: we cannot have any A(d) = thresh.

What is the right cutoff?
Need: sum_{i: deviant} loss_i <= (p-2)/2 - 0.58*sqrt(p) - O(log p)
For a single deviant with z_i = -C: loss_i = -log2(Phi(-C)) - 1
This must be <= (p-2)/2 - 0.58*sqrt(p) ~ p/2.

So we need: -log2(Phi(-C)) - 1 <= p/2
Which gives: Phi(-C) >= 2^{-p/2 - 1}
For C = sqrt(p): Phi(-sqrt(p)) ~ exp(-p/2)/sqrt(2*pi*p),
log2(Phi(-sqrt(p))) ~ -p/(2*ln(2)) - 0.5*log2(p) ~ -0.72*p

So a SINGLE T(d_i) = 0 gives loss ~ 0.72*p, which exceeds p/2.

Minimal requirement: T(d_i) >= E_B * (1 - epsilon) for all i.
Equivalently: A(d_i) <= thresh - E_B * (1 - epsilon) = (p-3)/4 * (1 + epsilon)

For the proof to work with a SINGLE deviant, we need:
  loss = -log2(Phi(z_i)) - 1 <= p/4  (say, half the budget)
  Phi(z_i) >= 2^{-p/4 - 1}
  z_i >= -Phi^{-1}(2^{-p/4-1})

For large p: Phi^{-1}(2^{-p/4}) ~ sqrt(p*ln(2)/2) ~ 0.59*sqrt(p)
So z_i >= -0.59*sqrt(p).

Since z_i = (T(d_i) - mu_cond) / sigma_cond,
T(d_i) >= mu_cond - 0.59*sqrt(p) * sigma_cond
       = (p-3)/4 - 0.59*sqrt(p) * sqrt(p)/4
       = (p-3)/4 - 0.15*p
       ~ 0.10*p

So we need T(d_i) >= 0.10*p for all i, which means
A(d_i) <= (p-3)/2 - 0.10*p ~ 0.40*p.

But E[A(d)] = (p+1)/4 ~ 0.25*p and the max A(d) for a random D11
is typically E[A] + O(sqrt(p)) ~ 0.25*p + O(sqrt(p)) << 0.40*p.

So the A-flat condition is VERY mild and holds for all random D11
with overwhelming probability.

But the RIGHT approach is to average over D11 as well, using the
existence argument: MOST D11 are A-flat, and for any A-flat D11,
the conditioning proof gives E[valid D12 | D11] >= 1.
"""

import numpy as np
from math import sqrt, log, log2
from scipy.stats import norm
from fractions import Fraction


def analyze_worst_case(p):
    """Analyze the worst-case product for different A-flatness cutoffs."""
    n = (p + 1) // 2
    r = n // 2
    E_B = (p - 3) / 4
    E_A = (p + 1) / 4
    thresh = (p - 3) // 2
    sigma_cond = sqrt(float(Fraction((p + 1) * (p - 1), 16 * (p - 2))) * (1 - Fraction(1, r)))

    print(f"\np = {p}, n = {n}, r = {r}")
    print(f"  E[A] = {E_A:.2f}, E[B] = {E_B:.2f}, thresh = {thresh}")
    print(f"  sigma_cond = {sigma_cond:.4f}")
    print(f"  Budget: -(p-1) = {-(p-1)}")

    # Baseline: all A(d) = E[A] => all T(d) = thresh - E[A] = (p-5)/4
    T_flat = thresh - E_A
    z_flat = (T_flat - E_B) / sigma_cond
    log2_prod_flat = r * log2(norm.cdf(z_flat))
    margin_flat = log2_prod_flat - (-(p - 1))

    print(f"\n  Baseline (all A = E[A] = {E_A:.1f}):")
    print(f"    T = {T_flat:.2f}, z = {z_flat:.4f}")
    print(f"    Phi(z) = {norm.cdf(z_flat):.6f}")
    print(f"    log2(prod) = {log2_prod_flat:.2f}")
    print(f"    margin = {margin_flat:.2f} bits")

    # What if ONE constraint has A(d) = E[A] + delta?
    print(f"\n  Effect of ONE deviant A(d) = E[A] + delta:")
    for delta in [1, 2, 3, 5, 10, int(sqrt(p)), int(2 * sqrt(p)), int(3 * sqrt(p))]:
        if E_A + delta > thresh:
            continue
        A_dev = E_A + delta
        T_dev = thresh - A_dev
        z_dev = (T_dev - E_B) / sigma_cond
        loss = -log2(norm.cdf(z_dev)) - 1 if norm.cdf(z_dev) > 0 else 999

        # Remaining r-1 constraints at baseline
        log2_prod_with_dev = (r - 1) * log2(norm.cdf(z_flat)) + log2(norm.cdf(z_dev))
        margin_with_dev = log2_prod_with_dev - (-(p - 1))

        print(f"    delta={delta:3d}: A={A_dev:.0f}, T={T_dev:.1f}, z={z_dev:.4f}, "
              f"Phi(z)={norm.cdf(z_dev):.6f}, extra_loss={loss:.2f} bits, "
              f"margin={margin_with_dev:.2f}")

    # The right condition: A(d) <= thresh - C*sqrt(p) for all d
    print(f"\n  A-flat cutoff analysis:")
    for C_factor in [0.5, 1.0, 1.5, 2.0, 3.0]:
        A_max = thresh - C_factor * sqrt(p)
        if A_max <= E_A:
            print(f"    C*sqrt(p) = {C_factor * sqrt(p):.1f}: A_max = {A_max:.1f} <= E[A], trivially ok")
            continue
        T_min = C_factor * sqrt(p)
        z_min = (T_min - E_B) / sigma_cond
        loss_per = -log2(norm.cdf(z_min)) - 1 if norm.cdf(z_min) > 0 else 999
        # Worst case: all r constraints at this level
        log2_prod_worst = r * log2(norm.cdf(z_min)) if norm.cdf(z_min) > 0 else -999
        margin_worst = log2_prod_worst - (-(p - 1))
        print(f"    A <= {A_max:.1f} (T >= {T_min:.1f}): z_min = {z_min:.4f}, "
              f"loss/const = {loss_per:.2f}, worst margin = {margin_worst:.2f}")


def compute_required_flatness(p):
    """Compute the minimum T(d) needed for the proof to work."""
    n = (p + 1) // 2
    r = n // 2
    E_B = (p - 3) / 4
    thresh = (p - 3) // 2
    sigma_cond = sqrt(float(Fraction((p + 1) * (p - 1), 16 * (p - 2))) * (1 - Fraction(1, r)))

    # We need: log2(prod Phi(z_i)) >= -(p-1) - O(log p)
    # With r constraints, average log2(Phi(z_i)) >= -((p-1) + C*log(p)) / r

    # For ALL constraints equal: log2(Phi(z)) >= -(p-1)/r = -(p-1)/((p+1)/4) ~ -4
    # Phi(z) >= 2^{-4} = 1/16
    # z >= Phi^{-1}(1/16) = -1.645 (approx)

    target_per_constraint = -(p - 1) / r - 0.5  # small margin
    Phi_min = 2 ** target_per_constraint
    z_min = norm.ppf(max(Phi_min, 1e-300))
    T_min = E_B + z_min * sigma_cond
    A_max = thresh - T_min

    return {
        'p': p,
        'target_per_constraint': target_per_constraint,
        'Phi_min': Phi_min,
        'z_min': z_min,
        'T_min': T_min,
        'A_max': A_max,
        'E_A': (p + 1) / 4,
        'excess_A': A_max - (p + 1) / 4,
        'excess_in_sqrt_p': (A_max - (p + 1) / 4) / sqrt(p),
    }


def main():
    print("=" * 70)
    print("DETAILED A-FLATNESS ANALYSIS")
    print("=" * 70)

    # First: analyze worst-case for several primes
    for p in [23, 43, 83, 127, 199, 499]:
        if p % 4 != 3:
            continue
        analyze_worst_case(p)

    # Required flatness
    print("\n\n" + "=" * 70)
    print("REQUIRED A-FLATNESS FOR PROOF TO WORK")
    print("=" * 70)
    print(f"{'p':>5} {'r':>4} {'target/r':>10} {'z_min':>8} {'T_min':>8} {'A_max':>8} {'E[A]':>8} {'excess':>8} {'excess/sp':>10}")
    print("-" * 80)

    for p in [11, 19, 23, 31, 43, 47, 59, 83, 127, 199, 499, 997]:
        if p % 4 != 3:
            continue
        r = compute_required_flatness(p)
        print(f"{r['p']:>5} {(p+1)//4:>4} {r['target_per_constraint']:>10.3f} "
              f"{r['z_min']:>8.4f} {r['T_min']:>8.2f} {r['A_max']:>8.2f} "
              f"{r['E_A']:>8.2f} {r['excess_A']:>8.2f} {r['excess_in_sqrt_p']:>10.4f}")

    print("""
INTERPRETATION:
  - 'excess' = A_max - E[A]: how much above average A(d) can be
  - 'excess/sqrt(p)': the normalized excess
  - For the proof to work: A(d) <= A_max for all d in D11

  As p -> infinity, excess/sqrt(p) -> constant (about 2-3 sqrt(p) units).
  This is easily satisfied by random D11 (by concentration of A(d)).

  The right proof approach:
  1. State theorem for "A-flat D11" where max A(d) <= E[A] + C*sqrt(p*log(p))
  2. Prove existence: random D11 is A-flat w.h.p. (Chebyshev + union bound)
  3. The conditioning proof then works for any A-flat D11.
  4. Since A-flat D11 exists, valid construction exists.

  Alternatively, the SIMPLEST approach:
  1. Average over BOTH D11 and D12
  2. E[valid pairs] = E_{D11}[E_{D12}[1{valid} | D11]]
  3. This avoids fixing D11 entirely and uses linearity of expectation
""")


if __name__ == '__main__':
    main()
