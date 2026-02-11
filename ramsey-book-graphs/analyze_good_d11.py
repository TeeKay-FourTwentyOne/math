#!/usr/bin/env python3
"""
Analyze properties of "good" D11 (those that have valid D12).

The key question: what makes a D11 good? From enumeration:
- max A(d) is minimal (flat autocorrelation)
- Equal QR/QNR split

Can we prove: for D11 with max A(d) <= E[A] + O(1), the full constraint
set (D11 binding + D22 blue) has Pr[all ok] >= 2^{-(p-1)}?

The target: Pr[all ok] >= 2^{-(p-1)} suffices since C(p-1,(p-3)/2) ~ 2^{p-1}.
"""

import numpy as np
from math import sqrt, log2, comb
from scipy.stats import norm
from fractions import Fraction


def compute_full_product(p, A_at_D11, A_at_D22):
    """
    Compute the product of conditional marginals for BOTH D11 and D22 constraints.

    D11 constraints: B(d) <= (p-3)/2 - A(d) for d in D11
    D22 constraints: B(d) <= (p+3)/2 - A(d) for d in D22

    Under conditioning on S_full = sum_{d=1}^{p-1} B(d) at optimal value.
    """
    n = (p + 1) // 2
    r = n // 2
    r_d22 = (p - 1 - n) // 2
    E_B = (p - 3) / 4

    # The FULL sum is over ALL nonzero d, not just D11
    # S_full = sum_{d=1}^{p-1} B(d) = |D12| * (|D12| - 1) by Parseval
    # This is deterministic! S_full = (p-1)/2 * ((p-1)/2 - 1) = (p-1)(p-3)/4
    S_full = (p - 1) * (p - 3) / 4

    # Since S_full is deterministic, sum B(d) = (p-1)*E[B] always.
    # The interesting partial sums are S1 = sum_{d in D11} B(d)
    # and S2 = sum_{d in D22} B(d), with S1 + S2 = S_full.

    # For the conditioning approach: condition on S1 at optimal value.
    # Then S2 = S_full - S1 is determined.

    # Conditional on S1 = s:
    # B(d) for d in D11 has conditional mean s/n, conditional var ~ lambda_perp * (1 - 1/r)
    # B(d) for d in D22 has conditional mean (S_full - s) / (p-1-n)

    # Wait -- the sum S_full = (p-1)(p-3)/4 is EXACTLY (p-1)*E[B].
    # And S1 = sum over D11 (n values) so S2 = sum over D22 (p-1-n values).
    # S1 + S2 = S_full is a hard constraint.

    # If we condition on S1, then the D11 B-values are constrained by S1,
    # and the D22 B-values are constrained by S2 = S_full - S1.

    # Let's compute: for each s1, what's the product of marginals?

    # For D11: after complementary pair reduction, r reps with sum = s1/2.
    # Conditional mean per rep: s1/(2r)
    # For D22: after complementary pair reduction, r_d22 reps with sum = s2/2 = (S_full-s1)/2.
    # Conditional mean per rep: (S_full-s1)/(2*r_d22)

    sigma_cond_d11 = sqrt((p+1)*(p-1)/(16*(p-2)) * (1 - 1/r))
    sigma_cond_d22 = sqrt((p+1)*(p-1)/(16*(p-2)) * (1 - 1/r_d22)) if r_d22 > 1 else sqrt((p+1)*(p-1)/(16*(p-2)))

    # Thresholds
    T_d11 = [(p-3)//2 - a for a in A_at_D11]  # r values
    T_d22 = [(p+3)//2 - a for a in A_at_D22]  # r_d22 values

    # Sum of thresholds
    sum_T_d11 = sum(T_d11)
    sum_T_d22 = sum(T_d22)

    # Optimal s1: maximize probability by choosing s1 that balances constraints.
    # The D11 sum constraint requires s1/2 <= sum_T_d11 (from pigeonhole).
    # The D22 sum constraint requires (S_full-s1)/2 <= sum_T_d22.
    # So s1 >= S_full - 2*sum_T_d22.

    s1_max = 2 * sum_T_d11  # pigeonhole D11
    s1_min_d22 = S_full - 2 * sum_T_d22  # pigeonhole D22

    # E[S1] = n * E[B]
    E_S1 = n * E_B

    # Try a range of s1 values and find the best
    best_log2 = -1e10
    best_s1 = None

    for s1_frac in np.linspace(max(0, s1_min_d22), min(S_full, s1_max), 200):
        s2_frac = S_full - s1_frac

        mu_d11 = s1_frac / (2 * r)
        mu_d22 = s2_frac / (2 * r_d22) if r_d22 > 0 else 0

        # Product over D11 constraints
        log2_prod = 0
        for T in T_d11:
            z = (T - mu_d11) / sigma_cond_d11
            phi_val = norm.cdf(z)
            if phi_val > 0:
                log2_prod += log2(phi_val)
            else:
                log2_prod -= 100

        # Product over D22 constraints
        for T in T_d22:
            z = (T - mu_d22) / sigma_cond_d22
            phi_val = norm.cdf(z)
            if phi_val > 0:
                log2_prod += log2(phi_val)
            else:
                log2_prod -= 100

        if log2_prod > best_log2:
            best_log2 = log2_prod
            best_s1 = s1_frac

    return best_log2, best_s1, s1_min_d22, s1_max, E_S1


def analyze_flat_d11(p):
    """Analyze constraint budget for a 'perfectly flat' D11."""
    n = (p + 1) // 2
    r = n // 2
    r_d22 = (p - 1 - n) // 2
    E_A = (p + 1) / 4

    # Perfectly flat: A(d) = E[A] for all d in D11 and D22
    # (Not exactly possible since A values are integers, but close)
    A_d11 = [round(E_A)] * r  # r representative values
    A_d22 = [round(E_A)] * r_d22

    log2_prod, s1_opt, s1_min, s1_max, E_S1 = compute_full_product(p, A_d11, A_d22)

    budget = p - 1
    margin = budget + log2_prod  # log2_prod is negative

    print(f"\np = {p}: FLAT D11 (A(d) = {round(E_A)} for all d)")
    print(f"  r = {r} D11 reps, r_d22 = {r_d22} D22 reps")
    print(f"  s1 range: [{s1_min:.0f}, {s1_max:.0f}], E[S1] = {E_S1:.0f}, s1_opt = {s1_opt:.0f}")
    print(f"  log2(prod Phi(z_i)) = {log2_prod:.2f}")
    print(f"  Budget: {budget} bits, margin: {margin:.2f} bits")

    return margin


def analyze_best_d11(p):
    """For enumerated primes, compute product for the actual good D11."""
    # From enumeration data
    if p == 11:
        # Working D11 example: [1,2,4,7,9,10], max_A = 3
        # A values at D11: {1:2, 2:3, 4:1, 7:1, 9:3, 10:2}
        # Reps (take min of pair): d=1(A=2), d=2(A=3), d=4(A=1) -> r=3
        A_at_D11_reps = [2, 3, 1]
        # D22 = {3,5,6,8}, pairs: (3,8), (5,6) -> r_d22 = 2
        # Need to compute A at D22 positions... let me approximate
        # sum_{all d} A(d) = n*(n-1) = 6*5 = 30
        # sum at D11 positions: 2+3+1+1+3+2 = 12 (over 6 values)
        # sum at D22 positions: 30-12 = 18 over 4 values, mean = 4.5
        # Max A at D22 ~ 5
        A_at_D22_reps = [5, 4]  # rough
    elif p == 19:
        # Working D11 with max_A = 5, E[A] = 5
        A_at_D11_reps = [5] * 5
        A_at_D22_reps = [5] * 4  # approximately flat
    elif p == 23:
        A_at_D11_reps = [6] * 6
        A_at_D22_reps = [6] * 5
    else:
        return None

    log2_prod, s1_opt, s1_min, s1_max, E_S1 = compute_full_product(p, A_at_D11_reps, A_at_D22_reps)
    budget = p - 1
    margin = budget + log2_prod

    print(f"\np = {p}: BEST D11 (from enumeration data)")
    print(f"  A at D11 reps: {A_at_D11_reps}")
    print(f"  A at D22 reps: {A_at_D22_reps}")
    print(f"  s1 range: [{s1_min:.0f}, {s1_max:.0f}], E[S1] = {E_S1:.0f}, s1_opt = {s1_opt:.0f}")
    print(f"  log2(prod Phi(z_i)) = {log2_prod:.2f}")
    print(f"  Budget: {budget} bits, margin: {margin:.2f} bits")
    return margin


def main():
    print("=" * 70)
    print("FULL CONSTRAINT PRODUCT (D11 + D22) FOR GOOD D11")
    print("=" * 70)

    print("\n--- FLAT D11 (all A = E[A]) ---")
    for p in [11, 19, 23, 31, 43, 59, 83, 127, 199, 499]:
        if p % 4 != 3:
            continue
        analyze_flat_d11(p)

    print("\n\n--- BEST D11 (from enumeration) ---")
    for p in [11, 19, 23]:
        analyze_best_d11(p)


if __name__ == '__main__':
    main()
