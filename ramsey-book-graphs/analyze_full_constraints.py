#!/usr/bin/env python3
"""
Analyze the FULL constraint structure for the conditioning proof.

The proof currently only handles D11-binding constraints:
  A(d) + B(d) <= (p-3)/2  for d in D11

But the valid construction must also satisfy:
  A(d) + B(d) <= (p+3)/2  for d NOT in D11  (V1V1 blue)
  A(d) + B(p-d) <= (p+3)/2  for d NOT in D11  (V2V2 red)

The V1V1 blue constraint at D22 positions can be binding when A(d) is large
at those positions.

Key questions:
1. What are A(d) values at D22 positions for good vs bad D11?
2. How tight are the D22 constraints compared to D11 constraints?
3. What's the effective number of binding constraints?
"""

import numpy as np
from math import sqrt, log2, comb
from scipy.stats import norm
from fractions import Fraction
from itertools import combinations


def get_symmetric_pairs(p):
    pairs = []
    seen = set()
    for d in range(1, p):
        if d not in seen:
            pairs.append((d, (-d) % p))
            seen.add(d)
            seen.add((-d) % p)
    return pairs


def random_symmetric_D11(p, rng):
    n = (p + 1) // 2
    r = n // 2
    half = list(range(1, (p + 1) // 2))
    chosen = rng.choice(half, size=r, replace=False)
    D11 = set()
    for x in chosen:
        D11.add(x)
        D11.add(p - x)
    return D11


def compute_autocorrelation(D11, p):
    """Compute A(d) = Delta(D11, D11, d) for ALL d in {1,...,p-1}."""
    A = {}
    D11_set = set(D11)
    for d in range(1, p):
        A[d] = sum(1 for a in D11 if (a - d) % p in D11_set)
    return A


def analyze_full_constraints(p, num_trials=5000):
    """Analyze both D11 and D22 constraints."""
    n = (p + 1) // 2
    r = n // 2
    d12_size = (p - 1) // 2
    thresh_red = (p - 3) // 2  # D11 binding
    thresh_blue = (p + 3) // 2  # D22 blue
    E_A = (p + 1) / 4
    E_B = (p - 3) / 4

    rng = np.random.default_rng(42)

    print(f"\np = {p}, n = {n}, r = {r}")
    print(f"  D11 binding threshold: {thresh_red}")
    print(f"  D22 blue threshold: {thresh_blue}")
    print(f"  E[A(d)] = {E_A:.2f}, E[B(d)] = {E_B:.2f}")
    print(f"  D11 slack: thresh_red - E[A] - E[B] = {thresh_red - E_A - E_B:.2f}")
    print(f"  D22 slack: thresh_blue - E[A] - E[B] = {thresh_blue - E_A - E_B:.2f}")

    d11_stats = []
    d22_stats = []

    for trial in range(num_trials):
        D11 = random_symmetric_D11(p, rng)
        D22 = set(range(1, p)) - D11
        A = compute_autocorrelation(D11, p)

        # D11 constraints: T_red(d) = thresh_red - A(d) for d in D11
        T_d11 = {d: thresh_red - A[d] for d in D11}
        max_A_d11 = max(A[d] for d in D11)
        min_T_d11 = min(T_d11.values())

        # D22 constraints (V1V1 blue): T_blue(d) = thresh_blue - A(d) for d in D22
        T_d22 = {d: thresh_blue - A[d] for d in D22}
        max_A_d22 = max(A[d] for d in D22)
        min_T_d22 = min(T_d22.values())

        # Effective constraint: for each d, B(d) <= T(d)
        # D11 constraints give T(d) = thresh_red - A(d) ~ (p-3)/4 - deviation
        # D22 constraints give T(d) = thresh_blue - A(d) ~ (p+3)/4 + deviation
        # The D22 constraint is LOOSER by (thresh_blue - thresh_red) = 3

        d11_stats.append({
            'max_A': max_A_d11,
            'min_T': min_T_d11,
            'mean_A': np.mean([A[d] for d in D11]),
            'num_tight': sum(1 for d in D11 if T_d11[d] <= E_B)
        })

        d22_stats.append({
            'max_A': max_A_d22,
            'min_T': min_T_d22,
            'mean_A': np.mean([A[d] for d in D22]),
            'num_tight': sum(1 for d in D22 if T_d22[d] <= E_B)
        })

    d11_max_A = np.array([s['max_A'] for s in d11_stats])
    d11_min_T = np.array([s['min_T'] for s in d11_stats])
    d22_max_A = np.array([s['max_A'] for s in d22_stats])
    d22_min_T = np.array([s['min_T'] for s in d22_stats])
    d22_num_tight = np.array([s['num_tight'] for s in d22_stats])

    print(f"\n  D11 positions (|D11| = {n}):")
    print(f"    max A(d): mean={np.mean(d11_max_A):.2f}, "
          f"P50={np.median(d11_max_A):.0f}, "
          f"P95={np.percentile(d11_max_A, 95):.0f}, "
          f"max={np.max(d11_max_A):.0f}")
    print(f"    min T_red(d): mean={np.mean(d11_min_T):.2f}, "
          f"P05={np.percentile(d11_min_T, 5):.0f}, "
          f"min={np.min(d11_min_T):.0f}")
    print(f"    Pr[feasible (all T_red >= 0)] = {np.mean(d11_min_T >= 0):.4f}")

    print(f"\n  D22 positions (|D22| = {p - 1 - n}):")
    print(f"    max A(d): mean={np.mean(d22_max_A):.2f}, "
          f"P50={np.median(d22_max_A):.0f}, "
          f"P95={np.percentile(d22_max_A, 95):.0f}, "
          f"max={np.max(d22_max_A):.0f}")
    print(f"    min T_blue(d): mean={np.mean(d22_min_T):.2f}, "
          f"P05={np.percentile(d22_min_T, 5):.0f}, "
          f"min={np.min(d22_min_T):.0f}")
    print(f"    Pr[feasible (all T_blue >= 0)] = {np.mean(d22_min_T >= 0):.4f}")
    print(f"    # D22 positions with T_blue <= E[B]: "
          f"mean={np.mean(d22_num_tight):.1f}")

    # For the D22 constraint to be binding, we need A(d) > thresh_blue - E[B]
    # thresh_blue - E[B] = (p+3)/2 - (p-3)/4 = (p+9)/4
    binding_level = (p + 9) / 4
    print(f"\n  D22 constraint binding level: A(d) > {binding_level:.1f}")
    print(f"    (compare E[A] = {E_A:.1f})")
    print(f"    Gap: {binding_level - E_A:.1f} above mean")
    print(f"    In std units: {(binding_level - E_A) / sqrt(E_A):.2f}")

    # The V2V2 constraint: for d in D22, A(d) + B(p-d) <= thresh_blue
    # B(p-d) is the autocorrelation of D12 at p-d. For d in D22, p-d in D11 (since D11 symmetric).
    # So we need B(p-d) <= thresh_blue - A(d) where p-d is in D11.
    # This is a CROSS constraint: for each d in D22, there's a constraint on B at the
    # complementary position p-d in D11.

    # Actually wait: if D11 is symmetric, d in D22 does NOT imply p-d in D11.
    # D11 is symmetric: d in D11 iff p-d in D11. So D22 is also symmetric.
    # d in D22 means p-d in D22 too.

    print(f"\n  V2V2 constraint: A(d) + B(p-d) <= {thresh_blue} for d in D22")
    print(f"  Since D22 is symmetric, p-d in D22 when d in D22.")
    print(f"  B(p-d) = B(d) when D12 is symmetric, but D12 is NOT symmetric in general.")
    print(f"  So the V2V2 constraints involve B evaluated at D22 positions.")

    # The key insight: D12 = {{0}} union S where S is a random (d12_size-1)-subset of {{1,...,p-1}}.
    # D12^T = {{0}} union {{p-s : s in S}}.
    # B(d) = Delta(D12, D12, d) counts pairs in D12 with difference d.
    # The V2V2 constraint uses Delta(D12^T, D12^T, d) = Delta(D12, D12, p-d) = B(p-d).
    #
    # For the V1V1 blue constraint at D22 positions, we need B(d) <= thresh_blue - A(d).
    # For the V2V2 red constraint at D22 positions, we need B(p-d) <= thresh_blue - A(d).
    #
    # Both involve B evaluated at D22 positions (for V1V1 blue, directly; for V2V2,
    # B(p-d) where d in D22 means p-d in D22 too since D22 is symmetric).

    # So the FULL set of B-constraints is:
    # For d in D11: B(d) <= thresh_red - A(d)  [binding, tight]
    # For d in D22: B(d) <= thresh_blue - A(d) [blue, potentially tight]
    # For d in D22: B(p-d) <= thresh_blue - A(d) [V2V2, = B(d') for d' = p-d in D22]
    # The V2V2 constraint at d is: B(p-d) <= thresh_blue - A(d).
    # Since A(d) = A(p-d) (D11 symmetric), this is B(p-d) <= thresh_blue - A(p-d).
    # So V2V2 at d and V1V1 blue at p-d give the SAME constraint.
    # Therefore V2V2 constraints are redundant with V1V1 blue constraints!

    print(f"\n  ** V2V2 at d gives B(p-d) <= thresh_blue - A(d) = thresh_blue - A(p-d) **")
    print(f"  ** V1V1 blue at p-d gives B(p-d) <= thresh_blue - A(p-d) **")
    print(f"  ** These are THE SAME constraint. V2V2 is redundant! **")
    print(f"\n  So the full constraint set is:")
    print(f"    For d in D11: B(d) <= {thresh_red} - A(d)  [red/binding]")
    print(f"    For d in D22: B(d) <= {thresh_blue} - A(d)  [blue]")
    print(f"  Total: {p-1} constraints on B(d) for d = 1,...,{p-1}")


def analyze_constraint_budget(p):
    """Analyze the log2 budget for both D11 and D22 constraints."""
    n = (p + 1) // 2
    r = n // 2
    thresh_red = (p - 3) // 2
    thresh_blue = (p + 3) // 2
    E_A = (p + 1) / 4
    E_B = (p - 3) / 4
    sigma_cond = sqrt((p + 1) * (p - 1) / (16 * (p - 2)) * (1 - 1 / r))

    print(f"\np = {p}: CONSTRAINT BUDGET ANALYSIS")
    print(f"  # D11 constraints (reps): {r}")
    print(f"  # D22 constraints (reps): {(p - 1 - n) // 2}")

    # D11: z_i = (thresh_red - A(d_i) - E[B]) / sigma_cond
    # For A(d_i) = E[A]: z_i = ((p-3)/2 - (p+1)/4 - (p-3)/4) / sigma_cond
    #                        = ((p-3)/4 - (p+1)/4 + (p-3)/4 - (p-3)/4) / sigma_cond
    # Wait let me just compute: thresh_red - E[A] - E[B] = (p-3)/2 - (p+1)/4 - (p-3)/4
    #  = (2(p-3) - (p+1) - (p-3)) / 4 = (2p-6-p-1-p+3)/4 = -4/4 = -1.

    slack_d11 = thresh_red - E_A - E_B  # = -1
    slack_d22 = thresh_blue - E_A - E_B  # = 2

    z_d11 = slack_d11 / sigma_cond  # z for typical D11 constraint
    z_d22 = slack_d22 / sigma_cond  # z for typical D22 constraint

    cost_d11 = -log2(norm.cdf(z_d11))  # bits per D11 constraint
    cost_d22 = -log2(norm.cdf(z_d22))  # bits per D22 constraint

    r_d22 = (p - 1 - n) // 2  # D22 reps

    total_cost = r * cost_d11 + r_d22 * cost_d22
    budget = p - 1  # log2(C(p-1, (p-3)/2)) ~ p-1

    print(f"\n  Slack (thresh - E[A] - E[B]):")
    print(f"    D11: {slack_d11:.2f}")
    print(f"    D22: {slack_d22:.2f}")
    print(f"  sigma_cond: {sigma_cond:.4f}")
    print(f"\n  z-scores (typical):")
    print(f"    D11: z = {z_d11:.4f}, Phi(z) = {norm.cdf(z_d11):.6f}, cost = {cost_d11:.4f} bits")
    print(f"    D22: z = {z_d22:.4f}, Phi(z) = {norm.cdf(z_d22):.6f}, cost = {cost_d22:.4f} bits")
    print(f"\n  Total cost:")
    print(f"    D11: {r} * {cost_d11:.4f} = {r * cost_d11:.2f} bits")
    print(f"    D22: {r_d22} * {cost_d22:.4f} = {r_d22 * cost_d22:.2f} bits")
    print(f"    Total: {total_cost:.2f} bits")
    print(f"    Budget: {budget} bits")
    print(f"    Margin: {budget - total_cost:.2f} bits")

    # Is the D22 cost significant?
    print(f"\n  D22 as fraction of total: {r_d22 * cost_d22 / total_cost:.1%}")

    return total_cost, budget


def main():
    print("=" * 70)
    print("FULL CONSTRAINT ANALYSIS (D11 + D22)")
    print("=" * 70)

    for p in [11, 19, 23, 31, 43, 59, 83, 127]:
        if p % 4 != 3:
            continue
        analyze_full_constraints(p, num_trials=2000)

    print("\n\n" + "=" * 70)
    print("CONSTRAINT BUDGET (D11 + D22 combined)")
    print("=" * 70)

    for p in [11, 19, 23, 31, 43, 59, 83, 127, 199, 499, 997]:
        if p % 4 != 3:
            continue
        analyze_constraint_budget(p)


if __name__ == '__main__':
    main()
