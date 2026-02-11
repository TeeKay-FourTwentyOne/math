#!/usr/bin/env python3
"""Exact moment computations for the first moment proof.

For a random D12 ⊂ Z_p with |D12| = s = (p-1)/2 and 0 ∈ D12:
- D12 = {0} ∪ S where S ⊂ {1,...,p-1} with |S| = s-1 = (p-3)/2

B(d) = Δ(D12, D12, d) = #{a ∈ D12 : (a-d) ∈ D12}

For d ≠ 0:
B(d) = Σ_{a ∈ D12} 1_{(a-d) mod p ∈ D12}
     = 1_{(-d mod p) ∈ D12}           [a=0 term]
       + 1_{d ∈ D12}                   [a=d, checking (a-d)=0 ∈ D12]
       + Σ_{a ∈ S, a≠d} 1_{(a-d) ∈ S} [other terms, both a and a-d nonzero]

Let's define:
  Y_a = 1_{a ∈ S} for a ∈ {1,...,p-1}

Then for d ≠ 0:
  B(d) = Y_{p-d}                      [a=0: is (-d) = (p-d) in S?]
       + Y_d                           [a=d: is d in D12? (a-d)=0 always in D12]
       + Σ_{a≠0,d} Y_a · Y_{(a-d) mod p}  [where (a-d) ≠ 0 since a≠d]

Wait, the last sum: for a ∈ S with a ≠ d, we check if (a-d) ∈ D12.
If (a-d) = 0, that's a=d which we excluded. If (a-d) ≠ 0, check if (a-d) ∈ S.
If (a-d) = -d (i.e., a=0), but a ∈ S means a ≠ 0. So all terms have a-d ≠ 0 and a ≠ 0.

So B(d) = Y_{p-d} + Y_d + Σ_{a ∈ {1,...,p-1}\\{d}} Y_a · Y_{(a-d) mod p}
         (where (a-d) mod p is never 0 since a ≠ d)

This is a sum of indicators of a random (s-1)-subset S of {1,...,p-1}.

KEY SIMPLIFICATION: If we don't condition on 0 ∈ D12, and instead take
D12 as a uniform random s-subset of Z_p, then B(d) has the standard
hypergeometric-type distribution. The 0-conditioning introduces mild
corrections of order O(1/p).

This script computes the exact E[B(d)], Var[B(d)], and verifies against
Monte Carlo.
"""

import numpy as np
from math import comb, factorial
import time
from fractions import Fraction


def exact_E_B(p):
    """Exact E[B(d)] for random D12 with |D12|=(p-1)/2 and 0 ∈ D12.

    S is a random (s-1)-subset of {1,...,p-1} where s = (p-1)/2.
    k = s - 1 = (p-3)/2 = |S|
    N = p - 1 = |{1,...,p-1}|

    B(d) = Y_{p-d} + Y_d + Q(d)
    where Q(d) = Σ_{a≠d, a∈{1..p-1}} Y_a · Y_{(a-d) mod p}

    E[Y_a] = k/N = (p-3)/(2(p-1))
    E[Y_a · Y_b] for a ≠ b = k(k-1)/(N(N-1)) = (p-3)(p-5)/(4(p-1)(p-2))

    For the Q(d) sum: there are (p-2) terms (a=1,...,p-1 except a=d).
    For each, (a-d) mod p ∈ {1,...,p-1}\\{0} and (a-d) ≠ a (since d≠0).
    Also (a-d) could equal d (when a = 2d mod p) or p-d (when a = 0 mod p, excluded).

    Actually, for a ∈ {1,...,p-1}\\{d}:
      - (a-d) mod p ∈ {1,...,p-1}\\{0}
      - (a-d) ≠ a (since d ≠ 0)
      - Is it possible that (a-d) = d? Only if a = 2d mod p.
      - Is it possible that a = p-d? Yes, then (a-d) = p-2d mod p.

    Key: for each a in the sum, we need E[Y_a · Y_{(a-d)}].
      - If a ≠ (a-d) mod p: E = k(k-1)/(N(N-1))
      - If a = (a-d): impossible since d ≠ 0

    So E[Q(d)] = (p-2) × k(k-1)/(N(N-1))

    Special terms Y_{p-d} and Y_d:
    If d = p-d (i.e., 2d = p, impossible since p odd): both are distinct.
    E[Y_{p-d}] = E[Y_d] = k/N
    If d ≠ p-d: E[Y_{p-d} + Y_d] = 2k/N
    """
    s = (p - 1) // 2  # |D12|
    k = s - 1         # |S| = |D12| - 1 (non-zero elements)
    N = p - 1         # |{1,...,p-1}|

    # E[Y_a] = k/N
    E_Y = Fraction(k, N)

    # E[Y_a · Y_b] for a ≠ b = k(k-1)/(N(N-1))
    E_YY = Fraction(k * (k - 1), N * (N - 1))

    # E[Q(d)] for d ≠ 0: (p-2) terms, each contributing E_YY
    # But wait: among the p-2 values of a ∈ {1,...,p-1}\\{d},
    # the map a → (a-d) mod p is a bijection from {1,...,p-1}\\{d} to {1,...,p-1}\\{0}
    # = {1,...,p-1}. But wait, 0 maps to p-d, and d maps to 0.
    # For a ∈ {1,...,p-1}\\{d}: (a-d) mod p takes values in {1,...,p-1}\\{0}
    # Wait no: when a ranges over {1,...,p-1}\\{d}, (a-d) mod p ranges over
    # {1-d, 2-d, ..., (p-1)-d}\\{0} mod p = {1,...,p-1}\\{-d} ... hmm no.
    # (a-d) for a = 1,...,p-1, a≠d: takes values p-d+1,...,p-1,1,...,p-d-1
    # (mod p), which is {1,...,p-1}\\{0} (since a≠d means a-d≠0).
    # Wait: a-d mod p = 0 iff a=d, which we excluded. So (a-d) ranges
    # over p-2 distinct values in {1,...,p-1}.

    # For each such pair (a, a-d): are a and (a-d) always distinct?
    # a = (a-d) mod p iff d = 0, which we excluded. So yes, always distinct.
    # Therefore E[Y_a · Y_{(a-d)}] = E_YY for ALL terms.

    E_Q = (p - 2) * E_YY

    # Y_{p-d} and Y_d: always distinct since d ≠ p-d (p odd, d ≠ 0 ≠ p/2)
    E_total = 2 * E_Y + E_Q

    return E_total


def exact_E_B_simplified(p):
    """Simplified formula for E[B(d)].

    E[B(d)] = 2(p-3)/(2(p-1)) + (p-2)(p-3)(p-5)/(4(p-1)(p-2))
            = (p-3)/(p-1) + (p-3)(p-5)/(4(p-1))
            = (p-3)/(p-1) × [1 + (p-5)/4]
            = (p-3)/(p-1) × (p-1)/4
            = (p-3)/4

    Wait, that simplifies nicely!
    """
    return Fraction(p - 3, 4)


def mc_E_B(p, num_trials=1000000):
    """Monte Carlo estimate of E[B(d)] for d=1."""
    s = (p - 1) // 2
    k = s - 1
    rng = np.random.default_rng(42)

    d = 1  # any nonzero d
    B_samples = []

    for _ in range(num_trials):
        S = set(rng.choice(range(1, p), size=k, replace=False).tolist())
        D12 = S | {0}

        B = 0
        for a in D12:
            if (a - d) % p in D12:
                B += 1
        B_samples.append(B)

    return np.mean(B_samples), np.var(B_samples)


def exact_Var_B(p):
    """Exact Var[B(d)] for random D12 with 0 ∈ D12.

    B(d) = Y_{p-d} + Y_d + Q(d)
    where Q(d) = Σ_{a∈T} Y_a · Y_{(a-d)}
    and T = {1,...,p-1}\\{d}, each (a-d) ∈ {1,...,p-1}.

    Let P_a = Y_a · Y_{(a-d) mod p} for a ∈ T.

    Var[B(d)] = Var[Y_{p-d}] + Var[Y_d] + Var[Q(d)]
              + 2Cov[Y_{p-d}, Y_d]
              + 2Cov[Y_{p-d}, Q(d)]
              + 2Cov[Y_d, Q(d)]

    This requires computing E[Y_a²], E[Y_a Y_b], E[Y_a Y_b Y_c],
    E[Y_a Y_b Y_c Y_d] for distinct indices.
    """
    s = (p - 1) // 2
    k = s - 1
    N = p - 1

    # Factorial moments of hypergeometric
    # E[Y_a] = k/N
    # E[Y_a²] = E[Y_a] = k/N (binary)
    # E[Y_a Y_b] = k(k-1)/(N(N-1))  for a≠b
    # E[Y_a Y_b Y_c] = k(k-1)(k-2)/(N(N-1)(N-2))  for distinct a,b,c
    # E[Y_a Y_b Y_c Y_d] = k(k-1)(k-2)(k-3)/(N(N-1)(N-2)(N-3))

    q1 = Fraction(k, N)
    q2 = Fraction(k * (k-1), N * (N-1))
    q3 = Fraction(k * (k-1) * (k-2), N * (N-1) * (N-2))
    q4 = Fraction(k * (k-1) * (k-2) * (k-3), N * (N-1) * (N-2) * (N-3))

    # Var[Y_a] = q1(1-q1) = q1 - q1²
    var_Y = q1 * (1 - q1)

    # Cov[Y_a, Y_b] = q2 - q1² for a ≠ b
    cov_YY = q2 - q1 * q1

    # Now decompose B(d) = Y_{p-d} + Y_d + Σ_{a∈T} P_a
    # where P_a = Y_a · Y_{(a-d)}

    # Let's denote the special indices:
    # α = p-d, β = d (the two standalone Y terms)
    # These are distinct since d ≠ p-d (p odd, d ≠ p/2).

    # The Q sum has p-2 terms with a ∈ T = {1,...,p-1}\\{d}.
    # P_a = Y_a · Y_{(a-d) mod p}

    # E[P_a] = q2 (since a ≠ (a-d) when d ≠ 0)
    E_P = q2

    # E[P_a²] = E[Y_a² · Y_{(a-d)}²] = E[Y_a · Y_{(a-d)}] = q2
    # (since Y's are binary: Y² = Y)
    E_P2 = q2

    # Var[P_a] = E_P2 - E_P² = q2 - q2²
    var_P = q2 - q2 * q2

    # Cov[P_a, P_b] for a ≠ b in T:
    # P_a P_b = Y_a Y_{a-d} Y_b Y_{b-d}
    # Need to check how many of {a, a-d, b, b-d} are distinct.

    # Case 1: all 4 distinct → E[P_a P_b] = q4
    # Case 2: exactly one collision → E[P_a P_b] = q3
    # Case 3: two collisions → E[P_a P_b] = q2

    # Possible collisions for a ≠ b, both in T, d ≠ 0:
    # a = b: excluded (a ≠ b)
    # a = b-d: i.e., b = a+d
    # a-d = b: i.e., b = a-d
    # a-d = b-d: i.e., a = b, excluded

    # So possible collisions: a = b-d OR a-d = b
    # a = b-d means b = a+d (one collision: a = b-d)
    # a-d = b (one collision)
    # Both: a = b-d AND a-d = b → b = a+d AND b = a-d → d = -d → 2d=0 → d=0 (excluded)

    # So:
    # For generic (a,b): 4 distinct indices → E[P_a P_b] = q4
    # For b = a+d mod p: {a, a-d, a+d, a} collision a appears twice...
    #   wait, b = a+d, so b-d = a. Indices: {a, a-d, a+d, a} → {a, a-d, a+d}, 3 distinct.
    #   E[P_a P_b] = E[Y_a · Y_{a-d} · Y_{a+d} · Y_a] = E[Y_a · Y_{a-d} · Y_{a+d}] = q3
    #   (since Y_a² = Y_a, and {a, a-d, a+d} are distinct iff d ≠ 0 and 2d ≠ 0)
    # For b = a-d mod p (i.e., a-d = b): Indices: {a, a-d, a-d, a-2d} → {a, a-d, a-2d}
    #   if a-2d ≠ a and a-2d ≠ a-d, i.e., 2d ≠ 0 and d ≠ 0 (both true).
    #   E[P_a P_b] = E[Y_a · Y_{a-d} · Y_{a-2d}] = q3
    #   (3 distinct indices)

    # But we also need to check:
    # When b = a+d: is b ∈ T = {1,...,p-1}\\{d}?
    #   b = a+d, need a+d ≠ d, i.e., a ≠ 0 (true since a ∈ T ⊂ {1,...,p-1})
    # When b = a-d: is b ∈ T? b = a-d, need a-d ∈ {1,...,p-1}\\{d},
    #   a-d ≠ 0 since a ≠ d (a ∈ T), a-d ≠ d iff a ≠ 2d.

    # Count of "collision" pairs:
    # For each a ∈ T, at most 2 values of b ∈ T give a collision: b=a+d and b=a-d.
    # Each collision pair is counted twice (once for (a,b) and once for (b,a)).
    # Number of collision pairs (unordered): at most p-2 (one per a, but each shared).
    # Actually, for b=a+d: the pair (a, a+d) is the same as (a+d, a+d+d) if considering
    # the collision b=a+d, so they form chains. Each collision pair is counted once for
    # each direction.

    # For the variance of Q(d), the key quantity is:
    # Var[Q] = Σ_a Var[P_a] + Σ_{a≠b} Cov[P_a, P_b]
    #        = (p-2) var_P + Σ_{a≠b} (E[P_a P_b] - E_P²)

    # Split Σ_{a≠b}:
    # Generic pairs (no collision): E[P_a P_b] = q4, contribution per pair: q4 - q2²
    # Collision pairs: E[P_a P_b] = q3, contribution per pair: q3 - q2²

    # Number of collision pairs: Need to count ORDERED pairs (a,b) with a≠b, both in T,
    # and (a=b-d OR a-d=b).

    # a = b-d (i.e., b = a+d mod p):
    # For each a ∈ T = {1,...,p-1}\\{d}, is b = a+d mod p also in T?
    # b = a+d mod p ∈ {1,...,p-1}? Yes, since a ≠ 0 and d ≠ 0, a+d = 0 iff a = -d = p-d.
    # So b = 0 when a = p-d. Is p-d ∈ T? T = {1,...,p-1}\\{d}, and p-d ≠ d (since p odd).
    # So p-d ∈ T, meaning a = p-d gives b = 0 ∉ T. Exclude this.
    # Also need b ≠ d: a+d ≡ d iff a ≡ 0 (impossible since a ∈ T ⊂ {1,...,p-1}).
    # So for a ∈ T, b = a+d is in T iff a ≠ p-d. That's (p-2) - 1 = p-3 values.

    # Similarly for a-d = b: same count by symmetry, p-3 ordered pairs.
    # But some pairs might be counted twice if a=b-d AND a-d=b, which requires 2d=0 (excluded).
    # So total collision ordered pairs: 2(p-3).

    num_ordered_pairs = (p - 2) * (p - 3)  # total ordered pairs in T
    num_collision_pairs = 2 * (p - 3)       # ordered collision pairs
    num_generic_pairs = num_ordered_pairs - num_collision_pairs

    sum_cov_PP = (num_generic_pairs * (q4 - q2 * q2)
                  + num_collision_pairs * (q3 - q2 * q2))

    var_Q = (p - 2) * var_P + sum_cov_PP

    # Now: Var[B] = Var[Y_α] + Var[Y_β] + Var[Q]
    #             + 2Cov[Y_α, Y_β] + 2Cov[Y_α, Q] + 2Cov[Y_β, Q]

    # Cov[Y_α, Y_β] = cov_YY (α = p-d, β = d, distinct)
    # Note: if α = β, then d = p-d, so p = 2d, impossible for p odd. So always distinct.

    # Cov[Y_α, Q] = Σ_{a∈T} Cov[Y_α, P_a]
    # = Σ_{a∈T} (E[Y_α P_a] - q1 q2)
    # P_a = Y_a Y_{(a-d)}
    # E[Y_α P_a] = E[Y_α Y_a Y_{(a-d)}]

    # Case analysis for α = p-d:
    # Indices in E[Y_{p-d} Y_a Y_{(a-d)}]:
    # - If p-d, a, (a-d) are all distinct: E = q3
    # - If p-d = a: then a = p-d ∈ T (yes since d ≠ p-d), gives E[Y_{p-d}² Y_{p-2d}] = E[Y_{p-d} Y_{p-2d}] = q2
    #   (2 distinct indices: p-d and p-2d; these are equal iff d=0, excluded)
    # - If p-d = (a-d): then a = p, but a ∈ {1,...,p-1}, and p mod p = 0, excluded.
    #   So this never happens.
    # So: 1 term with E = q2 (when a = p-d), rest have E = q3.
    # Number of generic terms: (p-2) - 1 = p-3.
    # Cov[Y_α, Q] = 1 × (q2 - q1 q2) + (p-3) × (q3 - q1 q2)

    cov_Ya_Q = (q2 - q1 * q2) + (p - 3) * (q3 - q1 * q2)

    # Similarly for Cov[Y_β, Q] with β = d:
    # E[Y_d P_a] = E[Y_d Y_a Y_{(a-d)}] for a ∈ T = {1,...,p-1}\\{d}
    # - If d, a, (a-d) all distinct: E = q3
    # - If d = a: excluded since a ∈ T = {1,...,p-1}\\{d}
    # - If d = (a-d): then a = 2d mod p ∈ T iff 2d ≠ d (always since d≠0, p>2)
    #   When a = 2d: indices are {d, 2d, d} → {d, 2d}, 2 distinct (since d≠0 ⟹ 2d≠d for p>2)
    #   E = E[Y_d · Y_{2d}] = q2
    # So: 1 term with E = q2 (when a = 2d mod p), rest have E = q3.

    cov_Yb_Q = (q2 - q1 * q2) + (p - 3) * (q3 - q1 * q2)

    # Total variance
    var_B = (2 * var_Y + var_Q
             + 2 * cov_YY
             + 2 * cov_Ya_Q
             + 2 * cov_Yb_Q)

    return var_B


def main():
    print("=" * 90)
    print("EXACT MOMENT COMPUTATIONS FOR FIRST MOMENT PROOF")
    print("=" * 90)

    print(f"\n{'p':>4s} {'n':>4s} {'E[B]':>10s} {'(p-3)/4':>10s} "
          f"{'Var[B]':>12s} {'Std[B]':>10s} "
          f"{'thresh':>8s} {'gap':>6s} {'gap/std':>8s}")
    print("  " + "-" * 80)

    for p in [7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83, 97, 103, 127]:
        n = (p + 1) // 2
        thresh = (p - 3) // 2  # binding threshold for A+B at D11 positions

        E_B = exact_E_B(p)
        E_B_simple = exact_E_B_simplified(p)
        assert E_B == E_B_simple, f"Mismatch at p={p}: {E_B} ≠ {E_B_simple}"

        var_B = exact_Var_B(p)

        # E[A(d)] for random symmetric D11 with |D11| = n = (p+1)/2
        # A(d) = Δ(D11,D11,d) for d ∈ D11 (will vary by d, but average):
        # E[A(d)] = |D11|(|D11|-1)/(p-1) = n(n-1)/(p-1) = (p+1)(p-1)/(4(p-1)) = (p+1)/4
        E_A = Fraction(p + 1, 4)

        E_AB = E_A + E_B  # = (p+1)/4 + (p-3)/4 = (2p-2)/4 = (p-1)/2
        gap = Fraction(thresh) - E_AB  # should be -1

        std_B_float = float(var_B) ** 0.5

        print(f"  {p:4d} {n:4d} {float(E_B):10.4f} {(p-3)/4:10.4f} "
              f"{float(var_B):12.4f} {std_B_float:10.4f} "
              f"{thresh:8d} {float(gap):6.2f} {float(gap)/std_B_float if std_B_float > 0 else 0:8.4f}")

    # Verify against MC for small primes
    print(f"\n{'='*90}")
    print("MC VERIFICATION")
    print(f"{'='*90}")

    for p in [7, 11, 19, 23]:
        n = (p + 1) // 2
        E_B_exact = float(exact_E_B(p))
        var_B_exact = float(exact_Var_B(p))

        mc_mean, mc_var = mc_E_B(p, num_trials=500000)

        print(f"\n  p={p}:")
        print(f"    E[B]: exact={E_B_exact:.6f}, MC={mc_mean:.6f}, "
              f"diff={abs(E_B_exact-mc_mean):.6f}")
        print(f"    Var[B]: exact={var_B_exact:.6f}, MC={mc_var:.6f}, "
              f"diff={abs(var_B_exact-mc_var):.6f}")

    # The key quantity: Pr[A(d)+B(d) ≤ thresh] for a single d
    # Under Gaussian approximation: this is Φ(gap/std)
    # where gap = thresh - E[A+B] and std = √(Var[A] + Var[B])
    print(f"\n{'='*90}")
    print("GAUSSIAN APPROXIMATION FOR PER-CONSTRAINT RATE")
    print(f"{'='*90}")

    from scipy.stats import norm

    print(f"\n  {'p':>4s} {'n':>4s} {'E[A+B]':>10s} {'thresh':>8s} "
          f"{'std(B)':>8s} {'Φ(gap/σ)':>10s} {'MC rate':>10s}")
    print("  " + "-" * 60)

    for p in [11, 19, 23, 31, 43, 47, 59]:
        n = (p + 1) // 2
        thresh = (p - 3) // 2

        E_B_f = float(exact_E_B(p))
        E_A_f = (p + 1) / 4
        E_AB = E_A_f + E_B_f
        gap = thresh - E_AB

        var_B_f = float(exact_Var_B(p))
        std_B = var_B_f ** 0.5

        # Note: Var[A] is harder to compute exactly because it involves
        # the structure of random symmetric D11. We approximate it.
        # For random D11 of size n = (p+1)/2 from (p-1)/2 pairs:
        # A(d) = #{(a,b) ∈ D11²: a-b=d}, which depends on the pair structure.
        # For now, use std_AB ≈ std_B as a rough estimate (since A and B have
        # similar distributions and are independent).
        std_AB = std_B * (2 ** 0.5)  # rough: std(A+B) ≈ √2 × std(B)

        phi = norm.cdf(gap / std_AB)

        # MC validation for small p
        if p <= 31:
            rng = np.random.default_rng(42)
            s = (p - 1) // 2
            k = s - 1
            hits = 0
            M = 200000
            for _ in range(M):
                # Random symmetric D11
                pairs = [(x, p-x) for x in range(1, (p+1)//2)]
                chosen = rng.choice(len(pairs), size=(p+1)//4, replace=False)
                D11 = set()
                for i in chosen:
                    D11.add(pairs[i][0])
                    D11.add(pairs[i][1])

                # Random D12
                others = rng.choice(range(1, p), size=k, replace=False)
                D12 = set([0]) | set(others.tolist())

                # Check one constraint: d=1
                d = 1
                A_d = sum(1 for a in D11 for b in D11 if (a-b) % p == d)
                B_d = sum(1 for a in D12 for b in D12 if (a-b) % p == d)
                if A_d + B_d <= thresh:
                    hits += 1
            mc_rate = hits / M
        else:
            mc_rate = float('nan')

        mc_str = f"{mc_rate:.4f}" if not np.isnan(mc_rate) else "  -"
        print(f"  {p:4d} {n:4d} {E_AB:10.4f} {thresh:8d} "
              f"{std_B:8.4f} {phi:10.4f} {mc_str:>10s}")

    # Asymptotic analysis
    print(f"\n{'='*90}")
    print("ASYMPTOTIC ANALYSIS")
    print(f"{'='*90}")
    print(f"\n  For large p:")
    print(f"    E[A(d)] = (p+1)/4")
    print(f"    E[B(d)] = (p-3)/4")
    print(f"    E[A+B] = (p-1)/2")
    print(f"    Threshold = (p-3)/2 = E[A+B] - 1")
    print(f"    Var[B(d)] = ?")

    # Compute Var[B] asymptotic
    # From the exact formula, for large p:
    # k = (p-3)/2, N = p-1
    # q1 = (p-3)/(2(p-1)) → 1/2
    # q2 = (p-3)(p-5)/(4(p-1)(p-2)) → 1/4
    # q3 → 1/8
    # q4 → 1/16
    # var_Y = q1(1-q1) → 1/4
    # cov_YY = q2 - q1² → 1/4 - 1/4 = 0 (leading order)

    for p in [31, 59, 83, 127, 199, 499, 997]:
        var_B = float(exact_Var_B(p))
        n = (p + 1) // 2
        print(f"    p={p:4d}: Var[B] = {var_B:.4f}, "
              f"Var[B]/p = {var_B/p:.6f}, "
              f"std = {var_B**0.5:.4f}")

    print(f"\n  Key finding: Var[B(d)] ≈ c·p for some constant c ≈ ?")
    print(f"  If Var[A+B] = O(p), then gap/std = -1/√(O(p)) = -O(1/√p)")
    print(f"  So Pr[single constraint ok] → 1/2 as p → ∞")
    print(f"  And Pr[all (p+1)/2 constraints ok] → (1/2)^{{(p+1)/2}}")
    print(f"  With #pairs ≈ 2^{{3(p-1)/2}}, E = 2^{{3(p-1)/2 - (p+1)/2}} = 2^{{p-2}} → ∞!")
    print(f"\n  THIS IS THE PROOF (modulo rigorous joint probability bound).")


if __name__ == '__main__':
    main()
