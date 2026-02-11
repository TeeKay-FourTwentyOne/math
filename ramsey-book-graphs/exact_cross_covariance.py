#!/usr/bin/env python3
"""Exact cross-covariance Cov[B(d1), B(d2)] for d1 != d2.

Setup (same as exact_moments.py):
  D12 = {0} ∪ S where S is a random k-subset of {1,...,p-1}
  k = (p-3)/2, N = p-1

  B(d) = Y_{p-d} + Y_d + Σ_{a ∈ T(d)} Y_a · Y_{(a-d) mod p}
  where T(d) = {1,...,p-1} \ {d}, and Y_a = 1[a ∈ S].

We compute Cov[B(d1), B(d2)] = E[B(d1)·B(d2)] - E[B(d1)]·E[B(d2)].

The product B(d1)·B(d2) decomposes into 9 types of cross-terms from
the 3+3 components of each B(d). For each term we count index collisions
to determine the correct hypergeometric moment (q1, q2, q3, or q4).
"""

import numpy as np
from fractions import Fraction
import sys
import time


def exact_cross_covariance(p, d1, d2):
    """Compute exact Cov[B(d1), B(d2)] for d1 != d2 (both nonzero, mod p).

    Returns (cov, E_B1B2, E_B1, E_B2) all as Fractions.
    """
    assert d1 != d2 and d1 != 0 and d2 != 0
    d1 = d1 % p
    d2 = d2 % p
    assert d1 != d2

    s = (p - 1) // 2
    k = s - 1  # = (p-3)/2
    N = p - 1

    # Hypergeometric moments for a random k-subset of {1,...,N}
    q1 = Fraction(k, N)
    q2 = Fraction(k * (k-1), N * (N-1))
    q3 = Fraction(k * (k-1) * (k-2), N * (N-1) * (N-2))
    q4 = Fraction(k * (k-1) * (k-2) * (k-3), N * (N-1) * (N-2) * (N-3))

    # E[B(d)] = (p-3)/4 for any d != 0
    E_B = Fraction(p - 3, 4)

    # B(d1) = L1 + L2 + Q1
    # B(d2) = M1 + M2 + Q2
    # where:
    #   L1 = Y_{p-d1}, L2 = Y_{d1}
    #   M1 = Y_{p-d2}, M2 = Y_{d2}
    #   Q1 = Σ_{a ∈ T1} Y_a · Y_{(a-d1)},  T1 = {1,...,p-1}\{d1}
    #   Q2 = Σ_{b ∈ T2} Y_b · Y_{(b-d2)},  T2 = {1,...,p-1}\{d2}
    #
    # B(d1)·B(d2) has 9 cross-terms:
    # (1) L1·M1 = Y_{p-d1} · Y_{p-d2}
    # (2) L1·M2 = Y_{p-d1} · Y_{d2}
    # (3) L2·M1 = Y_{d1} · Y_{p-d2}
    # (4) L2·M2 = Y_{d1} · Y_{d2}
    # (5) L1·Q2 = Y_{p-d1} · Σ_b Y_b · Y_{(b-d2)}
    # (6) L2·Q2 = Y_{d1} · Σ_b Y_b · Y_{(b-d2)}
    # (7) Q1·M1 = (Σ_a Y_a · Y_{(a-d1)}) · Y_{p-d2}
    # (8) Q1·M2 = (Σ_a Y_a · Y_{(a-d1)}) · Y_{d2}
    # (9) Q1·Q2 = (Σ_a Y_a · Y_{(a-d1)}) · (Σ_b Y_b · Y_{(b-d2)})

    # Helper: compute E[product of Y-indicators] given a multiset of indices
    # by determining the number of DISTINCT indices
    def E_Y_product(indices):
        """E[Y_{i1} · Y_{i2} · ...] where indices is a list of elements in {1..p-1}.
        Since Y_a^2 = Y_a, we just need the set of distinct indices."""
        distinct = len(set(indices))
        if distinct == 0:
            return Fraction(1)
        elif distinct == 1:
            return q1
        elif distinct == 2:
            return q2
        elif distinct == 3:
            return q3
        elif distinct == 4:
            return q4
        else:
            # For k >= 5 distinct indices: k(k-1)...(k-d+1) / (N(N-1)...(N-d+1))
            num = Fraction(1)
            for i in range(distinct):
                num *= Fraction(k - i, N - i)
            return num

    # Precompute the 4 "special" indices
    alpha1 = (p - d1) % p  # = p - d1
    beta1 = d1
    alpha2 = (p - d2) % p  # = p - d2
    beta2 = d2

    # =========================================================
    # Terms (1)-(4): products of two standalone Y's
    # =========================================================
    # (1) E[Y_{alpha1} · Y_{alpha2}]
    E_1 = E_Y_product([alpha1, alpha2])
    # (2) E[Y_{alpha1} · Y_{beta2}]
    E_2 = E_Y_product([alpha1, beta2])
    # (3) E[Y_{beta1} · Y_{alpha2}]
    E_3 = E_Y_product([beta1, alpha2])
    # (4) E[Y_{beta1} · Y_{beta2}]
    E_4 = E_Y_product([beta1, beta2])

    # =========================================================
    # Terms (5)-(8): standalone Y times a Q sum
    # =========================================================
    # (5) E[L1 · Q2] = Σ_{b ∈ T2} E[Y_{alpha1} · Y_b · Y_{(b-d2)}]
    #     T2 = {1,...,p-1}\{d2}
    E_5 = Fraction(0)
    for b in range(1, p):
        if b == d2:
            continue
        bd2 = (b - d2) % p
        E_5 += E_Y_product([alpha1, b, bd2])

    # (6) E[L2 · Q2] = Σ_{b ∈ T2} E[Y_{beta1} · Y_b · Y_{(b-d2)}]
    E_6 = Fraction(0)
    for b in range(1, p):
        if b == d2:
            continue
        bd2 = (b - d2) % p
        E_6 += E_Y_product([beta1, b, bd2])

    # (7) E[Q1 · M1] = Σ_{a ∈ T1} E[Y_a · Y_{(a-d1)} · Y_{alpha2}]
    E_7 = Fraction(0)
    for a in range(1, p):
        if a == d1:
            continue
        ad1 = (a - d1) % p
        E_7 += E_Y_product([a, ad1, alpha2])

    # (8) E[Q1 · M2] = Σ_{a ∈ T1} E[Y_a · Y_{(a-d1)} · Y_{beta2}]
    E_8 = Fraction(0)
    for a in range(1, p):
        if a == d1:
            continue
        ad1 = (a - d1) % p
        E_8 += E_Y_product([a, ad1, beta2])

    # =========================================================
    # Term (9): Q1 · Q2 (the big double sum)
    # =========================================================
    # E[Q1 · Q2] = Σ_{a ∈ T1} Σ_{b ∈ T2} E[Y_a · Y_{(a-d1)} · Y_b · Y_{(b-d2)}]
    E_9 = Fraction(0)
    for a in range(1, p):
        if a == d1:
            continue
        ad1 = (a - d1) % p
        for b in range(1, p):
            if b == d2:
                continue
            bd2 = (b - d2) % p
            E_9 += E_Y_product([a, ad1, b, bd2])

    # E[B(d1) · B(d2)] = sum of all 9 terms
    E_B1B2 = E_1 + E_2 + E_3 + E_4 + E_5 + E_6 + E_7 + E_8 + E_9

    cov = E_B1B2 - E_B * E_B

    return cov, E_B1B2, E_B, E_B


def exact_cross_covariance_fast(p, d1, d2):
    """Faster version that counts collision types analytically for the Q1*Q2 term.

    For large p, the brute-force double loop in term (9) is O(p^2).
    This version counts the number of (a,b) pairs with each collision pattern.
    """
    assert d1 != d2 and d1 != 0 and d2 != 0
    d1 = d1 % p
    d2 = d2 % p
    assert d1 != d2

    s = (p - 1) // 2
    k = s - 1
    N = p - 1

    q1 = Fraction(k, N)
    q2 = Fraction(k * (k-1), N * (N-1))
    q3 = Fraction(k * (k-1) * (k-2), N * (N-1) * (N-2))
    q4 = Fraction(k * (k-1) * (k-2) * (k-3), N * (N-1) * (N-2) * (N-3))

    E_B = Fraction(p - 3, 4)

    alpha1 = (p - d1) % p
    beta1 = d1
    alpha2 = (p - d2) % p
    beta2 = d2

    def E_Y_product(indices):
        distinct = len(set(indices))
        if distinct == 0:
            return Fraction(1)
        elif distinct == 1:
            return q1
        elif distinct == 2:
            return q2
        elif distinct == 3:
            return q3
        elif distinct == 4:
            return q4
        else:
            num = Fraction(1)
            for i in range(distinct):
                num *= Fraction(k - i, N - i)
            return num

    # Terms (1)-(4): same as brute force
    E_1 = E_Y_product([alpha1, alpha2])
    E_2 = E_Y_product([alpha1, beta2])
    E_3 = E_Y_product([beta1, alpha2])
    E_4 = E_Y_product([beta1, beta2])

    # Terms (5)-(8): standalone Y times Q sum
    # These are O(p) sums, manageable even for large p
    E_5 = Fraction(0)
    for b in range(1, p):
        if b == d2:
            continue
        E_5 += E_Y_product([alpha1, b, (b - d2) % p])

    E_6 = Fraction(0)
    for b in range(1, p):
        if b == d2:
            continue
        E_6 += E_Y_product([beta1, b, (b - d2) % p])

    E_7 = Fraction(0)
    for a in range(1, p):
        if a == d1:
            continue
        E_7 += E_Y_product([a, (a - d1) % p, alpha2])

    E_8 = Fraction(0)
    for a in range(1, p):
        if a == d1:
            continue
        E_8 += E_Y_product([a, (a - d1) % p, beta2])

    # Term (9): Q1 · Q2 -- analytical collision counting
    # For a ∈ T1 = {1,...,p-1}\{d1}, b ∈ T2 = {1,...,p-1}\{d2}:
    # Indices: {a, (a-d1) mod p, b, (b-d2) mod p}
    # Note: a != d1, b != d2, and a-d1 != 0 (since a != d1), b-d2 != 0 (since b != d2).
    # All four indices are in {1,...,p-1}.
    #
    # We need to count pairs (a,b) with various collision patterns.
    # Let ad1 = (a-d1) mod p, bd2 = (b-d2) mod p.
    #
    # Possible equalities among {a, ad1, b, bd2}:
    # (Note: a != ad1 since d1 != 0, b != bd2 since d2 != 0)
    # The possible cross-equalities are:
    #   (i)   a = b       iff a = b
    #   (ii)  a = bd2     iff a = b - d2  iff b = a + d2
    #   (iii) ad1 = b     iff a - d1 = b  iff b = a - d1
    #   (iv)  ad1 = bd2   iff a - d1 = b - d2  iff b = a - d1 + d2
    #
    # When d1 != d2, these four conditions define four different values of b for each a.
    # Two conditions can coincide:
    #   (i)=(ii): a = b and a = b - d2 => d2 = 0, impossible
    #   (i)=(iii): a = b and b = a - d1 => d1 = 0, impossible
    #   (i)=(iv): a = b and b = a - d1 + d2 => d1 = d2, impossible
    #   (ii)=(iii): b = a + d2 and b = a - d1 => d2 = -d1 => d1 + d2 = 0 mod p
    #   (ii)=(iv): b = a + d2 and b = a - d1 + d2 => d1 = 0, impossible
    #   (iii)=(iv): b = a - d1 and b = a - d1 + d2 => d2 = 0, impossible
    #
    # So the only possible double-collision is (ii)+(iii) when d1 + d2 ≡ 0 mod p.
    #
    # Triple/quadruple collisions would require additional equalities among the
    # already-collapsed set, which we handle below.

    # Define delta values that determine which b-values give collisions
    # For each a, collisions occur at:
    delta_i = 0          # b = a (collision a = b)
    delta_ii = d2        # b = a + d2 (collision a = bd2)
    delta_iii = (-d1) % p  # b = a - d1 (collision ad1 = b), i.e., b = a + (p-d1)
    delta_iv = (d2 - d1) % p  # b = a - d1 + d2 (collision ad1 = bd2)

    # These delta values (offsets for b relative to a) may coincide:
    deltas = {
        'i': delta_i,
        'ii': delta_ii,
        'iii': delta_iii,
        'iv': delta_iv,
    }

    # Group deltas by value to find coinciding collision types
    from collections import defaultdict
    delta_groups = defaultdict(list)
    for name, val in deltas.items():
        delta_groups[val].append(name)

    # For each group of coinciding deltas, count how many distinct indices
    # result when those collisions hold simultaneously.
    # We also need to determine: for each a, is b = a + delta valid?
    # Valid means b ∈ T2 = {1,...,p-1}\{d2}, and also a ∈ T1 = {1,...,p-1}\{d1}.

    # Strategy: for each distinct delta value, iterate over a to count valid pairs
    # and determine the number of distinct indices.
    # But actually, the number of distinct indices for a given set of collisions
    # is the same for all a (generically), except when a itself hits a special value.
    # Let's just enumerate.

    # For the GENERIC case (no collisions, 4 distinct indices):
    # Total pairs: |T1| * |T2| = (p-2)^2
    # Minus collision pairs.

    # Let me just carefully count by iterating over distinct delta values.

    # For each delta d (offset), count:
    #   C(d) = #{a ∈ T1 : b = (a + d) mod p ∈ T2}
    # where T1 = {1,...,p-1}\{d1}, T2 = {1,...,p-1}\{d2}.
    #
    # b = (a + d) mod p.
    # We need a ∈ {1,...,p-1}\{d1} and b ∈ {1,...,p-1}\{d2}.
    # b = 0 iff a = -d mod p = (p-d) mod p.
    # b = d2 iff a = (d2 - d) mod p.
    # a = d1 excluded.
    # a = 0 not in T1.
    #
    # Start with all a ∈ {1,...,p-1}\{d1}: that's p-2 values.
    # Remove a where b = 0: a = (p - d) mod p -- remove if this a is in T1
    #   i.e., if (p - d) mod p ∈ {1,...,p-1}\{d1}
    #   (p-d) mod p = 0 iff d = 0; if d != 0, it's in {1,...,p-1}
    #   So remove iff (p - d) mod p != d1, i.e., d != p - d1
    # Remove a where b = d2: a = (d2 - d) mod p -- remove if this a is in T1
    #   (d2 - d) mod p = 0 iff d = d2
    #   If d != d2 and (d2 - d) mod p != d1: remove 1 more
    #   If (d2 - d) mod p = d1: a = d1 already excluded, no extra removal

    def count_valid_pairs_with_offset(delta):
        """Count #{a in T1 : (a + delta) mod p in T2}."""
        cnt = 0
        for a in range(1, p):
            if a == d1:
                continue
            b = (a + delta) % p
            if b == 0 or b == d2:
                continue
            cnt += 1
        return cnt

    def count_and_classify_offset(delta, collision_names):
        """For pairs with b = a + delta, count pairs and determine
        the number of distinct indices {a, a-d1, b, b-d2}.

        With collisions in collision_names, some of these are equal.
        But the exact count may vary with a (due to additional coincidences
        like a = a-d1, which can't happen since d1 != 0, or a-d1 = b-d2
        additionally, etc.).

        Actually, the collisions from the delta already encode which of
        {a, ad1, b, bd2} are equal. But there could be EXTRA collisions
        (e.g., a might equal ad1, but that requires d1=0, impossible).
        Wait -- the only possible equalities are:
          a = ad1: d1 = 0, impossible
          b = bd2: d2 = 0, impossible
          a = b, a = bd2, ad1 = b, ad1 = bd2
        These are exactly our four collision types. So the delta determines
        ALL collisions. No extra ones are possible.

        However, we need to verify: a = ad1 can't happen, but what about
        a = b AND ad1 = bd2 simultaneously? That's (i) AND (iv):
          a = b => b = a
          ad1 = bd2 => a - d1 = a - d2 => d1 = d2, impossible.
        What about (i) AND (ii): a = b AND a = b - d2 => d2 = 0, impossible.
        What about (ii) AND (iii): a = bd2 AND ad1 = b
          => b = a + d2 and b = a - d1 => d1 + d2 = 0.
          In this case, delta_ii = d2 = delta_iii = p - d1, and these coincide.
          Distinct indices: {a, a-d1, a+d2, a} but d2 = -d1, so a+d2 = a-d1.
          So {a, a-d1} -- just 2 distinct indices.
        What about (ii) AND (iv): b = a+d2 AND b = a-d1+d2 => d1 = 0, impossible.
        """
        # The set of collision types that hold for this delta determines
        # which of {a, ad1, b, bd2} are identified.
        # Let's figure out the equivalence classes.

        # Start with 4 elements: a, ad1, b, bd2
        # Relations from collisions:
        uf = {0: 0, 1: 1, 2: 2, 3: 3}  # 0=a, 1=ad1, 2=b, 3=bd2

        def find(x):
            while uf[x] != x:
                x = uf[x]
            return x

        def union(x, y):
            rx, ry = find(x), find(y)
            if rx != ry:
                uf[rx] = ry

        for cn in collision_names:
            if cn == 'i':    # a = b
                union(0, 2)
            elif cn == 'ii':  # a = bd2
                union(0, 3)
            elif cn == 'iii': # ad1 = b
                union(1, 2)
            elif cn == 'iv':  # ad1 = bd2
                union(1, 3)

        num_distinct = len(set(find(x) for x in range(4)))
        return num_distinct

    # Build the contribution from each distinct delta offset
    # For offsets that DON'T correspond to any collision type: 4 distinct indices -> q4
    # For offsets that DO correspond to collision type(s): fewer distinct indices

    # Collect all collision delta values (may have repeats if deltas coincide)
    collision_deltas = set(deltas.values())

    # For the Q1*Q2 double sum, each (a, b) pair with a in T1, b in T2 contributes
    # E_Y_product({a, ad1, b, bd2}). The number of distinct indices depends on
    # whether b - a (mod p) equals one of the collision deltas.

    E_9 = Fraction(0)

    total_collision_count = 0
    for delta_val, coll_names in delta_groups.items():
        nd = count_and_classify_offset(delta_val, coll_names)
        cnt = count_valid_pairs_with_offset(delta_val)
        total_collision_count += cnt
        if nd == 1:
            E_9 += cnt * q1
        elif nd == 2:
            E_9 += cnt * q2
        elif nd == 3:
            E_9 += cnt * q3
        elif nd == 4:
            E_9 += cnt * q4

    # How many "generic" pairs are there (b - a not equal to any collision delta)?
    # Total valid pairs = |T1| * |T2| minus pairs where b = 0 or b = d2 for the
    # respective a values... Actually, total valid pairs:
    total_valid = 0
    for a in range(1, p):
        if a == d1:
            continue
        for b in range(1, p):
            if b == d2:
                continue
            total_valid += 1
    # = (p-2) * (p-2)  [since |T1| = |T2| = p-2... wait, T1 = {1..p-1}\{d1}, |T1| = p-2
    #  T2 = {1..p-1}\{d2}, |T2| = p-2. But these are independent loops, so total = (p-2)^2]
    assert total_valid == (p - 2) * (p - 2)

    # But from the collision deltas, we've already counted some.
    # The remaining pairs have all 4 indices distinct -> q4.
    # But wait, a collision delta corresponds to specific (a,b) pairs where b = a + delta.
    # For each delta_val in delta_groups, we counted count_valid_pairs_with_offset(delta_val).
    # But some delta_vals are NOT collision deltas. We need to add the generic ones.

    # Actually let's be more careful. The set of ALL delta values (b - a mod p) that can occur:
    # For a in T1, b in T2, b - a mod p can be anything in Z_p.
    # We've classified those deltas that appear in collision_deltas as "collision offsets".
    # All OTHER offsets give 4 distinct indices.

    # Count of collision-offset pairs:
    collision_pair_count = 0
    for delta_val in delta_groups:
        collision_pair_count += count_valid_pairs_with_offset(delta_val)

    generic_pair_count = total_valid - collision_pair_count
    E_9 += generic_pair_count * q4

    # Now we have all 9 terms
    E_B1B2 = E_1 + E_2 + E_3 + E_4 + E_5 + E_6 + E_7 + E_8 + E_9
    cov = E_B1B2 - E_B * E_B

    return cov, E_B1B2, E_B, E_B


def mc_cross_moment(p, d1, d2, num_trials=500000):
    """Monte Carlo estimate of E[B(d1)·B(d2)] and individual E[B(d)]."""
    s = (p - 1) // 2
    k = s - 1
    rng = np.random.default_rng(42)

    prod_samples = []
    B1_samples = []
    B2_samples = []

    for _ in range(num_trials):
        S = set(rng.choice(range(1, p), size=k, replace=False).tolist())
        D12 = S | {0}

        B1 = sum(1 for a in D12 if (a - d1) % p in D12)
        B2 = sum(1 for a in D12 if (a - d2) % p in D12)

        B1_samples.append(B1)
        B2_samples.append(B2)
        prod_samples.append(B1 * B2)

    E_B1 = np.mean(B1_samples)
    E_B2 = np.mean(B2_samples)
    E_prod = np.mean(prod_samples)
    cov_mc = E_prod - E_B1 * E_B2
    var_B1 = np.var(B1_samples)
    var_B2 = np.var(B2_samples)

    return E_prod, E_B1, E_B2, cov_mc, var_B1, var_B2


def mc_cross_moment_fft(p, d1, d2, num_trials=500000):
    """FFT-based MC for speed."""
    s = (p - 1) // 2
    k = s - 1
    rng = np.random.default_rng(42)

    prods = np.zeros(num_trials)
    B1s = np.zeros(num_trials)
    B2s = np.zeros(num_trials)

    for i in range(num_trials):
        S = rng.choice(range(1, p), size=k, replace=False)
        d12_ind = np.zeros(p)
        d12_ind[0] = 1
        d12_ind[S] = 1
        autocorr = np.fft.ifft(np.abs(np.fft.fft(d12_ind))**2).real
        b1 = round(autocorr[d1])
        b2 = round(autocorr[d2])
        B1s[i] = b1
        B2s[i] = b2
        prods[i] = b1 * b2

    E_B1 = np.mean(B1s)
    E_B2 = np.mean(B2s)
    E_prod = np.mean(prods)
    cov_mc = E_prod - E_B1 * E_B2
    var_B1 = np.var(B1s)
    var_B2 = np.var(B2s)

    return E_prod, E_B1, E_B2, cov_mc, var_B1, var_B2


def exact_Var_B(p):
    """Copy of the Var[B] computation from exact_moments.py for comparison."""
    s = (p - 1) // 2
    k = s - 1
    N = p - 1

    q1 = Fraction(k, N)
    q2 = Fraction(k * (k-1), N * (N-1))
    q3 = Fraction(k * (k-1) * (k-2), N * (N-1) * (N-2))
    q4 = Fraction(k * (k-1) * (k-2) * (k-3), N * (N-1) * (N-2) * (N-3))

    var_Y = q1 * (1 - q1)
    cov_YY = q2 - q1 * q1

    E_P = q2
    E_P2 = q2
    var_P = q2 - q2 * q2

    num_ordered_pairs = (p - 2) * (p - 3)
    num_collision_pairs = 2 * (p - 3)
    num_generic_pairs = num_ordered_pairs - num_collision_pairs

    sum_cov_PP = (num_generic_pairs * (q4 - q2 * q2)
                  + num_collision_pairs * (q3 - q2 * q2))

    var_Q = (p - 2) * var_P + sum_cov_PP

    cov_Ya_Q = (q2 - q1 * q2) + (p - 3) * (q3 - q1 * q2)
    cov_Yb_Q = (q2 - q1 * q2) + (p - 3) * (q3 - q1 * q2)

    var_B = (2 * var_Y + var_Q
             + 2 * cov_YY
             + 2 * cov_Ya_Q
             + 2 * cov_Yb_Q)

    return var_B


def main():
    print("=" * 90)
    print("EXACT CROSS-COVARIANCE Cov[B(d1), B(d2)] COMPUTATION")
    print("=" * 90)

    # =========================================================
    # 1. Verify brute-force vs fast for small p
    # =========================================================
    print("\n--- Verification: brute-force vs fast method ---")
    for p in [7, 11]:
        pairs_to_test = []
        # Collect several (d1, d2) pairs including special relationships
        for d1 in range(1, min(p, 6)):
            for d2 in range(d1 + 1, min(p, 6)):
                pairs_to_test.append((d1, d2))

        for d1, d2 in pairs_to_test:
            cov_bf, E_prod_bf, _, _ = exact_cross_covariance(p, d1, d2)
            cov_fast, E_prod_fast, _, _ = exact_cross_covariance_fast(p, d1, d2)
            match = "OK" if cov_bf == cov_fast else "MISMATCH"
            print(f"  p={p}, d1={d1}, d2={d2}: brute={float(cov_bf):.6f}, "
                  f"fast={float(cov_fast):.6f}  [{match}]")
            assert cov_bf == cov_fast, f"Mismatch at p={p}, d1={d1}, d2={d2}"
    print("  All brute-force vs fast checks PASSED.")

    # =========================================================
    # 2. Exact values and MC verification
    # =========================================================
    print("\n--- Exact Cov[B(d1),B(d2)] and MC verification ---")
    print(f"{'p':>4s} {'d1':>4s} {'d2':>4s} {'rel':>12s} "
          f"{'Cov(exact)':>14s} {'Cov(MC)':>12s} {'Var(B)':>12s} "
          f"{'Cov/Var':>10s}")
    print("  " + "-" * 85)

    for p in [7, 11, 19, 23]:
        var_B = exact_Var_B(p)
        var_B_f = float(var_B)

        # Generate pairs with various relationships
        test_pairs = []
        # Generic pair
        test_pairs.append((1, 2, "generic"))
        if p > 7:
            test_pairs.append((1, 3, "generic"))
            test_pairs.append((2, 5, "generic"))

        # d1 + d2 = p (i.e., d2 = p - d1)
        test_pairs.append((1, p - 1, "d1+d2=p"))
        if p > 11:
            test_pairs.append((2, p - 2, "d1+d2=p"))

        # d2 = 2*d1 mod p
        d2_double = (2 * 1) % p
        if d2_double != 1 and d2_double != 0:
            test_pairs.append((1, d2_double, "d2=2d1"))

        # d1 = 2*d2 mod p
        d1_double = (2 * 2) % p
        if d1_double != 2 and d1_double != 0:
            test_pairs.append((d1_double, 2, "d1=2d2"))

        # d2 - d1 special
        if p > 11:
            test_pairs.append((3, 6, "d2=2d1"))

        for d1, d2, rel in test_pairs:
            if d1 == d2 or d1 == 0 or d2 == 0 or d1 >= p or d2 >= p:
                continue

            cov_exact, E_prod_exact, _, _ = exact_cross_covariance_fast(p, d1, d2)
            cov_exact_f = float(cov_exact)

            # MC verification
            num_trials = 300000 if p <= 11 else 200000
            E_prod_mc, E_B1_mc, E_B2_mc, cov_mc, _, _ = mc_cross_moment_fft(
                p, d1, d2, num_trials=num_trials)

            ratio = cov_exact_f / var_B_f if var_B_f > 0 else 0

            print(f"  {p:4d} {d1:4d} {d2:4d} {rel:>12s} "
                  f"{cov_exact_f:14.6f} {cov_mc:12.6f} {var_B_f:12.6f} "
                  f"{ratio:10.6f}")

    # =========================================================
    # 3. Check if Cov depends on specific (d1, d2) values
    # =========================================================
    print("\n--- Does Cov depend on (d1,d2)? (all pairs for p=11) ---")
    p = 11
    var_B = exact_Var_B(p)
    covs = {}
    for d1 in range(1, p):
        for d2 in range(1, p):
            if d1 == d2:
                continue
            cov, _, _, _ = exact_cross_covariance_fast(p, d1, d2)
            # Classify the relationship
            rel = "generic"
            if (d1 + d2) % p == 0:
                rel = "d1+d2=0"
            elif (2 * d1) % p == d2:
                rel = "d2=2d1"
            elif (2 * d2) % p == d1:
                rel = "d1=2d2"
            key = (d1, d2)
            covs[key] = (cov, rel)

    # Group by value
    from collections import Counter
    val_counts = Counter()
    rel_vals = {}
    for (d1, d2), (cov, rel) in covs.items():
        val_counts[cov] += 1
        if rel not in rel_vals:
            rel_vals[rel] = set()
        rel_vals[rel].add(cov)

    print(f"  p={p}: {len(val_counts)} distinct covariance values among {len(covs)} pairs:")
    for val, cnt in sorted(val_counts.items(), key=lambda x: -x[1]):
        print(f"    Cov = {float(val):.8f} ({val})  -- {cnt} pairs")

    print(f"\n  By relationship type:")
    for rel, vals in sorted(rel_vals.items()):
        val_strs = [f"{float(v):.6f}" for v in sorted(vals)]
        print(f"    {rel}: {val_strs}")

    # =========================================================
    # 4. Scaling with p: is |Cov| = O(1)?
    # =========================================================
    print("\n--- Scaling of |Cov[B(d1),B(d2)]| with p ---")
    print(f"{'p':>6s} {'Cov(1,2)':>14s} {'Var(B)':>14s} "
          f"{'|Cov|/Var':>12s} {'Cov/1':>10s}")
    print("  " + "-" * 60)

    cov_values = []
    var_values = []
    p_values = []

    for p in [7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83]:
        d1, d2 = 1, 2
        cov, _, _, _ = exact_cross_covariance_fast(p, d1, d2)
        var_B = exact_Var_B(p)

        cov_f = float(cov)
        var_f = float(var_B)
        ratio = abs(cov_f) / var_f if var_f > 0 else 0

        cov_values.append(cov_f)
        var_values.append(var_f)
        p_values.append(p)

        print(f"  {p:6d} {cov_f:14.6f} {var_f:14.6f} "
              f"{ratio:12.6f} {cov_f:10.6f}")

    # =========================================================
    # 5. Also check other pair types at scale
    # =========================================================
    print("\n--- Scaling by relationship type ---")
    print(f"{'p':>6s} {'generic(1,2)':>14s} {'d1+d2=0(1,p-1)':>16s} "
          f"{'d2=2d1(1,2)':>14s} {'Var(B)':>12s}")
    print("  " + "-" * 65)

    for p in [7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83]:
        var_B = float(exact_Var_B(p))

        # Generic: (1, 3) for p > 7
        d1g, d2g = 1, 3
        if d1g == d2g or d2g >= p:
            d1g, d2g = 1, 2
        cov_gen, _, _, _ = exact_cross_covariance_fast(p, d1g, d2g)

        # d1 + d2 = 0 mod p: (1, p-1)
        cov_sum0, _, _, _ = exact_cross_covariance_fast(p, 1, p - 1)

        # d2 = 2*d1: (1, 2)
        cov_dbl, _, _, _ = exact_cross_covariance_fast(p, 1, 2)

        print(f"  {p:6d} {float(cov_gen):14.6f} {float(cov_sum0):16.6f} "
              f"{float(cov_dbl):14.6f} {var_B:12.4f}")

    # =========================================================
    # 6. Asymptotic formula
    # =========================================================
    print("\n--- Asymptotic analysis ---")
    print("\nFor large p, the hypergeometric moments converge:")
    print("  q1 -> 1/2, q2 -> 1/4, q3 -> 1/8, q4 -> 1/16")
    print("\nThe dominant term in Cov comes from the Q1*Q2 double sum.")
    print("Most pairs contribute q4 - q2^2 -> 1/16 - 1/16 = 0.")
    print("Collision pairs (O(p) of them) contribute q3 - q2^2 -> 1/8 - 1/16 = 1/16.")
    print("Plus O(p) terms from the Y*Q cross terms.")
    print("Plus O(1) terms from the Y*Y terms.")
    print()

    # Let's compute the leading-order contribution analytically
    # For the Q1*Q2 term:
    #   There are 4 collision offsets (some may coincide).
    #   For each offset, roughly p - O(1) valid pairs.
    #   Each collision pair contributes q3 - q2^2 extra (beyond generic q4).
    #   So the collision contribution is ~ 4 * p * (q3 - q4) = 4p * (q3 - q4)
    #   Actually, the Cov = E[B1*B2] - E[B1]*E[B2], and the dominant terms
    #   in E[B1*B2] from Q1*Q2 are:
    #     Σ (E[Y_a Y_{ad1} Y_b Y_{bd2}] - E[Y_a Y_{ad1}] * E[Y_b Y_{bd2}])
    #   For generic pairs: q4 - q2^2
    #   For collision pairs: q3 - q2^2 (if 1 collision) or q2 - q2^2 (if 2 collisions)

    # Let's extract the exact asymptotic coefficient
    print("Extracting Cov[B(d1),B(d2)] / 1 as p -> infinity:")
    for p in [199, 499, 997, 1999]:
        cov, _, _, _ = exact_cross_covariance_fast(p, 1, 2)
        cov_f = float(cov)
        var_f = float(exact_Var_B(p))
        print(f"  p={p:5d}: Cov = {cov_f:12.6f}, Var = {var_f:12.4f}, "
              f"|Cov|/Var = {abs(cov_f)/var_f:.8f}, Cov*16/p = {cov_f*16/p:.6f}")

    # Check if Cov converges to a constant
    print("\nChecking if Cov[B(d1),B(d2)] converges to a constant:")
    for p in [499, 997, 1999, 3989]:
        cov, _, _, _ = exact_cross_covariance_fast(p, 1, 2)
        cov_f = float(cov)
        print(f"  p={p:5d}: Cov(1,2) = {cov_f:.8f}")

    # Also check (1, p-1) case
    print("\nCov[B(1), B(p-1)] (complementary pair):")
    for p in [499, 997, 1999, 3989]:
        cov, _, _, _ = exact_cross_covariance_fast(p, 1, p - 1)
        cov_f = float(cov)
        print(f"  p={p:5d}: Cov(1,p-1) = {cov_f:.8f}")

    # =========================================================
    # 7. Full covariance matrix for small p
    # =========================================================
    print("\n--- Full covariance matrix structure (p=11) ---")
    p = 11
    var_B = exact_Var_B(p)
    print(f"  Var[B(d)] = {float(var_B):.6f} for all d != 0")

    dvals = list(range(1, p))
    n = len(dvals)
    cov_matrix = np.zeros((n, n))
    for i, d1 in enumerate(dvals):
        cov_matrix[i, i] = float(var_B)
        for j, d2 in enumerate(dvals):
            if i == j:
                continue
            cov, _, _, _ = exact_cross_covariance_fast(p, d1, d2)
            cov_matrix[i, j] = float(cov)

    # Eigenvalues of the covariance matrix
    eigvals = np.linalg.eigvalsh(cov_matrix)
    print(f"  Eigenvalues of {n}x{n} covariance matrix:")
    print(f"    min = {eigvals[0]:.6f}, max = {eigvals[-1]:.6f}")
    print(f"    all positive: {all(e > -1e-10 for e in eigvals)}")
    print(f"    condition number: {eigvals[-1] / max(eigvals[0], 1e-15):.2f}")

    # Correlation matrix
    corr_matrix = cov_matrix / float(var_B)
    print(f"\n  Correlation matrix Cov(d1,d2)/Var(B):")
    print(f"    Off-diagonal values (sample):")
    for i in range(min(5, n)):
        for j in range(i+1, min(5, n)):
            print(f"      corr({dvals[i]},{dvals[j]}) = {corr_matrix[i,j]:.6f}")

    # =========================================================
    # 8. Verify closed-form formulas
    # =========================================================
    print("\n--- CLOSED-FORM VERIFICATION ---")
    print()
    print("KEY IDENTITY: B(d) = B(p-d) always (autocorrelation symmetry).")
    print("Therefore Cov[B(d1), B(d2)] = Var[B(d)] when d1 + d2 = 0 mod p.")
    print()
    print("SUM CONSTRAINT: sum_{d=1}^{p-1} B(d) = s(s-1) = const (s = (p-1)/2).")
    print("Therefore Var(sum) = 0, which forces:")
    print("  (p-1) Var(B) + (p-1) [Var(B) + (p-3) c] = 0")
    print("  => c = -2 Var(B) / (p-3)")
    print("where c = Cov[B(d1),B(d2)] for d1+d2 != 0 mod p.")
    print()

    # Verify closed forms
    print("CLOSED FORMS:")
    print("  Var[B(d)]     = (p-3)(p+1) / (16(p-2))")
    print("  Cov (non-comp)= -(p+1) / (8(p-2))        [d1+d2 != 0 mod p]")
    print("  Cov (comp)    = Var[B(d)]                  [d1+d2 = 0 mod p]")
    print("  Correlation   = -2/(p-3)                   [non-complementary]")
    print()

    print("Verification:")
    all_ok = True
    for p in [7, 11, 19, 23, 31, 43, 47, 59, 67, 83, 199, 997]:
        var_B = exact_Var_B(p)
        var_formula = Fraction((p-3)*(p+1), 16*(p-2))
        cov_formula = Fraction(-(p+1), 8*(p-2))

        var_ok = (var_B == var_formula)
        if not var_ok:
            all_ok = False

        # Check non-complementary
        cov_nc, _, _, _ = exact_cross_covariance_fast(p, 1, 2)
        nc_ok = (cov_nc == cov_formula)
        if not nc_ok:
            all_ok = False

        # Check complementary
        cov_c, _, _, _ = exact_cross_covariance_fast(p, 1, p - 1)
        c_ok = (cov_c == var_B)
        if not c_ok:
            all_ok = False

        status = "OK" if (var_ok and nc_ok and c_ok) else "FAIL"
        print(f"  p={p:5d}: Var={float(var_B):10.6f}={float(var_formula):10.6f}, "
              f"Cov_nc={float(cov_nc):10.6f}={float(cov_formula):10.6f}, "
              f"Cov_c={float(cov_c):10.6f}=Var  [{status}]")

    print(f"\n  ALL CLOSED FORMS {'VERIFIED' if all_ok else 'FAILED'}.")

    # =========================================================
    # Summary
    # =========================================================
    print("\n" + "=" * 90)
    print("SUMMARY")
    print("=" * 90)

    print("""
EXACT CLOSED-FORM RESULTS:

  Var[B(d)] = (p-3)(p+1) / (16(p-2))  =  p/16 + O(1)

  Cov[B(d1), B(d2)] for d1+d2 != 0 mod p:
    = -(p+1) / (8(p-2))  =  -1/8 + O(1/p)

  Cov[B(d1), B(d2)] for d1+d2 = 0 mod p:
    = Var[B(d)]  (since B(d) = B(p-d) identically)

  Correlation (non-complementary):
    rho = Cov/Var = -2/(p-3)  ->  0 as p -> infinity

PROOF DERIVATION:
  1. B(d) = B(p-d) identically (autocorrelation symmetry), so Cov(d, p-d) = Var(d).
  2. sum_{d=1}^{p-1} B(d) = s(s-1) is constant => Var(sum) = 0.
  3. By symmetry, all non-complementary Cov values are equal (call it c).
  4. Expanding Var(sum) = 0 with (p-1)/2 complementary pairs:
       (p-1) Var + sum_{d1!=d2} Cov(d1,d2) = 0
     Each d has 1 complementary partner (Cov = Var) and (p-3) non-complementary (Cov = c):
       (p-1) Var + (p-1)[Var + (p-3)c] = 0
       => c = -2 Var/(p-3) = -(p+1)/(8(p-2))
  5. |Cov|/Var = 2/(p-3) = O(1/p).

IMPLICATIONS FOR THE PROOF:
  - For the multivariate CLT on B = (B(d1),...,B(d_m)):
    The covariance matrix Sigma has diagonal p/16 + O(1) and off-diagonal -1/8 + O(1/p).
    The matrix is (p/16) I - (1/8) J/(p-3) + lower order, where J is all-ones.
  - The eigenvalues are all Theta(p) except one that is 0 (from the sum constraint).
  - For the REDUCED vector (after identifying B(d) = B(p-d) and removing the
    constant-sum constraint), the covariance matrix is well-conditioned with
    eigenvalues Theta(p), supporting the multivariate CLT.
""")

    # Covariance matrix eigenvalue structure
    for p in [23, 47, 83]:
        n = p - 1
        var_f = float(Fraction((p-3)*(p+1), 16*(p-2)))
        cov_nc_f = float(Fraction(-(p+1), 8*(p-2)))
        # Build matrix
        M = np.full((n, n), cov_nc_f)
        np.fill_diagonal(M, var_f)
        # Complementary pairs: (d, p-d) have Cov = Var
        for i in range(n):
            d = i + 1
            j = (p - d) - 1  # index of p-d
            if i != j:
                M[i, j] = var_f
        eigvals = sorted(np.linalg.eigvalsh(M))
        print(f"  p={p}: eigenvalues of {n}x{n} Cov matrix:")
        print(f"    min={eigvals[0]:.4f}, 2nd={eigvals[1]:.4f}, "
              f"max={eigvals[-1]:.4f}, #near-zero={sum(1 for e in eigvals if abs(e) < 1e-6)}")


if __name__ == '__main__':
    main()
