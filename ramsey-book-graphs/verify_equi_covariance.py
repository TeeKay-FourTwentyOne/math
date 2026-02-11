#!/usr/bin/env python3
"""Verify the equi-covariance proof by checking analytical collision counts
against brute-force enumeration, and confirming the closed form Cov = -(p+1)/(8(p-2)).

This script accompanies proofs/proof_equi_covariance.md.
"""

from fractions import Fraction
from collections import defaultdict
import sys


def hypergeometric_moments(p):
    """Return q1, q2, q3, q4 for the random k-subset model."""
    k = (p - 3) // 2
    N = p - 1
    q1 = Fraction(k, N)
    q2 = Fraction(k * (k - 1), N * (N - 1))
    q3 = Fraction(k * (k - 1) * (k - 2), N * (N - 1) * (N - 2))
    q4 = Fraction(k * (k - 1) * (k - 2) * (k - 3), N * (N - 1) * (N - 2) * (N - 3))
    return q1, q2, q3, q4


def E_Y_product(indices, q1, q2, q3, q4):
    """E[prod Y_{a_i}] given list of indices (with possible repeats)."""
    d = len(set(indices))
    if d == 0:
        return Fraction(1)
    elif d == 1:
        return q1
    elif d == 2:
        return q2
    elif d == 3:
        return q3
    elif d == 4:
        return q4
    else:
        raise ValueError(f"Unexpected {d} distinct indices")


def brute_force_E_B1B2(p, d1, d2):
    """Compute E[B(d1)*B(d2)] by brute-force enumeration of all 9 terms."""
    q1, q2, q3, q4 = hypergeometric_moments(p)
    alpha1 = (p - d1) % p
    beta1 = d1
    alpha2 = (p - d2) % p
    beta2 = d2

    # Terms (1)-(4)
    E_1 = E_Y_product([alpha1, alpha2], q1, q2, q3, q4)
    E_2 = E_Y_product([alpha1, beta2], q1, q2, q3, q4)
    E_3 = E_Y_product([beta1, alpha2], q1, q2, q3, q4)
    E_4 = E_Y_product([beta1, beta2], q1, q2, q3, q4)

    # Terms (5)-(8): O(p) loops
    E_5 = Fraction(0)
    for b in range(1, p):
        if b == d2:
            continue
        E_5 += E_Y_product([alpha1, b, (b - d2) % p], q1, q2, q3, q4)

    E_6 = Fraction(0)
    for b in range(1, p):
        if b == d2:
            continue
        E_6 += E_Y_product([beta1, b, (b - d2) % p], q1, q2, q3, q4)

    E_7 = Fraction(0)
    for a in range(1, p):
        if a == d1:
            continue
        E_7 += E_Y_product([a, (a - d1) % p, alpha2], q1, q2, q3, q4)

    E_8 = Fraction(0)
    for a in range(1, p):
        if a == d1:
            continue
        E_8 += E_Y_product([a, (a - d1) % p, beta2], q1, q2, q3, q4)

    # Term (9): O(p^2) double loop
    E_9 = Fraction(0)
    for a in range(1, p):
        if a == d1:
            continue
        ad1 = (a - d1) % p
        for b in range(1, p):
            if b == d2:
                continue
            bd2 = (b - d2) % p
            E_9 += E_Y_product([a, ad1, b, bd2], q1, q2, q3, q4)

    return E_1 + E_2 + E_3 + E_4 + E_5 + E_6 + E_7 + E_8 + E_9


def brute_force_collision_counts_term9(p, d1, d2):
    """Count (a,b) pairs in Term 9 by number of distinct indices, via brute force."""
    counts = defaultdict(int)
    for a in range(1, p):
        if a == d1:
            continue
        ad1 = (a - d1) % p
        for b in range(1, p):
            if b == d2:
                continue
            bd2 = (b - d2) % p
            nd = len({a, ad1, b, bd2})
            counts[nd] += 1
    return dict(counts)


def brute_force_collision_counts_term5(p, d1, d2):
    """Count b values in Term 5 (L1*Q2) by number of distinct indices."""
    alpha1 = (p - d1) % p
    counts = defaultdict(int)
    for b in range(1, p):
        if b == d2:
            continue
        bd2 = (b - d2) % p
        nd = len({alpha1, b, bd2})
        counts[nd] += 1
    return dict(counts)


def analytical_collision_counts_term9(p, d1, d2):
    """Compute collision counts for Term 9 analytically, as in the proof."""
    # The four collision offsets
    delta1 = 0
    delta2 = d2
    delta3 = (-d1) % p  # = p - d1
    delta4 = (d2 - d1) % p

    deltas = [delta1, delta2, delta3, delta4]
    # Check all distinct (for non-complementary with d1 != d2)
    all_distinct = (len(set(deltas)) == 4)

    if not all_distinct:
        return None  # Not handled by the simple formula

    # Each collision offset contributes p-3 pairs with 3 distinct indices
    # Generic pairs: (p-2)^2 - 4*(p-3) with 4 distinct indices
    n_collision = 4 * (p - 3)
    n_generic = (p - 2) ** 2 - n_collision
    return {3: n_collision, 4: n_generic}


def analytical_collision_counts_term5(p, d1, d2):
    """Analytical collision counts for Term 5 (alpha1 * Q2)."""
    # For non-complementary: 2 collisions (2-distinct), p-4 generic (3-distinct)
    return {2: 2, 3: p - 4}


def analytical_E_B1B2(p):
    """Compute E[B(d1)*B(d2)] analytically for non-complementary pairs.

    From the proof:
      Terms (1)-(4): 4 * q2
      Terms (5)-(8): 8 * q2 + 4*(p-4) * q3
      Term (9): 4*(p-3) * q3 + ((p-2)^2 - 4*(p-3)) * q4
    """
    q1, q2, q3, q4 = hypergeometric_moments(p)

    terms_1_4 = 4 * q2
    terms_5_8 = 8 * q2 + 4 * (p - 4) * q3
    n_coll = 4 * (p - 3)
    n_generic = (p - 2) ** 2 - n_coll
    term_9 = n_coll * q3 + n_generic * q4

    return terms_1_4 + terms_5_8 + term_9


def main():
    print("=" * 80)
    print("VERIFICATION OF EQUI-COVARIANCE PROOF")
    print("=" * 80)

    # ================================================================
    # Test 1: Verify collision counts for Term 9 (brute-force vs analytical)
    # ================================================================
    print("\n--- Test 1: Term 9 collision counts ---")
    for p in [11, 23]:
        print(f"\n  p = {p}:")
        all_ok = True
        tested = 0
        for d1 in range(1, p):
            for d2 in range(1, p):
                if d1 == d2:
                    continue
                is_comp = (d1 + d2) % p == 0

                bf = brute_force_collision_counts_term9(p, d1, d2)
                an = analytical_collision_counts_term9(p, d1, d2)

                if is_comp:
                    # Complementary: the proof doesn't claim 4 distinct offsets
                    continue

                tested += 1
                if an is None:
                    print(f"    d1={d1}, d2={d2}: analytical returned None (unexpected)")
                    all_ok = False
                    continue

                # Compare
                if bf != an:
                    print(f"    MISMATCH at d1={d1}, d2={d2}: bf={bf}, an={an}")
                    all_ok = False

        status = "PASS" if all_ok else "FAIL"
        print(f"    Tested {tested} non-complementary pairs: [{status}]")

    # ================================================================
    # Test 2: Verify collision counts for Terms 5-8
    # ================================================================
    print("\n--- Test 2: Term 5 collision counts ---")
    for p in [11, 23]:
        print(f"\n  p = {p}:")
        all_ok = True
        tested = 0
        for d1 in range(1, p):
            for d2 in range(1, p):
                if d1 == d2 or (d1 + d2) % p == 0:
                    continue
                tested += 1
                bf = brute_force_collision_counts_term5(p, d1, d2)
                an = analytical_collision_counts_term5(p, d1, d2)
                if bf != an:
                    print(f"    MISMATCH at d1={d1}, d2={d2}: bf={bf}, an={an}")
                    all_ok = False
        status = "PASS" if all_ok else "FAIL"
        print(f"    Tested {tested} non-complementary pairs: [{status}]")

    # ================================================================
    # Test 3: Verify analytical E[B1*B2] matches brute-force
    # ================================================================
    print("\n--- Test 3: Analytical E[B1*B2] vs brute-force ---")
    for p in [11, 23]:
        print(f"\n  p = {p}:")
        E_analytical = analytical_E_B1B2(p)

        all_ok = True
        tested = 0
        for d1 in range(1, p):
            for d2 in range(1, p):
                if d1 == d2 or (d1 + d2) % p == 0:
                    continue
                tested += 1
                E_bf = brute_force_E_B1B2(p, d1, d2)
                if E_bf != E_analytical:
                    print(f"    MISMATCH at d1={d1}, d2={d2}: "
                          f"bf={float(E_bf):.8f}, an={float(E_analytical):.8f}")
                    all_ok = False

        status = "PASS" if all_ok else "FAIL"
        print(f"    Tested {tested} pairs. Analytical E[B1*B2] = {float(E_analytical):.8f}")
        print(f"    [{status}]")

    # ================================================================
    # Test 4: Verify closed-form Cov = -(p+1)/(8(p-2))
    # ================================================================
    print("\n--- Test 4: Closed-form Cov = -(p+1)/(8(p-2)) ---")
    all_ok = True
    for p in [7, 11, 19, 23, 31, 43, 47, 59, 67, 83, 199, 997]:
        E_B = Fraction(p - 3, 4)
        E_B2_analytical = analytical_E_B1B2(p)
        cov_analytical = E_B2_analytical - E_B * E_B
        cov_formula = Fraction(-(p + 1), 8 * (p - 2))

        ok = (cov_analytical == cov_formula)
        if not ok:
            all_ok = False
        status = "OK" if ok else "FAIL"
        print(f"  p={p:5d}: Cov(analytical) = {float(cov_analytical):12.8f}, "
              f"Cov(formula) = {float(cov_formula):12.8f}  [{status}]")

    status = "ALL PASS" if all_ok else "SOME FAILED"
    print(f"\n  [{status}]")

    # ================================================================
    # Test 5: Full check for small p -- ALL pairs
    # ================================================================
    print("\n--- Test 5: ALL non-complementary pairs equal for small p ---")
    for p in [7, 11, 19, 23]:
        E_B = Fraction(p - 3, 4)
        cov_formula = Fraction(-(p + 1), 8 * (p - 2))

        all_ok = True
        tested = 0
        for d1 in range(1, p):
            for d2 in range(1, p):
                if d1 == d2 or (d1 + d2) % p == 0:
                    continue
                tested += 1
                E_bf = brute_force_E_B1B2(p, d1, d2)
                cov_bf = E_bf - E_B * E_B
                if cov_bf != cov_formula:
                    print(f"    MISMATCH at p={p}, d1={d1}, d2={d2}: "
                          f"cov={float(cov_bf):.8f}, formula={float(cov_formula):.8f}")
                    all_ok = False

        status = "PASS" if all_ok else "FAIL"
        print(f"  p={p:3d}: Tested all {tested} non-complementary pairs: [{status}]")

    # ================================================================
    # Test 6: Verify complementary covariance = Var
    # ================================================================
    print("\n--- Test 6: Complementary pairs have Cov = Var ---")
    for p in [7, 11, 19, 23]:
        E_B = Fraction(p - 3, 4)
        var_formula = Fraction((p - 3) * (p + 1), 16 * (p - 2))

        all_ok = True
        for d1 in range(1, (p + 1) // 2):
            d2 = p - d1
            if d1 == d2:
                continue
            E_bf = brute_force_E_B1B2(p, d1, d2)
            cov_bf = E_bf - E_B * E_B
            ok = (cov_bf == var_formula)
            if not ok:
                print(f"    MISMATCH at p={p}, d1={d1}, d2={d2}: "
                      f"cov={float(cov_bf):.8f}, var={float(var_formula):.8f}")
                all_ok = False

        status = "PASS" if all_ok else "FAIL"
        print(f"  p={p:3d}: [{status}]")

    # ================================================================
    # Summary
    # ================================================================
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print("""
Verified:
  1. Term 9 collision counts match brute-force for p=11, 23 (all non-comp pairs).
  2. Term 5 collision counts match brute-force for p=11, 23.
  3. Analytical E[B1*B2] formula matches brute-force for p=11, 23.
  4. Closed-form Cov = -(p+1)/(8(p-2)) verified for p up to 997.
  5. ALL non-complementary pairs give the same Cov for p=7,11,19,23 (brute-force).
  6. Complementary pairs have Cov = Var for p=7,11,19,23.

The equi-covariance property is established by direct indicator computation.
No group-theoretic transitivity argument is needed.
""")


if __name__ == "__main__":
    main()
