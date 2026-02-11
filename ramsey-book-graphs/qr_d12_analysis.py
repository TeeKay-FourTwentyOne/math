"""
Analysis: Using D12 = QR(p) for primes p ≡ 3 mod 4.

Key hypothesis: For p ≡ 3 mod 4, setting D12 = QR(p) gives
B(d) = C(QR, QR; d) = (p-3)/4 for ALL nonzero d.

This would reduce the full constraint system to:
  A(d) <= (p-3)/4       for d in D11   (binding)
  A(d) <= (p+9)/4       for d in D22   (loose, +3 margin over avg)

where A(d) = C(D11, D11; d) with average (p+1)/4.

The binding constraint is: A(d) <= average - 1 for d in D11.

This is a pure autocorrelation problem on a symmetric subset of Z_p.
"""

import numpy as np
from collections import defaultdict


def quadratic_residues(p):
    """Compute QR(p) = {x^2 mod p : x in 1..p-1}"""
    return sorted(set((x * x) % p for x in range(1, p)))


def autocorrelation(S, m):
    """Compute C(S, S; d) for all d in Z_m."""
    indicator = np.zeros(m, dtype=int)
    for x in S:
        indicator[x % m] = 1
    result = np.zeros(m, dtype=int)
    for d in range(m):
        count = 0
        for x in range(m):
            if indicator[x] and indicator[(x + d) % m]:
                count += 1
        result[d] = count
    return result


def verify_qr_constant_autocorrelation():
    """Verify B(d) = (p-3)/4 for all nonzero d when D12 = QR(p), p ≡ 3 mod 4."""
    primes_3mod4 = [7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83]

    print("=" * 70)
    print("VERIFICATION: QR autocorrelation is constant for p ≡ 3 mod 4")
    print("=" * 70)

    all_constant = True
    for p in primes_3mod4:
        qr = quadratic_residues(p)
        B = autocorrelation(qr, p)

        expected = (p - 3) // 4
        nonzero_values = [B[d] for d in range(1, p)]
        is_constant = all(v == expected for v in nonzero_values)
        min_val = min(nonzero_values)
        max_val = max(nonzero_values)

        status = "CONSTANT ✓" if is_constant else f"NOT CONSTANT (range {min_val}-{max_val})"
        print(f"  p={p:3d}: |QR|={len(qr):2d}, B(d)={(p-3)//4} expected, "
              f"actual range [{min_val}, {max_val}] — {status}")

        if not is_constant:
            all_constant = False

    print()
    if all_constant:
        print(">>> CONFIRMED: B(d) = (p-3)/4 for ALL nonzero d, ALL primes p ≡ 3 mod 4 tested.")
    else:
        print(">>> FAILED: Some primes have non-constant QR autocorrelation!")

    return all_constant


def analyze_constraint_reduction():
    """
    With D12 = QR(p), the constraint system becomes:
    A(d) + (p-3)/4 <= n-2     for d in D11   =>  A(d) <= (p-3)/4
    A(d) + (p-3)/4 <= 2k-n+1  for d in D22   =>  A(d) <= 2k-n+1-(p-3)/4

    With k = |D11| = (p+1)/2 (standard size):
    D11 threshold: A(d) <= (p-3)/4
    D22 threshold: A(d) <= 2(p+1)/2 - (p+1)/2 + 1 - (p-3)/4 = (p+9)/4

    Average A(d) for d != 0: (p+1)/4
    D11 threshold = avg - 1
    D22 threshold = avg + 2
    """
    primes = [7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83]

    print("=" * 70)
    print("CONSTRAINT ANALYSIS with D12 = QR(p)")
    print("=" * 70)
    print(f"{'p':>4s} {'n':>4s} {'|D11|':>5s} {'avg A':>7s} {'D11 thr':>8s} "
          f"{'D22 thr':>8s} {'margin1':>8s} {'margin2':>8s} {'total slack':>12s}")
    print("-" * 70)

    for p in primes:
        n = (p + 1) // 2
        k = (p + 1) // 2  # |D11| = n (standard)

        avg_A = (p + 1) / 4
        d11_thresh = (p - 3) / 4  # = avg - 1
        d22_thresh = (p + 9) / 4  # = avg + 2

        margin1 = avg_A - d11_thresh  # how far below avg D11 needs
        margin2 = d22_thresh - avg_A  # how far above avg D22 allows

        # Total slack
        num_d11 = k  # (p+1)/2 positions
        num_d22 = p - 1 - k  # (p-3)/2 positions
        total_capacity = num_d11 * d11_thresh + num_d22 * d22_thresh
        total_needed = k * (k - 1) + (n - 1) * (n - 2)  # Σ A(d) = |D11|²-|D11|
        # Actually total A is |D11|(|D11|-1), B is constant
        total_A = k * (k - 1)
        total_slack = total_capacity - total_A

        print(f"{p:4d} {n:4d} {k:5d} {avg_A:7.2f} {d11_thresh:8.2f} "
              f"{d22_thresh:8.2f} {margin1:8.2f} {margin2:8.2f} {total_slack:12.1f}")


def test_d12_qr_on_known_solutions():
    """
    For known solutions, check: if we REPLACE D12 with QR(p), does the
    existing D11 still work? (I.e., is A(d) <= (p-3)/4 for all d in D11?)
    """
    # Known solutions from registry
    solutions = {
        22: {
            'D11': [1,2,5,10,11,13,16,17,18,19,20,23,24,25,26,27,30,32,33,38,41,42],
            'D12': [0,2,5,6,8,11,15,16,20,24,25,27,28,31,32,34,35,36,37,39,41],
        },
        24: {
            'D11': [2,6,8,10,11,13,14,16,19,20,21,26,27,28,31,33,34,36,37,39,41,45],
            'D12': [0,1,8,9,10,12,19,23,24,25,26,28,29,32,33,35,38,39,40,42,44,45,46],
        },
        30: {
            'D11': [1,2,3,5,6,8,9,10,11,17,19,25,27,29,30,32,34,40,42,48,49,50,51,53,54,56,57,58],
            'D12': [0,1,4,6,7,9,10,11,13,18,20,21,23,24,25,26,32,33,36,38,39,42,43,45,49,50,51,54,55],
        },
    }

    print("\n" + "=" * 70)
    print("TEST: Can known D11 work with D12 = QR(p)?")
    print("=" * 70)

    for n, sol in sorted(solutions.items()):
        p = 2 * n - 1
        D11 = set(sol['D11'])
        D11_actual_D12 = sol['D12']
        qr = set(quadratic_residues(p))

        k = len(D11)

        # Compute A(d) with actual D11
        A = autocorrelation(D11, p)

        # Compute B(d) with actual D12 and with QR
        B_actual = autocorrelation(D11_actual_D12, p)
        B_qr = autocorrelation(qr, p)

        d11_threshold = n - 2  # for A(d) + B(d)
        d22_threshold = 2 * k - n + 1

        # Check with actual D12
        violations_actual = 0
        for d in range(1, p):
            ab = A[d] + B_actual[d]
            if d in D11 and ab > d11_threshold:
                violations_actual += 1
            elif d not in D11 and d != 0 and ab > d22_threshold:
                violations_actual += 1

        # Check with QR D12
        violations_qr = 0
        violated_d11 = []
        violated_d22 = []
        for d in range(1, p):
            ab = A[d] + B_qr[d]
            if d in D11 and ab > d11_threshold:
                violations_qr += 1
                violated_d11.append((d, A[d], B_qr[d], ab))
            elif d not in D11 and d != 0 and ab > d22_threshold:
                violations_qr += 1
                violated_d22.append((d, A[d], B_qr[d], ab))

        print(f"\nn={n}, p={p}, |D11|={k}")
        print(f"  Actual D12: {violations_actual} violations (should be 0)")
        print(f"  D12 = QR:   {violations_qr} violations")

        if violations_qr > 0:
            A_at_d11 = [A[d] for d in range(1, p) if d in D11]
            qr_b = (p - 3) // 4
            a_threshold = d11_threshold - qr_b
            print(f"  A threshold for D11: A(d) <= {a_threshold} = {d11_threshold} - {qr_b}")
            print(f"  A(d) range at D11 positions: [{min(A_at_d11)}, {max(A_at_d11)}]")
            print(f"  A(d) mean at D11 positions: {np.mean(A_at_d11):.2f}")
            if violated_d11:
                print(f"  D11 violations (d, A, B, A+B):")
                for v in violated_d11[:5]:
                    print(f"    d={v[0]}: A={v[1]}, B={v[2]}, A+B={v[3]} > {d11_threshold}")
            if violated_d22:
                print(f"  D22 violations: {len(violated_d22)}")
        else:
            print(f"  >>> D12=QR WORKS with this D11!")

        # Also report: is the actual D12 close to QR?
        overlap = len(set(D11_actual_D12) & qr)
        print(f"  |actual_D12 ∩ QR| = {overlap} / {len(D11_actual_D12)}")


def search_d11_for_qr_d12():
    """
    For small primes p ≡ 3 mod 4, search for D11 (symmetric) such that
    A(d) <= (p-3)/4 for all d in D11, with |D11| = (p+1)/2.
    """
    print("\n" + "=" * 70)
    print("SEARCH: Find D11 with A(d) ≤ (p-3)/4 at D11 positions")
    print("=" * 70)

    import random

    for p in [7, 11, 19, 23, 31, 43]:
        n = (p + 1) // 2
        k = (p + 1) // 2  # target |D11|
        num_pairs_needed = k // 2

        # Build symmetric pairs
        pairs = []
        for x in range(1, (p + 1) // 2):
            pairs.append((x, p - x))

        total_pairs = len(pairs)
        a_threshold = (p - 3) // 4

        print(f"\np={p}, n={n}, |D11|={k}, pairs={total_pairs}, need={num_pairs_needed}")
        print(f"  A(d) threshold = {a_threshold}, average A = {(p+1)/4:.2f}")

        # Try random symmetric D11
        successes = 0
        trials = 10000 if p <= 23 else 1000

        for trial in range(trials):
            # Randomly select num_pairs_needed pairs
            selected = random.sample(range(total_pairs), num_pairs_needed)
            D11 = set()
            for i in selected:
                D11.add(pairs[i][0])
                D11.add(pairs[i][1])

            A = autocorrelation(D11, p)

            # Check: A(d) <= a_threshold for all d in D11
            ok = True
            for d in D11:
                if A[d] > a_threshold:
                    ok = False
                    break

            if ok:
                successes += 1

        rate = successes / trials * 100
        print(f"  Random search: {successes}/{trials} valid D11 found ({rate:.1f}%)")

        if successes > 0:
            # Find one and display
            for trial in range(100000):
                selected = random.sample(range(total_pairs), num_pairs_needed)
                D11 = set()
                for i in selected:
                    D11.add(pairs[i][0])
                    D11.add(pairs[i][1])

                A = autocorrelation(D11, p)
                ok = all(A[d] <= a_threshold for d in D11)

                if ok:
                    max_a_d11 = max(A[d] for d in D11)
                    max_a_d22 = max(A[d] for d in range(1, p) if d not in D11)
                    print(f"  Example D11: {sorted(D11)}")
                    print(f"  max A at D11: {max_a_d11} (threshold {a_threshold})")
                    print(f"  max A at D22: {max_a_d22} (threshold {(p+9)//4})")
                    break


def compare_actual_vs_qr_b():
    """Compare B(d) from actual solutions vs constant (p-3)/4 from QR."""
    solutions = {
        22: [0,2,5,6,8,11,15,16,20,24,25,27,28,31,32,34,35,36,37,39,41],
        24: [0,1,8,9,10,12,19,23,24,25,26,28,29,32,33,35,38,39,40,42,44,45,46],
        30: [0,1,4,6,7,9,10,11,13,18,20,21,23,24,25,26,32,33,36,38,39,42,43,45,49,50,51,54,55],
    }

    print("\n" + "=" * 70)
    print("COMPARISON: B(d) from actual D12 vs QR")
    print("=" * 70)

    for n, D12 in sorted(solutions.items()):
        p = 2 * n - 1
        qr = quadratic_residues(p)

        B_actual = autocorrelation(D12, p)
        B_qr = autocorrelation(qr, p)

        actual_vals = [B_actual[d] for d in range(1, p)]
        qr_vals = [B_qr[d] for d in range(1, p)]

        print(f"\nn={n}, p={p}:")
        print(f"  Actual B(d): range [{min(actual_vals)}, {max(actual_vals)}], "
              f"mean={np.mean(actual_vals):.2f}, std={np.std(actual_vals):.2f}")
        print(f"  QR B(d):     range [{min(qr_vals)}, {max(qr_vals)}], "
              f"mean={np.mean(qr_vals):.2f}, std={np.std(qr_vals):.2f}")
        print(f"  QR constant = {(p-3)//4}")

        # The actual D12 has varying B(d) - this variation is what allows
        # the SA to find solutions. But constant B means we need D11
        # to do ALL the work.


if __name__ == '__main__':
    # Step 1: Verify QR autocorrelation is constant for p ≡ 3 mod 4
    verify_qr_constant_autocorrelation()

    # Step 2: Constraint analysis
    analyze_constraint_reduction()

    # Step 3: Compare actual B(d) vs QR B(d)
    compare_actual_vs_qr_b()

    # Step 4: Test known D11 with QR D12
    test_d12_qr_on_known_solutions()

    # Step 5: Search for valid D11 with QR D12
    search_d11_for_qr_d12()
