#!/usr/bin/env python3
"""QR Perturbation Approach for R(B_{n-1}, B_n) = 4n-1

For primes p = 2n-1 ≡ 3 mod 4, test D12 = (QR ∪ {0}) \ {r} for r ∈ QR.
This introduces controlled variation in B(d) = C(D12, D12; d).

Key formula (for generic d ≠ 0, r, p-r):
  B_new(d) = (p-3)/4 - (chi(r+d) + chi(r-d)) / 2

where chi is the Legendre symbol. B_new(d) ∈ {(p-7)/4, (p-3)/4, (p+1)/4}.
Edge cases at d=r, d=p-r have correction -1/2 relative to naive formula.
"""

import json
import numpy as np
from itertools import combinations
from collections import Counter
import time
import os


def is_prime(n):
    if n < 2:
        return False
    for i in range(2, int(n**0.5) + 1):
        if n % i == 0:
            return False
    return True


def quadratic_residues(p):
    """Compute QR(p) = {x^2 mod p : x in 1..p-1}, sorted."""
    return sorted(set((x * x) % p for x in range(1, p)))


def legendre_symbol(a, p):
    """Compute (a/p) via Euler criterion."""
    a = a % p
    if a == 0:
        return 0
    val = pow(a, (p - 1) // 2, p)
    return 1 if val == 1 else -1


def autocorrelation(S, m):
    """Compute C(S, S; d) for all d in Z_m using FFT."""
    indicator = np.zeros(m)
    for x in S:
        indicator[x % m] = 1
    F = np.fft.fft(indicator)
    power = np.abs(F) ** 2
    result = np.real(np.fft.ifft(power))
    return np.rint(result).astype(int)


def autocorrelation_direct(S, m):
    """Direct computation of C(S, S; d) for verification."""
    S_set = frozenset(x % m for x in S)
    result = [0] * m
    for d in range(m):
        result[d] = sum(1 for x in S_set if (x + d) % m in S_set)
    return result


# ==================== Step 1: Formula Verification ====================

def verify_perturbation_formula():
    """Verify B_new(d) formula against direct computation for all primes."""
    print("=" * 80)
    print("STEP 1: Verify perturbation formula")
    print("=" * 80)
    print()
    print("Formula: B_new(d) = (p-3)/4 - (chi(r+d) + chi(r-d))/2")
    print("Valid for generic d. Edge cases (d=r, d=p-r) have correction -1/2.")
    print()

    primes = [7, 11, 19, 23, 31, 43, 47, 59]
    all_verified = True

    for p in primes:
        qr = quadratic_residues(p)
        qr_set = set(qr)
        base = (p - 3) // 4

        generic_ok = True
        edge_ok = True

        for r in qr:
            D12_new = sorted((qr_set | {0}) - {r})
            B_direct = autocorrelation_direct(D12_new, p)

            for d in range(1, p):
                chi_plus = legendre_symbol(r + d, p)
                chi_minus = legendre_symbol((r - d) % p, p)
                formula_val = base - (chi_plus + chi_minus) / 2
                is_edge = (d == r or d == (p - r) % p)

                if is_edge:
                    if abs(B_direct[d] - (formula_val - 0.5)) > 0.01:
                        edge_ok = False
                else:
                    if abs(B_direct[d] - formula_val) > 0.01:
                        generic_ok = False

        # Show B_new distribution for first r
        r0 = qr[0]
        D12_test = sorted((qr_set | {0}) - {r0})
        B_test = autocorrelation_direct(D12_test, p)
        dist = Counter(B_test[d] for d in range(1, p))
        dist_str = ", ".join(f"{v}:{c}" for v, c in sorted(dist.items()))

        status_g = "OK" if generic_ok else "FAIL"
        status_e = "OK" if edge_ok else "FAIL"
        print(f"  p={p:2d}: generic={status_g}, edge={status_e}, "
              f"base={base}, dist(r={r0}): {{{dist_str}}}")

        if not generic_ok:
            all_verified = False

    print()
    if all_verified:
        print("  >>> Formula VERIFIED for all generic positions.")
    else:
        print("  >>> Formula FAILED!")
    return all_verified


# ==================== Step 2: Sweep r for known D11 ====================

def sweep_all_r(D11, p, n, A_pre=None):
    """Try all r in QR(p) with D12 = (QR ∪ {0}) \\ {r}.

    Returns (successes, total, details).
    details: list of (r, valid, max_excess).
    """
    qr = quadratic_residues(p)
    qr_set = set(qr)
    D11_set = set(D11)
    k = len(D11_set)

    A = A_pre if A_pre is not None else autocorrelation(list(D11_set), p)

    d11_thresh = n - 2
    d22_thresh = 2 * k - n + 1

    successes = 0
    details = []

    for r in qr:
        D12_new = sorted((qr_set | {0}) - {r})
        B = autocorrelation(D12_new, p)

        max_excess = -999
        for d in range(1, p):
            ab = int(A[d] + B[d])
            thresh = d11_thresh if d in D11_set else d22_thresh
            excess = ab - thresh
            if excess > max_excess:
                max_excess = excess

        valid = (max_excess <= 0)
        if valid:
            successes += 1
        details.append((r, valid, max_excess))

    return successes, len(qr), details


def load_registry_solutions():
    """Load solutions for p ≡ 3 mod 4 primes from registry."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(script_dir, 'solutions_registry.json')
    try:
        with open(path, 'r') as f:
            data = json.load(f)
    except FileNotFoundError:
        print("  WARNING: solutions_registry.json not found")
        return {}

    solutions = {}
    for sol in data['solutions']:
        m = sol['m']
        if is_prime(m) and m % 4 == 3:
            solutions[m] = {
                'n': sol['n'],
                'D11': sol['D11'],
                'D12': sol['D12'],
                'solver': sol['solver'],
            }
    return solutions


def find_all_symmetric_d11(p, k):
    """Enumerate all symmetric D11 of size k in {1,...,p-1}."""
    pairs = [(x, p - x) for x in range(1, (p + 1) // 2)]
    num_pairs = k // 2
    result = []
    for combo in combinations(range(len(pairs)), num_pairs):
        D11 = []
        for i in combo:
            D11.append(pairs[i][0])
            D11.append(pairs[i][1])
        result.append(sorted(D11))
    return result


def run_step2():
    """Sweep all r for each prime p ≡ 3 mod 4."""
    print("\n" + "=" * 80)
    print("STEP 2: Sweep all r ∈ QR for each prime")
    print("=" * 80)

    primes = [7, 11, 19, 23, 31, 43, 47, 59]
    registry = load_registry_solutions()
    results = {}

    for p in primes:
        n = (p + 1) // 2
        t0 = time.time()

        if p in registry:
            D11 = registry[p]['D11']
            source = f"registry ({registry[p]['solver']})"
            successes, total, details = sweep_all_r(D11, p, n)
        else:
            # Exhaustive search over all symmetric D11 of standard size k=n
            all_d11 = find_all_symmetric_d11(p, n)
            best_d11 = None
            best_success = -1
            best_details = None

            for candidate in all_d11:
                s, t, det = sweep_all_r(candidate, p, n)
                if s > best_success:
                    best_success = s
                    best_d11 = candidate
                    best_details = det
                    if s == t:
                        break  # All r work, no need to continue

            D11 = best_d11
            successes = best_success
            total = (p - 1) // 2
            details = best_details
            source = f"best of {len(all_d11)} D11s"

        k = len(D11)
        elapsed = time.time() - t0
        pct = successes / total * 100 if total > 0 else 0

        working_r = [r for r, valid, _ in details if valid]
        failing = [(r, excess) for r, valid, excess in details if not valid]

        print(f"\n  p={p:3d}: n={n:2d}, k={k:2d}, |QR|={total:2d}, "
              f"success: {successes}/{total} ({pct:.1f}%)  "
              f"[{source}] ({elapsed:.1f}s)")
        print(f"    D11 = {D11}")
        if working_r:
            print(f"    Working r: {working_r}")
        if failing:
            excesses = [e for _, e in failing]
            print(f"    Failing: excess in [{min(excesses)}, {max(excesses)}]")

        # For one working pair, show A(d)+B(d) profile
        if working_r:
            r0 = working_r[0]
            qr_set = set(quadratic_residues(p))
            D12_new = sorted((qr_set | {0}) - {r0})
            A = autocorrelation(D11, p)
            B = autocorrelation(D12_new, p)
            D11_set = set(D11)
            d11_thresh = n - 2
            d22_thresh = 2 * k - n + 1
            d11_abs = [int(A[d] + B[d]) for d in range(1, p) if d in D11_set]
            d22_abs = [int(A[d] + B[d]) for d in range(1, p) if d not in D11_set]
            print(f"    r={r0}: A+B at D11 in [{min(d11_abs)},{max(d11_abs)}] "
                  f"(thresh {d11_thresh}), "
                  f"D22 in [{min(d22_abs)},{max(d22_abs)}] "
                  f"(thresh {d22_thresh})")

        results[p] = {
            'n': n, 'k': k, 'successes': successes, 'total': total,
            'rate': pct, 'working_r': working_r, 'D11': D11, 'source': source,
        }

    return results


# ==================== Step 3: Exhaustive D11 Search ====================

def run_step3():
    """Exhaustive D11 search for p=23 and p=31."""
    print("\n" + "=" * 80)
    print("STEP 3: Exhaustive D11 search (p=23, p=31)")
    print("=" * 80)

    results = {}

    for p in [23, 31]:
        n = (p + 1) // 2
        k = n
        all_d11 = find_all_symmetric_d11(p, k)
        num_d11 = len(all_d11)
        qr = quadratic_residues(p)
        num_qr = len(qr)

        print(f"\n  p={p}: n={n}, |D11 candidates|={num_d11}, |QR|={num_qr}")
        print(f"  Total (D11, r) pairs: {num_d11 * num_qr}")

        t0 = time.time()
        d11_with_working = 0
        total_working = 0
        total_pairs = 0
        best_d11 = None
        best_count = 0
        r_success_counts = Counter()

        for idx, D11 in enumerate(all_d11):
            A = autocorrelation(D11, p)
            s, t, details = sweep_all_r(D11, p, n, A_pre=A)
            total_pairs += t
            total_working += s

            if s > 0:
                d11_with_working += 1
            if s > best_count:
                best_count = s
                best_d11 = D11

            for r, valid, _ in details:
                if valid:
                    r_success_counts[r] += 1

            if (idx + 1) % 1000 == 0:
                elapsed = time.time() - t0
                print(f"    [{elapsed:.1f}s] {idx+1}/{num_d11}: "
                      f"{d11_with_working} D11s with working r, "
                      f"{total_working}/{total_pairs} pairs")

        elapsed = time.time() - t0
        d11_rate = d11_with_working / num_d11 * 100
        pair_rate = total_working / total_pairs * 100 if total_pairs > 0 else 0

        print(f"\n  Results for p={p} ({elapsed:.1f}s):")
        print(f"    D11s with >=1 working r: "
              f"{d11_with_working}/{num_d11} ({d11_rate:.1f}%)")
        print(f"    Total working pairs: "
              f"{total_working}/{total_pairs} ({pair_rate:.1f}%)")
        print(f"    Best D11: {best_d11} ({best_count}/{num_qr} r values)")

        print(f"    Per-r success rates:")
        for r in qr:
            cnt = r_success_counts.get(r, 0)
            print(f"      r={r:2d}: {cnt}/{num_d11} D11s ({cnt/num_d11*100:.1f}%)")

        results[p] = {
            'd11_with_working': d11_with_working,
            'total_d11': num_d11,
            'total_working': total_working,
            'total_pairs': total_pairs,
            'best_d11': best_d11,
            'best_count': best_count,
        }

    return results


# ==================== Cross-check ====================

def cross_check_small(D11, D12, m):
    """Full vertex-by-vertex verification for small m (≤ 15)."""
    if m > 15:
        return None
    n = (m + 1) // 2
    D11_set = set(D11)
    D12_set = set(D12)
    D22_set = set(range(1, m)) - D11_set

    def is_red(ub, ui, vb, vi):
        if ub == vb:
            d = (ui - vi) % m
            return d in (D11_set if ub == 0 else D22_set)
        # Cross-block: convention is d = (V1_index - V2_index)
        if ub == 0:
            d = (ui - vi) % m
        else:
            d = (vi - ui) % m
        return d in D12_set

    max_red_cn = 0
    max_blue_cn = 0
    verts = [(b, i) for b in range(2) for i in range(m)]

    for a in range(len(verts)):
        for b in range(a + 1, len(verts)):
            u, v = verts[a], verts[b]
            red_cn = blue_cn = 0
            for c in range(len(verts)):
                w = verts[c]
                if w == u or w == v:
                    continue
                r1 = is_red(u[0], u[1], w[0], w[1])
                r2 = is_red(v[0], v[1], w[0], w[1])
                if r1 and r2:
                    red_cn += 1
                elif not r1 and not r2:
                    blue_cn += 1
            if is_red(u[0], u[1], v[0], v[1]):
                max_red_cn = max(max_red_cn, red_cn)
            else:
                max_blue_cn = max(max_blue_cn, blue_cn)

    valid = (max_red_cn <= n - 2) and (max_blue_cn <= n - 1)
    return valid, max_red_cn, max_blue_cn


# ==================== Main ====================

def main():
    print("=" * 80)
    print("QR PERTURBATION APPROACH: D12 = (QR ∪ {0}) \\ {r}")
    print("for R(B_{n-1}, B_n) = 4n-1, primes p = 2n-1 ≡ 3 mod 4")
    print("=" * 80)

    # Step 1
    verify_perturbation_formula()

    # Step 2
    step2 = run_step2()

    # Cross-check one working pair for small p
    print("\n  Cross-check (full vertex verification for p <= 15):")
    for p in [7, 11]:
        if p in step2 and step2[p]['working_r']:
            D11 = step2[p]['D11']
            r0 = step2[p]['working_r'][0]
            qr_set = set(quadratic_residues(p))
            D12_new = sorted((qr_set | {0}) - {r0})
            result = cross_check_small(D11, D12_new, p)
            if result is not None:
                valid, max_r, max_b = result
                n = (p + 1) // 2
                print(f"    p={p}, r={r0}: valid={valid}, "
                      f"max_red_cn={max_r} (≤{n-2}), "
                      f"max_blue_cn={max_b} (≤{n-1})")

    # Step 3
    step3 = run_step3()

    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY TABLE")
    print("=" * 80)
    print(f"\n{'p':>4s} {'n':>4s} {'k':>4s} {'|QR|':>5s} "
          f"{'success':>10s} {'rate':>8s}")
    print("-" * 45)
    for p in sorted(step2):
        r = step2[p]
        print(f"{p:4d} {r['n']:4d} {r['k']:4d} {r['total']:5d} "
              f"{r['successes']:3d}/{r['total']:<3d}    {r['rate']:5.1f}%")

    # Step 4 analysis hints
    print("\n" + "=" * 80)
    print("STEP 4: Analysis")
    print("=" * 80)

    rates = [step2[p]['rate'] for p in sorted(step2)]
    primes_list = sorted(step2)

    if all(r > 0 for r in rates):
        print("\n  STRONG SIGNAL: All primes have success rate > 0%.")
        if min(rates) >= 10:
            print("  Success rate >= 10% everywhere — averaging argument viable.")
        else:
            print(f"  Minimum rate = {min(rates):.1f}% at p={primes_list[rates.index(min(rates))]}")
            print("  Rate may be decreasing — check trend.")
    else:
        zero_primes = [p for p, r in zip(primes_list, rates) if r == 0]
        print(f"\n  WEAK SIGNAL: Zero success at p = {zero_primes}")
        print("  Single-swap perturbation insufficient for these primes.")
        print("  Consider multi-swap: D12 = (QR ∪ {0} ∪ QNR_add) \\ QR_remove")

    # Rate trend
    print("\n  Rate trend:")
    for p in primes_list:
        r = step2[p]
        bar = "#" * int(r['rate'] / 2)
        print(f"    p={p:3d}: {r['rate']:5.1f}% {bar}")


if __name__ == '__main__':
    main()
