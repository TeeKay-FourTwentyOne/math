#!/usr/bin/env python3
"""
Exact coverage analysis: which n values have R(B_{n-1}, B_n) = 4n-1 proven?

Sources:
1. n ≤ 21: verified (Rousseau-Sheehan + computational verification)
2. q ≡ 1 mod 4 prime power: Paley construction (q = 2n-1)
3. p ≡ 3 mod 4 prime, p ≤ 227: hyperplane conditioning proof
"""

from math import gcd


def is_prime(n):
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0 or n % 3 == 0: return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0: return False
        i += 6
    return True


def factorize(n):
    """Return list of (prime, power) pairs."""
    factors = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            power = 0
            while n % d == 0:
                n //= d
                power += 1
            factors.append((d, power))
        d += 1
    if n > 1:
        factors.append((n, 1))
    return factors


def is_prime_power(n):
    """Check if n is a prime power (p^k for prime p, k ≥ 1)."""
    factors = factorize(n)
    return len(factors) == 1


def coverage_reason(n, m):
    """Determine why n is proven (or not)."""
    if n <= 21:
        return "verified (n ≤ 21)"

    # Check if m is a prime power ≡ 1 mod 4
    if is_prime_power(m) and m % 4 == 1:
        factors = factorize(m)
        p, k = factors[0]
        if k == 1:
            return f"Paley (p={m} ≡ 1 mod 4)"
        else:
            return f"Paley (q={p}^{k}={m} ≡ 1 mod 4)"

    # Check if m is a prime ≡ 3 mod 4 ≤ 227
    if is_prime(m) and m % 4 == 3 and m <= 227:
        return f"HC proof (p={m} ≡ 3 mod 4)"

    return None  # NOT proven


def main():
    max_n = 250
    proven = {}
    gaps = []

    for n in range(2, max_n + 1):
        m = 2 * n - 1
        reason = coverage_reason(n, m)
        if reason:
            proven[n] = reason
        else:
            factors = factorize(m)
            gaps.append((n, m, factors))

    total = max_n - 1  # n=2 to max_n
    print(f"COVERAGE: n = 2 to {max_n}")
    print(f"  Proven: {len(proven)}/{total} ({100*len(proven)/total:.1f}%)")
    print(f"  Gaps: {len(gaps)}/{total} ({100*len(gaps)/total:.1f}%)")
    print()

    # Breakdown by reason
    reason_counts = {}
    for n, reason in proven.items():
        key = reason.split("(")[0].strip()
        reason_counts[key] = reason_counts.get(key, 0) + 1

    print("BREAKDOWN:")
    for reason, count in sorted(reason_counts.items(), key=lambda x: -x[1]):
        print(f"  {reason}: {count}")
    print()

    # List gaps
    print(f"GAPS (first 50 of {len(gaps)}):")
    for n, m, factors in gaps[:50]:
        fact_str = " × ".join(f"{p}^{k}" if k > 1 else str(p) for p, k in factors)
        m_mod4 = m % 4
        # Check what blocks this m
        if is_prime(m):
            if m % 4 == 3:
                block = f"prime ≡ 3 mod 4, p={m} > 227"
            else:
                block = "BUG"  # should be covered by Paley
        elif is_prime_power(m):
            if m % 4 == 3:
                block = f"prime power ≡ 3 mod 4"
            else:
                block = "BUG"  # should be covered
        else:
            block = "composite (not prime power)"

        print(f"  n={n:3d}, m={m:3d} = {fact_str:20s} [{block}]")

    # Categorize gaps
    gap_types = {
        'composite_not_pp': [],
        'pp_3mod4': [],
        'prime_gt227': [],
    }
    for n, m, factors in gaps:
        if is_prime(m) and m % 4 == 3:
            gap_types['prime_gt227'].append((n, m))
        elif is_prime_power(m) and m % 4 == 3:
            gap_types['pp_3mod4'].append((n, m))
        else:
            gap_types['composite_not_pp'].append((n, m))

    print(f"\nGAP CATEGORIES:")
    print(f"  Composite (not prime power): {len(gap_types['composite_not_pp'])}")
    print(f"  Prime power ≡ 3 mod 4: {len(gap_types['pp_3mod4'])}")
    print(f"  Prime ≡ 3 mod 4, > 227: {len(gap_types['prime_gt227'])}")

    if gap_types['pp_3mod4']:
        print(f"\n  Prime power ≡ 3 mod 4 gaps:")
        for n, m in gap_types['pp_3mod4'][:10]:
            factors = factorize(m)
            fact_str = " × ".join(f"{p}^{k}" if k > 1 else str(p) for p, k in factors)
            print(f"    n={n}, m={m} = {fact_str}")

    if gap_types['prime_gt227']:
        print(f"\n  Prime ≡ 3 mod 4, > 227 gaps:")
        for n, m in gap_types['prime_gt227'][:5]:
            print(f"    n={n}, m={m}")

    # Density analysis
    print(f"\nDENSITY BY RANGE:")
    for lo in range(2, max_n, 50):
        hi = min(lo + 49, max_n)
        total_range = hi - lo + 1
        proven_range = sum(1 for n in range(lo, hi + 1) if n in proven)
        print(f"  n={lo:3d}-{hi:3d}: {proven_range:2d}/{total_range} "
              f"({100*proven_range/total_range:.0f}%)")


if __name__ == '__main__':
    main()
