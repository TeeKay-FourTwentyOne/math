#!/usr/bin/env python3
"""
Verify candidate separating invariants across ALL orbits at p=31, 43.
Not just a sample — compute for every orbit.

Key invariants to test:
1. DFT flatness (max|hat(1_{D11})(k)|^2 / mean)
2. (a_d11_var, a_d22_var)
3. (a_d22_max, a_d11_var)

The fundamental question: does any algebraic invariant separate working from non-working
orbits at ALL tested primes?
"""

import json
import cmath
import math
from collections import Counter


def primitive_root(p):
    phi = p - 1
    factors = set()
    n = phi
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors.add(d)
            n //= d
        d += 1
    if n > 1:
        factors.add(n)
    for g in range(2, p):
        if all(pow(g, phi // f, p) != 1 for f in factors):
            return g


def compute_A_values(d11_set, p):
    """Compute A(d) for all d in Z_p*."""
    a_vals = {}
    for d in range(1, p):
        count = sum(1 for x in d11_set if (x + d) % p in d11_set)
        a_vals[d] = count
    return a_vals


def compute_dft_flatness(d11_set, p):
    """Compute max|hat(1_{D11})(k)|^2 / mean over k=1..p-1."""
    mags = []
    for k in range(1, p):
        omega = cmath.exp(2j * cmath.pi * k / p)
        s = sum(omega ** d for d in d11_set)
        mags.append(abs(s) ** 2)
    mean_mag = sum(mags) / len(mags)
    max_mag = max(mags)
    return max_mag / mean_mag if mean_mag > 0 else 0, max_mag, mean_mag


def compute_orbit_invariants(d11_list, p):
    """Compute all candidate invariants for a D11."""
    d11_set = set(d11_list)
    k = len(d11_set)

    # A-values
    a_vals = compute_A_values(d11_set, p)

    # Split into D11 and D22 positions
    d22_set = set(range(1, p)) - d11_set
    a_d11 = [a_vals[d] for d in sorted(d11_set)]
    a_d22 = [a_vals[d] for d in sorted(d22_set)]

    # Variances
    mean_a = sum(a_d11) / len(a_d11)
    a_d11_var = sum((x - mean_a)**2 for x in a_d11) / len(a_d11)

    mean_a22 = sum(a_d22) / len(a_d22)
    a_d22_var = sum((x - mean_a22)**2 for x in a_d22) / len(a_d22)

    # Max values
    max_A = max(a_vals.values())
    a_d22_max = max(a_d22) if a_d22 else 0
    a_d11_max = max(a_d11) if a_d11 else 0
    a_d11_min = min(a_d11) if a_d11 else 0

    # DFT flatness
    dft_flat, dft_max, dft_mean = compute_dft_flatness(d11_set, p)

    # Sorted A-profiles (as tuples for hashing)
    a_d11_sorted = tuple(sorted(a_d11))
    a_d22_sorted = tuple(sorted(a_d22))

    # A-value range on D11
    a_d11_range = max(a_d11) - min(a_d11)

    # Median A on D11
    sa = sorted(a_d11)
    a_d11_median = (sa[len(sa)//2] + sa[(len(sa)-1)//2]) / 2

    return {
        'a_d11_var': round(a_d11_var, 4),
        'a_d22_var': round(a_d22_var, 4),
        'a_d22_max': a_d22_max,
        'a_d11_max': a_d11_max,
        'a_d11_min': a_d11_min,
        'a_d11_range': a_d11_range,
        'a_d11_median': a_d11_median,
        'max_A': max_A,
        'dft_flatness': round(dft_flat, 4),
        'dft_max': round(dft_max, 2),
        'a_d11_sorted': a_d11_sorted,
        'a_d22_sorted': a_d22_sorted,
    }


def get_all_orbits(p):
    """Generate ALL symmetric D11 orbits for prime p ≡ 3 mod 4."""
    from itertools import combinations

    n_pairs = (p - 1) // 2
    pairs = [(d, p - d) for d in range(1, (p + 1) // 2)]
    k_pairs = (p + 1) // 4  # number of pairs in D11

    g = primitive_root(p)

    # Generate all k_pairs-subsets of pairs
    seen_orbits = set()
    orbits = []

    for combo in combinations(range(n_pairs), k_pairs):
        # Convert to D11
        d11 = set()
        for idx in combo:
            d11.add(pairs[idx][0])
            d11.add(pairs[idx][1])
        d11_frozen = frozenset(d11)

        # Check if this orbit was already seen
        if d11_frozen in seen_orbits:
            continue

        # Generate full orbit under multiplication by g
        current = d11_frozen
        orbit_members = set()
        for _ in range(p - 1):
            orbit_members.add(current)
            current = frozenset((x * g) % p for x in current)

        seen_orbits.update(orbit_members)
        orbits.append(sorted(d11))

    return orbits


def main():
    # Load working orbit data
    with open('c0_estimation_results.json') as f:
        c0_data = json.load(f)

    with open('near_flat_p43_results.json') as f:
        p43_data = json.load(f)

    # ========= p = 31 =========
    print("=" * 70)
    print("p = 31: Checking ALL orbits")
    print("=" * 70)

    p = 31
    # Get working D11 sets
    working_d11_sets = set()
    for orbit in c0_data['p31']['orbit_results']:
        if orbit['num_valid'] > 0:
            working_d11_sets.add(frozenset(orbit['canonical']))

    print(f"Working orbits: {len(working_d11_sets)}")

    # Compute invariants for ALL orbits
    all_orbits = get_all_orbits(p)
    print(f"Total orbits: {len(all_orbits)}")

    working_invs = []
    nonworking_invs = []

    for i, d11 in enumerate(all_orbits):
        if (i + 1) % 50 == 0:
            print(f"  Processing orbit {i+1}/{len(all_orbits)}...")
        inv = compute_orbit_invariants(d11, p)
        is_working = frozenset(d11) in working_d11_sets
        inv['working'] = is_working
        inv['D11'] = d11

        if is_working:
            working_invs.append(inv)
        else:
            nonworking_invs.append(inv)

    print(f"\nWorking: {len(working_invs)}, Non-working: {len(nonworking_invs)}")

    # Test separators
    test_separators(working_invs, nonworking_invs, p)

    # ========= p = 43 =========
    print("\n" + "=" * 70)
    print("p = 43: Checking all orbits in dataset")
    print("=" * 70)

    p = 43
    # At p=43, full enumeration is expensive. Use the orbits from near_flat_p43_results.json
    # which covers all A-flat orbits.
    working_invs43 = []
    nonworking_invs43 = []

    total_orbits = len(p43_data['orbit_results'])
    print(f"Total orbits in dataset: {total_orbits}")

    for i, orbit in enumerate(p43_data['orbit_results']):
        if (i + 1) % 20 == 0:
            print(f"  Processing orbit {i+1}/{total_orbits}...")
        d11 = orbit['canonical']
        inv = compute_orbit_invariants(d11, p)
        is_working = orbit['working']
        inv['working'] = is_working
        inv['D11'] = d11

        if is_working:
            working_invs43.append(inv)
        else:
            nonworking_invs43.append(inv)

    print(f"\nWorking: {len(working_invs43)}, Non-working: {len(nonworking_invs43)}")

    test_separators(working_invs43, nonworking_invs43, p)

    # ========= Cross-prime summary =========
    print("\n" + "=" * 70)
    print("CROSS-PRIME SUMMARY")
    print("=" * 70)

    # Load p=23 exact data for comparison
    with open('schur_convexity_exact.json') as f:
        exact = json.load(f)

    for pkey in ['23', '19', '11']:
        p = int(pkey)
        working_invs_p = []
        nonworking_invs_p = []
        for orbit in exact[pkey]['orbit_summaries']:
            d11 = orbit['orbit_rep']
            inv = compute_orbit_invariants(d11, p)
            is_working = orbit['N'] > 0
            inv['working'] = is_working
            if is_working:
                working_invs_p.append(inv)
            else:
                nonworking_invs_p.append(inv)

        print(f"\np={p}: {len(working_invs_p)} working, {len(nonworking_invs_p)} non-working")
        test_separators(working_invs_p, nonworking_invs_p, p)


def test_separators(working, nonworking, p):
    """Test all candidate single and pair separators."""

    # Single invariants
    single_keys = ['dft_flatness', 'a_d11_var', 'a_d22_var', 'a_d22_max', 'max_A',
                   'a_d11_range', 'a_d11_min', 'a_d11_max', 'a_d11_median']

    print(f"\n--- Single invariant separation at p={p} ---")
    for key in single_keys:
        w_vals = [inv[key] for inv in working]
        nw_vals = [inv[key] for inv in nonworking]
        w_min, w_max = min(w_vals), max(w_vals)
        nw_min, nw_max = min(nw_vals), max(nw_vals)

        # Check if ranges are disjoint
        if w_max < nw_min or nw_max < w_min:
            print(f"  ** SEPARATES: {key}")
            print(f"     W=[{w_min}, {w_max}], NW=[{nw_min}, {nw_max}]")
        else:
            # Count misclassified
            # Try threshold: W <= T < NW or W >= T > NW
            # Best threshold minimizing errors
            all_vals = [(v, True) for v in w_vals] + [(v, False) for v in nw_vals]
            all_vals.sort(key=lambda x: x[0])

            # Count non-working in working range
            nw_in_w = sum(1 for v in nw_vals if w_min <= v <= w_max)
            w_in_nw = sum(1 for v in w_vals if nw_min <= v <= nw_max)

            overlap_frac_nw = nw_in_w / len(nw_vals) * 100
            print(f"  {key}: W=[{w_min}, {w_max}], NW=[{nw_min}, {nw_max}] -- "
                  f"{nw_in_w}/{len(nw_vals)} NW in W-range ({overlap_frac_nw:.1f}%)")

    # Pair separators
    print(f"\n--- Pair invariant separation at p={p} ---")
    pair_keys = [
        ('a_d11_var', 'a_d22_var'),
        ('a_d22_max', 'a_d11_var'),
        ('dft_flatness', 'a_d11_var'),
        ('dft_flatness', 'a_d22_var'),
        ('dft_flatness', 'a_d22_max'),
        ('a_d11_var', 'a_d22_max'),
        ('a_d11_range', 'a_d22_var'),
        ('a_d11_range', 'dft_flatness'),
        ('a_d11_min', 'a_d22_var'),
    ]

    for k1, k2 in pair_keys:
        w_pairs = set((inv[k1], inv[k2]) for inv in working)
        nw_pairs = set((inv[k1], inv[k2]) for inv in nonworking)
        overlap = w_pairs & nw_pairs
        if not overlap:
            # Check if a simple linear separator exists
            # (working on one side, non-working on other)
            print(f"  ** SEPARATES: ({k1}, {k2})")
            print(f"     W pairs: {sorted(w_pairs)}")
        else:
            print(f"  ({k1}, {k2}): {len(overlap)} overlapping pairs")
            if len(overlap) <= 3:
                print(f"     Overlap: {sorted(overlap)}")

    # Check sorted A-profile (fingerprint)
    w_profiles = set(inv['a_d11_sorted'] + inv['a_d22_sorted'] for inv in working)
    nw_profiles = set(inv['a_d11_sorted'] + inv['a_d22_sorted'] for inv in nonworking)
    profile_overlap = w_profiles & nw_profiles
    print(f"\n  Full A-profile: {len(profile_overlap)} overlapping profiles "
          f"(should be 0 since profiles determine orbit)")


if __name__ == '__main__':
    main()
