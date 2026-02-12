"""
Algebraic invariant analysis for D11 orbits.
Compute necklace representations, character sums, gap structures, and
find what separates working (N>0) from non-working (N=0) orbits.
"""

import json
import math
from collections import defaultdict
from itertools import groupby


def primitive_root(p):
    """Find smallest primitive root mod p."""
    if p == 2:
        return 1
    phi = p - 1
    # Factor phi
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
        ok = True
        for f in factors:
            if pow(g, phi // f, p) == 1:
                ok = False
                break
        if ok:
            return g
    return None


def make_symmetric(S, p):
    """Ensure S is symmetric mod p."""
    result = set(S)
    for x in S:
        result.add((-x) % p)
    return result


def pairs_from_d11(d11, p):
    """Convert D11 to set of pair indices.
    Pairs are {d, p-d} for d=1,...,(p-1)/2.
    Return frozenset of pair indices (1-indexed)."""
    pair_set = set()
    for d in d11:
        if d == 0:
            continue
        idx = min(d, p - d)
        pair_set.add(idx)
    return frozenset(pair_set)


def pair_permutation(g, p):
    """Compute how multiplication by g permutes the pairs {d, p-d}.
    Pairs are indexed by min(d, p-d) for d in 1...(p-1)/2.
    Returns dict: pair_idx -> pair_idx."""
    half = (p - 1) // 2
    perm = {}
    for d in range(1, half + 1):
        gd = (g * d) % p
        perm[d] = min(gd, p - gd)
    return perm


def necklace_from_pairs(pair_set, perm, half):
    """Compute necklace string from pair set and permutation.
    Follow the cycle of the permutation starting from pair 1,
    record 1 if pair is in set, 0 otherwise."""
    # Trace the cycle starting from pair index 1
    cycle = []
    current = 1
    for _ in range(half):
        cycle.append(current)
        current = perm[current]

    # Build binary string
    bits = ''.join('1' if c in pair_set else '0' for c in cycle)
    return bits, cycle


def canonical_necklace(bits):
    """Return lexicographically smallest rotation of binary string."""
    n = len(bits)
    doubled = bits + bits
    best = bits
    for i in range(1, n):
        rot = doubled[i:i+n]
        if rot < best:
            best = rot
    return best


def necklace_gaps_and_runs(bits):
    """Compute gap sizes (runs of 0s) and run sizes (runs of 1s), both sorted descending."""
    groups = [(k, len(list(g))) for k, g in groupby(bits)]

    # Handle wrap-around for circular string
    if len(groups) > 1 and groups[0][0] == groups[-1][0]:
        merged_len = groups[0][1] + groups[-1][1]
        merged_key = groups[0][0]
        groups = [(merged_key, merged_len)] + groups[1:-1]

    gaps = sorted([l for k, l in groups if k == '0'], reverse=True)
    runs = sorted([l for k, l in groups if k == '1'], reverse=True)
    return gaps, runs


def character_sum(d11, p, order):
    """Compute character sum S_chi = sum_{d in D11} chi(d) where chi has given order.
    chi(d) = exp(2*pi*i * ind_g(d) * (p-1)/(order * (p-1))) = exp(2*pi*i * ind_g(d) / order)
    Actually: chi of order k means chi = omega^((p-1)/k * .) where omega = exp(2*pi*i/(p-1))

    More precisely: let g be prim root, ind_g(d) is discrete log.
    Character of order k: chi_k(d) = exp(2*pi*i * ind_g(d) / k)

    But there are multiple characters of each order. We'll compute all of them.
    For each divisor k of p-1, characters of order k are chi^j where gcd(j, k) = 1 (up to
    grouping by order). Actually let's compute |sum chi^j(d)|^2 for all j.
    """
    pass  # Will implement differently below


def discrete_log_table(g, p):
    """Compute discrete log table: d -> ind_g(d) for d in 1..p-1."""
    table = {}
    val = 1
    for i in range(p - 1):
        table[val] = i
        val = (val * g) % p
    return table


def all_character_sums(d11, p, g):
    """Compute S_j = sum_{d in D11} omega^(j*ind_g(d)) for all j=0,...,p-2.
    where omega = exp(2*pi*i/(p-1)).

    Returns list of |S_j|^2 for j=0,...,p-2.
    Also returns the actual complex values.
    """
    dlog = discrete_log_table(g, p)
    n = p - 1  # = ord(g)
    omega = complex(math.cos(2 * math.pi / n), math.sin(2 * math.pi / n))

    sums = []
    abs_sq = []
    for j in range(n):
        s = complex(0, 0)
        for d in d11:
            if d == 0:
                continue
            s += omega ** (j * dlog[d])
        sums.append(s)
        abs_sq.append(s.real ** 2 + s.imag ** 2)

    return sums, abs_sq


def coset_intersections(d11, p, g, index):
    """For subgroup H of index `index` in Z_p*, compute |D11 ∩ C_k| for each coset C_k.
    H = {g^(index*t) : t} and cosets C_k = g^k * H for k=0,...,index-1."""
    dlog = discrete_log_table(g, p)
    counts = [0] * index
    for d in d11:
        if d == 0:
            continue
        k = dlog[d] % index
        counts[k] += 1
    return counts


def a_value(d11, d, p):
    """Compute A(d) = |{x in D11 : x-d in D11}| = Delta(D11, D11, d, p)."""
    count = 0
    for x in d11:
        if (x - d) % p in d11:
            count += 1
    return count


def compute_all_invariants(p, orbit_data, g):
    """Compute all algebraic invariants for all orbits at prime p."""
    half = (p - 1) // 2
    perm = pair_permutation(g, p)

    results = []
    for orb in orbit_data:
        d11_list = orb['orbit_rep']
        d11_set = make_symmetric(set(d11_list), p)
        d11_set.discard(0)
        N = orb['N']
        c0 = orb['c0']
        is_a_flat = orb['is_a_flat']
        max_A = orb['max_A_D11']

        # Pair representation
        pair_set = pairs_from_d11(d11_set, p)

        # Necklace
        bits, cycle = necklace_from_pairs(pair_set, perm, half)
        canon = canonical_necklace(bits)
        gaps, runs = necklace_gaps_and_runs(canon)

        # Character sums
        char_sums, char_abs_sq = all_character_sums(d11_set, p, g)

        # Group |S_j|^2 by the order of character j
        # Character j has order (p-1)/gcd(j, p-1)
        char_by_order = defaultdict(list)
        for j in range(p - 1):
            order = (p - 1) // math.gcd(j, p - 1) if j > 0 else 1
            char_by_order[order].append(char_abs_sq[j])

        # Coset intersections for various indices
        divisors = sorted(set(d for d in range(1, p) if (p - 1) % d == 0))
        coset_data = {}
        for idx in divisors:
            coset_data[idx] = sorted(coset_intersections(d11_set, p, g, idx))

        # QR/NQR intersection
        qr_count = sum(1 for d in d11_set if d != 0 and pow(d, (p - 1) // 2, p) == 1)
        nqr_count = len(d11_set) - qr_count

        # A-profile on D11 (re-sorted)
        a_prof_d11 = sorted([a_value(d11_set, d, p) for d in d11_set])

        # A-profile on D22 = complement
        d22_set = set(range(1, p)) - d11_set
        a_prof_d22 = sorted([a_value(d11_set, d, p) for d in d22_set])

        # Flatness measures
        a_range = max(a_prof_d11) - min(a_prof_d11) if a_prof_d11 else 0
        a_var = sum((x - sum(a_prof_d11) / len(a_prof_d11)) ** 2 for x in a_prof_d11) / len(a_prof_d11) if a_prof_d11 else 0

        # Power spectrum: |S_j|^2 for each j, grouped
        # For necklace: the DFT of the necklace indicator
        necklace_dft_sq = []
        n_bits = len(canon)
        for k in range(n_bits):
            omega = complex(math.cos(2 * math.pi * k / n_bits), math.sin(2 * math.pi * k / n_bits))
            s = sum(int(canon[i]) * omega ** (-i) for i in range(n_bits))
            necklace_dft_sq.append(round(s.real ** 2 + s.imag ** 2, 6))

        result = {
            'd11': sorted(d11_set),
            'N': N,
            'c0': c0,
            'is_a_flat': is_a_flat,
            'max_A_D11': max_A,
            'pair_set': sorted(pair_set),
            'necklace': canon,
            'necklace_raw': bits,
            'gaps': gaps,
            'runs': runs,
            'max_gap': max(gaps) if gaps else 0,
            'max_run': max(runs) if runs else 0,
            'num_gaps': len(gaps),
            'num_runs': len(runs),
            'qr_count': qr_count,
            'nqr_count': nqr_count,
            'char_abs_sq': [round(x, 6) for x in char_abs_sq],
            'char_by_order': {str(k): [round(x, 6) for x in v] for k, v in sorted(char_by_order.items())},
            'coset_data': {str(k): v for k, v in sorted(coset_data.items())},
            'a_profile_d11': a_prof_d11,
            'a_profile_d22': a_prof_d22,
            'a_range': a_range,
            'a_var': round(a_var, 6),
            'necklace_dft_sq': necklace_dft_sq,
        }
        results.append(result)

    return results


def analyze_separator(results, invariant_key, extractor=None):
    """Check if an invariant separates working from non-working orbits.
    Returns (separates: bool, working_vals, nonworking_vals)."""
    if extractor is None:
        extractor = lambda r: r[invariant_key]

    working = []
    nonworking = []
    for r in results:
        val = extractor(r)
        if r['N'] > 0:
            working.append(val)
        else:
            nonworking.append(val)

    # Check separation: is there a threshold that perfectly separates?
    # For numeric values
    try:
        w_set = set(working) if not isinstance(working[0], list) else set(tuple(v) for v in working)
        n_set = set(nonworking) if not isinstance(nonworking[0], list) else set(tuple(v) for v in nonworking)
        separates = len(w_set & n_set) == 0
    except (TypeError, IndexError):
        separates = False

    return separates, working, nonworking


def print_analysis(results, p):
    """Print comprehensive analysis."""
    working = [r for r in results if r['N'] > 0]
    nonworking = [r for r in results if r['N'] == 0]

    print(f"\n{'='*80}")
    print(f"ANALYSIS FOR p = {p}")
    print(f"{'='*80}")
    print(f"Total orbits: {len(results)}")
    print(f"Working orbits (N > 0): {len(working)}")
    print(f"Non-working orbits (N = 0): {len(nonworking)}")

    # ---- Necklace patterns ----
    print(f"\n--- NECKLACE PATTERNS ---")
    print(f"{'Necklace':<20s} {'N':>6s} {'c0':>8s} {'Flat':>5s} {'MaxA':>5s} {'Gaps':>15s} {'Runs':>15s} {'MaxGap':>6s} {'MaxRun':>6s}")
    for r in sorted(results, key=lambda x: (-x['N'], x['necklace'])):
        print(f"{r['necklace']:<20s} {r['N']:>6d} {r['c0']:>8.3f} {str(r['is_a_flat']):>5s} {r['max_A_D11']:>5d} "
              f"{str(r['gaps']):>15s} {str(r['runs']):>15s} {r['max_gap']:>6d} {r['max_run']:>6d}")

    # ---- Simple scalar invariants ----
    print(f"\n--- SCALAR INVARIANT SEPARATION TESTS ---")

    tests = [
        ('max_A_D11', 'Max A on D11'),
        ('is_a_flat', 'A-flat?'),
        ('max_gap', 'Max gap in necklace'),
        ('max_run', 'Max run in necklace'),
        ('num_gaps', 'Number of gaps'),
        ('num_runs', 'Number of runs'),
        ('a_range', 'A-profile range'),
        ('a_var', 'A-profile variance'),
        ('qr_count', 'QR count'),
        ('nqr_count', 'NQR count'),
    ]

    for key, label in tests:
        sep, wv, nv = analyze_separator(results, key)
        w_vals = sorted(set(wv))
        n_vals = sorted(set(nv))
        status = "SEPARATES" if sep else "FAILS"
        print(f"  {label:<30s} {status:<12s}  Working: {w_vals}  Non-working: {n_vals}")

    # ---- Gap/Run profile ----
    print(f"\n--- GAP PROFILE SEPARATION ---")
    sep, wv, nv = analyze_separator(results, 'gaps', lambda r: tuple(r['gaps']))
    print(f"  Gap profile: {'SEPARATES' if sep else 'FAILS'}")
    if sep or len(set(tuple(x) for x in wv)) <= 10:
        print(f"    Working gap profiles: {sorted(set(tuple(x) for x in wv))}")
    else:
        print(f"    Working gap profiles: {len(set(tuple(x) for x in wv))} distinct")
    if not sep:
        overlap = set(tuple(x) for x in wv) & set(tuple(x) for x in nv)
        print(f"    Overlapping profiles: {sorted(overlap)}")

    sep, wv, nv = analyze_separator(results, 'runs', lambda r: tuple(r['runs']))
    print(f"  Run profile: {'SEPARATES' if sep else 'FAILS'}")
    if sep or len(set(tuple(x) for x in wv)) <= 10:
        print(f"    Working run profiles: {sorted(set(tuple(x) for x in wv))}")

    # ---- Character sums by order ----
    print(f"\n--- CHARACTER SUM ANALYSIS ---")
    # Get all orders
    all_orders = set()
    for r in results:
        all_orders.update(int(k) for k in r['char_by_order'].keys())

    for order in sorted(all_orders):
        if order == 1:
            continue
        # For each order, compute max |S_chi|^2 and min |S_chi|^2
        key = str(order)
        def max_char(r, k=key):
            vals = r['char_by_order'].get(k, [0])
            return round(max(vals), 2)
        def min_char(r, k=key):
            vals = r['char_by_order'].get(k, [0])
            return round(min(vals), 2)
        def mean_char(r, k=key):
            vals = r['char_by_order'].get(k, [0])
            return round(sum(vals) / len(vals), 2)

        sep_max, wv_max, nv_max = analyze_separator(results, None, max_char)
        sep_min, wv_min, nv_min = analyze_separator(results, None, min_char)
        sep_mean, wv_mean, nv_mean = analyze_separator(results, None, mean_char)

        w_max_range = (min(wv_max), max(wv_max)) if wv_max else (0, 0)
        n_max_range = (min(nv_max), max(nv_max)) if nv_max else (0, 0)
        w_mean_range = (min(wv_mean), max(wv_mean)) if wv_mean else (0, 0)
        n_mean_range = (min(nv_mean), max(nv_mean)) if nv_mean else (0, 0)

        print(f"  Order {order}:")
        print(f"    Max |S|^2: Working [{w_max_range[0]:.1f}, {w_max_range[1]:.1f}]  "
              f"Non-working [{n_max_range[0]:.1f}, {n_max_range[1]:.1f}]  "
              f"{'SEPARATES' if sep_max else 'FAILS'}")
        print(f"    Mean |S|^2: Working [{w_mean_range[0]:.1f}, {w_mean_range[1]:.1f}]  "
              f"Non-working [{n_mean_range[0]:.1f}, {n_mean_range[1]:.1f}]  "
              f"{'SEPARATES' if sep_mean else 'FAILS'}")

    # ---- Coset intersections ----
    print(f"\n--- COSET INTERSECTION ANALYSIS ---")
    for r in results[:1]:
        divisors = sorted(int(k) for k in r['coset_data'].keys())

    for idx in divisors:
        if idx == 1 or idx == p - 1:
            continue
        key = str(idx)
        def coset_profile(r, k=key):
            return tuple(r['coset_data'].get(k, []))
        sep, wv, nv = analyze_separator(results, None, coset_profile)
        w_profiles = sorted(set(wv))
        n_profiles = sorted(set(nv))
        if sep:
            print(f"  Index {idx}: SEPARATES")
            print(f"    Working: {w_profiles}")
        else:
            overlap = set(wv) & set(nv)
            print(f"  Index {idx}: FAILS (overlap: {len(overlap)} profiles)")

    # ---- A-profile separation ----
    print(f"\n--- A-PROFILE SEPARATION ---")
    sep_d11, wv, nv = analyze_separator(results, 'a_profile_d11', lambda r: tuple(r['a_profile_d11']))
    sep_d22, wv2, nv2 = analyze_separator(results, 'a_profile_d22', lambda r: tuple(r['a_profile_d22']))
    print(f"  D11 A-profile: {'SEPARATES' if sep_d11 else 'FAILS'}")
    print(f"  D22 A-profile: {'SEPARATES' if sep_d22 else 'FAILS'}")

    # Combined D11+D22
    def combined_a(r):
        return (tuple(r['a_profile_d11']), tuple(r['a_profile_d22']))
    sep_both, wv3, nv3 = analyze_separator(results, None, combined_a)
    print(f"  Combined (D11, D22) A-profile: {'SEPARATES' if sep_both else 'FAILS'}")

    # ---- Necklace DFT ----
    print(f"\n--- NECKLACE DFT ANALYSIS ---")
    def necklace_max_dft(r):
        # Max non-DC component
        return round(max(r['necklace_dft_sq'][1:]), 2)
    def necklace_sum_dft(r):
        return round(sum(r['necklace_dft_sq'][1:]), 2)
    sep_dft_max, wv, nv = analyze_separator(results, None, necklace_max_dft)
    sep_dft_sum, wv2, nv2 = analyze_separator(results, None, necklace_sum_dft)
    print(f"  Max non-DC DFT: {'SEPARATES' if sep_dft_max else 'FAILS'}")
    print(f"    Working: {sorted(set(wv))}")
    print(f"    Non-working range: [{min(nv):.2f}, {max(nv):.2f}]")
    print(f"  Sum non-DC DFT: {'SEPARATES' if sep_dft_sum else 'FAILS'}")

    # ---- Combination tests ----
    print(f"\n--- COMBINATION TESTS ---")

    # is_a_flat AND max_gap <= threshold
    for max_gap_thresh in range(1, (p-1)//2 + 1):
        def combo1(r, t=max_gap_thresh):
            return r['is_a_flat'] and r['max_gap'] <= t
        w_pass = sum(1 for r in working if combo1(r))
        n_pass = sum(1 for r in nonworking if combo1(r))
        if w_pass == len(working) and n_pass == 0:
            print(f"  A-flat AND max_gap <= {max_gap_thresh}: SEPARATES (all {len(working)} working pass, 0 non-working pass)")
            break
        elif w_pass == len(working):
            print(f"  A-flat AND max_gap <= {max_gap_thresh}: all working pass but {n_pass} non-working also pass")

    # is_a_flat AND a_var <= threshold
    for thresh_100 in range(0, 500, 10):
        thresh = thresh_100 / 100.0
        def combo2(r, t=thresh):
            return r['is_a_flat'] and r['a_var'] <= t
        w_pass = sum(1 for r in working if combo2(r))
        n_pass = sum(1 for r in nonworking if combo2(r))
        if w_pass == len(working) and n_pass == 0:
            print(f"  A-flat AND a_var <= {thresh:.2f}: SEPARATES")
            break

    # max_A_D11 <= threshold
    for t in range(1, 15):
        w_pass = sum(1 for r in working if r['max_A_D11'] <= t)
        n_pass = sum(1 for r in nonworking if r['max_A_D11'] <= t)
        if w_pass == len(working) and n_pass == 0:
            print(f"  max_A_D11 <= {t}: SEPARATES")
            break

    # Necklace spread: max_gap + max_run <= threshold
    for t in range(2, 20):
        w_pass = sum(1 for r in working if r['max_gap'] + r['max_run'] <= t)
        n_pass = sum(1 for r in nonworking if r['max_gap'] + r['max_run'] <= t)
        if w_pass == len(working) and n_pass == 0:
            print(f"  max_gap + max_run <= {t}: SEPARATES")
            break

    # Character sum threshold tests
    for order in sorted(all_orders):
        if order == 1:
            continue
        key = str(order)
        for thresh_10 in range(0, 200):
            thresh = thresh_10 / 1.0
            def char_test(r, k=key, t=thresh):
                vals = r['char_by_order'].get(k, [0])
                return max(vals) <= t
            w_pass = sum(1 for r in working if char_test(r))
            n_pass = sum(1 for r in nonworking if char_test(r))
            if w_pass == len(working) and n_pass == 0:
                print(f"  max |S_chi|^2 (order {order}) <= {thresh}: SEPARATES")
                break

    # ---- Print detailed working orbit info ----
    print(f"\n--- DETAILED WORKING ORBITS ---")
    for r in sorted(working, key=lambda x: -x['N']):
        print(f"\n  D11: {r['d11']}")
        print(f"  N={r['N']}, c0={r['c0']:.3f}, A-flat={r['is_a_flat']}, maxA={r['max_A_D11']}")
        print(f"  Necklace: {r['necklace']}")
        print(f"  Gaps: {r['gaps']}, Runs: {r['runs']}")
        print(f"  QR: {r['qr_count']}, NQR: {r['nqr_count']}")
        print(f"  A-profile D11: {r['a_profile_d11']}")
        print(f"  A-profile D22: {r['a_profile_d22']}")
        # Character sums - just max for each order
        for order in sorted(int(k) for k in r['char_by_order'].keys()):
            if order == 1:
                continue
            vals = r['char_by_order'][str(order)]
            print(f"  Order {order}: max|S|^2={max(vals):.2f}, mean|S|^2={sum(vals)/len(vals):.2f}")

    return working, nonworking


def main():
    # Load orbit data
    with open('/Users/stephenpadgett/Projects/math/ramsey-book-graphs/schur_convexity_exact.json') as f:
        all_data = json.load(f)

    all_results = {}

    for p_str in ['11', '19', '23']:
        p = int(p_str)
        data = all_data[p_str]
        g = primitive_root(p)
        print(f"\nPrimitive root of {p}: {g}")

        results = compute_all_invariants(p, data['orbit_summaries'], g)
        all_results[p_str] = results

        working, nonworking = print_analysis(results, p)

    # ---- CROSS-VALIDATION ----
    print(f"\n{'='*80}")
    print("CROSS-VALIDATION SUMMARY")
    print(f"{'='*80}")

    # Collect candidate separators found at p=23 and test at p=11, p=19
    for p_str in ['11', '19', '23']:
        results = all_results[p_str]
        p = int(p_str)
        working = [r for r in results if r['N'] > 0]
        nonworking = [r for r in results if r['N'] == 0]

        print(f"\n  p={p}: {len(working)} working, {len(nonworking)} non-working")

        # Test: is_a_flat
        w_flat = sum(1 for r in working if r['is_a_flat'])
        n_flat = sum(1 for r in nonworking if r['is_a_flat'])
        print(f"    A-flat: {w_flat}/{len(working)} working, {n_flat}/{len(nonworking)} non-working")

        # Test: max_A_D11 <= k for various k
        for k in range(4, 10):
            w_pass = sum(1 for r in working if r['max_A_D11'] <= k)
            n_pass = sum(1 for r in nonworking if r['max_A_D11'] <= k)
            if w_pass == len(working):
                print(f"    max_A <= {k}: {w_pass}/{len(working)} working, {n_pass}/{len(nonworking)} non-working")
                break

        # Test: combined A-profile
        sep, _, _ = analyze_separator(results, None, lambda r: (tuple(r['a_profile_d11']), tuple(r['a_profile_d22'])))
        print(f"    Combined A-profile separates: {sep}")

        # Test: necklace pattern
        sep, _, _ = analyze_separator(results, 'necklace')
        print(f"    Necklace separates: {sep}")

        # Test various max_gap thresholds (combined with a_flat)
        for mg in range(1, (p-1)//2 + 1):
            def combo(r, t=mg):
                return r['is_a_flat'] and r['max_gap'] <= t
            w_p = sum(1 for r in working if combo(r))
            n_p = sum(1 for r in nonworking if combo(r))
            if w_p == len(working) and n_p == 0:
                print(f"    A-flat + max_gap <= {mg}: SEPARATES")
                break
            elif w_p == len(working):
                print(f"    A-flat + max_gap <= {mg}: all working, {n_p} non-working also pass")

    # ---- DEEPER ANALYSIS ----
    print(f"\n{'='*80}")
    print("DEEPER ANALYSIS")
    print(f"{'='*80}")

    for p_str in ['23', '19', '11']:
        p = int(p_str)
        results = all_results[p_str]
        working = [r for r in results if r['N'] > 0]
        nonworking = [r for r in results if r['N'] == 0]

        print(f"\n--- p={p} DEEPER ANALYSIS ---")

        # Test: A-profile range (max - min of A on D11) <= threshold
        for t in range(0, 10):
            w_pass = sum(1 for r in working if r['a_range'] <= t)
            n_pass = sum(1 for r in nonworking if r['a_range'] <= t)
            if w_pass > 0:
                print(f"  a_range <= {t}: {w_pass}/{len(working)} working, {n_pass}/{len(nonworking)} non-working pass")
            if w_pass == len(working) and n_pass == 0:
                print(f"    *** SEPARATES ***")
                break

        # Test: A-profile D22 range
        for r in results:
            r['a_d22_range'] = max(r['a_profile_d22']) - min(r['a_profile_d22']) if r['a_profile_d22'] else 0
            r['a_d22_max'] = max(r['a_profile_d22']) if r['a_profile_d22'] else 0
            r['a_d22_min'] = min(r['a_profile_d22']) if r['a_profile_d22'] else 0

        print(f"\n  D22 A-profile range:")
        for t in range(0, 12):
            w_pass = sum(1 for r in working if r['a_d22_range'] <= t)
            n_pass = sum(1 for r in nonworking if r['a_d22_range'] <= t)
            if w_pass > 0:
                print(f"    range <= {t}: {w_pass}/{len(working)} working, {n_pass}/{len(nonworking)} non-working pass")
            if w_pass == len(working) and n_pass == 0:
                print(f"    *** SEPARATES ***")
                break

        # Test: sorted D22 A-profile
        print(f"\n  Sorted D22 A-profiles of working orbits:")
        for r in working:
            print(f"    N={r['N']}: {r['a_profile_d22']}")
        sep_d22, _, _ = analyze_separator(results, 'a_profile_d22', lambda r: tuple(r['a_profile_d22']))
        print(f"  D22 A-profile separates: {sep_d22}")

        # Check D22 A-profile overlap
        w_d22 = set(tuple(r['a_profile_d22']) for r in working)
        n_d22 = set(tuple(r['a_profile_d22']) for r in nonworking)
        overlap = w_d22 & n_d22
        print(f"  D22 A-profile overlap: {overlap}")

        # Test: A-variance <= threshold (combined with a_flat)
        print(f"\n  A-variance tests:")
        for r in results:
            a_vals = r['a_profile_d11']
            mean = sum(a_vals) / len(a_vals)
            r['a_var_exact'] = sum((x - mean)**2 for x in a_vals) / len(a_vals)
        for thresh_10 in range(0, 50):
            thresh = thresh_10 / 10.0
            w_pass = sum(1 for r in working if r['a_var_exact'] <= thresh)
            n_pass = sum(1 for r in nonworking if r['a_var_exact'] <= thresh)
            if w_pass > 0 and (w_pass == len(working) or n_pass == 0):
                print(f"    a_var <= {thresh:.1f}: {w_pass}/{len(working)} working, {n_pass}/{len(nonworking)} non-working")
                if w_pass == len(working) and n_pass == 0:
                    print(f"    *** SEPARATES ***")
                    break

        # Test: max_A_D11 combined with a_var
        print(f"\n  Combined max_A + a_var tests:")
        for max_a in range(3, 10):
            for thresh_10 in range(0, 50):
                thresh = thresh_10 / 10.0
                w_pass = sum(1 for r in working if r['max_A_D11'] <= max_a and r['a_var_exact'] <= thresh)
                n_pass = sum(1 for r in nonworking if r['max_A_D11'] <= max_a and r['a_var_exact'] <= thresh)
                if w_pass == len(working) and n_pass == 0:
                    print(f"    max_A <= {max_a} AND a_var <= {thresh:.1f}: SEPARATES")
                    break

        # ---- Character sum: specific chi values ----
        print(f"\n  Individual character |S_j|^2 values:")
        # For p=23, p-1=22, chars of order 11 are j=2,4,6,...,20
        half_p = (p-1)//2
        # Look at specific character indices
        for j in range(1, p-1):
            order = (p-1) // math.gcd(j, p-1)
            if order <= 2:
                continue
            key_j = j
            w_vals_j = sorted(set(round(r['char_abs_sq'][j], 2) for r in working))
            n_vals_j = sorted(set(round(r['char_abs_sq'][j], 2) for r in nonworking))
            overlap_j = set(w_vals_j) & set(n_vals_j)
            if len(overlap_j) == 0 and len(w_vals_j) > 0:
                print(f"    j={j} (order {order}): SEPARATES! Working: {w_vals_j}, Non-working: {n_vals_j}")

        # ---- Check if all working orbits share same character power spectrum ----
        print(f"\n  Character power spectra (|S_j|^2 sorted) for working orbits:")
        for r in working:
            spectrum = tuple(sorted(round(x, 2) for x in r['char_abs_sq']))
            print(f"    N={r['N']}: {spectrum}")

        # Check if spectra are identical
        w_spectra = [tuple(sorted(round(x, 2) for x in r['char_abs_sq'])) for r in working]
        if len(set(w_spectra)) == 1:
            print(f"  ALL working orbits have IDENTICAL power spectrum!")
        else:
            print(f"  Working orbits have {len(set(w_spectra))} DISTINCT power spectra")

        # The key question: are there non-working orbits with the same spectrum?
        w_spectrum_set = set(w_spectra)
        n_with_same = sum(1 for r in nonworking if tuple(sorted(round(x, 2) for x in r['char_abs_sq'])) in w_spectrum_set)
        print(f"  Non-working orbits with same spectrum: {n_with_same}")

        # ---- Test: sum of A-values on specific subsets ----
        print(f"\n  Sum of A-values tests:")
        for r in results:
            r['a_sum_d11'] = sum(r['a_profile_d11'])
            r['a_sum_d22'] = sum(r['a_profile_d22'])
        sep_sum_d11, wv, nv = analyze_separator(results, 'a_sum_d11')
        sep_sum_d22, wv2, nv2 = analyze_separator(results, 'a_sum_d22')
        print(f"  Sum A on D11: {'SEP' if sep_sum_d11 else 'FAIL'} Working: {sorted(set(wv))} Non-working: {sorted(set(nv))}")
        print(f"  Sum A on D22: {'SEP' if sep_sum_d22 else 'FAIL'} Working: {sorted(set(wv2))} Non-working: {sorted(set(nv2))}")

        # ---- Necklace autocorrelation ----
        print(f"\n  Necklace autocorrelation:")
        for r in results:
            n_bits = len(r['necklace'])
            bits = [int(c) for c in r['necklace']]
            mean_b = sum(bits) / n_bits
            autocorr = []
            for lag in range(1, n_bits):
                c = sum((bits[i] - mean_b) * (bits[(i + lag) % n_bits] - mean_b) for i in range(n_bits)) / n_bits
                autocorr.append(round(c, 4))
            r['necklace_autocorr'] = autocorr
            r['necklace_autocorr_max'] = max(abs(c) for c in autocorr)
            r['necklace_autocorr_sum'] = sum(c**2 for c in autocorr)

        sep_ac, wv, nv = analyze_separator(results, 'necklace_autocorr_sum', lambda r: round(r['necklace_autocorr_sum'], 4))
        print(f"  Autocorr energy: {'SEP' if sep_ac else 'FAIL'}")
        print(f"    Working: {sorted(set(round(r['necklace_autocorr_sum'], 4) for r in working))}")
        print(f"    Non-working: [{min(round(r['necklace_autocorr_sum'], 4) for r in nonworking)}, {max(round(r['necklace_autocorr_sum'], 4) for r in nonworking)}]")

        # ---- A-profile L-infinity norm (deviation from flat) ----
        print(f"\n  A-profile L_inf deviation from flat value:")
        # The flat A-value for D11 would be |D11|(|D11|-1)/(p-1)
        for r in results:
            d11_size = len(r['d11'])
            expected = d11_size * (d11_size - 1) / (p - 1)
            r['a_linf'] = max(abs(a - expected) for a in r['a_profile_d11'])
            r['a_l2'] = math.sqrt(sum((a - expected)**2 for a in r['a_profile_d11']))

        for thresh_10 in range(0, 50):
            thresh = thresh_10 / 10.0
            w_pass = sum(1 for r in working if r['a_linf'] <= thresh)
            n_pass = sum(1 for r in nonworking if r['a_linf'] <= thresh)
            if w_pass > 0:
                print(f"    L_inf <= {thresh:.1f}: {w_pass}/{len(working)} working, {n_pass}/{len(nonworking)} non-working")
            if w_pass == len(working) and n_pass == 0:
                print(f"    *** SEPARATES ***")
                break

        for thresh_10 in range(0, 100):
            thresh = thresh_10 / 10.0
            w_pass = sum(1 for r in working if r['a_l2'] <= thresh)
            n_pass = sum(1 for r in nonworking if r['a_l2'] <= thresh)
            if w_pass > 0 and w_pass == len(working) and n_pass == 0:
                print(f"    L_2 <= {thresh:.1f}: SEPARATES")
                break

        # ---- Check the full character spectrum (not just grouped by order) ----
        print(f"\n  Full character spectrum analysis:")
        # For each j, compute the sorted list of |S_j|^2 across D11 elements
        # Specifically check: Σ_j |S_j|^4 (fourth moment of power spectrum)
        for r in results:
            r['spectral_fourth'] = sum(x**2 for x in r['char_abs_sq'])
            r['spectral_max'] = max(r['char_abs_sq'])
            # Entropy of normalized spectrum
            total = sum(r['char_abs_sq'])
            if total > 0:
                probs = [x / total for x in r['char_abs_sq'] if x > 0]
                r['spectral_entropy'] = -sum(p * math.log(p) for p in probs) if probs else 0
            else:
                r['spectral_entropy'] = 0

        sep_4th, wv, nv = analyze_separator(results, 'spectral_fourth', lambda r: round(r['spectral_fourth'], 1))
        print(f"  Spectral 4th moment: {'SEP' if sep_4th else 'FAIL'}")
        print(f"    Working: {sorted(set(round(r['spectral_fourth'], 1) for r in working))}")
        if len(nonworking) <= 20:
            print(f"    Non-working: {sorted(set(round(r['spectral_fourth'], 1) for r in nonworking))}")
        else:
            print(f"    Non-working range: [{min(round(r['spectral_fourth'], 1) for r in nonworking)}, {max(round(r['spectral_fourth'], 1) for r in nonworking)}]")

    # ---- EXHAUSTIVE PAIRWISE COMBINATION SEARCH for p=23 ----
    print(f"\n{'='*80}")
    print("EXHAUSTIVE TWO-INVARIANT COMBINATION SEARCH (p=23)")
    print(f"{'='*80}")

    results_23 = all_results['23']
    working_23 = [r for r in results_23 if r['N'] > 0]
    nonworking_23 = [r for r in results_23 if r['N'] == 0]

    # Build feature vectors for each orbit
    def get_features(r):
        feats = {
            'max_A': r['max_A_D11'],
            'a_flat': int(r['is_a_flat']),
            'max_gap': r['max_gap'],
            'max_run': r['max_run'],
            'num_gaps': r['num_gaps'],
            'num_runs': r['num_runs'],
            'a_range': r['a_range'],
            'a_var': round(r.get('a_var_exact', r['a_var']), 4),
        }
        # Add char sums
        for order_str, vals in r['char_by_order'].items():
            order = int(order_str)
            if order <= 2:
                continue
            feats[f'char_max_{order}'] = round(max(vals), 2)
            feats[f'char_min_{order}'] = round(min(vals), 2)
        # Necklace DFT
        feats['dft_max'] = round(max(r['necklace_dft_sq'][1:]), 2)
        feats['dft_sum'] = round(sum(r['necklace_dft_sq'][1:]), 2)
        # A profile stats
        d11_vals = r['a_profile_d11']
        feats['a_d11_min'] = min(d11_vals)
        feats['a_d11_max'] = max(d11_vals)
        d22_vals = r['a_profile_d22']
        feats['a_d22_min'] = min(d22_vals)
        feats['a_d22_max'] = max(d22_vals)
        feats['a_d22_range'] = max(d22_vals) - min(d22_vals)
        return feats

    all_feats = [(get_features(r), r['N'] > 0) for r in results_23]
    feat_names = sorted(all_feats[0][0].keys())

    # For each pair of features, try threshold combinations
    found_separators = []
    for i, f1 in enumerate(feat_names):
        for f2 in feat_names[i:]:
            w_pairs = set((feats[f1], feats[f2]) for feats, is_w in all_feats if is_w)
            n_pairs = set((feats[f1], feats[f2]) for feats, is_w in all_feats if not is_w)
            if len(w_pairs & n_pairs) == 0:
                found_separators.append((f1, f2, w_pairs, n_pairs))

    if found_separators:
        print(f"\nFound {len(found_separators)} separating feature pairs:")
        for f1, f2, w_pairs, n_pairs in found_separators[:20]:
            print(f"  ({f1}, {f2})")
            print(f"    Working values: {sorted(w_pairs)}")
    else:
        print(f"\nNo pair of scalar features perfectly separates working from non-working!")

    # Check sorted A-profile on D11
    print(f"\nSorted D11 A-profiles for all orbits:")
    for r in sorted(results_23, key=lambda x: (-x['N'], x['necklace'])):
        marker = "***" if r['N'] > 0 else "   "
        print(f"  {marker} N={r['N']:>4d} {r['necklace']} A_D11={r['a_profile_d11']} A_D22={r['a_profile_d22']}")

    # Check if D11 A-profile + D22 A-profile together separate
    # Already know combined separates. But what's the minimal distinguishing info?
    # Check: just the (min, max) of each profile
    def minmax_combo(r):
        return (min(r['a_profile_d11']), max(r['a_profile_d11']),
                min(r['a_profile_d22']), max(r['a_profile_d22']))
    sep_mm, wv, nv = analyze_separator(results_23, None, minmax_combo)
    print(f"\n  (min_D11, max_D11, min_D22, max_D22) separates: {sep_mm}")
    if sep_mm:
        print(f"    Working: {sorted(set(wv))}")

    # Check: second-smallest and second-largest A-values
    def profile_signature(r):
        d11 = r['a_profile_d11']
        d22 = r['a_profile_d22']
        return (d11[0], d11[1], d11[-2], d11[-1], d22[0], d22[1], d22[-2], d22[-1])
    sep_sig, wv, nv = analyze_separator(results_23, None, profile_signature)
    print(f"\n  Profile signature (first 2 + last 2 of each) separates: {sep_sig}")
    if sep_sig:
        print(f"    Working: {sorted(set(wv))}")

    # The KEY observation: look at D22 profile patterns of working vs non-working
    # Working D22 profiles all have specific structure
    print(f"\n  D22 A-profile sorted, working orbits:")
    for r in working_23:
        print(f"    N={r['N']:>4d}: {r['a_profile_d22']}")
    print(f"  D22 A-profile sorted, NON-working orbits that overlap:")
    w_d22 = set(tuple(r['a_profile_d22']) for r in working_23)
    for r in nonworking_23:
        if tuple(r['a_profile_d22']) in w_d22:
            print(f"    N={r['N']:>4d}: {r['a_profile_d22']}  necklace={r['necklace']}")

    # Save all results
    save_data = {}
    for p_str, results in all_results.items():
        save_data[p_str] = []
        for r in results:
            # Convert for JSON serialization
            entry = dict(r)
            entry['pair_set'] = list(entry['pair_set']) if isinstance(entry['pair_set'], (set, frozenset)) else entry['pair_set']
            # Remove non-serializable entries
            for key in ['a_var_exact', 'a_linf', 'a_l2', 'a_d22_range', 'a_d22_max', 'a_d22_min',
                        'a_sum_d11', 'a_sum_d22', 'necklace_autocorr', 'necklace_autocorr_max',
                        'necklace_autocorr_sum', 'spectral_fourth', 'spectral_max', 'spectral_entropy']:
                entry.pop(key, None)
            save_data[p_str].append(entry)

    with open('/Users/stephenpadgett/Projects/math/ramsey-book-graphs/algebraic_invariants.json', 'w') as f:
        json.dump(save_data, f, indent=2)
    print(f"\nData saved to algebraic_invariants.json")


if __name__ == '__main__':
    main()
