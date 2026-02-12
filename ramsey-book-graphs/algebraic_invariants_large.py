#!/usr/bin/env python3
"""
Compute algebraic invariants for working and non-working orbits at p=31 and p=43.

For each orbit:
1. Primitive root and pair permutation
2. Necklace representation in pair-cycle order
3. Gap/run structure
4. Character sums for all multiplicative characters
5. Subgroup intersection sizes
6. QR structure within necklace ordering

Saves results to algebraic_invariants_large.json
"""

import json
import math
import cmath
from collections import Counter
from itertools import groupby


def primitive_root(p):
    """Find smallest primitive root mod p."""
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
    return None


def discrete_log(a, g, p):
    """Compute discrete log of a base g mod p (brute force, fine for small p)."""
    val = 1
    for k in range(p - 1):
        if val == a % p:
            return k
        val = (val * g) % p
    return None


def pair_index(d, p):
    """Map element d in Z_p* to its pair {d, p-d}. Return min(d, p-d)."""
    return min(d, p - d)


def compute_pair_permutation(g, p):
    """
    Compute how multiplication by g permutes the (p-1)/2 pairs {d, p-d}.
    Returns list of cycles, where each pair is identified by min(d, p-d).
    """
    num_pairs = (p - 1) // 2
    pairs = sorted(set(pair_index(d, p) for d in range(1, p)))
    assert len(pairs) == num_pairs

    visited = set()
    cycles = []
    for start in pairs:
        if start in visited:
            continue
        cycle = []
        current = start
        while current not in visited:
            visited.add(current)
            cycle.append(current)
            # Multiply representative by g
            next_val = (current * g) % p
            current = pair_index(next_val, p)
        if cycle:
            cycles.append(cycle)
    return cycles


def necklace_from_d11(d11_set, pair_order, p):
    """
    Express D11 as a binary necklace in the given pair ordering.
    pair_order is a list of pair representatives (min(d, p-d)) in cycle order.
    Returns binary string: '1' if pair is in D11, '0' otherwise.
    """
    bits = []
    for rep in pair_order:
        if rep in d11_set or (p - rep) in d11_set:
            bits.append('1')
        else:
            bits.append('0')
    return ''.join(bits)


def gap_run_structure(necklace):
    """
    Compute gap sizes (runs of 0s) and run sizes (runs of 1s).
    Returns dict with gap_sizes, run_sizes, num_gaps, num_runs.
    """
    groups = [(k, len(list(g))) for k, g in groupby(necklace)]
    # Handle circular: if first and last have same character, merge them
    if len(groups) > 1 and groups[0][0] == groups[-1][0]:
        merged_len = groups[0][1] + groups[-1][1]
        groups = [(groups[0][0], merged_len)] + groups[1:-1]

    gaps = [l for c, l in groups if c == '0']
    runs = [l for c, l in groups if c == '1']
    return {
        'gap_sizes': sorted(gaps, reverse=True),
        'run_sizes': sorted(runs, reverse=True),
        'num_gaps': len(gaps),
        'num_runs': len(runs),
        'max_gap': max(gaps) if gaps else 0,
        'max_run': max(runs) if runs else 0,
    }


def multiplicative_characters(p):
    """
    Return all multiplicative characters of Z_p*.
    Character chi_k(g^j) = exp(2*pi*i*k*j/(p-1)) for k=0,...,p-2.
    Returns dict: order -> list of (k, chi_k) where chi_k is a dict {element: value}.
    """
    g = primitive_root(p)
    n = p - 1
    # Precompute discrete logs
    dlogs = {}
    val = 1
    for j in range(n):
        dlogs[val] = j
        val = (val * g) % p

    chars_by_order = {}
    for k in range(n):
        order = n // math.gcd(k, n) if k > 0 else 1
        omega = cmath.exp(2j * cmath.pi * k / n)
        chi = {}
        for a in range(1, p):
            chi[a] = omega ** dlogs[a]
        if order not in chars_by_order:
            chars_by_order[order] = []
        chars_by_order[order].append((k, chi))

    return chars_by_order


def character_sum(d11_set, chi):
    """Compute sum_{d in D11} chi(d)."""
    return sum(chi[d] for d in d11_set)


def character_sum_magnitudes(d11_set, p):
    """
    Compute |S_chi|^2 for all multiplicative characters.
    Returns dict: order -> list of |S_chi|^2 values.
    """
    g = primitive_root(p)
    n = p - 1

    dlogs = {}
    val = 1
    for j in range(n):
        dlogs[val] = j
        val = (val * g) % p

    results = {}
    for k in range(n):
        order = n // math.gcd(k, n) if k > 0 else 1
        omega = cmath.exp(2j * cmath.pi * k / n)
        s = sum(omega ** dlogs[d] for d in d11_set)
        mag_sq = abs(s) ** 2
        if order not in results:
            results[order] = []
        results[order].append(round(mag_sq, 6))

    return results


def subgroup_intersections(d11_set, p):
    """
    Compute |D11 ∩ C| for cosets of various subgroups of Z_p*.
    Returns dict with intersection data for index-2, index-3, index-6 subgroups.
    """
    g = primitive_root(p)
    n = p - 1
    results = {}

    for index in [2, 3, 5, 6, 7]:
        if n % index != 0:
            continue
        # Subgroup of index `index`: elements g^(index*j) for j=0,...,n/index-1
        sub_size = n // index
        subgroup = set()
        val = 1
        gen = pow(g, index, p)
        for j in range(sub_size):
            subgroup.add(val)
            val = (val * gen) % p

        # All cosets
        coset_intersections = []
        for c in range(index):
            coset = set((pow(g, c, p) * s) % p for s in subgroup)
            inter = len(d11_set & coset)
            coset_intersections.append(inter)

        results[f'index_{index}'] = {
            'coset_intersections': sorted(coset_intersections, reverse=True),
            'subgroup_size': sub_size,
            'balance': max(coset_intersections) - min(coset_intersections),
        }

    return results


def qr_structure(d11_set, p):
    """
    Compute quadratic residue structure of D11.
    Returns: number of QR and QNR elements, and their positions in D11.
    """
    qr = set()
    for x in range(1, p):
        qr.add(pow(x, 2, p))

    d11_qr = d11_set & qr
    d11_qnr = d11_set - qr

    return {
        'num_qr': len(d11_qr),
        'num_qnr': len(d11_qnr),
        'qr_elements': sorted(d11_qr),
        'qnr_elements': sorted(d11_qnr),
        'qr_fraction': len(d11_qr) / len(d11_set) if d11_set else 0,
    }


def compute_A_values(d11_set, p):
    """Compute A(d) = |{x in D11 : x+d in D11}| for d in Z_p*."""
    m = p  # since m = 2n-1 = p
    a_values = {}
    for d in range(1, p):
        count = sum(1 for x in d11_set if (x + d) % p in d11_set)
        a_values[d] = count
    return a_values


def compute_all_invariants(d11_list, p, label):
    """Compute all invariants for a D11 set."""
    d11_set = set(d11_list)
    g = primitive_root(p)
    n = p - 1
    num_pairs = n // 2

    # Pair permutation cycles
    pair_cycles = compute_pair_permutation(g, p)
    cycle_lengths = sorted([len(c) for c in pair_cycles], reverse=True)

    # Flatten cycles into a single pair ordering
    pair_order = []
    for cycle in pair_cycles:
        pair_order.extend(cycle)

    # Necklace
    necklace = necklace_from_d11(d11_set, pair_order, p)

    # Gap/run structure
    gap_run = gap_run_structure(necklace)

    # Character sum magnitudes
    char_sums = character_sum_magnitudes(d11_set, p)

    # Simplify: for each order, report sum and max of |S_chi|^2
    char_summary = {}
    for order, vals in char_sums.items():
        char_summary[str(order)] = {
            'count': len(vals),
            'values': vals,
            'sum': round(sum(vals), 4),
            'max': round(max(vals), 4),
            'min': round(min(vals), 4),
        }

    # Subgroup intersections
    sub_ints = subgroup_intersections(d11_set, p)

    # QR structure
    qr = qr_structure(d11_set, p)

    # A-values
    a_vals = compute_A_values(d11_set, p)
    max_a = max(a_vals.values())
    a_distribution = dict(Counter(a_vals.values()))

    # Compute z-sum: z_d = (A(d) - E[A]) / sigma
    k = len(d11_set)  # |D11|
    E_A = k * (k - 1) / (p - 1)
    # Variance from the known formula
    var_A = (k * (p - 1 - k) * (k - 1) * (p - k)) / ((p - 1)**2 * (p - 2))
    sigma = math.sqrt(var_A) if var_A > 0 else 1
    z_values = {d: (a_vals[d] - E_A) / sigma for d in range(1, p)}
    z_sum = sum(z_values.values())

    # DFT of indicator: hat(1_D11)(k) = sum_{d in D11} omega^{dk}
    # This is related to character sums
    dft_magnitudes = []
    for freq in range(1, p):
        omega = cmath.exp(2j * cmath.pi * freq / p)
        s = sum(omega ** d for d in d11_set)
        dft_magnitudes.append(round(abs(s) ** 2, 6))

    # Flat spectrum measure: how uniform are DFT magnitudes?
    dft_mean = sum(dft_magnitudes) / len(dft_magnitudes) if dft_magnitudes else 0
    dft_var = sum((x - dft_mean)**2 for x in dft_magnitudes) / len(dft_magnitudes) if dft_magnitudes else 0
    dft_max = max(dft_magnitudes) if dft_magnitudes else 0
    dft_min = min(dft_magnitudes) if dft_magnitudes else 0

    return {
        'label': label,
        'p': p,
        'D11': sorted(d11_set),
        'D11_size': len(d11_set),
        'max_A': max_a,
        'a_distribution': {str(k): v for k, v in sorted(a_distribution.items())},
        'z_sum': round(z_sum, 6),
        'primitive_root': g,
        'pair_cycle_lengths': cycle_lengths,
        'necklace': necklace,
        'gap_run': gap_run,
        'character_sums': char_summary,
        'subgroup_intersections': sub_ints,
        'qr_structure': qr,
        'dft_spectrum': {
            'mean_mag_sq': round(dft_mean, 6),
            'var_mag_sq': round(dft_var, 6),
            'max_mag_sq': round(dft_max, 6),
            'min_mag_sq': round(dft_min, 6),
            'flatness_ratio': round(dft_max / dft_mean, 6) if dft_mean > 0 else 0,
        },
    }


def load_working_orbits_p31():
    """Load working orbits at p=31 from c0_estimation_results.json."""
    with open('/Users/stephenpadgett/Projects/math/ramsey-book-graphs/c0_estimation_results.json') as f:
        data = json.load(f)

    working = []
    non_working = []
    for orbit in data['p31']['orbit_results']:
        d11 = orbit['canonical']
        info = {
            'D11': d11,
            'max_A': orbit['max_A'],
            'is_aflat': orbit['is_aflat'],
            'num_valid': orbit['num_valid'],
        }
        if orbit['num_valid'] > 0:
            working.append(info)
        else:
            non_working.append(info)

    return working, non_working


def load_working_orbits_p43():
    """Load working orbits at p=43 from near_flat_p43_results.json."""
    with open('/Users/stephenpadgett/Projects/math/ramsey-book-graphs/near_flat_p43_results.json') as f:
        data = json.load(f)

    working = []
    non_working = []
    for orbit in data['orbit_results']:
        d11 = orbit['canonical']
        info = {
            'D11': d11,
            'max_A': orbit['max_A'],
            'is_aflat': orbit['is_aflat'],
            'working': orbit['working'],
        }
        if orbit['working']:
            working.append(info)
        else:
            non_working.append(info)

    return working, non_working


def summarize_comparison(working_invs, nonworking_invs, p):
    """Compare invariants between working and non-working orbits."""
    summary = {'p': p}

    # Necklace comparison
    summary['working_necklaces'] = [inv['necklace'] for inv in working_invs]
    summary['working_max_A'] = [inv['max_A'] for inv in working_invs]

    # Gap/run comparison
    w_max_gaps = [inv['gap_run']['max_gap'] for inv in working_invs]
    w_max_runs = [inv['gap_run']['max_run'] for inv in working_invs]
    nw_max_gaps = [inv['gap_run']['max_gap'] for inv in nonworking_invs]
    nw_max_runs = [inv['gap_run']['max_run'] for inv in nonworking_invs]
    summary['gap_run'] = {
        'working_max_gap': {'min': min(w_max_gaps), 'max': max(w_max_gaps), 'mean': sum(w_max_gaps)/len(w_max_gaps)},
        'working_max_run': {'min': min(w_max_runs), 'max': max(w_max_runs), 'mean': sum(w_max_runs)/len(w_max_runs)},
        'nonworking_max_gap': {'min': min(nw_max_gaps), 'max': max(nw_max_gaps), 'mean': sum(nw_max_gaps)/len(nw_max_gaps)},
        'nonworking_max_run': {'min': min(nw_max_runs), 'max': max(nw_max_runs), 'mean': sum(nw_max_runs)/len(nw_max_runs)},
    }

    # DFT flatness comparison
    w_flat = [inv['dft_spectrum']['flatness_ratio'] for inv in working_invs]
    nw_flat = [inv['dft_spectrum']['flatness_ratio'] for inv in nonworking_invs]
    summary['dft_flatness'] = {
        'working': {'min': round(min(w_flat), 4), 'max': round(max(w_flat), 4), 'mean': round(sum(w_flat)/len(w_flat), 4)},
        'nonworking': {'min': round(min(nw_flat), 4), 'max': round(max(nw_flat), 4), 'mean': round(sum(nw_flat)/len(nw_flat), 4)},
    }

    # DFT max magnitude comparison
    w_dft_max = [inv['dft_spectrum']['max_mag_sq'] for inv in working_invs]
    nw_dft_max = [inv['dft_spectrum']['max_mag_sq'] for inv in nonworking_invs]
    summary['dft_max_mag_sq'] = {
        'working': {'min': round(min(w_dft_max), 4), 'max': round(max(w_dft_max), 4), 'mean': round(sum(w_dft_max)/len(w_dft_max), 4)},
        'nonworking': {'min': round(min(nw_dft_max), 4), 'max': round(max(nw_dft_max), 4), 'mean': round(sum(nw_dft_max)/len(nw_dft_max), 4)},
    }

    # QR balance comparison
    w_qr_frac = [inv['qr_structure']['qr_fraction'] for inv in working_invs]
    nw_qr_frac = [inv['qr_structure']['qr_fraction'] for inv in nonworking_invs]
    summary['qr_fraction'] = {
        'working': {'min': round(min(w_qr_frac), 4), 'max': round(max(w_qr_frac), 4), 'mean': round(sum(w_qr_frac)/len(w_qr_frac), 4)},
        'nonworking': {'min': round(min(nw_qr_frac), 4), 'max': round(max(nw_qr_frac), 4), 'mean': round(sum(nw_qr_frac)/len(nw_qr_frac), 4)},
    }

    # Subgroup balance comparison (index-2)
    def get_index2_balance(inv):
        si = inv.get('subgroup_intersections', {})
        if 'index_2' in si:
            return si['index_2']['balance']
        return None

    w_bal2 = [get_index2_balance(inv) for inv in working_invs if get_index2_balance(inv) is not None]
    nw_bal2 = [get_index2_balance(inv) for inv in nonworking_invs if get_index2_balance(inv) is not None]
    if w_bal2 and nw_bal2:
        summary['index2_balance'] = {
            'working': {'min': min(w_bal2), 'max': max(w_bal2), 'values': w_bal2},
            'nonworking': {'min': min(nw_bal2), 'max': max(nw_bal2), 'mean': sum(nw_bal2)/len(nw_bal2)},
        }

    # Character sum analysis: for each order, compare max |S_chi|^2
    # Focus on key orders
    key_orders_p31 = ['2', '3', '5', '6', '10', '15', '30']
    key_orders_p43 = ['2', '3', '6', '7', '14', '21', '42']
    key_orders = key_orders_p31 if p == 31 else key_orders_p43

    char_comparison = {}
    for order in key_orders:
        w_vals = []
        nw_vals = []
        for inv in working_invs:
            if order in inv['character_sums']:
                w_vals.append(inv['character_sums'][order]['max'])
        for inv in nonworking_invs:
            if order in inv['character_sums']:
                nw_vals.append(inv['character_sums'][order]['max'])
        if w_vals and nw_vals:
            char_comparison[f'order_{order}_max'] = {
                'working': {'min': round(min(w_vals), 4), 'max': round(max(w_vals), 4), 'mean': round(sum(w_vals)/len(w_vals), 4)},
                'nonworking': {'min': round(min(nw_vals), 4), 'max': round(max(nw_vals), 4), 'mean': round(sum(nw_vals)/len(nw_vals), 4)},
            }
    summary['character_sum_comparison'] = char_comparison

    # Z-sum comparison
    w_zsum = [inv['z_sum'] for inv in working_invs]
    nw_zsum = [inv['z_sum'] for inv in nonworking_invs]
    summary['z_sum'] = {
        'working': {'values': [round(z, 4) for z in w_zsum], 'mean': round(sum(w_zsum)/len(w_zsum), 4)},
        'nonworking': {'min': round(min(nw_zsum), 4), 'max': round(max(nw_zsum), 4), 'mean': round(sum(nw_zsum)/len(nw_zsum), 4)},
    }

    # Subgroup index-3, index-6 balance
    for idx_label in ['index_3', 'index_5', 'index_6', 'index_7']:
        def get_balance(inv, label=idx_label):
            si = inv.get('subgroup_intersections', {})
            if label in si:
                return si[label]['balance']
            return None

        w_b = [get_balance(inv) for inv in working_invs if get_balance(inv) is not None]
        nw_b = [get_balance(inv) for inv in nonworking_invs if get_balance(inv) is not None]
        if w_b and nw_b:
            summary[f'{idx_label}_balance'] = {
                'working': {'values': w_b, 'mean': round(sum(w_b)/len(w_b), 4)},
                'nonworking': {'min': min(nw_b), 'max': max(nw_b), 'mean': round(sum(nw_b)/len(nw_b), 4)},
            }

    return summary


def main():
    results = {}

    # =========== p = 31 ===========
    print("=" * 60)
    print("Processing p = 31")
    print("=" * 60)

    p31_working, p31_nonworking = load_working_orbits_p31()
    print(f"  Working orbits: {len(p31_working)}")
    print(f"  Non-working orbits: {len(p31_nonworking)}")

    g31 = primitive_root(31)
    print(f"  Primitive root: {g31}")
    pair_cycles_31 = compute_pair_permutation(g31, 31)
    print(f"  Pair cycle lengths: {sorted([len(c) for c in pair_cycles_31], reverse=True)}")

    # Compute invariants for working orbits
    p31_working_invs = []
    for i, orbit in enumerate(p31_working):
        print(f"  Computing working orbit {i+1}/{len(p31_working)}: D11={orbit['D11'][:5]}...")
        inv = compute_all_invariants(orbit['D11'], 31, f"p31_working_{i}")
        inv['num_valid'] = orbit['num_valid']
        inv['is_aflat'] = orbit['is_aflat']
        p31_working_invs.append(inv)

    # Sample of non-working orbits (at least 16: all 16 A-flat non-working + some NF)
    # Take all A-flat non-working and first few NF non-working
    p31_aflat_nw = [o for o in p31_nonworking if o['is_aflat']]
    p31_nf_nw = [o for o in p31_nonworking if not o['is_aflat']]
    sample_nw = p31_aflat_nw[:16] + p31_nf_nw[:10]

    p31_nonworking_invs = []
    for i, orbit in enumerate(sample_nw):
        print(f"  Computing non-working orbit {i+1}/{len(sample_nw)}: D11={orbit['D11'][:5]}...")
        inv = compute_all_invariants(orbit['D11'], 31, f"p31_nonworking_{i}")
        inv['is_aflat'] = orbit['is_aflat']
        p31_nonworking_invs.append(inv)

    # Summary comparison
    p31_summary = summarize_comparison(p31_working_invs, p31_nonworking_invs, 31)

    results['p31'] = {
        'working_orbits': p31_working_invs,
        'nonworking_sample': p31_nonworking_invs,
        'summary': p31_summary,
    }

    # Print key findings for p=31
    print("\n  --- p=31 Summary ---")
    print(f"  Working necklaces: {p31_summary['working_necklaces']}")
    print(f"  Working z-sums: {p31_summary['z_sum']['working']['values']}")
    print(f"  DFT flatness - working: {p31_summary['dft_flatness']['working']}")
    print(f"  DFT flatness - nonworking: {p31_summary['dft_flatness']['nonworking']}")
    print(f"  QR fraction - working: {p31_summary['qr_fraction']['working']}")
    print(f"  QR fraction - nonworking: {p31_summary['qr_fraction']['nonworking']}")

    # =========== p = 43 ===========
    print("\n" + "=" * 60)
    print("Processing p = 43")
    print("=" * 60)

    p43_working, p43_nonworking = load_working_orbits_p43()
    print(f"  Working orbits: {len(p43_working)}")
    print(f"  Non-working orbits: {len(p43_nonworking)}")

    g43 = primitive_root(43)
    print(f"  Primitive root: {g43}")
    pair_cycles_43 = compute_pair_permutation(g43, 43)
    print(f"  Pair cycle lengths: {sorted([len(c) for c in pair_cycles_43], reverse=True)}")

    # Compute invariants for working orbits
    p43_working_invs = []
    for i, orbit in enumerate(p43_working):
        print(f"  Computing working orbit {i+1}/{len(p43_working)}: D11={orbit['D11'][:5]}...")
        inv = compute_all_invariants(orbit['D11'], 43, f"p43_working_{i}")
        inv['is_aflat'] = orbit['is_aflat']
        p43_working_invs.append(inv)

    # Sample of non-working orbits (all A-flat + some NF)
    p43_aflat_nw = [o for o in p43_nonworking if o['is_aflat']]
    p43_nf_nw = [o for o in p43_nonworking if not o['is_aflat']]
    sample_nw_43 = p43_aflat_nw[:15] + p43_nf_nw[:10]

    p43_nonworking_invs = []
    for i, orbit in enumerate(sample_nw_43):
        print(f"  Computing non-working orbit {i+1}/{len(sample_nw_43)}: D11={orbit['D11'][:5]}...")
        inv = compute_all_invariants(orbit['D11'], 43, f"p43_nonworking_{i}")
        inv['is_aflat'] = orbit['is_aflat']
        p43_nonworking_invs.append(inv)

    # Summary comparison
    p43_summary = summarize_comparison(p43_working_invs, p43_nonworking_invs, 43)

    results['p43'] = {
        'working_orbits': p43_working_invs,
        'nonworking_sample': p43_nonworking_invs,
        'summary': p43_summary,
    }

    # Print key findings for p=43
    print("\n  --- p=43 Summary ---")
    print(f"  Working necklaces: {p43_summary['working_necklaces']}")
    print(f"  Working z-sums: {p43_summary['z_sum']['working']['values']}")
    print(f"  DFT flatness - working: {p43_summary['dft_flatness']['working']}")
    print(f"  DFT flatness - nonworking: {p43_summary['dft_flatness']['nonworking']}")
    print(f"  QR fraction - working: {p43_summary['qr_fraction']['working']}")
    print(f"  QR fraction - nonworking: {p43_summary['qr_fraction']['nonworking']}")

    # =========== Cross-prime comparison ===========
    print("\n" + "=" * 60)
    print("Cross-prime comparison")
    print("=" * 60)

    # Check: do working orbits tend to have lower DFT max?
    for label, w_invs, nw_invs in [("p31", p31_working_invs, p31_nonworking_invs),
                                     ("p43", p43_working_invs, p43_nonworking_invs)]:
        w_dft = sorted([inv['dft_spectrum']['max_mag_sq'] for inv in w_invs])
        nw_dft = sorted([inv['dft_spectrum']['max_mag_sq'] for inv in nw_invs])
        print(f"\n  {label} DFT max |hat(1)|^2:")
        print(f"    Working:     {w_dft}")
        print(f"    Non-working: {nw_dft[:10]}... (showing first 10)")

    # Check: subgroup intersection patterns
    for label, w_invs, nw_invs in [("p31", p31_working_invs, p31_nonworking_invs),
                                     ("p43", p43_working_invs, p43_nonworking_invs)]:
        print(f"\n  {label} Index-2 subgroup (QR/QNR) intersections:")
        for inv in w_invs:
            if 'index_2' in inv['subgroup_intersections']:
                ci = inv['subgroup_intersections']['index_2']['coset_intersections']
                print(f"    Working {inv['label']}: {ci}")
        for inv in nw_invs[:5]:
            if 'index_2' in inv['subgroup_intersections']:
                ci = inv['subgroup_intersections']['index_2']['coset_intersections']
                print(f"    Non-working {inv['label']}: {ci}")

    # Save results
    output_path = '/Users/stephenpadgett/Projects/math/ramsey-book-graphs/algebraic_invariants_large.json'
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nResults saved to {output_path}")

    # Print a final pattern analysis
    print("\n" + "=" * 60)
    print("PATTERN ANALYSIS")
    print("=" * 60)

    # Check for simple separating conditions
    print("\nChecking potential separating conditions...")

    for label, w_invs, nw_invs in [("p31", p31_working_invs, p31_nonworking_invs),
                                     ("p43", p43_working_invs, p43_nonworking_invs)]:
        p = w_invs[0]['p']
        print(f"\n  --- {label} (p={p}) ---")

        # 1. DFT flatness threshold
        w_flat = [inv['dft_spectrum']['flatness_ratio'] for inv in w_invs]
        nw_flat = [inv['dft_spectrum']['flatness_ratio'] for inv in nw_invs]
        print(f"  DFT flatness: W=[{min(w_flat):.3f}, {max(w_flat):.3f}], NW=[{min(nw_flat):.3f}, {max(nw_flat):.3f}]")

        # 2. Number of gaps/runs
        w_ngaps = [inv['gap_run']['num_gaps'] for inv in w_invs]
        nw_ngaps = [inv['gap_run']['num_gaps'] for inv in nw_invs]
        print(f"  Num gaps: W={sorted(w_ngaps)}, NW range=[{min(nw_ngaps)}, {max(nw_ngaps)}]")

        # 3. Max gap
        w_mgap = [inv['gap_run']['max_gap'] for inv in w_invs]
        nw_mgap = [inv['gap_run']['max_gap'] for inv in nw_invs]
        print(f"  Max gap: W={sorted(w_mgap)}, NW range=[{min(nw_mgap)}, {max(nw_mgap)}]")

        # 4. QR fraction
        w_qr = [inv['qr_structure']['qr_fraction'] for inv in w_invs]
        nw_qr = [inv['qr_structure']['qr_fraction'] for inv in nw_invs]
        print(f"  QR fraction: W=[{min(w_qr):.4f}, {max(w_qr):.4f}], NW=[{min(nw_qr):.4f}, {max(nw_qr):.4f}]")

        # 5. Quadratic character sum (order 2)
        w_quad = [inv['character_sums'].get('2', {}).get('max', 0) for inv in w_invs]
        nw_quad = [inv['character_sums'].get('2', {}).get('max', 0) for inv in nw_invs]
        print(f"  |S_chi2|^2 (Legendre): W={[round(x,2) for x in w_quad]}, NW range=[{min(nw_quad):.2f}, {max(nw_quad):.2f}]")

        # 6. Index-2 balance
        w_b2 = [inv['subgroup_intersections'].get('index_2', {}).get('balance', None) for inv in w_invs]
        nw_b2 = [inv['subgroup_intersections'].get('index_2', {}).get('balance', None) for inv in nw_invs]
        w_b2 = [x for x in w_b2 if x is not None]
        nw_b2 = [x for x in nw_b2 if x is not None]
        if w_b2 and nw_b2:
            print(f"  Index-2 balance: W={sorted(w_b2)}, NW range=[{min(nw_b2)}, {max(nw_b2)}]")

        # 7. Check all character orders for separation
        all_orders = set()
        for inv in w_invs + nw_invs:
            all_orders.update(inv['character_sums'].keys())

        for order in sorted(all_orders, key=lambda x: int(x)):
            w_max = [inv['character_sums'].get(order, {}).get('max', 0) for inv in w_invs]
            nw_max = [inv['character_sums'].get(order, {}).get('max', 0) for inv in nw_invs]
            # Check if there's separation
            if max(w_max) < min(nw_max) or min(w_max) > max(nw_max):
                print(f"  ** SEPARATION at order {order}: W=[{min(w_max):.4f}, {max(w_max):.4f}], NW=[{min(nw_max):.4f}, {max(nw_max):.4f}]")
            elif max(w_max) < sum(nw_max)/len(nw_max) * 0.7 or min(w_max) > sum(nw_max)/len(nw_max) * 1.3:
                print(f"  * Partial separation at order {order}: W=[{min(w_max):.4f}, {max(w_max):.4f}], NW mean={sum(nw_max)/len(nw_max):.4f}")


if __name__ == '__main__':
    main()
