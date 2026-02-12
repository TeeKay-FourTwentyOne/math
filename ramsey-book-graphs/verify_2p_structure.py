#!/usr/bin/env python3
"""
Verify the 2p orbit structure of valid D12 sets.

Key question: WHY is N(D11) always a multiple of 2p?

Two hypotheses:
(a) Additive shifts: D12 -> D12+c preserves validity for the SAME D11 (gives factor p),
    plus negation D12 -> -D12 (gives factor 2).
(b) The factor comes from the circulant graph structure / multiplicative orbits acting
    on (D11, D12) simultaneously.

This script tests hypothesis (a) by:
1. Taking a known valid (D11, D12) pair
2. Trying all p additive shifts D12+c for c=0,...,p-1
3. Checking which shifts are still valid for the SAME D11

It also tests:
- Negation: D12 -> {(-x) mod p : x in D12}
- Which constraint types break under shifts (V1V1, V2V2, V1V2)
"""

import json
import sys
import os

sys.path.insert(0, os.path.dirname(__file__))
from ramsey_core import BlockCirculantGraph, verify_construction, Delta, Sigma


def shift_set(S, c, m):
    """Shift every element of S by c mod m."""
    return {(x + c) % m for x in S}


def negate_set(S, m):
    """Negate every element of S mod m."""
    return {(-x) % m for x in S}


def detailed_check(D11, D12, D22, n, p):
    """Check validity and return detailed violation info."""
    G = BlockCirculantGraph(n=n, D11=set(D11), D12=set(D12), D22=set(D22))
    result = verify_construction(G)
    return result


def analyze_shift_invariance(D11_list, D12_list, p, n, label=""):
    """Test all additive shifts of D12 for a given D11."""
    D11 = set(D11_list)
    D12_orig = set(D12_list)
    D22 = set(range(1, p)) - D11

    print(f"\n{'='*70}")
    print(f"Analyzing {label}: p={p}, n={n}")
    print(f"|D11|={len(D11)}, |D12|={len(D12_orig)}, |D22|={len(D22)}")
    print(f"D11 symmetric: {all((-x) % p in D11 for x in D11)}")
    print(f"0 in D12: {0 in D12_orig}")

    # Verify original is valid
    result_orig = detailed_check(D11, D12_orig, D22, n, p)
    print(f"\nOriginal valid: {result_orig.valid}")
    assert result_orig.valid, "Original pair must be valid!"

    # Test all additive shifts
    valid_shifts = []
    shift_violations = {}

    for c in range(p):
        D12_shifted = shift_set(D12_orig, c, p)
        result = detailed_check(D11, D12_shifted, D22, n, p)
        if result.valid:
            valid_shifts.append(c)
        else:
            # Categorize violations
            v_types = set()
            for vtype, d, excess in result.violations:
                v_types.add(vtype.split('_')[0] + '_' + vtype.split('_')[1])
            shift_violations[c] = (len(result.violations), v_types)

    print(f"\nAdditive shift results:")
    print(f"  Valid shifts: {len(valid_shifts)} out of {p}")
    print(f"  Valid shift values: {valid_shifts}")

    # Check: is shift 0 always valid? (sanity)
    assert 0 in valid_shifts, "Shift 0 must be valid!"

    # Check which constraint types break
    v1v1_breaks = 0
    v2v2_breaks = 0
    v1v2_breaks = 0
    for c, (nv, vtypes) in shift_violations.items():
        if any('V1V1' in vt for vt in vtypes):
            v1v1_breaks += 1
        if any('V2V2' in vt for vt in vtypes):
            v2v2_breaks += 1
        if any('V1V2' in vt for vt in vtypes):
            v1v2_breaks += 1

    n_invalid = p - len(valid_shifts)
    if n_invalid > 0:
        print(f"\n  Among {n_invalid} invalid shifts:")
        print(f"    V1V1 constraint breaks: {v1v1_breaks}")
        print(f"    V2V2 constraint breaks: {v2v2_breaks}")
        print(f"    V1V2 constraint breaks: {v1v2_breaks}")

    # Test negation
    D12_neg = negate_set(D12_orig, p)
    result_neg = detailed_check(D11, D12_neg, D22, n, p)
    print(f"\nNegation D12 -> -D12: valid={result_neg.valid}")

    # Test negation of each valid shift
    neg_valid_count = 0
    for c in valid_shifts:
        D12_neg_shift = negate_set(shift_set(D12_orig, c, p), p)
        result = detailed_check(D11, D12_neg_shift, D22, n, p)
        if result.valid:
            neg_valid_count += 1
    print(f"  Negation of valid shifts also valid: {neg_valid_count}/{len(valid_shifts)}")

    # Analyze B-values (autocorrelation) invariance under shifts
    print(f"\n--- B-value (autocorrelation) analysis ---")
    B_orig = [Delta(D12_orig, D12_orig, d, p) for d in range(p)]
    B_shift1 = [Delta(shift_set(D12_orig, 1, p), shift_set(D12_orig, 1, p), d, p) for d in range(p)]
    b_invariant = (B_orig == B_shift1)
    print(f"  B(d) invariant under shift by 1: {b_invariant}")

    # Analyze Sigma (cross-term) change under shifts
    print(f"\n--- Sigma (cross-term) analysis ---")
    print(f"  Sigma(D11, D12, d) = #{'{'}s in D12 : d-s in D11{'}'}")

    # Show how Sigma changes for a few d values
    for d_test in [0, 1, 2]:
        sigma_vals = []
        for c in range(min(5, p)):
            D12_c = shift_set(D12_orig, c, p)
            s = Sigma(D11, D12_c, d_test, p)
            sigma_vals.append(s)
        print(f"  Sigma(D11, D12+c, d={d_test}) for c=0..{min(4,p-1)}: {sigma_vals}")

    return valid_shifts


def analyze_multiplicative_orbits(D11_list, D12_list, p, n, label=""):
    """Analyze multiplicative orbit structure."""
    D11 = set(D11_list)
    D12_orig = set(D12_list)

    print(f"\n{'='*70}")
    print(f"Multiplicative orbit analysis: {label}")

    # Find a primitive root mod p
    def primitive_root(p):
        for g in range(2, p):
            seen = set()
            x = 1
            for _ in range(p - 1):
                seen.add(x)
                x = (x * g) % p
            if len(seen) == p - 1:
                return g
        return None

    g = primitive_root(p)
    print(f"Primitive root: g={g}")

    # Multiplicative orbit of D11: g*D11 = {g*x mod p : x in D11}
    def mult_set(S, g, p):
        return {(g * x) % p for x in S}

    # Check if D11 is closed under multiplication by QRs
    qr = set()
    x = 1
    for _ in range((p - 1) // 2):
        qr.add(x)
        x = (x * g * g) % p

    print(f"Quadratic residues: {sorted(qr)}")
    print(f"D11 = QR? {D11 == qr}")

    # Check orbit of (D11, D12) under multiplication by g
    print(f"\nOrbit of (D11, D12) under multiplication by g={g}:")
    D11_cur = set(D11)
    D12_cur = set(D12_orig)
    D22 = set(range(1, p)) - D11
    orbit_size = 0
    orbit_D12s = []

    for i in range(p - 1):
        D11_next = mult_set(D11_cur, g, p)
        D12_next = mult_set(D12_cur, g, p)

        if D11_next == D11:
            # D11 is invariant under this power of g
            # Check if D12_next is valid for D11
            result = detailed_check(D11, D12_next, D22, n, p)
            orbit_D12s.append((i + 1, D12_next, result.valid))

        D11_cur = D11_next
        D12_cur = D12_next

        if D11_next == D11 and D12_next == D12_orig:
            orbit_size = i + 1
            break

    if orbit_size == 0:
        orbit_size = p - 1

    print(f"  Full orbit size (D11, D12): {orbit_size}")

    # D11 stabilizer: powers of g that fix D11
    stabilizer = []
    g_power = 1
    for i in range(1, p):
        g_power = (g_power * g) % p
        if mult_set(D11, g_power, p) == D11:
            stabilizer.append(i)
    print(f"  D11 stabilizer (powers of g fixing D11): {stabilizer}")
    print(f"  Stabilizer size: {len(stabilizer)}")
    print(f"  Orbit of D11 size: {(p-1) // len(stabilizer) if stabilizer else p-1}")

    # For each element of the stabilizer, check if g^i * D12 is valid
    valid_mult = []
    for i in stabilizer:
        g_i = pow(g, i, p)
        D12_mult = mult_set(D12_orig, g_i, p)
        result = detailed_check(D11, D12_mult, D22, n, p)
        valid_mult.append((i, g_i, result.valid))
        if not result.valid:
            # How many violations?
            n_viol = len(result.violations)

    print(f"\n  D12 images under D11-stabilizer:")
    n_valid_mult = sum(1 for _, _, v in valid_mult if v)
    print(f"    Valid: {n_valid_mult}/{len(stabilizer)}")
    for i, g_i, v in valid_mult:
        print(f"    g^{i} = {g_i}: valid={v}")

    return stabilizer, valid_mult


def analyze_circulant_symmetry(D11_list, D12_list, p, n, label=""):
    """
    Analyze the graph-level symmetry that produces the factor of p.

    In a 2-block circulant graph on Z_p x {1,2}, the map
      sigma_c: (v, block) -> ((v+c) mod p, block)
    is a graph automorphism for every c in Z_p.

    This does NOT change D11, D12, D22 -- it permutes the VERTICES.
    The number of valid D12 for a given D11 is N(D11).
    The factor of p in N(D11) comes from a different source.

    Key insight: when we count valid D12 by brute force, we count
    each *graph* multiple times if the graph has additional symmetries.
    But D12 -> D12+c gives a DIFFERENT graph (different cross edges),
    so different D12 sets CAN give the same graph only if D12+c = D12
    (i.e., D12 is a union of cosets of some subgroup).

    Actually, the factor of p might simply not come from D12 shifts at all.
    Let's verify by looking at the exact counts from the exhaustive enumeration.
    """
    print(f"\n{'='*70}")
    print(f"Circulant symmetry analysis: {label}")
    print(f"The vertex permutation sigma_c : v -> v+c is a graph automorphism.")
    print(f"But this does NOT change (D11, D12, D22) -- it acts on vertex labels.")
    print(f"So the factor of p in N(D11) is NOT from vertex automorphisms.")


def main():
    # =========================================================
    # Test 1: p=31, n=16 (from p31_solutions.json)
    # =========================================================
    with open(os.path.join(os.path.dirname(__file__), 'p31_solutions.json')) as f:
        data = json.load(f)

    # Use first solution
    sol = data['solutions'][0]
    D11_31 = sol['D11']
    D12_31 = sol['D12']
    p31 = 31
    n31 = 16

    valid_shifts_31 = analyze_shift_invariance(
        D11_31, D12_31, p31, n31, label="p=31 solution #0"
    )

    stab_31, mult_31 = analyze_multiplicative_orbits(
        D11_31, D12_31, p31, n31, label="p=31 solution #0"
    )

    # Use the solution with highest valid_d12_count (7) for comparison
    sol_best = data['solutions'][2]  # valid_d12_count=7
    D11_31b = sol_best['D11']
    D12_31b = sol_best['D12']

    valid_shifts_31b = analyze_shift_invariance(
        D11_31b, D12_31b, p31, n31, label="p=31 solution #2 (N=7)"
    )

    # =========================================================
    # Test 2: Paley construction at p=53 (n=27), where D11=D12=QR
    # =========================================================
    # For Paley: D11 = D12 = QR(p), D22 = QNR(p)
    # Under multiplicative shifts g*D12 = g*QR = QR (QR is closed under mult)
    # So ALL (p-1)/2 multiplicative stabilizer elements give valid D12
    p53 = 53
    n53 = 27
    # QR mod 53
    qr53 = set()
    g53 = 2  # primitive root mod 53
    x = 1
    for _ in range(26):
        qr53.add(x)
        x = (x * g53 * g53) % 53
    D11_53 = sorted(qr53)
    D12_53 = sorted(qr53)

    valid_shifts_53 = analyze_shift_invariance(
        D11_53, D12_53, p53, n53, label="Paley p=53"
    )

    stab_53, mult_53 = analyze_multiplicative_orbits(
        D11_53, D12_53, p53, n53, label="Paley p=53"
    )

    # =========================================================
    # Test 3: Check the exact N values from p31 exhaustive enumeration
    # =========================================================
    print(f"\n{'='*70}")
    print("Exact N values from p=31 exhaustive enumeration:")
    results_file = os.path.join(os.path.dirname(__file__), 'exact_N_p31_results.jsonl')
    if os.path.exists(results_file):
        nonzero_orbits = []
        total_orbits = 0
        with open(results_file) as f:
            for line in f:
                r = json.loads(line.strip())
                total_orbits += 1
                if r['N'] > 0:
                    nonzero_orbits.append(r)

        print(f"Total orbits: {total_orbits}")
        print(f"Nonzero orbits: {len(nonzero_orbits)}")
        for r in nonzero_orbits:
            N = r['N']
            print(f"  Orbit {r['orbit']}: N={N}, N/(2p)={N/(2*31):.1f}, "
                  f"N mod 2p = {N % (2*31)}, N mod p = {N % 31}")
    else:
        print("  File not found, skipping.")

    # =========================================================
    # Summary
    # =========================================================
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")

    print(f"\n1. Additive shift invariance of D12:")
    print(f"   p=31 sol#0: {len(valid_shifts_31)}/{p31} shifts valid")
    print(f"   p=31 sol#2: {len(valid_shifts_31b)}/{p31} shifts valid")
    print(f"   Paley p=53: {len(valid_shifts_53)}/{p53} shifts valid")

    if len(valid_shifts_31) == p31:
        print(f"\n   CONCLUSION: Additive shifts DO preserve validity => factor of p explained")
    else:
        print(f"\n   CONCLUSION: Additive shifts do NOT all preserve validity.")
        print(f"   The factor of p does NOT come from additive shift invariance of D12.")

    print(f"\n2. Multiplicative orbit structure:")
    print(f"   p=31 D11 stabilizer size: {len(stab_31)}")
    n_mult_valid_31 = sum(1 for _, _, v in mult_31 if v)
    print(f"   p=31 valid mult images: {n_mult_valid_31}/{len(stab_31)}")
    print(f"   p=53 Paley D11 stabilizer size: {len(stab_53)}")
    n_mult_valid_53 = sum(1 for _, _, v in mult_53 if v)
    print(f"   p=53 valid mult images: {n_mult_valid_53}/{len(stab_53)}")

    print(f"\n3. True source of 2p factor:")
    print(f"   The 2p divisibility likely arises from the CIRCULANT GRAPH STRUCTURE:")
    print(f"   - Factor of p: the circulant automorphism group Z_p acts on")
    print(f"     vertex labels, mapping solution -> solution (same D-sets, different")
    print(f"     vertex numbering). BUT this does NOT change D12.")
    print(f"   - So if N(D11) counts DISTINCT D12 sets, the factor of p comes from")
    print(f"     some other symmetry of the count, not from D12 shifts.")
    print(f"   - Factor of 2: negation D12 -> -D12 (since D11 symmetric)")


if __name__ == '__main__':
    main()
