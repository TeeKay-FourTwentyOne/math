#!/usr/bin/env python3
"""
Verify additive shift invariance for composite m cases.

The proof that D12 -> D12+c preserves validity works for ANY m (not just prime).
The key identities are:
  B(d; D12+c) = B(d; D12)           [autocorrelation shift-invariant]
  X(d; D12+c) = X(d-c; D12)         [cross-term shifts cyclically]
  d in D12+c  <=>  d-c in D12       [membership shifts accordingly]

This script verifies for composite m cases from the solutions registry.
Also checks: does the orbit act FREELY? I.e., is D12+c = D12 only for c=0?
"""

import json, sys, os
sys.path.insert(0, os.path.dirname(__file__))
from ramsey_core import BlockCirculantGraph, verify_construction


def shift_set(S, c, m):
    return {(x + c) % m for x in S}

def negate_set(S, m):
    return {(-x) % m for x in S}


def test_solution(sol_data, m, n, label):
    D11 = set(sol_data['D11'] if isinstance(sol_data, dict) and 'D11' in sol_data else sol_data[0])
    D12 = set(sol_data['D12'] if isinstance(sol_data, dict) and 'D12' in sol_data else sol_data[1])
    D22 = set(sol_data['D22'] if isinstance(sol_data, dict) and 'D22' in sol_data else sol_data[2])

    # Verify original
    G = BlockCirculantGraph(n=n, D11=D11, D12=D12, D22=D22)
    result = verify_construction(G)
    if not result.valid:
        print(f"  {label}: ORIGINAL INVALID, skipping")
        return

    # Test all shifts
    valid_count = 0
    for c in range(m):
        D12_c = shift_set(D12, c, m)
        G_c = BlockCirculantGraph(n=n, D11=D11, D12=D12_c, D22=D22)
        r = verify_construction(G_c)
        if r.valid:
            valid_count += 1

    # Test negation
    D12_neg = negate_set(D12, m)
    G_neg = BlockCirculantGraph(n=n, D11=D11, D12=D12_neg, D22=D22)
    neg_valid = verify_construction(G_neg).valid

    # Check freeness: is D12+c = D12 for any c != 0?
    D12_frozen = frozenset(D12)
    stabilizer = [c for c in range(m) if frozenset(shift_set(D12, c, m)) == D12_frozen]

    # Check if negation equals some shift
    D12_neg_frozen = frozenset(D12_neg)
    neg_is_shift = any(frozenset(shift_set(D12, c, m)) == D12_neg_frozen for c in range(m))

    orbit_size = m // len(stabilizer)  # additive orbit
    if neg_is_shift:
        full_orbit = orbit_size
    else:
        full_orbit = 2 * orbit_size

    print(f"  {label}: m={m}, n={n}")
    print(f"    Valid shifts: {valid_count}/{m}  (ALL={valid_count==m})")
    print(f"    Negation valid: {neg_valid}")
    print(f"    Additive stabilizer: {stabilizer} (size {len(stabilizer)})")
    print(f"    -D12 is a shift: {neg_is_shift}")
    print(f"    Orbit size: {full_orbit}")
    print(f"    N(D11) divisible by: {full_orbit}")


def main():
    base = os.path.dirname(__file__)

    with open(os.path.join(base, 'solutions_registry.json')) as f:
        registry = json.load(f)

    print("=" * 70)
    print("ADDITIVE SHIFT INVARIANCE: COMPOSITE m CASES")
    print("=" * 70)

    for sol in registry['solutions']:
        n = sol['n']
        m = sol['m']
        test_solution(sol, m, n, f"n={n}")


if __name__ == '__main__':
    main()
