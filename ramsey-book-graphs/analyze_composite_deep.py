"""
Deep analysis: investigate the n=35/m=69 blue slack issue,
and compare autocorrelation distributions across prime vs composite.
Also check the GF(49) autocorrelation issue.
"""

import json
from collections import defaultdict

def Delta(A, B, d, m):
    B_set = set(B)
    return sum(1 for a in A if (a - d) % m in B_set)

def Sigma(A, B, d, m):
    B_set = set(B)
    return sum(1 for a in A if (d - a) % m in B_set)

def is_symmetric(S, m):
    return all((-x) % m in S for x in S)

def full_validation(D11, D12, D22, m, n):
    """Full validation with detailed violation reporting."""
    D11_set = set(D11)
    D12_set = set(D12)
    D22_set = set(D22)
    D12T = {(-x) % m for x in D12}

    N = 2 * m
    d1 = len(D11) + len(D12)
    d2 = len(D22) + len(D12)

    red_threshold = n - 2
    blue_threshold = n - 1

    violations = []

    # V1V1
    for d in range(1, m):
        common_red = Delta(D11_set, D11_set, d, m) + Delta(D12_set, D12_set, d, m)
        if d in D11_set:
            if common_red > red_threshold:
                violations.append(('V1V1_red', d, common_red, red_threshold))
        else:
            common_blue = (N - 2) - d1 - d1 + common_red
            if common_blue > blue_threshold:
                violations.append(('V1V1_blue', d, common_blue, blue_threshold))

    # V2V2
    for d in range(1, m):
        common_red = Delta(D22_set, D22_set, d, m) + Delta(D12T, D12T, d, m)
        if d in D22_set:
            if common_red > red_threshold:
                violations.append(('V2V2_red', d, common_red, red_threshold))
        else:
            common_blue = (N - 2) - d2 - d2 + common_red
            if common_blue > blue_threshold:
                violations.append(('V2V2_blue', d, common_blue, blue_threshold))

    # V1V2
    for d in range(m):
        common_red = Sigma(D11_set, D12_set, d, m) + Delta(D12_set, D22_set, d, m)
        if d in D12_set:
            if common_red > red_threshold:
                violations.append(('V1V2_red', d, common_red, red_threshold))
        else:
            common_blue = (N - 2) - d1 - d2 + common_red
            if common_blue > blue_threshold:
                violations.append(('V1V2_blue', d, common_blue, blue_threshold))

    return violations

def main():
    with open('/Users/stephenpadgett/Projects/math/ramsey-book-graphs/solutions_registry.json') as f:
        registry = json.load(f)

    print("=" * 80)
    print("DEEP INVESTIGATION: VIOLATIONS AND GF(49)")
    print("=" * 80)

    for sol in registry['solutions']:
        m = sol['m']
        n = sol['n']
        D11 = sol['D11']
        D12 = sol['D12']
        D22 = sol.get('D22')
        if D22 is None:
            D22 = sorted(set(range(1, m)) - set(D11))

        violations = full_validation(D11, D12, D22, m, n)

        if violations:
            print(f"\nn={n}, m={m}: {len(violations)} VIOLATIONS found!")
            for vtype, d, val, threshold in violations:
                print(f"  {vtype} d={d}: {val} > {threshold}")
        else:
            print(f"n={n}, m={m}: OK (0 violations)")

    # Special investigation of GF(49)
    print("\n" + "=" * 80)
    print("GF(49) INVESTIGATION")
    print("=" * 80)

    sol49 = [s for s in registry['solutions'] if s['m'] == 49][0]
    D11 = sol49['D11']
    D12 = sol49['D12']
    D22 = sol49['D22']
    m = 49
    n = 25

    print(f"\nD11 = {D11}")
    print(f"D12 = {D12}")
    print(f"D22 = {D22}")
    print(f"|D11| = {len(D11)}, |D12| = {len(D12)}, |D22| = {len(D22)}")
    print(f"D11 == D12: {set(D11) == set(D12)}")

    # GF(49) uses different arithmetic (not Z_49)
    # The adjacency is defined by GF(49) arithmetic, not circulant over Z_49
    # So the circulant-based autocorrelation formulas don't apply directly
    print("\nNOTE: GF(49) uses finite field arithmetic, NOT Z_49 circulant.")
    print("The D11/D12 sets encode GF(49) elements as a+7b (0-indexed).")
    print("Circulant autocorrelation formulas are WRONG for GF(49).")
    print("The validator works because it builds the full adjacency matrix.")

    # For n=35 m=69
    print("\n" + "=" * 80)
    print("n=35, m=69 INVESTIGATION")
    print("=" * 80)

    sol69 = [s for s in registry['solutions'] if s['m'] == 69][0]
    D11 = sol69['D11']
    D12 = sol69['D12']
    D22 = sol69.get('D22')
    m = 69
    n = 35

    if D22 is None:
        D22 = sorted(set(range(1, m)) - set(D11))
        print(f"D22 computed as complement: {D22}")

    print(f"|D11| = {len(D11)}, |D12| = {len(D12)}, |D22| = {len(D22)}")
    print(f"D11 symmetric: {is_symmetric(set(D11), m)}")
    print(f"D22 symmetric: {is_symmetric(set(D22), m)}")

    # Check D12 size
    print(f"\n|D12| = {len(D12)}, expected n-1 = {n-1}")
    print(f"0 in D12: {0 in set(D12)}")

    # Check D22 = complement(D11)?
    expected_d22 = set(range(1, m)) - set(D11)
    print(f"D22 = complement(D11) in {{1,...,{m-1}}}: {set(D22) == expected_d22}")

    violations = full_validation(D11, D12, D22, m, n)
    if violations:
        print(f"\n{len(violations)} violations in circulant analysis:")
        for vtype, d, val, threshold in violations:
            print(f"  {vtype} d={d}: {val} > {threshold}")
    else:
        print("\n0 violations - construction is valid as circulant")

    # Check the registry verification data
    print(f"\nRegistry verification: {sol69.get('verification', {})}")

    # For n=35, let's check if D22 is explicitly in the registry or if
    # we need to compute it
    print(f"\nD22 in registry: {'D22' in sol69}")

    # Let's also check the standalone solution file
    with open('/Users/stephenpadgett/Projects/math/ramsey-book-graphs/solution_n35_sa.json') as f:
        sol35_file = json.load(f)
    print(f"\nSolution file data:")
    print(f"  D11_size: {sol35_file.get('D11_size')}")
    print(f"  D12_size: {sol35_file.get('D12_size')}")
    print(f"  'D22' key present: {'D22' in sol35_file}")
    if 'D22' in sol35_file:
        print(f"  D22 from file: {sol35_file['D22']}")
    else:
        print("  D22 not in file - computed as complement")

    # The issue: for n=35, m=69, m mod 4 = 1
    # Expected: |D11| = |D12| = |D22| = n-1 = 34
    # Actual: |D11| = 34, |D12| = 34, but we need to check |D22|
    print(f"\nm mod 4 = {m % 4}")
    print(f"If m≡1(4): expect |D11|=|D12|=|D22| = n-1 = {n-1}")
    print(f"If m≡3(4): expect |D11|=n={n}, |D12|=n-1={n-1}, |D22|=n-2={n-2}")

    # Actual sizes
    print(f"Actual: |D11|={len(D11)}, |D12|={len(D12)}, |D22|={len(D22)}")

    # For m≡1(4), d1 should equal d2 for balanced degrees
    d1 = len(D11) + len(D12)
    d2 = len(D22) + len(D12)
    print(f"d1 = {d1}, d2 = {d2}")

    # The issue: |D22| = m-1-|D11| = 68-34 = 34, not 33
    # Wait, let me recheck
    print(f"\nm-1-|D11| = {m-1-len(D11)}")

    # Hmm, D22 should be complement of D11 in {1,...,m-1}
    # If |D11|=34, then |D22| should be 68-34=34
    # But the registry says |D22|=33 in the notes. Let me check.

    # From registry: D22 is listed explicitly
    print(f"\nD22 from registry: {sol69.get('D22')}")
    if sol69.get('D22'):
        print(f"|D22 from registry| = {len(sol69['D22'])}")

if __name__ == '__main__':
    main()
