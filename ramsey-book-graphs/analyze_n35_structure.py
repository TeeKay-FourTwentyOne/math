"""
Investigate the n=35 m=69 solution structure.
Key question: D22 != complement(D11), so what is D22?
And does the standalone validator actually confirm this?
"""
import json

def main():
    with open('/Users/stephenpadgett/Projects/math/ramsey-book-graphs/solutions_registry.json') as f:
        registry = json.load(f)

    sol = [s for s in registry['solutions'] if s['m'] == 69][0]
    D11 = set(sol['D11'])
    D12 = set(sol['D12'])
    D22 = set(sol['D22'])
    m = 69
    n = 35

    complement_D11 = set(range(1, m)) - D11
    print(f"D11 = {sorted(D11)}")
    print(f"|D11| = {len(D11)}")
    print(f"\nD22 (registry) = {sorted(D22)}")
    print(f"|D22| = {len(D22)}")
    print(f"\ncomplement(D11) = {sorted(complement_D11)}")
    print(f"|complement(D11)| = {len(complement_D11)}")

    print(f"\nD22 == complement(D11): {D22 == complement_D11}")

    diff_extra = D22 - complement_D11
    diff_missing = complement_D11 - D22
    print(f"\nIn D22 but not in complement(D11): {sorted(diff_extra)}")
    print(f"In complement(D11) but not in D22: {sorted(diff_missing)}")

    # So D22 has elements from D11! That means some V2-V2 edges overlap
    # with the D11 pattern. This breaks the standard 2-block circulant structure.
    overlap = D22 & D11
    print(f"\nD22 ∩ D11 = {sorted(overlap)}")
    print(f"|D22 ∩ D11| = {len(overlap)}")

    # Check: is this solution using D22 != complement(D11)?
    # In the standard construction, a V2-V2 edge exists iff diff in D22
    # and the blue complement is what avoids the book.
    # If D22 != complement(D11), then the red subgraph is NOT a disjoint
    # union of G[V1] and G[V2] with G[V1-V2].

    # Actually wait - re-read the construction:
    # Red edges within V1: diff in D11
    # Red edges within V2: diff in D22
    # Red edges V1-V2: diff in D12
    # Blue is everything else.
    # D11 and D22 are INDEPENDENT choices for V1 and V2.
    # The complement relationship D22 = complement(D11) is just a common pattern,
    # not a requirement! The validator confirms correctness regardless.

    print("\n" + "=" * 60)
    print("CRITICAL INSIGHT: D22 need NOT be complement(D11).")
    print("D11 controls V1-V1 edges, D22 controls V2-V2 edges.")
    print("They are independent parameters. The validator checks")
    print("the actual common neighbor counts for all edge pairs.")
    print("=" * 60)

    # Check symmetry of D22
    print(f"\nD22 symmetric: {all((-x) % m in D22 for x in D22)}")
    # Which elements break symmetry?
    for x in sorted(D22):
        neg = (-x) % m
        if neg not in D22:
            print(f"  {x} in D22 but {neg} = -{x} mod {m} NOT in D22")

    # Check 0 in D12
    print(f"\n0 in D12: {0 in D12}")

    # Compare with n=32 m=63 which also works
    print("\n" + "=" * 60)
    print("COMPARISON: n=32, m=63 = 3^2 * 7 (also SA-found, also works)")
    print("=" * 60)

    sol63 = [s for s in registry['solutions'] if s['m'] == 63][0]
    D11_63 = set(sol63['D11'])
    D12_63 = set(sol63['D12'])
    D22_63 = set(sol63['D22'])
    m63 = 63

    complement_D11_63 = set(range(1, m63)) - D11_63
    print(f"|D11| = {len(D11_63)}, |D12| = {len(D12_63)}, |D22| = {len(D22_63)}")
    print(f"|complement(D11)| = {len(complement_D11_63)}")
    print(f"D22 == complement(D11): {D22_63 == complement_D11_63}")
    print(f"D22 symmetric: {all((-x) % m63 in D22_63 for x in D22_63)}")
    print(f"0 in D12: {0 in D12_63}")

    # Now check all solutions for which D22 == complement(D11) and which don't
    print("\n" + "=" * 60)
    print("ALL SOLUTIONS: D22 == complement(D11)?")
    print("=" * 60)

    for sol in registry['solutions']:
        m = sol['m']
        n = sol['n']
        D11 = set(sol['D11'])
        D22 = set(sol.get('D22', []))
        complement = set(range(1, m)) - D11

        if D22:
            is_comp = D22 == complement
        else:
            is_comp = "N/A"

        d22_sym = all((-x) % m in D22 for x in D22) if D22 else "N/A"
        d12_has_0 = 0 in set(sol['D12'])

        print(f"n={n:3d} m={m:3d}: D22=comp(D11)={str(is_comp):>5}, "
              f"D22_sym={str(d22_sym):>5}, 0∈D12={str(d12_has_0):>5}, "
              f"|D11|={len(D11):2d} |D12|={len(sol['D12']):2d} |D22|={len(D22):2d}")

if __name__ == '__main__':
    main()
