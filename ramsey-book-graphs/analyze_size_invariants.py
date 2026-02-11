"""
Analyze the true size invariants across all solutions.
"""
import json

def main():
    with open('/Users/stephenpadgett/Projects/math/ramsey-book-graphs/solutions_registry.json') as f:
        registry = json.load(f)

    print("TRUE SIZE INVARIANTS")
    print("=" * 120)
    print(f"{'n':>3} {'m':>3} {'m%4':>4} {'|D11|':>5} {'|D12|':>5} {'|D22|':>5} "
          f"{'|D11|+|D22|':>12} {'m-1':>5} {'comp?':>6} "
          f"{'d1':>4} {'d2':>4} {'N-2':>5} {'red_thr':>8} {'blue_thr':>9} "
          f"{'d1+d2':>6} {'N-2-d1':>7} {'N-2-d2':>7}")
    print("-" * 120)

    for sol in registry['solutions']:
        m = sol['m']
        n = sol['n']
        N = 2 * m
        D11 = sol['D11']
        D12 = sol['D12']
        D22 = sol.get('D22', [])

        is_comp = set(D22) == set(range(1, m)) - set(D11)

        d1 = len(D11) + len(D12)  # V1 degree
        d2 = len(D22) + len(D12)  # V2 degree

        red_thr = n - 2
        blue_thr = n - 1

        print(f"{n:>3} {m:>3} {m%4:>4} {len(D11):>5} {len(D12):>5} {len(D22):>5} "
              f"{len(D11)+len(D22):>12} {m-1:>5} {'Y' if is_comp else 'N':>6} "
              f"{d1:>4} {d2:>4} {N-2:>5} {red_thr:>8} {blue_thr:>9} "
              f"{d1+d2:>6} {N-2-d1:>7} {N-2-d2:>7}")

    print()
    print("OBSERVATIONS:")
    print("  d1 = |D11| + |D12|  (V1 vertex red degree)")
    print("  d2 = |D22| + |D12|  (V2 vertex red degree)")
    print("  Blue degree of V1 vertex: N-1-d1")
    print("  Blue degree of V2 vertex: N-1-d2")
    print()
    print("  For a balanced construction: d1 = d2 = N/2 = m")
    print("  This happens when |D11| = |D22|, i.e., |D11| = (m-1)/2")

    # Check which solutions have d1 = d2
    print()
    for sol in registry['solutions']:
        m = sol['m']
        n = sol['n']
        d1 = len(sol['D11']) + len(sol['D12'])
        d2 = len(sol.get('D22', [])) + len(sol['D12'])
        balanced = (d1 == d2)
        # Also check if d1, d2 are close to m
        print(f"  m={m:>3}: d1={d1:>3}, d2={d2:>3}, N-1={2*m-1:>3}, "
              f"balanced={'Y' if balanced else 'N'}, "
              f"d1/N={d1/(2*m):.3f}, d2/N={d2/(2*m):.3f}")

if __name__ == '__main__':
    main()
