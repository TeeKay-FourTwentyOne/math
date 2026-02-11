"""
Analyze the size pattern more carefully:
For which composite m do the sizes follow the prime pattern?
"""
import json

def factorize(n):
    factors = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            exp = 0
            while n % d == 0:
                exp += 1
                n //= d
            factors.append((d, exp))
        d += 1
    if n > 1:
        factors.append((n, 1))
    return factors

def is_prime(n):
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0: return False
    d = 3
    while d * d <= n:
        if n % d == 0: return False
        d += 2
    return True

def main():
    with open('/Users/stephenpadgett/Projects/math/ramsey-book-graphs/solutions_registry.json') as f:
        registry = json.load(f)

    print("SIZE PATTERN ANALYSIS")
    print("=" * 90)
    print(f"{'n':>3} {'m':>3} {'m%4':>4} {'type':>10} {'|D11|':>5} {'|D12|':>5} {'|D22|':>5} "
          f"{'expect':>20} {'matches':>8} {'0∈D12':>6} {'D11=D12':>8}")
    print("-" * 90)

    for sol in registry['solutions']:
        m = sol['m']
        n = sol['n']
        factors = factorize(m)

        D11 = sol['D11']
        D12 = sol['D12']
        D22 = sol.get('D22', [])

        if is_prime(m):
            mtype = 'prime'
        elif len(factors) == 1:
            p, a = factors[0]
            mtype = f'GF({p}^{a})'
        else:
            mtype = '*'.join(str(p) + ('^'+str(a) if a > 1 else '') for p, a in factors)

        # Expected sizes depend on m mod 4
        if m % 4 == 1:
            exp_d11 = n - 1
            exp_d12 = n - 1
            exp_d22 = n - 1
            expect_str = f"|all|={n-1}"
        else:  # m % 4 == 3
            exp_d11 = n
            exp_d12 = n - 1
            exp_d22 = n - 2
            expect_str = f"{n},{n-1},{n-2}"

        matches = (len(D11) == exp_d11 and len(D12) == exp_d12 and
                   (not D22 or len(D22) == exp_d22))

        # Alternative: does it match m≡3 pattern regardless of actual m%4?
        alt_match = (len(D11) == n and len(D12) == n-1 and (not D22 or len(D22) == n-2))
        # Or m≡1 pattern?
        alt1_match = (len(D11) == n-1 and len(D12) == n-1 and (not D22 or len(D22) == n-1))

        d11_eq_d12 = set(D11) == set(D12)
        has_0 = 0 in set(D12)

        note = ""
        if not matches:
            if alt_match:
                note = " (matches m≡3 pattern!)"
            elif alt1_match:
                note = " (matches m≡1 pattern!)"
            else:
                note = f" (ANOMALOUS: {len(D11)},{len(D12)},{len(D22)})"

        print(f"{n:>3} {m:>3} {m%4:>4} {mtype:>10} {len(D11):>5} {len(D12):>5} {len(D22):>5} "
              f"{expect_str:>20} {'YES' if matches else 'NO'+note:>8} "
              f"{'Y' if has_0 else 'N':>6} {'Y' if d11_eq_d12 else 'N':>8}")

    print()
    print("KEY FINDINGS:")
    print("  - m ≡ 1 (mod 4) composites: m=45(3^2*5), m=57(3*19), m=65(5*13), m=69(3*23)")
    print("  - m ≡ 3 (mod 4) composites: m=51(3*17), m=55(5*11), m=63(3^2*7)")
    print()

    # Check: for composites with m≡1(4), do we always get D11=D12?
    # For primes with m≡1(4), the Paley construction gives D11=D12=QR(m)
    print("\nFor m ≡ 1 (mod 4):")
    for sol in registry['solutions']:
        m = sol['m']
        if m % 4 != 1:
            continue
        d11_eq_d12 = set(sol['D11']) == set(sol['D12'])
        print(f"  m={m}: D11==D12={d11_eq_d12}, solver={sol.get('solver','?')}")

    print("\nFor m ≡ 3 (mod 4):")
    for sol in registry['solutions']:
        m = sol['m']
        if m % 4 != 3:
            continue
        d11_eq_d12 = set(sol['D11']) == set(sol['D12'])
        print(f"  m={m}: D11==D12={d11_eq_d12}, solver={sol.get('solver','?')}")

if __name__ == '__main__':
    main()
