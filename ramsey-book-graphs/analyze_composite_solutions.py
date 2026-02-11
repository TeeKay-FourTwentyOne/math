"""
Analyze composite m solutions from the solutions registry.
For each composite m, compute:
1. Factorization and algebraic structure
2. Symmetry of D11, D12, D22
3. Full autocorrelation profiles A(d), B(d)
4. Coset structure analysis (CRT decomposition)
5. Comparison with quadratic residue structure
"""

import json
import math
from collections import defaultdict
from itertools import product

def factorize(n):
    """Return prime factorization as list of (prime, exponent) pairs."""
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
    if n < 2:
        return False
    if n < 4:
        return True
    if n % 2 == 0:
        return False
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True

def is_prime_power(n):
    factors = factorize(n)
    return len(factors) == 1

def is_symmetric(S, m):
    """Check if S is closed under negation mod m."""
    return all((-x) % m in S for x in S)

def Delta(A, B, d, m):
    """Count pairs (a,b) in A x B where a - b = d (mod m)."""
    B_set = set(B)
    return sum(1 for a in A if (a - d) % m in B_set)

def Sigma(A, B, d, m):
    """Count pairs (a,b) in A x B where a + b = d (mod m)."""
    B_set = set(B)
    return sum(1 for a in A if (d - a) % m in B_set)

def compute_autocorrelation_profile(D11, D12, D22, m, n):
    """
    Compute A(d) and B(d) for all differences d.

    A(d) = autocorrelation count for red edges at difference d
    B(d) = blue common neighbors for blue edges at difference d

    Returns dict with all autocorrelation data.
    """
    D11_set = set(D11)
    D12_set = set(D12)
    D22_set = set(D22)
    D12T = {(-x) % m for x in D12}

    N = 2 * m
    d1 = len(D11) + len(D12)
    d2 = len(D22) + len(D12)

    red_threshold = n - 2
    blue_threshold = n - 1

    results = {
        'V1V1': {},
        'V2V2': {},
        'V1V2': {},
    }

    # V1-V1 edges
    for d in range(1, m):
        common_red = Delta(D11_set, D11_set, d, m) + Delta(D12_set, D12_set, d, m)
        if d in D11_set:
            # Red edge
            results['V1V1'][d] = {'color': 'red', 'common': common_red, 'ok': common_red <= red_threshold}
        else:
            # Blue edge
            common_blue = (N - 2) - d1 - d1 + common_red
            results['V1V1'][d] = {'color': 'blue', 'common': common_blue, 'ok': common_blue <= blue_threshold}

    # V2-V2 edges
    for d in range(1, m):
        common_red = Delta(D22_set, D22_set, d, m) + Delta(D12T, D12T, d, m)
        if d in D22_set:
            results['V2V2'][d] = {'color': 'red', 'common': common_red, 'ok': common_red <= red_threshold}
        else:
            common_blue = (N - 2) - d2 - d2 + common_red
            results['V2V2'][d] = {'color': 'blue', 'common': common_blue, 'ok': common_blue <= blue_threshold}

    # V1-V2 edges
    for d in range(m):
        common_red = Sigma(D11_set, D12_set, d, m) + Delta(D12_set, D22_set, d, m)
        if d in D12_set:
            results['V1V2'][d] = {'color': 'red', 'common': common_red, 'ok': common_red <= red_threshold}
        else:
            common_blue = (N - 2) - d1 - d2 + common_red
            results['V1V2'][d] = {'color': 'blue', 'common': common_blue, 'ok': common_blue <= blue_threshold}

    return results

def quadratic_residues_mod(m):
    """Compute quadratic residues mod m."""
    qr = set()
    for x in range(1, m):
        qr.add((x * x) % m)
    qr.discard(0)
    return qr

def crt_decompose(S, m, factors):
    """
    Decompose set S in Z_m into CRT components.
    Given m = p1^a1 * p2^a2 * ..., each element x in Z_m maps to
    (x mod p1^a1, x mod p2^a2, ...).
    Returns the projection of S onto each factor.
    """
    moduli = [p**a for p, a in factors]
    projections = {}
    for pi, mi in enumerate(moduli):
        proj = set()
        for x in S:
            proj.add(x % mi)
        projections[mi] = proj
    return projections

def is_coset(S, m):
    """
    Check if S is a coset of some subgroup of Z_m.
    Returns (True, subgroup, offset) or (False, None, None).
    """
    S_set = set(S)
    # Check all divisors of m as potential subgroup orders
    divisors = []
    for d in range(1, m + 1):
        if m % d == 0:
            divisors.append(d)

    for d in divisors:
        if d == 0:
            continue
        step = m // d  # elements of the subgroup {0, step, 2*step, ...}
        subgroup = set(range(0, m, step))
        if len(subgroup) != d:
            continue
        # Check if S is a coset of this subgroup
        for offset in range(step):
            coset = {(s + offset) % m for s in subgroup}
            if coset == S_set:
                return True, subgroup, offset

    return False, None, None

def union_of_cosets_analysis(S, m, factors):
    """
    Check if S is a union of cosets of subgroups from the CRT decomposition.
    For m = p * q, check if S is a union of residue classes mod p or mod q.
    """
    S_set = set(S)
    results = {}

    for p, a in factors:
        modulus = p ** a
        # Group elements by their residue mod modulus
        classes = defaultdict(set)
        for x in S:
            classes[x % modulus].add(x)

        # Check which residue classes are fully contained
        full_classes = []
        partial_classes = []
        total_elements_in_m = m // modulus  # expected elements per class

        for r, elements in sorted(classes.items()):
            all_in_class = {x for x in range(r, m, modulus)}
            if all_in_class.issubset(S_set):
                full_classes.append(r)
            else:
                partial_classes.append((r, len(elements), len(all_in_class)))

        results[modulus] = {
            'full_classes': full_classes,
            'partial_classes': partial_classes,
            'residues_hit': sorted(classes.keys()),
            'num_residues_hit': len(classes),
            'total_residues': modulus,
        }

    return results

def autocorrelation_stats(profile):
    """Compute summary statistics of autocorrelation profile."""
    red_commons = []
    blue_commons = []

    for block in ['V1V1', 'V2V2', 'V1V2']:
        for d, info in profile[block].items():
            if info['color'] == 'red':
                red_commons.append(info['common'])
            else:
                blue_commons.append(info['common'])

    def stats(vals):
        if not vals:
            return {}
        mean = sum(vals) / len(vals)
        var = sum((v - mean)**2 for v in vals) / len(vals)
        return {
            'min': min(vals),
            'max': max(vals),
            'mean': round(mean, 3),
            'var': round(var, 3),
            'count': len(vals),
        }

    return {
        'red': stats(red_commons),
        'blue': stats(blue_commons),
    }

def main():
    with open('/Users/stephenpadgett/Projects/math/ramsey-book-graphs/solutions_registry.json') as f:
        registry = json.load(f)

    print("=" * 80)
    print("ANALYSIS OF COMPOSITE m SOLUTIONS")
    print("=" * 80)

    # Classify all solutions
    prime_solutions = []
    composite_solutions = []
    prime_power_solutions = []

    for sol in registry['solutions']:
        m = sol['m']
        n = sol['n']
        factors = factorize(m)
        if is_prime(m):
            prime_solutions.append(sol)
        elif is_prime_power(m):
            prime_power_solutions.append(sol)
        else:
            composite_solutions.append(sol)

    print(f"\nTotal solutions: {len(registry['solutions'])}")
    print(f"  Prime m: {len(prime_solutions)} -> m = {[s['m'] for s in prime_solutions]}")
    print(f"  Prime power m: {len(prime_power_solutions)} -> m = {[s['m'] for s in prime_power_solutions]}")
    print(f"  Composite m: {len(composite_solutions)} -> m = {[s['m'] for s in composite_solutions]}")

    # Also analyze prime solutions for comparison
    all_solutions = composite_solutions + prime_power_solutions

    print("\n" + "=" * 80)
    print("DETAILED ANALYSIS OF EACH COMPOSITE/PRIME-POWER SOLUTION")
    print("=" * 80)

    for sol in all_solutions:
        m = sol['m']
        n = sol['n']
        factors = factorize(m)
        D11 = sol['D11']
        D12 = sol['D12']
        D22 = sol.get('D22')

        # Compute D22 if not provided
        if D22 is None:
            full = set(range(1, m))
            D22 = sorted(full - set(D11))

        print(f"\n{'=' * 70}")
        print(f"n = {n}, m = {m} = {'*'.join(str(p) + ('^'+str(a) if a>1 else '') for p,a in factors)}")
        print(f"N = {2*m}, R(B_{{{n-1}}}, B_{{{n}}}) = {4*n - 1}")
        print(f"Solver: {sol.get('solver', 'unknown')}")
        print(f"{'=' * 70}")

        # 1. Set sizes
        print(f"\n--- Set Sizes ---")
        print(f"|D11| = {len(D11)}, |D12| = {len(D12)}, |D22| = {len(D22)}")
        print(f"Expected for m≡1(4): |D11|=|D12|=|D22|=n-1={n-1}")
        print(f"Expected for m≡3(4): |D11|=n={n}, |D12|=n-1={n-1}, |D22|=n-2={n-2}")
        print(f"m mod 4 = {m % 4}")

        d1 = len(D11) + len(D12)
        d2 = len(D22) + len(D12)
        print(f"d1 (V1 degree) = {d1}")
        print(f"d2 (V2 degree) = {d2}")

        # 2. Symmetry
        print(f"\n--- Symmetry ---")
        d11_sym = is_symmetric(set(D11), m)
        d12_sym = is_symmetric(set(D12), m)
        d22_sym = is_symmetric(set(D22), m)
        print(f"D11 symmetric: {d11_sym}")
        print(f"D12 symmetric: {d12_sym}")
        print(f"D22 symmetric: {d22_sym}")

        # 3. D11 = D12 check (Paley-like)
        d11_eq_d12 = set(D11) == set(D12)
        print(f"D11 == D12: {d11_eq_d12}")

        # 4. Quadratic residue analysis
        print(f"\n--- Quadratic Residue Analysis ---")
        qr = quadratic_residues_mod(m)
        qnr = set(range(1, m)) - qr
        d11_set = set(D11)

        d11_in_qr = d11_set & qr
        d11_in_qnr = d11_set & qnr
        print(f"QR(Z_{m}): {len(qr)} elements")
        print(f"QNR(Z_{m}): {len(qnr)} elements")
        print(f"|D11 ∩ QR| = {len(d11_in_qr)}, |D11 ∩ QNR| = {len(d11_in_qnr)}")
        if len(D11) > 0:
            print(f"QR fraction in D11: {len(d11_in_qr)/len(D11):.3f}")

        if d11_set == qr:
            print("*** D11 = QR(Z_m) ***")
        elif d11_set == qnr:
            print("*** D11 = QNR(Z_m) ***")

        # 5. CRT decomposition analysis
        if len(factors) > 1:
            print(f"\n--- CRT Decomposition (Z_{m} ≅ " + " x ".join(f"Z_{p**a}" for p,a in factors) + ") ---")

            for name, S in [('D11', D11), ('D12', D12), ('D22', D22)]:
                print(f"\n  {name}:")
                crt = crt_decompose(set(S), m, factors)
                for mi, info in crt.items():
                    print(f"    mod {mi}: residues hit = {sorted(info)}, count = {len(info)}/{mi}")

                # Union of cosets analysis
                coset_info = union_of_cosets_analysis(set(S), m, factors)
                for mi, info in coset_info.items():
                    if info['full_classes']:
                        print(f"    Full classes mod {mi}: {info['full_classes']}")
                    n_partial = len(info['partial_classes'])
                    if n_partial > 0 and n_partial <= 10:
                        for r, count, total in info['partial_classes']:
                            print(f"    Partial class {r} mod {mi}: {count}/{total} elements")

        # Also check if it's prime power
        if is_prime_power(m):
            p, a = factors[0]
            print(f"\n--- Prime Power Structure (GF({m}) = GF({p}^{a})) ---")
            print(f"  This m is a prime power. Paley-type construction over GF({m}) possible.")

        # 6. Autocorrelation profile
        print(f"\n--- Autocorrelation Profile ---")
        profile = compute_autocorrelation_profile(D11, D12, D22, m, n)
        stats = autocorrelation_stats(profile)

        print(f"  Red edges:  min={stats['red']['min']}, max={stats['red']['max']}, "
              f"mean={stats['red']['mean']}, var={stats['red']['var']}, count={stats['red']['count']}")
        print(f"  Blue edges: min={stats['blue']['min']}, max={stats['blue']['max']}, "
              f"mean={stats['blue']['mean']}, var={stats['blue']['var']}, count={stats['blue']['count']}")
        print(f"  Red threshold: {n-2}, Blue threshold: {n-1}")

        # Check if all autocorrelations are exactly at threshold (SRG property)
        red_all_equal = (stats['red']['min'] == stats['red']['max'])
        blue_all_equal = (stats['blue']['min'] == stats['blue']['max'])
        print(f"  Red constant autocorrelation (SRG-like): {red_all_equal}" +
              (f" (value={stats['red']['min']})" if red_all_equal else ""))
        print(f"  Blue constant autocorrelation (SRG-like): {blue_all_equal}" +
              (f" (value={stats['blue']['min']})" if blue_all_equal else ""))

        # Slack analysis: how close to thresholds?
        red_slack = n - 2 - stats['red']['max']
        blue_slack = n - 1 - stats['blue']['max']
        print(f"  Red slack (threshold - max): {red_slack}")
        print(f"  Blue slack (threshold - max): {blue_slack}")

        # Detailed autocorrelation values
        if m <= 70:
            print(f"\n  V1V1 autocorrelation values:")
            red_vals_v11 = []
            blue_vals_v11 = []
            for d in range(1, m):
                info = profile['V1V1'][d]
                if info['color'] == 'red':
                    red_vals_v11.append(info['common'])
                else:
                    blue_vals_v11.append(info['common'])
            print(f"    Red: {sorted(set(red_vals_v11))} (distribution: {dict(sorted(defaultdict(int, {v: red_vals_v11.count(v) for v in set(red_vals_v11)}).items()))})")
            print(f"    Blue: {sorted(set(blue_vals_v11))} (distribution: {dict(sorted(defaultdict(int, {v: blue_vals_v11.count(v) for v in set(blue_vals_v11)}).items()))})")

    # Summary comparison table
    print("\n" + "=" * 80)
    print("SUMMARY COMPARISON TABLE")
    print("=" * 80)
    print(f"{'n':>4} {'m':>4} {'type':>12} {'|D11|':>6} {'|D12|':>6} {'|D22|':>6} {'D11sym':>7} {'D12sym':>7} {'D11=D12':>8} {'maxR':>5} {'maxB':>5} {'Rslack':>7} {'Bslack':>7}")
    print("-" * 100)

    for sol in registry['solutions']:
        m = sol['m']
        n = sol['n']
        factors = factorize(m)
        D11 = sol['D11']
        D12 = sol['D12']
        D22 = sol.get('D22')
        if D22 is None:
            D22 = sorted(set(range(1, m)) - set(D11))

        if is_prime(m):
            mtype = f"prime({m%4})"
        elif is_prime_power(m):
            p, a = factors[0]
            mtype = f"GF({p}^{a})"
        else:
            mtype = '*'.join(str(p) + ('^'+str(a) if a>1 else '') for p,a in factors)

        d11_sym = is_symmetric(set(D11), m)
        d12_sym = is_symmetric(set(D12), m)
        d11_eq_d12 = set(D11) == set(D12)

        profile = compute_autocorrelation_profile(D11, D12, D22, m, n)
        stats = autocorrelation_stats(profile)

        red_slack = n - 2 - stats['red']['max']
        blue_slack = n - 1 - stats['blue']['max']

        print(f"{n:>4} {m:>4} {mtype:>12} {len(D11):>6} {len(D12):>6} {len(D22):>6} "
              f"{'Y' if d11_sym else 'N':>7} {'Y' if d12_sym else 'N':>7} {'Y' if d11_eq_d12 else 'N':>8} "
              f"{stats['red']['max']:>5} {stats['blue']['max']:>5} {red_slack:>7} {blue_slack:>7}")

    # Key question: do composite solutions share features with prime solutions?
    print("\n" + "=" * 80)
    print("KEY OBSERVATIONS")
    print("=" * 80)

    print("""
1. SYMMETRY: Do composite D11 sets maintain symmetry?
2. SIZE PATTERN: Do |D11|, |D12|, |D22| follow the same formulas?
3. TIGHTNESS: Are composite solutions tight (slack = 0)?
4. D11 = D12: Only happens for Paley constructions (m ≡ 1 mod 4, prime or prime power).
5. AUTOCORRELATION VARIANCE: Higher variance for composite vs prime?
""")

if __name__ == '__main__':
    main()
