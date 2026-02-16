#!/usr/bin/env python3
"""
Verify SA solutions for composite m = 2n-1 values.

For each solution found in the log files, this script:
1. Parses D11 and D12 from the log
2. Derives D22 = {1,...,m-1} \ D11
3. Checks set sizes, symmetry, and range validity
4. Verifies ALL book graph constraints:
   - V1V1: For d in {1,...,m-1}, red common = A(d)+B(d), blue common via inclusion-exclusion
   - V2V2: For d in {1,...,m-1}, red common = C(d)+B'(d), blue common via inclusion-exclusion
   - V1V2: For d in {0,...,m-1}, red/blue common via Sigma+Delta
5. Optionally runs brute-force O(N^2) adjacency matrix verification

Uses the same constraint model as ramsey_core.py / validate_construction.py.
"""

import re
import sys
import time


# =====================================================================
# Core combinatorial functions (self-contained, no imports from ramsey_core)
# =====================================================================

def Delta(A, B, d, m):
    """Count a in A such that (a - d) mod m is in B."""
    return sum(1 for a in A if (a - d) % m in B)


def Sigma(A, B, d, m):
    """Count a in A such that (d - a) mod m is in B."""
    return sum(1 for a in A if (d - a) % m in B)


def is_symmetric(S, m):
    """Check if set S is closed under negation mod m."""
    return all((-x) % m in S for x in S)


# =====================================================================
# Log file parsing
# =====================================================================

def parse_log_file(filepath):
    """
    Parse a SA solver log file and extract solution(s).

    Returns a list of dicts with keys: n, m, D11, D12
    """
    with open(filepath, 'r') as f:
        text = f.read()

    solutions = []

    # Find all SOLUTION FOUND blocks
    # Pattern: lines after "*** SOLUTION FOUND! ***"
    blocks = text.split("*** SOLUTION FOUND! ***")

    for block in blocks[1:]:  # skip everything before the first solution
        lines = block.strip().split('\n')

        n_val = m_val = None
        d11_list = d12_list = None

        for line in lines:
            line = line.strip()

            # Parse n and m
            nm_match = re.match(r'n=(\d+),\s*m=(\d+)', line)
            if nm_match:
                n_val = int(nm_match.group(1))
                m_val = int(nm_match.group(2))

            # Parse D11
            d11_match = re.match(r'D11\s*=\s*\[([^\]]+)\]', line)
            if d11_match:
                d11_list = [int(x.strip()) for x in d11_match.group(1).split(',')]

            # Parse D12
            d12_match = re.match(r'D12\s*=\s*\[([^\]]+)\]', line)
            if d12_match:
                d12_list = [int(x.strip()) for x in d12_match.group(1).split(',')]

        if n_val is not None and d11_list is not None and d12_list is not None:
            solutions.append({
                'n': n_val,
                'm': m_val,
                'D11': d11_list,
                'D12': d12_list,
            })

    return solutions


# =====================================================================
# Verification using circulant constraint model
# =====================================================================

def verify_solution_circulant(n, m, D11_list, D12_list):
    """
    Verify using the circulant constraint model (same as ramsey_core.py).

    Returns (passed, report_lines, stats).
    """
    lines = []

    def log(msg=""):
        lines.append(msg)

    N = 2 * m  # Total vertices = 4n - 2
    red_threshold = n - 2    # Red edge: common_red <= n-2
    blue_threshold = n - 1   # Blue edge: common_blue <= n-1

    D11 = set(D11_list)
    D12 = set(D12_list)
    D22 = set(range(1, m)) - D11  # Complement in {1,...,m-1}

    log(f"  n={n}, m={m}, N={N}")
    log(f"  m = 2n-1 = {2*n-1} (composite: {not is_prime(m)})")
    log(f"  red_threshold={red_threshold}, blue_threshold={blue_threshold}")
    log()

    # --- Structural checks ---
    structural_ok = True

    log("  Structural checks:")

    # Size checks
    log(f"    |D11| = {len(D11)}  (expected ~(m-1)/2 = {(m-1)//2})")
    log(f"    |D12| = {len(D12)}  (expected n-1 = {n-1})")
    log(f"    |D22| = {len(D22)}  (= m-1-|D11| = {m - 1 - len(D11)})")

    if len(D12) != n - 1:
        log(f"    WARNING: |D12| = {len(D12)} != n-1 = {n-1}")

    # D11 should be a subset of {1,...,m-1}
    d11_range = all(1 <= d <= m - 1 for d in D11)
    log(f"    D11 in {{1,...,{m-1}}}: {'OK' if d11_range else 'FAIL'}")
    if not d11_range:
        structural_ok = False

    # D12 should be a subset of {0,...,m-1}
    d12_range = all(0 <= d <= m - 1 for d in D12)
    log(f"    D12 in {{0,...,{m-1}}}: {'OK' if d12_range else 'FAIL'}")
    if not d12_range:
        structural_ok = False

    # D11 symmetric?
    d11_sym = is_symmetric(D11, m)
    log(f"    D11 symmetric: {'YES' if d11_sym else 'NO'}")

    # D22 symmetric?
    d22_sym = is_symmetric(D22, m)
    log(f"    D22 symmetric: {'YES' if d22_sym else 'NO'}")

    # 0 in D12?
    log(f"    0 in D12: {'YES' if 0 in D12 else 'NO'}")

    # No duplicates
    d11_nodup = len(D11_list) == len(D11)
    d12_nodup = len(D12_list) == len(D12)
    log(f"    D11 no duplicates: {'OK' if d11_nodup else 'FAIL'}")
    log(f"    D12 no duplicates: {'OK' if d12_nodup else 'FAIL'}")
    if not (d11_nodup and d12_nodup):
        structural_ok = False

    # Degree computation
    d1 = len(D11) + len(D12)  # degree of V1 vertices
    d2 = len(D22) + len(D12)  # degree of V2 vertices
    log(f"    d1 = |D11|+|D12| = {len(D11)}+{len(D12)} = {d1}")
    log(f"    d2 = |D22|+|D12| = {len(D22)}+{len(D12)} = {d2}")

    # For the standard case: d1 = m-1+n-1 = m+n-2, d2 = m-1-|D11|+|D12|
    # The key requirement: d1 + d2 = N - 2 would give perfect balance
    log(f"    d1+d2 = {d1+d2}, N-2 = {N-2}")
    log()

    # --- Compute D12^T = {-x mod m : x in D12} ---
    D12T = {(-x) % m for x in D12}

    # --- V1V1 constraints ---
    log("  V1V1 constraints (d = 1,...,m-1):")
    v11_violations = []
    max_v11_red = 0
    max_v11_blue = 0

    for d in range(1, m):
        A_d = Delta(D11, D11, d, m)
        B_d = Delta(D12, D12, d, m)
        cv11 = A_d + B_d  # red common neighbors for V1V1 pair at difference d

        if d in D11:
            # Red edge: need cv11 <= red_threshold
            max_v11_red = max(max_v11_red, cv11)
            if cv11 > red_threshold:
                v11_violations.append(('V1V1_red', d, cv11, red_threshold))
        else:
            # Blue edge: blue_common = (N-2) - 2*d1 + cv11
            blue_common = (N - 2) - 2 * d1 + cv11
            max_v11_blue = max(max_v11_blue, blue_common)
            if blue_common > blue_threshold:
                v11_violations.append(('V1V1_blue', d, blue_common, blue_threshold))

    log(f"    max red common:  {max_v11_red} (threshold {red_threshold})")
    log(f"    max blue common: {max_v11_blue} (threshold {blue_threshold})")
    log(f"    violations: {len(v11_violations)}")
    if v11_violations:
        for vtype, d, val, thresh in v11_violations[:5]:
            log(f"      {vtype} d={d}: {val} > {thresh}")
        if len(v11_violations) > 5:
            log(f"      ... and {len(v11_violations) - 5} more")

    # --- V2V2 constraints ---
    log("  V2V2 constraints (d = 1,...,m-1):")
    v22_violations = []
    max_v22_red = 0
    max_v22_blue = 0

    for d in range(1, m):
        C_d = Delta(D22, D22, d, m)
        Bp_d = Delta(D12T, D12T, d, m)
        cv22 = C_d + Bp_d  # red common neighbors for V2V2 pair at difference d

        if d in D22:
            # Red edge in V2V2: need cv22 <= red_threshold
            max_v22_red = max(max_v22_red, cv22)
            if cv22 > red_threshold:
                v22_violations.append(('V2V2_red', d, cv22, red_threshold))
        else:
            # d in D11 means d NOT in D22, so V2V2 at difference d is blue
            blue_common = (N - 2) - 2 * d2 + cv22
            max_v22_blue = max(max_v22_blue, blue_common)
            if blue_common > blue_threshold:
                v22_violations.append(('V2V2_blue', d, blue_common, blue_threshold))

    log(f"    max red common:  {max_v22_red} (threshold {red_threshold})")
    log(f"    max blue common: {max_v22_blue} (threshold {blue_threshold})")
    log(f"    violations: {len(v22_violations)}")
    if v22_violations:
        for vtype, d, val, thresh in v22_violations[:5]:
            log(f"      {vtype} d={d}: {val} > {thresh}")
        if len(v22_violations) > 5:
            log(f"      ... and {len(v22_violations) - 5} more")

    # --- V1V2 constraints ---
    log("  V1V2 constraints (d = 0,...,m-1):")
    v12_violations = []
    max_v12_red = 0
    max_v12_blue = 0

    for d in range(m):
        sigma_val = Sigma(D11, D12, d, m)
        delta_val = Delta(D12, D22, d, m)
        common = sigma_val + delta_val

        if d in D12:
            # Red cross-edge: need common <= red_threshold
            max_v12_red = max(max_v12_red, common)
            if common > red_threshold:
                v12_violations.append(('V1V2_red', d, common, red_threshold))
        else:
            # Blue cross-edge: blue_common = (N-2) - d1 - d2 + common
            blue_common = (N - 2) - d1 - d2 + common
            max_v12_blue = max(max_v12_blue, blue_common)
            if blue_common > blue_threshold:
                v12_violations.append(('V1V2_blue', d, blue_common, blue_threshold))

    log(f"    max red common:  {max_v12_red} (threshold {red_threshold})")
    log(f"    max blue common: {max_v12_blue} (threshold {blue_threshold})")
    log(f"    violations: {len(v12_violations)}")
    if v12_violations:
        for vtype, d, val, thresh in v12_violations[:5]:
            log(f"      {vtype} d={d}: {val} > {thresh}")
        if len(v12_violations) > 5:
            log(f"      ... and {len(v12_violations) - 5} more")

    # --- Overall ---
    all_violations = v11_violations + v22_violations + v12_violations
    passed = structural_ok and len(all_violations) == 0

    log()
    max_red_overall = max(max_v11_red, max_v22_red, max_v12_red)
    max_blue_overall = max(max_v11_blue, max_v22_blue, max_v12_blue)
    log(f"  Overall: max_red_common={max_red_overall} (thresh {red_threshold}), "
        f"max_blue_common={max_blue_overall} (thresh {blue_threshold})")
    log(f"  Total violations: {len(all_violations)}")

    stats = {
        'max_red_common': max_red_overall,
        'max_blue_common': max_blue_overall,
        'red_threshold': red_threshold,
        'blue_threshold': blue_threshold,
        'num_violations': len(all_violations),
    }

    return passed, lines, stats


# =====================================================================
# Brute-force O(N^3) verification via adjacency matrix
# =====================================================================

def verify_brute_force(n, m, D11_list, D12_list):
    """
    Build the full N x N adjacency matrix and check every pair.

    This is the gold standard: no algebraic tricks, just raw counting.
    Returns (passed, report_lines, stats).
    """
    lines = []

    def log(msg=""):
        lines.append(msg)

    N = 2 * m
    red_threshold = n - 2
    blue_threshold = n - 1

    D11_set = set(D11_list)
    D12_set = set(D12_list)
    D22_set = set(range(1, m)) - D11_set

    t0 = time.time()

    # Build adjacency matrix
    # V1 = {0,...,m-1}, V2 = {m,...,2m-1}
    adj = [[0] * N for _ in range(N)]

    for u in range(N):
        for v in range(u + 1, N):
            u_in_V1 = u < m
            v_in_V1 = v < m

            if u_in_V1 and v_in_V1:
                diff = (v - u) % m
                if diff in D11_set:
                    adj[u][v] = adj[v][u] = 1
            elif not u_in_V1 and not v_in_V1:
                diff = (v - u) % m
                if diff in D22_set:
                    adj[u][v] = adj[v][u] = 1
            else:
                # Cross-block
                if u_in_V1:
                    diff = (v - m - u) % m
                else:
                    diff = (v - u + m) % m
                if diff in D12_set:
                    adj[u][v] = adj[v][u] = 1

    t1 = time.time()
    log(f"  Adjacency matrix built in {t1 - t0:.2f}s ({N}x{N})")

    # Verify degrees
    degs_v1 = set(sum(adj[u]) for u in range(m))
    degs_v2 = set(sum(adj[u]) for u in range(m, N))
    expected_d1 = len(D11_set) + len(D12_set)
    expected_d2 = len(D22_set) + len(D12_set)
    log(f"  V1 degrees: {degs_v1} (expected {{{expected_d1}}})")
    log(f"  V2 degrees: {degs_v2} (expected {{{expected_d2}}})")

    # Check all pairs
    max_red_cn = 0
    max_blue_cn = 0
    violations = []

    for u in range(N):
        for v in range(u + 1, N):
            if adj[u][v] == 1:
                # Red edge: count red common neighbors
                cn = 0
                for w in range(N):
                    if w != u and w != v and adj[u][w] == 1 and adj[v][w] == 1:
                        cn += 1
                if cn > max_red_cn:
                    max_red_cn = cn
                if cn > red_threshold:
                    violations.append(('RED', u, v, cn))
            else:
                # Blue edge: count blue common neighbors
                cn = 0
                for w in range(N):
                    if w != u and w != v and adj[u][w] == 0 and adj[v][w] == 0:
                        cn += 1
                if cn > max_blue_cn:
                    max_blue_cn = cn
                if cn > blue_threshold:
                    violations.append(('BLUE', u, v, cn))

    t2 = time.time()
    log(f"  Pair checking in {t2 - t1:.1f}s")
    log(f"  max_red_cn={max_red_cn} (thresh {red_threshold}), "
        f"max_blue_cn={max_blue_cn} (thresh {blue_threshold})")
    log(f"  violations: {len(violations)}")

    passed = len(violations) == 0

    stats = {
        'max_red_common': max_red_cn,
        'max_blue_common': max_blue_cn,
        'red_threshold': red_threshold,
        'blue_threshold': blue_threshold,
        'num_violations': len(violations),
    }

    return passed, lines, stats


def is_prime(n):
    """Simple primality test."""
    if n < 2:
        return False
    if n < 4:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True


def factorize(n):
    """Return prime factorization as a list of (prime, exponent) pairs."""
    factors = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            exp = 0
            while n % d == 0:
                n //= d
                exp += 1
            factors.append((d, exp))
        d += 1
    if n > 1:
        factors.append((n, 1))
    return factors


# =====================================================================
# Main
# =====================================================================

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Verify SA solutions for composite m = 2n-1"
    )
    parser.add_argument(
        '--brute-force', action='store_true',
        help="Also run O(N^3) brute-force verification (slow for large N)"
    )
    parser.add_argument(
        '--files', nargs='+', default=None,
        help="Log files to verify (default: standard composite solution files)"
    )
    args = parser.parse_args()

    # Default log files for composite m values
    if args.files is None:
        log_files = [
            '/tmp/sa_fast_n23.log',
            '/tmp/sa_fast_n26.log',
            '/tmp/sa_fast_n28.log',
            '/tmp/sa_fast_n29.log',
            '/tmp/sa_fast_n33.log',
        ]
    else:
        log_files = args.files

    print("=" * 70)
    print("COMPOSITE m SOLUTION VERIFICATION")
    print("Verifying R(B_{n-1}, B_n) = 4n-1 for composite m = 2n-1")
    print("=" * 70)
    print()

    all_passed = True
    summary = []

    for filepath in log_files:
        print(f"--- {filepath} ---")
        try:
            solutions = parse_log_file(filepath)
        except FileNotFoundError:
            print(f"  FILE NOT FOUND, skipping\n")
            all_passed = False
            summary.append((filepath, None, None, "FILE NOT FOUND"))
            continue

        if not solutions:
            print(f"  No solutions found in log file\n")
            summary.append((filepath, None, None, "NO SOLUTION"))
            continue

        for i, sol in enumerate(solutions):
            n = sol['n']
            m = sol['m']
            D11 = sol['D11']
            D12 = sol['D12']

            factors = factorize(m)
            factor_str = " * ".join(f"{p}^{e}" if e > 1 else str(p) for p, e in factors)

            print(f"  Solution {i+1}: n={n}, m={m} = {factor_str}")
            print(f"  Proving: R(B_{n-1}, B_{n}) = {4*n-1}")
            print()

            # Circulant model verification
            print("  [Circulant model verification]")
            passed, report, stats = verify_solution_circulant(n, m, D11, D12)

            for line in report:
                print(line)

            if passed:
                print(f"  >>> CIRCULANT VERIFICATION: PASS <<<")
            else:
                print(f"  >>> CIRCULANT VERIFICATION: FAIL <<<")
                all_passed = False

            # Brute force verification (optional)
            bf_passed = None
            if args.brute_force:
                print()
                print("  [Brute-force verification]")
                bf_passed, bf_report, bf_stats = verify_brute_force(n, m, D11, D12)
                for line in bf_report:
                    print(line)

                if bf_passed:
                    print(f"  >>> BRUTE-FORCE VERIFICATION: PASS <<<")
                else:
                    print(f"  >>> BRUTE-FORCE VERIFICATION: FAIL <<<")
                    all_passed = False

                # Cross-check stats
                if stats['max_red_common'] != bf_stats['max_red_common']:
                    print(f"  WARNING: max_red_common mismatch: "
                          f"circulant={stats['max_red_common']} vs brute={bf_stats['max_red_common']}")
                if stats['max_blue_common'] != bf_stats['max_blue_common']:
                    print(f"  WARNING: max_blue_common mismatch: "
                          f"circulant={stats['max_blue_common']} vs brute={bf_stats['max_blue_common']}")

            summary.append((filepath, n, m, "PASS" if passed else "FAIL"))
            print()

    # Final summary
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print()
    print(f"  {'File':<30} {'n':>4} {'m':>4}  {'m composite':>12}  {'Result':>8}")
    print(f"  {'-'*30} {'-'*4} {'-'*4}  {'-'*12}  {'-'*8}")
    for filepath, n, m, result in summary:
        fname = filepath.split('/')[-1]
        if n is not None:
            comp = "YES" if not is_prime(m) else "no (prime)"
            print(f"  {fname:<30} {n:>4} {m:>4}  {comp:>12}  {result:>8}")
        else:
            print(f"  {fname:<30} {'':>4} {'':>4}  {'':>12}  {result:>8}")

    print()
    if all_passed:
        print("ALL SOLUTIONS VERIFIED SUCCESSFULLY.")
        print()
        print("Each solution provides a 2-block circulant graph on N=4n-2 vertices")
        print("that avoids red B_{n-1} and blue B_n, proving R(B_{n-1}, B_n) >= 4n-1.")
        print("Combined with the known upper bound (Rousseau & Sheehan 1978),")
        print("this establishes R(B_{n-1}, B_n) = 4n-1 for these composite m values.")
    else:
        print("SOME VERIFICATIONS FAILED. See details above.")
        sys.exit(1)


if __name__ == "__main__":
    main()
