"""
Standalone brute-force validator for R(B_{n-1}, B_n) = 4n - 1.

Verifies a 2-block circulant construction on N = 4n - 2 vertices by:
1. Building the full N x N adjacency matrix from D11, D12, D22
2. Checking structural properties (symmetry, complement, cardinalities)
3. For every pair (u, v), counting common neighbors in the edge's color
4. Verifying no red edge has > n-2 red common neighbors (avoids B_{n-1})
5. Verifying no blue edge has > n-1 blue common neighbors (avoids B_n)

NO dependencies on ramsey_core.py - fully self-contained.

Usage:
  python validate_construction.py --json solution_n22.json
  python validate_construction.py --n 24 --d11 "1,2,5,..." --d12 "0,1,3,..." --d22 "3,4,..."
"""

import sys
import json
import argparse
import time
from itertools import combinations


def parse_int_list(s):
    """Parse a comma-separated string of integers."""
    return [int(x.strip()) for x in s.split(",") if x.strip()]


def load_from_json(filepath):
    """Load construction parameters from a JSON solution file."""
    with open(filepath, "r") as f:
        data = json.load(f)
    n = data["n"]
    params = data["parameters"]
    D11 = params["D11"]
    D12 = params["D12"]
    D22 = params["D22"]
    return n, D11, D12, D22


def validate(n, D11, D12, D22):
    """
    Validate a 2-block circulant construction for R(B_{n-1}, B_n) = 4n - 1.

    Returns (passed: bool, report: str).
    """
    lines = []

    def log(msg=""):
        lines.append(msg)

    m = 2 * n - 1
    N = 2 * m  # = 4n - 2

    red_threshold = n - 2  # max red common neighbors per red edge
    blue_threshold = n - 1  # max blue common neighbors per blue edge

    D11_set = set(D11)
    D12_set = set(D12)
    D22_set = set(D22)

    log("=" * 60)
    log(f"STANDALONE VALIDATOR: R(B_{{{n-1}}}, B_{{{n}}}) = {4*n - 1}")
    log(f"  n = {n}, m = {m}, N = {N}")
    log("=" * 60)
    log()

    # === STEP 1: Verify structural properties ===
    log("--- Structural Properties ---")

    structural_ok = True

    # Cardinality constraints:
    # D11 subset of {1,...,m-1}, symmetric. |D11| must be even (symmetric pairs).
    # D22 = complement of D11 in {1,...,m-1}, so |D22| = m - 1 - |D11|
    # D12 subset of {0,...,m-1}.
    # |D12| = n - 1 is the standard pattern (required for V1V2 auto-satisfaction).
    #
    # For m ≡ 3 (mod 4): |D11| = n, |D12| = n-1, |D22| = n-2
    # For m ≡ 1 (mod 4): |D11| = n-1 (since n is odd, need even |D11|),
    #   |D12| = n-1, |D22| = n-1. This gives symmetric degrees d1 = d2.

    d12_card_ok = len(D12) == n - 1
    d11_even_ok = len(D11) % 2 == 0
    d22_expected = m - 1 - len(D11)
    d22_card_ok = len(D22) == d22_expected

    log(f"|D11| = {len(D11)} (even: {'OK' if d11_even_ok else 'FAIL'})")
    log(f"|D12| = {len(D12)} (expected {n - 1}): {'OK' if d12_card_ok else 'FAIL'}")
    log(f"|D22| = {len(D22)} (expected {d22_expected} = m-1-|D11|): {'OK' if d22_card_ok else 'FAIL'}")

    if not (d12_card_ok and d11_even_ok and d22_card_ok):
        structural_ok = False

    # Check D11 is symmetric: d in D11 implies (m - d) in D11
    d11_sym = all((m - d) % m in D11_set for d in D11)
    log(f"D11 symmetric: {'OK' if d11_sym else 'FAIL'}")
    if not d11_sym:
        structural_ok = False

    # Check D22 is symmetric
    d22_sym = all((m - d) % m in D22_set for d in D22)
    log(f"D22 symmetric: {'OK' if d22_sym else 'FAIL'}")
    if not d22_sym:
        structural_ok = False

    # Check D22 = complement of D11 in {1,...,m-1}
    full_set = set(range(1, m))
    d22_is_complement = D22_set == full_set - D11_set
    log(f"D22 = complement(D11): {'OK' if d22_is_complement else 'FAIL'}")
    if not d22_is_complement:
        structural_ok = False

    # Check 0 in D12
    zero_in_d12 = 0 in D12_set
    log(f"0 in D12: {'OK' if zero_in_d12 else 'FAIL'}")
    if not zero_in_d12:
        structural_ok = False

    # Check no element in D11 or D22 is 0
    no_zero_d11 = 0 not in D11_set
    no_zero_d22 = 0 not in D22_set
    log(f"0 not in D11: {'OK' if no_zero_d11 else 'FAIL'}")
    log(f"0 not in D22: {'OK' if no_zero_d22 else 'FAIL'}")
    if not (no_zero_d11 and no_zero_d22):
        structural_ok = False

    # Check all elements are in valid range
    d11_range_ok = all(1 <= d <= m - 1 for d in D11)
    d12_range_ok = all(0 <= d <= m - 1 for d in D12)
    d22_range_ok = all(1 <= d <= m - 1 for d in D22)
    log(f"D11 range [1, {m-1}]: {'OK' if d11_range_ok else 'FAIL'}")
    log(f"D12 range [0, {m-1}]: {'OK' if d12_range_ok else 'FAIL'}")
    log(f"D22 range [1, {m-1}]: {'OK' if d22_range_ok else 'FAIL'}")
    if not (d11_range_ok and d12_range_ok and d22_range_ok):
        structural_ok = False

    # Check no duplicates
    d11_no_dup = len(D11) == len(D11_set)
    d12_no_dup = len(D12) == len(D12_set)
    d22_no_dup = len(D22) == len(D22_set)
    log(f"D11 no duplicates: {'OK' if d11_no_dup else 'FAIL'}")
    log(f"D12 no duplicates: {'OK' if d12_no_dup else 'FAIL'}")
    log(f"D22 no duplicates: {'OK' if d22_no_dup else 'FAIL'}")
    if not (d11_no_dup and d12_no_dup and d22_no_dup):
        structural_ok = False

    # Degrees
    d1 = len(D11) + len(D12)  # degree of V1 vertices
    d2 = len(D22) + len(D12)  # degree of V2 vertices
    log(f"d1 = |D11| + |D12| = {len(D11)} + {len(D12)} = {d1}")
    log(f"d2 = |D22| + |D12| = {len(D22)} + {len(D12)} = {d2}")
    log()

    # === STEP 2: Build the N x N adjacency matrix ===
    log("--- Building adjacency matrix ---")
    t0 = time.time()

    # adj[u][v] = 1 if red edge, 0 if blue/no edge (no self-loops)
    adj = [[0] * N for _ in range(N)]

    # V1 = {0, ..., m-1}, V2 = {m, ..., 2m-1}
    for u in range(N):
        for v in range(u + 1, N):
            u_in_V1 = u < m
            v_in_V1 = v < m

            if u_in_V1 and v_in_V1:
                # Both in V1: edge iff (v - u) mod m in D11
                diff = (v - u) % m
                if diff in D11_set:
                    adj[u][v] = 1
                    adj[v][u] = 1
            elif (not u_in_V1) and (not v_in_V1):
                # Both in V2: edge iff ((v-m) - (u-m)) mod m in D22
                diff = (v - u) % m
                if diff in D22_set:
                    adj[u][v] = 1
                    adj[v][u] = 1
            else:
                # One in V1, one in V2
                if u_in_V1:
                    # u in V1, v in V2: diff = (v - m - u) mod m
                    diff = (v - m - u) % m
                else:
                    # u in V2, v in V1: diff = (v - (u - m)) mod m = (v - u + m) mod m
                    diff = (v - u + m) % m
                if diff in D12_set:
                    adj[u][v] = 1
                    adj[v][u] = 1

    t1 = time.time()
    log(f"Adjacency matrix built in {t1 - t0:.2f}s")

    # Verify degree regularity
    deg = [sum(adj[u]) for u in range(N)]
    v1_degs = set(deg[:m])
    v2_degs = set(deg[m:])
    v1_deg_ok = v1_degs == {d1}
    v2_deg_ok = v2_degs == {d2}
    log(f"V1 degrees: {v1_degs} (expected {{{d1}}}): {'OK' if v1_deg_ok else 'FAIL'}")
    log(f"V2 degrees: {v2_degs} (expected {{{d2}}}): {'OK' if v2_deg_ok else 'FAIL'}")
    if not (v1_deg_ok and v2_deg_ok):
        structural_ok = False
    log()

    # === STEP 3: Check all C(N, 2) pairs ===
    total_pairs = N * (N - 1) // 2
    log(f"--- Checking all C({N}, 2) = {total_pairs} pairs ---")
    t2 = time.time()

    red_edges = 0
    blue_edges = 0
    max_red_cn = 0   # max red common neighbors over red edges
    max_blue_cn = 0  # max blue common neighbors over blue edges
    violations = []

    for u, v in combinations(range(N), 2):
        is_red = adj[u][v] == 1

        if is_red:
            red_edges += 1
            # Count red common neighbors: w such that adj[u][w]==1 and adj[v][w]==1
            cn = 0
            for w in range(N):
                if w != u and w != v and adj[u][w] == 1 and adj[v][w] == 1:
                    cn += 1
            if cn > max_red_cn:
                max_red_cn = cn
            if cn > red_threshold:
                violations.append(("RED", u, v, cn))
        else:
            blue_edges += 1
            # Count blue common neighbors: w such that adj[u][w]==0 and adj[v][w]==0
            cn = 0
            for w in range(N):
                if w != u and w != v and adj[u][w] == 0 and adj[v][w] == 0:
                    cn += 1
            if cn > max_blue_cn:
                max_blue_cn = cn
            if cn > blue_threshold:
                violations.append(("BLUE", u, v, cn))

    t3 = time.time()
    log(f"Pair checking completed in {t3 - t2:.1f}s")
    log()
    log(f"Total pairs checked: {total_pairs}")
    log(f"Red edges: {red_edges}")
    log(f"Blue edges: {blue_edges}")
    log(f"Max red common neighbors:  {max_red_cn} (threshold {red_threshold})")
    log(f"Max blue common neighbors: {max_blue_cn} (threshold {blue_threshold})")
    log(f"Violations: {len(violations)}")

    if violations:
        log()
        log("VIOLATIONS FOUND:")
        for color, u, v, cn in violations[:20]:
            log(f"  {color} edge ({u},{v}): {cn} common neighbors")
        if len(violations) > 20:
            log(f"  ... and {len(violations) - 20} more")

    # === STEP 4: Final verdict ===
    log()
    log("=" * 60)

    passed = structural_ok and len(violations) == 0

    if passed:
        log("RESULT: PASS")
        log()
        log(f"The 2-block circulant graph on {N} vertices contains")
        log(f"no red B_{{{n-1}}} and no blue B_{{{n}}}.")
        log(f"Therefore R(B_{{{n-1}}}, B_{{{n}}}) >= {4*n - 1}.")
        log(f"Combined with the upper bound R(B_{{{n-1}}}, B_{{{n}}}) <= {4*n - 1},")
        log(f"this proves R(B_{{{n-1}}}, B_{{{n}}}) = {4*n - 1}.")
    else:
        log("RESULT: FAIL")
        if not structural_ok:
            log("  Structural property check failed")
        if violations:
            log(f"  {len(violations)} constraint violation(s)")

    log("=" * 60)

    total_time = time.time() - t0
    log(f"Total validation time: {total_time:.1f}s")

    report = "\n".join(lines)

    return passed, report, {
        "max_red_common": max_red_cn,
        "max_blue_common": max_blue_cn,
        "red_threshold": red_threshold,
        "blue_threshold": blue_threshold,
        "red_edges": red_edges,
        "blue_edges": blue_edges,
        "total_pairs": total_pairs,
        "violations": len(violations),
    }


def main():
    parser = argparse.ArgumentParser(
        description="Standalone brute-force validator for R(B_{n-1}, B_n) = 4n-1"
    )
    parser.add_argument("--json", type=str, help="Path to JSON solution file")
    parser.add_argument("--n", type=int, help="Book parameter n")
    parser.add_argument("--d11", type=str, help="D11 as comma-separated integers")
    parser.add_argument("--d12", type=str, help="D12 as comma-separated integers")
    parser.add_argument("--d22", type=str, help="D22 as comma-separated integers")

    args = parser.parse_args()

    if args.json:
        n, D11, D12, D22 = load_from_json(args.json)
    elif args.n and args.d11 and args.d12 and args.d22:
        n = args.n
        D11 = parse_int_list(args.d11)
        D12 = parse_int_list(args.d12)
        D22 = parse_int_list(args.d22)
    else:
        parser.error("Provide either --json or all of --n, --d11, --d12, --d22")
        return

    passed, report, stats = validate(n, D11, D12, D22)
    print(report)

    if not passed:
        sys.exit(1)


if __name__ == "__main__":
    main()
