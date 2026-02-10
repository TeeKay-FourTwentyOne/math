"""
Computational verification of the V1V2 auto-satisfaction theorem.

Claim: When D22 = complement(D11) in {1,...,m-1}, D11 symmetric, and |D12| = n-1,
the V1V2 (cross-block) constraints are automatically satisfied:
  - Red V1V2 edges have <= n-2 red common neighbors
  - Blue V1V2 edges have <= n-1 blue common neighbors

This script checks the claim by brute-force on known constructions for
n = 6, 8, 10, 12, 14, 16, 18, 20, 22.

NO dependencies on ramsey_core.py - fully self-contained.
"""

from itertools import combinations

# Known constructions (from analyze_known.py and solution_n22.json)
KNOWN = {
    6: {
        "D11": [3, 5, 6, 8],
        "D12": [0, 1, 4, 6, 7],
        "D22": [1, 2, 4, 7, 9, 10],
    },
    8: {
        "D11": [3, 6, 7, 8, 9, 12],
        "D12": [0, 1, 4, 6, 8, 9, 13],
        "D22": [1, 2, 4, 5, 10, 11, 13, 14],
    },
    10: {
        "D11": [4, 5, 7, 9, 10, 12, 14, 15],
        "D12": [0, 1, 2, 6, 7, 10, 11, 13, 17],
        "D22": [1, 2, 3, 6, 8, 11, 13, 16, 17, 18],
    },
    12: {
        "D11": [5, 6, 7, 8, 9, 11, 12, 14, 15, 16, 17, 18],
        "D12": [0, 1, 2, 6, 10, 13, 14, 16, 18, 20, 21],
        "D22": [1, 2, 3, 4, 10, 13, 19, 20, 21, 22],
    },
    14: {
        "D11": [5, 7, 8, 9, 10, 11, 13, 14, 16, 17, 18, 19, 20, 22],
        "D12": [0, 1, 2, 7, 8, 10, 13, 14, 17, 18, 21, 23, 25],
        "D22": [1, 2, 3, 4, 6, 12, 15, 21, 23, 24, 25, 26],
    },
    16: {
        "D11": [6, 7, 8, 10, 11, 12, 14, 15, 16, 17, 19, 20, 21, 23, 24, 25],
        "D12": [0, 1, 2, 3, 8, 11, 12, 13, 15, 18, 20, 21, 24, 27, 29],
        "D22": [1, 2, 3, 4, 5, 9, 13, 18, 22, 26, 27, 28, 29, 30],
    },
    18: {
        "D11": [6, 8, 9, 10, 11, 13, 14, 15, 17, 18, 20, 21, 22, 24, 25, 26, 27, 29],
        "D12": [0, 1, 2, 3, 8, 9, 10, 12, 14, 15, 19, 20, 22, 24, 25, 28, 32],
        "D22": [1, 2, 3, 4, 5, 7, 12, 16, 19, 23, 28, 30, 31, 32, 33, 34],
    },
    20: {
        "D11": [7, 8, 9, 10, 12, 13, 14, 16, 18, 19, 20, 21, 23, 25, 26, 27, 29, 30, 31, 32],
        "D12": [0, 1, 2, 5, 6, 8, 12, 14, 15, 17, 18, 22, 23, 25, 27, 28, 33, 36, 37],
        "D22": [1, 2, 3, 4, 5, 6, 11, 15, 17, 22, 24, 28, 33, 34, 35, 36, 37, 38],
    },
    22: {
        "D11": [1, 2, 5, 10, 11, 13, 16, 17, 18, 19, 20, 23, 24, 25, 26, 27, 30, 32, 33, 38, 41, 42],
        "D12": [0, 2, 5, 6, 8, 11, 15, 16, 20, 24, 25, 27, 28, 31, 32, 34, 35, 36, 37, 39, 41],
        "D22": [3, 4, 6, 7, 8, 9, 12, 14, 15, 21, 22, 28, 29, 31, 34, 35, 36, 37, 39, 40],
    },
}


def verify_v1v2_theorem(n, D11, D12, D22):
    """
    For a given construction, check:
    1. Preconditions: D22=complement(D11), D11 symmetric, |D12|=n-1
    2. V1V2 constraints: brute-force check all cross-block pairs
    3. Report max red/blue common neighbors for V1V2 edges specifically
    Also report V1V1 and V2V2 for comparison.
    """
    m = 2 * n - 1
    N = 2 * m

    D11_set = set(D11)
    D12_set = set(D12)
    D22_set = set(D22)

    red_threshold = n - 2
    blue_threshold = n - 1

    # Check preconditions
    d11_sym = all((m - d) % m in D11_set for d in D11)
    full_set = set(range(1, m))
    d22_complement = D22_set == full_set - D11_set
    d12_size_ok = len(D12) == n - 1

    preconditions_met = d11_sym and d22_complement and d12_size_ok

    # Build adjacency matrix
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
                if u_in_V1:
                    diff = (v - m - u) % m
                else:
                    diff = (v - u + m) % m
                if diff in D12_set:
                    adj[u][v] = adj[v][u] = 1

    # Check V1V2 pairs only (u in V1, v in V2)
    max_red_cn_v1v2 = 0
    max_blue_cn_v1v2 = 0
    v1v2_red_violations = 0
    v1v2_blue_violations = 0
    v1v2_red_edges = 0
    v1v2_blue_edges = 0

    # Also track V1V1 and V2V2 for comparison
    max_red_cn_v1v1 = 0
    max_blue_cn_v1v1 = 0
    max_red_cn_v2v2 = 0
    max_blue_cn_v2v2 = 0
    v1v1_violations = 0
    v2v2_violations = 0

    for u in range(N):
        for v in range(u + 1, N):
            u_in_V1 = u < m
            v_in_V1 = v < m
            is_red = adj[u][v] == 1

            if is_red:
                cn = sum(1 for w in range(N) if w != u and w != v and adj[u][w] == 1 and adj[v][w] == 1)
            else:
                cn = sum(1 for w in range(N) if w != u and w != v and adj[u][w] == 0 and adj[v][w] == 0)

            if u_in_V1 and not v_in_V1:
                # V1V2 pair
                if is_red:
                    v1v2_red_edges += 1
                    max_red_cn_v1v2 = max(max_red_cn_v1v2, cn)
                    if cn > red_threshold:
                        v1v2_red_violations += 1
                else:
                    v1v2_blue_edges += 1
                    max_blue_cn_v1v2 = max(max_blue_cn_v1v2, cn)
                    if cn > blue_threshold:
                        v1v2_blue_violations += 1
            elif u_in_V1 and v_in_V1:
                # V1V1
                if is_red:
                    max_red_cn_v1v1 = max(max_red_cn_v1v1, cn)
                    if cn > red_threshold:
                        v1v1_violations += 1
                else:
                    max_blue_cn_v1v1 = max(max_blue_cn_v1v1, cn)
                    if cn > blue_threshold:
                        v1v1_violations += 1
            else:
                # V2V2
                if is_red:
                    max_red_cn_v2v2 = max(max_red_cn_v2v2, cn)
                    if cn > red_threshold:
                        v2v2_violations += 1
                else:
                    max_blue_cn_v2v2 = max(max_blue_cn_v2v2, cn)
                    if cn > blue_threshold:
                        v2v2_violations += 1

    total_v1v2_violations = v1v2_red_violations + v1v2_blue_violations
    total_other_violations = v1v1_violations + v2v2_violations

    return {
        "n": n,
        "m": m,
        "preconditions": {
            "D11_symmetric": d11_sym,
            "D22_is_complement": d22_complement,
            "|D12|_equals_n-1": d12_size_ok,
            "all_met": preconditions_met,
        },
        "v1v2": {
            "red_edges": v1v2_red_edges,
            "blue_edges": v1v2_blue_edges,
            "max_red_cn": max_red_cn_v1v2,
            "max_blue_cn": max_blue_cn_v1v2,
            "red_violations": v1v2_red_violations,
            "blue_violations": v1v2_blue_violations,
            "total_violations": total_v1v2_violations,
        },
        "v1v1": {
            "max_red_cn": max_red_cn_v1v1,
            "max_blue_cn": max_blue_cn_v1v1,
            "violations": v1v1_violations,
        },
        "v2v2": {
            "max_red_cn": max_red_cn_v2v2,
            "max_blue_cn": max_blue_cn_v2v2,
            "violations": v2v2_violations,
        },
        "red_threshold": red_threshold,
        "blue_threshold": blue_threshold,
    }


def main():
    print("=" * 70)
    print("COMPUTATIONAL VERIFICATION: V1V2 Auto-Satisfaction Theorem")
    print("=" * 70)
    print()
    print("Claim: When D22=complement(D11), D11 symmetric, |D12|=n-1,")
    print("       V1V2 constraints are automatically satisfied.")
    print()

    all_pass = True

    for n in sorted(KNOWN.keys()):
        data = KNOWN[n]

        # Handle set vs list (analyze_known.py uses sets)
        D11 = sorted(data["D11"]) if isinstance(data["D11"], set) else data["D11"]
        D12 = sorted(data["D12"]) if isinstance(data["D12"], set) else data["D12"]
        D22 = sorted(data["D22"]) if isinstance(data["D22"], set) else data["D22"]

        result = verify_v1v2_theorem(n, D11, D12, D22)

        pre = result["preconditions"]
        v12 = result["v1v2"]
        v11 = result["v1v1"]
        v22 = result["v2v2"]

        precond_str = "YES" if pre["all_met"] else "NO"
        v1v2_str = "PASS" if v12["total_violations"] == 0 else "FAIL"
        v1v1_str = "PASS" if v11["violations"] == 0 else "FAIL"
        v2v2_str = "PASS" if v22["violations"] == 0 else "FAIL"

        print(f"n={n:2d}  m={result['m']:2d}  "
              f"preconditions={precond_str:3s}  "
              f"V1V2={v1v2_str} (max_red={v12['max_red_cn']:2d}/{result['red_threshold']}, "
              f"max_blue={v12['max_blue_cn']:2d}/{result['blue_threshold']})  "
              f"V1V1={v1v1_str}  V2V2={v2v2_str}")

        if v12["total_violations"] > 0:
            print(f"       *** V1V2 VIOLATIONS: {v12['red_violations']} red, {v12['blue_violations']} blue")
            all_pass = False

        if not pre["all_met"]:
            details = []
            if not pre["D11_symmetric"]:
                details.append("D11 not symmetric")
            if not pre["D22_is_complement"]:
                details.append("D22 != complement(D11)")
            if not pre["|D12|_equals_n-1"]:
                details.append(f"|D12|={len(D12)} != {n-1}")
            print(f"       Precondition issues: {', '.join(details)}")

    print()
    print("=" * 70)
    if all_pass:
        print("THEOREM VERIFICATION: CONFIRMED")
        print("V1V2 constraints are satisfied for ALL tested constructions")
        print("where preconditions hold (n = 6, 8, 10, 12, 14, 16, 18, 20, 22).")
    else:
        print("THEOREM VERIFICATION: COUNTEREXAMPLE FOUND")
        print("At least one construction violates V1V2 constraints.")
    print("=" * 70)


if __name__ == "__main__":
    main()
