"""
Standalone brute-force validator for R(B_21, B_22) = 87.

Verifies the 2-block circulant construction on 86 vertices by:
1. Building the full 86x86 adjacency matrix from D11, D12, D22
2. Checking structural properties (symmetry, complement, cardinalities)
3. For every pair (u,v), counting common neighbors in the edge's color
4. Verifying no red edge has >20 red common neighbors (avoids B_21)
5. Verifying no blue edge has >21 blue common neighbors (avoids B_22)

NO dependencies on ramsey_core.py - fully self-contained.
"""

import sys
from itertools import combinations

# === THE SOLUTION ===
m = 43  # Z_m circulant, m = 2n-1 for n=22
N = 86  # total vertices = 2m

D11 = [1, 2, 5, 10, 11, 13, 16, 17, 18, 19, 20, 23, 24, 25, 26, 27, 30, 32, 33, 38, 41, 42]
D12 = [0, 2, 5, 6, 8, 11, 15, 16, 20, 24, 25, 27, 28, 31, 32, 34, 35, 36, 37, 39, 41]
D22 = [3, 4, 6, 7, 8, 9, 12, 14, 15, 21, 22, 28, 29, 31, 34, 35, 36, 37, 39, 40]

# Convert to sets for O(1) lookup
D11_set = set(D11)
D12_set = set(D12)
D22_set = set(D22)

# === THRESHOLDS ===
n = 22
red_threshold = n - 2   # 20: max red common neighbors per red edge
blue_threshold = n - 1   # 21: max blue common neighbors per blue edge

# === STEP 1: Verify structural properties ===
print("=" * 60)
print("STANDALONE VALIDATOR: R(B_21, B_22) = 87")
print("=" * 60)
print()
print("--- Structural Properties ---")

# Check cardinalities
print(f"|D11| = {len(D11)} (expected 22): {'OK' if len(D11) == 22 else 'FAIL'}")
print(f"|D12| = {len(D12)} (expected 21): {'OK' if len(D12) == 21 else 'FAIL'}")
print(f"|D22| = {len(D22)} (expected 20): {'OK' if len(D22) == 20 else 'FAIL'}")

# Check D11 is symmetric: d in D11 iff (m-d) in D11
d11_sym = all((m - d) % m in D11_set for d in D11)
print(f"D11 symmetric: {'OK' if d11_sym else 'FAIL'}")

# Check D22 is symmetric
d22_sym = all((m - d) % m in D22_set for d in D22)
print(f"D22 symmetric: {'OK' if d22_sym else 'FAIL'}")

# Check D22 = complement of D11 in {1,...,42}
full_set = set(range(1, m))
d22_is_complement = D22_set == full_set - D11_set
print(f"D22 = complement(D11): {'OK' if d22_is_complement else 'FAIL'}")

# Check 0 in D12
print(f"0 in D12: {'OK' if 0 in D12_set else 'FAIL'}")

# Degrees
d1 = len(D11) + len(D12)  # degree of V1 vertices
d2 = len(D22) + len(D12)  # degree of V2 vertices
print(f"d1 = |D11| + |D12| = {d1}")
print(f"d2 = |D22| + |D12| = {d2}")
print()

# === STEP 2: Build the 86x86 adjacency matrix ===
print("--- Building adjacency matrix ---")

# adj[u][v] = 1 if red edge, 0 if blue (no self-loops)
adj = [[0] * N for _ in range(N)]

# V1 = {0, ..., 42}, V2 = {43, ..., 85}
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
            # Both in V2: edge iff ((v-43) - (u-43)) mod m in D22
            diff = (v - u) % m  # same as ((v-43)-(u-43)) % m
            if diff in D22_set:
                adj[u][v] = 1
                adj[v][u] = 1
        else:
            # One in V1, one in V2
            if u_in_V1:
                # u in V1, v in V2: diff = (v - 43 - u) mod m
                diff = (v - m - u) % m
            else:
                # u in V2, v in V1: diff = (v - (u - 43)) mod m = (v - u + 43) mod m
                diff = (v - u + m) % m
            if diff in D12_set:
                adj[u][v] = 1
                adj[v][u] = 1

# Verify degrees
deg = [sum(adj[u]) for u in range(N)]
v1_degs = set(deg[:m])
v2_degs = set(deg[m:])
print(f"V1 degrees: {v1_degs} (expected {{{d1}}})")
print(f"V2 degrees: {v2_degs} (expected {{{d2}}})")
print()

# === STEP 3: Check all C(86,2) = 3655 pairs ===
print("--- Checking all pairs ---")

total_pairs = 0
red_edges = 0
blue_edges = 0
max_red_cn = 0   # max red common neighbors over red edges
max_blue_cn = 0  # max blue common neighbors over blue edges
violations = []

for u, v in combinations(range(N), 2):
    total_pairs += 1
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
            violations.append(('RED', u, v, cn))
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
            violations.append(('BLUE', u, v, cn))

print(f"Total pairs checked: {total_pairs}")
print(f"Red edges: {red_edges}")
print(f"Blue edges: {blue_edges}")
print(f"Max red common neighbors:  {max_red_cn} (threshold {red_threshold})")
print(f"Max blue common neighbors: {max_blue_cn} (threshold {blue_threshold})")
print(f"Violations: {len(violations)}")

if violations:
    print()
    print("VIOLATIONS FOUND:")
    for color, u, v, cn in violations[:20]:
        print(f"  {color} edge ({u},{v}): {cn} common neighbors")
    if len(violations) > 20:
        print(f"  ... and {len(violations) - 20} more")

# === STEP 4: Final verdict ===
print()
print("=" * 60)
structural_ok = (
    len(D11) == 22 and len(D12) == 21 and len(D22) == 20
    and d11_sym and d22_sym and d22_is_complement
    and 0 in D12_set
    and v1_degs == {d1} and v2_degs == {d2}
)

if structural_ok and len(violations) == 0:
    print("RESULT: PASS")
    print()
    print("The 2-block circulant graph on 86 vertices contains")
    print("no red B_21 and no blue B_22.")
    print("Therefore R(B_21, B_22) >= 87.")
    print("Combined with the upper bound R(B_21, B_22) <= 87,")
    print("this proves R(B_21, B_22) = 87.")
else:
    print("RESULT: FAIL")
    if not structural_ok:
        print("  Structural property check failed")
    if violations:
        print(f"  {len(violations)} constraint violation(s)")
    sys.exit(1)

print("=" * 60)
