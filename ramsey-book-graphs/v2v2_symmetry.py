"""
Investigate the relationship between V1V1 and V2V2 constraints.

Since D22 = complement(D11) and D12^T = {-x mod m : x in D12}:
- V1V1 at difference d: Delta(D11, D11, d) + Delta(D12, D12, d)
- V2V2 at difference d: Delta(D22, D22, d) + Delta(D12^T, D12^T, d)

Question: Is there a simple algebraic relation between these?
"""

import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ramsey_core import Delta, Sigma

KNOWN = {
    6: {"D11": {3, 5, 6, 8}, "D12": {0, 1, 4, 6, 7}, "D22": {1, 2, 4, 7, 9, 10}},
    8: {"D11": {3, 6, 7, 8, 9, 12}, "D12": {0, 1, 4, 6, 8, 9, 13}, "D22": {1, 2, 4, 5, 10, 11, 13, 14}},
    10: {"D11": {4, 5, 7, 9, 10, 12, 14, 15}, "D12": {0, 1, 2, 6, 7, 10, 11, 13, 17}, "D22": {1, 2, 3, 6, 8, 11, 13, 16, 17, 18}},
    12: {"D11": {5, 6, 7, 8, 9, 11, 12, 14, 15, 16, 17, 18}, "D12": {0, 1, 2, 6, 10, 13, 14, 16, 18, 20, 21}, "D22": {1, 2, 3, 4, 10, 13, 19, 20, 21, 22}},
    14: {"D11": {5, 7, 8, 9, 10, 11, 13, 14, 16, 17, 18, 19, 20, 22}, "D12": {0, 1, 2, 7, 8, 10, 13, 14, 17, 18, 21, 23, 25}, "D22": {1, 2, 3, 4, 6, 12, 15, 21, 23, 24, 25, 26}},
    16: {"D11": {6, 7, 8, 10, 11, 12, 14, 15, 16, 17, 19, 20, 21, 23, 24, 25}, "D12": {0, 1, 2, 3, 8, 11, 12, 13, 15, 18, 20, 21, 24, 27, 29}, "D22": {1, 2, 3, 4, 5, 9, 13, 18, 22, 26, 27, 28, 29, 30}},
    18: {"D11": {6, 8, 9, 10, 11, 13, 14, 15, 17, 18, 20, 21, 22, 24, 25, 26, 27, 29}, "D12": {0, 1, 2, 3, 8, 9, 10, 12, 14, 15, 19, 20, 22, 24, 25, 28, 32}, "D22": {1, 2, 3, 4, 5, 7, 12, 16, 19, 23, 28, 30, 31, 32, 33, 34}},
    20: {"D11": {7, 8, 9, 10, 12, 13, 14, 16, 18, 19, 20, 21, 23, 25, 26, 27, 29, 30, 31, 32}, "D12": {0, 1, 2, 5, 6, 8, 12, 14, 15, 17, 18, 22, 23, 25, 27, 28, 33, 36, 37}, "D22": {1, 2, 3, 4, 5, 6, 11, 15, 17, 22, 24, 28, 33, 34, 35, 36, 37, 38}},
}

n22_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "solution_n22.json")
if os.path.exists(n22_path):
    with open(n22_path) as f:
        data = json.load(f)
    KNOWN[22] = {"D11": set(data["parameters"]["D11"]), "D12": set(data["parameters"]["D12"]), "D22": set(data["parameters"]["D22"])}


print("V1V1 vs V2V2 RELATIONSHIP ANALYSIS")
print("=" * 80)

for n in sorted(KNOWN.keys()):
    m = 2 * n - 1
    D11 = KNOWN[n]["D11"]
    D12 = KNOWN[n]["D12"]
    D22 = KNOWN[n]["D22"]
    D12T = {(-x) % m for x in D12}

    print(f"\nn={n}, m={m}, |D11|={len(D11)}, |D12|={len(D12)}, |D22|={len(D22)}")

    # For each d, compute V1V1 and V2V2 lambda_red
    for d in range(1, m):
        v11 = Delta(D11, D11, d, m) + Delta(D12, D12, d, m)
        v22 = Delta(D22, D22, d, m) + Delta(D12T, D12T, d, m)

        # What's the relationship?
        # Delta(D22, D22, d) + Delta(D12^T, D12^T, d)
        # Since D22 = {1,..,m-1}\D11:
        # Delta(D22, D22, d) = |{a in D22 : (a-d) in D22}|
        #   = |{a not in D11, a!=0 : (a-d) not in D11, (a-d)!=0}|
        #
        # Total: |{a in {1,..,m-1} : (a-d) in {1,..,m-1}}| = m-2 (all except a=d give (a-d) in {1,..,m-1}\{0},
        #   but a=0 is excluded... hmm let me think more carefully)
        #
        # Actually {1,...,m-1} = D11 ∪ D22 (disjoint).
        # Delta over full set {1,...,m-1}:
        # |{a in {1,...,m-1} : (a-d) mod m in {1,...,m-1}}|
        # = |{a in {1,...,m-1} : a != d}| = m-2 (since d in {1,...,m-1})
        #
        # This equals:
        # Delta(D11, D11, d) + Delta(D11, D22, d) + Delta(D22, D11, d) + Delta(D22, D22, d)
        # = m - 2
        #
        # Similarly: Delta(D11, D22, d) = |{a in D11: (a-d) in D22}| = |D11| - Delta(D11, D11, d) - [d in D11 and 0 in... wait]
        # Actually (a-d) mod m could be 0 when a=d. But a is in D11 and 0 is not in D11 or D22.
        # Delta(D11, D11, d) + Delta(D11, D22, d) + |{a in D11: (a-d)=0}|
        #   = |{a in D11: (a-d) in D11}| + |{a in D11: (a-d) in D22}| + [d in D11]
        #   = |D11|   (since every a in D11 has (a-d) in D11 ∪ D22 ∪ {0})
        # So Delta(D11, D22, d) = |D11| - Delta(D11, D11, d) - [d in D11]
        pass

    # Let me just compute the algebraic relationship directly
    # V1V1(d) = Delta(D11,D11,d) + Delta(D12,D12,d) = A(d) + B(d)
    # V2V2(d) = Delta(D22,D22,d) + Delta(D12T,D12T,d) = C(d) + E(d)
    #
    # Note: Delta(D12T, D12T, d) = Delta(D12, D12, -d) (by substitution)
    # If D12 were symmetric, then B(d) = E(d). But D12 is NOT symmetric.
    #
    # Also: A(d) + C(d) + cross terms = m-2
    # Specifically: A(d) + Delta(D11,D22,d) + Delta(D22,D11,d) + C(d) = m-2
    # And Delta(D11,D22,d) = |D11| - A(d) - [d in D11]
    # And Delta(D22,D11,d) = |D22| - C(d) - [d in D22]  (by symmetric argument)
    # Wait, need to be careful. D22 doesn't include 0 either.
    # Delta(D22, D11, d) = |{a in D22 : (a-d) in D11}|
    # Delta(D22, D22, d) = |{a in D22 : (a-d) in D22}|
    # These + |{a in D22: (a-d) = 0}| = |D22|
    # So Delta(D22, D11, d) + C(d) + [d in D22] = |D22|
    # Delta(D22, D11, d) = |D22| - C(d) - [d in D22]
    #
    # Similarly: Delta(D11, D22, d) = |D11| - A(d) - [d in D11]
    #
    # Sum: A + (|D11|-A-[d in D11]) + (|D22|-C-[d in D22]) + C = m-2
    # |D11| + |D22| - [d in D11] - [d in D22] = m-2
    # Since d in {1,...,m-1} and D11 ∪ D22 = {1,...,m-1}: [d in D11] + [d in D22] = 1
    # So |D11| + |D22| - 1 = m - 2
    # |D11| + |D22| = m - 1 ✓ (since they partition {1,...,m-1})

    # So the constraint on C(d) = Delta(D22,D22,d) given A(d) = Delta(D11,D11,d) is:
    # From the partition identity, we get:
    # C(d) = (m-2) - A(d) - Delta(D11,D22,d) - Delta(D22,D11,d)
    # = (m-2) - A(d) - (|D11| - A(d) - [d in D11]) - (|D22| - C(d) - [d in D22])
    # This is circular. Let me do it differently.
    #
    # Direct identity:
    # For d in {1,...,m-1}, define 1_S(x) for S ⊂ {1,...,m-1}.
    # Delta(S,S,d) = sum_{x=1}^{m-1} 1_S(x) * 1_S(x-d)
    # where x-d is mod m, and we only count when x-d != 0 (but since S ⊂ {1,...,m-1}, 1_S(0)=0 anyway)
    #
    # Since D22 = {1,...,m-1}\D11, 1_{D22}(x) = 1 - 1_{D11}(x) for x in {1,...,m-1} and 1_{D22}(0) = 0.
    # So Delta(D22,D22,d) = sum_{x=1}^{m-1} (1-1_{D11}(x))(1-1_{D11}(x-d))  [for (x-d) mod m != 0]
    # = sum_{x, x!=d} [1 - 1_{D11}(x) - 1_{D11}(x-d) + 1_{D11}(x)*1_{D11}(x-d)]
    # = (m-2) - sum 1_{D11}(x) [x=1..m-1, x!=d] - sum 1_{D11}(x-d) [x=1..m-1, x!=d] + A(d)
    #
    # sum 1_{D11}(x) for x=1..m-1, x!=d = |D11| - [d in D11]
    # sum 1_{D11}(x-d) for x=1..m-1, x!=d = sum 1_{D11}(y) for y in {1-d,...,m-1-d}\{0} mod m
    #   = sum 1_{D11}(y) for y != 0 (as y ranges over all of Z_m \ {0}) = |D11|
    #   Wait: x ranges over {1,...,m-1}\{d}. So y = x-d ranges over {1-d,...,m-1-d}\{0} mod m.
    #   This is {0,1,...,m-1}\{0, (-d) mod m}... no.
    #   x in {1,...,m-1}\{d}, y = x-d mod m.
    #   When x = d, y = 0 (excluded). For other x in {1,...,m-1}, y ranges over {1,...,m-1}\{0} with one value missing.
    #   Actually x goes through all {1,...,m-1} except d, so y = x-d goes through all {1-d,...,-1,1,...,m-1-d} mod m,
    #   i.e., all of Z_m except y=0 (from x=d) AND we need to handle the wraparound.
    #   Since we exclude x=d (giving y=0), y ranges over Z_m\{0} MINUS nothing else = {1,...,m-1}.
    #   Wait, x ranges over {1,...,m-1}\{d}. There are m-2 values of x. y = x-d mod m gives m-2 distinct values.
    #   These are Z_m \ {0, something}? No: y=x-d, x in {1,...,m-1}\{d}. The missing y values are: y from x=0 (not included), y from x=d (excluded).
    #   x=0 -> y=-d, x=d -> y=0. So y ranges over Z_m \ {0, -d mod m}.
    #   So sum = |D11| - [(-d mod m) in D11] = |D11| - [d in D11] (since D11 is symmetric).
    #
    # Therefore:
    # C(d) = (m-2) - (|D11| - [d in D11]) - (|D11| - [d in D11]) + A(d)
    # = (m-2) - 2|D11| + 2[d in D11] + A(d)

    print(f"\n  Verifying: C(d) = (m-2) - 2|D11| + 2*[d in D11] + A(d)")
    constant = (m - 2) - 2 * len(D11)

    all_ok = True
    for d in range(1, m):
        A_d = Delta(D11, D11, d, m)
        C_d = Delta(D22, D22, d, m)
        predicted_C = constant + 2 * (1 if d in D11 else 0) + A_d
        if C_d != predicted_C:
            print(f"  FAILED at d={d}: C(d)={C_d}, predicted={predicted_C}")
            all_ok = False

    if all_ok:
        print(f"  VERIFIED for all d! C(d) = {constant} + 2*[d in D11] + A(d)")

    # Now V2V2(d) = C(d) + E(d) where E(d) = Delta(D12T, D12T, d)
    # V1V1(d) = A(d) + B(d) where B(d) = Delta(D12, D12, d)
    #
    # C(d) = A(d) + constant + 2[d in D11]
    # So V2V2(d) = A(d) + constant + 2[d in D11] + E(d)
    #            = V1V1(d) - B(d) + constant + 2[d in D11] + E(d)
    #            = V1V1(d) + (E(d) - B(d)) + constant + 2[d in D11]
    #
    # E(d) - B(d) = Delta(D12T, D12T, d) - Delta(D12, D12, d)
    # Note: Delta(D12T, D12T, d) = |{a in D12T : (a-d) in D12T}|
    #   = |{-a' for a' in D12 : (-a'-d) in D12T}| = |{a' in D12 : (a'+d) in D12}|
    #   Hmm wait. D12T = {-x : x in D12}. So a in D12T means a = -a' for some a' in D12.
    #   (a-d) in D12T means a-d = -b' for some b' in D12, i.e., b' = d-a = d+a'.
    #   So Delta(D12T, D12T, d) = |{a' in D12 : (d + a') mod m in D12}| = Sigma(D12, D12, d, m)
    # No wait: Delta(A,B,d) = |{a in A : (a-d) in B}|.
    # Delta(D12T, D12T, d) with A = B = D12T:
    #   = |{a in D12T : (a-d) in D12T}|
    #   Substituting a = -a', a' in D12:
    #   = |{a' in D12 : (-a' - d) in D12T}|
    #   (-a' - d) in D12T means -(- a' - d) = a' + d in D12
    #   = |{a' in D12 : (a' + d) mod m in D12}|
    #   This is Sigma(D12, D12, d, m)? No.
    #   Sigma(D12, D12, d) = |{a in D12 : (d - a) in D12}|
    #   What I have is |{a' in D12 : (a' + d) in D12}| = |{a' in D12: (d - (-a')) in D12}|
    #   = Sigma({-a' : a' in D12}, D12, d) = Sigma(D12T, D12, d)?
    #   Hmm, more directly: |{a' in D12 : (a'+d) in D12}| = Delta(D12, D12, -d, m)
    #   Because Delta(D12, D12, -d) = |{a in D12: (a - (-d)) in D12}| = |{a in D12: (a+d) in D12}|.
    #   YES! So Delta(D12T, D12T, d) = Delta(D12, D12, -d, m) = Delta(D12, D12, m-d, m).

    print(f"\n  Verifying: Delta(D12T, D12T, d) = Delta(D12, D12, m-d)")
    all_ok2 = True
    for d in range(1, m):
        E_d = Delta(D12T, D12T, d, m)
        B_neg_d = Delta(D12, D12, (m - d) % m, m)
        if E_d != B_neg_d:
            print(f"  FAILED at d={d}: E(d)={E_d}, B(m-d)={B_neg_d}")
            all_ok2 = False
    if all_ok2:
        print(f"  VERIFIED: Delta(D12T, D12T, d) = Delta(D12, D12, m-d) for all d")

    # So V2V2(d) = (m-2 - 2|D11| + 2[d in D11]) + A(d) + B(m-d)
    #            = V1V1(m-d) + (m-2 - 2|D11| + 2[d in D11]) + A(d) - A(m-d)
    # Hmm, this doesn't simplify as neatly. But we know:
    # V2V2(d) = A(d) + (m-2-2|D11|) + 2[d in D11] + B(m-d)

    # Since D11 is symmetric, d in D11 iff (m-d) in D11. And also A(d) = A(m-d) by symmetry of D11.
    # Wait, IS A(d) = A(m-d)? Delta(D11, D11, d) counts |{a in D11: (a-d) in D11}|.
    # Delta(D11, D11, m-d) = |{a in D11: (a-(m-d)) in D11}| = |{a in D11: (a+d) in D11}|
    # Since D11 symmetric: (a+d) in D11 iff -(a+d) in D11 iff (-a-d) in D11.
    # Substituting a' = -a (which is also in D11 since D11 symmetric):
    # = |{a' in D11: (a'-d) in D11}| = Delta(D11, D11, d). ✓
    # So YES, A(d) = A(m-d).

    # Also [d in D11] = [(m-d) in D11] by symmetry.
    # So V2V2(d) = A(d) + (m-2-2|D11|) + 2[d in D11] + B(m-d)
    # And V1V1(d) = A(d) + B(d)
    # V2V2(d) = V1V1(d) + (m-2-2|D11|) + 2[d in D11] + B(m-d) - B(d)

    print(f"\n  V2V2(d) = A(d) + B(m-d) + {(m-2-2*len(D11))} + 2*[d in D11]")
    print(f"  V1V1(d) = A(d) + B(d)")
    print(f"  V2V2(d) - V1V1(d) = B(m-d) - B(d) + {(m-2-2*len(D11))} + 2*[d in D11]")

    # Check V2V2 with the formula
    print(f"\n  Verifying combined formula:")
    all_ok3 = True
    for d in range(1, m):
        A_d = Delta(D11, D11, d, m)
        B_d = Delta(D12, D12, d, m)
        B_neg_d = Delta(D12, D12, (m - d) % m, m)
        v22_computed = A_d + B_neg_d + (m - 2 - 2 * len(D11)) + 2 * (1 if d in D11 else 0)
        v22_actual = Delta(D22, D22, d, m) + Delta(D12T, D12T, d, m)
        if v22_computed != v22_actual:
            print(f"  FAILED at d={d}: computed={v22_computed}, actual={v22_actual}")
            all_ok3 = False
    if all_ok3:
        print(f"  VERIFIED!")

    # Key insight: The V1V1 and V2V2 constraints are NOT independent.
    # V2V2(d) = A(d) + B(m-d) + const + 2*[d in D11]
    # For d in D11 (red in V1V1) <=> d in D11 <=> (m-d) in D11 <=> (m-d) not in D22
    # So d in D11 means d is RED in V1V1 and BLUE in V2V2.
    # And d not in D11 means d is BLUE in V1V1 and RED in V2V2.
    # The coloring is FLIPPED between V1V1 and V2V2!

    print(f"\n  KEY: d in D11 (red V1V1) implies d in D22 = False (blue V2V2)")
    print(f"  So V1V1 red edges = V2V2 blue edges and vice versa.")

print("\n" + "=" * 80)
print("SUMMARY OF V2V2 STRUCTURE")
print("=" * 80)
print("""
THEOREM: Given D22 = complement(D11) and both symmetric, with D12T = {-x : x in D12}:

  V2V2(d) = Delta(D11,D11,d) + Delta(D12,D12,m-d) + (m - 2 - 2|D11|) + 2*[d in D11]

Key properties:
1. V1V1 RED edges (d in D11) are V2V2 BLUE edges (d not in D22)
2. V1V1 BLUE edges (d not in D11) are V2V2 RED edges (d in D22)
3. The D12 contribution "flips" from B(d) in V1V1 to B(m-d) in V2V2

This means:
- V2V2 red constraint on d (d in D22, d not in D11):
  A(d) + B(m-d) + (m-2-2|D11|) <= n-2
  A(d) is already constrained by V1V1 blue on this same d

- If D12 is nearly symmetric, B(d) ≈ B(m-d) and the constraints are nearly symmetric too.

The construction problem reduces to finding D11 (symmetric) and D12 (|D12|=n-1, 0 in D12)
such that the V1V1 and V2V2 constraints are both satisfied.
""")
