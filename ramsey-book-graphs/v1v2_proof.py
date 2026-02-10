"""
Investigate WHY V1V2 common neighbors are perfectly flat.

Key discovery: For ALL known constructions, the V1V2 lambda values are:
  - Exactly n-2 for all red edges (d in D12)
  - Exactly n-1 for all blue edges (d not in D12)

This means V1V2 constraints are AUTOMATICALLY satisfied when:
  Sigma(D11, D12, d, m) + Delta(D12, D22, d, m) = n-2 for all d in D12
  Sigma(D11, D12, d, m) + Delta(D12, D22, d, m) = n-1 for all d not in D12

This is a bilinear condition on (D11, D12, D22) that should have an algebraic explanation.
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ramsey_core import BlockCirculantGraph, Delta, Sigma
import json

# Load all known constructions
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


def indicator_convolution(S, d, m):
    """Compute sum_{s in S} 1_{S}(d - s) = |{(a,b) in S x S : a + b = d}| = Sigma(S, S, d, m)."""
    return Sigma(S, S, d, m)


print("V1V2 FLATNESS PROOF INVESTIGATION")
print("=" * 80)

for n in sorted(KNOWN.keys()):
    m = 2 * n - 1
    D11 = KNOWN[n]["D11"]
    D12 = KNOWN[n]["D12"]
    D22 = KNOWN[n]["D22"]

    # V1V2 common neighbors for difference d:
    # f(d) = Sigma(D11, D12, d, m) + Delta(D12, D22, d, m)
    #
    # Sigma(D11, D12, d, m) = |{a in D11 : (d-a) in D12}|
    # Delta(D12, D22, d, m) = |{b in D12 : (b-d) in D22}|
    #
    # Since D22 = comp(D11) in {1,...,m-1}:
    # (b-d) in D22 iff (b-d) != 0 and (b-d) not in D11
    # i.e., b != d and (b-d) not in D11
    #
    # Delta(D12, D22, d, m) = |{b in D12 : b != d and (b-d) mod m not in D11}|
    #                       = |{b in D12 : b != d}| - |{b in D12 : b != d and (b-d) mod m in D11}|
    #
    # If d in D12: = (|D12| - 1) - Delta'(D12\{d}, D11, d, m)
    # If d not in D12: = |D12| - Delta(D12, D11, d, m)

    print(f"\nn={n}, m={m}")

    # Let's compute each term separately
    for d in range(m):
        sigma_term = Sigma(D11, D12, d, m)
        delta_term = Delta(D12, D22, d, m)
        total = sigma_term + delta_term

        # Also compute Sigma(D11, D12, d) and Delta(D12, comp(D11), d)
        # Delta(D12, D11, d) = |{b in D12 : (b-d) in D11}|
        delta_d12_d11 = Delta(D12, D11, d, m)

        # Delta(D12, D22, d) = |{b in D12: (b-d) in D22}|
        # Since D22 = {1,...,m-1}\D11:
        # = |{b in D12: (b-d) mod m != 0 and (b-d) mod m not in D11}|
        # = (number of b in D12 with b != d) - (number of b in D12, b!=d, with (b-d) in D11)
        # If d in D12: = (|D12|-1) - (delta_d12_d11 - (1 if d in D12 and 0 in D11 else 0))
        # Actually Delta(D12, D11, d) counts b in D12 such that (b-d) in D11
        # If b=d, then (b-d)=0, and 0 is NOT in D11, so we don't need to worry about the b=d case

        d_in_D12 = d in D12
        num_b_neq_d = (len(D12) - 1) if d_in_D12 else len(D12)
        # delta_d12_d11 already excludes b=d contribution since (0) not in D11

        reconstructed_delta_D22 = num_b_neq_d - delta_d12_d11
        assert reconstructed_delta_D22 == delta_term, f"Reconstruction failed at d={d}"

    # So f(d) = Sigma(D11, D12, d) + (|D12| - [d in D12]) - Delta(D12, D11, d)
    # = Sigma(D11, D12, d) - Delta(D12, D11, d) + |D12| - [d in D12]

    # Now: Sigma(D11, D12, d) = |{a in D11 : (d-a) in D12}| = sum_{a in D11} 1_{D12}(d-a)
    # Delta(D12, D11, d) = |{b in D12 : (b-d) in D11}| = sum_{b in D12} 1_{D11}(b-d)
    #                    = |{b in D12 : (b-d) in D11}|

    # Let's define g(d) = Sigma(D11, D12, d) - Delta(D12, D11, d)
    # = sum_{a in D11} 1_{D12}(d-a) - sum_{b in D12} 1_{D11}(b-d)
    # Note: Delta(D12, D11, d) = |{b in D12 : (b-d) in D11}|
    #     = sum_{b in D12} 1_{D11}(b - d)
    # And Sigma(D11, D12, d) = sum_{a in D11} 1_{D12}(d - a)
    #
    # These are related! If we let c = d-a for the first sum and c = b-d for the second:
    # Sigma = sum_{c} 1_{D11}(d-c) * 1_{D12}(c)  [c = d-a, a = d-c]
    # Delta = sum_{c} 1_{D12}(d+c) * 1_{D11}(c)   [c = b-d, b = d+c]
    #       = sum_{c} 1_{D12}(d+c) * 1_{D11}(c)

    # Hmm, let me just check that f(d) = n-2 + [d not in D12]
    vals = []
    for d in range(m):
        f = Sigma(D11, D12, d, m) + Delta(D12, D22, d, m)
        expected = (n - 2) + (0 if d in D12 else 1)
        if f != expected:
            print(f"  MISMATCH at d={d}: f={f}, expected={expected}")
        vals.append(f)

    unique_vals = set(vals)
    print(f"  f(d) values: {sorted(unique_vals)}")
    print(f"  Confirmed: f(d) = {n-2} for d in D12, f(d) = {n-1} for d not in D12: "
          f"{all(vals[d] == (n-2 if d in D12 else n-1) for d in range(m))}")

    # Now check if this is a consequence of size constraints + complement property
    # f(d) = Sigma(D11, D12, d) + |D12| - [d in D12] - Delta(D12, D11, d)
    #
    # So f(d) = n-2 + [d not in D12] means:
    # Sigma(D11, D12, d) - Delta(D12, D11, d) + |D12| - [d in D12] = n - 2 + [d not in D12]
    # Sigma(D11, D12, d) - Delta(D12, D11, d) = n - 2 + [d not in D12] - |D12| + [d in D12]
    # = n - 2 - |D12| + 1 = n - 1 - |D12|

    target = n - 1 - len(D12)
    print(f"  This requires: Sigma(D11, D12, d) - Delta(D12, D11, d) = {target} for all d")

    # Verify
    for d in range(m):
        s = Sigma(D11, D12, d, m)
        delta = Delta(D12, D11, d, m)
        diff = s - delta
        if diff != target:
            print(f"  FAILED at d={d}: Sigma={s}, Delta={delta}, diff={diff}")
            break
    else:
        print(f"  VERIFIED: Sigma(D11, D12, d) - Delta(D12, D11, d) = {target} for ALL d")

    # Now: Sigma(D11, D12, d) = sum_{a in D11} 1_{D12}(d - a)
    # Delta(D12, D11, d) = sum_{b in D12} 1_{D11}(b - d) = sum_{b in D12} 1_{D11}(-(d - b))
    # Since D11 is SYMMETRIC: 1_{D11}(-(d-b)) = 1_{D11}(d-b)
    # So Delta(D12, D11, d) = sum_{b in D12} 1_{D11}(d - b) = Sigma(D12, D11, d) [using sum convention]
    # Wait no. Let me be careful.
    # Delta(A, B, d, m) = |{a in A : (a - d) mod m in B}| = sum_{a in A} 1_B(a - d)
    # Sigma(A, B, d, m) = |{a in A : (d - a) mod m in B}| = sum_{a in A} 1_B(d - a)
    #
    # So Delta(D12, D11, d) = sum_{b in D12} 1_{D11}(b - d)
    # And Sigma(D11, D12, d) = sum_{a in D11} 1_{D12}(d - a)
    #
    # Since D11 is symmetric: 1_{D11}(x) = 1_{D11}(-x)
    # So 1_{D11}(b - d) = 1_{D11}(d - b)
    # Therefore Delta(D12, D11, d) = sum_{b in D12} 1_{D11}(d - b) = Sigma(D12, D11, d)
    #
    # So the condition becomes:
    # Sigma(D11, D12, d) - Sigma(D12, D11, d) = constant for all d
    #
    # In Fourier terms:
    # Sigma(A, B, d) = sum_{a in A} 1_B(d-a) is the convolution (1_A * 1_B)(d)
    # So Sigma(D11, D12, d) = (1_{D11} * 1_{D12})(d)
    # And Sigma(D12, D11, d) = (1_{D12} * 1_{D11})(d)
    #
    # CONVOLUTION IS COMMUTATIVE! So Sigma(D11, D12, d) = Sigma(D12, D11, d) for all d!
    # This means the condition is: 0 = constant = n - 1 - |D12|
    # i.e., |D12| = n - 1!

    print(f"  |D12| = {len(D12)}, n-1 = {n-1}, match: {len(D12) == n - 1}")

print("\n" + "=" * 80)
print("CONCLUSION: V1V2 FLATNESS IS AUTOMATIC!")
print("=" * 80)
print("""
The V1V2 common neighbor count f(d) = Sigma(D11, D12, d) + Delta(D12, D22, d)
simplifies when D22 = complement(D11):

  f(d) = Sigma(D11, D12, d) - Delta(D12, D11, d) + |D12| - [d in D12]

Since D11 is symmetric, Delta(D12, D11, d) = Sigma(D12, D11, d).
Since convolution is commutative, Sigma(D11, D12, d) = Sigma(D12, D11, d).
Therefore: Sigma(D11, D12, d) - Delta(D12, D11, d) = 0 for all d.

So f(d) = |D12| - [d in D12].

For f(d) = n-2 when d in D12: requires |D12| - 1 = n-2, i.e., |D12| = n-1.
For f(d) = n-1 when d not in D12: requires |D12| = n-1.

THEOREM: If D22 = {1,...,m-1} \\ D11, D11 is symmetric, and |D12| = n-1,
then V1V2 constraints are AUTOMATICALLY satisfied with:
  - Red common neighbors = n-2 (exactly at threshold)
  - Blue common neighbors = n-1 (exactly at threshold)

This reduces the problem to satisfying only V1V1 and V2V2 constraints!
""")

# Furthermore, since D22 = comp(D11), V2V2 is determined by D11 and D12^T
print("Since V2V2 uses D22 = comp(D11) and D12^T:")
print("The ENTIRE construction is determined by just D11 and D12.")
print("And the V1V2 constraints are free. The only work is V1V1 and V2V2.")
