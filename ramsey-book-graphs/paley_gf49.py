"""
Construct the Paley-type graph for n=25 (m=49=7^2).

Since 49 is a prime power ≡ 1 (mod 4), the Paley construction works:
- GF(49) = GF(7^2) = GF(7)[x]/(irreducible polynomial of degree 2)
- The quadratic residues in GF(49)* form a set of size (49-1)/2 = 24
- D11 = {elements of Z_49 corresponding to QR in GF(49)}, mapped via some bijection

For GF(p^2) where p is prime:
- Elements can be represented as a + b*alpha where alpha^2 = g for some non-residue g mod p
- Alternatively, use the polynomial ring GF(7)[x]/(x^2 + 1) if -1 is a non-residue mod 7
  (it is: 7 ≡ 3 mod 4, so -1 is a QNR mod 7)
- Elements: {a + b*i : a,b in GF(7)} where i^2 = -1, i.e., i^2 ≡ 6 mod 7

Actually for the Ramsey construction we need to work with the ADDITIVE group Z_49,
but the difference set D11 comes from the multiplicative structure of GF(49).

Key point: GF(49) ≠ Z_49 as a ring! GF(49) = GF(7)[x]/(x^2+1) has characteristic 7.
We need a way to embed the Paley construction into a circulant graph on Z_49.

One approach (from Lidicky et al.):
- For m = p^k with p prime, use the additive group of GF(p^k) ≅ (Z_p)^k
- But our circulant is on Z_m, which is NOT the same as (Z_p)^k when k > 1
- Z_49 ≅ Z_49 (cyclic), but GF(49)+ ≅ Z_7 x Z_7 (not cyclic!)

This means we can't directly use the Paley construction on Z_49.
We need a different approach.

Wait -- let me reconsider. The problem says "2-block circulant" graphs.
But for prime power m ≡ 1 mod 4, the known construction might use a
different graph structure (not necessarily circulant on Z_m).

Let me check: for m = p^2, the Paley graph on GF(p^2) has p^2 vertices
with adjacency defined by quadratic residues. This gives a strongly regular
graph. The Ramsey lower bound construction might use this directly.

Actually, re-reading the problem statement: "The lower bound colorings are
2-block circulant graphs and particular generalizations of Paley graphs."

So for prime powers ≡ 1 mod 4, the construction IS the Paley graph (or
a generalization), not necessarily a circulant on Z_m.

For our framework (which uses Z_m circulants), we need to find a way to
embed the GF(49) Paley construction into Z_49. One approach:

Since 49 = 7^2, we can try:
1. Use a bijection between Z_49 and GF(49) that preserves some structure
2. Or: find D11 ⊂ Z_49 such that the circulant graph "simulates" the Paley graph

Actually, there's a standard trick: for m = p^2, choose a primitive root g of Z_{p^2}*.
The units of Z_{p^2} have order phi(p^2) = p(p-1) = 42.
The quadratic residues in Z_{p^2}* are {g^{2k} : k = 0,...,20}, which form a subgroup of index 2.
These QR mod 49 can serve as D11!

Wait, but Z_49* is cyclic of order 42 = 2 * 3 * 7, while GF(49)* is cyclic of order 48 = 2^4 * 3.
These are different groups. The QR in Z_49* and QR in GF(49)* are different.

For the Paley construction to work for Ramsey book graphs, what matters is:
m ≡ 1 mod 4, and we need QR (or some difference set) such that the resulting
graph avoids the required book subgraphs.

Let me just try the Z_49 QR approach (quadratic residues mod 49) and see if it works.
"""

import sys, os, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ramsey_core import BlockCirculantGraph, verify_construction, Delta, Sigma

m = 49
n = 25
N = 4 * n - 2  # = 98

print(f"n={n}, m={m}, N={N}")
print(f"Need: R(B_24, B_25) >= 4*25-1 = 99, so graph on {N}=98 vertices")

# Z_49 structure
# Z_49* has order phi(49) = 42
# Find primitive root mod 49

def primitive_root_mod(m):
    """Find a primitive root mod m."""
    phi_m = m
    factors_m = []
    temp = m
    for p in range(2, int(temp**0.5) + 2):
        if temp % p == 0:
            factors_m.append(p)
            phi_m = phi_m * (p - 1) // p
            while temp % p == 0:
                temp //= p
    if temp > 1:
        factors_m.append(temp)
        phi_m = phi_m * (temp - 1) // temp

    phi_factors = []
    temp_phi = phi_m
    for p in range(2, int(temp_phi**0.5) + 2):
        if temp_phi % p == 0:
            phi_factors.append(p)
            while temp_phi % p == 0:
                temp_phi //= p
    if temp_phi > 1:
        phi_factors.append(temp_phi)

    for g in range(2, m):
        if math.gcd(g, m) != 1:
            continue
        ok = True
        for p in phi_factors:
            if pow(g, phi_m // p, m) == 1:
                ok = False
                break
        if ok:
            return g, phi_m
    return None, phi_m

g, phi = primitive_root_mod(m)
print(f"Primitive root mod {m}: g = {g}, phi({m}) = {phi}")

# Quadratic residues mod 49 (from units)
qr_49 = set()
for x in range(1, m):
    if math.gcd(x, m) == 1:
        qr_49.add(pow(x, 2, m))

print(f"QR mod {m}: {sorted(qr_49)} ({len(qr_49)} elements)")

# Non-residues (among units)
units = {x for x in range(1, m) if math.gcd(x, m) == 1}
qnr_49 = units - qr_49
print(f"QNR mod {m}: {sorted(qnr_49)} ({len(qnr_49)} elements)")

# Non-units (multiples of 7): {7, 14, 21, 28, 35, 42}
non_units = {x for x in range(1, m) if math.gcd(x, m) > 1}
print(f"Non-units: {sorted(non_units)} ({len(non_units)} elements)")

# Check: is -1 a QR mod 49?
# -1 mod 49 = 48. Is 48 in QR?
minus_one = m - 1
print(f"-1 mod {m} = {minus_one}, is QR: {minus_one in qr_49}")
# Since m = 49 ≡ 1 mod 4, -1 should be a QR.

# The Paley graph on Z_49 uses QR as the adjacency set.
# But QR only covers units. We need to also handle the 6 non-units.
# Total nonzero elements: 48 = 42 units + 6 non-units.

# For D11 with Paley-type: D11 = QR ∪ (some subset of non-units)
# Need D11 symmetric and |D11| = n = 25 (but wait, n=25 is odd -- same parity issue!)
# Actually m=49 ≡ 1 mod 4 is the "Paley works" case, so there should be a construction.

# For m ≡ 1 mod 4 (Paley case), the original result says the Paley graph itself works.
# The Paley graph uses D11 = QR. But |QR| = 21 (only units that are QR).
# We need |D11| = (m-1)/2 = 24 for the full circulant.
# 21 QR + how many of 6 non-units? Need 24 - 21 = 3 non-units.

# The non-units are {7, 14, 21, 28, 35, 42}.
# Their negations: -7=42, -14=35, -21=28. So pairs: {7,42}, {14,35}, {21,28}.
# Adding all 3 pairs = 6 elements -> 21 + 6 = 27, too many.
# Adding 0 pairs -> 21 elements, need 24.
# We need exactly 3 non-units in D11, forming 1.5 pairs -- impossible since D11 must be symmetric!

# Hmm. QR has 21 elements. Non-units: 6. Total nonzero: 48. Need |D11| = 24.
# D11 = QR ∪ (3 elements from non-units). But 3 non-units can't be symmetric
# (pairs are of size 2). So either 0, 2, 4, or 6 non-units.
# |D11| = 21 + 0 = 21, or 21 + 2 = 23, or 21 + 4 = 25, or 21 + 6 = 27.
# None gives 24! But wait -- we could also use QNR.

# Alternative: D11 = QNR ∪ (some non-units)
# |QNR| = 21. Same issue.

# Or: D11 takes some QR and some QNR, plus possibly non-units.
# Actually for the Paley approach with m = p^2, the right thing is probably
# to use the Paley graph on GF(p^2), not on Z_{p^2}.

# Let me try GF(49) approach.
# GF(49) = GF(7)[x]/(x^2 + 1) (since -1 is QNR mod 7)
# Elements: a + b*i where i^2 = -1 (mod 7), a,b in {0,1,...,6}
# Total: 49 elements. Nonzero: 48. QR in GF(49)*: 24 elements.
# GF(49)* is cyclic of order 48.

print("\n" + "=" * 80)
print("GF(49) PALEY CONSTRUCTION")
print("=" * 80)

# GF(49) = GF(7)[i] where i^2 = -1 (= 6 mod 7)
# Elements: (a, b) representing a + b*i, a,b in Z_7
# Multiplication: (a+bi)(c+di) = (ac - bd) + (ad + bc)i

class GF49:
    def __init__(self, a, b):
        self.a = a % 7
        self.b = b % 7

    def __mul__(self, other):
        return GF49(self.a * other.a - self.b * other.b,
                     self.a * other.b + self.b * other.a)

    def __eq__(self, other):
        return self.a == other.a and self.b == other.b

    def __hash__(self):
        return hash((self.a, self.b))

    def __repr__(self):
        return f"({self.a}+{self.b}i)"

    def is_zero(self):
        return self.a == 0 and self.b == 0

    def to_int(self):
        """Map to {0, ..., 48} via a + 7*b."""
        return self.a + 7 * self.b

    @staticmethod
    def from_int(n):
        return GF49(n % 7, n // 7)

# Enumerate all nonzero elements and find a generator
nonzero = [GF49(a, b) for a in range(7) for b in range(7) if not (a == 0 and b == 0)]
print(f"GF(49)* has {len(nonzero)} elements")

# Find a generator (primitive element)
def order_gf49(elem):
    """Compute multiplicative order of elem in GF(49)*."""
    if elem.is_zero():
        return 0
    power = elem
    for k in range(1, 49):
        if power.a == 1 and power.b == 0:
            return k
        power = power * elem
    return -1

gen = None
for elem in nonzero:
    if order_gf49(elem) == 48:
        gen = elem
        break

print(f"Generator of GF(49)*: {gen} with order {order_gf49(gen)}")

# Compute powers of generator
powers = [None] * 48
power = GF49(1, 0)
for k in range(48):
    powers[k] = power
    power = power * gen

# QR in GF(49)* = {g^{2k} : k = 0,...,23} (even powers)
qr_gf49 = {powers[2*k] for k in range(24)}
qnr_gf49 = {powers[2*k+1] for k in range(24)}

print(f"|QR in GF(49)*| = {len(qr_gf49)}")
print(f"|QNR in GF(49)*| = {len(qnr_gf49)}")

# Map GF(49) elements to integers 0..48
# Using the natural map: a + b*i -> a + 7*b
qr_ints = {e.to_int() for e in qr_gf49}
qnr_ints = {e.to_int() for e in qnr_gf49}

print(f"QR as integers: {sorted(qr_ints)}")
print(f"QNR as integers: {sorted(qnr_ints)}")

# For the Paley graph, D11 = QR (the 24 quadratic residues in GF(49)*)
# Since m = 49 ≡ 1 mod 4, -1 is a QR in GF(49)*, so QR is symmetric.
minus_one_gf = GF49(6, 0)  # -1 = 6 in GF(7)
print(f"-1 in GF(49): {minus_one_gf}, is QR: {minus_one_gf in qr_gf49}")

# Check symmetry: for each x in QR, is -x also in QR?
D11_paley = qr_ints
D11_paley_sym = all((-x % 49 if x != 0 else 0) in D11_paley or x == 0 for x in D11_paley)
print(f"D11 (QR) is symmetric under Z_49 negation: {D11_paley_sym}")

# WAIT: -x in Z_49 is NOT the same as -x in GF(49)!
# In GF(49): -x means negation in the field, which is (7-a, 7-b) for x = (a,b)
# In Z_49: -x means 49 - x
# The map a + 7b -> integer is additive over GF(7) addition but NOT over Z_49 addition!
# GF(49) has additive group Z_7 x Z_7, not Z_49.

# This is the fundamental issue: GF(49) addition != Z_49 addition.
# The Paley construction on GF(49) gives a graph whose vertex set is Z_7 x Z_7,
# NOT Z_49. So it's not directly a circulant on Z_49.

print("\n--- IMPORTANT: GF(49)+ ≅ Z_7 x Z_7 ≠ Z_49 ---")
print("The Paley graph on GF(49) is NOT a circulant on Z_49.")
print("It's a Cayley graph on the additive group Z_7 x Z_7.")
print("For the 2-block circulant framework (which uses Z_m), we need a different approach.")

# But the problem statement says m=49 ≡ 1 mod 4 is the "Paley works" case.
# This likely means the lower bound graph is the 2-block version of the
# Paley graph on Z_7 x Z_7, not on Z_49.

# Let's verify: construct the Paley-type graph on Z_7 x Z_7
# D11 = QR in GF(49)*
# The graph has edges between (a,b) and (c,d) if (c-a, d-b) [field subtraction] ∈ QR

# For the book graph Ramsey number, we need:
# Vertices: 2 * m = 2 * 49 = 98 = N = 4n-2 ✓ (n=25)
# Graph: 2-block structure with V1 = GF(49), V2 = GF(49)

# Actually wait. Re-reading the problem setup more carefully.
# The original approach uses a 2-block circulant on Z_m with m = 2n-1.
# For m = prime power ≡ 1 mod 4, the PALEY construction works differently:
# It uses the Paley graph/tournament to define the coloring.

# For the Paley construction on a prime p ≡ 1 mod 4:
# D11 = QR mod p. This works because Z_p = GF(p), so there's no issue.
# For p^k ≡ 1 mod 4 with k > 1: must use GF(p^k) Cayley graph, NOT Z_{p^k} circulant.

# So for n=25 (m=49), we can't use the Z_49 circulant framework directly.
# We need to either:
# 1. Work with the Z_7 x Z_7 Cayley graph (requires generalizing the validator)
# 2. Find an equivalent circulant on Z_49 (may not exist)
# 3. Use SA search on Z_49 instead

# Let's verify the Paley construction on GF(49) using brute force.
# This means: the graph on 2*49 = 98 vertices with:
# V1 = GF(49) vertices {0,...,48} (using a+7b encoding)
# V2 = GF(49) vertices {49,...,97}
# Edges within V1: u ~ v iff (v-u) in QR (field subtraction, i.e., GF(49) difference is QR)
# Edges within V2: same
# Cross edges: will need to determine D12

# For the Paley case (m ≡ 1 mod 4), the standard construction has:
# D11 = D22 = QR, and D12 = QR ∪ {0}
# (Since -1 is a QR, both D11 and D22 are symmetric)
# |D11| = |D22| = (m-1)/2 = 24, |D12| = (m-1)/2 + 1 = 25 = n
# d1 = d2 = 24 + 25 = 49 = m

# Wait, that's different from our usual setup where D22 = complement(D11)!
# For Paley: D11 = D22 = QR (same set), not complements.

# Let me check what sizes this gives:
print(f"\nPaley standard: D11 = D22 = QR, D12 = QR ∪ {{0}}")
print(f"|D11| = {len(qr_ints)}, |D22| = {len(qr_ints)}, |D12| = {len(qr_ints) + 1}")
print(f"d1 = {len(qr_ints) + len(qr_ints) + 1}, d2 = {len(qr_ints) + len(qr_ints) + 1}")

# Hmm, d1 = d2 = 24 + 25 = 49. For our framework: d1 = |D11| + |D12| = 24 + 25 = 49 = m.
# This matches the pattern d1 = m for n >= 12.

# But this uses D11 = D22, not D22 = complement(D11).
# For m ≡ 1 mod 4: the complement of QR is QNR, and |QNR| = (m-1)/2 = 24 too.
# So complement(D11) has the same size as D11, but is a DIFFERENT set.

# The Paley construction for m ≡ 1 mod 4 is fundamentally different from
# the m ≡ 3 mod 4 case! In the Paley case:
# - D11 = D22 (both are QR)
# - This is a SYMMETRIC construction where V1 and V2 have identical red neighborhoods
# vs. the m ≡ 3 mod 4 case where D22 = complement(D11).

# Let me verify this gives a valid construction.
# We need to work in GF(49) arithmetic, not Z_49 arithmetic.

print("\n" + "=" * 80)
print("VERIFYING PALEY CONSTRUCTION ON GF(49)")
print("=" * 80)

# Build lookup: integer -> GF49 element
int_to_gf = {i: GF49.from_int(i) for i in range(49)}
gf_to_int = {v: k for k, v in int_to_gf.items()}

# GF(49) subtraction: a - b in the field
def gf_sub(a_int, b_int):
    a = int_to_gf[a_int]
    b = int_to_gf[b_int]
    result = GF49((a.a - b.a) % 7, (a.b - b.b) % 7)
    return result.to_int()

def gf_add(a_int, b_int):
    a = int_to_gf[a_int]
    b = int_to_gf[b_int]
    result = GF49((a.a + b.a) % 7, (a.b + b.b) % 7)
    return result.to_int()

# D11 = D22 = QR in GF(49)*
D11 = qr_ints
D22 = qr_ints
D12 = qr_ints | {0}

print(f"D11 = QR: {sorted(D11)}")
print(f"|D11| = {len(D11)}, |D12| = {len(D12)}, |D22| = {len(D22)}")

# Verify construction using GF(49) arithmetic
# For each pair of vertices, compute common neighbors

red_threshold = n - 2  # = 23
blue_threshold = n - 1  # = 24

d1 = len(D11) + len(D12)  # degree of V1 vertices
d2 = len(D22) + len(D12)  # degree of V2 vertices

print(f"d1 = {d1}, d2 = {d2}")
print(f"Red threshold: {red_threshold}, Blue threshold: {blue_threshold}")

# V1V1: for each "difference" d in GF(49)\{0}, count common neighbors
# Common neighbors of u, v in V1 where v - u = d (field subtraction):
# w in V1: (w-u) in D11 and (w-v) in D11 -> (w-u) in QR and (w-u-d) in QR
#   count = |{x in QR : x-d in QR}| = autocorrelation of QR at d
# w in V2: (w-u) in D12 and (w-v) in D12 -> similar with D12

# But subtraction is in GF(49), not Z_49!
# So "Delta" must use GF(49) arithmetic.

def gf_delta(A, B, d, m_unused):
    """Count |{a in A : a - d in B}| using GF(49) arithmetic."""
    count = 0
    for a in A:
        diff = gf_sub(a, d)
        if diff in B:
            count += 1
    return count

def gf_sigma(A, B, d, m_unused):
    """Count |{a in A : d - a in B}| using GF(49) arithmetic."""
    count = 0
    for a in A:
        diff = gf_sub(d, a)
        if diff in B:
            count += 1
    return count

max_red = 0
max_blue = 0
violations = []

# V1V1
print("\nV1V1 check:")
for d_int in range(1, 49):
    # d must be nonzero in GF(49)
    if d_int == 0:
        continue
    common = gf_delta(D11, D11, d_int, 49) + gf_delta(D12, D12, d_int, 49)

    if d_int in D11:
        max_red = max(max_red, common)
        if common > red_threshold:
            violations.append(("V1V1_red", d_int, common))
    else:
        blue_common = (98 - 2) - d1 - d1 + common
        max_blue = max(max_blue, blue_common)
        if blue_common > blue_threshold:
            violations.append(("V1V1_blue", d_int, blue_common))

print(f"  Max red common (V1V1): {max_red} (threshold {red_threshold})")

# V2V2 (same as V1V1 since D11 = D22)
print("\nV2V2 check (same as V1V1 by symmetry)")
D12T = {gf_sub(0, x) for x in D12}  # D12^T using field negation
for d_int in range(1, 49):
    common = gf_delta(D22, D22, d_int, 49) + gf_delta(D12T, D12T, d_int, 49)

    if d_int in D22:
        max_red = max(max_red, common)
        if common > red_threshold:
            violations.append(("V2V2_red", d_int, common))
    else:
        blue_common = (98 - 2) - d2 - d2 + common
        max_blue = max(max_blue, blue_common)
        if blue_common > blue_threshold:
            violations.append(("V2V2_blue", d_int, blue_common))

print(f"  Max red common (V2V2): {max_red}")

# V1V2
print("\nV1V2 check:")
for d_int in range(49):
    common = gf_sigma(D11, D12, d_int, 49) + gf_delta(D12, D22, d_int, 49)

    if d_int in D12:
        max_red = max(max_red, common)
        if common > red_threshold:
            violations.append(("V1V2_red", d_int, common))
    else:
        blue_common = (98 - 2) - d1 - d2 + common
        max_blue = max(max_blue, blue_common)
        if blue_common > blue_threshold:
            violations.append(("V1V2_blue", d_int, blue_common))

print(f"  Max red common overall: {max_red}")
print(f"  Max blue common overall: {max_blue}")
print(f"  Violations: {len(violations)}")

if violations:
    print(f"  INVALID! Violations:")
    for vtype, d, val in violations[:20]:
        print(f"    {vtype} d={d} value={val}")
else:
    print(f"  VALID! R(B_24, B_25) >= 99 = 4*25-1")

    # Save the construction
    result = {
        "n": n,
        "target_vertices": N,
        "solver": "Paley-GF49",
        "field": "GF(49) = GF(7)[i], i^2=-1",
        "encoding": "a+bi -> a+7*b",
        "construction": "D11 = D22 = QR(GF(49)*), D12 = QR ∪ {0}",
        "parameters": {
            "m": m,
            "D11": sorted(D11),
            "D12": sorted(D12),
            "D22": sorted(D22)
        },
        "note": "Uses GF(49) arithmetic, NOT Z_49. Differences are field subtraction.",
        "verification": {
            "max_red_common": max_red,
            "max_blue_common": max_blue,
            "red_threshold": red_threshold,
            "blue_threshold": blue_threshold,
            "valid": len(violations) == 0
        }
    }

    outpath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "solution_n25_paley.json")
    with open(outpath, 'w') as f:
        json.dump(result, f, indent=2)
    print(f"  Saved to {outpath}")
