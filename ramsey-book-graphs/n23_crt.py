"""
CRT-based construction for n=23 (m=45 = 9 x 5).

Key idea: Z_45 ≅ Z_9 x Z_5 via CRT.
Element x in Z_45 corresponds to (x mod 9, x mod 5).

Strategy: Try to build D11 as a "product" or "semi-product" of sets in Z_9 and Z_5.
The Paley-like approach works for prime powers, so maybe we can combine
QR-like structures from Z_9 and Z_5.

Z_5 has QR = {1, 4} (since 1^2=1, 2^2=4, 3^2=4, 4^2=1)
Z_9 has: units = {1,2,4,5,7,8}, squares = {1,4,7} (since 1^2=1, 2^2=4, 4^2=7, 5^2=7, 7^2=4, 8^2=1)

So QR(Z_9) = {1, 4, 7} and QNR(Z_9) = {2, 5, 8} (among units)
Non-units in Z_9: {0, 3, 6} (divisible by 3)

This gives us a more nuanced partition of Z_45*:
- (QR_9, QR_5): elements that are QR in both components
- (QR_9, QNR_5): QR in Z_9, QNR in Z_5
- (QNR_9, QR_5): QNR in Z_9, QR in Z_5
- (QNR_9, QNR_5): QNR in both
- (non-unit_9, QR_5): non-unit in Z_9, QR in Z_5
- etc.
"""

import sys, os, math, json
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ramsey_core import Delta, BlockCirculantGraph, verify_construction, save_construction

m = 45
n = 23

def crt(a9, a5):
    """CRT: find x in Z_45 with x ≡ a9 (mod 9), x ≡ a5 (mod 5)."""
    # 45 = 9 * 5. 9^{-1} mod 5 = 4, 5^{-1} mod 9 = 2.
    return (a9 * 5 * 2 + a5 * 9 * 4) % 45

def to_crt(x):
    """Convert x in Z_45 to (x mod 9, x mod 5)."""
    return (x % 9, x % 5)

# Verify CRT
for x in range(45):
    a9, a5 = to_crt(x)
    assert crt(a9, a5) == x

# Classify elements of Z_45*
qr5 = {1, 4}  # QR in Z_5*
qnr5 = {2, 3}  # QNR in Z_5*

qr9 = {1, 4, 7}  # squares mod 9 (among units)
qnr9 = {2, 5, 8}  # non-squares mod 9 (among units)
non_unit9 = {0, 3, 6}  # elements of Z_9 divisible by 3

# Character classes in Z_45
classes = defaultdict(list)
for x in range(1, m):
    a9, a5 = to_crt(x)
    if a5 == 0:
        cls = f"(*,0)"
    elif a9 in non_unit9:
        if a5 in qr5:
            cls = "(nu,R)"
        else:
            cls = "(nu,N)"
    else:
        r9 = "R" if a9 in qr9 else "N"
        r5 = "R" if a5 in qr5 else "N"
        cls = f"({r9},{r5})"
    classes[cls].append(x)

print("Element classification in Z_45 = Z_9 x Z_5:")
for cls in sorted(classes.keys()):
    elts = sorted(classes[cls])
    print(f"  {cls}: {len(elts)} elements: {elts}")

# Symmetry check: is each class symmetric under negation?
print("\nSymmetry check (is -x in same class as x):")
for cls in sorted(classes.keys()):
    elts = set(classes[cls])
    negs = {(-x) % m for x in elts}
    print(f"  {cls}: symmetric = {elts == negs}")


# =============================================================================
# Try various set constructions
# =============================================================================
print("\n" + "=" * 80)
print("TRYING VARIOUS CRT-BASED D11 CONSTRUCTIONS")
print("=" * 80)

def try_construction(name, D11_set, D12_set):
    """Test a specific D11, D12 construction."""
    D11_set = set(D11_set) - {0}
    D22_set = set(range(1, m)) - D11_set

    G = BlockCirculantGraph(n=n, D11=D11_set, D12=D12_set, D22=D22_set)
    result = verify_construction(G)

    n_violations = len(result.violations)
    total_excess = sum(e for _, _, e in result.violations)

    print(f"\n  {name}:")
    print(f"    |D11|={len(D11_set)}, |D12|={len(D12_set)}, |D22|={len(D22_set)}")
    print(f"    valid={result.valid}, violations={n_violations}, total_excess={total_excess}")
    print(f"    max_red={result.max_red_common} (thresh {result.red_threshold})")
    print(f"    max_blue={result.max_blue_common} (thresh {result.blue_threshold})")

    if result.valid:
        print(f"    *** VALID CONSTRUCTION! ***")
        print(f"    D11 = {sorted(D11_set)}")
        print(f"    D12 = {sorted(D12_set)}")
        return True, result
    elif n_violations <= 3:
        print(f"    Violations: {result.violations}")

    return False, result


# Construction 1: D11 = (R,R) ∪ (N,N) (like Paley tensor)
print("\n--- Construction family: tensor-product-like ---")
rr = set(classes["(R,R)"])
nn = set(classes["(N,N)"])
rn = set(classes["(R,N)"])
nr = set(classes["(N,R)"])
nuR = set(classes["(nu,R)"])
nuN = set(classes["(nu,N)"])
star0 = set(classes.get("(*,0)", []))

# All "QR-like" elements
tensor_qr = rr | nn  # QR in both or QNR in both = "QR in product"
tensor_qnr = rn | nr  # mixed

print(f"\nTensor QR = (R,R)∪(N,N): {len(tensor_qr)} elements")
print(f"Tensor QNR = (R,N)∪(N,R): {len(tensor_qnr)} elements")
print(f"Non-units: (nu,R)={len(nuR)}, (nu,N)={len(nuN)}, (*,0)={len(star0)}")

# Construction 1a: D11 = tensor_QR, D12 = tensor_QR
try_construction("D11=D12=tensor_QR", tensor_qr, tensor_qr)

# Construction 1b: D11 = tensor_QR ∪ non-units, D12 = tensor_QR
try_construction("D11=tensor_QR∪nuR∪nuN, D12=tensor_QR", tensor_qr | nuR | nuN, tensor_qr)

# Construction 1c: D11 = tensor_QR ∪ (*,0), D12 = tensor_QR
try_construction("D11=tensor_QR∪(*,0), D12=tensor_QR", tensor_qr | star0, tensor_qr)

# Construction 2: Try D11 = D12 pattern (like Paley)
# D11 = (R,R) ∪ (N,N) ∪ some non-units
for nu_add in [set(), nuR, nuN, nuR | nuN, star0, nuR | star0, nuN | star0]:
    D11 = tensor_qr | nu_add
    if len(D11) % 2 == 0 and 10 <= len(D11) <= 30:  # reasonable size
        try_construction(f"D11=D12=tQR∪{len(nu_add)}nu ({len(D11)} elts)", D11, D11)

# Construction 3: "Additive character" approach
# Use the principal character of Z_5 to define D11
# D11 = {x : x mod 5 in QR_5} = all elements with QR residue mod 5
d11_r5 = {x for x in range(1, m) if x % 5 in qr5}
d11_n5 = {x for x in range(1, m) if x % 5 in qnr5}
d11_05 = {x for x in range(1, m) if x % 5 == 0}

print(f"\nElements by Z_5 residue:")
print(f"  QR_5: {len(d11_r5)} elements")
print(f"  QNR_5: {len(d11_n5)} elements")
print(f"  0 mod 5: {len(d11_05)} elements")

try_construction("D11=QR_5-residue, D12=QR_5-residue", d11_r5, d11_r5)

# Construction 4: Use Z_9 structure
d11_r9 = {x for x in range(1, m) if x % 9 in qr9}
d11_n9 = {x for x in range(1, m) if x % 9 in qnr9}

print(f"\nElements by Z_9 residue:")
print(f"  QR_9: {len(d11_r9)} elements")
print(f"  QNR_9: {len(d11_n9)} elements")

try_construction("D11=QR_9-residue, D12=QR_9-residue", d11_r9, d11_r9)

# Construction 5: Combined character
# D11 = {x : chi_9(x)*chi_5(x) = 1} where chi is Legendre symbol extended
def chi5(x):
    """Legendre-like symbol for Z_5."""
    a = x % 5
    if a == 0: return 0
    return 1 if a in qr5 else -1

def chi9(x):
    """Legendre-like symbol for Z_9."""
    a = x % 9
    if a in non_unit9: return 0
    return 1 if a in qr9 else -1

d11_product_char = {x for x in range(1, m) if chi9(x) * chi5(x) == 1}
d11_neg_product_char = {x for x in range(1, m) if chi9(x) * chi5(x) == -1}

print(f"\nProduct character D11 (chi9*chi5=+1): {len(d11_product_char)} elements")
print(f"Negative product char (chi9*chi5=-1): {len(d11_neg_product_char)} elements")

try_construction("D11=D12=chi9*chi5=+1", d11_product_char, d11_product_char)
try_construction("D11=D12=chi9*chi5=-1", d11_neg_product_char, d11_neg_product_char)

# Construction 6: SA-refined CRT construction
# Start with a CRT-based D11 and optimize D12 with SA
print("\n" + "=" * 80)
print("CRT-SEEDED SA SEARCH")
print("=" * 80)

import random

def sa_d12_for_fixed_d11(D11_set, max_iter=1000000, seed=0):
    """SA to find D12 given fixed D11."""
    random.seed(seed)
    D22_set = set(range(1, m)) - D11_set

    # Start with various D12 sizes
    for d12_size in [22, 21, 23]:
        if d12_size == 22:
            D12_list = [0] + random.sample(range(1, m), 21)
        elif d12_size == 21:
            D12_list = random.sample(range(1, m), 21)
        else:
            D12_list = [0] + random.sample(range(1, m), 22)
        D12_set = set(D12_list)

        G = BlockCirculantGraph(n=n, D11=D11_set, D12=D12_set, D22=D22_set)
        result = verify_construction(G)
        current_cost = sum(e for _, _, e in result.violations)
        best_cost = current_cost
        best_D12 = sorted(D12_set)

        T = 5.0
        alpha = 1 - 5.0 / max_iter

        for it in range(max_iter):
            if d12_size == 22:
                D12_others = [x for x in D12_set if x != 0]
            else:
                D12_others = list(D12_set)
            remove = random.choice(D12_others)
            not_in = list(set(range(m)) - D12_set)
            add = random.choice(not_in)

            D12_set.remove(remove)
            D12_set.add(add)

            G = BlockCirculantGraph(n=n, D11=D11_set, D12=D12_set, D22=D22_set)
            result = verify_construction(G)
            new_cost = sum(e for _, _, e in result.violations)

            delta = new_cost - current_cost
            if delta <= 0 or random.random() < math.exp(-delta / max(T, 0.001)):
                current_cost = new_cost
                if new_cost < best_cost:
                    best_cost = new_cost
                    best_D12 = sorted(D12_set)
                    if best_cost == 0:
                        return best_D12, best_cost, d12_size, result
            else:
                D12_set.remove(add)
                D12_set.add(remove)

            T *= alpha

        if best_cost <= 4:
            print(f"    |D12|={d12_size}: best cost = {best_cost}")

    return best_D12, best_cost, 22, None


# Try each CRT-based D11
crt_d11_candidates = [
    ("tensor_QR", tensor_qr),
    ("tensor_QR+nuR", tensor_qr | nuR),
    ("tensor_QR+nuN", tensor_qr | nuN),
    ("tensor_QR+star0", tensor_qr | star0),
    ("tensor_QR+nuR+nuN", tensor_qr | nuR | nuN),
    ("tensor_QNR", tensor_qnr),
    ("tensor_QNR+nuR", tensor_qnr | nuR),
    ("tensor_QNR+nuN", tensor_qnr | nuN),
    ("tensor_QNR+star0", tensor_qnr | star0),
    ("chi9*chi5=+1", d11_product_char),
    ("chi9*chi5=-1", d11_neg_product_char),
    ("QR_5-residue", d11_r5),
    ("QR_9-residue", d11_r9),
    ("QR_5-residue+0mod5", d11_r5 | d11_05),
]

for name, D11_set in crt_d11_candidates:
    D11_set = D11_set - {0}  # Remove 0
    if len(D11_set) < 10 or len(D11_set) > 34:
        continue
    # Check symmetry
    sym = all((-x) % m in D11_set for x in D11_set)
    if not sym:
        # Make symmetric
        D11_set = D11_set | {(-x) % m for x in D11_set}
        D11_set.discard(0)

    print(f"\n  CRT D11 = {name} ({len(D11_set)} elements, sym={sym})")
    d12, cost, d12_size, result = sa_d12_for_fixed_d11(D11_set, max_iter=500000, seed=42)
    print(f"    Best cost = {cost}")
    if cost == 0 and result is not None:
        print(f"    *** VALID! D12 = {d12} ***")
        G = BlockCirculantGraph(n=n, D11=D11_set, D12=set(d12), D22=set(range(1,m))-D11_set)
        filename = "/Users/stephenpadgett/Projects/math/ramsey-book-graphs/solution_n23.json"
        save_construction(G, result, filename, solver=f"CRT-{name}")
        print(f"    SAVED to {filename}")


print("\n" + "=" * 80)
print("DONE")
print("=" * 80)
