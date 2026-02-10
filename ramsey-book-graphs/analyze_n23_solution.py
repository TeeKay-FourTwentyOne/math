"""
Analyze the n=23 solution found by fast SA.

Understand:
1. CRT structure (Z_45 = Z_9 x Z_5 decomposition)
2. Delta value distribution
3. How tight the construction is
4. Whether it matches any algebraic pattern
"""

import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ramsey_core import Delta, Sigma, BlockCirculantGraph, verify_construction

m = 45
n = 23

# Load solution
with open("/Users/stephenpadgett/Projects/math/ramsey-book-graphs/solution_n23.json") as f:
    data = json.load(f)

D11 = set(data["parameters"]["D11"])
D12 = set(data["parameters"]["D12"])
D22 = set(data["parameters"]["D22"])

print(f"n={n}, m={m}, N={2*m}")
print(f"|D11|={len(D11)}, |D12|={len(D12)}, |D22|={len(D22)}")
print(f"D11 = {sorted(D11)}")
print(f"D12 = {sorted(D12)}")
print(f"D22 = {sorted(D22)}")

# Check properties
print(f"\n0 in D12: {0 in D12}")
print(f"D22 == complement(D11): {D22 == set(range(1,m)) - D11}")
print(f"D11 symmetric: {all((-x)%m in D11 for x in D11)}")
print(f"D12 symmetric: {all((-x)%m in D12 for x in D12)}")

# CRT decomposition
def to_crt(x):
    return (x % 9, x % 5)

print("\n" + "=" * 60)
print("CRT DECOMPOSITION (Z_45 = Z_9 x Z_5)")
print("=" * 60)

# QR structure
qr5 = {1, 4}
qnr5 = {2, 3}
qr9 = {1, 4, 7}
qnr9 = {2, 5, 8}
non_unit9 = {0, 3, 6}

print("\nD11 CRT coordinates:")
d11_crt = [(x, to_crt(x)) for x in sorted(D11)]
for x, (a, b) in d11_crt:
    r9 = "R" if a in qr9 else ("N" if a in qnr9 else "0")
    r5 = "R" if b in qr5 else ("N" if b in qnr5 else "0")
    print(f"  {x:2d} -> ({a},{b})  [{r9},{r5}]")

print("\nD12 CRT coordinates:")
d12_crt = [(x, to_crt(x)) for x in sorted(D12)]
for x, (a, b) in d12_crt:
    r9 = "R" if a in qr9 else ("N" if a in qnr9 else "0")
    r5 = "R" if b in qr5 else ("N" if b in qnr5 else "0")
    print(f"  {x:2d} -> ({a},{b})  [{r9},{r5}]")

# Count how many from each character class
from collections import Counter
def classify(x):
    a, b = to_crt(x)
    if b == 0:
        return "(*,0)"
    r5 = "R" if b in qr5 else "N"
    if a in non_unit9:
        return f"(nu,{r5})"
    r9 = "R" if a in qr9 else "N"
    return f"({r9},{r5})"

d11_classes = Counter(classify(x) for x in D11)
d12_classes = Counter(classify(x) for x in D12)
d22_classes = Counter(classify(x) for x in D22)

all_classes = sorted(set(list(d11_classes.keys()) + list(d12_classes.keys()) + list(d22_classes.keys())))

print(f"\nClass distribution:")
print(f"  {'Class':>8s}  {'D11':>4s}  {'D12':>4s}  {'D22':>4s}  {'Total':>5s}")
for cls in all_classes:
    total = sum(1 for x in range(1, m) if classify(x) == cls)
    print(f"  {cls:>8s}  {d11_classes.get(cls,0):4d}  {d12_classes.get(cls,0):4d}  {d22_classes.get(cls,0):4d}  {total:5d}")

# Delta value analysis
print("\n" + "=" * 60)
print("DELTA VALUE ANALYSIS")
print("=" * 60)

A = {}  # Delta(D11,D11,d)
B = {}  # Delta(D12,D12,d)
for d in range(1, m):
    A[d] = Delta(D11, D11, d, m)
    B[d] = Delta(D12, D12, d, m)

print("\nV1V1 common neighbors = A(d) + B(d):")
print(f"  {'d':>3s}  {'A(d)':>5s}  {'B(d)':>5s}  {'A+B':>5s}  {'Edge':>5s}  {'Margin':>7s}")
for d in range(1, m):
    ab = A[d] + B[d]
    if d in D11:
        edge = "RED"
        margin = (n - 2) - ab
    else:
        blue_common = (2*m - 2) - 2*(len(D11) + len(D12)) + ab
        edge = "BLUE"
        margin = (n - 1) - blue_common
    if margin <= 1:
        print(f"  {d:3d}  {A[d]:5d}  {B[d]:5d}  {ab:5d}  {edge:>5s}  {margin:7d} {'<-- TIGHT' if margin == 0 else ''}")

# Count tight constraints
tight_red = sum(1 for d in range(1,m) if d in D11 and A[d]+B[d] == n-2)
tight_blue = sum(1 for d in range(1,m) if d not in D11 and (2*m-2)-2*(len(D11)+len(D12))+A[d]+B[d] == n-1)
print(f"\nTight constraints: {tight_red} red (of {len(D11)}), {tight_blue} blue (of {len(D22)})")

# A(d) distribution by edge type
A_red = [A[d] for d in range(1,m) if d in D11]
A_blue = [A[d] for d in range(1,m) if d not in D11]
B_red = [B[d] for d in range(1,m) if d in D11]
B_blue = [B[d] for d in range(1,m) if d not in D11]

print(f"\nA(d) stats:")
print(f"  Red edges: min={min(A_red)}, max={max(A_red)}, avg={sum(A_red)/len(A_red):.2f}")
print(f"  Blue edges: min={min(A_blue)}, max={max(A_blue)}, avg={sum(A_blue)/len(A_blue):.2f}")
print(f"\nB(d) stats:")
print(f"  Red edges: min={min(B_red)}, max={max(B_red)}, avg={sum(B_red)/len(B_red):.2f}")
print(f"  Blue edges: min={min(B_blue)}, max={max(B_blue)}, avg={sum(B_blue)/len(B_blue):.2f}")

# Overall averages
avg_A = sum(A.values()) / (m-1)
avg_B = sum(B.values()) / (m-1)
print(f"\nOverall averages: A={avg_A:.2f}, B={avg_B:.2f}, A+B={avg_A+avg_B:.2f}")
print(f"Theoretical average: |D11|^2/m + |D12|^2/m = {len(D11)**2/m + len(D12)**2/m:.2f}")

# A(d) + B(d) histogram
from collections import Counter
ab_hist = Counter(A[d] + B[d] for d in range(1, m))
print(f"\nA(d)+B(d) histogram: {dict(sorted(ab_hist.items()))}")

print("\n" + "=" * 60)
print("ANALYSIS COMPLETE")
print("=" * 60)
