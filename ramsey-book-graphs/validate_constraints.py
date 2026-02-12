#!/usr/bin/env python3
"""Validate constraint model against ramsey_core.verify_construction."""
import json, sys
sys.path.insert(0, '.')
from ramsey_core import BlockCirculantGraph, verify_construction, Delta, Sigma

with open('p31_solutions.json') as f:
    data = json.load(f)

sol = data['solutions'][0]
p = 31
n = 16
D11 = set(sol['D11'])
D12 = set(sol['D12'])
D22 = set(range(1, p)) - D11

print(f"|D11|={len(D11)}, |D12|={len(D12)}, |D22|={len(D22)}")
d1 = len(D11) + len(D12)
d2 = len(D22) + len(D12)
N = 2 * p
print(f"d1={d1}, d2={d2}, N={N}")

# Ground truth
G = BlockCirculantGraph(n=n, D11=D11, D12=D12, D22=D22)
result = verify_construction(G)
print(f"\nramsey_core says: valid={result.valid}")
print(f"  max_red={result.max_red_common}, max_blue={result.max_blue_common}")
print(f"  red_thresh={result.red_threshold}, blue_thresh={result.blue_threshold}")

# My derived constraints
print("\n--- My constraint model ---")
# V1V1: common = A(d) + B(d) = Delta(D11,D11,d) + Delta(D12,D12,d)
# V2V2: common = C(d) + B'(d) where B'(d) = Delta(D12T,D12T,d) = B(p-d)
# V1V2: common = Sigma(D11,D12,d) + Delta(D12,D22,d)

# Blue common = (N-2) - deg_u - deg_v + red_common
# V1V1 blue: (N-2) - 2*d1 + (A+B) = A+B - 2  => need A+B-2 <= n-1=15 => A+B <= 17
# V2V2 blue: (N-2) - 2*d2 + (C+B') = C+B'+2  => need C+B'+2 <= n-1=15 => C+B' <= 13
# V1V2 blue: (N-2) - d1 - d2 + X = X          => need X <= n-1=15

D12T = {(-x) % p for x in D12}
all_ok = True

for d in range(1, p):
    A_d = Delta(D11, D11, d, p)
    B_d = Delta(D12, D12, d, p)
    C_d = Delta(D22, D22, d, p)
    Bp_d = Delta(D12T, D12T, d, p)  # = B(p-d)
    B_pd_check = Delta(D12, D12, p - d, p)

    assert Bp_d == B_pd_check, f"B'(d) != B(p-d) at d={d}: {Bp_d} vs {B_pd_check}"

    rc_v1v1 = A_d + B_d  # red common for V1V1
    rc_v2v2 = C_d + Bp_d  # red common for V2V2

    if d in D11:
        # V1V1 red: rc <= n-2=14
        if rc_v1v1 > n - 2:
            print(f"  FAIL V1V1 red d={d}: A+B={rc_v1v1} > {n-2}")
            all_ok = False
        # V2V2 blue: C+B' <= 13
        blue_v2v2 = rc_v2v2 + 2  # = (N-2)-2*d2 + rc_v2v2
        if blue_v2v2 > n - 1:
            print(f"  FAIL V2V2 blue d={d}: C+B'={rc_v2v2}, blue={blue_v2v2} > {n-1}")
            all_ok = False
    else:
        # d in D22
        # V2V2 red: rc <= n-2=14
        if rc_v2v2 > n - 2:
            print(f"  FAIL V2V2 red d={d}: C+B'={rc_v2v2} > {n-2}")
            all_ok = False
        # V1V1 blue: A+B <= 17
        blue_v1v1 = rc_v1v1 - 2
        if blue_v1v1 > n - 1:
            print(f"  FAIL V1V1 blue d={d}: A+B={rc_v1v1}, blue={blue_v1v1} > {n-1}")
            all_ok = False

# V1V2 cross
for d in range(p):
    sigma_val = Sigma(D11, D12, d, p)
    delta_val = Delta(D12, D22, d, p)
    X_d = sigma_val + delta_val

    if d in D12:
        if X_d > n - 2:
            print(f"  FAIL V1V2 red d={d}: X={X_d} > {n-2}")
            all_ok = False
    else:
        blue_v1v2 = X_d  # = (N-2)-d1-d2+X = X
        if blue_v1v2 > n - 1:
            print(f"  FAIL V1V2 blue d={d}: X={X_d} > {n-1}")
            all_ok = False

if all_ok:
    print("ALL CONSTRAINTS PASS - model is correct!")
else:
    print("SOME CONSTRAINTS FAILED")
