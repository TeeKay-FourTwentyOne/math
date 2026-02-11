"""
Product-based construction attempts for n=32, m=63 = 9*7.

Idea: Z_63 = Z_9 x Z_7 via CRT.
GF(9) and GF(7) both support Paley-like constructions.
Try to combine them into a Z_63 construction.

GF(7): q=7, QR = {1,2,4}, QNR = {3,5,6}. q ≡ 3 mod 4.
GF(9): q=9, QR = {1,2,4,5,7,8} (since x^4 = 1 for all x in GF(9)*).
  Wait, GF(9)* has order 8, QR = {x^2 : x in GF(9)*} has size 4.
  Actually for GF(9), we need proper field arithmetic.

Z_9: QR(Z_9*) = {1^2, 2^2, 3^2, 4^2, ...} mod 9 = {1,4,0,7,...}
  QR(Z_9*) = {1, 4, 7} (|QR|=3), QNR(Z_9*) = {2, 5, 8} (|QNR|=3)

Z_7: QR(Z_7*) = {1, 2, 4} (|QR|=3), QNR(Z_7*) = {3, 5, 6} (|QNR|=3)

CRT classes for Z_63* (excluding 0 and elements with component 0):
  (QR_9, QR_7): |QR_9|*|QR_7| = 3*3 = 9 elements
  (QR_9, QNR_7): 3*3 = 9
  (QNR_9, QR_7): 3*3 = 9
  (QNR_9, QNR_7): 3*3 = 9
  Total: 36 (plus zero divisors)

Elements with a=0 mod 9 (7 elements: 0,9,18,27,36,45,54), nonzero: {9,18,27,36,45,54}
  These are {0}*Z_7. The nonzero ones have b in {2,4,6,1,3,5} mod 7
  = {9:2, 18:4, 27:6, 36:1, 45:3, 54:5}
  (0,QR_7): 36,9,18 -> 3 elements
  (0,QNR_7): 27,45,54 -> 3 elements

Elements with b=0 mod 7 (9 elements: 0,7,14,21,28,35,42,49,56), nonzero: {7,14,21,28,35,42,49,56}
  a mod 9: {7:7, 14:5, 21:3, 28:1, 35:8, 42:6, 49:4, 56:2}
  (QR_9,0): 28,49,7 -> 3 elements  (a=1,4,7)
  (QNR_9,0): 14,35,56 -> 3 elements (a=5,8,2)
  There's also a=3,6: 21,42 -> these are divisible by 3 but not QR or QNR!
  Wait, 3^2=0 mod 9, 6^2=0 mod 9. So 3,6 are zero divisors in Z_9.
  QR(Z_9*) = {x^2 mod 9 : gcd(x,9)=1} = {1,4,7}
  Units of Z_9: {1,2,4,5,7,8}. So 3,6 are not units.
  21 mod 9 = 3, 42 mod 9 = 6 -> these are not in QR or QNR, they're non-units.

So the non-unit elements of Z_63 (those where gcd(x,63)>1 and x≠0) are:
  - Multiples of 7 with a not a unit of Z_9: 21 (a=3), 42 (a=6) -> 2 elements
  - Multiples of 9 with b not a unit of Z_7: only b=0, but that's 0 -> none extra
  - Actually, gcd(x,63)>1 iff 3|x or 7|x.
  - 3|x: a mod 3 = 0, i.e., a in {0,3,6}
    - a=0: 0,9,18,27,36,45,54 (7 elements)
    - a=3: 3,12,21,30,39,48,57 (7 elements)
    - a=6: 6,15,24,33,42,51,60 (7 elements)
  - 7|x: 0,7,14,21,28,35,42,49,56 (9 elements)
  - 3|x or 7|x = 7+7+7+9-3 = 27 (subtracting intersections: 0,21,42)
  - So 27 non-coprime elements. Units of Z_63: 63-27 = 36 elements (phi(63) = 63*(1-1/3)*(1-1/7) = 36).
  - Non-unit nonzero: 62-36 = 26 elements.

This is getting complex. Let me just code it up.
"""

import sys, os, math, random, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from n32_fast_numpy import compute_delta_fft, compute_cost_vectorized
from ramsey_core import BlockCirculantGraph, verify_construction, save_construction

M = 63
N_PARAM = 32


def crt_encode(a, b):
    """Map (a mod 9, b mod 7) to Z_63."""
    # inv(7,9) = 4 (since 7*4=28=1 mod 9)
    # inv(9,7) = 4 (since 9*4=36=1 mod 7)
    return (a * 7 * 4 + b * 9 * 4) % 63


def crt_decode(x):
    """Map Z_63 to (a mod 9, b mod 7)."""
    return (x % 9, x % 7)


def classify():
    """Classify Z_63 elements."""
    qr9 = {1, 4, 7}  # QR(Z_9*)
    qnr9 = {2, 5, 8}  # QNR(Z_9*)
    qr7 = {1, 2, 4}  # QR(Z_7*)
    qnr7 = {3, 5, 6}  # QNR(Z_7*)

    classes = {}
    for x in range(63):
        a, b = crt_decode(x)
        if a in qr9:
            atype = 'R'
        elif a in qnr9:
            atype = 'N'
        elif a == 0:
            atype = '0'
        else:  # a in {3, 6} - non-units of Z_9
            atype = 'Z'  # zero divisor

        if b in qr7:
            btype = 'R'
        elif b in qnr7:
            btype = 'N'
        elif b == 0:
            btype = '0'
        else:
            btype = '?'

        label = atype + btype
        classes.setdefault(label, []).append(x)

    return classes


def try_construction(D11_set, D12_set, label=""):
    """Test a construction and report."""
    m = M
    n = N_PARAM
    D11_ind = np.zeros(m, dtype=np.int64)
    D12_ind = np.zeros(m, dtype=np.int64)
    for x in D11_set:
        D11_ind[x] = 1
    for x in D12_set:
        D12_ind[x] = 1

    # Make D11 symmetric
    D11_sym = set()
    for x in D11_set:
        if x != 0:
            D11_sym.add(x)
            D11_sym.add((-x) % m)
    D11_ind = np.zeros(m, dtype=np.int64)
    for x in D11_sym:
        D11_ind[x] = 1

    delta_11 = compute_delta_fft(D11_ind, m)
    delta_12 = compute_delta_fft(D12_ind, m)
    D12T_ind = np.zeros(m, dtype=np.int64)
    for x in range(m):
        if D12_ind[x]:
            D12T_ind[(-x) % m] = 1
    delta_12T = compute_delta_fft(D12T_ind, m)

    cost = compute_cost_vectorized(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n)
    d11s = int(np.sum(D11_ind))
    d12s = int(np.sum(D12_ind))
    print(f"  {label}: |D11|={d11s}, |D12|={d12s}, cost={cost}")
    return cost


def main():
    classes = classify()
    print("CRT classes for Z_63 = Z_9 x Z_7:")
    for label in sorted(classes.keys()):
        elems = sorted(classes[label])
        print(f"  {label}: {len(elems)} elements = {elems}")

    print(f"\nElement 0 -> class: {crt_decode(0)}")

    # Verify CRT round-trip
    for x in range(63):
        a, b = crt_decode(x)
        assert crt_encode(a, b) == x

    # Symmetry analysis: x <-> -x mod 63
    print("\nSymmetry classes under negation:")
    for label in sorted(classes.keys()):
        neg_elements = [(-x) % M for x in classes[label]]
        neg_labels = set()
        for nx in neg_elements:
            for l2, e2 in classes.items():
                if nx in e2:
                    neg_labels.add(l2)
        print(f"  -{label} = {neg_labels}")

    # Try various constructions
    print("\n=== Construction Attempts ===\n")

    # Construction 1: "Paley-like" - put QR*QR and mixed in D11
    # For Paley on GF(q), D11 = QR. Here try combining QR structures.

    # The symmetric constraint means for d in D11, also -d in D11.
    # Need |D11| = 32 (or 30).

    # Strategy: D11 = union of classes that are closed under negation
    # Class sizes: RR=9, RN=9, NR=9, NN=9, R0=3, N0=3, 0R=3, 0N=3, ZR=6, ZN=6, Z0=2, 00=1
    # Negation: -RR: for (a,b) -> (-a mod 9, -b mod 7)
    # -1 mod 9 = 8 which is QNR. -1 mod 7 = 6 which is QNR.
    # So -(QR,QR) = (QNR*QR, QNR*QR) -- no wait.
    # -(a,b) = (-a, -b). If a in QR_9, is -a in QR_9?
    # -1 mod 9 = 8. 8 in QNR_9. So -QR_9 = {-1,-4,-7} = {8,5,2} = QNR_9.
    # -1 mod 7 = 6. 6 in QNR_7. So -QR_7 = {-1,-2,-4} = {6,5,3} = QNR_7.
    # Therefore: -RR = NN, -RN = NR, -R0 = N0, -0R = 0N, -ZR = ZN, -Z0 = Z0.

    print("Negation mapping of classes:")
    print("  -RR = NN (9+9=18)")
    print("  -RN = NR (9+9=18)")
    print("  -R0 = N0 (3+3=6)")
    print("  -0R = 0N (3+3=6)")
    print("  -ZR = ZN (6+6=12)")
    print("  -Z0 = Z0 (2, self-paired)")
    print()

    # So symmetric D11 must be union of:
    # - RR ∪ NN (size 18)
    # - RN ∪ NR (size 18)
    # - R0 ∪ N0 (size 6)
    # - 0R ∪ 0N (size 6)
    # - ZR ∪ ZN (size 12)
    # - Z0 (size 2)
    # Total: 18+18+6+6+12+2 = 62 (all nonzero elements)

    # Need |D11| = 32. Possible combinations:
    # RR∪NN (18) + RN∪NR (18) = 36 -> too big
    # RR∪NN (18) + R0∪N0 (6) + 0R∪0N (6) + Z0 (2) = 32! ✓
    # RR∪NN (18) + ZR∪ZN (12) + Z0 (2) = 32! ✓
    # RN∪NR (18) + R0∪N0 (6) + 0R∪0N (6) + Z0 (2) = 32! ✓
    # RN∪NR (18) + ZR∪ZN (12) + Z0 (2) = 32! ✓
    # Many other combinations...

    combos_32 = [
        ('RR∪NN + R0∪N0 + 0R∪0N + Z0', ['RR', 'NN', 'R0', 'N0', '0R', '0N', 'Z0']),
        ('RR∪NN + ZR∪ZN + Z0', ['RR', 'NN', 'ZR', 'ZN', 'Z0']),
        ('RN∪NR + R0∪N0 + 0R∪0N + Z0', ['RN', 'NR', 'R0', 'N0', '0R', '0N', 'Z0']),
        ('RN∪NR + ZR∪ZN + Z0', ['RN', 'NR', 'ZR', 'ZN', 'Z0']),
        ('RR∪NN + R0∪N0 + ZR∪ZN', ['RR', 'NN', 'R0', 'N0', 'ZR', 'ZN']),  # 18+6+12=36 too big
        ('RR∪NN + 0R∪0N + ZR∪ZN', ['RR', 'NN', '0R', '0N', 'ZR', 'ZN']),  # 18+6+12=36 too big
        ('RR∪NN + R0∪N0 + ZR∪ZN + Z0 - some', None),  # Can't easily subtract
    ]

    # Also for |D11|=30:
    combos_30 = [
        ('RR∪NN + ZR∪ZN', ['RR', 'NN', 'ZR', 'ZN']),  # 18+12=30
        ('RN∪NR + ZR∪ZN', ['RN', 'NR', 'ZR', 'ZN']),  # 18+12=30
        ('RR∪NN + R0∪N0 + 0R∪0N', ['RR', 'NN', 'R0', 'N0', '0R', '0N']),  # 18+6+6=30
        ('RN∪NR + R0∪N0 + 0R∪0N', ['RN', 'NR', 'R0', 'N0', '0R', '0N']),  # 18+6+6=30
    ]

    print("=== |D11| = 32 constructions ===\n")
    for desc, cls_list in combos_32:
        if cls_list is None:
            continue
        D11_set = set()
        for cls in cls_list:
            D11_set.update(classes[cls])
        D11_set.discard(0)
        if len(D11_set) != 32:
            print(f"  {desc}: |D11|={len(D11_set)} (skipping, want 32)")
            continue

        # Try D12 = {0} + various subsets
        # For D12, try QR-like sets and random
        for d12_desc, D12_set in [
            ("D12=class-based(RR+RN+R0+0R+ZR+half_Z0)", None),  # placeholder
        ]:
            pass

        # For now, just test D11 with random D12 (SA will optimize D12)
        random.seed(42)
        D12_set = {0}
        remaining = list(set(range(1, M)) - {0})
        random.shuffle(remaining)
        D12_set.update(remaining[:30])

        try_construction(D11_set, D12_set, desc)

    print("\n=== |D11| = 30 constructions ===\n")
    for desc, cls_list in combos_30:
        D11_set = set()
        for cls in cls_list:
            D11_set.update(classes[cls])
        D11_set.discard(0)
        if len(D11_set) != 30:
            print(f"  {desc}: |D11|={len(D11_set)} (skipping, want 30)")
            continue

        random.seed(42)
        D12_set = {0}
        remaining = list(set(range(1, M)) - {0})
        random.shuffle(remaining)
        D12_set.update(remaining[:30])

        try_construction(D11_set, D12_set, desc)

    # Now try SA-refining the best algebraic initializations
    print("\n=== SA refinement from algebraic starts ===\n")
    from n32_fast_numpy import compute_cost_vectorized as cost_fn

    best_overall = 999
    for desc, cls_list in combos_32 + combos_30:
        if cls_list is None:
            continue
        D11_set = set()
        for cls in cls_list:
            D11_set.update(classes[cls])
        D11_set.discard(0)
        d11_size = len(D11_set)
        if d11_size not in [30, 32]:
            continue

        # Run short SA from this D11 with random D12
        m = M
        n = N_PARAM
        pairs = []
        for x in range(1, m):
            neg_x = (-x) % m
            if x <= neg_x:
                pairs.append((x, neg_x))
        num_pairs = len(pairs)

        for seed in range(5):
            random.seed(seed)
            D11_ind = np.zeros(m, dtype=np.int64)
            for x in D11_set:
                D11_ind[x] = 1

            D12_ind = np.zeros(m, dtype=np.int64)
            D12_ind[0] = 1
            rest = random.sample(range(1, m), 30)
            for x in rest:
                D12_ind[x] = 1

            D12T_ind = np.zeros(m, dtype=np.int64)
            for x in range(m):
                if D12_ind[x]:
                    D12T_ind[(-x) % m] = 1

            delta_11 = compute_delta_fft(D11_ind, m)
            delta_12 = compute_delta_fft(D12_ind, m)
            delta_12T = compute_delta_fft(D12T_ind, m)

            current_cost = cost_fn(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n)
            best_cost = current_cost

            T = 8.0
            max_iter = 3000000
            alpha = 1 - 6.0 / max_iter

            for it in range(max_iter):
                r = random.random()
                if r < 0.15:
                    # D11 pair swap
                    in_pairs = [i for i in range(num_pairs) if D11_ind[pairs[i][0]]]
                    out_pairs = [i for i in range(num_pairs) if not D11_ind[pairs[i][0]]]
                    if not in_pairs or not out_pairs:
                        continue
                    rp = random.choice(in_pairs)
                    ap = random.choice(out_pairs)
                    old = delta_11.copy()
                    D11_ind[pairs[rp][0]] = 0
                    D11_ind[pairs[rp][1]] = 0
                    D11_ind[pairs[ap][0]] = 1
                    D11_ind[pairs[ap][1]] = 1
                    delta_11 = compute_delta_fft(D11_ind, m)
                    new_cost = cost_fn(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n)
                    dc = new_cost - current_cost
                    if dc <= 0 or random.random() < math.exp(-dc / max(T, 0.0001)):
                        current_cost = new_cost
                    else:
                        D11_ind[pairs[rp][0]] = 1
                        D11_ind[pairs[rp][1]] = 1
                        D11_ind[pairs[ap][0]] = 0
                        D11_ind[pairs[ap][1]] = 0
                        delta_11 = old
                else:
                    d12_list = np.where(D12_ind[1:] == 1)[0] + 1
                    not_in = np.where(D12_ind[1:] == 0)[0] + 1
                    if len(d12_list) == 0 or len(not_in) == 0:
                        continue
                    rem = d12_list[random.randint(0, len(d12_list)-1)]
                    add_el = not_in[random.randint(0, len(not_in)-1)]
                    old_d12 = delta_12.copy()
                    old_d12T = delta_12T.copy()
                    D12_ind[rem] = 0
                    D12_ind[add_el] = 1
                    D12T_ind[(-rem) % m] = 0
                    D12T_ind[(-add_el) % m] = 1
                    delta_12 = compute_delta_fft(D12_ind, m)
                    delta_12T = compute_delta_fft(D12T_ind, m)
                    new_cost = cost_fn(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n)
                    dc = new_cost - current_cost
                    if dc <= 0 or random.random() < math.exp(-dc / max(T, 0.0001)):
                        current_cost = new_cost
                    else:
                        D12_ind[rem] = 1
                        D12_ind[add_el] = 0
                        D12T_ind[(-rem) % m] = 1
                        D12T_ind[(-add_el) % m] = 0
                        delta_12 = old_d12
                        delta_12T = old_d12T

                if current_cost < best_cost:
                    best_cost = current_cost
                    if best_cost == 0:
                        break

                T = max(T * alpha, 0.0001)

            print(f"  {desc} (|D11|={d11_size}), seed={seed}: cost={best_cost}")
            if best_cost < best_overall:
                best_overall = best_cost
            if best_cost == 0:
                print("SOLUTION FOUND!")
                return

    print(f"\nBest overall from algebraic starts: {best_overall}")


if __name__ == "__main__":
    main()
