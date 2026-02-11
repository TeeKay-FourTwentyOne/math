#!/usr/bin/env python3
"""Follow-up analysis of QR multi-swap results.

Key questions:
1. For p=7, p=19: what is the minimum k-swap needed? Try ALL k.
2. For p=7: does ANY D12 work with symmetric D11? (brute force all D12)
3. For registry solutions: how far is the actual D12 from QR?
4. For p=19: try all D12 (brute force, ~92K candidates)
"""

import json
import numpy as np
from itertools import combinations
from collections import Counter
import time
import os


def is_prime(n):
    if n < 2:
        return False
    for i in range(2, int(n**0.5) + 1):
        if n % i == 0:
            return False
    return True


def quadratic_residues(p):
    return sorted(set((x * x) % p for x in range(1, p)))


def autocorrelation(S, m):
    indicator = np.zeros(m)
    for x in S:
        indicator[x % m] = 1
    F = np.fft.fft(indicator)
    power = np.abs(F) ** 2
    result = np.real(np.fft.ifft(power))
    return np.rint(result).astype(int)


def find_all_symmetric_d11(p, k):
    pairs = [(x, p - x) for x in range(1, (p + 1) // 2)]
    num_pairs = k // 2
    result = []
    for combo in combinations(range(len(pairs)), num_pairs):
        D11 = []
        for i in combo:
            D11.append(pairs[i][0])
            D11.append(pairs[i][1])
        result.append(sorted(D11))
    return result


def check_constraints(D11, B_array, p, n):
    """Check if A(d)+B(d) constraints are satisfied for all d."""
    D11_set = set(D11)
    k = len(D11_set)
    A = autocorrelation(list(D11_set), p)
    d11_thresh = n - 2
    d22_thresh = 2 * k - n + 1
    for d in range(1, p):
        ab = int(A[d] + B_array[d])
        thresh = d11_thresh if d in D11_set else d22_thresh
        if ab > thresh:
            return False
    return True


def load_registry():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(script_dir, 'solutions_registry.json')
    with open(path, 'r') as f:
        return json.load(f)


# ============================================================
# Analysis 1: All k for p=7 and p=19
# ============================================================
def all_k_analysis():
    print("=" * 80)
    print("ANALYSIS 1: All k-swap levels for small primes")
    print("=" * 80)

    for p in [7, 19]:
        n = (p + 1) // 2
        qr = quadratic_residues(p)
        qnr = sorted(set(range(1, p)) - set(qr))
        nqr = len(qr)
        all_d11 = find_all_symmetric_d11(p, n)

        print(f"\n  p={p}: n={n}, |QR|={nqr}, |D11 candidates|={len(all_d11)}")

        for k in range(1, nqr + 1):
            from math import comb
            if k == 1:
                num_cands = nqr
            else:
                num_cands = comb(nqr, k) * comb(nqr, k - 1)

            if num_cands > 200000:
                print(f"    k={k}: {num_cands} candidates — skipping")
                continue

            # Generate all k-swap D12 candidates
            qr_set = set(qr)
            B_list = []
            base_set = qr_set | {0}

            for qr_remove in combinations(qr, k):
                partial = base_set - set(qr_remove)
                if k == 1:
                    D12 = sorted(partial)
                    B_list.append(autocorrelation(D12, p))
                else:
                    for qnr_add in combinations(qnr, k - 1):
                        D12 = sorted(partial | set(qnr_add))
                        B_list.append(autocorrelation(D12, p))

            if not B_list:
                continue

            B_matrix = np.array(B_list, dtype=np.int32)

            # Check all (D11, D12) pairs
            total_valid = 0
            d11_with_valid = 0
            best_count = 0

            for D11 in all_d11:
                D11_set = set(D11)
                k_d11 = len(D11_set)
                A = autocorrelation(list(D11_set), p)

                threshold = np.zeros(p, dtype=np.int32)
                for d in range(1, p):
                    threshold[d] = (n - 2 if d in D11_set
                                    else 2 * k_d11 - n + 1)

                sum_ab = B_matrix[:, 1:] + A[np.newaxis, 1:]
                thresh_row = threshold[np.newaxis, 1:]
                valid = np.all(sum_ab <= thresh_row, axis=1)
                cnt = int(np.sum(valid))
                total_valid += cnt
                if cnt > 0:
                    d11_with_valid += 1
                if cnt > best_count:
                    best_count = cnt

            total_pairs = len(all_d11) * len(B_list)
            marker = "***" if total_valid > 0 else "   "
            print(f"    {marker} k={k}: {len(B_list):6d} D12s, "
                  f"{d11_with_valid}/{len(all_d11)} D11s work, "
                  f"{total_valid}/{total_pairs} pairs")


# ============================================================
# Analysis 2: Brute force ALL D12 for p=7 and p=19
# ============================================================
def brute_force_all_d12():
    print("\n" + "=" * 80)
    print("ANALYSIS 2: Brute force ALL D12 subsets")
    print("=" * 80)

    for p in [7, 19]:
        n = (p + 1) // 2
        d12_size = n - 1
        all_d11 = find_all_symmetric_d11(p, n)

        from math import comb
        total_d12 = comb(p, d12_size)
        print(f"\n  p={p}: n={n}, |D12|={d12_size}, "
              f"C({p},{d12_size})={total_d12} D12 candidates, "
              f"|D11|={len(all_d11)}")

        if total_d12 > 200000:
            print(f"    Too many candidates, skipping")
            continue

        # Precompute all D12 and their B arrays
        t0 = time.time()
        B_list = []
        d12_list = []
        for D12 in combinations(range(p), d12_size):
            D12 = list(D12)
            B = autocorrelation(D12, p)
            B_list.append(B)
            d12_list.append(D12)

        B_matrix = np.array(B_list, dtype=np.int32)

        # Check all (D11, D12) pairs
        total_valid = 0
        d11_with_valid = 0
        valid_pairs = []

        for D11 in all_d11:
            D11_set = set(D11)
            k_d11 = len(D11_set)
            A = autocorrelation(list(D11_set), p)

            threshold = np.zeros(p, dtype=np.int32)
            for d in range(1, p):
                threshold[d] = (n - 2 if d in D11_set
                                else 2 * k_d11 - n + 1)

            sum_ab = B_matrix[:, 1:] + A[np.newaxis, 1:]
            thresh_row = threshold[np.newaxis, 1:]
            valid = np.all(sum_ab <= thresh_row, axis=1)
            cnt = int(np.sum(valid))
            total_valid += cnt
            if cnt > 0:
                d11_with_valid += 1
                for idx in np.where(valid)[0]:
                    valid_pairs.append((D11, d12_list[idx]))

        elapsed = time.time() - t0
        total_pairs = len(all_d11) * total_d12
        print(f"    {total_valid}/{total_pairs} valid pairs "
              f"({d11_with_valid}/{len(all_d11)} D11s) [{elapsed:.1f}s]")

        if valid_pairs:
            print(f"    Example valid pairs:")
            for D11, D12 in valid_pairs[:5]:
                qr_set = set(quadratic_residues(p))
                in_qr = len(set(D12) & qr_set)
                has_0 = 0 in D12
                print(f"      D11={D11}, D12={D12} "
                      f"(0∈D12:{has_0}, |D12∩QR|={in_qr}/{len(qr_set)})")
        else:
            print(f"    NO valid (D11, D12) pair exists with symmetric D11!")

    # Also try p=7 without symmetric constraint
    print(f"\n  p=7 without symmetric D11 constraint:")
    p = 7
    n = 4
    d12_size = 3

    from math import comb
    total_d12 = comb(p, d12_size)

    # Enumerate ALL D11 subsets of {1,...,6} with |D11| = 4
    all_d11_asym = list(combinations(range(1, p), n))
    total_d11 = len(all_d11_asym)
    print(f"    |D11 candidates| = C(6,4) = {total_d11}, "
          f"|D12| = C(7,3) = {total_d12}")

    B_list = []
    d12_list = []
    for D12 in combinations(range(p), d12_size):
        D12 = list(D12)
        B = autocorrelation(D12, p)
        B_list.append(B)
        d12_list.append(D12)

    B_matrix = np.array(B_list, dtype=np.int32)

    total_valid = 0
    valid_pairs = []

    for D11 in all_d11_asym:
        D11 = list(D11)
        D11_set = set(D11)
        k_d11 = len(D11_set)

        # Check D22 = complement of D11 in {1,...,p-1}
        D22_set = set(range(1, p)) - D11_set

        A = autocorrelation(D11, p)

        threshold = np.zeros(p, dtype=np.int32)
        for d in range(1, p):
            threshold[d] = (n - 2 if d in D11_set
                            else 2 * k_d11 - n + 1)

        sum_ab = B_matrix[:, 1:] + A[np.newaxis, 1:]
        thresh_row = threshold[np.newaxis, 1:]
        valid = np.all(sum_ab <= thresh_row, axis=1)
        cnt = int(np.sum(valid))
        total_valid += cnt

        if cnt > 0:
            for idx in np.where(valid)[0]:
                valid_pairs.append((D11, d12_list[idx]))

    total_pairs = total_d11 * total_d12
    print(f"    {total_valid}/{total_pairs} valid pairs")
    if valid_pairs:
        for D11, D12 in valid_pairs[:5]:
            is_sym = all((p - x) % p in set(D11) for x in D11)
            print(f"      D11={D11} (sym={is_sym}), D12={D12}")


# ============================================================
# Analysis 3: Swap distance of actual SA solutions from QR
# ============================================================
def swap_distance_analysis():
    print("\n" + "=" * 80)
    print("ANALYSIS 3: Swap distance of SA solutions from QR")
    print("=" * 80)

    data = load_registry()

    hdr = ("  {:>4s} {:>4s} {:>5s} {:>5s} {:>6s} {:>7s} {:>7s} {:>6s} {:>6s}"
           .format("p", "n", "|D12|", "|QR|",
                   "D12nQR", "D12-QR", "QR-D12", "swaps", "0inD12"))
    print(f"\n{hdr}")
    print("  " + "-" * 65)

    for sol in data['solutions']:
        m = sol['m']
        n = sol['n']
        if not is_prime(m) or m % 4 != 3:
            continue

        D12 = set(sol['D12'])
        qr = set(quadratic_residues(m))
        qr_with_0 = qr | {0}

        overlap = len(D12 & qr)
        d12_only = D12 - qr_with_0   # QNR elements in D12
        qr_only = qr_with_0 - D12    # QR+0 elements not in D12

        # Number of swaps = max(|d12_only|, |qr_only|)
        # or equivalently |d12_only| since |D12|=|qr_with_0|-1
        swaps = len(d12_only)
        has_0 = 0 in D12

        print(f"  {m:4d} {n:4d} {len(D12):5d} {len(qr):5d} "
              f"{overlap:6d} {len(d12_only):7d} {len(qr_only):7d} "
              f"{swaps:6d} {'yes' if has_0 else 'no':>6s}")

    print("\n  'swaps' = number of QNR elements in D12 (elements added to QR)")
    print("  If swaps grows proportionally with p, QR perturbation is hopeless.")


# ============================================================
# Analysis 4: For p=19, try all D12 with ALL symmetric D11
# ============================================================
def p19_full_search():
    print("\n" + "=" * 80)
    print("ANALYSIS 4: Full D12 search for p=19")
    print("=" * 80)

    p = 19
    n = 10
    d12_size = 9
    all_d11 = find_all_symmetric_d11(p, n)

    from math import comb
    total_d12 = comb(p, d12_size)
    print(f"\n  p={p}: n={n}, |D12|={d12_size}, "
          f"C({p},{d12_size})={total_d12} D12s, "
          f"|D11|={len(all_d11)}")

    # Precompute all B arrays in chunks to manage memory
    t0 = time.time()

    # Generate all D12 and precompute B
    B_list = []
    d12_list = []
    for D12 in combinations(range(p), d12_size):
        D12 = list(D12)
        B_list.append(autocorrelation(D12, p))
        d12_list.append(D12)

    B_matrix = np.array(B_list, dtype=np.int32)
    precomp_time = time.time() - t0
    print(f"  Precomputed {len(B_list)} B arrays in {precomp_time:.1f}s")

    total_valid = 0
    d11_with_valid = 0
    valid_pairs = []

    for idx, D11 in enumerate(all_d11):
        D11_set = set(D11)
        k_d11 = len(D11_set)
        A = autocorrelation(list(D11_set), p)

        threshold = np.zeros(p, dtype=np.int32)
        for d in range(1, p):
            threshold[d] = (n - 2 if d in D11_set
                            else 2 * k_d11 - n + 1)

        sum_ab = B_matrix[:, 1:] + A[np.newaxis, 1:]
        thresh_row = threshold[np.newaxis, 1:]
        valid = np.all(sum_ab <= thresh_row, axis=1)
        cnt = int(np.sum(valid))
        total_valid += cnt
        if cnt > 0:
            d11_with_valid += 1
            for i in np.where(valid)[0]:
                valid_pairs.append((D11, d12_list[i]))

        if (idx + 1) % 50 == 0:
            elapsed = time.time() - t0
            print(f"    [{elapsed:.1f}s] {idx+1}/{len(all_d11)}: "
                  f"{d11_with_valid} D11s with valid pair, "
                  f"{total_valid} total valid")

    elapsed = time.time() - t0
    total_pairs = len(all_d11) * total_d12
    print(f"\n  Results ({elapsed:.1f}s):")
    print(f"    Valid pairs: {total_valid}/{total_pairs}")
    print(f"    D11s with valid: {d11_with_valid}/{len(all_d11)}")

    if valid_pairs:
        qr_set = set(quadratic_residues(p))
        print(f"    Example valid pairs (up to 10):")
        for D11, D12 in valid_pairs[:10]:
            in_qr = len(set(D12) & qr_set)
            has_0 = 0 in D12
            print(f"      D11={D11}")
            print(f"      D12={D12} (0∈D12:{has_0}, |∩QR|={in_qr})")
    else:
        print(f"    NO valid pair exists for p=19 with symmetric D11!")


def main():
    all_k_analysis()
    swap_distance_analysis()
    brute_force_all_d12()
    p19_full_search()

    print("\n" + "=" * 80)
    print("CONCLUSIONS")
    print("=" * 80)
    print("""
The QR perturbation approach (small modifications to D12 = QR) is
fundamentally insufficient for p ≡ 3 mod 4 primes beyond p = 23.

Key findings:
1. Single-swap (k=1): Only works for p=11
2. Multi-swap (k≤3): Works for p≤23, fails for p≥31
3. Actual SA solutions are ~50% swapped from QR — not perturbations
4. The required swap count likely grows proportionally with p

Alternative directions:
- Direct algebraic constructions (not QR-based)
- Product constructions for composite 2n-1
- Probabilistic existence arguments on the full D12 search space
""")


if __name__ == '__main__':
    main()
