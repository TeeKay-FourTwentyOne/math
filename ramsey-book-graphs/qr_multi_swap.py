#!/usr/bin/env python3
"""Multi-swap QR Perturbation for R(B_{n-1}, B_n) = 4n-1

Single-swap D12 = (QR ∪ {0}) \ {r} fails for p >= 19.
Extend to k-swap: D12 = (QR ∪ {0} ∪ QNR_add) \ QR_remove
with |QR_remove| = k, |QNR_add| = k-1 (maintains |D12| = n-1).

Search over k = 1, 2, 3 to find the minimum swap count needed.
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
    """FFT-based autocorrelation."""
    indicator = np.zeros(m)
    for x in S:
        indicator[x % m] = 1
    F = np.fft.fft(indicator)
    power = np.abs(F) ** 2
    result = np.real(np.fft.ifft(power))
    return np.rint(result).astype(int)


def find_all_symmetric_d11(p, k):
    """Enumerate all symmetric D11 of size k in {1,...,p-1}."""
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


def load_registry_solutions():
    """Load solutions for p ≡ 3 mod 4 primes from registry."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(script_dir, 'solutions_registry.json')
    try:
        with open(path, 'r') as f:
            data = json.load(f)
    except FileNotFoundError:
        return {}
    solutions = {}
    for sol in data['solutions']:
        m = sol['m']
        if is_prime(m) and m % 4 == 3:
            solutions[m] = {
                'n': sol['n'],
                'D11': sol['D11'],
                'D12': sol['D12'],
                'solver': sol['solver'],
            }
    return solutions


def precompute_b_matrix(p, k):
    """Precompute B(d) for all k-swap D12 candidates.

    D12 = (QR ∪ {0} ∪ QNR_add) \ QR_remove
    |QR_remove| = k, |QNR_add| = k-1

    Returns (B_matrix[N, p], candidate_info) where
    candidate_info[i] = (qr_remove_tuple, qnr_add_tuple).
    """
    qr = quadratic_residues(p)
    qr_set = set(qr)
    qnr = sorted(set(range(1, p)) - qr_set)

    candidates = []
    B_list = []

    base_set = qr_set | {0}

    for qr_remove in combinations(qr, k):
        qr_rem_set = set(qr_remove)
        partial = base_set - qr_rem_set
        if k == 1:
            # No QNR to add
            D12 = sorted(partial)
            B = autocorrelation(D12, p)
            B_list.append(B)
            candidates.append((qr_remove, ()))
        else:
            for qnr_add in combinations(qnr, k - 1):
                D12 = sorted(partial | set(qnr_add))
                B = autocorrelation(D12, p)
                B_list.append(B)
                candidates.append((qr_remove, qnr_add))

    if not B_list:
        return np.zeros((0, p), dtype=np.int32), []
    return np.array(B_list, dtype=np.int32), candidates


def check_d11_against_b_matrix(D11, B_matrix, p, n):
    """Check which D12 candidates (rows of B_matrix) satisfy all constraints.

    Returns boolean array of shape (N,).
    """
    D11_set = set(D11)
    k_d11 = len(D11_set)

    A = autocorrelation(list(D11_set), p)

    d11_thresh = n - 2
    d22_thresh = 2 * k_d11 - n + 1

    # Build threshold array for d=0,...,p-1
    threshold = np.zeros(p, dtype=np.int32)
    for d in range(1, p):
        threshold[d] = d11_thresh if d in D11_set else d22_thresh

    # excess[i, d] = A[d] + B_matrix[i, d] - threshold[d]
    # Valid if all excess[:, 1:] <= 0
    sum_ab = B_matrix[:, 1:] + A[np.newaxis, 1:]
    thresh_row = threshold[np.newaxis, 1:]
    valid = np.all(sum_ab <= thresh_row, axis=1)

    return valid, A


def sweep_k_swap(D11, p, n, k):
    """Try all k-swap D12 candidates for a given D11.

    Returns (num_valid, total_candidates, best_candidate_info).
    """
    B_matrix, candidates = precompute_b_matrix(p, k)
    if len(candidates) == 0:
        return 0, 0, None

    valid_mask, A = check_d11_against_b_matrix(D11, B_matrix, p, n)
    num_valid = int(np.sum(valid_mask))

    best_info = None
    if num_valid > 0:
        idx = np.argmax(valid_mask)
        best_info = candidates[idx]

    return num_valid, len(candidates), best_info


def count_candidates(p, k):
    """Count k-swap D12 candidates without enumerating."""
    nqr = (p - 1) // 2
    if k == 1:
        return nqr
    from math import comb
    return comb(nqr, k) * comb(nqr, k - 1)


def main():
    print("=" * 80)
    print("MULTI-SWAP QR PERTURBATION")
    print("D12 = (QR ∪ {0} ∪ QNR_add) \\ QR_remove")
    print("=" * 80)

    primes = [7, 11, 19, 23, 31, 43, 47, 59]
    registry = load_registry_solutions()

    # ========== Phase 1: Per-prime sweep with best D11 ==========
    print("\n" + "=" * 80)
    print("PHASE 1: Per-prime sweep (best D11, k=1,2,3)")
    print("=" * 80)

    # Precompute best D11 for small primes (from single-swap search)
    best_d11_cache = {}

    for p in primes:
        n = (p + 1) // 2
        t0 = time.time()

        print(f"\n  p={p} (n={n}):")
        print(f"    Candidates: k=1:{count_candidates(p,1)}, "
              f"k=2:{count_candidates(p,2)}, k=3:{count_candidates(p,3)}")

        if p in registry:
            D11 = registry[p]['D11']
            source = f"registry ({registry[p]['solver']})"
        else:
            # Find best D11: try all symmetric D11 of size n with k=2
            all_d11 = find_all_symmetric_d11(p, n)
            # Precompute B_matrix for k=2 (reuse for all D11)
            B_mat_k2, cands_k2 = precompute_b_matrix(p, 2)

            best = None
            best_count = -1
            for candidate in all_d11:
                valid_mask, _ = check_d11_against_b_matrix(
                    candidate, B_mat_k2, p, n)
                cnt = int(np.sum(valid_mask))
                if cnt > best_count:
                    best_count = cnt
                    best = candidate
                    if cnt == len(cands_k2):
                        break

            D11 = best
            source = f"best of {len(all_d11)} D11s (by k=2 success)"

        best_d11_cache[p] = D11
        k_d11 = len(D11)

        print(f"    D11 ({source}): {D11}")
        print(f"    |D11|={k_d11}, D11_thresh={n-2}, "
              f"D22_thresh={2*k_d11-n+1}")

        for k in [1, 2, 3]:
            tk = time.time()
            num_valid, total, best_info = sweep_k_swap(D11, p, n, k)
            elapsed_k = time.time() - tk
            pct = num_valid / total * 100 if total > 0 else 0

            status = "***" if num_valid > 0 else "   "
            print(f"    {status} k={k}: {num_valid:6d}/{total:<6d} "
                  f"({pct:5.1f}%) [{elapsed_k:.2f}s]", end="")

            if best_info and num_valid > 0:
                qr_rem, qnr_add = best_info
                print(f"  e.g. remove {qr_rem}, add {qnr_add}", end="")
            print()

        elapsed = time.time() - t0
        print(f"    Total: {elapsed:.1f}s")

    # ========== Phase 2: Exhaustive D11 search with k=2 ==========
    print("\n" + "=" * 80)
    print("PHASE 2: Exhaustive D11 search with k=2 (p=23, p=31)")
    print("=" * 80)

    for p in [23, 31]:
        n = (p + 1) // 2
        all_d11 = find_all_symmetric_d11(p, n)
        num_d11 = len(all_d11)

        # Precompute B_matrix for k=2 once
        t0 = time.time()
        B_matrix, candidates = precompute_b_matrix(p, 2)
        num_cands = len(candidates)
        print(f"\n  p={p}: n={n}, |D11|={num_d11}, "
              f"|D12 candidates (k=2)|={num_cands}")
        print(f"  Total pairs: {num_d11 * num_cands}")

        d11_with_valid = 0
        total_valid = 0
        best_d11 = None
        best_count = 0
        # Track per-candidate success
        cand_success = np.zeros(num_cands, dtype=np.int32)

        for idx, D11 in enumerate(all_d11):
            valid_mask, _ = check_d11_against_b_matrix(D11, B_matrix, p, n)
            cnt = int(np.sum(valid_mask))
            total_valid += cnt
            if cnt > 0:
                d11_with_valid += 1
                cand_success += valid_mask.astype(np.int32)
            if cnt > best_count:
                best_count = cnt
                best_d11 = D11

            if (idx + 1) % 2000 == 0:
                elapsed = time.time() - t0
                print(f"    [{elapsed:.1f}s] {idx+1}/{num_d11}: "
                      f"{d11_with_valid} D11s with valid pair")

        elapsed = time.time() - t0
        total_pairs = num_d11 * num_cands
        d11_rate = d11_with_valid / num_d11 * 100
        pair_rate = total_valid / total_pairs * 100

        print(f"\n  Results for p={p} ({elapsed:.1f}s):")
        print(f"    D11s with >=1 valid: "
              f"{d11_with_valid}/{num_d11} ({d11_rate:.1f}%)")
        print(f"    Total valid pairs: "
              f"{total_valid}/{total_pairs} ({pair_rate:.1f}%)")
        if best_d11:
            print(f"    Best D11: {best_d11} "
                  f"({best_count}/{num_cands} candidates work)")

        # Show top candidates
        if total_valid > 0:
            top_indices = np.argsort(cand_success)[::-1][:5]
            print(f"    Top D12 candidates (by # D11s they work with):")
            for i in top_indices:
                if cand_success[i] == 0:
                    break
                qr_rem, qnr_add = candidates[i]
                print(f"      remove {qr_rem}, add {qnr_add}: "
                      f"{cand_success[i]}/{num_d11} D11s "
                      f"({cand_success[i]/num_d11*100:.1f}%)")

    # ========== Phase 3: Exhaustive D11 search with k=3 for p=23 ==========
    print("\n" + "=" * 80)
    print("PHASE 3: Exhaustive D11 search with k=3 (p=23)")
    print("=" * 80)

    p = 23
    n = (p + 1) // 2
    all_d11 = find_all_symmetric_d11(p, n)
    num_d11 = len(all_d11)

    t0 = time.time()
    B_matrix, candidates = precompute_b_matrix(p, 3)
    num_cands = len(candidates)
    precomp_time = time.time() - t0
    print(f"\n  p={p}: |D11|={num_d11}, "
          f"|D12 candidates (k=3)|={num_cands} "
          f"[precomputed in {precomp_time:.1f}s]")
    print(f"  Total pairs: {num_d11 * num_cands}")

    d11_with_valid = 0
    total_valid = 0
    best_d11 = None
    best_count = 0
    cand_success = np.zeros(num_cands, dtype=np.int32)

    for idx, D11 in enumerate(all_d11):
        valid_mask, _ = check_d11_against_b_matrix(D11, B_matrix, p, n)
        cnt = int(np.sum(valid_mask))
        total_valid += cnt
        if cnt > 0:
            d11_with_valid += 1
            cand_success += valid_mask.astype(np.int32)
        if cnt > best_count:
            best_count = cnt
            best_d11 = D11

    elapsed = time.time() - t0
    total_pairs = num_d11 * num_cands
    d11_rate = d11_with_valid / num_d11 * 100
    pair_rate = total_valid / total_pairs * 100

    print(f"\n  Results for p={p} (k=3) ({elapsed:.1f}s):")
    print(f"    D11s with >=1 valid: "
          f"{d11_with_valid}/{num_d11} ({d11_rate:.1f}%)")
    print(f"    Total valid pairs: "
          f"{total_valid}/{total_pairs} ({pair_rate:.1f}%)")
    if best_d11:
        print(f"    Best D11: {best_d11} "
              f"({best_count}/{num_cands} candidates)")

    if total_valid > 0:
        top_indices = np.argsort(cand_success)[::-1][:5]
        print(f"    Top D12 candidates:")
        for i in top_indices:
            if cand_success[i] == 0:
                break
            qr_rem, qnr_add = candidates[i]
            print(f"      remove {qr_rem}, add {qnr_add}: "
                  f"{cand_success[i]}/{num_d11} D11s "
                  f"({cand_success[i]/num_d11*100:.1f}%)")

    # ========== Summary ==========
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print("\nThe question: how many QR-to-QNR swaps in D12 are needed")
    print("to make the perturbation approach work?")
    print("\nIf the required k stays bounded as p grows, the approach")
    print("is viable for a general proof. If k grows with p, it's not.")


if __name__ == '__main__':
    main()
