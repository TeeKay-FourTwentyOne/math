"""
Refine near-solutions (cost=4) by exhaustively trying all D12 modifications.

Given a D11 that achieves cost=4, try ALL possible D12 sets to find one with cost=0.

For m=63, |D12|=31, there are C(62,30) ≈ 10^17 choices. Too many for brute force.
But we can enumerate all SINGLE swaps (30 * 31 = 930 options) and all DOUBLE swaps
(C(30,2) * C(31,2) = 435 * 465 ≈ 202K options) from a near-optimal D12.

For triple swaps: C(30,3) * C(31,3) ≈ 4060 * 4495 ≈ 18M options - feasible!
"""

import sys, os, math, random, time
import numpy as np
from itertools import combinations
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ramsey_core import BlockCirculantGraph, verify_construction, save_construction


def compute_delta_fft(indicator, m):
    a = np.array(indicator, dtype=np.float64)
    fa = np.fft.rfft(a)
    return np.round(np.fft.irfft(np.conj(fa) * fa, n=m)).astype(np.int64)


def compute_cost(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n):
    d11_size = int(np.sum(D11_ind))
    d12_size = int(np.sum(D12_ind))
    d1 = d11_size + d12_size
    d22_size = m - 1 - d11_size
    d2 = d22_size + d12_size
    N = 2 * m
    red_thresh = n - 2
    blue_thresh = n - 1

    common_v11 = delta_11[1:] + delta_12[1:]
    d11_mask = D11_ind[1:].astype(bool)
    cost = int(np.sum(np.maximum(common_v11[d11_mask] - red_thresh, 0)))
    blue_common = (N - 2) - 2 * d1 + common_v11[~d11_mask]
    cost += int(np.sum(np.maximum(blue_common - blue_thresh, 0)))

    v22_const = m - 2 - 2 * d11_size
    d22_delta = delta_11[1:] + v22_const + 2 * D11_ind[1:]
    common_v22 = d22_delta + delta_12T[1:]
    cost += int(np.sum(np.maximum(common_v22[~d11_mask] - red_thresh, 0)))
    blue_common_v22 = (N - 2) - 2 * d2 + common_v22[d11_mask]
    cost += int(np.sum(np.maximum(blue_common_v22 - blue_thresh, 0)))
    return cost


def find_near_solution(n, d11_size, max_iter=100000000):
    """Use SA to find a near-solution (cost <= 4) and return the D11, D12."""
    m = 2 * n - 1
    d12_size = n - 1

    pairs = []
    for x in range(1, m):
        neg_x = (-x) % m
        if x <= neg_x:
            pairs.append((x, neg_x))
    num_pairs = len(pairs)

    best_overall = 999
    best_D11 = None
    best_D12 = None

    for seed in range(500):
        random.seed(seed)
        selected = random.sample(range(num_pairs), d11_size // 2)
        D11_ind = np.zeros(m, dtype=np.int64)
        for i in selected:
            D11_ind[pairs[i][0]] = 1
            D11_ind[pairs[i][1]] = 1

        D12_ind = np.zeros(m, dtype=np.int64)
        D12_ind[0] = 1
        rest = random.sample(range(1, m), d12_size - 1)
        for x in rest:
            D12_ind[x] = 1

        D12T_ind = np.zeros(m, dtype=np.int64)
        for x in range(m):
            if D12_ind[x]:
                D12T_ind[(-x) % m] = 1

        delta_11 = compute_delta_fft(D11_ind, m)
        delta_12 = compute_delta_fft(D12_ind, m)
        delta_12T = compute_delta_fft(D12T_ind, m)

        current_cost = compute_cost(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n)
        best_cost = current_cost
        b_D11 = np.copy(D11_ind)
        b_D12 = np.copy(D12_ind)

        T = 15.0
        T_min = 0.0001
        alpha = math.exp(math.log(T_min / T) / max_iter)

        for it in range(max_iter):
            r = random.random()
            if r < 0.22:
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
                new_cost = compute_cost(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n)
                dc = new_cost - current_cost
                if dc <= 0 or random.random() < math.exp(-dc / T):
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
                new_cost = compute_cost(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n)
                dc = new_cost - current_cost
                if dc <= 0 or random.random() < math.exp(-dc / T):
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
                b_D11 = np.copy(D11_ind)
                b_D12 = np.copy(D12_ind)
                if best_cost == 0:
                    return b_D11, b_D12, 0

            T *= alpha

        print(f"  seed={seed}: best={best_cost}", flush=True)

        if best_cost < best_overall:
            best_overall = best_cost
            best_D11 = np.copy(b_D11)
            best_D12 = np.copy(b_D12)
            if best_cost <= 4:
                # Found a near-solution, try to refine it
                result = refine_d12(n, b_D11, b_D12)
                if result is not None:
                    return result[0], result[1], 0

    return best_D11, best_D12, best_overall


def refine_d12(n, D11_ind, D12_start, max_swaps=3):
    """Given a fixed D11, try to find a D12 with cost=0 by local search around D12_start."""
    m = 2 * n - 1
    delta_11 = compute_delta_fft(D11_ind, m)

    d12_in = [x for x in range(1, m) if D12_start[x]]
    d12_out = [x for x in range(1, m) if not D12_start[x]]

    total_tried = 0

    # Try all single swaps
    print(f"  Trying single swaps ({len(d12_in)} * {len(d12_out)} = {len(d12_in)*len(d12_out)})...", flush=True)
    for rem in d12_in:
        for add_el in d12_out:
            D12_ind = np.copy(D12_start)
            D12_ind[rem] = 0
            D12_ind[add_el] = 1
            D12T_ind = np.zeros(m, dtype=np.int64)
            for x in range(m):
                if D12_ind[x]:
                    D12T_ind[(-x) % m] = 1
            delta_12 = compute_delta_fft(D12_ind, m)
            delta_12T = compute_delta_fft(D12T_ind, m)
            c = compute_cost(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n)
            total_tried += 1
            if c == 0:
                print(f"  FOUND at single swap: remove {rem}, add {add_el}")
                return D11_ind, D12_ind
    print(f"  Single swaps: {total_tried} tried, none found", flush=True)

    if max_swaps < 2:
        return None

    # Try all double swaps
    n_double = len(list(combinations(range(len(d12_in)), 2))) * len(list(combinations(range(len(d12_out)), 2)))
    print(f"  Trying double swaps (~{n_double})...", flush=True)
    for r1, r2 in combinations(d12_in, 2):
        for a1, a2 in combinations(d12_out, 2):
            D12_ind = np.copy(D12_start)
            D12_ind[r1] = 0
            D12_ind[r2] = 0
            D12_ind[a1] = 1
            D12_ind[a2] = 1
            D12T_ind = np.zeros(m, dtype=np.int64)
            for x in range(m):
                if D12_ind[x]:
                    D12T_ind[(-x) % m] = 1
            delta_12 = compute_delta_fft(D12_ind, m)
            delta_12T = compute_delta_fft(D12T_ind, m)
            c = compute_cost(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n)
            total_tried += 1
            if c == 0:
                print(f"  FOUND at double swap: remove {r1},{r2}, add {a1},{a2}")
                return D11_ind, D12_ind
            if total_tried % 100000 == 0:
                print(f"    {total_tried} tried...", flush=True)
    print(f"  Double swaps done: {total_tried} total tried", flush=True)

    if max_swaps < 3:
        return None

    # Try all triple swaps
    n_triple = len(list(combinations(range(len(d12_in)), 3))) * len(list(combinations(range(len(d12_out)), 3)))
    print(f"  Trying triple swaps (~{n_triple})...", flush=True)
    for r1, r2, r3 in combinations(d12_in, 3):
        for a1, a2, a3 in combinations(d12_out, 3):
            D12_ind = np.copy(D12_start)
            D12_ind[r1] = 0
            D12_ind[r2] = 0
            D12_ind[r3] = 0
            D12_ind[a1] = 1
            D12_ind[a2] = 1
            D12_ind[a3] = 1
            D12T_ind = np.zeros(m, dtype=np.int64)
            for x in range(m):
                if D12_ind[x]:
                    D12T_ind[(-x) % m] = 1
            delta_12 = compute_delta_fft(D12_ind, m)
            delta_12T = compute_delta_fft(D12T_ind, m)
            c = compute_cost(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n)
            total_tried += 1
            if c == 0:
                print(f"  FOUND at triple swap: remove {r1},{r2},{r3}, add {a1},{a2},{a3}")
                return D11_ind, D12_ind
            if total_tried % 500000 == 0:
                print(f"    {total_tried} tried...", flush=True)
    print(f"  Triple swaps done: {total_tried} total tried", flush=True)
    return None


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 32
    d11_size = int(sys.argv[2]) if len(sys.argv) > 2 else 32
    max_sa_iter = int(sys.argv[3]) if len(sys.argv) > 3 else 50000000

    m = 2 * n - 1
    print(f"n={n}, m={m}")
    print(f"|D11|={d11_size}, |D12|={n-1}")
    print(f"SA iterations per seed: {max_sa_iter}")
    print(flush=True)

    D11, D12, cost = find_near_solution(n, d11_size, max_sa_iter)

    if cost == 0:
        D11_list = sorted(int(i) for i in range(m) if D11[i])
        D12_list = sorted(int(i) for i in range(m) if D12[i])
        print(f"\n*** SOLUTION FOUND! ***")
        print(f"D11 = {D11_list}")
        print(f"D12 = {D12_list}")

        D11_s = set(D11_list)
        D12_s = set(D12_list)
        D22_s = set(range(1, m)) - D11_s
        G = BlockCirculantGraph(n=n, D11=D11_s, D12=D12_s, D22=D22_s)
        result = verify_construction(G)
        print(f"Full verification: valid={result.valid}")
        if result.valid:
            fn = f"/Users/stephenpadgett/Projects/math/ramsey-book-graphs/solution_n{n}.json"
            save_construction(G, result, fn, solver="refine-SA")
            print(f"SAVED to {fn}")
    else:
        print(f"\nBest cost = {cost}")
