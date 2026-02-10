"""
Analyze near-solutions for n=32 (m=63).

Find a near-solution with SA, then analyze where violations occur.
This helps understand what makes n=32 hard.
"""

import sys, os, math, random, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ramsey_core import BlockCirculantGraph, verify_construction, save_construction
from n32_fast_numpy import compute_delta_fft, compute_cost_vectorized, M, N_PARAM, D12_SIZE


def get_violations(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n):
    """Return list of (type, d, excess) violations."""
    d11_size = int(np.sum(D11_ind))
    d12_size = int(np.sum(D12_ind))
    d22_size = m - 1 - d11_size
    d1 = d11_size + d12_size
    d2 = d22_size + d12_size
    N = 2 * m
    red_thresh = n - 2
    blue_thresh = n - 1
    violations = []

    for d in range(1, m):
        common = int(delta_11[d] + delta_12[d])
        if D11_ind[d]:
            if common > red_thresh:
                violations.append(('V1V1_red', d, common - red_thresh, common))
        else:
            blue_common = (N - 2) - 2 * d1 + common
            if blue_common > blue_thresh:
                violations.append(('V1V1_blue', d, blue_common - blue_thresh, blue_common))

    v22_const = m - 2 - 2 * d11_size
    for d in range(1, m):
        d22_delta = int(delta_11[d]) + v22_const + (2 if D11_ind[d] else 0)
        common = d22_delta + int(delta_12T[d])
        if not D11_ind[d]:
            if common > red_thresh:
                violations.append(('V2V2_red', d, common - red_thresh, common))
        else:
            blue_common = (N - 2) - 2 * d2 + common
            if blue_common > blue_thresh:
                violations.append(('V2V2_blue', d, blue_common - blue_thresh, blue_common))

    return violations


def find_near_solution(d11_size, seed, max_iter=8000000):
    """Run SA and return best near-solution."""
    m = M
    n = N_PARAM
    rng = random.Random(seed)

    pairs = []
    for x in range(1, m):
        neg_x = (-x) % m
        if x <= neg_x:
            pairs.append((x, neg_x))
    num_pairs = len(pairs)

    num_selected = d11_size // 2
    selected = rng.sample(range(num_pairs), num_selected)
    D11_ind = np.zeros(m, dtype=np.int64)
    for i in selected:
        D11_ind[pairs[i][0]] = 1
        D11_ind[pairs[i][1]] = 1

    d12_rest = rng.sample(range(1, m), D12_SIZE - 1)
    D12_ind = np.zeros(m, dtype=np.int64)
    D12_ind[0] = 1
    for x in d12_rest:
        D12_ind[x] = 1

    delta_11 = compute_delta_fft(D11_ind, m)
    delta_12 = compute_delta_fft(D12_ind, m)
    D12T_ind = np.zeros(m, dtype=np.int64)
    for x in range(m):
        if D12_ind[x]:
            D12T_ind[(-x) % m] = 1
    delta_12T = compute_delta_fft(D12T_ind, m)

    current_cost = compute_cost_vectorized(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n)
    best_cost = current_cost
    best_D11 = np.copy(D11_ind)
    best_D12 = np.copy(D12_ind)

    T = 10.0
    alpha = 1 - 6.0 / max_iter
    last_imp = 0

    for it in range(max_iter):
        if random.random() < 0.22:
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
            new_cost = compute_cost_vectorized(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n)
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
            d12_list = list(np.where(D12_ind[1:] == 1)[0] + 1)
            not_in = list(np.where(D12_ind[1:] == 0)[0] + 1)
            if not d12_list or not not_in:
                continue
            rem = random.choice(d12_list)
            add_el = random.choice(not_in)
            old_d12 = delta_12.copy()
            old_d12T = delta_12T.copy()
            D12_ind[rem] = 0
            D12_ind[add_el] = 1
            D12T_ind[(-rem) % m] = 0
            D12T_ind[(-add_el) % m] = 1
            delta_12 = compute_delta_fft(D12_ind, m)
            delta_12T = compute_delta_fft(D12T_ind, m)
            new_cost = compute_cost_vectorized(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n)
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
            best_D11 = np.copy(D11_ind)
            best_D12 = np.copy(D12_ind)
            last_imp = it
            if best_cost == 0:
                break

        if it - last_imp > 500000 and best_cost > 0 and best_cost <= 12:
            T = 4.0
            last_imp = it
            D11_ind = np.copy(best_D11)
            D12_ind = np.copy(best_D12)
            D12T_ind = np.zeros(m, dtype=np.int64)
            for x in range(m):
                if D12_ind[x]:
                    D12T_ind[(-x) % m] = 1
            delta_11 = compute_delta_fft(D11_ind, m)
            delta_12 = compute_delta_fft(D12_ind, m)
            delta_12T = compute_delta_fft(D12T_ind, m)
            current_cost = best_cost

        T = max(T * alpha, 0.0001)

    return best_cost, best_D11, best_D12


def exhaustive_d12_fix(D11_ind, D12_ind, m, n):
    """Given D11 fixed, try to find D12 that gives cost=0 by exhaustive single-swap search."""
    delta_11 = compute_delta_fft(D11_ind, m)
    D12T_ind = np.zeros(m, dtype=np.int64)
    for x in range(m):
        if D12_ind[x]:
            D12T_ind[(-x) % m] = 1
    delta_12 = compute_delta_fft(D12_ind, m)
    delta_12T = compute_delta_fft(D12T_ind, m)

    best_cost = compute_cost_vectorized(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n)

    d12_list = list(np.where(D12_ind[1:] == 1)[0] + 1)
    not_in = list(np.where(D12_ind[1:] == 0)[0] + 1)

    improved = True
    while improved:
        improved = False
        for rem in d12_list:
            for add_el in not_in:
                D12_ind[rem] = 0
                D12_ind[add_el] = 1
                D12T_ind[(-rem) % m] = 0
                D12T_ind[(-add_el) % m] = 1
                delta_12 = compute_delta_fft(D12_ind, m)
                delta_12T = compute_delta_fft(D12T_ind, m)
                cost = compute_cost_vectorized(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n)
                if cost < best_cost:
                    best_cost = cost
                    d12_list = list(np.where(D12_ind[1:] == 1)[0] + 1)
                    not_in = list(np.where(D12_ind[1:] == 0)[0] + 1)
                    improved = True
                    print(f"    Improved to cost={cost}")
                    if cost == 0:
                        return cost, D12_ind
                    break
                else:
                    D12_ind[rem] = 1
                    D12_ind[add_el] = 0
                    D12T_ind[(-rem) % m] = 1
                    D12T_ind[(-add_el) % m] = 0
            if improved:
                break

    delta_12 = compute_delta_fft(D12_ind, m)
    delta_12T = compute_delta_fft(D12T_ind, m)
    return best_cost, D12_ind


def main():
    m = M
    n = N_PARAM

    print(f"Finding near-solution for n={n}, m={m}")
    print(f"Trying seeds 0-4 with |D11|=32, 8M iterations each...")
    print(flush=True)

    best_overall = 999
    best_data = None

    for seed in range(5):
        t0 = time.time()
        cost, D11, D12 = find_near_solution(32, seed, 8000000)
        elapsed = time.time() - t0
        print(f"\nSeed {seed}: best cost = {cost} ({elapsed:.0f}s)")

        if cost < best_overall:
            best_overall = cost
            best_data = (np.copy(D11), np.copy(D12), seed)

        if cost <= 8:
            # Analyze violations
            delta_11 = compute_delta_fft(D11, m)
            delta_12 = compute_delta_fft(D12, m)
            D12T_ind = np.zeros(m, dtype=np.int64)
            for x in range(m):
                if D12[x]:
                    D12T_ind[(-x) % m] = 1
            delta_12T = compute_delta_fft(D12T_ind, m)

            violations = get_violations(D11, D12, delta_11, delta_12, delta_12T, m, n)
            print(f"  Violations: {len(violations)}")
            for vtype, d, excess, val in violations:
                print(f"    {vtype} at d={d}: excess={excess}, value={val}, thresh={30 if 'red' in vtype else 31}")

            D11_list = sorted(int(i) for i in range(m) if D11[i])
            D12_list = sorted(int(i) for i in range(m) if D12[i])
            print(f"  D11 = {D11_list}")
            print(f"  D12 = {D12_list}")

            # CRT analysis of violated elements
            for vtype, d, excess, val in violations:
                a, b = d % 9, d % 7
                print(f"    d={d} -> ({a},{b}) mod (9,7)")

            # Try exhaustive D12 fix
            print(f"\n  Trying exhaustive D12 fix from this near-solution...")
            fixed_cost, fixed_D12 = exhaustive_d12_fix(np.copy(D11), np.copy(D12), m, n)
            print(f"  After fix: cost = {fixed_cost}")

        if cost == 0:
            print("SOLUTION FOUND!")
            break

    if best_overall > 0:
        print(f"\nBest overall cost: {best_overall}")

        # Try more D11 sizes
        print(f"\nTrying alternative |D11| sizes...")
        for d11_size in [28, 34, 26, 36]:
            print(f"\n--- |D11| = {d11_size} ---")
            for seed in range(3):
                t0 = time.time()
                cost, D11, D12 = find_near_solution(d11_size, seed, 5000000)
                elapsed = time.time() - t0
                print(f"  seed={seed}: cost={cost} ({elapsed:.0f}s)", flush=True)

                if cost < best_overall:
                    best_overall = cost
                    print(f"  ** NEW BEST: {cost} **")

                if cost == 0:
                    print("SOLUTION FOUND!")
                    D11_list = sorted(int(i) for i in range(m) if D11[i])
                    D12_list = sorted(int(i) for i in range(m) if D12[i])
                    print(f"  D11 = {D11_list}")
                    print(f"  D12 = {D12_list}")
                    sys.exit(0)


if __name__ == "__main__":
    main()
