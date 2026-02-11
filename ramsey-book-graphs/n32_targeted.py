"""
Targeted SA for n=32 that exploits the violation structure.

Key insight: Violations come in complementary pairs (d, m-d) with excess exactly 1.
The constraint A(d)+B(d) <= 30 fails at exactly 2 positions (giving cost=4).

Strategy:
1. Run standard SA until cost <= 4
2. Switch to targeted mode: identify the 2 violated positions
3. Try systematic perturbations that specifically reduce A(d)+B(d) at those positions
"""

import sys, os, math, random, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ramsey_core import BlockCirculantGraph, verify_construction, save_construction

N_PARAM = 32
M = 63
D12_SIZE = 31


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


def get_violated_diffs(D11_ind, D12_ind, delta_11, delta_12, m, n):
    """Get list of d values where A(d)+B(d) > n-2 and d in D11."""
    red_thresh = n - 2
    violated = []
    for d in range(1, m):
        if D11_ind[d]:
            val = int(delta_11[d] + delta_12[d])
            if val > red_thresh:
                violated.append((d, val))
    return violated


def targeted_search(D11_ind, D12_ind, max_iter=5000000):
    """Targeted search around a near-solution to eliminate remaining violations."""
    m = M
    n = N_PARAM

    pairs = []
    for x in range(1, m):
        neg_x = (-x) % m
        if x <= neg_x:
            pairs.append((x, neg_x))
    num_pairs = len(pairs)

    D12T_ind = np.zeros(m, dtype=np.int64)
    for x in range(m):
        if D12_ind[x]:
            D12T_ind[(-x) % m] = 1

    delta_11 = compute_delta_fft(D11_ind, m)
    delta_12 = compute_delta_fft(D12_ind, m)
    delta_12T = compute_delta_fft(D12T_ind, m)

    best_cost = compute_cost(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n)
    best_D11 = np.copy(D11_ind)
    best_D12 = np.copy(D12_ind)

    if best_cost == 0:
        return 0, best_D11, best_D12

    # SA with warm restarts
    T = 3.0
    T_min = 0.00001
    alpha = 1 - 4.0 / max_iter
    current_cost = best_cost
    last_imp = 0

    for it in range(max_iter):
        r = random.random()

        if r < 0.35:
            # D11 pair swap
            in_pairs = [i for i in range(num_pairs) if D11_ind[pairs[i][0]]]
            out_pairs = [i for i in range(num_pairs) if not D11_ind[pairs[i][0]]]
            if not in_pairs or not out_pairs:
                continue

            # Bias towards swapping pairs that are involved in violations
            violated = get_violated_diffs(D11_ind, D12_ind, delta_11, delta_12, m, n)
            if violated and random.random() < 0.5:
                # Try to swap out a pair containing a violated difference
                v_diffs = {d for d, _ in violated}
                relevant_in = [i for i in in_pairs if pairs[i][0] in v_diffs or pairs[i][1] in v_diffs]
                if relevant_in:
                    rp = random.choice(relevant_in)
                else:
                    rp = random.choice(in_pairs)
            else:
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
            if dc <= 0 or random.random() < math.exp(-dc / max(T, T_min)):
                current_cost = new_cost
            else:
                D11_ind[pairs[rp][0]] = 1
                D11_ind[pairs[rp][1]] = 1
                D11_ind[pairs[ap][0]] = 0
                D11_ind[pairs[ap][1]] = 0
                delta_11 = old

        elif r < 0.75:
            # D12 swap
            d12_list = [x for x in range(1, m) if D12_ind[x]]
            not_in = [x for x in range(1, m) if not D12_ind[x]]
            if not d12_list or not not_in:
                continue

            # Bias: try removing D12 elements that contribute to violations
            violated = get_violated_diffs(D11_ind, D12_ind, delta_11, delta_12, m, n)
            if violated and random.random() < 0.3:
                v_diffs = {d for d, _ in violated}
                # Elements x in D12 contribute to Delta(D12,D12,d) when (x-d) in D12
                # Removing x from D12 reduces Delta(D12,D12,d) by at most 2
                # Prefer removing x that contributes to violated d
                contributors = set()
                for d in v_diffs:
                    for x in d12_list:
                        if D12_ind[(x - d) % m]:
                            contributors.add(x)
                if contributors:
                    rem = random.choice(list(contributors & set(d12_list)))
                else:
                    rem = random.choice(d12_list)
            else:
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
            new_cost = compute_cost(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n)
            dc = new_cost - current_cost
            if dc <= 0 or random.random() < math.exp(-dc / max(T, T_min)):
                current_cost = new_cost
            else:
                D12_ind[rem] = 1
                D12_ind[add_el] = 0
                D12T_ind[(-rem) % m] = 1
                D12T_ind[(-add_el) % m] = 0
                delta_12 = old_d12
                delta_12T = old_d12T

        elif r < 0.90:
            # Coupled D11+D12 swap
            in_pairs = [i for i in range(num_pairs) if D11_ind[pairs[i][0]]]
            out_pairs = [i for i in range(num_pairs) if not D11_ind[pairs[i][0]]]
            d12_list = [x for x in range(1, m) if D12_ind[x]]
            d12_not_in = [x for x in range(1, m) if not D12_ind[x]]
            if not in_pairs or not out_pairs or not d12_list or not d12_not_in:
                continue
            rp = random.choice(in_pairs)
            ap = random.choice(out_pairs)
            d12_rem = random.choice(d12_list)
            d12_add = random.choice(d12_not_in)

            old_d11 = delta_11.copy()
            old_d12 = delta_12.copy()
            old_d12T = delta_12T.copy()
            D11_ind[pairs[rp][0]] = 0
            D11_ind[pairs[rp][1]] = 0
            D11_ind[pairs[ap][0]] = 1
            D11_ind[pairs[ap][1]] = 1
            D12_ind[d12_rem] = 0
            D12_ind[d12_add] = 1
            D12T_ind[(-d12_rem) % m] = 0
            D12T_ind[(-d12_add) % m] = 1
            delta_11 = compute_delta_fft(D11_ind, m)
            delta_12 = compute_delta_fft(D12_ind, m)
            delta_12T = compute_delta_fft(D12T_ind, m)
            new_cost = compute_cost(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n)
            dc = new_cost - current_cost
            if dc <= 0 or random.random() < math.exp(-dc / max(T, T_min)):
                current_cost = new_cost
            else:
                D11_ind[pairs[rp][0]] = 1
                D11_ind[pairs[rp][1]] = 1
                D11_ind[pairs[ap][0]] = 0
                D11_ind[pairs[ap][1]] = 0
                D12_ind[d12_rem] = 1
                D12_ind[d12_add] = 0
                D12T_ind[(-d12_rem) % m] = 1
                D12T_ind[(-d12_add) % m] = 0
                delta_11 = old_d11
                delta_12 = old_d12
                delta_12T = old_d12T
        else:
            # 2x D12 swap
            d12_list = [x for x in range(1, m) if D12_ind[x]]
            d12_not_in = [x for x in range(1, m) if not D12_ind[x]]
            if len(d12_list) < 2 or len(d12_not_in) < 2:
                continue
            rems = random.sample(d12_list, 2)
            adds = random.sample(d12_not_in, 2)
            old_d12 = delta_12.copy()
            old_d12T = delta_12T.copy()
            for rv in rems:
                D12_ind[rv] = 0
                D12T_ind[(-rv) % m] = 0
            for av in adds:
                D12_ind[av] = 1
                D12T_ind[(-av) % m] = 1
            delta_12 = compute_delta_fft(D12_ind, m)
            delta_12T = compute_delta_fft(D12T_ind, m)
            new_cost = compute_cost(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n)
            dc = new_cost - current_cost
            if dc <= 0 or random.random() < math.exp(-dc / max(T, T_min)):
                current_cost = new_cost
            else:
                for rv in rems:
                    D12_ind[rv] = 1
                    D12T_ind[(-rv) % m] = 1
                for av in adds:
                    D12_ind[av] = 0
                    D12T_ind[(-av) % m] = 0
                delta_12 = old_d12
                delta_12T = old_d12T

        if current_cost < best_cost:
            best_cost = current_cost
            best_D11 = np.copy(D11_ind)
            best_D12 = np.copy(D12_ind)
            last_imp = it
            if best_cost == 0:
                break

        if it - last_imp > 400000 and best_cost > 0:
            T = 3.0
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

        T = max(T * alpha, T_min)

    return best_cost, best_D11, best_D12


def main():
    m = M
    n = N_PARAM
    d11_size = 32

    max_seeds = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    max_iter = int(sys.argv[2]) if len(sys.argv) > 2 else 15000000

    pairs = []
    for x in range(1, m):
        neg_x = (-x) % m
        if x <= neg_x:
            pairs.append((x, neg_x))
    num_pairs = len(pairs)

    print(f"n={n}, m={m}, |D11|={d11_size}")
    print(f"Seeds: {max_seeds}, Iter per seed: {max_iter}")
    print(f"Using violation-biased moves", flush=True)

    t0 = time.time()
    global_best = 999

    for seed in range(max_seeds):
        random.seed(seed)
        selected = random.sample(range(num_pairs), d11_size // 2)
        D11_ind = np.zeros(m, dtype=np.int64)
        for i in selected:
            D11_ind[pairs[i][0]] = 1
            D11_ind[pairs[i][1]] = 1

        D12_ind = np.zeros(m, dtype=np.int64)
        D12_ind[0] = 1
        rest = random.sample(range(1, m), D12_SIZE - 1)
        for x in rest:
            D12_ind[x] = 1

        best, D11, D12 = targeted_search(D11_ind, D12_ind, max_iter)
        total = time.time() - t0
        print(f"  seed={seed}: best={best} ({total:.0f}s)", flush=True)

        if best < global_best:
            global_best = best

        if best == 0:
            D11_list = sorted(int(i) for i in range(M) if D11[i])
            D12_list = sorted(int(i) for i in range(M) if D12[i])
            print(f"\n  *** SOLUTION FOUND! ***")
            print(f"  D11 = {D11_list}")
            print(f"  D12 = {D12_list}")
            D11_s = set(D11_list)
            D12_s = set(D12_list)
            D22_s = set(range(1, M)) - D11_s
            G = BlockCirculantGraph(n=N_PARAM, D11=D11_s, D12=D12_s, D22=D22_s)
            result = verify_construction(G)
            print(f"  Full verification: valid={result.valid}")
            print(f"  max_red={result.max_red_common} (thresh {result.red_threshold})")
            print(f"  max_blue={result.max_blue_common} (thresh {result.blue_threshold})")
            if result.valid:
                fn = f"/Users/stephenpadgett/Projects/math/ramsey-book-graphs/solution_n{N_PARAM}.json"
                save_construction(G, result, fn, solver="targeted-SA")
                print(f"  SAVED to {fn}")
                sys.exit(0)

    print(f"\nNo solution found. Best = {global_best}")


if __name__ == "__main__":
    main()
