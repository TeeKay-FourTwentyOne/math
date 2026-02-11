"""
Massive SA run for n=32 with tuned parameters.

Key insights from analysis:
1. |D11|=30 or 32 are the only viable sizes (slack -0.03 and -1.00)
2. V1V1 and V2V2 are identical - only need to satisfy V1V1
3. Violations always come in complementary pairs (d, m-d)
4. All violations have excess exactly 1 (one unit over threshold)
5. D12 fix alone can't resolve violations - need joint D11/D12 changes

Strategy:
- Longer runs (20M+ iterations per seed)
- More aggressive D11 mutations (30% instead of 22%)
- Coupled D11+D12 moves to escape paired violations
- Slower cooling for better exploration
- Try removing the 0-in-D12 constraint
"""

import sys, os, math, random, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ramsey_core import BlockCirculantGraph, verify_construction, save_construction

N_PARAM = 32
M = 63
D12_SIZE = N_PARAM - 1  # 31


def compute_delta_fft(indicator, m):
    a = np.array(indicator, dtype=np.float64)
    fa = np.fft.rfft(a)
    result = np.fft.irfft(np.conj(fa) * fa, n=m)
    return np.round(result).astype(np.int64)


def compute_cost(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n):
    d11_size = int(np.sum(D11_ind))
    d12_size = int(np.sum(D12_ind))
    d22_size = m - 1 - d11_size
    d1 = d11_size + d12_size
    d2 = d22_size + d12_size
    N = 2 * m
    red_thresh = n - 2
    blue_thresh = n - 1

    common_v11 = delta_11[1:] + delta_12[1:]
    d11_mask = D11_ind[1:].astype(bool)
    red_excess = np.maximum(common_v11[d11_mask] - red_thresh, 0)
    blue_common = (N - 2) - 2 * d1 + common_v11[~d11_mask]
    blue_excess = np.maximum(blue_common - blue_thresh, 0)

    v22_const = m - 2 - 2 * d11_size
    d22_delta = delta_11[1:] + v22_const + 2 * D11_ind[1:]
    common_v22 = d22_delta + delta_12T[1:]
    red_excess_v22 = np.maximum(common_v22[~d11_mask] - red_thresh, 0)
    blue_common_v22 = (N - 2) - 2 * d2 + common_v22[d11_mask]
    blue_excess_v22 = np.maximum(blue_common_v22 - blue_thresh, 0)

    return (int(np.sum(red_excess)) + int(np.sum(blue_excess)) +
            int(np.sum(red_excess_v22)) + int(np.sum(blue_excess_v22)))


def solve(d11_size, seed, max_iter=20000000, require_0_in_d12=True):
    m = M
    n = N_PARAM
    random.seed(seed)

    pairs = []
    for x in range(1, m):
        neg_x = (-x) % m
        if x <= neg_x:
            pairs.append((x, neg_x))
    num_pairs = len(pairs)

    # Init D11
    selected = random.sample(range(num_pairs), d11_size // 2)
    D11_ind = np.zeros(m, dtype=np.int64)
    for i in selected:
        D11_ind[pairs[i][0]] = 1
        D11_ind[pairs[i][1]] = 1

    # Init D12
    D12_ind = np.zeros(m, dtype=np.int64)
    if require_0_in_d12:
        D12_ind[0] = 1
        rest = random.sample(range(1, m), D12_SIZE - 1)
        for x in rest:
            D12_ind[x] = 1
    else:
        rest = random.sample(range(0, m), D12_SIZE)
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
    best_D11 = np.copy(D11_ind)
    best_D12 = np.copy(D12_ind)

    T = 12.0
    T_min = 0.00003
    alpha = 1 - 5.0 / max_iter  # Slower cooling
    last_imp = 0
    reheat_count = 0
    t0 = time.time()

    # Determine which D12 elements are protected from swaps
    protected = set()
    if require_0_in_d12:
        protected = {0}

    for it in range(max_iter):
        r = random.random()

        if r < 0.30:
            # D11 pair swap (more frequent than before)
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
            if dc <= 0 or random.random() < math.exp(-dc / max(T, T_min)):
                current_cost = new_cost
            else:
                D11_ind[pairs[rp][0]] = 1
                D11_ind[pairs[rp][1]] = 1
                D11_ind[pairs[ap][0]] = 0
                D11_ind[pairs[ap][1]] = 0
                delta_11 = old

        elif r < 0.80:
            # D12 single swap
            d12_swappable = [x for x in range(m) if D12_ind[x] and x not in protected]
            not_in = [x for x in range(m) if not D12_ind[x]]
            if not d12_swappable or not not_in:
                continue
            rem = random.choice(d12_swappable)
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
            # Coupled move: swap D11 pair AND D12 element simultaneously
            in_pairs = [i for i in range(num_pairs) if D11_ind[pairs[i][0]]]
            out_pairs = [i for i in range(num_pairs) if not D11_ind[pairs[i][0]]]
            d12_swappable = [x for x in range(m) if D12_ind[x] and x not in protected]
            d12_not_in = [x for x in range(m) if not D12_ind[x]]
            if not in_pairs or not out_pairs or not d12_swappable or not d12_not_in:
                continue

            rp = random.choice(in_pairs)
            ap = random.choice(out_pairs)
            d12_rem = random.choice(d12_swappable)
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
            # Multi-swap D12 (2 elements)
            d12_swappable = [x for x in range(m) if D12_ind[x] and x not in protected]
            d12_not_in = [x for x in range(m) if not D12_ind[x]]
            if len(d12_swappable) < 2 or len(d12_not_in) < 2:
                continue
            rems = random.sample(d12_swappable, 2)
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

        # Reheat if stuck
        if it - last_imp > 700000 and best_cost > 0 and best_cost <= 16:
            T = 5.0
            last_imp = it
            reheat_count += 1
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

        if it % 2000000 == 0 and it > 0:
            elapsed = time.time() - t0
            its = it / elapsed
            print(f"    {it//1000000}M: cost={current_cost} best={best_cost} T={T:.5f} rh={reheat_count} ({its:.0f} it/s)", flush=True)

    elapsed = time.time() - t0
    return best_cost, best_D11, best_D12, elapsed, reheat_count


def main():
    max_seeds = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    max_iter = int(sys.argv[2]) if len(sys.argv) > 2 else 20000000

    d11_sizes = [30, 32]  # 30 first (better slack)

    print(f"n={N_PARAM}, m={M}")
    print(f"Max seeds per size: {max_seeds}")
    print(f"Max iter per seed: {max_iter}")
    print(f"D11 sizes: {d11_sizes}")
    print(flush=True)

    t0 = time.time()
    global_best = 999

    for d11_size in d11_sizes:
        print(f"\n=== |D11| = {d11_size} ===", flush=True)
        for seed in range(max_seeds):
            best, D11, D12, elapsed, reheats = solve(d11_size, seed, max_iter)
            total = time.time() - t0
            print(f"  seed={seed}: best={best} rh={reheats} ({total:.0f}s)", flush=True)

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
                    save_construction(G, result, fn, solver="massive-SA")
                    print(f"  SAVED to {fn}")
                    sys.exit(0)

    print(f"\nNo solution found. Best = {global_best}")


if __name__ == "__main__":
    main()
