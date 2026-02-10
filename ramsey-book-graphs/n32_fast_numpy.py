"""
Ultra-fast SA solver for n=32 (m=63) using fully vectorized numpy.

Key: Vectorize EVERYTHING. No Python loops in the hot path.
FFT-based delta computation + vectorized cost function.
"""

import sys, os, math, random, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ramsey_core import BlockCirculantGraph, verify_construction, save_construction

N_PARAM = 32
M = 63
N_VERTS = 126
RED_THRESH = N_PARAM - 2  # 30
BLUE_THRESH = N_PARAM - 1  # 31
D12_SIZE = N_PARAM - 1  # 31


def compute_delta_fft(indicator, m):
    """Autocorrelation via FFT. Returns array of length m."""
    a = np.array(indicator, dtype=np.float64)
    fa = np.fft.rfft(a)
    result = np.fft.irfft(np.conj(fa) * fa, n=m)
    return np.round(result).astype(np.int64)


def compute_cost_vectorized(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n):
    """Fully vectorized cost computation. No Python loops."""
    d11_size = int(np.sum(D11_ind))
    d12_size = int(np.sum(D12_ind))
    d22_size = m - 1 - d11_size
    d1 = d11_size + d12_size
    d2 = d22_size + d12_size
    N = 2 * m

    red_thresh = n - 2
    blue_thresh = n - 1

    # V1V1 constraints (d=1..m-1)
    common_v11 = delta_11[1:] + delta_12[1:]
    d11_mask = D11_ind[1:].astype(bool)

    # Red V1V1: d in D11
    red_excess_v11 = np.maximum(common_v11[d11_mask] - red_thresh, 0)

    # Blue V1V1: d not in D11
    blue_common_v11 = (N - 2) - 2 * d1 + common_v11[~d11_mask]
    blue_excess_v11 = np.maximum(blue_common_v11 - blue_thresh, 0)

    # V2V2 constraints (d=1..m-1)
    v22_const = m - 2 - 2 * d11_size
    d22_delta = delta_11[1:] + v22_const + 2 * D11_ind[1:]
    common_v22 = d22_delta + delta_12T[1:]

    # Red V2V2: d not in D11 (d in D22)
    red_excess_v22 = np.maximum(common_v22[~d11_mask] - red_thresh, 0)

    # Blue V2V2: d in D11
    blue_common_v22 = (N - 2) - 2 * d2 + common_v22[d11_mask]
    blue_excess_v22 = np.maximum(blue_common_v22 - blue_thresh, 0)

    cost = (int(np.sum(red_excess_v11)) + int(np.sum(blue_excess_v11)) +
            int(np.sum(red_excess_v22)) + int(np.sum(blue_excess_v22)))

    return cost


def solve(d11_size, seed, max_iter=15000000):
    """Single SA run."""
    m = M
    n = N_PARAM

    rng = random.Random(seed)
    np_rng = np.random.RandomState(seed)

    # Build symmetric pairs
    pairs = []
    for x in range(1, m):
        neg_x = (-x) % m
        if x <= neg_x:
            pairs.append((x, neg_x))
    num_pairs = len(pairs)

    # Initialize D11
    num_selected = d11_size // 2
    selected = rng.sample(range(num_pairs), num_selected)
    D11_ind = np.zeros(m, dtype=np.int64)
    for i in selected:
        D11_ind[pairs[i][0]] = 1
        D11_ind[pairs[i][1]] = 1

    # Initialize D12
    d12_rest = rng.sample(range(1, m), D12_SIZE - 1)
    D12_ind = np.zeros(m, dtype=np.int64)
    D12_ind[0] = 1
    for x in d12_rest:
        D12_ind[x] = 1

    # Precompute pair membership for fast lookup
    pair_of = {}
    for i, (a, b) in enumerate(pairs):
        pair_of[a] = i
        pair_of[b] = i

    # Compute initial deltas
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
    T_min = 0.00005
    alpha = 1 - 6.0 / max_iter

    last_improvement = 0
    reheat_count = 0
    t0 = time.time()

    # Pre-generate random numbers in batches for speed
    batch_size = 100000
    rand_batch = np_rng.random(batch_size)
    rand_idx = 0

    def next_rand():
        nonlocal rand_batch, rand_idx
        if rand_idx >= batch_size:
            rand_batch = np_rng.random(batch_size)
            rand_idx = 0
        v = rand_batch[rand_idx]
        rand_idx += 1
        return v

    for it in range(max_iter):
        r = next_rand()

        if r < 0.22:
            # D11 pair swap
            in_pairs = [i for i in range(num_pairs) if D11_ind[pairs[i][0]]]
            out_pairs = [i for i in range(num_pairs) if not D11_ind[pairs[i][0]]]
            if not in_pairs or not out_pairs:
                continue
            rp = in_pairs[int(next_rand() * len(in_pairs))]
            ap = out_pairs[int(next_rand() * len(out_pairs))]

            old_delta = delta_11.copy()
            D11_ind[pairs[rp][0]] = 0
            D11_ind[pairs[rp][1]] = 0
            D11_ind[pairs[ap][0]] = 1
            D11_ind[pairs[ap][1]] = 1
            delta_11 = compute_delta_fft(D11_ind, m)

            new_cost = compute_cost_vectorized(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n)
            dc = new_cost - current_cost

            if dc <= 0 or next_rand() < math.exp(-dc / max(T, T_min)):
                current_cost = new_cost
            else:
                D11_ind[pairs[rp][0]] = 1
                D11_ind[pairs[rp][1]] = 1
                D11_ind[pairs[ap][0]] = 0
                D11_ind[pairs[ap][1]] = 0
                delta_11 = old_delta

        elif r < 0.82:
            # D12 single swap
            d12_list = np.where(D12_ind[1:] == 1)[0] + 1
            not_in_list = np.where(D12_ind[1:] == 0)[0] + 1
            if len(d12_list) == 0 or len(not_in_list) == 0:
                continue
            rem = d12_list[int(next_rand() * len(d12_list))]
            add_el = not_in_list[int(next_rand() * len(not_in_list))]

            old_d12 = delta_12.copy()
            old_d12T = delta_12T.copy()
            D12_ind[rem] = 0
            D12_ind[add_el] = 1
            neg_rem = (-rem) % m
            neg_add = (-add_el) % m
            D12T_ind[neg_rem] = 0
            D12T_ind[neg_add] = 1

            delta_12 = compute_delta_fft(D12_ind, m)
            delta_12T = compute_delta_fft(D12T_ind, m)

            new_cost = compute_cost_vectorized(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n)
            dc = new_cost - current_cost

            if dc <= 0 or next_rand() < math.exp(-dc / max(T, T_min)):
                current_cost = new_cost
            else:
                D12_ind[rem] = 1
                D12_ind[add_el] = 0
                D12T_ind[neg_rem] = 1
                D12T_ind[neg_add] = 0
                delta_12 = old_d12
                delta_12T = old_d12T

        else:
            # Multi-swap: 2 D12 elements
            d12_list = np.where(D12_ind[1:] == 1)[0] + 1
            not_in_list = np.where(D12_ind[1:] == 0)[0] + 1
            if len(d12_list) < 2 or len(not_in_list) < 2:
                continue
            ri = np_rng.choice(len(d12_list), 2, replace=False)
            ai = np_rng.choice(len(not_in_list), 2, replace=False)
            rems = d12_list[ri]
            adds = not_in_list[ai]

            old_d12 = delta_12.copy()
            old_d12T = delta_12T.copy()
            for r in rems:
                D12_ind[r] = 0
                D12T_ind[(-r) % m] = 0
            for a in adds:
                D12_ind[a] = 1
                D12T_ind[(-a) % m] = 1

            delta_12 = compute_delta_fft(D12_ind, m)
            delta_12T = compute_delta_fft(D12T_ind, m)

            new_cost = compute_cost_vectorized(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n)
            dc = new_cost - current_cost

            if dc <= 0 or next_rand() < math.exp(-dc / max(T, T_min)):
                current_cost = new_cost
            else:
                for r in rems:
                    D12_ind[r] = 1
                    D12T_ind[(-r) % m] = 1
                for a in adds:
                    D12_ind[a] = 0
                    D12T_ind[(-a) % m] = 0
                delta_12 = old_d12
                delta_12T = old_d12T

        if current_cost < best_cost:
            best_cost = current_cost
            best_D11 = np.copy(D11_ind)
            best_D12 = np.copy(D12_ind)
            last_improvement = it
            if best_cost == 0:
                break

        # Reheat if stuck for 600K iterations at low cost
        if it - last_improvement > 600000 and best_cost > 0 and best_cost <= 12:
            T = 5.0
            last_improvement = it
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

        if it % 1000000 == 0 and it > 0:
            elapsed = time.time() - t0
            iters_per_sec = it / elapsed
            print(f"    {it//1000000}M: cost={current_cost} best={best_cost} T={T:.5f} rh={reheat_count} ({elapsed:.0f}s, {iters_per_sec:.0f} it/s)", flush=True)

    elapsed = time.time() - t0
    return best_cost, best_D11, best_D12, elapsed, reheat_count


def main():
    max_seeds = int(sys.argv[1]) if len(sys.argv) > 1 else 50
    max_iter = int(sys.argv[2]) if len(sys.argv) > 2 else 15000000

    d11_sizes = [32, 30]

    print(f"n={N_PARAM}, m={M}, N={N_VERTS}")
    print(f"|D12| = {D12_SIZE}")
    print(f"D11 sizes: {d11_sizes}")
    print(f"Seeds per size: {max_seeds}, Iterations: {max_iter}")
    print(flush=True)

    global_best = 999
    t0 = time.time()

    for d11_size in d11_sizes:
        print(f"\n=== |D11| = {d11_size} ===", flush=True)
        for seed in range(max_seeds):
            best, best_D11, best_D12, elapsed, reheats = solve(d11_size, seed, max_iter)
            total = time.time() - t0
            print(f"  seed={seed}: best={best} reheats={reheats} ({total:.0f}s)", flush=True)

            if best < global_best:
                global_best = best

            if best == 0:
                D11_list = sorted(int(i) for i in range(M) if best_D11[i])
                D12_list = sorted(int(i) for i in range(M) if best_D12[i])
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
                    save_construction(G, result, fn, solver="fast-numpy-SA")
                    print(f"  SAVED to {fn}")
                    sys.exit(0)

    print(f"\nNo solution found. Best cost = {global_best}")


if __name__ == "__main__":
    main()
