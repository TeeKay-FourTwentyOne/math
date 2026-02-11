"""
Fast general SA solver using NumPy FFT.

Usage: python fast_sa_numpy.py <n> [max_seeds] [max_iter]

Key improvements over fast_sa_general.py:
- FFT-based autocorrelation: O(m log m) vs O(m^2)
- Vectorized cost function: no Python loops
- ~4x faster per iteration
"""

import sys, os, math, random, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ramsey_core import BlockCirculantGraph, verify_construction, save_construction


def compute_delta_fft(indicator, m):
    """Autocorrelation via FFT."""
    a = np.array(indicator, dtype=np.float64)
    fa = np.fft.rfft(a)
    result = np.fft.irfft(np.conj(fa) * fa, n=m)
    return np.round(result).astype(np.int64)


def compute_cost_vec(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n):
    """Vectorized cost computation."""
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

    red_excess_v11 = np.maximum(common_v11[d11_mask] - red_thresh, 0)
    blue_common_v11 = (N - 2) - 2 * d1 + common_v11[~d11_mask]
    blue_excess_v11 = np.maximum(blue_common_v11 - blue_thresh, 0)

    v22_const = m - 2 - 2 * d11_size
    d22_delta = delta_11[1:] + v22_const + 2 * D11_ind[1:]
    common_v22 = d22_delta + delta_12T[1:]

    red_excess_v22 = np.maximum(common_v22[~d11_mask] - red_thresh, 0)
    blue_common_v22 = (N - 2) - 2 * d2 + common_v22[d11_mask]
    blue_excess_v22 = np.maximum(blue_common_v22 - blue_thresh, 0)

    return (int(np.sum(red_excess_v11)) + int(np.sum(blue_excess_v11)) +
            int(np.sum(red_excess_v22)) + int(np.sum(blue_excess_v22)))


def solve_n(n_param, max_seeds=20, max_iter=10000000):
    m = 2 * n_param - 1
    N = 2 * m
    d12_size = n_param - 1

    pairs = []
    for x in range(1, m):
        neg_x = (-x) % m
        if x <= neg_x:
            pairs.append((x, neg_x))
    num_pairs = len(pairs)

    half = (m - 1) // 2
    d11_sizes = []
    for s in [half - 1, half, half + 1, half - 2, half + 2]:
        if s > 0 and s < m - 1 and s % 2 == 0 and (s // 2) <= num_pairs:
            d11_sizes.append(s)
    d11_sizes = sorted(set(d11_sizes), key=lambda s: abs(s - half))

    print(f"n={n_param}, m={m}, N={N}")
    print(f"Symmetric pairs: {num_pairs}")
    print(f"|D12| = {d12_size}")
    print(f"D11 sizes to try: {d11_sizes}")
    print(flush=True)

    start_time = time.time()

    for d11_size in d11_sizes:
        print(f"\n--- |D11| = {d11_size} ---", flush=True)
        for seed in range(max_seeds):
            random.seed(seed)
            np_rng = np.random.RandomState(seed)

            num_selected = d11_size // 2
            selected = random.sample(range(num_pairs), num_selected)
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

            current_cost = compute_cost_vec(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n_param)
            best_cost = current_cost
            best_D11 = np.copy(D11_ind)
            best_D12 = np.copy(D12_ind)

            T = 10.0
            T_min = 0.00005
            alpha = 1 - 6.0 / max_iter
            last_imp = 0
            reheat_count = 0

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
                    new_cost = compute_cost_vec(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n_param)
                    dc = new_cost - current_cost
                    if dc <= 0 or random.random() < math.exp(-dc / max(T, T_min)):
                        current_cost = new_cost
                    else:
                        D11_ind[pairs[rp][0]] = 1
                        D11_ind[pairs[rp][1]] = 1
                        D11_ind[pairs[ap][0]] = 0
                        D11_ind[pairs[ap][1]] = 0
                        delta_11 = old

                elif r < 0.82:
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
                    new_cost = compute_cost_vec(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n_param)
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

                else:
                    d12_list = np.where(D12_ind[1:] == 1)[0] + 1
                    not_in = np.where(D12_ind[1:] == 0)[0] + 1
                    if len(d12_list) < 2 or len(not_in) < 2:
                        continue
                    ri = np_rng.choice(len(d12_list), 2, replace=False)
                    ai = np_rng.choice(len(not_in), 2, replace=False)
                    rems = d12_list[ri]
                    adds = not_in[ai]
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
                    new_cost = compute_cost_vec(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n_param)
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

                if it - last_imp > 500000 and best_cost > 0 and best_cost <= 12:
                    T = 4.0
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

                if it % 1000000 == 0 and it > 0:
                    elapsed = time.time() - start_time
                    its = it / (time.time() - start_time + 0.001)
                    print(f"    {it//1000000}M: cost={current_cost} best={best_cost} T={T:.5f} rh={reheat_count} ({its:.0f} it/s)", flush=True)

            elapsed = time.time() - start_time
            print(f"  seed={seed}: best cost = {best_cost} ({elapsed:.1f}s)", flush=True)

            if best_cost == 0:
                D11_list = sorted(int(i) for i in range(m) if best_D11[i])
                D12_list = sorted(int(i) for i in range(m) if best_D12[i])
                print(f"\n  *** SOLUTION FOUND! ***")
                print(f"  D11 = {D11_list}")
                print(f"  D12 = {D12_list}")

                D11_s = set(D11_list)
                D12_s = set(D12_list)
                D22_s = set(range(1, m)) - D11_s
                G = BlockCirculantGraph(n=n_param, D11=D11_s, D12=D12_s, D22=D22_s)
                result = verify_construction(G)
                print(f"  Full verification: valid={result.valid}")
                print(f"  max_red={result.max_red_common} (thresh {result.red_threshold})")
                print(f"  max_blue={result.max_blue_common} (thresh {result.blue_threshold})")

                if result.valid:
                    fn = f"/Users/stephenpadgett/Projects/math/ramsey-book-graphs/solution_n{n_param}.json"
                    save_construction(G, result, fn, solver="fast-numpy-SA")
                    print(f"  SAVED to {fn}")
                    return True

    return False


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python fast_sa_numpy.py <n> [max_seeds] [max_iter]")
        sys.exit(1)

    n_param = int(sys.argv[1])
    max_seeds = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    max_iter = int(sys.argv[3]) if len(sys.argv) > 3 else 10000000

    solved = solve_n(n_param, max_seeds, max_iter)
    if not solved:
        print(f"\nNo solution found for n={n_param}")
