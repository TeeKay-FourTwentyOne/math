"""
n=32 solver with very slow cooling and NO reheat.

The reheat strategy may be hurting by constantly resetting from cost=8.
Let's try a single long SA run with very slow cooling.
"""

import sys, os, math, random, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ramsey_core import BlockCirculantGraph, verify_construction, save_construction

M = 63
N_PARAM = 32
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


def solve(d11_size, seed, max_iter=30000000):
    m = M
    n = N_PARAM
    random.seed(seed)

    pairs = []
    for x in range(1, m):
        neg_x = (-x) % m
        if x <= neg_x:
            pairs.append((x, neg_x))
    num_pairs = len(pairs)

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

    # Very slow cooling, no reheat
    T = 15.0
    T_min = 0.00001
    # Cool from 15 to 0.00001 over 30M iterations
    # alpha = exp(ln(T_min/T_start) / max_iter) = exp(ln(0.00001/15) / 30M)
    alpha = math.exp(math.log(T_min / T) / max_iter)

    t0 = time.time()

    for it in range(max_iter):
        r = random.random()

        if r < 0.25:
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

        else:
            # Multi-swap D12
            d12_list = np.where(D12_ind[1:] == 1)[0] + 1
            not_in = np.where(D12_ind[1:] == 0)[0] + 1
            if len(d12_list) < 2 or len(not_in) < 2:
                continue
            ri = [random.randint(0, len(d12_list)-1)]
            ri.append(random.randint(0, len(d12_list)-2))
            if ri[1] >= ri[0]:
                ri[1] += 1
            ai = [random.randint(0, len(not_in)-1)]
            ai.append(random.randint(0, len(not_in)-2))
            if ai[1] >= ai[0]:
                ai[1] += 1
            rems = [d12_list[i] for i in ri]
            adds = [not_in[i] for i in ai]
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
            if dc <= 0 or random.random() < math.exp(-dc / T):
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
            if best_cost == 0:
                break

        T *= alpha

        if it % 2000000 == 0 and it > 0:
            elapsed = time.time() - t0
            its = it / elapsed
            print(f"    {it//1000000}M: cost={current_cost} best={best_cost} T={T:.5f} ({its:.0f} it/s)", flush=True)

    elapsed = time.time() - t0
    return best_cost, best_D11, best_D12, elapsed


def main():
    max_seeds = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    max_iter = int(sys.argv[2]) if len(sys.argv) > 2 else 30000000

    print(f"n={N_PARAM}, m={M}")
    print(f"Very slow cooling, NO reheat")
    print(f"Seeds: {max_seeds}, Iter: {max_iter}", flush=True)

    t0 = time.time()
    global_best = 999

    for seed in range(max_seeds):
        best, D11, D12, elapsed = solve(32, seed, max_iter)
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
            if result.valid:
                fn = f"/Users/stephenpadgett/Projects/math/ramsey-book-graphs/solution_n{N_PARAM}.json"
                save_construction(G, result, fn, solver="slow-cool-SA")
                print(f"  SAVED to {fn}")
                sys.exit(0)

    print(f"\nBest = {global_best}")


if __name__ == "__main__":
    main()
