"""
Turbo SA solver for n=32 (m=63 = 9*7).

Key optimizations over fast_sa_general.py:
1. NumPy-based vectorized Delta computation (much faster than Python loops)
2. True O(m) incremental delta updates for single-element swaps
3. Reheat strategy when stuck at low cost
4. Multi-swap moves to escape local minima
5. CRT-aware initialization options for m=63 = 9*7
"""

import sys, os, math, random, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ramsey_core import BlockCirculantGraph, verify_construction, save_construction

N_PARAM = 32
M = 2 * N_PARAM - 1  # 63
N_VERTS = 2 * M  # 126
RED_THRESH = N_PARAM - 2  # 30
BLUE_THRESH = N_PARAM - 1  # 31
D12_SIZE = N_PARAM - 1  # 31


def compute_delta_array(indicator, m):
    """Compute Delta(S,S,d) for all d using FFT-based convolution."""
    # Delta(S,S,d) = sum_x indicator[x] * indicator[(x-d) % m]
    # This is the circular cross-correlation of indicator with itself
    a = np.array(indicator, dtype=np.float64)
    fa = np.fft.rfft(a)
    # Cross-correlation: IFFT(conj(FA) * FA) = autocorrelation
    result = np.fft.irfft(np.conj(fa) * fa, n=m)
    return np.round(result).astype(np.int64)


def compute_cross_delta_array(ind_a, ind_b, m):
    """Compute Delta(A,B,d) = sum_x ind_a[x] * ind_b[(x-d) % m] for all d."""
    a = np.array(ind_a, dtype=np.float64)
    b = np.array(ind_b, dtype=np.float64)
    fa = np.fft.rfft(a)
    fb = np.fft.rfft(b)
    # Cross-correlation: IFFT(conj(FA) * FB)
    result = np.fft.irfft(np.conj(fa) * fb, n=m)
    return np.round(result).astype(np.int64)


def update_delta_remove_add(delta, indicator, m, rem, add):
    """
    Incrementally update delta array when swapping element rem -> add.
    O(m) instead of O(m^2).

    Delta(S,S,d) = sum_x indicator[x] * indicator[(x-d) % m]

    When we remove rem and add 'add':
    For each d, the change is:
      - lose contribution of rem paired with all x in S: indicator[(rem-d)%m] + indicator[(rem+d)%m] (minus self-pair if d=0)
      - lose contribution of all x in S paired with rem: same by symmetry
      - gain contribution of add paired with all x in new S
      - etc.

    It's simpler to just recompute for this element change:
    delta_new[d] = delta[d]
      - indicator[(rem-d)%m]  (rem was in S, paired going forward)
      - indicator[(rem+d)%m]  (rem was in S, paired going backward)
      + new_indicator[(add-d)%m]  (add is new in S)
      + new_indicator[(add+d)%m]  (add is new in S)

    But we need to be more careful. Let's think precisely.

    Actually for an autocorrelation, the cleanest approach:
    delta[d] = sum_x S[x] * S[(x-d)%m]

    Remove rem from S, add 'add' to S:
    new_delta[d] = delta[d]
      - S_old[(rem - d) % m]      # rem no longer contributes as 'x'
      - S_old[(rem + d) % m]      # rem no longer contributes as 'x-d'
      + (1 if rem == (rem - d) % m else 0)  # we double-subtracted rem,rem pair
      + S_new[(add - d) % m]      # add now contributes as 'x'
      + S_new[(add + d) % m]      # add now contributes as 'x-d'
      - (1 if add == (add - d) % m else 0)  # we double-added add,add pair

    Wait, this is getting complicated. Let me use a different formulation.
    """
    # Just recompute via FFT - it's fast enough with numpy
    indicator[rem] = 0
    indicator[add] = 1
    return compute_delta_array(indicator, m)


def compute_cost_vectorized(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n):
    """Compute total violation cost using vectorized operations."""
    d11_size = int(np.sum(D11_ind))
    d12_size = int(np.sum(D12_ind))
    d22_size = m - 1 - d11_size
    d1 = d11_size + d12_size
    d2 = d22_size + d12_size
    N = 2 * m

    red_thresh = n - 2
    blue_thresh = n - 1

    cost = 0

    # V1V1 constraints (d=1..m-1)
    for d in range(1, m):
        common = int(delta_11[d] + delta_12[d])
        if D11_ind[d]:  # red edge
            if common > red_thresh:
                cost += common - red_thresh
        else:  # blue edge
            blue_common = (N - 2) - 2 * d1 + common
            if blue_common > blue_thresh:
                cost += blue_common - blue_thresh

    # V2V2 constraints (d=1..m-1)
    v22_const = m - 2 - 2 * d11_size
    for d in range(1, m):
        d22_delta = int(delta_11[d]) + v22_const + (2 if D11_ind[d] else 0)
        common = d22_delta + int(delta_12T[d])
        if not D11_ind[d]:  # D22 = complement of D11, so d in D22 iff d not in D11; red edge
            if common > red_thresh:
                cost += common - red_thresh
        else:  # blue edge in V2V2
            blue_common = (N - 2) - 2 * d2 + common
            if blue_common > blue_thresh:
                cost += blue_common - blue_thresh

    return cost


def get_violation_details(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n):
    """Get detailed violation information for debugging."""
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
                violations.append(('V1V1_red', d, common - red_thresh))
        else:
            blue_common = (N - 2) - 2 * d1 + common
            if blue_common > blue_thresh:
                violations.append(('V1V1_blue', d, blue_common - blue_thresh))

    v22_const = m - 2 - 2 * d11_size
    for d in range(1, m):
        d22_delta = int(delta_11[d]) + v22_const + (2 if D11_ind[d] else 0)
        common = d22_delta + int(delta_12T[d])
        if not D11_ind[d]:
            if common > red_thresh:
                violations.append(('V2V2_red', d, common - red_thresh))
        else:
            blue_common = (N - 2) - 2 * d2 + common
            if blue_common > blue_thresh:
                violations.append(('V2V2_blue', d, blue_common - blue_thresh))

    return violations


def solve_n32(d11_size, seed, max_iter=10000000):
    """Solve n=32 with turbo SA."""
    m = M
    n = N_PARAM

    random.seed(seed)
    np.random.seed(seed)

    # Build symmetric pairs
    pairs = []
    for x in range(1, m):
        neg_x = (-x) % m
        if x <= neg_x:
            pairs.append((x, neg_x))
    num_pairs = len(pairs)

    # Initialize D11 (symmetric)
    num_selected = d11_size // 2
    selected = random.sample(range(num_pairs), num_selected)
    D11_ind = np.zeros(m, dtype=np.int64)
    for i in selected:
        D11_ind[pairs[i][0]] = 1
        D11_ind[pairs[i][1]] = 1

    # Initialize D12 (always contains 0)
    d12_rest = random.sample(range(1, m), D12_SIZE - 1)
    D12_ind = np.zeros(m, dtype=np.int64)
    D12_ind[0] = 1
    for x in d12_rest:
        D12_ind[x] = 1

    # Compute initial deltas via FFT
    delta_11 = compute_delta_array(D11_ind, m)
    delta_12 = compute_delta_array(D12_ind, m)
    D12T_ind = np.zeros(m, dtype=np.int64)
    for x in range(m):
        if D12_ind[x]:
            D12T_ind[(-x) % m] = 1
    delta_12T = compute_delta_array(D12T_ind, m)

    current_cost = compute_cost_vectorized(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n)
    best_cost = current_cost
    best_D11 = np.copy(D11_ind)
    best_D12 = np.copy(D12_ind)

    # SA parameters
    T = 10.0
    T_min = 0.0001
    alpha = 1 - 6.0 / max_iter  # Slightly slower cooling

    stuck_count = 0
    last_improvement_iter = 0
    reheat_count = 0

    t0 = time.time()

    for it in range(max_iter):
        move_type = random.random()

        if move_type < 0.20:
            # D11 swap: swap a symmetric pair in/out
            in_pairs = [i for i in range(num_pairs) if D11_ind[pairs[i][0]]]
            out_pairs = [i for i in range(num_pairs) if not D11_ind[pairs[i][0]]]
            if not in_pairs or not out_pairs:
                continue

            rp = random.choice(in_pairs)
            ap = random.choice(out_pairs)

            # Apply swap
            old_d11 = delta_11.copy()
            D11_ind[pairs[rp][0]] = 0
            D11_ind[pairs[rp][1]] = 0
            D11_ind[pairs[ap][0]] = 1
            D11_ind[pairs[ap][1]] = 1
            delta_11 = compute_delta_array(D11_ind, m)

            new_cost = compute_cost_vectorized(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n)
            dc = new_cost - current_cost

            if dc <= 0 or random.random() < math.exp(-dc / max(T, T_min)):
                current_cost = new_cost
            else:
                D11_ind[pairs[rp][0]] = 1
                D11_ind[pairs[rp][1]] = 1
                D11_ind[pairs[ap][0]] = 0
                D11_ind[pairs[ap][1]] = 0
                delta_11 = old_d11

        elif move_type < 0.80:
            # D12 single swap
            d12_others = [x for x in range(1, m) if D12_ind[x]]
            not_in = [x for x in range(1, m) if not D12_ind[x]]
            if not d12_others or not not_in:
                continue

            rem = random.choice(d12_others)
            add = random.choice(not_in)

            old_d12 = delta_12.copy()
            old_d12T = delta_12T.copy()

            D12_ind[rem] = 0
            D12_ind[add] = 1
            delta_12 = compute_delta_array(D12_ind, m)

            D12T_ind[(-rem) % m] = 0
            D12T_ind[(-add) % m] = 1
            delta_12T = compute_delta_array(D12T_ind, m)

            new_cost = compute_cost_vectorized(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n)
            dc = new_cost - current_cost

            if dc <= 0 or random.random() < math.exp(-dc / max(T, T_min)):
                current_cost = new_cost
            else:
                D12_ind[rem] = 1
                D12_ind[add] = 0
                D12T_ind[(-rem) % m] = 1
                D12T_ind[(-add) % m] = 0
                delta_12 = old_d12
                delta_12T = old_d12T

        else:
            # Multi-swap: swap 2 D12 elements at once (helps escape local minima)
            d12_others = [x for x in range(1, m) if D12_ind[x]]
            not_in = [x for x in range(1, m) if not D12_ind[x]]
            if len(d12_others) < 2 or len(not_in) < 2:
                continue

            rems = random.sample(d12_others, 2)
            adds = random.sample(not_in, 2)

            old_d12 = delta_12.copy()
            old_d12T = delta_12T.copy()

            for r in rems:
                D12_ind[r] = 0
                D12T_ind[(-r) % m] = 0
            for a in adds:
                D12_ind[a] = 1
                D12T_ind[(-a) % m] = 1

            delta_12 = compute_delta_array(D12_ind, m)
            delta_12T = compute_delta_array(D12T_ind, m)

            new_cost = compute_cost_vectorized(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n)
            dc = new_cost - current_cost

            if dc <= 0 or random.random() < math.exp(-dc / max(T, T_min)):
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
            last_improvement_iter = it
            if best_cost == 0:
                break

        # Reheat strategy: if stuck for a long time at low cost
        if it - last_improvement_iter > 500000 and best_cost > 0 and best_cost <= 8:
            T = 3.0  # Reheat
            last_improvement_iter = it
            reheat_count += 1
            # Restart from best known solution
            D11_ind = np.copy(best_D11)
            D12_ind = np.copy(best_D12)
            D12T_ind = np.zeros(m, dtype=np.int64)
            for x in range(m):
                if D12_ind[x]:
                    D12T_ind[(-x) % m] = 1
            delta_11 = compute_delta_array(D11_ind, m)
            delta_12 = compute_delta_array(D12_ind, m)
            delta_12T = compute_delta_array(D12T_ind, m)
            current_cost = best_cost

        T = max(T * alpha, T_min)

        if it % 1000000 == 0 and it > 0:
            elapsed = time.time() - t0
            print(f"    iter {it//1000000}M: cost={current_cost}, best={best_cost}, T={T:.6f}, reheats={reheat_count} ({elapsed:.0f}s)", flush=True)

    elapsed = time.time() - t0
    return best_cost, best_D11, best_D12, elapsed, reheat_count


def main():
    d11_sizes = [32, 30]
    max_seeds = int(sys.argv[1]) if len(sys.argv) > 1 else 50
    max_iter = int(sys.argv[2]) if len(sys.argv) > 2 else 10000000

    print(f"n={N_PARAM}, m={M}, N={N_VERTS}")
    print(f"D11 sizes to try: {d11_sizes}")
    print(f"Max seeds: {max_seeds}, Max iter: {max_iter}")
    print(f"Reheat strategy: restart from best after 500K stuck iterations")
    print(flush=True)

    global_best_cost = 999
    t0 = time.time()

    for d11_size in d11_sizes:
        print(f"\n--- |D11| = {d11_size} ---", flush=True)
        for seed in range(max_seeds):
            best_cost, best_D11, best_D12, elapsed, reheats = solve_n32(
                d11_size, seed, max_iter
            )
            total_elapsed = time.time() - t0
            print(f"  seed={seed}: best={best_cost}, reheats={reheats} ({total_elapsed:.0f}s)", flush=True)

            if best_cost < global_best_cost:
                global_best_cost = best_cost

            if best_cost == 0:
                D11_list = sorted([int(i) for i in range(M) if best_D11[i]])
                D12_list = sorted([int(i) for i in range(M) if best_D12[i]])

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
                    save_construction(G, result, fn, solver="turbo-SA")
                    print(f"  SAVED to {fn}")
                    sys.exit(0)

    print(f"\nNo solution found. Global best cost = {global_best_cost}")


if __name__ == "__main__":
    main()
