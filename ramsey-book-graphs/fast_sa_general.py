"""
General fast SA solver for R(B_{n-1}, B_n) = 4n-1.

Usage: python fast_sa_general.py <n> [max_seeds] [max_iter]
"""

import sys, os, math, random, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ramsey_core import Delta, BlockCirculantGraph, verify_construction, save_construction


def compute_all_deltas(S_set, m):
    indicator = [0] * m
    for x in S_set:
        indicator[x % m] = 1
    deltas = [0] * m
    for d in range(1, m):
        count = 0
        for x in range(m):
            if indicator[x] and indicator[(x - d) % m]:
                count += 1
        deltas[d] = count
    return deltas


def fast_cost(D11_set, D12_set, delta_11, delta_12, delta_12T, n, m):
    d11_size = len(D11_set)
    d12_size = len(D12_set)
    d22_size = m - 1 - d11_size
    d1 = d11_size + d12_size
    d2 = d22_size + d12_size
    N = 2 * m

    red_thresh = n - 2
    blue_thresh = n - 1

    cost = 0

    for d in range(1, m):
        common = delta_11[d] + delta_12[d]
        if d in D11_set:
            if common > red_thresh:
                cost += common - red_thresh
        else:
            blue_common = (N - 2) - 2 * d1 + common
            if blue_common > blue_thresh:
                cost += blue_common - blue_thresh

    v22_const = m - 2 - 2 * d11_size
    for d in range(1, m):
        d22_delta = delta_11[d] + v22_const + (2 if d in D11_set else 0)
        common = d22_delta + delta_12T[d]
        if d not in D11_set:
            if common > red_thresh:
                cost += common - red_thresh
        else:
            blue_common = (N - 2) - 2 * d2 + common
            if blue_common > blue_thresh:
                cost += blue_common - blue_thresh

    # V1V2 check
    v12_blue = (N - 2) - d1 - d2 + d12_size
    if v12_blue > blue_thresh:
        cost += (v12_blue - blue_thresh) * m

    return cost


def solve_n(n_param, max_seeds=20, max_iter=5000000):
    m = 2 * n_param - 1
    N = 2 * m

    # Symmetric pairs
    pairs = []
    for x in range(1, m):
        neg_x = (-x) % m
        if x <= neg_x:
            pairs.append((x, neg_x))

    num_pairs = len(pairs)
    d12_size = n_param - 1

    # Determine valid D11 sizes (must be even for symmetry)
    # Standard: |D11| ≈ (m-1)/2
    half = (m - 1) // 2
    d11_sizes = []
    for s in [half - 1, half, half + 1, half - 2, half + 2]:
        if s > 0 and s < m - 1 and s % 2 == 0 and (s // 2) <= num_pairs:
            d11_sizes.append(s)
    # Remove duplicates and sort by closeness to half
    d11_sizes = sorted(set(d11_sizes), key=lambda s: abs(s - half))

    print(f"n={n_param}, m={m}, N={N}")
    print(f"Symmetric pairs: {num_pairs}")
    print(f"|D12| = {d12_size}")
    print(f"D11 sizes to try: {d11_sizes}")

    start_time = time.time()

    for d11_size in d11_sizes:
        print(f"\n--- |D11| = {d11_size} ---")
        for seed in range(max_seeds):
            random.seed(seed)

            num_selected = d11_size // 2
            selected_pairs = random.sample(range(num_pairs), num_selected)
            D11_set = set()
            for i in selected_pairs:
                D11_set.add(pairs[i][0])
                D11_set.add(pairs[i][1])

            D12_rest = random.sample(range(1, m), d12_size - 1)
            D12_set = set([0] + D12_rest)

            delta_11 = compute_all_deltas(D11_set, m)
            delta_12 = compute_all_deltas(D12_set, m)
            D12T_set = {(-x) % m for x in D12_set}
            delta_12T = compute_all_deltas(D12T_set, m)

            current_cost = fast_cost(D11_set, D12_set, delta_11, delta_12, delta_12T, n_param, m)
            best_cost = current_cost
            best_D11 = sorted(D11_set)
            best_D12 = sorted(D12_set)

            T = 8.0
            alpha = 1 - 5.0 / max_iter

            for it in range(max_iter):
                if random.random() < 0.25:
                    in_pairs = [i for i in range(num_pairs) if pairs[i][0] in D11_set]
                    out_pairs = [i for i in range(num_pairs) if pairs[i][0] not in D11_set]
                    if not in_pairs or not out_pairs:
                        continue

                    rp = random.choice(in_pairs)
                    ap = random.choice(out_pairs)
                    D11_set.discard(pairs[rp][0])
                    D11_set.discard(pairs[rp][1])
                    D11_set.add(pairs[ap][0])
                    D11_set.add(pairs[ap][1])

                    old_d11 = delta_11
                    delta_11 = compute_all_deltas(D11_set, m)

                    new_cost = fast_cost(D11_set, D12_set, delta_11, delta_12, delta_12T, n_param, m)
                    dc = new_cost - current_cost

                    if dc <= 0 or random.random() < math.exp(-dc / max(T, 0.001)):
                        current_cost = new_cost
                    else:
                        D11_set.discard(pairs[ap][0])
                        D11_set.discard(pairs[ap][1])
                        D11_set.add(pairs[rp][0])
                        D11_set.add(pairs[rp][1])
                        delta_11 = old_d11
                else:
                    D12_others = [x for x in D12_set if x != 0]
                    rem = random.choice(D12_others)
                    not_in = list(set(range(1, m)) - D12_set)
                    add = random.choice(not_in)

                    D12_set.remove(rem)
                    D12_set.add(add)

                    old_d12 = delta_12
                    old_d12T = delta_12T
                    delta_12 = compute_all_deltas(D12_set, m)
                    D12T_set = {(-x) % m for x in D12_set}
                    delta_12T = compute_all_deltas(D12T_set, m)

                    new_cost = fast_cost(D11_set, D12_set, delta_11, delta_12, delta_12T, n_param, m)
                    dc = new_cost - current_cost

                    if dc <= 0 or random.random() < math.exp(-dc / max(T, 0.001)):
                        current_cost = new_cost
                    else:
                        D12_set.remove(add)
                        D12_set.add(rem)
                        delta_12 = old_d12
                        delta_12T = old_d12T

                if current_cost < best_cost:
                    best_cost = current_cost
                    best_D11 = sorted(D11_set)
                    best_D12 = sorted(D12_set)
                    if best_cost == 0:
                        break

                T = max(T * alpha, 0.0001)

            elapsed = time.time() - start_time
            print(f"  seed={seed}: best cost = {best_cost} ({elapsed:.1f}s)")

            if best_cost == 0:
                print(f"\n  *** SOLUTION FOUND! ***")
                print(f"  D11 = {best_D11}")
                print(f"  D12 = {best_D12}")

                D11_s = set(best_D11)
                D12_s = set(best_D12)
                D22_s = set(range(1, m)) - D11_s
                G = BlockCirculantGraph(n=n_param, D11=D11_s, D12=D12_s, D22=D22_s)
                result = verify_construction(G)
                print(f"  Full verification: valid={result.valid}")
                print(f"  max_red={result.max_red_common} (thresh {result.red_threshold})")
                print(f"  max_blue={result.max_blue_common} (thresh {result.blue_threshold})")

                if result.valid:
                    fn = f"/Users/stephenpadgett/Projects/math/ramsey-book-graphs/solution_n{n_param}.json"
                    save_construction(G, result, fn, solver="fast-joint-SA")
                    print(f"  SAVED to {fn}")
                    return True

    return False


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python fast_sa_general.py <n> [max_seeds] [max_iter]")
        sys.exit(1)

    n_param = int(sys.argv[1])
    max_seeds = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    max_iter = int(sys.argv[3]) if len(sys.argv) > 3 else 5000000

    solved = solve_n(n_param, max_seeds, max_iter)
    if not solved:
        print(f"\nNo solution found for n={n_param}")
