"""
Fast SA for n=24 (m=47, prime ≡ 3 mod 4).

m=47 is prime, 47 ≡ 3 mod 4, so -1 is a QNR.
This means QR is NOT symmetric under negation: -QR = QNR.

For primes ≡ 3 mod 4, the known constructions (n=6,8,10,...,22) use:
  |D11| = n-1 (= m-1)/2 = 23 for n=24)
  |D12| = n-1 = 23
  0 in D12
  D22 = complement(D11) in {1,...,m-1}, so |D22| = (m-1) - |D11|

D11 symmetric with |D11| = (m-1)/2 = 23 is the standard size.
For m ≡ 3 mod 4 with m prime: since m-1 = 46 is even, (m-1)/2 = 23 is odd.
Symmetric pairs: 23 pairs, so |D11| must be even. But 23 is odd!
Wait: m=47, pairs = {x, 47-x} for x=1..23 → 23 pairs.
For |D11| to be symmetric and have 23 elements, we'd need 11 full pairs + 1 self-paired.
But no x satisfies x ≡ -x mod 47 (would need 2x ≡ 0, x=0 which is excluded).
So |D11| must be even: 22 or 24.

For n=24: red_thresh = 22, blue_thresh = 23.
Try |D11| = 22 and |D11| = 24.
"""

import sys, os, math, random, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ramsey_core import Delta, BlockCirculantGraph, verify_construction, save_construction

m = 47
n = 24
N = 2 * m  # = 94


def compute_all_deltas(S_set):
    """Compute Delta(S, S, d) for all d in {0,...,m-1}."""
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


def fast_cost(D11_set, D12_set, delta_11, delta_12, delta_12T):
    """Compute violation cost from precomputed Delta arrays."""
    d11_size = len(D11_set)
    d12_size = len(D12_set)
    d22_size = m - 1 - d11_size
    d1 = d11_size + d12_size
    d2 = d22_size + d12_size

    red_thresh = n - 2  # 22
    blue_thresh = n - 1  # 23

    cost = 0

    # V1V1
    for d in range(1, m):
        common = delta_11[d] + delta_12[d]
        if d in D11_set:
            if common > red_thresh:
                cost += common - red_thresh
        else:
            blue_common = (N - 2) - 2 * d1 + common
            if blue_common > blue_thresh:
                cost += blue_common - blue_thresh

    # V2V2
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

    # V1V2: automatic when |D12| = n-1 and D22 = complement(D11)
    # V1V2 red: |D12| - 1 should be <= red_thresh = 22 → |D12|-1 = 22 ✓
    # V1V2 blue: (N-2) - d1 - d2 + |D12| should be <= blue_thresh
    v12_blue = (N - 2) - d1 - d2 + d12_size
    if v12_blue > blue_thresh:
        cost += (v12_blue - blue_thresh) * m  # applies to all blue V1V2 edges

    return cost


def update_deltas(D12_set):
    """Recompute D12 and D12T deltas."""
    delta_12 = compute_all_deltas(D12_set)
    D12T_set = {(-x) % m for x in D12_set}
    delta_12T = compute_all_deltas(D12T_set)
    return delta_12, delta_12T


# Compute symmetric pairs
pairs = []
for x in range(1, m):
    neg_x = (-x) % m
    if x <= neg_x:
        pairs.append((x, neg_x))

print(f"Symmetric pairs in Z_{m}: {len(pairs)}")


def fast_joint_sa(d11_size, d12_size=23, max_iter=5000000, seed=0, verbose=True):
    """Fast joint SA optimizing D11 and D12."""
    random.seed(seed)

    num_pairs = d11_size // 2
    selected_pairs = random.sample(range(len(pairs)), num_pairs)
    D11_set = set()
    for i in selected_pairs:
        D11_set.add(pairs[i][0])
        D11_set.add(pairs[i][1])

    # Initialize D12
    D12_rest = random.sample(range(1, m), d12_size - 1)
    D12_set = set([0] + D12_rest)

    delta_11 = compute_all_deltas(D11_set)
    delta_12, delta_12T = update_deltas(D12_set)

    current_cost = fast_cost(D11_set, D12_set, delta_11, delta_12, delta_12T)
    best_cost = current_cost
    best_D11 = sorted(D11_set)
    best_D12 = sorted(D12_set)

    T = 8.0
    alpha = 1 - 5.0 / max_iter

    for it in range(max_iter):
        if random.random() < 0.25:
            # Swap a pair in D11
            in_pairs = [i for i in range(len(pairs)) if pairs[i][0] in D11_set]
            out_pairs = [i for i in range(len(pairs)) if pairs[i][0] not in D11_set]
            if not in_pairs or not out_pairs:
                continue

            remove_pair = random.choice(in_pairs)
            add_pair = random.choice(out_pairs)

            D11_set.discard(pairs[remove_pair][0])
            D11_set.discard(pairs[remove_pair][1])
            D11_set.add(pairs[add_pair][0])
            D11_set.add(pairs[add_pair][1])

            old_delta_11 = delta_11
            delta_11 = compute_all_deltas(D11_set)

            new_cost = fast_cost(D11_set, D12_set, delta_11, delta_12, delta_12T)
            delta_c = new_cost - current_cost

            if delta_c <= 0 or random.random() < math.exp(-delta_c / max(T, 0.001)):
                current_cost = new_cost
            else:
                D11_set.discard(pairs[add_pair][0])
                D11_set.discard(pairs[add_pair][1])
                D11_set.add(pairs[remove_pair][0])
                D11_set.add(pairs[remove_pair][1])
                delta_11 = old_delta_11
        else:
            # Swap an element in D12 (keep 0)
            D12_others = [x for x in D12_set if x != 0]
            remove = random.choice(D12_others)
            not_in = list(set(range(1, m)) - D12_set)
            add = random.choice(not_in)

            D12_set.remove(remove)
            D12_set.add(add)

            old_delta_12 = delta_12
            old_delta_12T = delta_12T
            delta_12, delta_12T = update_deltas(D12_set)

            new_cost = fast_cost(D11_set, D12_set, delta_11, delta_12, delta_12T)
            delta_c = new_cost - current_cost

            if delta_c <= 0 or random.random() < math.exp(-delta_c / max(T, 0.001)):
                current_cost = new_cost
            else:
                D12_set.remove(add)
                D12_set.add(remove)
                delta_12 = old_delta_12
                delta_12T = old_delta_12T

        if current_cost < best_cost:
            best_cost = current_cost
            best_D11 = sorted(D11_set)
            best_D12 = sorted(D12_set)
            if best_cost == 0:
                if verbose:
                    print(f"  SOLUTION FOUND at iter {it}!")
                return best_D11, best_D12, best_cost

        T = max(T * alpha, 0.0001)

        if verbose and it % 500000 == 0:
            print(f"    iter {it}: cost={current_cost}, best={best_cost}, T={T:.6f}")

    return best_D11, best_D12, best_cost


print("=" * 80)
print(f"FAST JOINT SA FOR n={n} (m={m})")
print("=" * 80)

start_time = time.time()

found = False
for d11_size in [22, 24]:
    if found:
        break
    print(f"\n--- |D11| = {d11_size} ---")
    for seed in range(20):
        d11, d12, cost = fast_joint_sa(d11_size, d12_size=23, max_iter=3000000, seed=seed)
        elapsed = time.time() - start_time
        print(f"  seed={seed}: best cost = {cost} (elapsed: {elapsed:.1f}s)")
        if cost == 0:
            print(f"\n  *** SOLUTION FOUND! ***")
            print(f"  D11 = {d11}")
            print(f"  D12 = {d12}")

            D11_set = set(d11)
            D12_set = set(d12)
            D22_set = set(range(1, m)) - D11_set
            G = BlockCirculantGraph(n=n, D11=D11_set, D12=D12_set, D22=D22_set)
            result = verify_construction(G)
            print(f"  Full verification: valid={result.valid}")
            print(f"  max_red={result.max_red_common} (thresh {result.red_threshold})")
            print(f"  max_blue={result.max_blue_common} (thresh {result.blue_threshold})")

            if result.valid:
                filename = "/Users/stephenpadgett/Projects/math/ramsey-book-graphs/solution_n24.json"
                save_construction(G, result, filename, solver="fast-joint-SA")
                print(f"  SAVED to {filename}")
                found = True
                break

total_time = time.time() - start_time
print(f"\nTotal time: {total_time:.1f}s")
print("DONE")
