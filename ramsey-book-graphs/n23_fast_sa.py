"""
Fast SA for n=23 (m=45) with incremental Delta updates.

The key optimization: when swapping one element in D12, we can update
Delta(D12, D12, d) incrementally instead of recomputing from scratch.

When we remove element 'r' and add element 'a':
  For each d: Delta_new(d) = Delta_old(d)
    - [r-d in D12_old] - [d in D12_old and r = (d+something) in old]  (removing r)
    + [a-d in D12_new] + [d in D12_new and a = ...] (adding a)

More precisely:
  Delta(S, S, d) = |{s in S : s-d in S}|
  When S changes from S to S' = S - {r} + {a}:
    delta_change(d) = -[(r-d) % m in S'] - [(r) in S' and (r-d = something)]
    Actually, it's easier to just recount for the affected elements.

Simpler approach: maintain a characteristic vector and use convolution.
Delta(S, S, d) = sum_{x in Z_m} 1_S(x) * 1_S(x-d)

When we change D12 by removing r and adding a, the change to Delta for each d is:
  Delta_new(d) - Delta_old(d) =
    1_{new}(a) * 1_{new}(a-d) + 1_{new}(a+d) * 1_{new}(a) [adding a: a contributes as "x" and as "x-d"]
    - 1_{old}(r) * 1_{old}(r-d) - 1_{old}(r+d) * 1_{old}(r) [removing r]

Wait, let me be more careful. Let D12_new = D12_old - {r} + {a}.

Delta(D12_new, D12_new, d) - Delta(D12_old, D12_old, d)
  = sum_x [1_new(x)*1_new(x-d) - 1_old(x)*1_old(x-d)]

The only x values where the indicator changes are x=a (added), x=r (removed),
x-d=a (i.e., x=a+d), and x-d=r (i.e., x=r+d).

So we need to check x in {a, r, a+d, r+d} (all mod m).

This gives O(1) update per d, and O(m) total per swap. Much faster than O(m^2).
"""

import sys, os, math, random, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ramsey_core import Delta, BlockCirculantGraph, verify_construction, save_construction

m = 45
n = 23
N = 2 * m  # = 90


def compute_all_deltas(S_set):
    """Compute Delta(S, S, d) for all d in {1,...,m-1}."""
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

    red_thresh = n - 2  # 21
    blue_thresh = n - 1  # 22

    cost = 0

    # V1V1: common = delta_11[d] + delta_12[d]
    for d in range(1, m):
        common = delta_11[d] + delta_12[d]
        if d in D11_set:
            # Red edge
            if common > red_thresh:
                cost += common - red_thresh
        else:
            # Blue edge
            blue_common = (N - 2) - 2 * d1 + common
            if blue_common > blue_thresh:
                cost += blue_common - blue_thresh

    # V2V2: need Delta(D22,D22,d) + Delta(D12T,D12T,d)
    # Delta(D22,D22,d) = Delta(D11,D11,d) + (m-2-2*d11_size) + 2*[d in D11]
    v22_const = m - 2 - 2 * d11_size
    for d in range(1, m):
        d22_delta = delta_11[d] + v22_const + (2 if d in D11_set else 0)
        common = d22_delta + delta_12T[d]
        if d not in D11_set:
            # d in D22 = red edge
            if common > red_thresh:
                cost += common - red_thresh
        else:
            # d in D11 -> blue edge in V2V2
            blue_common = (N - 2) - 2 * d2 + common
            if blue_common > blue_thresh:
                cost += blue_common - blue_thresh

    # V1V2: by theorem, V1V2(d) = |D12| - [d in D12], automatically OK
    # if |D12| = n-1 = 22 and red_thresh = 21
    # V1V2 red: |D12| - 1 = 21 = red_thresh (OK)
    # V1V2 blue: (N-2) - d1 - d2 + |D12| = 88 - d1 - d2 + d12_size
    #          = 88 - (d11_size + d12_size) - (d22_size + d12_size) + d12_size
    #          = 88 - d11_size - d22_size - d12_size
    #          = 88 - (m-1) - d12_size = 88 - 44 - d12_size = 44 - d12_size
    # For d12_size = 22: blue_common = 22 = blue_thresh (OK)
    # So V1V2 is automatically satisfied when |D12| = 22 and D22=complement(D11)

    return cost


def update_delta_12_swap(delta_12, delta_12T, D12_set, remove, add):
    """Update delta_12 and delta_12T when swapping remove->add in D12_set.
    D12_set should already have the swap applied."""
    # D12_set is the NEW set (already has 'add', doesn't have 'remove')

    for d in range(1, m):
        # Update delta_12[d] = Delta(D12, D12, d) = |{x in D12 : (x-d) in D12}|
        # Removing r: lose contributions where x=r or (x-d)=r
        # Adding a: gain contributions where x=a or (x-d)=a

        change = 0

        # x = remove was removed:
        # Old contribution: 1_{old}(remove) * 1_{old}((remove-d)%m) = 1 * [((remove-d)%m) in D12_old]
        # But remove is NOT in D12_new. Was (remove-d)%m in D12_old?
        rd = (remove - d) % m
        if rd != remove:  # if rd == remove, it was in old but now removed
            if rd in D12_set:  # rd is still in new set
                change -= 1
        # else rd == remove: it was in old (so 1*1=1) but now it's removed; but since
        # remove is gone, this pair doesn't exist in new. Handle below.

        # x-d = remove, i.e., x = (remove+d)%m was a pair contributor:
        rpd = (remove + d) % m
        if rpd != remove:  # if rpd == remove, already counted
            if rpd in D12_set:  # rpd still in new set
                change -= 1
        # else rpd == remove: same element, skip

        # Special case: if rd == remove, the pair (remove, remove-d) was counted once
        # and we already handled it in the first block... actually let me re-derive.

        # Actually, the clean way: Delta changes only at positions involving remove or add.
        # Let me just recount the affected terms.

        # In D12_old: count pairs (x, x-d) where both in D12_old
        # In D12_new: count pairs (x, x-d) where both in D12_new
        # D12_new = D12_old - {remove} + {add}

        # Pairs lost (in old but not new):
        # Any pair involving 'remove': either x=remove or x-d=remove
        # Pair (remove, (remove-d)%m): was counted if (remove-d)%m in D12_old
        #   = 1 if (remove-d)%m in (D12_new + {remove} - {add}) = D12_set ∪ {remove} - {add}
        # Pair ((remove+d)%m, remove): counted if (remove+d)%m in D12_old
        # But these might double-count if d=0 (not possible since d>=1)
        # or if (remove+d)%m == remove (d=0, not possible)

        # Let me just use the simple approach: recompute delta for this d directly.
        pass

    # Actually the incremental approach is error-prone. Let me use a hybrid:
    # Precompute indicator arrays and use them directly.
    # With m=45, recomputing all deltas is O(m^2) = O(2025) per swap, which is fine.

    indicator = [0] * m
    for x in D12_set:
        indicator[x % m] = 1

    for d in range(1, m):
        count = 0
        for x in range(m):
            if indicator[x] and indicator[(x - d) % m]:
                count += 1
        delta_12[d] = count

    # D12T = {(-x)%m : x in D12}
    indicator_T = [0] * m
    for x in D12_set:
        indicator_T[(-x) % m] = 1

    for d in range(1, m):
        count = 0
        for x in range(m):
            if indicator_T[x] and indicator_T[(x - d) % m]:
                count += 1
        delta_12T[d] = count


# =============================================================================
# Compute symmetric pairs
# =============================================================================
pairs = []
for x in range(1, m):
    neg_x = (-x) % m
    if x <= neg_x:
        pairs.append((x, neg_x))


# =============================================================================
# SA with fast cost function
# =============================================================================
def fast_joint_sa(d11_size, max_iter=5000000, seed=0, verbose=True):
    """Fast joint SA optimizing D11 and D12."""
    random.seed(seed)

    # Initialize D11 from random pairs
    num_pairs = d11_size // 2
    selected_pairs = random.sample(range(22), num_pairs)
    D11_set = set()
    for i in selected_pairs:
        D11_set.add(pairs[i][0])
        D11_set.add(pairs[i][1])

    # Initialize D12: 0 + random 21 elements
    D12_rest = random.sample(range(1, m), 21)
    D12_set = set([0] + D12_rest)

    # Precompute Delta arrays
    delta_11 = compute_all_deltas(D11_set)
    delta_12 = compute_all_deltas(D12_set)
    D12T_set = {(-x) % m for x in D12_set}
    delta_12T = compute_all_deltas(D12T_set)

    current_cost = fast_cost(D11_set, D12_set, delta_11, delta_12, delta_12T)
    best_cost = current_cost
    best_D11 = sorted(D11_set)
    best_D12 = sorted(D12_set)

    T = 8.0
    T_min = 0.0001
    alpha = 1 - 5.0 / max_iter

    for it in range(max_iter):
        if random.random() < 0.25:
            # Swap a pair in D11
            in_pairs = [i for i in range(22) if pairs[i][0] in D11_set]
            out_pairs = [i for i in range(22) if pairs[i][0] not in D11_set]
            if not in_pairs or not out_pairs:
                continue

            remove_pair = random.choice(in_pairs)
            add_pair = random.choice(out_pairs)

            # Apply swap
            D11_set.discard(pairs[remove_pair][0])
            D11_set.discard(pairs[remove_pair][1])
            D11_set.add(pairs[add_pair][0])
            D11_set.add(pairs[add_pair][1])

            # Recompute delta_11 (D11 changed)
            old_delta_11 = delta_11[:]
            delta_11 = compute_all_deltas(D11_set)

            new_cost = fast_cost(D11_set, D12_set, delta_11, delta_12, delta_12T)
            delta_c = new_cost - current_cost

            if delta_c <= 0 or random.random() < math.exp(-delta_c / max(T, 0.001)):
                current_cost = new_cost
            else:
                # Undo
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

            # Save old deltas
            old_delta_12 = delta_12[:]
            old_delta_12T = delta_12T[:]

            # Update D12T
            D12T_set.discard((-remove) % m)
            D12T_set.add((-add) % m)

            # Recompute deltas (O(m^2) ~ 2025 ops, fast)
            update_delta_12_swap(delta_12, delta_12T, D12_set, remove, add)

            new_cost = fast_cost(D11_set, D12_set, delta_11, delta_12, delta_12T)
            delta_c = new_cost - current_cost

            if delta_c <= 0 or random.random() < math.exp(-delta_c / max(T, 0.001)):
                current_cost = new_cost
            else:
                # Undo
                D12_set.remove(add)
                D12_set.add(remove)
                D12T_set.discard((-add) % m)
                D12T_set.add((-remove) % m)
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

        T = max(T * alpha, T_min)

        if verbose and it % 500000 == 0:
            print(f"    iter {it}: cost={current_cost}, best={best_cost}, T={T:.6f}")

    return best_D11, best_D12, best_cost


# =============================================================================
# Run searches
# =============================================================================
print("=" * 80)
print("FAST JOINT SA FOR n=23 (m=45)")
print("=" * 80)

start_time = time.time()

found = False
for d11_size in [22, 24, 20]:
    if found:
        break
    print(f"\n--- |D11| = {d11_size} ---")
    for seed in range(20):
        d11, d12, cost = fast_joint_sa(d11_size, max_iter=3000000, seed=seed)
        elapsed = time.time() - start_time
        print(f"  seed={seed}: best cost = {cost} (elapsed: {elapsed:.1f}s)")
        if cost == 0:
            print(f"\n  *** SOLUTION FOUND! ***")
            print(f"  D11 = {d11}")
            print(f"  D12 = {d12}")

            # Full verification
            D11_set = set(d11)
            D12_set = set(d12)
            D22_set = set(range(1, m)) - D11_set
            G = BlockCirculantGraph(n=n, D11=D11_set, D12=D12_set, D22=D22_set)
            result = verify_construction(G)
            print(f"  Full verification: valid={result.valid}")
            print(f"  max_red={result.max_red_common} (thresh {result.red_threshold})")
            print(f"  max_blue={result.max_blue_common} (thresh {result.blue_threshold})")

            if result.valid:
                filename = "/Users/stephenpadgett/Projects/math/ramsey-book-graphs/solution_n23.json"
                save_construction(G, result, filename, solver="fast-joint-SA")
                print(f"  SAVED to {filename}")
                found = True
                break

total_time = time.time() - start_time
print(f"\nTotal time: {total_time:.1f}s")
print("DONE")
