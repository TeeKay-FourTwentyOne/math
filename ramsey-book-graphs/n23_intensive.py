"""
Intensive SA search for n=23 (m=45) with longer runs and more restarts.

Key insight from exhaustive search: best budget = 9 for both |D11|=22 and |D11|=24.
SA gets very close (cost=2) but never hits 0. Need longer/smarter search.

Strategy:
1. Try |D11|=24 candidates (only 60, more special structure)
2. Use longer SA runs with adaptive cooling
3. Try joint optimization of D11 and D12 simultaneously
4. Also try non-standard sizes (|D11| = 20 or 26)
"""

import sys, os, math, random, json
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ramsey_core import Delta, BlockCirculantGraph, verify_construction, save_construction

m = 45
n = 23

# Compute symmetric pairs
pairs = []
for x in range(1, m):
    neg_x = (-x) % m
    if x <= neg_x:
        pairs.append((x, neg_x))


def full_cost(D11_set, D12_set):
    """Compute total violation cost using full verification."""
    D22_set = set(range(1, m)) - D11_set
    G = BlockCirculantGraph(n=n, D11=D11_set, D12=D12_set, D22=D22_set)
    result = verify_construction(G)
    return sum(excess for _, _, excess in result.violations), result


def joint_sa(d11_size, d12_size=22, max_iter=2000000, seed=0, verbose=True):
    """Joint SA optimizing both D11 and D12 simultaneously."""
    random.seed(seed)

    # Initialize D11 from random pairs
    num_pairs = d11_size // 2
    selected_pairs = random.sample(range(22), num_pairs)
    D11_set = set()
    for i in selected_pairs:
        D11_set.add(pairs[i][0])
        D11_set.add(pairs[i][1])

    # Initialize D12: 0 + random 21 elements
    D12_rest = random.sample(range(1, m), d12_size - 1)
    D12_set = set([0] + D12_rest)

    current_cost, _ = full_cost(D11_set, D12_set)
    best_cost = current_cost
    best_D11 = sorted(D11_set)
    best_D12 = sorted(D12_set)
    best_result = None

    T = 10.0
    T_min = 0.0001
    # Slower cooling for longer exploration
    alpha = 1 - 3.0 / max_iter

    stale_count = 0
    last_best = best_cost

    for it in range(max_iter):
        # Decide whether to modify D11 or D12
        if random.random() < 0.3:
            # Modify D11: swap a pair in/out
            in_pairs = [i for i in range(22)
                        if pairs[i][0] in D11_set]
            out_pairs = [i for i in range(22)
                         if pairs[i][0] not in D11_set]
            if not in_pairs or not out_pairs:
                continue
            remove_pair = random.choice(in_pairs)
            add_pair = random.choice(out_pairs)

            # Apply swap
            D11_set.discard(pairs[remove_pair][0])
            D11_set.discard(pairs[remove_pair][1])
            D11_set.add(pairs[add_pair][0])
            D11_set.add(pairs[add_pair][1])

            new_cost, result = full_cost(D11_set, D12_set)
            delta = new_cost - current_cost

            if delta <= 0 or random.random() < math.exp(-delta / max(T, 0.001)):
                current_cost = new_cost
            else:
                # Undo
                D11_set.discard(pairs[add_pair][0])
                D11_set.discard(pairs[add_pair][1])
                D11_set.add(pairs[remove_pair][0])
                D11_set.add(pairs[remove_pair][1])
                result = None
        else:
            # Modify D12: swap an element in/out (keep 0)
            D12_others = [x for x in D12_set if x != 0]
            if not D12_others:
                continue
            remove = random.choice(D12_others)
            not_in = list(set(range(1, m)) - D12_set)
            if not not_in:
                continue
            add = random.choice(not_in)

            D12_set.remove(remove)
            D12_set.add(add)

            new_cost, result = full_cost(D11_set, D12_set)
            delta = new_cost - current_cost

            if delta <= 0 or random.random() < math.exp(-delta / max(T, 0.001)):
                current_cost = new_cost
            else:
                D12_set.remove(add)
                D12_set.add(remove)
                result = None

        if current_cost < best_cost:
            best_cost = current_cost
            best_D11 = sorted(D11_set)
            best_D12 = sorted(D12_set)
            best_result = result
            if best_cost == 0:
                if verbose:
                    print(f"  SOLUTION FOUND at iter {it}!")
                return best_D11, best_D12, best_cost, best_result

        T *= alpha

        # Reheat if stuck
        stale_count += 1
        if stale_count > 50000:
            if best_cost != last_best:
                stale_count = 0
                last_best = best_cost
            else:
                T = max(T, 2.0)
                stale_count = 0

        if verbose and it % 200000 == 0:
            print(f"    iter {it}: cost={current_cost}, best={best_cost}, T={T:.5f}")

    return best_D11, best_D12, best_cost, best_result


# =============================================================================
# APPROACH 1: Joint SA with |D11|=22
# =============================================================================
print("=" * 80)
print("APPROACH 1: Joint SA with |D11|=22")
print("=" * 80)

for seed in range(10):
    print(f"\n--- Run seed={seed} ---")
    d11, d12, cost, result = joint_sa(22, 22, max_iter=1000000, seed=seed)
    print(f"  Best cost = {cost}")
    if cost == 0:
        print(f"  D11 = {d11}")
        print(f"  D12 = {d12}")
        print(f"  max_red={result.max_red_common}, max_blue={result.max_blue_common}")
        D11_set = set(d11)
        D12_set = set(d12)
        D22_set = set(range(1, m)) - D11_set
        G = BlockCirculantGraph(n=n, D11=D11_set, D12=D12_set, D22=D22_set)
        filename = "/Users/stephenpadgett/Projects/math/ramsey-book-graphs/solution_n23.json"
        save_construction(G, result, filename, solver="joint-SA-n23")
        print(f"  *** SAVED to {filename} ***")
        break
    else:
        print(f"  Best D11 = {d11[:8]}...")
        print(f"  Best D12 = {d12[:8]}...")

# =============================================================================
# APPROACH 2: Joint SA with |D11|=24
# =============================================================================
print("\n" + "=" * 80)
print("APPROACH 2: Joint SA with |D11|=24")
print("=" * 80)

for seed in range(10):
    print(f"\n--- Run seed={seed} ---")
    d11, d12, cost, result = joint_sa(24, 22, max_iter=1000000, seed=100+seed)
    print(f"  Best cost = {cost}")
    if cost == 0:
        print(f"  D11 = {d11}")
        print(f"  D12 = {d12}")
        print(f"  max_red={result.max_red_common}, max_blue={result.max_blue_common}")
        D11_set = set(d11)
        D12_set = set(d12)
        D22_set = set(range(1, m)) - D11_set
        G = BlockCirculantGraph(n=n, D11=D11_set, D12=D12_set, D22=D22_set)
        filename = "/Users/stephenpadgett/Projects/math/ramsey-book-graphs/solution_n23.json"
        save_construction(G, result, filename, solver="joint-SA-n23")
        print(f"  *** SAVED to {filename} ***")
        break
    else:
        print(f"  Best D11 = {d11[:8]}...")
        print(f"  Best D12 = {d12[:8]}...")

# =============================================================================
# APPROACH 3: Try |D11|=20 (10 pairs, smaller D11)
# =============================================================================
print("\n" + "=" * 80)
print("APPROACH 3: Joint SA with |D11|=20")
print("=" * 80)

for seed in range(5):
    print(f"\n--- Run seed={seed} ---")
    d11, d12, cost, result = joint_sa(20, 22, max_iter=1000000, seed=200+seed)
    print(f"  Best cost = {cost}")
    if cost == 0:
        print(f"  D11 = {d11}")
        print(f"  D12 = {d12}")
        D11_set = set(d11)
        D12_set = set(d12)
        D22_set = set(range(1, m)) - D11_set
        G = BlockCirculantGraph(n=n, D11=D11_set, D12=D12_set, D22=D22_set)
        filename = "/Users/stephenpadgett/Projects/math/ramsey-book-graphs/solution_n23.json"
        save_construction(G, result, filename, solver="joint-SA-n23")
        print(f"  *** SAVED to {filename} ***")
        break

# =============================================================================
# APPROACH 4: Try non-standard |D12| values
# =============================================================================
print("\n" + "=" * 80)
print("APPROACH 4: Try |D12|=21 (0 not in D12)")
print("=" * 80)

def joint_sa_no_zero(d11_size, d12_size=21, max_iter=1000000, seed=0, verbose=True):
    """SA where D12 does NOT contain 0."""
    random.seed(seed)

    num_pairs = d11_size // 2
    selected_pairs = random.sample(range(22), num_pairs)
    D11_set = set()
    for i in selected_pairs:
        D11_set.add(pairs[i][0])
        D11_set.add(pairs[i][1])

    # D12 without 0
    D12_rest = random.sample(range(1, m), d12_size)
    D12_set = set(D12_rest)

    current_cost, _ = full_cost(D11_set, D12_set)
    best_cost = current_cost
    best_D11 = sorted(D11_set)
    best_D12 = sorted(D12_set)
    best_result = None

    T = 10.0
    alpha = 1 - 3.0 / max_iter

    for it in range(max_iter):
        if random.random() < 0.3:
            in_pairs = [i for i in range(22) if pairs[i][0] in D11_set]
            out_pairs = [i for i in range(22) if pairs[i][0] not in D11_set]
            if not in_pairs or not out_pairs:
                continue
            remove_pair = random.choice(in_pairs)
            add_pair = random.choice(out_pairs)
            D11_set.discard(pairs[remove_pair][0])
            D11_set.discard(pairs[remove_pair][1])
            D11_set.add(pairs[add_pair][0])
            D11_set.add(pairs[add_pair][1])
            new_cost, result = full_cost(D11_set, D12_set)
            delta = new_cost - current_cost
            if delta <= 0 or random.random() < math.exp(-delta / max(T, 0.001)):
                current_cost = new_cost
            else:
                D11_set.discard(pairs[add_pair][0])
                D11_set.discard(pairs[add_pair][1])
                D11_set.add(pairs[remove_pair][0])
                D11_set.add(pairs[remove_pair][1])
                result = None
        else:
            D12_list = list(D12_set)
            remove = random.choice(D12_list)
            candidates = list(set(range(m)) - D12_set)
            add = random.choice(candidates)
            D12_set.remove(remove)
            D12_set.add(add)
            new_cost, result = full_cost(D11_set, D12_set)
            delta = new_cost - current_cost
            if delta <= 0 or random.random() < math.exp(-delta / max(T, 0.001)):
                current_cost = new_cost
            else:
                D12_set.remove(add)
                D12_set.add(remove)
                result = None

        if current_cost < best_cost:
            best_cost = current_cost
            best_D11 = sorted(D11_set)
            best_D12 = sorted(D12_set)
            best_result = result
            if best_cost == 0:
                if verbose:
                    print(f"  SOLUTION FOUND at iter {it}!")
                return best_D11, best_D12, best_cost, best_result

        T *= alpha
        if verbose and it % 200000 == 0:
            print(f"    iter {it}: cost={current_cost}, best={best_cost}, T={T:.5f}")

    return best_D11, best_D12, best_cost, best_result


# Try with D11=D12=QR-like pattern (like Paley) but adapted
for d11_size in [22, 24]:
    for d12_size in [21, 22, 23]:
        print(f"\n--- |D11|={d11_size}, |D12|={d12_size}, 0 not forced in D12 ---")
        d11, d12, cost, result = joint_sa_no_zero(d11_size, d12_size, max_iter=500000, seed=300)
        print(f"  Best cost = {cost}")
        if cost == 0:
            print(f"  D11 = {d11}")
            print(f"  D12 = {d12}")
            D11_set = set(d11)
            D12_set = set(d12)
            D22_set = set(range(1, m)) - D11_set
            G = BlockCirculantGraph(n=n, D11=D11_set, D12=D12_set, D22=D22_set)
            filename = "/Users/stephenpadgett/Projects/math/ramsey-book-graphs/solution_n23.json"
            save_construction(G, result, filename, solver="joint-SA-n23-free")
            print(f"  *** SAVED to {filename} ***")

print("\n" + "=" * 80)
print("SEARCH COMPLETE")
print("=" * 80)
