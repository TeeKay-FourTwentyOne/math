"""Search for n=22 with restarts and slower cooling."""

import time
import random
import math
from ramsey_core import BlockCirculantGraph, verify_construction, compute_violation_cost

def quadratic_residues(m):
    return {(x * x) % m for x in range(1, m)} - {0}

n = 22
m = 43
half_m = m // 2

QR43 = quadratic_residues(43)

def get_D22_half(D11_half):
    return set(range(1, half_m + 1)) - D11_half

def make_graph(D11_half, D12):
    return BlockCirculantGraph(n=n, D11=D11_half, D12=D12, D22=get_D22_half(D11_half))

def get_cost(D11_half, D12):
    return compute_violation_cost(make_graph(D11_half, D12))

def random_move(D11_half, D12):
    D11_new = D11_half.copy()
    D12_new = D12.copy()

    if random.random() < 0.33:
        d = random.randint(1, half_m)
        if d in D11_new:
            D11_new.remove(d)
        else:
            D11_new.add(d)
    else:
        d = random.randint(0, m - 1)
        if d in D12_new:
            D12_new.remove(d)
        else:
            D12_new.add(d)

    return D11_new, D12_new

def sa_search(D11_init, D12_init, max_iter=200000, seed=None):
    """Run SA from given initial state."""
    if seed is not None:
        random.seed(seed)

    current_D11 = D11_init.copy()
    current_D12 = D12_init.copy()
    current_cost = get_cost(current_D11, current_D12)

    best_D11 = current_D11.copy()
    best_D12 = current_D12.copy()
    best_cost = current_cost

    temp = 50.0
    cooling = 0.99998

    for it in range(max_iter):
        if current_cost == 0:
            return best_D11, best_D12, best_cost, it

        new_D11, new_D12 = random_move(current_D11, current_D12)
        new_cost = get_cost(new_D11, new_D12)

        delta = new_cost - current_cost
        if delta < 0 or (temp > 0.001 and random.random() < math.exp(-delta / temp)):
            current_D11 = new_D11
            current_D12 = new_D12
            current_cost = new_cost

            if current_cost < best_cost:
                best_cost = current_cost
                best_D11 = current_D11.copy()
                best_D12 = current_D12.copy()

        temp *= cooling

    return best_D11, best_D12, best_cost, max_iter

# Initial seeds to try
D11_seeds = [
    {1, 4, 6, 9, 10, 11, 13, 14, 15, 16, 17},  # First 11 from QR43
    set(range(1, 12)),  # {1..11}
    {1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21},  # odd numbers
    {2, 4, 6, 8, 10, 12, 14, 16, 18, 20},  # even numbers (10 elements)
]

D12_seeds = [
    QR43.copy(),
    QR43 | {0},
    set(range(43)) - QR43,  # Non-residues
]

print("Searching for R(B_21, B_22) construction")
print(f"n={n}, m={m}, N={4*n-2}")
print("=" * 60)

global_best_cost = float('inf')
global_best_D11 = None
global_best_D12 = None
total_start = time.time()

# Try all combinations with multiple random seeds
for d11_idx, D11_init in enumerate(D11_seeds):
    for d12_idx, D12_init in enumerate(D12_seeds):
        init_cost = get_cost(D11_init, D12_init)
        print(f"\nConfig D11[{d11_idx}], D12[{d12_idx}]: initial cost = {init_cost}")

        for seed in range(3):  # 3 random restarts per config
            start = time.time()
            best_D11, best_D12, best_cost, iters = sa_search(
                D11_init, D12_init, max_iter=100000, seed=seed*1000 + d11_idx*100 + d12_idx
            )
            elapsed = time.time() - start

            if best_cost < global_best_cost:
                global_best_cost = best_cost
                global_best_D11 = best_D11
                global_best_D12 = best_D12
                print(f"  Seed {seed}: cost={best_cost} (NEW BEST!) in {elapsed:.1f}s")
            else:
                print(f"  Seed {seed}: cost={best_cost} in {elapsed:.1f}s")

            if best_cost == 0:
                break

        if global_best_cost == 0:
            break
    if global_best_cost == 0:
        break

# If not found, do intensive search from best found
if global_best_cost > 0:
    print(f"\nIntensive search from best (cost={global_best_cost})...")
    for seed in range(5):
        start = time.time()
        best_D11, best_D12, best_cost, _ = sa_search(
            global_best_D11, global_best_D12, max_iter=200000, seed=seed*7777
        )
        elapsed = time.time() - start

        if best_cost < global_best_cost:
            global_best_cost = best_cost
            global_best_D11 = best_D11
            global_best_D12 = best_D12
            print(f"  Attempt {seed}: cost={best_cost} (NEW BEST!) in {elapsed:.1f}s")
        else:
            print(f"  Attempt {seed}: cost={best_cost} in {elapsed:.1f}s")

        if best_cost == 0:
            break

total_elapsed = time.time() - total_start
print()
print("=" * 60)
print("FINAL RESULT")
print("=" * 60)
print(f"Best violation count: {global_best_cost}")
print(f"Total time: {total_elapsed:.1f}s")

G = make_graph(global_best_D11, global_best_D12)
print(f"\nD11 = {sorted(G.D11)}")
print(f"D12 = {sorted(G.D12)}")
print(f"D22 = {sorted(G.D22)}")

result = verify_construction(G)
if global_best_cost > 0:
    print(f"\nViolations ({len(result.violations)}):")
    for v in result.violations[:20]:
        print(f"  {v}")
else:
    print(f"\nVerification: valid={result.valid}")
    print(f"Max red common: {result.max_red_common} (threshold: {result.red_threshold})")
    print(f"Max blue common: {result.max_blue_common} (threshold: {result.blue_threshold})")
