"""
Search for n=22 with D22 as independent variable (not constrained to complement of D11).
"""

import time
import random
import math
from ramsey_core import BlockCirculantGraph, verify_construction, compute_violation_cost

n = 22
m = 43
half_m = m // 2

def make_graph(D11_half, D12, D22_half):
    """Create graph with independent D11, D12, D22."""
    return BlockCirculantGraph(n=n, D11=D11_half, D12=D12, D22=D22_half)

def get_cost(D11_half, D12, D22_half):
    return compute_violation_cost(make_graph(D11_half, D12, D22_half))

def random_move(D11_half, D12, D22_half):
    """Apply random move to one of the three sets."""
    D11_new = D11_half.copy()
    D12_new = D12.copy()
    D22_new = D22_half.copy()

    r = random.random()
    if r < 0.25:  # D11 move
        d = random.randint(1, half_m)
        if d in D11_new:
            D11_new.remove(d)
        else:
            D11_new.add(d)
    elif r < 0.6:  # D12 move (larger set, more moves)
        d = random.randint(0, m - 1)
        if d in D12_new:
            D12_new.remove(d)
        else:
            D12_new.add(d)
    else:  # D22 move
        d = random.randint(1, half_m)
        if d in D22_new:
            D22_new.remove(d)
        else:
            D22_new.add(d)

    return D11_new, D12_new, D22_new

# Best found configuration (cost 8) - convert to half representatives
# D11 full = [2, 3, 4, 7, 9, 11, 12, 15, 19, 21, 22, 24, 28, 31, 32, 34, 36, 39, 40, 41]
D11_best = {2, 3, 4, 7, 9, 11, 12, 15, 19, 21}  # half reps (min(d, 43-d))

# D12 = [6, 10, 11, 12, 13, 14, 15, 16, 17, 18, 21, 24, 26, 27, 28, 30, 33, 35, 37, 38, 39]
D12_best = {6, 10, 11, 12, 13, 14, 15, 16, 17, 18, 21, 24, 26, 27, 28, 30, 33, 35, 37, 38, 39}

# D22 full = [1, 5, 6, 8, 10, 13, 14, 16, 17, 18, 20, 23, 25, 26, 27, 29, 30, 33, 35, 37, 38, 42]
D22_best = {1, 5, 6, 8, 10, 13, 14, 16, 17, 18, 20}  # half reps

print("Search with RELAXED D22 constraint (D22 independent of D11)", flush=True)
print(f"n={n}, m={m}, N={4*n-2}", flush=True)
print(f"Starting D11_half: {sorted(D11_best)}", flush=True)
print(f"Starting D12: {sorted(D12_best)}", flush=True)
print(f"Starting D22_half: {sorted(D22_best)}", flush=True)

init_cost = get_cost(D11_best, D12_best, D22_best)
print(f"Initial cost: {init_cost}", flush=True)
print("=" * 60, flush=True)

global_best_cost = init_cost
global_best_D11 = D11_best.copy()
global_best_D12 = D12_best.copy()
global_best_D22 = D22_best.copy()

total_start = time.time()

for restart in range(40):
    random.seed(restart * 123456789)

    # Start from global best with small perturbation
    current_D11 = global_best_D11.copy()
    current_D12 = global_best_D12.copy()
    current_D22 = global_best_D22.copy()

    # Small perturbation
    for _ in range(restart % 5):
        current_D11, current_D12, current_D22 = random_move(current_D11, current_D12, current_D22)

    current_cost = get_cost(current_D11, current_D12, current_D22)

    best_D11 = current_D11.copy()
    best_D12 = current_D12.copy()
    best_D22 = current_D22.copy()
    best_cost = current_cost

    temp = 30.0
    cooling = 0.999985

    for it in range(200000):
        if current_cost == 0:
            break

        new_D11, new_D12, new_D22 = random_move(current_D11, current_D12, current_D22)
        new_cost = get_cost(new_D11, new_D12, new_D22)

        delta = new_cost - current_cost
        if delta < 0 or (temp > 0.0001 and random.random() < math.exp(-delta / temp)):
            current_D11 = new_D11
            current_D12 = new_D12
            current_D22 = new_D22
            current_cost = new_cost

            if current_cost < best_cost:
                best_cost = current_cost
                best_D11 = current_D11.copy()
                best_D12 = current_D12.copy()
                best_D22 = current_D22.copy()

        temp *= cooling

    if best_cost < global_best_cost:
        global_best_cost = best_cost
        global_best_D11 = best_D11
        global_best_D12 = best_D12
        global_best_D22 = best_D22
        elapsed = time.time() - total_start
        print(f"Restart {restart}: cost={best_cost} (NEW BEST!) [{elapsed:.1f}s]", flush=True)

        if global_best_cost == 0:
            break
    elif restart % 5 == 0:
        print(f"Restart {restart}: cost={best_cost}", flush=True)

total_elapsed = time.time() - total_start
print(flush=True)
print("=" * 60, flush=True)
print("FINAL RESULT", flush=True)
print("=" * 60, flush=True)
print(f"Best cost: {global_best_cost}", flush=True)
print(f"Total time: {total_elapsed:.1f}s", flush=True)

G = make_graph(global_best_D11, global_best_D12, global_best_D22)
print(f"D11 = {sorted(G.D11)}", flush=True)
print(f"D12 = {sorted(G.D12)}", flush=True)
print(f"D22 = {sorted(G.D22)}", flush=True)

# Check if D22 is still complement of D11
D11_full = set(G.D11)
D22_full = set(G.D22)
complement_D11 = set(range(1, m)) - D11_full
is_complement = (D22_full == complement_D11)
print(f"D22 == complement(D11): {is_complement}", flush=True)

result = verify_construction(G)
if global_best_cost > 0:
    print(f"\nViolations ({len(result.violations)}):", flush=True)
    for v in result.violations:
        print(f"  {v}", flush=True)
else:
    print(f"\nVALID: {result.valid}", flush=True)
    print(f"Max red: {result.max_red_common}/{result.red_threshold}", flush=True)
    print(f"Max blue: {result.max_blue_common}/{result.blue_threshold}", flush=True)

    # Save solution
    from ramsey_core import save_construction
    save_construction(G, result, "solution_n22.json", solver="SA_relaxed")
    print("Saved to solution_n22.json", flush=True)
