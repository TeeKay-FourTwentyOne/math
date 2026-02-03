"""Quick search for n=22 using simulated annealing."""

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
D11_half = {1, 4, 6, 9, 10, 11, 13, 14, 15, 16, 17}  # First 11 from QR43
D12 = QR43.copy()

def get_D22_half(D11_half):
    return set(range(1, half_m + 1)) - D11_half

def make_graph(D11_half, D12):
    return BlockCirculantGraph(n=n, D11=D11_half, D12=D12, D22=get_D22_half(D11_half))

def get_cost(D11_half, D12):
    return compute_violation_cost(make_graph(D11_half, D12))

def random_move(D11_half, D12):
    """Apply random move and return new sets."""
    D11_new = D11_half.copy()
    D12_new = D12.copy()

    if random.random() < 0.33:  # D11 move
        d = random.randint(1, half_m)
        if d in D11_new:
            D11_new.remove(d)
        else:
            D11_new.add(d)
    else:  # D12 move
        d = random.randint(0, m - 1)
        if d in D12_new:
            D12_new.remove(d)
        else:
            D12_new.add(d)

    return D11_new, D12_new

# Initial state
current_D11 = D11_half.copy()
current_D12 = D12.copy()
current_cost = get_cost(current_D11, current_D12)

best_D11 = current_D11.copy()
best_D12 = current_D12.copy()
best_cost = current_cost

print(f"Initial cost: {current_cost}")

# SA parameters
temp = 30.0
cooling = 0.99995
max_iter = 500000

start = time.time()
random.seed(42)

for it in range(max_iter):
    if current_cost == 0:
        print(f"SOLUTION FOUND at iteration {it}!")
        break

    # Try random move
    new_D11, new_D12 = random_move(current_D11, current_D12)
    new_cost = get_cost(new_D11, new_D12)

    delta = new_cost - current_cost
    if delta < 0 or random.random() < math.exp(-delta / max(temp, 0.01)):
        current_D11 = new_D11
        current_D12 = new_D12
        current_cost = new_cost

        if current_cost < best_cost:
            best_cost = current_cost
            best_D11 = current_D11.copy()
            best_D12 = current_D12.copy()
            elapsed = time.time() - start
            print(f"Iter {it}: cost={best_cost}, temp={temp:.4f} ({elapsed:.1f}s)")

    temp *= cooling

    if it % 50000 == 0 and it > 0:
        elapsed = time.time() - start
        print(f"Progress: iter {it}, best={best_cost}, temp={temp:.4f} ({elapsed:.1f}s)")

elapsed = time.time() - start
print()
print("=" * 60)
print(f"FINAL RESULT after {elapsed:.1f}s")
print(f"Best cost: {best_cost}")

G = make_graph(best_D11, best_D12)
print(f"D11 = {sorted(G.D11)}")
print(f"D12 = {sorted(G.D12)}")
print(f"D22 = {sorted(G.D22)}")

if best_cost > 0:
    result = verify_construction(G)
    print(f"\nViolations ({len(result.violations)}):")
    for v in result.violations[:15]:
        print(f"  {v}")
    if len(result.violations) > 15:
        print(f"  ... and {len(result.violations) - 15} more")
elif best_cost == 0:
    result = verify_construction(G)
    print(f"\nVerification: valid={result.valid}")
    print(f"Max red common: {result.max_red_common} (threshold: {result.red_threshold})")
    print(f"Max blue common: {result.max_blue_common} (threshold: {result.blue_threshold})")
