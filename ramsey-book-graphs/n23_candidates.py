"""
Generate candidate D11 sets for n=23 (m=45) using structural constraints.

Key constraints:
- D11 must be symmetric (x in D11 iff -x mod 45 in D11)
- |D11| must be even (22 or 24)
- D22 = complement(D11) in {1,...,44}
- |D12| = 22 (= n-1), with 0 in D12

Z_45 = Z_9 x Z_5. The negation map x -> -x becomes (a,b) -> (-a mod 9, -b mod 5).

Symmetric pairs under negation in Z_45:
{x, 45-x} for x = 1,...,22. There are 22 pairs.

For |D11| = 22: pick 11 of 22 symmetric pairs.
For |D11| = 24: pick 12 of 22 symmetric pairs.

We'll evaluate each candidate by computing Delta(D11,D11,d) and checking whether
a valid D12 can potentially exist (using average-case analysis).
"""

import sys, os, json, math, itertools, random
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ramsey_core import Delta, Sigma, BlockCirculantGraph, verify_construction

m = 45
n = 23

# Compute symmetric pairs
pairs = []
for x in range(1, m):
    neg_x = (-x) % m
    if x <= neg_x:
        pairs.append((x, neg_x))

print(f"Symmetric pairs in Z_{m}: {len(pairs)}")
for i, (a, b) in enumerate(pairs):
    crt_a = (a % 9, a % 5)
    crt_b = (b % 9, b % 5)
    print(f"  Pair {i}: {{{a}, {b}}} -> CRT ({crt_a}, {crt_b})")

# Multiplier orbits
def units_mod(m):
    return {x for x in range(1, m) if math.gcd(x, m) == 1}

units = units_mod(m)
print(f"\nUnits of Z_{m}: {len(units)} elements")

# Orbit structure on pairs
# Two pairs are in the same orbit if one can be mapped to the other by a unit multiplier
pair_visited = [False] * len(pairs)
pair_orbits = []
for i in range(len(pairs)):
    if pair_visited[i]:
        continue
    orbit = set()
    for u in units:
        mapped = tuple(sorted([(u * pairs[i][0]) % m, (u * pairs[i][1]) % m]))
        # Find which pair this is
        for j, (a, b) in enumerate(pairs):
            if (a, b) == mapped or (b, a) == mapped:
                orbit.add(j)
                pair_visited[j] = True
                break
    pair_orbits.append(sorted(orbit))

print(f"\nPair orbits under unit multipliers: {len(pair_orbits)}")
for i, orb in enumerate(pair_orbits):
    print(f"  Orbit {i} (size {len(orb)}): pairs {orb}")
    for j in orb:
        print(f"    {pairs[j]}")


# For each D11 size, evaluate candidates
def evaluate_d11(D11_set, verbose=False):
    """Evaluate a D11 candidate by computing Delta(D11,D11,d) for all d."""
    d11_size = len(D11_set)
    d22_size = m - 1 - d11_size

    # For this to work, we need |D12| = n-1 = 22
    d12_size = n - 1  # = 22

    # Compute A(d) = Delta(D11, D11, d) for all d
    A = {}
    for d in range(1, m):
        A[d] = Delta(D11_set, D11_set, d, m)

    # Average A(d)
    avg_A = sum(A.values()) / (m - 1)

    # For V1V1: need A(d) + B(d) <= n-2 = 21 for d in D11
    # For V1V1 blue: need A(d) + B(d) <= blue_max for d not in D11
    # Average B(d) = |D12|(|D12|-1)/(m-1) = 22*21/44 = 10.5
    avg_B = d12_size * (d12_size - 1) / (m - 1)
    avg_AB = avg_A + avg_B

    # P(0)/m = (|D11|^2 + |D12|^2)/m = (d11_size^2 + d12_size^2)/45
    P0_m = (d11_size ** 2 + d12_size ** 2) / m

    # Red threshold
    red_thresh = n - 2  # = 21

    # A(d) for d in D11 vs not in D11
    A_red = [A[d] for d in range(1, m) if d in D11_set]
    A_blue = [A[d] for d in range(1, m) if d not in D11_set]

    max_A_red = max(A_red)
    avg_A_red = sum(A_red) / len(A_red) if A_red else 0

    # The max A(d) for red edges determines the "budget" for B(d)
    # Need B(d) <= red_thresh - A(d) for each d in D11
    # The tightest constraint is at d with max A(d) in D11: need B(d) <= 21 - max_A_red
    max_B_allowed_red = red_thresh - max_A_red

    # For V2V2: need A(d) + B(m-d) + (m-2-2|D11|) + 2*[d in D11] <= 21 for d in D22
    # D22 = {d not in D11}, so [d in D11] = 0 for these.
    # const = m - 2 - 2*d11_size = 45 - 2 - 2*d11_size = 43 - 2*d11_size
    v22_const = m - 2 - 2 * d11_size

    # For d in D22 (d not in D11):
    # V2V2(d) = A(d) + B(m-d) + v22_const <= 21
    # So B(m-d) <= 21 - A(d) - v22_const = 21 - A(d) - (43 - 2*d11_size) = 2*d11_size - 22 - A(d)
    max_A_blue = max(A_blue) if A_blue else 0
    max_B_allowed_v22 = 2 * d11_size - 22 - max_A_blue

    if verbose:
        print(f"\n  |D11| = {d11_size}, |D22| = {d22_size}, |D12| = {d12_size}")
        print(f"  avg A(d) = {avg_A:.3f}, P(0)/m = {P0_m:.3f}")
        print(f"  A(d) for d in D11: min={min(A_red)}, max={max_A_red}, avg={avg_A_red:.3f}")
        print(f"  A(d) for d not in D11: min={min(A_blue)}, max={max_A_blue}, avg={sum(A_blue)/len(A_blue):.3f}")
        print(f"  V1V1 red: need A(d)+B(d) <= {red_thresh}, max_A_red = {max_A_red}, max B allowed = {max_B_allowed_red}")
        print(f"  V2V2 red const = {v22_const}")
        print(f"  V2V2 red: need A(d)+B(m-d) <= {red_thresh - v22_const} for d in D22")
        print(f"  max A(d) for blue = {max_A_blue}, max B allowed for V2V2 = {max_B_allowed_v22}")

    return {
        "D11_size": d11_size,
        "max_A_red": max_A_red,
        "max_A_blue": max_A_blue,
        "avg_A_red": avg_A_red,
        "max_B_allowed_red": max_B_allowed_red,
        "max_B_allowed_v22": max_B_allowed_v22,
        "v22_const": v22_const,
        "A": A,
    }


# Exhaustive search over pair selections for |D11| = 22 (11 pairs from 22)
# C(22,11) = 705,432 -- feasible!
print("\n" + "=" * 80)
print("EVALUATING D11 CANDIDATES WITH |D11| = 22 (picking 11 of 22 pairs)")
print("=" * 80)

best_candidates = []
total_checked = 0

# Instead of exhaustive (too slow in single run), let's use orbit structure
# and random sampling with scoring

def make_d11_from_pairs(selected_pair_indices):
    """Create D11 from selected pair indices."""
    D11 = set()
    for i in selected_pair_indices:
        D11.add(pairs[i][0])
        D11.add(pairs[i][1])
    return D11


def score_d11(D11_set, d11_size):
    """Score a D11 set: lower is better. Returns (max_A_red, feasibility_score)."""
    A_red_vals = []
    A_blue_vals = []
    for d in range(1, m):
        a = Delta(D11_set, D11_set, d, m)
        if d in D11_set:
            A_red_vals.append(a)
        else:
            A_blue_vals.append(a)

    max_A_red = max(A_red_vals)
    max_A_blue = max(A_blue_vals)

    # Score: lower max_A_red is better (more room for B(d))
    # Also want low max_A_blue (for V2V2)
    red_thresh = n - 2  # 21
    v22_const = m - 2 - 2 * d11_size

    # Budget for B at tightest red constraint
    max_B_red = red_thresh - max_A_red
    # Budget for B at tightest V2V2 constraint
    max_B_v22 = red_thresh - v22_const - max_A_blue

    # Need both budgets > 0 (and ideally large)
    return max_A_red, max_A_blue, min(max_B_red, max_B_v22)


# Random search with scoring
random.seed(42)
print("\nRandom sampling 50000 D11 with |D11|=22...")
best_score = -999
best_d11_22 = None

for trial in range(50000):
    selected = random.sample(range(22), 11)
    D11 = make_d11_from_pairs(selected)
    max_A_red, max_A_blue, min_budget = score_d11(D11, 22)

    if min_budget > best_score:
        best_score = min_budget
        best_d11_22 = (D11.copy(), selected, max_A_red, max_A_blue, min_budget)
        if trial % 5000 == 0 or min_budget > 3:
            print(f"  Trial {trial}: max_A_red={max_A_red}, max_A_blue={max_A_blue}, min_budget={min_budget}")

if best_d11_22:
    D11_best, pairs_best, mar, mab, mb = best_d11_22
    print(f"\nBest |D11|=22 candidate:")
    print(f"  D11 = {sorted(D11_best)}")
    print(f"  max_A_red = {mar}, max_A_blue = {mab}, min_budget = {mb}")
    evaluate_d11(D11_best, verbose=True)

# Now try |D11| = 24 (12 of 22 pairs)
print("\n" + "=" * 80)
print("EVALUATING D11 CANDIDATES WITH |D11| = 24 (picking 12 of 22 pairs)")
print("=" * 80)

random.seed(43)
print("\nRandom sampling 50000 D11 with |D11|=24...")
best_score_24 = -999
best_d11_24 = None

for trial in range(50000):
    selected = random.sample(range(22), 12)
    D11 = make_d11_from_pairs(selected)
    max_A_red, max_A_blue, min_budget = score_d11(D11, 24)

    if min_budget > best_score_24:
        best_score_24 = min_budget
        best_d11_24 = (D11.copy(), selected, max_A_red, max_A_blue, min_budget)
        if trial % 5000 == 0 or min_budget > 3:
            print(f"  Trial {trial}: max_A_red={max_A_red}, max_A_blue={max_A_blue}, min_budget={min_budget}")

if best_d11_24:
    D11_best, pairs_best, mar, mab, mb = best_d11_24
    print(f"\nBest |D11|=24 candidate:")
    print(f"  D11 = {sorted(D11_best)}")
    print(f"  max_A_red = {mar}, max_A_blue = {mab}, min_budget = {mb}")
    evaluate_d11(D11_best, verbose=True)

print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)
print(f"Best |D11|=22: min_budget = {best_score}")
print(f"Best |D11|=24: min_budget = {best_score_24}")
print("The configuration with higher min_budget has more room for D12 to satisfy constraints.")
