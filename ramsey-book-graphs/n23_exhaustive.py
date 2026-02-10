"""
Exhaustive search for n=23 (m=45) construction.

Phase 1: Exhaustively search ALL D11 candidates (C(22,11) = 705,432 for |D11|=22).
Phase 2: For the best D11 candidates, exhaustively search for compatible D12.

Key constraints (from pattern analysis):
- D11 symmetric, D22 = complement(D11) in {1,...,44}
- |D12| = 22, 0 in D12
- V1V1 red: Delta(D11,D11,d) + Delta(D12,D12,d) <= 21 for d in D11
- V1V1 blue: blue_common <= 22 for d not in D11
- V2V2 red: Delta(D22,D22,d) + Delta(D12T,D12T,d) <= 21 for d in D22
- V1V2: AUTOMATIC (proven)

Since D22 = complement(D11):
  Delta(D22,D22,d) = Delta(D11,D11,d) + (m-2-2|D11|) + 2*[d in D11]
  V2V2(d) = Delta(D11,D11,d) + Delta(D12,D12,m-d) + (m-2-2|D11|) + 2*[d in D11]

So the constraints on D12 given D11 are:
  For d in D11: A(d) + B(d) <= 21     (V1V1 red)
  For d not in D11: A(d) + B(45-d) + (43-2|D11|) <= 21   (V2V2 red)
  Blue constraints are implied by totals.

where A(d) = Delta(D11,D11,d), B(d) = Delta(D12,D12,d).
"""

import sys, os, math, itertools, json
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ramsey_core import Delta, BlockCirculantGraph, verify_construction, save_construction

m = 45
n = 23

# Compute symmetric pairs: {x, m-x} for x = 1,...,22
pairs = []
for x in range(1, m):
    neg_x = (-x) % m
    if x <= neg_x:
        pairs.append((x, neg_x))

assert len(pairs) == 22

def make_d11_from_pairs(selected_pair_indices):
    D11 = set()
    for i in selected_pair_indices:
        D11.add(pairs[i][0])
        D11.add(pairs[i][1])
    return D11

def compute_A(D11_set):
    """Compute A(d) = Delta(D11,D11,d) for all d in {1,...,m-1}."""
    A = {}
    for d in range(1, m):
        A[d] = Delta(D11_set, D11_set, d, m)
    return A

def score_d11_fast(D11_set, d11_size):
    """Fast scoring: returns (min_budget, max_A_red, max_A_blue)."""
    A_red_max = 0
    A_blue_max = 0
    for d in range(1, m):
        a = Delta(D11_set, D11_set, d, m)
        if d in D11_set:
            A_red_max = max(A_red_max, a)
        else:
            A_blue_max = max(A_blue_max, a)

    red_thresh = 21  # n - 2
    v22_const = m - 2 - 2 * d11_size

    max_B_red = red_thresh - A_red_max
    max_B_v22 = red_thresh - v22_const - A_blue_max

    return min(max_B_red, max_B_v22), A_red_max, A_blue_max


# =============================================================================
# PHASE 1: Exhaustive D11 search for |D11| = 22 (11 of 22 pairs)
# =============================================================================
print("=" * 80)
print("PHASE 1: Exhaustive D11 search for |D11| = 22")
print(f"Total candidates: C(22,11) = {math.comb(22, 11)}")
print("=" * 80)

best_budget = -999
best_candidates_22 = []
count = 0

for selected in itertools.combinations(range(22), 11):
    D11 = make_d11_from_pairs(selected)
    budget, mar, mab = score_d11_fast(D11, 22)
    count += 1

    if budget > best_budget:
        best_budget = budget
        best_candidates_22 = [(selected, sorted(D11), mar, mab, budget)]
        if count % 100000 == 0 or budget > 8:
            print(f"  [{count}] New best budget={budget}, max_A_red={mar}, max_A_blue={mab}")
    elif budget == best_budget:
        best_candidates_22.append((selected, sorted(D11), mar, mab, budget))

print(f"\nChecked {count} candidates")
print(f"Best min_budget = {best_budget}")
print(f"Number of D11 with best budget: {len(best_candidates_22)}")

# Show top candidates (up to 20)
for i, (sel, d11, mar, mab, bud) in enumerate(best_candidates_22[:20]):
    print(f"  #{i}: D11 = {d11}, max_A_red={mar}, max_A_blue={mab}, budget={bud}")


# =============================================================================
# PHASE 1b: Exhaustive D11 search for |D11| = 24 (12 of 22 pairs)
# =============================================================================
print("\n" + "=" * 80)
print("PHASE 1b: Exhaustive D11 search for |D11| = 24")
print(f"Total candidates: C(22,12) = {math.comb(22, 12)}")
print("=" * 80)

best_budget_24 = -999
best_candidates_24 = []
count = 0

for selected in itertools.combinations(range(22), 12):
    D11 = make_d11_from_pairs(selected)
    budget, mar, mab = score_d11_fast(D11, 24)
    count += 1

    if budget > best_budget_24:
        best_budget_24 = budget
        best_candidates_24 = [(selected, sorted(D11), mar, mab, budget)]
        if count % 100000 == 0 or budget > 8:
            print(f"  [{count}] New best budget={budget}, max_A_red={mar}, max_A_blue={mab}")
    elif budget == best_budget_24:
        best_candidates_24.append((selected, sorted(D11), mar, mab, budget))

print(f"\nChecked {count} candidates")
print(f"Best min_budget = {best_budget_24}")
print(f"Number of D11 with best budget: {len(best_candidates_24)}")

for i, (sel, d11, mar, mab, bud) in enumerate(best_candidates_24[:20]):
    print(f"  #{i}: D11 = {d11}, max_A_red={mar}, max_A_blue={mab}, budget={bud}")


# =============================================================================
# PHASE 2: For best D11, search for compatible D12
# =============================================================================
print("\n" + "=" * 80)
print("PHASE 2: D12 search for best D11 candidates")
print("=" * 80)

# Pick the best D11 size
if best_budget >= best_budget_24:
    print(f"Using |D11|=22 candidates (budget={best_budget} >= {best_budget_24})")
    candidates = best_candidates_22
    d11_size = 22
else:
    print(f"Using |D11|=24 candidates (budget={best_budget_24} > {best_budget})")
    candidates = best_candidates_24
    d11_size = 24

# For each D11, compute per-d budgets for B(d)
# Then try to find D12 satisfying all budgets simultaneously

def compute_b_budgets(D11_set, d11_size):
    """Compute the maximum allowed B(d) = Delta(D12,D12,d) for each d."""
    red_thresh = 21
    v22_const = m - 2 - 2 * d11_size

    A = compute_A(D11_set)

    # V1V1 red: A(d) + B(d) <= 21 for d in D11
    # So B(d) <= 21 - A(d) for d in D11
    budgets = {}
    for d in range(1, m):
        if d in D11_set:
            # V1V1 red constraint on B(d)
            budgets[d] = red_thresh - A[d]
        # V2V2 red: for each d' in D22, A(d') + B(m-d') + v22_const <= 21
        # i.e., B(m-d') <= 21 - v22_const - A(d')
        # Let d = m - d', then B(d) <= 21 - v22_const - A(m-d) for (m-d) in D22

    for d_prime in range(1, m):
        if d_prime not in D11_set:  # d' in D22
            d = (-d_prime) % m  # d = m - d'
            budget_v22 = red_thresh - v22_const - A[d_prime]
            if d in budgets:
                budgets[d] = min(budgets[d], budget_v22)
            else:
                budgets[d] = budget_v22

    # Blue constraints:
    # V1V1 blue: for d not in D11, blue_common = (N-2) - 2*d1 + (A(d)+B(d))
    # Need blue_common <= 22, so A(d)+B(d) <= 22 - (N-2) + 2*d1 = 22 - 88 + 2*(22+22) = 22 - 88 + 88 = 22
    # Wait: N = 2*45 = 90, d1 = |D11| + |D12| = d11_size + 22
    d1 = d11_size + 22
    N = 90
    blue_limit_v11 = 22 - (N - 2) + 2 * d1  # = 22 - 88 + 2*d1 = 2*d1 - 66
    for d in range(1, m):
        if d not in D11_set:
            budget_blue_v11 = blue_limit_v11 - A[d]
            if d in budgets:
                budgets[d] = min(budgets[d], budget_blue_v11)
            else:
                budgets[d] = budget_blue_v11

    # V2V2 blue: for d in D22 (d not in D11):
    # blue_common = (N-2) - 2*d2 + V2V2(d)
    # V2V2(d) = A(d) + B(m-d) + v22_const + 2*[d in D11]
    # For d in D11: blue_common = (N-2) - 2*d2 + A(d) + B(m-d) + v22_const + 2
    # Need <= 22
    d2 = (m - 1 - d11_size) + 22  # |D22| + |D12|
    for d in range(1, m):
        if d in D11_set:  # V2V2 blue for d in D11 (d in QR for V2V2)
            d_neg = (-d) % m
            budget_blue_v22 = 22 - (N - 2) + 2 * d2 - v22_const - 2 - A[d]
            if d_neg in budgets:
                budgets[d_neg] = min(budgets[d_neg], budget_blue_v22)
            else:
                budgets[d_neg] = budget_blue_v22

    return budgets, A


# Try an LP-relaxation approach: for each D11, determine if budgets are satisfiable
print("\nAnalyzing budget constraints for top candidates...")

for idx, (sel, d11_list, mar, mab, bud) in enumerate(candidates[:5]):
    D11_set = set(d11_list)
    budgets, A = compute_b_budgets(D11_set, d11_size)

    print(f"\n--- Candidate #{idx}: D11 = {d11_list} ---")
    print(f"  max_A_red={mar}, max_A_blue={mab}, overall_budget={bud}")

    # Show budgets
    all_positive = all(b >= 0 for b in budgets.values())
    min_b = min(budgets.values())
    print(f"  B(d) budgets: min={min_b}, all positive: {all_positive}")

    if not all_positive:
        neg_ds = [d for d, b in budgets.items() if b < 0]
        print(f"  INFEASIBLE: negative budgets at d = {neg_ds}")
        continue

    # Show A(d) distribution for d in D11
    A_vals_red = [(d, A[d]) for d in range(1, m) if d in D11_set]
    A_vals_red.sort(key=lambda x: -x[1])
    print(f"  A(d) for d in D11 (sorted desc): {[(d, a) for d, a in A_vals_red[:5]]}...")
    print(f"  B budget for d in D11: {[(d, budgets.get(d, '?')) for d, a in A_vals_red[:5]]}...")

    # Average B(d) budget
    avg_budget = sum(budgets.values()) / len(budgets)
    print(f"  Average B budget: {avg_budget:.2f}")

    # What B(d) values can a D12 with |D12|=22 produce?
    # B(d) = Delta(D12,D12,d), average = 22*21/44 = 10.5
    avg_B = 22 * 21 / 44
    print(f"  Expected average B(d) = {avg_B:.2f}")

    # Check: can we find D12 with B(d) <= budget(d) for all d?
    # D12 has 22 elements from {0,...,44}, with 0 in D12 (by our convention)
    # Try random sampling of D12
    import random
    random.seed(100 + idx)

    best_violations = 999
    best_d12 = None
    for trial in range(100000):
        # Pick 21 more elements from {1,...,44}
        D12_rest = random.sample(range(1, m), 21)
        D12_set = set([0] + D12_rest)

        # Check budgets
        violations = 0
        for d in range(1, m):
            B_d = Delta(D12_set, D12_set, d, m)
            if d in budgets and B_d > budgets[d]:
                violations += B_d - budgets[d]

        if violations < best_violations:
            best_violations = violations
            best_d12 = sorted(D12_set)
            if violations == 0:
                print(f"  FOUND VALID D12 at trial {trial}!")
                print(f"  D12 = {best_d12}")
                # Full verification
                D22_set = set(range(1, m)) - D11_set
                G = BlockCirculantGraph(n=n, D11=D11_set, D12=D12_set, D22=D22_set)
                result = verify_construction(G)
                print(f"  Full verification: valid={result.valid}")
                print(f"  max_red={result.max_red_common}, max_blue={result.max_blue_common}")
                if result.valid:
                    filename = f"/Users/stephenpadgett/Projects/math/ramsey-book-graphs/solution_n23.json"
                    save_construction(G, result, filename, solver="exhaustive-n23")
                    print(f"  SAVED to {filename}")
                break

        if trial % 20000 == 0 and trial > 0:
            print(f"  Trial {trial}: best violations = {best_violations}")

    if best_violations > 0:
        print(f"  No valid D12 found (best violations = {best_violations})")
        print(f"  Best D12 attempt: {best_d12}")


# =============================================================================
# PHASE 3: SA-based D12 search for best D11s
# =============================================================================
print("\n" + "=" * 80)
print("PHASE 3: SA-based joint D12 search")
print("=" * 80)

import random, math as mth

def sa_search_d12(D11_set, d11_size, max_iter=500000, seed=0):
    """Simulated annealing to find D12 compatible with given D11."""
    random.seed(seed)
    D22_set = set(range(1, m)) - D11_set

    budgets, A = compute_b_budgets(D11_set, d11_size)

    # Initialize D12 randomly: 0 plus 21 elements from {1,...,44}
    D12_list = [0] + random.sample(range(1, m), 21)
    D12_set = set(D12_list)
    not_in_D12 = set(range(1, m)) - (D12_set - {0})

    def cost(D12_s):
        """Cost = sum of budget violations."""
        total = 0
        for d in range(1, m):
            B_d = Delta(D12_s, D12_s, d, m)
            if d in budgets and B_d > budgets[d]:
                total += B_d - budgets[d]
        return total

    current_cost = cost(D12_set)
    best_cost = current_cost
    best_D12 = sorted(D12_set)

    T = 5.0
    T_min = 0.001
    alpha = 0.99999

    for it in range(max_iter):
        # Swap: remove one element (not 0), add one not in D12
        D12_others = [x for x in D12_set if x != 0]
        remove = random.choice(D12_others)
        not_in = list(set(range(1, m)) - D12_set)
        add = random.choice(not_in)

        D12_set.remove(remove)
        D12_set.add(add)

        new_cost = cost(D12_set)
        delta = new_cost - current_cost

        if delta <= 0 or random.random() < mth.exp(-delta / T):
            current_cost = new_cost
            if new_cost < best_cost:
                best_cost = new_cost
                best_D12 = sorted(D12_set)
                if best_cost == 0:
                    return best_D12, best_cost
        else:
            D12_set.remove(add)
            D12_set.add(remove)

        T = max(T * alpha, T_min)

        if it % 100000 == 0:
            print(f"    SA iter {it}: cost={current_cost}, best={best_cost}, T={T:.4f}")

    return best_D12, best_cost


for idx, (sel, d11_list, mar, mab, bud) in enumerate(candidates[:10]):
    D11_set = set(d11_list)
    print(f"\n--- SA for candidate #{idx}: D11 = {d11_list[:6]}... ---")

    for seed in range(3):
        d12, cost = sa_search_d12(D11_set, d11_size, max_iter=300000, seed=seed*42)
        print(f"  Seed {seed}: best cost = {cost}")
        if cost == 0:
            D12_set = set(d12)
            D22_set = set(range(1, m)) - D11_set
            G = BlockCirculantGraph(n=n, D11=D11_set, D12=D12_set, D22=D22_set)
            result = verify_construction(G)
            print(f"  Full verification: valid={result.valid}")
            print(f"  max_red={result.max_red_common}, max_blue={result.max_blue_common}")
            print(f"  D12 = {d12}")
            if result.valid:
                filename = "/Users/stephenpadgett/Projects/math/ramsey-book-graphs/solution_n23.json"
                save_construction(G, result, filename, solver="exhaustive+SA-n23")
                print(f"  *** SOLUTION FOUND AND SAVED! ***")
                sys.exit(0)
            break
