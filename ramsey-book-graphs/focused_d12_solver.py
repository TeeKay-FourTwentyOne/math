"""
Focused D12 solver: fixes a good D11 and searches only over D12.

For primes p ≡ 3 mod 4, the joint D11+D12 SA search gets stuck because
it needs to simultaneously find a good D11 AND a matching D12. This solver
separates the problem:
  Phase A: Find LP-feasible D11 candidates (fast, spectral optimization)
  Phase B: For each D11, search for D12 using SA (focused, smaller search space)

Usage: python focused_d12_solver.py <n> [max_seeds] [max_iter]
"""

import sys, os, math, random, time
import numpy as np
from scipy.optimize import linprog

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ramsey_core import BlockCirculantGraph, verify_construction, save_construction


def compute_delta_fft(indicator, m):
    a = np.array(indicator, dtype=np.float64)
    fa = np.fft.rfft(a)
    result = np.fft.irfft(np.conj(fa) * fa, n=m)
    return np.round(result).astype(np.int64)


def compute_cost(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n):
    d11_size = int(np.sum(D11_ind))
    d12_size = int(np.sum(D12_ind))
    d22_size = m - 1 - d11_size
    d1 = d11_size + d12_size
    d2 = d22_size + d12_size
    N = 2 * m
    red_thresh = n - 2
    blue_thresh = n - 1

    common_v11 = delta_11[1:] + delta_12[1:]
    d11_mask = D11_ind[1:].astype(bool)

    red_excess_v11 = np.maximum(common_v11[d11_mask] - red_thresh, 0)
    blue_common_v11 = (N - 2) - 2 * d1 + common_v11[~d11_mask]
    blue_excess_v11 = np.maximum(blue_common_v11 - blue_thresh, 0)

    v22_const = m - 2 - 2 * d11_size
    d22_delta = delta_11[1:] + v22_const + 2 * D11_ind[1:]
    common_v22 = d22_delta + delta_12T[1:]

    red_excess_v22 = np.maximum(common_v22[~d11_mask] - red_thresh, 0)
    blue_common_v22 = (N - 2) - 2 * d2 + common_v22[d11_mask]
    blue_excess_v22 = np.maximum(blue_common_v22 - blue_thresh, 0)

    return (int(np.sum(red_excess_v11)) + int(np.sum(blue_excess_v11)) +
            int(np.sum(red_excess_v22)) + int(np.sum(blue_excess_v22)))


def get_pairs(m):
    pairs = []
    seen = set()
    for x in range(1, m):
        if x not in seen:
            pairs.append((x, (-x) % m))
            seen.add(x)
            seen.add((-x) % m)
    return pairs


def lp_feasible_d11(p, D11):
    """Check LP feasibility and return margin."""
    n = (p + 1) // 2
    F = (p - 1) // 2
    D11_set = set(D11)
    d11_size = len(D11)
    d12_size = (p - 1) // 2

    ind11 = np.zeros(p, dtype=np.float64)
    for j in D11:
        ind11[j] = 1.0
    dft11 = np.fft.fft(ind11)
    power11 = np.abs(dft11) ** 2

    P0 = d11_size ** 2 + d12_size ** 2
    S = p * (d11_size + d12_size) - P0

    c = np.zeros(F + 1)
    c[F] = 1.0

    A_ub_list = []
    b_ub_list = []

    for d in sorted(D11):
        row = np.zeros(F + 1)
        for k in range(1, F + 1):
            row[k - 1] = (2.0 / p) * math.cos(2 * math.pi * k * d / p)
        row[F] = -1.0
        A_ub_list.append(row)
        b_ub_list.append(-P0 / p)

    row_t = np.zeros(F + 1)
    row_t[F] = 1.0
    A_ub_list.append(row_t)
    b_ub_list.append(n - 2)

    D22 = set(range(1, p)) - D11_set
    for d in sorted(D22):
        row = np.zeros(F + 1)
        for k in range(1, F + 1):
            row[k - 1] = (2.0 / p) * math.cos(2 * math.pi * k * d / p)
        A_ub_list.append(row)
        b_ub_list.append(n + 1 - P0 / p)

    A_eq = np.zeros((1, F + 1))
    for k in range(F):
        A_eq[0, k] = 2.0
    b_eq = np.array([S])
    bounds = [(power11[k], None) for k in range(1, F + 1)] + [(None, None)]

    result = linprog(c, A_ub=np.array(A_ub_list), b_ub=np.array(b_ub_list),
                     A_eq=A_eq, b_eq=b_eq, bounds=bounds, method='highs',
                     options={'presolve': True})

    if result.success:
        return (n - 2) - result.x[F]
    return -999


def find_good_d11(m, num_candidates=500, verbose=True):
    """Find LP-feasible D11 candidates for prime m ≡ 3 mod 4."""
    pairs = get_pairs(m)
    num_pairs = len(pairs)
    num_selected = (m + 1) // 4  # |D11| = (m+1)/2, need (m+1)/4 pairs
    rng = np.random.default_rng(42)

    candidates = []
    for trial in range(num_candidates):
        chosen = rng.choice(num_pairs, size=num_selected, replace=False)
        D11 = set()
        for i in chosen:
            D11.add(pairs[i][0])
            D11.add(pairs[i][1])

        margin = lp_feasible_d11(m, D11)
        if margin > -900:
            candidates.append((margin, D11))

    candidates.sort(key=lambda x: -x[0])  # best margin first
    if verbose:
        print(f"  Found {len(candidates)}/{num_candidates} LP-feasible D11", flush=True)
        if candidates:
            print(f"  Best margin: {candidates[0][0]:.4f}", flush=True)
    return candidates


def find_good_d11_composite(m, num_candidates=200, verbose=True):
    """Find good D11 for composite m using autocorrelation analysis."""
    pairs = get_pairs(m)
    num_pairs = len(pairs)
    n = (m + 1) // 2

    # For composite m, try several D11 sizes
    half = (m - 1) // 2
    d11_sizes = []
    for s in [half - 1, half, half + 1, half - 2, half + 2]:
        if 0 < s < m - 1 and s % 2 == 0 and (s // 2) <= num_pairs:
            d11_sizes.append(s)
    d11_sizes = sorted(set(d11_sizes), key=lambda s: abs(s - half))

    rng = np.random.default_rng(42)
    candidates = []

    for d11_size in d11_sizes:
        num_pair_sel = d11_size // 2
        for trial in range(num_candidates // len(d11_sizes)):
            chosen = rng.choice(num_pairs, size=num_pair_sel, replace=False)
            D11_ind = np.zeros(m, dtype=np.int64)
            for i in chosen:
                D11_ind[pairs[i][0]] = 1
                D11_ind[pairs[i][1]] = 1

            delta_11 = compute_delta_fft(D11_ind, m)
            d11_mask = D11_ind[1:].astype(bool)

            # For a good D11, max A(d) at D11 positions should be low
            A_d11 = delta_11[1:][d11_mask]
            max_A = int(np.max(A_d11))

            # Also check: average A(d) at D11 positions vs threshold
            avg_A = float(np.mean(A_d11))
            thresh = n - 2

            # Score: lower max_A and lower spread are better
            score = -(max_A + 0.5 * float(np.std(A_d11)))

            D11 = set(int(i) for i in range(m) if D11_ind[i])
            candidates.append((score, D11, d11_size, max_A))

    candidates.sort(key=lambda x: -x[0])
    if verbose:
        print(f"  Tested {len(candidates)} D11 candidates", flush=True)
        if candidates:
            print(f"  Best: |D11|={candidates[0][2]}, max_A={candidates[0][3]}", flush=True)
    return [(c[0], c[1]) for c in candidates[:20]]  # top 20


def sa_d12_only(n, D11, max_iter=20000000, verbose=True):
    """SA search for D12 with D11 fixed."""
    m = 2 * n - 1
    d12_size = n - 1

    D11_ind = np.zeros(m, dtype=np.int64)
    for x in D11:
        D11_ind[x] = 1
    delta_11 = compute_delta_fft(D11_ind, m)

    best_cost = float('inf')
    best_D12 = None

    # Random D12 initialization
    D12_ind = np.zeros(m, dtype=np.int64)
    D12_ind[0] = 1  # 0 always in D12
    rest = random.sample(range(1, m), d12_size - 1)
    for x in rest:
        D12_ind[x] = 1

    D12T_ind = np.zeros(m, dtype=np.int64)
    for x in range(m):
        if D12_ind[x]:
            D12T_ind[(-x) % m] = 1

    delta_12 = compute_delta_fft(D12_ind, m)
    delta_12T = compute_delta_fft(D12T_ind, m)

    current_cost = compute_cost(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n)
    best_cost = current_cost
    best_D12 = np.copy(D12_ind)

    T = 8.0
    T_min = 0.00005
    alpha = 1 - 5.0 / max_iter
    last_imp = 0
    reheat_count = 0

    for it in range(max_iter):
        # Only D12 moves: single swap (70%) or double swap (30%)
        if random.random() < 0.7:
            # Single swap
            d12_list = np.where(D12_ind[1:] == 1)[0] + 1
            not_in = np.where(D12_ind[1:] == 0)[0] + 1
            if len(d12_list) == 0 or len(not_in) == 0:
                continue
            rem = d12_list[random.randint(0, len(d12_list) - 1)]
            add = not_in[random.randint(0, len(not_in) - 1)]

            D12_ind[rem] = 0
            D12_ind[add] = 1
            D12T_ind[(-rem) % m] = 0
            D12T_ind[(-add) % m] = 1

            delta_12_new = compute_delta_fft(D12_ind, m)
            delta_12T_new = compute_delta_fft(D12T_ind, m)

            new_cost = compute_cost(D11_ind, D12_ind, delta_11, delta_12_new, delta_12T_new, m, n)
            dc = new_cost - current_cost

            if dc <= 0 or random.random() < math.exp(-dc / max(T, T_min)):
                current_cost = new_cost
                delta_12 = delta_12_new
                delta_12T = delta_12T_new
            else:
                D12_ind[rem] = 1
                D12_ind[add] = 0
                D12T_ind[(-rem) % m] = 1
                D12T_ind[(-add) % m] = 0
        else:
            # Double swap
            d12_list = np.where(D12_ind[1:] == 1)[0] + 1
            not_in = np.where(D12_ind[1:] == 0)[0] + 1
            if len(d12_list) < 2 or len(not_in) < 2:
                continue
            ri = np.random.choice(len(d12_list), 2, replace=False)
            ai = np.random.choice(len(not_in), 2, replace=False)
            rems = d12_list[ri]
            adds = not_in[ai]

            for rv in rems:
                D12_ind[rv] = 0
                D12T_ind[(-rv) % m] = 0
            for av in adds:
                D12_ind[av] = 1
                D12T_ind[(-av) % m] = 1

            delta_12_new = compute_delta_fft(D12_ind, m)
            delta_12T_new = compute_delta_fft(D12T_ind, m)

            new_cost = compute_cost(D11_ind, D12_ind, delta_11, delta_12_new, delta_12T_new, m, n)
            dc = new_cost - current_cost

            if dc <= 0 or random.random() < math.exp(-dc / max(T, T_min)):
                current_cost = new_cost
                delta_12 = delta_12_new
                delta_12T = delta_12T_new
            else:
                for rv in rems:
                    D12_ind[rv] = 1
                    D12T_ind[(-rv) % m] = 1
                for av in adds:
                    D12_ind[av] = 0
                    D12T_ind[(-av) % m] = 0

        if current_cost < best_cost:
            best_cost = current_cost
            best_D12 = np.copy(D12_ind)
            last_imp = it
            if best_cost == 0:
                break

        # Aggressive reheating for D12-only search
        if it - last_imp > 300000 and best_cost > 0 and best_cost <= 20:
            T = 5.0
            last_imp = it
            reheat_count += 1
            D12_ind = np.copy(best_D12)
            D12T_ind = np.zeros(m, dtype=np.int64)
            for x in range(m):
                if D12_ind[x]:
                    D12T_ind[(-x) % m] = 1
            delta_12 = compute_delta_fft(D12_ind, m)
            delta_12T = compute_delta_fft(D12T_ind, m)
            current_cost = best_cost

        T = max(T * alpha, T_min)

        if it % 2000000 == 0 and it > 0 and verbose:
            print(f"      {it // 1000000}M: cost={current_cost} best={best_cost} T={T:.4f} rh={reheat_count}", flush=True)

    return best_cost, best_D12


def solve(n, max_seeds=30, max_iter=20000000):
    m = 2 * n - 1
    print(f"=== Focused D12 Solver: n={n}, m={m} ===", flush=True)

    is_prime = all(m % i != 0 for i in range(2, int(m**0.5) + 1))

    start = time.time()

    # Phase A: Find good D11 candidates
    print(f"\nPhase A: Finding good D11 candidates...", flush=True)
    if is_prime and m % 4 == 3:
        candidates = find_good_d11(m, num_candidates=1000)
    else:
        candidates = find_good_d11_composite(m, num_candidates=500)

    if not candidates:
        print("No good D11 found!", flush=True)
        return False

    # Phase B: For each D11, search for D12
    print(f"\nPhase B: Searching for D12 (max {max_seeds} seeds per D11)...", flush=True)

    num_d11_to_try = min(10, len(candidates))
    for d11_idx in range(num_d11_to_try):
        _, D11 = candidates[d11_idx]
        D11_sorted = sorted(D11)
        print(f"\n  D11 #{d11_idx} (|D11|={len(D11)}): {D11_sorted[:6]}...", flush=True)

        for seed in range(max_seeds):
            random.seed(seed * 1000 + d11_idx)
            np.random.seed(seed * 1000 + d11_idx)

            best_cost, best_D12 = sa_d12_only(n, D11, max_iter=max_iter)
            elapsed = time.time() - start
            print(f"    seed={seed}: best={best_cost} ({elapsed:.1f}s)", flush=True)

            if best_cost == 0:
                D11_list = sorted(D11)
                D12_list = sorted(int(i) for i in range(m) if best_D12[i])
                D22_list = sorted(set(range(1, m)) - set(D11_list))

                print(f"\n  *** SOLUTION FOUND! ***", flush=True)
                print(f"  D11 = {D11_list}", flush=True)
                print(f"  D12 = {D12_list}", flush=True)

                G = BlockCirculantGraph(n=n, D11=set(D11_list), D12=set(D12_list), D22=set(D22_list))
                result = verify_construction(G)
                print(f"  Verification: valid={result.valid}", flush=True)
                print(f"  max_red={result.max_red_common} (thresh {result.red_threshold})", flush=True)
                print(f"  max_blue={result.max_blue_common} (thresh {result.blue_threshold})", flush=True)

                if result.valid:
                    fn = f"/Users/stephenpadgett/Projects/math/ramsey-book-graphs/solution_n{n}.json"
                    save_construction(G, result, fn, solver="focused-d12-SA")
                    print(f"  SAVED to {fn}", flush=True)
                    return True

    print(f"\nNo solution found after trying {num_d11_to_try} D11 candidates.", flush=True)
    return False


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python focused_d12_solver.py <n> [max_seeds] [max_iter]")
        sys.exit(1)

    n = int(sys.argv[1])
    max_seeds = int(sys.argv[2]) if len(sys.argv) > 2 else 30
    max_iter = int(sys.argv[3]) if len(sys.argv) > 3 else 20000000

    solved = solve(n, max_seeds, max_iter)
    if not solved:
        print(f"\nNo solution found for n={n}")
