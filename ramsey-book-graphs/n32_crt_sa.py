"""
CRT-aware SA solver for n=32, m=63 = 9*7.

Uses the Chinese Remainder Theorem decomposition Z_63 = Z_9 x Z_7.
Elements of Z_63 map to pairs (a mod 9, a mod 7).

Key insight: For composite m, the CRT decomposition reveals
multiplicative structure that can guide construction.

QR in Z_9: {0,1,2,4,5,7,8} (squares mod 9: 0^2=0, 1^2=1, 2^2=4, 3^2=0, 4^2=7, 5^2=7, 6^2=0, 7^2=4, 8^2=1)
Actually, QR(Z_9*) = {1,4,7} (nonzero squares), QNR(Z_9*) = {2,5,8}
QR(Z_7*) = {1,2,4}, QNR(Z_7*) = {3,5,6}

CRT classes (ignoring 0):
- (QR_9, QR_7): both quadratic residues
- (QR_9, QNR_7): QR in Z_9, QNR in Z_7
- (QNR_9, QR_7): QNR in Z_9, QR in Z_7
- (QNR_9, QNR_7): both non-residues

Also need to handle (0, *) and (*, 0) elements.
"""

import sys, os, math, random, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ramsey_core import BlockCirculantGraph, verify_construction, save_construction

N_PARAM = 32
M = 63
N_VERTS = 126


def crt_decompose(x, m1=9, m2=7):
    """Map x in Z_63 to (x mod 9, x mod 7)."""
    return (x % m1, x % m2)


def crt_compose(a, b, m1=9, m2=7):
    """Map (a mod 9, b mod 7) back to Z_63."""
    # 63 = 9 * 7. Need inv(9, 7) and inv(7, 9).
    # 9 mod 7 = 2, inv(2,7) = 4 (since 2*4=8=1 mod 7)
    # 7 mod 9 = 7, inv(7,9) = 4 (since 7*4=28=1 mod 9)
    return (a * 7 * 4 + b * 9 * 4) % (m1 * m2)


def classify_elements():
    """Classify Z_63 elements by CRT quadratic residue structure."""
    # QR mod 9: squares of {1,2,3,4,5,6,7,8} mod 9 = {1,4,0,7,7,0,4,1} -> QR* = {1,4,7}
    qr9 = {1, 4, 7}
    qnr9 = {2, 5, 8}
    zero9 = {0}
    # QR mod 7: squares of {1,2,3,4,5,6} mod 7 = {1,4,2,2,4,1} -> QR* = {1,2,4}
    qr7 = {1, 2, 4}
    qnr7 = {3, 5, 6}
    zero7 = {0}

    classes = {}
    for x in range(63):
        a, b = crt_decompose(x)
        if a == 0 and b == 0:
            classes.setdefault('00', []).append(x)
        elif a == 0:
            if b in qr7:
                classes.setdefault('0R', []).append(x)
            else:
                classes.setdefault('0N', []).append(x)
        elif b == 0:
            if a in qr9:
                classes.setdefault('R0', []).append(x)
            else:
                classes.setdefault('N0', []).append(x)
        else:
            if a in qr9 and b in qr7:
                classes.setdefault('RR', []).append(x)
            elif a in qr9 and b in qnr7:
                classes.setdefault('RN', []).append(x)
            elif a in qnr9 and b in qr7:
                classes.setdefault('NR', []).append(x)
            else:
                classes.setdefault('NN', []).append(x)
    return classes


def compute_delta_array(indicator, m):
    """Compute autocorrelation via FFT."""
    a = np.array(indicator, dtype=np.float64)
    fa = np.fft.rfft(a)
    result = np.fft.irfft(np.conj(fa) * fa, n=m)
    return np.round(result).astype(np.int64)


def compute_cost(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n):
    """Compute total violation cost."""
    d11_size = int(np.sum(D11_ind))
    d12_size = int(np.sum(D12_ind))
    d22_size = m - 1 - d11_size
    d1 = d11_size + d12_size
    d2 = d22_size + d12_size
    N = 2 * m

    red_thresh = n - 2
    blue_thresh = n - 1
    cost = 0

    for d in range(1, m):
        common = int(delta_11[d] + delta_12[d])
        if D11_ind[d]:
            if common > red_thresh:
                cost += common - red_thresh
        else:
            blue_common = (N - 2) - 2 * d1 + common
            if blue_common > blue_thresh:
                cost += blue_common - blue_thresh

    v22_const = m - 2 - 2 * d11_size
    for d in range(1, m):
        d22_delta = int(delta_11[d]) + v22_const + (2 if D11_ind[d] else 0)
        common = d22_delta + int(delta_12T[d])
        if not D11_ind[d]:
            if common > red_thresh:
                cost += common - red_thresh
        else:
            blue_common = (N - 2) - 2 * d2 + common
            if blue_common > blue_thresh:
                cost += blue_common - blue_thresh

    return cost


def init_from_crt_class(classes, d11_size, d12_size, strategy, rng):
    """Initialize D11 and D12 using CRT class structure."""
    m = M

    # Build symmetric pairs
    pairs = []
    for x in range(1, m):
        neg_x = (-x) % m
        if x <= neg_x:
            pairs.append((x, neg_x))

    if strategy == 'paley_like':
        # Mimic Paley: put "QR-like" elements in D11
        # RR, RN, 0R, R0 -> 9 + 9 + 3 + 3 = 24 elements
        # NR, NN, 0N, N0 -> 9 + 9 + 3 + 3 = 24 elements
        # Need to pick symmetric subsets of these
        base = set()
        for cls in ['RR', 'RN', 'R0', '0R']:
            base.update(classes.get(cls, []))
        # Make it symmetric and the right size
        D11 = set()
        for x in base:
            if x != 0:
                D11.add(x)
                D11.add((-x) % m)
        # Trim to d11_size
        pair_list = [(min(x, (-x)%m), max(x, (-x)%m)) for x in D11 if x != 0]
        pair_set = set(pair_list)
        pair_list = list(pair_set)
        rng.shuffle(pair_list)
        D11 = set()
        for p in pair_list[:d11_size // 2]:
            D11.add(p[0])
            D11.add(p[1])

    elif strategy == 'anti_paley':
        # Put "QNR-like" elements in D11
        base = set()
        for cls in ['NR', 'NN', 'N0', '0N']:
            base.update(classes.get(cls, []))
        D11 = set()
        for x in base:
            if x != 0:
                D11.add(x)
                D11.add((-x) % m)
        pair_list = [(min(x, (-x)%m), max(x, (-x)%m)) for x in D11 if x != 0]
        pair_set = set(pair_list)
        pair_list = list(pair_set)
        rng.shuffle(pair_list)
        D11 = set()
        for p in pair_list[:d11_size // 2]:
            D11.add(p[0])
            D11.add(p[1])

    elif strategy == 'mixed':
        # Mix classes randomly
        all_nonzero = list(range(1, m))
        rng.shuffle(all_nonzero)
        D11 = set()
        for x in all_nonzero:
            if len(D11) >= d11_size:
                break
            neg_x = (-x) % m
            if x not in D11 and neg_x not in D11:
                D11.add(x)
                D11.add(neg_x)

    else:  # 'random'
        num_selected = d11_size // 2
        selected = rng.sample(range(len(pairs)), num_selected)
        D11 = set()
        for i in selected:
            D11.add(pairs[i][0])
            D11.add(pairs[i][1])

    # Initialize D12
    D12 = {0}
    remaining = list(set(range(1, m)))
    rng.shuffle(remaining)
    for x in remaining:
        if len(D12) >= d12_size:
            break
        D12.add(x)

    D11_ind = np.zeros(m, dtype=np.int64)
    for x in D11:
        D11_ind[x] = 1
    D12_ind = np.zeros(m, dtype=np.int64)
    for x in D12:
        D12_ind[x] = 1

    return D11_ind, D12_ind


def run_sa(D11_ind, D12_ind, max_iter, seed_label=""):
    """Run SA from given initial state."""
    m = M
    n = N_PARAM

    pairs = []
    for x in range(1, m):
        neg_x = (-x) % m
        if x <= neg_x:
            pairs.append((x, neg_x))
    num_pairs = len(pairs)

    D12T_ind = np.zeros(m, dtype=np.int64)
    for x in range(m):
        if D12_ind[x]:
            D12T_ind[(-x) % m] = 1

    delta_11 = compute_delta_array(D11_ind, m)
    delta_12 = compute_delta_array(D12_ind, m)
    delta_12T = compute_delta_array(D12T_ind, m)

    current_cost = compute_cost(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n)
    best_cost = current_cost
    best_D11 = np.copy(D11_ind)
    best_D12 = np.copy(D12_ind)

    T = 10.0
    T_min = 0.0001
    alpha = 1 - 6.0 / max_iter

    stuck_count = 0
    last_improvement = 0
    reheat_count = 0
    t0 = time.time()

    for it in range(max_iter):
        move_type = random.random()

        if move_type < 0.20:
            # D11 pair swap
            in_pairs = [i for i in range(num_pairs) if D11_ind[pairs[i][0]]]
            out_pairs = [i for i in range(num_pairs) if not D11_ind[pairs[i][0]]]
            if not in_pairs or not out_pairs:
                continue
            rp = random.choice(in_pairs)
            ap = random.choice(out_pairs)
            old = delta_11.copy()
            D11_ind[pairs[rp][0]] = 0
            D11_ind[pairs[rp][1]] = 0
            D11_ind[pairs[ap][0]] = 1
            D11_ind[pairs[ap][1]] = 1
            delta_11 = compute_delta_array(D11_ind, m)
            new_cost = compute_cost(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n)
            if new_cost - current_cost <= 0 or random.random() < math.exp(-(new_cost - current_cost) / max(T, T_min)):
                current_cost = new_cost
            else:
                D11_ind[pairs[rp][0]] = 1
                D11_ind[pairs[rp][1]] = 1
                D11_ind[pairs[ap][0]] = 0
                D11_ind[pairs[ap][1]] = 0
                delta_11 = old

        elif move_type < 0.80:
            # D12 single swap
            d12_others = [x for x in range(1, m) if D12_ind[x]]
            not_in = [x for x in range(1, m) if not D12_ind[x]]
            if not d12_others or not not_in:
                continue
            rem = random.choice(d12_others)
            add_el = random.choice(not_in)
            old_d12 = delta_12.copy()
            old_d12T = delta_12T.copy()
            D12_ind[rem] = 0
            D12_ind[add_el] = 1
            D12T_ind[(-rem) % m] = 0
            D12T_ind[(-add_el) % m] = 1
            delta_12 = compute_delta_array(D12_ind, m)
            delta_12T = compute_delta_array(D12T_ind, m)
            new_cost = compute_cost(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n)
            if new_cost - current_cost <= 0 or random.random() < math.exp(-(new_cost - current_cost) / max(T, T_min)):
                current_cost = new_cost
            else:
                D12_ind[rem] = 1
                D12_ind[add_el] = 0
                D12T_ind[(-rem) % m] = 1
                D12T_ind[(-add_el) % m] = 0
                delta_12 = old_d12
                delta_12T = old_d12T

        else:
            # Multi-swap: 2 D12 elements
            d12_others = [x for x in range(1, m) if D12_ind[x]]
            not_in = [x for x in range(1, m) if not D12_ind[x]]
            if len(d12_others) < 2 or len(not_in) < 2:
                continue
            rems = random.sample(d12_others, 2)
            adds = random.sample(not_in, 2)
            old_d12 = delta_12.copy()
            old_d12T = delta_12T.copy()
            for r in rems:
                D12_ind[r] = 0
                D12T_ind[(-r) % m] = 0
            for a in adds:
                D12_ind[a] = 1
                D12T_ind[(-a) % m] = 1
            delta_12 = compute_delta_array(D12_ind, m)
            delta_12T = compute_delta_array(D12T_ind, m)
            new_cost = compute_cost(D11_ind, D12_ind, delta_11, delta_12, delta_12T, m, n)
            if new_cost - current_cost <= 0 or random.random() < math.exp(-(new_cost - current_cost) / max(T, T_min)):
                current_cost = new_cost
            else:
                for r in rems:
                    D12_ind[r] = 1
                    D12T_ind[(-r) % m] = 1
                for a in adds:
                    D12_ind[a] = 0
                    D12T_ind[(-a) % m] = 0
                delta_12 = old_d12
                delta_12T = old_d12T

        if current_cost < best_cost:
            best_cost = current_cost
            best_D11 = np.copy(D11_ind)
            best_D12 = np.copy(D12_ind)
            last_improvement = it
            if best_cost == 0:
                break

        # Reheat if stuck
        if it - last_improvement > 400000 and best_cost > 0 and best_cost <= 12:
            T = 4.0
            last_improvement = it
            reheat_count += 1
            D11_ind = np.copy(best_D11)
            D12_ind = np.copy(best_D12)
            D12T_ind = np.zeros(m, dtype=np.int64)
            for x in range(m):
                if D12_ind[x]:
                    D12T_ind[(-x) % m] = 1
            delta_11 = compute_delta_array(D11_ind, m)
            delta_12 = compute_delta_array(D12_ind, m)
            delta_12T = compute_delta_array(D12T_ind, m)
            current_cost = best_cost

        T = max(T * alpha, T_min)

        if it % 2000000 == 0 and it > 0:
            elapsed = time.time() - t0
            print(f"    iter {it//1000000}M: cost={current_cost}, best={best_cost}, T={T:.6f}, reheats={reheat_count} ({elapsed:.0f}s)", flush=True)

    elapsed = time.time() - t0
    return best_cost, best_D11, best_D12, elapsed, reheat_count


def main():
    max_seeds = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    max_iter = int(sys.argv[2]) if len(sys.argv) > 2 else 10000000

    classes = classify_elements()
    print(f"CRT classes for Z_63 = Z_9 x Z_7:")
    for cls_name, elements in sorted(classes.items()):
        print(f"  {cls_name}: {len(elements)} elements = {sorted(elements)}")

    d12_size = N_PARAM - 1  # 31

    # Verify CRT compose/decompose
    for x in range(63):
        a, b = crt_decompose(x)
        assert crt_compose(a, b) == x, f"CRT failed for {x}"

    print(f"\nn={N_PARAM}, m={M}, N={N_VERTS}")
    print(f"Max seeds: {max_seeds}, Max iter: {max_iter}")
    print(flush=True)

    strategies = ['paley_like', 'anti_paley', 'mixed', 'random']
    d11_sizes = [32, 30]

    global_best_cost = 999
    t0 = time.time()

    for d11_size in d11_sizes:
        for strat in strategies:
            print(f"\n--- |D11|={d11_size}, strategy={strat} ---", flush=True)
            seeds_per_strat = max_seeds // (len(strategies) * len(d11_sizes))

            for seed in range(seeds_per_strat):
                rng = random.Random(seed * 1000 + hash(strat) % 1000)
                D11_ind, D12_ind = init_from_crt_class(classes, d11_size, d12_size, strat, rng)
                random.seed(seed + hash(strat))

                best_cost, best_D11, best_D12, elapsed, reheats = run_sa(
                    D11_ind, D12_ind, max_iter
                )
                total_elapsed = time.time() - t0
                print(f"  seed={seed}: best={best_cost}, reheats={reheats} ({total_elapsed:.0f}s)", flush=True)

                if best_cost < global_best_cost:
                    global_best_cost = best_cost

                if best_cost == 0:
                    D11_list = sorted([int(i) for i in range(M) if best_D11[i]])
                    D12_list = sorted([int(i) for i in range(M) if best_D12[i]])
                    print(f"\n  *** SOLUTION FOUND! ***")
                    print(f"  D11 = {D11_list}")
                    print(f"  D12 = {D12_list}")
                    D11_s = set(D11_list)
                    D12_s = set(D12_list)
                    D22_s = set(range(1, M)) - D11_s
                    G = BlockCirculantGraph(n=N_PARAM, D11=D11_s, D12=D12_s, D22=D22_s)
                    result = verify_construction(G)
                    print(f"  Full verification: valid={result.valid}")
                    print(f"  max_red={result.max_red_common} (thresh {result.red_threshold})")
                    print(f"  max_blue={result.max_blue_common} (thresh {result.blue_threshold})")
                    if result.valid:
                        fn = f"/Users/stephenpadgett/Projects/math/ramsey-book-graphs/solution_n{N_PARAM}.json"
                        save_construction(G, result, fn, solver="CRT-SA")
                        print(f"  SAVED to {fn}")
                        sys.exit(0)

    print(f"\nNo solution found. Global best cost = {global_best_cost}")


if __name__ == "__main__":
    main()
