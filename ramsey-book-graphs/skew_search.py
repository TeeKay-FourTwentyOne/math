"""
Skew-Hadamard Search for R(B_21, B_22) = 87

Exploits the Paley tournament structure for p = 43 = 3 (mod 4).
D12 (cross-block edges) doesn't need symmetry, so QR_43 is a natural fit.

Key algebraic fact: Delta(QR_43, QR_43, d) = (p-3)/4 = 10 for all d != 0.
"""

import sys, os, time, random, math
from typing import Set, Tuple, Optional

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ramsey_core import BlockCirculantGraph, verify_construction, compute_violation_cost

# === Constants ===
N_PARAM = 22
M = 2 * N_PARAM - 1  # 43
N_VERTS = 4 * N_PARAM - 2  # 86
HALF_M = M // 2  # 21

# === Quadratic residues mod 43 ===
QR = {pow(x, 2, M) for x in range(1, M)}
QNR = set(range(1, M)) - QR

# Symmetric pairs {d, M-d} for d = 1,...,21
PAIRS = [(d, M - d) for d in range(1, HALF_M + 1)]


def pairs_to_set(pair_indices):
    result = set()
    for i in pair_indices:
        result.add(PAIRS[i][0])
        result.add(PAIRS[i][1])
    return result


def evaluate(D11, D12, D22):
    G = BlockCirculantGraph(n=N_PARAM, D11=D11, D12=D12, D22=D22)
    result = verify_construction(G)
    cost = sum(e for _, _, e in result.violations)
    return len(result.violations), cost, result


def print_result(label, D11, D12, D22):
    nviol, cost, result = evaluate(D11, D12, D22)
    d1 = len(D11) + len(D12)
    d2 = len(D22) + len(D12)
    print(f"  {label}: {nviol} violations, cost={cost}, "
          f"|D11|={len(D11)}, |D12|={len(D12)}, |D22|={len(D22)}, "
          f"d1={d1}, d2={d2}")
    if 0 < nviol <= 16:
        for vtype, d, excess in result.violations:
            print(f"    {vtype} d={d} excess={excess}")
    return nviol, cost


# =====================================================================
# PART 1: Quick algebraic exploration
# =====================================================================
def explore_algebraic():
    print("=" * 70)
    print("PART 1: Algebraic Exploration")
    print("=" * 70)

    print(f"\nQR_43 = {sorted(QR)}")
    print(f"|QR| = {len(QR)}, |QNR| = {len(QNR)}")

    # Pairs where the small element (d <= 21) is in QR
    qr_aligned = [i for i in range(HALF_M) if PAIRS[i][0] in QR]
    qnr_aligned = [i for i in range(HALF_M) if PAIRS[i][0] in QNR]
    print(f"QR-aligned pairs: {len(qr_aligned)}, QNR-aligned: {len(qnr_aligned)}")

    best_nviol = 999
    best_config = None
    best_label = ""

    def try_config(label, D11, D12, D22):
        nonlocal best_nviol, best_config, best_label
        nv, cost = print_result(label, D11, D12, D22)
        if nv < best_nviol or (nv == best_nviol and cost < best_config[1]):
            best_nviol = nv
            best_config = (D11.copy(), D12.copy(), D22.copy(), cost)
            best_label = label
        return nv

    # --- A: D12 = QR, D22 = complement(D11) ---
    print("\n--- A: D12=QR, D22=complement(D11) ---")
    for np11 in range(7, 16):
        sel = qr_aligned[:min(np11, len(qr_aligned))]
        if len(sel) < np11:
            sel = sel + qnr_aligned[:np11 - len(sel)]
        D11 = pairs_to_set(sel)
        D22 = set(range(1, M)) - D11
        try_config(f"A np={np11:2d}", D11, QR.copy(), D22)

    # --- B: D12 = QR, D11 = D22 ---
    print("\n--- B: D12=QR, D11=D22 ---")
    for np in range(7, 16):
        D11 = pairs_to_set(list(range(np)))
        try_config(f"B np={np:2d}", D11, QR.copy(), D11.copy())

    # --- C: D12 = QR, independent D11/D22, various densities ---
    print("\n--- C: D12=QR, independent D11/D22 ---")
    for np11 in [9, 10, 11, 12]:
        for np22 in [9, 10, 11, 12]:
            D11 = pairs_to_set(qr_aligned[:np11] if np11 <= len(qr_aligned)
                               else list(range(np11)))
            D22 = pairs_to_set(qnr_aligned[:np22] if np22 <= len(qnr_aligned)
                               else list(range(np22)))
            try_config(f"C {np11:2d}/{np22:2d}", D11, QR.copy(), D22)

    # --- D: D12 = QNR ---
    print("\n--- D: D12=QNR, complement D22 ---")
    for np11 in range(8, 14):
        sel = qr_aligned[:min(np11, len(qr_aligned))]
        if len(sel) < np11:
            sel = sel + qnr_aligned[:np11 - len(sel)]
        D11 = pairs_to_set(sel)
        D22 = set(range(1, M)) - D11
        try_config(f"D np={np11:2d}", D11, QNR.copy(), D22)

    # --- E: D12 = QR + {0} ---
    print("\n--- E: D12=QR+{0} ---")
    D12e = QR | {0}
    for np11 in range(8, 14):
        D11 = pairs_to_set(list(range(np11)))
        D22 = set(range(1, M)) - D11
        try_config(f"E np={np11:2d}", D11, D12e, D22)

    # --- F: D12 = QNR + {0} ---
    print("\n--- F: D12=QNR+{0} ---")
    D12f = QNR | {0}
    for np11 in range(8, 14):
        D11 = pairs_to_set(list(range(np11)))
        D22 = set(range(1, M)) - D11
        try_config(f"F np={np11:2d}", D11, D12f, D22)

    print(f"\n*** Best algebraic: {best_label} with {best_nviol} violations ***")
    return best_config


# =====================================================================
# PART 2: Local search with fixed D12
# =====================================================================
def local_search_fixed_d12(D12, D11_init, D22_init,
                           max_iter=200000, n_trials=10, label=""):
    print(f"\n{'='*70}")
    print(f"PART 2: Local search, D12 fixed ({label}), |D12|={len(D12)}")
    print(f"{'='*70}")

    best_global = float('inf')
    best_config = None

    for trial in range(n_trials):
        if trial == 0:
            D11 = D11_init.copy()
            D22 = D22_init.copy()
        else:
            # Random perturbation of init
            p11 = list(range(HALF_M))
            p22 = list(range(HALF_M))
            random.shuffle(p11)
            random.shuffle(p22)
            n11 = random.randint(8, 13)
            n22 = random.randint(8, 13)
            D11 = pairs_to_set(p11[:n11])
            D22 = pairs_to_set(p22[:n22])

        G = BlockCirculantGraph(n=N_PARAM, D11=D11, D12=D12, D22=D22)
        current_cost = compute_violation_cost(G)
        best_cost = current_cost
        best_D11 = D11.copy()
        best_D22 = D22.copy()

        temp = 5.0
        cooling = 0.99997

        for it in range(max_iter):
            if current_cost == 0:
                print(f"  SOLUTION FOUND! trial={trial} iter={it}")
                return D11, D12, D22

            # Flip one pair in D11 or D22
            if random.random() < 0.5:
                pi = random.randint(0, HALF_M - 1)
                d, md = PAIRS[pi]
                new_D11 = D11.copy()
                new_D11.symmetric_difference_update({d, md})
                new_D22 = D22
            else:
                pi = random.randint(0, HALF_M - 1)
                d, md = PAIRS[pi]
                new_D22 = D22.copy()
                new_D22.symmetric_difference_update({d, md})
                new_D11 = D11

            G2 = BlockCirculantGraph(n=N_PARAM, D11=new_D11, D12=D12, D22=new_D22)
            new_cost = compute_violation_cost(G2)

            delta = new_cost - current_cost
            if delta < 0 or random.random() < math.exp(-delta / max(temp, 0.01)):
                D11 = new_D11
                D22 = new_D22
                current_cost = new_cost
                if current_cost < best_cost:
                    best_cost = current_cost
                    best_D11 = D11.copy()
                    best_D22 = D22.copy()

            temp *= cooling

            if it % 50000 == 0:
                print(f"  t{trial} i{it}: cur={current_cost} best={best_cost} T={temp:.3f}")

        print(f"  Trial {trial}: best={best_cost}")
        if best_cost < best_global:
            best_global = best_cost
            best_config = (best_D11.copy(), D12.copy(), best_D22.copy())

    if best_config:
        print(f"\nBest: cost={best_global}")
        print_result("Best fixed-D12", *best_config)
    return best_config


# =====================================================================
# PART 3: Full local search (D12 also free)
# =====================================================================
def local_search_full(max_iter=300000, n_trials=8):
    print(f"\n{'='*70}")
    print(f"PART 3: Full local search (D11, D12, D22 all free)")
    print(f"{'='*70}")

    best_global = float('inf')
    best_config = None

    for trial in range(n_trials):
        # Initialize
        if trial < 2:
            D12 = QR.copy()
        elif trial < 4:
            D12 = QNR.copy()
        elif trial < 6:
            D12 = QR | {0}
        else:
            D12 = set(random.sample(range(M), random.randint(19, 24)))

        n11 = random.randint(8, 13)
        n22 = random.randint(8, 13)
        p11 = random.sample(range(HALF_M), n11)
        p22 = random.sample(range(HALF_M), n22)
        D11 = pairs_to_set(p11)
        D22 = pairs_to_set(p22)

        G = BlockCirculantGraph(n=N_PARAM, D11=D11, D12=D12, D22=D22)
        current_cost = compute_violation_cost(G)
        best_cost = current_cost
        best_cfg = (D11.copy(), D12.copy(), D22.copy())

        temp = 8.0
        cooling = 0.99998

        for it in range(max_iter):
            if current_cost == 0:
                print(f"  SOLUTION FOUND! trial={trial} iter={it}")
                return D11, D12, D22

            r = random.random()
            new_D11, new_D12, new_D22 = D11, D12, D22

            if r < 0.35:
                pi = random.randint(0, HALF_M - 1)
                d, md = PAIRS[pi]
                new_D11 = D11.copy()
                new_D11.symmetric_difference_update({d, md})
            elif r < 0.70:
                pi = random.randint(0, HALF_M - 1)
                d, md = PAIRS[pi]
                new_D22 = D22.copy()
                new_D22.symmetric_difference_update({d, md})
            else:
                d = random.randint(0, M - 1)
                new_D12 = D12.copy()
                new_D12.symmetric_difference_update({d})

            G2 = BlockCirculantGraph(n=N_PARAM, D11=new_D11, D12=new_D12, D22=new_D22)
            new_cost = compute_violation_cost(G2)

            delta = new_cost - current_cost
            if delta < 0 or random.random() < math.exp(-delta / max(temp, 0.01)):
                D11, D12, D22 = new_D11, new_D12, new_D22
                current_cost = new_cost
                if current_cost < best_cost:
                    best_cost = current_cost
                    best_cfg = (D11.copy(), D12.copy(), D22.copy())

            temp *= cooling

            if it % 100000 == 0:
                print(f"  t{trial} i{it}: cur={current_cost} best={best_cost} T={temp:.3f}")

        print(f"  Trial {trial}: best={best_cost}")
        if best_cost < best_global:
            best_global = best_cost
            best_config = best_cfg

    if best_config:
        print(f"\nBest full search: cost={best_global}")
        print_result("Best full", *best_config)
    return best_config


# =====================================================================
# MAIN
# =====================================================================
if __name__ == "__main__":
    random.seed(42)
    start = time.time()

    # Part 1: Quick algebraic scan
    alg_result = explore_algebraic()
    print(f"\nPart 1 elapsed: {time.time()-start:.1f}s")

    # Part 2: Local search with best fixed D12
    if alg_result:
        D11, D12, D22, _ = alg_result
        t2 = time.time()
        ls_result = local_search_fixed_d12(
            D12, D11, D22, max_iter=200000, n_trials=8, label="best-alg")
        print(f"Part 2 elapsed: {time.time()-t2:.1f}s")

    # Also try QR and QNR as D12
    for d12_set, d12_label in [(QR, "QR"), (QNR, "QNR")]:
        t2 = time.time()
        init_D11 = pairs_to_set(list(range(10)))
        init_D22 = pairs_to_set(list(range(10)))
        local_search_fixed_d12(
            d12_set.copy(), init_D11, init_D22,
            max_iter=200000, n_trials=5, label=d12_label)
        print(f"  {d12_label} search elapsed: {time.time()-t2:.1f}s")

    # Part 3: Full search
    t3 = time.time()
    full_result = local_search_full(max_iter=300000, n_trials=8)
    print(f"Part 3 elapsed: {time.time()-t3:.1f}s")

    print(f"\nTotal elapsed: {time.time()-start:.1f}s")
