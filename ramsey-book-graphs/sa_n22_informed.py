"""
Informed SA search for R(B_21, B_22) = 87.

Uses structural insights from known m = 3 (mod 4) constructions:
- D22 = complement(D11) (universal)
- 0 in D12 (universal)
- |D11| = 22, |D12| = 21 => d1 = 43, d2 = 41 (primary)
- OR |D11| = 20, |D12| = 21 => d1 = 41, d2 = 43 (flipped)

Improvements over prior SA:
1. Fix cardinalities (swap moves only, maintaining |D11| and |D12|)
2. Fix 0 in D12
3. Track per-constraint costs for smarter moves
4. Multi-temperature parallel runs
"""

import sys, os, time, random, math
from typing import Set, List, Tuple, Dict
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ramsey_core import (
    BlockCirculantGraph, verify_construction, Delta, Sigma,
    make_symmetric
)

M = 43
N = 86
HALF_M = M // 2  # 21
N_PARAM = 22
RED_THRESH = 20
BLUE_THRESH = 21

# Symmetric pairs
PAIRS = [(d, M - d) for d in range(1, HALF_M + 1)]


class FastEvaluator:
    """Fast incremental evaluation of violation cost."""

    def __init__(self, D11_half: Set[int], D12: Set[int]):
        """Initialize with D11 half-set (pair indices 0..20) and D12."""
        self.D11_half = set(D11_half)
        self.D12 = set(D12)

        # Build full sets
        self._rebuild()

    def _rebuild(self):
        """Rebuild full D11, D22 and precompute all lambda values."""
        self.D11 = set()
        for i in self.D11_half:
            d, md = PAIRS[i]
            self.D11.add(d)
            self.D11.add(md)

        self.D22 = set(range(1, M)) - self.D11
        self.D12T = {(-x) % M for x in self.D12}

        self.d1 = len(self.D11) + len(self.D12)
        self.d2 = len(self.D22) + len(self.D12)

        self.blue_V1V1 = 2 * self.d1 - (N - 2) + BLUE_THRESH
        self.blue_V2V2 = 2 * self.d2 - (N - 2) + BLUE_THRESH
        self.blue_V1V2 = self.d1 + self.d2 - (N - 2) + BLUE_THRESH

        # Precompute all lambda_red values
        self.lam_V1V1 = {}
        for d in range(1, M):
            self.lam_V1V1[d] = (
                Delta(self.D11, self.D11, d, M) +
                Delta(self.D12, self.D12, d, M)
            )

        self.lam_V2V2 = {}
        for d in range(1, M):
            self.lam_V2V2[d] = (
                Delta(self.D22, self.D22, d, M) +
                Delta(self.D12T, self.D12T, d, M)
            )

        self.lam_V1V2 = {}
        for d in range(M):
            self.lam_V1V2[d] = (
                Sigma(self.D11, self.D12, d, M) +
                Delta(self.D12, self.D22, d, M)
            )

    def compute_cost(self) -> Tuple[int, int, List]:
        """Compute total cost, violation count, and violation list."""
        cost = 0
        viols = []

        for d in range(1, M):
            lam = self.lam_V1V1[d]
            if d in self.D11:
                if lam > RED_THRESH:
                    excess = lam - RED_THRESH
                    cost += excess
                    viols.append(("V1V1_red", d, excess))
            else:
                lam_blue = (N - 2) - 2 * self.d1 + lam
                if lam_blue > BLUE_THRESH:
                    excess = lam_blue - BLUE_THRESH
                    cost += excess
                    viols.append(("V1V1_blue", d, excess))

        for d in range(1, M):
            lam = self.lam_V2V2[d]
            if d in self.D22:
                if lam > RED_THRESH:
                    excess = lam - RED_THRESH
                    cost += excess
                    viols.append(("V2V2_red", d, excess))
            else:
                lam_blue = (N - 2) - 2 * self.d2 + lam
                if lam_blue > BLUE_THRESH:
                    excess = lam_blue - BLUE_THRESH
                    cost += excess
                    viols.append(("V2V2_blue", d, excess))

        for d in range(M):
            lam = self.lam_V1V2[d]
            if d in self.D12:
                if lam > RED_THRESH:
                    excess = lam - RED_THRESH
                    cost += excess
                    viols.append(("V1V2_red", d, excess))
            else:
                lam_blue = (N - 2) - self.d1 - self.d2 + lam
                if lam_blue > BLUE_THRESH:
                    excess = lam_blue - BLUE_THRESH
                    cost += excess
                    viols.append(("V1V2_blue", d, excess))

        return cost, len(viols), viols


def sa_search(d11_target: int, d12_target: int,
              max_iter: int = 2000000, n_trials: int = 16,
              T_init: float = 5.0, cooling: float = 0.999995):
    """
    Run SA with fixed cardinalities using swap moves.

    Swap moves maintain set sizes:
    - D11 swap: move one pair in, one pair out
    - D12 swap: move one element in, one element out
    """
    half_d11 = d11_target // 2  # pairs in D11
    best_global_cost = float('inf')
    best_global_config = None
    best_global_viols = None

    for trial in range(n_trials):
        # Random initialization with exact cardinalities
        all_pairs = list(range(HALF_M))
        random.shuffle(all_pairs)
        D11_half = set(all_pairs[:half_d11])

        # D12: always include 0, then random selection for rest
        d12_rest = list(range(1, M))
        random.shuffle(d12_rest)
        D12 = {0} | set(d12_rest[:d12_target - 1])

        ev = FastEvaluator(D11_half, D12)
        cost, nviols, viols = ev.compute_cost()
        best_cost = cost
        best_D11_half = D11_half.copy()
        best_D12 = D12.copy()
        best_viols_local = viols

        temp = T_init

        for it in range(max_iter):
            if cost == 0:
                print(f"  SOLUTION FOUND at trial={trial}, iter={it}")
                ev._rebuild()
                return ev.D11, ev.D12, ev.D22, viols

            # Choose move type
            r = random.random()
            if r < 0.5:
                # D11 swap: move one pair out, one pair in
                out_list = list(D11_half)
                in_list = list(set(range(HALF_M)) - D11_half)
                if not out_list or not in_list:
                    continue
                pi_out = random.choice(out_list)
                pi_in = random.choice(in_list)

                D11_half.discard(pi_out)
                D11_half.add(pi_in)
                ev.D11_half = D11_half
                ev._rebuild()
                new_cost, new_nviols, new_viols = ev.compute_cost()

                delta = new_cost - cost
                if delta < 0 or random.random() < math.exp(-delta / max(temp, 0.001)):
                    cost = new_cost
                    if cost < best_cost:
                        best_cost = cost
                        best_D11_half = D11_half.copy()
                        best_D12 = D12.copy()
                        best_viols_local = new_viols
                else:
                    # Undo
                    D11_half.discard(pi_in)
                    D11_half.add(pi_out)
                    ev.D11_half = D11_half
                    ev._rebuild()

            else:
                # D12 swap: move one element out, one element in (keep 0)
                out_list = [x for x in D12 if x != 0]
                in_list = list(set(range(1, M)) - D12)
                if not out_list or not in_list:
                    continue
                d_out = random.choice(out_list)
                d_in = random.choice(in_list)

                D12.discard(d_out)
                D12.add(d_in)
                ev.D12 = D12
                ev._rebuild()
                new_cost, new_nviols, new_viols = ev.compute_cost()

                delta = new_cost - cost
                if delta < 0 or random.random() < math.exp(-delta / max(temp, 0.001)):
                    cost = new_cost
                    if cost < best_cost:
                        best_cost = cost
                        best_D11_half = D11_half.copy()
                        best_D12 = D12.copy()
                        best_viols_local = new_viols
                else:
                    # Undo
                    D12.discard(d_in)
                    D12.add(d_out)
                    ev.D12 = D12
                    ev._rebuild()

            temp *= cooling

            if it % 200000 == 0:
                print(f"  t{trial} i{it}: cost={cost} best={best_cost} T={temp:.4f} "
                      f"d1={ev.d1} d2={ev.d2}")
                sys.stdout.flush()

        print(f"  Trial {trial}: best={best_cost} ({len(best_viols_local)} viols)")
        sys.stdout.flush()

        if best_cost < best_global_cost:
            best_global_cost = best_cost
            # Restore best configuration
            D11_half = best_D11_half.copy()
            D12 = best_D12.copy()
            ev = FastEvaluator(D11_half, D12)
            best_global_config = (ev.D11.copy(), ev.D12.copy(), ev.D22.copy())
            best_global_viols = best_viols_local

    return best_global_config, best_global_cost, best_global_viols


def main():
    random.seed(42)

    print("=" * 70)
    print("Informed SA for R(B_21, B_22) = 87")
    print("Based on universal pattern: |D11|=22, |D12|=21, D22=complement")
    print("=" * 70)

    # Primary config: |D11|=22, |D12|=21
    print("\n--- Config A: |D11|=22, |D12|=21 (d1=43, d2=41) ---")
    start = time.time()
    result = sa_search(22, 21, max_iter=2000000, n_trials=16,
                       T_init=5.0, cooling=0.999997)

    if isinstance(result, tuple) and len(result) == 4:
        D11, D12, D22, viols = result
        print(f"\n*** SOLUTION FOUND! ***")
        G = BlockCirculantGraph(n=N_PARAM, D11=D11, D12=D12, D22=D22)
        v = verify_construction(G)
        print(f"Valid: {v.valid}")
        print(f"D11={sorted(D11)}")
        print(f"D12={sorted(D12)}")
    else:
        config, cost, viols = result
        print(f"\nBest config A: cost={cost}")
        if config:
            D11, D12, D22 = config
            print(f"|D11|={len(D11)}, |D12|={len(D12)}, |D22|={len(D22)}")
            print(f"d1={len(D11)+len(D12)}, d2={len(D22)+len(D12)}")
        if viols:
            print(f"Violations ({len(viols)}):")
            for vt, d, ex in viols[:20]:
                print(f"  {vt} d={d} excess={ex}")

    print(f"Config A elapsed: {time.time()-start:.1f}s")

    # Flipped config: |D11|=20, |D12|=21
    print("\n--- Config B: |D11|=20, |D12|=21 (d1=41, d2=43) ---")
    start = time.time()
    result = sa_search(20, 21, max_iter=2000000, n_trials=16,
                       T_init=5.0, cooling=0.999997)

    if isinstance(result, tuple) and len(result) == 4:
        D11, D12, D22, viols = result
        print(f"\n*** SOLUTION FOUND! ***")
        G = BlockCirculantGraph(n=N_PARAM, D11=D11, D12=D12, D22=D22)
        v = verify_construction(G)
        print(f"Valid: {v.valid}")
        print(f"D11={sorted(D11)}")
        print(f"D12={sorted(D12)}")
    else:
        config, cost, viols = result
        print(f"\nBest config B: cost={cost}")
        if config:
            D11, D12, D22 = config
            print(f"|D11|={len(D11)}, |D12|={len(D12)}, |D22|={len(D22)}")
            print(f"d1={len(D11)+len(D12)}, d2={len(D22)+len(D12)}")
        if viols:
            print(f"Violations ({len(viols)}):")
            for vt, d, ex in viols[:20]:
                print(f"  {vt} d={d} excess={ex}")

    print(f"Config B elapsed: {time.time()-start:.1f}s")


if __name__ == "__main__":
    main()
