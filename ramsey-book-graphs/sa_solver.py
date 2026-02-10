"""
Generalized Simulated Annealing solver for R(B_{n-1}, B_n) = 4n - 1.

Parameterized for arbitrary n. Uses structural insights:
- D22 = complement(D11) (universal)
- 0 in D12 (universal)
- Fixed cardinalities with swap moves
- Multi-trial parallel runs

Based on sa_n22_informed.py which successfully solved n=22.

Usage:
    python sa_solver.py --n 24 --trials 16 --max-iter 2000000
    python sa_solver.py --n 24 --config both --trials 32
"""

import sys, os, time, random, math, json, argparse
from typing import Set, List, Tuple, Optional

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ramsey_core import (
    BlockCirculantGraph, verify_construction, save_construction,
    Delta, Sigma, make_symmetric
)


class FastEvaluator:
    """Fast incremental evaluation of violation cost for arbitrary n."""

    def __init__(self, n: int, D11_half: Set[int], D12: Set[int]):
        self.n = n
        self.m = 2 * n - 1
        self.N = 4 * n - 2
        self.half_m = self.m // 2
        self.red_thresh = n - 2
        self.blue_thresh = n - 1
        self.pairs = [(d, self.m - d) for d in range(1, self.half_m + 1)]

        self.D11_half = set(D11_half)
        self.D12 = set(D12)
        self._rebuild()

    def _rebuild(self):
        """Rebuild full D11, D22 and precompute all lambda values."""
        self.D11 = set()
        for i in self.D11_half:
            d, md = self.pairs[i]
            self.D11.add(d)
            self.D11.add(md)

        self.D22 = set(range(1, self.m)) - self.D11
        self.D12T = {(-x) % self.m for x in self.D12}

        self.d1 = len(self.D11) + len(self.D12)
        self.d2 = len(self.D22) + len(self.D12)

        m = self.m
        N = self.N

        # Precompute all lambda_red values
        self.lam_V1V1 = {}
        for d in range(1, m):
            self.lam_V1V1[d] = (
                Delta(self.D11, self.D11, d, m) +
                Delta(self.D12, self.D12, d, m)
            )

        self.lam_V2V2 = {}
        for d in range(1, m):
            self.lam_V2V2[d] = (
                Delta(self.D22, self.D22, d, m) +
                Delta(self.D12T, self.D12T, d, m)
            )

        self.lam_V1V2 = {}
        for d in range(m):
            self.lam_V1V2[d] = (
                Sigma(self.D11, self.D12, d, m) +
                Delta(self.D12, self.D22, d, m)
            )

    def compute_cost(self) -> Tuple[int, int, List]:
        """Compute total cost, violation count, and violation list."""
        cost = 0
        viols = []
        m = self.m
        N = self.N
        red_thresh = self.red_thresh
        blue_thresh = self.blue_thresh

        for d in range(1, m):
            lam = self.lam_V1V1[d]
            if d in self.D11:
                if lam > red_thresh:
                    excess = lam - red_thresh
                    cost += excess
                    viols.append(("V1V1_red", d, excess))
            else:
                lam_blue = (N - 2) - 2 * self.d1 + lam
                if lam_blue > blue_thresh:
                    excess = lam_blue - blue_thresh
                    cost += excess
                    viols.append(("V1V1_blue", d, excess))

        for d in range(1, m):
            lam = self.lam_V2V2[d]
            if d in self.D22:
                if lam > red_thresh:
                    excess = lam - red_thresh
                    cost += excess
                    viols.append(("V2V2_red", d, excess))
            else:
                lam_blue = (N - 2) - 2 * self.d2 + lam
                if lam_blue > blue_thresh:
                    excess = lam_blue - blue_thresh
                    cost += excess
                    viols.append(("V2V2_blue", d, excess))

        for d in range(m):
            lam = self.lam_V1V2[d]
            if d in self.D12:
                if lam > red_thresh:
                    excess = lam - red_thresh
                    cost += excess
                    viols.append(("V1V2_red", d, excess))
            else:
                lam_blue = (N - 2) - self.d1 - self.d2 + lam
                if lam_blue > blue_thresh:
                    excess = lam_blue - blue_thresh
                    cost += excess
                    viols.append(("V1V2_blue", d, excess))

        return cost, len(viols), viols


def sa_search(n: int, d11_target: int, d12_target: int,
              max_iter: int = 2000000, n_trials: int = 16,
              T_init: float = 5.0, cooling: float = 0.999995,
              seed: Optional[int] = None, verbose: bool = True,
              reheat_interval: int = 0, reheat_factor: float = 0.5) -> Tuple:
    """
    Run SA with fixed cardinalities using swap moves.

    Args:
        reheat_interval: If > 0, reheat every this many iterations when stuck.
        reheat_factor: Fraction of T_init to reheat to.

    Returns:
        On success: (D11, D12, D22, []) as a 4-tuple
        On failure: ((D11, D12, D22), best_cost, viols) as a 3-tuple
    """
    m = 2 * n - 1
    half_m = m // 2

    if seed is not None:
        random.seed(seed)

    half_d11 = d11_target // 2  # pairs in D11

    best_global_cost = float('inf')
    best_global_config = None
    best_global_viols = None

    for trial in range(n_trials):
        # Random initialization with exact cardinalities
        all_pairs = list(range(half_m))
        random.shuffle(all_pairs)
        D11_half = set(all_pairs[:half_d11])

        # D12: always include 0, then random selection for rest
        d12_rest = list(range(1, m))
        random.shuffle(d12_rest)
        D12 = {0} | set(d12_rest[:d12_target - 1])

        ev = FastEvaluator(n, D11_half, D12)
        cost, nviols, viols = ev.compute_cost()
        best_cost = cost
        best_D11_half = D11_half.copy()
        best_D12 = D12.copy()
        best_viols_local = viols

        temp = T_init
        stall_counter = 0
        last_best = best_cost

        for it in range(max_iter):
            if cost == 0:
                if verbose:
                    print(f"  SOLUTION FOUND at trial={trial}, iter={it}")
                ev._rebuild()
                return ev.D11, ev.D12, ev.D22, viols

            # Reheat if stuck
            if reheat_interval > 0 and it > 0 and it % reheat_interval == 0:
                if best_cost >= last_best:
                    stall_counter += 1
                    reheat_temp = T_init * reheat_factor ** stall_counter
                    if reheat_temp > temp:
                        temp = reheat_temp
                        if verbose:
                            print(f"  t{trial} i{it}: REHEAT to T={temp:.4f} (stall#{stall_counter})")
                else:
                    stall_counter = 0
                last_best = best_cost

            # Choose move type
            r = random.random()
            if r < 0.5:
                # D11 swap: move one pair out, one pair in
                out_list = list(D11_half)
                in_list = list(set(range(half_m)) - D11_half)
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
                    D11_half.discard(pi_in)
                    D11_half.add(pi_out)
                    ev.D11_half = D11_half
                    ev._rebuild()

            else:
                # D12 swap: move one element out, one element in (keep 0)
                out_list = [x for x in D12 if x != 0]
                in_list = list(set(range(1, m)) - D12)
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
                    D12.discard(d_in)
                    D12.add(d_out)
                    ev.D12 = D12
                    ev._rebuild()

            temp *= cooling

            if verbose and it % 200000 == 0:
                print(f"  t{trial} i{it}: cost={cost} best={best_cost} T={temp:.4f} "
                      f"d1={ev.d1} d2={ev.d2}")
                sys.stdout.flush()

        if verbose:
            print(f"  Trial {trial}: best={best_cost} ({len(best_viols_local)} viols)")
            sys.stdout.flush()

        if best_cost < best_global_cost:
            best_global_cost = best_cost
            D11_half = best_D11_half.copy()
            D12 = best_D12.copy()
            ev = FastEvaluator(n, D11_half, D12)
            best_global_config = (ev.D11.copy(), ev.D12.copy(), ev.D22.copy())
            best_global_viols = best_viols_local

    return best_global_config, best_global_cost, best_global_viols


def compute_cardinalities(n: int):
    """
    Compute primary and flipped cardinalities for given n.

    Primary config: |D11| = (m+1)/2, |D12| = (m-1)/2
    Flipped config: |D11| = (m-1)/2, |D12| = (m-1)/2 with swapped degrees

    Returns list of (label, d11_size, d12_size) configurations to try.
    """
    m = 2 * n - 1
    configs = []

    if m % 2 == 1:
        # m is always odd since m = 2n-1
        primary_d11 = (m + 1) // 2  # e.g., 22 for m=43
        primary_d12 = (m - 1) // 2  # e.g., 21 for m=43
        flipped_d11 = (m - 1) // 2 - 1  # e.g., 20 for m=43: |D11| symmetric so must be even
        flipped_d12 = (m - 1) // 2  # same |D12|

        # Primary: d1 = d11 + d12 = m, d2 = (m-1-d11) + d12 = m - 2
        configs.append(("primary", primary_d11, primary_d12))

        # Flipped: d1 = d11 + d12 = m - 2, d2 = (m-1-d11) + d12 = m
        configs.append(("flipped", flipped_d11, flipped_d12))

    return configs


def run_solver(n: int, config: str = "both", trials: int = 16,
               max_iter: int = 2000000, T_init: float = 5.0,
               cooling: float = 0.999997, seed: int = 42,
               reheat_interval: int = 0, reheat_factor: float = 0.5):
    """
    Main solver entry point.

    Args:
        n: Book graph parameter (solving R(B_{n-1}, B_n) = 4n-1)
        config: "primary", "flipped", or "both"
        trials: Number of SA trials per configuration
        max_iter: Maximum iterations per trial
        T_init: Initial temperature
        cooling: Cooling rate
        seed: Random seed
        reheat_interval: If > 0, reheat every this many iterations when stuck
        reheat_factor: Fraction of T_init for reheat
    """
    m = 2 * n - 1
    N = 4 * n - 2

    random.seed(seed)

    print("=" * 70)
    print(f"SA Solver for R(B_{{{n-1}}}, B_{{{n}}}) = {4*n - 1}")
    print(f"  m={m}, N={N}, red_thresh={n-2}, blue_thresh={n-1}")
    print(f"  trials={trials}, max_iter={max_iter}, T_init={T_init}, cooling={cooling}")
    if reheat_interval > 0:
        print(f"  reheat_interval={reheat_interval}, reheat_factor={reheat_factor}")
    print("=" * 70)

    all_configs = compute_cardinalities(n)
    if config == "primary":
        all_configs = [c for c in all_configs if c[0] == "primary"]
    elif config == "flipped":
        all_configs = [c for c in all_configs if c[0] == "flipped"]

    for label, d11_size, d12_size in all_configs:
        d1 = d11_size + d12_size
        d22_size = (m - 1) - d11_size
        d2 = d22_size + d12_size

        print(f"\n--- Config {label}: |D11|={d11_size}, |D12|={d12_size} (d1={d1}, d2={d2}) ---")
        start = time.time()

        result = sa_search(n, d11_size, d12_size,
                           max_iter=max_iter, n_trials=trials,
                           T_init=T_init, cooling=cooling,
                           reheat_interval=reheat_interval,
                           reheat_factor=reheat_factor)

        if isinstance(result, tuple) and len(result) == 4:
            D11, D12, D22, viols = result
            print(f"\n*** SOLUTION FOUND! ***")

            G = BlockCirculantGraph(n=n, D11=D11, D12=D12, D22=D22)
            v = verify_construction(G)
            print(f"Valid: {v.valid}")
            print(f"Max red common: {v.max_red_common} (threshold: {v.red_threshold})")
            print(f"Max blue common: {v.max_blue_common} (threshold: {v.blue_threshold})")
            print(f"D11={sorted(D11)}")
            print(f"D12={sorted(D12)}")
            print(f"D22={sorted(D22)}")

            # Save solution
            outfile = os.path.join(
                os.path.dirname(os.path.abspath(__file__)),
                f"solution_n{n}.json"
            )
            save_construction(G, v, outfile, solver=f"SA-informed-{label}")
            print(f"Solution saved to {outfile}")

            elapsed = time.time() - start
            print(f"Elapsed: {elapsed:.1f}s")
            return True, G, v
        else:
            cfg, cost, viols = result
            print(f"\nBest {label}: cost={cost}")
            if cfg:
                D11, D12, D22 = cfg
                print(f"  |D11|={len(D11)}, |D12|={len(D12)}, |D22|={len(D22)}")
                print(f"  d1={len(D11)+len(D12)}, d2={len(D22)+len(D12)}")
            if viols:
                print(f"  Violations ({len(viols)}):")
                for vt, d, ex in viols[:10]:
                    print(f"    {vt} d={d} excess={ex}")

            elapsed = time.time() - start
            print(f"Elapsed: {elapsed:.1f}s")

    print("\nNo solution found in any configuration.")
    return False, None, None


def main():
    parser = argparse.ArgumentParser(
        description="SA solver for R(B_{n-1}, B_n) = 4n-1"
    )
    parser.add_argument("--n", type=int, required=True,
                        help="Book graph parameter n")
    parser.add_argument("--config", choices=["primary", "flipped", "both"],
                        default="both",
                        help="Which cardinality config to try (default: both)")
    parser.add_argument("--trials", type=int, default=16,
                        help="Number of SA trials per config (default: 16)")
    parser.add_argument("--max-iter", type=int, default=2000000,
                        help="Max iterations per trial (default: 2000000)")
    parser.add_argument("--T-init", type=float, default=5.0,
                        help="Initial temperature (default: 5.0)")
    parser.add_argument("--cooling", type=float, default=0.999997,
                        help="Cooling rate (default: 0.999997)")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed (default: 42)")
    args = parser.parse_args()

    run_solver(
        n=args.n,
        config=args.config,
        trials=args.trials,
        max_iter=args.max_iter,
        T_init=args.T_init,
        cooling=args.cooling,
        seed=args.seed,
    )


if __name__ == "__main__":
    main()
