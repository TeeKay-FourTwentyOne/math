"""
Search for R(B_21, B_22) = 87 construction.

Uses quadratic residues mod 43 as seed and D22 = complement(D11).
"""

import time
import random
import math
from typing import Set, Tuple, Optional
from dataclasses import dataclass

from ramsey_core import (
    BlockCirculantGraph, verify_construction, compute_violation_cost,
    save_construction, make_symmetric
)


def quadratic_residues(m: int) -> Set[int]:
    """Compute quadratic residues mod m (excluding 0)."""
    return {(x * x) % m for x in range(1, m)} - {0}


@dataclass
class SearchStateComplement:
    """Search state with D22 = complement(D11)."""
    n: int
    m: int
    D11_half: Set[int]  # Representatives for D11 (1 to m//2)
    D12: Set[int]       # Full D12 (0 to m-1)
    cost: int = 0

    def get_D22_half(self) -> Set[int]:
        """D22 = complement of D11 in {1, ..., m-1}."""
        # D11 is symmetric, so we work with half representatives
        # D22_half = {1..m//2} - D11_half
        return set(range(1, self.m // 2 + 1)) - self.D11_half

    def to_graph(self) -> BlockCirculantGraph:
        """Convert state to graph."""
        return BlockCirculantGraph(
            n=self.n,
            D11=self.D11_half.copy(),
            D12=self.D12.copy(),
            D22=self.get_D22_half()
        )

    def compute_cost(self) -> int:
        """Compute violation cost."""
        G = self.to_graph()
        self.cost = compute_violation_cost(G)
        return self.cost

    def copy(self) -> "SearchStateComplement":
        """Create a copy of this state."""
        return SearchStateComplement(
            n=self.n,
            m=self.m,
            D11_half=self.D11_half.copy(),
            D12=self.D12.copy(),
            cost=self.cost
        )


def get_neighbors_complement(state: SearchStateComplement):
    """Generate all possible single-flip moves (D11 and D12 only, D22 is derived)."""
    moves = []
    half_m = state.m // 2

    # D11 moves: flip elements in range 1 to m//2
    for d in range(1, half_m + 1):
        moves.append(("D11", d))

    # D12 moves: flip elements in range 0 to m-1
    for d in range(state.m):
        moves.append(("D12", d))

    # No D22 moves - it's derived from D11

    return moves


def apply_move_complement(state: SearchStateComplement, move: Tuple[str, int]) -> SearchStateComplement:
    """Apply a move and return new state."""
    new_state = state.copy()
    set_name, element = move

    if set_name == "D11":
        if element in new_state.D11_half:
            new_state.D11_half.remove(element)
        else:
            new_state.D11_half.add(element)
    elif set_name == "D12":
        if element in new_state.D12:
            new_state.D12.remove(element)
        else:
            new_state.D12.add(element)

    new_state.compute_cost()
    return new_state


def tabu_search_complement(
    n: int,
    initial_D11_half: Set[int],
    initial_D12: Set[int],
    max_iter: int = 500000,
    tabu_tenure: int = 100,
    restart_threshold: int = 20000,
    verbose: bool = True
) -> Tuple[Optional[BlockCirculantGraph], int, float]:
    """
    Tabu search with D22 = complement(D11).

    Returns:
        (best_graph, best_cost, elapsed_time)
    """
    m = 2 * n - 1

    # Convert initial D11 to half representatives
    D11_half = {min(d, m - d) for d in initial_D11_half if d != 0}

    state = SearchStateComplement(
        n=n, m=m,
        D11_half=D11_half,
        D12=initial_D12.copy()
    )
    state.compute_cost()

    best_state = state.copy()
    best_cost = state.cost

    tabu_list = {}
    no_improve_count = 0

    start_time = time.time()

    if verbose:
        G = state.to_graph()
        print(f"Initial: cost={state.cost}, |D11|={len(G.D11)}, |D12|={len(G.D12)}, |D22|={len(G.D22)}")

    for iteration in range(max_iter):
        if state.cost == 0:
            elapsed = time.time() - start_time
            if verbose:
                print(f"SOLUTION FOUND at iteration {iteration}!")
            return state.to_graph(), 0, elapsed

        # Get all moves
        moves = get_neighbors_complement(state)
        random.shuffle(moves)

        # Find best non-tabu move (or aspiration)
        best_move = None
        best_move_cost = float('inf')
        best_new_state = None

        for move in moves:
            new_state = apply_move_complement(state, move)

            is_tabu = tabu_list.get(move, 0) > iteration
            aspiration = new_state.cost < best_cost

            if not is_tabu or aspiration:
                if new_state.cost < best_move_cost:
                    best_move = move
                    best_move_cost = new_state.cost
                    best_new_state = new_state

        if best_move is None:
            move = random.choice(moves)
            best_new_state = apply_move_complement(state, move)
            best_move = move

        state = best_new_state
        tabu_list[best_move] = iteration + tabu_tenure

        if state.cost < best_cost:
            best_cost = state.cost
            best_state = state.copy()
            no_improve_count = 0
            if verbose:
                elapsed = time.time() - start_time
                print(f"Iter {iteration}: new best cost = {best_cost} (elapsed: {elapsed:.1f}s)")
        else:
            no_improve_count += 1

        if no_improve_count >= restart_threshold:
            if verbose:
                print(f"Iter {iteration}: restarting from best (no improvement for {restart_threshold} iters)")
            # Restart from best with small perturbation
            state = best_state.copy()
            # Apply random perturbation
            moves = get_neighbors_complement(state)
            for _ in range(3):
                move = random.choice(moves)
                state = apply_move_complement(state, move)
            no_improve_count = 0

    elapsed = time.time() - start_time
    return best_state.to_graph(), best_cost, elapsed


def simulated_annealing_complement(
    n: int,
    initial_D11_half: Set[int],
    initial_D12: Set[int],
    max_iter: int = 1000000,
    initial_temp: float = 50.0,
    cooling_rate: float = 0.99999,
    verbose: bool = True
) -> Tuple[Optional[BlockCirculantGraph], int, float]:
    """
    Simulated annealing with D22 = complement(D11).
    """
    m = 2 * n - 1

    D11_half = {min(d, m - d) for d in initial_D11_half if d != 0}

    state = SearchStateComplement(
        n=n, m=m,
        D11_half=D11_half,
        D12=initial_D12.copy()
    )
    state.compute_cost()

    best_state = state.copy()
    best_cost = state.cost
    temp = initial_temp

    start_time = time.time()

    if verbose:
        G = state.to_graph()
        print(f"Initial: cost={state.cost}, |D11|={len(G.D11)}, |D12|={len(G.D12)}, |D22|={len(G.D22)}")

    for iteration in range(max_iter):
        if state.cost == 0:
            elapsed = time.time() - start_time
            if verbose:
                print(f"SOLUTION FOUND at iteration {iteration}!")
            return state.to_graph(), 0, elapsed

        moves = get_neighbors_complement(state)
        move = random.choice(moves)
        new_state = apply_move_complement(state, move)

        delta = new_state.cost - state.cost
        if delta < 0 or random.random() < math.exp(-delta / temp):
            state = new_state

            if state.cost < best_cost:
                best_cost = state.cost
                best_state = state.copy()
                if verbose:
                    elapsed = time.time() - start_time
                    print(f"Iter {iteration}: new best cost = {best_cost}, temp = {temp:.4f} (elapsed: {elapsed:.1f}s)")

        temp *= cooling_rate

        if iteration % 100000 == 0 and iteration > 0 and verbose:
            elapsed = time.time() - start_time
            print(f"Iter {iteration}: current cost = {state.cost}, best = {best_cost}, temp = {temp:.4f} (elapsed: {elapsed:.1f}s)")

    elapsed = time.time() - start_time
    return best_state.to_graph(), best_cost, elapsed


def get_violations(G: BlockCirculantGraph):
    """Get list of constraint violations."""
    result = verify_construction(G)
    return result.violations


def main():
    n = 22
    m = 43
    N = 86

    print(f"Searching for R(B_21, B_22) construction")
    print(f"n={n}, m={m}, N={N}")
    print("=" * 60)

    # Quadratic residues mod 43
    QR43 = quadratic_residues(43)
    print(f"QR_43 = {sorted(QR43)}")
    print(f"|QR_43| = {len(QR43)}")

    # Non-residues
    NQR43 = set(range(1, 43)) - QR43
    print(f"NQR_43 = {sorted(NQR43)}")
    print(f"|NQR_43| = {len(NQR43)}")

    # NOTE: For m=43 (3 mod 4), QR43 is NOT symmetric under negation.
    # The half-representatives of QR43 cover all of {1,...,21}.
    # We need to construct symmetric D11 sets differently.

    # Strategy: Use QR43 directly for D11 (it has 21 elements).
    # When we expand to symmetric, we'd get 42 elements (too many).
    # Instead, we use QR43 as D11 WITHOUT symmetry enforcement,
    # OR we pick a subset of half-representatives.

    # Approach: Since the constraint D22 = complement(D11), and
    # we want |D11| ~ 21, let's pick ~10-11 half-representatives
    # so |D11| ~ 20-22 after symmetry.

    # Try different D11 seeds based on QR43:
    # - Take every other element of QR43's half-reps
    # - Random subsets of size ~10-11
    half_m = m // 2  # 21

    # Seed 1: First 11 elements of {1..21}
    D11_seed1 = set(range(1, 12))  # {1,...,11} -> |D11| = 22

    # Seed 2: Every other element
    D11_seed2 = set(range(1, 22, 2))  # {1,3,5,...,21} -> |D11| = 22

    # Seed 3: QR43 intersected with {1..21} (already all of them, so pick half)
    D11_seed3 = {1, 4, 6, 9, 10, 11, 13, 14, 15, 16, 17}  # First 11 from QR43

    # Initial configurations to try (ordered by initial cost, best first)
    configs = [
        ("D11=QR43[:11], D12=QR43", D11_seed3, QR43.copy()),  # Cost 56
        ("D11=QR43[:11], D12=QR43+{0}", D11_seed3, QR43 | {0}),  # Cost 86
        ("D11=odd, D12=QR43", D11_seed2, QR43.copy()),  # Cost 148
    ]

    best_overall_cost = float('inf')
    best_overall_graph = None
    best_config_name = None

    total_start = time.time()

    for config_name, D11_init, D12_init in configs:
        print()
        print("=" * 60)
        print(f"Trying configuration: D11={config_name.split(',')[0]}, D12={config_name.split(',')[-1].strip()}")
        print(f"Initial |D11|={len(D11_init)}, |D12|={len(D12_init)}")

        # Run tabu search first
        print("\n--- Tabu Search ---")
        G, cost, elapsed = tabu_search_complement(
            n=n,
            initial_D11_half=D11_init,
            initial_D12=D12_init,
            max_iter=200000,
            tabu_tenure=100,
            restart_threshold=15000,
            verbose=True
        )

        print(f"Tabu result: cost={cost}, time={elapsed:.1f}s")

        if cost < best_overall_cost:
            best_overall_cost = cost
            best_overall_graph = G
            best_config_name = config_name

        if cost == 0:
            print("SOLUTION FOUND!")
            break

        # If tabu didn't find solution, try SA from best tabu result
        if cost > 0:
            print("\n--- Simulated Annealing (from tabu best) ---")
            G_sa, cost_sa, elapsed_sa = simulated_annealing_complement(
                n=n,
                initial_D11_half=G.D11,
                initial_D12=G.D12,
                max_iter=300000,
                initial_temp=20.0,
                cooling_rate=0.99998,
                verbose=True
            )

            print(f"SA result: cost={cost_sa}, time={elapsed_sa:.1f}s")

            if cost_sa < best_overall_cost:
                best_overall_cost = cost_sa
                best_overall_graph = G_sa
                best_config_name = config_name + " (SA)"

            if cost_sa == 0:
                print("SOLUTION FOUND!")
                break

    total_elapsed = time.time() - total_start

    print()
    print("=" * 60)
    print("FINAL RESULTS")
    print("=" * 60)
    print(f"Best violation count: {best_overall_cost}")
    print(f"Best configuration: {best_config_name}")
    print(f"Total time elapsed: {total_elapsed:.1f}s")

    if best_overall_graph:
        print(f"\nBest (D11, D12) configuration:")
        print(f"D11 = {sorted(best_overall_graph.D11)}")
        print(f"D12 = {sorted(best_overall_graph.D12)}")
        print(f"D22 = {sorted(best_overall_graph.D22)}")

        if best_overall_cost > 0:
            violations = get_violations(best_overall_graph)
            print(f"\nViolations ({len(violations)}):")
            for v in violations[:20]:  # Show first 20
                print(f"  {v}")
            if len(violations) > 20:
                print(f"  ... and {len(violations) - 20} more")

        if best_overall_cost == 0:
            print("\nVerifying solution...")
            result = verify_construction(best_overall_graph)
            print(f"Valid: {result.valid}")
            print(f"Max red common: {result.max_red_common} (threshold: {result.red_threshold})")
            print(f"Max blue common: {result.max_blue_common} (threshold: {result.blue_threshold})")

            save_construction(best_overall_graph, result, "solution_n22.json", solver="tabu+SA")
            print("Saved to solution_n22.json")


if __name__ == "__main__":
    random.seed(42)
    main()
