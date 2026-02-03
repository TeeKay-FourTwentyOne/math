"""
Heuristic solver for Ramsey book graph constructions.

Implements tabu search and simulated annealing to find valid constructions.
"""

import sys
import time
import random
import math
from typing import Optional, Set, Tuple, List
from dataclasses import dataclass, field
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

from ramsey_core import (
    BlockCirculantGraph, verify_construction, compute_violation_cost,
    save_construction, make_symmetric
)


@dataclass
class SearchState:
    """State for heuristic search."""
    n: int
    m: int
    D11_half: Set[int]  # Representatives for D11 (1 to m//2)
    D12: Set[int]       # Full D12 (0 to m-1)
    D22_half: Set[int]  # Representatives for D22 (1 to m//2)
    cost: int = 0

    def to_graph(self) -> BlockCirculantGraph:
        """Convert state to graph."""
        return BlockCirculantGraph(
            n=self.n,
            D11=self.D11_half.copy(),
            D12=self.D12.copy(),
            D22=self.D22_half.copy()
        )

    def compute_cost(self) -> int:
        """Compute violation cost."""
        G = self.to_graph()
        self.cost = compute_violation_cost(G)
        return self.cost

    def copy(self) -> "SearchState":
        """Create a copy of this state."""
        return SearchState(
            n=self.n,
            m=self.m,
            D11_half=self.D11_half.copy(),
            D12=self.D12.copy(),
            D22_half=self.D22_half.copy(),
            cost=self.cost
        )


def random_state(n: int, density: float = 0.3) -> SearchState:
    """Generate a random initial state."""
    m = 2 * n - 1
    half_m = m // 2

    D11_half = set(random.sample(range(1, half_m + 1),
                                  int(half_m * density)))
    D12 = set(random.sample(range(m), int(m * density)))
    D22_half = set(random.sample(range(1, half_m + 1),
                                  int(half_m * density)))

    state = SearchState(n=n, m=m, D11_half=D11_half, D12=D12, D22_half=D22_half)
    state.compute_cost()
    return state


def get_neighbors(state: SearchState) -> List[Tuple[str, int]]:
    """
    Generate all possible single-flip moves.

    Returns list of (set_name, element) tuples.
    """
    moves = []
    half_m = state.m // 2

    # D11 moves: flip elements in range 1 to m//2
    for d in range(1, half_m + 1):
        moves.append(("D11", d))

    # D12 moves: flip elements in range 0 to m-1
    for d in range(state.m):
        moves.append(("D12", d))

    # D22 moves: flip elements in range 1 to m//2
    for d in range(1, half_m + 1):
        moves.append(("D22", d))

    return moves


def apply_move(state: SearchState, move: Tuple[str, int]) -> SearchState:
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
    elif set_name == "D22":
        if element in new_state.D22_half:
            new_state.D22_half.remove(element)
        else:
            new_state.D22_half.add(element)

    new_state.compute_cost()
    return new_state


def tabu_search(
    n: int,
    max_iter: int = 100000,
    tabu_tenure: int = 50,
    restart_threshold: int = 5000,
    seed: Optional[int] = None,
    verbose: bool = False
) -> Optional[BlockCirculantGraph]:
    """
    Tabu search for valid construction.

    Args:
        n: Book parameter
        max_iter: Maximum iterations
        tabu_tenure: How long moves stay tabu
        restart_threshold: Restart after this many iterations without improvement
        seed: Random seed
        verbose: Print progress

    Returns:
        Valid BlockCirculantGraph or None
    """
    if seed is not None:
        random.seed(seed)

    state = random_state(n)
    best_state = state.copy()
    best_cost = state.cost

    tabu_list = {}  # move -> iteration when it becomes non-tabu
    no_improve_count = 0

    if verbose:
        print(f"Initial cost: {state.cost}")

    for iteration in range(max_iter):
        if state.cost == 0:
            if verbose:
                print(f"Found valid construction at iteration {iteration}")
            return state.to_graph()

        # Get all moves
        moves = get_neighbors(state)
        random.shuffle(moves)

        # Find best non-tabu move (or aspiration)
        best_move = None
        best_move_cost = float('inf')

        for move in moves:
            new_state = apply_move(state, move)

            # Accept if: not tabu OR aspiration (better than global best)
            is_tabu = tabu_list.get(move, 0) > iteration
            aspiration = new_state.cost < best_cost

            if not is_tabu or aspiration:
                if new_state.cost < best_move_cost:
                    best_move = move
                    best_move_cost = new_state.cost
                    best_new_state = new_state

        if best_move is None:
            # All moves are tabu, pick random
            move = random.choice(moves)
            best_new_state = apply_move(state, move)
            best_move = move

        # Apply move
        state = best_new_state
        tabu_list[best_move] = iteration + tabu_tenure

        # Update global best
        if state.cost < best_cost:
            best_cost = state.cost
            best_state = state.copy()
            no_improve_count = 0
            if verbose and iteration % 1000 == 0:
                print(f"Iter {iteration}: new best cost = {best_cost}")
        else:
            no_improve_count += 1

        # Restart if stuck
        if no_improve_count >= restart_threshold:
            if verbose:
                print(f"Iter {iteration}: restarting (no improvement for {restart_threshold} iters)")
            state = random_state(n)
            no_improve_count = 0

    if verbose:
        print(f"Reached max iterations. Best cost: {best_cost}")

    if best_cost == 0:
        return best_state.to_graph()
    return None


def simulated_annealing(
    n: int,
    max_iter: int = 100000,
    initial_temp: float = 10.0,
    cooling_rate: float = 0.9999,
    seed: Optional[int] = None,
    verbose: bool = False
) -> Optional[BlockCirculantGraph]:
    """
    Simulated annealing for valid construction.

    Args:
        n: Book parameter
        max_iter: Maximum iterations
        initial_temp: Initial temperature
        cooling_rate: Temperature decay per iteration
        seed: Random seed
        verbose: Print progress

    Returns:
        Valid BlockCirculantGraph or None
    """
    if seed is not None:
        random.seed(seed)

    state = random_state(n)
    best_state = state.copy()
    best_cost = state.cost
    temp = initial_temp

    if verbose:
        print(f"Initial cost: {state.cost}")

    for iteration in range(max_iter):
        if state.cost == 0:
            if verbose:
                print(f"Found valid construction at iteration {iteration}")
            return state.to_graph()

        # Pick random move
        moves = get_neighbors(state)
        move = random.choice(moves)
        new_state = apply_move(state, move)

        # Accept or reject
        delta = new_state.cost - state.cost
        if delta < 0 or random.random() < math.exp(-delta / temp):
            state = new_state

            if state.cost < best_cost:
                best_cost = state.cost
                best_state = state.copy()
                if verbose and iteration % 10000 == 0:
                    print(f"Iter {iteration}: new best cost = {best_cost}, temp = {temp:.4f}")

        temp *= cooling_rate

    if verbose:
        print(f"Reached max iterations. Best cost: {best_cost}")

    if best_cost == 0:
        return best_state.to_graph()
    return None


def worker_tabu(args):
    """Worker function for parallel tabu search."""
    n, max_iter, seed, verbose = args
    return tabu_search(n, max_iter=max_iter, seed=seed, verbose=verbose)


def worker_sa(args):
    """Worker function for parallel simulated annealing."""
    n, max_iter, seed, verbose = args
    return simulated_annealing(n, max_iter=max_iter, seed=seed, verbose=verbose)


def parallel_search(
    n: int,
    method: str = "tabu",
    num_workers: int = 10,
    max_iter: int = 100000,
    verbose: bool = False
) -> Optional[BlockCirculantGraph]:
    """
    Run parallel search with multiple random seeds.

    Args:
        n: Book parameter
        method: "tabu" or "sa" (simulated annealing)
        num_workers: Number of parallel workers
        max_iter: Max iterations per worker
        verbose: Print progress

    Returns:
        Valid BlockCirculantGraph or None
    """
    worker = worker_tabu if method == "tabu" else worker_sa
    args = [(n, max_iter, seed, verbose) for seed in range(num_workers)]

    print(f"Starting parallel {method} search with {num_workers} workers...")

    # Use process pool for true parallelism
    with ProcessPoolExecutor(max_workers=min(num_workers, multiprocessing.cpu_count())) as executor:
        futures = [executor.submit(worker, arg) for arg in args]

        for future in as_completed(futures):
            result = future.result()
            if result is not None:
                # Cancel remaining futures
                for f in futures:
                    f.cancel()
                return result

    return None


def main():
    """Main entry point."""
    import argparse
    parser = argparse.ArgumentParser(description="Heuristic solver for Ramsey book graphs")
    parser.add_argument("n", type=int, help="Book parameter n")
    parser.add_argument("--method", choices=["tabu", "sa", "parallel"],
                       default="tabu", help="Search method")
    parser.add_argument("--max-iter", type=int, default=100000,
                       help="Maximum iterations")
    parser.add_argument("--workers", type=int, default=10,
                       help="Number of parallel workers")
    parser.add_argument("--seed", type=int, help="Random seed (single-threaded only)")
    parser.add_argument("--verbose", "-v", action="store_true",
                       help="Verbose output")
    parser.add_argument("--output", type=str, help="Output JSON file")
    args = parser.parse_args()

    start_time = time.time()

    if args.method == "parallel":
        result = parallel_search(args.n, num_workers=args.workers,
                                max_iter=args.max_iter, verbose=args.verbose)
    elif args.method == "tabu":
        result = tabu_search(args.n, max_iter=args.max_iter,
                            seed=args.seed, verbose=args.verbose)
    else:
        result = simulated_annealing(args.n, max_iter=args.max_iter,
                                     seed=args.seed, verbose=args.verbose)

    elapsed = time.time() - start_time

    if result:
        verification = verify_construction(result)
        print(f"\nFound construction for n={args.n} in {elapsed:.2f}s!")
        print(f"D11 = {sorted(result.D11)}")
        print(f"D12 = {sorted(result.D12)}")
        print(f"D22 = {sorted(result.D22)}")
        print(f"Valid: {verification.valid}")
        print(f"Max red common: {verification.max_red_common}")
        print(f"Max blue common: {verification.max_blue_common}")

        if args.output:
            save_construction(result, verification, args.output,
                            solver=f"heuristic_{args.method}")
            print(f"Saved to {args.output}")
    else:
        print(f"\nNo construction found for n={args.n} after {elapsed:.2f}s")


if __name__ == "__main__":
    main()
