"""
Cayley graph search over the Dihedral group D_43 for R(B_21, B_22) = 87.

The dihedral group D_m has order 2m = 86, matching our vertex count N = 86.
Elements: (r, s) where r in Z_m (rotation) and s in {0, 1} (reflection).
Group operation: (r1, s1) * (r2, s2) = (r1 + (-1)^s1 * r2 mod m, s1 XOR s2).

A Cayley graph Cay(D_m, S) has:
- Vertices: all 2m group elements
- Edges: (g, h) iff g^{-1} h in S (or equivalently h g^{-1} in S for right Cayley)

Connection set S must be:
- Inverse-closed: s in S => s^{-1} in S
- Not contain identity

This gives a vertex-transitive graph, so all vertices have the same degree |S|.

Key difference from 2-block circulant: the non-abelian structure means the
adjacency pattern is NOT determined by simple difference sets.
"""

import sys, os, time, random, math
from typing import Set, Tuple, List, Dict, Optional
from itertools import combinations

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

M = 43
N = 2 * M  # 86
N_PARAM = 22
RED_THRESH = 20  # n - 2
BLUE_THRESH = 21  # n - 1


# === Dihedral Group D_43 ===

def d_inv(r: int, s: int) -> Tuple[int, int]:
    """Inverse of (r, s) in D_m.
    (r, 0)^{-1} = (-r, 0)
    (r, 1)^{-1} = (r, 1)  [reflections are self-inverse]
    """
    if s == 0:
        return ((-r) % M, 0)
    else:
        return (r, 1)


def d_mul(r1: int, s1: int, r2: int, s2: int) -> Tuple[int, int]:
    """Multiply (r1, s1) * (r2, s2) in D_m."""
    if s1 == 0:
        return ((r1 + r2) % M, s2)
    else:
        return ((r1 - r2) % M, 1 - s2)


def d_diff(r1: int, s1: int, r2: int, s2: int) -> Tuple[int, int]:
    """Compute g^{-1} h where g = (r1,s1), h = (r2,s2)."""
    gi = d_inv(r1, s1)
    return d_mul(gi[0], gi[1], r2, s2)


# Enumerate all non-identity elements
ALL_ELEMENTS = []
for r in range(M):
    for s in range(2):
        if r == 0 and s == 0:
            continue
        ALL_ELEMENTS.append((r, s))

# Index mapping: element -> index (0..84)
ELEM_TO_IDX = {}
IDX_TO_ELEM = {}
for i, e in enumerate(ALL_ELEMENTS):
    ELEM_TO_IDX[e] = i
    IDX_TO_ELEM[i] = e

# Vertex indexing: (r, s) -> vertex number
def vertex_idx(r: int, s: int) -> int:
    return r + s * M

def vertex_elem(v: int) -> Tuple[int, int]:
    return (v % M, v // M)


# === Inverse-closed connection sets ===

def find_inverse_pairs():
    """Find all {g, g^{-1}} pairs among non-identity elements.
    Self-inverse elements form singleton 'pairs'."""
    pairs = []
    seen = set()

    for e in ALL_ELEMENTS:
        if e in seen:
            continue
        inv_e = d_inv(e[0], e[1])
        if inv_e == e:
            # Self-inverse
            pairs.append((e,))
            seen.add(e)
        else:
            pairs.append((e, inv_e))
            seen.add(e)
            seen.add(inv_e)

    return pairs


INVERSE_PAIRS = find_inverse_pairs()
# Each pair contributes 1 (self-inverse) or 2 elements to S


def connection_set_from_pair_selection(selected_pairs: Set[int]) -> Set[Tuple[int, int]]:
    """Build connection set S from selected pair indices."""
    S = set()
    for i in selected_pairs:
        for e in INVERSE_PAIRS[i]:
            S.add(e)
    return S


# === Common neighbor counting ===

def count_common_neighbors(S: Set[Tuple[int, int]], g_diff: Tuple[int, int]) -> int:
    """
    Count common neighbors of vertices u and v where g^{-1}h = g_diff.

    In Cay(D_m, S), common neighbors of u,v are:
    w such that u^{-1}w in S and v^{-1}w in S.

    Let g = u^{-1}v = g_diff. Then:
    u^{-1}w = s1 in S, and v^{-1}w = (u^{-1}v)^{-1} (u^{-1}w) = g^{-1} s1.
    So we need s1 in S and g^{-1} s1 in S.

    lambda_red(g_diff) = |{s in S : g_diff^{-1} s in S}|
                       = |S intersect (g_diff * S)|
    """
    g_inv = d_inv(g_diff[0], g_diff[1])
    count = 0
    for s in S:
        # Check if g_inv * s is in S
        prod = d_mul(g_inv[0], g_inv[1], s[0], s[1])
        if prod in S:
            count += 1
    return count


def evaluate_connection_set(S: Set[Tuple[int, int]]) -> Tuple[int, int, List]:
    """
    Evaluate a connection set S.

    Returns (cost, n_violations, violation_list).
    """
    degree = len(S)
    cost = 0
    viols = []

    # For each non-identity element g, check the constraint
    for g in ALL_ELEMENTS:
        lam_red = count_common_neighbors(S, g)

        if g in S:
            # Red edge
            if lam_red > RED_THRESH:
                excess = lam_red - RED_THRESH
                cost += excess
                viols.append(("red", g, excess))
        else:
            # Blue edge: both endpoints have degree = |S|
            lam_blue = (N - 2) - 2 * degree + lam_red
            if lam_blue > BLUE_THRESH:
                excess = lam_blue - BLUE_THRESH
                cost += excess
                viols.append(("blue", g, excess))

    return cost, len(viols), viols


# === SA Search ===

def sa_search_dihedral(target_degree: int, max_iter: int = 1000000,
                       n_trials: int = 8, T_init: float = 5.0,
                       cooling: float = 0.999995):
    """SA search over inverse-closed connection sets in D_43."""

    # Classify pairs by size
    pair_sizes = [len(p) for p in INVERSE_PAIRS]
    n_pairs = len(INVERSE_PAIRS)

    print(f"Dihedral D_{M}: {len(ALL_ELEMENTS)} non-identity elements, "
          f"{n_pairs} inverse pairs")
    print(f"  Self-inverse (reflections): {sum(1 for s in pair_sizes if s == 1)}")
    print(f"  Rotation pairs: {sum(1 for s in pair_sizes if s == 2)}")
    print(f"Target degree: {target_degree}")

    best_global_cost = float('inf')
    best_global_S = None
    best_global_viols = None

    for trial in range(n_trials):
        # Random initial selection achieving target degree
        indices = list(range(n_pairs))
        random.shuffle(indices)
        selected = set()
        current_size = 0
        for i in indices:
            if current_size + pair_sizes[i] <= target_degree:
                selected.add(i)
                current_size += pair_sizes[i]
            if current_size == target_degree:
                break

        if current_size != target_degree:
            # Try to adjust
            pass

        S = connection_set_from_pair_selection(selected)
        cost, nviols, viols = evaluate_connection_set(S)
        best_cost = cost
        best_selected = selected.copy()
        best_viols_local = viols

        temp = T_init

        for it in range(max_iter):
            if cost == 0:
                print(f"  SOLUTION FOUND at trial={trial}, iter={it}")
                return S, viols

            # Swap move: remove one pair, add another
            in_list = list(selected)
            out_list = list(set(range(n_pairs)) - selected)
            if not in_list or not out_list:
                continue

            # Try to maintain degree
            pi_out = random.choice(in_list)
            # Find a replacement with same size
            size_needed = pair_sizes[pi_out]
            candidates = [j for j in out_list if pair_sizes[j] == size_needed]
            if not candidates:
                # Allow size change
                candidates = out_list
            pi_in = random.choice(candidates)

            selected.discard(pi_out)
            selected.add(pi_in)
            new_S = connection_set_from_pair_selection(selected)
            new_cost, new_nviols, new_viols = evaluate_connection_set(new_S)

            delta = new_cost - cost
            if delta < 0 or random.random() < math.exp(-delta / max(temp, 0.001)):
                S = new_S
                cost = new_cost
                viols = new_viols
                if cost < best_cost:
                    best_cost = cost
                    best_selected = selected.copy()
                    best_viols_local = viols
            else:
                # Undo
                selected.discard(pi_in)
                selected.add(pi_out)

            temp *= cooling

            if it % 100000 == 0:
                print(f"  t{trial} i{it}: cost={cost} best={best_cost} T={temp:.4f} "
                      f"|S|={len(S)}")
                sys.stdout.flush()

        print(f"  Trial {trial}: best={best_cost} ({len(best_viols_local)} viols, "
              f"|S|={sum(pair_sizes[i] for i in best_selected)})")
        sys.stdout.flush()

        if best_cost < best_global_cost:
            best_global_cost = best_cost
            best_global_S = connection_set_from_pair_selection(best_selected)
            best_global_viols = best_viols_local

    return best_global_S, best_global_cost, best_global_viols


def main():
    random.seed(42)

    print("=" * 70)
    print("Cayley Graph Search on Dihedral Group D_43")
    print(f"Looking for Cay(D_43, S) avoiding red B_21 and blue B_22")
    print("=" * 70)

    # Inverse pairs analysis
    pair_sizes = [len(p) for p in INVERSE_PAIRS]
    print(f"\n{len(INVERSE_PAIRS)} inverse pairs:")
    print(f"  {sum(1 for s in pair_sizes if s == 1)} self-inverse (reflections)")
    print(f"  {sum(1 for s in pair_sizes if s == 2)} rotation pairs")

    # Try degrees around 42-43 (matching the 2-block circulant sweet spot)
    for target_deg in [42, 43, 41, 44]:
        print(f"\n{'='*60}")
        print(f"Target degree: {target_deg}")
        print("=" * 60)

        result = sa_search_dihedral(
            target_degree=target_deg,
            max_iter=500000,
            n_trials=4,
            T_init=5.0,
            cooling=0.999995
        )

        if isinstance(result, tuple) and len(result) == 2:
            S, viols = result
            print(f"\nSOLUTION FOUND! |S|={len(S)}")
        else:
            S, cost, viols = result
            print(f"\nBest for deg={target_deg}: cost={cost}")
            if viols:
                print(f"Violations ({len(viols)}):")
                for vt, g, ex in viols[:10]:
                    print(f"  {vt} g={g} excess={ex}")

    print(f"\nDone.")


if __name__ == "__main__":
    main()
