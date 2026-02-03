"""
Core library for Ramsey number search on book graphs.

Book graph B_n: K_2 + K̄_n (n triangles sharing a common edge)
Goal: Find constructions proving R(B_{n-1}, B_n) ≥ 4n - 1

Uses 2-block circulant graphs parameterized by difference sets.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Set, Dict, Tuple, List, Optional
import json


def make_symmetric(S: Set[int], m: int) -> Set[int]:
    """Ensure set S is symmetric: x ∈ S ⟺ -x (mod m) ∈ S."""
    result = set(S)
    for x in S:
        result.add((-x) % m)
    return result


def is_symmetric(S: Set[int], m: int) -> bool:
    """Check if set S is symmetric under negation mod m."""
    return all((-x) % m in S for x in S)


@dataclass
class BlockCirculantGraph:
    """
    2-block circulant graph on 2m = 4n-2 vertices.

    Vertices: V = V1 ∪ V2 where V1 = {0,...,m-1}, V2 = {m,...,2m-1}

    Adjacency defined by difference sets:
    - D11: edges within V1 (must be symmetric)
    - D12: edges from V1 to V2 (not necessarily symmetric)
    - D22: edges within V2 (must be symmetric)

    Edge (u,v) exists iff:
    - u,v ∈ V1: (v-u) mod m ∈ D11
    - u,v ∈ V2: (v-u) mod m ∈ D22
    - u ∈ V1, v ∈ V2: (v-u) mod m ∈ D12
    """
    n: int                      # Book parameter
    D11: Set[int] = field(default_factory=set)
    D12: Set[int] = field(default_factory=set)
    D22: Set[int] = field(default_factory=set)

    def __post_init__(self):
        self.m = 2 * self.n - 1       # Block size
        self.N = 4 * self.n - 2       # Total vertices

        # Ensure D11 and D22 are symmetric
        self.D11 = make_symmetric(self.D11, self.m)
        self.D22 = make_symmetric(self.D22, self.m)

        # Remove 0 from D11 and D22 (no self-loops)
        self.D11.discard(0)
        self.D22.discard(0)

        # Precompute D12 transpose: {-x mod m : x ∈ D12}
        self._D12T = {(-x) % self.m for x in self.D12}

    def adjacent(self, u: int, v: int) -> bool:
        """Check if vertices u and v are adjacent."""
        if u == v:
            return False

        m = self.m
        if u > v:
            u, v = v, u

        if u < m and v < m:
            # Both in V1
            return (v - u) % m in self.D11
        elif u >= m and v >= m:
            # Both in V2
            return (v - u) % m in self.D22
        else:
            # u in V1, v in V2
            return (v - u) % m in self.D12

    def to_adjacency_matrix(self) -> List[List[int]]:
        """Convert to full adjacency matrix (for verification)."""
        adj = [[0] * self.N for _ in range(self.N)]
        for u in range(self.N):
            for v in range(u + 1, self.N):
                if self.adjacent(u, v):
                    adj[u][v] = adj[v][u] = 1
        return adj


def Delta(A: Set[int], B: Set[int], d: int, m: int) -> int:
    """
    Count pairs (a,b) ∈ A×B where a - b ≡ d (mod m).
    Equivalently, count a ∈ A such that (a - d) mod m ∈ B.
    """
    return sum(1 for a in A if (a - d) % m in B)


def Sigma(A: Set[int], B: Set[int], d: int, m: int) -> int:
    """
    Count pairs (a,b) ∈ A×B where a + b ≡ d (mod m).
    Equivalently, count a ∈ A such that (d - a) mod m ∈ B.
    """
    return sum(1 for a in A if (d - a) % m in B)


def common_neighbors_V1V1(G: BlockCirculantGraph, d: int) -> int:
    """
    Count common neighbors for edge within V1 with difference d.

    For edge (u, v) with v - u ≡ d (mod m), common neighbors are:
    - w ∈ V1: (w-u) ∈ D11 and (w-v) ∈ D11, i.e., pairs where a-b = d
    - w ∈ V2: (w-u) ∈ D12 and (w-v) ∈ D12, i.e., pairs where a-b = d
    """
    return Delta(G.D11, G.D11, d, G.m) + Delta(G.D12, G.D12, d, G.m)


def common_neighbors_V2V2(G: BlockCirculantGraph, d: int) -> int:
    """
    Count common neighbors for edge within V2 with difference d.

    For edge (u, v) in V2 with v - u ≡ d (mod m), common neighbors are:
    - w ∈ V2: pairs in D22×D22 with difference d
    - w ∈ V1: need (u-w) ∈ D12 and (v-w) ∈ D12, i.e., pairs in D12^T × D12^T
    """
    return Delta(G.D22, G.D22, d, G.m) + Delta(G._D12T, G._D12T, d, G.m)


def common_neighbors_V1V2(G: BlockCirculantGraph, d: int) -> int:
    """
    Count common neighbors for cross-block edge.

    For edge (u, v) with u ∈ V1, v ∈ V2, and v - u ≡ d (mod m):
    - w ∈ V1: need (w-u) ∈ D11 and (v-w) ∈ D12
      This means (w-u) + (v-w) = v-u = d, so sum of elements = d
    - w ∈ V2: need (u-w) ∈ D12^T (i.e., w-u ∈ D12) and (w-v) ∈ D22
      Let w-u = a ∈ D12, w-v = b ∈ D22, then a-b = d
    """
    return Sigma(G.D11, G.D12, d, G.m) + Delta(G.D12, G.D22, d, G.m)


def common_neighbors_brute_force(G: BlockCirculantGraph, u: int, v: int) -> int:
    """Brute-force count of common neighbors (for verification)."""
    count = 0
    for w in range(G.N):
        if w != u and w != v and G.adjacent(u, w) and G.adjacent(v, w):
            count += 1
    return count


@dataclass
class VerificationResult:
    """Result of verifying a graph construction."""
    valid: bool
    max_red_common: int         # Max common neighbors for red edges
    max_blue_common: int        # Max common neighbors for blue edges
    red_threshold: int          # n - 2 (must have < n-1 = threshold+1)
    blue_threshold: int         # n - 1 (must have < n = threshold+1)
    violations: List[Tuple[str, int, int]]  # (edge_type, d, excess)

    def to_dict(self) -> Dict:
        return {
            "valid": self.valid,
            "max_red_common": self.max_red_common,
            "max_blue_common": self.max_blue_common,
            "red_threshold": self.red_threshold,
            "blue_threshold": self.blue_threshold,
            "num_violations": len(self.violations)
        }


def verify_construction(G: BlockCirculantGraph) -> VerificationResult:
    """
    Verify that the graph avoids red B_{n-1} and blue B_n.

    For a book B_k, we need edge + k common neighbors.
    - Red edges: must have < n-1 common neighbors (avoid B_{n-1})
    - Blue edges: must have < n common neighbors (avoid B_n)

    Blue common neighbors by inclusion-exclusion:
        blue_common = (N - 2) - deg_u - deg_v + red_common
    where deg_u, deg_v are the red degrees of the endpoints.
    """
    n = G.n
    m = G.m
    N = G.N

    red_threshold = n - 2    # Red edge OK if common_red <= n-2
    blue_threshold = n - 1   # Blue edge OK if common_blue <= n-1

    max_red = 0
    max_blue = 0
    violations = []

    # Compute vertex degrees (same for all vertices in each block by circulant structure)
    d1 = len(G.D11) + len(G.D12)  # Degree of vertices in V1
    d2 = len(G.D22) + len(G.D12)  # Degree of vertices in V2

    # Check edges within V1 (differences 1 to m-1)
    for d in range(1, m):
        common = common_neighbors_V1V1(G, d)
        if d in G.D11:
            # Red edge
            max_red = max(max_red, common)
            if common > red_threshold:
                violations.append(("V1V1_red", d, common - red_threshold))
        else:
            # Blue edge: both endpoints in V1 have degree d1
            blue_common = (N - 2) - d1 - d1 + common
            max_blue = max(max_blue, blue_common)
            if blue_common > blue_threshold:
                violations.append(("V1V1_blue", d, blue_common - blue_threshold))

    # Check edges within V2 (differences 1 to m-1)
    for d in range(1, m):
        common = common_neighbors_V2V2(G, d)
        if d in G.D22:
            # Red edge
            max_red = max(max_red, common)
            if common > red_threshold:
                violations.append(("V2V2_red", d, common - red_threshold))
        else:
            # Blue edge: both endpoints in V2 have degree d2
            blue_common = (N - 2) - d2 - d2 + common
            max_blue = max(max_blue, blue_common)
            if blue_common > blue_threshold:
                violations.append(("V2V2_blue", d, blue_common - blue_threshold))

    # Check cross-block edges (differences 0 to m-1)
    for d in range(m):
        common = common_neighbors_V1V2(G, d)
        if d in G.D12:
            # Red edge
            max_red = max(max_red, common)
            if common > red_threshold:
                violations.append(("V1V2_red", d, common - red_threshold))
        else:
            # Blue edge: one endpoint in V1 (degree d1), one in V2 (degree d2)
            blue_common = (N - 2) - d1 - d2 + common
            max_blue = max(max_blue, blue_common)
            if blue_common > blue_threshold:
                violations.append(("V1V2_blue", d, blue_common - blue_threshold))

    return VerificationResult(
        valid=len(violations) == 0,
        max_red_common=max_red,
        max_blue_common=max_blue,
        red_threshold=red_threshold,
        blue_threshold=blue_threshold,
        violations=violations
    )


def compute_violation_cost(G: BlockCirculantGraph) -> int:
    """
    Compute total violation cost (sum of excesses).
    Used as objective function for heuristic search.
    """
    result = verify_construction(G)
    return sum(excess for _, _, excess in result.violations)


def save_construction(G: BlockCirculantGraph, result: VerificationResult,
                      filename: str, solver: str = "unknown"):
    """Save a valid construction to JSON file."""
    data = {
        "n": G.n,
        "target_vertices": G.N,
        "solver": solver,
        "parameters": {
            "m": G.m,
            "D11": sorted(G.D11),
            "D12": sorted(G.D12),
            "D22": sorted(G.D22)
        },
        "verification": {
            "max_red_common": result.max_red_common,
            "max_blue_common": result.max_blue_common,
            "red_threshold": result.red_threshold,
            "blue_threshold": result.blue_threshold,
            "valid": result.valid
        }
    }
    with open(filename, 'w') as f:
        json.dump(data, f, indent=2)


def load_construction(filename: str) -> BlockCirculantGraph:
    """Load a construction from JSON file."""
    with open(filename, 'r') as f:
        data = json.load(f)
    return BlockCirculantGraph(
        n=data["n"],
        D11=set(data["parameters"]["D11"]),
        D12=set(data["parameters"]["D12"]),
        D22=set(data["parameters"]["D22"])
    )


if __name__ == "__main__":
    # Quick sanity check
    G = BlockCirculantGraph(n=6, D11={1, 2}, D12={0, 1, 2}, D22={1, 3})
    print(f"n={G.n}, m={G.m}, N={G.N}")
    print(f"D11={sorted(G.D11)}, D12={sorted(G.D12)}, D22={sorted(G.D22)}")
    result = verify_construction(G)
    print(f"Valid: {result.valid}")
    print(f"Max red common: {result.max_red_common}, threshold: {result.red_threshold}")
    print(f"Max blue common: {result.max_blue_common}, threshold: {result.blue_threshold}")
