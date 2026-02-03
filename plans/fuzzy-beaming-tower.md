# Implementation Plan: Ramsey Number Search Pipeline

## Goal
Build a computational pipeline to find lower bound constructions for R(B_{n-1}, B_n) ≥ 4n-1.

**Target**: n=18, 22, 27 (composite values where 2n-1 is not a prime power ≡ 1 mod 4)

## Approach: 2-Block Circulant Graphs

Parameterize graphs by three difference sets (D11, D12, D22) ⊂ Z_m instead of full adjacency matrices.

---

## Phase 1: Core Library (`ramsey_core.py`)

### 1.1 Data Structure

```python
class BlockCirculantGraph:
    def __init__(self, n: int, D11: set, D12: set, D22: set):
        self.n = n           # Book parameter
        self.m = 2*n - 1     # Block size
        self.N = 4*n - 2     # Total vertices = 2m
        self.D11 = D11       # Edges within V1 (symmetric)
        self.D12 = D12       # Edges between V1 and V2
        self.D22 = D22       # Edges within V2 (symmetric)
```

**Vertices**: V = V1 ∪ V2 where V1 = {0, ..., m-1}, V2 = {m, ..., 2m-1}

**Adjacency**:
- u, v ∈ V1: adjacent iff (v - u) mod m ∈ D11
- u, v ∈ V2: adjacent iff (v - u) mod m ∈ D22
- u ∈ V1, v ∈ V2: adjacent iff (v - u) mod m ∈ D12

**Symmetry constraints**:
- D11 and D22 must be symmetric: x ∈ D ⟺ -x ∈ D (mod m)
- D12 is NOT necessarily symmetric

### 1.2 Fast Common Neighbor Counting

**Helper functions** (O(|A| · |B|) but can precompute):

```python
def Delta(A, B, d, m):
    """Count pairs (a,b) ∈ A×B where a - b ≡ d (mod m)"""
    return sum(1 for a in A if (a - d) % m in B)

def Sigma(A, B, d, m):
    """Count pairs (a,b) ∈ A×B where a + b ≡ d (mod m)"""
    return sum(1 for a in A if (d - a) % m in B)
```

**Common neighbor count λ(u, v)** for edge (u, v):

1. **Both in V1** (let d = v - u mod m):
   ```
   λ = Δ(D11, D11, d) + Δ(D12, D12, d)
   ```

2. **Both in V2** (let d = v - u mod m):
   ```
   λ = Δ(D22, D22, d) + Δ(D12^T, D12^T, d)
   ```
   where D12^T = {(-x) mod m : x ∈ D12}

3. **Cross-block** u ∈ V1, v ∈ V2 (let d = v - u mod m):
   ```
   λ = Σ(D11, D12, d) + Δ(D12, D22, d)
   ```

### 1.3 Verification Logic

For each edge type and difference value d:
- Compute λ_red = common neighbors in red (edge present)
- Compute λ_blue = (N - 2) - λ_red (common neighbors in complement)

**Constraints**:
- RED edge (d ∈ D): require λ_red < n - 1 (no B_{n-1})
- BLUE edge (d ∉ D): require λ_blue < n (no B_n)

---

## Phase 2: SAT Solver (`solver_sat.py`)

### 2.1 Variables

- `x[i]` for i ∈ {1, ..., ⌊m/2⌋}: represents i ∈ D11 (symmetry enforced)
- `y[i]` for i ∈ {0, ..., m-1}: represents i ∈ D12
- `z[i]` for i ∈ {1, ..., ⌊m/2⌋}: represents i ∈ D22

### 2.2 Constraints

For each representative edge (d = 1, ..., m-1 for intra-block, d = 0, ..., m-1 for cross-block):
- Express λ as sum of indicator variables
- Use `pysat.card.CardEnc.atmost(literals, bound)` for cardinality constraints

### 2.3 Solver Config
- Use Glucose4 backend
- Timeout: 10 minutes per instance

---

## Phase 3: Heuristic Solver (`solver_heuristic.py`)

### 3.1 Tabu Search

```python
def tabu_search(n, max_iter=100000, tabu_tenure=50):
    # Initialize random D11, D12, D22 (respecting symmetry)
    # Cost = Σ max(0, λ_actual - λ_allowed) over all edges
    # Move: flip single element in one set
    # Tabu: don't reverse recent flips for T steps
    # Restart: perturb 10% if stuck
```

### 3.2 Parallel Execution
- Run 10 seeds in parallel
- Stop when any seed finds valid construction

---

## Phase 4: Output

### JSON Format
```json
{
  "n": 18,
  "target_vertices": 70,
  "parameters": {
    "m": 35,
    "D11": [1, 3, 4, ...],
    "D12": [0, 2, 5, ...],
    "D22": [2, 6, 7, ...]
  },
  "verification": {
    "max_red_book_size": 17,
    "max_blue_book_size": 18,
    "valid": true
  }
}
```

---

## Implementation Order

1. `ramsey_core.py` - BlockCirculantGraph class with Δ, Σ helpers and verify()
2. `test_core.py` - Unit test with n=6 against brute-force matrix check
3. `solver_sat.py` - SAT solver using PySAT
4. `solver_heuristic.py` - Tabu search / simulated annealing
5. Run for n=18

---

## Files to Create

```
ramsey-book-graphs/
├── ramsey_core.py       # Core library
├── solver_sat.py        # SAT-based solver
├── solver_heuristic.py  # Tabu/SA solver
├── test_core.py         # Unit tests
└── results/             # Output JSON files
```

---

## Verification Strategy

1. **Unit test (n=6)**: Compare fast formula results against brute-force adjacency matrix
2. **Sanity check**: Random D11/D12/D22 should almost always fail verification
3. **End-to-end**: Run SAT for small n, confirm valid construction
4. **Production**: Run for n=18 with SAT, fall back to heuristic if timeout
