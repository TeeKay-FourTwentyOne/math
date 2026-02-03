# Implementation Instructions: Computational Search for R(B₂₁, B₂₂) = 87

## Executive Summary

**Objective**: Find a 2-block circulant graph on N = 86 vertices that avoids:
- Red B₂₁ (book with 21 triangles sharing a common edge)
- Blue B₂₂ (book with 22 triangles sharing a common edge)

**Success Criterion**: Discovery of connection sets (D₁₁, D₁₂) that satisfy all book-avoidance constraints.

**Why n=22**: This is the **first genuinely open case**:
- n ≤ 21: Solved by Wesley (2024) and Lidický et al. (2024)
- n=27: Trivial (2n-1=53 ≡ 1 mod 4, Paley construction guaranteed)
- n=22: 2n-1=43 is prime but ≡ 3 mod 4, no theoretical guarantee exists

---

## Phase 0: Data Acquisition and Validation

### 0.1 Clone Repositories

```bash
# Create working directory
mkdir -p ~/ramsey_books && cd ~/ramsey_books

# Clone known solution repositories
git clone https://github.com/gwen-mckinley/ramsey-books-wheels.git
git clone https://github.com/Steven-VO/circulant-Ramsey.git
```

### 0.2 Parse and Validate Known Solutions

Create `parse_known_solutions.py`:

```python
"""
Parse all known 2-block circulant constructions for n ≤ 21.
Verify the D₂₂ = complement(D₁₁) hypothesis.
Extract structural patterns.
"""

import os
import json
from pathlib import Path

def parse_solution_file(filepath):
    """
    Parse a solution file and extract D₁₁, D₁₂, D₂₂.
    Format varies by repository - implement parsers for each.
    """
    # Implementation depends on file format
    pass

def verify_complement_hypothesis(D11, D22, m):
    """
    Verify that D₂₂ = Z_m \ {0} \ D₁₁
    Returns True if hypothesis holds, False otherwise.
    """
    full_set = set(range(1, m))  # Z_m \ {0}
    expected_D22 = full_set - set(D11)
    return set(D22) == expected_D22

def analyze_D12_structure(D11, D12, m):
    """
    Analyze relationship between D₁₁ and D₁₂.
    Report: |D₁₂|/m ratio, overlap with D₁₁, algebraic structure.
    """
    overlap = set(D11) & set(D12)
    return {
        'size_ratio': len(D12) / m,
        'overlap_size': len(overlap),
        'overlap_ratio': len(overlap) / len(D11) if D11 else 0,
        'D12_subset_D11': set(D12).issubset(set(D11)),
        'D11_subset_D12': set(D11).issubset(set(D12)),
    }

def main():
    results = []
    
    # Process all known solutions
    # ... implementation ...
    
    # Output summary
    print("=== Complement Hypothesis Verification ===")
    # All should be True
    
    print("\n=== D₁₂ Structure Analysis ===")
    # Report patterns
    
    # Save to JSON for later use
    with open('known_solutions_analysis.json', 'w') as f:
        json.dump(results, f, indent=2)

if __name__ == '__main__':
    main()
```

**Expected Outcome**: Confirm D₂₂ = complement(D₁₁) holds universally for n ≤ 21.

---

## Phase 1: Core Implementation

### 1.1 File Structure

```
ramsey_books/
├── core/
│   ├── __init__.py
│   ├── circulant.py      # Circulant matrix operations
│   ├── constraints.py    # Book-avoidance constraints (CORRECTED formulas)
│   ├── verification.py   # Solution verification
│   └── utils.py          # Helper functions
├── solvers/
│   ├── __init__.py
│   ├── sat_encoder.py    # SAT encoding with symmetry breaking
│   ├── heuristic.py      # Local search / simulated annealing
│   └── hybrid.py         # Combined approach
├── analysis/
│   ├── parse_known.py    # Parse known solutions
│   └── patterns.py       # Pattern extraction
├── tests/
│   ├── test_constraints.py
│   └── test_verification.py
└── main.py
```

### 1.2 Core Constraints Module (CRITICAL)

Create `core/constraints.py`:

```python
"""
Book-avoidance constraints for 2-block circulant graphs.

CRITICAL: This module implements the CORRECTED common neighbor formulas.
The original spec contained a mathematical error.

Correct formula (Inclusion-Exclusion):
    λ_blue(u,v) = (N - 2) - deg(u) - deg(v) + λ_red(u,v)

For 2-block circulant with vertex partition V₁, V₂ (each size m = (N-2)/2 + 1 = N/2):
    - If u,v ∈ V₁: λ_blue = N - 2 - 2*d₁ + λ_red
    - If u,v ∈ V₂: λ_blue = N - 2 - 2*d₂ + λ_red  
    - If u ∈ V₁, v ∈ V₂: λ_blue = N - 2 - d₁ - d₂ + λ_red

Where d₁ = |D₁₁| + |D₁₂| and d₂ = |D₂₂| + |D₁₂| are vertex degrees.
"""

import numpy as np
from typing import Set, Tuple, List
from functools import lru_cache


def delta(A: Set[int], B: Set[int], d: int, m: int) -> int:
    """
    Compute Δ(A, B, d) = |{(a,b) ∈ A×B : a - b ≡ d (mod m)}|
    
    This counts ordered pairs where the difference is exactly d.
    Efficient O(|A|) implementation using indicator arrays.
    """
    if not A or not B:
        return 0
    
    # Create indicator array for B
    B_indicator = np.zeros(m, dtype=np.int32)
    for b in B:
        B_indicator[b % m] = 1
    
    count = 0
    for a in A:
        # We want a - b ≡ d (mod m), so b ≡ a - d (mod m)
        target = (a - d) % m
        count += B_indicator[target]
    
    return count


def sigma(A: Set[int], B: Set[int], d: int, m: int) -> int:
    """
    Compute Σ(A, B, d) = |{(a,b) ∈ A×B : a + b ≡ d (mod m)}|
    
    This counts ordered pairs where the sum is exactly d.
    """
    if not A or not B:
        return 0
    
    # Create indicator array for B
    B_indicator = np.zeros(m, dtype=np.int32)
    for b in B:
        B_indicator[b % m] = 1
    
    count = 0
    for a in A:
        # We want a + b ≡ d (mod m), so b ≡ d - a (mod m)
        target = (d - a) % m
        count += B_indicator[target]
    
    return count


def compute_lambda_red_V1V1(D11: Set[int], D12: Set[int], d: int, m: int) -> int:
    """
    Compute λ_red(u, v) for u, v ∈ V₁ with difference d.
    
    Per Wesley (2024) Lemma:
        λ_red = Δ(D₁₁, D₁₁, d) + Δ(D₁₂, D₁₂, d)
    """
    return delta(D11, D11, d, m) + delta(D12, D12, d, m)


def compute_lambda_red_V2V2(D22: Set[int], D12: Set[int], d: int, m: int) -> int:
    """
    Compute λ_red(u, v) for u, v ∈ V₂ with difference d.
    
    Per Wesley (2024) Lemma:
        λ_red = Δ(D₂₂, D₂₂, d) + Δ(D₁₂, D₁₂, d)
    
    Note: D₁₂ appears because edges from V₂ to V₁ use the same connection set.
    """
    return delta(D22, D22, d, m) + delta(D12, D12, d, m)


def compute_lambda_red_V1V2(D11: Set[int], D12: Set[int], D22: Set[int], 
                             d: int, m: int) -> int:
    """
    Compute λ_red(u, v) for u ∈ V₁, v ∈ V₂ with difference d.
    
    Per Wesley (2024) Lemma:
        λ_red = Σ(D₁₁, D₁₂, d) + Δ(D₁₂, D₂₂, d)
    
    Note the use of Σ (sum) for the first term, not Δ (difference).
    """
    return sigma(D11, D12, d, m) + delta(D12, D22, d, m)


def compute_lambda_blue(N: int, deg_u: int, deg_v: int, lambda_red: int) -> int:
    """
    Compute λ_blue(u, v) using the CORRECTED inclusion-exclusion formula.
    
    CRITICAL: This is the corrected formula. The original spec was WRONG.
    
    Correct: λ_blue = (N - 2) - deg(u) - deg(v) + λ_red
    Wrong:   λ_blue = (N - 2) - λ_red  [MISSING DEGREE TERMS]
    """
    return (N - 2) - deg_u - deg_v + lambda_red


class BookConstraintChecker:
    """
    Verifies book-avoidance constraints for a 2-block circulant graph.
    
    A graph contains B_k (book with k triangles) iff there exists an edge (u,v)
    with at least k common neighbors in the same color.
    
    For R(B_{n-1}, B_n):
        - Red edges must have < n-1 common red neighbors
        - Blue edges (non-edges) must have < n common blue neighbors
    """
    
    def __init__(self, n: int, D11: Set[int], D12: Set[int], D22: Set[int] = None):
        """
        Initialize constraint checker.
        
        Args:
            n: Parameter for R(B_{n-1}, B_n)
            D11: Connection set within V₁ (elements of Z_m \ {0})
            D12: Connection set between V₁ and V₂ (elements of Z_m)
            D22: Connection set within V₂ (default: complement of D11)
        """
        self.n = n
        self.m = 2 * n - 1  # |V₁| = |V₂| = m
        self.N = 2 * self.m  # Total vertices = 4n - 2
        
        self.D11 = set(D11)
        self.D12 = set(D12)
        
        # Enforce complement constraint if D22 not provided
        if D22 is None:
            full_set = set(range(1, self.m))
            self.D22 = full_set - self.D11
        else:
            self.D22 = set(D22)
        
        # Compute vertex degrees
        self.d1 = len(self.D11) + len(self.D12)  # degree of vertices in V₁
        self.d2 = len(self.D22) + len(self.D12)  # degree of vertices in V₂
        
        # Thresholds
        self.red_threshold = n - 2   # Red edges must have λ_red < n-1, so ≤ n-2
        self.blue_threshold = n - 1  # Blue edges must have λ_blue < n, so ≤ n-1
        
        # Cache for computed values
        self._cache = {}
    
    def _symmetrize_D11(self) -> Set[int]:
        """Ensure D11 is symmetric: d ∈ D11 ⟺ (m-d) ∈ D11"""
        sym = set()
        for d in self.D11:
            sym.add(d)
            sym.add(self.m - d)
        return sym
    
    def check_edge_V1V1(self, d: int) -> Tuple[bool, dict]:
        """
        Check constraints for edge between vertices in V₁ at difference d.
        
        Returns:
            (is_valid, details) where details contains λ_red, λ_blue, etc.
        """
        if d == 0:
            return True, {'skip': 'self-loop'}
        
        # Determine edge color
        is_red = d in self.D11
        
        # Compute λ_red
        lambda_red = compute_lambda_red_V1V1(self.D11, self.D12, d, self.m)
        
        # Compute λ_blue (CORRECTED FORMULA)
        lambda_blue = compute_lambda_blue(self.N, self.d1, self.d1, lambda_red)
        
        details = {
            'd': d,
            'block': 'V1V1',
            'is_red': is_red,
            'lambda_red': lambda_red,
            'lambda_blue': lambda_blue,
            'deg_u': self.d1,
            'deg_v': self.d1,
        }
        
        if is_red:
            is_valid = lambda_red <= self.red_threshold
            details['constraint'] = f'λ_red ≤ {self.red_threshold}'
            details['satisfied'] = is_valid
        else:
            is_valid = lambda_blue <= self.blue_threshold
            details['constraint'] = f'λ_blue ≤ {self.blue_threshold}'
            details['satisfied'] = is_valid
        
        return is_valid, details
    
    def check_edge_V2V2(self, d: int) -> Tuple[bool, dict]:
        """
        Check constraints for edge between vertices in V₂ at difference d.
        """
        if d == 0:
            return True, {'skip': 'self-loop'}
        
        is_red = d in self.D22
        
        lambda_red = compute_lambda_red_V2V2(self.D22, self.D12, d, self.m)
        lambda_blue = compute_lambda_blue(self.N, self.d2, self.d2, lambda_red)
        
        details = {
            'd': d,
            'block': 'V2V2',
            'is_red': is_red,
            'lambda_red': lambda_red,
            'lambda_blue': lambda_blue,
            'deg_u': self.d2,
            'deg_v': self.d2,
        }
        
        if is_red:
            is_valid = lambda_red <= self.red_threshold
            details['constraint'] = f'λ_red ≤ {self.red_threshold}'
            details['satisfied'] = is_valid
        else:
            is_valid = lambda_blue <= self.blue_threshold
            details['constraint'] = f'λ_blue ≤ {self.blue_threshold}'
            details['satisfied'] = is_valid
        
        return is_valid, details
    
    def check_edge_V1V2(self, d: int) -> Tuple[bool, dict]:
        """
        Check constraints for edge between V₁ and V₂ at difference d.
        """
        is_red = d in self.D12
        
        lambda_red = compute_lambda_red_V1V2(self.D11, self.D12, self.D22, d, self.m)
        lambda_blue = compute_lambda_blue(self.N, self.d1, self.d2, lambda_red)
        
        details = {
            'd': d,
            'block': 'V1V2',
            'is_red': is_red,
            'lambda_red': lambda_red,
            'lambda_blue': lambda_blue,
            'deg_u': self.d1,
            'deg_v': self.d2,
        }
        
        if is_red:
            is_valid = lambda_red <= self.red_threshold
            details['constraint'] = f'λ_red ≤ {self.red_threshold}'
            details['satisfied'] = is_valid
        else:
            is_valid = lambda_blue <= self.blue_threshold
            details['constraint'] = f'λ_blue ≤ {self.blue_threshold}'
            details['satisfied'] = is_valid
        
        return is_valid, details
    
    def check_all_constraints(self) -> Tuple[bool, List[dict]]:
        """
        Check all book-avoidance constraints.
        
        Returns:
            (all_valid, list of violation details)
        """
        violations = []
        
        # Check V₁-V₁ edges (d from 1 to m-1, symmetric so only need 1 to m//2)
        for d in range(1, self.m):
            valid, details = self.check_edge_V1V1(d)
            if not valid:
                violations.append(details)
        
        # Check V₂-V₂ edges
        for d in range(1, self.m):
            valid, details = self.check_edge_V2V2(d)
            if not valid:
                violations.append(details)
        
        # Check V₁-V₂ edges (d from 0 to m-1)
        for d in range(self.m):
            valid, details = self.check_edge_V1V2(d)
            if not valid:
                violations.append(details)
        
        return len(violations) == 0, violations
    
    def get_violation_count(self) -> int:
        """Quick count of constraint violations (for optimization)."""
        _, violations = self.check_all_constraints()
        return len(violations)


def verify_known_solution(n: int, D11: Set[int], D12: Set[int], 
                          D22: Set[int] = None) -> bool:
    """
    Verify a known solution satisfies all constraints.
    Use this to validate implementation against archived solutions.
    """
    checker = BookConstraintChecker(n, D11, D12, D22)
    is_valid, violations = checker.check_all_constraints()
    
    if not is_valid:
        print(f"VERIFICATION FAILED for n={n}")
        print(f"Found {len(violations)} violations:")
        for v in violations[:5]:  # Show first 5
            print(f"  {v}")
        if len(violations) > 5:
            print(f"  ... and {len(violations) - 5} more")
    
    return is_valid
```

### 1.3 SAT Encoder with Symmetry Breaking

Create `solvers/sat_encoder.py`:

```python
"""
SAT encoding for 2-block circulant book Ramsey search.

Key optimizations:
1. Enforce D₂₂ = complement(D₁₁) as hard constraint (eliminates D₂₂ variables)
2. Symmetry breaking: fix x₁ = True, enforce x_d → x_{m-d} for d > m/2
3. Incremental solving: warm-start from smaller n

Variable scheme for n=22 (m=43):
- x_1, x_2, ..., x_21: D₁₁ membership (x_d ⟺ d ∈ D₁₁, only d ≤ m/2 due to symmetry)
- y_0, y_1, ..., y_42: D₁₂ membership (y_d ⟺ d ∈ D₁₂)

Total: 21 + 43 = 64 variables (down from ~107 without optimizations)
"""

from typing import Set, List, Tuple, Optional
from pysat.solvers import Glucose4
from pysat.formula import CNF


class BookRamseySATEncoder:
    """
    Encodes book Ramsey constraints as SAT clauses.
    """
    
    def __init__(self, n: int):
        self.n = n
        self.m = 2 * n - 1
        self.N = 2 * self.m
        
        # Variable indices
        # D₁₁: variables 1 to m//2 (symmetry: x_d ⟺ x_{m-d})
        # D₁₂: variables (m//2 + 1) to (m//2 + m)
        self.d11_offset = 0
        self.d12_offset = self.m // 2
        
        self.cnf = CNF()
        self.var_count = 0
    
    def d11_var(self, d: int) -> int:
        """Get variable index for D₁₁ membership of d."""
        if d == 0:
            raise ValueError("0 is never in D₁₁")
        # Use symmetry: variable for min(d, m-d)
        d_sym = min(d, self.m - d)
        return self.d11_offset + d_sym
    
    def d12_var(self, d: int) -> int:
        """Get variable index for D₁₂ membership of d."""
        return self.d12_offset + d + 1
    
    def d22_literal(self, d: int) -> int:
        """
        Get literal for D₂₂ membership of d.
        Since D₂₂ = complement(D₁₁), d ∈ D₂₂ ⟺ d ∉ D₁₁ ⟺ ¬x_d
        """
        if d == 0:
            raise ValueError("0 is never in D₂₂")
        return -self.d11_var(d)
    
    def add_symmetry_breaking(self):
        """
        Add symmetry breaking clauses:
        1. Fix x_1 = True (break cyclic rotation symmetry)
        2. For d > m/2: enforce x_d ↔ x_{m-d} (already implicit in variable scheme)
        """
        # Fix first element
        self.cnf.append([self.d11_var(1)])  # x_1 = True
    
    def encode_V1V1_constraints(self):
        """
        Encode constraints for V₁-V₁ edges.
        
        For each difference d ∈ [1, m-1]:
        - If d ∈ D₁₁ (red edge): λ_red ≤ n-2
        - If d ∉ D₁₁ (blue edge): λ_blue ≤ n-1
        
        λ_red = Δ(D₁₁, D₁₁, d) + Δ(D₁₂, D₁₂, d)
        λ_blue = (N-2) - 2*d₁ + λ_red
        """
        # This requires cardinality constraints
        # Use pseudo-boolean encoding or sequential counter
        # For simplicity, enumerate and add implications
        
        # ... complex encoding logic ...
        pass
    
    def encode_all_constraints(self):
        """Build complete SAT formula."""
        self.add_symmetry_breaking()
        self.encode_V1V1_constraints()
        # ... V2V2 and V1V2 constraints ...
    
    def solve(self, timeout: int = None) -> Optional[Tuple[Set[int], Set[int]]]:
        """
        Solve the SAT instance.
        
        Returns:
            (D₁₁, D₁₂) if satisfiable, None otherwise
        """
        with Glucose4(bootstrap_with=self.cnf) as solver:
            if solver.solve():
                model = solver.get_model()
                D11, D12 = self._decode_model(model)
                return D11, D12
            return None
    
    def _decode_model(self, model: List[int]) -> Tuple[Set[int], Set[int]]:
        """Extract D₁₁ and D₁₂ from SAT model."""
        D11 = set()
        D12 = set()
        
        for d in range(1, self.m // 2 + 1):
            var = self.d11_var(d)
            if model[var - 1] > 0:
                D11.add(d)
                D11.add(self.m - d)  # Add symmetric element
        
        for d in range(self.m):
            var = self.d12_var(d)
            if model[var - 1] > 0:
                D12.add(d)
        
        return D11, D12
```

### 1.4 Heuristic Solver with Algebraic Seeding

Create `solvers/heuristic.py`:

```python
"""
Local search / simulated annealing solver with algebraic seeding.

Key insight: Seed with Paley-like structures even though 43 ≡ 3 (mod 4).
The quadratic residues of F₄₃ may provide a good starting point for perturbation.
"""

import numpy as np
from typing import Set, Tuple, Optional
import random

from core.constraints import BookConstraintChecker


def quadratic_residues(p: int) -> Set[int]:
    """
    Compute quadratic residues mod p.
    QR_p = {x² mod p : x ∈ {1, ..., p-1}}
    """
    return {(x * x) % p for x in range(1, p)}


def algebraic_seed_D11(m: int) -> Set[int]:
    """
    Generate initial D₁₁ using quadratic residues.
    
    For m=43 (prime ≡ 3 mod 4):
    QR₄₃ = {1, 4, 6, 9, 10, 11, 13, 14, 15, 16, 17, 21, 
            23, 24, 25, 31, 35, 36, 38, 40, 41}
    
    Note: Since 43 ≡ 3 (mod 4), -1 is NOT a QR, so the Paley graph
    would be directed. We use QR anyway as a starting point.
    """
    if m % 4 == 3:
        # Paley doesn't work directly, but QRs still provide structure
        return quadratic_residues(m)
    else:
        # Standard Paley seeding for m ≡ 1 (mod 4)
        return quadratic_residues(m)


def algebraic_seed_D12(m: int, D11: Set[int]) -> Set[int]:
    """
    Generate initial D₁₂ based on D₁₁.
    
    Heuristic from known solutions: D₁₂ often equals or is close to D₁₁.
    """
    # Start with D₁₂ = D₁₁ ∪ {0}
    D12 = D11.copy()
    D12.add(0)
    return D12


class LocalSearchSolver:
    """
    Simulated annealing solver for book Ramsey search.
    """
    
    def __init__(self, n: int, seed: int = None):
        self.n = n
        self.m = 2 * n - 1
        
        if seed is not None:
            random.seed(seed)
            np.random.seed(seed)
    
    def objective(self, D11: Set[int], D12: Set[int]) -> int:
        """
        Objective: number of constraint violations.
        Goal: minimize to 0.
        """
        checker = BookConstraintChecker(self.n, D11, D12)
        return checker.get_violation_count()
    
    def random_neighbor(self, D11: Set[int], D12: Set[int]) -> Tuple[Set[int], Set[int]]:
        """
        Generate random neighbor by flipping one element.
        
        Moves:
        1. Add/remove element from D₁₁ (maintaining symmetry)
        2. Add/remove element from D₁₂
        """
        D11_new = D11.copy()
        D12_new = D12.copy()
        
        if random.random() < 0.5:
            # Modify D₁₁
            candidates = list(range(1, self.m // 2 + 1))
            d = random.choice(candidates)
            if d in D11_new:
                D11_new.discard(d)
                D11_new.discard(self.m - d)
            else:
                D11_new.add(d)
                D11_new.add(self.m - d)
        else:
            # Modify D₁₂
            d = random.randint(0, self.m - 1)
            if d in D12_new:
                D12_new.discard(d)
            else:
                D12_new.add(d)
        
        return D11_new, D12_new
    
    def solve(self, max_iter: int = 1000000, 
              initial_temp: float = 1.0,
              cooling_rate: float = 0.99999) -> Optional[Tuple[Set[int], Set[int]]]:
        """
        Run simulated annealing.
        """
        # Initialize with algebraic seed
        D11 = algebraic_seed_D11(self.m)
        D12 = algebraic_seed_D12(self.m, D11)
        
        current_obj = self.objective(D11, D12)
        best_obj = current_obj
        best_D11, best_D12 = D11.copy(), D12.copy()
        
        temp = initial_temp
        
        for iteration in range(max_iter):
            if current_obj == 0:
                print(f"Solution found at iteration {iteration}!")
                return D11, D12
            
            # Generate neighbor
            D11_new, D12_new = self.random_neighbor(D11, D12)
            new_obj = self.objective(D11_new, D12_new)
            
            # Accept/reject
            delta = new_obj - current_obj
            if delta < 0 or random.random() < np.exp(-delta / temp):
                D11, D12 = D11_new, D12_new
                current_obj = new_obj
                
                if current_obj < best_obj:
                    best_obj = current_obj
                    best_D11, best_D12 = D11.copy(), D12.copy()
                    print(f"Iteration {iteration}: new best = {best_obj}")
            
            # Cool down
            temp *= cooling_rate
            
            if iteration % 10000 == 0:
                print(f"Iteration {iteration}, temp={temp:.6f}, "
                      f"current={current_obj}, best={best_obj}")
        
        print(f"No solution found. Best objective: {best_obj}")
        return None
```

---

## Phase 2: Validation

### 2.1 Test Against Known Solutions

Create `tests/test_constraints.py`:

```python
"""
Test constraint checker against known solutions.

CRITICAL: All known solutions for n ≤ 21 must pass validation.
If any fail, the implementation is WRONG.
"""

import pytest
from core.constraints import BookConstraintChecker, verify_known_solution


# Known solution for n=6 (from Wesley 2024)
# R(B_5, B_6) = 23, so N = 22, m = 11
SOLUTION_N6 = {
    'D11': {1, 3, 4, 5, 9},  # Example - replace with actual
    'D12': {0, 1, 3, 4, 5, 9},
}

# Known solution for n=10 (from Wesley 2024)
# R(B_9, B_10) = 39, so N = 38, m = 19
SOLUTION_N10 = {
    'D11': {1, 4, 5, 6, 7, 9, 11, 16, 17},  # Example - replace with actual
    'D12': {0, 1, 4, 5, 6, 7, 9, 11, 16, 17},
}


def test_n6_solution():
    """Validate n=6 solution."""
    assert verify_known_solution(6, 
        SOLUTION_N6['D11'], 
        SOLUTION_N6['D12']
    )


def test_n10_solution():
    """Validate n=10 solution."""
    assert verify_known_solution(10,
        SOLUTION_N10['D11'],
        SOLUTION_N10['D12']
    )


def test_complement_constraint():
    """Verify D₂₂ = complement(D₁₁) is enforced."""
    n = 10
    m = 2 * n - 1
    D11 = {1, 4, 5, 6, 7, 9, 11, 16, 17}
    D12 = {0, 1, 4, 5, 6, 7, 9, 11, 16, 17}
    
    checker = BookConstraintChecker(n, D11, D12)
    
    # D₂₂ should be automatically computed as complement
    expected_D22 = set(range(1, m)) - D11
    assert checker.D22 == expected_D22


def test_blue_formula_correctness():
    """
    CRITICAL TEST: Verify blue neighbor formula includes degree terms.
    
    λ_blue = (N - 2) - deg(u) - deg(v) + λ_red  [CORRECT]
    λ_blue = (N - 2) - λ_red                    [WRONG - old spec]
    """
    from core.constraints import compute_lambda_blue
    
    N = 86  # n=22
    deg_u = 40
    deg_v = 40
    lambda_red = 10
    
    correct = compute_lambda_blue(N, deg_u, deg_v, lambda_red)
    
    # Correct formula: 84 - 40 - 40 + 10 = 14
    assert correct == 14
    
    # Wrong formula would give: 84 - 10 = 74 (MASSIVELY WRONG)
    wrong = (N - 2) - lambda_red
    assert wrong == 74
    assert correct != wrong


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
```

### 2.2 Run Validation Suite

```bash
# Run all tests
cd ~/ramsey_books
python -m pytest tests/ -v

# All tests MUST pass before proceeding to search
```

---

## Phase 3: Search Execution

### 3.1 Main Search Script

Create `main.py`:

```python
"""
Main entry point for R(B₂₁, B₂₂) search.
"""

import argparse
import json
import time
from datetime import datetime

from solvers.heuristic import LocalSearchSolver
from solvers.sat_encoder import BookRamseySATEncoder
from core.constraints import verify_known_solution


def run_heuristic_search(n: int, max_iter: int = 10000000):
    """Run heuristic search."""
    print(f"=== Heuristic Search for n={n} ===")
    print(f"Target: R(B_{n-1}, B_{n}) = {4*n - 1}")
    print(f"Graph order: N = {4*n - 2}")
    print(f"Block size: m = {2*n - 1}")
    print()
    
    solver = LocalSearchSolver(n)
    start = time.time()
    result = solver.solve(max_iter=max_iter)
    elapsed = time.time() - start
    
    if result:
        D11, D12 = result
        print(f"\n=== SOLUTION FOUND in {elapsed:.2f}s ===")
        print(f"D₁₁ = {sorted(D11)}")
        print(f"D₁₂ = {sorted(D12)}")
        
        # Verify
        if verify_known_solution(n, D11, D12):
            print("✓ Solution VERIFIED")
            
            # Save
            solution = {
                'n': n,
                'timestamp': datetime.now().isoformat(),
                'elapsed_seconds': elapsed,
                'D11': sorted(D11),
                'D12': sorted(D12),
            }
            filename = f'solution_n{n}_{datetime.now().strftime("%Y%m%d_%H%M%S")}.json'
            with open(filename, 'w') as f:
                json.dump(solution, f, indent=2)
            print(f"Saved to {filename}")
        else:
            print("✗ Verification FAILED - BUG IN SOLVER")
    else:
        print(f"\nNo solution found in {elapsed:.2f}s")


def run_sat_search(n: int, timeout: int = 3600):
    """Run SAT-based search."""
    print(f"=== SAT Search for n={n} ===")
    
    encoder = BookRamseySATEncoder(n)
    encoder.encode_all_constraints()
    
    print(f"Variables: {encoder.var_count}")
    print(f"Clauses: {len(encoder.cnf.clauses)}")
    
    result = encoder.solve(timeout=timeout)
    
    if result:
        D11, D12 = result
        print(f"D₁₁ = {sorted(D11)}")
        print(f"D₁₂ = {sorted(D12)}")
    else:
        print("UNSAT or timeout")


def main():
    parser = argparse.ArgumentParser(description='Book Ramsey Number Search')
    parser.add_argument('--n', type=int, default=22, help='Parameter n for R(B_{n-1}, B_n)')
    parser.add_argument('--method', choices=['heuristic', 'sat', 'both'], default='heuristic')
    parser.add_argument('--max-iter', type=int, default=10000000)
    parser.add_argument('--timeout', type=int, default=3600)
    
    args = parser.parse_args()
    
    if args.method in ['heuristic', 'both']:
        run_heuristic_search(args.n, args.max_iter)
    
    if args.method in ['sat', 'both']:
        run_sat_search(args.n, args.timeout)


if __name__ == '__main__':
    main()
```

### 3.2 Execution Commands

```bash
# First run: heuristic with algebraic seeding
python main.py --n 22 --method heuristic --max-iter 50000000

# Alternative: SAT solver (may be slower but complete)
python main.py --n 22 --method sat --timeout 86400

# Parallel runs with different seeds
for seed in 1 2 3 4 5 6 7 8; do
    python main.py --n 22 --method heuristic --seed $seed &
done
wait
```

---

## Phase 4: Post-Discovery Analysis

If a solution is found:

1. **Verify independently**: Re-check all constraints manually
2. **Generate full adjacency matrix**: Construct the 86×86 matrix
3. **Check for algebraic structure**: Is D₁₁ a union of cyclotomic cosets?
4. **Test generalizability**: Does the pattern extend to n=23?
5. **Document**: Prepare for submission to House of Graphs

---

## Summary of Critical Fixes

| Original Spec | Correction | Impact |
|---------------|------------|--------|
| Target n=18 | n=22 | n=18 already solved |
| Target n=27 | Skip | Paley construction works |
| λ_blue = (N-2) - λ_red | λ_blue = (N-2) - deg(u) - deg(v) + λ_red | Catastrophic error |
| Search over D₁₁, D₁₂, D₂₂ | Enforce D₂₂ = complement(D₁₁) | ~2^21 search space reduction |
| No symmetry breaking | Fix x₁ = True, use symmetric variables | Further reduction |

---

## Expected Resource Requirements

- **Variables**: 64 (21 for D₁₁ with symmetry, 43 for D₁₂)
- **Constraints**: O(m³) ≈ 80,000 clauses for SAT
- **Heuristic runtime**: Hours to days
- **SAT runtime**: Days to weeks (with proper solver)

Good luck! This is genuine frontier mathematics.
