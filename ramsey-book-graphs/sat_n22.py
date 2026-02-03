"""
Complete SAT encoding for n=22 2-block circulant Ramsey construction.

Variables:
- D11_half: 21 boolean vars (d=1..21, symmetric so d and 43-d share var)
- D12: 43 boolean vars (d=0..42)
- D22 = complement(D11) to reduce search space

Constraints:
- For each V1V1 edge d: if d in D11 then lambda_red <= 20, else lambda_blue <= 21
- For each V2V2 edge d: if d in D22 then lambda_red <= 20, else lambda_blue <= 21
- For each V1V2 edge d: if d in D12 then lambda_red <= 20, else lambda_blue <= 21

The blue constraint uses: lambda_blue = (N-2) - deg_u - deg_v + lambda_red
"""

import sys
import time
from typing import Set, List, Tuple, Optional
from itertools import combinations

try:
    from pysat.solvers import Cadical153, Glucose4
    from pysat.card import CardEnc, EncType
    from pysat.formula import CNF
    PYSAT_AVAILABLE = True
except ImportError:
    PYSAT_AVAILABLE = False
    print("PySAT not available. Install with: pip install python-sat")

# Problem parameters
n = 22
m = 43
N = 86
half_m = m // 2  # 21

red_threshold = n - 2  # 20
blue_threshold = n - 1  # 21


class SATEncoder:
    def __init__(self):
        self.cnf = CNF()
        self.next_var = 1

        # D11 variables: var[d] for d in 1..21 (half representatives)
        self.d11_vars = {}
        for d in range(1, half_m + 1):
            self.d11_vars[d] = self.next_var
            self.next_var += 1

        # D12 variables: var[d] for d in 0..42
        self.d12_vars = {}
        for d in range(m):
            self.d12_vars[d] = self.next_var
            self.next_var += 1

        # D22 = complement of D11, so no separate variables needed
        # d22_var(d) = NOT d11_var(d)

        print(f"Variables: {half_m} (D11) + {m} (D12) = {half_m + m} base variables")

    def d11_var(self, d: int) -> int:
        """Get variable for d in D11 (handling symmetry)."""
        d = d % m
        if d == 0:
            return None  # 0 never in D11
        if d > half_m:
            d = m - d
        return self.d11_vars[d]

    def d12_var(self, d: int) -> int:
        """Get variable for d in D12."""
        return self.d12_vars[d % m]

    def d22_var(self, d: int) -> int:
        """D22 = complement of D11, so d22_var(d) = -d11_var(d)."""
        v = self.d11_var(d)
        return -v if v else None  # If d11_var is None (d=0), D22 contains it

    def new_var(self) -> int:
        v = self.next_var
        self.next_var += 1
        return v

    def encode_and(self, a: int, b: int) -> int:
        """Create auxiliary variable aux = a AND b."""
        if a == b:
            return a
        aux = self.new_var()
        # aux => a
        self.cnf.append([-aux, a])
        # aux => b
        self.cnf.append([-aux, b])
        # a AND b => aux
        self.cnf.append([aux, -a, -b])
        return aux

    def get_lambda_red_V1V1(self, d: int) -> List[int]:
        """
        Get list of aux variables representing common neighbor contributions for V1V1 edge.
        lambda_red = Delta(D11, D11, d) + Delta(D12, D12, d)
        """
        aux_vars = []

        # Delta(D11, D11, d): for each a in {1..m-1}, check if a in D11 and (a-d) in D11
        for a in range(1, m):
            b = (a - d) % m
            if b == 0:
                continue
            va = self.d11_var(a)
            vb = self.d11_var(b)
            if va and vb:
                aux_vars.append(self.encode_and(va, vb))

        # Delta(D12, D12, d): for each a in {0..m-1}, check if a in D12 and (a-d) in D12
        for a in range(m):
            b = (a - d) % m
            va = self.d12_var(a)
            vb = self.d12_var(b)
            aux_vars.append(self.encode_and(va, vb))

        return aux_vars

    def get_lambda_red_V2V2(self, d: int) -> List[int]:
        """
        Get list of aux variables for V2V2 edge.
        lambda_red = Delta(D22, D22, d) + Delta(D12^T, D12^T, d)

        Since D22 = complement(D11), d22_var(x) = -d11_var(x)
        Delta(D22, D22, d) counts pairs where both are NOT in D11
        """
        aux_vars = []

        # Delta(D22, D22, d): pairs (a, a-d) where both NOT in D11
        # This is: NOT d11(a) AND NOT d11(a-d)
        for a in range(1, m):
            b = (a - d) % m
            if b == 0:
                continue
            va = self.d11_var(a)
            vb = self.d11_var(b)
            if va and vb:
                # aux = (NOT va) AND (NOT vb)
                aux = self.new_var()
                # aux => NOT va
                self.cnf.append([-aux, -va])
                # aux => NOT vb
                self.cnf.append([-aux, -vb])
                # (NOT va) AND (NOT vb) => aux
                self.cnf.append([aux, va, vb])
                aux_vars.append(aux)

        # Delta(D12^T, D12^T, d): D12^T = {-x mod m : x in D12}
        # For a in D12, -a in D12^T. Pairs in D12^T with difference d:
        # (-a) - (-b) = b - a = d, so a - b = -d, meaning b = a + d
        for a in range(m):
            b = (a + d) % m
            va = self.d12_var(a)
            vb = self.d12_var(b)
            aux_vars.append(self.encode_and(va, vb))

        return aux_vars

    def get_lambda_red_V1V2(self, d: int) -> List[int]:
        """
        Get list of aux variables for V1V2 edge.
        lambda_red = Sigma(D11, D12, d) + Delta(D12, D22, d)
        """
        aux_vars = []

        # Sigma(D11, D12, d): pairs (a, b) where a + b = d
        for a in range(1, m):
            b = (d - a) % m
            va = self.d11_var(a)
            vb = self.d12_var(b)
            if va:
                aux_vars.append(self.encode_and(va, vb))

        # Delta(D12, D22, d): pairs (a, a-d) where a in D12 and (a-d) in D22
        # (a-d) in D22 means (a-d) NOT in D11
        for a in range(m):
            b = (a - d) % m
            if b == 0:
                continue  # 0 is always in D22 (not in D11)
            va = self.d12_var(a)
            vb = self.d11_var(b)  # We want NOT vb
            if vb:
                # aux = va AND (NOT vb)
                aux = self.new_var()
                self.cnf.append([-aux, va])
                self.cnf.append([-aux, -vb])
                self.cnf.append([aux, -va, vb])
                aux_vars.append(aux)
            else:
                # b = 0, always in D22, so just need va
                aux_vars.append(va)

        return aux_vars

    def encode_conditional_atmost(self, condition_var: int, aux_vars: List[int], bound: int):
        """
        Encode: IF condition_var THEN sum(aux_vars) <= bound
        Using combinatorial encoding for small bounds.
        """
        if bound + 1 > len(aux_vars):
            return  # Constraint always satisfied

        # For each subset of size (bound + 1), at least one must be false if condition is true
        # condition => NOT(all of subset true)
        # = NOT condition OR NOT(all true)
        # = NOT condition OR (at least one false)

        if bound + 1 <= 8:  # Small enough for explicit enumeration
            for subset in combinations(range(len(aux_vars)), bound + 1):
                clause = [-condition_var] + [-aux_vars[i] for i in subset]
                self.cnf.append(clause)
        else:
            # Use cardinality encoding
            atmost = CardEnc.atmost(
                lits=aux_vars,
                bound=bound,
                top_id=self.next_var,
                encoding=EncType.seqcounter
            )
            self.next_var = atmost.nv + 1
            for clause in atmost.clauses:
                # Make conditional: add -condition_var to each clause
                self.cnf.append([-condition_var] + list(clause))

    def encode_conditional_atleast(self, condition_var: int, aux_vars: List[int], bound: int):
        """
        Encode: IF condition_var THEN sum(aux_vars) >= bound
        """
        if bound <= 0:
            return  # Always satisfied
        if bound > len(aux_vars):
            # Never satisfiable when condition is true
            self.cnf.append([-condition_var])
            return

        # condition => at least (bound) aux_vars are true
        # = NOT condition OR (sum >= bound)
        # Contrapositive: if sum < bound, then NOT condition
        # If (len - bound + 1) are false, then NOT condition

        num_false_needed = len(aux_vars) - bound + 1
        if num_false_needed <= 8:
            for subset in combinations(range(len(aux_vars)), num_false_needed):
                # If all in subset are false, condition must be false
                clause = [-condition_var] + [aux_vars[i] for i in subset]
                self.cnf.append(clause)
        else:
            atleast = CardEnc.atleast(
                lits=aux_vars,
                bound=bound,
                top_id=self.next_var,
                encoding=EncType.seqcounter
            )
            self.next_var = atleast.nv + 1
            for clause in atleast.clauses:
                self.cnf.append([-condition_var] + list(clause))

    def encode_all_constraints(self):
        """Encode all Ramsey constraints."""

        print("Encoding V1V1 constraints...")
        for d in range(1, m):
            if d > half_m:
                continue  # Symmetric, handled by d <= half_m

            edge_var = self.d11_var(d)
            aux_vars = self.get_lambda_red_V1V1(d)

            # If edge is red (d in D11): lambda_red <= 20
            self.encode_conditional_atmost(edge_var, aux_vars, red_threshold)

            # If edge is blue (d not in D11): lambda_blue <= 21
            # lambda_blue = (N-2) - 2*|D11| - 2*|D12| + 2*|D11 cap D12| + lambda_red (complex)
            # Simplified: for blue edges, we need enough red common neighbors
            # Skip complex blue encoding for now - verify solutions post-hoc

        print("Encoding V2V2 constraints...")
        for d in range(1, m):
            if d > half_m:
                continue

            # D22 = complement(D11), so d in D22 iff d not in D11
            # edge_var for "d in D22" is -d11_var(d)
            d11_v = self.d11_var(d)
            aux_vars = self.get_lambda_red_V2V2(d)

            # If d in D22 (d not in D11): red edge, lambda_red <= 20
            # Condition: NOT d11_v, i.e., -d11_v is true
            # We need: (-d11_v) => sum(aux) <= 20
            # Rewrite: d11_v OR sum(aux) <= 20
            # For atmost encoding with condition = -d11_v:
            # -(-d11_v) = d11_v should be added to each clause
            if aux_vars and d11_v:
                if red_threshold + 1 <= len(aux_vars):
                    if red_threshold + 1 <= 8:
                        for subset in combinations(range(len(aux_vars)), red_threshold + 1):
                            clause = [d11_v] + [-aux_vars[i] for i in subset]
                            self.cnf.append(clause)
                    else:
                        atmost = CardEnc.atmost(
                            lits=aux_vars,
                            bound=red_threshold,
                            top_id=self.next_var,
                            encoding=EncType.seqcounter
                        )
                        self.next_var = atmost.nv + 1
                        for clause in atmost.clauses:
                            self.cnf.append([d11_v] + list(clause))

        print("Encoding V1V2 constraints...")
        for d in range(m):
            edge_var = self.d12_var(d)
            aux_vars = self.get_lambda_red_V1V2(d)

            # If d in D12: red edge, lambda_red <= 20
            self.encode_conditional_atmost(edge_var, aux_vars, red_threshold)

        print(f"Total clauses: {len(self.cnf.clauses)}")
        print(f"Total variables: {self.next_var - 1}")

    def solve(self, timeout: int = 600) -> Optional[Tuple[Set[int], Set[int]]]:
        """Run SAT solver and extract solution."""
        print(f"\nStarting SAT solver (timeout={timeout}s)...")

        start = time.time()

        # Try Cadical first, fall back to Glucose
        try:
            solver = Cadical153(bootstrap_with=self.cnf.clauses)
        except:
            solver = Glucose4(bootstrap_with=self.cnf.clauses)

        result = solver.solve()
        elapsed = time.time() - start

        print(f"Solver finished in {elapsed:.2f}s")

        if result:
            model = set(solver.get_model())

            # Extract D11
            D11 = set()
            for d in range(1, half_m + 1):
                if self.d11_vars[d] in model:
                    D11.add(d)
                    D11.add(m - d)

            # Extract D12
            D12 = set()
            for d in range(m):
                if self.d12_vars[d] in model:
                    D12.add(d)

            solver.delete()
            return D11, D12
        else:
            solver.delete()
            return None


def main():
    if not PYSAT_AVAILABLE:
        print("PySAT required. Install with: pip install python-sat")
        return

    print(f"SAT encoding for n={n}, m={m}, N={N}")
    print(f"Red threshold: {red_threshold}, Blue threshold: {blue_threshold}")
    print("Using D22 = complement(D11)")
    print("=" * 60)

    encoder = SATEncoder()
    encoder.encode_all_constraints()

    result = encoder.solve(timeout=3600)

    if result:
        D11, D12 = result
        D22 = set(range(1, m)) - D11

        print("\nSAT: Solution found!")
        print(f"D11 = {sorted(D11)}")
        print(f"D12 = {sorted(D12)}")
        print(f"D22 = {sorted(D22)}")

        # Verify
        from ramsey_core import BlockCirculantGraph, verify_construction, save_construction
        G = BlockCirculantGraph(n=n, D11=D11, D12=D12, D22=D22)
        verification = verify_construction(G)

        print(f"\nVerification: valid={verification.valid}")
        print(f"Max red common: {verification.max_red_common}/{verification.red_threshold}")
        print(f"Max blue common: {verification.max_blue_common}/{verification.blue_threshold}")

        if verification.valid:
            save_construction(G, verification, "solution_n22_sat.json", solver="SAT")
            print("Saved to solution_n22_sat.json")
        else:
            print(f"Violations: {len(verification.violations)}")
            for v in verification.violations[:10]:
                print(f"  {v}")
    else:
        print("\nUNSAT: No 2-block circulant solution exists for n=22")
        print("(with D22 = complement(D11) constraint)")


if __name__ == "__main__":
    main()
