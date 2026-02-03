"""
Complete SAT encoding for n=22 with both red and blue constraints.

Key insight: Blue constraint depends on vertex degrees.
For V1V1 blue edge: lambda_blue = (N-2) - 2*d1 + lambda_red where d1 = |D11| + |D12|
For V2V2 blue edge: lambda_blue = (N-2) - 2*d2 + lambda_red where d2 = |D22| + |D12|
For V1V2 blue edge: lambda_blue = (N-2) - d1 - d2 + lambda_red

Strategy: Fix cardinalities based on heuristic best solution and encode.
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
    print("PySAT required. Install with: pip install python-sat")

n = 22
m = 43
N = 86
half_m = m // 2

red_threshold = n - 2  # 20
blue_threshold = n - 1  # 21


class SATEncoderV2:
    def __init__(self, target_d11_size: int, target_d12_size: int):
        """
        target_d11_size: target for |D11| (full, not half)
        target_d12_size: target for |D12|
        """
        self.cnf = CNF()
        self.next_var = 1
        self.target_d11_size = target_d11_size
        self.target_d12_size = target_d12_size

        # With D22 = complement(D11):
        # |D22| = (m-1) - |D11| = 42 - target_d11_size

        # Degrees
        self.d1 = target_d11_size + target_d12_size
        self.d2 = (42 - target_d11_size) + target_d12_size

        # D11 variables
        self.d11_vars = {}
        for d in range(1, half_m + 1):
            self.d11_vars[d] = self.next_var
            self.next_var += 1

        # D12 variables
        self.d12_vars = {}
        for d in range(m):
            self.d12_vars[d] = self.next_var
            self.next_var += 1

        print(f"Target |D11|={target_d11_size}, |D12|={target_d12_size}")
        print(f"Degrees: d1={self.d1}, d2={self.d2}")
        print(f"Blue constraint for V1V1: lambda_red >= {2*self.d1 - 84 + blue_threshold} = {2*self.d1 - 63}")
        print(f"Blue constraint for V2V2: lambda_red >= {2*self.d2 - 63}")
        print(f"Blue constraint for V1V2: lambda_red >= {self.d1 + self.d2 - 63}")

    def d11_var(self, d: int) -> int:
        d = d % m
        if d == 0:
            return None
        if d > half_m:
            d = m - d
        return self.d11_vars[d]

    def d12_var(self, d: int) -> int:
        return self.d12_vars[d % m]

    def new_var(self) -> int:
        v = self.next_var
        self.next_var += 1
        return v

    def encode_and(self, a: int, b: int) -> int:
        if a == b:
            return a
        aux = self.new_var()
        self.cnf.append([-aux, a])
        self.cnf.append([-aux, b])
        self.cnf.append([aux, -a, -b])
        return aux

    def encode_and_not(self, a: int, b: int) -> int:
        """Encode a AND (NOT b)"""
        aux = self.new_var()
        self.cnf.append([-aux, a])
        self.cnf.append([-aux, -b])
        self.cnf.append([aux, -a, b])
        return aux

    def encode_not_and_not(self, a: int, b: int) -> int:
        """Encode (NOT a) AND (NOT b)"""
        aux = self.new_var()
        self.cnf.append([-aux, -a])
        self.cnf.append([-aux, -b])
        self.cnf.append([aux, a, b])
        return aux

    def get_lambda_red_V1V1(self, d: int) -> List[int]:
        aux_vars = []
        for a in range(1, m):
            b = (a - d) % m
            if b == 0:
                continue
            va = self.d11_var(a)
            vb = self.d11_var(b)
            if va and vb:
                aux_vars.append(self.encode_and(va, vb))
        for a in range(m):
            b = (a - d) % m
            va = self.d12_var(a)
            vb = self.d12_var(b)
            aux_vars.append(self.encode_and(va, vb))
        return aux_vars

    def get_lambda_red_V2V2(self, d: int) -> List[int]:
        aux_vars = []
        for a in range(1, m):
            b = (a - d) % m
            if b == 0:
                continue
            va = self.d11_var(a)
            vb = self.d11_var(b)
            if va and vb:
                aux_vars.append(self.encode_not_and_not(va, vb))
        for a in range(m):
            b = (a + d) % m
            va = self.d12_var(a)
            vb = self.d12_var(b)
            aux_vars.append(self.encode_and(va, vb))
        return aux_vars

    def get_lambda_red_V1V2(self, d: int) -> List[int]:
        aux_vars = []
        for a in range(1, m):
            b = (d - a) % m
            va = self.d11_var(a)
            vb = self.d12_var(b)
            if va:
                aux_vars.append(self.encode_and(va, vb))
        for a in range(m):
            b = (a - d) % m
            if b == 0:
                aux_vars.append(self.d12_var(a))  # 0 always in D22
            else:
                va = self.d12_var(a)
                vb = self.d11_var(b)
                if vb:
                    aux_vars.append(self.encode_and_not(va, vb))
        return aux_vars

    def encode_conditional_atmost(self, cond: int, aux_vars: List[int], bound: int):
        if bound + 1 > len(aux_vars):
            return
        if bound + 1 <= 10:
            for subset in combinations(range(len(aux_vars)), bound + 1):
                clause = [-cond] + [-aux_vars[i] for i in subset]
                self.cnf.append(clause)
        else:
            atmost = CardEnc.atmost(lits=aux_vars, bound=bound, top_id=self.next_var, encoding=EncType.seqcounter)
            self.next_var = atmost.nv + 1
            for clause in atmost.clauses:
                self.cnf.append([-cond] + list(clause))

    def encode_conditional_atleast(self, cond: int, aux_vars: List[int], bound: int):
        if bound <= 0:
            return
        if bound > len(aux_vars):
            self.cnf.append([-cond])
            return
        num_false = len(aux_vars) - bound + 1
        if num_false <= 10:
            for subset in combinations(range(len(aux_vars)), num_false):
                clause = [-cond] + [aux_vars[i] for i in subset]
                self.cnf.append(clause)
        else:
            atleast = CardEnc.atleast(lits=aux_vars, bound=bound, top_id=self.next_var, encoding=EncType.seqcounter)
            self.next_var = atleast.nv + 1
            for clause in atleast.clauses:
                self.cnf.append([-cond] + list(clause))

    def encode_all_constraints(self):
        # Cardinality constraints on D11 and D12
        d11_list = list(self.d11_vars.values())
        d12_list = list(self.d12_vars.values())

        # |D11| = target means exactly target_d11_size / 2 half-reps (since symmetric)
        half_target = self.target_d11_size // 2
        print(f"Constraining |D11_half| = {half_target}")

        eq_d11 = CardEnc.equals(lits=d11_list, bound=half_target, top_id=self.next_var, encoding=EncType.seqcounter)
        self.next_var = eq_d11.nv + 1
        for clause in eq_d11.clauses:
            self.cnf.append(clause)

        print(f"Constraining |D12| = {self.target_d12_size}")
        eq_d12 = CardEnc.equals(lits=d12_list, bound=self.target_d12_size, top_id=self.next_var, encoding=EncType.seqcounter)
        self.next_var = eq_d12.nv + 1
        for clause in eq_d12.clauses:
            self.cnf.append(clause)

        # V1V1 constraints
        print("Encoding V1V1 constraints...")
        min_red_V1V1 = max(0, 2 * self.d1 - 63)
        print(f"  V1V1 blue edges need lambda_red >= {min_red_V1V1}")

        for d in range(1, half_m + 1):
            edge_var = self.d11_var(d)
            aux_vars = self.get_lambda_red_V1V1(d)

            # Red constraint: edge_var => lambda_red <= 20
            self.encode_conditional_atmost(edge_var, aux_vars, red_threshold)

            # Blue constraint: NOT edge_var => lambda_red >= min_red_V1V1
            # Equivalent: edge_var OR (lambda_red >= min_red)
            if min_red_V1V1 > 0:
                self.encode_conditional_atleast(-edge_var, aux_vars, min_red_V1V1)

        # V2V2 constraints
        print("Encoding V2V2 constraints...")
        min_red_V2V2 = max(0, 2 * self.d2 - 63)
        print(f"  V2V2 blue edges need lambda_red >= {min_red_V2V2}")

        for d in range(1, half_m + 1):
            d11_v = self.d11_var(d)
            aux_vars = self.get_lambda_red_V2V2(d)

            # d in D22 iff d not in D11, i.e., -d11_v
            # Red constraint: NOT d11_v => lambda_red <= 20
            self.encode_conditional_atmost(-d11_v, aux_vars, red_threshold)

            # Blue constraint: d11_v => lambda_red >= min_red_V2V2
            if min_red_V2V2 > 0:
                self.encode_conditional_atleast(d11_v, aux_vars, min_red_V2V2)

        # V1V2 constraints
        print("Encoding V1V2 constraints...")
        min_red_V1V2 = max(0, self.d1 + self.d2 - 63)
        print(f"  V1V2 blue edges need lambda_red >= {min_red_V1V2}")

        for d in range(m):
            edge_var = self.d12_var(d)
            aux_vars = self.get_lambda_red_V1V2(d)

            # Red constraint: edge_var => lambda_red <= 20
            self.encode_conditional_atmost(edge_var, aux_vars, red_threshold)

            # Blue constraint: NOT edge_var => lambda_red >= min_red_V1V2
            if min_red_V1V2 > 0:
                self.encode_conditional_atleast(-edge_var, aux_vars, min_red_V1V2)

        print(f"Total clauses: {len(self.cnf.clauses)}")
        print(f"Total variables: {self.next_var - 1}")

    def solve(self, timeout: int = 3600) -> Optional[Tuple[Set[int], Set[int]]]:
        print(f"\nStarting SAT solver...")
        start = time.time()

        try:
            solver = Cadical153(bootstrap_with=self.cnf.clauses)
        except:
            solver = Glucose4(bootstrap_with=self.cnf.clauses)

        result = solver.solve()
        elapsed = time.time() - start
        print(f"Solver finished in {elapsed:.2f}s")

        if result:
            model = set(solver.get_model())
            D11 = set()
            for d in range(1, half_m + 1):
                if self.d11_vars[d] in model:
                    D11.add(d)
                    D11.add(m - d)
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
        print("PySAT required")
        return

    print(f"SAT encoding for n={n}")
    print("=" * 60)

    # Try cardinalities near the cost-8 solution
    cardinalities = [
        (20, 21),  # Best heuristic solution
        (20, 20),
        (20, 22),
        (22, 21),
        (18, 21),
    ]

    for d11_size, d12_size in cardinalities:
        print(f"\n{'='*60}")
        print(f"Trying |D11|={d11_size}, |D12|={d12_size}")
        print("=" * 60)

        encoder = SATEncoderV2(d11_size, d12_size)
        encoder.encode_all_constraints()

        result = encoder.solve(timeout=300)

        if result:
            D11, D12 = result
            D22 = set(range(1, m)) - D11

            print("\nSAT!")
            print(f"D11 = {sorted(D11)}")
            print(f"D12 = {sorted(D12)}")

            from ramsey_core import BlockCirculantGraph, verify_construction, save_construction
            G = BlockCirculantGraph(n=n, D11=D11, D12=D12, D22=D22)
            v = verify_construction(G)

            print(f"Valid: {v.valid}")
            print(f"Max red: {v.max_red_common}/{v.red_threshold}")
            print(f"Max blue: {v.max_blue_common}/{v.blue_threshold}")

            if v.valid:
                save_construction(G, v, "solution_n22_sat.json", solver="SAT")
                print("SOLUTION FOUND AND SAVED!")
                return
            else:
                print(f"Violations: {len(v.violations)}")
                for viol in v.violations[:5]:
                    print(f"  {viol}")
        else:
            print("\nUNSAT for these cardinalities")

    print("\nNo valid solution found for tested cardinalities")


if __name__ == "__main__":
    main()
