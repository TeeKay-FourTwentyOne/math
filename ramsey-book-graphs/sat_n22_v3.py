"""
Correct SAT encoding for n=22.

Blue constraint: lambda_blue <= 21
lambda_blue = (N-2) - deg_u - deg_v + lambda_red
So: lambda_red <= deg_u + deg_v - 63

For V1V1: lambda_red <= 2*d1 - 63 (blue edges)
For V2V2: lambda_red <= 2*d2 - 63 (blue edges)
For V1V2: lambda_red <= d1 + d2 - 63 (blue edges)

Red constraint: lambda_red <= 20 (for red edges)
"""

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

n = 22
m = 43
N = 86
half_m = m // 2
red_threshold = 20
blue_threshold = 21


class SATEncoderV3:
    def __init__(self, d11_size: int, d12_size: int):
        self.cnf = CNF()
        self.next_var = 1
        self.d11_size = d11_size
        self.d12_size = d12_size
        self.d1 = d11_size + d12_size
        self.d2 = (42 - d11_size) + d12_size

        self.d11_vars = {d: self.alloc_var() for d in range(1, half_m + 1)}
        self.d12_vars = {d: self.alloc_var() for d in range(m)}

        # Blue thresholds (max lambda_red for blue edges)
        self.blue_V1V1 = 2 * self.d1 - 63
        self.blue_V2V2 = 2 * self.d2 - 63
        self.blue_V1V2 = self.d1 + self.d2 - 63

        print(f"|D11|={d11_size}, |D12|={d12_size}, d1={self.d1}, d2={self.d2}")
        print(f"Red edges: lambda_red <= {red_threshold}")
        print(f"Blue V1V1: lambda_red <= {self.blue_V1V1}")
        print(f"Blue V2V2: lambda_red <= {self.blue_V2V2}")
        print(f"Blue V1V2: lambda_red <= {self.blue_V1V2}")

    def alloc_var(self) -> int:
        v = self.next_var
        self.next_var += 1
        return v

    def d11_var(self, d: int) -> Optional[int]:
        d = d % m
        if d == 0:
            return None
        if d > half_m:
            d = m - d
        return self.d11_vars[d]

    def d12_var(self, d: int) -> int:
        return self.d12_vars[d % m]

    def encode_and(self, a: int, b: int) -> int:
        if a == b:
            return a
        aux = self.alloc_var()
        self.cnf.append([-aux, a])
        self.cnf.append([-aux, b])
        self.cnf.append([aux, -a, -b])
        return aux

    def encode_not_and_not(self, a: int, b: int) -> int:
        aux = self.alloc_var()
        self.cnf.append([-aux, -a])
        self.cnf.append([-aux, -b])
        self.cnf.append([aux, a, b])
        return aux

    def encode_and_not(self, a: int, b: int) -> int:
        aux = self.alloc_var()
        self.cnf.append([-aux, a])
        self.cnf.append([-aux, -b])
        self.cnf.append([aux, -a, b])
        return aux

    def get_aux_V1V1(self, d: int) -> List[int]:
        aux = []
        for a in range(1, m):
            b = (a - d) % m
            if b == 0:
                continue
            va, vb = self.d11_var(a), self.d11_var(b)
            if va and vb:
                aux.append(self.encode_and(va, vb))
        for a in range(m):
            b = (a - d) % m
            aux.append(self.encode_and(self.d12_var(a), self.d12_var(b)))
        return aux

    def get_aux_V2V2(self, d: int) -> List[int]:
        aux = []
        for a in range(1, m):
            b = (a - d) % m
            if b == 0:
                continue
            va, vb = self.d11_var(a), self.d11_var(b)
            if va and vb:
                aux.append(self.encode_not_and_not(va, vb))
        for a in range(m):
            b = (a + d) % m
            aux.append(self.encode_and(self.d12_var(a), self.d12_var(b)))
        return aux

    def get_aux_V1V2(self, d: int) -> List[int]:
        aux = []
        for a in range(1, m):
            b = (d - a) % m
            va = self.d11_var(a)
            if va:
                aux.append(self.encode_and(va, self.d12_var(b)))
        for a in range(m):
            b = (a - d) % m
            if b == 0:
                aux.append(self.d12_var(a))
            else:
                vb = self.d11_var(b)
                if vb:
                    aux.append(self.encode_and_not(self.d12_var(a), vb))
        return aux

    def add_atmost(self, cond: int, aux: List[int], bound: int):
        """Add: cond => sum(aux) <= bound"""
        if bound < 0:
            self.cnf.append([-cond])
            return
        if bound >= len(aux):
            return
        if bound + 1 <= 10:
            for subset in combinations(range(len(aux)), bound + 1):
                self.cnf.append([-cond] + [-aux[i] for i in subset])
        else:
            enc = CardEnc.atmost(lits=aux, bound=bound, top_id=self.next_var, encoding=EncType.seqcounter)
            self.next_var = enc.nv + 1
            for c in enc.clauses:
                self.cnf.append([-cond] + list(c))

    def encode(self):
        # Cardinality constraints
        d11_list = list(self.d11_vars.values())
        d12_list = list(self.d12_vars.values())

        half_d11 = self.d11_size // 2
        print(f"Constraining |D11_half|={half_d11}, |D12|={self.d12_size}")

        enc = CardEnc.equals(lits=d11_list, bound=half_d11, top_id=self.next_var, encoding=EncType.seqcounter)
        self.next_var = enc.nv + 1
        self.cnf.extend(enc.clauses)

        enc = CardEnc.equals(lits=d12_list, bound=self.d12_size, top_id=self.next_var, encoding=EncType.seqcounter)
        self.next_var = enc.nv + 1
        self.cnf.extend(enc.clauses)

        # V1V1 constraints
        print("Encoding V1V1...")
        for d in range(1, half_m + 1):
            ev = self.d11_var(d)
            aux = self.get_aux_V1V1(d)
            # Red: ev => lambda <= 20
            self.add_atmost(ev, aux, red_threshold)
            # Blue: NOT ev => lambda <= blue_V1V1
            self.add_atmost(-ev, aux, self.blue_V1V1)

        # V2V2 constraints
        print("Encoding V2V2...")
        for d in range(1, half_m + 1):
            d11v = self.d11_var(d)
            aux = self.get_aux_V2V2(d)
            # Red (d in D22 = d not in D11): NOT d11v => lambda <= 20
            self.add_atmost(-d11v, aux, red_threshold)
            # Blue (d not in D22 = d in D11): d11v => lambda <= blue_V2V2
            self.add_atmost(d11v, aux, self.blue_V2V2)

        # V1V2 constraints
        print("Encoding V1V2...")
        for d in range(m):
            ev = self.d12_var(d)
            aux = self.get_aux_V1V2(d)
            # Red: ev => lambda <= 20
            self.add_atmost(ev, aux, red_threshold)
            # Blue: NOT ev => lambda <= blue_V1V2
            self.add_atmost(-ev, aux, self.blue_V1V2)

        print(f"Clauses: {len(self.cnf.clauses)}, Variables: {self.next_var - 1}")

    def solve(self) -> Optional[Tuple[Set[int], Set[int]]]:
        print("Solving...")
        start = time.time()
        try:
            solver = Cadical153(bootstrap_with=self.cnf.clauses)
        except:
            solver = Glucose4(bootstrap_with=self.cnf.clauses)

        result = solver.solve()
        print(f"Time: {time.time() - start:.2f}s")

        if result:
            model = set(solver.get_model())
            D11 = set()
            for d in range(1, half_m + 1):
                if self.d11_vars[d] in model:
                    D11.add(d)
                    D11.add(m - d)
            D12 = {d for d in range(m) if self.d12_vars[d] in model}
            solver.delete()
            return D11, D12
        solver.delete()
        return None


def main():
    if not PYSAT_AVAILABLE:
        print("Need PySAT")
        return

    from ramsey_core import BlockCirculantGraph, verify_construction, save_construction

    # Test multiple cardinalities
    configs = [(20, 21), (20, 20), (20, 22), (18, 21), (22, 21), (22, 20), (22, 22)]

    for d11, d12 in configs:
        print(f"\n{'='*60}")
        print(f"Testing |D11|={d11}, |D12|={d12}")
        print("=" * 60)

        enc = SATEncoderV3(d11, d12)
        enc.encode()
        result = enc.solve()

        if result:
            D11, D12 = result
            D22 = set(range(1, m)) - D11

            print(f"SAT! D11={sorted(D11)}")
            print(f"D12={sorted(D12)}")

            G = BlockCirculantGraph(n=n, D11=D11, D12=D12, D22=D22)
            v = verify_construction(G)

            print(f"Valid: {v.valid}, Red: {v.max_red_common}/{v.red_threshold}, Blue: {v.max_blue_common}/{v.blue_threshold}")

            if v.valid:
                save_construction(G, v, "solution_n22.json", "SAT")
                print("VALID SOLUTION FOUND!")
                return
            else:
                print(f"Violations: {len(v.violations)}")
                for viol in v.violations[:5]:
                    print(f"  {viol}")
        else:
            print("UNSAT")

    print("\nNo valid solution found")


if __name__ == "__main__":
    main()
