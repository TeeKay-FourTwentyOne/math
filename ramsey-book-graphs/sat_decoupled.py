"""
SAT solver for R(B_21, B_22) with independent D11, D12, D22.

Key change from sat_n22_v3.py: D22 is NOT derived as complement(D11).
Instead, D22 has its own variables, allowing asymmetric degree configurations.

Variables:
  - p_i (i=1..21): pair i in D11 (symmetric: d and m-d both in/out)
  - y_d (d=0..42): d in D12
  - q_i (i=1..21): pair i in D22 (symmetric)
  Total: 21 + 43 + 21 = 85 primary variables

For each (|D11|, |D12|, |D22|) triple, degrees are fixed:
  d1 = |D11| + |D12|, d2 = |D22| + |D12|
This makes blue thresholds constants, enabling clean SAT encoding.
"""

import sys
import os
import time
from typing import Set, List, Tuple, Optional
from itertools import combinations

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    from pysat.solvers import Cadical153, Glucose4
    from pysat.card import CardEnc, EncType
    from pysat.formula import CNF
    PYSAT_AVAILABLE = True
except ImportError:
    PYSAT_AVAILABLE = False
    print("PySAT not available. Install: pip install python-sat")

from ramsey_core import BlockCirculantGraph, verify_construction

n = 22
m = 43
N = 86
half_m = m // 2  # 21
red_threshold = 20  # n - 2
blue_threshold = 21  # n - 1


class DecoupledSATEncoder:
    """SAT encoder with independent D11, D12, D22."""

    def __init__(self, d11_size: int, d12_size: int, d22_size: int):
        self.cnf = CNF()
        self.next_var = 1
        self.d11_size = d11_size
        self.d12_size = d12_size
        self.d22_size = d22_size
        self.d1 = d11_size + d12_size  # degree of V1 vertices
        self.d2 = d22_size + d12_size  # degree of V2 vertices

        # Allocate variables
        # D11: pair variables (1..half_m) - symmetric, so one var per pair
        self.d11_vars = {d: self._alloc() for d in range(1, half_m + 1)}
        # D12: individual variables (0..m-1)
        self.d12_vars = {d: self._alloc() for d in range(m)}
        # D22: pair variables (1..half_m) - symmetric
        self.d22_vars = {d: self._alloc() for d in range(1, half_m + 1)}

        # Blue thresholds (lambda_red upper bound for blue edges)
        self.blue_V1V1 = 2 * self.d1 - (N - 2) + blue_threshold
        self.blue_V2V2 = 2 * self.d2 - (N - 2) + blue_threshold
        self.blue_V1V2 = self.d1 + self.d2 - (N - 2) + blue_threshold

        print(f"|D11|={d11_size}, |D12|={d12_size}, |D22|={d22_size}")
        print(f"d1={self.d1}, d2={self.d2}")
        print(f"Red: lambda_red <= {red_threshold}")
        print(f"Blue V1V1: lambda_red <= {self.blue_V1V1}")
        print(f"Blue V2V2: lambda_red <= {self.blue_V2V2}")
        print(f"Blue V1V2: lambda_red <= {self.blue_V1V2}")

        # Check feasibility
        if self.blue_V1V1 < 0 or self.blue_V2V2 < 0 or self.blue_V1V2 < 0:
            print("  -> INFEASIBLE (negative blue threshold)")
            self.feasible = False
        else:
            self.feasible = True

    def _alloc(self) -> int:
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

    def d22_var(self, d: int) -> Optional[int]:
        d = d % m
        if d == 0:
            return None
        if d > half_m:
            d = m - d
        return self.d22_vars[d]

    def encode_and(self, a: int, b: int) -> int:
        """aux <=> a AND b"""
        if a == b:
            return a
        aux = self._alloc()
        self.cnf.append([-aux, a])
        self.cnf.append([-aux, b])
        self.cnf.append([aux, -a, -b])
        return aux

    def get_aux_V1V1(self, d: int) -> List[int]:
        """Auxiliary vars for common neighbors of V1-V1 edge at difference d.
        lambda_red = Delta(D11,D11,d) + Delta(D12,D12,d)"""
        aux = []
        # Delta(D11, D11, d): pairs (a, a-d) both in D11
        for a in range(1, m):
            b = (a - d) % m
            if b == 0:
                continue
            va = self.d11_var(a)
            vb = self.d11_var(b)
            if va and vb:
                aux.append(self.encode_and(va, vb))
        # Delta(D12, D12, d): pairs (a, a-d) both in D12
        for a in range(m):
            b = (a - d) % m
            aux.append(self.encode_and(self.d12_var(a), self.d12_var(b)))
        return aux

    def get_aux_V2V2(self, d: int) -> List[int]:
        """Auxiliary vars for common neighbors of V2-V2 edge at difference d.
        lambda_red = Delta(D22,D22,d) + Delta(D12^T,D12^T,d)"""
        aux = []
        # Delta(D22, D22, d)
        for a in range(1, m):
            b = (a - d) % m
            if b == 0:
                continue
            va = self.d22_var(a)
            vb = self.d22_var(b)
            if va and vb:
                aux.append(self.encode_and(va, vb))
        # Delta(D12^T, D12^T, d) where D12^T = {-x : x in D12}
        # pairs (-a, -b) with (-a)-(-b) = b-a = -d+2b-2b... no.
        # D12^T diff d: need -a - (-b) = b - a = d, so a = b - d
        # Equivalently: a such that (-a) and (-(a+d)) both in D12,
        # i.e., D12 var at (-a)%m and D12 var at (-(a+d))%m
        # Simpler: Delta(D12^T, D12^T, d) counts x in D12^T where (x-d) in D12^T
        # x in D12^T means -x in D12. So count -x in D12 where -(x-d) in D12
        # Let a = -x, then a in D12 and a+d in D12
        for a in range(m):
            b = (a + d) % m
            aux.append(self.encode_and(self.d12_var(a), self.d12_var(b)))
        return aux

    def get_aux_V1V2(self, d: int) -> List[int]:
        """Auxiliary vars for common neighbors of V1-V2 edge at difference d.
        lambda_red = Sigma(D11,D12,d) + Delta(D12,D22,d)"""
        aux = []
        # Sigma(D11, D12, d): count a in D11 where (d-a) in D12
        for a in range(1, m):
            b = (d - a) % m
            va = self.d11_var(a)
            if va:
                aux.append(self.encode_and(va, self.d12_var(b)))
        # Delta(D12, D22, d): count a in D12 where (a-d) in D22
        for a in range(m):
            b = (a - d) % m
            if b == 0:
                continue
            vb = self.d22_var(b)
            if vb:
                aux.append(self.encode_and(self.d12_var(a), vb))
        return aux

    def add_atmost(self, cond: int, aux: List[int], bound: int):
        """Add conditional constraint: cond => sum(aux) <= bound"""
        if bound < 0:
            self.cnf.append([-cond])
            return
        if bound >= len(aux):
            return
        if bound + 1 <= 8:
            for subset in combinations(range(len(aux)), bound + 1):
                self.cnf.append([-cond] + [-aux[i] for i in subset])
        else:
            enc = CardEnc.atmost(lits=aux, bound=bound,
                                 top_id=self.next_var,
                                 encoding=EncType.seqcounter)
            self.next_var = enc.nv + 1
            for c in enc.clauses:
                self.cnf.append([-cond] + list(c))

    def encode(self):
        if not self.feasible:
            return

        # Cardinality constraints on set sizes
        d11_list = list(self.d11_vars.values())
        d12_list = list(self.d12_vars.values())
        d22_list = list(self.d22_vars.values())

        half_d11 = self.d11_size // 2
        half_d22 = self.d22_size // 2

        print(f"Constraining |D11_half|={half_d11}, |D12|={self.d12_size}, "
              f"|D22_half|={half_d22}")

        enc = CardEnc.equals(lits=d11_list, bound=half_d11,
                             top_id=self.next_var, encoding=EncType.seqcounter)
        self.next_var = enc.nv + 1
        self.cnf.extend(enc.clauses)

        enc = CardEnc.equals(lits=d12_list, bound=self.d12_size,
                             top_id=self.next_var, encoding=EncType.seqcounter)
        self.next_var = enc.nv + 1
        self.cnf.extend(enc.clauses)

        enc = CardEnc.equals(lits=d22_list, bound=half_d22,
                             top_id=self.next_var, encoding=EncType.seqcounter)
        self.next_var = enc.nv + 1
        self.cnf.extend(enc.clauses)

        # V1V1 constraints (d = 1 to half_m, symmetric)
        print("Encoding V1V1...")
        for d in range(1, half_m + 1):
            ev = self.d11_var(d)
            aux = self.get_aux_V1V1(d)
            # Red: d in D11 => lambda_red <= 20
            self.add_atmost(ev, aux, red_threshold)
            # Blue: d not in D11 => lambda_red <= blue_V1V1
            self.add_atmost(-ev, aux, self.blue_V1V1)

        # V2V2 constraints
        print("Encoding V2V2...")
        for d in range(1, half_m + 1):
            ev = self.d22_var(d)
            aux = self.get_aux_V2V2(d)
            # Red: d in D22 => lambda_red <= 20
            self.add_atmost(ev, aux, red_threshold)
            # Blue: d not in D22 => lambda_red <= blue_V2V2
            self.add_atmost(-ev, aux, self.blue_V2V2)

        # V1V2 constraints
        print("Encoding V1V2...")
        for d in range(m):
            ev = self.d12_var(d)
            aux = self.get_aux_V1V2(d)
            # Red: d in D12 => lambda_red <= 20
            self.add_atmost(ev, aux, red_threshold)
            # Blue: d not in D12 => lambda_red <= blue_V1V2
            self.add_atmost(-ev, aux, self.blue_V1V2)

        print(f"Clauses: {len(self.cnf.clauses)}, Variables: {self.next_var - 1}")

    def solve(self) -> Optional[Tuple[Set[int], Set[int], Set[int]]]:
        if not self.feasible:
            print("Skipping (infeasible)")
            return None

        print("Solving...")
        start = time.time()
        try:
            solver = Cadical153(bootstrap_with=self.cnf.clauses)
        except Exception:
            solver = Glucose4(bootstrap_with=self.cnf.clauses)

        result = solver.solve()
        elapsed = time.time() - start
        print(f"Time: {elapsed:.2f}s")

        if result:
            model = set(solver.get_model())
            D11 = set()
            for d in range(1, half_m + 1):
                if self.d11_vars[d] in model:
                    D11.add(d)
                    D11.add(m - d)
            D12 = {d for d in range(m) if self.d12_vars[d] in model}
            D22 = set()
            for d in range(1, half_m + 1):
                if self.d22_vars[d] in model:
                    D22.add(d)
                    D22.add(m - d)
            solver.delete()
            return D11, D12, D22
        solver.delete()
        return None


def main():
    if not PYSAT_AVAILABLE:
        print("Need PySAT: pip install python-sat")
        return

    from ramsey_core import BlockCirculantGraph, verify_construction, save_construction

    # Degree configurations to try
    # d1 = |D11| + |D12|, d2 = |D22| + |D12|
    # Need: d1, d2 >= 32 (blue feasibility), |D11|,|D22| even
    #
    # Previous best with complement: d1=43, d2=41 (|D11|=22, |D12|=21, |D22|=20)
    # Density trap suggestion: try asymmetric, d1 > d2
    #
    # Strategy: sweep around the sweet spot

    configs = []

    # Near-symmetric configs (close to complement baseline)
    for nd12 in [20, 21, 22]:
        for nd11 in [18, 20, 22, 24]:
            for nd22 in [18, 20, 22, 24]:
                d1 = nd11 + nd12
                d2 = nd22 + nd12
                # Blue feasibility
                if 2 * d1 - 84 + blue_threshold < 0:
                    continue
                if 2 * d2 - 84 + blue_threshold < 0:
                    continue
                if d1 + d2 - 84 + blue_threshold < 0:
                    continue
                # Skip extreme configs
                if d1 < 38 or d1 > 48 or d2 < 38 or d2 > 48:
                    continue
                configs.append((nd11, nd12, nd22))

    # Remove duplicates and sort by distance from (20, 21, 20)
    configs = sorted(set(configs),
                     key=lambda c: abs(c[0]-20) + abs(c[1]-21) + abs(c[2]-20))

    print(f"Testing {len(configs)} configurations")
    print(f"{'|D11|':>5} {'|D12|':>5} {'|D22|':>5} {'d1':>4} {'d2':>4} "
          f"{'blue11':>6} {'blue22':>6} {'blue12':>6} {'Result':>10}")
    print("-" * 65)

    for nd11, nd12, nd22 in configs:
        d1 = nd11 + nd12
        d2 = nd22 + nd12
        b11 = 2 * d1 - 84 + blue_threshold
        b22 = 2 * d2 - 84 + blue_threshold
        b12 = d1 + d2 - 84 + blue_threshold

        print(f"\n{'='*60}")
        print(f"|D11|={nd11}, |D12|={nd12}, |D22|={nd22} "
              f"(d1={d1}, d2={d2}, bV1V1={b11}, bV2V2={b22}, bV1V2={b12})")
        print("=" * 60)

        enc = DecoupledSATEncoder(nd11, nd12, nd22)
        if not enc.feasible:
            print(f"{'SKIP':>10} (infeasible)")
            continue

        enc.encode()
        result = enc.solve()

        if result:
            D11, D12, D22 = result
            G = BlockCirculantGraph(n=n, D11=D11, D12=D12, D22=D22)
            v = verify_construction(G)
            status = "VALID" if v.valid else f"FAIL({len(v.violations)}v)"
            print(f"SAT! {status}")
            print(f"D11={sorted(D11)}")
            print(f"D12={sorted(D12)}")
            print(f"D22={sorted(D22)}")
            print(f"Red max={v.max_red_common}/{v.red_threshold}, "
                  f"Blue max={v.max_blue_common}/{v.blue_threshold}")

            if v.valid:
                save_construction(G, v, "solution_n22_decoupled.json",
                                  "SAT-decoupled")
                print("\n*** VALID SOLUTION FOUND! ***")
                return
            else:
                print(f"Violations ({len(v.violations)}):")
                for vt, d, ex in v.violations[:10]:
                    print(f"  {vt} d={d} excess={ex}")
        else:
            print("UNSAT")

    print("\nNo valid solution found across all configurations")


if __name__ == "__main__":
    main()
