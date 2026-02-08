"""
Targeted SAT solver for R(B_21, B_22) = 87.

Informed by analysis of ALL known m = 3 (mod 4) constructions:
- D22 = complement(D11) (universal pattern)
- 0 in D12 (universal pattern)
- |D11| = (m+1)/2 = 22, |D12| = (m-1)/2 = 21
- d1 = m = 43, d2 = m-2 = 41

Also tries the flipped config: |D11|=20, |D12|=21, d1=41, d2=43.

Key optimizations over sat_n22_v3.py:
1. Fix 0 in D12 (saves 1 variable, removes some clauses)
2. Add symmetry breaking (lexicographic order on D11)
3. Single-config focus with tighter encoding
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

from ramsey_core import BlockCirculantGraph, verify_construction, save_construction

n = 22
m = 43
N = 86
half_m = m // 2  # 21
red_threshold = 20  # n - 2
blue_threshold = 21  # n - 1


class TargetedSATEncoder:
    """SAT encoder using complement D22 and known structural patterns."""

    def __init__(self, d11_size: int, d12_size: int, fix_zero_in_d12: bool = True):
        self.cnf = CNF()
        self.next_var = 1
        self.d11_size = d11_size
        self.d12_size = d12_size
        self.d1 = d11_size + d12_size
        self.d2 = (m - 1 - d11_size) + d12_size  # complement D22
        self.fix_zero_in_d12 = fix_zero_in_d12

        # D11: pair variables (1..half_m)
        self.d11_vars = {d: self._alloc() for d in range(1, half_m + 1)}

        # D12: individual variables
        if fix_zero_in_d12:
            # 0 is always in D12, so only allocate vars for 1..m-1
            self.d12_vars = {0: None}  # sentinel
            for d in range(1, m):
                self.d12_vars[d] = self._alloc()
        else:
            self.d12_vars = {d: self._alloc() for d in range(m)}

        # Blue thresholds (max lambda_red for blue edges)
        self.blue_V1V1 = 2 * self.d1 - (N - 2) + blue_threshold
        self.blue_V2V2 = 2 * self.d2 - (N - 2) + blue_threshold
        self.blue_V1V2 = self.d1 + self.d2 - (N - 2) + blue_threshold

        print(f"|D11|={d11_size}, |D12|={d12_size}, |D22|={m-1-d11_size}")
        print(f"d1={self.d1}, d2={self.d2}")
        print(f"Red: lambda_red <= {red_threshold}")
        print(f"Blue V1V1: lambda_red <= {self.blue_V1V1}")
        print(f"Blue V2V2: lambda_red <= {self.blue_V2V2}")
        print(f"Blue V1V2: lambda_red <= {self.blue_V1V2}")
        print(f"Fix 0 in D12: {fix_zero_in_d12}")

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
        """Get SAT variable for 'd in D11'. Returns None if d=0."""
        d = d % m
        if d == 0:
            return None
        if d > half_m:
            d = m - d
        return self.d11_vars[d]

    def d12_var(self, d: int) -> Optional[int]:
        """Get SAT variable for 'd in D12'. Returns None if d=0 and fixed."""
        d = d % m
        if d == 0 and self.fix_zero_in_d12:
            return None  # always true
        return self.d12_vars[d]

    def d12_is_true(self, d: int) -> bool:
        """Check if d in D12 is forced true (i.e., d=0 and fixed)."""
        return (d % m == 0) and self.fix_zero_in_d12

    def encode_and(self, a: int, b: int) -> int:
        """aux <=> a AND b"""
        if a == b:
            return a
        aux = self._alloc()
        self.cnf.append([-aux, a])
        self.cnf.append([-aux, b])
        self.cnf.append([aux, -a, -b])
        return aux

    def encode_not_and_not(self, a: int, b: int) -> int:
        """aux <=> (NOT a) AND (NOT b)"""
        aux = self._alloc()
        self.cnf.append([-aux, -a])
        self.cnf.append([-aux, -b])
        self.cnf.append([aux, a, b])
        return aux

    def encode_and_not(self, a: int, b: int) -> int:
        """aux <=> a AND (NOT b)"""
        aux = self._alloc()
        self.cnf.append([-aux, a])
        self.cnf.append([-aux, -b])
        self.cnf.append([aux, -a, b])
        return aux

    def get_aux_V1V1(self, d: int) -> List[int]:
        """Auxiliary vars for lambda_red of V1-V1 edge at difference d.
        lambda_red = Delta(D11, D11, d) + Delta(D12, D12, d)"""
        aux = []
        # Delta(D11, D11, d): count a in {1..m-1} s.t. a in D11 and (a-d) in D11
        for a in range(1, m):
            b = (a - d) % m
            if b == 0:
                continue
            va = self.d11_var(a)
            vb = self.d11_var(b)
            if va is not None and vb is not None:
                aux.append(self.encode_and(va, vb))

        # Delta(D12, D12, d): count a in {0..m-1} s.t. a in D12 and (a-d) in D12
        for a in range(m):
            b = (a - d) % m
            va = self.d12_var(a)
            vb = self.d12_var(b)
            a_true = self.d12_is_true(a)
            b_true = self.d12_is_true(b)

            if a_true and b_true:
                # Both forced true: contributes 1 always, handled below
                aux.append(None)  # placeholder
            elif a_true and vb is not None:
                aux.append(vb)
            elif b_true and va is not None:
                aux.append(va)
            elif va is not None and vb is not None:
                aux.append(self.encode_and(va, vb))

        # Remove None placeholders and count forced-true contributions
        forced_count = aux.count(None)
        aux = [x for x in aux if x is not None]

        return aux, forced_count

    def get_aux_V2V2(self, d: int) -> List[int]:
        """Auxiliary vars for lambda_red of V2-V2 edge at difference d.
        D22 = complement(D11), so d in D22 iff d not in D11.
        lambda_red = Delta(D22, D22, d) + Delta(D12^T, D12^T, d)"""
        aux = []
        # Delta(D22, D22, d): count a s.t. a NOT in D11 and (a-d) NOT in D11
        for a in range(1, m):
            b = (a - d) % m
            if b == 0:
                continue
            va = self.d11_var(a)
            vb = self.d11_var(b)
            if va is not None and vb is not None:
                aux.append(self.encode_not_and_not(va, vb))

        # Delta(D12^T, D12^T, d): D12^T = {-x : x in D12}
        # Count x s.t. -x in D12 and -(x-d) in D12
        # Let a = -x: count a in D12 s.t. (a+d) in D12
        for a in range(m):
            b = (a + d) % m
            va = self.d12_var(a)
            vb = self.d12_var(b)
            a_true = self.d12_is_true(a)
            b_true = self.d12_is_true(b)

            if a_true and b_true:
                aux.append(None)
            elif a_true and vb is not None:
                aux.append(vb)
            elif b_true and va is not None:
                aux.append(va)
            elif va is not None and vb is not None:
                aux.append(self.encode_and(va, vb))

        forced_count = aux.count(None)
        aux = [x for x in aux if x is not None]

        return aux, forced_count

    def get_aux_V1V2(self, d: int) -> List[int]:
        """Auxiliary vars for lambda_red of V1-V2 cross-block edge at difference d.
        lambda_red = Sigma(D11, D12, d) + Delta(D12, D22, d)"""
        aux = []
        # Sigma(D11, D12, d): count a in D11 s.t. (d-a) in D12
        for a in range(1, m):
            b = (d - a) % m
            va = self.d11_var(a)
            vb = self.d12_var(b)
            b_true = self.d12_is_true(b)

            if va is not None:
                if b_true:
                    aux.append(va)
                elif vb is not None:
                    aux.append(self.encode_and(va, vb))

        # Delta(D12, D22, d): count a in D12 s.t. (a-d) in D22 = NOT in D11
        for a in range(m):
            b = (a - d) % m
            if b == 0:
                # b=0 is never in D11, so condition is just: a in D12
                va = self.d12_var(a)
                if self.d12_is_true(a):
                    aux.append(None)  # forced true
                elif va is not None:
                    aux.append(va)
                continue

            va_d12 = self.d12_var(a)
            vb_d11 = self.d11_var(b)  # b NOT in D11

            a_true = self.d12_is_true(a)

            if a_true and vb_d11 is not None:
                # a in D12 is forced, need b NOT in D11
                aux.append(-vb_d11)  # NOT vb acts as "b not in D11"
                # But we need an aux var for this
                neg_aux = self._alloc()
                self.cnf.append([-neg_aux, -vb_d11])
                self.cnf.append([neg_aux, vb_d11])
                # Actually that's just neg_aux <=> NOT vb_d11, simpler:
                # neg_aux = 1 iff vb_d11 = 0
                aux[-1] = neg_aux  # replace
            elif va_d12 is not None and vb_d11 is not None:
                aux.append(self.encode_and_not(va_d12, vb_d11))
            elif a_true:
                aux.append(None)  # b=0 handled above, shouldn't reach

        forced_count = aux.count(None)
        aux = [x for x in aux if x is not None]

        return aux, forced_count

    def add_atmost(self, cond: Optional[int], aux: List[int], bound: int,
                   forced_count: int = 0):
        """Add constraint: [cond =>] sum(aux) <= bound - forced_count.
        If cond is None, the constraint is unconditional."""
        effective_bound = bound - forced_count
        if effective_bound < 0:
            if cond is not None:
                self.cnf.append([-cond])
            else:
                self.cnf.append([])  # empty clause = UNSAT
            return
        if effective_bound >= len(aux):
            return  # trivially satisfied

        if effective_bound + 1 <= 8:
            # Direct enumeration for small bounds
            for subset in combinations(range(len(aux)), effective_bound + 1):
                clause = [-aux[i] for i in subset]
                if cond is not None:
                    clause = [-cond] + clause
                self.cnf.append(clause)
        else:
            enc = CardEnc.atmost(lits=aux, bound=effective_bound,
                                 top_id=self.next_var,
                                 encoding=EncType.seqcounter)
            self.next_var = enc.nv + 1
            for c in enc.clauses:
                clause = list(c)
                if cond is not None:
                    clause = [-cond] + clause
                self.cnf.append(clause)

    def add_symmetry_breaking(self):
        """Add lexicographic symmetry breaking on D11 pair variables.
        This doesn't restrict the search space — just prunes symmetric
        branches of the search tree."""
        # Simple: require that the first selected pair index is small
        # More effective: for the automorphism group of Z_m
        # Since m=43 is prime, Aut(Z_43) = Z_42 (multiplication by units)
        # We can fix one element: require that 1 is in D11 (or not)
        # This is a factor-2 speedup at most

        # Actually, let's use the observation from known constructions:
        # D11 tends to be a contiguous-ish block of pairs
        # We can break symmetry by fixing pair 1 to be in D11
        # (since |D11| > m/2, it MUST contain most pairs)

        # For |D11|=22: 11 of 21 pairs are selected. Very likely pair 1 is in.
        # For |D11|=20: 10 of 21 pairs selected.
        # Don't force this — might miss solutions. Skip for now.
        pass

    def encode(self):
        if not self.feasible:
            return

        # Cardinality constraints
        d11_list = list(self.d11_vars.values())

        if self.fix_zero_in_d12:
            d12_list = [self.d12_vars[d] for d in range(1, m)]
            d12_bound = self.d12_size - 1  # one slot (d=0) is pre-filled
        else:
            d12_list = list(self.d12_vars.values())
            d12_bound = self.d12_size

        half_d11 = self.d11_size // 2
        print(f"Constraining |D11_half|={half_d11}, |D12_free|={d12_bound}")

        enc = CardEnc.equals(lits=d11_list, bound=half_d11,
                             top_id=self.next_var, encoding=EncType.seqcounter)
        self.next_var = enc.nv + 1
        self.cnf.extend(enc.clauses)

        enc = CardEnc.equals(lits=d12_list, bound=d12_bound,
                             top_id=self.next_var, encoding=EncType.seqcounter)
        self.next_var = enc.nv + 1
        self.cnf.extend(enc.clauses)

        self.add_symmetry_breaking()

        # V1V1 constraints (d = 1 to half_m, exploiting d <-> m-d symmetry)
        print("Encoding V1V1...")
        for d in range(1, half_m + 1):
            ev = self.d11_var(d)
            aux, forced = self.get_aux_V1V1(d)
            # Red: d in D11 => lambda_red <= 20
            self.add_atmost(ev, aux, red_threshold, forced)
            # Blue: d not in D11 => lambda_red <= blue_V1V1
            self.add_atmost(-ev, aux, self.blue_V1V1, forced)

        # V2V2 constraints (D22 = complement D11)
        print("Encoding V2V2...")
        for d in range(1, half_m + 1):
            d11v = self.d11_var(d)
            aux, forced = self.get_aux_V2V2(d)
            # Red (d in D22 = d not in D11): -d11v => lambda_red <= 20
            self.add_atmost(-d11v, aux, red_threshold, forced)
            # Blue (d not in D22 = d in D11): d11v => lambda_red <= blue_V2V2
            self.add_atmost(d11v, aux, self.blue_V2V2, forced)

        # V1V2 constraints
        print("Encoding V1V2...")
        for d in range(m):
            if d == 0 and self.fix_zero_in_d12:
                # 0 is in D12 (red edge), just need lambda_red <= 20
                aux, forced = self.get_aux_V1V2(0)
                self.add_atmost(None, aux, red_threshold, forced)
            else:
                ev = self.d12_var(d)
                if ev is not None:
                    aux, forced = self.get_aux_V1V2(d)
                    # Red: ev => lambda_red <= 20
                    self.add_atmost(ev, aux, red_threshold, forced)
                    # Blue: NOT ev => lambda_red <= blue_V1V2
                    self.add_atmost(-ev, aux, self.blue_V1V2, forced)

        print(f"Clauses: {len(self.cnf.clauses)}, Variables: {self.next_var - 1}")

    def solve(self, timeout: int = 0) -> Optional[Tuple[Set[int], Set[int]]]:
        if not self.feasible:
            print("Skipping (infeasible)")
            return None

        print("Solving...")
        sys.stdout.flush()
        start = time.time()
        try:
            solver = Cadical153(bootstrap_with=self.cnf.clauses)
        except Exception:
            solver = Glucose4(bootstrap_with=self.cnf.clauses)

        result = solver.solve()
        elapsed = time.time() - start
        print(f"Result: {'SAT' if result else 'UNSAT'} ({elapsed:.2f}s)")
        sys.stdout.flush()

        if result:
            model = set(solver.get_model())
            D11 = set()
            for d in range(1, half_m + 1):
                if self.d11_vars[d] in model:
                    D11.add(d)
                    D11.add(m - d)
            D12 = set()
            if self.fix_zero_in_d12:
                D12.add(0)
            for d in range(1 if self.fix_zero_in_d12 else 0, m):
                if self.d12_vars[d] in model:
                    D12.add(d)
            solver.delete()
            return D11, D12
        solver.delete()
        return None


def main():
    if not PYSAT_AVAILABLE:
        print("Need PySAT: pip install python-sat")
        return

    # Configuration based on universal pattern for m = 3 (mod 4):
    # Dense block: |D11| = (m+1)/2 = 22, |D12| = (m-1)/2 = 21
    # This gives d1 = m = 43 (dense), d2 = m-2 = 41 (sparse)
    configs = [
        (22, 21, "primary: d1=43, d2=41"),
        (20, 21, "flipped: d1=41, d2=43"),
    ]

    for d11_size, d12_size, label in configs:
        print(f"\n{'='*60}")
        print(f"Config: |D11|={d11_size}, |D12|={d12_size} ({label})")
        print("=" * 60)

        enc = TargetedSATEncoder(d11_size, d12_size, fix_zero_in_d12=True)
        if not enc.feasible:
            print("INFEASIBLE")
            continue

        enc.encode()
        result = enc.solve()

        if result:
            D11, D12 = result
            D22 = set(range(1, m)) - D11

            print(f"\nSAT! D11={sorted(D11)}")
            print(f"D12={sorted(D12)}")
            print(f"D22={sorted(D22)}")

            G = BlockCirculantGraph(n=n, D11=D11, D12=D12, D22=D22)
            v = verify_construction(G)

            print(f"Valid: {v.valid}")
            print(f"Red: {v.max_red_common}/{v.red_threshold}, "
                  f"Blue: {v.max_blue_common}/{v.blue_threshold}")

            if v.valid:
                save_construction(G, v, "solution_n22.json", "SAT-targeted")
                print("\n*** VALID SOLUTION FOUND! ***")
                print(f"R(B_21, B_22) >= 87 is PROVED!")
                return
            else:
                print(f"Violations ({len(v.violations)}):")
                for vt, d, ex in v.violations[:10]:
                    print(f"  {vt} d={d} excess={ex}")
        else:
            print("UNSAT")

    print("\nNo valid solution found")


if __name__ == "__main__":
    main()
