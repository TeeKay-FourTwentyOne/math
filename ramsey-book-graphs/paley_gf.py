"""
Paley-type constructions for R(B_{n-1}, B_n) over Galois fields.

When m = 2n-1 is a prime power with m ≡ 1 (mod 4), the quadratic residues
of GF(m) give a construction proving R(B_{n-1}, B_n) >= 4n - 1.

The Paley construction:
  D11 = QR(m) (quadratic residues of GF(m))
  D22 = QNR(m) = complement of QR in {1,...,m-1} (quadratic non-residues)
  D12 = {0} ∪ {x : some condition based on a fixed non-residue}

For prime p ≡ 1 (mod 4): -1 is a QR, so QR is symmetric.
For prime power q = p^k ≡ 1 (mod 4): work in GF(q).

Usage:
    python paley_gf.py --n 25   # m=49=7^2, GF(49)
    python paley_gf.py --n 14   # m=27=3^3, GF(27)
"""

import sys, os, json, argparse
from typing import Set, List, Tuple, Optional, Dict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ramsey_core import (
    BlockCirculantGraph, verify_construction, save_construction
)


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n < 4:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True


def prime_power_factorization(n: int) -> Optional[Tuple[int, int]]:
    """If n = p^k for prime p, return (p, k). Otherwise return None."""
    for p in range(2, int(n**0.5) + 2):
        if n % p == 0:
            k = 0
            val = 1
            while val < n:
                val *= p
                k += 1
            if val == n and is_prime(p):
                return (p, k)
            return None
    if is_prime(n):
        return (n, 1)
    return None


class GaloisField:
    """
    Arithmetic in GF(p^k) represented as polynomials over GF(p)
    modulo an irreducible polynomial.

    Elements are stored as tuples of length k (coefficients of polynomial).
    """

    def __init__(self, p: int, k: int):
        self.p = p
        self.k = k
        self.q = p ** k
        self.irred = self._find_irreducible()
        self._elements = None
        self._log_table = None
        self._exp_table = None

    def _find_irreducible(self) -> Tuple[int, ...]:
        """Find an irreducible polynomial of degree k over GF(p)."""
        if self.k == 1:
            return (0, 1)  # x (trivial, not actually used for prime fields)

        # For common cases, use known irreducible polynomials
        p, k = self.p, self.k
        if p == 7 and k == 2:
            # x^2 + 1 is irreducible over GF(7) since -1 is not a QR mod 7
            # Check: 0^2=0, 1^2=1, 2^2=4, 3^2=2, 4^2=2, 5^2=4, 6^2=1 mod 7
            # -1 = 6 mod 7, not in {0,1,2,4} => x^2+1 irreducible
            return (1, 0, 1)  # 1 + 0*x + 1*x^2
        if p == 3 and k == 3:
            # x^3 + 2x + 1 is irreducible over GF(3)
            return (1, 2, 0, 1)  # 1 + 2x + 0x^2 + x^3
        if p == 5 and k == 2:
            # x^2 + 2 is irreducible over GF(5) since -2 = 3 is not a QR mod 5
            return (2, 0, 1)
        if p == 2 and k == 2:
            return (1, 1, 1)  # x^2 + x + 1
        if p == 2 and k == 3:
            return (1, 1, 0, 1)  # x^3 + x + 1

        # General search for irreducible polynomial
        return self._search_irreducible()

    def _search_irreducible(self) -> Tuple[int, ...]:
        """Search for an irreducible polynomial of degree k over GF(p)."""
        p, k = self.p, self.k

        def poly_to_tuple(coeffs):
            return tuple(c % p for c in coeffs)

        def poly_eval_mod(poly, x, mod_poly):
            """Evaluate polynomial at x modulo mod_poly, all over GF(p)."""
            # Not needed for irreducibility check, using different approach
            pass

        # Enumerate monic polynomials of degree k
        def next_poly(coeffs):
            """Generate next polynomial (increment least significant coeff)."""
            result = list(coeffs)
            for i in range(len(result) - 1):  # skip leading coeff
                result[i] = (result[i] + 1) % p
                if result[i] != 0:
                    return tuple(result)
            return None

        # Start from x^k + 0*x^(k-1) + ... + 0*x + 1
        coeffs = [0] * (k + 1)
        coeffs[k] = 1  # monic
        coeffs[0] = 1  # constant term nonzero for irreducibility

        for _ in range(p ** k):
            poly = tuple(coeffs)
            if self._is_irreducible_poly(poly):
                return poly
            # Increment
            for i in range(k):
                coeffs[i] = (coeffs[i] + 1) % p
                if coeffs[i] != 0:
                    break

        raise RuntimeError(f"No irreducible polynomial found for GF({p}^{k})")

    def _is_irreducible_poly(self, poly: Tuple[int, ...]) -> bool:
        """Check if polynomial is irreducible over GF(p) using brute force for small fields."""
        p, k = self.p, self.k

        # A polynomial of degree k is irreducible iff it has no roots in GF(p)
        # AND is not a product of irreducible polys of smaller degree.
        # For degree 2 and 3, checking no roots suffices.
        if k <= 3:
            for x in range(p):
                val = 0
                xpow = 1
                for c in poly:
                    val = (val + c * xpow) % p
                    xpow = (xpow * x) % p
                if val == 0:
                    return False
            return True

        # For higher degree, use the standard criterion:
        # f is irreducible iff gcd(f, x^(p^i) - x) = 1 for i = 1..k//2
        # and f | (x^(p^k) - x)
        # This is more complex; for now use brute force factoring
        return self._brute_force_irreducibility(poly)

    def _brute_force_irreducibility(self, poly):
        """Check irreducibility by trying all factor pairs."""
        p, k = self.p, self.k
        # Try dividing by all monic polynomials of degree 1 to k//2
        for d in range(1, k // 2 + 1):
            for divisor in self._all_monic_polys(d):
                q, r = self._poly_divmod(poly, divisor)
                if all(c == 0 for c in r):
                    return False
        return True

    def _all_monic_polys(self, degree):
        """Generate all monic polynomials of given degree over GF(p)."""
        p = self.p
        count = p ** degree
        for i in range(count):
            coeffs = []
            val = i
            for _ in range(degree):
                coeffs.append(val % p)
                val //= p
            coeffs.append(1)  # monic
            yield tuple(coeffs)

    def _poly_divmod(self, a, b):
        """Polynomial division a / b over GF(p). Returns (quotient, remainder)."""
        p = self.p
        a = list(a)
        b = list(b)
        # Remove trailing zeros
        while len(a) > 1 and a[-1] == 0:
            a.pop()
        while len(b) > 1 and b[-1] == 0:
            b.pop()

        if len(a) < len(b):
            return (0,), tuple(a)

        inv_lead = pow(b[-1], p - 2, p)  # modular inverse of leading coeff of b
        q = [0] * (len(a) - len(b) + 1)
        r = list(a)

        for i in range(len(q) - 1, -1, -1):
            if len(r) >= len(b) + i:
                coeff = (r[len(b) + i - 1] * inv_lead) % p
                q[i] = coeff
                for j in range(len(b)):
                    r[i + j] = (r[i + j] - coeff * b[j]) % p

        # Trim remainder
        while len(r) > 1 and r[-1] == 0:
            r.pop()
        return tuple(q), tuple(r)

    def _build_tables(self):
        """Build log/exp tables using a primitive element."""
        if self._elements is not None:
            return

        q = self.q
        p, k = self.p, self.k

        if k == 1:
            # Prime field: elements are 0..p-1
            self._elements = list(range(p))
            # Find primitive root
            g = self._find_primitive_root_prime()
            self._exp_table = [0] * (q - 1)
            self._log_table = [0] * q
            val = 1
            for i in range(q - 1):
                self._exp_table[i] = val
                self._log_table[val] = i
                val = (val * g) % p
            return

        # Extension field: elements as integer-encoded tuples
        # Element (a0, a1, ..., a_{k-1}) encoded as a0 + a1*p + a2*p^2 + ...
        self._elements = list(range(q))

        # Find primitive element by trying generator x, x+1, etc.
        g = self._find_primitive_element()
        self._exp_table = [0] * (q - 1)
        self._log_table = [0] * q
        val = 1  # encoded as integer
        for i in range(q - 1):
            self._exp_table[i] = val
            self._log_table[val] = i
            val = self._gf_mul(val, g)

    def _encode(self, coeffs: Tuple[int, ...]) -> int:
        """Encode polynomial coefficients as integer."""
        val = 0
        mult = 1
        for c in coeffs:
            val += (c % self.p) * mult
            mult *= self.p
        return val

    def _decode(self, val: int) -> Tuple[int, ...]:
        """Decode integer to polynomial coefficients."""
        coeffs = []
        for _ in range(self.k):
            coeffs.append(val % self.p)
            val //= self.p
        return tuple(coeffs)

    def _gf_add(self, a: int, b: int) -> int:
        """Add two elements in GF(p^k)."""
        ac = self._decode(a)
        bc = self._decode(b)
        return self._encode(tuple((ac[i] + bc[i]) % self.p for i in range(self.k)))

    def _gf_mul(self, a: int, b: int) -> int:
        """Multiply two elements in GF(p^k)."""
        p, k = self.p, self.k
        ac = self._decode(a)
        bc = self._decode(b)

        # Polynomial multiplication
        prod = [0] * (2 * k - 1)
        for i in range(k):
            for j in range(k):
                prod[i + j] = (prod[i + j] + ac[i] * bc[j]) % p

        # Reduce modulo irreducible polynomial
        irred = self.irred
        for i in range(2 * k - 2, k - 1, -1):
            if prod[i] != 0:
                # Subtract prod[i] * irred * x^(i-k)
                coeff = prod[i]
                inv_lead = pow(irred[k], p - 2, p)
                factor = (coeff * inv_lead) % p
                for j in range(k + 1):
                    prod[i - k + j] = (prod[i - k + j] - factor * irred[j]) % p

        return self._encode(tuple(prod[:k]))

    def _find_primitive_root_prime(self) -> int:
        """Find primitive root of GF(p)."""
        p = self.p
        if p == 2:
            return 1
        # Factor p-1
        factors = set()
        n = p - 1
        d = 2
        while d * d <= n:
            while n % d == 0:
                factors.add(d)
                n //= d
            d += 1
        if n > 1:
            factors.add(n)

        for g in range(2, p):
            if all(pow(g, (p - 1) // f, p) != 1 for f in factors):
                return g
        raise RuntimeError(f"No primitive root found for GF({p})")

    def _find_primitive_element(self) -> int:
        """Find a primitive element of GF(p^k)."""
        q = self.q
        # Factor q-1
        factors = set()
        n = q - 1
        d = 2
        while d * d <= n:
            while n % d == 0:
                factors.add(d)
                n //= d
            d += 1
        if n > 1:
            factors.add(n)

        # Try candidates: x, x+1, x+2, ...
        for candidate_int in range(2, q):
            if self._is_primitive(candidate_int, factors):
                return candidate_int
        raise RuntimeError(f"No primitive element found for GF({q})")

    def _is_primitive(self, g: int, factors: Set[int]) -> bool:
        """Check if g is a primitive element (generates multiplicative group)."""
        q = self.q
        for f in factors:
            # Compute g^((q-1)/f)
            exp = (q - 1) // f
            val = self._gf_pow(g, exp)
            if val == 1:  # 1 encoded as (1, 0, ..., 0)
                return False
        return True

    def _gf_pow(self, base: int, exp: int) -> int:
        """Compute base^exp in GF(p^k)."""
        if exp == 0:
            return 1  # multiplicative identity
        result = 1
        b = base
        while exp > 0:
            if exp % 2 == 1:
                result = self._gf_mul(result, b)
            b = self._gf_mul(b, b)
            exp //= 2
        return result

    def quadratic_residues(self) -> Set[int]:
        """
        Return the set of quadratic residues in GF(q)*.

        x is a QR if x = y^2 for some y in GF(q)*.
        Equivalently, x is a QR iff x^((q-1)/2) = 1.
        """
        self._build_tables()
        q = self.q
        qr = set()
        for i in range(q - 1):
            elem = self._exp_table[i]
            if i % 2 == 0:  # even power of generator => QR
                qr.add(elem)
        return qr

    def quadratic_non_residues(self) -> Set[int]:
        """Return the set of quadratic non-residues in GF(q)*."""
        self._build_tables()
        q = self.q
        qnr = set()
        for i in range(q - 1):
            elem = self._exp_table[i]
            if i % 2 == 1:  # odd power of generator => QNR
                qnr.add(elem)
        return qnr

    def element_to_int(self, elem: int) -> int:
        """Convert GF element (encoded as int) to an integer in {0, ..., q-1}.

        For the circulant construction, we need a bijection from GF(q) to Z_q.
        We use the natural encoding: element a0 + a1*x + ... maps to
        a0 + a1*p + a2*p^2 + ...
        """
        return elem  # already encoded this way

    def int_to_element(self, i: int) -> int:
        """Convert integer in {0, ..., q-1} to GF element."""
        return i


def paley_construction(n: int) -> Optional[BlockCirculantGraph]:
    """
    Construct a Paley-type graph for R(B_{n-1}, B_n).

    For m = 2n-1 a prime power with m ≡ 1 (mod 4):
    - D11 = QR(m) (quadratic residues in GF(m)*, mapped to Z_m)
    - D22 = complement(D11) = QNR(m)
    - D12 = {0} ∪ QR(m) or {0} ∪ QNR(m) (try both)

    The key property of Paley graphs when m ≡ 1 (mod 4):
    -1 is a QR, so D11 = QR is symmetric under negation.
    """
    m = 2 * n - 1
    N = 4 * n - 2

    factorization = prime_power_factorization(m)
    if factorization is None:
        print(f"m={m} is not a prime power, Paley construction not applicable.")
        return None

    p, k = factorization
    if m % 4 != 1:
        print(f"m={m} ≡ {m % 4} (mod 4), need m ≡ 1 (mod 4) for Paley.")
        return None

    print(f"Constructing Paley graph for n={n}, m={m}={p}^{k}")

    gf = GaloisField(p, k)

    # Get QR and QNR
    qr_elements = gf.quadratic_residues()
    qnr_elements = gf.quadratic_non_residues()

    # Map to Z_m indices
    qr_set = {gf.element_to_int(e) for e in qr_elements}
    qnr_set = {gf.element_to_int(e) for e in qnr_elements}

    # Remove 0 (not in QR or QNR by definition, but ensure)
    qr_set.discard(0)
    qnr_set.discard(0)

    print(f"  |QR| = {len(qr_set)}, |QNR| = {len(qnr_set)}")
    print(f"  QR ∩ QNR = {len(qr_set & qnr_set)} (should be 0)")
    print(f"  |QR| + |QNR| = {len(qr_set) + len(qnr_set)} (should be {m-1})")

    # Check symmetry: -x ∈ QR iff x ∈ QR (when m ≡ 1 mod 4)
    sym_violations = sum(1 for x in qr_set if (m - x) % m not in qr_set)
    print(f"  QR symmetry violations: {sym_violations}")

    # Try different D12 constructions
    best_graph = None
    best_result = None

    # Construction 1: D11 = QR, D12 = {0} ∪ (first (m-3)/2 of QR)
    # Construction 2: D11 = QR, D12 = {0} ∪ (first (m-3)/2 of QNR)
    # Construction 3: Standard Paley: D12 = {0} ∪ QR[:(m-1)/2 - 1]
    # Construction 4: D12 based on half-QR + half-QNR

    d12_size = (m - 1) // 2  # = n - 1

    # The Paley construction for book graphs:
    # D11 = QR (size (m-1)/2)
    # D22 = QNR (size (m-1)/2)
    # D12 needs size (m-1)/2 and must include 0

    # Primary: |D11| = (m+1)/2, but QR has only (m-1)/2 elements.
    # So for the standard Paley, |D11| = (m-1)/2 = n-1, which doesn't match
    # the primary config (needs (m+1)/2 = n).
    # Instead, this matches the "equal" case: |D11| = |D22| = (m-1)/2.

    # For the book graph Paley construction (Rousseau-Sheehan style):
    # We actually need d1 = d2 = m-1 for the Paley tournament approach,
    # OR we use a modified construction.

    # Standard approach: try D11 = QR, D22 = QNR, D12 = {0} ∪ S
    # where S ⊂ {1,...,m-1} with |S| = d12_size - 1

    configs_to_try = []

    # Config A: D12 = {0} ∪ first (d12_size-1) QR elements
    qr_sorted = sorted(qr_set)
    if len(qr_sorted) >= d12_size - 1:
        d12_a = {0} | set(qr_sorted[:d12_size - 1])
        configs_to_try.append(("QR-based D12", d12_a))

    # Config B: D12 = {0} ∪ first (d12_size-1) QNR elements
    qnr_sorted = sorted(qnr_set)
    if len(qnr_sorted) >= d12_size - 1:
        d12_b = {0} | set(qnr_sorted[:d12_size - 1])
        configs_to_try.append(("QNR-based D12", d12_b))

    # Config C: D12 = {0} ∪ mixed (half QR, half QNR)
    half = (d12_size - 1) // 2
    d12_c = {0} | set(qr_sorted[:half]) | set(qnr_sorted[:d12_size - 1 - half])
    configs_to_try.append(("Mixed D12", d12_c))

    # Config D: D12 = {0, 1, 2, ..., d12_size-1} (consecutive)
    d12_d = set(range(d12_size))
    configs_to_try.append(("Consecutive D12", d12_d))

    for label, D12 in configs_to_try:
        G = BlockCirculantGraph(n=n, D11=qr_set, D12=D12, D22=qnr_set)
        result = verify_construction(G)
        cost = sum(e for _, _, e in result.violations)
        print(f"\n  {label}: |D12|={len(D12)}, valid={result.valid}, cost={cost}")
        print(f"    d1={len(G.D11)+len(G.D12)}, d2={len(G.D22)+len(G.D12)}")
        print(f"    max_red={result.max_red_common}, max_blue={result.max_blue_common}")

        if result.valid:
            if best_result is None or cost < sum(e for _, _, e in best_result.violations):
                best_graph = G
                best_result = result
        elif best_result is None or (not best_result.valid and cost < sum(e for _, _, e in best_result.violations)):
            best_graph = G
            best_result = result

    return best_graph, best_result


def paley_sa_hybrid(n: int, max_iter: int = 2000000, n_trials: int = 8):
    """
    Use Paley QR as D11 (fixed) and search for D12 using SA.

    This combines the algebraic structure of Paley with heuristic search
    for the cross-block edges.
    """
    m = 2 * n - 1
    factorization = prime_power_factorization(m)
    if factorization is None:
        return None

    p, k = factorization
    gf = GaloisField(p, k)
    qr_set = {gf.element_to_int(e) for e in gf.quadratic_residues()}
    qr_set.discard(0)

    print(f"\nPaley-SA hybrid: n={n}, m={m}, |QR|={len(qr_set)}")
    print(f"  D11 = QR (fixed), searching for optimal D12")

    # Import SA search machinery
    from sa_solver import FastEvaluator
    import random, math

    half_m = m // 2
    pairs = [(d, m - d) for d in range(1, half_m + 1)]

    # Build D11_half from QR
    D11_half = set()
    for i, (d, md) in enumerate(pairs):
        if d in qr_set:
            D11_half.add(i)

    d12_target = (m - 1) // 2  # standard size

    best_global_cost = float('inf')
    best_D12 = None

    for trial in range(n_trials):
        # Random D12 initialization
        d12_rest = list(range(1, m))
        random.shuffle(d12_rest)
        D12 = {0} | set(d12_rest[:d12_target - 1])

        ev = FastEvaluator(n, D11_half, D12)
        cost, nviols, viols = ev.compute_cost()
        best_cost = cost
        best_trial_D12 = D12.copy()

        temp = 5.0
        cooling = 0.999997

        for it in range(max_iter):
            if cost == 0:
                print(f"  SOLUTION FOUND at trial={trial}, iter={it}")
                ev._rebuild()
                return ev.D11, ev.D12, ev.D22

            # Only D12 swaps (D11 fixed as QR)
            out_list = [x for x in D12 if x != 0]
            in_list = list(set(range(1, m)) - D12)
            if not out_list or not in_list:
                continue
            d_out = random.choice(out_list)
            d_in = random.choice(in_list)

            D12.discard(d_out)
            D12.add(d_in)
            ev.D12 = D12
            ev._rebuild()
            new_cost, _, _ = ev.compute_cost()

            delta = new_cost - cost
            if delta < 0 or random.random() < math.exp(-delta / max(temp, 0.001)):
                cost = new_cost
                if cost < best_cost:
                    best_cost = cost
                    best_trial_D12 = D12.copy()
            else:
                D12.discard(d_in)
                D12.add(d_out)
                ev.D12 = D12
                ev._rebuild()

            temp *= cooling

            if it % 200000 == 0:
                print(f"  t{trial} i{it}: cost={cost} best={best_cost} T={temp:.4f}")
                sys.stdout.flush()

        print(f"  Trial {trial}: best={best_cost}")
        if best_cost < best_global_cost:
            best_global_cost = best_cost
            best_D12 = best_trial_D12.copy()

    if best_global_cost == 0:
        ev = FastEvaluator(n, D11_half, best_D12)
        return ev.D11, ev.D12, ev.D22
    else:
        print(f"\n  Best cost: {best_global_cost}")
        return None


def main():
    parser = argparse.ArgumentParser(
        description="Paley construction for R(B_{n-1}, B_n)"
    )
    parser.add_argument("--n", type=int, required=True,
                        help="Book graph parameter n")
    parser.add_argument("--sa-hybrid", action="store_true",
                        help="Use Paley-SA hybrid (fix D11=QR, search D12)")
    parser.add_argument("--trials", type=int, default=8,
                        help="Trials for SA hybrid (default: 8)")
    parser.add_argument("--max-iter", type=int, default=2000000,
                        help="Max iterations for SA hybrid (default: 2000000)")
    args = parser.parse_args()

    n = args.n
    m = 2 * n - 1

    factorization = prime_power_factorization(m)
    if factorization is None:
        print(f"m={m} is not a prime power. Paley construction not applicable.")
        sys.exit(1)

    p, k = factorization
    print(f"n={n}, m={m}={p}^{k}, m mod 4 = {m % 4}")

    if m % 4 != 1:
        print(f"WARNING: m ≡ {m % 4} (mod 4), standard Paley requires m ≡ 1 (mod 4)")

    if args.sa_hybrid:
        result = paley_sa_hybrid(n, max_iter=args.max_iter, n_trials=args.trials)
        if result is not None:
            D11, D12, D22 = result
            G = BlockCirculantGraph(n=n, D11=D11, D12=D12, D22=D22)
            v = verify_construction(G)
            print(f"\n*** SOLUTION: valid={v.valid} ***")
            print(f"D11={sorted(D11)}")
            print(f"D12={sorted(D12)}")
            print(f"D22={sorted(D22)}")
            if v.valid:
                outfile = os.path.join(
                    os.path.dirname(os.path.abspath(__file__)),
                    f"solution_n{n}.json"
                )
                save_construction(G, v, outfile, solver="Paley-SA-hybrid")
                print(f"Saved to {outfile}")
    else:
        result = paley_construction(n)
        if result is not None:
            G, v = result
            if v.valid:
                print(f"\n*** VALID CONSTRUCTION FOUND ***")
                outfile = os.path.join(
                    os.path.dirname(os.path.abspath(__file__)),
                    f"solution_n{n}.json"
                )
                save_construction(G, v, outfile, solver="Paley-GF")
                print(f"Saved to {outfile}")
            else:
                print(f"\nNo valid pure Paley construction found.")
                print(f"Try --sa-hybrid to search for D12 with SA.")


if __name__ == "__main__":
    main()
