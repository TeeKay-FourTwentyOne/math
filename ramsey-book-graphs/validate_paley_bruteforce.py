"""
Standalone brute-force validator for Paley constructions on GF(q), q = p^k, q ≡ 1 (mod 4).

For given q, constructs the 2-block graph with D11 = D12 = QR(GF(q)*), D22 = QNR(GF(q)*),
builds the full 2q x 2q adjacency matrix, and checks all C(2q, 2) pairs.

Includes full GF(p^k) arithmetic from scratch. NO external dependencies.

Usage:
  python validate_paley_bruteforce.py --q 81
  python validate_paley_bruteforce.py --q 97
"""

import sys
import argparse
import time
from itertools import combinations, product as iprod


# === GF(p^k) ARITHMETIC ===

def find_irreducible(p, k):
    """Find an irreducible polynomial of degree k over GF(p).
    Returns list of coefficients [a0, a1, ..., ak] (monic, so ak=1)."""
    if k == 1:
        return [0, 1]  # x is irreducible

    for coeffs in iprod(range(p), repeat=k):
        poly = list(coeffs) + [1]  # monic
        if is_irreducible(poly, p, k):
            return poly
    return None


def is_irreducible(poly, p, k):
    """Check if monic polynomial of degree k is irreducible over GF(p)."""
    if k <= 0:
        return False
    if k == 1:
        return True

    # Check no roots in GF(p)
    for x in range(p):
        val = 0
        xpow = 1
        for c in poly:
            val = (val + c * xpow) % p
            xpow = (xpow * x) % p
        if val == 0:
            return False

    if k == 2:
        return True  # No roots => irreducible for degree 2

    # For higher degrees, use: f is irreducible iff
    # f divides x^{p^k} - x AND gcd(f, x^{p^{k/r}} - x) = 1 for prime divisors r of k
    prime_divisors = set()
    temp = k
    for d in range(2, k + 1):
        if d * d > temp:
            break
        while temp % d == 0:
            prime_divisors.add(d)
            temp //= d
    if temp > 1:
        prime_divisors.add(temp)

    x_poly = [0, 1]

    for r in prime_divisors:
        exp = p ** (k // r)
        xp = poly_pow_mod(x_poly, exp, poly, p)
        diff = list(xp)
        while len(diff) < 2:
            diff.append(0)
        diff[1] = (diff[1] - 1) % p
        g = poly_gcd(poly, diff, p)
        while len(g) > 1 and g[-1] == 0:
            g.pop()
        if len(g) > 1:
            return False

    return True


def poly_mod(a, b, p):
    """Compute a mod b over GF(p)[x]."""
    a = list(a)
    b = list(b)
    while len(b) > 1 and b[-1] % p == 0:
        b.pop()
    db = len(b) - 1
    if db == 0 and b[0] % p == 0:
        return a
    if len(a) - 1 < db:
        result = [c % p for c in a]
        while len(result) > 1 and result[-1] == 0:
            result.pop()
        return result if result else [0]
    lead_inv = pow(b[db] % p, p - 2, p)
    for i in range(len(a) - 1, db - 1, -1):
        if a[i] % p != 0:
            c = (a[i] * lead_inv) % p
            for j in range(db + 1):
                a[i - db + j] = (a[i - db + j] - c * b[j]) % p
    result = [a[i] % p for i in range(min(db, len(a)))]
    while len(result) < db:
        result.append(0)
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    if not result:
        result = [0]
    return result


def poly_mul_mod(a, b, mod, p):
    """Multiply polynomials a*b mod 'mod' over GF(p)."""
    prod = [0] * (len(a) + len(b) - 1)
    for i in range(len(a)):
        for j in range(len(b)):
            prod[i + j] = (prod[i + j] + a[i] * b[j]) % p
    return poly_mod(prod, mod, p)


def poly_pow_mod(base, exp, mod, p):
    """Compute base^exp mod 'mod' over GF(p) using repeated squaring."""
    result = [1]
    base = poly_mod(list(base), mod, p)
    while exp > 0:
        if exp % 2 == 1:
            result = poly_mul_mod(result, base, mod, p)
        base = poly_mul_mod(base, base, mod, p)
        exp //= 2
    return result


def poly_gcd(a, b, p):
    """Compute GCD of polynomials over GF(p)."""
    a = [x % p for x in a]
    b = [x % p for x in b]
    while True:
        while len(b) > 1 and b[-1] == 0:
            b.pop()
        if len(b) == 1 and b[0] == 0:
            break
        a, b = b, poly_mod(a, b, p)
    return a


class GF:
    """GF(p^k) with elements as tuples of length k."""

    def __init__(self, p, k):
        self.p = p
        self.k = k
        self.q = p ** k
        self.modpoly = find_irreducible(p, k) if k > 1 else None

    def zero(self):
        return tuple([0] * self.k)

    def one(self):
        r = [0] * self.k
        r[0] = 1
        return tuple(r)

    def add(self, a, b):
        return tuple((a[i] + b[i]) % self.p for i in range(self.k))

    def sub(self, a, b):
        return tuple((a[i] - b[i]) % self.p for i in range(self.k))

    def neg(self, a):
        return tuple((-a[i]) % self.p for i in range(self.k))

    def mul(self, a, b):
        if self.k == 1:
            return ((a[0] * b[0]) % self.p,)
        prod = [0] * (2 * self.k - 1)
        for i in range(self.k):
            for j in range(self.k):
                prod[i + j] = (prod[i + j] + a[i] * b[j]) % self.p
        for i in range(2 * self.k - 2, self.k - 1, -1):
            if prod[i] != 0:
                c = prod[i]
                for j in range(self.k + 1):
                    prod[i - self.k + j] = (prod[i - self.k + j] - c * self.modpoly[j]) % self.p
        return tuple(prod[i] % self.p for i in range(self.k))

    def is_zero(self, a):
        return all(x == 0 for x in a)

    def to_int(self, elem):
        val = 0
        for i in range(self.k):
            val += elem[i] * (self.p ** i)
        return val

    def from_int(self, n):
        coeffs = []
        for _ in range(self.k):
            coeffs.append(n % self.p)
            n //= self.p
        return tuple(coeffs)

    def order(self, a):
        if self.is_zero(a):
            return 0
        power = a
        for k in range(1, self.q):
            if power == self.one():
                return k
            power = self.mul(power, a)
        return -1

    def find_generator(self):
        """Find a primitive element of GF(q)*."""
        for n in range(2, self.q):
            elem = self.from_int(n)
            if self.order(elem) == self.q - 1:
                return elem
        return None

    def compute_qr(self):
        """Compute quadratic residues as set of ints."""
        gen = self.find_generator()
        assert gen is not None
        qr = set()
        power = self.one()
        for k in range(self.q - 1):
            if k % 2 == 0:
                qr.add(self.to_int(power))
            power = self.mul(power, gen)
        return qr


def is_prime_power(n):
    """Return (p, k) if n = p^k, else None."""
    if n < 2:
        return None
    for p in range(2, int(n ** 0.5) + 2):
        if n % p == 0:
            k, val = 0, n
            while val % p == 0:
                val //= p
                k += 1
            return (p, k) if val == 1 else None
    return (n, 1)


def validate_paley(q):
    """Full brute-force validation of Paley construction on GF(q)."""
    pp = is_prime_power(q)
    if pp is None:
        print(f"ERROR: {q} is not a prime power")
        return False
    if q % 4 != 1:
        print(f"ERROR: {q} is not ≡ 1 (mod 4)")
        return False

    p, k = pp
    n = (q + 1) // 2
    m = q
    N = 2 * m

    red_threshold = n - 2
    blue_threshold = n - 1

    print("=" * 60)
    print(f"BRUTE-FORCE PALEY VALIDATOR: R(B_{{{n-1}}}, B_{{{n}}}) = {4*n - 1}")
    print(f"  q = {q} = {p}^{k}, n = {n}, m = {m}, N = {N}")
    print("=" * 60)
    print()

    # Build GF(q)
    t0 = time.time()
    gf = GF(p, k)
    if k > 1:
        print(f"GF({q}) irreducible polynomial: {gf.modpoly}")

    # Compute QR
    print("Computing quadratic residues...")
    qr_ints = gf.compute_qr()
    all_nonzero = set(range(1, q))
    qnr_ints = all_nonzero - qr_ints

    D11_set = qr_ints
    D12_set = qr_ints
    D22_set = qnr_ints

    print(f"|QR| = {len(qr_ints)}, |QNR| = {len(qnr_ints)} (expected {(q-1)//2})")

    # Verify -1 is QR (needed for symmetry when q ≡ 1 mod 4)
    neg_one = gf.to_int(gf.neg(gf.one()))
    neg_one_is_qr = neg_one in qr_ints
    print(f"-1 = {neg_one}, is QR: {neg_one_is_qr} (expected True for q ≡ 1 mod 4)")

    # Check QR symmetry under negation
    qr_sym = all(gf.to_int(gf.neg(gf.from_int(x))) in qr_ints for x in qr_ints)
    print(f"QR symmetric under negation: {qr_sym}")

    # Check D22 = complement
    complement_ok = D22_set == all_nonzero - D11_set
    print(f"D22 = complement(D11): {complement_ok}")

    # Verify SRG property: Delta(QR,QR,d) for all nonzero d
    print("Checking SRG property...")
    expected_lambda = (q - 5) // 4
    expected_mu = (q - 1) // 4
    delta_qr = set()
    delta_qnr = set()
    for d in range(1, q):
        d_elem = gf.from_int(d)
        count = 0
        for a in qr_ints:
            a_elem = gf.from_int(a)
            diff = gf.sub(a_elem, d_elem)
            diff_int = gf.to_int(diff)
            if diff_int in qr_ints:
                count += 1
        if d in qr_ints:
            delta_qr.add(count)
        else:
            delta_qnr.add(count)

    srg_ok = delta_qr == {expected_lambda} and delta_qnr == {expected_mu}
    print(f"  Delta(QR,QR,d) for d in QR: {sorted(delta_qr)} (expected {{{expected_lambda}}}): "
          f"{'OK' if delta_qr == {expected_lambda} else 'FAIL'}")
    print(f"  Delta(QR,QR,d) for d in QNR: {sorted(delta_qnr)} (expected {{{expected_mu}}}): "
          f"{'OK' if delta_qnr == {expected_mu} else 'FAIL'}")

    d1 = len(D11_set) + len(D12_set)
    d2 = len(D22_set) + len(D12_set)
    print(f"d1 = {d1}, d2 = {d2}")
    print()

    # Build full adjacency matrix
    print(f"Building {N}x{N} adjacency matrix...")
    t1 = time.time()

    adj = [[0] * N for _ in range(N)]

    for u in range(N):
        for v in range(u + 1, N):
            u_in_V1 = u < m
            v_in_V1 = v < m

            if u_in_V1 and v_in_V1:
                diff = gf.to_int(gf.sub(gf.from_int(v), gf.from_int(u)))
                if diff in D11_set:
                    adj[u][v] = adj[v][u] = 1
            elif not u_in_V1 and not v_in_V1:
                diff = gf.to_int(gf.sub(gf.from_int(v - m), gf.from_int(u - m)))
                if diff in D22_set:
                    adj[u][v] = adj[v][u] = 1
            else:
                if u_in_V1:
                    diff = gf.to_int(gf.sub(gf.from_int(v - m), gf.from_int(u)))
                else:
                    diff = gf.to_int(gf.sub(gf.from_int(v), gf.from_int(u - m)))
                if diff in D12_set:
                    adj[u][v] = adj[v][u] = 1

    t2 = time.time()
    print(f"Adjacency matrix built in {t2 - t1:.1f}s")

    # Verify degrees
    deg = [sum(adj[u]) for u in range(N)]
    v1_degs = set(deg[:m])
    v2_degs = set(deg[m:])
    print(f"V1 degrees: {v1_degs} (expected {{{d1}}}): {'OK' if v1_degs == {d1} else 'FAIL'}")
    print(f"V2 degrees: {v2_degs} (expected {{{d2}}}): {'OK' if v2_degs == {d2} else 'FAIL'}")
    print()

    # Check ALL pairs
    total_pairs = N * (N - 1) // 2
    print(f"Checking all C({N}, 2) = {total_pairs} pairs...")
    t3 = time.time()

    red_edges = 0
    blue_edges = 0
    max_red_cn = 0
    max_blue_cn = 0
    violations = []

    for u, v in combinations(range(N), 2):
        is_red = adj[u][v] == 1

        if is_red:
            red_edges += 1
            cn = 0
            for w in range(N):
                if w != u and w != v and adj[u][w] == 1 and adj[v][w] == 1:
                    cn += 1
            if cn > max_red_cn:
                max_red_cn = cn
            if cn > red_threshold:
                violations.append(("RED", u, v, cn))
        else:
            blue_edges += 1
            cn = 0
            for w in range(N):
                if w != u and w != v and adj[u][w] == 0 and adj[v][w] == 0:
                    cn += 1
            if cn > max_blue_cn:
                max_blue_cn = cn
            if cn > blue_threshold:
                violations.append(("BLUE", u, v, cn))

    t4 = time.time()
    print(f"Pair checking completed in {t4 - t3:.1f}s")
    print()
    print(f"Total pairs: {total_pairs}")
    print(f"Red edges: {red_edges}, Blue edges: {blue_edges}")
    print(f"Max red common neighbors:  {max_red_cn} (threshold {red_threshold})")
    print(f"Max blue common neighbors: {max_blue_cn} (threshold {blue_threshold})")
    print(f"Violations: {len(violations)}")

    if violations:
        print()
        print("VIOLATIONS:")
        for color, u, v, cn in violations[:20]:
            print(f"  {color} edge ({u},{v}): {cn} common neighbors")
        if len(violations) > 20:
            print(f"  ... and {len(violations) - 20} more")

    print()
    print("=" * 60)
    structural_ok = (neg_one_is_qr and qr_sym and complement_ok and srg_ok
                     and v1_degs == {d1} and v2_degs == {d2})
    passed = structural_ok and len(violations) == 0

    if passed:
        print("RESULT: PASS")
        print()
        print(f"The Paley 2-block graph over GF({q}) on {N} vertices contains")
        print(f"no red B_{{{n-1}}} and no blue B_{{{n}}}.")
        print(f"Therefore R(B_{{{n-1}}}, B_{{{n}}}) >= {4*n - 1}.")
        print(f"Combined with the upper bound R(B_{{{n-1}}}, B_{{{n}}}) <= {4*n - 1},")
        print(f"this proves R(B_{{{n-1}}}, B_{{{n}}}) = {4*n - 1}.")
    else:
        print("RESULT: FAIL")
        if not structural_ok:
            print("  Structural check failed")
        if violations:
            print(f"  {len(violations)} constraint violation(s)")

    print("=" * 60)
    total_time = time.time() - t0
    print(f"Total time: {total_time:.1f}s")

    return passed


def main():
    parser = argparse.ArgumentParser(
        description="Brute-force Paley construction validator for GF(q)"
    )
    parser.add_argument("--q", type=int, required=True,
                        help="Prime power q ≡ 1 (mod 4)")
    args = parser.parse_args()

    ok = validate_paley(args.q)
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
