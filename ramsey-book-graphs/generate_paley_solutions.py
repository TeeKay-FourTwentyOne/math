"""
Generate explicit Paley constructions for prime q ≡ 1 (mod 4).

For prime q, QR = {x^2 mod q : x in 1..q-1} (quadratic residues).
D11 = D12 = QR, D22 = QNR = {1..q-1} - QR.
"""

import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ramsey_core import BlockCirculantGraph, verify_construction, save_construction


def quadratic_residues_prime(q):
    """Compute QR for prime q."""
    qr = set()
    for x in range(1, q):
        qr.add(pow(x, 2, q))
    return qr


def generate_paley(q):
    """Generate and verify Paley construction for prime q ≡ 1 mod 4."""
    assert q % 4 == 1, f"q={q} not ≡ 1 mod 4"

    n = (q + 1) // 2
    m = q

    qr = quadratic_residues_prime(q)
    qnr = set(range(1, q)) - qr

    assert len(qr) == (q - 1) // 2
    assert len(qnr) == (q - 1) // 2

    # D11 = QR, D12 = QR (no 0), D22 = QNR
    # But wait - for the Paley construction, 0 is NOT in D12.
    # The standard Paley has D12 = QR (without 0).
    D11 = qr
    D12 = qr  # 0 not included
    D22 = qnr

    print(f"\nq={q}, n={n}, m={m}, N={2*m}")
    print(f"|D11|={len(D11)}, |D12|={len(D12)}, |D22|={len(D22)}")

    G = BlockCirculantGraph(n=n, D11=D11, D12=D12, D22=D22)
    result = verify_construction(G)

    print(f"Valid: {result.valid}")
    print(f"Max red common: {result.max_red_common}, threshold: {result.red_threshold}")
    print(f"Max blue common: {result.max_blue_common}, threshold: {result.blue_threshold}")
    if result.violations:
        print(f"Violations: {result.violations[:5]}")

    return G, result


def main():
    # Generate for all prime q ≡ 1 mod 4 that we're missing
    targets = [53, 61]  # n=27 (q=53), n=31 (q=61)

    for q in targets:
        G, result = generate_paley(q)
        n = (q + 1) // 2

        if result.valid:
            fn = f"/Users/stephenpadgett/Projects/math/ramsey-book-graphs/solution_n{n}_paley.json"
            save_construction(G, result, fn, solver="Paley-prime")
            print(f"SAVED to {fn}")
        else:
            print(f"FAILED for q={q}")


if __name__ == "__main__":
    main()
