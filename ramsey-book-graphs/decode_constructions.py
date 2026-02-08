"""
Decode graph6 strings from the Lidicky et al. / Steven-VO repository
to extract 2-block circulant structure for R(B_{n-1}, B_n) constructions.

Focus: cases where m = 2n-1 is congruent to 3 (mod 4).
"""

import urllib.request
import sys


def graph6_to_adjacency(g6):
    """Decode a graph6 string to adjacency matrix."""
    g6 = g6.strip()
    raw = [ord(c) for c in g6]

    if raw[0] == 126:  # '~' prefix: large graph (n >= 63)
        # Next 3 bytes encode n (for n < 258048)
        n = ((raw[1] - 63) << 12) | ((raw[2] - 63) << 6) | (raw[3] - 63)
        data_start = 4
    else:
        n = raw[0] - 63
        data_start = 1

    bits = []
    for i in range(data_start, len(raw)):
        d = raw[i] - 63
        for j in range(5, -1, -1):
            bits.append((d >> j) & 1)

    adj = [[0] * n for _ in range(n)]
    bit_idx = 0
    for j in range(1, n):
        for i in range(j):
            if bit_idx < len(bits):
                adj[i][j] = adj[j][i] = bits[bit_idx]
                bit_idx += 1
    return n, adj


def extract_and_verify(n_verts, adj):
    """Extract circulant structure and verify."""
    m = n_verts // 2
    D11 = set()
    for j in range(1, m):
        if adj[0][j]:
            D11.add(j)
    D12 = set()
    for j in range(m, 2 * m):
        if adj[0][j]:
            D12.add(j % m)
    D22 = set()
    for j in range(m + 1, 2 * m):
        if adj[m][j]:
            D22.add((j - m) % m)

    # Full circulant verification
    is_circ = True
    for u in range(m):
        for v in range(u + 1, m):
            d = (v - u) % m
            if adj[u][v] != (d in D11):
                is_circ = False
                return D11, D12, D22, is_circ
    for u in range(m, 2 * m):
        for v in range(u + 1, 2 * m):
            d = (v - u) % m
            if adj[u][v] != (d in D22):
                is_circ = False
                return D11, D12, D22, is_circ
    for u in range(m):
        for v in range(m, 2 * m):
            d = (v - u) % m
            if adj[u][v] != (d in D12):
                is_circ = False
                return D11, D12, D22, is_circ
    return D11, D12, D22, is_circ


def analyze_construction(n_val, g6):
    """Analyze a single construction."""
    m = 2 * n_val - 1
    mod4 = m % 4
    is_target = (mod4 == 3)

    print(f"\n{'=' * 70}")
    print(f"n={n_val}, m={m}, m mod 4 = {mod4}", end="")
    if is_target:
        print(" *** 3 (mod 4) ***")
    else:
        print(" [1 (mod 4) - Paley applies]")
    print(f"{'=' * 70}")

    nv, adj = graph6_to_adjacency(g6)
    expected = 4 * n_val - 2
    print(f"Vertices: {nv} (expected {expected})")

    if nv != expected:
        print("  VERTEX COUNT MISMATCH!")
        return

    D11, D12, D22, is_circ = extract_and_verify(nv, adj)
    print(f"2-block circulant: {is_circ}")
    print(f"D11 = {sorted(D11)} (size {len(D11)})")
    print(f"D12 = {sorted(D12)} (size {len(D12)})")
    print(f"D22 = {sorted(D22)} (size {len(D22)})")

    complement = set(range(1, m)) - D11
    print(f"D22 == complement(D11): {D22 == complement}")

    D11_sym = all((-d % m) in D11 for d in D11)
    D22_sym = all((-d % m) in D22 for d in D22)
    print(f"D11 symmetric: {D11_sym}, D22 symmetric: {D22_sym}")

    d1 = len(D11) + len(D12)
    d2 = len(D22) + len(D12)
    print(f"Degrees: d1={d1}, d2={d2}")

    # Quadratic residue analysis
    if m > 2 and is_prime(m):
        QR = {pow(x, 2, m) for x in range(1, m)}
        QNR = set(range(1, m)) - QR
        print(f"QR_{m} = {sorted(QR)} (size {len(QR)})")

        d11_qr = len(D11 & QR)
        d11_qnr = len(D11 & QNR)
        print(f"D11: {d11_qr} in QR, {d11_qnr} in QNR")

        d12_nz = D12 - {0}
        d12_qr = len(d12_nz & QR)
        d12_qnr = len(d12_nz & QNR)
        print(f"D12 (nonzero): {d12_qr} in QR, {d12_qnr} in QNR, 0 in D12: {0 in D12}")

        # Check if D12 is QR, QNR, or close
        if d12_nz == QR:
            print("  >>> D12 \\ {0} == QR <<<")
        elif d12_nz == QNR:
            print("  >>> D12 \\ {0} == QNR <<<")

        if D11 == QR:
            print("  >>> D11 == QR <<<")

        # Compute Delta(D11, D11, d) for a few d values
        print(f"\nDelta(D11, D11, d) values:")
        for d in range(1, min(m, 10)):
            delta = sum(1 for a in D11 if (a - d) % m in D11)
            print(f"  d={d}: {delta}", end="")
        print()

        # Compute Delta(D12, D12, d) for a few d values
        print(f"Delta(D12, D12, d) values:")
        for d in range(1, min(m, 10)):
            delta = sum(1 for a in D12 if (a - d) % m in D12)
            print(f"  d={d}: {delta}", end="")
        print()


def is_prime(n):
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


def main():
    # Try to fetch the data from the repository
    url = "https://raw.githubusercontent.com/Steven-VO/circulant-Ramsey/master/RamseyGraphs/Books/Bn-1Bn.txt"

    print("Fetching graph6 data from repository...")
    try:
        req = urllib.request.urlopen(url, timeout=15)
        raw_data = req.read().decode('utf-8')
        print(f"Fetched {len(raw_data)} bytes")
    except Exception as e:
        print(f"Could not fetch: {e}")
        print("Using hardcoded data for available cases...")
        raw_data = None

    if raw_data:
        # Parse the file: lines alternate between "R(Bx,By) > z" and graph6 strings
        lines = [l.strip() for l in raw_data.split('\n') if l.strip()]
        constructions = {}
        i = 0
        while i < len(lines) - 1:
            header = lines[i]
            g6 = lines[i + 1]
            # Parse n from "R(B_{n-1}, B_n) > 4n-2"
            # Format: R(Bx,By) > z
            if header.startswith('R(B'):
                parts = header.split(')')
                inner = parts[0][2:]  # "Bx,By"
                bvals = inner.split(',')
                b1 = int(bvals[0][1:])  # x from "Bx"
                b2 = int(bvals[1][1:])  # y from "By"
                n_val = b2  # For R(B_{n-1}, B_n), n = b2
                constructions[n_val] = g6
            i += 2

        print(f"\nParsed {len(constructions)} constructions")
        print(f"n values: {sorted(constructions.keys())}")

        # Analyze target cases (m = 3 mod 4)
        target_n = []
        for n_val in sorted(constructions.keys()):
            m = 2 * n_val - 1
            if m % 4 == 3:
                target_n.append(n_val)

        print(f"\nTarget n values (m = 3 mod 4): {target_n}")

        # Analyze ALL constructions
        for n_val in sorted(constructions.keys()):
            m = 2 * n_val - 1
            if m % 4 == 3:  # Only the target cases
                try:
                    analyze_construction(n_val, constructions[n_val])
                except Exception as e:
                    print(f"Error for n={n_val}: {e}")

        # Also show summary for m=1 mod 4 cases for comparison
        print(f"\n\n{'#' * 70}")
        print("COMPARISON: m = 1 (mod 4) cases (Paley construction applies)")
        print(f"{'#' * 70}")
        for n_val in sorted(constructions.keys()):
            m = 2 * n_val - 1
            if m % 4 == 1:
                try:
                    analyze_construction(n_val, constructions[n_val])
                except Exception as e:
                    print(f"Error for n={n_val}: {e}")


if __name__ == "__main__":
    main()
