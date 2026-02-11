#!/usr/bin/env python3
"""Analyze which D11 sets have valid D12, and their structural properties.

Key questions:
1. Are the QR-based D11 among those with valid D12?
2. What fraction of D11 have valid D12 (does it grow with p)?
3. What structural property distinguishes working D11 from non-working?
"""

import numpy as np
import json
import os
from math import comb


def quadratic_residues(p):
    """Return the set of quadratic residues mod p."""
    return set((x * x) % p for x in range(1, p))


def is_qr_based(D11, p):
    """Check if D11 is the QR set or its complement."""
    QR = quadratic_residues(p)
    QR.discard(0)
    QNR = set(range(1, p)) - QR
    return set(D11) == QR or set(D11) == QNR


def max_autocorrelation(D11, p):
    """Compute max A(d) for d in D11."""
    indicator = np.zeros(p, dtype=np.float64)
    for j in D11:
        indicator[j] = 1.0
    fft_val = np.fft.fft(indicator)
    autocorr = np.fft.ifft(np.abs(fft_val) ** 2).real
    A = np.round(autocorr).astype(int)
    return max(int(A[d]) for d in D11)


def analyze_prime(p):
    """Load enumeration data and analyze which D11 work."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    data_path = os.path.join(script_dir, "..", "enumeration_data", f"enumeration_p{p}.json")

    with open(data_path, "r") as f:
        data = json.load(f)

    summary = data["summary"]
    details = data["per_d11"]
    n = summary["n"]
    threshold = summary["threshold_binding"]

    QR = quadratic_residues(p)
    QR.discard(0)
    QNR = set(range(1, p)) - QR

    print(f"\n{'=' * 80}")
    print(f"ANALYSIS: p = {p}, n = {n}")
    print(f"{'=' * 80}")
    print(f"QR({p}) = {sorted(QR)}")
    print(f"QNR({p}) = {sorted(QNR)}")
    print(f"Binding threshold: A(d)+B(d) <= {threshold}")

    working = []
    non_working = []

    for entry in details:
        D11 = set(entry["D11"])
        count = entry["num_valid_d12"]
        if count > 0:
            working.append(entry)
        else:
            non_working.append(entry)

    print(f"\nWorking D11: {len(working)}/{len(details)}")
    print(f"Non-working D11: {len(non_working)}/{len(details)}")

    # Check QR
    qr_found = False
    qnr_found = False
    for entry in details:
        D11 = set(entry["D11"])
        if D11 == QR:
            print(f"\nQR D11 = {sorted(QR)}: {entry['num_valid_d12']} valid D12")
            qr_found = True
        if D11 == QNR:
            print(f"QNR D11 = {sorted(QNR)}: {entry['num_valid_d12']} valid D12")
            qnr_found = True

    if not qr_found:
        print(f"\nWARNING: QR set not found among symmetric D11!")
        # QR might not be symmetric for p = 3 mod 4
        # Check: d in QR iff p-d in QR? For p = 3 mod 4, -1 is a QNR, so
        # p-d is a QR iff -d is a QR iff d is a QNR. So QR is NOT symmetric!
        print(f"  (-1) is a QR mod {p}? {((p-1)//2) % 2 == 0}")
        print(f"  For p = 3 mod 4, -1 is QNR, so QR is NOT symmetric.")
        print(f"  This means the Paley construction uses D11 that is NOT purely QR or QNR.")

    # Structural analysis of working D11
    print(f"\nWorking D11 details:")
    for entry in working:
        D11 = set(entry["D11"])
        qr_count = len(D11 & QR)
        qnr_count = len(D11 & QNR)
        max_A = max_autocorrelation(D11, p)
        A_values = entry.get("A_values", {})
        print(f"  {entry['D11']}: "
              f"{entry['num_valid_d12']} valid, "
              f"|QR|={qr_count}, |QNR|={qnr_count}, "
              f"max_A={max_A}, "
              f"A_vals={A_values}")

    # What about A-flatness? Compute max A(d) for all D11
    print(f"\nA-flatness analysis (max A(d) for d in D11):")
    max_A_working = []
    max_A_nonworking = []
    for entry in details:
        D11 = set(entry["D11"])
        mA = max_autocorrelation(D11, p)
        if entry["num_valid_d12"] > 0:
            max_A_working.append(mA)
        else:
            max_A_nonworking.append(mA)

    if max_A_working:
        print(f"  Working D11: max_A in {sorted(set(max_A_working))}")
    if max_A_nonworking:
        print(f"  Non-working D11: max_A in {sorted(set(max_A_nonworking))}")

    # Is there a clear max_A threshold?
    if max_A_working:
        cutoff = max(max_A_working)
        could_work = sum(1 for entry in details
                         if max_autocorrelation(set(entry["D11"]), p) <= cutoff)
        print(f"  Max A among working: {cutoff}")
        print(f"  D11 with max_A <= {cutoff}: {could_work}/{len(details)}")

    # For p=3 mod 4, the Paley-type construction uses D11 based on
    # quadratic residues of a MODIFIED form. Let's check:
    # The standard Paley approach for p = 3 mod 4 takes:
    # D11 = {d : chi(d) = 1} where chi is the Legendre symbol,
    # but this isn't symmetric since chi(-1) = -1.
    # Instead, one uses D11 = {d : chi(d) = 1} ∪ {d : chi(-d) = 1}
    # or some other symmetrization.

    # Actually, for the proof, we don't need ALL D11 to work.
    # We just need the expected number of valid (D11, D12) pairs to be > 0.
    # Since working D11 exist and each has multiple valid D12, this holds.

    return {
        "p": p,
        "num_working": len(working),
        "num_total": len(details),
        "fraction_working": len(working) / len(details),
    }


def main():
    print("=" * 80)
    print("ANALYSIS OF WHICH D11 SETS ADMIT VALID D12")
    print("=" * 80)

    results = {}
    for p in [11, 19, 23]:
        results[p] = analyze_prime(p)

    print(f"\n{'=' * 80}")
    print("SUMMARY")
    print(f"{'=' * 80}")
    print(f"{'p':>4s} {'working':>8s} {'total':>8s} {'fraction':>10s}")
    for p, r in results.items():
        print(f"  {p:4d} {r['num_working']:8d} {r['num_total']:8d} "
              f"{r['fraction_working']:10.4f}")

    print(f"\nKey insight: For p = 3 mod 4, QR is NOT symmetric (since -1 is QNR).")
    print(f"The Paley construction must use a different D11 (not pure QR/QNR).")
    print(f"Not all symmetric D11 have valid D12, but a significant fraction do.")
    print(f"The proof only needs E[#valid pairs] > 0, which holds since working D11 exist.")


if __name__ == "__main__":
    main()
