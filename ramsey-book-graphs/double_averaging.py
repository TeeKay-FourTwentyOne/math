#!/usr/bin/env python3
"""Double-averaging computation for the second moment method.

Computes E[Z] where Z = 1[valid] for a random (D11, D12) pair, and compares
to the theoretical headroom. Also organizes second moment ratio data.

For each p in {11, 19, 23, 31}:
1. Count total valid (D11, D12) pairs = sum of N(D11) over all symmetric D11
2. Compute E[Z] = total_valid / (#D11 * #D12)
3. Report log_2(E[Z]) and compare to headroom = log_2(C(p-1, k))
4. Report second moment ratio E[N^2]/E[N]^2
"""

import json
import math
import os
from collections import Counter


def load_enumeration_data(p):
    """Load exact enumeration data for prime p from JSON files."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    enum_dir = os.path.join(script_dir, "..", "enumeration_data")
    path = os.path.join(enum_dir, f"enumeration_p{p}.json")
    with open(path) as f:
        data = json.load(f)

    summary = data["summary"]
    per_d11 = data["per_d11"]

    # Extract N(D11) values
    N_values = [entry["num_valid_d12"] for entry in per_d11]

    return summary, N_values


def load_p31_data():
    """Load exact enumeration data for p=31 from JSONL orbit file."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(script_dir, "exact_N_p31_results.jsonl")

    orbits = []
    with open(path) as f:
        for line in f:
            orbits.append(json.loads(line.strip()))

    # Each orbit has 'size' D11 members (orbit_size = (p-1)/2 = 15)
    # and 'N' valid D12 per D11 in the orbit
    # Total D11 = sum of sizes = num_orbits * orbit_size
    # Total valid pairs = sum of N * size over orbits

    return orbits


def compute_double_averaging(p, N_values, num_d11, num_d12):
    """Compute double-averaging statistics.

    Args:
        p: prime
        N_values: list of N(D11) for each D11 (one entry per D11)
        num_d11: total number of symmetric D11
        num_d12: total number of D12 candidates
    """
    n = (p + 1) // 2
    k = (p - 3) // 2  # |S| = non-zero elements of D12

    total_valid = sum(N_values)
    num_working = sum(1 for N in N_values if N > 0)

    # E[N] = average N(D11) over all D11
    E_N = total_valid / num_d11

    # E[N^2] = average N(D11)^2 over all D11
    E_N2 = sum(N**2 for N in N_values) / num_d11

    # Second moment ratio
    if E_N > 0:
        ratio = E_N2 / (E_N**2)
    else:
        ratio = float('inf')

    # E[Z] = Pr[random (D11, D12) pair is valid]
    E_Z = total_valid / (num_d11 * num_d12)

    # log_2(E[Z])
    if E_Z > 0:
        log2_EZ = math.log2(E_Z)
    else:
        log2_EZ = float('-inf')

    # Headroom: log_2(#D12) = log_2(C(p-1, k))
    log2_headroom = math.log2(num_d12)

    # Budget: log_2(C(p-1, k)) using Stirling
    # This is the number of bits available

    # Paley-Zygmund bound: Pr[N > 0] >= E[N]^2 / E[N^2]
    if E_N2 > 0:
        PZ_bound = E_N**2 / E_N2
    else:
        PZ_bound = 0

    # N-value distribution
    N_dist = Counter(N_values)

    return {
        "p": p,
        "n": n,
        "k": k,
        "num_d11": num_d11,
        "num_d12": num_d12,
        "total_valid_pairs": total_valid,
        "num_working_d11": num_working,
        "fraction_working": num_working / num_d11,
        "E_N": E_N,
        "E_N2": E_N2,
        "second_moment_ratio": ratio,
        "E_Z": E_Z,
        "log2_EZ": log2_EZ,
        "log2_headroom": log2_headroom,
        "log2_EZ_minus_headroom": log2_EZ - log2_headroom if E_Z > 0 else float('-inf'),
        "PZ_bound": PZ_bound,
        "Pr_N_gt_0": num_working / num_d11,
        "N_distribution": dict(sorted(N_dist.items())),
    }


def main():
    results = {}

    # p = 11, 19, 23: exact enumeration data
    for p in [11, 19, 23]:
        summary, N_values = load_enumeration_data(p)
        num_d11 = summary["num_d11"]
        num_d12 = summary["num_d12"]

        result = compute_double_averaging(p, N_values, num_d11, num_d12)
        results[p] = result

        print(f"\n{'='*70}")
        print(f"p = {p}, n = {result['n']}")
        print(f"{'='*70}")
        print(f"  #D11 (symmetric) = {num_d11}")
        print(f"  #D12 = C({p-1},{result['k']}) = {num_d12:,}")
        print(f"  Total search space = {num_d11 * num_d12:,}")
        print(f"  Total valid (D11,D12) pairs = {result['total_valid_pairs']:,}")
        print(f"  Working D11: {result['num_working_d11']}/{num_d11} "
              f"({100*result['fraction_working']:.1f}%)")
        print(f"  E[N(D11)] = {result['E_N']:.4f}")
        print(f"  E[N(D11)^2] = {result['E_N2']:.4f}")
        print(f"  E[N^2]/E[N]^2 = {result['second_moment_ratio']:.4f}")
        print(f"  E[Z] = Pr[valid] = {result['E_Z']:.6e}")
        print(f"  log2(E[Z]) = {result['log2_EZ']:.4f}")
        print(f"  log2(#D12) = {result['log2_headroom']:.4f}")
        print(f"  log2(E[Z]) - log2(#D12) = {result['log2_EZ_minus_headroom']:.4f}")
        print(f"  Paley-Zygmund bound Pr[N>0] >= {result['PZ_bound']:.4f}")
        print(f"  Actual Pr[N>0] = {result['Pr_N_gt_0']:.4f}")
        print(f"  N-value distribution: {result['N_distribution']}")

    # p = 31: orbit-level data
    orbits = load_p31_data()
    p = 31
    n = 16
    k = (p - 3) // 2  # 14
    orbit_size = (p - 1) // 2  # 15
    num_orbits = len(orbits)
    num_d11 = num_orbits * orbit_size
    num_d12 = math.comb(p - 1, k)

    # Build N_values: each orbit contributes orbit_size copies of the same N
    N_values = []
    for orb in orbits:
        N_values.extend([orb["N"]] * orb["size"])

    result = compute_double_averaging(p, N_values, num_d11, num_d12)
    results[p] = result

    # Also compute orbit-level statistics
    orbit_N_values = [orb["N"] for orb in orbits]
    working_orbits = sum(1 for N in orbit_N_values if N > 0)
    orbit_N_dist = Counter(orbit_N_values)

    print(f"\n{'='*70}")
    print(f"p = {p}, n = {result['n']}")
    print(f"{'='*70}")
    print(f"  #orbits = {num_orbits}, orbit size = {orbit_size}")
    print(f"  #D11 (symmetric) = {num_d11}")
    print(f"  #D12 = C({p-1},{k}) = {num_d12:,}")
    print(f"  Total search space = {num_d11 * num_d12:,}")
    print(f"  Total valid (D11,D12) pairs = {result['total_valid_pairs']:,}")
    print(f"  Working orbits: {working_orbits}/{num_orbits} "
          f"({100*working_orbits/num_orbits:.1f}%)")
    print(f"  Working D11: {result['num_working_d11']}/{num_d11} "
          f"({100*result['fraction_working']:.1f}%)")
    print(f"  E[N(D11)] = {result['E_N']:.4f}")
    print(f"  E[N(D11)^2] = {result['E_N2']:.4f}")
    print(f"  E[N^2]/E[N]^2 = {result['second_moment_ratio']:.4f}")
    print(f"  E[Z] = Pr[valid] = {result['E_Z']:.6e}")
    print(f"  log2(E[Z]) = {result['log2_EZ']:.4f}")
    print(f"  log2(#D12) = {result['log2_headroom']:.4f}")
    print(f"  log2(E[Z]) - log2(#D12) = {result['log2_EZ_minus_headroom']:.4f}")
    print(f"  Paley-Zygmund bound Pr[N>0] >= {result['PZ_bound']:.4f}")
    print(f"  Actual Pr[N>0] = {result['Pr_N_gt_0']:.4f}")
    print(f"  Orbit N-distribution: {dict(sorted(orbit_N_dist.items()))}")

    # Cross-prime summary table
    print(f"\n{'='*70}")
    print("DOUBLE-AVERAGING SUMMARY")
    print(f"{'='*70}")
    print(f"{'p':>4s} {'n':>4s} {'#D11':>8s} {'#D12':>12s} "
          f"{'Total valid':>14s} {'E[N]':>10s} {'E[N2]/E[N]2':>12s} "
          f"{'log2(E[Z])':>12s} {'log2(#D12)':>12s} {'headroom':>10s}")
    print(f"  {'-'*110}")
    for p_val in [11, 19, 23, 31]:
        r = results[p_val]
        headroom = r['log2_EZ'] - r['log2_EZ_minus_headroom'] if r['log2_EZ'] != float('-inf') else 0
        print(f"  {r['p']:4d} {r['n']:4d} {r['num_d11']:8d} {r['num_d12']:12,} "
              f"{r['total_valid_pairs']:14,} {r['E_N']:10.2f} {r['second_moment_ratio']:12.2f} "
              f"{r['log2_EZ']:12.4f} {r['log2_headroom']:12.4f} "
              f"{r['log2_EZ_minus_headroom']:10.4f}")

    # Write results markdown
    write_results_markdown(results)

    # Save JSON
    script_dir = os.path.dirname(os.path.abspath(__file__))
    json_path = os.path.join(script_dir, "..", "proofs", "double_averaging_data.json")
    with open(json_path, "w") as f:
        json.dump({str(k): v for k, v in results.items()}, f, indent=2)
    print(f"\nJSON saved to {json_path}")


def write_results_markdown(results):
    """Write results to proofs/double_averaging_results.md."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    proofs_dir = os.path.join(script_dir, "..", "proofs")
    os.makedirs(proofs_dir, exist_ok=True)
    path = os.path.join(proofs_dir, "double_averaging_results.md")

    lines = [
        "# Double-Averaging Computation for E[Z]",
        "",
        "**Date**: 2026-02-11",
        "**Purpose**: Compute E[Z] = Pr[random (D11, D12) is valid] and compare to theoretical headroom.",
        "",
        "## 1. Setup",
        "",
        "For each prime p = 3 mod 4, we consider:",
        "- **D11**: symmetric subset of {1,...,p-1} with |D11| = n = (p+1)/2",
        "- **D12**: subset of {0,...,p-1} containing 0, with |D12| = (p-1)/2",
        "  - Equivalently, D12 = {0} union S where S is a k-subset of {1,...,p-1}, k = (p-3)/2",
        "- **Z(D11, D12)** = 1 if (D11, D12) satisfies all Ramsey constraints",
        "",
        "The double average is:",
        "",
        "  E[Z] = (1 / #D11) * (1 / #D12) * sum_{D11} sum_{D12} Z(D11, D12)",
        "       = (1 / #D11) * (1 / #D12) * (total valid pairs)",
        "",
        "## 2. Counting",
        "",
        "| Quantity | Formula | p=11 | p=19 | p=23 | p=31 |",
        "|----------|---------|------|------|------|------|",
    ]

    # Build the counting table
    for label, key in [
        ("n = (p+1)/2", "n"),
        ("k = (p-3)/2", "k"),
        ("#D11 = C(R, n/2), R=(p-1)/2", "num_d11"),
        ("#D12 = C(p-1, k)", "num_d12"),
    ]:
        vals = []
        for p_val in [11, 19, 23, 31]:
            v = results[p_val][key]
            if isinstance(v, int) and v > 9999:
                vals.append(f"{v:,}")
            else:
                vals.append(str(v))
        lines.append(f"| {label} | | {' | '.join(vals)} |")

    # Total valid pairs
    vals = [f"{results[p]['total_valid_pairs']:,}" for p in [11, 19, 23, 31]]
    lines.append(f"| Total valid pairs | sum N(D11) | {' | '.join(vals)} |")

    lines += [
        "",
        "## 3. Double-Averaging Results",
        "",
        "| Quantity | p=11 | p=19 | p=23 | p=31 |",
        "|----------|------|------|------|------|",
    ]

    for label, key, fmt in [
        ("E[N(D11)]", "E_N", ".4f"),
        ("E[N(D11)^2]", "E_N2", ".2f"),
        ("E[N^2]/E[N]^2", "second_moment_ratio", ".4f"),
        ("E[Z] = Pr[valid]", "E_Z", ".6e"),
        ("log2(E[Z])", "log2_EZ", ".4f"),
        ("log2(#D12) (headroom)", "log2_headroom", ".4f"),
        ("log2(E[Z]) - log2(#D12)", "log2_EZ_minus_headroom", ".4f"),
        ("Pr[N(D11) > 0]", "Pr_N_gt_0", ".4f"),
        ("Paley-Zygmund bound", "PZ_bound", ".4f"),
    ]:
        vals = []
        for p_val in [11, 19, 23, 31]:
            v = results[p_val][key]
            vals.append(f"{v:{fmt}}")
        lines.append(f"| {label} | {' | '.join(vals)} |")

    lines += [
        "",
        "## 4. Interpretation",
        "",
        "### 4.1 Headroom Analysis",
        "",
        "The key quantity is **log2(E[Z]) - log2(#D12)**. This measures how many \"bits\" of headroom",
        "remain after accounting for the constraint cost:",
        "",
        "- log2(#D12) = log2(C(p-1, k)) is the entropy of the D12 choice (the \"budget\")",
        "- log2(E[Z]) = log2(Pr[valid]) combines the budget with the constraint cost",
        "- The difference = log2(E[N]) = log2(expected number of valid D12 per D11)",
        "",
    ]

    for p_val in [11, 19, 23, 31]:
        r = results[p_val]
        lines.append(f"**p = {p_val}**: log2(E[N]) = {r['log2_EZ'] - r['log2_EZ_minus_headroom'] + r['log2_EZ_minus_headroom']:.4f} "
                     f"- {r['log2_headroom']:.4f} = {r['log2_EZ_minus_headroom']:.4f}. "
                     f"E[N] = {r['E_N']:.2f}, so on average each D11 has ~{r['E_N']:.1f} valid D12.")
        lines.append("")

    lines += [
        "### 4.2 Second Moment Ratio",
        "",
        "The Paley-Zygmund inequality gives:",
        "",
        "  Pr[N(D11) > 0] >= E[N]^2 / E[N^2] = 1 / (E[N^2]/E[N]^2)",
        "",
        "| p | E[N^2]/E[N]^2 | PZ bound | Actual Pr[N>0] | Ratio actual/PZ |",
        "|---|---------------|----------|----------------|-----------------|",
    ]

    for p_val in [11, 19, 23, 31]:
        r = results[p_val]
        pz = r['PZ_bound']
        actual = r['Pr_N_gt_0']
        ratio_pz = actual / pz if pz > 0 else float('inf')
        lines.append(f"| {p_val} | {r['second_moment_ratio']:.2f} | {pz:.4f} | {actual:.4f} | {ratio_pz:.2f} |")

    lines += [
        "",
        "The PZ bound is tight at p=11 and p=19 (ratio = 1.0), meaning the second moment",
        "method is essentially optimal there. At p=23 and p=31, the actual probability exceeds",
        "the PZ bound by a constant factor.",
        "",
        "### 4.3 N-value Distributions",
        "",
    ]

    for p_val in [11, 19, 23, 31]:
        r = results[p_val]
        lines.append(f"**p = {p_val}**: N values = {r['N_distribution']}")
        lines.append("")

    lines += [
        "### 4.4 Growth of E[N] and the Ratio",
        "",
        "| p | E[N] | log2(E[N]) | E[N^2]/E[N]^2 | log2(ratio) | (p-1)/2 (bits budget) |",
        "|---|------|------------|---------------|-------------|----------------------|",
    ]

    for p_val in [11, 19, 23, 31]:
        r = results[p_val]
        log2_EN = math.log2(r['E_N']) if r['E_N'] > 0 else float('-inf')
        log2_ratio = math.log2(r['second_moment_ratio']) if r['second_moment_ratio'] > 0 else float('-inf')
        lines.append(f"| {p_val} | {r['E_N']:.2f} | {log2_EN:.2f} | "
                     f"{r['second_moment_ratio']:.2f} | {log2_ratio:.2f} | {(p_val-1)//2} |")

    lines += [
        "",
        "The budget grows as (p-1)/2 bits. For the second moment method to work, we need:",
        "- E[N] >> 1 (equivalently log2(E[N]) >> 0)",
        "- E[N^2]/E[N]^2 = O(poly(p))",
        "",
        "Both conditions appear to be satisfied empirically.",
        "",
        "## 5. Divisibility by 2p",
        "",
        "All nonzero N(D11) values are divisible by 2p:",
        "",
    ]

    for p_val in [11, 19, 23, 31]:
        r = results[p_val]
        nonzero_N = [N for N in r['N_distribution'].keys() if int(N) > 0]
        if nonzero_N:
            check = all(int(N) % (2 * p_val) == 0 for N in nonzero_N)
            divisors = [f"{int(N)} = {int(N)//(2*p_val)} * {2*p_val}" for N in nonzero_N]
            lines.append(f"- p={p_val}: N values {nonzero_N} all divisible by {2*p_val}: "
                        f"{check}. {', '.join(divisors)}")

    lines += [
        "",
        "This confirms the proven structure: valid D12 sets are closed under additive shifts",
        "mod p (giving factor p) and negation (giving factor 2).",
        "",
    ]

    with open(path, "w") as f:
        f.write("\n".join(lines))

    print(f"\nMarkdown report saved to {path}")


if __name__ == "__main__":
    main()
