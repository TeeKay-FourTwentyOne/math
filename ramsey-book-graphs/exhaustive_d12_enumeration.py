#!/usr/bin/env python3
"""Exhaustive D12 enumeration for L6 verification at small primes.

For each prime p in {11, 19, 23}, with m = p, n = (p+1)/2:
1. Enumerate all symmetric D11 subsets of {1,...,p-1} with |D11| = n.
2. For each D11, enumerate ALL D12 = {0} ∪ S where S is a k-subset of {1,...,p-1},
   k = (p-3)/2.
3. Check validity: A(d)+B(d) <= (p-3)/2 for all d in D11 (binding),
                   A(d)+B(d) <= (p+3)/2 for all d in D22 (loose).
4. Report per-D11 valid counts, minimum across all D11, second moment ratio.

Optimized: precompute all D12 autocorrelations via batch FFT, then vectorized
constraint checking per D11.
"""

import numpy as np
import json
import time
import os
import sys
from itertools import combinations
from math import comb
from fractions import Fraction
from collections import defaultdict


def get_symmetric_pairs(p):
    """Get all negation pairs {d, p-d} from {1,...,p-1}."""
    pairs = []
    seen = set()
    for d in range(1, p):
        if d not in seen:
            pairs.append((d, (-d) % p))
            seen.add(d)
            seen.add((-d) % p)
    return pairs


def enumerate_symmetric_D11(p, size):
    """Enumerate all symmetric subsets of {1,...,p-1} with given even size.

    Elements come in pairs {d, p-d}. We choose size//2 pairs.
    """
    pairs = get_symmetric_pairs(p)
    num_pairs = size // 2
    for chosen_pairs in combinations(range(len(pairs)), num_pairs):
        D11 = set()
        for i in chosen_pairs:
            d, neg_d = pairs[i]
            D11.add(d)
            D11.add(neg_d)
        yield frozenset(D11)


def batch_autocorrelation(indicator_matrix, p):
    """Compute Delta(S,S,d) for all sets and all d using FFT.

    indicator_matrix: (num_sets, p) where each row is an indicator vector.
    Returns: (num_sets, p) integer matrix where [i,d] = Delta(S_i, S_i, d).
    """
    fft_vals = np.fft.fft(indicator_matrix, axis=1)
    autocorr = np.fft.ifft(np.abs(fft_vals) ** 2, axis=1).real
    return np.round(autocorr).astype(np.int32)


def autocorrelation_single(D11_set, p):
    """Compute A(d) = Delta(D11,D11,d) for a single set via FFT."""
    indicator = np.zeros(p, dtype=np.float64)
    for j in D11_set:
        indicator[j] = 1.0
    fft_val = np.fft.fft(indicator)
    autocorr = np.fft.ifft(np.abs(fft_val) ** 2).real
    return np.round(autocorr).astype(np.int32)


def build_d12_candidates(p, k):
    """Build indicator matrix for all D12 = {0} ∪ S, S a k-subset of {1,...,p-1}.

    Returns indicator_matrix of shape (C(p-1, k), p).
    """
    candidates = list(range(1, p))
    num_d12 = comb(p - 1, k)
    matrix = np.zeros((num_d12, p), dtype=np.float64)
    matrix[:, 0] = 1.0  # 0 is always in D12

    for i, chosen in enumerate(combinations(candidates, k)):
        for j in chosen:
            matrix[i, j] = 1.0

    return matrix


def run_prime(p, verbose=True):
    """Run exhaustive enumeration for prime p.

    Returns a results dict with per-D11 statistics.
    """
    n = (p + 1) // 2
    d11_size = n          # |D11| = (p+1)/2
    d12_size = n - 1      # |D12| = (p-1)/2
    k = d12_size - 1      # |S| = (p-3)/2 (non-zero elements of D12)
    threshold_binding = (p - 3) // 2   # A(d)+B(d) <= this for d in D11
    threshold_loose = (p + 3) // 2     # A(d)+B(d) <= this for d in D22

    pairs = get_symmetric_pairs(p)
    num_pairs = len(pairs)
    num_d11 = comb(num_pairs, d11_size // 2)
    num_d12 = comb(p - 1, k)
    total = num_d11 * num_d12

    if verbose:
        print(f"\n{'=' * 80}")
        print(f"PRIME p = {p}, n = {n}, m = {p}")
        print(f"|D11| = {d11_size}, |D12| = {d12_size}, |S| = k = {k}")
        print(f"Binding: A(d)+B(d) <= {threshold_binding} for d in D11")
        print(f"Loose:   A(d)+B(d) <= {threshold_loose} for d in D22")
        print(f"Search space: {num_d11} D11 x {num_d12} D12 = {total:,}")
        print(f"{'=' * 80}")

    t0 = time.time()

    # Step 1: Precompute all D12 autocorrelations via batch FFT
    if verbose:
        print(f"\nPrecomputing {num_d12:,} D12 autocorrelations via batch FFT...")
    d12_matrix = build_d12_candidates(p, k)
    B_matrix = batch_autocorrelation(d12_matrix, p)  # (num_d12, p)
    t1 = time.time()
    if verbose:
        print(f"  Done in {t1 - t0:.1f}s. B_matrix shape: {B_matrix.shape}")

    # Step 2: Enumerate symmetric D11, check constraints vectorized
    d11_list = list(enumerate_symmetric_D11(p, d11_size))
    assert len(d11_list) == num_d11, f"Expected {num_d11} D11, got {len(d11_list)}"

    results_per_d11 = []

    for d11_idx, D11 in enumerate(d11_list):
        # Compute A(d) for this D11
        A = autocorrelation_single(D11, p)

        # Build masks for D11 and D22 positions (excluding d=0)
        d11_indices = np.array(sorted(D11), dtype=np.int32)
        d22_set = set(range(1, p)) - D11
        d22_indices = np.array(sorted(d22_set), dtype=np.int32)

        # F[i, d] = A(d) + B_i(d)
        F_binding = A[d11_indices] + B_matrix[:, d11_indices]  # (num_d12, |D11|)
        valid_binding = np.all(F_binding <= threshold_binding, axis=1)

        F_loose = A[d22_indices] + B_matrix[:, d22_indices]    # (num_d12, |D22|)
        valid_loose = np.all(F_loose <= threshold_loose, axis=1)

        valid_mask = valid_binding & valid_loose
        count = int(valid_mask.sum())

        # Per-D11 statistics
        entry = {
            "D11": sorted(D11),
            "A_values": {int(d): int(A[d]) for d in sorted(D11)},
            "num_valid_d12": count,
            "fraction_valid": count / num_d12,
        }

        # Per-constraint marginal counts (how many D12 satisfy each constraint independently)
        marginal_counts = {}
        for j, d in enumerate(d11_indices):
            marginal_counts[int(d)] = int(np.sum(
                A[d] + B_matrix[:, d] <= threshold_binding
            ))
        entry["marginal_counts"] = marginal_counts

        # Product of marginals / num_d12^{|D11|-1}
        # joint_prob = count / num_d12
        # product_of_marginals = prod(marginal_counts[d]/num_d12 for d in D11)
        prod_marginals = 1.0
        for d in sorted(D11):
            prod_marginals *= marginal_counts[int(d)] / num_d12
        entry["product_of_marginals"] = prod_marginals
        entry["joint_prob"] = count / num_d12
        if prod_marginals > 0:
            entry["ratio_joint_over_product"] = (count / num_d12) / prod_marginals
        else:
            entry["ratio_joint_over_product"] = float('inf') if count > 0 else float('nan')

        # S1 statistics: for valid D12, compute S1 = sum_{d in D11} B(d)
        if count > 0:
            valid_indices = np.where(valid_mask)[0]
            S1_values = B_matrix[valid_indices][:, d11_indices].sum(axis=1)
            entry["S1_stats"] = {
                "mean": float(S1_values.mean()),
                "std": float(S1_values.std()),
                "min": int(S1_values.min()),
                "max": int(S1_values.max()),
            }

        # For p=11: pairwise overlap |D12 ∩ D12'| among valid D12 sets
        if p <= 19 and count > 0 and count <= 5000:
            valid_indices = np.where(valid_mask)[0]
            valid_d12_rows = d12_matrix[valid_indices]  # (count, p) float indicators
            # Pairwise overlap = dot product of indicator vectors
            overlap_matrix = valid_d12_rows @ valid_d12_rows.T
            # Extract upper triangle (excluding diagonal)
            iu = np.triu_indices(count, k=1)
            overlaps = overlap_matrix[iu].astype(int)
            overlap_dist = defaultdict(int)
            for ov in overlaps:
                overlap_dist[int(ov)] += 1
            entry["pairwise_overlap_distribution"] = dict(sorted(overlap_dist.items()))

        results_per_d11.append(entry)

        if verbose and (d11_idx + 1) % max(1, num_d11 // 10) == 0:
            elapsed = time.time() - t0
            print(f"  D11 {d11_idx + 1}/{num_d11}: "
                  f"valid={count}, frac={count/num_d12:.4f}, "
                  f"ratio={entry.get('ratio_joint_over_product', 0):.3f}, "
                  f"{elapsed:.1f}s")

    elapsed = time.time() - t0

    # Aggregate statistics
    valid_counts = [r["num_valid_d12"] for r in results_per_d11]
    min_count = min(valid_counts)
    max_count = max(valid_counts)
    all_have_valid = all(c > 0 for c in valid_counts)
    ratios = [r["ratio_joint_over_product"] for r in results_per_d11
              if r["ratio_joint_over_product"] != float('inf')
              and not (isinstance(r["ratio_joint_over_product"], float)
                       and r["ratio_joint_over_product"] != r["ratio_joint_over_product"])]

    summary = {
        "p": p,
        "n": n,
        "d11_size": d11_size,
        "d12_size": d12_size,
        "k": k,
        "threshold_binding": threshold_binding,
        "threshold_loose": threshold_loose,
        "num_d11": num_d11,
        "num_d12": num_d12,
        "total_pairs_checked": total,
        "min_valid_d12": min_count,
        "max_valid_d12": max_count,
        "all_d11_have_valid_d12": all_have_valid,
        "elapsed_seconds": elapsed,
        "ratio_stats": {
            "min": min(ratios) if ratios else None,
            "max": max(ratios) if ratios else None,
            "mean": sum(ratios) / len(ratios) if ratios else None,
        },
    }

    if verbose:
        print(f"\n{'=' * 80}")
        print(f"RESULTS FOR p = {p}")
        print(f"{'=' * 80}")
        print(f"Total D11: {num_d11}")
        print(f"Total D12: {num_d12}")
        print(f"Total pairs checked: {total:,}")
        print(f"Elapsed: {elapsed:.1f}s")
        print(f"")
        print(f"Min valid D12 across all D11: {min_count}")
        print(f"Max valid D12 across all D11: {max_count}")
        print(f"ALL D11 have at least 1 valid D12: {all_have_valid}")
        print(f"")
        if ratios:
            print(f"Ratio (joint/product of marginals):")
            print(f"  min = {min(ratios):.4f}")
            print(f"  max = {max(ratios):.4f}")
            print(f"  mean = {sum(ratios)/len(ratios):.4f}")

        print(f"\nPer-D11 details:")
        print(f"{'D11':>40s} {'#valid':>8s} {'frac':>8s} {'ratio':>8s}")
        print(f"  {'-'*70}")
        for r in results_per_d11:
            d11_str = str(r["D11"])
            if len(d11_str) > 38:
                d11_str = d11_str[:35] + "..."
            ratio_str = f"{r['ratio_joint_over_product']:.4f}" if isinstance(
                r['ratio_joint_over_product'], (int, float)) and r[
                    'ratio_joint_over_product'] == r['ratio_joint_over_product'] else "N/A"
            print(f"  {d11_str:>38s} {r['num_valid_d12']:8d} "
                  f"{r['fraction_valid']:8.4f} {ratio_str:>8s}")

    return summary, results_per_d11


def main():
    print("=" * 80)
    print("EXHAUSTIVE D12 ENUMERATION FOR L6 VERIFICATION")
    print("=" * 80)

    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.join(script_dir, "..", "enumeration_data")
    os.makedirs(output_dir, exist_ok=True)
    proofs_dir = os.path.join(script_dir, "..", "proofs")

    all_summaries = {}
    all_details = {}

    for p in [11, 19, 23]:
        summary, details = run_prime(p)
        all_summaries[str(p)] = summary
        all_details[str(p)] = details

        # Save per-prime raw data
        outpath = os.path.join(output_dir, f"enumeration_p{p}.json")
        with open(outpath, "w") as f:
            json.dump({
                "summary": summary,
                "per_d11": details,
            }, f, indent=2, default=str)
        print(f"\nRaw data saved to {outpath}")

    # Cross-prime summary
    print(f"\n{'=' * 80}")
    print("CROSS-PRIME SUMMARY")
    print(f"{'=' * 80}")
    print(f"{'p':>4s} {'n':>4s} {'#D11':>6s} {'#D12':>10s} "
          f"{'min_valid':>10s} {'max_valid':>10s} {'all_ok':>7s} "
          f"{'min_ratio':>10s} {'time':>8s}")
    print(f"  {'-'*78}")
    for p_str, s in all_summaries.items():
        min_r = s["ratio_stats"]["min"]
        min_r_str = f"{min_r:.4f}" if min_r is not None else "N/A"
        print(f"  {s['p']:4d} {s['n']:4d} {s['num_d11']:6d} {s['num_d12']:10d} "
              f"{s['min_valid_d12']:10d} {s['max_valid_d12']:10d} "
              f"{'YES' if s['all_d11_have_valid_d12'] else 'NO':>7s} "
              f"{min_r_str:>10s} {s['elapsed_seconds']:8.1f}s")

    # Write summary markdown
    md_path = os.path.join(proofs_dir, "enumeration_results.md")
    write_markdown_report(md_path, all_summaries, all_details)
    print(f"\nMarkdown report saved to {md_path}")


def write_markdown_report(path, summaries, all_details):
    """Write a markdown summary of the enumeration results."""
    lines = [
        "# Exhaustive D12 Enumeration Results",
        "",
        f"**Date**: {time.strftime('%Y-%m-%d')}",
        "**Purpose**: Verify that every symmetric D11 admits a valid D12 for small primes.",
        "",
        "## Summary",
        "",
        "| p | n | #D11 | #D12 | Min valid | Max valid | All D11 ok? | Min ratio | Time |",
        "|---|---|------|------|-----------|-----------|-------------|-----------|------|",
    ]

    for p_str, s in summaries.items():
        min_r = s["ratio_stats"]["min"]
        min_r_str = f"{min_r:.4f}" if min_r is not None else "N/A"
        lines.append(
            f"| {s['p']} | {s['n']} | {s['num_d11']} | {s['num_d12']:,} "
            f"| {s['min_valid_d12']:,} | {s['max_valid_d12']:,} "
            f"| {'YES' if s['all_d11_have_valid_d12'] else 'NO'} "
            f"| {min_r_str} | {s['elapsed_seconds']:.1f}s |"
        )

    lines += [
        "",
        "## Interpretation",
        "",
        "- **Min valid**: The minimum number of valid D12 sets across all symmetric D11.",
        "- **All D11 ok?**: YES means every symmetric D11 has at least one valid D12.",
        "- **Ratio**: joint probability / product of marginals. Ratio >= 1 means "
        "positive association; ratio < 1 but > 0 means the joint event is less likely "
        "than independence would predict, but still occurs.",
        "",
        "If 'All D11 ok' = YES for all primes, this confirms the existence claim in the "
        "proof of R(B_{n-1}, B_n) = 4n-1 for these values of n.",
        "",
        "## Constraints",
        "",
        "For each (D11, D12) pair with m = p:",
        "- **Binding** (d in D11): A(d) + B(d) <= (p-3)/2",
        "- **Loose** (d in D22):   A(d) + B(d) <= (p+3)/2",
        "",
        "where A(d) = Delta(D11,D11,d), B(d) = Delta(D12,D12,d).",
        "",
    ]

    # Per-prime details
    for p_str, details in all_details.items():
        p = int(p_str)
        s = summaries[p_str]
        lines += [
            f"## p = {p} (n = {s['n']})",
            "",
            f"- {s['num_d11']} symmetric D11, each checked against {s['num_d12']:,} D12 candidates",
            f"- Binding threshold: A(d)+B(d) <= {s['threshold_binding']}",
            f"- Loose threshold: A(d)+B(d) <= {s['threshold_loose']}",
            "",
            "| D11 | #Valid D12 | Fraction | Ratio (joint/prod) |",
            "|-----|------------|----------|---------------------|",
        ]
        for r in details:
            d11_str = str(r["D11"])
            ratio = r["ratio_joint_over_product"]
            ratio_str = f"{ratio:.4f}" if isinstance(ratio, (int, float)) and ratio == ratio else "N/A"
            lines.append(
                f"| {d11_str} | {r['num_valid_d12']:,} | {r['fraction_valid']:.4f} | {ratio_str} |"
            )
        lines.append("")

    with open(path, "w") as f:
        f.write("\n".join(lines))


if __name__ == "__main__":
    main()
