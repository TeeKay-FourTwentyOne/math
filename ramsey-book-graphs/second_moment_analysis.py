#!/usr/bin/env python3
"""Second moment analysis: E[N^2] for valid D12 count.

For each working D11 at p=11,19,23:
- Enumerate all valid D12
- Compute E[N^2] = |{(D12, D12') : both valid}| (ordered pairs)
- Compute E[N]^2 = (number of valid D12)^2
- Ratio E[N^2] / E[N]^2  (should be O(poly(p)) for Paley-Zygmund)
- Breakdown by overlap |D12 ∩ D12'|

Paley-Zygmund: Pr[N > 0] >= E[N]^2 / E[N^2]
If ratio = E[N^2]/E[N]^2 = O(poly(p)), then Pr[N>0] = Omega(1/poly(p)) > 0.
"""

import numpy as np
import json
import os
import sys
from itertools import combinations
from math import comb, log2
from collections import defaultdict, Counter


def autocorrelation_full(D11, p):
    """Compute A(d) for all d."""
    indicator = np.zeros(p, dtype=np.float64)
    for j in D11:
        indicator[j] = 1.0
    fft_val = np.fft.fft(indicator)
    autocorr = np.fft.ifft(np.abs(fft_val) ** 2).real
    return np.round(autocorr).astype(int)


def autocorrelation_single_set(S, p):
    """Compute B(d) for a single set S."""
    indicator = np.zeros(p, dtype=np.float64)
    for j in S:
        indicator[j] = 1.0
    fft_val = np.fft.fft(indicator)
    autocorr = np.fft.ifft(np.abs(fft_val) ** 2).real
    return np.round(autocorr).astype(int)


def is_valid_d12(D12, A, D11, D22, threshold_binding, threshold_loose, p):
    """Check if D12 is valid given D11's autocorrelation A."""
    B = autocorrelation_single_set(D12, p)
    for d in D11:
        if A[d] + B[d] > threshold_binding:
            return False
    for d in D22:
        if A[d] + B[d] > threshold_loose:
            return False
    return True


def enumerate_valid_d12(D11, p):
    """Enumerate all valid D12 for a given D11."""
    n = (p + 1) // 2
    k = (p - 3) // 2
    threshold_binding = (p - 3) // 2
    threshold_loose = (p + 3) // 2
    D22 = set(range(1, p)) - D11

    A = autocorrelation_full(D11, p)

    # Batch approach: generate all D12 indicator matrices, FFT, check
    candidates = list(range(1, p))
    valid_d12_list = []

    # For small p, enumerate all D12
    batch = []
    batch_sets = []
    for combo in combinations(candidates, k):
        D12 = frozenset({0} | set(combo))
        batch_sets.append(D12)
        ind = np.zeros(p, dtype=np.float64)
        for j in D12:
            ind[j] = 1.0
        batch.append(ind)

        if len(batch) >= 10000:
            # Process batch
            mat = np.array(batch)
            fft_vals = np.fft.fft(mat, axis=1)
            B_batch = np.round(np.fft.ifft(np.abs(fft_vals) ** 2, axis=1).real).astype(int)

            d11_arr = np.array(sorted(D11), dtype=int)
            d22_arr = np.array(sorted(D22), dtype=int)

            for i in range(len(batch)):
                B = B_batch[i]
                ok = True
                for d in d11_arr:
                    if A[d] + B[d] > threshold_binding:
                        ok = False
                        break
                if ok:
                    for d in d22_arr:
                        if A[d] + B[d] > threshold_loose:
                            ok = False
                            break
                if ok:
                    valid_d12_list.append(batch_sets[i])

            batch = []
            batch_sets = []

    # Process remaining
    if batch:
        mat = np.array(batch)
        fft_vals = np.fft.fft(mat, axis=1)
        B_batch = np.round(np.fft.ifft(np.abs(fft_vals) ** 2, axis=1).real).astype(int)

        d11_arr = np.array(sorted(D11), dtype=int)
        d22_arr = np.array(sorted(D22), dtype=int)

        for i in range(len(batch)):
            B = B_batch[i]
            ok = True
            for d in d11_arr:
                if A[d] + B[d] > threshold_binding:
                    ok = False
                    break
            if ok:
                for d in d22_arr:
                    if A[d] + B[d] > threshold_loose:
                        ok = False
                        break
            if ok:
                valid_d12_list.append(batch_sets[i])

    return valid_d12_list


def compute_second_moment(valid_d12_list, p):
    """Compute E[N^2] and overlap distribution.

    E[N^2] = sum over all ordered pairs (D12, D12') where both are valid
           = N^2  (since we have all valid D12, E[N^2] = N^2 if deterministic)

    But the interesting quantity is: over random (D11, D12), what is E[N^2]?
    Here N = #{valid D12 for this D11}, and we compute E[N^2] = N^2.

    For the Paley-Zygmund argument over random D12:
    N(D12) = 1 if D12 is valid, 0 otherwise.
    E[N] = |valid| / |total D12|
    E[N^2] = E[N] (since N^2 = N for indicator)

    Actually, for Paley-Zygmund we want:
    N = |{valid D12}| (count, not indicator)
    E_D11[N] = average over D11 of |valid D12 for D11|
    E_D11[N^2] = average over D11 of |valid D12 for D11|^2

    But the task says "for each pair of valid D12, record overlap".
    This computes the CONDITIONAL second moment: given D11 is good,
    what is E_{D12,D12'}[1(both valid)] where D12, D12' are independent?

    Let me compute:
    - For a FIXED D11, N = |valid D12|
    - E[N^2] comes from summing over all ordered pairs of D12 that are valid
    - The pairwise overlap |D12 ∩ D12'| tells us the correlation structure
    """
    N = len(valid_d12_list)
    if N == 0:
        return {"N": 0, "N_sq": 0, "ratio": 0}

    # Compute pairwise overlaps
    overlap_counts = Counter()
    for i in range(N):
        for j in range(N):
            s = len(valid_d12_list[i] & valid_d12_list[j])
            overlap_counts[s] += 1

    # E[N^2] = sum over all valid pairs = N^2 (trivially, but let's verify)
    total_pairs = sum(overlap_counts.values())
    assert total_pairs == N * N, f"Expected {N*N}, got {total_pairs}"

    return {
        "N": N,
        "N_sq": N * N,
        "overlap_distribution": dict(sorted(overlap_counts.items())),
    }


def compute_random_overlap_distribution(p, k, num_samples=100000, rng=None):
    """Compute the overlap distribution for two random D12 = {0} ∪ S.

    This is the baseline: what overlap do two random D12 have?
    """
    if rng is None:
        rng = np.random.default_rng(42)

    candidates = np.arange(1, p)
    overlap_counts = Counter()

    for _ in range(num_samples):
        S1 = set(rng.choice(candidates, size=k, replace=False))
        S2 = set(rng.choice(candidates, size=k, replace=False))
        D12_1 = frozenset({0} | S1)
        D12_2 = frozenset({0} | S2)
        s = len(D12_1 & D12_2)
        overlap_counts[s] += 1

    return {s: c / num_samples for s, c in sorted(overlap_counts.items())}


def analyze_prime(p):
    """Full second moment analysis for one prime."""
    n = (p + 1) // 2
    k = (p - 3) // 2
    total_d12 = comb(p - 1, k)

    print(f"\n{'='*80}")
    print(f"p = {p}, n = {n}, |D12| = {n-1}, k = {k}")
    print(f"Total D12: {total_d12}")
    print(f"{'='*80}")
    sys.stdout.flush()

    # Load enumeration data to find working D11
    script_dir = os.path.dirname(os.path.abspath(__file__))
    data_path = os.path.join(script_dir, "..", "enumeration_data", f"enumeration_p{p}.json")
    with open(data_path) as f:
        data = json.load(f)

    working_d11 = [entry for entry in data["per_d11"] if entry["num_valid_d12"] > 0]
    print(f"Working D11: {len(working_d11)}")

    # Group by num_valid (since many D11 have the same count)
    groups = defaultdict(list)
    for entry in working_d11:
        groups[entry["num_valid_d12"]].append(entry)

    results_by_group = {}
    all_overlap_dist = Counter()
    total_N = 0
    total_N_sq = 0
    num_d11_analyzed = 0

    for nv in sorted(groups.keys(), reverse=True):
        entries = groups[nv]
        # Analyze one representative from each group
        entry = entries[0]
        D11 = set(entry["D11"])
        max_A = max(entry["A_values"].values())

        print(f"\n--- D11 group: {nv} valid D12, {len(entries)} D11 with this count, max_A(D11)={max_A} ---")
        sys.stdout.flush()

        valid_d12 = enumerate_valid_d12(D11, p)
        assert len(valid_d12) == nv, f"Expected {nv} valid, got {len(valid_d12)}"

        sm = compute_second_moment(valid_d12, p)

        print(f"  N = {sm['N']}")
        print(f"  N^2 = {sm['N_sq']}")
        print(f"  Overlap distribution:")
        for s, count in sorted(sm["overlap_distribution"].items()):
            pct = 100 * count / sm["N_sq"]
            print(f"    |D12 ∩ D12'| = {s}: {count} pairs ({pct:.1f}%)")

        # For ALL D11 in this group
        group_N = nv * len(entries)
        group_N_sq = nv * nv * len(entries)  # N^2 for each D11

        total_N += group_N
        total_N_sq += group_N_sq
        num_d11_analyzed += len(entries)

        # Accumulate overlap distribution (weighted by group size)
        for s, count in sm["overlap_distribution"].items():
            all_overlap_dist[s] += count * len(entries)

        results_by_group[nv] = {
            "num_d11": len(entries),
            "N": nv,
            "overlap_distribution": sm["overlap_distribution"],
        }

    # Now compute the Paley-Zygmund quantities
    # Over random D11 (uniform over all symmetric D11):
    total_d11 = data["summary"]["num_d11"]
    E_N = total_N / total_d11  # average valid D12 per D11
    E_N_sq = total_N_sq / total_d11  # average N^2 per D11
    ratio = E_N_sq / (E_N * E_N) if E_N > 0 else float('inf')

    print(f"\n{'='*60}")
    print(f"PALEY-ZYGMUND ANALYSIS for p = {p}")
    print(f"{'='*60}")
    print(f"Total D11: {total_d11}")
    print(f"Working D11: {num_d11_analyzed}")
    print(f"Sum of N over all D11: {total_N}")
    print(f"Sum of N^2 over all D11: {total_N_sq}")
    print(f"E[N] = {E_N:.4f}")
    print(f"E[N^2] = {E_N_sq:.4f}")
    print(f"E[N]^2 = {E_N*E_N:.4f}")
    print(f"E[N^2] / E[N]^2 = {ratio:.4f}")
    print(f"Pr[N > 0] >= E[N]^2 / E[N^2] = {1/ratio:.6f}" if ratio > 0 and ratio < float('inf') else "")

    # Now the more interesting analysis: for a FIXED good D11, analyze over random D12
    # Here N is NOT over D11, but the count of valid D12 for a fixed D11.
    # For Paley-Zygmund over (D11, D12) jointly:
    # Let X(D11, D12) = 1{D12 valid for D11}
    # N(D11) = sum_{D12} X(D11, D12)
    # E_{D11}[N] = sum_D11 N(D11) / |all D11| = total_N / total_d11
    # E_{D11}[N^2] = sum_D11 N(D11)^2 / |all D11| = total_N_sq / total_d11
    # By PZ: Pr[N > 0] >= E[N]^2 / E[N^2]
    # Note: E[N]^2 = (total_N / total_d11)^2

    print(f"\nPaley-Zygmund bound over RANDOM D11:")
    print(f"  Pr[N(D11) > 0] >= E[N]^2 / E[N^2] = {E_N*E_N / E_N_sq:.6f}")
    print(f"  (This is the fraction of D11 with at least one valid D12)")
    print(f"  Actual fraction: {num_d11_analyzed}/{total_d11} = {num_d11_analyzed/total_d11:.6f}")
    print(f"  PZ bound is {'tight' if abs(E_N*E_N/E_N_sq - num_d11_analyzed/total_d11) < 0.01 else 'not tight'}")

    # For a FIXED good D11, the second moment over D12 pairs:
    # Let Y(D12) = 1{D12 valid for this D11}
    # N = sum_{D12} Y(D12)
    # E[N] = N / total_d12  (fraction of D12 that are valid)
    # E[N^2] = sum_{D12, D12'} Pr[both valid]
    #        = N^2 / total_d12^2   (since both D12 are chosen independently)
    #
    # Wait, this is just (E[N])^2 if D12 is drawn randomly and independently.
    # The second moment is trivially N^2 for ordered pairs of valid D12.
    # But E[Y(D12) * Y(D12')] = Pr[D12 valid AND D12' valid]
    #   for INDEPENDENT random D12, D12'  = Pr[D12 valid]^2 = (N/total_d12)^2
    # So E[N^2] = total_d12^2 * (N/total_d12)^2 = N^2. Trivially exact.
    #
    # The interesting PZ direction is: fix D12 and randomize D11.
    # Or better: let N = #{(D11, D12) : valid pair}, randomize both.

    # Actually the MOST useful for the proof is:
    # Fix "good" D11 class. For random D12:
    # N = 1{valid}. Then Pr[N=1] = |valid| / total_d12.
    # But the theorist wants E[N^2] over ALL (D11, D12) pairs for Paley-Zygmund.

    # The useful version: N = number of valid D12 for random D11
    # E[N] = total_valid / total_d11 (for each D11, count valid D12)
    # E[N^2] = sum over D11 of N(D11)^2 / total_d11
    # Already computed above.

    # Aggregate overlap distribution
    print(f"\n--- AGGREGATE OVERLAP DISTRIBUTION ---")
    total_ordered = sum(all_overlap_dist.values())
    for s in sorted(all_overlap_dist.keys()):
        count = all_overlap_dist[s]
        pct = 100 * count / total_ordered
        print(f"  |D12 ∩ D12'| = {s}: {count} pairs ({pct:.1f}%)")

    # Random baseline overlap
    print(f"\n--- RANDOM OVERLAP BASELINE ---")
    rng = np.random.default_rng(42)
    random_dist = compute_random_overlap_distribution(p, k, num_samples=100000, rng=rng)
    for s, prob in sorted(random_dist.items()):
        print(f"  |D12 ∩ D12'| = {s}: {prob:.4f}")

    sys.stdout.flush()

    return {
        "p": p,
        "total_d11": total_d11,
        "working_d11": num_d11_analyzed,
        "total_N": total_N,
        "total_N_sq": total_N_sq,
        "E_N": E_N,
        "E_N_sq": E_N_sq,
        "ratio": ratio,
        "PZ_bound": E_N * E_N / E_N_sq if E_N_sq > 0 else 0,
        "actual_fraction": num_d11_analyzed / total_d11,
        "results_by_group": results_by_group,
        "overlap_distribution": dict(all_overlap_dist),
    }


def main():
    print("=" * 80)
    print("SECOND MOMENT ANALYSIS: E[N^2] FOR VALID D12 COUNT")
    print("=" * 80)
    sys.stdout.flush()

    all_results = {}
    for p in [11, 19, 23]:
        all_results[p] = analyze_prime(p)

    # Cross-prime summary
    print(f"\n{'='*80}")
    print("CROSS-PRIME PALEY-ZYGMUND SUMMARY")
    print(f"{'='*80}")
    print(f"{'p':>4s} {'E[N]':>10s} {'E[N^2]':>12s} {'E[N]^2':>12s} {'ratio':>10s} {'PZ_bound':>10s} {'actual':>10s}")
    for p, r in all_results.items():
        print(f"  {p:4d} {r['E_N']:10.4f} {r['E_N_sq']:12.2f} {r['E_N']*r['E_N']:12.4f} "
              f"{r['ratio']:10.4f} {r['PZ_bound']:10.6f} {r['actual_fraction']:10.6f}")

    # Does ratio grow as poly(p)?
    print(f"\nRatio E[N^2]/E[N]^2 growth:")
    for p, r in all_results.items():
        ratio = r["ratio"]
        print(f"  p={p}: ratio = {ratio:.4f}, log(ratio)/log(p) = {np.log(ratio)/np.log(p):.4f}")

    # Save results
    script_dir = os.path.dirname(os.path.abspath(__file__))
    outpath = os.path.join(script_dir, "..", "enumeration_data", "second_moment_analysis.json")
    save_data = {}
    for p, r in all_results.items():
        save_data[str(p)] = {k: v for k, v in r.items() if k != "results_by_group"}
    with open(outpath, "w") as f:
        json.dump(save_data, f, indent=2, default=str)
    print(f"\nResults saved to {outpath}")


if __name__ == "__main__":
    main()
