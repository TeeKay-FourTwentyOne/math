#!/usr/bin/env python3
"""Extended second moment computation for Ramsey book graph construction.

Computes E[N^2]/E[N]^2 for primes p=43 and p=47,
where N(D11) = number of valid D12 for a given symmetric D11.

KEY INSIGHT: For p >= 43, uniform random D12 sampling finds 0 valid D12
even with 10M samples (rate ~1e-7 for the known solution). Instead, we
use a MARGINAL PROBABILITY PRODUCT proxy:

  log N(D11) ~ sum_d log Pr[B(d) <= slack(d)]

where slack(d) = threshold - A(d), computed over random D12. This product
treats the constraints as approximately independent (which is justified by
the low correlation rho = -2/(p-3)).

For the second moment, the ratio E[N^2]/E[N]^2 measures concentration.
We compute the proxy log-N for many D11 and use the distribution to
estimate the ratio.

Usage: python -u second_moment_extended.py
"""

import numpy as np
import json
import time
import sys
import os
from math import comb, log2, log, sqrt, exp, lgamma
from itertools import combinations
from collections import defaultdict


def symmetric_pairs(p):
    return [(d, p - d) for d in range(1, (p + 1) // 2)]


def enumerate_all_symmetric_d11(p):
    n = (p + 1) // 2
    num_choose = n // 2
    pairs = symmetric_pairs(p)
    all_d11 = []
    for chosen in combinations(range(len(pairs)), num_choose):
        d11 = set()
        for idx in chosen:
            d, comp = pairs[idx]
            d11.add(d)
            d11.add(comp)
        all_d11.append(frozenset(d11))
    return all_d11, pairs


def batch_autocorrelation_fft(all_d11, p, batch_size=5000):
    results = {}
    d11_list = list(all_d11)
    for start in range(0, len(d11_list), batch_size):
        end = min(start + batch_size, len(d11_list))
        batch = d11_list[start:end]
        bs = len(batch)
        indicator_matrix = np.zeros((bs, p), dtype=np.float64)
        for i, d11 in enumerate(batch):
            for j in d11:
                indicator_matrix[i, j] = 1.0
        fft_vals = np.fft.fft(indicator_matrix, axis=1)
        autocorr = np.fft.ifft(np.abs(fft_vals) ** 2, axis=1).real
        A_batch = np.round(autocorr).astype(np.int32)
        for i, d11 in enumerate(batch):
            results[d11] = A_batch[i]
    return results


def estimate_marginal_rates(p, A, d11_arr, d22_arr, num_samples, batch_size, rng):
    """Estimate marginal Pr[B(d) <= threshold - A(d)] for each constraint d.

    Returns: (binding_rates, loose_rates, log2_prod_binding, log2_prod_loose,
              full_valid_count, total_sampled)
    """
    k = (p - 3) // 2
    threshold_binding = (p - 3) // 2
    threshold_loose = (p + 3) // 2
    candidates = np.arange(1, p)

    # Accumulators for marginal rates
    binding_ok_counts = np.zeros(len(d11_arr), dtype=np.int64)
    loose_ok_counts = np.zeros(len(d22_arr), dtype=np.int64)
    full_valid_count = 0
    total_sampled = 0

    while total_sampled < num_samples:
        bs = min(batch_size, num_samples - total_sampled)
        indicator_matrix = np.zeros((bs, p), dtype=np.float64)
        indicator_matrix[:, 0] = 1.0
        for i in range(bs):
            chosen = rng.choice(candidates, size=k, replace=False)
            indicator_matrix[i, chosen] = 1.0
        fft_vals = np.fft.fft(indicator_matrix, axis=1)
        B_batch = np.round(np.fft.ifft(np.abs(fft_vals) ** 2, axis=1).real).astype(np.int32)

        # Binding constraints: A(d) + B(d) <= threshold_binding for d in D11
        F_binding = A[d11_arr][np.newaxis, :] + B_batch[:, d11_arr]
        binding_ok = (F_binding <= threshold_binding)  # (bs, len(d11_arr))
        binding_ok_counts += binding_ok.sum(axis=0)

        # Loose constraints: A(d) + B(d) <= threshold_loose for d in D22
        if d22_arr.size > 0:
            F_loose = A[d22_arr][np.newaxis, :] + B_batch[:, d22_arr]
            loose_ok = (F_loose <= threshold_loose)
            loose_ok_counts += loose_ok.sum(axis=0)

            # Full validity check
            all_ok = np.all(binding_ok, axis=1) & np.all(loose_ok, axis=1)
            full_valid_count += int(all_ok.sum())
        else:
            all_ok = np.all(binding_ok, axis=1)
            full_valid_count += int(all_ok.sum())

        total_sampled += bs

    # Compute marginal rates
    binding_rates = binding_ok_counts / total_sampled
    loose_rates = loose_ok_counts / total_sampled if d22_arr.size > 0 else np.array([])

    # Log2 product of marginals (proxy for log N)
    eps = 1e-30
    log2_prod_binding = float(np.sum(np.log2(np.maximum(binding_rates, eps))))
    log2_prod_loose = float(np.sum(np.log2(np.maximum(loose_rates, eps)))) if d22_arr.size > 0 else 0.0

    return (binding_rates, loose_rates, log2_prod_binding, log2_prod_loose,
            full_valid_count, total_sampled)


def load_known_solutions(p):
    """Load known D11 solutions for a given prime."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    solutions = {}
    known = {
        31: {6, 7, 8, 10, 11, 12, 14, 15, 16, 17, 19, 20, 21, 23, 24, 25},
        43: {1, 2, 5, 10, 11, 13, 16, 17, 18, 19, 20, 23, 24, 25, 26, 27,
             30, 32, 33, 38, 41, 42},
    }
    if p in known:
        solutions["hardcoded"] = known[p]

    # Check registry
    registry_path = os.path.join(script_dir, "solutions_registry.json")
    if os.path.exists(registry_path):
        try:
            with open(registry_path) as f:
                registry = json.load(f)
            for sol in registry.get("solutions", []):
                if sol["m"] == p:
                    D11_raw = set(sol["D11"])
                    target = (p + 1) // 2
                    if len(D11_raw) != target:
                        D11_raw = set(range(1, p)) - D11_raw
                    if len(D11_raw) == target:
                        solutions["registry"] = D11_raw
                    break
        except:
            pass

    # Check SA solution files
    for fname in [f"solution_n{(p+1)//2}_sa.json", f"p{p}_solutions.json"]:
        fpath = os.path.join(script_dir, fname)
        if os.path.exists(fpath):
            try:
                with open(fpath) as f:
                    data = json.load(f)
                if isinstance(data, dict) and "D11" in data:
                    D11_raw = set(data["D11"])
                    target = (p + 1) // 2
                    if len(D11_raw) != target:
                        D11_raw = set(range(1, p)) - D11_raw
                    solutions[fname] = D11_raw
                elif isinstance(data, list):
                    for sol in data[:5]:
                        if "D11" in sol:
                            D11_raw = set(sol["D11"])
                            target = (p + 1) // 2
                            if len(D11_raw) != target:
                                D11_raw = set(range(1, p)) - D11_raw
                            solutions[f"{fname}_0"] = D11_raw
                            break
            except:
                pass

    return solutions


def analyze_prime(p, d12_samples_per_d11, batch_size, num_d11_reps):
    """Analyze a single prime using marginal probability product proxy.

    For each sampled D11:
    1. Compute autocorrelation A(d)
    2. Sample D12 and measure marginal Pr[B(d) <= slack(d)] for each d
    3. Compute log2_N_proxy = sum_d log2(marginal_rate(d)) + log2(|D12_space|)
    4. Compute N_proxy = 2^(log2_N_proxy)

    Then compute E[N_proxy], E[N_proxy^2], ratio.
    """
    t_start = time.time()
    n = (p + 1) // 2
    k = (p - 3) // 2
    total_d12 = comb(p - 1, k)
    log2_total_d12 = log2(total_d12)
    a_flat_threshold = (p + 1) // 4
    threshold_binding = (p - 3) // 2

    num_pairs = (p - 1) // 2
    num_choose = n // 2
    total_d11 = comb(num_pairs, num_choose)

    print(f"\n{'='*80}")
    print(f"SECOND MOMENT ANALYSIS FOR p = {p}")
    print(f"{'='*80}")
    print(f"  n={n}, |D11|={n}, |D12|={k+1}, A-flat threshold={a_flat_threshold}")
    print(f"  C({num_pairs},{num_choose}) = {total_d11:,} symmetric D11")
    print(f"  C({p-1},{k}) = {total_d12:,} total D12")
    print(f"  log2(|D12|) = {log2_total_d12:.2f}")
    print(f"  Binding threshold: {threshold_binding}")
    sys.stdout.flush()

    # Step 1: Enumerate and compute A
    print(f"\n  Enumerating D11...")
    t0 = time.time()
    all_d11, pairs = enumerate_all_symmetric_d11(p)
    print(f"  {len(all_d11):,} D11 in {time.time()-t0:.1f}s")
    sys.stdout.flush()

    print(f"  Computing autocorrelation...")
    t0 = time.time()
    A_dict = batch_autocorrelation_fft(all_d11, p, batch_size=5000)
    print(f"  Done in {time.time()-t0:.1f}s")
    sys.stdout.flush()

    # Step 2: Classify by max_A
    maxA_to_d11 = defaultdict(list)
    for d11 in all_d11:
        A = A_dict[d11]
        max_A = max(int(A[d]) for d in d11)
        maxA_to_d11[max_A].append(d11)

    max_A_dist = {m: len(lst) for m, lst in maxA_to_d11.items()}
    num_aflat = sum(cnt for m, cnt in max_A_dist.items() if m <= a_flat_threshold)

    print(f"\n  max_A distribution:")
    for val in sorted(max_A_dist.keys()):
        tag = ""
        if val <= a_flat_threshold:
            tag = " [A-flat]"
        elif val > threshold_binding:
            tag = " [INFEASIBLE]"
        print(f"    max_A={val:3d}: {max_A_dist[val]:6d} D11{tag}")
    sys.stdout.flush()

    # Step 3: Known solutions
    print(f"\n  Loading known solutions...")
    known_sols = load_known_solutions(p)
    known_d11_fs = None
    known_maxA = None

    for name, D11_set in known_sols.items():
        D11_fs = frozenset(D11_set)
        if D11_fs in A_dict:
            A = A_dict[D11_fs]
            mA = max(int(A[d]) for d in D11_set)
            print(f"    {name}: max_A={mA}, A-flat={mA <= a_flat_threshold}")
            known_d11_fs = D11_fs
            known_maxA = mA
        else:
            print(f"    {name}: NOT in enumerated D11")
    sys.stdout.flush()

    # Step 4: Compute marginal proxy for stratified sample of D11
    rng = np.random.default_rng(seed=2024)

    print(f"\n  Computing marginal proxy for D11 sample "
          f"({num_d11_reps} reps/class, {d12_samples_per_d11:,} D12/D11)...")
    sys.stdout.flush()

    # Collect proxy values
    all_log2_N_proxy = []  # (weight, log2_N_proxy, max_A, is_aflat)
    class_results = {}

    for max_A in sorted(maxA_to_d11.keys()):
        d11_list = maxA_to_d11[max_A]
        class_size = len(d11_list)

        if max_A > threshold_binding:
            # Infeasible: N=0
            class_results[max_A] = {
                "count": class_size, "n_sampled": 0,
                "log2_proxy_mean": float('-inf'),
                "log2_proxy_std": 0.0,
            }
            continue

        # Sample representatives
        n_reps = min(num_d11_reps, class_size)

        # For known class, sample more
        if known_maxA is not None and max_A == known_maxA:
            n_reps = min(num_d11_reps * 3, class_size)

        if n_reps < class_size:
            indices = rng.choice(class_size, size=n_reps, replace=False)
            reps = [d11_list[i] for i in indices]
        else:
            reps = list(d11_list)

        # Include known D11 if in this class
        if known_d11_fs is not None and max_A == known_maxA:
            if known_d11_fs not in reps:
                reps.append(known_d11_fs)

        log2_proxies = []
        t0 = time.time()

        for d11 in reps:
            A = A_dict[d11]
            d11_arr = np.array(sorted(d11), dtype=np.int32)
            d22_arr = np.array(sorted(set(range(1, p)) - d11), dtype=np.int32)

            (br, lr, l2pb, l2pl, fvc, ts) = estimate_marginal_rates(
                p, A, d11_arr, d22_arr, d12_samples_per_d11, batch_size, rng
            )

            log2_N = l2pb + l2pl + log2_total_d12
            log2_proxies.append(log2_N)

            # Store per-D11 data
            is_known = (d11 == known_d11_fs)
            is_aflat = (max_A <= a_flat_threshold)
            weight = class_size / n_reps  # importance weight
            all_log2_N_proxy.append({
                "log2_N": log2_N,
                "log2_binding": l2pb,
                "log2_loose": l2pl,
                "max_A": max_A,
                "is_aflat": is_aflat,
                "is_known": is_known,
                "weight": weight,
                "full_valid_count": fvc,
                "total_sampled": ts,
                "binding_min": float(br.min()) if len(br) > 0 else 0,
                "binding_max": float(br.max()) if len(br) > 0 else 0,
                "binding_mean": float(br.mean()) if len(br) > 0 else 0,
                "loose_min": float(lr.min()) if len(lr) > 0 else 0,
                "loose_max": float(lr.max()) if len(lr) > 0 else 0,
                "loose_mean": float(lr.mean()) if len(lr) > 0 else 0,
            })

        elapsed = time.time() - t0

        arr = np.array(log2_proxies)
        tag = " <--KNOWN" if (known_maxA is not None and max_A == known_maxA) else ""
        aflat_tag = " [A-flat]" if max_A <= a_flat_threshold else ""

        cr = {
            "count": class_size, "n_sampled": len(reps),
            "log2_proxy_mean": float(np.mean(arr)),
            "log2_proxy_std": float(np.std(arr)),
            "log2_proxy_min": float(np.min(arr)),
            "log2_proxy_max": float(np.max(arr)),
        }
        class_results[max_A] = cr

        print(f"    max_A={max_A:2d}: {len(reps):4d} reps, "
              f"log2(N)={cr['log2_proxy_mean']:8.1f} +/- {cr['log2_proxy_std']:5.1f} "
              f"[{cr['log2_proxy_min']:8.1f}, {cr['log2_proxy_max']:8.1f}] "
              f"({elapsed:.1f}s){tag}{aflat_tag}")
        sys.stdout.flush()

    # Step 5: Compute ratio using proxy
    print(f"\n  Computing second moment ratio from marginal proxy...")

    # The proxy gives log2(N) for each D11. To compute E[N^2]/E[N]^2,
    # we use: N_proxy = 2^log2_N, and compute the ratio.
    #
    # For numerical stability, we work in log space:
    # E[N] = sum(w_i * 2^x_i) / sum(w_i)
    # E[N^2] = sum(w_i * 2^(2*x_i)) / sum(w_i)
    #
    # To avoid overflow, use log-sum-exp trick.

    # Separate by all and A-flat
    def compute_ratio(entries, label):
        if not entries:
            print(f"    {label}: no entries")
            return float('inf'), {}

        log2_Ns = np.array([e["log2_N"] for e in entries])
        weights = np.array([e["weight"] for e in entries])

        # Filter out -inf entries
        finite = np.isfinite(log2_Ns)
        if not np.any(finite):
            print(f"    {label}: all log2(N) = -inf")
            return float('inf'), {}

        log2_Ns_f = log2_Ns[finite]
        weights_f = weights[finite]

        # Convert to natural log for logsumexp
        ln2 = np.log(2)
        ln_Ns = log2_Ns_f * ln2
        ln_weights = np.log(weights_f)

        # E[N] using log-sum-exp
        # E[N] = sum(w * exp(ln_N)) / sum(w)
        ln_wN = ln_weights + ln_Ns
        max_ln_wN = np.max(ln_wN)
        ln_sum_wN = max_ln_wN + np.log(np.sum(np.exp(ln_wN - max_ln_wN)))

        ln_sum_w = np.log(np.sum(weights_f))
        ln_E_N = ln_sum_wN - ln_sum_w

        # E[N^2] using log-sum-exp
        ln_wN2 = ln_weights + 2 * ln_Ns
        max_ln_wN2 = np.max(ln_wN2)
        ln_sum_wN2 = max_ln_wN2 + np.log(np.sum(np.exp(ln_wN2 - max_ln_wN2)))
        ln_E_N2 = ln_sum_wN2 - ln_sum_w

        # Ratio = E[N^2] / E[N]^2 = exp(ln_E_N2 - 2*ln_E_N)
        ln_ratio = ln_E_N2 - 2 * ln_E_N
        ratio = np.exp(ln_ratio)

        # Also compute log2 versions
        log2_E_N = ln_E_N / ln2
        log2_E_N2 = ln_E_N2 / ln2
        log2_ratio = ln_ratio / ln2

        stats = {
            "log2_E_N": float(log2_E_N),
            "log2_E_N2": float(log2_E_N2),
            "ratio": float(ratio),
            "log2_ratio": float(log2_ratio),
            "n_finite": int(np.sum(finite)),
            "n_total": len(entries),
            "log2_N_mean": float(np.mean(log2_Ns_f)),
            "log2_N_std": float(np.std(log2_Ns_f)),
            "log2_N_min": float(np.min(log2_Ns_f)),
            "log2_N_max": float(np.max(log2_Ns_f)),
        }

        print(f"    {label}:")
        print(f"      n_finite = {stats['n_finite']} / {stats['n_total']}")
        print(f"      log2(N) mean={stats['log2_N_mean']:.1f}, "
              f"std={stats['log2_N_std']:.1f}, "
              f"range=[{stats['log2_N_min']:.1f}, {stats['log2_N_max']:.1f}]")
        print(f"      log2(E[N])   = {stats['log2_E_N']:.2f}")
        print(f"      log2(E[N^2]) = {stats['log2_E_N2']:.2f}")
        print(f"      ratio = E[N^2]/E[N]^2 = {stats['ratio']:.4f}")
        print(f"      log2(ratio) = {stats['log2_ratio']:.2f}")

        return float(ratio), stats

    # All D11
    ratio_all, stats_all = compute_ratio(all_log2_N_proxy, "ALL D11")

    # A-flat D11
    aflat_entries = [e for e in all_log2_N_proxy if e["is_aflat"]]
    ratio_aflat, stats_aflat = compute_ratio(aflat_entries, "A-FLAT D11")

    # Known D11
    known_entries = [e for e in all_log2_N_proxy if e.get("is_known")]
    if known_entries:
        ke = known_entries[0]
        print(f"\n    KNOWN D11:")
        print(f"      log2(N_proxy) = {ke['log2_N']:.2f}")
        print(f"      log2_binding  = {ke['log2_binding']:.2f}")
        print(f"      log2_loose    = {ke['log2_loose']:.2f}")
        print(f"      binding_rate  = [{ke['binding_min']:.4f}, {ke['binding_max']:.4f}]"
              f" mean={ke['binding_mean']:.4f}")
        print(f"      loose_rate    = [{ke['loose_min']:.4f}, {ke['loose_max']:.4f}]"
              f" mean={ke['loose_mean']:.4f}")
        print(f"      full_valid    = {ke['full_valid_count']}/{ke['total_sampled']}")

    elapsed_total = time.time() - t_start
    print(f"\n    Total time: {elapsed_total:.1f}s")
    sys.stdout.flush()

    return {
        "p": p, "n": n,
        "total_d11": total_d11,
        "total_d12": total_d12,
        "num_aflat": num_aflat,
        "ratio_all": ratio_all,
        "ratio_aflat": ratio_aflat,
        "stats_all": stats_all,
        "stats_aflat": stats_aflat,
        "class_results": {str(k): v for k, v in class_results.items()},
        "max_A_distribution": {str(k): v for k, v in sorted(max_A_dist.items())},
        "known_maxA": known_maxA,
        "d11_proxy_data": all_log2_N_proxy,
        "elapsed": elapsed_total,
    }


def growth_rate_analysis(all_ratios):
    """Analyze growth rate of the second moment ratio vs p."""
    print(f"\n{'='*80}")
    print("GROWTH RATE ANALYSIS")
    print(f"{'='*80}")

    print(f"\n  {'p':>4s} {'ratio_all':>12s} {'ratio_aflat':>12s} {'source':>16s}")
    print(f"  {'-'*48}")
    for e in all_ratios:
        r_all = e.get("ratio_all", None)
        r_af = e.get("ratio_aflat", None)
        r_all_s = f"{r_all:.2f}" if isinstance(r_all, (int, float)) and r_all < 1e10 else "inf"
        r_af_s = f"{r_af:.2f}" if isinstance(r_af, (int, float)) and r_af < 1e10 else "inf"
        print(f"  {e['p']:4d} {r_all_s:>12s} {r_af_s:>12s} {e.get('source',''):>16s}")

    # Fit to ratio_all
    valid = [e for e in all_ratios
             if isinstance(e.get("ratio_all"), (int, float)) and 0 < e["ratio_all"] < 1e10]
    ps = np.array([e["p"] for e in valid])
    rs = np.array([e["ratio_all"] for e in valid])

    if len(ps) < 3:
        print("\n  Insufficient data for fitting.")
        return

    print(f"\n  Model fits for ratio_all ({len(ps)} points):")

    fits = {}
    # O(1)
    m = float(np.mean(rs))
    ssr = float(np.sum((rs - m) ** 2))
    print(f"    O(1):      mean={m:.1f}, SSR={ssr:.1f}")
    fits["O(1)"] = ssr

    # O(log p)
    c = np.polyfit(np.log(ps), rs, 1)
    ssr = float(np.sum((rs - np.polyval(c, np.log(ps))) ** 2))
    print(f"    O(log p):  {c[0]:.1f}*ln(p)+{c[1]:.1f}, SSR={ssr:.1f}")
    fits["O(log p)"] = ssr

    # O(sqrt p)
    c = np.polyfit(np.sqrt(ps), rs, 1)
    ssr = float(np.sum((rs - np.polyval(c, np.sqrt(ps))) ** 2))
    print(f"    O(sqrt p): {c[0]:.1f}*sqrt(p)+{c[1]:.1f}, SSR={ssr:.1f}")
    fits["O(sqrt p)"] = ssr

    # O(p)
    c = np.polyfit(ps, rs, 1)
    ssr = float(np.sum((rs - np.polyval(c, ps)) ** 2))
    print(f"    O(p):      {c[0]:.2f}*p+{c[1]:.1f}, SSR={ssr:.1f}")
    fits["O(p)"] = ssr

    # O(p^2)
    c = np.polyfit(ps ** 2, rs, 1)
    ssr = float(np.sum((rs - np.polyval(c, ps ** 2)) ** 2))
    print(f"    O(p^2):    {c[0]:.4f}*p^2+{c[1]:.1f}, SSR={ssr:.1f}")
    fits["O(p^2)"] = ssr

    # Power law
    try:
        from scipy.optimize import curve_fit
        popt, _ = curve_fit(lambda x, a, c: a * x ** c, ps, rs, p0=[0.01, 2], maxfev=10000)
        ssr = float(np.sum((rs - popt[0] * ps ** popt[1]) ** 2))
        print(f"    a*p^c:     {popt[0]:.6f}*p^{popt[1]:.2f}, SSR={ssr:.1f}")
        fits["a*p^c"] = ssr
    except:
        pass

    best = min(fits, key=fits.get)
    print(f"\n  Best fit: {best}")

    # Log-log analysis
    print(f"\n  Log-log:")
    for i in range(len(ps)):
        print(f"    p={ps[i]:.0f}: log2(ratio)={log2(rs[i]):.2f}, log2(p)={log2(ps[i]):.2f}")

    if len(ps) >= 2:
        print(f"  Successive log-log slopes:")
        for i in range(1, len(ps)):
            sl = (log2(rs[i]) - log2(rs[i - 1])) / (log2(ps[i]) - log2(ps[i - 1]))
            print(f"    p={ps[i - 1]:.0f}->{ps[i]:.0f}: {sl:.2f}")

    # A-flat
    valid_af = [e for e in all_ratios
                if isinstance(e.get("ratio_aflat"), (int, float)) and 0 < e["ratio_aflat"] < 1e10]
    if len(valid_af) >= 3:
        ps2 = np.array([e["p"] for e in valid_af])
        rs2 = np.array([e["ratio_aflat"] for e in valid_af])
        print(f"\n  ratio_aflat log-log:")
        for i in range(len(ps2)):
            print(f"    p={ps2[i]:.0f}: log2={log2(rs2[i]):.2f}")
        if len(ps2) >= 2:
            print(f"  Successive slopes:")
            for i in range(1, len(ps2)):
                sl = (log2(rs2[i]) - log2(rs2[i - 1])) / (log2(ps2[i]) - log2(ps2[i - 1]))
                print(f"    p={ps2[i - 1]:.0f}->{ps2[i]:.0f}: {sl:.2f}")

    print(f"\n  CONCLUSION:")
    print(f"  ratio_all best model: {best}")
    print(f"  If ratio = O(poly(p)), Paley-Zygmund gives Pr[N>0] = Omega(1/poly(p)),")
    print(f"  sufficient for the existence proof.")
    sys.stdout.flush()


def main():
    print("=" * 80)
    print("EXTENDED SECOND MOMENT COMPUTATION (MARGINAL PROXY)")
    print("R(B_{n-1}, B_n) = 4n-1 for p = 3 mod 4 primes")
    print("=" * 80)
    print(f"Date: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print()
    print("METHOD: Instead of counting exact valid D12 (infeasible for p >= 43),")
    print("we estimate log N(D11) via product of marginal constraint probabilities.")
    print("This proxy treats B(d) values at different d as approximately independent")
    print("(justified by correlation rho = -2/(p-3) -> 0).")
    sys.stdout.flush()

    # Known exact results from exhaustive enumeration
    known = [
        {"p": 11, "ratio_all": 2.0, "ratio_aflat": 1.0, "source": "exact_enum"},
        {"p": 19, "ratio_all": 14.0, "ratio_aflat": 3.0, "source": "exact_enum"},
        {"p": 23, "ratio_all": 14.5, "ratio_aflat": 2.95, "source": "exact_enum"},
        {"p": 31, "ratio_all": 134.0, "ratio_aflat": 6.9, "source": "sampled"},
    ]

    results = {}

    # p = 43: 352K D11
    # Use 100K D12 per D11 for marginal estimation, 30 reps per class
    r43 = analyze_prime(
        p=43,
        d12_samples_per_d11=100_000,
        batch_size=10000,
        num_d11_reps=30,
    )
    results["43"] = r43

    # p = 47: 1.35M D11
    # Use 100K D12 per D11, 20 reps per class
    r47 = analyze_prime(
        p=47,
        d12_samples_per_d11=100_000,
        batch_size=10000,
        num_d11_reps=20,
    )
    results["47"] = r47

    # Combined summary
    print(f"\n{'='*80}")
    print("COMBINED SUMMARY")
    print(f"{'='*80}")

    all_ratios = list(known)
    for p_str, r in sorted(results.items(), key=lambda x: int(x[0])):
        all_ratios.append({
            "p": r["p"],
            "ratio_all": r["ratio_all"],
            "ratio_aflat": r["ratio_aflat"],
            "source": "marginal_proxy",
        })
    all_ratios.sort(key=lambda x: x["p"])

    print(f"\n  {'p':>4s} {'total_D11':>12s} {'ratio_all':>12s} {'ratio_aflat':>12s}")
    print(f"  {'-'*45}")
    for e in known:
        print(f"  {e['p']:4d} {'(exact)':>12s} {e['ratio_all']:12.2f} {e['ratio_aflat']:12.2f}")
    for p_str, r in sorted(results.items(), key=lambda x: int(x[0])):
        r_all = f"{r['ratio_all']:.2f}" if r['ratio_all'] < 1e10 else "inf"
        r_af = f"{r['ratio_aflat']:.2f}" if r['ratio_aflat'] < 1e10 else "inf"
        print(f"  {r['p']:4d} {r['total_d11']:12,d} {r_all:>12s} {r_af:>12s}")

    growth_rate_analysis(all_ratios)

    # Save (strip large proxy data for JSON)
    output_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "second_moment_extended_results.json"
    )
    save_results = {}
    for p_str, r in results.items():
        r_save = dict(r)
        # Keep only summary of proxy data
        proxy_summary = []
        for e in r.get("d11_proxy_data", []):
            proxy_summary.append({
                "log2_N": e["log2_N"],
                "log2_binding": e["log2_binding"],
                "log2_loose": e["log2_loose"],
                "max_A": e["max_A"],
                "is_aflat": e["is_aflat"],
                "is_known": e["is_known"],
                "weight": e["weight"],
                "binding_mean": e["binding_mean"],
                "loose_mean": e["loose_mean"],
            })
        r_save["d11_proxy_data"] = proxy_summary
        save_results[p_str] = r_save

    save_data = {
        "known": known,
        "computed": save_results,
        "all_ratios": all_ratios,
        "timestamp": time.strftime('%Y-%m-%d %H:%M:%S'),
        "method": "marginal_probability_product_proxy",
    }
    with open(output_path, "w") as f:
        json.dump(save_data, f, indent=2, default=str)
    print(f"\n  Results saved to {output_path}")
    sys.stdout.flush()


if __name__ == "__main__":
    main()
