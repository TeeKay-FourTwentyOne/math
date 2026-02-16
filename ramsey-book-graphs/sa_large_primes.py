#!/usr/bin/env python3
"""
SA search for valid 2-block circulant constructions at larger primes p ≡ 3 mod 4.
Jointly optimizes D11 and D12. Reports first solution found for each prime.
"""

import numpy as np
import time
import json
import sys


def is_prime(n):
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0: return False
    for d in range(3, int(n**0.5) + 1, 2):
        if n % d == 0: return False
    return True


def compute_cost(D11_set, D12_set, p):
    """Compute total constraint violation cost."""
    n = (p + 1) // 2
    D22_set = set(range(1, p)) - D11_set
    D12T_set = {(-x) % p for x in D12_set}

    cost = 0
    for d in range(1, p):
        # A(d) = |{x in D11 : (x+d)%p in D11}|
        A_d = sum(1 for x in D11_set if (x + d) % p in D11_set)
        # B(d) = |{x in D12 : (x+d)%p in D12}|
        B_d = sum(1 for x in D12_set if (x + d) % p in D12_set)
        # C(d) = |{x in D22 : (x+d)%p in D22}|
        C_d = sum(1 for x in D22_set if (x + d) % p in D22_set)
        # Bp(d) = B(p-d)
        Bp_d = sum(1 for x in D12_set if (x + p - d) % p in D12_set)

        if d in D11_set:
            # V1V1 red: A+B <= n-2
            cost += max(0, A_d + B_d - (n - 2))
            # V2V2 blue: C+Bp+2 <= n-1 => C+Bp <= n-3
            cost += max(0, C_d + Bp_d - (n - 3))
        else:
            # V2V2 red: C+Bp <= n-2
            cost += max(0, C_d + Bp_d - (n - 2))
            # V1V1 blue: A+B-2 <= n-1 => A+B <= n+1
            cost += max(0, A_d + B_d - (n + 1))

    # V1V2 constraints
    for d in range(p):
        sigma = sum(1 for x in D11_set if (x + d) % p in D12_set)
        sigma += sum(1 for x in D12_set if (x + d) % p in D11_set)
        delta_12_22 = sum(1 for x in D12_set if (x + d) % p in D22_set)
        X_d = sigma + delta_12_22

        if d in D12_set:
            cost += max(0, X_d - (n - 2))
        else:
            cost += max(0, X_d - (n - 1))

    return cost


def sa_search(p, max_iter=5000000, T_init=15.0, T_min=0.001, seed=None):
    """Simulated annealing search for valid (D11, D12) at prime p."""
    rng = np.random.RandomState(seed)
    n = (p + 1) // 2
    k = (p - 1) // 2

    # Initialize random symmetric D11
    half = (p - 1) // 2
    pairs = [(i + 1, p - (i + 1)) for i in range(half)]
    n_choose = (p + 1) // 4
    chosen = rng.choice(half, size=n_choose, replace=False)
    D11_set = set()
    for i in chosen:
        D11_set.add(pairs[i][0])
        D11_set.add(pairs[i][1])

    # Initialize random D12 containing 0
    others = list(range(1, p))
    rng.shuffle(others)
    D12_set = {0} | set(others[:k - 1])

    # Compute initial cost
    cost = compute_cost(D11_set, D12_set, p)
    best_cost = cost
    best_D11 = set(D11_set)
    best_D12 = set(D12_set)

    alpha = (T_min / T_init) ** (1.0 / max_iter)
    T = T_init

    d12_list = sorted(D12_set - {0})
    d12_non = sorted(set(range(1, p)) - D12_set)

    for it in range(max_iter):
        r = rng.random()

        if r < 0.25 and len(pairs) > 0:
            # D11 pair swap
            # Remove a pair, add another
            d11_pairs = [i for i in range(half) if pairs[i][0] in D11_set]
            d22_pairs = [i for i in range(half) if pairs[i][0] not in D11_set]
            if len(d11_pairs) > 0 and len(d22_pairs) > 0:
                remove_idx = d11_pairs[rng.randint(len(d11_pairs))]
                add_idx = d22_pairs[rng.randint(len(d22_pairs))]
                D11_set.discard(pairs[remove_idx][0])
                D11_set.discard(pairs[remove_idx][1])
                D11_set.add(pairs[add_idx][0])
                D11_set.add(pairs[add_idx][1])

                new_cost = compute_cost(D11_set, D12_set, p)
                delta = new_cost - cost
                if delta <= 0 or rng.random() < np.exp(-delta / T):
                    cost = new_cost
                else:
                    D11_set.discard(pairs[add_idx][0])
                    D11_set.discard(pairs[add_idx][1])
                    D11_set.add(pairs[remove_idx][0])
                    D11_set.add(pairs[remove_idx][1])
        else:
            # D12 swap
            if len(d12_list) > 0 and len(d12_non) > 0:
                rm_idx = rng.randint(len(d12_list))
                add_idx = rng.randint(len(d12_non))
                old_elem = d12_list[rm_idx]
                new_elem = d12_non[add_idx]

                D12_set.discard(old_elem)
                D12_set.add(new_elem)

                new_cost = compute_cost(D11_set, D12_set, p)
                delta = new_cost - cost
                if delta <= 0 or rng.random() < np.exp(-delta / T):
                    cost = new_cost
                    d12_list[rm_idx] = new_elem
                    d12_non[add_idx] = old_elem
                else:
                    D12_set.discard(new_elem)
                    D12_set.add(old_elem)

        if cost < best_cost:
            best_cost = cost
            best_D11 = set(D11_set)
            best_D12 = set(D12_set)
            if cost == 0:
                return True, best_D11, best_D12, it

        T *= alpha

    return False, best_D11, best_D12, max_iter


def main():
    primes = []
    for n in range(43, 500, 2):
        if n % 4 == 3 and is_prime(n):
            primes.append(n)

    max_trials = 20
    max_iter = 2000000
    results = []

    print(f"{'p':>5} {'found':>6} {'trial':>6} {'iter':>10} {'time':>8} {'best_cost':>10}")

    for p in primes:
        t0 = time.time()
        found = False
        best_overall = 999999

        for trial in range(max_trials):
            ok, D11, D12, iters = sa_search(p, max_iter=max_iter,
                                             T_init=20.0, T_min=0.001,
                                             seed=trial * 1000 + p)
            cost = 0 if ok else compute_cost(D11, D12, p)
            best_overall = min(best_overall, cost)

            if ok:
                elapsed = time.time() - t0
                print(f"{p:5d} {'YES':>6} {trial+1:6d} {iters:10d} {elapsed:8.1f}s {0:10d}")
                results.append({
                    'p': p,
                    'found': True,
                    'trial': trial + 1,
                    'time': elapsed,
                    'D11': sorted(D11),
                    'D12': sorted(D12),
                })
                found = True
                break

        if not found:
            elapsed = time.time() - t0
            print(f"{p:5d} {'NO':>6} {max_trials:6d} {'--':>10} {elapsed:8.1f}s {best_overall:10d}")
            results.append({
                'p': p,
                'found': False,
                'time': elapsed,
                'best_cost': best_overall,
            })

    # Summary
    print("\n" + "=" * 60)
    found_count = sum(1 for r in results if r['found'])
    total = len(results)
    print(f"Found solutions: {found_count}/{total}")
    not_found = [r for r in results if not r['found']]
    if not_found:
        print(f"Primes without solution: {[r['p'] for r in not_found]}")

    with open('ramsey-book-graphs/sa_large_primes_results.json', 'w') as f:
        json.dump(results, f, indent=2, default=int)


if __name__ == '__main__':
    main()
