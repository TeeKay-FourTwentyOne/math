# p=43 Exhaustive Orbit Computation: Specification

## Goal

Determine exact p_working at p=43: the fraction of multiplicative orbits of symmetric D11 that admit at least one valid D12. This is the single most important missing data point for the paper.

## Parameters

- p = 43, n = (p+1)/2 = 22, k = (p-3)/2 = 20, m = p = 43
- Symmetric D11 ⊂ {1,...,42} with |D11| = 22 (each complementary pair {d, 43-d} either both in or both out)
- D12 = {0} ∪ S where S is a 20-subset of {1,...,42}
- Total symmetric D11: C(21, 11) = 352,716
- Orbits under Z_43* (multiplicative group): 352,716 / 21 = 16,796 orbits (each orbit has exactly 21 elements since 43 is prime)

## What to compute

For each of the 16,796 orbits:
1. Pick one representative D11 from the orbit
2. Run simulated annealing to search for a valid D12
3. Record: orbit ID, representative D11, whether a valid D12 was found, the D12 if found, and the SA cost (number of violated constraints at termination)

A valid (D11, D12) pair satisfies: A(d) + B(d) ≤ T(d) for all d = 1,...,42, where:
- A(d) = |{(a,b) ∈ D11 × D11 : a-b ≡ d mod 43, a≠b}|
- B(d) = |{(a,b) ∈ D12 × D12 : a-b ≡ d mod 43, a≠b}|
- T(d) = (p-3)/2 = 20 if d ∈ D11
- T(d) = (p+3)/2 = 23 if d ∈ D22 = {1,...,42} \ D11

## SA Design

### Cost function
cost(D12) = #{d : A(d) + B(d) > T(d)}

A valid D12 has cost = 0.

### Move operator
Swap one element of S = D12 \ {0} with one element of {1,...,42} \ S. This preserves |D12| = 21.

### Temperature schedule
Standard geometric cooling. Start temperature ~5, cooling rate 0.9995-0.9999, ~80K-100K iterations per trial.

### Trials per orbit
Run 5 independent SA trials per orbit representative. If ANY trial finds cost=0, the orbit is working. This gives false negative rate < 1% for orbits where the valid D12 fraction is ≥ 1e-7 (based on the p=43 correction data showing SA with 80K × 2 trials found 5/168 A-flat orbits).

### Early termination
If cost=0 found, stop remaining trials for this orbit immediately.

## Critical implementation details

### Generating orbit representatives
1. Enumerate all C(21,11) = 352,716 symmetric D11 by choosing 11 elements from {1,...,21} as "pair representatives" (each chosen element d includes both d and 43-d in D11)
2. For each D11, compute its canonical form under Z_43*: apply all 21 multipliers g ∈ {1,2,...,42}/±1, compute sorted(g·D11 mod 43) for each, take the lexicographic minimum as the canonical representative
3. Group by canonical form to get 16,796 orbits
4. Store one representative per orbit

### Precomputing A-values
For each orbit representative D11, precompute A(d) for all d = 1,...,42. This is fixed for the entire SA search and should not be recomputed on each iteration.

### Efficient B-value updates
When swapping element x_out for x_in in D12:
- B(d) changes only for d values related to x_out and x_in
- Specifically: for each existing element y ∈ D12, B(y - x_out) decreases by 1 and B(x_out - y) decreases by 1, and B(y - x_in) increases by 1 and B(x_in - y) increases by 1
- Also handle the 0 element: B(x_out) and B(43-x_out) decrease, B(x_in) and B(43-x_in) increase
- This gives O(|D12|) = O(p) update per move, not O(p²) full recomputation

### Verification
When SA reports cost=0, verify the solution independently by recomputing ALL A(d) + B(d) from scratch. Do not trust incremental updates alone.

## Parallelization

The 16,796 orbits are completely independent. The computation is embarrassingly parallel.

Suggested structure:
- Package as a single Python script that takes command-line arguments: `--start_orbit N --end_orbit M`
- Each invocation processes orbits N through M
- Output results to a JSON file per chunk: `results_N_M.json`
- A merge script combines all chunk results into the final table

This allows running on 1 core (all orbits sequentially), 8-16 cores (MacBook), or 100+ cores (cloud) with no code changes.

### Estimated runtime per orbit
~5 trials × 100K iterations × O(p) per iteration = ~2.5M × 43 ≈ 100M operations per orbit. At ~100M ops/sec in Python: ~1 second per orbit. Total: ~17K seconds ≈ 5 hours single-core.

With NumPy optimization of the inner loop, 10-50x speedup is achievable, bringing it to 5-30 minutes single-core. With 8 cores: under 5 minutes.

## Output format

```json
{
  "p": 43,
  "total_orbits": 16796,
  "orbits_tested": 16796,
  "working_orbits": <count>,
  "p_working": <fraction>,
  "p_times_p_working": <p * fraction>,
  "sa_config": {
    "trials_per_orbit": 5,
    "iterations_per_trial": 100000,
    "start_temp": 5.0,
    "cooling_rate": 0.9998
  },
  "results": [
    {
      "orbit_id": 0,
      "d11_representative": [list of elements],
      "working": true/false,
      "best_cost": 0,
      "valid_d12": [list of elements] or null,
      "trials_run": 1-5
    },
    ...
  ]
}
```

## Success criterion

The paper needs: p × p_working at p=43 with tight confidence bounds. If p × p_working > 0.5, the Θ(1/p) conjecture is well-supported. If p × p_working < 0.3, we have a problem that needs discussion in the paper.

## Existing code to reference

The project already has:
- `ramsey_core.py` or similar with `verify_construction()` function
- `near_flat_analysis.py` with SA implementation
- `invariant_hunt.py` with orbit enumeration logic
- `solutions_registry.json` with known valid pairs

Reuse the existing SA and verification code where possible. The main new work is the orbit enumeration and parallelization wrapper.
