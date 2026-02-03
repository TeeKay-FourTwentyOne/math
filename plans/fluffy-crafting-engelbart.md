# Implementation Plan: R(B₂₁, B₂₂) Ramsey Search - Round 2

## Summary

Fix a critical bug in the blue neighbor formula and prepare the codebase for searching n=22.

## Critical Bug

**Location**: `ramsey-book-graphs/ramsey_core.py` lines 207, 224, 237

**Current (WRONG)**:
```python
blue_common = (N - 2) - common
```

**Correct (inclusion-exclusion)**:
```python
blue_common = (N - 2) - deg_u - deg_v + common
```

Where:
- `d1 = |D11| + |D12|` (degree of V1 vertices)
- `d2 = |D22| + |D12|` (degree of V2 vertices)

This bug makes blue constraints ~3x stricter than they should be. All prior searches used incorrect constraints.

---

## Implementation Steps

### Step 1: Fix `verify_construction()` in `ramsey_core.py`

1. Compute degrees at start of function:
   ```python
   d1 = len(G.D11) + len(G.D12)
   d2 = len(G.D22) + len(G.D12)
   ```

2. Update line 207 (V1-V1 blue edges):
   ```python
   blue_common = (N - 2) - d1 - d1 + common
   ```

3. Update line 224 (V2-V2 blue edges):
   ```python
   blue_common = (N - 2) - d2 - d2 + common
   ```

4. Update line 237 (V1-V2 blue edges):
   ```python
   blue_common = (N - 2) - d1 - d2 + common
   ```

### Step 2: Add Blue Formula Tests in `test_core.py`

Add brute-force blue neighbor counting:
```python
def blue_common_brute(G, u, v):
    """Count vertices that are NOT adjacent to either u or v."""
    count = 0
    for w in range(G.N):
        if w != u and w != v:
            if not G.adjacent(u, w) and not G.adjacent(v, w):
                count += 1
    return count
```

Add test class `TestBlueNeighborFormula`:
- Test that formula matches brute force for random graphs
- Test the mathematical identity: `red_common + blue_common + deg_u + deg_v - 2*adjacent(u,v) = N - 2`

### Step 3: Add D22 = Complement(D11) Mode (Optional Enhancement)

Add `use_complement_D22` parameter to `BlockCirculantGraph`:
- When enabled, D22 = Z_m \ {0} \ D11 (reduces search space by ~2^21)
- Update `SearchState` in `solver_heuristic.py` to support this mode
- Remove D22 moves from neighbor generation when in complement mode

### Step 4: Add Algebraic Seeding (Optional Enhancement)

Add to `solver_heuristic.py`:
```python
def quadratic_residues(m):
    return {(x * x) % m for x in range(1, m)}
```

Use QR(43) as initial seed for n=22 search even though 43 ≡ 3 (mod 4).

---

## Files to Modify

| File | Changes |
|------|---------|
| `ramsey_core.py` | Fix blue formula in `verify_construction()` |
| `test_core.py` | Add `TestBlueNeighborFormula` tests |
| `solver_heuristic.py` | Add complement mode + algebraic seeding (optional) |
| `solver_sat.py` | Fix blue encoding (lower priority, same bug) |

---

## Verification

1. Run `python -m pytest test_core.py -v` - all tests should pass
2. Run blue formula sanity check with small n (n=4,5,6)
3. Verify known solutions for n ≤ 21 still work (if we have solution data)
4. Run short heuristic search for n=6 to validate end-to-end

---

## Priority

1. **Critical**: Fix blue formula bug (Step 1)
2. **High**: Add tests to catch this bug (Step 2)
3. **Medium**: Complement mode (Step 3)
4. **Low**: Algebraic seeding, SAT solver fixes (Step 4)

After the fix, the system will be ready to run searches for n=22.
