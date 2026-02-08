# Plan: Final Documentation & Validator Script for R(B_21, B_22) = 87

## Context

We found a valid 2-block circulant graph on 86 vertices proving R(B_21, B_22) >= 87. The solution has been verified by the circulant-formula checker and by targeted checks (V2V2 blue, V1V1 red, all 6 category maxima). The Validator now requires a standalone O(n^3) brute-force script that builds the full 86x86 adjacency matrix and checks every pair directly.

### The Solution
```
D11 = [1,2,5,10,11,13,16,17,18,19,20,23,24,25,26,27,30,32,33,38,41,42]
D12 = [0,2,5,6,8,11,15,16,20,24,25,27,28,31,32,34,35,36,37,39,41]
D22 = complement(D11) = [3,4,6,7,8,9,12,14,15,21,22,28,29,31,34,35,36,37,39,40]
```

## Step 1: Write standalone validator script

Create `validate_n22_full.py` that does the following with NO dependency on `ramsey_core.py` (fully self-contained):

1. **Build the 86x86 adjacency matrix** from D11, D12, D22 directly:
   - Vertices 0-42 = V1, vertices 43-85 = V2
   - Edge (u,v) where u,v in V1: `(v-u) % 43 in D11`
   - Edge (u,v) where u,v in V2: `((v-43)-(u-43)) % 43 in D22`
   - Edge (u,v) where u in V1, v in V2: `(v-43-u) % 43 in D12` (for u < v with u in V1, v in V2)
   - Symmetric: adj[u][v] = adj[v][u]

2. **Verify structural properties**:
   - D11 is symmetric: d in D11 iff (43-d) in D11
   - D22 = {1..42} \ D11
   - 0 in D12
   - |D11|=22, |D12|=21, |D22|=20

3. **For each of C(86,2) = 3655 pairs (u,v)**:
   - Determine color: red if adj[u][v]==1, blue if adj[u][v]==0
   - Count common neighbors in that color by iterating all 84 other vertices w:
     - Red common neighbor: adj[u][w]==1 AND adj[v][w]==1
     - Blue common neighbor: adj[u][w]==0 AND adj[v][w]==0
   - For red edges: count red common neighbors, verify <= 20
   - For blue edges: count blue common neighbors, verify <= 21

4. **Report**:
   - Total red edges, total blue edges
   - Max red common neighbors (across all red edges)
   - Max blue common neighbors (across all blue edges)
   - Any violations found
   - Final PASS/FAIL

Key point: This directly counts blue common neighbors (vertices adjacent to neither endpoint) rather than using the inclusion-exclusion formula. This is the most transparent verification possible.

### File to create
- `ramsey-book-graphs/validate_n22_full.py` — standalone validator

### File to modify
- `ramsey-book-graphs/docs/EXPERIMENT_LOG.md` — add final validation result

## Step 2: Run the validator

Run `python -u validate_n22_full.py` and capture full output. Expected runtime: a few seconds (86^3 ~ 636K operations).

## Step 3: Update documentation

Append the validator output to `EXPERIMENT_LOG.md` under a new "Final Validation" section.

## Verification

The script output should show:
- 3655 total pairs checked
- Max red common neighbors = 20 (threshold 20)
- Max blue common neighbors = 21 (threshold 21)
- 0 violations
- PASS
