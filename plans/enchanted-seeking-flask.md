# Agent Team Design: Solving R(B_{n-1}, B_n) = 4n-1

## Context

The goal is to prove R(B_{n-1}, B_n) = 4n-1 for all n, where B_n = K_2 + K̄_n (book graph). The upper bound is proven (Rousseau & Sheehan 1978). The lower bound is verified for n ≤ 22 and when 2n-1 is a prime power ≡ 1 (mod 4). The n=22 case was recently cracked in this repo using informed simulated annealing after discovering universal structural patterns in known constructions.

**Next frontier cases:**
- n=24 (m=47, prime ≡ 3 mod 4) — same type as n=22; informed SA likely works
- n=25 (m=49=7², prime power ≡ 1 mod 4) — Paley over GF(49) should work
- n=23 (m=45=9×5, composite) — hardest near-term target; may need new construction family
- General n — the ultimate goal; requires a universal construction or new proof technique

**Key insight from n=22:** The breakthrough came from analyzing known solutions to discover universal patterns (D22=complement(D11), 0∈D12, degrees {m, m-2}), then constraining the SA search accordingly. The solution was found in ~370K iterations.

---

## Team: `ramsey-attack`

### Team Lead (coordinator)
Me. Assigns tasks, routes results between agents, makes strategic decisions, synthesizes findings.

### Agent 1: `researcher`
**Role:** Pattern analysis, algebraic theory, and proof strategy.
**subagent_type:** `general-purpose`

Combines the analytical and theoretical workstreams. Responsible for:
- Deep structural analysis of all known constructions (n=6..22), including spectral properties, cyclotomic class membership, and D12 structure
- Formulating algebraic conjectures about what determines D11 and D12
- Studying composite m cases (CRT decomposition over Z_9 × Z_5 for n=23, etc.)
- Pursuing the "Delta fluctuation" proof — can we show that for symmetric subsets of Z_p with |D11|=(p+1)/2, the Delta function always fluctuates enough?
- Literature search for recent progress (Wesley 2024, Lidicky et al. 2024 updates)
- Writing a formal proof outline partitioning all n by proof method
- Reporting structural constraints to `solver` to inform searches

**Key files:** `analyze_known.py`, `ramsey_core.py`, `solution_n22.json`, `docs/EXPERIMENT_LOG.md`, `ramsey-book-graphs.pdf`

### Agent 2: `solver`
**Role:** Computational construction search for specific n values.
**subagent_type:** `general-purpose`

Responsible for:
- Generalizing `sa_n22_informed.py` into a parameterized solver for arbitrary n
- Solving n=24 (m=47) using informed SA with universal pattern constraints
- Implementing GF(49) Paley construction for n=25
- Attacking n=23 (m=45) via informed SA, CRT-aware search, and dihedral Cayley alternatives
- Running long SA searches as background tasks
- Incorporating new constraints from `researcher` findings

**Key files:** `sa_n22_informed.py`, `solver_heuristic.py`, `ramsey_core.py`, `cayley_dihedral.py`

### Agent 3: `validator`
**Role:** Independent verification and correctness guarantees.
**subagent_type:** `general-purpose`

Responsible for:
- Generalizing `validate_n22_full.py` into a standalone validator for arbitrary n
- Independently verifying every claimed construction (brute-force O(N³) check, zero dependency on ramsey_core.py)
- Cross-validating formula implementations (Delta/Sigma fast formulas vs adjacency matrix)
- Computationally verifying theoretical claims from `researcher`
- Maintaining a solutions registry (JSON) tracking all verified results
- Running `test_core.py` and extending it for new n values

**Key files:** `validate_n22_full.py`, `test_core.py`, `ramsey_core.py`

---

## Task Phases

### Phase 1: Foundation (all parallel)

| Task | Agent | Description |
|------|-------|-------------|
| Generalize solver | solver | Parameterize `sa_n22_informed.py` for arbitrary n with auto-detected pattern constraints |
| Generalize validator | validator | Parameterize `validate_n22_full.py` for arbitrary n, standalone |
| Deep pattern analysis | researcher | Analyze all known constructions (n=6..22): spectral properties, D12 structure, cyclotomic classes. Categorize all n≤50 by proof method |

### Phase 2: Extend Verified Range (sequential pipeline per n)

**n=24 (highest priority — same type as n=22):**
1. `solver` runs informed SA for n=24
2. `validator` verifies the solution
3. `researcher` incorporates into pattern analysis

**n=25 (Paley — should be straightforward):**
1. `solver` implements GF(49) construction
2. `validator` verifies

**n=23 (hard — composite m):**
1. `researcher` analyzes Z_45 algebraic structure, sends constraints to solver
2. `solver` tries informed SA, CRT decomposition, dihedral Cayley
3. `validator` verifies (if found)

### Phase 3: General Proof (parallel with Phase 2)

| Task | Agent | Description |
|------|-------|-------------|
| Prime p≡3(4) proof | researcher | Formalize the Delta fluctuation argument: prove existence of good D11/D12 for all such primes |
| Composite construction | researcher | Develop construction family for composite m using CRT or extension fields |
| Proof outline | researcher | Write structured proof covering all cases, identifying proven vs open lemmas |
| Verify theory claims | validator | Computationally check conjectures for n≤30 |

### Phase 4: Synthesis

Team lead combines computational results (extended verified range) and theoretical progress (proof outline, new infinite families) into a coherent output.

---

## Coordination Patterns

1. **Solver → Validator pipeline:** Every claimed solution goes to `validator` before being reported as confirmed
2. **Researcher → Solver feedback:** New structural constraints discovered by `researcher` are sent to `solver` to narrow searches (this was the n=22 breakthrough pattern)
3. **Cross-validation:** Mathematical claims from `researcher` are sent to `validator` for computational checking — critical since all reasoning is AI-generated
4. **Background computation:** Long SA/SAT runs launched as background tasks by `solver`; solver works on other cases while waiting

---

## Verification

- Every construction verified by standalone brute-force validator (no shared code dependency)
- `test_core.py` run to confirm core library integrity
- Theoretical claims checked computationally where possible
- Solutions stored in JSON with provenance metadata
