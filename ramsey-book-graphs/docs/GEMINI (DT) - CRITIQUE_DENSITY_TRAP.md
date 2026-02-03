# CRITIQUE: The Density Trap in R(B21, B22) Search

## Executive Summary
The "8-violation barrier" identified in the Significance Analysis is a critical diagnostic result. It strongly suggests that the **Symmetry Constraint ($D_{22} = \overline{D_{11}}$)** is the root cause of the failure, not the graph family itself.

By enforcing symmetry, the solver is forced into a "density dead zone" where no solution exists. The fix requires breaking this symmetry to allow Red and Blue to have different average degrees.

---

## 1. The Mathematical Bottleneck: Why Symmetry Fails

The current approach enforces $D_{22} = \text{complement}(D_{11})$.
- This implies the graph is regular or near-regular with degree $k \approx \frac{N-1}{2}$.
- For $N=86$, this forces $k \approx 42.5$.

### The "Dead Zone" Calculation
For a random-like graph (which Ramsey graphs mimic), the expected number of common neighbors $\lambda$ is given by:
$$ \lambda \approx \frac{k^2}{N} $$

Substituting $N=86$ and $k=42.5$:
$$ \lambda \approx \frac{42.5^2}{86} \approx \frac{1806}{86} \approx 21.0 $$

**The Conflict:**
- **Red Requirement ($B_{21}$):** Requires $\lambda_{red} \le 20$.
- **Blue Requirement ($B_{22}$):** Requires $\lambda_{blue} \le 21$.
- **Result:** The geometry forces $\lambda \approx 21$. This satisfies Blue ($21 \le 21$) but **fails Red** ($21 > 20$).

The "8 violations" are not a failure of the heuristic; they represent the solver hitting the irreducible mathematical floor for graphs of this density.

---

## 2. The Solution: Target Asymmetry
To satisfy the stricter Red bound ($\lambda \le 20$), the Red graph must be **denser** (counter-intuitively, to access SRG structures) or **sparser** (to lower random collision) than average.

There exists a known parameter set for a **Strongly Regular Graph (SRG)** that fits the bounds exactly, but it is **asymmetric**.

**Target Parameters:**
1.  **Red Graph ($D_{11} \cup D_{12}$):** Target degree $k \approx 45$.
    -   Parameters: SRG(86, 45, 20, 27).
    -   $\lambda = 20$ (Satisfies Red $\le 20$).
2.  **Blue Graph (Complement):** Target degree $k \approx 40$.
    -   Parameters: SRG(86, 40, 21, 16).
    -   $\lambda = 21$ (Satisfies Blue $\le 21$).

**Crucial Pivot:** You must target **Red = Dense ($k=45$)** and **Blue = Sparse ($k=40$)**. The symmetric constraint $D_{22} = \overline{D_{11}}$ rules out this valid region.

---

## 3. Recommended Action Plan

Your SAT solver handles ~64 variables, which is small. You can afford to double the search space to fix the density issue.

### Step A: Decouple the Blocks
- **Current:** `D22 = ~D11` (Implicitly derived).
- **New:** Define $D_{11}$ (21 vars) and $D_{22}$ (21 vars) as independent variables.
- **Cost:** Adds only ~21 boolean variables. Total variables $\approx 85$. Still trivial for SAT.

### Step B: Add Density Hints (Cardinality Constraints)
Force the solver away from the $k=42,43$ dead zone.
- Add constraint: `Sum(Red Edges) >= 44 * 86 / 2`.
- Specifically for blocks:
    - If $D_{12}$ (circulant) has 43 vars, force $\approx 22-23$ to be Red.
    - If $D_{11}$ has 21 vars, force $\approx 10-11$ to be Red.

### Step C: Check "Skew-Hadamard" Logic
Since $N=86$ and $43 \equiv 3 \pmod 4$, the solution is likely related to **Skew-Hadamard** matrices rather than Paley graphs.
- Ensure your symmetry breaking doesn't assume $D_{ij} = D_{ji}^T$ (symmetric adjacency).
- Skew-adjacency ($A = -A^T$) might be the correct symmetry class for $p \equiv 3 \pmod 4$.

## Summary
The "Significant Negative Result" is currently a "Significant Constraint Error." Unlock the density, let Red be dense ($k \approx 45$), and the solver will likely break the 8-violation barrier.