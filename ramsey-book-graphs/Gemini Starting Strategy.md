Here is the technical specification document designed for Claude Code. It condenses the "2-Block Circulant" approach into an actionable implementation plan.

# ---

**Technical Specification: Ramsey Number Search Agent**

**Objective:**

Develop a computational pipeline to find lower bound constructions for the Book Ramsey Number $R(B\_{n-1}, B\_n)$.

**Target:** Prove $R(B\_{n-1}, B\_n) \\ge 4n \- 1$ for composite $n$ (specifically $n=18, 22, 27$).

**Method:** Search for "2-Block Circulant" graphs on $4n-2$ vertices that avoid Red $B\_{n-1}$ and Blue $B\_n$.

## **1\. Mathematical Data Structure**

The search space is too large for raw adjacency matrices ($2^{2415}$ for $n=18$). You must implement the **2-Block Circulant** structure, which reduces the problem to finding three sets of integers $D\_{11}, D\_{12}, D\_{22} \\subset \\mathbb{Z}\_m$.

### **1.1 Definitions**

* **Order:** $N \= 4n \- 2$.  
* **Block Size:** $m \= N / 2 \= 2n \- 1$.  
* **Vertices:** $V \= V\_1 \\cup V\_2$, where $V\_1 \= \\{0, \\dots, m-1\\}$ and $V\_2 \= \\{m, \\dots, 2m-1\\}$.  
* **Defining Sets:**  
  * $D\_{11} \\subseteq \\{1, \\dots, m-1\\}$ (symmetric: $x \\in D\_{11} \\iff \-x \\in D\_{11}$). Represents edges inside $V\_1$.  
  * $D\_{22} \\subseteq \\{1, \\dots, m-1\\}$ (symmetric). Represents edges inside $V\_2$.  
  * $D\_{12} \\subseteq \\{0, \\dots, m-1\\}$ (not necessarily symmetric). Represents edges between $V\_1$ and $V\_2$.

### **1.2 Adjacency Logic**

A graph $G$ is defined by $(D\_{11}, D\_{12}, D\_{22})$. Vertices $u, v$ are adjacent if:

1. $u, v \\in V\_1$ and $(v \- u) \\pmod m \\in D\_{11}$.  
2. $u, v \\in V\_2$ and $(v \- u) \\pmod m \\in D\_{22}$.  
3. $u \\in V\_1, v \\in V\_2$ and $(v \- u) \\pmod m \\in D\_{12}$.

## ---

**2\. Core Library Implementation (ramsey\_core.py)**

Implement a Python class BlockCirculantGraph with the following optimized methods.

### **2.1 Fast Common Neighbor Counting (The "Wesley Lemma")**

Do **not** build the full adjacency matrix and square it (too slow, $O(N^3)$). Use the arithmetic properties of the sets to count common neighbors in $O(m)$ time.

Implement these two helper functions:

* $\\Delta(A, B, d) \= |\\{(a, b) \\in A \\times B \\mid a \- b \\equiv d \\pmod m\\}|$  
* $\\Sigma(A, B, d) \= |\\{(a, b) \\in A \\times B \\mid a \+ b \\equiv d \\pmod m\\}|$

**The Count Formula:**

For any edge $(u, v)$, the number of common neighbors $\\lambda(u, v)$ is:

1. **If $u, v \\in V\_1$** (let $d \= v \- u$):  
   $$\\lambda \= \\Delta(D\_{11}, D\_{11}, d) \+ \\Delta(D\_{12}, D\_{12}, d)$$  
2. **If $u, v \\in V\_2$** (let $d \= v \- u$):  
   $$\\lambda \= \\Delta(D\_{22}, D\_{22}, d) \+ \\Delta(D\_{12}^T, D\_{12}^T, d)$$  
   *(Note: $D\_{12}^T$ corresponds to edges from $V\_2$ to $V\_1$. If using undirected graphs, ensure logic handles the symmetry).*  
3. **If $u \\in V\_1, v \\in V\_2$** (let $d \= v \- u$):  
   $$\\lambda \= \\Sigma(D\_{11}, D\_{12}, d) \+ \\Delta(D\_{12}, D\_{22}, d)$$

### **2.2 Verification**

* **Input:** $n$, sets $D\_{11}, D\_{12}, D\_{22}$.  
* **Logic:**  
  1. Iterate over all possible difference values $d \\in \\{1, \\dots, m-1\\}$ for intra-block edges.  
  2. Iterate over all $d \\in \\{0, \\dots, m-1\\}$ for inter-block edges.  
  3. Calculate $\\lambda$ using the formulas above.  
  4. **Constraint Check:**  
     * If edge is present (RED), require $\\lambda \< n \- 1$ (strictly less than $n-1$ common neighbors implies no $B\_{n-1}$).  
     * If edge is absent (BLUE), require $\\lambda\_{blue} \< n$. Note: $\\lambda\_{blue} \= (N \- 2\) \- \\lambda\_{red}$.

## ---

**3\. Solver Implementation (Two Engines)**

Create two separate solver scripts using the core library.

### **3.1 Engine A: SAT Solver (solver\_sat.py)**

Use pysat (PySAT) to find exact solutions for small $n$ ($n \\le 25$).

**Encoding:**

1. **Variables:** Create Boolean variables for each potential member of the sets.  
   * $x\_i$ for $i \\in \\{1, \\dots, \\lfloor m/2 \\rfloor\\}$ representing $D\_{11}$ (enforce symmetry).  
   * $y\_i$ for $i \\in \\{0, \\dots, m-1\\}$ representing $D\_{12}$.  
   * $z\_i$ for $i \\in \\{1, \\dots, \\lfloor m/2 \\rfloor\\}$ representing $D\_{22}$.  
2. **Constraints:**  
   * For every pair $(u, v)$, the "Book Constraint" is a cardinality constraint on the sum of active variables in the neighbor sets.  
   * Use pysat.card.CardEnc.atmost to encode: $\\sum (\\text{neighbors}) \\le k$.  
   * *Optimization:* Only encode constraints for "representative" edges ($d=1, \\dots$ in the cyclic group) due to symmetry.

### **3.2 Engine B: Tabu Search / Simulated Annealing (solver\_heuristic.py)**

Use this for larger $n$ or if SAT times out ($n=18$ might be on the boundary).

**Algorithm:**

1. **Initialize:** Random sets $D\_{11}, D\_{12}, D\_{22}$.  
2. **Objective Function:** $Cost \= \\sum\_{\\text{edges}} \\max(0, \\lambda\_{actual} \- \\lambda\_{allowed})$.  
3. **Move:** Flip membership of a random element $k$ in one of the sets.  
4. **Tabu Mechanism:** Keep a list of recent flips; do not reverse them for $T$ steps.  
5. **Restart:** If stuck in local minima for $X$ steps, do a "soft restart" (perturb 10% of bits).

## ---

**4\. Verification & Output Artifacts**

The agent must produce a final JSON artifact for any valid construction found:

JSON

{  
  "n": 18,  
  "target\_vertices": 70,  
  "parameters": {  
    "m": 35,  
    "D11": \[1, 3, 4,...\],  
    "D12": \[0, 2, 5,...\],  
    "D22": \[2, 6, 7,...\]  
  },  
  "verification": {  
    "max\_red\_book\_size": 17,  
    "max\_blue\_book\_size": 18,  
    "valid": true  
  }  
}

## **5\. Execution Strategy for Claude Code**

**Phase 1: Tooling**

1. Implement ramsey\_core.py with the fast counting formulas.  
2. Write a unit test with $n=6$ (known solution exists) to verify the formulas against a brute-force matrix check.

**Phase 2: Discovery ($n=18$)**

1. Run solver\_sat.py with glucose4. Set timeout to 10 minutes.  
2. If SAT fails, run solver\_heuristic.py with 10 parallel seeds.

**Phase 3: Analysis**

1. Once a graph is found, verify it using the ramsey\_core.py logic.  
2. Output the JSON structure.