# Algebraic Reduction of 2-Block Circulant Ramsey Book Graphs

---

## 1. Setup and Notation

Let $m \geq 3$ be odd. We consider a 2-coloring (red/blue) of the complete graph $K_{2m}$ using a **2-block circulant** structure:

- **Vertex set**: $V = V_1 \cup V_2$, where $V_1 = \{0, 1, \dots, m-1\}$ and $V_2 = \{0', 1', \dots, (m-1)'\}$.
- **Difference sets**: $D_{11}, D_{12}, D_{22} \subseteq \{0, 1, \dots, m-1\}$, with arithmetic in $\mathbb{Z}_m$.
- **Edge coloring**:
  - $V_1$–$V_1$: edge $\{i, j\}$ is **red** iff $(i - j) \bmod m \in D_{11}$
  - $V_1$–$V_2$: edge $\{i, j'\}$ is **red** iff $(i - j) \bmod m \in D_{12}$
  - $V_2$–$V_2$: edge $\{i', j'\}$ is **red** iff $(i - j) \bmod m \in D_{22}$
  - All other edges are **blue**.

The **book graph** $B_n = K_2 + \bar{K}_n$ consists of $n$ triangles sharing a common edge ($n + 2$ vertices). To show $R(B_{n-1}, B_n) \geq 4n - 1$, we exhibit a 2-coloring of $K_{4n-2}$ containing no red $B_{n-1}$ and no blue $B_n$. The validity conditions are:

> **(R)** Every red edge has at most $n - 2$ red common neighbors.
>
> **(B)** Every blue edge has at most $n - 1$ blue common neighbors.

We impose the following **structural constraints** on the difference sets, with $m = 2n - 1$:

| Label | Constraint | Consequence |
|-------|-----------|-------------|
| (S1) | $D_{11} \subseteq \{1, \dots, m-1\}$ is symmetric: $x \in D_{11} \iff -x \in D_{11}$ | Within-block adjacency is undirected |
| (S2) | $D_{22} = \{1, \dots, m-1\} \setminus D_{11}$ | $D_{22}$ is the complement of $D_{11}$ among nonzero elements |
| (S3) | $|D_{12}| = n - 1$ | Cross-block red degree is controlled |
| (S4) | $0 \in D_{12}$ | Convention (applicable when $m \equiv 3 \pmod{4}$) |

Note that $0 \notin D_{11}$ and $0 \notin D_{22}$ since vertices are not self-adjacent.

---

## 2. Autocorrelation and the Universal Symmetry Lemma

Define the **autocorrelation** of a set $S \subseteq \mathbb{Z}_m$ at shift $d$:

$$C(S, S;\, d) \;=\; \bigl|\{x \in S : x + d \in S\}\bigr|$$

We write:

$$A(d) = C(D_{11}, D_{11};\, d), \qquad B(d) = C(D_{12}, D_{12};\, d)$$

**Lemma (Universal Symmetry).** *For any set $S \subseteq \mathbb{Z}_m$ and any $d$, $C(S, S;\, d) = C(S, S;\, -d)$.*

*Proof.* Substitute $y = x + d$. Then $x = y - d$, and $\{x \in S : x + d \in S\} \to \{y - d : y \in S,\; y - d \in S\} = \{y \in S : y - d \in S\}$. Hence $C(S, S;\, d) = |\{y \in S : y - d \in S\}| = C(S, S;\, -d)$. ∎

**Consequence.** Both $A(d) = A(-d)$ and $B(d) = B(-d)$ hold, even though $D_{12}$ is not necessarily symmetric. This is critical: it means $B(-d) = B(d)$ throughout, which will cause constraint pairs to merge.

---

## 3. Theorem 1: V₁V₂ Auto-Satisfaction

**Under constraints (S1)–(S3), every $V_1$–$V_2$ edge automatically satisfies both validity conditions (R) and (B). Specifically:**

**(a)** For any $d \in \mathbb{Z}_m$, the $V_1$–$V_2$ edge at difference $d$ has exactly $|D_{12}| - \mathbb{1}[d \in D_{12}]$ red common neighbors.

**(b)** Consequently:
  - Every **red** $V_1$–$V_2$ edge ($d \in D_{12}$) has exactly $n - 2$ red common neighbors. ✓ (R)
  - Every **blue** $V_1$–$V_2$ edge ($d \notin D_{12}$) has exactly $n - 1$ blue common neighbors. ✓ (B)

**Therefore, the problem of constructing a valid 2-block circulant $R(B_{n-1}, B_n)$ lower bound reduces entirely to satisfying (R) and (B) on $V_1$–$V_1$ and $V_2$–$V_2$ edges.**

### Proof

**Part (a): Red common neighbor count.** Fix $u = a \in V_1$ and $v = b' \in V_2$ with $d = (a - b) \bmod m$. Red common neighbors come from two blocks.

*From $V_1$:* Vertex $w \in V_1$ ($w \neq u$) contributes iff $(w - a) \in D_{11}$ and $(w - b) \in D_{12}$. Setting $s = w - a$: require $s \in D_{11}$ and $s + d \in D_{12}$. Define:

$$\sigma_1(d) = \bigl|\{s \in D_{11} : s + d \in D_{12}\}\bigr|$$

(The exclusions $s \neq 0$ and $s \neq -d$ are automatic since $0 \notin D_{11}$.)

*From $V_2$:* Vertex $w' \in V_2$ ($w' \neq v$) contributes iff $(a - w) \in D_{12}$ and $(w - b) \in D_{22}$. Setting $t = a - w$: require $t \in D_{12}$ and $d - t \in D_{22}$. Define:

$$\sigma_2(d) = \bigl|\{t \in D_{12} : d - t \in D_{22}\}\bigr|$$

**Step 1: Expand $\sigma_2$ via the complement (S2).** Since $D_{22} = \{1, \dots, m{-}1\} \setminus D_{11}$, the condition $d - t \in D_{22}$ requires $d - t \neq 0$ (i.e., $t \neq d$) and $d - t \notin D_{11}$. Therefore:

$$\sigma_2(d) = \bigl|\{t \in D_{12} : t \neq d\}\bigr| - \bigl|\{t \in D_{12} : t \neq d,\; d - t \in D_{11}\}\bigr|$$

The first term is $|D_{12}| - \mathbb{1}[d \in D_{12}]$.

For the second term: since $0 \notin D_{11}$, the condition $d - t \in D_{11}$ already forces $d - t \neq 0$, hence $t \neq d$. The constraint $t \neq d$ is therefore redundant, giving:

$$\bigl|\{t \in D_{12} : t \neq d,\; d - t \in D_{11}\}\bigr| = \bigl|\{t \in D_{12} : d - t \in D_{11}\}\bigr|$$

**Step 2: Apply symmetry (S1).** Since $D_{11} = -D_{11}$, we have $d - t \in D_{11} \iff t - d \in D_{11}$. Substituting $s = t - d$:

$$\bigl|\{t \in D_{12} : d - t \in D_{11}\}\bigr| = \bigl|\{s \in D_{11} : s + d \in D_{12}\}\bigr| = \sigma_1(d)$$

**Step 3: Combine.**

$$\sigma_1(d) + \sigma_2(d) = \bigl(|D_{12}| - \mathbb{1}[d \in D_{12}]\bigr) - \sigma_1(d) + \sigma_1(d) = |D_{12}| - \mathbb{1}[d \in D_{12}]$$

$$\boxed{\text{Red common neighbors of } V_1\text{–}V_2 \text{ edge at } d \;=\; |D_{12}| - \mathbb{1}[d \in D_{12}]}$$

**Part (b): Validity of red edges.** If $d \in D_{12}$: red CN $= |D_{12}| - 1 = (n-1) - 1 = n - 2 \leq n - 2$. ✓ (R)

**Part (c): Validity of blue edges.** If $d \notin D_{12}$, we count blue common neighbors $\beta_1(d) + \beta_2(d)$.

*From $V_1$:* $\beta_1(d) = |D_{22}| - |\{s \in D_{22} : s + d \in D_{12}\}|$

*From $V_2$:* $\beta_2(d) = |D_{11}| - \sigma_1(d)$

The key observation: $\sigma_1(d) + |\{s \in D_{22} : s + d \in D_{12}\}|$ counts all nonzero $s$ with $s + d \in D_{12}$. Since $d \notin D_{12}$, the element $s = 0$ does not contribute ($0 + d = d \notin D_{12}$), so this sum equals $|\{s \in \mathbb{Z}_m : s + d \in D_{12}\}| = |D_{12}|$. Therefore $|\{s \in D_{22} : s + d \in D_{12}\}| = |D_{12}| - \sigma_1(d)$.

$$\beta_1(d) + \beta_2(d) = |D_{22}| - |D_{12}| + \sigma_1(d) + |D_{11}| - \sigma_1(d) = |D_{11}| + |D_{22}| - |D_{12}|$$

$$= (m - 1) - (n - 1) = n - 1 \leq n - 1 \quad \checkmark\text{ (B)}$$

This completes the proof. ∎

---

## 4. Theorem 2: V₂V₂ Algebraic Determination

### Lemma (Complement Autocorrelation)

*For $D_{22} = \{1, \dots, m-1\} \setminus D_{11}$ and any nonzero $d$:*

$$C(D_{22}, D_{22};\, d) = A(d) + (m - 2) - 2|D_{11}| + 2\cdot\mathbb{1}[d \in D_{11}]$$

*Proof.* Define $f(x) = \mathbb{1}[x \in D_{11}]$ and $g(x) = \mathbb{1}[x \in D_{22}] = \mathbb{1}[x \neq 0] - f(x)$. Expand:

$$C(D_{22}, D_{22};\, d) = \sum_{s} g(s)\,g(s+d) = T_1 - T_2 - T_3 + T_4$$

where:

$T_1 = \sum_s \mathbb{1}[s \neq 0]\cdot\mathbb{1}[s+d \neq 0]$. For $d \neq 0$, we exclude $s = 0$ and $s = -d$ (distinct), giving $T_1 = m - 2$.

$T_2 = \sum_s \mathbb{1}[s \neq 0]\cdot f(s+d)$. Substituting $t = s+d$: as $s$ ranges over $\mathbb{Z}_m \setminus \{0\}$, $t$ ranges over $\mathbb{Z}_m \setminus \{d\}$. So $T_2 = |D_{11}| - \mathbb{1}[d \in D_{11}]$.

$T_3 = \sum_s f(s)\cdot\mathbb{1}[s+d \neq 0] = |D_{11}| - f(-d) = |D_{11}| - \mathbb{1}[d \in D_{11}]$, using symmetry (S1).

$T_4 = \sum_s f(s)\,f(s+d) = A(d)$.

Combining: $C(D_{22}, D_{22};\, d) = (m-2) - 2(|D_{11}| - \mathbb{1}[d \in D_{11}]) + A(d)$. ∎

### Theorem

*For any nonzero $d \in \mathbb{Z}_m$, the red common neighbor count for a $V_2$–$V_2$ edge is:*

$$R_{22}(d) = A(d) + B(d) + (m - 2 - 2|D_{11}|) + 2\cdot\mathbb{1}[d \in D_{11}]$$

*Proof.* Fix $u = a', v = b' \in V_2$ with $d = (a - b) \bmod m$.

*Red CN from $V_1$:* Vertex $w \in V_1$ contributes iff $(w - a) \in D_{12}$ and $(w - b) \in D_{12}$. Setting $s = w - a$: require $s \in D_{12}$ and $s + d \in D_{12}$.

$$\text{Count} = C(D_{12}, D_{12};\, d) = B(d)$$

*Red CN from $V_2$:* Vertex $w' \in V_2$ ($w' \neq u, v$) contributes iff $(w - a) \in D_{22}$ and $(w - b) \in D_{22}$. Setting $s = w - a$: require $s \in D_{22}$ and $s + d \in D_{22}$. (Exclusions $s \neq 0$ and $s + d \neq 0$ are automatic since $0 \notin D_{22}$.)

$$\text{Count} = C(D_{22}, D_{22};\, d)$$

Summing and applying the Complement Autocorrelation Lemma yields the result. ∎

---

## 5. V₁V₁ Common Neighbor Counts

For a **red** $V_1$–$V_1$ edge ($d \in D_{11}$):

*Red CN from $V_1$:* Vertex $w \in V_1$ ($w \neq u, v$) contributes iff $s = w - a \in D_{11}$ and $s + d \in D_{11}$. (Exclusions automatic since $0 \notin D_{11}$.)

$$\text{Count} = A(d)$$

*Red CN from $V_2$:* Vertex $w' \in V_2$ contributes iff $(a - w) \in D_{12}$ and $(b - w) \in D_{12}$. Setting $t = a - w$: require $t \in D_{12}$ and $t - d \in D_{12}$.

$$\text{Count} = C(D_{12}, D_{12};\, -d) = B(-d) = B(d)$$

where the last equality is by the Universal Symmetry Lemma.

$$\boxed{R_{11}(d) = A(d) + B(d)}$$

---

## 6. Theorem 3: The Complete Structural Reduction

**Under constraints (S1)–(S3), the four validity constraint classes (red/blue on $V_1V_1$ and $V_2V_2$) collapse into two inequalities governing $A(d) + B(d)$, indexed by whether $d \in D_{11}$ or $d \in D_{22}$.**

**Theorem.** *A 2-block circulant construction satisfying (S1)–(S3) is valid if and only if for all $d \in \{1, \dots, m-1\}$:*

$$\boxed{A(d) + B(d) \;\leq\; \begin{cases} n - 2 & \text{if } d \in D_{11} \\[4pt] 2|D_{11}| - n + 1 & \text{if } d \in D_{22} \end{cases}}$$

### Proof

We verify that the four constraint classes yield exactly these two conditions.

**Case 1: Red $V_1V_1$ ($d \in D_{11}$).** Requirement: $R_{11}(d) \leq n - 2$.

$$A(d) + B(d) \leq n - 2$$

**Case 2: Blue $V_2V_2$ ($d \in D_{11}$).** The edge at difference $d$ is blue in $V_2$ (since $d \in D_{11}$ means $d \notin D_{22}$). We count blue common neighbors.

*Blue CN from $V_1$:* Vertex $w$ contributes iff $(w - a) \notin D_{12}$ and $(w - b) \notin D_{12}$. By inclusion–exclusion:

$$\text{Count} = m - 2|D_{12}| + B(d) = (2n-1) - 2(n-1) + B(d) = 1 + B(d)$$

*Blue CN from $V_2$:* Vertex $w'$ ($w' \neq u, v$) contributes iff $(w - a) \notin D_{22}$ and $(w - b) \notin D_{22}$. For $w \neq a, b$: both $(w-a)$ and $(w-b)$ are nonzero, so $\notin D_{22}$ means $\in D_{11}$.

$$\text{Count} = A(d)$$

Total: $1 + A(d) + B(d) \leq n - 1$, giving $A(d) + B(d) \leq n - 2$. **Identical to Case 1.** ✓

**Case 3: Red $V_2V_2$ ($d \in D_{22}$).** Requirement: $R_{22}(d) \leq n - 2$. Since $d \in D_{22}$ means $d \notin D_{11}$, the indicator $\mathbb{1}[d \in D_{11}] = 0$. Substituting from Theorem 2 with $m = 2n - 1$ and $k = |D_{11}|$:

$$A(d) + B(d) + (2n - 3 - 2k) \leq n - 2$$

$$A(d) + B(d) \leq 2k - n + 1$$

**Case 4: Blue $V_1V_1$ ($d \in D_{22}$).** Blue common neighbors:

*From $V_1$:* $C(D_{22}, D_{22};\, d) = A(d) + (2n - 3 - 2k)$ (with $\mathbb{1}[d \in D_{11}] = 0$).

*From $V_2$:* $m - 2|D_{12}| + C(D_{12}, D_{12};\, -d) = 1 + B(-d) = 1 + B(d)$.

Total: $A(d) + B(d) + (2n - 2 - 2k) \leq n - 1$, giving $A(d) + B(d) \leq 2k - n + 1$. **Identical to Case 3.** ✓

This completes the proof. ∎

---

## 7. Structural Consequences

### 7.1 The reduction to additive combinatorics

The search for valid constructions is now entirely a problem of finding:

1. A **symmetric** set $D_{11} \subseteq \{1, \dots, m-1\}$ (with $|D_{11}| = k$ even)
2. A set $D_{12} \subseteq \mathbb{Z}_m$ (with $|D_{12}| = n - 1$)

such that the combined autocorrelation $A(d) + B(d)$ respects a "waterlevel" threshold that toggles based on membership in $D_{11}$. This connects directly to the theory of almost-difference sets and low-autocorrelation sequences.

### 7.2 Extremal tightness

Every verified construction satisfies the $D_{11}$ constraint with equality: $\max_{d \in D_{11}} [A(d) + B(d)] = n - 2$. The V₁V₂ edges also saturate their thresholds exactly. This zero-slack property suggests constructions are essentially unique up to automorphism.

### 7.3 Constraint asymmetry and the Paley case

When $|D_{11}| = n - 1$ (as in the Paley construction, where $D_{11} = \text{QR}(\mathbb{F}_q^*)$), the $D_{22}$ threshold becomes $2(n-1) - n + 1 = n - 1$, which is strictly larger than the $D_{11}$ threshold of $n - 2$. The binding constraints are therefore always those on $D_{11}$ differences. The $D_{22}$ constraints provide a unit of slack — a structural buffer built into the complement architecture.

### 7.4 Algebraic core

The entire reduction rests on three facts:

1. **Complement partition**: $D_{11} \sqcup D_{22} = \{1, \dots, m{-}1\}$ converts summations over $D_{22}$ to summations over all nonzero elements minus $D_{11}$.

2. **Symmetry of $D_{11}$**: $D_{11} = -D_{11}$ enables the cross-correlation identity that drives the V₁V₂ cancellation.

3. **Universal autocorrelation symmetry**: $C(S, S;\, d) = C(S, S;\, -d)$ for any set $S$, which causes the red/blue constraint pairs on the same edge type to merge, halving the system from four independent conditions to two.
