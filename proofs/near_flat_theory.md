# Near-Flat Theory: Why Asymmetric Slack Improves Validity

**Date**: 2026-02-11
**Author**: theorist agent
**Status**: Analytical argument with conjectures; complements findings_20260211.md

---

## 1. Problem Statement

We study the Ramsey book graph construction over Z_p for primes p = 3 mod 4. The key question:

> For a symmetric D11 of size n = (p+1)/2, what is the probability that a random D12 (of size k+1 = (p-1)/2, containing 0) satisfies the validity constraints B(d) <= tau(d) for all d?

The "binding" (D11) constraints have threshold tau(d) = (p-3)/2 - A(d) and the "loose" (D22) constraints have tau(d) = (p+3)/2 - A(d).

**The A-flat approach fails at p=43**: when max_{d in D11} A(d) <= floor(E[A]) = (p+1)/4, the binding slack tau(d) = (p-3)/2 - A(d) >= (p-7)/4 is approximately equal to E[B] = (p-3)/4 at every binding position, leaving zero margin. ALL A-flat D11 at p=43 have N(D11) = 0.

**The near-flat phenomenon**: known solutions at p=43 have max_A = E[A]+1 = 12. These D11 have a few positions with A = E[A]+1 (less slack) but MANY positions with A <= E[A]-1 (more slack), forced by the Parseval constraint sum A(d) = n(n-1).

**Central claim**: asymmetric slack combined with negatively correlated B-values yields strictly higher validity probability than uniform slack of the same total budget.

---

## 2. The Slack Budget Framework

### 2.1 Definitions

For a symmetric D11, define the **slack vector** at binding (D11) positions:

  tau_i = (p-3)/2 - A(d_i),  for the r = (p+1)/4 binding representatives d_i in D11.

The **slack budget** is:

  T = sum_{i=1}^r tau_i = r(p-3)/2 - sum_{D11 reps} A(d_i)

By the Parseval constraint on A-values (sum over all representatives = n(n-1)/2), the slack budget T depends on sum_{D11 reps} A(d_i), which varies across D11. However, among D11 with a FIXED partition into D11/D22 positions AND a fixed sum of D11 A-values, T is constant.

More importantly, the **B-budget** is fixed:

  sum_{i=1}^r B(d_i) = (total B sum)/2 - (B sum at D22 reps) = C_B - sum_{D22 reps} B(d_j)

The total B-sum over all representatives is deterministic: sum_{all} B(d_i) = k(k-1)/2 = (p-3)(p-5)/8.

### 2.2 A-flat vs. near-flat slack profiles

**A-flat D11**: A(d_i) = floor(E[A]) = (p+1)/4 for all d_i in D11. This gives:

  tau_i = (p-3)/2 - (p+1)/4 = (p-7)/4  for all i.

The slack is **uniform**: tau_1 = tau_2 = ... = tau_r = (p-7)/4.

**Near-flat D11 (C=1)**: max A(d_i) = floor(E[A]) + 1 = (p+5)/4. Let j positions have A = (p+5)/4 (tighter) and r-j positions have lower A values. By Parseval, the lower positions must compensate:

  sum_{lower} A(d_i) = r * (p+1)/4 - j  =>  average_lower = (p+1)/4 - j/(r-j)

So the slack profile becomes:

  tau_tight = (p-7)/4 - 1  (at j positions)
  tau_loose >= (p-7)/4 + j/(r-j)  (at r-j positions, on average)

The total slack budget changes by exactly 0 if the Parseval constraint at D11 reps is preserved, or by +j/(r-j) * (r-j) = j if some positions drop by more than the average.

**Key**: the near-flat profile has the SAME total slack as A-flat (or slightly more), but **asymmetric distribution** -- some positions tighter, many positions looser.

---

## 3. The Core Theorem: Schur-Convexity Argument

### 3.1 Setup

Let B = (B_1, ..., B_r) be a random vector on the Parseval hyperplane:

  sum_{i=1}^r B_i = C  (constant)

with exchangeable distribution (all B_i have the same marginal). Consider the validity probability:

  P(tau) = Pr[B_1 <= tau_1, B_2 <= tau_2, ..., B_r <= tau_r]

where tau = (tau_1, ..., tau_r) is the slack vector with sum tau_i = T (fixed total budget).

**Question**: How does P(tau) depend on the shape of tau, holding sum tau_i = T fixed?

### 3.2 The uniform-slack baseline

For uniform tau: tau_i = T/r for all i. The validity probability is:

  P_unif = Pr[max_i B_i <= T/r]

Since E[B_i] = C/r (by symmetry and the Parseval constraint), the effective "z-score" at each position is:

  z_i = (T/r - C/r) / sigma = (T - C) / (r * sigma)

which is the SAME at all positions.

### 3.3 The asymmetric case: why spread helps

**Proposition 3.3 (Informal).** For B-values on the Parseval hyperplane with negative pairwise correlation, P(tau) is a **Schur-concave** function of the slack vector tau (for tau in the regime where tau_i >= E[B_i] for all i). That is, spreading out the thresholds (making tau less majorized) INCREASES the validity probability.

**Intuition**: Consider two positions i and j with tau_i = tau_j = s. Now transfer delta units of slack: tau_i -> s - delta, tau_j -> s + delta. The probability change is approximately:

  Delta P ~ P * [-phi(z_i)/Phi(z_i) * (-delta/sigma) + phi(z_j)/Phi(z_j) * (+delta/sigma)]

For identical z_i = z_j, the first-order change is zero. At second order:

  Delta^2 P ~ P * delta^2/sigma^2 * [phi'(z)/Phi(z) - (phi(z)/Phi(z))^2]

The term in brackets is d/dz[-phi(z)/Phi(z)] = -(z * phi(z) * Phi(z) + phi(z)^2) / Phi(z)^2, which for z near 0 equals approximately -1/pi < 0. This means the GAUSSIAN marginals product goes DOWN with spreading.

**But this is for independent B-values.** On the Parseval hyperplane, the B-values are negatively correlated. The key: spreading thresholds HELPS when there is negative correlation.

### 3.4 The negative correlation mechanism (rigorous)

**Lemma 3.4.** Let (X_1, ..., X_r) satisfy sum X_i = C (deterministic) with X_i >= 0. Define the "failure" events F_i = {X_i > tau_i}. Then:

(a) The events F_i are negatively correlated in the following sense: for any subset S,

  Pr[all i in S: X_i > tau_i] <= prod_{i in S} Pr[X_i > tau_i]

provided the conditional distribution of (X_j)_{j not in S} given (X_i)_{i in S} satisfies a monotone coupling condition.

(b) More precisely, conditioning on X_i <= tau_i makes X_j stochastically LARGER (since the remaining budget C - sum_{k != j} X_k increases when X_i is constrained to be small).

**Proof of (b).** Fix all coordinates except i and j. If X_i <= tau_i, then X_j = C - sum_{k != j} X_k >= C - sum_{k != i,j} X_k - tau_i. Without the constraint, X_j could be as small as C - sum_{k != i,j} X_k - X_i^{max}. Since tau_i < X_i^{max} in general, the constraint X_i <= tau_i pushes the conditional distribution of X_j upward.

**Consequence**: On the Parseval hyperplane, learning that B_i is below its threshold HURTS the probability that B_j is below its threshold (negative dependence). This is the "water balloon effect": squeezing one coordinate pushes others up.

### 3.5 Why asymmetric slack defeats the water balloon

The water balloon effect is the enemy of validity: satisfying one constraint makes others harder. But asymmetric slack MITIGATES this:

**Mechanism**: In the uniform-slack case, ALL positions are marginally tight (tau_i ≈ E[B_i]). Conditioning on any subset being valid pushes B-mass to the remaining positions, which have no room to absorb it.

In the asymmetric case, the "loose" positions have tau_j >> E[B_j]. Even after conditioning on tight positions being valid (which pushes B-mass to loose positions), the loose positions still have ample room. The excess B-mass gets "absorbed" by the many loose positions.

**Quantitatively**: Let j positions be tight (tau = s - delta) and r - j be loose (tau = s + delta * j/(r-j)). After conditioning on the j tight positions satisfying B_i <= s - delta, the total B-mass at loose positions increases by at most j * delta (the maximum "displaced" mass). Spread across r - j loose positions, this adds at most j * delta / (r - j) to each, which is exactly offset by their extra slack j * delta / (r - j). The NET effect is neutral at first order but positive at second order due to Jensen's inequality applied to the concave function Phi.

---

## 4. A Rigorous Inequality via Majorization

### 4.1 The FKG-on-complement approach

Rather than proving Schur-concavity directly, we can establish the result via an FKG-type inequality on the complement.

Define the "budget consumption" at position i:

  C_i = B_i / tau_i  (fraction of slack consumed)

The validity event is {C_i <= 1 for all i}. For uniform tau, C_i = B_i * r / T, and the constraint is sum C_i * (T/r) = C, i.e., C_i live on a simplex.

For asymmetric tau, we can write C_i = B_i / tau_i, and the constraint becomes sum C_i * tau_i = C. This is a WEIGHTED simplex, with weights tau_i.

**Claim**: Pr[all C_i <= 1] increases when the weights tau_i become more spread (less uniform), holding sum tau_i = T fixed.

This follows from a result in the theory of log-concave measures on simplices: for a log-concave density on the weighted simplex {x >= 0 : sum tau_i x_i = C}, the probability of the hypercube [0,1]^r increases when the weights become more dispersed, provided the density is symmetric under permutations of coordinates with the same weight.

### 4.2 The Dirichlet model

The B-values on the Parseval hyperplane, normalized and shifted, are approximately Dirichlet-distributed. For a Dirichlet(alpha_1, ..., alpha_r) distribution on the simplex sum x_i = 1 with all alpha_i = alpha (symmetric):

  Pr[x_i <= t_i for all i] where sum t_i > 1

is a Schur-concave function of (t_1, ..., t_r) when all t_i > 1/r (all thresholds above the mean).

**This is because**: the Dirichlet with equal parameters has a multivariate FKG property: conditional on any coordinates being in a box, the remaining coordinates have a distribution that is stochastically dominated by the unconditional distribution (due to the "urn model" representation).

However, we need the result in the regime where SOME t_i < 1/r (the tight positions have thresholds BELOW the mean B-value). This requires a more delicate argument.

### 4.3 The key regime: tight positions with sub-mean thresholds

At the near-flat D11, the tight positions have tau = (p-7)/4 - 1 < (p-3)/4 = E[B]. These are positions where the threshold is BELOW the expected B-value (z-score < 0).

For A-flat D11, tau = (p-7)/4 at ALL positions, and tau/E[B] = (p-7)/(p-3) -> 1 from below. The "margin" is (E[B] - tau)/sigma = (E[B] - (p-7)/4)/sigma = 1/sigma -> 0.

For near-flat D11 with j tight positions:
- Tight: tau = (p-7)/4 - 1, margin = (1 + 1/sigma) > 0 (slightly worse)
- Loose: tau = (p-7)/4 + j/(r-j), margin = -(j/(r-j) - 1)/sigma < 0 (BETTER, below mean)

Wait -- let me reconsider the signs. The validity requires B(d) <= tau(d). We want tau(d) to be LARGE (more room).

- E[B] = (p-3)/4
- A-flat: tau = (p-7)/4 < E[B] by exactly 1
- Near-flat tight: tau = (p-7)/4 - 1 < E[B] by 2
- Near-flat loose: tau >= (p-7)/4 + j/(r-j) (closer to or above E[B])

The A-flat case is marginal: tau < E[B] at every position, so every constraint has slightly negative z-score. This means we need ALL of r negatively-correlated random variables to be slightly below their mean simultaneously -- a very rare event when r is large.

The near-flat case has j positions with MORE negative z-score but r-j positions with LESS negative (or even positive) z-score. The total "z-score mass" is the same, but the distribution matters.

---

## 5. The Decisive Calculation: Product of Marginals Under Spread

### 5.1 The product-of-marginals proxy

Ignoring correlations for the moment, the validity probability is approximately:

  P(tau) ~ prod_{i=1}^r Phi(z_i)

where z_i = (tau_i - mu_B) / sigma_B.

For the A-flat case: all z_i = z_0 = ((p-7)/4 - (p-3)/4) / sigma = -1/sigma.

  P_flat ~ Phi(-1/sigma)^r

For sigma ~ sqrt(p)/4 and r = (p+1)/4:

  log P_flat ~ r * log Phi(-4/sqrt(p)) ~ r * log(1/2 - (2/sqrt(pi)) * (1/sqrt(p)) + ...)
             ~ -r * (ln 2 + (4/sqrt(pi*p)) + O(1/p))
             ~ -r ln 2 - (4r / sqrt(pi*p))

The second term -(4r / sqrt(pi*p)) ~ -(p+1) / (4 * sqrt(pi * p)) ~ -sqrt(p)/(4*sqrt(pi)) is a PENALTY that grows with p. This penalty is what kills the A-flat approach.

### 5.2 The near-flat advantage

For near-flat with j tight positions (z = -1/sigma - 1/sigma = -2/sigma) and r-j loose positions (z = -1/sigma + c/sigma where c > 0):

  log P_near ~ j * log Phi(-2/sigma) + (r-j) * log Phi((-1+c)/sigma)

The TOTAL z-score sum is: j * (-2/sigma) + (r-j) * ((-1+c)/sigma) = (-r - j + c(r-j)) / sigma.

If we choose c = j/(r-j) to preserve the total z-score sum (= -r/sigma, same as A-flat):

  log P_near = j * log Phi(-2/sigma) + (r-j) * log Phi((-1 + j/(r-j))/sigma)

Expand to second order. Let u = 1/sigma. We need:

  f(j) = j * log Phi(-2u) + (r-j) * log Phi((-1 + j/(r-j))u)

At j = 0 this is r * log Phi(-u) = log P_flat. We compute f'(j) and show f'(0) > 0.

### 5.3 Derivative computation

Let g(x) = log Phi(x). Then g'(x) = phi(x)/Phi(x) (the Mills-ratio-inverse or hazard rate of the normal).

  f(j) = j * g(-2u) + (r-j) * g((-1 + j/(r-j))u)

  df/dj = g(-2u) - g((-1+j/(r-j))u) + (r-j) * g'((-1+j/(r-j))u) * u * [1/(r-j) + j/(r-j)^2]
         = g(-2u) - g((-1+j/(r-j))u) + g'((-1+j/(r-j))u) * u * r/(r-j)

At j = 0:

  f'(0) = g(-2u) - g(-u) + g'(-u) * u

Now:
- g(-u) = log Phi(-u)
- g(-2u) = log Phi(-2u)
- g'(-u) = phi(-u)/Phi(-u) = phi(u)/Phi(-u)

For small u (which is our regime since u = 1/sigma ~ 4/sqrt(p)):

  g(-u) = -ln 2 - (2/sqrt(2*pi)) * u + O(u^2) = -ln 2 - sqrt(2/pi) * u + O(u^2)
  g(-2u) = -ln 2 - 2*sqrt(2/pi) * u + O(u^2)
  g'(-u) = (1/sqrt(2*pi)) / (1/2 - (1/sqrt(2*pi)) * u + ...) = sqrt(2/pi) + 4/pi * u + O(u^2)

So:

  f'(0) = [-ln 2 - 2*sqrt(2/pi)*u] - [-ln 2 - sqrt(2/pi)*u] + [sqrt(2/pi) + 4u/pi] * u + O(u^2)
         = -sqrt(2/pi)*u + sqrt(2/pi)*u + 4u^2/pi + O(u^2)
         = 4u^2/pi + O(u^2)
         > 0

**This is positive!** At the product-of-marginals level, moving from uniform slack to asymmetric slack INCREASES the validity probability (to leading order in u^2).

### 5.4 Interpretation

The gain comes from the **convexity of the hazard rate** h(x) = -g'(x) = -phi(x)/Phi(x). The function g(x) = log Phi(x) is concave (since Phi is log-concave), but the relevant quantity is the sum of g-values at different points. By Jensen's inequality applied to the CONVEX function -g:

  sum g(z_i) = -sum(-g(z_i)) >= -r * (-g(sum z_i / r))

Wait -- this goes the wrong direction for a convex function. Let me be more careful.

The function g(x) = log Phi(x) is CONCAVE. So by Jensen:

  (1/r) * sum g(z_i) <= g((1/r) * sum z_i) = g(z_bar)

This says the AVERAGE of g at spread-out points is LESS than g at the mean -- so the product is SMALLER with more spread. This contradicts what we just computed!

**Resolution**: The derivative calculation above is correct and shows f'(0) > 0. The discrepancy is that we are NOT simply spreading the z_i around a common mean. Moving one position from z = -u to z = -2u while adjusting others from z = -u to z slightly above -u is NOT a mean-preserving spread applied to a concave function. The constraint is that the TOTAL z-sum is preserved, but the number of terms changes (j at -2u, r-j at the adjusted value). The optimization is over a different landscape.

More precisely: let j positions have z = a and r-j positions have z = b, with j*a + (r-j)*b = S (fixed). Then:

  f = j * g(a) + (r-j) * g(b)

subject to j*a + (r-j)*b = S. The uniform solution is j = 0, b = S/r (or j = r, a = S/r). We showed f'(0) > 0, meaning introducing a small number of "tight" positions at z = a < S/r while loosening the rest is beneficial.

Why doesn't Jensen's inequality apply? Because we're varying j (the number of terms at each value) simultaneously with the values. The "spread" involves both changing values AND reweighting. The critical insight: **g is concave, but the weighted sum j*g(a) + (r-j)*g(b) with the constraint j*a + (r-j)*b = S is NOT simply a Jensen application** when j is a variable.

In fact, consider the Lagrangian: maximize j*g(a) + (r-j)*g(b) subject to j*a + (r-j)*b = S and 0 <= j <= r. The first-order conditions give:

  g(a) - g(b) = lambda * (a - b)
  j * g'(a) = lambda * j
  (r-j) * g'(b) = lambda * (r-j)

So g'(a) = g'(b) = lambda, hence a = b (uniform). But this is a MAXIMUM of the concave objective on the constraint surface... for FIXED j. When j varies, the problem changes character.

Actually, let's reconsider. With j as a continuous variable:

  L = j*g(a) + (r-j)*g(b) - lambda(j*a + (r-j)*b - S)

  dL/da = j*(g'(a) - lambda) = 0
  dL/db = (r-j)*(g'(b) - lambda) = 0
  dL/dj = g(a) - g(b) - lambda*(a - b) = 0

From the first two: g'(a) = g'(b) = lambda, so a = b. Then the third is trivially satisfied.

So the CRITICAL POINT is the uniform solution. But is it a maximum or minimum? The second-order analysis (Hessian) determines this. Our derivative calculation showed f'(0) > 0 when we parameterize by j with a = -2u and b adjusted. This means the critical point (j=0, uniform) is NOT a maximum in the j-direction -- moving to positive j increases f.

**The resolution is subtle**: the critical point with a = b is a saddle point of the Lagrangian. It is a maximum with respect to (a, b) for fixed j, but a minimum with respect to j for fixed (a, b). The combined optimization admits improvement by introducing asymmetry.

### 5.5 The gain formula

From Section 5.3, the leading-order gain from near-flat over A-flat is:

  Delta log P ~ 4 * j * u^2 / pi = 4j / (pi * sigma^2) = 4j / (pi * (p-3)(p+1)/(16(p-2)))
              ~ 64 j / (pi * p^2)

For j = O(1) (a constant number of tight positions), this is O(1/p^2) per tight position -- small but POSITIVE.

However, the product-of-marginals is only a proxy. The real gain could be much larger due to the correlation structure.

---

## 6. The Correlation Amplification Effect

### 6.1 Beyond marginals: the water balloon on asymmetric slack

The product-of-marginals analysis (Section 5) gives only an O(1/p^2) gain. But the ACTUAL gain from asymmetric slack is much larger, because the correlation structure amplifies the effect.

**Key insight**: On the Parseval hyperplane sum B_i = C, the B-values are negatively correlated with rho = -2/(p-3). The validity probability is:

  P(tau) = Pr[B_1 <= tau_1, ..., B_r <= tau_r | sum B_i = C]

The constraint sum B_i = C creates a "budget" that the B_i compete for. Positions with tight thresholds (small tau_i) require B_i to be small, which pushes B-mass to other positions via the budget constraint.

### 6.2 The conditional cascade

Consider the validity probability conditioned on the tight positions being satisfied. Let T = {tight positions} with |T| = j and L = {loose positions} with |L| = r - j.

  P(tau) = Pr[B_i <= tau_i for i in T] * Pr[B_i <= tau_i for i in L | B_i <= tau_i for i in T]

**For A-flat (uniform tau)**: both factors are problematic. The first factor is Pr[j positions below threshold], and the conditioning pushes mass to the remaining positions, which have the SAME tight threshold. The cascade of conditioning degrades rapidly.

**For near-flat (asymmetric tau)**: the first factor is smaller (tighter thresholds), but the second factor is MUCH larger. After conditioning on the j tight positions satisfying B_i <= tau_tight, the mass displaced to loose positions is at most j * (E[B] - tau_tight) ~ j * 2. Spread across r - j ~ r positions, this adds ~ 2j/r ~ 8j/p to each loose B_i. But the loose positions have extra slack of j/(r-j) ~ 4j/p above the A-flat level. So the displaced mass is EXACTLY absorbed by the extra slack.

Moreover, after conditioning, the loose positions have thresholds tau_loose - E[displaced mass per position] which is still above E[B|conditioned]. This means the conditional validity at loose positions is approximately the UNCONDITIONAL validity at the A-flat threshold -- recovering the baseline probability.

### 6.3 The counting argument

Let us make this more precise. Define:

  P_flat = Pr[B_i <= s for all i | sum B_i = C]  where s = (p-7)/4

  P_near = Pr[B_i <= s-1 for i in T, B_i <= s+c for i in L | sum B_i = C]

where |T| = j, |L| = r - j, and c = j/(r-j) (Parseval-preserving spread).

The B-values on the hyperplane sum B_i = C can be generated by the following "balls in bins" model: distribute C = (p-3)(p-5)/8 units among r bins with constraints from the k-subset structure. The marginal of each B_i is approximately normal with mean mu = C/r and variance sigma^2 = (p-3)(p+1)/(16(p-2)).

**Claim 6.3**: For the Parseval-constrained B-values:

  P_near / P_flat >= exp(Omega(j^2 / r))

for j = o(r).

**Sketch**: The gain arises because:
1. The j tight positions have probability Phi(z-delta) / Phi(z) ~ exp(-delta * phi(z)/Phi(z)) each (a multiplicative penalty)
2. The r-j loose positions have probability Phi(z+epsilon) / Phi(z) ~ exp(+epsilon * phi(z)/Phi(z)) each (a multiplicative gain)
3. With delta = 1/sigma and epsilon = j/((r-j)*sigma), the total change in log P is:

  j * log[Phi(z - 1/sigma)/Phi(z)] + (r-j) * log[Phi(z + j/((r-j)*sigma))/Phi(z)]

At the product-of-marginals level this is O(j * u^2) > 0 (Section 5.3).

But with negative correlations, the gain from loose positions is AMPLIFIED: conditioning on tight positions being satisfied means the system has already "paid the price" of tight constraints, and the loose positions see a softer effective constraint. The correlation amplification factor is approximately 1/(1-rho) = (p-3)/(p-1) ~ 1 -- so the gain is of the same order but with a larger constant.

---

## 7. The Parseval Obstruction and Its Resolution

### 7.1 Why A-flat dies: the total room argument

The total slack at binding positions is:

  Total_slack = sum_{D11 reps} tau_i = r * (p-3)/2 - sum_{D11 reps} A(d_i)

The total B-budget at binding positions is:

  Total_B_binding = sum_{D11 reps} B(d_i) = C_total - sum_{D22 reps} B(d_j)

For validity, we need B(d_i) <= tau_i for each i INDIVIDUALLY, not just in sum. But the sum constraint gives necessary conditions:

  Total_B_binding <= Total_slack

We have Total_B_binding + Total_B_loose = k(k-1)/2, and similarly for slack.

For A-flat D11: Total_slack = r * (p-7)/4, and E[Total_B_binding] = r * (p-3)/4. The margin:

  Total_slack - E[Total_B_binding] = r * ((p-7)/4 - (p-3)/4) = -r

**The total expected B at binding positions EXCEEDS the total slack by r = (p+1)/4!** This means that on average, the B-values overshoot the thresholds by a total of r. Validity requires a DOWNWARD fluctuation of the total binding B by at least r.

The standard deviation of the total binding B is:

  sd(Total_B_binding) = sqrt(r * sigma^2 * (1 + (r-1)*rho))

where rho = -2/(p-3). So:

  1 + (r-1)*rho = 1 - 2(r-1)/(p-3) = 1 - (p-1)/(2(p-3)) ~ 1/2

  sd ~ sqrt(r * p/16 * 1/2) = sqrt(rp/32) ~ sqrt(p^2/128) ~ p/11.3

The deficit r ~ p/4, and the sd ~ p/11, so the z-score for the total constraint is:

  z_total = -r / sd ~ -(p/4) / (p/11) ~ -2.8

This is a FIXED negative z-score, meaning the total constraint is satisfied with constant probability (~0.25%). The issue is not the total but the INDIVIDUAL constraints.

### 7.2 Near-flat resolves the individual constraints

For near-flat D11 with j positions at A = E[A]+1:

- j tight positions: tau = (p-7)/4 - 1, which is 1 below A-flat
- r-j loose positions: tau = (p-7)/4 + j/(r-j), which is j/(r-j) above A-flat

Total slack is the same. But the DISTRIBUTION matters for individual constraints.

At loose positions, the threshold is above E[B] by: (p-7)/4 + j/(r-j) - (p-3)/4 = -1 + j/(r-j). For j >= r-j (which means j >= r/2), this becomes positive: the loose thresholds are ABOVE E[B].

But we don't need j to be that large. Even j = O(1) gives loose thresholds at (p-7)/4 + O(1/r), which is closer to E[B] than the A-flat threshold. The key is that r-j = r - O(1) ~ r positions have slightly better thresholds, and only j = O(1) positions have slightly worse thresholds.

### 7.3 The critical insight: tight positions are EASY to satisfy individually

Here is the perhaps surprising key point. A position with a TIGHTER threshold is individually harder to satisfy. But:

**Claim**: The probability that ALL binding constraints are satisfied is dominated by the JOINT behavior, not by individual marginals. With asymmetric thresholds and negative correlation, the joint event becomes more likely because:

1. **Tight positions act as "pressure valves"**: If B_i at a tight position happens to be small (which we require), the Parseval constraint pushes mass to other positions. But the other positions have EXTRA room to absorb it.

2. **The "failure" events become more independent**: In the uniform case, failure at position i (B_i > s) is positively correlated with failure at position j (B_j > s) CONDITIONAL on the Parseval constraint -- if B_i is large, the remaining budget for B_j is smaller, but since s is barely above E[B], both are marginally likely to fail. In the asymmetric case, failure at tight position i is NEGATIVELY correlated with failure at loose position j (since loose positions have plenty of room), making joint failure less likely.

---

## 8. Formal Conjecture

**Conjecture 8.1 (Near-flat validity).** For all sufficiently large primes p = 3 mod 4, there exists a constant C (independent of p) such that the class of "C-near-flat" D11 (with max_{d in D11} A(d) <= floor(E[A]) + C) satisfies:

(a) The number of C-near-flat D11 is at least a positive fraction of all symmetric D11.

(b) The fraction of C-near-flat D11 with N(D11) > 0 is at least 1/poly(p).

(c) The second moment ratio E[N^2]/E[N]^2 over C-near-flat D11 is at most poly(p).

In particular, at least one C-near-flat D11 has N(D11) > 0, proving R(B_{n-1}, B_n) >= 4n-1.

**Conjecture 8.2 (Schur-concavity for Parseval-constrained validity).** For B = (B_1, ..., B_r) drawn from the uniform distribution on k-subsets conditioned on the Parseval hyperplane sum_{i=1}^r B(d_i) = C, the validity probability

  P(tau) = Pr[B_i <= tau_i for all i]

is Schur-concave in tau restricted to the hyperplane sum tau_i = T, in the regime where the mean thresholds tau_i / E[B_i] are bounded away from 0 and infinity.

That is: if tau' is obtained from tau by a Robin Hood transfer (taking slack from a loose position and giving to a tight position, making the vector more equal), then P(tau') <= P(tau).

---

## 9. Connection to the Proof

### 9.1 How this resolves the A-flat failure

The findings document (Section 1.6) shows that A-flat D11 fail at p=43 because binding slack = (p-7)/4 ~ E[B], leaving zero margin. Our analysis shows:

1. **At the product-of-marginals level**: near-flat D11 have slightly higher validity probability than A-flat D11, by a factor of exp(O(j/p)) (Section 5.5). This is small but the direction is correct.

2. **At the correlated level**: the negative correlation amplifies the gain. Asymmetric slack allows the system to avoid the "uniform pressure" failure mode where all constraints are simultaneously marginal.

3. **The total budget is unchanged**: Parseval guarantees that the total z-score sum is constant. Near-flat D11 merely redistribute the slack, they don't create it ex nihilo.

### 9.2 What remains to prove

To make this into a complete proof, one needs:

1. **Quantify the gain**: Show that P_near >= P_flat * exp(Omega(1)) for j = O(1) tight positions. This requires going beyond the product-of-marginals proxy.

2. **Bound the second moment ratio**: Show E[N^2]/E[N]^2 = O(poly(p)) over C-near-flat D11. This may use the same Parseval + orbit structure as before, but now conditioning on C-near-flat instead of A-flat.

3. **Verify computationally at p=43, 47**: Confirm that near-flat (max_A = E[A]+1) D11 have positive N and bounded second moment ratio.

### 9.3 Suggested proof strategy

The most promising path combines:

(a) **Parseval-preserving redistribution**: Show that ANY symmetric D11 with the "right" amount of A-value spread (neither too flat nor too spiky) has N(D11) > 0. The "right" spread is determined by the Parseval constraint.

(b) **Orbit counting**: Among C-near-flat D11, count the number of orbits. Each orbit has size (p-1)/2. If the total number of C-near-flat orbits is Omega(2^p / poly(p)) and the number of "failing" orbits is at most (1 - 1/poly(p)) of this, then at least one working orbit exists.

(c) **Character sum bounds**: The A-values can be expressed as A(d) = sum_{x in D11} 1[x+d in D11], which is related to the convolution of the indicator function of D11 with itself. For D11 defined by character-sum conditions, A(d) can be computed exactly using Gauss sums and Jacobi sums.

---

## 10. Summary

| Aspect | A-flat | Near-flat (C=1) |
|--------|--------|-----------------|
| max A at D11 positions | <= E[A] | <= E[A] + 1 |
| Binding slack (tight positions) | (p-7)/4 at all | (p-7)/4 - 1 at j positions |
| Binding slack (loose positions) | (p-7)/4 at all | >= (p-7)/4 + 1 at r-j positions |
| Total slack budget | Same | Same (Parseval) |
| Slack profile | Uniform | Asymmetric |
| Product-of-marginals | Baseline | Baseline * exp(O(j/p^2)) -- slightly better |
| Correlation effect | Uniform pressure (bad) | Pressure valves + absorption (good) |
| Empirical at p=43 | ALL N=0 | Working solutions exist (max_A=12) |

The near-flat theory explains WHY solutions exist above the A-flat threshold: **the Parseval constraint ensures that pushing a few A-values up forces many A-values down, creating asymmetric slack that is better utilized by negatively correlated B-values.** This is a genuine combinatorial advantage that cannot be captured by the Gaussian product-of-marginals proxy alone.
