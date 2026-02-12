# Literature Review: Schur-Convexity of Orthant Probabilities Under Negative Dependence

**Date**: 2026-02-11
**Purpose**: Support L6 proof (second moment method) for R(B_{n-1}, B_n) = 4n-1

---

## 1. The Conjecture

On the hyperplane Sigma B_i = C (deterministic Parseval constraint), for B drawn from the uniform k-subset measure on Z_p\{0}:

> Pr[B_1 <= T_1, ..., B_R <= T_R] is **Schur-convex** in (T_1,...,T_R).

That is, if T' majorizes T (T' has more spread, same sum), then the joint probability under T' is >= that under T. This would imply that the "correlation bonus" c_0 = Pr[all B_i <= T_i] / Prod Pr[B_i <= T_i] is maximized when thresholds are most spread out -- precisely the structure we see in the D11 constraints.

---

## 2. Directly Applicable Results

### 2.1. Tong (1982): Rectangular Probabilities for Schur-Concave Densities

**Reference**: Y.L. Tong, "Rectangular and Elliptical Probability Inequalities for Schur-Concave Random Variables," *Annals of Statistics* 10(2), 1982.

**Result**: If the density f(x) of X = (X_1,...,X_n) is Schur-concave, then P(|X_i| <= a_i, i=1,...,n) is **Schur-concave** in (a_1,...,a_n).

**Applicability to our problem**: This goes in the **wrong direction** for us. Tong's result says that for Schur-concave densities (including exchangeable normals), rectangular probabilities are Schur-CONCAVE in the thresholds -- meaning UNIFORM thresholds MAXIMIZE the probability. We need Schur-CONVEXITY.

**Key gap**: Tong's result applies to the density itself being Schur-concave. Our distribution on the hyperplane Sigma B_i = C is NOT a standard exchangeable distribution with Schur-concave density. The k-subset measure is fundamentally discrete and has strong algebraic structure.

### 2.2. Joag-Dev and Proschan (1983): Negative Association

**Reference**: K. Joag-Dev and F. Proschan, "Negative Association of Random Variables with Applications," *Annals of Statistics* 11(1), 286-295, 1983.

**Key results**:
1. **Definition**: Random variables X_1,...,X_k are negatively associated (NA) if Cov(f(X_i, i in A), g(X_j, j in B)) <= 0 for all disjoint subsets A, B and all nondecreasing functions f, g.
2. **Orthant inequality for NA**: P(X_i <= x_i for all i) <= Prod P(X_i <= x_i). Dually, P(X_i >= x_i for all i) <= Prod P(X_i >= x_i).
3. **Closure**: Nondecreasing functions of mutually exclusive subsets of NA random variables are NA.
4. **Negatively correlated normals are NA**.

**Applicability**: Confirms that for our k-subset B-values, P(all B_i <= T_i) <= Prod P(B_i <= T_i). But this gives the OPPOSITE of what we want for the "correlation bonus" -- NA says the joint is BELOW the product. However, our problem involves B-values on a FIXED-SUM HYPERPLANE, and the B-values conditional on the hyperplane may NOT be NA (they are constrained, which can create positive-like effects for certain threshold configurations).

### 2.3. Proschan and Sethuraman (1977): Schur Functions in Statistics

**Reference**: F. Proschan and J. Sethuraman, "Schur Functions in Statistics I. The Preservation Theorem," *Annals of Statistics* 5(2), 1977.

**Result**: If f(x) is Schur-concave and phi(lambda, x) is TP2 satisfying the semigroup property, then the integral h(lambda_1,...,lambda_n) = integral f(x) Prod phi(lambda_i, x_i) dx is Schur-concave.

**Applicability**: This is a powerful preservation result but applies to CONTINUOUS integral transforms. The k-subset measure does not have this product structure.

---

## 3. Near-Miss Results

### 3.1. Slepian's Inequality and Its Majorization Extensions

**Classical Slepian (1962)**: If X, Y are centered Gaussian vectors with Var(X_i) = Var(Y_i) and Cov(X_i,X_j) <= Cov(Y_i,Y_j) for i != j, then P(X_i <= u_i for all i) >= P(Y_i <= u_i for all i).

**Fang and Zhang (2010)**: "Slepian's inequality with respect to majorization," *Linear Algebra and Applications* 434, 1107-1118. Extends Slepian to majorization: if the covariance matrix of Y majorizes that of X (in some sense), the orthant probabilities can be compared.

**Fernandez et al. (2013)**: "Slepian's inequality for Gaussian processes with respect to weak majorization," *Journal of Inequalities and Applications*, 2013:5. Further extends to weak majorization of covariance structures.

**Why these almost but don't quite work**: These results compare orthant probabilities as the COVARIANCE STRUCTURE changes, not as the THRESHOLDS change. Our problem fixes the covariance structure (equicorrelated with rho = -2/(p-3)) and varies the thresholds. The classical Slepian inequality also goes in the wrong direction for negative correlations -- it says that MORE positive correlation gives LARGER joint probability. For negative correlations, we are outside the regime of these results.

### 3.2. Gaussian Correlation Inequality (Royen 2014)

**Reference**: T. Royen, proved the full Gaussian Correlation Inequality.

**Result**: For any centered Gaussian measure and symmetric convex sets A, B: mu(A intersect B) >= mu(A) * mu(B).

**Why it doesn't apply**: This concerns intersection of TWO symmetric convex sets, not the relationship between orthant probabilities and threshold spread. Also, our sets (half-spaces B_i <= T_i) are not symmetric.

### 3.3. Borcea, Branden, Liggett (2009): Strong Rayleigh Property

**Reference**: J. Borcea, P. Branden, T.M. Liggett, "Negative Dependence and the Geometry of Polynomials," *JAMS* 22(2), 521-567, 2009.

**Key results**:
1. Introduces **strongly Rayleigh** (SR) measures via stability of generating polynomials.
2. SR implies conditional negative association (CNA+), which implies NA.
3. The uniform measure on k-element subsets of [n] IS strongly Rayleigh.
4. SR measures enjoy concentration (Lipschitz functions concentrate).

**Why it almost works but fails**: The k-subset measure IS strongly Rayleigh, so our B-value indicators are SR. SR implies NA, which gives P(all B_i <= T_i) <= Prod P(B_i <= T_i). But this gives an UPPER bound on the joint probability, not the Schur-convexity of the joint probability as a function of thresholds.

**Potential connection**: The CNA+ property (conditional negative association closed under external fields) means that even after conditioning on some coordinates, the remaining ones are NA. This is relevant because our Parseval constraint Sigma B_i = C is a conditioning event. If we could show that CONDITIONAL on Sigma B_i = C, the joint probability P(all B_i <= T_i | Sigma B_i = C) has Schur-convexity in T, this would suffice. But CNA+ only gives conditional NA, which bounds the ratio below 1, not Schur-convexity of the ratio.

### 3.4. Recent: Concentration for Sampling Without Replacement via Majorization (2025)

**Reference**: arXiv:2503.20473, "Concentration inequalities for the sum in sampling without replacement: an approach via majorization," 2025.

**Key results**:
1. Uses Schur-convexity as a proof technique for concentration of hypergeometric sums.
2. Proves that the mean absolute deviation E[|X_P|] is Schur-convex in the population P.
3. Uses majorization to reduce arbitrary populations to structured benchmark cases.

**Why it's relevant but doesn't directly apply**: This paper uses Schur-convexity of a DIFFERENT quantity (mean absolute deviation) with respect to the POPULATION, not with respect to THRESHOLDS. However, the proof technique -- using majorization to reduce to structured cases -- is potentially transferable to our setting.

---

## 4. State of Knowledge: Negative Association + Majorization

### 4.1. The Central Tension

The literature reveals a fundamental tension:

**For independent random variables**: P(all X_i <= T_i) = Prod P(X_i <= T_i) = Prod Phi(T_i). Since log Phi is concave, this product is Schur-CONCAVE by the Schur-Ostrowski criterion. Uniform thresholds maximize the product.

**For negatively associated variables**: P(all X_i <= T_i) <= Prod P(X_i <= T_i), so the joint probability is BELOW the product. The joint probability is not simply a product, and its Schur-convexity/concavity depends on the dependence structure.

**For variables on a fixed-sum hyperplane**: Conditioning on Sigma X_i = C introduces strong negative correlations. The conditional joint probability can be ABOVE or BELOW the product of conditional marginals, depending on the threshold configuration. Our conjecture says this conditional probability is Schur-convex in thresholds.

### 4.2. What Is Known

1. **NA implies NOD**: P(all X_i <= T_i) <= Prod P(X_i <= T_i). This is a UNIFORM bound, independent of T.

2. **Exchangeable distributions with Schur-concave density**: Rectangular probabilities are Schur-concave (Tong 1982). This includes exchangeable normals.

3. **Equicorrelated Gaussian orthant probabilities**: Extensively studied (Steck 1962, many others). For equicorrelated normals with rho < 0, P(all X_i <= t_i) can be computed via differential equations. The behavior as thresholds spread is NOT well-characterized in the literature.

4. **k-subset measures are strongly Rayleigh**: Borcea-Branden-Liggett 2009. This gives CNA+ but not Schur-convexity of orthant probabilities.

5. **Concentration for sums**: Majorization has been used for concentration of sampling-without-replacement sums (2025), but not for orthant probabilities of these sums.

### 4.3. What Is NOT Known (Open Problems)

1. **Schur-convexity of orthant probabilities for NA/SR variables**: NO general result exists showing that P(all X_i <= T_i) is Schur-convex in T for negatively associated or strongly Rayleigh variables. This appears to be OPEN.

2. **Conditional orthant probabilities on hyperplane**: NO result in the literature addresses the Schur-convexity of P(all X_i <= T_i | Sigma X_i = C) as a function of T for k-subset measures.

3. **The "correlation bonus" ratio**: The ratio c_0 = P(all X_i <= T_i) / Prod P(X_i <= T_i) as a function of T, for negatively dependent variables, has not been systematically studied for its majorization properties.

4. **Slepian for thresholds (not covariances)**: The majorization extensions of Slepian (Fang-Zhang 2010, 2013) concern majorization of COVARIANCE MATRICES, not of THRESHOLDS. A threshold-majorization version of Slepian for equicorrelated variables appears to be unstudied.

---

## 5. Suggested Proof Techniques

Based on the literature survey, several approaches emerge:

### 5.1. Direct Combinatorial / Algebraic Approach (Most Promising)

Since the k-subset measure has algebraic structure (especially in Z_p), one could try to:
- Express P(all B_i <= T_i) as a sum over k-subsets satisfying the constraints
- Use the character sum / Fourier structure of Z_p to factor or bound this sum
- Show directly that transferring threshold from a loose constraint to a tight one (one step of majorization via T-transforms) increases the count

This bypasses the need for general probabilistic results about Schur-convexity under negative dependence.

### 5.2. Coupling / T-Transform Argument

By the T-transform characterization of majorization, it suffices to show: if T' differs from T by transferring epsilon from coordinate i to coordinate j (with T_i > T_j), then P(all B <= T') >= P(all B <= T). This reduces to showing that "loosening a tight constraint while tightening a loose one" increases probability. Under the fixed-sum Parseval constraint, this is plausible because:
- The tight constraint is likely binding (many B-values close to T_j)
- The loose constraint has margin (few B-values near T_i)
- Tightening i by epsilon removes few subsets; loosening j by epsilon adds many

### 5.3. Conditional Gaussianization

On the hyperplane Sigma B_i = C, the B-values have an approximate multivariate normal distribution with equicorrelation rho < 0. For this Gaussian proxy, one could try to directly compute d/depsilon P(all X_i <= T_i + epsilon e_ij) where e_ij transfers from i to j, and show it is positive when T_i > T_j.

For equicorrelated Gaussian with rho < 0 on the hyperplane, this reduces to a 2-dimensional integral involving the conditional bivariate density of (B_i, B_j) given the other constraints. The known differential equations for equicorrelated normal orthant probabilities (Steck 1962) might be useful here.

### 5.4. Entropy / Log-Sobolev Approach

The number of k-subsets satisfying threshold constraints can be written as the permanent of a 0-1 matrix or as a partition function. Showing Schur-convexity of this partition function might be approachable via log-Sobolev inequalities for strongly log-concave distributions (Anari et al. 2021).

### 5.5. Avoid the Conjecture Entirely

Rather than proving Schur-convexity, one could try to bound c_0 directly for the SPECIFIC threshold vectors arising from D11 (which have specific algebraic structure). The thresholds tau(d) = T(d) - A(d) have structure determined by the autocorrelation A(d) of D11. For the Paley-type D11, these have extra symmetry that might directly give c_0 = O(1).

---

## 6. Key References (Full Citations)

1. **Borcea, J., Branden, P., Liggett, T.M.** (2009). Negative dependence and the geometry of polynomials. *JAMS* 22(2), 521-567. [arXiv:0707.2340](https://arxiv.org/abs/0707.2340)

2. **Dubhashi, D., Ranjan, D.** (1998). Balls and bins: A study in negative dependence. *Random Structures & Algorithms* 13(2), 99-124. [Wiley](https://onlinelibrary.wiley.com/doi/10.1002/(SICI)1098-2418(199809)13:2<99::AID-RSA1>3.0.CO;2-M)

3. **Fang, L., Zhang, X.** (2010). Slepian's inequality with respect to majorization. *Linear Algebra and Applications* 434, 1107-1118. [ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0024379510005446)

4. **Joag-Dev, K., Proschan, F.** (1983). Negative association of random variables with applications. *Annals of Statistics* 11(1), 286-295. [Project Euclid](https://projecteuclid.org/journals/annals-of-statistics/volume-11/issue-1/Negative-Association-of-Random-Variables-with-Applications/10.1214/aos/1176346079.full)

5. **Marshall, A.W., Olkin, I., Arnold, B.C.** (2011). *Inequalities: Theory of Majorization and Its Applications* (2nd ed.). Springer. [Springer](https://link.springer.com/book/10.1007/978-0-387-68276-1)

6. **Pemantle, R.** (2000). Towards a theory of negative dependence. *J. Math. Physics* 41, 1371-1390. [arXiv:math/0404095](https://arxiv.org/abs/math/0404095)

7. **Proschan, F., Sethuraman, J.** (1977). Schur functions in statistics I. The preservation theorem. *Annals of Statistics* 5(2), 256-262. [Project Euclid](https://projecteuclid.org/journals/annals-of-statistics/volume-5/issue-2)

8. **Rinott, Y.** (1973). Multivariate majorization and rearrangement inequalities with some applications to probability and statistics. *Israel J. Math.* 15, 60-77. [Springer](https://link.springer.com/article/10.1007/BF02771774)

9. **Royen, T.** (2014). A simple proof of the Gaussian correlation conjecture extended to multivariate gamma distributions. *Far East J. Theor. Stat.* 48(2), 139-145. [arXiv:1408.1028](https://arxiv.org/abs/1408.1028)

10. **Shaked, M., Shanthikumar, J.G.** (2007). *Stochastic Orders*. Springer. [Springer](https://link.springer.com/book/10.1007/978-0-387-34675-5)

11. **Steck, G.P.** (1962). Orthant probabilities for the equicorrelated multivariate normal distribution. *Biometrika* 49(3-4), 433-445. [Oxford Academic](https://academic.oup.com/biomet/article/49/3-4/433/223641)

12. **Tong, Y.L.** (1980). *Probability Inequalities in Multivariate Distributions*. Academic Press.

13. **Tong, Y.L.** (1982). Rectangular and elliptical probability inequalities for Schur-concave random variables. *Annals of Statistics* 10(2), 637-642. [Project Euclid](https://projecteuclid.org/euclid.aos/1176345807)

14. **Tong, Y.L.** (1990). *The Multivariate Normal Distribution*. Springer. [Springer](https://link.springer.com/book/10.1007/978-1-4613-9655-0)

15. **arXiv:2503.20473** (2025). Concentration inequalities for the sum in sampling without replacement: an approach via majorization. [arXiv](https://arxiv.org/abs/2503.20473)

16. **Fernandez, L., Sempi, C., Valdes, J.E.** (2013). Slepian's inequality for Gaussian processes with respect to weak majorization. *J. Inequal. Appl.* 2013:5. [Springer](https://link.springer.com/article/10.1186/1029-242X-2013-5)

17. **Pinasco, D., Smucler, E.** (2020). Orthant probabilities and the attainment of maxima on a vertex of a simplex. [arXiv:2004.04682](https://arxiv.org/abs/2004.04682)

---

## 7. Summary Assessment

**The conjecture appears to be genuinely new.** No result in the existing literature proves or disproves Schur-convexity of orthant probabilities as a function of thresholds for negatively dependent random variables (whether NA, SR, or conditional on a fixed sum).

The closest results are:
- **Tong (1982)**: Proves Schur-CONCAVITY for Schur-concave densities (OPPOSITE direction)
- **Borcea-Branden-Liggett (2009)**: Establishes that k-subset measures are SR/CNA+, giving P(joint) <= Product, but says nothing about Schur-convexity in thresholds
- **Fang-Zhang (2010)**: Slepian majorization, but for COVARIANCE comparison, not THRESHOLD comparison

**The most promising approaches** are:
1. **Direct combinatorial/algebraic** (Section 5.1): Use the specific structure of B-values in Z_p
2. **T-transform coupling** (Section 5.2): Show single-step majorization increases probability
3. **Avoid the conjecture** (Section 5.5): Bound c_0 directly for the specific D11 structures that arise, rather than proving a general Schur-convexity result
