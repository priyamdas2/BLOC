# BLOC: Black-box Optimization over Correlation Matrices

**BLOC** (*Black-box Optimization over Correlation matrices*) is a general-purpose framework for **sparse correlation and covariance matrix estimation** under **arbitrary (possibly non-convex, non-differentiable, or black-box) objective functions**.

BLOC is designed to operate **directly on the space of valid correlation matrices**, guaranteeing **positive definiteness and unit diagonals at every iteration**, while enabling **global, gradient-free optimization**. The method couples a **geometric reparameterization of the correlation matrix manifold** with a powerful **recursive pattern search algorithm**, making it both flexible and scalable in high dimensions.

---

## 🔑 Key Features

- **Guaranteed validity**  
  Every iterate produced by BLOC is a valid correlation matrix (symmetric, positive definite, unit diagonal).

- **Black-box optimization**  
  The objective function only needs to be *evaluated*—no gradients, Hessians, or likelihood structure required.

- **Penalty-agnostic**  
  Supports convex and non-convex penalties (LASSO, Ridge, Elastic Net, SCAD, MCP, capped-ℓ₁, or user-defined penalties).

- **Global exploration**  
  Built-in restart and step-size reset mechanisms enable systematic escape from poor local minima.

- **Parallelizable**  
  Coordinate evaluations can be run in parallel; a \( d \times d \) correlation matrix admits up to \( d(d-1)/2 \) simultaneous coordinate polls.

- **Scalable to high dimensions**  
  Suitable for \( p < n \) and \( p > n \) regimes.

---

## 📌 Problem Setting

Let Γ₀ denote a sparse *d × d* correlation matrix.  
BLOC targets optimization problems of the form:

<pre>
minimize   hₙ(Γ) + ∑<sub>i ≠ j</sub> p<sub>λ</sub>(|γ<sub>ij</sub>|)
subject to Γ ∈ 𝒞<sub>d</sub>
</pre>

where:

- **𝒞<sub>d</sub>** denotes the space of full-rank correlation matrices (symmetric, positive definite, unit diagonal).
- **h<sub>n</sub>(·)** is an arbitrary empirical loss (e.g., Gaussian negative log-likelihood, Frobenius-norm loss, robust loss, depth-based loss).
- **p<sub>λ</sub>(·)** is a sparsity-inducing penalty applied to off-diagonal entries (e.g., LASSO, SCAD, MCP, capped-ℓ₁, or user-defined penalties).

Direct optimization over **𝒞<sub>d</sub>** is challenging due to the positive-definiteness and unit-diagonal constraints.  
BLOC resolves this by reparameterizing the correlation-matrix manifold via an angular Cholesky transformation, converting the problem into an unconstrained Euclidean optimization while ensuring every iterate remains a valid correlation matrix.

---

## 🧠 Method Overview

### 1. Angular Cholesky Reparameterization

Every correlation matrix \( \Gamma \in \mathcal{C}_d \) can be uniquely written as
\[
\Gamma = L L^\top,
\]
where each row of \( L \) lies on a unit hypersphere. BLOC represents these rows using **hyperspherical angles**, yielding a **smooth bijection** between correlation matrices and an open hyperrectangle in \( \mathbb{R}^{d(d-1)/2} \).

<p align="center">
  <img src="images/BLOC_diagram.jpg" width="85%">
</p>

This mapping ensures:
- unconstrained optimization in Euclidean space,
- automatic enforcement of correlation-matrix constraints.

---

### 2. Recursive Modified Pattern Search (RMPS)

On the transformed space, BLOC applies **Recursive Modified Pattern Search (RMPS)**, a derivative-free global optimization algorithm featuring:

- coordinate-wise polling in \( \pm \) directions,
- adaptive step-size reduction,
- run-wise restarts for global exploration,
- parallel evaluation of candidate points.

<p align="center">
  <img src="images/BLOC_concept_v2.png" width="90%">
</p>

---

## 📐 Theoretical Guarantees

BLOC is supported by rigorous theory, including:

- **Statistical guarantees**  
  - Frobenius-norm convergence rates for local minimizers  
  - Sparsistency (exact support recovery) under standard non-convex penalties  

- **Algorithmic guarantees**  
  - Stationarity when no improving coordinate direction exists  
  - Probabilistic reachability of neighborhoods of the global minimizer  
  - Global convergence in probability under mild regularity conditions  
  - Sublinear convergence rates for smooth convex objectives (within a run)

These results hold for **general losses**, extending classical theory beyond Gaussian likelihoods.

---

## 🧪 Empirical Performance

Extensive simulations and benchmark studies show that BLOC:

- outperforms state-of-the-art sparse covariance estimators under non-convex penalties,
- achieves lower Frobenius and spectral norm errors,
- delivers superior sparsity recovery,
- remains stable in high dimensions.

A real-data application to **pan-gynecologic proteomics networks** demonstrates BLOC’s ability to produce interpretable, pathway-informed correlation structures.

---

## 🗂 Repository Structure

```text
BLOC/
├── BLOC/                  # Core MATLAB source code
├── DEMO/                  # Reproducible demo scripts (serial & parallel)
├── Supp functions/        # Supporting utilities
├── images/                # Figures used in paper and README
├── manopt/                # Placeholder (see note below)
├── README.md              # This file
```
