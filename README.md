# BLOC: Black-box Optimization over Correlation Matrices

**BLOC** is a general-purpose, derivative-free optimization framework for estimating **sparse covariance and correlation matrices** under **arbitrary user-defined loss and penalty functions**. The method operates directly on the nonlinear manifold of correlation matrices—symmetric, positive definite matrices with unit diagonal—while supporting **non-convex**, **non-differentiable**, and **black-box** objectives.

This repository accompanies the paper:

> **BLOC: A Flexible Framework for Sparse Covariance Estimation with Existing and User-Defined Penalties**

and provides fully reproducible code for benchmark studies, simulation experiments, and real-data analyses reported in the manuscript.

---

## Motivation

Many existing covariance and correlation estimation methods are penalty-specific, rely on convex relaxations or local approximations, and lack mechanisms to escape local optima in non-convex settings. **BLOC** addresses these limitations by combining (i) a smooth angular reparameterization of correlation matrices via Cholesky decomposition, which converts constrained matrix optimization into unconstrained Euclidean optimization, and (ii) a global, derivative-free optimization engine based on **Recursive Modified Pattern Search (RMPS)** with adaptive step-size control and restart mechanisms. As a result, BLOC enables sparse covariance estimation with **any penalty of choice** (e.g., Lasso, Elastic Net, SCAD, MCP, Ridge, or user-defined penalties) while guaranteeing that all iterates remain valid correlation matrices.

---

## Key Contributions

- Penalty-agnostic framework: plug in any sparsity-inducing penalty or loss function  
- Global optimization: explicit restart strategy to escape local minima in non-convex landscapes  
- Geometry-aware design: optimization over a bijective angular representation of correlation matrices  
- Theoretical guarantees: stationarity, reachability, global convergence in probability, and sublinear convergence rates under mild conditions  
- Scalable implementation: supports parallel evaluations and high-dimensional settings  

---

## Method Overview

BLOC solves optimization problems of the form  
\[
\min_{\mathbf{C} \in \mathcal{C}_M} h(\mathbf{C}), \quad 
\mathcal{C}_M = \{\mathbf{C} \succ 0,\; \mathbf{C} = \mathbf{C}^\top,\; \mathrm{diag}(\mathbf{C}) = 1\},
\]
by mapping \(\mathbf{C}\) to angular parameters via a bijective Cholesky–hyperspherical transformation, optimizing the transformed objective using RMPS in unconstrained Euclidean space, and mapping the solution back to a valid correlation matrix.

A schematic illustration of the search principle (Fermi’s principle) underlying RMPS is shown below.

![Fermi principle illustration](image/Fermi.png)

---

## Repository Structure

```text
BLOC/
├── Benchmark/                 # Global optimization benchmark experiments
├── Simulation/                # Simulation study (covariance estimation)
├── Simulation_Frobenius/      # p > n Frobenius norm simulations
├── Real_data_analysis/        # Pan-gynecologic proteomics application
├── image/                     # Figures used in the manuscript
├── utils/                     # Helper functions
├── README.md                  # This file
```
