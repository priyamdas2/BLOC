# BLOC: Black-box Optimization over Correlation Matrices

**BLOC** (*Black-box Optimization over Correlation matrices*) is a general-purpose framework for **sparse correlation and covariance matrix estimation** under **arbitrary (possibly non-convex, non-differentiable, or black-box) objective functions**.

BLOC is designed to operate **directly on the space of valid correlation matrices**, guaranteeing **positive definiteness and unit diagonals at every iteration**, while enabling **global, gradient-free optimization**. The method couples a **geometric reparameterization of the correlation matrix manifold** with a powerful **recursive pattern search algorithm**, making it both flexible and scalable in high dimensions.

---

### Notation note

Throughout the code, comments, and documentation, the symbols **`p`** and **`d`** may occasionally be used interchangeably.  
Both symbols refer to the **same quantity**: the **dimension of the correlation matrix of interest**.

This interchangeable usage reflects differences in notation commonly adopted across optimization and statistical literature.  
For clarity, readers may treat **`p ≡ d`** everywhere in this repository.

---

## 📂 Repository Structure

The current top-level layout of this repository (branch: `main`) is:

```text
BLOC/
├── Benchmark/                 # Benchmark experiments and comparisons
├── DEMO/                      # Demo scripts illustrating how to use BLOC (serial + parallel examples)
├── Real data analysis/        # Real-data application workflows (proteomics network analysis)
├── Simulation study/          # Simulation experiments for baseline / moderate dimension settings, using Gaussian likelihood
├── Simulation study large p/  # Simulation experiments for large-d (or p) / high-dimensional regimes, using Frobenius norm
├── images/                    # Figures used in the paper/README (e.g., diagrams, flowcharts)
└── README.md                  # Main repository README
```

## 🔑 Key Features

- **Guaranteed validity**  
  Every iterate produced by BLOC is a valid correlation matrix: symmetric, positive definite, and with unit diagonal entries.

- **Black-box optimization**  
  The objective function only needs to be *evaluated*. No gradients, Hessians, likelihood structure, or smoothness assumptions are required.

- **Penalty-agnostic**  
  Supports both convex and non-convex penalties, including LASSO, Ridge, Elastic Net, SCAD, MCP, capped-ℓ₁, as well as fully user-defined penalties.

- **Global exploration**  
  Built-in restart and step-size reset mechanisms allow the algorithm to systematically escape poor local minima and explore the objective landscape.

- **Parallelizable**  
  Coordinate-wise objective evaluations can be executed in parallel.  
  For a d-by-d correlation matrix, up to d(d−1)/2 coordinate directions can be evaluated simultaneously.

- **Scalable to high dimensions**  
  Designed to perform reliably in both low-dimensional (d < n) and high-dimensional (d > n) regimes.

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

Every correlation matrix Γ belonging to the space 𝒞_d admits a unique Cholesky decomposition of the form

Γ = L · Lᵀ,

where L is a lower-triangular matrix with unit-norm rows and positive diagonal entries.  
Each row of L therefore lies on a unit hypersphere.

BLOC represents these rows using **hyperspherical (angular) coordinates**, which induces a **smooth, bijective mapping** between the space of valid correlation matrices and an open hyperrectangle in Euclidean space of dimension d(d−1)/2. This transformation allows unconstrained optimization while guaranteeing that every iterate corresponds to a valid correlation matrix.


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
  
---


## 🧬 Pathway-informed correlation estimation for pan-gynecologic proteomics data

We illustrate the practical utility of **BLOC** through an application to **reverse-phase protein array (RPPA)** data from *The Cancer Genome Atlas (TCGA)*, focusing on five pan-gynecologic cancers:

- Breast carcinoma (**BRCA**)
- Cervical squamous cell carcinoma (**CESC**)
- Ovarian serous cystadenocarcinoma (**OV**)
- Uterine corpus endometrial carcinoma (**UCEC**)
- Uterine carcinosarcoma (**UCS**)

These cancers share common hormonal drivers and signaling dysregulation, yet differ markedly in clinical presentation and disease progression. Understanding **protein–protein correlation structure across key signaling pathways** is therefore of strong biological and translational interest.

The proteomics data are obtained from *The Cancer Proteome Atlas (TCPA)* and consist of RPPA-based protein abundance measurements spanning diverse cellular processes, including cell cycle regulation, hormone receptor activity, PI3K/AKT signaling, apoptosis, metabolism, immune response, and epithelial–mesenchymal transition.

---

### Pathway-aware modeling via penalty covers

Following prior biological literature, we focus on a curated panel of **27 proteins**, grouped into **five non-overlapping functional pathways**:

- **Cell Cycle**
- **Hormone Receptor**
- **Hormone Signaling**
- **PI3K/AKT**
- **Breast Reactive**

A key strength of BLOC is that it can optimize **any user-defined objective function over the space of correlation matrices**. We leverage this flexibility to incorporate biological prior knowledge through a **structured penalty cover**.

Specifically, correlations **within the same pathway** are *not penalized*, while correlations **across different pathways** are penalized using a non-convex SCAD penalty. This design reflects the expectation that proteins within a pathway are biologically coordinated, whereas cross-pathway correlations should only persist if strongly supported by data.

As a result, the induced sparsity pattern is **block-structured**, preserving within-pathway coherence while adaptively shrinking weaker cross-pathway associations.

---

### Estimated correlation structures across cancers

Below we display the estimated sparse correlation heatmaps obtained using **BLOC with a SCAD penalty and pathway-based penalty cover** for each cancer type.

<p align="center">
  <img src="images/plot_corr_BRCA_cropped.jpg" width="45%">
  <img src="images/plot_corr_CESC_cropped.jpg" width="45%">
</p>

<p align="center">
  <img src="images/plot_corr_OV_cropped.jpg" width="45%">
  <img src="images/plot_corr_UCEC_cropped.jpg" width="45%">
</p>

<p align="center">
  <img src="images/plot_corr_UCS_cropped.jpg" width="45%">
</p>

**Figure:** Estimated sparse correlation heatmaps for five pan-gynecologic cancers using BLOC with SCAD and a pathway-informed penalty cover. Within-pathway structure is preserved, while cross-pathway associations reveal tumor-specific differences in signaling integration.

---

### Biological interpretation

Several biologically meaningful patterns emerge:

- **Within-pathway coherence** is consistently preserved, with strong positive correlations observed among Cell Cycle proteins (e.g., FOXM1, PCNA, CYCLINB1) and within the PI3K/AKT module (e.g., AKT phosphorylation sites, PRAS40, GSK3 isoforms, PTEN).

- **BRCA and UCEC** exhibit pronounced cross-talk between Hormone Receptor and Hormone Signaling pathways, consistent with estrogen-driven biology.

- **OV** shows partial integration between PI3K/AKT and Cell Cycle modules, suggesting coordinated proliferative signaling.

- **CESC** displays comparatively weaker and more diffuse cross-pathway correlations, reflecting a less hormone-dependent disease mechanism.

- **UCS** presents a strikingly modular structure, with near-independent Cell Cycle and PI3K/AKT blocks. While this may reflect genuine biological heterogeneity, it should be interpreted cautiously due to limited sample size.

Sample sizes differ substantially across cancers (BRCA: 879, CESC: 171, OV: 428, UCEC: 404, UCS: 48), and sparsity levels—particularly for UCS—are influenced in part by statistical power.

---

### Summary

This real-data application highlights the strength of BLOC in **embedding biological prior knowledge directly into correlation estimation**. By preserving expected within-pathway associations and selectively shrinking cross-pathway edges, BLOC yields interpretable, tumor-specific network structures that align closely with known signaling biology and may inform downstream therapeutic insights.



