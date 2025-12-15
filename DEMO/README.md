# DEMO — BLOC Usage Examples (Non-parallel and Parallel)

This folder contains a **fully self-contained demonstration** showing how to use the **BLOC** optimization framework for constrained optimization over the space of **valid correlation matrices**.

The demo illustrates both:
- **BLOC** (non-parallel solver), and  
- **BLOCparallel** (parallel solver),

through simple, reproducible benchmark examples.

The main demo script is:

**`BLOC_demo.m`**

All files required to run this demo are located **within this `DEMO/` folder**.

---

## Overview of the demo

The script `BLOC_demo.m` demonstrates two representative optimization tasks.

### Example 1: Non-parallel BLOC

- **Matrix dimension:** 10 × 10  
- **Decision variable:** Correlation matrix  
- **Objective function:** `modified_ackley(C)`  
- **Solver:** `BLOC` (non-parallel)  

This example illustrates the **basic usage** of BLOC with minimal solver options.

### Example 2: Parallel BLOC

- **Matrix dimension:** 20 × 20  
- **Decision variable:** Correlation matrix  
- **Objective function:** `modified_ackley(C)`  
- **Solver:** `BLOCparallel`  

This example demonstrates the **parallel version** of BLOC on a higher-dimensional problem and highlights additional solver controls.

In both examples, the workflow is:
1. Generate a **random feasible correlation matrix** as the starting point
2. Run the corresponding BLOC solver
3. Report the optimized solution, objective value, and computation time

---

## Folder structure

The `DEMO/` folder is designed to be **self-sufficient**. A typical structure is:

- `BLOC_demo.m`  
  Main demo script illustrating non-parallel and parallel usage.

- `BLOC/`  
  Folder containing the core BLOC solver implementations and internal utilities.

- `modified_ackley.m`  
  Benchmark objective function used in the demo.

- `randCorrMatrix.m`  
  Utility function for generating valid random correlation matrices.

- Additional helper files required internally by BLOC.

No external dependencies are required beyond MATLAB itself (and the Parallel Computing Toolbox for the parallel example).

---

## Requirements

- **MATLAB:** R2018b or newer is recommended.
- **Parallel Computing Toolbox:** Required **only** for running the `BLOCparallel` example.  
  The non-parallel example runs without any additional toolboxes.

---

## How to run the demo

### Recommended method

1. Open MATLAB.
2. Set the current folder to `DEMO/`.
3. Run:
   `BLOC_demo`

### Alternative method

You may also run the demo from any location by adding the folder to the MATLAB path:

`addpath('path/to/your/repository/DEMO');`  
`BLOC_demo`

---

## Expected output

Running `BLOC_demo.m` produces two clearly labeled result sections in the MATLAB Command Window:

1. **BLOC (Non-parallel) Optimization Results**
2. **BLOCparallel Optimization Results**

Each section reports:
- Problem dimension (size of the correlation matrix)
- Final minimized objective value
- Total wall-clock computation time (in seconds)
- The optimized correlation matrix
- A consistency check by re-evaluating the objective at the final solution

This output format is intended to be **self-explanatory** and suitable for both users and reviewers.

---

## Starting point and feasibility

The demo uses the function:

`C0 = randCorrMatrix(M, rand_seed);`

This guarantees that the initial point:
- is symmetric
- has unit diagonal
- is positive semi-definite

Starting from a valid correlation matrix is essential, since BLOC performs optimization **directly on the constrained correlation matrix space**.

---

## Solver options used in the demo

The demo uses modest option settings for clarity and speed.

### BLOC (non-parallel)

`options.MaxRuns = 5;`  
`options.MaxTime = 20;`

- `MaxRuns`: maximum number of outer runs (restarts)
- `MaxTime`: maximum allowed wall-clock time (seconds)

### BLOCparallel

`options.MaxRuns      = 10;`  
`options.phi          = 1e-18;`  
`options.MaxIter      = 5000;`  
`options.DisplayEvery = 5;`

- `MaxRuns`: maximum number of outer runs
- `phi`: minimum allowed step-size
- `MaxIter`: maximum number of iterations per run
- `DisplayEvery`: every that many seconds update is displayed (only if `DisplayUpdate` is set to 1; that's true by default)

These values can be modified directly in `BLOC_demo.m` to trade off speed versus accuracy.

---

## Using BLOC with your own objective function

To use BLOC with a custom objective, define a function handle of the form:

`objFun = @(C) your_objective(C);`

where `C` is an `M × M` correlation matrix and the function returns a scalar.

Example:

`M = 15;`  
`objFun = @(C) norm(C - eye(M), 'fro');`  
`C0 = randCorrMatrix(M, 123);`

`options.MaxRuns = 5;`  
`options.MaxTime = 30;`

`[C_hat, fval, comp_time] = BLOC(objFun, C0, options);`

---

## Troubleshooting

- **Undefined function errors:**  
  Ensure MATLAB is running inside the `DEMO/` folder or that it has been added to the path.

- **Parallel solver not working:**  
  Verify that the Parallel Computing Toolbox is installed. The non-parallel example does not require it.

- **Path issues:**  
  Check the following line in `BLOC_demo.m` and update it if needed:  
  `addpath('./BLOC/');`

---

## Citation

If you use BLOC in academic work, please cite the software or accompanying paper as indicated in the main repository.

---
## 💬 Contact

For questions, please contact:  
**Priyam Das**  
[dasp4@vcu.edu](mailto:dasp4@vcu.edu)

