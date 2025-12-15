## Manopt Toolbox Requirement

For copyright reasons, the **Manopt** source files are not redistributed in this repository.  
Please download the Manopt MATLAB toolbox directly from the official website:

👉 https://www.manopt.org/

After downloading and extracting the toolbox, the `manopt` directory should have the following structure:

```text
manopt/
├── autodiff
├── core
├── lifts
├── manifolds
├── solvers
└── tools
```
After downloading, please replace the existing `manopt` folder in this repository (which contains only instructions and no source files) with the downloaded `manopt` folder.

## Benchmark study code reproducibility

To reproduce all numerical results reported in the main manuscript, please follow the steps below.

### Step 1: Low- and Moderate-Dimensional Benchmark Experiments
Run the script `Benchmark_study.m` for the following test functions:

- `input_vals = 1` : Ackley  
- `input_vals = 2` : Griewank  
- `input_vals = 3` : Rosenbrock  
- `input_vals = 4` : Rastrigin  

For each function, set the problem dimension to:


### Step 2: High-Dimensional Benchmark Experiments
Run the script `Benchmark_study_high_dim.m` using the same set of test functions:
- Ackley, Griewank, Rosenbrock, and Rastrigin (`input_vals = 1–4`)

with the high-dimensional setting:


### Step 3: Results Aggregation
Finally, execute the script `Summary.m` to compile and summarize all benchmark results.

All final outputs, as reported in the main draft, are automatically generated and saved in the directory:


These files contain the complete set of tables and summaries used in the manuscript.


