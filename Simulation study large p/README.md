## Simulation Study (Frobenius Norm Scenario)

This folder contains scripts to reproduce the simulation results under the Frobenius norm scenario (`p > n`; here `p` is equivalent to `d` as defined in the manuscript). The simulations are organized by method and followed by a summary step.

### Method: BLOC
Run `simulation_BLOC_p_100_n_50.m`, `simulation_BLOC_p_100_n_100.m` and `simulation_BLOC_p_200_n_100.m` for problem sizes `(p,n) = (100,50), (100,100), (200,100)`, with covariance structures specified by `method_num = 1` (Block matrix), `method_num = 2` (Toeplitz), and `method_num = 3` (Banded), and penalty choices given by `penalty = 1` (SCAD) and `penalty = 2` (MCP).

### Method: Wei-Zhao
Run `Summary_WeiZhaoMCP_COV.m` and `Summary_WeiZhaoMCP.m` for problem sizes `(p,n) = (100,50), (100,100), (200,100)`, with covariance structures specified by `method_num = 1` (Block matrix), `method_num = 2` (Toeplitz), and `method_num = 3` (Banded), and penalty choices given by `penalty = 1` (SCAD) and `penalty = 2` (MCP).

### Other Methods
Run `simulation_other_methods.m` for `p = 100` with `Ns = [50, 100]`, and for `p = 200` with `Ns = [100]`, using covariance structures defined by `method_num = 1` (Block matrix), `method_num = 2` (Toeplitz), and `method_num = 3` (Banded).

### Results Aggregation
Finally, run `Summary_ALL_methods_v2.m` for problem sizes `(p,n) = (100,50), (100,100), (200,100)` and covariance structures `method_num = 1` (Block matrix), `method_num = 2` (Toeplitz), and `method_num = 3` (Banded). This script generates all reported mean and standard error summaries corresponding to this simulation scenario, as presented in the manuscript.


