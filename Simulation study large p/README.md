## Simulation Study (Frobenius Norm Scenario)

### Notation note

Throughout the code, comments, and documentation, the symbols **`p`** and **`d`** may occasionally be used interchangeably.  
Both symbols refer to the **same quantity**: the **dimension of the correlation matrix of interest** (i.e., the number of variables).

This interchangeable usage reflects differences in notation commonly adopted across optimization and statistical literature.  
For clarity, readers may treat **`p ≡ d`** everywhere in this repository.


This folder contains scripts to reproduce the simulation results under the Frobenius norm scenario (`p > n`; here `p` is equivalent to `d` as defined in the manuscript). The simulations are organized by method and followed by a summary step.

#### Method: BLOC
Run `simulation_BLOC.m` for problem sizes `(p,n) = (50,50), (50,100), (100,100)`, with covariance structures specified by `method_num = 1` (Block matrix), `method_num = 2` (Toeplitz), and `method_num = 3` (Banded), and penalty choices given by `penalty = 1` (SCAD) and `penalty = 2` (MCP).

#### Other Methods
Run `simulation_other_methods.m` for `p = 50` with `Ns = [50]`, and for `p = 100` with `Ns = [50, 100]`, using covariance structures defined by `method_num = 1` (Block matrix), `method_num = 2` (Toeplitz), and `method_num = 3` (Banded).

#### Results Aggregation
Finally, run `Summary_all_methods.m` for problem sizes `(p,n) = (50,50), (50,100), (100,100)` and covariance structures `method_num = 1` (Block matrix), `method_num = 2` (Toeplitz), and `method_num = 3` (Banded). This script generates all reported mean and standard error summaries corresponding to this simulation scenario, as presented in the manuscript.


