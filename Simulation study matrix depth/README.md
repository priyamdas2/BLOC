## Simulation Study (Matrix depth Scenario)

This folder contains scripts to reproduce the simulation results under the matrix depth problem. The simulations are organized by method and followed by a summary step.

### Method: BLOC
Run `simulation_BLOC.m` with covariance structures specified by `method_num = 1` (Block matrix) and `method_num = 3` (Banded), and penalty choices given by `penalty = 1` (SCAD) and `penalty = 2` (MCP).

### Other Methods
Run `simulation_other_methods.m` using covariance structures defined by `method_num = 1` (Block matrix) and `method_num = 3` (Banded), and penalty choices given by `penalty = 1` (SCAD) and `penalty = 2` (MCP).

### Results Aggregation
Finally, run `Summary_ALL.m` for covariance structures `method_num = 1` (Block matrix) and `method_num = 3` (Banded). This script generates all reported mean and standard error summaries corresponding to this simulation scenario, as presented in the manuscript.



