%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEMO: BLOC Optimization on Modified Benchmark Functions
%
% This script demonstrates the usage of:
%   (i)  BLOC        : Non-parallel optimization
%   (ii) BLOCparallel: Parallel optimization
%
% Both solvers optimize a modified benchmark function defined
% over the space of valid correlation matrices.
%
% This file is intended as a FINAL DEMO script.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars;
warning('off', 'all');

% Add BLOC source code to MATLAB path
addpath('./BLOC/');

%% ============================================================
%  DEMO 1: BLOC (Non-parallel) on a 10 x 10 Modified Ackley Function
%  ============================================================

% Problem dimension (correlation matrix size)
M = 10;

% Objective function handle
% Input:  M x M correlation matrix
% Output: scalar objective value
objFun = @(C) modified_ackley(C);

% ------------------------------------------------------------
% Generate a random feasible starting point
% ------------------------------------------------------------
% randCorrMatrix ensures:
%   - symmetry
%   - unit diagonal
%   - positive semi-definiteness
rand_seed = 2025;
C0 = randCorrMatrix(M, rand_seed);

% ------------------------------------------------------------
% Set BLOC optimization options
% ------------------------------------------------------------
options.MaxRuns = 5;     % Maximum number of outer runs
options.MaxTime = 20;    % Maximum wall-clock time (seconds)

% ------------------------------------------------------------
% Run BLOC
% ------------------------------------------------------------
[C_BLOC, fval_BLOC, comp_time_BLOC] = BLOC(objFun, C0, options);

% ------------------------------------------------------------
% Print results (professionally formatted)
% ------------------------------------------------------------
fprintf('\n============================================================\n');
fprintf('BLOC (Non-parallel) Optimization Results\n');
fprintf('------------------------------------------------------------\n');
fprintf('Problem dimension           : %d x %d correlation matrix\n', M, M);
fprintf('Objective function value    : %.6f\n', fval_BLOC);
fprintf('Total computation time (s)  : %.2f\n', comp_time_BLOC);
fprintf('------------------------------------------------------------\n');
fprintf('Optimized correlation matrix C_BLOC:\n');
disp(C_BLOC);

fprintf('Objective value re-evaluated at C_BLOC (consistency check): %.6f\n', ...
        objFun(C_BLOC));
fprintf('============================================================\n');


%% ============================================================
%  DEMO 2: BLOCparallel on a 20 x 20 Modified Ackley Function
%  ============================================================

% Problem dimension (larger-scale example)
M = 20;

% Objective function handle
objFun = @(C) modified_ackley(C);

% ------------------------------------------------------------
% Generate a random feasible starting point
% ------------------------------------------------------------
rand_seed = 1;
C0 = randCorrMatrix(M, rand_seed);

% ------------------------------------------------------------
% Set BLOCparallel optimization options
% ------------------------------------------------------------
options.MaxRuns      = 10;        % Maximum number of outer runs
options.phi          = 1e-18;     % Numerical tolerance for step refinement
options.MaxIter      = 5000;      % Maximum iterations per run
options.DisplayEvery = 5;         % Display progress every given time window

% ------------------------------------------------------------
% Run BLOCparallel
% ------------------------------------------------------------
[C_BLOCparallel, fval_BLOCparallel, comp_time_BLOCparallel] = ...
    BLOCparallel(objFun, C0, options);

% ------------------------------------------------------------
% Print results (professionally formatted)
% ------------------------------------------------------------
fprintf('\n============================================================\n');
fprintf('BLOCparallel Optimization Results\n');
fprintf('------------------------------------------------------------\n');
fprintf('Problem dimension           : %d x %d correlation matrix\n', M, M);
fprintf('Objective function value    : %.6f\n', fval_BLOCparallel);
fprintf('Total computation time (s)  : %.2f\n', comp_time_BLOCparallel);
fprintf('------------------------------------------------------------\n');
fprintf('Optimized correlation matrix C_BLOCparallel:\n');
disp(C_BLOCparallel);

fprintf('Objective value re-evaluated at C_BLOCparallel (consistency check): %.6f\n', ...
        objFun(C_BLOCparallel));
fprintf('============================================================\n');
