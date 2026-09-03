function f = obj_depth_scad_corr_avg(C, proj2, U, beta, lambda, a)
% ================================================================
% Average directional-depth objective + SCAD penalty
%
% Objective:
%
%   f(C) = - D_avg(C)
%          + (1/q) sum_{i<j} p_lambda(|C_ij|)
%
% where
%
%   D_avg(C) = (1/K) sum_k d_k(C)
%
% and
%
%   d_k(C) =
%     min{ P_n[(u_k'X)^2 <= beta*u_k'C*u_k],
%          P_n[(u_k'X)^2 >= beta*u_k'C*u_k] }.
%
% C      : p x p candidate correlation matrix
% proj2  : n x K matrix, (X*U).^2
% U      : p x K fixed unit directions
% beta   : [Phi^{-1}(0.75)]^2
% lambda : SCAD tuning parameter
% a      : SCAD shape parameter, typically 3.7
% ================================================================

p = size(C,1);

% ------------------------------------------------
% Average directional depth
% ------------------------------------------------

CU = C * U;

% u_k' C u_k for every direction
directional_scale = sum(U .* CU, 1);

% Gaussian median threshold
threshold = beta * directional_scale;

% Empirical directional probabilities
prop_below = mean(proj2 <= threshold, 1);
prop_above = mean(proj2 >= threshold, 1);

% Depth separately for every direction
direction_depth = min(prop_below, prop_above);

% IMPORTANT CHANGE:
% average instead of minimum across directions
D_avg = mean(direction_depth);


% ------------------------------------------------
% SCAD penalty on off-diagonal entries
% ------------------------------------------------

idx = triu(true(p),1);

theta = abs(C(idx));

P = zeros(size(theta));

% Region 1: 0 <= |theta| <= lambda
I1 = (theta <= lambda);

P(I1) = lambda .* theta(I1);

% Region 2: lambda < |theta| <= a*lambda
I2 = (theta > lambda) & ...
     (theta <= a*lambda);

P(I2) = ...
    (-theta(I2).^2 ...
     + 2*a*lambda.*theta(I2) ...
     - lambda^2) ...
    ./ (2*(a-1));

% Region 3: |theta| > a*lambda
I3 = (theta > a*lambda);

P(I3) = ...
    ((a+1)*lambda^2)/2;


% ------------------------------------------------
% Normalize penalty by number of off-diagonals
% ------------------------------------------------

num_offdiag = p*(p-1)/2;

penalty_term = sum(P) / num_offdiag;


% ------------------------------------------------
% Final objective
% ------------------------------------------------

f = -D_avg + penalty_term;

end