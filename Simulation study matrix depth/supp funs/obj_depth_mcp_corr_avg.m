function f = obj_depth_mcp_corr_avg(C, proj2, U, beta, lambda, gamma)
% ================================================================
% Average directional-depth objective + MCP penalty
%
% Objective:
%
%   f(C) = - D_avg(C)
%          + (1/q) sum_{i<j} p_lambda(|C_ij|)
%
% where D_avg is the average empirical directional depth.
%
% C      : p x p candidate correlation matrix
% proj2  : n x K matrix, (X*U).^2
% U      : p x K fixed unit directions
% beta   : [Phi^{-1}(0.75)]^2
% lambda : MCP tuning parameter
% gamma  : MCP shape parameter
% ================================================================

p = size(C,1);

% ------------------------------------------------
% Average directional depth
% ------------------------------------------------

CU = C * U;

% u_k' C u_k for every direction
directional_scale = sum(U .* CU, 1);

threshold = beta * directional_scale;

prop_below = mean(proj2 <= threshold, 1);
prop_above = mean(proj2 >= threshold, 1);

direction_depth = min(prop_below, prop_above);

% IMPORTANT CHANGE:
% average rather than minimum
D_avg = mean(direction_depth);


% ------------------------------------------------
% MCP penalty on off-diagonal entries
% ------------------------------------------------

idx = triu(true(p),1);

theta = abs(C(idx));

P = zeros(size(theta));

% MCP:
%
% p_lambda(theta)
% = lambda*theta - theta^2/(2*gamma),
%       theta <= gamma*lambda
%
% = gamma*lambda^2/2,
%       theta > gamma*lambda

I1 = (theta <= gamma*lambda);

P(I1) = ...
    lambda .* theta(I1) ...
    - theta(I1).^2 ./ (2*gamma);

I2 = (theta > gamma*lambda);

P(I2) = ...
    gamma * lambda^2 / 2;


% ------------------------------------------------
% Normalize penalty
% ------------------------------------------------

num_offdiag = p*(p-1)/2;

penalty_term = sum(P) / num_offdiag;


% ------------------------------------------------
% Final objective
% ------------------------------------------------

f = -D_avg + penalty_term;

end