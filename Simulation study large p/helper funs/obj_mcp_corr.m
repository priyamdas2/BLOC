function f = obj_mcp_corr(Theta, S, lambda, gamma)
% MCP objective for correlation estimation with eigenvalue lower bound.
% Inputs:
%   Theta  : candidate correlation matrix (d x d)
%   S      : sample correlation matrix (d x d)
%   lambda : main regularization parameter (>0)
%   gamma  : MCP parameter (>1), default 3
%
% Output:
%   f      : objective value (Inf if constraints violated)

    if nargin < 4 || isempty(gamma), gamma = 3; end

    % Constraint checks
    if max(abs(diag(Theta)-1)) > 1e-6, f = Inf; return; end


    % Data-fit term
    fit = 0.5 * norm(Theta - S, 'fro')^2;

    % MCP penalty on off-diagonals
    A = abs(Theta);
    d = size(A,1);
    A(1:d+1:end) = 0;

    t = A(:);
    P = zeros(size(t));

    % MCP: p(t) = λ t - t^2/(2γ) for t <= γλ ; else = (γ λ^2)/2
    thresh = gamma*lambda;
    idx1 = t <= thresh;
    idx2 = ~idx1;

    P(idx1) = lambda.*t(idx1) - (t(idx1).^2) / (2*gamma);
    P(idx2) = 0.5 * gamma * lambda^2;

    pen = sum(P);

    f = fit + pen;
end
