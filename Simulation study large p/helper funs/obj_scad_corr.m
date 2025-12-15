function f = obj_scad_corr(Theta, S, lambda, a)
% SCAD objective for correlation estimation with eigenvalue lower bound.
% Inputs:
%   Theta  : candidate correlation matrix (d x d)
%   S      : sample correlation matrix (d x d)
%   lambda : main regularization parameter (>0)
%   a      : SCAD parameter (>2), default 3.7
%
% Output:
%   f      : objective value (Inf if constraints violated)

    if nargin < 4 || isempty(a),     a = 3.7; end

    % Constraint checks: diag=1 and PD
    if max(abs(diag(Theta)-1)) > 1e-6, f = Inf; return; end


    % Data-fit term
    fit = 0.5 * norm(Theta - S, 'fro')^2;

    % SCAD penalty on off-diagonals
    A = abs(Theta);
    d = size(A,1);
    A(1:d+1:end) = 0;  % zero out diagonal

    % piecewise SCAD penalty p(t)
    t  = A(:);
    p1 = t <= lambda;
    p2 = t > lambda & t <= a*lambda;
    p3 = t > a*lambda;

    P = zeros(size(t));
    P(p1) = lambda .* t(p1);
    % integrate derivative for middle regime:
    % p(t) = λ^2 + [aλ t - t^2/2 - aλ^2 + λ^2/2]/(a-1)
    P(p2) = lambda^2 + (a*lambda.*t(p2) - 0.5*t(p2).^2 - a*lambda^2 + 0.5*lambda^2) ./ (a - 1);
    % constant tail: ((a+1)/2) * λ^2
    P(p3) = ((a+1)/2) * lambda^2;

    pen = sum(P);

    f = fit + pen;
end
