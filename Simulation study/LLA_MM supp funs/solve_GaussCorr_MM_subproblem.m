function [Cnew,status,optval] = ...
    solve_GaussCorr_MM_subproblem( ...
    Ct,B,W,eig_floor)
% ================================================================
% One convex LLA/MM subproblem.
%
% At current iterate Ct, minimize
%
%       trace(Ct^{-1} C)
%       + trace(S C^{-1})
%       + sum_ij W_ij |C_ij|
%
% subject to
%
%       diag(C) = 1
%       C >= eig_floor * I
%
% where S = B*B'.
%
%
% Why this is MM:
%
% log(det(C)) is concave, hence
%
% log(det(C))
% <=
% log(det(Ct))
% + trace(Ct^{-1}(C-Ct)).
%
% Ignoring constants leaves trace(Ct^{-1}C).
%
%
% The inverse term is represented using:
%
%       [ C    B  ]
%       [ B'   Y  ] >= 0
%
% which implies
%
%       Y >= B' C^{-1} B.
%
% Minimizing trace(Y) therefore gives
%
%       trace(Y) = trace(S C^{-1}).
%
% Requires CVX.
% ================================================================

p = size(Ct,1);


% ================================================================
% Tangent coefficient for logdet
% ================================================================

A = Ct \ eye(p);

A = (A+A')/2;


% ================================================================
% CVX convex SDP
% ================================================================

cvx_begin sdp quiet

    variable Cnew(p,p) symmetric
    variable Y(p,p) symmetric

    expression penalty_term

    % W already includes the BLOC mask P.
    %
    % Because P = ones(p)-eye(p), this counts BOTH
    % C_ij and C_ji exactly as BLOC does.
    penalty_term = ...
        sum(sum(W .* abs(Cnew)));

    minimize( ...
        trace(A*Cnew) ...
        + trace(Y) ...
        + penalty_term ...
        )

    subject to

        % Correlation-matrix diagonal
        diag(Cnew) == ones(p,1);

        % Positive definiteness
        Cnew - eig_floor*eye(p) >= 0;

        % Schur complement for trace(S*C^{-1})
        [Cnew, B; ...
         B',   Y] >= 0;

cvx_end


status = cvx_status;

optval = cvx_optval;


% ================================================================
% Fail explicitly if CVX did not solve the subproblem
% ================================================================

if ~contains(cvx_status,'Solved')

    error( ...
        'CVX subproblem failed. Status: %s', ...
        cvx_status);

end


% Numerical symmetrization
Cnew = (Cnew+Cnew')/2;

Cnew(1:p+1:end) = 1;

end