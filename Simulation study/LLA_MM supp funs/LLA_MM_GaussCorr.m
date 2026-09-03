function [C_final, out] = LLA_MM_GaussCorr( ...
    X, lambda, penalty, shape, P, C0, opts)
% ================================================================
% LLA_MM_GaussCorr
%
% LLA/MM comparator for the SAME penalized Gaussian correlation
% objective used by BLOC:
%
%   min_C
%
%       log(det(C))
%       + mean_i x_i' C^{-1} x_i
%       + sum_{(i,j):P_ij=1} p_lambda(|C_ij|)
%
% subject to
%
%       C = C'
%       diag(C) = 1
%       C >= eig_floor * I
%
% penalty:
%       'SCAD'
%       'MCP'
%
% shape:
%       SCAD : a     (use 3.7)
%       MCP  : gamma (use 3)
%
% P must be EXACTLY the same penalty mask as BLOC:
%
%       P = ones(p,p) - eye(p);
%
% Requires CVX.
%
% Stage 1:
%       L1/MM contraction initialization.
%
% Stage 2:
%       SCAD/MCP LLA tightening.
% ================================================================


% ================================================================
% Defaults
% ================================================================

if nargin < 7 || isempty(opts)
    opts = struct();
end

if ~isfield(opts,'maxOuter')
    opts.maxOuter = 20;
end

if ~isfield(opts,'maxInit')
    opts.maxInit = 10;
end

if ~isfield(opts,'tol')
    opts.tol = 1e-5;
end

% Match the numerical PD threshold used in the BLOC likelihood
if ~isfield(opts,'eig_floor')
    opts.eig_floor = 1e-4;
end

if ~isfield(opts,'doL1Init')
    opts.doL1Init = true;
end

if ~isfield(opts,'verbose')
    opts.verbose = true;
end


% ================================================================
% Input checks
% ================================================================

[n,p] = size(X);

if ~isequal(size(P),[p,p])
    error('P must be p x p.');
end

if ~all(ismember(P(:),[0 1]))
    error('P must contain only zeros and ones.');
end

if norm(P-P','fro') > 1e-12
    error('P must be symmetric.');
end

if nargin < 6 || isempty(C0)
    C0 = eye(p);
end

C0 = (C0+C0')/2;
C0(1:p+1:end) = 1;


% ================================================================
% Sample second moment
%
% BLOC uses
%
%       mean_i x_i' C^{-1} x_i
%
% which equals
%
%       trace((X'X/n) C^{-1}).
% ================================================================

S = (X'*X)/n;

S = (S+S')/2;


% ================================================================
% Factor S = B*B'
%
% Needed for the Schur-complement representation of
% trace(S*C^{-1}).
% ================================================================

[L,flag] = chol(S,'lower');

if flag == 0

    B = L;

else

    % Numerical PSD fallback
    [V,D] = eig(S);

    dd = real(diag(D));

    dd(dd < 0 & dd > -1e-10) = 0;

    if any(dd < 0)
        error('X''*X/n is not positive semidefinite.');
    end

    B = V*diag(sqrt(dd));

end


% ================================================================
% STAGE 1:
% L1/MM contraction initialization
%
% At zero, both SCAD and MCP have derivative lambda.
% Thus this gives the natural weighted-L1 contraction stage.
% ================================================================

C = C0;

if opts.doL1Init

    if opts.verbose

        fprintf('\n');
        fprintf('=============================================\n');
        fprintf('        L1/MM INITIALIZATION\n');
        fprintf('=============================================\n');

    end

    % Same full off-diagonal counting convention as BLOC
    W = lambda*P;

    for it0 = 1:opts.maxInit

        C_old = C;

        [C,cvx_status,~] = ...
            solve_GaussCorr_MM_subproblem( ...
            C_old, ...
            B, ...
            W, ...
            opts.eig_floor);

        relchg = ...
            norm(C-C_old,'fro') / ...
            max(1,norm(C_old,'fro'));

        if opts.verbose

            fprintf( ...
                'L1 iter %2d | rel.change = %.4e | CVX = %s\n', ...
                it0,relchg,cvx_status);

        end

        if relchg < opts.tol
            break;
        end

    end

end


% ================================================================
% STAGE 2:
% SCAD/MCP LLA tightening
% ================================================================

obj_history = nan(opts.maxOuter+1,1);

obj_history(1) = ...
    obj_GaussCorr_nonconvex_matched( ...
    X, ...
    C, ...
    lambda, ...
    penalty, ...
    shape, ...
    P);


if opts.verbose

    fprintf('\n');
    fprintf('=============================================\n');
    fprintf('        %s LLA/MM TIGHTENING\n',upper(penalty));
    fprintf('=============================================\n');

    fprintf( ...
        'Initial true objective = %.12f\n', ...
        obj_history(1));

end


for it = 1:opts.maxOuter

    C_old = C;

    % ============================================================
    % LLA weights:
    %
    %       w_ij = p'_lambda(|C_ij|)
    %
    % using EXACTLY the BLOC mask P.
    % ============================================================

    W = ...
        folded_penalty_weights_masked( ...
        C_old, ...
        lambda, ...
        penalty, ...
        shape, ...
        P);


    % ============================================================
    % Solve convex MM/LLA surrogate
    % ============================================================

    [C,cvx_status,subproblem_value] = ...
        solve_GaussCorr_MM_subproblem( ...
        C_old, ...
        B, ...
        W, ...
        opts.eig_floor);


    % ============================================================
    % Evaluate ORIGINAL nonconvex BLOC objective
    % ============================================================

    obj_history(it+1) = ...
        obj_GaussCorr_nonconvex_matched( ...
        X, ...
        C, ...
        lambda, ...
        penalty, ...
        shape, ...
        P);


    relchg = ...
        norm(C-C_old,'fro') / ...
        max(1,norm(C_old,'fro'));

    objchg = ...
        abs(obj_history(it+1)-obj_history(it)) / ...
        max(1,abs(obj_history(it)));


    if opts.verbose

        fprintf( ...
            ['Iter %2d | Obj = %.12f | ' ...
             'rel.C = %.4e | rel.obj = %.4e | CVX = %s\n'], ...
             it, ...
             obj_history(it+1), ...
             relchg, ...
             objchg, ...
             cvx_status);

    end


    % ============================================================
    % MM diagnostic
    %
    % Up to numerical CVX error, original objective should
    % be non-increasing.
    % ============================================================

    if obj_history(it+1) > obj_history(it) + 1e-6

        warning( ...
            ['Original objective increased from %.12f ' ...
             'to %.12f. Check CVX accuracy.'], ...
             obj_history(it), ...
             obj_history(it+1));

    end


    % Require BOTH matrix and objective convergence
    if max(relchg,objchg) < opts.tol
        break;
    end

end


% ================================================================
% Final cleanup
% ================================================================

C_final = (C+C')/2;

C_final(1:p+1:end) = 1;


% ================================================================
% Output information
% ================================================================

out.iter = it;

out.obj_history = ...
    obj_history(1:it+1);

out.final_objective = ...
    obj_history(it+1);

out.lambda = lambda;

out.penalty = penalty;

out.shape = shape;

out.subproblem_value = ...
    subproblem_value;

out.min_eig = ...
    min(real(eig(C_final)));

end