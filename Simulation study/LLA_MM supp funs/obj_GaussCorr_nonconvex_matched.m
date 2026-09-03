function f = ...
    obj_GaussCorr_nonconvex_matched( ...
    X,C,lambda,penalty,shape,P)
% ================================================================
% EXACT matched negative penalized log-likelihood corresponding
% to the BLOC Gaussian SCAD/MCP functions.
%
% Your BLOC function computes
%
%       ll =
%       -log(det(C))
%       -mean_i x_i'C^{-1}x_i
%       -penalty.
%
% This function returns
%
%       f = -ll,
%
% so SMALLER is better:
%
%       f =
%       log(det(C))
%       +mean_i x_i'C^{-1}x_i
%       +penalty.
%
% P is the SAME BLOC mask:
%
%       P = ones(p,p)-eye(p).
% ================================================================

[~,p] = size(X);

C = (C+C')/2;


% ================================================================
% Cholesky decomposition
% ================================================================

[U,flag] = chol(C);

if flag ~= 0

    f = Inf;

    return;

end


% ================================================================
% Gaussian negative log-likelihood
% ================================================================

log_det_C = ...
    2*sum(log(diag(U)));

% This is exactly the operation used in your BLOC code
Y = X/U;

quad_form = ...
    sum(Y.^2,2);

f = ...
    log_det_C ...
    + mean(quad_form);


% ================================================================
% EXACT SAME penalized entries as BLOC
% ================================================================

t = ...
    abs(C(P==1));

pen = ...
    zeros(size(t));


switch upper(penalty)

    % ============================================================
    % SCAD
    % ============================================================

    case 'SCAD'

        a = shape;

        if a <= 2
            error('SCAD parameter a must be > 2.');
        end

        I1 = ...
            (t <= lambda);

        I2 = ...
            (t > lambda) & ...
            (t <= a*lambda);

        I3 = ...
            (t > a*lambda);

        pen(I1) = ...
            lambda.*t(I1);

        pen(I2) = ...
            ( ...
            -t(I2).^2 ...
            + 2*a*lambda.*t(I2) ...
            - lambda^2 ...
            ) ...
            ./ ...
            (2*(a-1));

        pen(I3) = ...
            0.5*(a+1)*lambda^2;


    % ============================================================
    % MCP / MC+
    % ============================================================

    case 'MCP'

        gamma = shape;

        if gamma <= 1
            error('MCP gamma must be > 1.');
        end

        I1 = ...
            (t <= gamma*lambda);

        I2 = ...
            (t > gamma*lambda);

        pen(I1) = ...
            lambda.*t(I1) ...
            - t(I1).^2/(2*gamma);

        pen(I2) = ...
            0.5*gamma*lambda^2;


    otherwise

        error( ...
            'penalty must be ''SCAD'' or ''MCP''.');

end


% ================================================================
% Add nonconvex penalty
% ================================================================

f = ...
    f + sum(pen);

end