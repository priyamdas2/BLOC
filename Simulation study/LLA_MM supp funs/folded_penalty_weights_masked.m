function W = ...
    folded_penalty_weights_masked( ...
    C,lambda,penalty,shape,P)
% ================================================================
% LLA derivative weights
%
%       W_ij = P_ij * p'_lambda(|C_ij|)
%
% Uses exactly the same full off-diagonal mask as BLOC.
%
% SCAD:
%
% p'(t) =
%
%   lambda,                    0 <= t <= lambda
%
%   (a*lambda-t)/(a-1),        lambda < t < a*lambda
%
%   0,                         t >= a*lambda
%
%
% MCP:
%
% p'(t) =
%
%   lambda - t/gamma,          t <= gamma*lambda
%
%   0,                         t > gamma*lambda
% ================================================================

T = abs(C);

W = zeros(size(C));


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
            (T <= lambda);

        I2 = ...
            (T > lambda) & ...
            (T < a*lambda);

        W(I1) = lambda;

        W(I2) = ...
            (a*lambda-T(I2)) ...
            ./ ...
            (a-1);

        % Entries with T >= a*lambda have weight zero


    % ============================================================
    % MCP / MC+
    % ============================================================

    case 'MCP'

        gamma = shape;

        if gamma <= 1
            error('MCP gamma must be > 1.');
        end

        W = ...
            max( ...
            lambda-T/gamma, ...
            0);


    otherwise

        error( ...
            'penalty must be ''SCAD'' or ''MCP''.');

end


% ================================================================
% EXACT same penalty mask used by BLOC
% ================================================================

W = W .* P;

end