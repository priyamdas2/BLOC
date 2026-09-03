function X = WeiZhaoMCP(S, tau, phi0, gamma_bt, lambda_val, a_val)
% WeiZhaoMCP
%
% Correlation-based Wei-Zhao MCP estimator using the MM +
% proximal-gradient/backtracking algorithm.
%
% This is a MATLAB translation of the supplied Python PDTE_FC algorithm,
% with the correlation-matrix modification from Wei & Zhao (2023):
% after each proximal step the diagonal is set to 1.
%
% Inputs
%   S          : sample correlation matrix
%   tau        : log-determinant barrier parameter
%   phi0       : initial quadratic coefficient / step-size denominator
%   gamma_bt   : backtracking multiplier (>1)
%   lambda_val : MCP regularization parameter
%   a_val      : MCP shape parameter (Wei-Zhao simulations use a = 2)
%
% Output
%   X          : estimated positive-definite correlation matrix
%
% Algorithmic constants below match the supplied Python code:
%   maximum MM stages = 100
%   relative inner stopping tolerance = 1e-5
%   relative outer stopping tolerance = 1e-5

max_outer = 100;
tol_inner = 1e-5;
tol_outer = 1e-5;

d = size(S,1);

% Python class receives x externally. For the correlation estimator,
% use the identity initialization, consistent with the paper.
x = eye(d);

value_list = [];

for outer = 1:max_outer

    % The supplied Python code resets the inner iteration counter
    % at every MM stage.
    iteration = 0;

    x_outer_old = x;

    % MM weights: derivative of MCP at the previous outer iterate.
    w = mcp_derivative(x_outer_old, lambda_val, a_val);

    % Do not penalize diagonal entries.
    w(1:d+1:end) = 0;

    while true

        iteration = iteration + 1;
        x_old = x;

        if iteration == 1
            phi_t = phi0;
        else
            phi_t = max(phi0, phi_t/gamma_bt);
        end

        while true

            % Gradient of
            %   0.5*||X-S||_F^2 - tau*log(det(X))
            grad_old = gradient_function(x_old, S, tau);

            % Proximal-gradient step.
            x_trial = x_old - (1/phi_t)*grad_old;

            % Weighted elementwise soft thresholding.
            x_trial = softthresholding(x_trial, w/phi_t);

            % Correlation-based modification in Wei-Zhao:
            % enforce unit diagonal in the proximal step.
            x_trial(1:d+1:end) = 1;

            % Python evaluates the smooth function and then checks PD.
            % In MATLAB, checking PD first avoids complex log(det(.))
            % for a non-PD trial point. This does not change whether
            % such a trial is accepted.
            % np.linalg.eigh() in the Python code assumes a symmetric
            % matrix. Use the symmetric part only for this PD check; do
            % not alter the iterate itself.
            eig_trial = eig((x_trial + x_trial')/2);
            is_pd = min(real(eig_trial)) > 0 && ...
                    max(abs(imag(eig_trial))) < 1e-10;

            if is_pd

                a1 = f_function(x_old, S, tau);
                a2 = sum(grad_old .* (x_trial - x_old), 'all');
                a3 = norm(x_trial - x_old, 'fro')^2;

                g_value = a1 + a2 + 0.5*phi_t*a3;
                f_value = f_function(x_trial, S, tau);

                if g_value >= f_value
                    x_new = x_trial;
                    break;
                end
            end

            phi_t = gamma_bt*phi_t;
        end

        % Same objective-value tracking as the supplied Python code.
        value_list(end+1,1) = function_value( ...
            x_new, S, tau, lambda_val, a_val); %#ok<AGROW>

        x = x_new;

        % NumPy's matrix norm defaults to Frobenius norm.
        rel_inner = norm(x_new - x_old, 'fro') / norm(x_old, 'fro');

        if rel_inner < tol_inner
            break;
        end
    end

    rel_outer = norm(x_new - x_outer_old, 'fro') / ...
                norm(x_outer_old, 'fro');

    if rel_outer < tol_outer
        fprintf('=> WeiZhaoMCP converged.\n');
        X = x_new;
        return;
    end
end

X = x_new;

end


% ========================================================================
% Smooth objective:
% 0.5 ||X-S||_F^2 - tau log det(X)
% ========================================================================
function val = f_function(X, S, tau)

quad = 0.5*norm(X - S, 'fro')^2;

if tau == 0
    % Mathematically identical, while avoiding 0*log(det(X)) numerical
    % problems when det(X) underflows in moderate/high dimension.
    val = quad;
else
    % Stable log-determinant for a positive-definite matrix.
    R = chol((X + X')/2);
    logdetX = 2*sum(log(diag(R)));
    val = quad - tau*logdetX;
end

end


% ========================================================================
% Gradient of smooth objective
% ========================================================================
function G = gradient_function(X, S, tau)

if tau == 0
    G = X - S;
else
    invX = inv(X);
    invX = 0.5*(invX + invX');
    G = X - S - tau*invX;
end

end


% ========================================================================
% Weighted elementwise soft thresholding
% ========================================================================
function Xout = softthresholding(B, Lam)

Xout = sign(B) .* max(abs(B) - Lam, 0);

end


% ========================================================================
% MCP derivative:
% p'_lambda(|x|) = (lambda - |x|/a) I{|x| <= a lambda}
% ========================================================================
function W = mcp_derivative(X, lambda_val, a_val)

AX = abs(X);
W = (lambda_val - AX/a_val) .* ...
    (AX <= lambda_val*a_val);

end


% ========================================================================
% MCP penalty
% ========================================================================
function P = mcp_penalty(X, lambda_val, a_val)

AX = abs(X);

is_linear = (AX <= lambda_val*a_val);
is_constant = (a_val*lambda_val < AX);

linear_part = ...
    (lambda_val*AX - (X.^2)/(2*a_val)) .* is_linear;

constant_part = ...
    (0.5*(lambda_val^2)*a_val) .* is_constant;

P = linear_part + constant_part;

end


% ========================================================================
% Full nonconvex objective used only for value tracking, as in Python
% ========================================================================
function val = function_value(X, S, tau, lambda_val, a_val)

P = mcp_penalty(X, lambda_val, a_val);

quad = 0.5*norm(X - S, 'fro')^2;

if tau == 0
    smooth_val = quad;
else
    R = chol((X + X')/2);
    logdetX = 2*sum(log(diag(R)));
    smooth_val = quad - tau*logdetX;
end

val = smooth_val - trace(P) + sum(P, 'all');

end
