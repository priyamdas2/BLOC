function ll = log_likelihood_gaussian(X, Sigma)
% Computes the log-likelihood of multivariate normal observations X
% assuming zero mean and covariance matrix Sigma

    [n, d] = size(X);

    % Validate Sigma
    if size(Sigma,1) ~= d || size(Sigma,2) ~= d
        error('Dimension mismatch: Sigma must be a d x d matrix where d = size(X,2)');
    end
    if norm(Sigma - Sigma') > 10^-10
        error('Sigma must be symmetric');
    end

    % Cholesky decomposition
    [U, p] = chol(Sigma);
    if p > 0
        error('Sigma must be positive definite for log-likelihood computation');
    end

    % Log determinant using Cholesky
    log_det_Sigma = 2 * sum(log(diag(U)));

    % Quadratic term using backsolve (more stable than inv)
    Y = X / U;
    quad_form = sum(Y.^2, 2);  % length-n vector

    % Log-likelihood
    ll = -0.5 * n * d * log(2*pi) ...
         - 0.5 * n * log_det_Sigma ...
         - 0.5 * sum(quad_form);
end