function ll = log_likelihood_gaussian_RIDGE(X, Sigma, lambda, P)
% Computes penalized log-likelihood of multivariate Gaussian data X
% using covariance matrix Sigma, with L2 penalty weighted by mask matrix P.
%
% Inputs:
%   X      : n x d data matrix (assumed zero-mean)
%   Sigma  : d x d positive definite covariance matrix
%   lambda : scalar penalty coefficient
%   P      : d x d binary mask (1 = penalize, 0 = ignore)

    [n, d] = size(X);

    % Validate Sigma and P
    if size(Sigma,1) ~= d || size(Sigma,2) ~= d
        error('Sigma must be a d x d matrix where d = size(X,2)');
    end
    if ~isequal(size(P), [d, d])
        error('Mask matrix P must be of size d x d');
    end
    if norm(Sigma - Sigma', 'fro') > 1e-10
        error('Sigma must be symmetric');
    end
    if ~all(ismember(P(:), [0 1]))
        error('P must be a binary matrix (only 0s and 1s)');
    end
    if min(eig(Sigma)) < 1e-4
        Sigma = make_pd(Sigma);
    end

    % Cholesky decomposition
    [U, p_chol] = chol(Sigma);
    if p_chol > 0
        disp(min(eig(Sigma)));
        error('Sigma must be positive definite');
    end

    % Log determinant using Cholesky
    log_det_Sigma = 2 * sum(log(diag(U)));

    % Quadratic form using Cholesky
    Y = X / U;
    quad_form = sum(Y.^2, 2);

    % Base log-likelihood
    ll = - log_det_Sigma - mean(quad_form);

    % L2 (Ridge) penalty
    l2_penalty = lambda * sum(Sigma(P == 1).^2);

    % Penalized log-likelihood
    ll = ll - l2_penalty;
end
