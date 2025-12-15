function ll_ENet = log_likelihood_gaussian_ElasticNet(X, Sigma, lambda, alpha, P)
% Computes penalized log-likelihood of multivariate Gaussian data X
% using covariance matrix Sigma, with Elastic Net penalty weighted by mask matrix P.
%
% Inputs:
%   X      : n x d data matrix (assumed zero-mean)
%   Sigma  : d x d positive definite covariance matrix
%   lambda : scalar penalty coefficient
%   alpha  : mixing parameter (0 = Ridge, 1 = Lasso)
%   P      : d x d binary mask (1 = penalize, 0 = ignore)

    [n, d] = size(X);

    % Validate inputs
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
    if alpha < 0 || alpha > 1
        error('alpha must be in [0, 1]');
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

    % Elastic Net penalty (L1 + L2)
    penalized_elements = Sigma(P == 1);
    l1_term = alpha * sum(abs(penalized_elements));
    l2_term = (1 - alpha) / 2 * sum(penalized_elements .^ 2);

    penalty = lambda * (l1_term + l2_term);

    % Penalized log-likelihood
    ll_ENet = ll - penalty;
end
