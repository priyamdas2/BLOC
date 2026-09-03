function ll = log_likelihood_gaussian_SCAD(X, Sigma, lambda, a, P)
% Computes penalized log-likelihood of multivariate Gaussian data X
% using covariance matrix Sigma, with SCAD penalty weighted by mask matrix P.
%
% Inputs:
%   X      : n x d data matrix (assumed zero-mean)
%   Sigma  : d x d positive definite covariance matrix
%   lambda : scalar penalty coefficient
%   a      : SCAD tuning parameter (>2)
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
    if a <= 2
        error('SCAD parameter a must be greater than 2');
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

    % SCAD penalty
    penalized_elements = Sigma(P == 1);
    abs_vals = abs(penalized_elements);
    scad_penalty = sum(arrayfun(@(t) SCAD_scalar(t, lambda, a), abs_vals));

    % Penalized log-likelihood
    ll = ll - scad_penalty;
end

function val = SCAD_scalar(theta, lambda, a)
% SCAD penalty for a single scalar theta
    t = abs(theta);
    if t <= lambda
        val = lambda * t;
    elseif t <= a * lambda
        val = (-t^2 + 2 * a * lambda * t - lambda^2) / (2 * (a - 1));
    else
        val = 0.5 * (a + 1) * lambda^2;
    end
end