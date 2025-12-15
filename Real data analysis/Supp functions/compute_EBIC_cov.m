function ebic = compute_EBIC_cov(Sigma_est, X, threshold, gamma)
% COMPUTE_EBIC_COV computes the EBIC score for a given estimated covariance matrix
% under the multivariate normal model.
%
% Inputs:
%   Sigma_est : Estimated covariance matrix (p x p)
%   X         : Data matrix (n x p), rows are observations
%   gamma     : EBIC sparsity penalty parameter (default: 0.5)
%
% Output:
%   ebic : EBIC score

    if nargin < 3
        threshold = 1e-03;
    end
    if nargin < 4
        gamma = 0.5;
    end

    [n, p] = size(X);

    % Sample covariance from data
    S = cov(X, 1);  % Use 1/n normalization

    % Compute log-likelihood under multivariate normal
    % ell = -n/2 * (log(det(Sigma_est)) + trace(S * inv(Sigma_est)))
    % or more stable via Cholesky
    try
        R = chol(Sigma_est);
        log_det = 2 * sum(log(diag(R)));
        Sigma_inv = inv(Sigma_est);
        log_likelihood = -n/2 * (log_det + trace(Sigma_inv * S));
    catch
        ebic = Inf;  % Not PD or ill-conditioned
        return;
    end

    % Count number of non-zero off-diagonal entries
    off_diag_mask = ~eye(p);
    num_edges = nnz(abs(Sigma_est(off_diag_mask)) > threshold) / 2;  % symmetric count

    % EBIC formula
    ebic = -2 * log_likelihood + num_edges * log(n) + 4 * gamma * log(nchoosek(p, 2));
end
