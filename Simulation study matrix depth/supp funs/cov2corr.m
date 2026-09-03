function X_corr = cov2corr(X)
% cov2corr - Convert a covariance (or general square) matrix to correlation matrix
%
% Input:
%   X       : square variance/covariance matrix (p x p)
%
% Output:
%   X_corr  : correlation matrix (p x p)

    % Ensure symmetry
    X = (X + X')/2;

    % Extract standard deviations
    d = sqrt(diag(X));

    % Avoid division by zero
    D_inv = diag(1 ./ d);

    % Compute correlation matrix
    X_corr = D_inv * X * D_inv;
end
