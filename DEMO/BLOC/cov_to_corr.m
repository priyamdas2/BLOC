function R = cov_to_corr(Sigma)
% Converts a covariance matrix to a correlation matrix

    d = sqrt(diag(Sigma));         % standard deviations
    D_inv = diag(1 ./ d);          % inverse of std devs
    R = D_inv * Sigma * D_inv;     % normalize by std devs
end