function C_pd = make_pd(C)
% Ensures a symmetric matrix C is positive definite by eigenvalue adjustment
    eps_val = 10^-3;
    [V, D] = eig((C + C') / 2);  % ensure symmetry
    D = diag(D);
    D(D < eps_val) = eps_val;   % adjust small/negative eigenvalues
    C_pd = V * diag(D) * V';

    % Normalize to get back a correlation matrix (1s on diagonal)
    d = sqrt(diag(C_pd));
    C_pd = C_pd ./ (d * d');  % re-standardize
end