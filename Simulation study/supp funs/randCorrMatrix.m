function C = randCorrMatrix(M, rand_seed)
    % Generate a random matrix
    rng(rand_seed)
    A = randn(M, M);

    % Make it symmetric positive definite
    Sigma = A' * A;

    % Convert covariance matrix to correlation matrix
    D = sqrt(diag(Sigma));
    C = Sigma ./ (D * D');
end