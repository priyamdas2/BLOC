function L = Corr2L(C)
    % C must be a positive definite correlation matrix of size M x M
    L = chol(C, 'lower');
end