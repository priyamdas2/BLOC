function x = vecFromCorrMatrix(R)
% ================================================================
% Convert a p x p correlation matrix to a vector containing
% the upper-triangular off-diagonal entries.
% ================================================================

p = size(R,1);

idx = triu(true(p),1);

x = R(idx);

x = x(:);

end