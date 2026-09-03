function R = corrMatrixFromVec(x,p)
% ================================================================
% Reconstruct a p x p correlation matrix from its
% upper-triangular off-diagonal elements.
%
% Diagonal is fixed at 1.
% ================================================================

expected_length = ...
    p*(p-1)/2;

if length(x) ~= expected_length

    error( ...
        ['Length of x is %d but expected ' ...
         'p*(p-1)/2 = %d.'], ...
        length(x), ...
        expected_length);

end

R = eye(p);

idx = triu(true(p),1);

R(idx) = x(:);

% Copy upper triangular entries to lower triangle
R = R + triu(R,1)';

end