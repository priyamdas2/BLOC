function val = modified_griewank(C)
% MODIFIED_GRIEWANK_CORR applies Griewank function to non-diagonal elements of a pd matrix.
% Input:
%   C : N x N pd matrix
% Output:
%   val : scalar Griewank function value of non-diagonal elements

    if size(C,1) ~= size(C,2)
        error('Input matrix must be square.');
    end

    % Extract non-diagonal elements
    N = size(C,1);
    non_diag_mask = ~eye(N);
    x = 100*C(non_diag_mask);  % Column vector of length N^2 - N
    d = numel(x);

    % Griewank function components
    sum_term = sum(x.^2) / 4000;
    i = 1:d;
    prod_term = prod(cos(x' ./ sqrt(i)));

    % Final Griewank value
    val = 1 + sum_term - prod_term;
end