function val = modified_rastrigin(C)
% MODIFIED_RASTRIGIN applies the Rastrigin function to non-diagonal elements of a PD matrix.
% Input:
%   C : N x N positive definite matrix
% Output:
%   val : scalar Rastrigin function value of non-diagonal elements

    if size(C,1) ~= size(C,2)
        error('Input matrix must be square.');
    end

    % Extract non-diagonal elements
    N = size(C,1);
    non_diag_mask = ~eye(N);
    x = 10 * C(non_diag_mask);  % Scale for variability like in benchmark tests
    d = numel(x);

    % Rastrigin function components
    A = 10;
    val = A * d + sum(x.^2 - A * cos(2 * pi * x));
end