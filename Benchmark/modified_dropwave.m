function val = modified_dropwave(C)
% MODIFIED_DROPWAVE_CORR applies Drop-Wave function to non-diagonal elements of a pd matrix.
% Input:
%   C : N x N pd matrix 
% Output:
%   val : scalar Drop-Wave function value of non-diagonal elements

    if size(C,1) ~= size(C,2)
        error('Input matrix must be square.');
    end

    % Extract and scale non-diagonal elements
    N = size(C,1);
    non_diag_mask = ~eye(N);
    x = 5 * C(non_diag_mask);  % scale like original Drop-Wave domain [-5,5]

    r2 = sum(x.^2);
    numerator = 1 + cos(12 * sqrt(r2));
    denominator = 0.5 * r2 + 2;

    val = - numerator / denominator;
end