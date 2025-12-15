function val = modified_ackley(C)
% MODIFIED_ACKLEY_CORR applies Ackley function to non-diagonal elements of a correlation matrix.
% Input:
%   C : N x N pd matrix 
% Output:
%   val : scalar Ackley function value of non-diagonal elements

    if size(C,1) ~= size(C,2)
        error('Input matrix must be square.');
    end

    % Extract non-diagonal elements
    N = size(C,1);
    non_diag_mask = ~eye(N);
    x = 10*C(non_diag_mask);  % Column vector of length N^2 - N

    % Ackley parameters
    a = 20;
    b = 0.2;
    c = 2 * pi;
    d = numel(x);

    % Ackley function
    sum_sq = sum(x.^2);
    sum_cos = sum(cos(c * x));
    term1 = -a * exp(-b * sqrt(sum_sq / d));
    term2 = -exp(sum_cos / d);
    val = term1 + term2 + a + exp(1);
end